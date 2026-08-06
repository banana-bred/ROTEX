! ================================================================================================================================ !
program rotex
  !! The main program

  use rotex__kinds,      only: dp
  use rotex__globals,    only: G, read_namelists
  use rotex__drivers,    only: make_grid, diagonalize_rotational_hamiltonian &
                             , make_output_directories, do_coulomb_born_approx, do_kmat_xs, combine_cb_smat_xs
  use rotex__types,      only: eigenh_type, n_states_type, asymtop_rot_channel_l_type &
                             , asymtop_rot_channel_l_vector_type, cmatrix_type, rvector_type, asymtop_rot_transition_type
  use rotex__system,     only: mkdir, die, stdout
  use rotex__hamilton,   only: h_asym, assign_projections, rotate_eigvecs, resolve_c2prime_phi, wangify_symtop_eigvecs
  use rotex__characters, only: lower

  implicit none (type, external)

  integer :: jmin, jmax
    !! The min/max values of J = N + l to consider
  integer :: nmin_ft, nmax_ft
    !! The min/max values of N to consider for the frame transformation.
  integer :: i_n, num_n
  integer, allocatable :: N_values(:)
    !! Contains the values of N to consider

  real(dp) :: start, finish, Eground, c2phi
  real(dp), allocatable :: egrid_tot_smat(:)
    !! The total energy grid for the S-matrix calculations

  character(:), allocatable :: pcb_output_directory, tcb_output_directory, smat_output_directory

  type(n_states_type), allocatable :: n_states(:)
    !! The states of the system for a given N

  type(asymtop_rot_transition_type), allocatable :: transitions_cb(:), transitions_smat(:)

  type(rvector_type), allocatable :: Eel_grid_cb(:)
    !! Array of arrays of electron/collision energies for each transition
  type(rvector_type), allocatable :: xs_xcite_tcb(:), xs_xcite_pcb(:)
    !! Array of arrays of the Total and Partial Coulomb-Born cross sections for each transition
  type(rvector_type), allocatable :: xs_xcite_smat(:), xs_dxcite_smat(:)

  ! GO
  start = time()

  call print_header()
  call read_namelists()
  ! call validate_user_input(G)
  call make_output_directories(pcb_output_directory, tcb_output_directory, smat_output_directory)

  ! -- determine the (number of) N values
  if(G%USE_KMAT .eqv. .true.) then
    ! -- if we  use the K-matrix, then the S-matrix will be block
    !    diagonal in J = [ max(0, Nmin - lmax_kmat), Nmax + lmax_kmat ], but for
    !    each J-block of S to be unitary, we must calculate extra rotational wavefunctions
    !    N = [ max(0, Jmin - lmax_kmat), Jmax + lmax  ]
    jmin = max(0, G%NMIN - G%LMAX_KMAT)
    jmax = G%NMAX + G%LMAX_KMAT
    nmin_ft = max(0, jmin - G%LMAX_KMAT)
    nmax_ft = jmax + G%LMAX_KMAT
    num_n    = nmax_ft - nmin_ft + 1
    n_values = [(nmin_ft + i_n, i_n=0, nmax_ft - nmin_ft, 1)]
    write(stdout, '(a)') "using a K-matrix to get cross sections, so we need to calculate more states &
      & so that each J-block of the S-matrix will be unitary."
    write(stdout, '(a, i0, "/", i0)') "smallest/largest value for n considered: ", nmin_ft, nmax_ft
    write(stdout, *)
  else
    ! -- when not using a K-matrix (i.e. just getting Coulomb-Born
    !    cross sections) we only need the N values that we want for
    !    (de-)excitation
    num_n    = G%NMAX - G%NMIN + 1
    n_values = [(G%NMIN + i_n, i_n=0, G%NMAX - G%NMIN, 1)]
    write(stdout, '(a, i0, "/", i0)') "not using a k-matrix to get cross sections, so we only need the user-supplied&
      & values for nmin/nmax: "&
      , G%NMIN, G%NMAX
    write(stdout, *)
  endif

  ! -- this array holds the information on the rotational states of the target
  allocate(N_states(num_N))

  select case(G%ROTOR_KIND)
  case("l","L")
    call print_dipole(G%DIPOLE_ABC(1))
  case("a","A","s","S")
    call print_dipoles
  case default
    call die("Unacceptable G%ROTOR_KIND: "//G%ROTOR_KIND)
  end select

  ! -- diagonalize the hamiltonian and assign state labels for the rotational states that
  !    will be involved in the transitions/collisions
  call diagonalize_rotational_hamiltonian(num_N, N_values, N_states)

  ! -- SYMTOPS ONLY: rotate exactly degenerate ±K pairs into a Wang basis and tag
  !    each state with its C2' character. If USE_KMAT is false, φ is set to 0.0
  select case(G%ROTOR_KIND)
  case("s", "S")
    if(lower(G%RR_DIAG_AXIS) .ne. lower(G%SYMAXIS)) then
      call die("Wangify requires eigvecs that are already in the G%SYMAXIS frame.&
        & G%RR_DIAG_AXIS does not agree ("//G%RR_DIAG_AXIS//"). Not carrying out&
        & the RFT on this.")
    endif
    c2phi = resolve_c2prime_phi(G%USE_KMAT, G%C2PRIME_PHI_DEG)
    write(stdout, '("Wang basis: rotating degenerate ±K pairs, C2prime azimuthal φ = ", F8.5, " rad")') c2phi
    do i_n = 1, num_n
      call wangify_symtop_eigvecs(               &
          N       = N_states(i_n)%N              &
        , eigvecs = N_states(i_n)%eigenH%eigvecs &
        , irchar  = N_states(i_n)%rchar          &
        , phi     = c2phi                        &
      )
      ! -- check rchar. The N=0 rotational function is constant, so it's
      !    invariant under any rotation, so its C₂' character must be +1 regardless of
      !    which C₂' axis (any φ) was chosen. Failure implies that the wangify routine
      !    is not consistent and that nothing that follows is coherent
      if(N_states(i_n)%N .ne. 0) cycle
      if(N_states(i_n)%rchar(1) .ne. +1) call die("rchar convention check failed: N=0 must have rchar=+1")
    enddo
  end select

  if(G%use_CB .eqv. .true.) then
    call do_coulomb_born_approx( &
        n_states                 &
      , Eel_grid_cb              &
      , transitions_cb           &
      , xs_xcite_pcb             &
      , xs_xcite_tcb             &
      , pcb_output_directory     &
      , tcb_output_directory)
  endif

  ! -- do this only AFTER we have called DO_COULOMB_BORN_APPROX because it may use
  !    CDMS energies which will change the energies of our rotational levels (but not the eigenvectors)
  if(G%PRINT_ROT_STATES) call print_rot_targ_states(n_states)

  usingkmat: if(G%use_kmat) then

    write(stdout, *)
    write(stdout, '(A)') "--------------------------"
    write(stdout, '(A)') "Switching to K-matrix data"
    write(stdout, '(A)') "--------------------------"
    write(stdout, *)

    ! -- ensure that the rigid rotor eigenvectors are in the SYMAXIS frame
    if(lower(G%RR_DIAG_AXIS) .eq. lower(G%SYMAXIS)) then
      write(stdout, '("No eigenvector rotation needed. SYMAXIS AND RR_DIAG_AXIS are the same: ",A)') G%SYMAXIS
    else
      write(stdout, '("Rotating rigid-rotor eigenvectors from RR_DIAG_AXIS to SYMAXIS", A, " --> ", A)') &
        G%RR_DIAG_AXIS, G%SYMAXIS
      do i_n=1,  num_n
        call rotate_eigvecs(                         &
            N         = N_states(i_n)%N              &
          , from_axis = G%RR_DIAG_AXIS               &
          , to_axis   = G%SYMAXIS                    &
          , eigvecs   = N_states(i_n)%eigenH%eigvecs &
        )
      enddo
    endif

    ! -- make the total energy grid for the K/S-matrix cross sections
    eground = n_states(1) % eigenh % eigvals(1)
    call make_grid(egrid_tot_smat, eground)

    call do_kmat_xs(          &
        n_states              &
      , egrid_tot_smat        &
      , smat_output_directory &
      , transitions_smat      &
      , xs_xcite_smat         &
      , xs_dxcite_smat        &
    )

    if(G%USE_CB .eqv. .false.) exit usingkmat

    ! -- now we combine the S-matrix cross sections (which are spin-averaged) with the
    !    Coulom-Born correcion (σTCB - σPCB) which are spin-independent. Each CB σ has its own
    !    energy grid, while the σSmat cross sections were all calculated on the same grid of total
    !    energies. The CB cross sections are resonance free, so we'll have to interpolate these
    !    to match the σSmat energy grid before we can do σTot = σSmat + σTCB - σPCB
    write(stdout, *)
    write(stdout, '("Combining Coulomb-Born and S-matrix cross sections.")')
    write(stdout, '("⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻")')

    call combine_cb_smat_xs( &
        Eel_grid_cb          &
      , egrid_tot_smat       &
      , transitions_cb       &
      , xs_xcite_pcb         &
      , xs_xcite_tcb         &
      , transitions_smat     &
      , xs_xcite_smat        &
      , xs_dxcite_smat       &
    )

  endif usingkmat

  finish = time()
  call print_footer(start, finish)

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_header()
    use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
    use rotex__system, only: stdout, determine_system_properties
#ifdef USE_OPENMP
    use omp_lib,       only: omp_get_num_threads, omp_get_thread_num
#endif
    implicit none
    write(stdout, *)
    write(stdout, *)
    write(stdout, '(A)') "        ___           ___           ___           ___     e⁻    ___"
    write(stdout, '(A)') "       /\  \         /\  \         /\  \         /\  \         |\__\"
    write(stdout, '(A)') "      /::\  \       /::\  \        \:\  \       /::\  \        |:|  |"
    write(stdout, '(A)') "     /:/\:\  \     /:/\:\  \  e⁻    \:\  \     /:/\:\  \       |:|  |        e⁻"
    write(stdout, '(A)') "    /::\~\:\  \   /:/  \:\  \       /::\  \   /::\~\:\  \      |:|__|__"
    write(stdout, '(A)') "   /:/\:\ \:\__\ /:/__/ \:\__\     /:/\:\__\ /:/\:\ \:\__\ ____/::::\__\"
    write(stdout, '(A)') "   \/_|::\/:/  / \:\  \ /:/  /    /:/  \/__/ \:\~\:\ \/__/ \::::/~~/~"
    write(stdout, '(A)') "      |:|::/  /   \:\  /:/  /    /:/  /       \:\ \:\__\    ~~|:|~~|"
    write(stdout, '(A)') "      |:|\/__/     \:\/:/  /    /:/  /         \:\ \/__/      |:|  |"
    write(stdout, '(A)') "   e⁻ |:|  |        \::/  /     \/__/   e⁻      \:\__\        |:|  |"
    write(stdout, '(A)') "       \|__|         \/__/                       \/__/         \|__|"
    write(stdout, *)
    write(stdout, *)
    write(stdout, '(A)') "              ROTational (de-)EXcitation by electron impact"
    write(stdout, *)
    write(stdout, *)
    ! -- determine system / environment properies for system interaction later on
    call determine_system_properties
    write(stdout, "(2A)") "Fortran compiler and version :: ", compiler_version()
    write(stdout, *)
    write(stdout, "(2A)") "Fortran compiler options :: ", compiler_options()
    write(stdout, *)
#ifdef USE_OPENMP
    !$omp parallel
    if(omp_get_thread_num() .eq. 0) write(stdout, '("Number of available threads via OpenMP: ", I0)') &
      omp_get_num_threads()
    !$omp end parallel
#else
    write(stdout, '("OpenMP has *not* been explicitly requested")')
#endif
  end subroutine print_header

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_footer(time_start, time_end)
    use rotex__kinds, only: dp
    use rotex__system,     only: stdout
    use rotex__characters, only: s2hms
    implicit none
    real(dp) :: time_start
    real(dp) :: time_end
    write(stdout, *)
    write(stdout, '(A)') "========================================================================================================"
    write(stdout,'(A)') "Program complete (^:"
    write(stdout, *)
    write(stdout,'("Elapsed time: ",A)') s2hms(time_end - time_start)
    write(stdout, '(A)') "========================================================================================================"
  end subroutine print_footer

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_rot_targ_states(n_states)
    use rotex__types,      only: N_states_type
    use rotex__system,     only: stdout
    use rotex__characters, only: lower
    use rotex__constants,  only: au2ev
    implicit none
    type(N_states_type), intent(in) :: n_states(:)
    integer  :: i,j,n,ka,kc,k,jstart
    real(dp) :: e
    character(:), allocatable :: fmt
    write(stdout,*)
    write(stdout, '(A)') "Rotational target states"
    write(stdout, '(A)') "⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻"
    select case(G%ROTOR_KIND)
    case("s","S")
      write(stdout, '(4X, 2A5, A14)') "N", "K", "E (meV)"
      do i=1, size(N_states, 1)
        n = n_states(i)%n
        do j=1,2*n+1
          select case(G%SYMAXIS)
          case("a", "A")
            k = n_states(i)%ka(j)
          case("c", "C")
            k = n_states(i)%kc(j)
          case default
            call die("Unacceptable G%SYMAXIS in print_rot_targ_states")
          end select
          ! -- print |K| only
          if(k.lt.0) cycle
          e  = n_states(i)%eigenh%eigvals(j)*au2ev*1000
          ! sym = n_states(i)%eigenh%sym(j)
          if(abs(e) .lt. 0.001_dp) then
            fmt =  '(4X, 2I5, E14.5)'
          else
            fmt =  '(4X, 2I5, F14.5)'
          endif
          write(stdout, fmt) n, k, e
        enddo
      enddo
    case("a","A")
      write(stdout, '(4X, 3A5, A14)') "N", "Ka", "Kc", "E (meV)"
      do i=1, size(N_states, 1)
        n = n_states(i)%n
        do j=1,2*n+1
          ka = n_states(i)%ka(j)
          kc = n_states(i)%kc(j)
          e  = n_states(i)%eigenh%eigvals(j)*au2ev*1000
          ! sym = n_states(i)%eigenh%sym(j)
          if(abs(e) .lt. 0.001_dp) then
            fmt =  '(4X, 3I5, E14.5)'
          else
            fmt =  '(4X, 3I5, F14.5)'
          endif
          write(stdout, fmt) n, ka, kc, e
        enddo
      enddo
    case("l")
      write(stdout, '(4X, A5, A14)') "N", "E (meV)"
      do i=1, size(N_states, 1)
        n = n_states(i)%n
        e = n_states(i)%eigenh%eigvals(1)*au2ev*1000
        if(abs(e) .lt. 0.001_dp) then
          fmt = '(4X, I5, E14.5)'
        else
          fmt = '(4X, I5, F14.5)'
        endif
        write(stdout, fmt) n, e
      enddo
    case default
      call die("Unacceptable G%ROTOR_KIND "// G%ROTOR_KIND // " in print_rot_targ_states")
    end select
    write(stdout, *)
  end subroutine print_rot_targ_states

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_dipole(dipole)
    !! Print the dipole components in the determined ABC frame
    use rotex__system,    only: stdout
    use rotex__constants, only: au2deb
    implicit none
    real(dp),     intent(in) :: dipole
    write(stdout, '("Permanent dipole moment: ")')
    write(stdout, '("μ: ", F7.4, " Debye")') dipole*au2deb
  end subroutine print_dipole
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_dipoles
    !! Print the dipole components in the determined ABC frame
    use rotex__system,    only: stdout
    use rotex__constants, only: au2deb
    implicit none
    write(stdout, '("Cartesian dipole moments in the inertial frame ABC:")')
    write(stdout, '("μ(A): ", F7.4, " Debye")') G%DIPOLE_ABC(1)*au2deb
    write(stdout, '("μ(B): ", F7.4, " Debye")') G%DIPOLE_ABC(2)*au2deb
    write(stdout, '("μ(C): ", F7.4, " Debye")') G%DIPOLE_ABC(3)*au2deb
  end subroutine print_dipoles

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_datetime(funit)
    implicit none(type, external)
    integer, intent(in) :: funit
    character(5)  :: zone
    integer :: values(8)
    call date_and_time(zone=zone, values=values)
    write(funit,*)
    write(funit, '("Date [ymd]: ", I4.4,"-",I2.2,"-",I2.2)') values(1:3)
    write(funit, '("Time [hms]: ", I4.2,"-",I2.2,"-",I2.2)') values(5:7)
    write(funit, '("Time zone: ", A)') zone
    write(funit,*)
  end subroutine print_datetime


  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine validate_user_input
  !   !! Validates the  G global config
  !   implicit none(type, external)
  !   stop "VALIDATE ROUTINE"
  !     ! if(is_supported_pg(target_point_group) .eqv. .false.) call die("Invalid TARGET_POINT_GROUP: "//target_point_group)
  !     ! if(is_supported_pg(scattering_point_group) .eqv. .false.) &
  !     !   call die("Invalid SCATTERING_POINT_GROUP: "//target_point_group)
  !     ! if(is_abelian_pg(scattering_point_group) .eqv. .false.) &
  !     !   call die("SCATTERING POINT GROUP "//scattering_point_group//" is not Abelian. I don't believe you !&
  !     ! & (please supply the scattering point group that was used for the scattering calculation)")
  ! end subroutine validate_user_input

! ================================================================================================================================ !
end program rotex
! ================================================================================================================================ !
