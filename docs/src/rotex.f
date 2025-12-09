! ================================================================================================================================ !
program rotex
  !! The main program

  use rotex__drivers,    only: make_grid, diagonalize_rotational_hamiltonian &
    , make_output_directories, do_coulomb_born_approx, do_kmat_xs, combine_cb_smat_xs
  use rotex__reading,    only: read_namelists
  use rotex__types,      only: dp, eigenh_type, n_states_type, config_type, asymtop_rot_channel_l_type &
    , asymtop_rot_channel_l_vector_type, cmatrix_type, rvector_type, asymtop_rot_transition_type
  use rotex__rft,        only: rft_nonlinear
  use rotex__system,     only: mkdir, die, stdout
  use rotex__hamilton,   only: h_asym, assign_projections

  implicit none


  integer :: jmin, jmax
    !! The min/max values of J = N + l to consider
  integer :: nmin_ft, nmax_ft
    !! The min/max values of N to consider for the frame transformation.
  integer :: i_n, num_n
  integer, allocatable :: N_values(:)
    !! Contains the values of N to consider

  real(dp) :: start, finish, Eground
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

  type(config_type):: cfg

  ! GO
  start = time()

  call print_header()
  call read_namelists(cfg)
  call make_output_directories( cfg%output_directory, cfg%use_CB, cfg%spinmults, cfg%use_kmat &
                              , pcb_output_directory, tcb_output_directory, smat_output_directory)

  ! -- determine the (number of) N values
  if(cfg%use_kmat .eqv. .true.) then
    ! -- if we  use the K-matrix, then the S-matrix will be block
    !    diagonal in J = [ max(0, Nmin - lmax_kmat), Nmax + lmax_kmat ], but for
    !    each J-block of S to be unitary, we must calculate extra rotational wavefunctions
    !    N = [ max(0, Jmin - lmax_kmat), Jmax + lmax  ]
    jmin = max(0, cfg%nmin - cfg%lmax_kmat)
    jmax = cfg%nmax + cfg%lmax_kmat
    nmin_ft = max(0, jmin - cfg%lmax_kmat)
    nmax_ft = jmax + cfg%lmax_kmat
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
    num_n    = cfg%nmax - cfg%nmin + 1
    n_values = [(cfg%nmin + i_n, i_n=0, cfg%nmax - cfg%nmin, 1)]
    write(stdout, '(a, i0, "/", i0)') "not using a k-matrix to get cross sections, so we only need the user-supplied&
      & values for nmin/nmax: "&
      , cfg%nmin, cfg%nmax
    write(stdout, *)
  endif

  ! -- this array holds the information on the rotational states of the target
  allocate(N_states(num_N))

  ! -- diagonalize the hamiltonian and assign state labels for the rotational states that
  !    will be involved in the transitions/collisions
  select case(cfg%zaxis)
  case("a") ; call print_dipoles(cfg%cartesian_dipole_moments, x = "B", y = "C", z = "A")
  case("b") ; call print_dipoles(cfg%cartesian_dipole_moments, x = "C", y = "A", z = "B")
  case("c") ; call print_dipoles(cfg%cartesian_dipole_moments, x = "A", y = "B", z = "C")
  end select

  call diagonalize_rotational_hamiltonian(cfg, num_N, N_values, N_states)

  if(cfg%use_CB .eqv. .true.) then
    call do_coulomb_born_approx( &
        cfg                      &
      , n_states                 &
      , Eel_grid_cb              &
      , transitions_cb           &
      , xs_xcite_pcb             &
      , xs_xcite_tcb             &
      , pcb_output_directory     &
      , tcb_output_directory)
  endif

  ! -- do this only AFTER we have called DO_COULOMB_BORN_APPROX because it may use
  !    CDMS energies which will change the energies of our rotational levels (but not the eigenvectors)
  call print_rot_targ_states(n_states)

  usingkmat: if(cfg%use_kmat) then

    write(stdout, *)
    write(stdout, '(A)') "--------------------------"
    write(stdout, '(A)') "Switching to K-matrix data"
    write(stdout, '(A)') "--------------------------"
    write(stdout, *)

    ! -- make the total energy grid for the K/S-matrix cross sections
    eground = n_states(1) % eigenh % eigvals(1)
    call make_grid(egrid_tot_smat, eground, cfg%num_egrid_segs, cfg%egrid_segs, cfg%num_egrid, cfg%egrid_spacing)

    call do_kmat_xs(          &
        cfg                   &
      , n_states              &
      , egrid_tot_smat        &
      , smat_output_directory &
      , transitions_smat      &
      , xs_xcite_smat         &
      , xs_dxcite_smat)

    if(cfg%use_cb .eqv. .false.) exit usingkmat

    ! -- now we combine the S-matrix cross sections (which are spin-averaged) with the
    !    Coulom-Born correcion (σTCB - σPCB) which are spin-independent. Each CB σ has its own
    !    energy grid, while the σSmat cross sections were all calculated on the same grid of total
    !    energies. The CB cross sections are resonance free, so we'll have to interpolate these
    !    to match the σSmat energy grid before we can do σTot = σSmat + σTCB - σPCB
    write(stdout, *)
    write(stdout, '("Combining Coulomb-Born and S-matrix cross sections.")')
    write(stdout, '("⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻")')

    call combine_cb_smat_xs( &
        cfg &
      , Eel_grid_cb          &
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
    use rotex__system,   only: stdout, determine_system_properties
    use iso_fortran_env, only: compiler_version, compiler_options
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
  end subroutine print_header

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_footer(time_start, time_end)
    use rotex__types,      only: dp
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
    use rotex__types,     only: N_states_type
    use rotex__system,    only: stdout
    use rotex__constants, only: au2ev
    implicit none
    type(N_states_type), intent(in) :: n_states(:)
    integer  :: i,j,n,ka,kc
    real(dp) :: e
    character(:), allocatable :: fmt
    write(stdout,*)
    write(stdout, '(A)') "Rotational target states"
    write(stdout, '(A)') "⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻"
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
    write(stdout, *)
  end subroutine print_rot_targ_states

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_dipoles(dipole_xyz, x, y, z)
    !! Print the dipole components in the determined ABC frame
    use rotex__system,    only: stdout
    use rotex__constants, only: au2deb
    implicit none
    real(dp),     intent(in) :: dipole_xyz(3)
    character(1), intent(in) :: x, y, z
    write(stdout, '("Cartesian dipole moments  in the inertial frame ABC:")')
    write(stdout, '("μ(", A1, "): ", F7.4, " Debye")') x, dipole_xyz(1)*au2deb
    write(stdout, '("μ(", A1, "): ", F7.4, " Debye")') y, dipole_xyz(2)*au2deb
    write(stdout, '("μ(", A1, "): ", F7.4, " Debye")') z, dipole_xyz(3)*au2deb
  end subroutine print_dipoles

! ================================================================================================================================ !
end program rotex
! ================================================================================================================================ !
