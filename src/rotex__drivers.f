! ================================================================================================================================ !
module rotex__drivers
  !! Driver used by the PROGRAM in main

  use rotex__kinds,     only: dp
  use rotex__globals,   only: G
  use rotex__constants, only: UKRMOLX, MQDTR2K, DEFAULT_CHAR1

  implicit none (type, external)

  private

  public :: do_coulomb_born_approx
  public :: do_kmat_xs
  public :: make_grid
  public :: convert_multipoles
#ifdef USE_CDMSREADER
  public :: get_cdms_data
  public :: get_cdms_einsta
  public :: get_cdms_state_energies
#endif
  public :: diagonalize_rotational_hamiltonian
  public :: make_output_directories
  public :: combine_cb_smat_xs

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine do_coulomb_born_approx( &
        n_states                            &
      , egrid_elec_cb                       &
      , transitions_cb                      &
      , xs_xcite_pcb                        &
      , xs_xcite_tcb                        &
      , pcb_output_directory                &
      , tcb_output_directory                &
    )
    !! Use the Coulomb-Born approximation to get scattering cross sections for e⁻ + target.
    !! Also determines Einstein A coefficients and excited state average lifetimes within
    !! !! the radiative dipole/multipole approximation. The target charge Z is supplied and
    !! used by the routine that are called here, so various target charges can be used, including Z = 0 (neutral targets).
    !! This routine loops over the different pairs  Nτ —→ N'τ'

    use rotex__kinds,      only: dp
    use rotex__arrays,     only: append
    use rotex__symmetry,   only: is_spin_forbidden, spin_symmetry, symtop_rotstate_is_allowed
    use rotex__system,     only: stdout, die
    use rotex__functions,  only: istriangle, logrange
    use rotex__types,      only: asymtop_rot_channel_type, asymtop_rot_transition_type, n_states_type, rvector_type
    use rotex__channel_ops, only: operator(.eq.)
    use rotex__writing,    only: write_CB_xs_to_file, write_lifetimes_to_file
    use rotex__cbxs,       only: get_cb_xs_asym, get_einsta_only, xtrapolate_cb_xs
    use rotex__characters, only: i2c => int2char
#ifdef USE_CDMSREADER
    use CDMSreader__types, only: asymtop_state_type      => asymtop_state_nohfs &
                               , asymtop_transition_type => asymtop_transition_nohfs
#endif

    implicit none (type, external)

    type(n_states_type), intent(inout), allocatable :: n_states(:)
    type(rvector_type),  intent(out), allocatable :: egrid_elec_cb(:)
      !! Array of arrays of electron/collision energies for each transition
    type(asymtop_rot_transition_type), intent(out), allocatable :: transitions_cb(:)
      !! Array containing transitions between rotational states for the long-range CB excitations
    type(rvector_type),  intent(out), allocatable :: xs_xcite_pcb(:), xs_xcite_tcb(:)
      !! Array of arrays of the Partial and Total Coulomb-Born cross sections for each transition
    character(*),        intent(in) :: pcb_output_directory
    character(*),        intent(in) :: tcb_output_directory

    integer :: ipair, inlo, inup, itaulo, itauup, itrans
    integer :: nlo, nup, kalo, kaup, kclo, kcup, neleclo, nelecup
    integer :: num_klo, num_kup
    integer :: num_n
    integer :: sym
    integer, allocatable :: lambdas(:)

    real(dp) :: elo, eup, de, estart, eend, ei
    real(dp) :: einsta
    real(dp), allocatable :: Eel(:)
    real(dp), allocatable :: sigma_pcb(:), sigma_tcb(:)
    complex(dp), allocatable :: eigveclo(:), eigvecup(:)

    complex(dp) :: spherical_dipole_moments(3)
    complex(dp) :: spherical_quadrupole_moments(5)

    type(asymtop_rot_channel_type) :: lo, up
      !! Rotational channels of the target for an excitation pair
    type(asymtop_rot_transition_type) :: transition
      !! A single rotational transition
#ifdef USE_CDMSREADER
    type(asymtop_state_type),          allocatable :: cdms_states(:)
      !! Array of states involved in the CDMS transitions, averaged over hyperfine contributions. Output
      !! of CDMSreader
    type(asymtop_transition_type),     allocatable :: CDMS_transitions(:)
      !! Array of CDMS transitions, averaged over hyperfine contributions. Output of CDMSreader
#endif

    num_n = size(n_states, 1)

    ! -- cartesian multipoles → spherical multipoles.
    if(G%DO_DIPOLE     .eqv. .true.) &
      call convert_multipoles(G%CARTESIAN_DIPOLE_MOMENTS,     spherical_dipole_moments)
    if(G%DO_QUADRUPOLE .eqv. .true.) &
      call convert_multipoles(G%CARTESIAN_QUADRUPOLE_MOMENTS, spherical_quadrupole_moments)

    lambdas = pack([1, 2],  [G%DO_DIPOLE, G%DO_QUADRUPOLE])

#ifdef USE_CDMSREADER
    ! -- whether to use the CDMS einstein A coefficients in our Coulomb-Born approximation
    if(G%USE_CDMS_EINSTA .eqv. .true.) then
      call get_cdms_data(G%CDMS_FILE, G%OUTPUT_DIRECTORY, cdms_states, cdms_transitions)
      call get_cdms_state_energies(n_states, cdms_states) ! n_states <— cdms data
    endif
#endif

    write(stdout, *)
    write(stdout, '(A)') "--------------------------------------------------"
    write(stdout, '(A)') "Looping over allowed Coulomb-Born excitation pairs"
    write(stdout, '(A)') "⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻"

    ipair =  0
    nlo_loop: do inlo = 1, num_n

      nlo = n_states(inlo) % n
      if(nlo .lt. G%NMIN) cycle nlo_loop
      if(nlo .gt. G%NMAX) cycle nlo_loop

      ! -- only consider 1 electronic state for now
      neleclo = 1

      select case(G%ROTOR_KIND)
      case("l")      ; num_klo = 1
      case("a", "s") ; num_klo = 2*nlo+1
      case default   ; call die("ROTOR_KIND ( "// G%ROTOR_KIND //" )is neither 'l', 'a', or 's'")
      end select

      ! -- only consider excitation pairs; de-excitation is handled symmetrically
      nup_loop: do inup = inlo, num_n

        nup = n_states(inup) % n
        if(nup .lt. G%NMIN) cycle nup_loop
        if(nup .gt. G%NMAX) cycle nup_loop

        ! -- only consider 1 electronic state for now
        nelecup = 1

        ! -- skip if the pair Nlo, Nup don't satisfy the triangle inequality
        !    |Nlo-Nup| <= λ <= Nlo+Nup for any λ
        if(all( istriangle(nlo, nup, lambdas) .eqv. .false. )) cycle

        write(stdout, '("N : ", I0, " —> ", I0)') nlo, nup

        select case(G%ROTOR_KIND)
        case("l")      ; num_kup = 1
        case("a", "s") ; num_kup = 2*nup+1
        end select

        taulo_loop: do itaulo=1, num_klo

          select case(G%ROTOR_KIND)
          case("l")

            kalo = 0
            kclo = nlo

          case("a")

            kalo = n_states(inlo) % ka(itaulo)
            kclo = n_states(inlo) % kc(itaulo)

          case("s")

            select case(G%ROTOR_ZAXIS)
            case("A","a")

              kalo = n_states(inlo) % ka(itaulo)
              kclo = 0

              ! -- skip forbidden state
              if(symtop_rotstate_is_allowed(nlo, Kalo) .eqv. .false.) cycle taulo_loop

            case("C", "c")

              kalo = 0
              kclo = n_states(inlo) % kc(itaulo)

              ! -- skip forbidden state
              if(symtop_rotstate_is_allowed(nlo, Kclo) .eqv. .false.) cycle taulo_loop

            case default
              call die("G%ROTOR_ZAXIS must be A or C")
            end select

          end select
          elo  = n_states(inlo) % eigenh % eigvals(itaulo)
          sym  = spin_symmetry(nlo, kalo, kclo)
          lo = asymtop_rot_channel_type(nelec= neleclo, n=nlo, ka=kalo, kc=kclo, e=elo, sym=sym)

          tauup_loop: do itauup=1, num_kup

            select case(G%ROTOR_KIND)
            case("l")

              kaup = 0
              kcup = nup

            case("a")

              kaup = n_states(inup) % ka(itauup)
              kcup = n_states(inup) % kc(itauup)

            case("s")

              select case(G%ROTOR_ZAXIS)
              case("A","a")

                kaup = n_states(inup) % ka(itauup)
                kcup = 0

                ! -- skip forbidden state
                if(symtop_rotstate_is_allowed(nup, Kcup) .eqv. .false.) cycle tauup_loop

              case("C", "c")

                kaup = 0
                kcup = n_states(inup) % kc(itauup)

                ! -- skip forbidden state
                if(symtop_rotstate_is_allowed(nup, Kcup) .eqv. .false.) cycle tauup_loop

              case default
                call die("G%ROTOR_ZAXIS must be A or C")
              end select
            end select

            eup  = n_states(inup) % eigenh % eigvals(itauup)
            sym  = spin_symmetry(nup, kaup, kcup)
            up = asymtop_rot_channel_type(nelec= nelecup, n=nup, ka=kaup, kc=kcup, e=eup, sym=sym)

            ! -- ignore elastic collisions and de-excitation pairs
            if(lo  .eq. up)  cycle tauup_loop
            if(Eup .lt. Elo) cycle

            ! -- add this pair of states to the list of transitions
            ipair = ipair + 1
            transition = asymtop_rot_transition_type(lo = lo, up = up)

            select case(G%ROTOR_KIND)
            case("a", "s")
              write(stdout, '(2X, "(Ka,Kc) : ", I0,",",I0, " --> ", I0,",",I0, " ... ")', advance = "no") &
                kalo, kclo, kaup, kcup
              eigveclo = n_states(inlo) % eigenh % eigvecs(:, itaulo)
              eigvecup = n_states(inup) % eigenh % eigvecs(:, itauup)
            end select

            ! -- check if we need to respect ortho-para symmetry
            if( is_spin_forbidden(nlo, kalo, kclo, nup, kaup, kcup, G%SPIN_ISOMER_KIND, G%ROTOR_ZAXIS) ) then
              write(stdout, "(A)") "is forbidden (ortho-para violation) !"
              cycle tauup_loop
            endif

            ! -- coulomb-born energy grid starts at the excitation threshold + Ei, and goes to
            !    some absolute electron energy ef that does not depend on the threshold
            de = eup - elo
            ! -- the starting value of the calculation grid based on the maximum allowed value of
            !    η' for an excitation.
            ei = 1.0_dp / (2.0_dp*G%ETA_THRESH**2)
            estart = de + ei
            eend = G%EF
            if(eend .le. estart) call die("The starting energy of this coulomb-born energy grid is higher than the ending energy")
            Eel = logrange(estart, eend, G%NE)

            einsta = 0

#ifdef USE_CDMSREADER
            ! -- get CDMS Einstein A coefficients
            if(G%USE_CDMS_EINSTA .eqv. .true.) then
              call get_cdms_einsta(nlo, kalo, kclo, nup, kaup, kcup, cdms_transitions, einsta)
              select case(G%ROTOR_KIND)
              case("a") ; continue
              case("s") ; call die("Can't get CDMS for a symmetric top yet")
              case("l") ; call die("Can't get CDMS for linear rotor yet")
              end select
            endif
#endif

            if(G%ONLY_EINSTA .eqv. .false.) then

              ! -- get partial CB cross sections Nτ —> N'τ'
              call get_cb_xs_asym(             &
                  Eel                          &
                , sigma_pcb                    &
                , nlo                          &
                , nup                          &
                , elo                          &
                , eup                          &
                , eigveclo                     &
                , eigvecup                     &
                , einsta                       &
                , spherical_dipole_moments     &
                , spherical_quadrupole_moments &
                , [.false., .false.]           & ! explicitly request partial xs
                , G%LMAX_PARTIAL               &
              )

            elseif(einsta .eq. 0) then

              ! -- get calculate only Einstein A coefficients if they're not found in the CDMS
              call get_einsta_only(            &
                  einsta                       &
                , nlo                          &
                , nup                          &
                , elo                          &
                , eup                          &
                , eigveclo                     &
                , eigvecup                     &
                , spherical_dipole_moments     &
                , spherical_quadrupole_moments &
              )

            endif

            ! -- add the einstein coefficients to the upper state for this transition
            n_states(inup) % einsta(itauup) = n_states(inup) % einsta(itauup) + einsta

            if(G%ONLY_EINSTA .eqv. .true.) then
              if(EinstA .eq. 0) then
                write(stdout, "(A)") "is forbidden (Einstein A coefficient is 0)!"
              else
                write(stdout, "(A)") "done !"
              endif
              cycle tauup_loop
            endif

            if(all(sigma_pcb .eq. 0)) then
              write(stdout, "(A)") "is forbidden (all cross sections are 0)!"
              cycle tauup_loop
            endif

            ! -- append to the list of transitions, electron energies, and PCB cross sections
            call append(transitions_CB, transition)

            ! -- only include the grid elements where we successfully got a well-behaved
            !    cross section. This defines the shared Eel Coulomb-Born grid that will be
            !    used for the TCB cross sections as well
            Eel       = pack(Eel,       sigma_pcb .gt. 0.0_dp)
            sigma_pcb = pack(sigma_pcb, sigma_pcb .gt. 0.0_dp)

            ! -- get total CB cross sections Nτ —> N'τ'
            call get_cb_xs_asym(             &
                Eel                          &
              , sigma_tcb                    &
              , nlo                          &
              , nup                          &
              , elo                          &
              , eup                          &
              , eigveclo                     &
              , eigvecup                     &
              , einsta                       &
              , spherical_dipole_moments     &
              , spherical_quadrupole_moments &
              , G%ANALYTIC_TOTAL_CB          & ! as per user request
              , G%LMAX_TOTAL                 &
            )

            ! -- extrapolate excitation CB cross sections as 1/E to threshold ?
            ! if(G%do_xtrap) call xtrapolate_cb_xs(G%Ei_xtrap, dE, G%nE_xtrap, Eel, sigma_pcb, sigma_tcb)
            if(G%DO_XTRAP) call xtrapolate_cb_xs(G%EI_XTRAP, dE, G%NE_XTRAP, Eel, sigma_pcb, sigma_tcb)

            ! -- append energy grid, including extrapolated energies if that happened
            call append(egrid_elec_cb, Eel)

            ! -- append PCB and TCB excitation to corresponding arrays, including
            !    extrapolated cross sections if that happened
            call append(xs_xcite_pcb, sigma_pcb)
            call append(xs_xcite_tcb, sigma_tcb)

            write(stdout, "(A)") "done !"

          enddo tauup_loop

        enddo taulo_loop

      enddo nup_loop

    enddo Nlo_loop

    call write_lifetimes_to_file(n_states)

    ! -- exit if we don't have to write cross sections
    if(G%ONLY_EINSTA) return

    ! -- reduce the resolution on all states to just have |K|
    if(G%symtop_reduce_projection .AND. G%rotor_kind .eq. "s") &
      call reduce_symtop_ksign(transitions_CB, xs_xcite_pcb, xs_xcite_tcb)

    ! -- write the cross (potentially reduced) cross sections to disk
    do itrans = 1, size(transitions_CB)

        transition = transitions_CB(itrans)
        sigma_pcb  = xs_xcite_pcb(itrans)%vec(:)
        sigma_tcb  = xs_xcite_tcb(itrans)%vec(:)
        Eel        = egrid_elec_cb(itrans)%vec(:)

        ! -- write this excitation to disk
        call write_CB_xs_to_file(    &
            "PCB"                 &
          , pcb_output_directory  &
          , Eel                   &
          , sigma_pcb             &
          , transition%lo         &
          , transition%up         &
          , i2c(G%LMAX_PARTIAL) &
        )
        call write_CB_xs_to_file(    &
            "TCB"                 &
          , tcb_output_directory  &
          , Eel                   &
          , sigma_tcb             &
          , transition%lo         &
          , transition%up         &
          , i2c(G%LMAX_PARTIAL) &
        )

        ! -- get corresponding de-excitation cross section
        call convert_xcite2dxcite(Eel, sigma_pcb, transition)
        call convert_xcite2dxcite(Eel, sigma_tcb, transition)

        dE = transition % up % E - transition % lo % E

        ! -- write this de-excitation to disk
        call write_CB_xs_to_file( &
            "PCB"                 &
          , pcb_output_directory  &
          , Eel - dE              &
          , sigma_pcb             &
          , transition%up         &
          , transition%lo         &
          , i2c(G%LMAX_PARTIAL) &
        )
        call write_CB_xs_to_file( &
            "TCB"                 &
          , tcb_output_directory  &
          , Eel - dE              &
          , sigma_tcb             &
          , transition%up         &
          , transition%lo         &
          , i2c(G%LMAX_PARTIAL) &
        )

    enddo

  end subroutine do_coulomb_born_approx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine do_kmat_xs( &
        n_states                &
      , egrid_tot_smat          &
      , smat_output_directory   &
      , transitions             &
      , xs_xcite_spinavg        &
      , xs_dxcite_spinavg       &
    )
    !! Read K-matrices from an electron-molecule scattering calculation, get S-matrices, add the rotation via
    !! the rotational frame transformation for asymmetric tops, then use the MQDT channel elimination to
    !! get electron-impact excitation cross sections entirely from the K-matrix scattering data.

    use rotex__kinds,     only: dp
    use rotex__channel_ops, only: findloc_transitions
    use rotex__types,     only: n_states_type, cmatrix_type, elec_channel_type &
                           , asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type &
                           , asymtop_rot_transition_type, rvector_type
    use rotex__system,    only: die, DS => DIRECTORY_SEPARATOR, stdout
    use rotex__constants, only: SPINMULT_NAMES
    use rotex__rft,       only: rft_nonlinear
    use rotex__arrays,    only: append_uniq
    use rotex__mqdtxs,    only: get_smat_probs
    use rotex__characters, only: i2c => int2char
    use rotex__reading,   only: read_kmats
    use rotex__writing,   only: write_smat_xs_to_file, write_channels_to_file

    implicit none (type, external)

    type(n_states_type), intent(in) :: n_states(:)
    real(dp), intent(in) :: egrid_tot_smat(:)
      !! Total energy grid for the S-matrix cross sections
    character(*),        intent(in) :: smat_output_directory
    type(asymtop_rot_transition_type), intent(out), allocatable :: transitions(:)
    type(rvector_type),  intent(out), allocatable :: xs_xcite_spinavg(:)
      !! Array of arrays of excitation cross sections for all spin multiplicities
    type(rvector_type),  intent(out), allocatable :: xs_dxcite_spinavg(:)
      !! Array of arrays of de-excitation cross sections for all spin multiplicities

    ! -- default values of the channel and K-matrix energies
    character(1), parameter :: UKRMOLX_CHANNEL_ENERGY_UNITS = "r"!ydberg
    character(1), parameter :: MQDTR2K_CHANNEL_ENERGY_UNITS = "e"!ydberg
    character(1), parameter :: UKRMOLX_KMAT_ENERGY_UNITS    = "r"!ydberg
    character(1), parameter :: MQDTR2K_KMAT_ENERGY_UNITS    = "e"!lectron-Volts

    integer :: ispin, nspins, jmin, jmax, itrans, ntrans, ne
    integer, allocatable :: idxmap(:)

    real(dp), allocatable :: kmat(:,:)

    character(1) :: kmat_eval_E_units, channel_e_units
    character(:), allocatable :: smat_output_directory_this_spin, smat_output_directory_all_spins
    character(:), allocatable :: channels_file_this_spin

    type(cmatrix_type),                      allocatable :: smat_j(:)
    type(elec_channel_type),                 allocatable :: elec_channels(:)
    type(rvector_type),                      allocatable :: prob_smat(:), xs_xcite(:), xs_dxcite(:)
    type(asymtop_rot_channel_l_type),        allocatable :: asymtop_rot_channels_l(:)
    type(asymtop_rot_channel_l_vector_type), allocatable :: asymtop_rot_channels_l_j(:)
    type(asymtop_rot_transition_type),       allocatable :: transitions_this_spin(:)

    nspins = size(G%SPINMULTS, 1)
    spinsdo: do ispin= 1, nspins

      ! -- determine energy and channel units
      kmat_eval_e_units = G%KMAT_ENERGY_UNITS_OVERRIDE
      channel_e_units   = G%CHANNEL_ENERGY_UNITS_OVERRIDE
      ! -- check if we need to use default values
      select case(G%KMAT_OUTPUT_TYPE)
      case(UKRMOLX)
        if(kmat_eval_e_units .eq. DEFAULT_CHAR1) kmat_eval_e_units = UKRMOLX_KMAT_ENERGY_UNITS
        if(channel_e_units   .eq. DEFAULT_CHAR1) channel_e_units   = UKRMOLX_CHANNEL_ENERGY_UNITS
      case(MQDTR2K)
        if(kmat_eval_e_units .eq. DEFAULT_CHAR1) kmat_eval_e_units = MQDTR2K_KMAT_ENERGY_UNITS
        if(channel_e_units   .eq. DEFAULT_CHAR1) channel_e_units   = MQDTR2K_CHANNEL_ENERGY_UNITS
      case default
        call die("KMAT_OUTPUT_TYPE ("//G%KMAT_OUTPUT_TYPE//") must be one of "//UKRMOLX//" or "//MQDTR2K)
      end select

      call read_kmats(                          &
          kmat              = kmat              &
        , spinmult = G%SPINMULTS(ispin) &
        , elec_channels     = elec_channels     &
        , channel_e_units   = channel_e_units   &
        , kmat_eval_E_units = kmat_eval_e_units &
        )

      if(maxval(elec_channels % l) .gt. G%LMAX_KMAT) call die("K-matrix has at least one channel with&
        & l > LMAX_KMAT: " // i2c(maxval(elec_channels % l)) // " > " // i2c(G%LMAX_KMAT))

      write(stdout, '(A)') "Performing the rotational frame transformation"
      write(stdout, '(A)') "----------------------------------------------"
      write(stdout, *)

      select case(G%ROTOR_KIND)
      case("l")
        call die("Linear RFT not programmed yet. Just use ABC with large A (hundreds+)")
      case("s", "a")
        continue
      case default
        call die("ROTOR_KIND  must be (l)inear, (s)ymmetric top, or (a)symmetric top")
      end select

      ! -- min and max values of total J
      jmin = max(0, G%NMIN - G%LMAX_KMAT)
      jmax = abs(G%NMAX + G%LMAX_KMAT)
      allocate(smat_j(jmin:jmax))
      allocate(asymtop_rot_channels_l_j(jmin:jmax))
      call rft_nonlinear( kmat                     &
                        , jmin, jmax               &
                        , smat_j(jmin:jmax)        &
                        , elec_channels            &
                        , n_states                 &
                        , asymtop_rot_channels_l   &
                        , asymtop_rot_channels_l_j(jmin:jmax) &
      )
      deallocate(elec_channels)


      channels_file_this_spin = G%OUTPUT_DIRECTORY // SPINMULT_NAMES(G%SPINMULTS(ispin)) // ".channels"

      call write_channels_to_file(            &
          channels_file_this_spin             &
        , jmin                                &
        , jmax                                &
        , n_states                            &
        , asymtop_rot_channels_l              &
        , asymtop_rot_channels_l_j(jmin:jmax) &
      )

      write(stdout, *)
      write(stdout, '(A)') "--------------------------"
      write(stdout, '(A)') "Calculating cross sections"
      write(stdout, '(A)') "⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻"

      call get_smat_probs(         &
          egrid_tot_smat           &
        , prob_smat                &
        , transitions_this_spin    &
        , smat_j(jmin:jmax)        &
        , jmin, jmax               &
        , asymtop_rot_channels_l_j &
        , asymtop_rot_channels_l   &
        )


      deallocate(smat_J)
      deallocate(asymtop_rot_channels_l)
      deallocate(asymtop_rot_channels_l_j)

      ! -- this routine also remove transitions whose cross sectiosn are too small
      call get_xs_from_smat(    &
          prob_smat             &
        , transitions_this_spin &
        , egrid_tot_smat        &
        , xs_xcite              &
        , xs_dxcite             &
      )

      write(stdout, *)
      write(stdout, '(A)') "Writing " // SPINMULT_NAMES(G%SPINMULTS(ispin)) // " S-matrix cross sections to disk"
      write(stdout, *)

      smat_output_directory_this_spin = &
        & smat_output_directory // SPINMULT_NAMES(G%SPINMULTS(ispin)) // DS

      do itrans = 1, size(transitions_this_spin, 1)
        call write_smat_xs_to_file(              &
            "Smat"                          &
          , smat_output_directory_this_spin &
          , egrid_tot_smat                  &
          , transitions_this_spin(itrans)   &
          , xs_xcite(itrans)%vec            &
          , xs_dxcite(itrans)%vec           &
          , G%LMAX_KMAT                   &
        )
      enddo

      ! -- append the transitions for this spin to the array of transitions
      !    for all spins uniquely (don't double-count existing transitions).
      call append_uniq(transitions, transitions_this_spin)

      ! -- if there's only one spin, don't worry about spin averaging
      if(nspins .eq. 1) then
        call move_alloc(xs_xcite,  xs_xcite_spinavg)
        call move_alloc(xs_dxcite, xs_dxcite_spinavg)
        return
      endif

      ! -- xs_xcite, xs_dxcite,and transitions_this_spin are always indexed in the same basis.
      !    This is probably true for transitions, but not necessarily.Figure out how
      !    they map to transitions so that we ensure contributions are added properly
      idxmap = findloc_transitions(transitions_this_spin, transitions)

      ! -- allocate the spinavg cross section arrays if necessary before adding this spin's contribution
      ne     = size(egrid_tot_smat, 1)
      ntrans = size(transitions_this_spin, 1)
      if(.not. allocated(xs_xcite_spinavg)) then
        allocate(xs_xcite_spinavg(ntrans))
        allocate(xs_dxcite_spinavg(ntrans))
        do itrans=1, ntrans ; allocate(xs_xcite_spinavg(itrans)%vec(ne),  source = 0.0_dp) ; enddo
        do itrans=1, ntrans ; allocate(xs_dxcite_spinavg(itrans)%vec(ne), source = 0.0_dp) ; enddo
      endif

      ! -- σ_allspins += σ_thisspin * (2S+1) / Σ(2S+1)
      call add_xs_contrib_this_spin( &
          transitions                &
        , idxmap                     &
        , ispin                      &
        , xs_xcite                   &
        , xs_dxcite                  &
        , xs_xcite_spinavg           &
        , xs_dxcite_spinavg)

      if(allocated(prob_smat))   deallocate(prob_smat)

    enddo spinsdo

    ! -- write spin-averaged S-matrix cross sections to disk
    write(stdout, '(A)') "Writing  S-matrix cross sections ∀ spins to disk"

    smat_output_directory_all_spins = &
      & smat_output_directory // "all_spins" // DS

    do itrans = 1, size(transitions, 1)
      call write_smat_xs_to_file(              &
          "Smat"                          &
        , smat_output_directory_all_spins &
        , egrid_tot_smat                  &
        , transitions(itrans)             &
        , xs_xcite_spinavg(itrans)%vec    &
        , xs_dxcite_spinavg(itrans)%vec   &
        , G%LMAX_KMAT                   &
        )
    enddo

  end subroutine do_kmat_xs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine convert_multipoles(cartesian_moments_array, spherical_moments_array)
    !! Convert the supplied array of multipole moments from cartesian, obtained as typical output from
    !! quantum chemistry codes, to spherical multipole tensors
    use rotex__kinds,     only: dp
    use rotex__system,    only: die
    use rotex__constants, only: im, pi

    implicit none (type, external)

    real(dp), intent(in) :: cartesian_moments_array(:)
      !! Array containing cartesian multipole moments.
    complex(dp), intent(out) :: spherical_moments_array(:)
      !! Array containing spherical multipole moments.

    integer :: nelements
    integer :: ilambda_start
    real(dp) :: dx, dy, dz
    real(dp) :: Qxx, Qxy, Qxz, Qyy, Qyz, Qzz

    ilambda_start = 1

    nelements = size(cartesian_moments_array)

    select case(nelements)
    ! -- dipole
    case(3)
      dx = cartesian_moments_array(1)
      dy = cartesian_moments_array(2)
      dz = cartesian_moments_array(3)

      spherical_moments_array = [                            & !  μ
                                  (dx + im*dy)/sqrt(2.0_dp)  & ! -1
                                , cmplx(dz, kind = dp)       & !  0
                                , -(dx - im*dy)/sqrt(2.0_dp) & !  1
                                ]

    ! -- quadrupole
    case(6)

      ! Qxx = cartesian_moments_array(ilambda_start)
      ! Qxy = cartesian_moments_array(ilambda_start+1)
      ! Qxz = cartesian_moments_array(ilambda_start+2)
      ! Qyy = cartesian_moments_array(ilambda_start+3)
      ! Qyz = cartesian_moments_array(ilambda_start+4)
      ! Qzz = cartesian_moments_array(ilambda_start+5)

      ! spherical_moments_array = [                      & !  μ
      !          (Qxx - Qyy - 2*im*Qxy)/2                & ! -2
      !        , (Qxz - im*Qyz)/sqrt(2.0_dp)             & ! -1
      !        , (2*Qzz - Qxy - Qyy)/sqrt(6.0_dp) + 0*im & !  0
      !        ,-(Qxz + im*Qyz)/sqrt(2.0_dp)             & !  1
      !        , (Qxx - Qyy + 2*im*Qxy)/2                & !  2
      !        ]

    case default
      call die("Passed a cartesian multipole moment array of unacceptable size !")

    end select

  end subroutine convert_multipoles

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine make_grid(grid, E0)
    !! Make a segmented grid starting at E0
    !! Example with
    !!   G%NUM_EGRID_SEGS = 3, G%EGRID_SEGS = [1e-3, 1e-2, 1e-1, 1], G%NUM_EGRID = [1000,1000, 100]
    !! E0+1e-3           E0+1e-2              E0+1e-1           E0+1.0
    !!   !------------------!-------------------!- - - - - - - - -!
    !!      1000 energies     1000 energies      100 energies
    use rotex__system, only: die
    use rotex__arrays, only: realloc
    implicit none (type, external)
    real(dp), intent(inout), allocatable :: grid(:)
      !! The energy grid
    real(dp), intent(in) :: E0
      !! The lowest energy
    integer :: ie, iseg
    real(dp), allocatable :: dE(:)
    allocate(dE(G%NUM_EGRID_SEGS), source = 0.0_dp)
    select case(G%EGRID_SPACING)
    case("lin")
      dE = [(                                                                      &
          (G%EGRID_SEGS(iseg+1) - G%EGRID_SEGS(iseg))/(G%NUM_EGRID(iseg)-1) &
        , iseg=1, G%NUM_EGRID_SEGS                                                     &
      )]
      grid = [                                                                                               &
          E0                                                                                                 &
        , E0 + G%EGRID_SEGS(1)                                                                              &
        , ( (E0 + G%EGRID_SEGS(iseg) + ie*dE(iseg), ie=1, G%NUM_EGRID(iseg)-1), iseg=1, G%NUM_EGRID_SEGS ) &
      ]
    case("log")
      dE = [(                                                                                       &
          (G%EGRID_SEGS(iseg+1)/G%EGRID_SEGS(iseg))**(1/real(G%NUM_EGRID(iseg)-1,  kind=dp)) &
        , iseg=1, G%NUM_EGRID_SEGS                                                                      &
      )]
      grid = [                                                                                                &
          E0                                                                                                  &
        , E0 + G%EGRID_SEGS(1)                                                                               &
        , ( (E0 + G%EGRID_SEGS(iseg) * dE(iseg)**ie, ie=1, G%NUM_EGRID(iseg)-1 ), iseg=1, G%NUM_EGRID_SEGS ) &
      ]
    case default
      call die("EGRID_G%EGRID_SPACING (" // G%EGRID_SPACING // ") must be LIN or LOG")
    end select
  end subroutine make_grid

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine make_output_directories(pcb_output_directory, tcb_output_directory, smat_output_directory)
    use rotex__system,    only: DS => DIRECTORY_SEPARATOR, mkdir
    use rotex__constants, only: SPINMULT_NAMES
    implicit none (type, external)
    character(:), intent(out), allocatable :: pcb_output_directory, tcb_output_directory, smat_output_directory
    integer :: ispin
    call mkdir(G%OUTPUT_DIRECTORY)
    if(G%USE_CB .eqv. .true.) then
      pcb_output_directory = G%OUTPUT_DIRECTORY // "PCB" // DS
      tcb_output_directory = G%OUTPUT_DIRECTORY // "TCB" // DS
      call mkdir(PCB_output_directory)
      call mkdir(TCB_output_directory)
    endif
    if(G%USE_KMAT .eqv. .true.) then
      smat_output_directory = G%OUTPUT_DIRECTORY // "Smat" // DS
      call mkdir(smat_output_directory)
      do ispin = 1, size(G%SPINMULTS, 1)
        call mkdir(smat_output_directory // SPINMULT_NAMES(G%SPINMULTS(ispin))  // DS)
      enddo
    endif
  end subroutine make_output_directories

#ifdef USE_CDMSREADER
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_CDMS_data(filename, output_directory, CDMS_states, CDMS_transitions)
    !! Read a file from the CDMS search to get Einstein A coefficients that will be used in determining
    !! Coulomb-Born cross sections.
    use CDMSreader__readwrite, only: CDMS_readfile_nohfs
    use CDMSreader__types,     only: asymtop_state_type => asymtop_state_nohfs &
                                   , asymtop_transition_type => asymtop_transition_nohfs
    implicit none (type, external)
    character(*), intent(in) :: filename
    character(*), intent(in) :: output_directory
    integer :: funit_in, funit_out
    type(asymtop_state_type),       allocatable, intent(out) :: CDMS_states(:)
    type(asymtop_transition_type),  allocatable, intent(out) :: CDMS_transitions(:)
    open(newunit = funit_in,  file = filename)
    open(newunit = funit_out, file = output_directory // "CDMSreader.output.dat")
    call CDMS_readfile_nohfs( funit_in              &
                            , funit_out                   &
                            , states_nohfs = CDMS_states &
                            , transitions_nohfs = CDMS_transitions)
    close(funit_in)
    close(funit_out)
  end subroutine get_CDMS_data

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_cdms_state_energies(n_states, cdms_states)
    !! Update the energy of our rotational states with those from the CDMS
    use rotex__types,      only: n_states_type
    use rotex__constants,  only: au2invcm
    use cdmsreader__types, only: asymtop_state_type => asymtop_state
    implicit none (type, external)
    type(n_states_type),       intent(inout) :: n_states(:)
    class(asymtop_state_type), intent(in)    :: cdms_states(:)
    integer :: in, itau, icdms
    integer :: n, n_cdms
    integer :: ka, ka_cdms
    integer :: kc, kc_cdms
    real(dp) :: e_cdms
    do in = 1, size(n_states, 1)
      n = n_states(in) % n
      do itau = 1, 2*n+1
        ka = n_states(in) % ka(itau)
        kc = n_states(in) % kc(itau)
        do icdms = 1, size(cdms_states, 1)
          n_cdms = cdms_states(icdms) % dn / 2
          ka_cdms = cdms_states(icdms) % dka / 2
          kc_cdms = cdms_states(icdms) % dkc / 2
          if(n  .ne. n_cdms)  cycle
          if(ka .ne. ka_cdms) cycle
          if(kc .ne. kc_cdms) cycle
          ! n_states E <— CDMS E
          e_cdms      = cdms_states(icdms) % e / au2invcm ! -- CDMS energies in cm⁻¹
          n_states(in) % eigenh % eigvals(itau) = e_cdms
        enddo
      enddo
    enddo
  end subroutine get_CDMS_state_energies

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_CDMS_einstA(Nlo, Kalo, Kclo, Nup, Kaup, Kcup, CDMS_transitions, EinstA)
    !! Given the quantum numbers of a rotational transition, find the matching CDMS transition
    !! and get the corresponding Einstein A coefficient
    use rotex__system,     only: stderr, die
    use rotex__constants,  only: au2sec
    use CDMSreader__types, only: asymtop_transition_type       => asymtop_transition &
                               , asymtop_transition_nohfs_type => asymtop_transition_nohfs
    implicit none (type, external)
    integer,                        intent(in)  :: Nlo, Kalo, Kclo, Nup, Kaup, Kcup
    class(asymtop_transition_type), intent(in)  :: CDMS_transitions(:)
    real(dp),                       intent(out) :: EinstA
    integer :: i
    EinstA = 0
    select type(transitions => CDMS_transitions)
    type is(asymtop_transition_nohfs_type)
      do i=1, size(transitions, 1)
        if(Nlo  .ne. transitions(i) % lo % dN  / 2) cycle
        if(Nup  .ne. transitions(i) % up % dN  / 2) cycle
        if(Kalo .ne. transitions(i) % lo % dKa / 2) cycle
        if(Kclo .ne. transitions(i) % lo % dKc / 2) cycle
        if(Kaup .ne. transitions(i) % up % dKa / 2) cycle
        if(Kcup .ne. transitions(i) % up % dKc / 2) cycle
        ! -- EinstA coeffs are in 1/s
        EinstA = transitions(i) % EinstA * au2sec
        return
      enddo
    class default
      call die("CDMS_TRANSITIONS is not of type ASYMTOP_TRANSITION_NOHFS_TYPE")
    end select
  end subroutine get_CDMS_einstA
#endif

  ! -------------------------------------------------------------------------------------------------------------------------------
  module subroutine diagonalize_rotational_hamiltonian(num_n, n_values, n_states)
    !! Build the rigid-rotor hamiltonian for each N and diagonalize it. Keep eigenenergies and
    !! eigenvectors, stored in the eigenH type of n_states
    use rotex__kinds,      only: dp
    use rotex__types,      only: n_states_type, eigenh_type
    use rotex__system,     only: die
    use rotex__arrays,     only: size_check
    use rotex__hamilton,   only: h_asym, assign_projections, rotate_eigvecs
    use rotex__characters, only: lower
    implicit none (type, external)
    integer,             intent(in)  :: num_n, n_values(:)
    type(n_states_type), intent(out) :: n_states(:)
    integer  :: i_n, n, K
    real(dp) :: E_rot
    type(eigenh_type) :: hka, hkb, hkc
    complex(dp), allocatable :: eigvecs(:,:)
    character(1) :: current_axis
    call size_check(n_values, num_n, "N_VALUES")
    call size_check(n_states, num_n, "N_STATES")
    do i_n = 1, num_n

      n = n_values(i_n)
      n_states(i_n) % n = n

      select case(G%ROTOR_KIND)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!! asymmetric rotors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      case("a", "A")

        allocate(n_states(i_n) % einsta(2*n+1), source = 0.0_dp)

        associate(a => G%ABC(1), b => G%ABC(2), c => G%ABC(3))

          ! -- diagonalize in z=A frame so that we can use the CD coefficients and get Ka
          current_axis = "a"
          call asym_rigid_rotor(n, hka, b, c, a, G%CD4, G%CD6) ! <-- A basis, Ka = Kz, energies
          N_states(i_N) % eigenH = HKa
          eigvecs = HKa % eigvecs
          call assign_projections(N, eigvecs, N_states(i_N) % Ka) ! Ka labels

          ! -- diagonalize (without CD) in z=C basis to get Kc labels
          call asym_rigid_rotor(n, hkc, a, b, c) ! <-- C basis, Kc = Kz
          call assign_projections(N, hkc%eigvecs, N_states(i_N) % Kc) ! Kc labels

          ! -- rotate eigenvectors if needed to that xyz align with scattering calculations
          select case(G%ROTOR_ZAXIS)
          case("a","A")
            continue
          case("b","B","c","C")
            ! call rigid_rotor(n, hkb, c, a, b) ! <-- C basis, Kc = Kz
            ! N_states(i_N) % eigenH = hkb
            call rotate_eigvecs(N, current_axis, G%ROTOR_ZAXIS, eigvecs)
            N_states(i_N) % eigenH % eigvecs = eigvecs
          end select

        end associate

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!! symmetric rotors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      case("s", "S")

        allocate(n_states(i_n) % einsta(2*n+1), source=0.0_dp)

        ! -- appropriate energies whether z is aligned with A or C
        associate(a=>G%ABC(1), c=>G%ABC(3))
          select case(G%ROTOR_ZAXIS)
          case("a", "A")
            call sym_rigid_rotor(N, N_states(i_N)%eigenH, a, c, G%ADD_CD4, G%ADD_CD6, G%CD4, G%CD6)
            N_states(i_N)%Ka = [(K, K=-N, N)]
          case("c", "C")
            call sym_rigid_rotor(N, N_states(i_N)%eigenH, c, a, G%ADD_CD4, G%ADD_CD6, G%CD4, G%CD6)
            N_states(i_N)%Kc = [(K, K=-N, N)]
          case default
            call die("ZAXIS must be 'a' or 'c' for symmetric tops. Got "//G%ROTOR_ZAXIS)
          end select
        end associate

        ! -- if C₂ axis is the z-axis of the scattering calculations and is different than
        !    the rotational z-axis, rotate eigenvectors so that the z-axis lines up with the C₂ axis
        if(lower(G%ROTOR_C2AXIS) .eq. lower(G%ROTOR_ZAXIS)) cycle
        current_axis = G%ROTOR_ZAXIS
        call rotate_eigvecs(N, current_axis, G%ROTOR_C2AXIS, N_states(i_N)%eigenH%eigvecs)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!! linear rotors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      case("l")

        E_rot = G%B_ROT * (N*(N+1)) - G%D_ROT * (N*(N+1))**2 + G%H_ROT * (N*(N+1))**3
        allocate(n_states(i_n) % einsta(1),             source = 0.0_dp)
        allocate(n_states(i_n) % eigenH % eigvals(1),   source = E_rot )
        ! allocate(n_states(i_n) % Kc(1),               source = 0.0_dp)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      case default

        call die("Undefined value for namelist variable ROTOR_KIND: " // G%ROTOR_KIND)

      end select

    enddo

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  contains
  ! ------------------------------------------------------------------------------------------------------------------------------ !
    subroutine sym_rigid_rotor(nn, ham, bpara, bperp, add_cd4, add_cd6, cd4, cd6)
      !! Wrapper for calling the hamiltonian routine
      use rotex__hamilton, only: h_sym
      use rotex__types,    only: cd4_type, cd6_type
      implicit none (type, external)
      integer,           intent(in)  :: nn
      type(eigenh_type), intent(out) :: ham
      real(dp),          intent(in)  :: bpara,  bperp
      logical,           intent(in)  :: add_cd4, add_cd6
      type(cd4_type),    intent(in)  :: cd4
      type(cd6_type),    intent(in)  :: cd6
      integer :: test
      ! -- add cd4 ?
      test = merge(1, 0, add_cd4 .eqv. .true.)
      ! -- add cd4 & cd6 ?
      test = merge(2, test, (add_cd6 .eqv. .true.) .AND. (add_cd4 .eqv. .true.))
      select case(test)
        case(2) ; call h_sym(nn, ham, bpara, bperp, cd4, cd6)
        case(1) ; call h_sym(nn, ham, bpara, bperp, cd4)
        case(0) ; call h_sym(nn, ham, bpara, bperp)
        case default
          call die("Somehow got something other than 0,1,2 !")
      end select
    end subroutine sym_rigid_rotor
    subroutine asym_rigid_rotor(nn, ham, bx, by, bz, cd4, cd6)
      !! Wrapper for calling the hamiltonian routine
      use rotex__types, only: cd4_type, cd6_type
      implicit none (type, external)
      integer,           intent(in)  :: nn
      type(eigenh_type), intent(out) :: ham
      real(dp),          intent(in)  :: bx, by, bz
      type(cd4_type), intent(in), optional :: cd4
      type(cd6_type), intent(in), optional :: cd6
      logical :: add_cd4, add_cd6
      integer :: test
      add_cd4 = present(cd4)
      add_cd6 = present(cd6)
      ! -- add cd4 ?
      test = merge(1, 0, add_cd4 .eqv. .true.)
      ! -- add cd4 & cd6 ?
      test = merge(2, test, (add_cd6 .eqv. .true.) .AND. (add_cd4 .eqv. .true.))
      select case(test)
        case(2) ; call h_asym(nn, ham, bx, by, bz, cd4, cd6)
        case(1) ; call h_asym(nn, ham, bx, by, bz, cd4)
        case(0) ; call h_asym(nn, ham, bx, by, bz)
        case default
          call die("Somehow got something other than 0,1,2 !")
      end select
    end subroutine asym_rigid_rotor
  end subroutine diagonalize_rotational_hamiltonian

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine add_xs_contrib_this_spin( &
      transitions                           &
    , idxmap                                &
    , ispin                                 &
    , xs_xcite                              &
    , xs_dxcite                             &
    , xs_xcite_spinavg                      &
    , xs_dxcite_spinavg)
    !! This subrotine is called when calculating cross sections for a particular spin multiplicity.
    !! It adds the contribution of XS_(D)XCITE for this spin to the total XS_(D)XCITE_SPINAVG as
    !!   XS_SPINAVG += spinmul(ispin)*XS(ispin) / sum(spinmults)

    use rotex__types,  only: asymtop_rot_transition_type, rvector_type
    use rotex__utils,  only: assert
    use rotex__arrays, only: size_check

    implicit none (type, external)

    type(asymtop_rot_transition_type), intent(in) :: transitions(:)
      !! List of transitions between states lo and up
    integer, intent(in) :: idxmap(:)
      !! Mapping from XS_(D)XCITE to the arrays {TRANSITIONS, XS_(D)XCITE_SPINAVG}
    integer, intent(in) :: ispin
      !! The current index of spin multiplicities
    type(rvector_type), intent(in) :: xs_xcite(:), xs_dxcite(:)
      !! Cross sections for this spin multiplicity
    type(rvector_type), intent(inout) :: xs_dxcite_spinavg(:), xs_xcite_spinavg(:)

    integer :: itrans, ntrans
    real(dp) :: sum_spinmults

    ntrans = size(transitions, 1)
    call size_check(xs_xcite,          ntrans, "XS_XCITE")
    call size_check(xs_dxcite,         ntrans, "XS_DXCITE")
    call assert(size(xs_xcite_spinavg, 1)  .eq. ntrans, "XS_XCITE and XS_XCITE_SPINAVG do not have&
      & the same number of elements, which means that they do not have the same number of transitions.&
      & This should not be the case if only one electronic state was included, which is the only case that&
      & this code is currently expected to handle.")
    call assert(size(xs_dxcite_spinavg, 1) .eq. ntrans, "XS_DXCITE_SPINAVG and XS_DXCITE do not have&
      & the same number of elements, which means that they do not have the same number of transitions.&
      & However, this XS_XCITE and XS_DXCITE_SPINAVAG have the same number of elements, so something&
      & very unexpected has happened")

    sum_spinmults = real(sum(G%SPINMULTS), kind = dp)

    ! -- σ_allspins += + σ(ispin) * (2S+1)/Σ(2S+1)
    do itrans=1, ntrans
      xs_xcite_spinavg(itrans)%vec  = xs_xcite_spinavg(itrans)%vec  &
                                    + xs_xcite(idxmap(itrans))%vec  * G%SPINMULTS(ispin) / sum_spinmults
      xs_dxcite_spinavg(itrans)%vec = xs_dxcite_spinavg(itrans)%vec &
                                    + xs_dxcite(idxmap(itrans))%vec * G%SPINMULTS(ispin) / sum_spinmults
    enddo

  end subroutine add_xs_contrib_this_spin

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_xs_from_smat( &
        prob                          &
      , transitions                   &
      , egrid_tot                     &
      , xs_xcite                      &
      , xs_dxcite                     &
    )
    !! Calculate excitation and de-excitation cross sections from probabilities

    use rotex__kinds,       only: dp
    use rotex__types,       only: rvector_type, asymtop_rot_transition_type, asymtop_rot_channel_type
    use rotex__arrays,      only: size_check
    use rotex__symmetry,    only: rotstate_is_allowed
    use rotex__constants,   only: pi

    implicit none (type, external)

    type(rvector_type), intent(in) :: prob(:)
      !! Array of arrays of probabilities P
    type(asymtop_rot_transition_type), intent(inout), allocatable :: transitions(:)
      !! Array of transitions between states lo -> up
    real(dp), intent(in) :: egrid_tot(:)
      !! Array of total energies on which the calculations were performed
    type(rvector_type), intent(out), allocatable :: xs_xcite(:)
      !! Array of arrays of excitation cross sections σ ~ P/E
    type(rvector_type), intent(out), allocatable :: xs_dxcite(:)
      !! Array of arrays of de-excitation cross sections σ ~ P/E

    logical, allocatable :: keeptrans(:)
    integer :: ntrans, itrans, iemin, nlo, nup, ne, kaup, kcup, kalo, kclo, ksymlo, ksymup
    real(dp) :: Elo, Eup
    real(dp), allocatable :: Eel_ex(:), Eel_dex(:)
    type(asymtop_rot_channel_type) :: lo, up

    ntrans = size(prob, 1)
    call size_check(transitions, ntrans, "TRANSITIONS")

    ne = size(egrid_tot, 1)

    ! -- allocate cross sections to have the same size as egrid
    allocate(xs_xcite(ntrans), xs_dxcite(ntrans))
    do itrans = 1, ntrans
      allocate(xs_xcite(itrans)%vec(ne),  source=0.0_dp)
      allocate(xs_dxcite(itrans)%vec(ne), source=0.0_dp)
    enddo

    allocate(keeptrans(ntrans), source = .false.)

    do itrans = 1, ntrans

      ! -- get state info for this transition
      lo   = transitions(itrans) % lo
      up   = transitions(itrans) % up
      nlo  = lo % n
      kalo = lo%ka
      kclo = lo%kc
      nup  = up % n
      kaup = up%ka
      kcup = up%kc
      Elo  = lo % E
      Eup  = up % E

      if(rotstate_is_allowed(nlo, kalo, kclo) .eqv. .false.) cycle
      if(rotstate_is_allowed(nup, kaup, kcup) .eqv. .false.) cycle

      ! -- find the starting point of our energy grid
      iemin = findloc(prob(itrans) % vec(:) .gt. 0.0_dp, .true., 1)

      ! -- electron energy grids
      Eel_ex  = egrid_tot(iemin:) - Elo
      Eel_dex = egrid_tot(iemin:) - Eup

      ! -- cross sections
      xs_xcite(itrans)  % vec(iemin:) = prob(itrans)%vec(iemin:) * pi/(2*Eel_ex(:))  / (2*Nlo+1)
      xs_dxcite(itrans) % vec(iemin:) = prob(itrans)%vec(iemin:) * pi/(2*Eel_dex(:)) / (2*Nup+1)

      ! -- filter out tiny cross sections
      if( all(xs_xcite(itrans)  % vec(iemin:) .lt. G%XS_ZERO_THRESHOLD) ) cycle
      if( all(xs_dxcite(itrans) % vec(iemin:) .lt. G%XS_ZERO_THRESHOLD) ) cycle

      keeptrans(itrans) = .true.

    enddo

    xs_xcite    = pack(xs_xcite,    keeptrans)
    xs_dxcite   = pack(xs_dxcite,   keeptrans)
    transitions = pack(transitions, keeptrans)

    if(G%ROTOR_KIND .ne. "s") return

    ! -- reduce the resolution on all states to just have |K|
    if(G%symtop_reduce_projection) call reduce_symtop_ksign(transitions, xs_xcite, xs_dxcite)

  end subroutine get_xs_from_smat

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine combine_cb_smat_xs( &
      egrid_cb                                 &
    , egrid_tot_smat                           &
    , transitions_cb                           &
    , xs_xcite_pcb                             &
    , xs_xcite_tcb                             &
    , transitions_smat                         &
    , xs_xcite_smat                            &
    , xs_dxcite_smat                           &
    )
    !! Combine Coulomb-Born and S-matrix cross sections to be on the same electron energy grid.
    !! In general, a different number of transitions will exist for the Coulomb-Born cross sections
    !! and for the S-matrix cross sections. This routine matches the transitions, interpolates the
    !! CB cross sections to the total energy grid used by the S-matrix routines, and adds them
    !! together:
    !!   σ(tot) = σ(S-mat) + σ(TCB) - σ(PCB)

    use rotex__channel_ops, only: findloc_transitions
    use rotex__types,     only: rvector_type, art_type => asymtop_rot_transition_type, n_states_type
    use rotex__system,    only: stdout, stderr, die, DS => DIRECTORY_SEPARATOR
    use rotex__arrays,    only: append_uniq, size_check
    use rotex__splines,   only: interpolate_replace
    use rotex__constants, only: au2ev
    use rotex__writing,   only: write_cb_xs_to_file, write_smat_xs_to_file, write_total_xs_to_file

    implicit none (type, external)

    type(rvector_type), intent(in) :: egrid_cb(:)
      !! Array of arrays of electron energy grids for the CB cross sections for each transition
    real(dp),           intent(in) :: egrid_tot_smat(:)
      !! The total energy grid on which the S-matrix cross sections were evaluated
    type(art_type),     intent(in) :: transitions_cb(:)
      !! Array of arrays of transition between two rotational states lo and up in the CB approx
    type(rvector_type), intent(inout) :: xs_xcite_pcb(:)
      !! Array of arrays of partial CB cross sections for each transition (excitation only, de-excitation handled by detailed balance)
    type(rvector_type), intent(inout) :: xs_xcite_tcb(:)
      !! Array of arrays of total CB cross sections for each transition (excitation only, de-excitation handled by detailed balance)
    type(art_type),     intent(in) :: transitions_smat(:)
      !! Array of arrays of transition between two rotational states lo and up using the S-matrix
    type(rvector_type), intent(inout) :: xs_xcite_smat(:)
      !! Array of arrays of S-matrix cross sections for each transition; excitation
    type(rvector_type), intent(inout) :: xs_dxcite_smat(:)
      !! Array of arrays of S-matrix cross sections for each transition; de-excitation

    integer :: kalo, kclo, kaup, kcup
    integer :: itrans_cb, itrans_smat, itrans_all
    integer :: ntrans_cb, ntrans_smat, ntrans_all, nlo, nup
    integer, allocatable :: idx_all2smat(:), idx_all2cb(:)
      !! Mappings between various transition arrays
    integer, allocatable :: idx_smat(:)
      !! Indices for the S-matrix energy grid after interpolation (some energies might have been removed)

    logical :: warned_smat
    real(dp) :: Elo, Eup, dE, Eground
    real(dp), allocatable :: xs_pcb(:), xs_tcb(:), Egrid_xcite(:), Egrid_dxcite(:)
    real(dp), allocatable :: egrid_tot_cb(:)
    real(dp), allocatable :: xs_xcite_combined(:)
      !! Array of combined S-matrix + TCB - PCB cross sections; excitation
    real(dp), allocatable :: xs_dxcite_combined(:)
      !! Array of combined S-matrix + TCB - PCB cross sections; de-excitation

    character(:), allocatable :: total_xs_output_dir

    type(art_type), allocatable :: transitions_all(:)

    ! -- size checks on the arrays
    ntrans_cb = size(transitions_cb, 1)
    call size_check(egrid_cb,     ntrans_cb, "EGRID_CB")
    call size_check(xs_xcite_tcb, ntrans_cb, "XS_XCITE_TCB")
    call size_check(xs_xcite_pcb, ntrans_cb, "XS_XCITE_PCB")
    write(stdout, '("Detected ", I0, " Coulomb-Born excitations")') ntrans_cb
    ntrans_smat = size(transitions_smat, 1)
    call size_check(xs_xcite_smat,  ntrans_smat, "XS_XCITE_SMAT")
    call size_check(xs_dxcite_smat, ntrans_smat, "XS_DXCITE_SMAT")
    write(stdout, '("Detected ", I0, " S-matrix excitations")') ntrans_smat

    ! -- create a unified transitions array
    call append_uniq(transitions_all, transitions_smat)
    call append_uniq(transitions_all, transitions_cb)
    ntrans_all = size(transitions_all, 1)
    write(stdout, '("Detected ", I0, " unique excitations between the Coulomb-Born and S-matrix&
      & transitions")') ntrans_all
    write(stdout, *)

    ! -- mapping from total transitions to CB and S-matrix transitions
    idx_all2cb   = findloc_transitions(transitions_all, transitions_cb)
    idx_all2smat = findloc_transitions(transitions_all, transitions_smat)

    warned_smat = .false.

    Eground = egrid_tot_smat(1)

    ! -- make sure energies line up
    if( Eground .ne. minval(transitions_all(:)%lo%E) ) then
      write(stderr, '("EGROUND: ", E30.20)') Eground
      write(stderr, '("MINVAL(TRANSITIONS_ALL(:)%LO%E): ", E30.20)') minval(transitions_all(:)%lo%E)
      call die("Ground state energy mismatch ! Cannot align collision energy grids.")
    endif


    ! -- loop over all transitions
    transloop: do itrans_all = 1, ntrans_all

      itrans_smat = idx_all2smat(itrans_all)
      itrans_cb = idx_all2cb(itrans_all)

      total_xs_output_dir = G%OUTPUT_DIRECTORY // "Total" // DS

      if(itrans_cb .eq. 0) then

        if(itrans_smat .eq. 0) call die("Found a transition that is not in either the S-matrix&
          & or the CB transitions arrays ! Sounds like a bug.")

        ! -- no CB transition, but we do have an S-matrix. Print the corresponding S-matrix
        !    cross sections as the Total cross sections

        call write_smat_xs_to_file(         &
            "Total"                         &
          , total_xs_output_dir             &
          , egrid_tot_smat                  &
          , transitions_all(itrans_all)     &
          , xs_xcite_smat(itrans_smat)%vec  &
          , xs_dxcite_smat(itrans_smat)%vec &
          , G%LMAX_KMAT                     &
        )

        ! -- transitions are unique, so it should be safe to deallocate here
        deallocate(xs_xcite_smat(itrans_smat)%vec)
        deallocate(xs_dxcite_smat(itrans_smat)%vec)

        ! -- move on to next transition
        cycle transloop

      elseif(itrans_smat .eq. 0) then

        ! -- no S-matrix transition, but we do have an Coulomb-Born transition.
        if(warned_smat .eqv. .false.) then
          warned_smat = .true.
          write(stderr, '(A)') &
            "At least one transition has been detected that is available in the CB approx but not with&
            & the supplied S-matrix. This may or may not be an error."
          write(stderr, '(A)') &
            & "This should ONLY happen if XS_ZERO_THRESHOLD was set to some positive value, in which&
            & case the S-matrix amplitudes probably experience some destructive interference."
          write(stderr, '(A)') &
            & "Otherwise, the S-matrix should cover at least all dipole-allowed transitions."
        endif

        ! -- write CB excitation as Total excitation
        call write_CB_xs_to_file(             &
            "Total"                        &
          , total_xs_output_dir            &
          , egrid_cb(itrans_cb)%vec        &
          , xs_xcite_tcb(itrans_cb)%vec     &
          , transitions_all(itrans_all)%lo &
          , transitions_all(itrans_all)%up &
          , "inf"                          &
          , "TCB"                          &
        )

        ! -- convert σ excite  to σ de-excite
        nlo = transitions_all(itrans_all)%lo%n
        nup = transitions_all(itrans_all)%up%n
        Elo = transitions_all(itrans_all)%lo%E
        Eup = transitions_all(itrans_all)%up%E
        dE = Eup - Elo
        Egrid_dxcite = egrid_cb(itrans_cb)%vec - dE
        xs_tcb = xs_xcite_tcb(itrans_cb)%vec * egrid_cb(itrans_cb)%vec / Egrid_dxcite &
               * real(2*nlo+1, kind=dp) / real(2*nup+1, kind=dp)
        ! -- write CB de-excitation as Total de-excitation
        call write_CB_xs_to_file(             &
            "Total"                        &
          , total_xs_output_dir            &
          , egrid_cb(itrans_cb)%vec - dE    &
          , xs_tcb    &
          , transitions_all(itrans_all)%up &
          , transitions_all(itrans_all)%lo &
          , "inf"                          &
        )

        deallocate(xs_tcb)

        ! -- move to next transition
        cycle transloop

      endif

      ! -- At this part in the loop, we have both CB and S-matrix cross sections that
      !    need to be combined. Interpolate the CB cross sections (they have no resonances)
      !    to match the S-matrix energy grid

      nlo = transitions_all(itrans_all)%lo%n
      nup = transitions_all(itrans_all)%up%n
      Elo = transitions_all(itrans_all)%lo%E
      Eup = transitions_all(itrans_all)%up%E
      dE  = Eup - Elo

      ! -- CB grid is electron energy, convert to total energy and ensure
      !    that they start at the same energy relative to 0 (which may or may
      !    not be the ground state energy because of CDMS averaging)
      egrid_tot_cb = egrid_cb(itrans_cb) % vec + Elo

      ! -- PCB
      call interpolate_replace(       &
          egrid_tot_cb                &
        , egrid_tot_smat              &
        , xs_xcite_pcb(itrans_cb)%vec &
        , idx_smat                    &
      )
      ! -- TCB
      call interpolate_replace(       &
          egrid_tot_cb                &
        , egrid_tot_smat              &
        , xs_xcite_tcb(itrans_cb)%vec &
      )

      ! -- σTot = σSmat + σTCB - σPCB (excitation)
      xs_xcite_combined = xs_xcite_smat(itrans_smat)%vec(idx_smat(:)) &
                        + xs_xcite_tcb(itrans_cb)%vec(:)              &
                        - xs_xcite_pcb(itrans_cb)%vec(:)


      ! -- detailed balance, get de-excitation CB cross sections because we haven't stored
      !    those explicitly

      ! -- electron energy grid for de-excitation
      Egrid_xcite  = egrid_tot_smat(idx_smat(:)) - Elo
      Egrid_dxcite = egrid_tot_smat(idx_smat(:)) - Eup

      ! -- check for negative (de-)excitation energies
      if(any(Egrid_xcite .le. 0.0_dp)) then
        kaup = transitions_all(itrans_all)%up%Ka
        kcup = transitions_all(itrans_all)%up%Kc
        kalo = transitions_all(itrans_all)%lo%Ka
        kclo = transitions_all(itrans_all)%lo%Kc
        write(stderr, '("(N,Ka,Kc) -> (N,Ka,Kc): ", "(",2(I0,","),I0,")", " -> ", "(",2(I0,","),I0,")")') &
          nlo, kalo, kclo, nup, kaup, kcup
        write(stderr, '("Egrid_xcite(1) (eV): ", E30.20)') Egrid_xcite(1)*au2ev
        call die("Excitation electron energy grid has nonpositive&
        & energies, which means there was an issue converting the CB electron-energy grid to&
        & a total-energy grid")
      endif
      if(any(Egrid_dxcite .le. 0.0_dp)) then
        kaup = transitions_all(itrans_all)%up%Ka
        kcup = transitions_all(itrans_all)%up%Kc
        kalo = transitions_all(itrans_all)%lo%Ka
        kclo = transitions_all(itrans_all)%lo%Kc
        write(stderr, '("(N,Ka,Kc) -> (N,Ka,Kc): ", "(",2(I0,","),I0,")", " -> ", "(",2(I0,","),I0,")")') &
          nlo, kalo, kclo, nup, kaup, kcup
        write(stderr, '("Egrid_dxcite(1) (eV): ", E30.20)') Egrid_dxcite(1)*au2ev
        call die("De-excitation electron energy grid has nonpositive&
        & energies, which means there was an issue converting the CB electron-energy grid to&
        & a total-energy grid")
      endif

      ! -- σ(E1) * E1 = σ(E2) * E2
      xs_pcb = xs_xcite_pcb(itrans_cb)%vec(:)           &
             * Egrid_xcite(:)         / Egrid_dxcite(:) &
             * real(2*nlo+1, kind=dp) / real(2*nup+1, kind=dp)
      xs_tcb = xs_xcite_tcb(itrans_cb)%vec(:)           &
             * Egrid_xcite(:)         / Egrid_dxcite(:) &
             * real(2*nlo+1, kind=dp) / real(2*nup+1, kind=dp)

      ! -- σTot = σSmat + σTCB - σPCB (de-excitation)
      xs_dxcite_combined = xs_dxcite_smat(itrans_smat)%vec(idx_smat(:)) &
                         + xs_tcb(:)                                    &
                         - xs_pcb(:)

      ! -- no longer needed, deallocate
      deallocate(xs_xcite_smat(itrans_smat)%vec)
      deallocate(xs_xcite_pcb(itrans_cb)%vec)
      deallocate(xs_xcite_tcb(itrans_cb)%vec)
      deallocate(xs_dxcite_smat(itrans_smat)%vec)

      ! -- Write excitation to disk
      call write_total_xs_to_file(            &
          "Total"                             &
        , total_xs_output_dir                 &
        , Egrid_xcite                         &
        , Egrid_dxcite                        &
        , transitions_all(itrans_all)         &
        , xs_xcite_combined                   &
        , xs_dxcite_combined                  &
        , G%LMAX_KMAT                       &
      )

    enddo transloop

  end subroutine combine_cb_smat_xs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine convert_xcite2dxcite(Eel_xcite, xs, transition)
    !! Convert XS from excitation cross sections to de-excitation cross sections
    use rotex__kinds, only: dp
    use rotex__types, only: asymtop_rot_transition_type
    use rotex__system, only: die
    implicit none (type, external)
    real(dp), intent(in) :: Eel_xcite(:)
    real(dp), intent(inout) :: xs(:)
    type(asymtop_rot_transition_type), intent(in) :: transition
    integer :: Nlo, Nup
    real(dp) :: Elo, Eup, dE
    real(dp), allocatable :: Eel_dxcite(:)
    Nlo = transition % lo % N
    Nup = transition % up % N
    Elo = transition % lo % E
    Eup = transition % up % E
    dE = Eup - Elo
    Eel_dxcite = Eel_xcite(:) - dE
    if(any(Eel_dxcite .le. 0)) call die("Non-positive electron energies when calculating&
      & de-excitation cross sections from excitation cross sections ! The excitation electron&
      & energy grid has values that are below threshold, which is not expected behavior !")
    xs = xs * Eel_xcite / Eel_dxcite * real(2*nlo+1, kind=dp) / real(2*nup+1, kind=dp)
  end subroutine convert_xcite2dxcite

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure subroutine reduce_symtop_ksign(transitions, xs_xcite, xs_dxcite)
    !! Given arrays of transitions between states, excitation cross sections, and de-excitation cross sections,
    !! average over the different ±K for each state, .e.g,
    !!  (N,K)
    !!    (1-1) -> (2-1)    +> (1 1) -> (2 1)
    !!    (1-1) -> (2 1)   /
    !!    (1 1) -> (2-1)  /
    !!    (1 1) -> (2 1) /

    use rotex__types,  only: asymtop_rot_transition_type, rvector_type, asymtop_rot_channel_type
    use rotex__arrays, only: size_check
    use rotex__channel_ops, only: findloc_transitions, operator(.ne.)

    implicit none (type, external)

    type(asymtop_rot_transition_type), intent(inout), allocatable :: transitions(:)
      !! Array of state transitions
    type(rvector_type), intent(inout), allocatable :: xs_xcite(:), xs_dxcite(:)
      !! Excitation and de-excitation cross section arrays

    integer :: n, itrans
    integer :: kalo, kclo
    integer :: kaup, kcup
    integer, allocatable :: idx(:)
    logical, allocatable :: mask(:)
    type(asymtop_rot_transition_type) :: transabs
    type(asymtop_rot_channel_type) :: lo, up, loabs, upabs

    n = size(transitions, 1)
    call size_check(xs_xcite,  n, "XS_XCITE")
    call size_check(xs_dxcite, n, "XS_DXCITE")

    allocate(mask(n), source = .true.)
    idx = [(itrans, itrans=1, n)]

    ! -- initial marking
    do itrans=1,n

      ! -- N, Ka, Kc
      up   = transitions(itrans)%up
      Kaup = up % Ka
      Kcup = up % Kc
      lo   = transitions(itrans)%lo
      Kalo = lo % Ka
      Kclo = lo % Kc

      ! -- N, |Ka|, |Kc|
      upabs = up
      upabs%Ka = abs(upabs%Ka)
      upabs%Kc = abs(upabs%Kc)
      loabs = lo
      loabs%Ka = abs(loabs%Ka)
      loabs%Kc = abs(loabs%Kc)

      transabs%lo = loabs
      transabs%up = upabs

      xs_xcite(itrans)%vec  = xs_xcite(itrans)%vec  / symtop_degen(Kalo, Kclo)
      xs_dxcite(itrans)%vec = xs_dxcite(itrans)%vec / symtop_degen(Kaup, Kcup)

      ! -- is NKaKc -> N'Ka'Kc' = N|Ka||Kc| -> N|Ka'||Kc'| ?
      idx(itrans) = findloc_transitions(transabs, transitions)
      if(idx(itrans) .eq. itrans) cycle

      ! -- if not, remove it and add its contribution to the corresponding transition
      mask(itrans) = .false.

      xs_xcite(idx(itrans))%vec  = xs_xcite(idx(itrans))%vec  + xs_xcite(itrans)%vec
      xs_dxcite(idx(itrans))%vec = xs_dxcite(idx(itrans))%vec + xs_dxcite(itrans)%vec

    enddo

    ! -- filter out transitions with K<0
    transitions = pack(transitions, mask)
    xs_xcite    = pack(xs_xcite, mask)
    xs_dxcite   = pack(xs_dxcite, mask)

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  contains
  ! ------------------------------------------------------------------------------------------------------------------------------- !

    pure elemental function symtop_degen(ka, kc) result(res)
      !! Returns the degeneracy of a transition:
      !! K=0: 1
      !! K≠0: 2
      use rotex__globals, only: G
      use rotex__system, only: die
      implicit none (type, external)
      integer, intent(in) :: Ka, Kc
      integer :: res
      integer :: Ksym
      if(G%ROTOR_KIND .ne. "s") call die("SYMTOP_DEGEN found a non-symtop rotor")
      if(all(G%ROTOR_ZAXIS .ne. ["a", "c"])) call die("SYMTOP_DEGEN needs ROTOR_ZAXIS to be A or C")
      Ksym = merge(Ka, Kc, G%ROTOR_ZAXIS .eq. "a")
      res = merge(1, 2, Ksym .eq. 0)
    end function symtop_degen

  end subroutine reduce_symtop_ksign

! ================================================================================================================================ !
end module rotex__drivers
! ================================================================================================================================ !
