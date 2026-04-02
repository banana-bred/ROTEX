! ================================================================================================================================ !
module rotex__RFT
  !! Procedures used to carry out the rotational frame transformation
  use rotex__globals, only: G
  use rotex__kinds, only: dp
  use rotex__system, only: stderr, die, stdout

  implicit none (type, external)

  private

  ! public :: RFT_linear
  public :: RFT_nonlinear

  interface real2complex_ylm
    module procedure :: real2complex_ylm_r
    module procedure :: real2complex_ylm_c
  end interface real2complex_ylm

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine RFT_nonlinear( Kmat                     &
                                 , kmat_eval_energies       &
                                 , spinmult                 &
                                 , Jmin, Jmax               &
                                 , Smat_J                   &
                                 , elec_channels            &
                                 , N_states                 &
                                 , asymtop_rot_channels_l   &
                                 , asymtop_rot_channels_l_J &
                                 )
    !! Build the electronic S-matrix from the electronic K-matrix, then perform the rotational frame transformation on the S-matrix

    use rotex__kinds,     only: dp
    use rotex__types,     only: elec_channel_type, asymtop_rot_channel_l_type, N_states_type &
                              , asymtop_rot_channel_l_vector_type, r3carr_type
    use rotex__channel_ops, only: sort_channels_by_energy
    use rotex__system,    only: stdout, die
    use rotex__symmetry,  only: possible_spin_symmetries, spin_symmetry
    use rotex__constants, only: im
    use rotex__arrays,    only: append, is_unitary
    use rotex__writing,   only: write_Smat_J_elems_to_file

    implicit none (type, external)

    real(dp), intent(inout), allocatable :: Kmat(:,:,:)
      !! The K-matrix, needed as input for the RFT. nchans × nchans × nE
    real(dp), intent(in) :: kmat_eval_energies(:)
      !! The evaluation energies of the K-matrix
    integer, intent(in) :: spinmult
      !! Current spin multiplicity
    integer, intent(in) :: Jmin
      !!  The lowest value of J = N + l
    integer, intent(in) :: Jmax
      !!  The largest value of J = N + l
    type(r3carr_type), intent(out) :: Smat_J(Jmin:Jmax)
      !! The rotational S-matrices \(S^J\), produced by the RFT, for each J
    type(elec_channel_type),        intent(in)                 :: elec_channels(:)
      !! The array of electronic channels (n, l, ml), needed as input for the RFT
    type(N_states_type),            intent(in)                 :: N_states(:)
      !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
    type(asymtop_rot_channel_l_type), intent(out),   allocatable :: asymtop_rot_channels_l(:)
      !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix
    type(asymtop_rot_channel_l_vector_type), intent(out) :: asymtop_rot_channels_l_J(Jmin:Jmax)
      !! The array of arrays of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblocks at each J

    integer :: nelec, l, ie, ml
    integer :: ichan, nchans_rot, nchans_elec, ne
    real(dp), allocatable :: sin_elec(:,:,:), cos_elec(:,:,:)
      !! nchan × nchan × ne arrays of electronic sine/cosine matrix elements used in the EDFT
    complex(dp), allocatable :: csin_elec(:,:,:), ccos_elec(:,:,:)
      !! nchan × nchan × ne arrays of electronic sine/cosine matrix elements used in the EDFT, complex
    complex(dp), allocatable :: smat_elec(:,:)
      !! nchan × nchan array of S-matrix elements. Not used in EDFT because
      !! we do that on the sine and cosine matrices

    ne = size(kmat_eval_energies, 1)
    nchans_elec = size(elec_channels, 1)

    call build_rotational_channels( &
        n_states                    &
      , elec_channels               &
      , asymtop_rot_channels_l      &
    )

    ! -- channels are (probably) sorted by contsruction, but we should still sort them by energy here in case
    !    anything above changes in the future
    call sort_channels_by_energy(asymtop_rot_channels_l)

    nchans_rot = size(asymtop_rot_channels_l, 1)
    if(nchans_rot .lt. 1) call die("Number of rotational channels must not be less than 1 !")

    allocate(sin_elec(nchans_elec, nchans_elec, ne), source=0._dp)
    allocate(cos_elec(nchans_elec, nchans_elec, ne), source=0._dp)

    if(G%EDFT) then

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! ENERGY DEPENDENT FRAME TRANSFORMATION !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! -- correct for branch cuts to get a smooth S-matrix
      write(stdout, '("K -> sin, cos ... ")', advance="no")
      call K2sincos(kmat, kmat_eval_energies, sin_elec, cos_elec, elec_channels, spinmult)
      write(stdout, '("done !")')

      csin_elec = cmplx(sin_elec, kind=dp) ; deallocate(sin_elec)
      ccos_elec = cmplx(cos_elec, kind=dp) ; deallocate(cos_elec)

      if(G%REAL_SPHERICAL_HARMONICS) then
        write(stdout, '("X_{lλ}(θ,φ) -> Y_l^λ(θ,φ)... ")', advance="no")
        do concurrent(ie=1:ne)
          call real2complex_ylm(csin_elec(:,:,ie), elec_channels)
          call real2complex_ylm(ccos_elec(:,:,ie), elec_channels)
        enddo
        write(stdout, '("done !")')
      endif

      write(stdout, '(A)') "Frame transformation: sin_elec, cos_elec -> S^J"
      call do_edrft(               &
          csin_elec                &
        , ccos_elec                &
        , kmat_eval_energies       &
        , smat_J                   &
        , Jmin                     &
        , Jmax                     &
        , N_states                 &
        , elec_channels            &
        , asymtop_rot_channels_l   &
        , asymtop_rot_channels_l_J &
      )

      ! -- write S^J(e) elements and channels to file
      call write_Smat_J_elems_to_file(spinmult, Jmin, Jmin, kmat_eval_energies, asymtop_rot_channels_l_J, smat_J)
      ! call write_Smat_J_elems_to_file(Jmin, Jmax, kmat_eval_energies, asymtop_rot_channels_l_J, smat_J)

    else

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! ENERGY INDEPENDENT FRAME TRANSFORMATION !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! -- just do the Cayley transform and get the S-matrix directly
      allocate(smat_elec, mold=cmplx(kmat(:,:,1), kind=dp)) ; smat_elec = 0.0_dp
      call K2S_cayley(kmat(:,:,1), smat_elec)

      ! -- make sure we're in the basis of COMPLEX spherical harmonics
      if(G%REAL_SPHERICAL_HARMONICS) call real2complex_ylm(smat_elec, elec_channels)

      write(stdout, '(A)') "Frame transformation: S_elec -> S^J"

      call do_eirft(               &
          smat_elec                &
        , smat_J                   &
        , Jmin                     &
        , Jmax                     &
        , N_states                 &
        , elec_channels            &
        , asymtop_rot_channels_l   &
        , asymtop_rot_channels_l_J &
      )

    endif

  end subroutine RFT_nonlinear

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine build_U_RFT( &
        J                      &
      , Omega                  &
      , N_states               &
      , elec_channels          &
      , rot_channels           &
      , U                      &
    )
    !! Build the transformation matrix U for the RFT for a given J, Ω

    use rotex__types,     only: N_states_type, elec_channel_type, asymtop_rot_channel_l_type
    use rotex__arrays,    only: size_check, adjoint
    use rotex__wigner,    only: clebsch
    use rotex__functions, only: neg

    implicit none(type, external)

    integer,                          intent(in)  :: J, Omega
    type(N_states_type),              intent(in)  :: N_states(:)
    type(elec_channel_type),          intent(in)  :: elec_channels(:)
    type(asymtop_rot_channel_l_type), intent(in)  :: rot_channels(:)
    complex(dp),                      intent(out) :: U(:,:)

    integer :: nrot, nelec
    integer :: i_N, ik, itau
    integer :: Ni, Kai, Kci, li, neleci, irot, ielec, Ki
    integer :: nelecj, lj, lambdaj, jelec
    logical, allocatable :: mask(:)

    nrot  = size(rot_channels, 1)
    nelec = size(elec_channels, 1)

    call size_check(U, [nrot, nelec], "U_RFT")

    U = (0.0_dp, 0.0_dp)

    do irot=1, nrot
      call get_channel_qnums_rot(rot_channels, irot, neleci, Ni, Kai, Kci, li)

      i_N = findloc(N_states % N, value = Ni, dim = 1)
      if(i_N .eq. 0) cycle

      select case(G%ROTOR_KIND)
      case("a", "A")
        mask = (n_states(i_N) % ka(:) .eq. kai) .AND. (n_states(i_N) % kc(:) .eq. kci)
      case("s", "S")
        select case(G%ROTOR_ZAXIS)
        case("a", "A")
          mask = (n_states(i_N) % ka(:) .eq. kai)
        case("c", "C")
          mask = (n_states(i_N) % kc(:) .eq. kci)
        case default
          call die("Somehow got a symmetric top with a G%ROTOR_ZAXIS " // G%ROTOR_ZAXIS // " that is neither A nor C")
        end select
      case default
        call die("ROTOR_KIND " // G%ROTOR_KIND // " not allowed in RFT")
      end select

      itau  = findloc(mask, value = .true., dim = 1)

      do jelec = 1, nelec

        nelecj  = elec_channels(jelec) % nelec
        lj      = elec_channels(jelec) % l
        lambdaj = elec_channels(jelec) % ml

        ! -- enforce transformation between the same electronic state n and partial wave l
        if(neleci .ne. nelecj) cycle
        if(li .ne. lj) cycle
        Ki = Omega - lambdaj
        if(abs(Ki) .gt. Ni) cycle

        ik = Ki + Ni + 1

        U(irot, jelec) = neg(lj + lambdaj)                          &
                       * N_states(i_N) % eigenH % eigvecs(ik, itau) &
                       * clebsch(lj, -lambdaj, J, Omega, Ni, Ki)

      enddo
    enddo

  end subroutine build_U_RFT

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine do_eirft(           &
        smat_elec                &
      , smat_J                   &
      , Jmin                     &
      , Jmax                     &
      , N_states                 &
      , elec_channels            &
      , asymtop_rot_channels_l   &
      , asymtop_rot_channels_l_J &
    )
    !! Perform the Energy Independent Rotational Frame Transformation on the electronic S-matrix

    use rotex__types,  only: r3carr_type, n_states_type, asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type&
                           , elec_channel_type
    use rotex__arrays, only: realloc, is_symmetric, is_unitary, adjoint, unitary_defect

    implicit none (type, external)

    complex(dp),                             intent(in) :: smat_elec(:,:)
      !! Electronic S-matrix: n×n
    type(r3carr_type),                      intent(out):: smat_J(Jmin:Jmax)
      !! Array of S^J sub blocks
    integer,                                 intent(in) :: Jmin, Jmax
      !! Minimum and maximum value of J to consider
    type(n_states_type),                     intent(in) :: N_states(:)
      !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
    type(elec_channel_type),                 intent(in) :: elec_channels(:)
      !! The array of electronic channels (n, l, ml), needed as input for the RFT
    type(asymtop_rot_channel_l_type),        intent(in) :: asymtop_rot_channels_l(:)
      !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix
    type(asymtop_rot_channel_l_vector_type), intent(out):: asymtop_rot_channels_l_J(Jmin:Jmax)
      !! The array of arrays of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblocks at each J

    logical :: bad_unitarity
    integer :: J, Omega, nchans_elec, nchans_rot, ne
    integer, allocatable :: idx(:)
    type(asymtop_rot_channel_l_type), allocatable :: rot_channels(:)
    complex(dp), allocatable :: U(:,:), smat_rot(:,:)
    complex(dp), allocatable :: C(:,:)

    bad_unitarity = .false.
    nchans_elec = size(elec_channels, 1)

    do J = Jmin, Jmax

      call collect_j_channels_indices(J, asymtop_rot_channels_l, idx)
      asymtop_rot_channels_l_J(J) % channels = asymtop_rot_channels_l(idx)
      rot_channels = asymtop_rot_channels_l_J(J) % channels

      nchans_rot = size(rot_channels, 1)

      call realloc(U,        nchans_rot, nchans_elec)
      call realloc(smat_rot, nchans_rot, nchans_rot)
      smat_rot = (0.0_dp, 0.0_dp)

      do Omega = -J, J
        call build_U_RFT(J, Omega, N_states, elec_channels, rot_channels, U)
        ! smat_rot = smat_rot + matmul(U, matmul(smat_elec, adjoint(U)))
        smat_rot = smat_rot + matmul(U, matmul(smat_elec, adjoint(U)))
        ! C = C + matmul(u, adjoint(u))
      enddo

      if(.not. allocated(smat_J(J) % arr)) allocate(smat_J(J) % arr(nchans_rot, nchans_rot, 1))
      smat_J(J) % arr(:,:,1) = smat_rot

      if(.not. is_symmetric(smat_rot)) call die("The S-matrix is not symmetric after the energy-independent RFT !")
      if(is_unitary(smat_rot)) cycle

      bad_unitarity = .true.
      write(stderr, '("WARN: The S-matrix for J = ", I0, " is nonunitary with unitary defect ", ES12.4)') &
        J, unitary_defect(smat_rot)

    enddo

    if(bad_unitarity .eqv. .false.) return

    call die("At least one J-block of the S-matrix is non-uitary after the energy-independent RFT. Aborting !")

  end subroutine do_eirft

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine do_edrft(           &
        sin_elec                 &
      , cos_elec                 &
      , kmat_eval_energies       &
      , smat_J                   &
      , Jmin                     &
      , Jmax                     &
      , N_states                 &
      , elec_channels            &
      , asymtop_rot_channels_l   &
      , asymtop_rot_channels_l_J &
    )
    !! Perform the Energy Dependent Rotational Frame Transformation on the electronic Sine and Cosine matrices

    use rotex__types,  only: r3carr_type, N_states_type, elec_channel_type, &
                             asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type
    use rotex__arrays, only: realloc, adjoint, is_unitary, unitary_defect, is_symmetric &
                           , interp_matrix_at_energy

    implicit none (type, external)

    complex(dp),                              intent(in)    :: sin_elec(:,:,:), cos_elec(:,:,:)
      !! Electronic sine and cosine matrices: n×n
    real(dp),                                 intent(in)    :: Kmat_eval_energies(:)
      !! Array of K-matrix evaluation energies (therefore the evaluation energies of the other
      !! electronic matrices as well)
    type(r3carr_type),                       intent(inout) :: smat_J(Jmin:Jmax)
      !! Array of S^J sub blocks
    integer,                                  intent(in)    :: Jmin, Jmax
      !! Smallest and largest values of J to consider
    type(N_states_type),                      intent(in)    :: N_states(:)
      !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
    type(elec_channel_type),                  intent(in)    :: elec_channels(:)
      !! The array of electronic channels (n, l, ml), needed as input for the RFT
    type(asymtop_rot_channel_l_type),         intent(in)    :: asymtop_rot_channels_l(:)
      !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix
    type(asymtop_rot_channel_l_vector_type),  intent(out)    :: asymtop_rot_channels_l_J(Jmin:Jmax)
      !! The array of arrays of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblocks at each J

    logical :: bad_unitarity, bad_symmetry, did_symmetrize
    integer :: J, Omega, ie, jrot, ne, nchans_elec, nchans_rot
    integer, allocatable :: idx(:)
    real(dp) :: E_elec_rhs
    type(asymtop_rot_channel_l_type), allocatable :: rot_channels(:)
    complex(dp), allocatable :: U(:,:), Udagg(:,:), smat_rot(:,:), sin_rot(:,:), cos_rot(:,:)
    complex(dp), allocatable :: sinE(:,:), cosE(:,:), tmp(:)

    bad_unitarity = .false.
    bad_symmetry = .false.
    ne            = size(kmat_eval_energies, 1)
    nchans_elec   = size(elec_channels, 1)

    write(stdout, '(2X, 2(A5, " /"), A5)') "Jmin", "J", "Jmax"
    Jloop: do J = Jmin, Jmax

      write(stdout, '(2X, 2(I5," /"),I5,X)') Jmin, J, Jmax
      call collect_j_channels_indices(J, asymtop_rot_channels_l, idx)
      asymtop_rot_channels_l_J(J) % channels = asymtop_rot_channels_l(idx)
      rot_channels = asymtop_rot_channels_l_J(J) % channels

      nchans_rot = size(rot_channels, 1)

      call realloc(U,        nchans_rot,  nchans_elec)
      call realloc(Udagg,    nchans_elec, nchans_rot)
      call realloc(smat_rot, nchans_rot,  nchans_rot)
      call realloc(sin_rot,  nchans_rot,  nchans_rot)
      call realloc(cos_rot,  nchans_rot,  nchans_rot)
      call realloc(sinE,     nchans_elec, nchans_elec)
      call realloc(cosE,     nchans_elec, nchans_elec)
      call realloc(tmp,      nchans_elec)

      if(.not. allocated(smat_J(J) % arr)) allocate(smat_J(J) % arr(nchans_rot, nchans_rot, ne))

      did_symmetrize = .false.

      !TODO: OMP maybe ?
      nrg: do ie=1, ne

        smat_rot = (0.0_dp, 0.0_dp)
        sin_rot  = (0.0_dp, 0.0_dp)
        cos_rot  = (0.0_dp, 0.0_dp)

        do Omega = -J, J

          call build_U_RFT(J, Omega, N_states, elec_channels, rot_channels, U)
          Udagg = adjoint(U)

          ! -- loop over RHS channels because we need to interpolate/pick matrices at a specific energy
          do jrot = 1, nchans_rot

            ! -- take the evaluation energy of the right-hand channel, and interpolate the
            !    matrices to get the corresponding value
            E_elec_rhs = kmat_eval_energies(ie) - rot_channels(jrot) % E
            call interp_matrix_at_energy(E_elec_rhs, kmat_eval_energies, sin_elec, sinE)
            call interp_matrix_at_energy(E_elec_rhs, kmat_eval_energies, cos_elec, cosE)

            ! -- U {sin(E),cos(E)} U⁺
            sin_rot(:, jrot) = sin_rot(:, jrot) + matmul( U, matmul(sinE, Udagg(:, jrot)) )
            cos_rot(:, jrot) = cos_rot(:, jrot) + matmul( U, matmul(cosE, Udagg(:, jrot)) )
          enddo

        enddo

        call sincos2S(sin_rot, cos_rot, smat_rot, did_symmetrize)
        smat_J(J) % arr(:,:,ie) = smat_rot

        if(is_unitary(smat_rot) .eqv. .false.) then
          bad_unitarity = .true.
          write(stderr, '("WARN: The EDFT S-matrix for J= ", I0, ", ie = ", I0, ", is notunitary with&
            &  unitary defect ", ES12.4)') J, ie, unitary_defect(smat_rot)
        endif

        if(is_symmetric(smat_rot) .eqv. .false.) then
          bad_symmetry = .true.
          write(stderr, '("WARN: The EDFT S-matrix for J = ", I0, ", ie = ", I0, " is not symmetric")') J, ie
        endif

      enddo nrg

      if(did_symmetrize) write(stderr, '("INFO: intermediate K was detected as non-symmetric and was symmetrized (J = ",I0,")")') J

    enddo Jloop

    if(bad_unitarity .OR. bad_symmetry) call die("At least one J-block of the EDFT S-matrix is&
      & non-unitary or non-symmetric after asymmetrization !")

  end subroutine do_edrft

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine build_rotational_channels(n_states, elec_channels, rot_channels)
    !! Build rotational+electronic channels channels: (N Ka Kc)+(l λ) = (N Ka Kc l λ)
    use rotex__kinds,    only: dp
    use rotex__types,    only: n_states_type, elec_channel_type, asymtop_rot_channel_l_type
    use rotex__channel_ops, only: operator(.eq.)
    use rotex__arrays,   only: append
    use rotex__symmetry, only: spin_symmetry
    use rotex__system,   only: die
    implicit none (type, external)
    type(N_states_type),            intent(in)               :: n_states(:)
    type(elec_channel_type),        intent(in)               :: elec_channels(:)
    type(asymtop_rot_channel_l_type), intent(out), allocatable :: rot_channels(:)
    integer  :: i_N_state, i_tau, i_elec_channel
    integer  :: n, ka, kc, ksym, nelec, iq, sym
    integer  :: l
    real(dp) :: e, e_elec, e_rot
    type(asymtop_rot_channel_l_type) :: channel
    ka = 0; kc = 0
    ! -- build rotational channels from elec_channels and N_states
    do i_n_state = 1, size(n_states, 1)

      n  = n_states(i_n_state) % n

      do i_tau = 1, 2*n + 1

        select case(G%ROTOR_KIND)
        case("a", "A")
          ka = n_states(i_n_state) % ka(i_tau)
          kc = n_states(i_n_state) % kc(i_tau)
        case("s", "S")
          select case(G%ROTOR_ZAXIS)
          case("a", "A")
            ksym = n_states(i_n_state) % ka(i_tau)
            ka = ksym
            kc = 0
          case("c", "C")
            ksym = n_states(i_n_state) % kc(i_tau)
            ka = 0
            kc = ksym
          case default
            call die("Symtop rotational G%ROTOR_ZAXIS must be A or C")
          end select

          ! ! -- skip forbidden channels
          ! if(symtop_rotstate_is_allowed(N, Ksym) .eqv. .false.) cycle

        end select
        do i_elec_channel = 1, size(elec_channels, 1)
          nelec = elec_channels(i_elec_channel) % nelec
          l     = elec_channels(i_elec_channel) % l
          iq    = elec_channels(i_elec_channel) % iq
          ! -- get the channel energy (rotational + electronic)
          e_elec = elec_channels(i_elec_channel) % e
          e_rot  = n_states(i_n_state) % eigenh % eigvals(i_tau)
          e      = e_rot + e_elec
          ! -- for now, enforce ground state RE only
          if(nelec .ne. 1) cycle
          sym = spin_symmetry(n, ka, kc)
          channel = asymtop_rot_channel_l_type(nelec=nelec, l=l, iq=iq, n=n, ka=ka, kc=kc, e=e, sym=sym)
          if(allocated(rot_channels)) then
            if(any(channel .eq. rot_channels)) cycle
          endif
          call append(rot_channels, channel)
        enddo
      enddo
    enddo
  end subroutine build_rotational_channels

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function get_eigenvector_permutation_idx(U1, U0, eigenphases1, eigenphases0) result(res)
    !! For eigenvectors U1 and U0 computed at energy E1 and E0 with eigenphases EIGENPHASES1 and EIGENPHASES0,
    !! return the index permutation of matrices 0 such that their inner products with 1 are the best.
    !! This is essentially trying to find the "best" permutation of eigenvectors to get the "same"
    !! eigenphases in the same order across energies
    implicit none (type, external)
    real(dp), intent(in) :: U1(:,:)
      !! Eigenvectors for ie-1
    real(dp), intent(in) :: U0(:,:)
      !! Eigenvectors for ie
    real(dp), intent(in) :: eigenphases1(:)
      !! Eigenphases for ie-1
    real(dp), intent(in) :: eigenphases0(:)
      !! Eigenphases for ie

    integer, allocatable :: res(:)

    integer :: i, n
    integer, allocatable :: idx(:)
    real(dp), allocatable :: inner(:)

    n = size(U1, 1)

    ! -- NOTE right here we could check the signs of the vectors, but we don't have to do that
    !    when we're looking at energies because they come from the same UKRmol+ calculation.
    !    We would use the sine and cosine matrices, but for now we don't need to do this
    !    call CHECK_TARGET_SIGN_FLIP .....

    allocate(idx(n))
    ! -- compare inner products..
    do i=1,n
      inner = matmul(U1(:, i), U0(:,:))
      idx(i) = maxloc(abs(inner), 1)
    enddo

    call move_alloc(idx, res)

  end function get_eigenvector_permutation_idx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine K2sincos(kmat, kmat_eval_energies, sine, cosine, elec_channels, spinmult)
    !! Take a K-matrix that is defined for several energies, and
    !! transform it into the S-matrix at those energies via diagonalization:
    !!   1) Diagonalize K(E) ~ tan(δ(E))
    !!   2) Identify eigenphases across (E)
    !!   3) Correct for ±π shifts due to branch cuts
    !!   4) Build the sine and cosine matrices. These will be interpolated linearly later in the RFT and
    !!      closed channel elimination procedure.
    !!   5) Write eigenphase info to file

    use rotex__types,      only: elec_channel_type
    use rotex__linalg,     only: dsyev
    use rotex__arrays,     only: size_check, is_unitary, unitary_defect, adjoint, is_symmetric
    use rotex__characters, only: i2c => int2char
    use rotex__constants,  only: spinmult_names, au2ev, pi
    use rotex__writing,    only: write_elec_mat_elems_to_file

    implicit none (type, external)

    real(dp),                intent(in)  :: kmat(:,:,:)
    real(dp),                intent(in)  :: kmat_eval_energies(:)
    real(dp),                intent(out) :: sine(:,:,:), cosine(:,:,:)
    type(elec_channel_type), intent(in)  :: elec_channels(:)
    integer,                 intent(in)  :: spinmult
      !! The current spin multiplicity

    logical :: add_energy_pre, add_energy_post
    integer :: ne, n, i, j, ie
    integer :: funit_eigenphases, funit_sine, funit_cosine
    integer, allocatable :: idx(:)
    real(dp) :: D1, D2, D3, D4, h1, h2
    real(dp) :: E_pre, E_post, E1, E2
    real(dp), allocatable :: sin_pre(:), sin_post(:), cos_pre(:), cos_post(:)
    real(dp), allocatable :: tmp(:)
    real(dp), allocatable :: eval_E_pre(:), eval_E_post(:)
    real(dp), allocatable :: U(:,:,:), eigenphases(:,:), eigenphases_unwrapped(:,:)
    character(:), allocatable :: eigenphases_dir
    character(:), allocatable :: eigenphases_file, sin_file, cos_file

    ! -- lapack variables
    integer :: info, lwork
    real(dp), allocatable :: w(:), work(:)
    character(1), parameter :: UPLO = 'U'

    n  = size(kmat, 1)
    ne = size(kmat, 3)

    call size_check(kmat_eval_energies, ne,         "KMAT_EVAL_ENERGIES")
    call size_check(sine,               [n, n, ne], "SINE_FLAT")
    call size_check(cosine,             [n, n, ne], "COSINE_FLAT")
    call size_check(elec_channels,      n,          "ELEC_CHANNELS")

    lwork = 3*n+1
    allocate(U(n,n,ne), source=0._dp)
    allocate(eigenphases(n, ne))
    allocate(work(lwork))

    ! -- 1) diagonalize K, get eigenphases, permute them to ensure consistent ordering across geometries
    do ie=1, ne

      U(:,:,ie) = kmat(:,:,ie)

      call dsyev('V', UPLO, n, U(:,:,ie), n, eigenphases(:, ie), work, lwork, info)
      if(info .ne. 0) call die("DSYEV exited with nonzero INFO = " // i2c(info))

      ! -- δ <- tan(δ)
      eigenphases(:, ie) = atan(eigenphases(:, ie))

      ! -- need two energies for comparison
      if(ie .eq. 1) cycle

      ! -- 2)  identify eigenvalues by inner product of eigenvectors to permute columns of
      !        the matrix U
      ! -- compare with the previous energy to identify eigenvectors/eigenphases. This assumes
      !    that the inner product doesn't change too much between energies, so a relatively dense
      !    grid may be necessary near steep resonances. .
      idx = get_eigenvector_permutation_idx(U(:, :, ie-1), U(:, :, ie), eigenphases(:, ie-1), eigenphases(:, ie))
      U(:,:,ie) = U(:,idx,ie)
      eigenphases(:, ie) = eigenphases(idx, ie)

    enddo

    ! -- 3) correct jumps by ±π
    allocate(eigenphases_unwrapped, source=eigenphases)
    call unwrap_eigenphases(eigenphases_unwrapped)

    sine(:,:,:)   = 0._dp
    cosine(:,:,:) = 0._dp

    ! -- 4) construct the sine and cosine matrices
    do ie=1,ne

      do concurrent (i=1:n)
        sine(i,i,ie)   = sin(eigenphases_unwrapped(i, ie))
        cosine(i,i,ie) = cos(eigenphases_unwrapped(i, ie))
      enddo

      ! -- Eigenphases δ -> sin(δ), cos(δ)
      sine(:,:,ie)   = matmul(U(:,:,ie), matmul(sine(:,:,ie),   adjoint(U(:,:,ie))))
      cosine(:,:,ie) = matmul(U(:,:,ie), matmul(cosine(:,:,ie), adjoint(U(:,:,ie))))

      if(is_symmetric(sine(:,:,ie)) .AND. is_symmetric(cosine(:,:,ie))) cycle

      ! -- error
      write(stderr, '("The Sine/Cosine matrices are not symmetric for energy ", I0, ": ", E15.7, " eV")') &
        kmat_eval_energies(ie)*au2ev
      write(stderr, '("Maxval( |sin-transpose(sine)| ): ", E15.7)') maxval(abs(sine(:,:,ie)-transpose(sine(:,:,ie))))
      write(stderr, '("Maxval( |cosin-transpose(cosine)| ): ", E15.7)') maxval(abs(cosine(:,:,ie)-transpose(cosine(:,:,ie))))
      call die("Non-symmetric electronic sine/cosine matrix detected")

    enddo

    ! -- 5) write to file
    call write_elec_mat_elems_to_file(kmat_eval_energies, elec_channels, eigenphases, eigenphases_unwrapped &
      , sine, cosine, spinmult_names(spinmult))

  end subroutine K2sincos

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine K2S_cayley(K, S)
    !! electronic Kmat -> electronic Smat via Cayley transform S = (I + ik) / (I - iK).
    !! Also ensures that the S-matrix is in the basis of complex-valued spherical harmonics.
    !! This is to be used when G%EDFT is false, it assumes only one evaluation energy
    use rotex__linalg, only: zgesv
    use rotex__arrays, only: eye
    implicit none (type, external)
    real(dp),    intent(in)  :: K(:,:)
    complex(dp), intent(out) :: S(:,:)
    complex(dp), allocatable :: A(:,:)
    integer :: n, info
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: I(:,:)
    n = size(K, 1)
    allocate(ipiv(n))
    I = real(eye(n), kind=dp)
    A = cmplx(I, -K, kind=dp)
    S = cmplx(I,  K, kind=dp)
    call zgesv(n, n, A, n, ipiv, S, n, info)
    if(info .eq. 0) return
    write(stderr, '("INFO = ", I0)') INFO
    call die("ZGESV exited with nonzero INFO")
  end subroutine K2S_cayley

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine setup_RFT_block(J, channels_l, channels_l_J, idx, rot_channels, nchans_J)
  !   implicit none (type, external)
  !   ! @@@
  ! end subroutine setup_RFT_block

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  !subroutine do_rft(             &
  !      smat_elec                &
  !    , sin_elec                 &
  !    , cos_elec                 &
  !    , kmat_eval_energies       &
  !    , Smat_j                   &
  !    , jmin                     &
  !    , jmax                     &
  !    , n_states                 &
  !    , elec_channels            &
  !    , asymtop_rot_channels_l   &
  !    , asymtop_rot_channels_l_j &
  !  )
  !  !! Perform the rotational frame transformation on the electronic S-matrix

  !  use rotex__kinds, only: dp
  !  use rotex__types, only: cmatrix_type, asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type &
  !                        , elec_channel_type, N_states_type, cvector_type
  !  use rotex__wigner,     only: clebsch
  !  use rotex__system,     only: stdout, die
  !  use rotex__arrays,     only: realloc, is_unitary, is_symmetric, uniq
  !  use rotex__functions,  only: neg
  !  use rotex__characters, only: i2c => int2char

  !  implicit none (type, external)

  !  complex(dp),        intent(in)  :: smat_elec(:,:)
  !  real(dp),        intent(in)     :: sin_elec(:,:,:), cos_elec(:,:,:)
  !  real(dp),        intent(in)     :: kmat_eval_energies(:)
  !  type(cmatrix_type), intent(out) :: Smat_j(jmin:jmax, size(kmat_eval_energies, 1))
  !    !! Rotationally resolved S-matrix at each J
  !  integer, intent(in) :: jmin, jmax
  !    !! Min/max values of the total angular momentum J
  !  type(n_states_type), intent(in) :: n_states(:)
  !  type(elec_channel_type), intent(in) :: elec_channels(:)
  !  type(asymtop_rot_channel_l_type), intent(in) :: asymtop_rot_channels_l(:)
  !  type(asymtop_rot_channel_l_vector_type), intent(out) :: asymtop_rot_channels_l_j(jmin:jmax)

  !  type(asymtop_rot_channel_l_type), allocatable :: rot_channels(:)
  !  logical :: flag
  !  logical, allocatable :: mask(:)
  !  integer :: nsyms, nchans_elec, nchans_elec_flat
  !  integer :: J, i, isym, sym, nchans_sym, ie
  !  integer :: nchans_J, ne
  !  integer, allocatable :: idx(:), uniq_syms(:)
  !  complex(dp), allocatable :: Smat_rot(:,:), Smat_rot_sym(:,:)
  !  complex(dp), allocatable :: cos_rot(:,:), sin_rot(:,:), sin_rot_sym(:,:), cos_rot_sym(:,:)
  !  complex(dp), allocatable :: U(:,:)

  !  flag = .false.
  !  nchans_elec = size(elec_channels, 1)

  !  ne = size(kmat_eval_energies, 1)

  !  ! -- loop over different values of the agular momentum J
  !  jloop: do J=Jmin,Jmax

  !    ! -- determine number of rotational channels in this block of total J
  !    call collect_J_channels_indices(j, asymtop_rot_channels_l, idx)
  !    asymtop_rot_channels_l_j(j) % channels = asymtop_rot_channels_l(idx)

  !    nchans_J = size(asymtop_rot_channels_l_j(j) % channels, 1)

  !    ! -- the total S-matrix for this J
  !    call realloc(U,        nchans_J, nchans_elec)
  !    call realloc(smat_rot, nchans_J, nchans_J)
  !    U = 0
  !    smat_rot = 0

  !    ! -- get the unique symmetry elements. We do one frame transformation per
  !    !    symmetry element
  !    uniq_syms = asymtop_rot_channels_l_j(j) % channels % sym
  !    uniq_syms = uniq(uniq_syms)
  !    nsyms = size(uniq_syms, 1)

  !    ! call do_rft_no_sym(J, n_states, elec_channels, asymtop_rot_channels_l_j(j)%channels, smat_elec, smat_rot, U, point_group)

  !    if(G%EDFT) then
  !      allocate(sin_elec(nchans_elec, nchans_elec, ne), source=0._dp)
  !      allocate(cos_elec(nchans_elec, nchans_elec, ne), source=0._dp)
  !    else
  !      allocate(smat_elec(nchans_elec, nchans_elec), source=(0._dp, 0._dp))
  !    endif

  !    !TODO omp
  !    nrg: do ie=1, ne


  !      ! -- loop over symmetries
  !      ! do isym=1, nsyms

  !      !TODO: symmetry loop. loop over symmetries in the RFT. This is not needed if the electronic
  !      ! calculations are in the full symmetry. expand to full symmetry beforehand if needed
  !      ! TODO: get rid of symmetry stuff
  !      sym = uniq_syms(isym)

  !      ! -- map all J -> this sym
  !      mask = asymtop_rot_channels_l_j(j)%channels%sym .eq. sym
  !      idx = pack([(i,i=1,nchans_j)], mask)
  !      rot_channels = asymtop_rot_channels_l_j(j) % channels(idx)

  !      ! -- allocate U, Smat_rot/Sine_rot/Cosine_rot for this sym
  !      nchans_sym = size(idx, 1)
  !      call realloc(U,            nchans_sym, nchans_elec)
  !      if(G%EDFT) then
  !        call realloc(sin_rot_sym, nchans_sym, nchans_sym)
  !        call realloc(cos_rot_sym, nchans_sym, nchans_sym)
  !        sin_rot_sym = 0
  !        cos_rot_sym = 0
  !      else
  !        call realloc(smat_rot_sym, nchans_sym, nchans_sym)
  !        smat_rot_sym = 0
  !      endif
  !      U = 0

  !      call do_rft_this_sym( &
  !          j                 &
  !        , sym               &
  !        , n_states          &
  !        , elec_channels     &
  !        , rot_channels      &
  !        , smat_elec         &
  !        , sin_elec          &
  !        , cos_elec          &
  !        , smat_rot_sym      &
  !        , sin_rot_sym       &
  !        , cos_rot_sym       &
  !        , U                 &
  !      )

  !      ! -- add this contribution back to the total S/Sin/Cos-matrix for this J
  !      if(G%EDFT) then
  !        sin_rot(idx, idx) = sin_rot_sym(:,:)
  !        cos_rot(idx, idx) = cos_rot_sym(:,:)
  !      else
  !        smat_rot(idx, idx) = smat_rot_sym(:,:)
  !      endif

  !      ! -- export this S^J(E)
  !      if(G%EDFT) call sincos2s(sin_rot, cos_rot, smat_rot)
  !      if(is_symmetric(smat_rot) .eqv. .false.) call die("The S-matrix is not symmetric after the RFT")

  !      if(is_unitary(Smat_rot) .eqv. .true.) cycle

  !      flag = .true.

  !      ! -- warn about nonunitarity
  !      block
  !        use rotex__system, only: stderr
  !        use rotex__arrays, only: eye, adjoint, norm_frob, unitary_defect
  !        associate(S => smat_rot)
  !          write(stderr, '(A, I0, A, F7.5)') &
  !            "WARN: The S-matrix for J = ", J, " is nonunitary with unitary defect ", unitary_defect(S)
  !        end associate
  !      end block

  !    enddo nrg

  !  enddo jloop

  !  if(flag .eqv. .false.) return

  !  call die("At least one J-block of the S-matrix is non-unitary. This may cause some issues in the&
  !    & ensuing MQDT closed-channel elimination procedure which takes the closed channels into account for each J-block.&
  !    & Therefore, each J-block should be unitary, even if they involve states with N< N_min or N > N_max. It is probably&
  !    & worth noting that the S-matrix at this point was detected to be non-unitary, but each symmetry sub-block was unitary.")

  !end subroutine do_rft

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine sincos2s(sinec, cosinec, s, did_symmetrize)
    !! Sine,Cosine -> K=SC⁻¹ -> S
    use rotex__system, only: die
    use rotex__arrays, only: is_symmetric
    use rotex__linalg, only: right_divide
    implicit none (type, external)
    complex(dp), intent(in)  :: sinec(:,:), cosinec(:,:)
    complex(dp), intent(out) :: s(:,:)
    logical,     intent(out) :: did_symmetrize
    real(dp) :: TOL
    real(dp), allocatable :: sine(:,:), cosine(:,:), K(:,:)
    TOL = G%POST_RFT_SINCOS2S_IMAG_TOL
    if(any(abs(sinec%im)   .gt. TOL) .OR. any(abs(cosinec%im) .gt. TOL)) then
      write(stderr, '("maxval(|sin%im|): ", ES20.12)') maxval(abs(sinec%im))
      write(stderr, '("maxval(|cos%im|): ", ES20.12)') maxval(abs(cosinec%im))
      call die("SINE/COSINE matrix is non-real after frame transformation")
    endif
    sine   = sinec % re
    cosine = cosinec % re
    ! -- K=Sin*Cos⁻¹
    K = right_divide(sine, cosine)
    ! -- symmetrize K if needed
    if(is_symmetric(K) .eqv. .false.) then
      K = (K + transpose(K)) * 0.5_dp
      did_symmetrize = .true.
    endif
    ! -- S = (I+iK)/(I-iK)
    call K2S_cayley(K, S)
  end subroutine sincos2s

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine real2complex_ylm_r(M, chans)
    !! Real-valued version of the complex-valued equivalent
    use rotex__types,  only: elec_channel_type
    use rotex__system, only: die
    implicit none (type, external)
    real(dp), intent(inout) :: M(:,:)
      !! The S/K-matrix
    type(elec_channel_type), intent(in) :: chans(:)
    complex(dp), allocatable :: MC(:,:)
    MC = cmplx(M, 0.0_dp, kind = dp)
    call real2complex_ylm_c(MC, chans)
    if(maxval(abs(MC%im)) .gt. 1e-12) call die("Nonzero imaginary values detected in&
      & transformed S/K matrix (complex spherical harmonics basis)")
    M = MC % re
  end subroutine real2complex_ylm_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine real2complex_ylm_c(M, chans)
    !! Take an S/K-matrix that is in a basis of electronic channels for exactly
    !! one electronic state and a basis of real-valued spherical harmonics, and
    !! transform it to an S/K-matrix in a basis of the same electronic state but
    !! complex-valued spherical harmonics for the scattering electron

    use rotex__kinds,  only: dp
    use rotex__types,  only: elec_channel_type
    use rotex__arrays, only: adjoint, size_check, is_unitary
    use rotex__system, only: die

    implicit none (type, external)

    complex(dp), intent(inout) :: M(:,:)
      !! The S/K-matrix
    type(elec_channel_type), intent(in) :: chans(:)
      !! The electronic channels

    integer :: n, i, j
    integer :: eleci, li, lambi, elecj, lj, lambj
    complex(dp), allocatable :: U(:,:)

    n = size(chans, 1)
    call size_check(M, [n,n], "M")

    allocate(U(n,n)) ; U = 0

    ! -- Build the ℭ → ℛ  spherical harmonics transformatices
    !    (incident electron partial waves)
    do j=1, n
      elecj = chans(j) % nelec
      lj    = chans(j) % l
      lambj = chans(j) % ml
      do i=1, n
        eleci = chans(i) % nelec
        li    = chans(i) % l
        lambi = chans(i) % ml
        if(elecj .ne. eleci) call die("There are >1 electronic states detected in the S/K-matrix")
        if(li    .ne. lj) cycle
        ! -- transform incident electron partial waves
        ! select case(point_group)
        ! case("c2v", "c2", "d2", "d2h")
          U(i,j) = ylm_r2c(lambi, lambj)
        ! case("cs")
        !   U(i,j) = ylm_r2c(lambi, lambj)
          ! U(i,j) = ylm_r2c_cs(lambi, lambj)
          ! U(i,j) = ylm_r2c(lambi, lambj)
          ! U(i,j) = ylm_c2r_cs(lambi, lambj)
        ! case default
        !   call die("Unsupported point group in REAL2COMPLEX_YLM transformation: " // point_group)
        ! end select
      enddo
    enddo

    ! -- unitary checks
    if( is_unitary(U) .eqv. .false.) call die("The transformation matrix for the spherical harmonics is not unitary")

    ! -- ℛ  → ℭ
    ! select case(point_group)
    ! case("c2v", "c2", "d2", "d2h")
      M = matmul(adjoint(U), matmul(M, U))
    ! case("cs")
    !   M = matmul(adjoint(U), matmul(M, U))
    ! end select

    deallocate(U)

  end subroutine real2complex_ylm_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function ylm_r2c(mr, mc) result(res)
    !! Determine the coefficient of the real spherical harmonics for expressing the
    !! complex spherical harmonics in terms of the real spherical harmonics, e.g.,
    !!   \(Y_{l}^m = c_1 X_{l|m|} + c_2X_{l-|m|}\),
    !! given mc=m and one of mr=|m| or mr=-|m|. Assumes that the degree (l) of the spherical harmonics
    !! is the same
    use rotex__kinds, only: dp
    use rotex__functions, only: neg
    implicit none (type, external)
    integer, intent(in) :: mr
      !! The order (m) for the real spherical harmonic
    integer, intent(in) :: mc
      !! The order (±|m|) for the complex spherical harmonic
    complex(dp), parameter :: one = cmplx(1, 0, kind=dp)
    real(dp),    parameter :: invsq2 = 1.0_dp/sqrt(real(2, kind=dp))
    complex(dp), parameter :: im  = cmplx(0, 1, kind=dp)
    logical :: is_cos
    complex(dp) :: coef1, coef2
    complex(dp) :: res
    res = 0
    if(abs(mc) .ne. abs(mr)) return
    is_cos = mr .gt. 0
    select case(mc)
    case(:-1)
      coef1 = 1
      coef2 = merge(one, -im, is_cos)
    case(1:)
      coef1 = neg(mc)
      coef2 = merge(one,  im, is_cos)
    case(0)
      res = 1
      return
    end select
    res = coef1 * coef2 * invsq2
  end function ylm_r2c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine get_channel_qnums_rot(channels, irot, nelec, N, Ka, Kc, l, sym)
    !! Given an array of channels with quantum numbers, extract the quantum numbers at index irot
    use rotex__types,  only: asymtop_rot_channel_l_type
    use rotex__system, only: die
    implicit none (type, external)
    type(asymtop_rot_channel_l_type), intent(in)  :: channels(:)
    integer,                          intent(in)  :: irot
    integer,                          intent(out) :: nelec, N, Ka, Kc, l
    integer, optional,                intent(out) :: sym
    if(irot .gt. ubound(channels, 1)) call die("Targeted channel > channel upper bound")
    if(irot .lt. lbound(channels, 1)) call die("Targeted channel < channel lower bound")
    nelec  = channels(irot) % nelec
    N      = channels(irot) % N
    Ka     = channels(irot) % Ka
    Kc     = channels(irot) % Kc
    l      = channels(irot) % l
    if(present(sym) .eqv. .false.) return
    sym    = channels(irot) % sym
  end subroutine get_channel_qnums_rot

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine collect_j_channels_indices(j, channels_l, idx)
    !! Go through channels_l and add the channels to channels_l_j that
    !! obey the degenerate triangle inequality for N, l, J
    use rotex__types,     only: asymtop_rot_channel_l_type
    use rotex__functions, only: istriangle
    implicit none (type, external)
    integer, intent(in) :: j
      !! Total angular momentum J
    type(asymtop_rot_channel_l_type), intent(in) :: channels_l(:)
      !! All rotatinal channels
    integer, intent(out), allocatable :: idx(:)
      !! Rotational channels for this J will be channels_l(idx)
    integer :: nchans_rot
    integer :: n, l
    integer :: ichan
    integer :: count
    nchans_rot = size(channels_l, 1)
    count = 0
    ! -- count number of channels
    do ichan=1, nchans_rot
      n     = channels_l(ichan) % n
      l     = channels_l(ichan) % l
      if(istriangle(n, l, j) .eqv. .false.) cycle
      count = count + 1
    enddo
    allocate(idx(count), source=0)
    if(count .eq. 0) return
    ! -- fill idx if count > 0
    count = 0
    do ichan=1, nchans_rot
      n     = channels_l(ichan) % n
      l     = channels_l(ichan) % l
      if(istriangle(n, l, j) .eqv. .false.) cycle
      count = count + 1
      idx(count) = ichan
    enddo
  end subroutine collect_j_channels_indices

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine do_rft_no_sym(j, n_states, elec_channels, rot_channels, Smat_elec, Smat_rot, U, point_group)
  !   !! Do the rotational frame transformation for a specific symmetry
  !   use rotex__kinds,      only: dp
  !   use rotex__types,      only: elec_channel_type, asymtop_rot_channel_l_type, n_states_type
  !   use rotex__arrays,     only: size_check, is_unitary, is_symmetric, adjoint
  !   use rotex__wigner,     only: clebsch
  !   use rotex__system,     only: die, stderr, stdout
  !   use rotex__functions,  only: neg
  !   use rotex__characters, only: i2c => int2char
  !   implicit none (type, external)
  !   integer,                          intent(in)  :: j
  !     !! Total angular momentum quantum number J
  !   type(n_states_type),              intent(in)  :: n_states(:)
  !     !! N, Ka, and Kc for each N
  !   type(elec_channel_type),          intent(in)  :: elec_channels(:)
  !     !! Electronic channel basis for Smat_elec
  !   type(asymtop_rot_channel_l_type), intent(in)  :: rot_channels(:)
  !     !! Rotational channel basis for Smat_rot (this symmetry)
  !   complex(dp),                      intent(in)  :: smat_elec(:,:)
  !     !! Electronic S-matrix
  !   complex(dp),                      intent(out) :: smat_rot(:,:)
  !     !! Rotatinal S-matrix
  !   complex(dp),                         intent(out) :: U(:,:)
  !     !! Unitary transformation matrix
  !   character(1),                     intent(in) :: point_group
  !     !! The point group of the scattering calculations
  !   integer :: irot
  !   integer :: nchans_elec, nchans_rot
  !   integer :: ni, kai, kci, li, lj, ki, lambdaj, symchan
  !   integer :: neleci, nelecj
  !   integer :: in, itau, iK, jelec
  !   integer :: Omega
  !   logical, allocatable :: mask(:)
  !   complex(dp), allocatable :: C(:, :)
  !   nchans_rot  = size(rot_channels, 1)
  !   nchans_elec = size(elec_channels, 1)
  !   Smat_rot = 0
  !   ! -- build the rectangular transformation matrix U <LF|BF> for each Ω
  !   allocate(C(nchans_rot, nchans_rot))
  !   C = 0
  !   do Omega = -J,J
  !     U = 0
  !     do irot = 1, nchans_rot
  !       call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li, symchan)
  !       in    = findloc(n_states % n, value = ni, dim = 1)
  !       mask = (n_states(in) % ka(:) .eq. kai) .AND. (n_states(in) % kc(:) .eq. kci)
  !       itau  = findloc(mask, value = .true., dim = 1)
  !       do jelec = 1, nchans_elec
  !         nelecj  = elec_channels(jelec) % nelec
  !         lj      = elec_channels(jelec) % l
  !         lambdaj = elec_channels(jelec) % ml
  !         ! -- enforce transformation between the same electronic state n and partial wave l
  !         if (neleci .ne. nelecj) cycle
  !         if (li     .ne. lj) cycle
  !         Ki  = Omega - lambdaj
  !         if(abs(Ki) .gt. Ni) cycle
  !         ik = Ki + Ni + 1
  !         U(irot, jelec) = neg(lj + lambdaj)             &
  !             * N_states(in) % eigenH % eigvecs(ik, itau)&
  !             * clebsch(lj, -lambdaj, J, Omega, Ni, Ki)
  !       enddo
  !     enddo
  !     ! -- S^J = Σ_Ω USU+ (for each Ω)
  !     Smat_rot = Smat_rot + matmul( U, matmul(Smat_elec, adjoint(U)) )
  !     ! -- diagnostics if we fail to produce a unitary S
  !     C = C + matmul(U, adjoint(U))
  !   enddo
  !
  !   if(is_symmetric(Smat_rot) .eqv. .false.) then
  !     write(stderr, '("Maxval S-S+: ", E30.20)') maxval(abs(Smat_rot-adjoint(Smat_rot)))
  !     call die("The S-matrix is not symmetric for J = " // i2c(j) // " ❌")
  !   endif
  !
  !   if(is_unitary(Smat_rot)   .eqv. .true.) then
  !     write(stdout, '("S-matrix is unitary for J = ", I0, " ✔️")') J
  !     return
  !   endif
  !
  !   err: block
  !     use rotex__arrays, only: unitary_defect
  !     write(stderr, '("Channels: ", 6(A5,X), A20)') "i", "nelec", "N", "Ka", "Kc", "l", "Σ|S(i,:)|²"
  !     do irot=1, nchans_rot
  !       call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li)
  !       write(stderr, '(10X, 6(I5,X), E20.10)', advance = "no") irot, neleci, ni, kai, kci, li &
  !         , sum(abs(Smat_rot(irot,:))**2)
  !       if(all(U(irot,:) .eq. 0._dp)) write(stderr, '(" <-- ", A)', advance = "no") "Does not couple to any electronic channels !"
  !       write(stderr, *)
  !     enddo
  !     write(stderr, *)
  !     write(stderr, '("Rank of UU+: ", I0)') rank(C)
  !     write(stderr, '("Unitary defect in UUT:  ", F15.9)') unitary_defect(C)
  !     write(stderr, '("Unitary defect in USUT: ", F15.9)') unitary_defect(Smat_rot)
  !     call die("The S-matrix is not unitary for J = " // i2c(J) // " ❌")
  !   end block err
  !
  ! end subroutine do_rft_no_sym

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine do_rft_this_sym( &
  !       j                     &
  !     , sym                   &
  !     , n_states              &
  !     , elec_channels         &
  !     , rot_channels          &
  !     , whichmat              &
  !     , Smat_elec             &
  !     , sin_elec              &
  !     , cos_elec              &
  !     , Smat_rot              &
  !     , sin_rot               &
  !     , cos_rot               &
  !     , U                     &
  !   )
  !   !! Do the rotational frame transformation for a specific symmetry.
  !   !!   Energy independent -> S-matrix
  !   !!   Energy dependent ---> Sine and Cosine matrices, evaluation energy matching the outgoing channel
  !   use rotex__kinds,      only: dp
  !   use rotex__types,      only: elec_channel_type, asymtop_rot_channel_l_type, n_states_type
  !   use rotex__arrays,     only: size_check, is_unitary, is_symmetric, adjoint, linear_interpolation, idx_binsearch
  !   use rotex__wigner,     only: clebsch
  !   use rotex__system,     only: die, stderr, stdout
  !   use rotex__functions,  only: neg
  !   use rotex__characters, only: i2c => int2char

  !   implicit none (type, external)

  !   integer,                          intent(in)  :: j
  !     !! Total angular momentum quantum number J
  !   integer,                          intent(in)  :: sym
  !     !! The current symmetry
  !   type(n_states_type),              intent(in)  :: n_states(:)
  !     !! N, Ka, and Kc for each N
  !   type(elec_channel_type),          intent(in)  :: elec_channels(:)
  !     !! Electronic channel basis for Smat_elec
  !   type(asymtop_rot_channel_l_type), intent(in)  :: rot_channels(:)
  !     !! Rotational channel basis for Smat_rot (this symmetry)
  !   character(4),                     intent(in)  :: whichmat
  !     !! 'SMAT': use the S-matrix only
  !     !! 'SICO': use the Sine and Cosine matrices
  !   complex(dp),                      intent(in)  :: smat_elec(:,:)
  !     !! Electronic S-matrix
  !   real(dp),                         intent(inout)  :: sin_elec(:,:,:)
  !     !! Electronic Sine-matrix
  !   real(dp),                         intent(inout)  :: cos_elec(:,:,:)
  !     !! Electronic Cosine-matrix
  !   complex(dp),                      intent(out) :: sin_rot(:,:)
  !     !! Rotatoinal Sine-matrix
  !   complex(dp),                      intent(out) :: cos_rot(:,:)
  !     !! Rotatoinal Cosine-matrix
  !   complex(dp),                      intent(out) :: smat_rot(:,:)
  !     !! Rotatoinal S-matrix
  !   complex(dp),                      intent(inout) :: U(:,:)
  !     !! Unitary transformation matrix

  !   integer :: irot
  !   integer :: nchans_elec, nchans_rot
  !   integer :: ni, kai, kci, li, nj, kaj, kcj, lj, ki, lambdaj, symchan
  !   integer :: neleci, nelecj
  !   integer :: in, itau, iK, jelec
  !   integer :: Omega
  !   logical, allocatable :: mask(:)
  !   complex(dp), allocatable :: C(:,:)

  !   nchans_rot  = size(rot_channels, 1)
  !   nchans_elec = size(elec_channels, 1)

  !   select case(whichmat)
  !   case("SMAT")
  !     Smat_rot = 0
  !   case("SICO")
  !     sin_rot = 0
  !     cos_rot = 0
  !   case default
  !     call die("Unacceptable value of WHICHMAT: "//whichmat//". Must be 'SMAT' or 'SICO'")
  !   end select

  !   allocate(C, source=U) ; C = 0

  !   ! -- build the rectangular transformation matrix U <LF|BF> for each Ω
  !   do Omega = -J,J

  !     U = 0

  !     do irot = 1, nchans_rot

  !       call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li, symchan)

  !       ! if(sym .ne. symchan) call die("Channel symmetry does not match transformation symmetry !")

  !       ! -- get the corresponding eigenvector for this state, indexed by itau.
  !       in    = findloc(n_states % n, value = ni, dim = 1)
  !       select case(G%ROTOR_KIND)
  !       case("a", "A")
  !         mask = (n_states(in) % ka(:) .eq. kai) .AND. (n_states(in) % kc(:) .eq. kci)
  !       case("s", "S")
  !         select case(G%ROTOR_ZAXIS)
  !         case("a", "A")
  !           mask = (n_states(in) % ka(:) .eq. kai)
  !         case("c", "C")
  !           mask = (n_states(in) % kc(:) .eq. kci)
  !         case default
  !           call die("Somehow got a symmetric top with a G%ROTOR_ZAXIS " // G%ROTOR_ZAXIS // " that is neither A nor C")
  !         end select
  !       case default
  !         call die("ROTOR_KIND " // G%ROTOR_KIND // " not allowed in RFT")
  !       end select
  !       itau  = findloc(mask, value = .true., dim = 1)

  !       do jelec = 1, nchans_elec
  !         nelecj  = elec_channels(jelec) % nelec
  !         lj      = elec_channels(jelec) % l
  !         lambdaj = elec_channels(jelec) % ml

  !         ! -- enforce transformation between the same electronic state n and partial wave l
  !         if (neleci .ne. nelecj) cycle
  !         if (li     .ne. lj) cycle
  !         Ki  = Omega - lambdaj
  !         if(abs(Ki) .gt. Ni) cycle

  !         ik = Ki + Ni + 1

  !         U(irot, jelec) = neg(lj + lambdaj)              &
  !             * N_states(in) % eigenH % eigvecs(ik, itau) &
  !             * clebsch(lj, -lambdaj, J, Omega, Ni, Ki)

  !       enddo
  !     enddo

  !     @@@@@
  !     - we should pass the energy in from the outer loop ? an duse that to search and get indicde, then get our Sine cosien
  !       mats.
  !     - pass this routine all the enegies/matrices that it will need
  !     - elementiwise interpolation based on rhs ? use array linear_interpolate with channel energies
  !       - array slices
  !     figure out edft here. shoudl be energy of adjoint U; do an interpolation ?


  !     ! -- M^J = Σ_Ω UMU⁺ (for each Ω), where M is some matrix in the basis of channels
  !     select case(whichmat)
  !     case("SMAT")
  !       Smat_rot = Smat_rot + matmul( U, matmul(Smat_elec, adjoint(U)) )
  !     case("SICO")
  !       ie_elec = [( idx_binserach )]
  !       call linear_interpolation()
  !       sin_rot = sin_rot + matmul( U, matmul(sin_elec, adjoint(U)) )
  !       cos_rot = cos_rot + matmul( U, matmul(cos_elec, adjoint(U)) )
  !     end select

  !     C = C + matmul(U, adjoint(U))

  !   enddo

  !   error_checks: block
  !     use rotex__utils,  only: printmat
  !     use rotex__system, only: warn
  !     use rotex__arrays, only: unitary_defect, eye, norm_frob
  !     logical :: symflag = .false.
  !     logical :: unitaryflag = .false.

  !     select case(whichmat)
  !     case("SMAT")

  !       if(is_symmetric(Smat_rot) .eqv. .false.) then
  !         ! -- not symmetric
  !         symflag = .true.
  !         call warn("The S-matrix is not symmetric for symmetry " // i2c(sym) // " ❌")
  !       endif

  !       if(is_unitary(Smat_rot, 1e-7_dp)   .eqv. .true.) then
  !         write(stdout, '("S-matrix is unitary for J = ", I0, ", symmetry ", I0, " ✔️")') J, sym
  !         if(symflag .eqv. .false.) return
  !       else
  !         unitaryflag = .true.
  !       endif


  !       write(stderr, '("Symmetry: ", I0)') sym
  !       write(stderr, '("Channels: ", 6(A5,X), A20)') "i", "nelec", "N", "Ka", "Kc", "l", "Σ|S(i,:)|²"
  !       do irot=1, nchans_rot
  !         call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li)
  !         write(stderr, '(10X, 6(I5,X), E20.10)', advance = "no") irot, neleci, ni, kai, kci, li &
  !           , sum(abs(Smat_rot(irot,:))**2)
  !         if(all(U(irot,:) .eq. 0._dp)) write(stderr, '(" <-- ", A)', advance = "no") "Does not couple to any electronic channels !"
  !         write(stderr, *)
  !       enddo
  !       write(stderr, *)
  !       write(stderr, '("This is symmetry ", I0, ", J = ", I0)') sym, J
  !       write(stderr, '(A30, F15.9)') "Unitary defect in UU⁺: ", unitary_defect(C)
  !       write(stderr, '(A30, F15.9)') "Unitary defect in USU⁺: ", unitary_defect(Smat_rot)
  !       if(unitaryflag) call warn("The S-matrix is not unitary for symmetry " // i2c(sym) // " ❌")
  !       if(symflag) call warn("The S-matrix is not symmetric for symmetry " // i2c(sym) // " ❌")
  !       if(unitaryflag .or. symflag) error stop

  !     case("SICO")

  !       if(is_symmetric(sin_rot) .eqv. .false.) then
  !         ! -- not symmetric
  !         symflag = .true.
  !         call warn("The Sine-matrix is not symmetric for symmetry " // i2c(sym) // " ❌")
  !       endif
  !       if(is_symmetric(cos_rot) .eqv. .false.) then
  !         ! -- not symmetric
  !         symflag = .true.
  !         call warn("The Cosine-matrix is not symmetric for symmetry " // i2c(sym) // " ❌")
  !       endif

  !       if(symflag .eqv. .false.) exit error_checks

  !       write(stderr, '("Symmetry: ", I0)') sym
  !       write(stderr, '("Channels: ", 6(A5,X), A20)') "i", "nelec", "N", "Ka", "Kc", "l", "Σ|S(i,:)|²"
  !       do irot=1, nchans_rot
  !         call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li)
  !         write(stderr, '(10X, 6(I5,X), E20.10)', advance = "no") irot, neleci, ni, kai, kci, li &
  !           , sum(abs(Smat_rot(irot,:))**2)
  !         if(all(U(irot,:) .eq. 0._dp)) write(stderr, '(" <-- ", A)', advance = "no") "Does not couple to any electronic channels !"
  !         write(stderr, *)
  !       enddo
  !       write(stderr, *)
  !       write(stderr, '("This is symmetry ", I0, ", J = ", I0)') sym, J
  !       write(stderr, '(A30, F15.9)') "Unitary defect in UU⁺: ", unitary_defect(C)

  !     end select

  !   end block error_checks

  ! end subroutine do_rft_this_sym

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine unwrap_eigenphases(eigenphases)
    !! Unwraps each phase index of eigenphases so that, at consecutive energy steps, the values change by at most π/2
    use rotex__constants, only: pi
    implicit none (type, external)
    real(dp), intent(inout) :: eigenphases(:,:)
    real(dp), parameter :: PERIOD = pi
    integer :: ne, nphases, ie, iphase
    real(dp) :: halfperiod, prev, dphase
    nphases = size(eigenphases, 1)
    ne      = size(eigenphases, 2)
    if(ne .eq. 1) return
    halfperiod = PERIOD / 2
    do concurrent (iphase=1:nphases)
      prev = eigenphases(iphase, 1)
      do ie=2, ne
        dphase = eigenphases(iphase, ie) - prev
        dphase = modulo(dphase + halfperiod, period) - halfperiod
        eigenphases(iphase, ie) = prev + dphase
        prev = eigenphases(iphase, ie)
      enddo
    enddo
  end subroutine unwrap_eigenphases

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! pure elemental function ylm_r2c_cs(mr,mc) result(res)
  !   !! Determine the coefficient for the complex complex spherical harmonic for expressing the
  !   !! real spherical harmonics in terms of the comple spherical harmonics, except that in Cs symmetry,
  !   !! we work with linear combinations that are even or odd under mirror plane reflection, so it's
  !   !! actually something like
  !   !!   \(Y_{l}^m = c_1 X_{l|m|}^{even} + c_2 X_{l-|m|}^{odd}\),
  !   !!   \(Y_{l}^m = c_1 X_{l|m|}^{even} - c_2 X_{l-|m|}^{odd}\),
  !   !! UKRmol+ uses the convention of having the YZ plane be the mirror plane for Cs symmetry; this affects
  !   !! which values of m belong to the irreps A' (even) and A'' (odd). For the YZ plane,
  !   !! $x \to -x$ corresponds to $\phi \to \pi-\phi$, under which cos is odd (A'') and sin is even (A').
  !   !! However, instead of assuming that cos is odd, this routine will just query some symmetry arrays
  !   !! that were determined when reading channel files
  !   use rotex__kinds,    only: dp
  !   use rotex__functions,  only: neg, iseven
  !   use rotex__system,   only: die
  !   use rotex__symmetry, only: m_parity, even, odd
  !   implicit none (type, external)
  !   integer, intent(in) :: mr
  !     !! The order (m) for the real spherical harmonic
  !   integer, intent(in) :: mc
  !     !! The ordedr (±|m|) for the complex spherical harmonic
  !   complex(dp), parameter :: one = cmplx(1, 0, kind=dp)
  !   real(dp),    parameter :: invsq2 = 1.0_dp/sqrt(real(2, kind=dp))
  !   complex(dp), parameter :: im  = cmplx(0, 1, kind=dp)
  !   logical :: is_sin
  !   complex(dp) :: coef1, coef2
  !   complex(dp) :: res
  !   if(allocated(m_parity) .eqv. .false.) call die("M_PARITY array is not allocated, but is needed")
  !   res = 0
  !   if(abs(mc) .ne. abs(mr)) return
  !   ! -- the ±i term goes to sin, but sin is not necessarily determined by the sign of mr
  !   is_sin = mr .lt. 0
  !   ! is_sin = iseven(mr) .eqv. (m_parity(mr) .eq. odd)
  !   select case(mc)
  !   case(:-1)
  !     coef1 = 1
  !     coef2 = merge(-im, one, is_sin)
  !   case(1:)
  !     coef1 = neg(mc)
  !     coef2 = merge(im, one, is_sin)
  !   case(0)
  !     res = 1
  !     return
  !   end select
  !   res = coef1 * coef2 * invsq2
  ! end function ylm_r2c_cs


  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! pure elemental function ylm_c2r(mr, mc) result(res)
  !   !! Determine the coefficient of the complex spherical harmonics for expressing the
  !   !! real spherical harmonics in terms of the complex spherical harmonics, e.g.,
  !   !!   \(X_{lm} = c_1 Y_l^{|m|} + c_2Y_l^{-|m|}\),
  !   !! given mc=m and one of mr=|m| or mr=-|m|. Assumes that the degree (l) of the spherical harmonics
  !   !! is the same
  !   use rotex__kinds, only: dp
  !   use rotex__functions, only: neg
  !   implicit none (type, external)
  !   integer, intent(in) :: mr
  !     !! The order (m) for the real spherical harmonic
  !   integer, intent(in) :: mc
  !     !! The order (±|m|) for the complex spherical harmonic
  !   complex(dp), parameter :: one = cmplx(1, 0, kind=dp)
  !   real(dp),    parameter :: invsq2 = 1.0_dp/sqrt(real(2, kind=dp))
  !   complex(dp), parameter :: im  = cmplx(0, 1, kind=dp)
  !   complex(dp) :: coef1, coef2
  !   complex(dp) :: res
  !   res = 0
  !   if(abs(mc) .ne. abs(mr)) return
  !   select case(mr)
  !   case(:-1)
  !     coef1 = im
  !     coef2 = merge(neg(mc+1), 1, mc .gt. 0)
  !   case(1:)
  !     coef1 = 1
  !     coef2 = merge(neg(mc), 1, mc .gt. 0)
  !   case(0)
  !     res = 1
  !     return
  !   end select
  !   res = coef1 * coef2 * invsq2
  ! end function ylm_c2r

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! pure elemental function ylm_c2r_cs(mr,mc) result(res)
  !   !! Determine the coefficient for the complex complex spherical harmonic for expressing the
  !   !! real spherical harmonics in terms of the comple spherical harmonics, except that in Cs symmetry,
  !   !! we work with linear combinations that are even or odd under mirror plane reflection, so it's
  !   !! actually something like
  !   !!   \(X_{lm}^{even} = c_1 Y_l^{|m|} + c_2Y_l^{-|m|}\),
  !   !!   \(X_{lm}^{odd}  = c_1 Y_l^{|m|} - c_2Y_l^{-|m|}\),
  !   !! UKRmol+ uses the convention of having the YZ plane be the mirror plane for Cs symmetry; this affects
  !   !! which values of m belong to the irreps A' (even) and A'' (odd). For the YZ plane,
  !   !! $x \to -x$ corresponds to $\phi \to \pi-\phi$, under which cos is odd (A'') and sin is even (A').
  !   !! However, instead of assuming that cos is odd, this routine will just query some symmetry arrays
  !   !! that were determined when reading channel files
  !   use rotex__kinds,    only: dp
  !   use rotex__system,   only: die
  !   use rotex__symmetry, only: m_parity, even, odd
  !   use rotex__functions, only: neg
  !   implicit none (type, external)
  !   integer, intent(in) :: mr
  !     !! The order (m) for the real spherical harmonic
  !   integer, intent(in) :: mc
  !     !! The ordedr (±|m|) for the complex spherical harmonic
  !   complex(dp), parameter :: one = cmplx(1, 0, kind=dp)
  !   real(dp),    parameter :: invsq2 = 1.0_dp/sqrt(real(2, kind=dp))
  !   complex(dp), parameter :: im  = cmplx(0, 1, kind=dp)
  !   complex(dp) :: coef1, coef2
  !   complex(dp) :: res
  !   if(allocated(m_parity) .eqv. .false.) call die("M_PARITY array is not allocated, but is needed")
  !   res = 0
  !   if(abs(mc) .ne. abs(mr)) return
  !   select case(m_parity(mr))
  !   case(odd)
  !     coef1 = -im
  !     coef2 = merge(1, neg(mc+1), mc .gt. 0)
  !   case(even)
  !     coef1 = 1
  !     coef2 = merge(1, neg(mc), mc .gt. 0)
  !   case(0)
  !     res = 1
  !     return
  !   end select
  !   res = coef1 * coef2 * invsq2
  ! end function ylm_c2r_cs


! ================================================================================================================================ !
end module rotex__RFT
! ================================================================================================================================ !
