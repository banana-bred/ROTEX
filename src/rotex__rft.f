! ================================================================================================================================ !
module rotex__RFT
  !! Procedures used to carry out the rotational frame transformation
  use rotex__globals, only: G
  use rotex__kinds,   only: dp
  use rotex__system,  only: stderr, die, stdout

  implicit none (type, external)

  private

  public :: do_eirft
  ! public :: do_edrft
  public :: do_edrft_chunk
  public :: K2sincos
  public :: K2S_cayley

  public :: real2complex_ylm

  interface real2complex_ylm
    module procedure :: real2complex_ylm_r
    module procedure :: real2complex_ylm_c
  end interface real2complex_ylm

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine do_eirft(           &
        smat_elec                &
      , smat_rot                 &
      , J                        &
      , N_states                 &
      , elec_channels            &
      , rot_channels &
    )
    !! Perform the Energy Independent Rotational Frame Transformation on the electronic S-matrix

    use rotex__types,      only: r3carr_type, n_states_type, asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type&
                               , elec_channel_type
    use rotex__arrays,     only: realloc, is_symmetric, is_unitary, adjoint, unitary_defect
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    complex(dp),                             intent(in) :: smat_elec(:,:)
      !! Electronic S-matrix: n×n
    complex(dp),                      intent(out) :: smat_rot(:,:)
      !! Array of S^J sub blocks
    integer,                                 intent(in) :: J
      !! The current value of J
    type(n_states_type),                     intent(in) :: N_states(:)
      !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
    type(elec_channel_type),                 intent(in) :: elec_channels(:)
      !! The array of electronic channels (n, l, ml), needed as input for the RFT
    type(asymtop_rot_channel_l_type), intent(out):: rot_channels(:)
      !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblock for this J

    logical :: bad_unitarity
    integer :: Omega, nchans_elec, nchans_rot
    complex(dp), allocatable :: U(:,:)
    ! complex(dp), allocatable :: C(:,:)

    bad_unitarity = .false.
    nchans_elec = size(elec_channels, 1)
    nchans_rot = size(rot_channels, 1)

    call realloc(U,        nchans_rot, nchans_elec)
    smat_rot = (0.0_dp, 0.0_dp)

    do Omega = -J, J
      call build_U_RFT(J, Omega, N_states, elec_channels, rot_channels, U)
      ! smat_rot = smat_rot + matmul(U, matmul(smat_elec, adjoint(U)))
      smat_rot = smat_rot + matmul(U, matmul(smat_elec, adjoint(U)))
      ! C = C + matmul(u, adjoint(u))
    enddo

    if(.not. is_symmetric(smat_rot)) call die("The S-matrix is not symmetric after the energy-independent RFT !")
    if(is_unitary(smat_rot)) return

    bad_unitarity = .true.
    write(stderr, '("WARN: The S-matrix for J = ", I0, " is nonunitary with unitary defect ", ES12.4)') &
      J, unitary_defect(smat_rot)

    write(stderr, '("Unitary defect in S: ", ES20.12)') unitary_defect(smat_rot)
    call die("J = "//i2c(J)//"subblock of the S-matrix is non-unitary after the energy-independent RFT. Aborting !")

  end subroutine do_eirft

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  !module subroutine do_edrft(           &
  !      sin_elec                 &
  !    , cos_elec                 &
  !    , kmat_eval_energies       &
  !    , smat_rot_flat            &
  !    , Jmin, Jmax, J            &
  !    , N_states                 &
  !    , elec_channels            &
  !    , rot_channels &
  !  )
  !  !! Perform the Energy Dependent Rotational Frame Transformation on the electronic Sine and Cosine matrices

  !  use rotex__types,  only: r3carr_type, N_states_type, elec_channel_type, &
  !                           asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type
  !  use rotex__arrays, only: realloc, adjoint, is_unitary, unitary_defect, is_symmetric &
  !                         , interp_array_at_energy, packmat, size_check
!#ifdef USE_FORBEAR
  !  use rotex__progress,   only: progressbar_type
!#endif
!#ifdef USE_OPENMP
  !  use omp_lib, only: omp_get_thread_num
!#endif

  !  implicit none (type, external)

  !  complex(dp),                              intent(in)    :: sin_elec(:,:,:), cos_elec(:,:,:)
  !    !! Electronic sine and cosine matrices: n×n×nE
  !  real(dp),                                 intent(in)    :: Kmat_eval_energies(:)
  !    !! Array of K-matrix evaluation energies (therefore the evaluation energies of the other
  !    !! electronic matrices as well)
  !  complex(dp),                       intent(out) :: smat_rot_flat(:,:)
  !    !! The flattened n(n+1)/2 × nE S-matrix subblock for this J
  !  integer,                                  intent(in)    :: Jmin, Jmax, J
  !    !! The min, max, and current values of J
  !  type(N_states_type),                      intent(in)    :: N_states(:)
  !    !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
  !  type(elec_channel_type),                  intent(in)    :: elec_channels(:)
  !    !! The array of electronic channels (n, l, ml), needed as input for the RFT
  !  type(asymtop_rot_channel_l_type),  intent(in)    :: rot_channels(:)
  !    !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblocks at each J

  !  logical :: bad_unitarity, bad_symmetry
  !  ! logical :: did_symmetrize
  !  logical :: bad_unitarity_local, bad_symmetry_local
  !  ! logical :: did_symmetrize_local
  !  integer :: Omega, ie, jrot, ne, nchans_elec, nchans_rot
  !  integer :: ithread
  !  real(dp) :: E_elec_rhs
  !  complex(dp), allocatable :: U(:,:), Udagg(:,:), sin_rot(:,:), cos_rot(:,:)
  !  complex(dp), allocatable :: sinE(:,:), cosE(:,:), smat_rot(:,:)
  !  character(27) :: prefix_string

!#ifdef USE_FORBEAR
  !  integer  :: iprogress, iprogress_last
  !  real(dp) :: rprogress, rprogress_inc
  !  type(progressbar_type) progressbar
!#endif

  !  ne          = size(kmat_eval_energies, 1)
  !  nchans_rot  = size(rot_channels,       1)
  !  nchans_elec = size(elec_channels,      1)

  !  call size_check(sin_elec, [nchans_elec, nchans_elec, ne], "SIN ELEC")
  !  call size_check(cos_elec, [nchans_elec, nchans_elec, ne], "COS ELEC")
  !  call size_check(smat_rot_flat, [(nchans_rot*(nchans_rot+1))/2, ne], "SMAT ROT")

  !  write(prefix_string, '("EDRFT", 2X, 2(I5," /"),I5,X)') Jmin, J, Jmax

!#ifdef USE_FORBEAR
  !  call progressbar % initialize( &
  !      filled_char_string = "|" &
  !    , empty_char_string = " " &
  !    , bracket_left_string = "[" &
  !    , prefix_string = prefix_string &
  !    , suffix_string = "] " &
  !    , add_progress_percent = .true. &
  !  )
  !  call progressbar % start
  !  call progressbar % update(current = 0.0_dp)
  !  rprogress = 0.0_dp
  !  iprogress_last = 0
  !  rprogress_inc = 1.0_dp / real(ne, kind=dp)
!#else
  !  write(stdout, '(A)') prefix_string//".."
!#endif

  !  bad_unitarity  = .false.
  !  bad_symmetry   = .false.
  !  ! did_symmetrize = .false.

  !  !$omp parallel default(none) &
  !  !$omp& shared(ne, J, nchans_rot, nchans_elec, N_states, elec_channels, rot_channels&
  !  !$omp&   , kmat_eval_energies, sin_elec, cos_elec, bad_unitarity, bad_symmetry &
  !  !$omp&   , smat_rot_flat &
!#ifdef USE_FORBEAR
  !  !$omp&   , progressbar, rprogress, rprogress_inc, iprogress, iprogress_last) &
!#else
  !  !$omp& ) &
!#endif
  !  !$omp& private(U, Udagg, sin_rot, sinE, cos_rot, cosE, ie, Omega &
  !  !$omp&   , bad_symmetry_local, bad_unitarity_local, E_elec_RHS, jrot, ithread, smat_rot)
!#ifdef USE_OPENMP
  !  ithread = omp_get_thread_num()
!#else
  !  ithread = 0
!#endif

  !  bad_symmetry_local   = .false.
  !  bad_unitarity_local  = .false.
  !  ! did_symmetrize_local = .false.

  !  call realloc(U,        nchans_rot,  nchans_elec)
  !  call realloc(Udagg,    nchans_elec, nchans_rot)
  !  call realloc(sin_rot,  nchans_rot,  nchans_rot)
  !  call realloc(cos_rot,  nchans_rot,  nchans_rot)
  !  call realloc(smat_rot, nchans_rot,  nchans_rot)
  !  call realloc(sinE,     nchans_elec, nchans_elec)
  !  call realloc(cosE,     nchans_elec, nchans_elec)

  !  !$omp do schedule(static)
  !  nrg: do ie=1, ne

  !    sin_rot  = (0.0_dp, 0.0_dp)
  !    cos_rot  = (0.0_dp, 0.0_dp)

  !    do Omega = -J, J

  !      call build_U_RFT(J, Omega, N_states, elec_channels, rot_channels, U)
  !      Udagg = adjoint(U)

  !      ! -- loop over RHS channels because we need to interpolate/pick matrices at a specific energy
  !      do jrot = 1, nchans_rot

  !        ! -- take the evaluation energy of the right-hand channel, and interpolate the
  !        !    matrices to get the corresponding value
  !        E_elec_rhs = kmat_eval_energies(ie) - rot_channels(jrot) % E
  !        call interp_array_at_energy(E_elec_rhs, kmat_eval_energies, sin_elec, sinE)
  !        call interp_array_at_energy(E_elec_rhs, kmat_eval_energies, cos_elec, cosE)

  !        ! -- U {sin(E),cos(E)} U⁺
  !        sin_rot(:, jrot) = sin_rot(:, jrot) + matmul( U, matmul(sinE, Udagg(:, jrot)) )
  !        cos_rot(:, jrot) = cos_rot(:, jrot) + matmul( U, matmul(cosE, Udagg(:, jrot)) )
  !      enddo

  !    enddo

  !    ! call sincos2S(sin_rot, cos_rot, smat_rot, did_symmetrize_local)
  !    call sincos2S(sin_rot, cos_rot, smat_rot)
  !    call packmat(smat_rot, smat_rot_flat(:, ie), "U")

!#ifdef USE_FORBEAR
  !    !$omp atomic
  !    rprogress = rprogress + rprogress_inc
  !    update_progress: if(ithread .eq. 0) then
  !      iprogress = floor(rprogress * 100)
  !      if(iprogress .eq. iprogress_last) exit update_progress
  !      if(iprogress .eq. 100) exit update_progress
  !      iprogress_last = iprogress
  !      call progressbar % update(current = rprogress)
  !    endif update_progress
!#endif

!! print*, ie, unitary_defect(smat_rot(:,:,ie))
  !    if(is_unitary(smat_rot) .eqv. .false.) then
  !      bad_unitarity_local = .true.
  !      write(stderr, '("WARN: The EDFT S-matrix for J= ", I0, ", ie = ", I0, ", is nonunitary  with&
  !        & unitary defect ", ES12.4)') J, ie, unitary_defect(smat_rot)
  !    endif

  !    if(is_symmetric(smat_rot) .eqv. .false.) then
  !      bad_symmetry_local = .true.
  !      write(stderr, '("WARN: The EDFT S-matrix for J = ", I0, ", ie = ", I0, " is not symmetric")') J, ie
  !    endif

  !  enddo nrg
  !  !$omp end do
  !  ! -- collect issues across threads
  !  !$omp critical
  !    bad_unitarity  = bad_unitarity  .OR. bad_unitarity_local
  !    bad_symmetry   = bad_symmetry   .OR. bad_symmetry_local
  !  !$omp end critical
  !  !$omp end parallel

!#ifdef USE_FORBEAR
  !  call progressbar % update(current = 1.0_dp)
!#else
  !  write(stdout, '(" done !")')
!#endif

  !  if(bad_unitarity .OR. bad_symmetry) call die("At least one J-block of the EDFT S-matrix is&
  !    & non-unitary or non-symmetric after asymmetrization !")

  !end subroutine do_edrft

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine do_edrft_chunk(     &
        sin_elec                 &
      , cos_elec                 &
      , kmat_eval_energies       &
      , ie0, ie1                 &
      , smat_rot_flat_chunk      &
      , J                        &
      , N_states                 &
      , elec_channels            &
      , rot_channels             &
    )
    !! Perform the Energy Dependent Rotational Frame Transformation on the electronic Sine and Cosine matrices

    use rotex__types,  only: N_states_type, elec_channel_type, asymtop_rot_channel_l_type
    use rotex__arrays, only: realloc, adjoint, is_unitary, unitary_defect, is_symmetric &
                           , interp_array_at_energy, packmat, size_check

    implicit none (type, external)

    complex(dp),                      intent(in) :: sin_elec(:,:,:), cos_elec(:,:,:)
      !! Electronic sine and cosine matrices: n×n×nE
    complex(dp),                      intent(out), contiguous :: smat_rot_flat_chunk(:,:)
      !! The flattened n(n+1)/2 × nE S-matrix subblock for this J
    real(dp),                         intent(in) :: kmat_eval_energies(:)
      !! Grid of electronic matrix evaluation energies
    integer,                          intent(in)  :: ie0, ie1
      !! Chunk energy grid bounds
    integer,                          intent(in)  :: J
      !! The current value of J
    type(N_states_type),              intent(in)  :: N_states(:)
      !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
    type(elec_channel_type),          intent(in)  :: elec_channels(:)
      !! The array of electronic channels (n, l, ml), needed as input for the RFT
    type(asymtop_rot_channel_l_type), intent(in)  :: rot_channels(:)
      !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblocks at each J

    logical :: bad_unitarity, bad_symmetry
    logical :: bad_unitarity_local, bad_symmetry_local
    integer :: Omega, ie, ie_loc, jrot, ne, nchans_elec, nchans_rot
    integer :: ne_chunk, nchans_rot_flat
    real(dp) :: E_elec_rhs
    complex(dp), allocatable :: U(:,:), Udagg(:,:), sin_rot(:,:), cos_rot(:,:)
    complex(dp), allocatable :: sinE(:,:), cosE(:,:), smat_rot(:,:)

    ne          = size(kmat_eval_energies, 1)
    ne_chunk    = ie1 - ie0 + 1
    nchans_rot  = size(rot_channels,       1)
    nchans_rot_flat = (nchans_rot*(nchans_rot+1))/2
    nchans_elec = size(elec_channels,      1)

    ! -- bound checks on energy grid
    if(ie0 .lt. 1)   call die("Cannot process a chunk with lower bound index IE0 < 1")
    if(ie1 .gt. ne)  call die("Cannot process a chunk with upper bound index IE1 > NE")
    if(ie0 .gt. ie1) call die("Cannoe process a chunk with IE1 < IE0 (ubound < lbound)")

    ! -- size checks on arrays
    call size_check(sin_elec,            [nchans_elec,     nchans_elec, ne], "SIN_ELEC")
    call size_check(cos_elec,            [nchans_elec,     nchans_elec, ne], "COS_ELEC")
    call size_check(smat_rot_flat_chunk, [nchans_rot_flat, ne_chunk],        "SMAT_ROT_FLAT_CHUNK")

    bad_unitarity  = .false.
    bad_symmetry   = .false.

    !$omp parallel default(none) &
    !$omp& shared(ne, J, nchans_rot, nchans_elec, N_states, elec_channels, rot_channels&
    !$omp&   , kmat_eval_energies, sin_elec, cos_elec, bad_unitarity, bad_symmetry &
    !$omp&   , smat_rot_flat_chunk, ie0, ie1) &
    !$omp& private(U, Udagg, sin_rot, sinE, cos_rot, cosE, ie, ie_loc, Omega &
    !$omp&   , bad_symmetry_local, bad_unitarity_local, E_elec_RHS, jrot, smat_rot)

    bad_symmetry_local   = .false.
    bad_unitarity_local  = .false.

    call realloc(U,        nchans_rot,  nchans_elec)
    call realloc(Udagg,    nchans_elec, nchans_rot)
    call realloc(sin_rot,  nchans_rot,  nchans_rot)
    call realloc(cos_rot,  nchans_rot,  nchans_rot)
    call realloc(smat_rot, nchans_rot,  nchans_rot)
    call realloc(sinE,     nchans_elec, nchans_elec)
    call realloc(cosE,     nchans_elec, nchans_elec)

    !$omp do schedule(static)
    !TODO: flattened electronic matrices here and before, take up less space
    nrg: do ie=ie0, ie1

      ie_loc = ie - ie0 + 1

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
          call interp_array_at_energy(E_elec_rhs, kmat_eval_energies, sin_elec, sinE, .true., .true.)
          call interp_array_at_energy(E_elec_rhs, kmat_eval_energies, cos_elec, cosE, .true., .true.)

          ! -- U {sin(E),cos(E)} U⁺
          sin_rot(:, jrot) = sin_rot(:, jrot) + matmul( U, matmul(sinE, Udagg(:, jrot)) )
          cos_rot(:, jrot) = cos_rot(:, jrot) + matmul( U, matmul(cosE, Udagg(:, jrot)) )
        enddo

      enddo

      call sincos2S(sin_rot, cos_rot, smat_rot)
      call packmat(smat_rot, smat_rot_flat_chunk(:, ie_loc), "U")

      if(is_unitary(smat_rot) .eqv. .false.) then
        bad_unitarity_local = .true.
        write(stderr, '("WARN: The EDFT S-matrix for J= ", I0, ", ie = ", I0, ", is nonunitary  with&
          & unitary defect ", ES12.4)') J, ie, unitary_defect(smat_rot)
      endif

      if(is_symmetric(smat_rot) .eqv. .false.) then
        bad_symmetry_local = .true.
        write(stderr, '("WARN: The EDFT S-matrix for J = ", I0, ", ie = ", I0, " is not symmetric")') J, ie
      endif

    enddo nrg
    !$omp end do
    ! -- collect issues across threads
    !$omp critical
      bad_unitarity  = bad_unitarity  .OR. bad_unitarity_local
      bad_symmetry   = bad_symmetry   .OR. bad_symmetry_local
    !$omp end critical
    !$omp end parallel

    if(bad_unitarity .OR. bad_symmetry) call die("At least one J-block of the EDFT S-matrix is&
      & non-unitary or non-symmetric after asymmetrization !")

  end subroutine do_edrft_chunk

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
    integer :: Ni, Kai, Kci, li, neleci, irot, Ki
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
  pure function get_eigenvector_permutation_idx(U1, U0) result(res)
    !! For eigenvectors U1 and U0 computed at energy E1 and E0 with eigenphases EIGENPHASES1 and EIGENPHASES0,
    !! return the index permutation of matrices 0 such that their inner products with 1 are the best.
    !! This is essentially trying to find the "best" permutation of eigenvectors to get the "same"
    !! eigenphases in the same order across energies
    implicit none (type, external)
    real(dp), intent(in) :: U1(:,:)
      !! Eigenvectors for ie-1
    real(dp), intent(in) :: U0(:,:)
      !! Eigenvectors for ie

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
    use rotex__globals,    only: spinmult_names
    use rotex__constants,  only: au2ev
    use rotex__writing,    only: write_elec_mat_elems_to_file

    implicit none (type, external)

    real(dp),                intent(in)  :: kmat(:,:,:)
    real(dp),                intent(in)  :: kmat_eval_energies(:)
    real(dp),                intent(out) :: sine(:,:,:), cosine(:,:,:)
    type(elec_channel_type), intent(in)  :: elec_channels(:)
    integer,                 intent(in)  :: spinmult
      !! The current spin multiplicity

    logical :: bad_sym, bad_sym_local
    integer :: ne, n, i, ie
    integer :: ithread
    integer, allocatable :: idx(:)
    real(dp), allocatable :: U(:,:,:), eigenphases(:,:), eigenphases_unwrapped(:,:)
    real(dp), allocatable :: tmp(:,:), Udagg(:,:)

    ! -- lapack variables
    integer :: info, lwork
    real(dp), allocatable :: work(:)
    character(1), parameter :: UPLO = 'U'

    n  = size(kmat, 1)
    ne = size(kmat, 3)

    call size_check(kmat_eval_energies, ne,         "KMAT_EVAL_ENERGIES")
    call size_check(sine,               [n, n, ne], "SINE")
    call size_check(cosine,             [n, n, ne], "COSINE")
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
      idx = get_eigenvector_permutation_idx(U(:, :, ie-1), U(:, :, ie))
      U(:,:,ie) = U(:,idx,ie)
      eigenphases(:, ie) = eigenphases(idx, ie)

    enddo

    ! -- 3) correct jumps by ±π
    allocate(eigenphases_unwrapped, source=eigenphases)
    call unwrap_eigenphases(eigenphases_unwrapped)

    bad_sym  = .false.

    ! -- 4) construct the sine and cosine matrices
    !$omp parallel default(none) &
    !$omp& shared(ne, n, sine, cosine, U, eigenphases_unwrapped, kmat_eval_energies, bad_sym )&
    !$omp& private(ie, i, ithread, tmp, Udagg, bad_sym_local)

    allocate(tmp(n,n))
    allocate(Udagg(n,n))
    bad_sym_local = .false.

    !$omp do schedule(static)
    do ie=1,ne

      sine(:,:,ie) = 0._dp
      cosine(:,:,ie) = 0._dp
      do concurrent (i=1:n)
        sine(i,i,ie)   = sin(eigenphases_unwrapped(i, ie))
        cosine(i,i,ie) = cos(eigenphases_unwrapped(i, ie))
      enddo

      ! -- Eigenphases δ -> sin(δ), cos(δ)
      Udagg = adjoint(U(:,:,ie))
      tmp = matmul(sine(:,:,ie), Udagg)
      sine(:,:,ie)   = matmul(U(:,:,ie), tmp)
      tmp = matmul(cosine(:,:,ie), Udagg)
      cosine(:,:,ie) = matmul(U(:,:,ie), tmp)

      if(is_symmetric(sine(:,:,ie)) .AND. is_symmetric(cosine(:,:,ie))) cycle

      ! -- error
      bad_sym_local = .true.

      !$omp critical
      bad_sym = bad_sym .OR. bad_sym_local
      !$omp end critical

      write(stderr, '("The Sine/Cosine matrices are not symmetric for energy ", I0, ": ", E15.7, " eV")') &
        kmat_eval_energies(ie)*au2ev

    enddo
    !$omp end do
    !$omp end parallel

    if (bad_sym) call die("Non-symmetric electronic sine/cosine matrix detected")

    ! -- 5) write to file
    write(stdout, '("Writing sin(E), cos(E), eigenphases(E) to disk..")')
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

  !    nrg: do ie=1, ne


  !      ! -- loop over symmetries
  !      ! do isym=1, nsyms

  !      ! calculations are in the full symmetry. expand to full symmetry beforehand if needed
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
  subroutine sincos2s(sinec, cosinec, s)
    !! Sine,Cosine -> K=SC⁻¹ -> S
    use rotex__system, only: die
    use rotex__arrays, only: is_symmetric
    use rotex__linalg, only: right_divide
    implicit none (type, external)
    complex(dp), intent(in)  :: sinec(:,:), cosinec(:,:)
    complex(dp), intent(out) :: s(:,:)
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
