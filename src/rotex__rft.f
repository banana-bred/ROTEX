! ================================================================================================================================ !
module rotex__RFT
  !! Procedures used to carry out the rotational frame transformation
  use rotex__globals, only: G
  use rotex__kinds, only: dp

  implicit none (type, external)

  private

  ! public :: RFT_linear
  public :: RFT_nonlinear

  interface real2complex_ylm
    ! module procedure :: real2complex_ylm_r
    module procedure :: real2complex_ylm_c
  end interface real2complex_ylm

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine RFT_nonlinear( Kmat_flat                &
                                 , Jmin, Jmax               &
                                 , Smat_J                   &
                                 , elec_channels            &
                                 , N_states                 &
                                 , asymtop_rot_channels_l   &
                                 , asymtop_rot_channels_l_J &
                                 )
    !! Build the electronic S-matrix from the electronic K-matrix, then perform the rotational frame transformation on the S-matrix

    use rotex__kinds,     only: dp
    use rotex__types,     only: elec_channel_type, asymtop_rot_channel_l_type, N_states_type, cmatrix_type &
                              , asymtop_rot_channel_l_vector_type
    use rotex__channel_ops, only: sort_channels_by_energy
    use rotex__system,    only: stdout, die
    use rotex__symmetry,  only: possible_spin_symmetries, spin_symmetry
    use rotex__constants, only: im
    use rotex__arrays,    only: append, is_unitary

    implicit none (type, external)

    real(dp),                       intent(inout), allocatable :: Kmat_flat(:,:)
      !! The flattened K-matrix, needed as input for the RFT. nchans(nchans+1)/2 × nE
    integer, intent(in) :: Jmin
      !!  The lowest value of J = N + l
    integer, intent(in) :: Jmax
      !!  The largest value of J = N + l
    type(cmatrix_type), intent(out) :: Smat_J(Jmin:Jmax)
      !! The rotational S-matrices \(S^J\), produced by the RFT
    type(elec_channel_type),        intent(in)                 :: elec_channels(:)
      !! The array of electronic channels (n, l, ml), needed as input for the RFT
    type(N_states_type),            intent(in)                 :: N_states(:)
      !! The array of rotational states of the target (N, Ka, Kc), needed as input for the RFT
    type(asymtop_rot_channel_l_type), intent(out),   allocatable :: asymtop_rot_channels_l(:)
      !! The array of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix
    type(asymtop_rot_channel_l_vector_type), intent(out) :: asymtop_rot_channels_l_J(Jmin:Jmax)
      !! The array of arrays of rotational channels (N, Ka, Kc, l) that make up the basis of the S-matrix subblocks at each J

    integer :: nelec, l
    integer :: ml
    integer :: ichan, nchans_rot, nchans_elec
    complex(dp), allocatable :: Smat_elec(:,:)

    nchans_elec = size(elec_channels, 1)

    call build_rotational_channels( n_states         &
                                  , elec_channels    &
                                  , asymtop_rot_channels_l)

    ! -- channels are sorted by contsruction, but we should still sort them by energy here in case
    !    anything above changes in the future
    call sort_channels_by_energy(asymtop_rot_channels_l)

    nchans_rot = size(asymtop_rot_channels_l, 1)
    if(nchans_rot .lt. 1) call die("Number of rotational channels must not be less than 1 !")

    allocate(smat_elec(nchans_elec, nchans_elec))
    smat_elec = 0
    ! call K2S(Kmat, smat_elec, elec_channels, real_spherical_harmonics, point_group)

    if(G%EDFT) then
      ! -- correct for branch cuts to get a smooth S-matrix
      @@@@@
      call K2S_smooth(kmat_flat, smat_elec_flat, elec_channels)
    else
      ! -- just do the Cayley transform
      call K2S_cayley(kmat_flat, smat_elec_flat, elec_channels)
    endif

    write(stdout, '(A)') "Channel-by-channel unitarity of the electronic S-matrix:"
    write(stdout, '(8X, 4A4, A15)') "i", "n", "l", "ml", "||S(:,i)||₂"
    do ichan=1, nchans_elec
      nelec = elec_channels(ichan) % nelec
      l     = elec_channels(ichan) % l
      ml    = elec_channels(ichan) % ml
      write(stdout, '(8X, 4I4, 2X, F8.4)') ichan, nelec, l, ml, sum(abs(Smat_elec(ichan,:))**2)
    enddo
    write(stdout, *)
    if(is_unitary(Smat_elec) .eqv. .false.) call die("Electronic S-matrix is not unitary !")

    ! -- Smat_elec -> Smat (RFT)
    ! (n) l λ -> N τ(Ka,Kc) l

    write(stdout, '(A)') "Frame transformation: S_elec -> S^J"

    call do_rft(                 &
        smat_elec                &
      , smat_j                   &
      , jmin                     &
      , jmax                     &
      , n_states                 &
      , elec_channels            &
      , asymtop_rot_channels_l   &
      , asymtop_rot_channels_l_j &
    )

  end subroutine RFT_nonlinear

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine build_rotational_channels(n_states, elec_channels, rot_channels)
    !! Build rotational+electronic channels channels: (N Ka Kc)+(l λ) = (N Ka Kc l λ)
    use rotex__kinds,    only: dp
    use rotex__types,    only: n_states_type, elec_channel_type, asymtop_rot_channel_l_type
    use rotex__arrays,   only: append
    use rotex__symmetry, only: spin_symmetry
    use rotex__system,   only: die
    implicit none (type, external)
    type(N_states_type),            intent(in)               :: n_states(:)
    type(elec_channel_type),        intent(in)               :: elec_channels(:)
    type(asymtop_rot_channel_l_type), intent(out), allocatable :: rot_channels(:)
    integer  :: i_N_state, i_tau, i_elec_channel
    integer  :: n, ka, kc, ksym, nelec, iq, sym
    integer  :: l, lprev
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
        lprev = elec_channels(1) % l
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
          ! -- don't worry about the different projections of ml for the final channels
          if(l .eq. lprev .and. l .ne. 0) cycle
          lprev = l
          sym = spin_symmetry(n, ka, kc)
          channel = asymtop_rot_channel_l_type(nelec=nelec, l=l, iq=iq, n=n, ka=ka, kc=kc, e=e, sym=sym)
          call append(rot_channels, channel)
        enddo
      enddo
    enddo
  end subroutine build_rotational_channels

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine K2S_smooth(kmat_flat, smat_flat, elec_channels)
    !! Take a flattened K-matrix that is defined for several energies, and
    !! transform it into the S-matrix. Along the way, eigenphases will be put in the same branch
    !! so as to avoid jumps in S-matrix elements as a function of evaluation energy

    use rotex__types,  only: elec_channel_type
    use rotex__linalg, only: dsyev
    use rotex__arrays, only: nflat2n
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    real(dp),                intent(in)  :: kmat_flat(:,:)
    complex(dp),             intent(out) :: smat_flat(:,:)
    type(elec_channel_type), intent(in)  :: elec_channels(:)

    integer :: ne, nflat, n
    real(dp), allocatable :: U(:,:,:)

    ! -- lapack variables
    integer :: info, lwork
    real(dp), allocatable :: w(:), work(:)
    character(1), parameter :: UPLO = 'U'

    nflat = size(kmat_flat, 1)
    ne = size(kmat_flat, 2)

    n = nflat2n(nchans_flat)

    allocate(U(n,n,ne), source=0._dp)
    lwork = 3*n+1
    allocate(w(n))
    allocate(work(lwork))

    ! -- 1) diagonalize K, get eigenphases, print eigenphases as functin of energy
    do concurrent(ie=1, ne)

      call unpackmat(kmat_flat, U(:,:,ie), UPLO)

      call dsyev('V', UPLO, n, U, n, w, work, lwork, info)
      if(info .ne. 0) call die("DSYEV exited with nonzero INFO = " // i2c(info))

      ! -- δ <- tan(δ)
      w = atan(w)

      @@@@@
      !now do the drip thing. compare inner products, etc etc, permute columns, etc etc,
      ! maybe make some array helpers for that

    enddo


    ! -- 2)  identify eigenvalues by inner product of eigenvectors to permute columns of
    !        the matrix U at each geometry

    ! -- 3) now that we can identify eigenphases across energies, fix their continuity
    !      between branches by jumps of  +/- π

    ! -- 4) maybe print some stuff

  end subroutine K2S_smooth

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine K2S_cayley(kmat_flat, smat_flat, elec_channels)
    !! electronic Kmat -> electronic Smat via Cayley transform. Also ensures that the
    !! S-matrix is in the basis of complex-valued spherical harmonics

    use rotex__types,      only: elec_channel_type
    use rotex__arrays,     only: adjoint, eye, size_check, unpackmat, packmat
    use rotex__constants,  only: im
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    use rotex__linalg,     only: dsyev, zgesv

    implicit none (type, external)

    real(dp),    intent(in)  :: kmat_flat(:,:)
    complex(dp), intent(out) :: smat_flat(:,:)
    type(elec_channel_type), intent(in) :: elec_channels(:)

    character(1), parameter :: jobz = "V"
    character(1), parameter :: uplo = "U"
    integer :: n, info, ne
    integer, allocatable :: ipiv(:)
    real(dp),    allocatable :: I(:,:)
    real(dp),    allocatable :: kmat(:,:)
    complex(dp), allocatable :: A(:,:)

    ! -- array sizes
    n = size(elec_channels,1)
    ne = size(kmat_flat, 2)
    call size_check(Kmat_flat, [n,ne], "KMAT_FLAT")
    call size_check(Smat_flat, [n,ne], "SMAT_FLAT")

    ! -- S = (I + iK) / (I - iK)
    allocate(ipiv(n))
    I    = real(eye(n), kind=dp)

    do concurrent(ie=1:ne)

      ! -- unpack kmat_flat -> Kmat
      call unpackmat(kmat_flat, Kmat)

      A    = cmplx(I(:,:), -Kmat(:,:), kind=dp)
      smat = cmplx(I(:,:),  Kmat(:,:), kind=dp)
      call zgesv(n, n, A, n, ipiv, smat, n, info)

      ! -- transform real-valued Xlm basis to complex-valued Ylm if desired
      if(G%REAL_SPHERICAL_HARMONICS .eqv. .false.) return

      call real2complex_ylm(smat, elec_channels)

      ! -- pack the Smat -> Smat_flat
      call packmat(Smat, smat_flat)

    enddo

  end subroutine K2S_cayley

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine K2S(Kmat, Smat, elec_channels)
    !! electronic Kmat -> electronic Smat
    use rotex__kinds,      only: dp
    use rotex__types,      only: elec_channel_type
    use rotex__arrays,     only: adjoint, eye
    use rotex__constants,  only: im
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    use rotex__linalg,     only: dsyev, zgesv
    implicit none (type, external)
    real(dp),    intent(in)  :: Kmat(:,:)
    complex(dp), intent(out) :: Smat(:,:)
    type(elec_channel_type), intent(in) :: elec_channels(:)
    character(1), parameter :: jobz = "V"
    character(1), parameter :: uplo = "U"
    integer :: ichan
    integer :: nchans_elec
    integer :: nn
    integer :: lda
    integer :: ldwork
    integer :: info
    real(dp), allocatable :: w(:)
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: U(:,:)
    nchans_elec = size(elec_channels,1)
    smat = 0
    ! -- ensure the K-matrix is represented in a basis of complex spherical harmonic
    U = Kmat(:,:)
    ! if(real_spherical_harmonics .eqv. .true.) call real2complex_ylm(U, elec_channels)
    nn = size(U, 1)
    lda = nn
    ldwork = 3*nn-1
    allocate(w(nn))
    allocate(work(ldwork))
    ! -- diagonalize K = U tan(δ) UT
    call dsyev(jobz, uplo, nn, U, lda, w, work, ldwork, info)
    if(info .ne. 0) call die("Error in diagonalizaion of the electronic K-matrix ! INFO return as " // i2c(info))
    ! -- build diagonal e^{2iδ}
    do concurrent (ichan=1:nn)
      smat(ichan,ichan) = exp(2*im*atan(w(ichan)))
    enddo
    ! -- S = U e^{2iδ} UT
    smat = matmul( U, matmul(smat , adjoint(U)) )

    ! -- transform real-valued Ylm basis to complex-valued Ylm if desired
    if(G%REAL_SPHERICAL_HARMONICS .eqv. .false.) return

    call real2complex_ylm(smat, elec_channels)

  end subroutine K2S

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine do_rft(             &
        Smat_elec                &
      , Smat_j                   &
      , jmin                     &
      , jmax                     &
      , n_states                 &
      , elec_channels            &
      , asymtop_rot_channels_l   &
      , asymtop_rot_channels_l_j &
    )
    !! Perform the rotational frame transformation on the electronic S-matrix

    use rotex__kinds, only: dp
    use rotex__types, only: cmatrix_type, asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type &
                          , elec_channel_type, N_states_type
    use rotex__wigner,     only: clebsch
    use rotex__system,     only: stdout, die
    use rotex__arrays,     only: realloc, is_unitary, uniq
    use rotex__functions,  only: neg
    use rotex__characters, only: i2c => int2char
    implicit none (type, external)
    complex(dp),        intent(in)  :: Smat_elec(:,:)
    type(cmatrix_type), intent(out) :: Smat_j(jmin:jmax)
      !! Rotationally resolved S-matrix at each J
    integer, intent(in) :: jmin, jmax
      !! Min/max values of the total angular momentum J
    type(n_states_type), intent(in) :: n_states(:)
    type(elec_channel_type), intent(in) :: elec_channels(:)
    type(asymtop_rot_channel_l_type), intent(in) :: asymtop_rot_channels_l(:)
    type(asymtop_rot_channel_l_vector_type), intent(out) :: asymtop_rot_channels_l_j(jmin:jmax)
    type(asymtop_rot_channel_l_type), allocatable :: rot_channels(:)
    logical :: flag
    logical, allocatable :: mask(:)
    integer :: nsyms, nchans_elec
    integer :: J, i, isym, sym, nchans_sym
    integer :: nchans_J
    integer, allocatable :: idx(:), uniq_syms(:)
    complex(dp), allocatable :: Smat_rot(:,:), Smat_rot_sym(:,:)
    complex(dp), allocatable :: U(:,:)

    flag = .false.
    nchans_elec = size(smat_elec, 1)

    ! -- loop over different values of the agular momentum J
    jloop: do J=Jmin,Jmax

      ! -- determine number of rotational channels in this block of total J
      call collect_J_channels_indices(j, asymtop_rot_channels_l, idx)
      asymtop_rot_channels_l_j(j) % channels = asymtop_rot_channels_l(idx)

      nchans_J = size(asymtop_rot_channels_l_j(j) % channels, 1)

      ! -- the total S-matrix for this J
      call realloc(U,        nchans_J, nchans_elec)
      call realloc(smat_rot, nchans_J, nchans_J)
      U = 0
      smat_rot = 0

      ! -- get the unique symmetry elements. We do one frame transformation per
      !    symmetry element
      uniq_syms = asymtop_rot_channels_l_j(j) % channels % sym
      uniq_syms = uniq(uniq_syms)
      nsyms = size(uniq_syms, 1)

      ! call do_rft_no_sym(J, n_states, elec_channels, asymtop_rot_channels_l_j(j)%channels, smat_elec, smat_rot, U, point_group)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! For symmetry enforcement, we need to consider the
      ! electronic parity too
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! -- loop over symmetries
      do isym=1, nsyms
        sym = uniq_syms(isym)
        ! -- map all J -> this sym
        mask = asymtop_rot_channels_l_j(j)%channels%sym .eq. sym
        ! mask = asymtop_rot_channels_l_j(j)%channels%sym .eq. asymtop_rot_channels_l_j(j)%channels%sym
        idx = pack([(i,i=1,nchans_j)], mask)
        rot_channels = asymtop_rot_channels_l_j(j) % channels(idx)
        ! -- allocate U, Smat_rot for this sym
        nchans_sym = size(idx, 1)
        call realloc(U,            nchans_sym, nchans_elec)
        call realloc(smat_rot_sym, nchans_sym, nchans_sym)
        U = 0
        smat_rot_sym = 0
        call do_rft_this_sym(j, sym, n_states, elec_channels, rot_channels, smat_elec, smat_rot_sym, U)
        ! -- add this contribution back to the total S-matrix for this J
        smat_rot(idx, idx) = smat_rot_sym(:,:)
      enddo

      ! -- export this S^J
      smat_j(j) % mtrx = smat_rot(:,:)

      ! if(is_unitary(Smat_rot) .eqv. .true.) cycle

      ! flag = .true.

      ! ! -- warn about nonunitarity
      ! block
      !   use rotex__system, only: stderr
      !   use rotex__arrays, only: eye, adjoint, norm_frob, unitary_defect
      !   associate(S => Smat_J(J)%mtrx)
      !     write(stderr, '(A, I0, A, F7.5)') &
      !       "WARN: The S-matrix for J = ", J, " is nonunitary with unitary defect ", unitary_defect(S)
      !   end associate
      ! end block

    enddo jloop

    ! if(flag .eqv. .false.) return

    ! call die("At least one J-block of the S-matrix is non-unitary. This may cause some issues in the&
    !   & ensuing MQDT closed-channel elimination procedure which takes the closed channels into account for each J-block.&
    !   & Therefore, each J-block should be unitary, even if they involve states with N< N_min or N > N_max. It is probably&
    !   & worth noting that the S-matrix at this point was detected to be non-unitary, but each symmetry sub-block was unitary.")

  end subroutine do_rft

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! pure subroutine real2complex_ylm_r(M, chans)
  !   !! Real-valued version of the complex-valued equivalent
  !   use rotex__types,  only: dp, elec_channel_type
  !   use rotex__system, only: die
  !   implicit none (type, external)
  !   real(dp), intent(inout) :: M(:,:)
  !     !! The S/K-matrix
  !   type(elec_channel_type), intent(in) :: chans(:)
  !   complex(dp), allocatable :: MC(:,:)
  !   MC = cmplx(M, 0.0_dp, kind = dp)
  !   call real2complex_ylm_c(MC, chans)
  !   if(maxval(abs(MC%im)) .gt. 1e-12) call die("Nonzero imaginary values detected in&
  !     & transformed S/K matrix (complex spherical harmonics basis)")
  !   M = MC % re
  ! end subroutine real2complex_ylm_r

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
    allocate(idx(count))
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

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine do_rft_this_sym(j, sym, n_states, elec_channels, rot_channels, Smat_elec, Smat_rot, U)
    !! Do the rotational frame transformation for a specific symmetry
    use rotex__kinds,      only: dp
    use rotex__types,      only: elec_channel_type, asymtop_rot_channel_l_type, n_states_type
    use rotex__arrays,     only: size_check, is_unitary, is_symmetric, adjoint
    use rotex__wigner,     only: clebsch
    use rotex__system,     only: die, stderr, stdout
    use rotex__functions,  only: neg
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    integer,                          intent(in)  :: j
      !! Total angular momentum quantum number J
    integer,                          intent(in)  :: sym
      !! The current symmetry
    type(n_states_type),              intent(in)  :: n_states(:)
      !! N, Ka, and Kc for each N
    type(elec_channel_type),          intent(in)  :: elec_channels(:)
      !! Electronic channel basis for Smat_elec
    type(asymtop_rot_channel_l_type), intent(in)  :: rot_channels(:)
      !! Rotational channel basis for Smat_rot (this symmetry)
    complex(dp),                      intent(in)  :: smat_elec(:,:)
      !! Electronic S-matrix
    complex(dp),                      intent(out) :: smat_rot(:,:)
      !! Rotatinal S-matrix
    complex(dp),                      intent(inout) :: U(:,:)
      !! Unitary transformation matrix

    integer :: irot
    integer :: nchans_elec, nchans_rot
    integer :: ni, kai, kci, li, nj, kaj, kcj, lj, ki, lambdaj, symchan
    integer :: neleci, nelecj
    integer :: in, itau, iK, jelec
    integer :: Omega
    logical, allocatable :: mask(:)
    complex(dp), allocatable :: C(:,:)

    nchans_rot  = size(rot_channels, 1)
    nchans_elec = size(elec_channels, 1)
    Smat_rot = 0

    allocate(C, source=U) ; C = 0

    ! -- build the rectangular transformation matrix U <LF|BF> for each Ω
    do Omega = -J,J

      U = 0

      do irot = 1, nchans_rot

        call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li, symchan)

        ! if(sym .ne. symchan) call die("Channel symmetry does not match transformation symmetry !")

        ! -- get the corresponding eigenvector for this state, indexed by itau.
        in    = findloc(n_states % n, value = ni, dim = 1)
        select case(G%ROTOR_KIND)
        case("a", "A")
          mask = (n_states(in) % ka(:) .eq. kai) .AND. (n_states(in) % kc(:) .eq. kci)
        case("s", "S")
          select case(G%ROTOR_ZAXIS)
          case("a", "A")
            mask = (n_states(in) % ka(:) .eq. kai)
          case("c", "C")
            mask = (n_states(in) % kc(:) .eq. kci)
          case default
            call die("Somehow got a symmetric top with a G%ROTOR_ZAXIS " // G%ROTOR_ZAXIS // " that is neither A nor C")
          end select
        case default
          call die("ROTOR_KIND " // G%ROTOR_KIND // " not allowed in RFT")
        end select
        itau  = findloc(mask, value = .true., dim = 1)

        do jelec = 1, nchans_elec
          nelecj  = elec_channels(jelec) % nelec
          lj      = elec_channels(jelec) % l
          lambdaj = elec_channels(jelec) % ml

          ! -- enforce transformation between the same electronic state n and partial wave l
          if (neleci .ne. nelecj) cycle
          if (li     .ne. lj) cycle
          Ki  = Omega - lambdaj
          if(abs(Ki) .gt. Ni) cycle

          ik = Ki + Ni + 1

          U(irot, jelec) = neg(lj + lambdaj)              &
              * N_states(in) % eigenH % eigvecs(ik, itau) &
              * clebsch(lj, -lambdaj, J, Omega, Ni, Ki)

        enddo
      enddo

      ! -- S^J = Σ_Ω USU⁺ (for each Ω)
      Smat_rot = Smat_rot + matmul( U, matmul(Smat_elec, adjoint(U)) )
      C = C + matmul(U, adjoint(U))

    enddo


    error_checks: block
      use rotex__utils,  only: printmat
      use rotex__system, only: warn
      use rotex__arrays, only: unitary_defect, eye, norm_frob
      logical :: symflag = .false.
      logical :: unitaryflag = .false.

      if(is_symmetric(Smat_rot) .eqv. .false.) then
        ! -- not symmetric
        symflag = .true.
        call warn("The S-matrix is not symmetric for symmetry " // i2c(sym) // " ❌")
      endif

      if(is_unitary(Smat_rot, 1e-7_dp)   .eqv. .true.) then
        write(stdout, '("S-matrix is unitary for J = ", I0, ", symmetry ", I0, " ✔️")') J, sym
        if(symflag .eqv. .false.) return
      else
        unitaryflag = .true.
      endif


      write(stderr, '("Symmetry: ", I0)') sym
      write(stderr, '("Channels: ", 6(A5,X), A20)') "i", "nelec", "N", "Ka", "Kc", "l", "Σ|S(i,:)|²"
      do irot=1, nchans_rot
        call get_channel_qnums_rot(rot_channels, irot, neleci, ni, kai, kci, li)
        write(stderr, '(10X, 6(I5,X), E20.10)', advance = "no") irot, neleci, ni, kai, kci, li &
          , sum(abs(Smat_rot(irot,:))**2)
        if(all(U(irot,:) .eq. 0._dp)) write(stderr, '(" <-- ", A)', advance = "no") "Does not couple to any electronic channels !"
        write(stderr, *)
      enddo
      write(stderr, *)
      write(stderr, '("This is symmetry ", I0, ", J = ", I0)') sym, J
      write(stderr, '(A30, F15.9)') "Unitary defect in UU⁺: ", unitary_defect(C)
      write(stderr, '(A30, F15.9)') "Unitary defect in USU⁺: ", unitary_defect(Smat_rot)
      if(unitaryflag) call warn("The S-matrix is not unitary for symmetry " // i2c(sym) // " ❌")
      if(symflag) call warn("The S-matrix is not symmetric for symmetry " // i2c(sym) // " ❌")
      if(unitaryflag .or. symflag) error stop
    end block error_checks

  end subroutine do_rft_this_sym

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
