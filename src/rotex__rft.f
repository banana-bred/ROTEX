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
                                 , kmat_eval_energies       &
                                 , spinmult                 &
                                 , Jmin, Jmax               &
                                 , Smat_J_flat              &
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
    real(dp), intent(in) :: kmat_eval_energies(:)
      !! The evaluation energies of the K-matrix
    integer, intent(in) :: spinmult
      !! Current spin multiplicity
    integer, intent(in) :: Jmin
      !!  The lowest value of J = N + l
    integer, intent(in) :: Jmax
      !!  The largest value of J = N + l
    type(cmatrix_type), intent(out) :: Smat_J_flat(Jmin:Jmax,1:size(kmat_flat, 2))
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
    integer :: ichan, nchans_rot, nchans_elec, nchans_elec_flat, ne
    real(dp), allocatable :: sin_elec_flat(:,:), cos_elec_flat(:,:)
      !! nchan(nchan+1)/2 × ne arrays of flattened sine/cosine matrix elements used in the EDFT
    complex(dp), allocatable :: smat_elec_flat(:)
      !! nchan(nchan+1)/2 array of flattened S-matrix elements. Not used in EDFT because
      !! we do that on the sine and cosine matrices

    ne = size(kmat_eval_energies, 1)
    nchans_elec = size(elec_channels, 1)
    nchans_elec_flat = (nchans_elec*(nchans_elec+1))/2

    call build_rotational_channels( n_states         &
                                  , elec_channels    &
                                  , asymtop_rot_channels_l)

    ! -- channels are sorted by contsruction, but we should still sort them by energy here in case
    !    anything above changes in the future
    call sort_channels_by_energy(asymtop_rot_channels_l)

    nchans_rot = size(asymtop_rot_channels_l, 1)
    if(nchans_rot .lt. 1) call die("Number of rotational channels must not be less than 1 !")

    allocate(sin_elec_flat(nchans_elec_flat, ne), source=0._dp)
    allocate(cos_elec_flat(nchans_elec_flat, ne), source=0._dp)

    if(G%EDFT) then
      ! -- correct for branch cuts to get a smooth S-matrix
      call K2sincos(kmat_flat, kmat_eval_energies, sin_elec_flat, cos_elec_flat, elec_channels, spinmult)
    else
      ! -- just do the Cayley transform and get the S-matrix directly
      call K2S_cayley(kmat_flat, smat_elec_flat, elec_channels)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! ROTATIONAL FRAME TRANSFORMATION !!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(stdout, '(A)') "Frame transformation: S_elec -> S^J"
    call do_rft(                 &
        smat_elec_flat           &
      , sin_elec_flat            &
      , cos_elec_flat            &
      , smat_j_flat              &
      , ne                       &
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
    real(dp), intent(in) :: eigenphases1(:,:)
      !! Eigenphases for ie-1
    real(dp), intent(in) :: eigenphases0(:,:)
      !! Eigenphases for ie

    integer, allocatable :: res(:)

    integer :: i, n
    integer, allocatable :: idx(:)
    real(dp), allocatable :: inner(:)

    n = size(U1, 1)

    ! -- NOTE right here we could chekc the signs of the vectors, but we don't have to do that
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
  subroutine K2sincos(kmat_flat, kmat_eval_energies, sin_flat, cos_flat, elec_channels, spinmult)
    !! Take a flattened K-matrix that is defined for several energies, and
    !! transform it into the sine and cosine matrices, interpolated in energy.
    !!   1) Diagonalize K(E) ~ tan(δ(E))
    !!   2) Identify eigenphases across (E)
    !!   3) Build sine/cosine matrices. These will be interpolated linearly later in the RFT and
    !!      closed channel elimination procedure. Sine and Cosine -> K=Sin*Cos⁻¹ -> S
    !!   4) Write eigenphase info to file
    !! Along the way, eigenphases will be put in the same branch
    !! so as to avoid jumps in matrix elements as a function of evaluation energy.
    !! Extrapolation beyond the available evaluation energy grid will be handled in the RFT and cross sections
    !! as needed via simple linear or constant inteprolation
    !! TODO: can be more efficient. use kmat_flat as eigenphases, etc

    use rotex__types,  only: elec_channel_type
    use rotex__linalg, only: dsyev
    use rotex__arrays, only: nflat2n, size_check, is_unitary, unitary_defect, adjoint, is_symmetric, packmat
    use rotex__characters, only: i2c => int2char
    use rotex__constants, only: SPINMULT_NAMES, au2ev, pi

    implicit none (type, external)

    real(dp),                intent(in)  :: kmat_flat(:,:)
    real(dp),                intent(in)  :: kmat_eval_energies(:)
    real(dp),                intent(out) :: sin_flat(:,:), cos_flat(:,:)
    type(elec_channel_type), intent(in)  :: elec_channels(:)
    integer,                 intent(in)  :: spinmult
      !! The current spin multiplicity

    logical :: add_energy_pre, add_energy_post
    integer :: ne, nflat, n, i, j, ie
    integer :: funit_eigenphases, funit_sine, funit_cosine
    integer, allocatable :: idx(:)
    real(dp) :: D1, D2, D3, D4, h1, h2
    real(dp) :: E_pre, E_post, E1, E2
    real(dp), allocatable :: sin_pre(:), sin_post(:), cos_pre(:), cos_post(:)
    real(dp), allocatable :: tmp(:)
    real(dp), allocatable :: eval_E_pre(:), eval_E_post(:)
    real(dp), allocatable :: U(:,:,:), sine(:,:), cosine(:,:), eigenphases(:,:)
    complex(dp), allocatable :: Smat(:,:)
    character(:), allocatable :: eigenphases_dir
    character(:), allocatable :: eigenphases_file, smat_file

    ! -- lapack variables
    integer :: info, lwork
    real(dp), allocatable :: w(:), work(:)
    character(1), parameter :: UPLO = 'U'

    nflat = size(kmat_flat, 1)
    ne = size(kmat_flat, 2)

    n = nflat2n(nflat)

    call size_check(kmat_eval_energies, ne, "KMAT_EVAL_ENERGIES")
    call size_check(sin_flat, [nflat, ne], "SINE_FLAT")
    call size_check(cos_flat, [nflat, ne], "COSINE_FLAT")
    call size_check(elec_channels, n, "ELEC_CHANNELS")

    allocate(U(n,n,ne), source=0._dp)
    allocate(Smat(n, n), source=0._dp)
    lwork = 3*n+1
    allocate(sine(n,n), source=0._dp)
    allocate(cosine(n,n), source=0._dp)
    allocate(eigenphases(n, ne))
    allocate(work(lwork))

    ! -- 1) diagonalize K, get eigenphases, permute them to ensure consistent ordering across geometries
    do ie=1, ne

      call unpackmat(kmat_flat, U(:,:,ie), UPLO)

      call dsyev('V', UPLO, n, U(:,:,ie), n, eigenphases(:, ie), work, lwork, info)
      if(info .ne. 0) call die("DSYEV exited with nonzero INFO = " // i2c(info))

      ! -- δ <- tan(δ)
      eigenphases = atan(eigenphases)

      ! -- need two energies for comparison
      if(ie .eq. 1) cycle

      ! -- compare with the previous energy to identify eigenvectors/eigenphases. This assumes
      !    that the inner product doesn't change too much between energies, so a relatively dense
      !    grid may be necessary near steep resonances. .
      idx = get_eigenvector_permutation_idx(U(:, :, ie-1), U(:, :, ie), eigenphases(:, ie-1), eigenphases(:, ie))
      U(:,:,ie) = U(idx,idx,ie)

    enddo

    ! -- 2)  identify eigenvalues by inner product of eigenvectors to permute columns of
    !        the matrix U at each energy
    do concurrent (i=1:n)
      do ie=1, ne-2

        ! -- minimize the absolute value of the first derivative between the first two points
        if(ie .eq. 1) then
          D1 =       eigenphases(i, ie+1) - eigenphases(i, ie)
          D2 =  pi + eigenphases(i, ie+1) - eigenphases(i, ie)
          D3 = -pi + eigenphases(i, ie+1) - eigenphases(i, ie)
          select case( minloc(abs([D1,D2,D3]), 1) )
          case(1) ; continue
          case(2)
            ! -- shift by +π
            if(sgn(D2) .ne. sgn(D1)) cycle
            eigenphases(i, ie+1) = eigenphases(i, ie+1) + pi
          case(3)
            ! -- shift by -π
            if(sgn(D3) .ne. sgn(D1)) cycle
            eigenphases(i, ie+1) = eigenphases(i, ie+1) - pi
          case default
            call die("Somehow, minloc of a 3-element vector return something other than 1,2,3")
          end select
        endif

        ! -- minimize the difference in first derivatives between points [ie,ie+1] and [ie,ie+2]
        h1 = kmat_eval_energies(ie+1) - kmat_eval_energies(ie)
        h2 = kmat_eval_energies(ie+2) - kmat_eval_energies(ie+1)
        D1 = ( eigenphases(i, ie+1) - eigenphases(i, ie) ) / h1
        do
          D2 = (    eigenphases(i, ie+2) - eigenphases(i, ie+1) ) / h2
          D3 = ( pi+eigenphases(i, ie+2) - eigenphases(i, ie+1) ) / h2
          D4 = (-pi+eigenphases(i, ie+2) - eigenphases(i, ie+1) ) / h2
          select case(minloc(abs([D2-D1,D3-D1,D4-D1]),1))
          case(1)
            exit
          case(2)
            eigenphases(i, ie+2:) = eigenphases(i, ie+2:) + pi
          case(3)
            eigenphases(i, ie+2:) = eigenphases(i, ie+2:) - pi
          case default
            call die("Somehow, minloc of a 3-element vector return something other than 1,2,3")
          end select
        enddo

      enddo
    enddo

    ! -- 3) construct Sine and Cosine matrices
    do ie=1,ne

      sine   = 0._dp
      cosine = 0._dp
      do concurrent (i=1:n)
        sine(i,i)   = sin(eigenphases(i))
        cosine(i,i) = cos(eigenphases(i))
      enddo

      ! -- Eigenphases δ -> sin(δ), cos(δ)
      sine   = matmul(U, matmul(sine,   adjoint(U)))
      cosine = matmul(U, matmul(cosine, adjoint(U)))
      call packmat(sine,   sin_flat(:, ie))
      call packmat(cosine, cos_flat(:, ie))

      if(is_symmetric(sine) .AND. is_symmetric(cosine)) cycle

      ! -- error
      write(stderr, '("The Sine/Cosine matrices are not symmetric for energy ", I0, ": ", E15.7, " eV")') &
        kmat_eval_energies(ie)*au2ev
      write(stderr, '("Maxval( |sin-transpose(sine)| ): ", E15.7)') maxval(abs(sine-transpose(sine)))
      write(stderr, '("Maxval( |cosin-transpose(cosine)| ): ", E15.7)') maxval(abs(cosine-transpose(cosine)))
      call die("Non-symmetric electronic sine/cosine matrix detected")

    enddo


    ! -- 4) write to file
    !!!!!!!!!!!!!!!!!!!!!
    eigenphases_dir  = G%OUTPUT_DIRECTORY // "eigenphases/"
    eigenphases_file = eigenphases_dir // spinmult_names(spinmult) // "eigenphases.dat"
    sin_file = eigenphases_dir // spinmult_names(spinmult) // "sine_elements.dat"
    cos_file = eigenphases_dir // spinmult_names(spinmult) // "cosine_elements.dat"
    open(newunit=funit_eigenphases, file=eigenphases_file)
    open(newunit=funit_sine,        file=sin_file)
    open(newunit=funit_cosine,      file=cos_file)
    ! -- channel header
    write(funit_eigenphases, '("# ", '//i2c(n)//'(I0,",",I0,",",I0,2X))') &
      ( elec_channels(i) % nelec                                          &
      , elec_channels(i) % l                                              &
      , elec_channels(i) % ml, i=1,n)
    write(funit_sine, '("# ", '//i2c(nflat)//'(I0,",",I0,",",I0," <-> "I0,",",I0,",",I0,2X))') &
      ((elec_channels(i) % nelec                                                               &
      , elec_channels(i) % l                                                                   &
      , elec_channels(i) % ml                                                                  &
      , elec_channels(j) % nelec                                                               &
      , elec_channels(j) % l                                                                   &
      , elec_channels(j) % ml                                                                  &
      , i=1, n), j=1,i)
    write(funit_cosine, '("# ", '//i2c(nflat)//'(I0,",",I0,",",I0," <-> "I0,",",I0,",",I0,2X))') &
      ((elec_channels(i) % nelec                                                               &
      , elec_channels(i) % l                                                                   &
      , elec_channels(i) % ml                                                                  &
      , elec_channels(j) % nelec                                                               &
      , elec_channels(j) % l                                                                   &
      , elec_channels(j) % ml                                                                  &
      , i=1, n), j=1,i)
    ! -- data
    do ie=1, ne
      write(funit_eigenphases, '(E15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      write(funit_sine,        '(E15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      write(funit_cosine,      '(E15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      do i=1,nflat
        write(funit_eigenphases, '(E15.7,X)', advance='no') eigenphases(i)
        write(funit_sine,        '(E15.7,X)', advance='no') sin_flat(i, ie)
        write(funit_cosine,      '(E15.7,X)', advance='no') cosin_flat(i, ie)
      enddo
      write(funit_eigenphases,*)
      write(funit_sine,*)
      write(funit_cosine,*)
    enddo
    close(funit_eigenphases)
    close(funit_sine)
    close(funit_cosine)


  end subroutine K2sincos

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine K2S_cayley(kmat_flat, smat_flat, elec_channels)
    !! electronic Kmat -> electronic Smat via Cayley transform. Also ensures that the
    !! S-matrix is in the basis of complex-valued spherical harmonics

    use rotex__types,      only: elec_channel_type
    use rotex__arrays,     only: adjoint, eye, size_check, unpackmat, packmat, is_unitary, unitary_defect
    use rotex__constants,  only: im, au2ev
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    use rotex__linalg,     only: dsyev, zgesv

    implicit none (type, external)

    real(dp),    intent(in)  :: kmat_flat(:,:)
    complex(dp), intent(out) :: smat_flat(:,:)
    type(elec_channel_type), intent(in) :: elec_channels(:)

    character(1), parameter :: jobz = "V"
    character(1), parameter :: uplo = "U"
    integer :: n, info, ne, ie
    integer, allocatable :: ipiv(:)
    real(dp),    allocatable :: I(:,:)
    real(dp),    allocatable :: kmat(:,:)
    complex(dp), allocatable :: A(:,:), smat(:,:)

    ! -- array sizes
    n = size(elec_channels,1)
    ne = size(kmat_flat, 2)
    call size_check(Kmat_flat, [n,ne], "KMAT_FLAT")
    call size_check(Smat_flat, [n,ne], "SMAT_FLAT")

    allocate(A(n,n))
    allocate(smat(n,n))

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
      if(G%REAL_SPHERICAL_HARMONICS) call real2complex_ylm(smat, elec_channels)

      ! -- pack the Smat -> Smat_flat
      call packmat(Smat, smat_flat)

      if(is_unitary(smat)) cycle

      write(stderr, '("The S-matrix is not unitary for energy ", I0, ": ", E15.7, " eV")') kmat_eval_energies(ie)*au2ev
      write(stderr, '("Unitary defect in S: ", E15.7)') unitary_defect(Smat)
      call die("Non-unitary electronic S-matrix detected")

    enddo

  end subroutine K2S_cayley

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine do_rft(             &
        smat_elec_flat           &
      , sin_elec_flat            &
      , cos_elec_flat            &
      , Smat_j_flat              &
      , ne                       &
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
    use rotex__arrays,     only: realloc, is_unitary, uniq, nflat2n
    use rotex__functions,  only: neg
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    complex(dp),        intent(in)  :: smat_elec_flat(:)
    real(dp),        intent(in)  :: sin_elec_flat(:,:), cos_elec_flat(:,:)
    type(cmatrix_type), intent(out) :: Smat_j_flat(jmin:jmax, ne)
      !! Rotationally resolved S-matrix at each J
    integer, intent(in) :: ne
      !! Number of electronicK/S-matrix evaluation energies
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

    !@@@

    flag = .false.
    nchans_elec = size(elec_channels, 1)

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
        idx = pack([(i,i=1,nchans_j)], mask)
        rot_channels = asymtop_rot_channels_l_j(j) % channels(idx)

        ! -- allocate U, Smat_rot for this sym
        nchans_sym = size(idx, 1)
        call realloc(U,            nchans_sym, nchans_elec)
        call realloc(smat_rot_sym, nchans_sym, nchans_sym)
        U = 0
        smat_rot_sym = 0

        @@@@@@@@@@
        where to put energy loop ? out here ? in there ? omp where ? this should be parallelized over energies probably, even though
          we will only have a few hundred or whatever; just in case
        if(G%EDFT) then
          call unpack(sin_elec_flat, sin_elec)
          call unpack(cos_elec_flat, cos_elec)
        else
          call unpack(smat_elec_flat, smat_elec)
        endif

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
