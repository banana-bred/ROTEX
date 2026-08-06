! ================================================================================================================================ !
module rotex__hamilton
  !! Module containing procedures to construct and diagonalize rotational Hamiltonians
  use rotex__kinds,  only: dp
  use rotex__types,  only: eigenH_type, N_states_type
  use rotex__system, only: stderr, die

  implicit none (type, external)

  private

  ! public :: H_linear
  public :: H_asym
  public :: H_sym
  public :: assign_projections
  public :: rotate_eigvecs
  public :: wangify_symtop_eigvecs
  public :: resolve_c2prime_phi

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine H_sym(N, eigenH, Bpara, Bperp, cd4, cd6)
    !! Get the 2N+1 rotational states for a symmetric top, optionally adding
    !! diagonal centrifugal distortion (CD) terms
    !!   E(N,K) = Bperp*N*(N+1) + (Bpara-Bperp)*K² + centrifugal distortion terms

    use rotex__types,    only: cd4_type, cd6_type

    implicit none (type, external)

    integer, intent(in) :: N
      !! The rotational quantum number \(N\)
    type(eigenH_type), intent(out) :: eigenH
      !! The eigenvectors and eigenvalues of \(H\)
      !! The angular momentum number \(N\)
    real(dp), intent(in) :: Bpara
      !! The non-degenerate rotational constant (parallel to symmetry axis)
    real(dp), intent(in) :: Bperp
      !! The degenerate rotational constants (perpendicular to symmetry axis)
    type(cd4_type), intent(in), optional :: cd4
      !! The quartic centrifugal distortion parameters
    type(cd6_type), intent(in), optional :: cd6
      !! The sextic centrifugal distortion parameters

    integer :: numK, ik, K
    real(dp) :: NNp1, KK, E

    NNp1 = real(N*(N+1), kind=dp)

    numK = 2*N+1

    allocate(eigenH%eigvals(numK))
    allocate(eigenH%eigvecs(numK,numK), source=(0.0_dp, 0.0_dp))

    do K=-N,N

      iK = K+N+1

      ! -- diagonal eigvecs for symmetric top
      eigenH%eigvecs(iK, iK) = (1.0_dp, 0.0_dp)

      KK = real(K*K, kind=dp)
      E = Bperp*NNp1 + (Bpara-Bperp)*KK
      cd: if(present(cd4)) then
        ! -- 4th order CD
        E = E                            &
          - cd4%dn  * NNp1 * NNp1        &
          - cd4%dnk * NNp1 * KK          &
          - cd4%dk  * KK   * KK
        if(.not. present(cd6)) exit cd
        ! -- 6th order CD
        E = E                            &
          - cd6%hn  * NNp1 * NNp1 * NNp1 &
          - cd6%hnk * NNp1 * NNp1 * KK   &
          - cd6%hkn * NNp1 * KK   * KK   &
          - cd6%hk  * KK   * KK   * KK
      endif cd
      eigenH%eigvals(iK) = E
    enddo

  end subroutine H_sym

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine H_asym(N, eigenH, Bx, By, Bz, cd4, cd6)
    !! Construct the \(2N + 1 \times 2N + 1\) symmetric top rigid-rotor Hamiltonian
    !!   \(H = Bx N_x^2 + By N_y^2 + Bz N_z^2\)

    use rotex__utils,     only: assert
    use rotex__arrays,    only: realloc, is_symmetric
    use rotex__types,     only: cd4_type, cd6_type
    use rotex__arrays,    only: eye
    use rotex__constants, only: zero, two, four

    implicit none (type, external)

    integer, intent(in) :: N
      !! The rotational quantum number \(N\)
    type(eigenH_type), intent(out) :: eigenH
      !! The eigenvectors and eigenvalues of \(H\)
      !! The angular momentum number \(N\)
    real(dp), intent(in) :: Bx, By, Bz
      !! Rotational constants
    type(cd4_type), intent(in), optional :: cd4
      !! The quartic centrifugal distortion parameters
    type(cd6_type), intent(in), optional :: cd6
      !! The sextic centrifugal distortion parameters

    integer :: num_K
    integer :: k1, k2
    integer :: i1, i2
    integer, allocatable :: basis(:)

    real(dp) :: k
    real(dp), allocatable :: H(:,:)

    num_K = 2*N + 1

    allocate(H(num_K, num_K))

    basis = [(k1, k1=-N, N)]

    ! -- check allocation on eigen energies and vectors, avoid reallocation
    call realloc(eigenH%eigvals, num_K)
    call realloc(eigenH%eigvecs, num_K, num_K)

    ! -- fill the H, loop over projections k, k + 2, k - 2. See the matrix elements of, e.g.,
    !    "Molecular Symmetry and Spectroscopy", 2nd edition, by P.R. Bunker and P. Jensen, 11.2.4
    do concurrent(k2 = -N:N, k1 = -N:N)
      k  = k2
      i1 = (k1 + N) + 1
      i2 = (k2 + N) + 1
      if(k1 .eq. k2) then
        H(i1, i2) = (Bx+By)/2 * N*(N+1) + (Bz - (Bx+By)/2)*k*k ! 11-55: 11-56, 11-57
      elseif(k1 .eq. k2 + 2) then
        H(i1, i2) = (Bx-By)/4 * sqrt( ( N*(N+1) - (k+1)*(k+2) ) * ( N*(N+1) - k*(k+1) ) ) ! 11-55: 11-59
      elseif(k1 .eq. k2 - 2) then
        H(i1, i2) = (Bx-By)/4 * sqrt( ( N*(N+1) - (k-1)*(k-2) ) * ( N*(N+1) - k*(k-1) ) ) ! 11-55: 11-58
      else
        H(i1, i2) = zero
      endif
    enddo

    ! -- centrifugal distortion if supplied to the subroutine
    if(present(cd4)) call add_cd4(H, cd4, basis)
    if(present(cd6)) call add_cd6(H, cd6, basis)

    call assert(is_symmetric(H), "The rotational Hamiltonian is not symmetric")

    ! -- diagonalize H
    diag: block
      use rotex__linalg,     only: dsyev
      use rotex__characters, only: int2char
      integer :: lda
      integer :: lwork
      integer :: info
      character :: jobz
      character :: uplo
      real(dp), allocatable :: work(:)
      jobz    = "V" ! -- return eigenvectors
      uplo    = "U" ! -- use upper triangle "ads"
      lda     = max(1, num_K)
      lwork   = max(1, 3*num_K - 1)
      allocate(work(lwork))
      call dsyev(jobz, uplo, num_K, H, lda, eigenH % eigvals, work, lwork, info)
      if(info .ne. 0) call die("Procedure DSYEV returned with INFO = " // int2char(info))
      eigenH%eigvecs = cmplx(H, 0.0_dp, kind=dp)
    end block diag

  end subroutine H_asym

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure subroutine add_cd4(H, cd4, basis)
    !! Add the quartic centrifugal distortion effects to the Hamiltonian
    use rotex__types,     only: cd4_type
    use rotex__utils,     only: assert
    use rotex__arrays,    only: size_check
    use rotex__functions, only: isodd
    implicit none (type, external)
    real(dp),       intent(inout) :: H(:,:)
      !! Hamiltonian (already filled with A, B, C parameters)
    type(cd4_type), intent(in)    :: cd4
      !! ΔN, ΔNK, ΔK, δN, δK
    integer,        intent(in)    :: basis(:)
      !! Array of projections K of N that define the basis
    integer :: i1, i2, k1, k2
    integer :: n, nn, numk, k, kk, kp, km
    real(dp) :: dn, dnk, dk, deltan, deltak, s
    numk = size(basis, 1)
    call size_check(h, [numk, numk], "H in ADD_CD4")
    call assert(isodd(numk), "2N+1 must be odd")
    ! -- N and the eigenvalues of N² → N(N+1)
    n  = (numk-1)/2
    nn = n*(n+1)
    dn     = cd4%dn
    dnk    = cd4%dnk
    dk     = cd4%dk
    deltan = cd4%deltan
    deltak = cd4%deltak
    ! n2p1 = real(n*(n+1), kind=dp)
    do concurrent(i2=1:numk, i1=1:numk)
      k1 = basis(i1)
      k2 = basis(i2)
      k  = k2
      kk = k*k
      select case(k1 - k2)
      case(0)
        ! -- <k|h|k>
        h(i1,i2) = h(i1,i2) &
                 - dn  * (nn*nn) & ! ΔN  N⁴
                 - dnk * (nn*kk) & ! ΔNK N²Nz²
                 - dk  * (kk*kk)   ! ΔK  Nz⁴
      case(2)
        ! -- <k-2|h|k>
        s = sqrt(real( (nn-(k+1)*(k+2)) * (nn-k*(k+1)), kind = dp )) ! <N,K+2|(N⁺)²|N,K>
        kp = k + 2 ! k1
        h(i1,i2) = h(i1,i2) &
                 -          deltan * nn*s &          ! ½ δN < k+2 | [N²,  (N⁺)²]₊ |k > → δN N(N+1) <k+2|(N⁺)²|k>
                 - 0.5_dp * deltak * (kk + kp*kp)*s  ! ½ δK < k+2 | [Nz², (N⁺)²]₊ |k > → ½ δK [ k² + (k+2)² ] <k+2|(N⁺)²|k>
      case(-2)
        ! -- <k+2|h|k>
        s = sqrt(real( (nn-(k-1)*(k-2)) * (nn-k*(k-1)), kind = dp )) ! <N,K-2|(N⁺)²|N,K>
        km = k - 2 ! k1
        h(i1,i2) = h(i1,i2) &
                 -          deltan * nn*s &           ! ½ δN [N²,  (N⁻)²]₊ → δN N(N+1) <k-2|(N⁻)²|k>
                 - 0.5_dp * deltak * (kk + km*km) * s ! ½ δK [Nz², (N⁻)²]₊ → ½ δK [ k² + (k-2)² ] <k-2|(N⁻)²|k>
      end select
    enddo
  end subroutine add_cd4

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine add_cd6(H, cd6, basis)
    !! Add the quartic centrifugal distortion effects to the Hamiltonian
    use rotex__types,     only: cd6_type
    use rotex__utils,     only: assert
    use rotex__arrays,    only: size_check
    use rotex__functions, only: isodd
    implicit none (type, external)
    real(dp),       intent(inout) :: H(:,:)
      !! Hamiltonian (already filled with A, B, C parameters)
    type(cd6_type), intent(in)    :: cd6
      !! HN, HNK, HKN, HK, ηN, ηNK, ηK
    integer,        intent(in)    :: basis(:)
      !! Array of projections K of N that define the basis
    integer :: i1, i2, k1, k2
    integer :: n, nn, nn2, nn3, numk
    integer :: k, kk, kk2, kk3
    integer :: km, km2, km4, kp, kp2, kp4
    real(dp) :: hn, hnk, hkn, hk, etan, etank, etak, s
    numk = size(basis, 1)
    call size_check(h, [numk, numk], "H in ADD_CD6")
    call assert(isodd(numk), "2N+1 must be odd")
    ! -- N and the eigenvalues of N² → N(N+1)
    n  = (numk-1)/2
    nn = n*(n+1)
    nn2 = nn*nn
    nn3 = nn2*nn
    hn    = cd6%hn
    hnk   = cd6%hnk
    hkn   = cd6%hkn
    hk    = cd6%hk
    etan  = cd6%etan
    etank = cd6%etank
    etak  = cd6%etak
    ! n2p1 = real(n*(n+1), kind=dp)
    do concurrent(i2=1:numk, i1=1:numk)
      k1 = basis(i1)
      k2 = basis(i2)
      k  = k2
      kk = k*k
      kk2 = kk*kk
      kk3 = kk2*kk
      select case(k1 - k2)
      case(0)
        ! -- <k|h|k>
        h(i1,i2) = h(i1,i2) &
                 - hn  * nn3 & ! <-------- HN  (N²)³
                 - hnk * nn2 * kk  & ! <-- HNK (N²)² Nz²
                 - hkn * nn  * kk2 & ! <-- HNK  N²   Nz⁴
                 - hk  *       kk3 ! <---- HK        Nz⁶

      case(2)
        ! -- <k-2|h|k>
        s = sqrt(real( (nn-(k+1)*(k+2)) * (nn-k*(k+1)), kind = dp )) ! <N,K+2|(N⁺)²|N,K>
        kp = k + 2 ! k1
        kp2 = kp*kp
        kp4 = kp2*kp2
        h(i1,i2) = h(i1,i2) &
                 -          etan  * nn2 *             s & ! <-- ½ ηN  [(N²)²,  (N₊)²]₊
                 - 0.5_dp * etank * nn  * (kk +kp2) * s & ! <-- ½ ηNK [ N²Nz², (N₊)²]₊
                 - 0.5_dp * etak  *       (kk2+kp4) * s ! <---- ½ ηNK [   Nz⁴, (N₊)²]₊
      case(-2)
        ! -- <k+2|h|k>
        s = sqrt(real( (nn-(k-1)*(k-2)) * (nn-k*(k-1)), kind = dp )) ! <N,K-2|(N⁺)²|N,K>
        km = k - 2 ! k1
        km2 = km*km
        km4 = km2*km2
        h(i1,i2) = h(i1,i2) &
                 -          etan  * nn2 *             s & ! <-- ½ ηN  [(N²)²,  (N₋)²]₊
                 - 0.5_dp * etank * nn  * (kk +km2) * s & ! <-- ½ ηNK [ N²Nz², (N₋)²]₊
                 - 0.5_dp * etak  *       (kk2+km4) * s ! <---- ½ ηNK [   Nz⁴, (N₋)²]₊
      end select
    enddo
  end subroutine add_cd6

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine assign_projections(N, eigvecs, absKvals, sort_eigvecs)
    use rotex__arrays, only: realloc
    !! Using the eigenvectors and energies from a diagonalized rotational Hamiltonian,
    !! determine which projection is maximal. The eigenvectors can be in the Ka or Kc basis.
    !! This routine will return the array Kvals which indicats the absolute value of the projection that
    !! contributes the most to a particular eigenvector

    implicit none (type, external)

    integer, intent(in) :: N
      !! The rotational quantum number \(N\)
    integer, intent(out), allocatable :: absKvals(:)
      !! Array of the absolte value of |K| that contributes the most to a particular eigenvector
    complex(dp), intent(in) :: eigvecs(:,:)
      !! Eigenvectors
    logical, intent(in), optional :: sort_eigvecs
      !! Sort the eigenvectors ?

    logical :: sort_eigvecs_local
    integer ::  i
    integer :: num_K
    integer :: K

    integer, allocatable :: Kvals(:)

    num_K = 2*N + 1
    Kvals = [(K, K = -N, N)]
    sort_eigvecs_local = .true. ; if(present(sort_eigvecs)) sort_eigvecs_local = sort_eigvecs

    call realloc(absKvals, num_k)

    ! -- determine the projections
    do concurrent (i=1:num_K)
      absKvals(i) = abs( Kvals( maxloc(abs(eigvecs(:,i))**2, 1) ) )
    enddo

  end subroutine assign_projections

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine rotate_eigvecs(N, from_axis, to_axis, eigvecs)
    !! Rotate the rigid rotor eigenvectors from one of the principal axes A,B,C to another
    !! principal axis A,B,C using the Wigner D-matrix, while ensuring that the coordinate
    !! system remains right-handed and that each of A,B,C get one of x,y,z.
    !! The three coordinate systems are:
    !!   ABC = zxy
    !!   ABC = xyz
    !!   ABC = yzx
    !! where A,B,C is defined to be right-handed as well.
    use rotex__arrays,     only: is_unitary, unitary_defect
    use rotex__kinds,      only: dp
    use rotex__arrays,     only: adjoint
    use rotex__characters, only: lower
    use rotex__constants,  only: pi, im
    use rotex__wigner,     only: wigner_big_D, wigner_little_d
    implicit none (type, external)
    integer,      intent(in)    :: N
      !! The rotational angular moment quantum number
    character(1), intent(in) :: from_axis
      !! The starting z-axis
    character(1), intent(in)    :: to_axis
      !! The target z-axis to which we rotate
    complex(dp),  intent(inout) :: eigvecs(:,:)
      !! The eigenvectors

    real(dp) :: a, b, g
      !! Euler angles α β γ
    real(dp) :: R(3,3)
    complex(dp), allocatable :: D(:,:)

    if(lower(from_axis) .eq. lower(to_axis)) return ! no rotation needed

    ! -- determine the Euler angles α,β,γ for D=Rz(α)Ry(β)Rz(γ). The following are INTRINSIC
    !    rotations, i.e. body frame rotations where A,B,C stay fixed
    select case(from_axis)
    case("a","A")
      select case(to_axis)
      case("b", "B") ; a = pi/2 ; b = pi/2 ; g = 0
      case("c", "C") ; a = pi   ; b = pi/2 ; g = pi/2
      case default
        call die("Unacceptable TO_AXIS (" // to_axis //") provided")
      end select
    case("b","B")
      select case(to_axis)
      case("c", "C") ; a = pi/2 ; b = pi/2 ; g = 0
      case("a", "A") ; a = pi   ; b = pi/2 ; g = pi/2
      case default
        call die("Unacceptable TO_AXIS (" // to_axis //") provided")
      end select
    case("c", "C")
      select case(to_axis)
      case("a", "A") ; a = pi/2 ; b = pi/2 ; g = 0
      case("b", "B") ; a = pi   ; b = pi/2 ; g = pi/2
      case default
        call die("Unacceptable TO_AXIS (" // to_axis //") provided")
      end select
    case default
      call die("FROM_AXIS (" // ") can only be 'A'")
    end select

    ! -- swap α and γ for intrinsic -> extrinsic rotation
    D = wigner_big_D(N, g, b, a, use_analytic = .false.)
    ! D = wigner_little_d(N, b, use_analytic = .false.)

    ! -- unitarity check on D
    if(is_unitary(D) .eqv. .false.) then
      write(stderr, '("Unitary defect in D: ", F7.5)') unitary_defect(D)
      call die("Wigner D-matrix is not unitary !")
    endif

    eigvecs = matmul(adjoint(D), eigvecs)

    if(is_unitary(eigvecs)) return

    write(stderr, '("From axis: ", A)') from_axis
    write(stderr, '("To axis: ", A)') to_axis
    write(stderr, '("Unitary defect in rotated eigenvectors: ", F7.5)') unitary_defect(eigvecs)
    call die("Rotated eigenvectors are not unitary !")

  end subroutine rotate_eigvecs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function resolve_c2prime_phi(use_kmat, phi_override_deg) result(phi)
    !! Get the azimuthal angle φ of the C2' axis in the frame the rotor eigenvectors live in
    !! (after any RR_DIAG_AXIS -> SYMAXIS stuff). This is here to fix Wang phases (exp[-2iK*phi])
    !! and Wigner D-matrix rotations used to calculate rchar.
    !! If the S-matrix afterwards had forbidden elements that are allowed, this φ is probably wrong.
    !! The C2' axis must be the C2 axis that, e.g., the C2v scattering frame contains.
    !! For equilibrium H₃⁺ this is φ=π/2
    use rotex__constants, only: pi
    use rotex__globals,   only: G, is_unset
    implicit none(type, external)
    logical,  intent(in)           :: use_kmat
    real(dp), intent(in), optional :: phi_override_deg !! If present and set, override the inferred value
    real(dp) :: phi
    over: if(present(phi_override_deg)) then
      if(is_unset(phi_override_deg)) exit over
      phi = phi_override_deg * pi / 180.0_dp
      return
    endif over
    if(use_kmat) then
      phi = 0.5_dp * pi
    else
      phi = 0.0_dp ! don't do any extra rotation; not needed when using only multipoles
    endif
  end function resolve_c2prime_phi

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine wangify_symtop_eigvecs(N, eigvecs, irchar, phi)
    !! Rotate the symmetric-top rigid-rotor eigenvectors into the Wang basis, tagging each state
    !! with its C₂' character IRCHAR.
    !!
    !! H_SYM is diagonal in the signed-K basis -> its eigenvectors are the identity.
    !! Within each exactly degenerate ±K pair, any unitary mixture is an equally valid eigenbasis.
    !! The Wang combinations are the ones that ALSO diagonalize C₂', the operation carrying the
    !! information on nuclear permutation symmetry.
    !!
    !! Row ordering of the primitive |N,K⟩ basis (K = -N..N) is unchanged, so the RFT and CB
    !! routines work as expected, given that they were designed with asymtops in mind and, from that
    !! perspective, the only thing that really changes are the eigenvector expansion coefficients.

    use rotex__arrays,    only: eye, is_unitary, unitary_defect
    use rotex__constants, only: pi
    use wignerd,          only: wigner_big_D

    implicit none(type, external)

    integer,      intent(in)               :: N            !! Rotational quantum number
    complex(dp),  intent(inout)            :: eigvecs(:,:) !! Eigenvectors for this N
    integer,      intent(out), allocatable :: irchar(:)    !! Rotational character ±1
    real(dp),     intent(in)               :: phi          !! Azimuthal angle φ of the C₂' axis in the SYMAXIS frame

    real(dp), parameter :: MAX_IDENTITY_DEFECT = 1e-10_dp
    real(dp), parameter :: invsq2 = 1.0_dp / sqrt(2.0_dp)

    integer :: K, ip, im, i, numk
    real(dp) :: identity_defect, r
    complex(dp) :: ph
    complex(dp), allocatable :: W(:,:) !! Wang matrix
    complex(dp), allocatable :: RC2(:,:) !! C₂' rotation matrix
    complex(dp), allocatable :: tmp(:,:)

    numk = 2*N + 1

    ! -- make sure the eigenvectors are unit vectors
    identity_defect = maxval(abs(eigvecs - cmplx(eye(numk), kind=dp)))
    if(identity_defect .gt. MAX_IDENTITY_DEFECT) then
      write(stderr, '("identity_defect: ", ES20.12)') identity_defect
      call die("Attempt to Wangify non-unit eigenvectors")
    endif

    allocate(W(numk, numk), source=(0._dp,0._dp))
    ! -- K=0
    W(N+1, N+1) = (1.0_dp, 0.0_dp)
    ! -- K≠0
    do K=1, N
      ip = K + N + 1
      im =-K + N + 1
      ph = exp(cmplx(0, -2*K*phi, kind=dp))
      w(ip, ip) = invsq2
      w(im, ip) = invsq2 * ph
      w(ip, im) = invsq2
      w(im, im) =-invsq2 * ph
    enddo

    eigvecs = w

    if(is_unitary(eigvecs) .eqv. .false.) then
      write(stderr, '("Unitary defect in Wang eigenvectors: ", ES20.12)') unitary_defect(eigvecs)
      call die("Wang eigenvectors are not unitary !")
    endif

    ! -- C₂'(φ) = Rz(φ) Ry(π) Rz(-φ)
    RC2 = wigner_big_D(N, phi, pi, -phi, use_analytic=.false.)

    ! -- irchar(i) = <W_i|Rc2|W_i>
    allocate(irchar(numk), source=0)
    tmp = matmul(Rc2, eigvecs)
    do i=1, numk
      r = real(dot_product(eigvecs(:, i), tmp(:, i)), kind=dp)
      if(abs(abs(r) - 1.0_dp) .gt. MAX_IDENTITY_DEFECT) then
        write(stderr, '("N = ", I0, ", column ", I0, ": <w|C2p|w> = ", ES12.5)') N, i, r
        call die("Wang state is not a C2' eigenstate. Check G%C2PRIME_AZIMUTH or D-matrix conventions.")
      endif
      irchar(i) = nint(sign(1.0_dp, r))
    enddo

  end subroutine wangify_symtop_eigvecs

! ================================================================================================================================ !
end module rotex__hamilton
! ================================================================================================================================ !
