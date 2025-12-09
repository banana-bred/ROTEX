! ================================================================================================================================ !
module rotex__hamilton
  !! Module containing procedures to construct and diagonalize rotational Hamiltonians
  use rotex__types, only: dp, eigenH_type, N_states_type

  implicit none

  private

  ! public :: H_linear
  public :: H_asym
  public :: assign_projections
  ! public :: get_different_K_projections

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine H_asym(N, eigenH, Bx, By, Bz, cd4, cd6)
    !! Construct the \(2N + 1 \times 2N + 1\) symmetric top rigid-rotor Hamiltonian
    !!   \(H = Bx N_x^2 + By N_y^2 + Bz N_z^2\)

    use rotex__utils,     only: assert
    use rotex__arrays,    only: realloc, is_symmetric
    use rotex__types,     only: cd4_type, cd6_type
    use rotex__system,    only: die, stderr
    use rotex__arrays,    only: eye
    use rotex__constants, only: zero, two, four

    implicit none

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
      eigenH % eigvecs = H
      allocate(work(lwork))
      call dsyev(jobz, uplo, num_K, eigenH % eigvecs, lda, eigenH % eigvals, work, lwork, info)
      if(info .ne. 0) call die("Procedure DSYEV returned with INFO = " // int2char(info))
    end block diag

  end subroutine H_asym

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure subroutine add_cd4(H, cd4, basis)
    !! Add the quartic centrifugal distortion effects to the Hamiltonian
    use rotex__types,     only: cd4_type
    use rotex__utils,     only: assert
    use rotex__arrays,    only: size_check
    use rotex__functions, only: isodd
    implicit none
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
    implicit none
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
  module subroutine assign_projections(N, eigenH, absKvals, sort_eigvecs)
    use rotex__arrays, only: realloc
    !! Using the eigenvectors and energies from a diagonalized rotational Hamiltonian,
    !! determine which projection is maximal. The eigenvectors can be in the Ka or Kc basis.
    !! This routine will return the array Kvals which indicats the absolute value of the projection that
    !! contributes the most to a particular eigenvector

    implicit none

    integer, intent(in) :: N
      !! The rotational quantum number \(N\)
    integer, intent(out), allocatable :: absKvals(:)
      !! Array of the absolte value of |K| that contributes the most to a particular eigenvector
    type(eigenH_type), intent(in) :: eigenH
      !! Contains the eigenvalues and eigenvectors of the rotational Hamiltonian for the given value of N, digaonalized in
      !! the Ka or the Kc symmetric top basis
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
      absKvals(i) = abs( Kvals( maxloc(abs(eigenH % eigvecs(:,i))**2, 1) ) )
    enddo

  end subroutine assign_projections

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine sort_eigenvectors(eigenH)
  !   !! Sort the eigenvectors and energies so that the energies are in increasing order
  !   use rotex__arrays, only: sort_index
  !   implicit none
  !   type(eigenH_type), intent(inout) :: eigenH
  !   integer :: n
  !   integer, allocatable :: indices(:)
  !   n = size(eigenH % eigvals, 1)
  !   allocate(indices(n))
  !   ! -- sort energies, lowest first
  !   call sort_index(eigenH % eigvals, indices)
  !   ! -- swap columns based on sorted energies
  !   eigenH % eigvals = eigenH % eigvals(indices)
  !   eigenH % eigvecs(:,:) = eigenH % eigvecs(:, indices)
  ! end subroutine sort_eigenvectors

! ================================================================================================================================ !
end module rotex__hamilton
! ================================================================================================================================ !
