! ================================================================================================================================ !
module rotex__arrays
  !! Various routines for arrays
  use rotex__types, only: dp

  implicit none

  private

  public :: append
  public :: append_uniq
  public :: uniq
  public :: remove_value
  public :: eye
  public :: adjoint
  public :: size_check
  public :: is_unitary
  public :: is_symmetric
  public :: norm_frob
  public :: realloc
  public :: unitary_defect
  public :: sort_index
  ! public :: vec2diag

  interface append_uniq
    module procedure :: append_uniq_i
    module procedure :: append_uniq_transition
  end interface append_uniq
  interface append
    module procedure :: append_i
    module procedure :: append_r
    module procedure :: append_rvector
    module procedure :: append_asymtop_transition
    module procedure :: append_elec_channel
    module procedure :: append_elec_channels
    module procedure :: append_asymtop_rot_channel
    module procedure :: append_asymtop_rot_channel_l
  end interface append

  interface adjoint
    module procedure :: adjoint_i
    module procedure :: adjoint_r
    module procedure :: adjoint_c
  end interface adjoint

  interface size_check
    module procedure :: size_check_1d
    module procedure :: size_check_2d
  end interface size_check

  interface norm_frob
    module procedure :: norm_frob_i
    module procedure :: norm_frob_r
    module procedure :: norm_frob_c
  end interface norm_frob

  interface realloc
    module procedure :: realloc_1d_int
    module procedure :: realloc_1d_real
    module procedure :: realloc_1d_cmplx
    module procedure :: realloc_1d_elec_channel
    module procedure :: realloc_2d_real
    module procedure :: realloc_2d_cmplx
  end interface realloc

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_1d_int(arr, n)
    implicit none
    integer, intent(inout), allocatable :: arr(:)
    integer,  intent(in)                 :: n
    if(allocated(arr)) then
      if(size(arr, 1) .eq. n) return
      deallocate(arr)
      allocate(arr(n))
      return
    endif
    allocate(arr(n))
  end subroutine realloc_1d_int
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_1d_real(arr, n)
    implicit none
    real(dp), intent(inout), allocatable :: arr(:)
    integer,  intent(in)                 :: n
    if(allocated(arr)) then
      if(size(arr, 1) .eq. n) return
      deallocate(arr)
      allocate(arr(n))
      return
    endif
    allocate(arr(n))
  end subroutine realloc_1d_real
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_1d_cmplx(arr, n)
    implicit none
    complex(dp), intent(inout), allocatable :: arr(:)
    integer,  intent(in)                 :: n
    if(allocated(arr)) then
      if(size(arr, 1) .eq. n) return
      deallocate(arr)
      allocate(arr(n))
      return
    endif
    allocate(arr(n))
  end subroutine realloc_1d_cmplx
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_1d_elec_channel(arr, n)
    use rotex__types, only: elec_channel_type
    implicit none
    type(elec_channel_type), intent(inout), allocatable :: arr(:)
    integer,  intent(in)                 :: n
    if(allocated(arr)) then
      if(size(arr, 1) .eq. n) return
      deallocate(arr)
      allocate(arr(n))
      return
    endif
    allocate(arr(n))
  end subroutine realloc_1d_elec_channel
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_2d_real(arr, n, m)
    implicit none
    real(dp), intent(inout), allocatable :: arr(:,:)
    integer,  intent(in)                 :: n, m
    if(allocated(arr)) then
      if(all(shape(arr) .eq. [n,m])) return
      deallocate(arr)
      allocate(arr(n,m))
      return
    endif
    allocate(arr(n,m))
  end subroutine realloc_2d_real
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine realloc_2d_cmplx(arr, n, m)
    implicit none
    complex(dp), intent(inout), allocatable :: arr(:,:)
    integer,  intent(in)                 :: n, m
    if(allocated(arr)) then
      if(all(shape(arr) .eq. [n,m])) return
      deallocate(arr)
      allocate(arr(n,m))
      return
    endif
    allocate(arr(n,m))
  end subroutine realloc_2d_cmplx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function norm_frob_i(A) result(res)
    !! Returns the Frobenius norm for a matrix A
    implicit none
    integer, intent(in) :: A(:,:)
    real(dp) :: res
    res = sqrt(real(sum(A*A), kind=dp))
  end function norm_frob_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function norm_frob_r(A) result(res)
    !! Returns the Frobenius norm for a matrix A
    implicit none
    real(dp), intent(in) :: A(:,:)
    real(dp) :: res
    res = sqrt(sum(A*A))
  end function norm_frob_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function norm_frob_c(A) result(res)
    !! Returns the Frobenius norm for a matrix A
    implicit none
    complex(dp), intent(in) :: A(:,:)
    real(dp) :: res
    res = sqrt(sum(abs(A)**2))
  end function norm_frob_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function adjoint_i(A) result(res)
    !! Returns the adjoint of an integer-valued matrix
    implicit none
    integer, intent(in) :: A(:,:)
    integer :: res(size(A, 2), size(A, 1))
    res = transpose(A)
  end function adjoint_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function adjoint_r(A) result(res)
    !! Returns the adjoint of a real-valued matrix
    implicit none
    real(dp), intent(in) :: A(:,:)
    real(dp) :: res(size(A, 2), size(A, 1))
    res = transpose(A)
  end function adjoint_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function adjoint_c(A) result(res)
    !! Returns the adjoint of a complex-valued matrix
    implicit none
    complex(dp), intent(in) :: A(:,:)
    complex(dp) :: res(size(A, 2), size(A, 1))
    res = conjg(transpose(A))
  end function adjoint_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function unitary_defect(A) result(rF)
    !! Return the unitary defect with respect to the Frobenius norm
    !! \(rF = ||A^{\dagger}A-I||_F / sqrt{n}\)
    use rotex__system, only: die
    implicit none
    class(*), intent(in) :: A(:,:)
    real(dp) :: rF
    integer :: n, m
    integer :: i, j
    real(dp) :: err2

    n = size(A, 1)
    m = size(A, 2)

    if(n .ne. m) call die("Cannot determine the unitarity of a nonsquare matrix")

    err2 = 0

    select type(A)

    type is (real(dp))

    realmat: block
      real(dp), allocatable :: G(:,:)
      real(dp) :: d
      allocate(G(n,m))
      ! -- Gram matrix
      G = matmul(adjoint(A), A)
      do j=1,n
        do i=1,n
          ! -- subtract identity if i==j
          d = G(i,j) - merge(1, 0, i .eq. j)
          ! -- accumulate norm squared
          err2 = err2 + d*d
        enddo
      enddo
    end block realmat

    type is (complex(dp))

    cmplxmat: block
      complex(dp), allocatable :: G(:,:)
      complex(dp) :: d
      allocate(G(n,m))
      ! -- Gram matrix
      G = matmul(adjoint(A), A)
      do j=1,n
        do i=1,n
          ! -- subtract identity if i==j
          d = G(i,j) - merge(1, 0, i .eq. j)
          ! -- accumulate norm squared
          err2 = err2 + abs(d)**2
        enddo
      enddo
    end block cmplxmat

    class default

      call die("Trying to determine the unitary of a matrix that is neither real nor complex")

    end select

    ! -- relative size-aware defect
    rF = sqrt(err2) / sqrt(real(n, kind = dp))

  end function unitary_defect

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_unitary(A, rtol) result(res)
    !! Determine whether a matrix is unitary w.r.t the Frobenius norm
    use rotex__constants, only: macheps => macheps_dp
    class(*), intent(in) :: A(:,:)
      !! The matrix
    real(dp), intent(in), optional :: rtol
      !! The optional relative tolerance, default 1e-10
    logical :: res
    integer :: n, m
    real(dp) :: rF, tol
    tol = 1e-10_dp ; if(present(rtol)) tol = rtol
    res = .false.
    n = size(A, 1)
    m = size(A, 2)
    if(n .ne. m) return
    rF = unitary_defect(A)
    tol = max(tol, 100*macheps*sqrt(real(n, kind=dp)))
    res = rF .lt. tol
  end function is_unitary

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_symmetric(A, rtol) result(res)
    !! Determine whether a matrix is symmetric
    use rotex__constants, only: macheps => macheps_dp
    class(*), intent(in) :: A(:,:)
      !! The matrix
    real(dp), intent(in), optional :: rtol
      !! The optional relative tolerance, default 1e-10
    logical :: res
    integer :: n, m
    real(dp) :: tol, diff, denom
    tol = 1e-10_dp ; if(present(rtol)) tol = rtol
    res = .false.
    n = size(A, 1)
    m = size(A, 2)
    if(n .ne. m) return
    select type(AA => A)
    type is (integer)
      diff = norm_frob(AA - transpose(AA))
      denom = max(1._dp, norm_frob(AA))
    type is (real(dp))
      diff = norm_frob(AA - transpose(AA))
      denom = max(1._dp, norm_frob(AA))
    type is (complex(dp))
      diff = norm_frob(AA - transpose(AA))
      denom = max(1._dp, norm_frob(AA))
    end select
    res = diff .le. tol*denom + sqrt(macheps)
  end function is_symmetric

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function eye(n) result(res)
    !! Return an n x n identity matrix
    use rotex__system, only: die
    integer, intent(in) :: n
    integer :: res(n,n)
    integer :: i
    if(n .lt. 1) call die("Attempt to make an identity matrix with dims < 1")
    res = 0
    do concurrent(i=1:n)
      res(i,i) = 1
    enddo
  end function eye

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine remove_value(arr, val)
    !! Remove all instances of the value val from the array arr
    implicit none
    integer, intent(inout), allocatable :: arr(:)
    integer, intent(in) :: val
    integer, allocatable :: tmp(:)
    integer :: i, n
    if(.not. allocated(arr)) return
    n = size(arr, 1)
    do i=1, n
      if(val .ne. arr(i)) call append(tmp, arr(i))
    enddo
    if(.not. allocated(tmp)) then
      deallocate(arr)
      return
    endif
    call move_alloc(tmp, arr)
  end subroutine remove_value

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine append_uniq_i(arr, new)
    !! Append unique element "new" to array "arr"
    implicit none
    integer, intent(in)                 :: new
    integer, intent(inout), allocatable :: arr(:)
    select case(allocated(arr))
    case(.true.)
      if(any(arr .eq. new)) return
      arr = [arr, new]
    case(.false.)
      arr = [new]
    end select
  end subroutine append_uniq_i

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine append_uniq_transition(old, new)
    !! Append unique element "new" to array "old"
    use rotex__types, only: asymtop_rot_transition_type, operator(.ne.), operator(.isin.)
    implicit none
    type(asymtop_rot_transition_type), intent(in)                 :: new(:)
    type(asymtop_rot_transition_type), intent(inout), allocatable :: old(:)

    logical, allocatable :: mask(:)
    integer :: itrans, nold, num_new_uniq, nnew
    type(asymtop_rot_transition_type), allocatable :: tmp(:)

    select case(allocated(old))

    case(.false.)

      old = new

    case(.true.)

      ! -- count number of uniq elements
      nnew = size(new, 1)
      nold = size(old, 1)
      allocate(mask(nnew), source=.false.)
      if((new(1) .isin. old(:)) .eqv. .false.) mask(1) = .true. ! check first element
      ! get uniq elemnt indices
      do itrans = 2, nnew
        if(new(itrans) .isin. new(:itrans-1)) cycle ! skip dupicate elements in new
        if(new(itrans) .isin. old(:)) cycle
        mask(itrans) = .true.
      enddo
      num_new_uniq = count(mask)

      ! -- return if nothing to append
      if(num_new_uniq .eq. 0) return

      allocate(tmp(nold + num_new_uniq))
      tmp(1:nold)  = old(:)
      tmp(nold+1:) = pack(new, mask)
      call move_alloc(tmp, old)

    end select

  end subroutine append_uniq_transition

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine append_i(arr, val)
    !! Append the value val to the array arr
    implicit none
    integer, intent(inout), allocatable :: arr(:)
    integer, intent(in) :: val
    select case(allocated(arr))
    case(.true.)
      arr = [arr, val]
    case(.false.)
      arr = [val]
    end select
  end subroutine append_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine append_r(arr, val)
    !! Append the value val to the array arr
    implicit none
    real(dp), intent(inout), allocatable :: arr(:)
    real(dp), intent(in) :: val
    select case(allocated(arr))
    case(.true.)
      arr = [arr, val]
    case(.false.)
      arr = [val]
    end select
  end subroutine append_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine append_rvector(arr, val)
    !! Append the value val to the array arr
    use rotex__types, only: rvector_type
    implicit none
    type(rvector_type), intent(inout), allocatable :: arr(:)
    real(dp), intent(in) :: val(:)
    type(rvector_type) :: elem
    allocate(elem%vec(size(val, 1)))
    elem % vec = val
    select case(allocated(arr))
    case(.true.)
      arr = [arr, elem]
    case(.false.)
      arr = [elem]
    end select
  end subroutine append_rvector
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine append_asymtop_transition(arr, val)
    !! Append the value val to the array arr
    use rotex__types, only: asymtop_rot_transition_type
    implicit none
    type(asymtop_rot_transition_type), intent(inout), allocatable :: arr(:)
    type(asymtop_rot_transition_type), intent(in) :: val
    select case(allocated(arr))
    case(.true.)
      arr = [arr, val]
    case(.false.)
      arr = [val]
    end select
  end subroutine append_asymtop_transition
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module subroutine append_elec_channel(channels, channel)
    !! Append elec_channel to the array elec_channels
    use rotex__types, only: elec_channel_type
    implicit none
    type(elec_channel_type), intent(inout), allocatable :: channels(:)
    type(elec_channel_type), intent(in) :: channel
    select case(allocated(channels))
    case(.true.)  ; channels = [channels, channel]
    case(.false.) ; channels = [channel]
    end select
  end subroutine append_elec_channel
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module subroutine append_elec_channels(channels, channels2)
    !! Append channels2 to the array channels
    use rotex__types, only: elec_channel_type
    implicit none
    type(elec_channel_type), intent(inout), allocatable :: channels(:)
    type(elec_channel_type), intent(in) :: channels2(:)
    select case(allocated(channels))
    case(.true.)  ; channels = [channels, channels2]
    case(.false.) ; channels = [channels2]
    end select
  end subroutine append_elec_channels
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module subroutine append_asymtop_rot_channel(channels, channel)
    !! Append asymtop_rot_channel to the array asymtop_rot_channels
    use rotex__types, only: asymtop_rot_channel_type
    implicit none
    type(asymtop_rot_channel_type), intent(inout), allocatable :: channels(:)
    type(asymtop_rot_channel_type), intent(in) :: channel
    select case(allocated(channels))
    case(.true.)  ; channels = [channels, channel]
    case(.false.) ; channels = [channel]
    end select
  end subroutine append_asymtop_rot_channel
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module subroutine append_asymtop_rot_channel_l(channels, channel)
    !! Append asymtop_rot_channel to the array asymtop_rot_channels
    use rotex__types, only: asymtop_rot_channel_l_type
    implicit none
    type(asymtop_rot_channel_l_type), intent(inout), allocatable :: channels(:)
    type(asymtop_rot_channel_l_type), intent(in) :: channel
    select case(allocated(channels))
    case(.true.)  ; channels = [channels, channel]
    case(.false.) ; channels = [channel]
    end select
  end subroutine append_asymtop_rot_channel_l

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine size_check_1d(arr, larr, name)
    !! Check that the size of the array arr is of length larr
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    implicit none
    class(*), intent(in) :: arr(:)
    integer, intent(in) :: larr
    character(*), intent(in) :: name
    if(size(arr, 1) .ne. larr) &
      call die("Array " // name // "(:) " // i2c(shape(arr)) // " must have the shape " // i2c([larr]))
  end subroutine size_check_1d
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine size_check_2d(arr, larr, name)
    !! Check that the size of the array arr is of length larr
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    implicit none
    class(*), intent(in) :: arr(:,:)
    integer, intent(in) :: larr(:)
    character(*), intent(in) :: name
    if(any(shape(arr) .ne. larr)) &
      call die("Array " // name // "(:,:) " // i2c(shape(arr)) // " must have the shape " // i2c(larr))
  end subroutine size_check_2d

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function uniq(arr) result(res)
    !! Returns the unique elements of arr
    implicit none
    integer, intent(in)  :: arr(:)
    integer, allocatable :: res(:)
    integer :: n, i, k
    n = size(arr, 1)
    k = 0
    allocate(res(n))
    do i=1, n
      if(any(res(1:k) .eq. arr(i))) cycle
      k = k + 1
      res(k) = arr(i)
    enddo
    if(k .eq. 0) then
      call realloc(res, 0)
      return
    endif
    res = res(1:k)
  end function uniq

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine sort_index(vals, idx)
    !! Sort the array vals and return the permutation indices
    implicit none
    real(dp),    intent(in)  :: vals(:)
    integer,     intent(out) :: idx(:)
    integer :: i, j, n
    integer :: key
    real(dp) :: targ
    n = size(vals, 1)
    call size_check(idx, n, "IDX")
    do concurrent(i=1:n) ; idx(i) = i ; enddo
    if(n .le. 1) return
    do i = 2, n
      key = idx(i)
      j = i - 1
      targ = vals(key)
      do
        if(j .lt. 1) exit
        if(vals(idx(j)) .le. targ) exit
        idx(j+1) = idx(j)
        j = j - 1
      enddo
      idx(j+1) = key
    enddo
  end subroutine sort_index

! ================================================================================================================================ !
end module rotex__arrays
! ================================================================================================================================ !
