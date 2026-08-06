! ================================================================================================================================ !
module rotex__utils
  !! Some small utilities
  use rotex__kinds, only: dp, qp

  implicit none (type, external)

  private

  public :: read_blank
  public :: assert
  public :: isint
  public :: kbn_sum
  public :: upcast
  public :: downcast
  public :: printmat
  public :: isin
  public :: estimate_total_storage_size
  public :: bytes2human
  public :: halfint_float_to_rational

  interface isint
    module procedure :: isint_r
    module procedure :: isint_c
  end interface isint

  interface kbn_sum
    module procedure :: kbn_sum_rdp
    module procedure :: kbn_sum_rqp
    module procedure :: kbn_sum_cdp
    module procedure :: kbn_sum_cqp
  end interface kbn_sum

  interface downcast
    module procedure :: downcast_r
    module procedure :: downcast_c
  end interface downcast

  interface upcast
    module procedure :: upcast_r
    module procedure :: upcast_c
  end interface upcast

  interface printmat
    module procedure :: printmat_i
    module procedure :: printmat_r
    module procedure :: printmat_c
  end interface printmat

  interface bytes2human
    module procedure :: bytes2human_int32
    module procedure :: bytes2human_int64
  end interface bytes2human

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine assert(test, message)
    use rotex__system, only: die
    implicit none (type, external)
    logical,      intent(in) :: test
    character(*), intent(in) :: message
    if(test .eqv. .true.) return
    call die(message)
  end subroutine assert

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_blank(read_unit, num_read)
    !! Reads num_read lines from unit read_unit, not storing any information. If num_read is not supplied, read one line.
    implicit none (type, external)
    integer, intent(in)           :: read_unit
    integer, intent(in), optional :: num_read
    integer :: k, n
    n = 1 ; if(present(num_read)) n = num_read
    do k = 1, n ; read(read_unit,*) ; enddo
  end subroutine read_blank

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isint_r(x) result(res)
    implicit none (type, external)
    real(dp), intent(in) :: x
    logical :: res
    real(dp) :: tol
    tol = 8*spacing(x)
    res = .false.
    if(abs(x - anint(x)) .gt. tol) return
    res = .true.
  end function isint_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isint_c(z) result(res)
    implicit none (type, external)
    complex(dp), intent(in) :: z
    logical :: res
    real(dp) :: tol
    real(dp) :: a
    a = z%re
    tol = 8*spacing(a)
    res = .false.
    if(abs(z%im) .gt. tol) return
    if(abs(a - anint(a)) .gt. tol) return
    res = .true.
  end function isint_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine kbn_sum_rdp(summation, c, input)
    !! Improved Kahan-Babuška algorithm accumulation for summations
    implicit none (type, external)
    real(dp), intent(inout) :: summation, c
    real(dp), intent(in)    :: input
    real(dp) :: t
    t = summation + input
    if(abs(summation) .ge. abs(input)) then
      c = c + (summation-t) + input
    else
      c = c + (input - t) + summation
    endif
    summation = t
  end subroutine kbn_sum_rdp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine kbn_sum_rqp(summation, c, input)
    !! Improved Kahan-Babuška algorithm accumulation for summations
    implicit none (type, external)
    real(qp), intent(inout) :: summation, c
    real(qp), intent(in)    :: input
    real(qp) :: t
    t = summation + input
    if(abs(summation) .ge. abs(input)) then
      c = c + (summation-t) + input
    else
      c = c + (input - t) + summation
    endif
    summation = t
  end subroutine kbn_sum_rqp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine kbn_sum_cdp(summation, c, input)
    !! Improved Kahan-Babuška algorithm accumulation for summations
    implicit none (type, external)
    complex(dp), intent(inout) :: summation, c
    complex(dp), intent(in)    :: input
    complex(dp) :: t
    t = summation + input
    if(abs(summation) .ge. abs(input)) then
      c = c + (summation-t) + input
    else
      c = c + (input - t) + summation
    endif
    summation = t
  end subroutine kbn_sum_cdp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine kbn_sum_cqp(summation, c, input)
    !! Improved Kahan-Babuška algorithm accumulation for summations
    implicit none (type, external)
    complex(qp), intent(inout) :: summation, c
    complex(qp), intent(in)    :: input
    complex(qp) :: t
    t = summation + input
    if(abs(summation) .ge. abs(input)) then
      c = c + (summation-t) + input
    else
      c = c + (input - t) + summation
    endif
    summation = t
  end subroutine kbn_sum_cqp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module subroutine downcast_r(hi, lo)
    !! Send the value of hi to lo, respecting the kind of the types
    use rotex__kinds, only: dp, qp
    implicit none (type, external)
    real(qp), intent(in)  :: hi
    real(dp), intent(out) :: lo
    lo = real(hi, kind = dp)
  end subroutine downcast_r
  pure elemental module subroutine downcast_c(hi, lo)
    !! Send the value of hi to lo, respecting the kind of the types
    use rotex__kinds, only: dp, qp
    implicit none (type, external)
    complex(qp), intent(in)  :: hi
    complex(dp), intent(out) :: lo
    lo = cmplx(hi%re, hi%im, kind = dp)
  end subroutine downcast_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module subroutine upcast_r(lo, hi)
    !! Send the value of lo to hi, respecting the kind of the types
    use rotex__kinds, only: dp, qp
    implicit none (type, external)
    real(dp), intent(in)  :: lo
    real(qp), intent(out) :: hi
    hi = real(lo, kind = qp)
  end subroutine upcast_r
  pure elemental module subroutine upcast_c(lo, hi)
    !! Send the value of lo to hi, respecting the kind of the types
    use rotex__kinds, only: dp, qp
    implicit none (type, external)
    complex(dp), intent(in)  :: lo
    complex(qp), intent(out) :: hi
    hi = cmplx(lo%re, lo%im, kind = qp)
  end subroutine upcast_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine printmat_i(M, funit, header)
    !! Prints a matrix to the supplied funit, otherwise print to stdout
    use rotex__kinds,  only: dp
    use rotex__system, only: stdout
    integer,      intent(in)           :: M(:,:)
    integer,      intent(in), optional :: funit
    character(*), intent(in), optional :: header
    character(9), parameter :: fmt = '(X,I7)'
    integer :: nr,nc, i, j, funit_local
    funit_local = stdout ; if(present(funit)) funit_local = funit
    nr = size(M, 1)
    nc = size(M, 2)
    write(funit_local, *)
    if(present(header)) write(funit_local, '(A)') header
    do i=1, nr
      do j=1, nc
        write(funit_local, fmt, advance = "no") M(i,j)
      enddo
      write(funit_local, *)
    enddo
  end subroutine printmat_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine printmat_r(M, funit, header)
    !! Prints a matrix to the supplied funit, otherwise print to stdout
    use rotex__kinds,  only: dp
    use rotex__system, only: stdout
    real(dp),     intent(in)           :: M(:,:)
    integer,      intent(in), optional :: funit
    character(*), intent(in), optional :: header
    real(dp), parameter :: absmin = 1e-5_dp
    real(dp), parameter :: absmax = 1e2_dp
    character(9) :: exp_fmt = '(X,E11.4)'
    character(9) :: flt_fmt = '(X,F11.8)'
    integer :: nr,nc, i, j, funit_local
    real(dp) :: absm
    character(:), allocatable :: fmt
    funit_local = stdout ; if(present(funit)) funit_local = funit
    nr = size(M, 1)
    nc = size(M, 2)
    write(funit_local, *)
    if(present(header)) write(funit_local, '(A)') header
    do i=1, nr
      do j=1, nc
        absm = M(i,j)
        if( (absm .le. absmin .AND. absm .ne. 0.0_dp) .OR. absm .ge. absmax ) then
          fmt = exp_fmt
        else
          fmt = flt_fmt
        endif
        write(funit_local, fmt, advance = "no") M(i,j)
      enddo
      write(funit_local, *)
    enddo
  end subroutine printmat_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine printmat_c(M, funit, header)
    !! Prints a matrix to the supplied funit, otherwise print to stdout
    use rotex__kinds,    only: dp
    use rotex__system,   only: stdout
    use ieee_arithmetic, only: copysign => ieee_copy_sign
    complex(dp),     intent(in)           :: M(:,:)
    integer,      intent(in), optional :: funit
    character(*), intent(in), optional :: header
    real(dp), parameter :: absmin = 1e-5_dp
    real(dp), parameter :: absmax = 1e2_dp
    character(9) :: exp_fmt = '(X,E11.4)'
    character(9) :: flt_fmt = '(X,F11.8)'
    integer :: nr,nc, i, j, funit_local
    integer :: signc
    real(dp) :: absmr, absmc
    character(:), allocatable :: fmtr, fmtc
    funit_local = stdout ; if(present(funit)) funit_local = funit
    nr = size(M, 1)
    nc = size(M, 2)
    write(funit_local, *)
    if(present(header)) write(funit_local, '(A)') header
    do i=1, nr
      do j=1, nc
        absmr = abs(M(i,j)%re)
        absmc = abs(M(i,j)%im)
        if( (absmr .le. absmin .AND. absmr .ne. 0._dp) .OR. absmr .ge. absmax ) then
          fmtr = exp_fmt
        else
          fmtr = flt_fmt
        endif
        if( (absmc .le. absmin .AND. absmc .ne. 0._dp) .OR. absmc .ge. absmax ) then
          fmtc = exp_fmt
        else
          fmtc = flt_fmt
        endif
        signc = nint(copysign(1._dp, M(i,j)%re))
        write(funit_local, fmtr,  advance = "no") M(i,j)%re
        write(funit_local, '(A)', advance = "no") merge(" +", " -", M(i,j)%im .ge. 0.0_dp)
        write(funit_local, fmtc,  advance = "no") absmc
        write(funit_local, '(A)', advance = "no") " im,"
      enddo
      write(funit_local, *)
    enddo
  end subroutine printmat_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isin(x, xl, xr, lclosed, rclosed) result(res)
    !! Test whether x is in the interval spanned by x1,x2
    !! l/rclosed if true include xl and xr, respectively. They are true by default
    implicit none (type, external)
    real(dp), intent(in) :: x, xl, xr
    logical, intent(in), optional :: lclosed, rclosed
    logical :: res
    logical :: lclosed_, rclosed_
    lclosed_ = .true. ; if(present(lclosed)) lclosed_ = lclosed
    rclosed_ = .true. ; if(present(rclosed)) rclosed_ = rclosed
    res = .false.
    if(x .lt. xl) return
    if(x .gt. xr) return
    if(lclosed_ .eqv. .false.) then
      if(x .eq. xl) return
    endif
    if(rclosed_ .eqv. .false.) then
      if(x .eq. xr) return
    endif
    res = .true.
  end function isin

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine estimate_total_storage_size(obj, dims, storage, units, multby)
    !! Estimate how much space will be taken up by obj
    use rotex__kinds, only: int64, prob_rk
    use rotex__types, only: prob_vector_type
    implicit none (type, external)
    ! class(*),     intent(in)  :: obj(..)
    class(*),     intent(in)  :: obj
      !! Some object that represents the type that will take up space
    integer,      intent(in)  :: dims(:)
      !! Dimensions of the object
    real(dp),     intent(out) :: storage
      !! Number of `units` that will be taken up, approximately`
    character(2), intent(out) :: units
      !! The units of `storage. One of: B, KB, MB, GB
    integer, intent(in), optional :: multby
      !! Multiply the storage size by this amount
    integer(int64) :: nbytes
    select type(obj)
    type is (prob_vector_type)
      nbytes = int(storage_size(0.0_prob_rk)/8, kind=int64) * int(product(dims), kind=int64)
    class default
      nbytes = int(storage_size(obj)/8, kind=int64) * product(int(dims, kind=int64))
    end select
    if(present(multby)) nbytes = nbytes * multby
    call bytes2human(nbytes, storage, units)
  end subroutine estimate_total_storage_size

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module subroutine bytes2human_int32(nbytes, storage, units)
    !! Convert an integer number of bytes to something like KB, MB, GB
    use rotex__kinds, only: int32, int64
    implicit none (type, external)
    integer(int32), intent(in) :: nbytes
    real(dp),       intent(out) :: storage
    character(2),   intent(out) :: units
    call bytes2human_int64(int(nbytes, kind=int64), storage, units)
  end subroutine bytes2human_int32

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module subroutine bytes2human_int64(nbytes, storage, units)
    !! Convert an integer number of bytes to something like KB, MB, GB
    use rotex__kinds, only: int64
    implicit none (type, external)
    integer(int64), intent(in) :: nbytes
    real(dp),       intent(out) :: storage
    character(2),   intent(out) :: units
    integer(int64), parameter :: NBYTES_PER_KB = 1000_int64
    integer(int64), parameter :: NBYTES_PER_MB = 1000_int64 * NBYTES_PER_KB
    integer(int64), parameter :: NBYTES_PER_GB = 1000_int64 * NBYTES_PER_MB
    real(dp) :: nbytes_r
    nbytes_r = real(nbytes, kind=dp)
    if(nbytes .lt. NBYTES_PER_KB) then
      storage = nbytes_r
      units = " B"
    elseif(nbytes .lt. NBYTES_PER_MB) then
      storage = nbytes_r / NBYTES_PER_KB
      units = "KB"
    elseif(nbytes .lt. NBYTES_PER_GB) then
      storage = nbytes_r / NBYTES_PER_MB
      units = "MB"
    else
      storage = nbytes_r / NBYTES_PER_GB
      units = "GB"
    endif
  end subroutine bytes2human_int64

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine halfint_float_to_rational(x, numer, denom)
    implicit none(type, external)
    real(dp), intent(in) :: x
    integer, intent(out) :: numer, denom
    numer = nint(2.0_dp * x)
    denom = 2
    if(mod(numer, 2) .ne. 0) return
    numer = numer / 2
    denom = 1
  end subroutine halfint_float_to_rational

! ================================================================================================================================ !
end module rotex__utils
! ================================================================================================================================ !
