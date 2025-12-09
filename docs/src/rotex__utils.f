! ================================================================================================================================ !
module rotex__utils
  !! Some small utilities
  use rotex__kinds, only: dp, qp

  implicit none

  private

  public :: read_blank
  public :: assert
  public :: isint
  public :: kbn_sum
  public :: upcast
  public :: downcast
  public :: printmat
  public :: isin

  interface isint
    module procedure :: isint_r
    module procedure :: isint_c
  end interface isint

  interface kbn_sum
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

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine assert(test, message)
    use rotex__system, only: die
    implicit none
    logical,      intent(in) :: test
    character(*), intent(in) :: message
    if(test .eqv. .true.) return
    call die(message)
  end subroutine assert

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_blank(read_unit, num_read)
    !! Reads num_read lines from unit read_unit, not storing any information. If num_read is not supplied, read one line.
    implicit none
    integer, intent(in)           :: read_unit
    integer, intent(in), optional :: num_read
    integer :: k, n
    n = 1 ; if(present(num_read)) n = num_read
    do k = 1, n ; read(read_unit,*) ; enddo
  end subroutine read_blank

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isint_r(x) result(res)
    implicit none
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
    implicit none
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
  pure elemental subroutine kbn_sum_rqp(summation, c, input)
    !! Improved Kahan-Babuška algorithm accumulation for summations
    implicit none
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
    implicit none
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
    implicit none
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
    use rotex__types, only: dp, qp
    implicit none
    real(qp), intent(in)  :: hi
    real(dp), intent(out) :: lo
    lo = real(hi, kind = dp)
  end subroutine downcast_r
  pure elemental module subroutine downcast_c(hi, lo)
    !! Send the value of hi to lo, respecting the kind of the types
    use rotex__types, only: dp, qp
    implicit none
    complex(qp), intent(in)  :: hi
    complex(dp), intent(out) :: lo
    lo = cmplx(hi%re, hi%im, kind = dp)
  end subroutine downcast_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module subroutine upcast_r(lo, hi)
    !! Send the value of lo to hi, respecting the kind of the types
    use rotex__types, only: dp, qp
    implicit none
    real(dp), intent(in)  :: lo
    real(qp), intent(out) :: hi
    hi = real(lo, kind = qp)
  end subroutine upcast_r
  pure elemental module subroutine upcast_c(lo, hi)
    !! Send the value of lo to hi, respecting the kind of the types
    use rotex__types, only: dp, qp
    implicit none
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
        write(funit_local, '(A)', advance = "no") " +"
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
    implicit none
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

! ================================================================================================================================ !
end module rotex__utils
! ================================================================================================================================ !
