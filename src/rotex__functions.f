! ================================================================================================================================ !
module rotex__functions
  !! Contains various small functions used in the code

  use rotex__kinds,  only: dp

  private

  public :: delta
  public :: factorial
  public :: log_factorial
  public :: are_approx_eq
  public :: expm1
  public :: logp1
  public :: isnatural
  public :: isinteger
  public :: arg
  public :: inv
  public :: iseven
  public :: isodd
  public :: istriangle
  public :: logrange
  public :: logb
  public :: neg

  interface are_approx_eq
    !! Compare two a and b and see if the magnitude of their difference is smaller than a tolerance,
    !! taking machine epsilon*max(|a|,|b|) for their precision as the default value
    module procedure :: are_approx_eqz
    module procedure :: are_approx_eqr
  end interface are_approx_eq

  interface factorial
    !! !n
    module procedure :: factorial_int
    module procedure :: factorial_real
  end interface factorial

  interface log_factorial
    !! ln(!n)
    module procedure :: log_factorial_i
    module procedure :: log_factorial_r
  end interface log_factorial

  interface expm1
    module procedure :: expm1r
    module procedure :: expm1z
  end interface expm1

  interface logp1
    module procedure :: logp1r
    module procedure :: logp1z
  end interface logp1

  interface logb
    module procedure :: logb_ii
    module procedure :: logb_rr
    module procedure :: logb_cc
  end interface logb

  interface logrange
    module procedure :: logrange_ib
    module procedure :: logrange_rb
  end interface logrange

  interface
    !! Interface to the C functions expm1 for real x
    function c_expm1(x) bind(C, name="expm1") result(res)
      use, intrinsic :: iso_c_binding, only: c_double
      implicit none
      real(c_double), value :: x
      real(c_double) :: res
    end function c_expm1
  end interface

  interface
    !! Interface to the C functions logp1 for real x
    function c_logp1(x) bind(C, name="logp1") result(res)
      use, intrinsic :: iso_c_binding, only: c_double
      implicit none
      real(c_double), value :: x
      real(c_double) :: res
    end function c_logp1
  end interface

  interface isinteger
    !! Check if a real/complex number is an integer
    module procedure :: isintegerr
    module procedure :: isintegerz
  end interface isinteger

  interface isnatural
    !! Check if a real/complex number is a natrual number
    module procedure :: isnaturalr
    module procedure :: isnaturalz
  end interface isnatural

  interface arg
    !! Return the phase of a complex number in (-π,π]
    module procedure argi
    module procedure argr
    module procedure argc
  end interface arg

  interface inv
    !! Compute 1/z
    module procedure invi
    module procedure invr
    module procedure invz
  end interface inv

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elementaL module function are_approx_eqz(a, b, tol) result(res)
    !! Comparre two a and b and see if magnitude of their difference is smaller than a tolerance,
    !! taking machine epsilon for their precision as the default value
    use rotex__constants, only: macheps => macheps_dp

    implicit none

    complex(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    logical :: res

    real(dp) :: tol_local

    tol_local = macheps*max(abs(a), abs(b)) ; if(present(tol)) tol_local = tol

    res = .false.
    if(abs(a - b) .gt. tol_local) return
    res = .true.

  end function are_approx_eqz
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elementaL module function are_approx_eqr(a, b, tol) result(res)
    !! Comparre two a and b and see if magnitude of their difference is smaller than a tolerance,
    !! taking machine epsilon for their precision as the default value
    use rotex__constants, only: macheps => macheps_dp

    implicit none

    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: tol
    logical :: res

    real(dp) :: tol_local

    tol_local = macheps*max(abs(a), abs(b)) ; if(present(tol)) tol_local = tol

    res = .false.
    if(abs(a - b) .gt. tol_local) return
    res = .true.

  end function are_approx_eqr

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function delta(m, n) result(res)
    !! Return the Kronecker delta function \(\delta_{m,n}\)
    use rotex__constants, only: zero, one
    implicit none
    integer, intent(in) :: m, n
    real(dp) :: res
    res = zero ; if(m .eq. n) res = one
  end function delta

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function factorial_int(n) result(res)
    !! !n
    implicit none
    integer, intent(in) :: n
    real(dp) :: res
    res = gamma(real(n + 1, kind = dp))
  end function factorial_int

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function factorial_real(n) result(res)
    !! !n
    use rotex__constants, only: one
    implicit none
    real(dp), intent(in) :: n
    real(dp) :: res
    res = gamma(n + one)
  end function factorial_real

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function expm1r(x) result(res)
    !! Returns \(e^x - 1\). TODO, replace with a fast and accurate version
    !! coded natively in fortran !

    implicit none

    real(dp), intent(in), value :: x
    real(dp) :: res
    res = c_expm1(x)
  end function expm1r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function expm1z(z) result(res)
    !! Returns \(e^x - 1\). TODO, replace with a fast and accurate version
    !! coded natively in fortran !
    implicit none
    complex(dp), intent(in) :: z ! = a + ib
    complex(dp) :: res
    real(dp) :: sin_bover2
    sin_bover2 = sin(z%im / 2)
    res % re = c_expm1(z%re)*cos(z%im) - 2*sin_bover2*sin_bover2
    res % im = exp(z%re)*sin(z%im)
  end function expm1z

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function logp1r(x) result(res)
    !! Returns \(\log(x+1)\) for real z. TODO, replace with a fast and accurate version
    !! coded natively in fortran !
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: res
    res = c_logp1(x)
  end function logp1r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function logp1z(z) result(res)
    !! Returns \(\log(z+1)\) for complex z. TODO, replace with a fast and accurate version ?
    !! coded natively in fortran !
    use rotex__constants, only: one
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res
    complex(dp) :: zp1
    zp1 = z + one
    res = log(zp1) - (zp1 - one - z)/zp1
  end function logp1z

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isintegerr(x) result(res)
    !! Check if x is an integer
    implicit none
    real(dp), intent(in) :: x
    logical :: res
    res = .false.
    if(x .ne. nint(x)) return
    res = .true.
  end function isintegerr
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isintegerz(z) result(res)
    !! Check if z is an integer
    use rotex__constants, only: zero
    implicit none
    complex(dp), intent(in) :: z
    logical :: res
    res = .false.
    if(z%im .ne. zero) return
    if(z%re .ne. nint(z%re)) return
    res = .true.
  end function isintegerz

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isnaturalr(x) result(res)
    !! Check if x is a natural number (1, 2, 3, ...)
    use rotex__constants, only: zero
    implicit none
    real(dp), intent(in) :: x
    logical :: res
    integer :: nintx
    res = .false.
    nintx = nint(x)
    if(x .ne. nintx) return
    if(nintx .le. 0) return
    res = .true.
  end function isnaturalr
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isnaturalz(z) result(res)
    !! Check if z is a natural number (1, 2, 3, ...)
    use rotex__constants, only: zero
    implicit none
    complex(dp), intent(in) :: z
    logical :: res
    integer :: nintz
    res = .false.
    nintz = nint(z%re)
    if(z%im .ne. zero) return
    if(z%re .ne. nintz) return
    if(nintz .le. 0) return
    res = .true.
  end function isnaturalz

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function argi(z) result(res)
    !! Return the phase of an integer (as a complex number) in (-π,π]
    use rotex__constants, only: pi
    implicit none
    integer, intent(in) :: z
    real(dp) :: res
    res = 0.0_dp
    if(z .ge. 0) return
    res = pi
  end function argi
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function argr(z) result(res)
    !! Return the phase of a real number (as a complex number) in (-π,π]
    use rotex__constants, only: pi
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: res
    res = 0.0_dp
    if(z .ge. 0._dp) return
    res = pi
  end function argr
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function argc(z) result(res)
    !! Return the phase of a complex number in (-π,π]
    use rotex__constants, only: pi, zero
    implicit none
    complex(dp), intent(in) :: z
    real(dp) :: res
    res = zero
    if(res .eq. (zero, zero)) return
    res = atan2(z%im, z%re)
  end function argc

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function invi(i) result(res)
    !! Return the real 1/i for integer i
    implicit none
    integer, intent(in) :: i
    real(dp) :: res
    res = 1/real(i, kind = dp)
  end function invi
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function invr(x) result(res)
    !! Return 1/x
    use rotex__constants, only: zero
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: res
    res = 1/x
  end function invr
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function invz(z) result(res)
    !! Return 1/z for complex z. If the imaginary part of z is 0, flip its sign
    use ieee_arithmetic, only: ieee_is_negative, ieee_copy_sign
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res
    res = 1/z
    if(z%im .ne. 0) return
    if(ieee_is_negative(z%im) .eqv. .true.) return
    res%im = -res%im
  end function invz

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function iseven(n) result(res)
    integer, intent(in) :: n
    logical :: res
    res = .false.
    if(iand(n, 1) .eq. 1) return
    res = .true.
  end function iseven
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isodd(n) result(res)
    integer, intent(in) :: n
    logical :: res
    res = .not. iseven(n)
  end function isodd

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function istriangle(a, b, c) result(res)
    !! returns whether the arguments satisfy the triangle inequality,
    !! assuming that they're positive quantities

    implicit none

    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    logical :: res

    res = .false.

    if(c .gt. a+b) return
    if(a .gt. b+c) return
    if(b .gt. a+c) return

    res = .true.

  end function istriangle

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function log_factorial_i(n) result(res)
    !! returns ln(n!) for integer n
    use rotex__system, only: die
    implicit none
    integer, intent(in) :: n
    real(dp) :: res
    select case(n)
    case(0)  ; res = 0.0_dp
    case(1)  ; res = 1.0_dp
    case(2:) ; res = log_gamma(n + 1.0_dp)
    case default
      call die("ln(n!) not defined for negative integer n")
    end select
  end function log_factorial_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function log_factorial_r(n) result(res)
    !! returns ln(n!) for real n
    use rotex__system, only: die
    implicit none
    real(dp), intent(in) :: n
    real(dp) :: res
    if(n .lt. 0) call die("ln(n!) not defined for negative real n")
    res = log_gamma(n+1)
  end function log_factorial_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function logb_ii(b, x) result(res)
    !! Returns the logarithm of x in the base b
    implicit none
    integer, intent(in) :: b
    integer, intent(in) :: x
    real(dp) :: res
    res = log(real(x, kind = dp))/log(real(b, kind = dp))
  end function logb_ii
  pure elemental function logb_rr(b, x) result(res)
    !! Returns the logarithm of x in the base b
    implicit none
    real(dp), intent(in) :: b
    real(dp), intent(in) :: x
    real(dp) :: res
    res = log(x)/log(b)
  end function logb_rr
  pure elemental function logb_cc(b, x) result(res)
    !! Returns the logarithm of x in the base b
    implicit none
    complex(dp), intent(in) :: b
    complex(dp), intent(in) :: x
    complex(dp) :: res
    res = log(x)/log(b)
  end function logb_cc

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function logrange_ib(a, b, n, base, inclast) result(res)
    !! Returns an array of n log-spaced values from a to b. By default, the base is 10
    !! if omitted but can be changed by the user
    implicit none
    integer(dp), intent(in) :: a, b
    integer,  intent(in) :: n
    integer,  intent(in), optional :: base
      !! Logarithm base
    logical, intent(in), optional :: inclast
      !! Whether to include the last value b
    real(dp) :: res(n)
    logical :: inclast_
    integer :: base_
    inclast_ = .true. ; if(present(inclast)) inclast_ = inclast
    base_ = 10 ; if(present(base)) base_ = base
    res = logrange_rb(real(a, kind = dp), real(b, kind = dp), n, real(base_, kind = dp), inclast_)
  end function logrange_ib
  pure function logrange_rb(a, b, n, base, inclast) result(res)
    !! Returns an array of n log-spaced values from a to b. By default, the base is 10
    !! if omitted but can be changed by the user
    implicit none
    real(dp), intent(in) :: a, b
    integer,  intent(in) :: n
    real(dp), intent(in), optional :: base
      !! Logarithm base
    logical, intent(in), optional :: inclast
      !! Whether to include the last value b
    real(dp) :: res(n)
    logical :: inclast_
    integer  :: i, m
    real(dp) :: base_
    real(dp) :: start, finish
    real(dp) :: dx
    real(dp) :: logx(n)
    base_    = 10._dp ; if(present(base))    base_    = base
    inclast_ = .true. ; if(present(inclast)) inclast_ = inclast
    if(inclast_) then
      m = n
    else
      m = n+1
    endif
    start     = logb(base_, a)
    finish    = logb(base_, b)
    dx        = (finish - start) / (m-1)
    logx      = [( start + i*dx, i=0, m-1 )]
    res       = base_**logx(1:n)
  end function logrange_rb

  ! -------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function neg(i) result(res)
    !! Returns the integer \((-1)^{i)\)}
    implicit none
    integer, intent(in) :: i
    integer :: res
    res = 1
    if(iand(i, 1) .eq. 0) return
    res = -1
  end function neg

! ================================================================================================================================ !
end module rotex__functions
! ================================================================================================================================ !
