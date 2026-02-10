
module test_utils

  use rotex__kinds,  only: dp, qp
  use rotex__system, only: die
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit

  implicit none

  private

  public :: randc
  public :: randr
  public :: isnear
  public :: tolerance
  public :: print_params
  public :: print_z

  character(14), parameter, public :: inum_fmt = '(A, " = ", I0)'
  character(18), parameter, public :: rnum_fmt = '(A, " = ", E25.16)'
  character(42), parameter, public :: cnum_fmt = '(A, " = ", E25.16, " + ", E25.16, " * im")'

  interface print_params
    module procedure :: print_z
    module procedure :: print_abcz_rr
    module procedure :: print_abcz_cc
    module procedure :: print_abcz_cr
  end interface print_params

  interface randr
    module procedure :: randr_dp
  end interface randr

  interface randc
    module procedure :: randc_dp
    module procedure :: randc_qp
  end interface randc

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function isnear(x, y, atol, rtol) result(res)
    !! Test if x ≈ y
    implicit none
    complex(dp), intent(in) :: x, y
    real(dp),    intent(in) :: atol, rtol
    logical :: res
    res = abs(x-y) .le. max(atol, rtol*max(abs(x), abs(y)))
  end function isnear

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental module subroutine randr_dp(x, scl)
    !! Generate a random number x with optional scaling: [-1,1]*scl
    !! By default, scaling is just 1
    implicit none
    real(dp), intent(out) :: x
    real(dp), intent(in), optional :: scl
    real(dp) :: a, sr
    sr = 1 ; if(present(scl)) sr = scl
    call random_number(a)
    x = (2*a-1)*sr
  end subroutine randr_dp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental module subroutine randc_dp(z, scale_re, scale_im)
    !! Generate a random number z with optional real and imaginary scaling
    implicit none
    complex(dp), intent(out) :: z
    real(dp), intent(in), optional :: scale_re, scale_im
    real(dp) :: a, b, sr, si
    sr = 1 ; if(present(scale_re)) sr = scale_re
    si = 1 ; if(present(scale_im)) si = scale_im
    call random_number(a)
    call random_number(b)
    z = cmplx((2*a-1)*sr, (2*b-1)*si, kind = dp)
  end subroutine randc_dp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental module subroutine randc_qp(z, scale_re, scale_im)
    !! Generate a random number z with optional real and imaginary scaling
    implicit none
    complex(qp), intent(out) :: z
    real(qp), intent(in), optional :: scale_re, scale_im
    real(qp) :: a, b, sr, si
    sr = 1 ; if(present(scale_re)) sr = scale_re
    si = 1 ; if(present(scale_im)) si = scale_im
    call random_number(a)
    call random_number(b)
    z = cmplx((2*a-1)*sr, (2*b-1)*si, kind = qp)
  end subroutine randc_qp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function tolerance() result(res)
    !! Returns 100 ε
    use rotex__constants, only: macheps => macheps_dp
    implicit none
    real(dp) :: res
    res = 100*macheps
  end function tolerance

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine print_abcz_rr(a, b, c, z)
    !! Prints the parameters a, b, c, z, and info about them
    implicit none
    real(dp), intent(in) :: a, b, c, z
    write(stderr, rnum_fmt) "a", a
    write(stderr, rnum_fmt) "b", b
    write(stderr, rnum_fmt) "c", c
    write(stderr, rnum_fmt) "z", z
  end subroutine print_abcz_rr
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine print_abcz_cr(a, b, c, z)
    !! Prints the parameters a, b, c, z, and info about them
    implicit none
    complex(dp), intent(in) :: a, b, c
    real(dp), intent(in) :: z
    write(stderr, cnum_fmt) "a", a
    write(stderr, cnum_fmt) "b", b
    write(stderr, cnum_fmt) "c", c
    write(stderr, rnum_fmt) "z", z
    write(stderr, rnum_fmt) "|a|", abs(a)
    write(stderr, rnum_fmt) "|b|", abs(b)
    write(stderr, rnum_fmt) "|c|", abs(c)
  end subroutine print_abcz_cr
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine print_abcz_cc(a, b, c, z)
    !! Prints the parameters a, b, c, z, and info about them
    implicit none
    complex(dp), intent(in) :: a, b, c, z
    write(stderr, cnum_fmt) "a", a
    write(stderr, cnum_fmt) "b", b
    write(stderr, cnum_fmt) "c", c
    write(stderr, cnum_fmt) "z", z
    write(stderr, rnum_fmt) "|a|", abs(a)
    write(stderr, rnum_fmt) "|b|", abs(b)
    write(stderr, rnum_fmt) "|c|", abs(c)
    write(stderr, rnum_fmt) "|z|", abs(z)
  end subroutine print_abcz_cc

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine print_z(z)
    !! Prints the parameter z, and info about it
    implicit none
    class(*), intent(in) :: z
    select type(z)
    type is (integer)
      write(stderr, inum_fmt) "z", z
      write(stderr, inum_fmt) "|z|", abs(z)
    type is (real(dp))
      write(stderr, rnum_fmt) "z", z
      write(stderr, rnum_fmt) "|z|", abs(z)
    type is (real(qp))
      write(stderr, rnum_fmt) "z", z
      write(stderr, rnum_fmt) "|z|", abs(z)
    type is (complex(dp))
      write(stderr, cnum_fmt) "z", z
      write(stderr, rnum_fmt) "|z|", abs(z)
    type is (complex(qp))
      write(stderr, cnum_fmt) "z", z
      write(stderr, rnum_fmt) "|z|", abs(z)
    class default
      call die("Z is not an integer, real, or complex or any expected kind")
    end select
  end subroutine print_z


! ================================================================================================================================ !
end module test_utils
! ================================================================================================================================ !
