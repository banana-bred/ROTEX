! ================================================================================================================================ !
module test_hypergeometric
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use rotex__kinds,     only: dp
  use rotex__constants, only: zero, one, two, three, four, pi, macheps => macheps_dp
  use rotex__hypergeometric, only: f21
  use test_utils, only: randr, randc, print_params

  implicit none

  private

  public :: collect_hypergeometric_2f1

  real(dp), parameter :: ATOL = 1e-14_dp
  real(dp), parameter :: RTOL = 5e-13_dp
  real(dp), parameter :: ABC_SCALE = 30._dp
  integer, parameter :: ntests_per_test = 1000

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine collect_hypergeometric_2f1(testsuites)
    !! Collect all unit tests to be exported

    implicit none

    type(unittest_type), allocatable, intent(out) :: testsuites(:)
      !! collection of tests

    testsuites = [                                    &
                   new_unittest("z = 0 (real)", test_z_eq_0_r) &
                 , new_unittest("z = 0 (cmplx)", test_z_eq_0_c) &
                 , new_unittest("F( 1,   1  ;   2; z ) = -z⁻¹ln(1-z)", test_dlmf_15_4_1) &
                 , new_unittest("F( 1/2, 1  ; 3/2; z²) = ln([1+z]/[1-z])/2z", test_dlmf_15_4_2) &
                 , new_unittest("F( 1/2, 1  ; 3/2;-z²) = atan(z)/z", test_dlmf_15_4_3) &
                 , new_unittest("F( 1/2, 1/2; 3/2; z²) = asin(z)/z", test_dlmf_15_4_4) &
                 , new_unittest("F( 1/2, 1/2; 3/2; z ) = ln(z+√(1+z²))/z", test_dlmf_15_4_5) &
                 , new_unittest("F( a  , b  ; a  ; z ) = (1-z)^(-b)", test_dlmf_15_4_6a) &
                 , new_unittest("F( a  , b  ; b  ; z ) = (1-z)^(-a)", test_dlmf_15_4_6b) &
                 ! , new_unittest("test_2f1_finite")
      ]

  end subroutine collect_hypergeometric_2f1

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_z_eq_0_c(error)
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    complex(dp) :: a, b, c, ans, tol
    real(dp) :: z
    z = 0._dp
    ans = (1._dp, 0._dp)
    tol = atol + rtol*abs(ans)
    do i = 1, ntests_per_test
      call randc(a, 100._dp, 100._dp)
      call randc(b, 100._dp, 100._dp)
      call randc(c, 100._dp, 100._dp)
      call check(error, f21(a, b, c, z), ans, thr = macheps, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_z_eq_0_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_z_eq_0_r(error)
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, c, z, ans, tol
    z = 0._dp
    ans = 1._dp
    tol = atol + rtol*abs(ans)
    do i = 1, ntests_per_test
      call randr(a, 100._dp)
      call randr(b, 100._dp)
      call randr(c, 100._dp)
      call check(error, f21(a, b, c, z), ans, thr = macheps, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_z_eq_0_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_1(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, c, z, ans, tol
    a=1;b=1;c=2
    do i = 1, ntests_per_test
      call randr(z)
      ans = -logp1(-z)/z
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, c, z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_dlmf_15_4_1

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_2(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, c, z, ans, tol
    a=1._dp/2._dp;b=1;c=3._dp/2._dp
    do i = 1, ntests_per_test
      call randr(z)
      ans = log((1+z)/(1-z))/(2*z)
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, c, z*z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_dlmf_15_4_2

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_3(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, c, z, ans, tol
    a=1._dp/2._dp;b=1;c=3._dp/2._dp
    do i = 1, ntests_per_test
      call randr(z)
      ans = atan(z)/z
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, c, -z*z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_dlmf_15_4_3

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_4(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, c, z, ans, tol
    a=1._dp/2._dp;b=1._dp/2._dp;c=3._dp/2._dp
    do i = 1, ntests_per_test
      call randr(z)
      ans = asin(z)/z
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, c, z*z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_dlmf_15_4_4

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_5(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, c, z, ans, tol
    a=1._dp/2._dp;b=1._dp/2._dp;c=3._dp/2._dp
    do i = 1, ntests_per_test
      call randr(z)
      ans = log(z + sqrt(1._dp+z*z))/z
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, c, -z*z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, c, z)
      return
    enddo
  end subroutine test_dlmf_15_4_5

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_6a(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, z, ans, tol
    do i = 1, ntests_per_test
      call randr(z)
      call randr(a, ABC_SCALE)
      call randr(b, ABC_SCALE)
      ans = (1-z)**(-b)
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, a, z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, a, z)
      return
    enddo
  end subroutine test_dlmf_15_4_6a

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_dlmf_15_4_6b(error)
    use rotex__functions, only: logp1
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp) :: a, b, z, ans, tol
    do i = 1, ntests_per_test
      call randr(z)
      call randr(a, ABC_SCALE)
      call randr(b, ABC_SCALE)
      ans = (1-z)**(-a)
      tol = atol + rtol*abs(ans)
      call check(error, f21(a, b, b, z), ans, thr = tol, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(a, b, b, z)
      return
    enddo
  end subroutine test_dlmf_15_4_6b


! ================================================================================================================================ !
end module test_hypergeometric
! ================================================================================================================================ !
