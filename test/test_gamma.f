! ================================================================================================================================ !
module test_gamma
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use test_utils,       only: randc, print_params
  use rotex__kinds,     only: dp
  use rotex__polygamma, only: log_gamma
  use rotex__constants, only: macheps_dp

  implicit none

  private

  ! public :: collect_gamma
  public :: collect_log_gamma

  ! -- these (+implementation) need to be refined but it's fine for now
  real(dp), parameter :: rtol_dp = 5e-14_dp
  real(dp), parameter :: atol_dp = 1e-14_dp
  integer, parameter :: ntests_per_test = 100

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine collect_log_gamma(testsuites)
    !! Collect all unit tests to be exported
    implicit none
    type(unittest_type), allocatable, intent(out) :: testsuites(:)
      !! collection of tests
    testsuites = [                                    &
        new_unittest("Log Gamma: z = 1", test_z_eq_1) &
      , new_unittest("Log Gamma: z = n + ½", test_z_is_halfint) &
      ]
  end subroutine collect_log_gamma
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------------------------------------------------------------ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_z_eq_1(error)
    implicit none
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: lhs, rhs
    real(dp), parameter :: z(*) = [1._dp, 2._dp]
    integer :: i
    real(dp), parameter :: atol = atol_dp
    real(dp), parameter :: rtol = rtol_dp
    real(dp) :: thr
    real(dp) :: diff
    do i=1, size(z, 1)
      lhs = log_gamma(z(i))
      rhs = 0._dp
      thr = max(atol, rtol*max(abs(lhs), abs(rhs)))
      call check(error, lhs, rhs, thr = thr, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(z(i))
    enddo
  end subroutine test_z_eq_1

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_z_is_halfint(error)
    use rotex__constants, only: pi
    implicit none
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: lhs, rhs
    real(dp) :: z
    integer :: n, k, nmax
    real(dp), parameter :: atol = atol_dp * 100
    real(dp), parameter :: rtol = rtol_dp  * 100
    real(dp) :: thr
    real(dp), parameter :: sqpi = sqrt(pi)
    real(dp), parameter :: ans(*) = [ &
      log(                           sqpi) &
    , log(1._dp         / 2._dp    * sqpi) &
    , log(3._dp         / 4._dp    * sqpi) &
    , log(15._dp        / 8._dp    * sqpi) &
    , log(105._dp       / 16._dp   * sqpi) &
    , log(945._dp       / 32._dp   * sqpi) &
    , log(10395._dp     / 64._dp   * sqpi) &
    , log(135135._dp    / 128._dp  * sqpi) &
    , log(2027025._dp   / 256._dp  * sqpi) &
    , log(34459425._dp  / 512._dp  * sqpi) &
    , log(654729075._dp / 1024._dp * sqpi) &
    ]
    nmax = size(ans, 1) -1
    do n = 0, nmax
      z = real(n, kind = dp) + 0.5_dp
      lhs = log_gamma(z)
      rhs = ans(n+1)
      thr = max(atol, rtol*max(abs(lhs), abs(rhs)))
      call check(error, lhs, rhs, thr = thr, rel = .false.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(z)
      return
    enddo
  end subroutine test_z_is_halfint

! ================================================================================================================================ !
end module test_gamma
! ================================================================================================================================ !
