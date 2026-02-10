! ================================================================================================================================ !
module test_digamma
  use testdrive,        only: new_unittest, unittest_type, error_type, check
  use test_utils,       only: randr, print_params
  use rotex__kinds,     only: dp
  use rotex__polygamma, only: digamma
  use rotex__constants, only: macheps_dp
  ! use hypergeometric__constants, only: zero, one, two, three, four, pi

  implicit none

  private

  public :: collect_digamma

  real(dp), parameter :: rtol_dp = sqrt(macheps_dp)
  ! real(qp), parameter :: rtol_qp = 2000 * macheps_qp
  ! real(dp), parameter :: atol_dp = 10*macheps_dp
  ! real(qp), parameter :: atol_qp = 10*macheps_qp
  integer, parameter :: ntests_per_test = 100

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine collect_digamma(testsuites)
    !! Collect all unit tests to be exported
    implicit none
    type(unittest_type), allocatable, intent(out) :: testsuites(:)
      !! collection of tests
    testsuites = [                                                 &
        new_unittest("Digamma : z = 1", test_digamma_z_eq_1) &
      , new_unittest("Digamma : z = n + 1/2", test_digamma_z_is_halfint) &
      , new_unittest("Digamma : z → z + 1 recurrence", test_digamma_recurr) &
      , new_unittest("Digamma : z → 1 - z reflection", test_digamma_reflect) &
      , new_unittest("Digamma : ψ(nz) multiplication", test_digamma_multiplication) &
      ! , new_unittest("Digamma : ψ(z*) = ψ(z)* symmetry", test_digamma_conjg) &
                 ! , new_unittest("gamma",   test_gamma) &
                 ! ,   new_unittest("real-valued gamma", test_gamma_r) &
                 ! , new_unittest("complex-valued gamma", test_gamma_c) &
                 ! , new_unittest("log gamma", test_log_gamma) &
      ]
  end subroutine collect_digamma
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------------------------------------------------------------ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_digamma_z_eq_1(error)
    !! Test the identity ψ(1) = -γ
    use rotex__constants, only: euler_mascheroni
    implicit none
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: z, lhs, rhs
    z = real(1, kind = dp)
    lhs = digamma(z)
    rhs = real(-euler_mascheroni, kind = dp)
    call check(error, lhs, rhs, thr = rtol_dp, rel = .true.)
    if(allocated(error) .eqv. .false.) return
    call print_params(z)
  end subroutine test_digamma_z_eq_1

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_digamma_z_is_halfint(error)
    !! Test the identity ψ(n+½)
    use rotex__constants, only: euler_mascheroni
    implicit none
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: z, lhs, rhs
    integer, parameter :: nmax = 15
    integer, parameter :: nstep = 3
    integer :: n, k
    real(dp) :: nre
    do n = 0, nmax, nstep
      nre = real(n, kind = dp)
      z = nre + 0.5_dp
      lhs = digamma(z)
      rhs = -euler_mascheroni - 2._dp*log(2._dp)
      do k = 1, n
        rhs = rhs + 2._dp / real(2*k-1, kind = dp)
      enddo
      call check(error, lhs, rhs, thr = rtol_dp, rel = .true.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(z)
      return
    enddo
  end subroutine test_digamma_z_is_halfint

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_digamma_recurr(error)
    !! Test the recursion relations ψ(z+z) = ψ(z) + 1/z
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: k
    real(dp) :: z, lhs, rhs
    do k = 1, ntests_per_test
      call randr(z, 100._dp)
      lhs = digamma(z + 1._dp)
      rhs = digamma(z) + 1._dp/z
      call check(error, lhs, rhs, thr = rtol_dp, rel = .true.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(z)
      return
    enddo
  end subroutine test_digamma_recurr

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_digamma_reflect(error)
    !! Test the reflection relation ψ(1-z)- ψ(z) = πcot(πz)
    use rotex__constants, only: pi
    use rotex__functions, only: cotpi
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer :: k
    real(dp) :: z, lhs, rhs
    do k = 1, ntests_per_test
      call randr(z, 0.6_dp)
      lhs = digamma(1._dp - z) - digamma(z)
      rhs = pi * cotpi(z)
      call check(error, lhs, rhs, thr = rtol_dp, rel = .true.)
      if(allocated(error) .eqv. .false.) cycle
      call print_params(z)
      return
    enddo
  end subroutine test_digamma_reflect

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine test_digamma_multiplication(error)
    !! Test the multiplication theorem
    implicit none
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: nmax = 10
    integer :: k, n, i
    real(dp) :: z, lhs, rhs
    real(dp) :: nc, ic
    do n = 1, nmax
      nc = real(n, kind = dp)
      do k = 1, ntests_per_test
        call randr(z, 0.6_dp)
        lhs = digamma(nc*z)
        rhs = log(nc)
        do i=0, n-1
          ic = real(i, kind = dp)
          rhs = rhs + digamma(z + ic/nc) / nc
        enddo
        call check(error, lhs, rhs, thr = rtol_dp, rel = .true.)
        if(allocated(error) .eqv. .false.) cycle
        call print_params(z)
        return
      enddo
    enddo
  end subroutine test_digamma_multiplication

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine test_digamma_conjg(error)
  !   !! Test the conjugation symmetry ψ(z)* = ψ(z*)
  !   implicit none
  !   type(error_type), allocatable, intent(out) :: error
  !   complex(dp) :: z, lhs, rhs
  !   integer :: k
  !   do k = 1, ntests_per_test
  !     call randc(z, 0.6_dp, 0.6_dp)
  !     lhs = digamma(conjg(z))
  !     rhs = conjg(digamma(z))
  !     call check(error, lhs, rhs, thr = rtol_dp, rel = .true.)
  !     if(allocated(error) .eqv. .false.) cycle
  !     call print_params(z)
  !     return
  !   enddo
  ! end subroutine test_digamma_conjg

! ================================================================================================================================ !
end module test_digamma
! ================================================================================================================================ !
