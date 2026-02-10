! ================================================================================================================================ !
program tester

  use testdrive,    only: run_testsuite, new_testsuite, testsuite_type
  use test_digamma, only: collect_digamma
  use test_gamma,   only: collect_log_gamma!, collect_gamma
  use test_hypergeometric, only: collect_hypergeometric_2f1
  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
                 ! new_testsuite("Confluent Hypergeometric Function", collect_hypergeometric_1f1) &
               ! , new_testsuite("Gauss Hypergeometric Function", collect_hypergeometric_2f1) &
      new_testsuite("Digamma function ψ(z)", collect_digamma)             &
    , new_testsuite("Log Gamma function ln(Γ(z))", collect_log_gamma)             &
    ! , new_testsuite("Gamma function Γ(z)", collect_gamma)             &
    , new_testsuite("₂F₁ hypergeometric function", collect_hypergeometric_2f1) &
    ]

  do is = 1, size(testsuites)
    write(error_unit, *)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat .gt. 0) then
    write(error_unit, *)
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

! ================================================================================================================================ !
end program tester
! ================================================================================================================================ !
