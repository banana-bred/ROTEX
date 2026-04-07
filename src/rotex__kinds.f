! ================================================================================================================================ !
module rotex__kinds
  !! Contains kind definitions (double, xtended double, quad precision)
  use, intrinsic :: iso_fortran_env, only: int32, int64
  implicit none (type, external)
  private
  public :: int32
  public :: int64
  integer, parameter, public :: sp  = selected_real_kind(6)
  integer, parameter, public :: dp  = selected_real_kind(15)
  integer, parameter, public :: xdp = selected_real_kind(18)
  integer, parameter, public :: qp  = selected_real_kind(33)
#if defined(PROB_SINGLE) && defined(PROB_DOUBLE)
#error Define only one of PROB_SINGLE or PROB_DOUBLE
#elif defined(PROB_SINGLE)
  integer, parameter, public :: prob_rk = sp
#else
  integer, parameter, public :: prob_rk = dp
#endif
! ================================================================================================================================ !
end module rotex__kinds
! ================================================================================================================================ !
