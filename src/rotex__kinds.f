! ================================================================================================================================ !
module rotex__kinds
  !! Contains kind definitions (double, xtended double, quad precision)
  use, intrinsic :: iso_fortran_env, only: int32, int64
  implicit none (type, external)
  private
  public :: int32
  public :: int64
  integer, parameter, public :: dp  = selected_real_kind(15)
  integer, parameter, public :: xdp = selected_real_kind(18)
  integer, parameter, public :: qp  = selected_real_kind(33)
! ================================================================================================================================ !
end module rotex__kinds
! ================================================================================================================================ !
