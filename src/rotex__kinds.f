! ================================================================================================================================ !
module rotex__kinds
  !! Contains kind definitions (double, xtended double, quad precision)
  implicit none
  private
  integer, parameter, public :: dp  = selected_real_kind(15)
  integer, parameter, public :: xdp = selected_real_kind(18)
  integer, parameter, public :: qp  = selected_real_kind(33)
! ================================================================================================================================ !
end module rotex__kinds
! ================================================================================================================================ !
