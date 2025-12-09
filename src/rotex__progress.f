! ================================================================================================================================ !
module rotex__progress
  !! Interface to [forbear](https://github.com/szaghi/forbear), the Fortran (progress) B(e)ar environment
#ifdef USE_FORBEAR
  use forbear, only: progressbar_type => bar_object
#endif

  implicit none

  public

! ================================================================================================================================ !
end module rotex__progress
! ================================================================================================================================ !
