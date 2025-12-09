! ================================================================================================================================ !
module rotex__splines
  !! Wrapper routines for the spline fitting procedures

  implicit none

  private

  public :: interpolate_replace

  integer, parameter :: SPLINE_ORDER_KX = 3
  integer, parameter :: DB1INK_IKNOT = 0
  integer, parameter :: DB1VAL_IDX = 0

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure module subroutine interpolate_replace(xold, xnew, fx, idxx)
    !! Interpolate f(xold) -> f(xnew). Only consider xnew values that are contained within xold for now.
    !! Xnew returns untouched, but fx is overwritten with interpolated values on the
    !! grid xnew(idxx)

    use rotex__kinds,   only: dp
    use rotex__arrays,  only: size_check
    use bspline_module, only: db1ink, db1val

    implicit none

    real(dp), intent(in) :: xold(:)
      !! The grid of values on which our function has been evaluated
    real(dp), intent(in) :: xnew(:)
      !! The new grid of values that we want our function to be evaluated on.
    real(dp), intent(inout), allocatable :: fx(:)
      !! On input:  evaluated function f(xold)
      !! On output: the evaluated function f(xnew)
    integer, intent(out), allocatable, optional :: idxx(:)
      !! The indices of values of xnew that are used to evaluate fx

    logical, allocatable :: mask(:)
    integer :: nxold, nxnew_chopped, iflag, inbvx, ix
    real(dp), allocatable :: tx(:), bcoef(:), w0(:), fx_copy(:), xnew_chopped(:)

    nxold = size(xold, 1)
    call size_check(fx, nxold, "FX")

    mask = xnew .ge. xold(1) .AND. xnew .le. xold(nxold)

    ! -- get x-subgrid
    xnew_chopped  = pack(xnew, mask)
    nxnew_chopped = size(xnew_chopped, 1)

    if(present(idxx)) then
      idxx = [(ix, ix=1, size(xnew, 1))]
      idxx = pack(idxx, mask)
    endif

    allocate(tx(nxold+SPLINE_ORDER_KX), source = 0.0_dp)
    allocate(bcoef(nxold),              source = 0.0_dp)
    allocate(fx_copy(nxnew_chopped),    source = 0.0_dp)
    allocate(w0(3*SPLINE_ORDER_KX),     source = 0.0_dp)

    ! -- interpolate
    call db1ink(xold, nxold, fx, SPLINE_ORDER_KX, DB1INK_IKNOT, tx, bcoef, iflag)

    ! -- evaluate
    inbvx = 1
    do ix = 1, nxnew_chopped
      call db1val(         &
          xnew_chopped(ix) &
        , DB1VAL_IDX       &
        , tx               &
        , nxold            &
        , SPLINE_ORDER_KX  &
        , bcoef            &
        , fx_copy(ix)      &
        , iflag            &
        , inbvx            &
        , w0               &
        , extrap = .false. &
      )
    enddo

    call move_alloc(fx_copy, fx)

  end subroutine interpolate_replace

! ================================================================================================================================ !
end module rotex__splines
! ================================================================================================================================ !
