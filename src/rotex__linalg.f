! ================================================================================================================================ !
module rotex__linalg
  !! Linear algebra interfaces to LAPACK routines and other linear algebra stuff

  use rotex__kinds, only: dp
! #ifdef WITH_STDLIB
  ! use stdlib_linalg_lapack, only: zgesv  => stdlib_zgesv,  dsyev  => stdlib_dsyev  &
  !                               , zgetrs => stdlib_zgetrs, zgetrf => stdlib_zgetrf
! #endif

  implicit none

  private

  public :: dsyev
  public :: zgesv
  public :: operator(.matmul.)

  public :: right_divide

#ifndef WITH_STDLIB
  interface
    subroutine zgesv(n, nrhs, a, lda, ipiv,b, ldb, info)
      import dp
      implicit none
      integer,     intent(in)    :: lda, ldb, n, nrhs
      integer,     intent(out)   :: info,ipiv(*)
      complex(dp), intent(inout) :: a(lda,*),b(ldb,*)
    end subroutine zgesv
  end interface
  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      implicit none
      character(1), intent(in)    :: jobz, uplo
      real(dp),     intent(inout) :: a(lda, n)
      integer,      intent(in)    :: lda, lwork, n
      integer,      intent(out)   :: info
      real(dp),     intent(out)   :: w(*), work(*)
    end subroutine dsyev
  end interface
  interface
    subroutine zgetrf(m, n, a, lda, ipiv, info)
      import dp
      implicit none
      integer,     intent(in)    :: m, n, lda
      integer,     intent(out)   :: ipiv(*), info
      complex(dp), intent(inout) :: a(lda, *)
    end subroutine zgetrf
  end interface
  interface zgetrs
    subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info )
      import dp
      implicit none
      character(1), intent(in)    :: trans
      integer,      intent(out)   :: info
      integer,      intent(in)    :: lda, ldb, n, nrhs, ipiv(*)
      complex(dp),  intent(in)    :: a(lda,*)
      complex(dp),  intent(inout) :: b(ldb,*)
    end subroutine zgetrs
  end interface zgetrs
#endif

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function right_divide(A, B) result(X)
    !! Returns X = AB⁻¹ without evaluating B⁻¹
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    implicit none
    complex(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp) :: X(size(A, 1), size(A, 2))
    complex(dp), allocatable :: BT(:,:), AT(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n, m, info
    n = size(B, 1)
    m = size(A, 1)
    if(size(B,2) .ne. n) call die("Trying to invert a nonsquare matrix B !")
    if(size(A,2) .ne. n) call die("Can't form AB⁻¹ because the dimensions of A are wrong !")
    AT = transpose(A)
    BT = transpose(B)
    allocate(ipiv(n))
    call zgetrf(n, n, BT, n, ipiv, info)
    if(info .ne. 0) call die("ZGETRF exited with INFO = "//i2c(info))
    ! -- solve X = AB⁻¹ is BT XT = AT, which is solved by ZGETRS (AX=B)
    call zgetrs('T', n, m, BT, n, ipiv, AT, n, info)
    if(info .ne. 0) call die("ZGETRS exited with INFO = "//i2c(info))
    X = transpose(AT)
  end function right_divide

! ================================================================================================================================ !
end module rotex__linalg
! ================================================================================================================================ !
