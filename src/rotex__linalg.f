! ================================================================================================================================ !
module rotex__linalg
  !! Linear algebra interfaces to LAPACK routines and other linear algebra stuff

  use rotex__kinds, only: dp
! #ifdef WITH_STDLIB
  ! use stdlib_linalg_lapack, only: zgesv  => stdlib_zgesv,  dsyev  => stdlib_dsyev  &
  !                               , zgetrs => stdlib_zgetrs, zgetrf => stdlib_zgetrf
! #endif

  implicit none (type, external)

  private

  public :: dsyev
  public :: zgesv
  public :: zheev
  public :: operator(.matmul.)

  public :: right_divide

  interface right_divide
    module procedure :: right_divide_r
    module procedure :: right_divide_c
  end interface right_divide

#ifndef WITH_STDLIB
  interface
    subroutine dpotrf(uplo, n, a, lda, info)
      !! Computes Cholesky factorization of real SPD A
      import dp
      implicit none (type, external)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda
      real(dp),     intent(in)    :: a(lda, *)
      integer,      intent(out)   :: info
    end subroutine dpotrf
  end interface

  interface
    subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
      !! solves AX=B with SPD A using Cholesky factorization
      import dp
      implicit none (type, external)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, nrhs, lda, ldb
      real(dp),     intent(in)    :: a(lda, *)
      real(dp),     intent(inout) :: b(ldb,*)
      integer,      intent(out)   :: info
    end subroutine dpotrs
  end interface

  interface
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      !! LU factorization
      import dp
      implicit none (type, external)
      integer,  intent(in)    :: m, n, lda
      real(dp), intent(inout) :: a(lda, *)
      integer,  intent(out)   :: info, ipiv(*)
    end subroutine dgetrf
  end interface

  interface
    subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
      !! Solves A X = B or AT X = B using LU factorization from dgetrf
      import dp
      implicit none (type, external)
      character(1), intent(in)  :: trans
      integer,      intent(in)  :: n, nrhs, lda, ldb, ipiv(*)
      real(dp),     intent(in)  :: a(lda,*)
      real(dp),     intent(out) :: b(ldb,*)
      integer,      intent(out) :: info
    end subroutine dgetrs
  end interface

  interface
    subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
      !! Factorization of realsymmetric A
      import dp
      implicit none (type, external)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, lda, lwork
      real(dp),     intent(inout) :: a(lda,*)
      real(dp),     intent(out)   :: work(*)
      integer,      intent(out)   :: ipiv(*), info
    end subroutine dsytrf
  end interface

  interface
    subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
      !! Solve AX=B with real symmetric A using the dsytrf factorization
      import dp
      implicit none (type, external)
      character(1), intent(in)    :: uplo
      integer,      intent(in)    :: n, nrhs, lda, ldb
      real(dp),     intent(in)    :: a(lda,*)
      real(dp),     intent(inout) :: b(ldb, *)
      integer,      intent(out)   :: ipiv(*), info
    end subroutine dsytrs
  end interface

  interface
    subroutine zgesv(n, nrhs, a, lda, ipiv,b, ldb, info)
      import dp
      implicit none (type, external)
      integer,     intent(in)    :: lda, ldb, n, nrhs
      integer,     intent(out)   :: info,ipiv(*)
      complex(dp), intent(inout) :: a(lda,*),b(ldb,*)
    end subroutine zgesv
  end interface

  interface
    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      import dp
      implicit none (type, external)
      character(1), intent(in)    :: jobz, uplo
      real(dp),     intent(inout) :: a(lda, n)
      integer,      intent(in)    :: lda, lwork, n
      integer,      intent(out)   :: info
      real(dp),     intent(out)   :: w(*), work(*)
    end subroutine dsyev
  end interface

  interface
    subroutine zheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
      import dp
      implicit none(type, external)
      character(1), intent(in)    :: jobz, uplo
      integer,      intent(in)    :: n
      complex(dp),  intent(inout) :: A(lda, *)
      integer,      intent(in)    :: lda
      real(dp),     intent(out)   :: w(*)
      complex(dp),  intent(out)   :: work(*)
      integer,      intent(in)    :: lwork
      real(dp),     intent(out)   :: rwork(*)
      integer,      intent(out)   :: info
    end subroutine zheev
  end interface

  interface
    subroutine zgetrf(m, n, a, lda, ipiv, info)
      import dp
      implicit none (type, external)
      integer,     intent(in)    :: m, n, lda
      integer,     intent(out)   :: ipiv(*), info
      complex(dp), intent(inout) :: a(lda, *)
    end subroutine zgetrf
  end interface

  interface zgetrs
    subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info )
      import dp
      implicit none (type, external)
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
  module function right_divide_r(A, B) result(X)
    !! Returns X = AB⁻¹ without evaluating B⁻¹ for real matrices
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    use rotex__arrays,     only: is_symmetric
    implicit none (type, external)
    real(dp), intent(in) :: A(:,:), B(:,:)
    real(dp) :: X(size(A, 1), size(A, 2))
    real(dp), allocatable :: BT(:,:), AT(:,:)
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: work(:)
    integer :: n, m, info, lwork

    n = size(B, 1)
    m = size(A, 1)

    if(size(B,2) .ne. n) call die("Trying to invert a nonsquare matrix B !")
    if(size(A,2) .ne. n) call die("Can't form AB⁻¹ because the dimensions of A are wrong !")

    AT = transpose(A)

    if(is_symmetric(B)) then
      BT = B
      ! -- try SPD, cholesky
      call dpotrf('L', n, BT, n, info)
      if(info .eq. 0) then
        call dpotrs('L', n, m, BT, n, AT, n, info)
        if(info .ne. 0) call die("DPOTRS exited with nonzero INFO")
      ! -- not SPD
      else
        allocate(ipiv(n))

        ! -- workspace query
        lwork = -1
        allocate(work(1))
        call dsytrf('L', n, BT, n, ipiv, work, lwork, info)
        if(info .ne. 0) call die("DSYTRF exited with nonzero INFO on workspace query")

        lwork = nint(work(1))
        deallocate(work) ; allocate(work(lwork))

        call dsytrf('L', n, BT, n, ipiv, work, lwork, info)
        if(info .ne. 0) call die("DSYTRF exited with nonzero INFO: "//i2c(INFO))
        call dsytrs('L', n, m, BT, n, ipiv, AT, n, info)
        if(info .ne. 0) call die("DSYTRS exited with nonzero INFO: "//i2c(INFO))
      endif

    else
      ! -- general B
      BT = transpose(B)
      allocate(ipiv(n))

      call dgetrf(n, n, BT, n, ipiv, info)
      if(info .ne. 0) call die("DGETRF exited with nonzero info")
      call dgetrs('N', n, m, BT, n, ipiv, AT, n, info)
      if(info .ne. 0) call die("DGETRS exited with nonzero info")

    endif

    X = transpose(AT)

  end function right_divide_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module function right_divide_c(A, B) result(X)
    !! Returns X = AB⁻¹ without evaluating B⁻¹ for complex matrices
    use rotex__system,     only: die
    use rotex__characters, only: i2c => int2char
    implicit none (type, external)
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
  end function right_divide_c

! ================================================================================================================================ !
end module rotex__linalg
! ================================================================================================================================ !
