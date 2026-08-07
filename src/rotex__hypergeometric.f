! ================================================================================================================================ !
module rotex__hypergeometric
  !! For calculating the Gauss hypergeometric function ₂F₁(a,b;c;z) and Appel hypergeometric function F₂(a,b₁,b₂;c₁,c₂,x,y)
  use rotex__kinds,     only: dp, qp
  use rotex__constants, only: macheps_dp
  use rotex__system,    only: stderr, die

  implicit none (type, external)

  private

  public :: f21
  public :: appell_f2_xy_ac1

  interface f21
    !! Interface for calculating the Gauss hypergeometric function ₂F₁(a,b;c;z)
    !! for real and complex parameters
    module procedure :: f21_r
    module procedure :: f21_c
  end interface f21

  real(dp), parameter :: TS_DEFAULT_TOLERANCE = macheps_dp

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function f21_r(a, b, c, x) result(res)
    !! Returns one of the following transforms for real a, b, c, x. These transformations are not valid for
    !! integral b-a or c-a-b
    !!   1. ₂F₁(a,b;c;x) = (1-x)^{-b} ₂F₁(b,c-a;c;x/(x-1))
    !!   2. ₂F₁(a,b;c;x) = (1-x)^{-a} ₂F₁(a,c-b;c;x/(x-1))
    !!   3. ₂F₁(a,b;c;x) = (1-x)^{-a} Γ(c)Γ(b-a)/(Γ(b)Γ(c-a)) ₂F₁(a,c-b;a-b+1;1/(1-x))
    !!                   + (1-x)^{-b} Γ(c)Γ(a-b)/(Γ(a)Γ(c-b)) ₂F₁(b,c-a;b-a+1;1/(1-x))
    !!   4. ₂F₁(a,b;c;x) =               Γ(c)Γ(c-a-b)/Γ(c-a)Γ(c-b) ₂F₁(a, b, a+b-c+1, 1-x)
    !!                   + (1-x)^(c-a-b) Γ(c)Γ(a+b-c)/Γ(a)Γ(b)     ₂F₁(c-a,c-b;c-a-b+1;1-x)
    !!   5. ₂F₁(a,b;c;x) = (x)^(-a)              Γ(c)Γ(c-a-b)/Γ(c-a)Γ(c-b) ₂F₁(a,a-c+1;a+b-c+1;1-1/x)
    !!                   + x^(a-c) (1-x)^(c-a-b) Γ(c)Γ(a+b-c)/Γ(a)Γ(b)     ₂F₁(c-a,1-a;c-a-b+1;1-1/x)
    !!   6. ₂F₁(a,b;c;x) = (-x)^(-a) Γ(c)Γ(b-a)/Γ(b)Γ(c-a) ₂F₁(a,a-c+1;a-b+1;1/x)
    !!                   + (-x)^(-b) Γ(c)Γ(a-b)/Γ(a)Γ(c-b) ₂F₁(b-c+1,b;b-a+1;1/x)
    !! Regions of validity:
    !!   1. |a| < |b|, -1 ≤ x < 0
    !!   2. |a| > |b|, -1 ≤ x < 0
    !!   3. -∞ < x < -1
    !!   4. ½ < x < 1
    !!   5. 1 < x ≤ 2
    !!   6. 2 < x < ∞

    use rotex__utils,     only: isint
    use rotex__system,    only: die
    use rotex__functions, only: inv, isinteger, arg
    use rotex__constants, only: pi
    use rotex__polygamma, only: lgamma => log_gamma

    implicit none (type, external)
    real(dp), intent(in) :: a, b, c, x
    real(dp) :: res

    integer :: m
    real(dp) :: wx

    if(isint(c) .eqv. .true.) then
      if(anint(c) .lt. 1._dp) call die("Hypergeometric function not not defined for c = 0, -1, -2, ..")
    endif

    ! -- transform closer to 0 for better convergence
    if(x .gt. 2.0_dp) then
      ! -- 2 < x < ∞
      if(abs(arg(1-x)) .ge. pi) call die("The transformation for 2<x<∞ is not defined for |arg(1-x)| ≥ π")
      wx = 1._dp/x
      res = (-x)**(-a) * gamma(c)*gamma(b-a)/gamma(b)/gamma(c-a) * f21(a,a-c+1,a-b+1,wx) &
          + (-x)**(-b) * gamma(c)*gamma(a-b)/gamma(a)/gamma(c-b) * f21(b-c+1,b,b-a+1,wx)
      return
    elseif(x .gt. 1.0_dp) then
      ! -- 1 < x ≤ 2
      if(abs(arg(1-x)) .ge. pi) call die("The transformation for 1<x≤2 is not defined for |arg(1-x)| ≥ π")
      wx = 1._dp - 1._dp/x
      res = x**(-a)*gamma(c)*gamma(c-a-b)/gamma(c-a)/gamma(c-b)                 * f21_dispatch_r(a,a-c+1,a+b-c+1,wx) &
          + x**(a-c)*(1._dp-x)**(c-a-b)*gamma(c)*gamma(a+b-c)/gamma(a)/gamma(b) * f21_dispatch_r(c-a,1-a,c-a-b+1,wx)
      return
    elseif(x .gt. 0.5_dp) then
      ! -- ½ < x ≤ 1
      if(isinteger(c-a-b)) then
        m = nint(c-a-b)
        if(c-a-b .ge. 0) then
          ! -- DLMF 15.8.10
          res = gamma(c) * olver_2f1_dlmf_15_8_10(a, b, m, x)
          return
        else
          res = gamma(c) * olver_2f1_dlmf_15_8_12(a, b, abs(m), x)
          return
        endif
      endif
      wx = 1._dp - x
      res =                gamma(c)*gamma(c-a-b)/gamma(c-a)/gamma(c-b) * f21_dispatch_r(a,  b,  a+b-c+1,    wx) &
          + (1-x)**(c-a-b)*gamma(c)*gamma(a+b-c)/gamma(a)/gamma(b)     * f21_dispatch_r(c-a,c-b,c-a-b+1._dp,wx)
      return
    elseif(x .ge. 0.0_dp) then
      ! -- 0 < x ≤ ½
      res = f21_dispatch_r(a, b, c, x)
      return
    endif

    if(x .ge. -1) then
      ! -- -1 ≤ x < 0
      wx = x/(x-1)
      if(abs(a) .lt. abs(b)) then
        res = (1-x)**(-b) * f21_dispatch_r(b, c-a, c, wx)
      else
        res = (1-x)**(-a) * f21_dispatch_r(a, c-b, c, wx)
      endif
      return
    endif

    ! -- -∞ < x < -1
    wx = inv(1._dp-x)
    res = wx**(a)*gamma(c)*gamma(b-a)/gamma(b)/gamma(c-a)*f21_dispatch_r(a,c-b,a-b+1,wx) &
        + wx**(b)*gamma(c)*gamma(a-b)/gamma(a)/gamma(c-b)*f21_dispatch_r(b,c-a,b-a+1,wx)

  end function f21_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function f21_c(a, b, c, x) result(res)
    !! Returns one of the following transforms for complex a, b, c; real x
    !!   1. ₂F₁(a,b;c;x) = (1-x)^{-b} ₂F₁(b,c-a;c;x/(x-1))
    !!   2. ₂F₁(a,b;c;x) = (1-x)^{-a} ₂F₁(a,c-b;c;x/(x-1))
    !!   3. ₂F₁(a,b;c;x) = (1-x)^{-a} Γ(c)Γ(b-a)/(Γ(b)Γ(c-a)) ₂F₁(a,c-b;a-b+1;1/(1-x))
    !!                   + (1-x)^{-b} Γ(c)Γ(a-b)/(Γ(a)Γ(c-b)) ₂F₁(b,c-a;b-a+1;1/(1-x))
    !!   4. ₂F₁(a,b;c;x) =               Γ(c)Γ(c-a-b)/Γ(c-a)Γ(c-b) ₂F₁(a, b, a+b-c+1, 1-x)
    !!                   + (1-x)^(c-a-b) Γ(c)Γ(a+b-c)/Γ(a)Γ(b)     ₂F₁(c-a,c-b;c-a-b+1;1-x)
    !! Regions of validity:
    !!   1. |a| < |b|, -1 ≤ x < 0
    !!   2. |a| > |b|, -1 ≤ x < 0
    !!   3. -∞ < x < -1
    !!   4. ½ < x < 1

    use rotex__utils,     only: isint
    use rotex__system,    only: die
    use rotex__functions, only: inv
    use rotex__polygamma, only: lgamma => log_gamma

    implicit none (type, external)
    complex(dp), intent(in) :: a, b, c
    real(dp),    intent(in) :: x
    complex(dp) :: res
    complex(dp) :: wx
    complex(dp) :: zx

    if(isint(c) .eqv. .true.) then
      if(nint(c%re) .lt. 1) call die("Hypergeometric function not not defined for c = 0, -1, -2, ..")
    endif

    ! if(x .ge. 0) call die("Hypergeometric function got x > 0, which shouldn't happen when using complex parameters.&
    !   & x>0 should only happen in e-neutral scattering.")

    zx = cmplx(x, kind = dp)

    ! -- transform closer to 0 for better convergence
    if(x .gt. 2._dp) then
      wx = 1._dp / zx
      res = exp((-a)*log(-zx)+lgamma(c)+lgamma(b-a)-lgamma(b)-lgamma(c-a)) * f21_dispatch_c(a,a-c+1,a-b+1,wx) &
          + exp((-b)*log(-zx)+lgamma(c)+lgamma(a-b)-lgamma(a)-lgamma(c-b)) * f21_dispatch_c(b-c+1,b,b-a+1,wx)
      return
    elseif(x .gt. 1) then
      wx = 1._dp - 1._dp/x
      res = exp((-a)*log(zx)+lgamma(c)+lgamma(c-a-b)-lgamma(c-a)-lgamma(c-b)) * f21_dispatch_c(a,a-c+1,a+b-c+1,wx) &
          + exp((a-c)*log(zx)+(c-a-b)*log(1-zx)+lgamma(c)+lgamma(a+b-c)-lgamma(a)-lgamma(b)) * f21_dispatch_c(c-a,1-a,c-a-b+1,wx)
      return
    elseif(x .gt. 0.5_dp) then
      wx = 1._dp - zx
      res = exp(lgamma(c) + lgamma(c-a-b) - lgamma(c-a) - lgamma(c-b))           * f21_dispatch_c(a,b,a+b-c+1, wx) &
          + exp((c-a-b)*log(1-x)+lgamma(c)+lgamma(a+b-c)-lgamma(a)-lgamma(b)) * f21_dispatch_c(c-a,c-b,c-a-b+1._dp,wx)
      return
    ! -- 0 < x ≤ ½
    elseif(x .ge. 0.0_dp) then
      ! res = michelf21(a, b, c, zx)
      res = f21_dispatch_c(a, b, c, zx)
      return
    endif

    ! -- -1 ≤ x < 0
    if(x .ge. -1) then
      wx = cmplx(x/(x-1), 0.0_dp, kind=dp)
      if(abs(a) .lt. abs(b)) then
        ! res = (1-zx)**(-b) * michelf21(b, c-a, c, wx)
        res = (1-zx)**(-b) * f21_dispatch_c(b, c-a, c, wx)
      else
        ! res = (1-zx)**(-a) * michelf21(a, c-b, c, wx)
        res = (1-zx)**(-a) * f21_dispatch_c(a, c-b, c, wx)
      endif
      return
    endif

    ! -- -∞ < x < -1
    wx = cmplx(inv(1-x), 0.0_dp, kind=dp)
    res = exp(-a*log(1-x)+lgamma(c)+lgamma(b-a)-lgamma(b)-lgamma(c-a))*f21_dispatch_c(a,c-b,a-b+1,wx) &
        + exp(-b*log(1-x)+lgamma(c)+lgamma(a-b)-lgamma(a)-lgamma(c-b))*f21_dispatch_c(b,c-a,b-a+1,wx)

  end function f21_c


  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function f21_dispatch_r(a, b, c, x, ts_tol) result(res)
    !! Checks if x is indeed in (0,1/2), and then makes a choice of evaluating the ODE (large a,b,c)
    !! or defaulting to the usual Taylor series

    use rotex__utils,     only: isin
    use rotex__globals,   only: ABC_THRESHOLD => HYPGEO_ABC_THRESHOLD

    implicit none (type, external)

    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: ts_tol
    real(dp) :: res

    logical :: bigabc
    real(dp) :: ts_tol_
      !! The tolerance for which \( \frac{ \left\lvert S_{N+1} - S_{N} \right\rvert }{ \left\lvert S_N \right\rvert }\)
      !! must be met for the series to be considered converged. If this is not supplied, this value will be taken
      !! to be machine epsilon `macheps_dp` from the `hypergeometric__constants` module.

    ! -- argument in range
    if(isin(x, 0._dp, 0.5_dp, lclosed=.true., rclosed=.true.) .eqv. .false.) then
      write(stderr, '("Re(x): ", e20.10)') x
      call die("Re(x) must be between 0 and 1/2 in F21_DISPATCH")
    endif

    ! -- size check on a, b, c
    bigabc = abs(a) .ge. ABC_THRESHOLD .OR. abs(b) .ge. ABC_THRESHOLD .OR. abs(c) .ge. ABC_THRESHOLD

    ! -- ODE if a, b, c too big
    if(bigabc) then
      write(stderr, '("WARN: Large value of a parameter detected ! The electron energy&
      & is probably very close to a threshold, resulting in very large η=-Z/k.")')
      write(stderr, '("      |A|: ", F7.3)') abs(a)
      write(stderr, '("      |B|: ", F7.3)') abs(b)
      write(stderr, '("      |C|: ", F7.3)') abs(c)
      ts_tol_ = TS_DEFAULT_TOLERANCE ; if(present(ts_tol)) ts_tol_ = ts_tol
      res = f21_ts_r(a, b, c, x, ts_tol_)
      return
    endif

    ! -- Taylor series otherwise
    ts_tol_ = TS_DEFAULT_TOLERANCE ; if(present(ts_tol)) ts_tol_ = ts_tol
    res = f21_ts_r(a, b, c, x, ts_tol_)

  end function f21_dispatch_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function f21_dispatch_c(za, zb, zc, z, ts_tol) result(res)
    !! Checks if x is indeed in (0,1/2), and then makes a choice of evaluating the ODE (large a,b,c)
    !! or defaulting to the usual Taylor series

    use rotex__utils,  only: isin
    use rotex__globals, only: ABC_THRESHOLD => HYPGEO_ABC_THRESHOLD
    use rotex__system, only: die, stderr

    implicit none (type, external)

    complex(dp), intent(in) :: za
    complex(dp), intent(in) :: zb
    complex(dp), intent(in) :: zc
    complex(dp), intent(in) :: z
    real(dp), intent(in), optional :: ts_tol
    complex(dp) :: res

    logical :: bigabc
    real(dp) :: x, ts_tol_
      !! The tolerance for which \( \frac{ \left\lvert S_{N+1} - S_{N} \right\rvert }{ \left\lvert S_N \right\rvert }\)
      !! must be met for the series to be considered converged. If this is not supplied, this value will be taken
      !! to be machine epsilon `macheps_dp` from the `hypergeometric__constants` module.

    ! -- argument in range
    if(z%im .ne. 0._dp) then
      write(stderr, '("Im(Z): ", e20.10)') z%im
      call die("Z is nonreal in F21_DISPATCH !")
    endif
    x = z%re
    if(isin(x, 0._dp, 0.5_dp, lclosed=.true., rclosed=.true.) .eqv. .false.) then
      write(stderr, '("Re(Z): ", e20.10)') x
      call die("Re(Z) must be between 0 and 1/2 in F21_DISPATCH")
    endif

    ! -- size check on a, b, c
    bigabc = abs(za) .ge. ABC_THRESHOLD .OR. abs(zb) .ge. ABC_THRESHOLD .OR. abs(zc) .ge. ABC_THRESHOLD

    ! -- ODE if a, b, c too big
    if(bigabc) then
      ! res = f21_ode_eval(za, zb, zc, x)
      write(stderr, '("WARN: Large value of a parameter detected ! The electron energy&
      & is probably very close to a threshold, resulting in very large η=-Z/k.")')
      write(stderr, '("      |A|: ", F7.3)') abs(za)
      write(stderr, '("      |B|: ", F7.3)') abs(zb)
      write(stderr, '("      |C|: ", F7.3)') abs(zc)
      ! res = f21_ode_eval(za, zb, zc, x)
      ts_tol_ = TS_DEFAULT_TOLERANCE ; if(present(ts_tol)) ts_tol_ = ts_tol
      res = f21_ts_c(za, zb, zc, z, ts_tol_)
      return
    endif

    ! -- Taylor series otherwise
    ts_tol_ = TS_DEFAULT_TOLERANCE ; if(present(ts_tol)) ts_tol_ = ts_tol
    res = f21_ts_c(za, zb, zc, z, ts_tol_)

  end function f21_dispatch_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function f21_ts_r(a, b, c, x, tol) result(res)
  ! pure elemental function f21_ts_r(a, b, c, x, tol) result(res)
    !! Returns the Gauss hypergeometric function ₂F₁(a,b,;c;z\) via a Taylor series method, with quad precision
    use rotex__utils,     only: downcast, upcast, kbn_sum, isint
    use rotex__system,    only: die, stderr
    use rotex__constants, only: macheps => macheps_dp, zero, one
    implicit none (type, external)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: tol
      !! The tolerance for which \( \frac{ \left\lvert S_{N+1} - S_{N} \right\rvert }{ \left\lvert S_N \right\rvert }\)
      !! must be met for the series to be considered converged. If this is not supplied, this value will be taken
      !! to be machine epsilon `macheps_dp` from the `hypergeometric__constants` module.
    real(dp) :: res
    integer, parameter :: kmax = 20000
    real(dp) :: tol_local
    real(qp) :: tol_local_qp
    integer :: k
    real(qp) :: kq
    real(qp) :: r
    real(qp) :: aa, bb, cc, xx
    real(qp) :: numer, denom
    real(qp) :: sumq, diff, comp
    tol_local = TS_DEFAULT_TOLERANCE ; if(present(tol)) tol_local = tol
    call upcast(tol_local, tol_local_qp)

    ! -- terminating seris -> calculate exactly
    if(isint(a) .AND. nint(a) .lt. 0) then
      res = f21_finite_r(a, b, c, x)
      return
    elseif(isint(b) .AND. nint(b) .lt. 0) then
      res = f21_finite_r(a, b, c, x)
      return
    endif

    ! -- upcast to quad precision
    call upcast(a, aa)
    call upcast(b, bb)
    call upcast(c, cc)
    call upcast(x, xx)

    sumq = 1
    comp = 0
    diff = 1
    k = 0

    do
      k = k + 1
      kq = real(k, kind = qp)
      numer = (aa  + kq - 1._qp ) * (bb + kq - 1._qp)
      denom = kq * (cc + kq - 1_qp )
      r = numer / denom
      diff = diff * r * xx
      call kbn_sum(sumq, comp, diff)
      ! write(6,*) numer, denom
      ! write(6, *) sumq, comp, diff
      ! write(6,*) ""
      if( abs(diff) .le. tol_local_qp * abs(sumq) ) exit
      if(k .lt. kmax) cycle
      ! -- non-convergence
      write(stderr, '("KMAX: ", I0)') kmax
      call die("k = kmax has been achieved without convergence in f21_ts_r")
    enddo

    ! -- downcast to double precision for return value
    call downcast(sumq + comp, res)

  end function f21_ts_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function f21_ts_c(a, b, c, z, tol) result(res)
    !! Returns the Gauss hypergeometric function ₂F₁(a,b,;c;z\) via a Taylor series method, with quad precision
    use rotex__utils,     only: downcast, upcast, kbn_sum, isint
    use rotex__system,    only: die
    use rotex__constants, only: macheps => macheps_dp, zero, one
    implicit none (type, external)
    complex(dp), intent(in) :: a
    complex(dp), intent(in) :: b
    complex(dp), intent(in) :: c
    complex(dp), intent(in) :: z
    real(dp), intent(in), optional :: tol
      !! The tolerance for which \( \frac{ \left\lvert S_{N+1} - S_{N} \right\rvert }{ \left\lvert S_N \right\rvert }\)
      !! must be met for the series to be considered converged. If this is not supplied, this value will be taken
      !! to be machine epsilon `macheps_dp` from the `hypergeometric__constants` module.
    complex(dp) :: res
    integer, parameter :: kmax = 20000
    real(dp) :: tol_local
    real(qp) :: tol_local_qp
    integer :: k
    real(qp) :: kq
    complex(qp) :: r
    complex(qp) :: aa, bb, cc, zz
    complex(qp) :: numer, denom
    complex(qp) :: sumq, diff, comp
    tol_local = TS_DEFAULT_TOLERANCE ; if(present(tol)) tol_local = tol
    call upcast(tol_local, tol_local_qp)

    ! -- terminating seris -> calculate exactly
    if(isint(a) .AND. nint(a%re) .lt. 0) then
      res = f21_finite_c(a, b, c, z)
      return
    elseif(isint(b) .AND. nint(b%re) .lt. 0) then
      res = f21_finite_c(a, b, c, z)
      return
    endif

    ! -- upcast to quad precision
    call upcast(a, aa)
    call upcast(b, bb)
    call upcast(c, cc)
    call upcast(z, zz)

    sumq = 1
    comp = 0
    diff = 1
    k = 0

    do
      k = k + 1
      kq = real(k, kind = qp)
      numer = (aa  + kq - 1._qp ) * (bb + kq - 1._qp)
      denom = kq * (cc + kq - 1_qp )
      r = numer / denom
      diff = diff * r * zz
      call kbn_sum(sumq, comp, diff)
      if( abs(diff) .le. tol_local_qp * abs(sumq) ) exit
      if(k .lt. kmax) cycle
      call die("k = kmax has been achieved without convergence in gauss_2f1_ts")
    enddo
    ! -- downcast to double precision for return value
    call downcast(sumq + comp, res)
  end function f21_ts_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function f21_finite_r(a,b,c,x) result(res)
    !! Calculate the finite sum of ₂F₁(a,b;c;x) when a or b is a negative integer because the
    !! rising factorial will eventually be 0
    use rotex__utils,  only: isint, kbn_sum
    use rotex__system, only: die, stderr
    implicit none (type, external)
    real(dp), intent(in) :: a,b,c,x
    real(dp) :: res, Sk, comp
    integer :: n, m, nc
    integer :: k
    n = 1
    m = 1
    if(isint(a)) n = nint(a)
    if(isint(b)) m = nint(b)
    if(n .lt. 0 .AND. m .lt. 0) then
      n = max(n,m)
    elseif(n.lt.0 .neqv. m.lt.0) then
      n = min(n,m)
    else
      call die("Finite 2F1 will not be finite because n and m are both positive")
    endif
    ! -- poch(c) might terminate before n terms, so guard against that
    ccheck: if(isint(c)) then
      nc = nint(c)
      if(nc .gt. 0) exit ccheck
      if(abs(n) .gt. abs(nc)) &
        call die("Finite ₂F₁ hits a pole of poch(c), because c is a negative integer close to 0 than a or b")
    endif ccheck
    Sk = 1
    res  = Sk
    comp = 0
    do k=1,abs(n)
      Sk = Sk * (a+k-1)*(b+k-1) / ((c+k-1)*real(k, kind=dp)) * x
      call kbn_sum(res, comp, Sk)
    end do
    res = res + comp
  end function f21_finite_r

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function f21_finite_c(a,b,c,z) result(res)
    !! Calculate the finite sum of ₂F₁(a,b;c;z) when a or b is a negative integer because the
    !! rising factorial will eventually be 0
    use rotex__utils,  only: isint, kbn_sum
    use rotex__system, only: die, stderr
    implicit none (type, external)
    complex(dp), intent(in) :: a,b,c,z
    complex(dp) :: res, Sk, comp
    integer :: n, m, nc
    integer :: k
    n = 1
    m = 1
    if(isint(a)) n = nint(a%re)
    if(isint(b)) m = nint(b%re)
    if(n .lt. 0 .AND. m .lt. 0) then
      n = max(n,m)
    elseif(n.lt.0 .neqv. m.lt.0) then
      n = min(n,m)
    else
      call die("Finite 2F1 will not be finite because n and m are both positive")
    endif
    ! -- poch(c) might terminate before n terms, so guard against that
    ccheck: if(isint(c)) then
      nc = nint(c%re)
      if(nc .gt. 0) exit ccheck
      if(abs(n) .gt. abs(nc)) &
        call die("Finite ₂F₁ hits a pole of poch(c), because c is a negative integer close to 0 than a or b")
    endif ccheck
    Sk = 1
    res  = Sk
    comp = 0
    do k=1,abs(n)
      Sk = Sk * (a+k-1)*(b+k-1) / ((c+k-1)*real(k, kind=dp)) * z
      call kbn_sum(res, comp, Sk)
    end do
    res = res + comp
  end function f21_finite_c

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function olver_2f1_dlmf_15_8_10(a, b, m, x, tol) result(res)
    !! Returns Olver's hypergeometric function F(a,b;a+b+m;x)
    !! for integral a+b+m,m > 0,0

    use rotex__constants, only: macheps => macheps_dp
    use rotex__polygamma, only: lgamma => log_gamma, digamma
    use rotex__functions, only: logp1, factorial
    use rotex__system,    only: stderr

    implicit none (type, external)

    real(dp), intent(in) :: a, b
    integer,  intent(in) :: m
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: tol
    real(dp) :: res

    integer, parameter :: KMAX = 20000

    logical :: continueloop
    integer :: k
    real(dp) :: tol_
    real(dp) :: Sk, rk, termk, termkm1, sumk, sumkm1

    tol_ = macheps ; if(present(tol)) tol_ = tol
    res = 0._dp

    if(m .gt. 0) then
      k = 0
      Sk = factorial(m-k-1)
      res = res + Sk

      do k=1, m-1
        Sk = Sk * (a+(k-1)) * (b+(k-1)) / (m-k) / k * (x-1)
        res = res + Sk
      enddo

      ! -- multiplicative prefactor
      if(res .ne. 0) res = res * exp(-lgamma(a+m) - lgamma(b+m))

    endif

    ! -- initial terms. The sum is calculated as Σk Sk*rk where Sk and rk
    !    are defined recursively
    k = 0
    rk = logp1(-x) - digamma(k+1) - digamma(k+m+1) + digamma(a+k+m) + digamma(b+k+m)
    Sk = 1._dp / factorial(k+m)
    termk   = Sk*rk
    termkm1 = termk
    sumk    = termk
    sumkm1  = termk

    continueloop = .true.
    do while(continueloop)
      k = k+1

      ! -- digamma recurrence
      rk = rk                                 &
         - 1._dp/real(k,         kind=dp)     &
         - 1._dp/real(k+m,       kind=dp)     &
         + 1._dp/( a+real((k-1)+m, kind=dp) ) &
         + 1._dp/( b+real((k-1)+m, kind=dp) )
      Sk = ( a+real(m+(k-1), kind=dp) ) &
         * ( b+real(m+(k-1), kind=dp) ) &
         / real(k*(k+m), kind=dp)       &
         * (1._dp - x) * Sk

      termk  = Sk*rk
      sumkm1 = sumk
      sumk   = sumk + termk

      continueloop = k .lt. 2 &
        .OR. (abs(termk)   .gt. tol_*max(abs(sumk),   1._dp)) &
        .OR. (abs(termkm1) .gt. tol_*max(abs(sumkm1), 1._dp))

      termkm1 = termk

      ! -- loop bound
      if(k .lt. KMAX) cycle

      write(stderr, '("Max loop interation of ", I0, " achieved.")') KMAX
      error stop

    enddo

    res = res - sumk * (x-1._dp)**m * exp(-lgamma(a) - lgamma(b))

  end function olver_2f1_dlmf_15_8_10

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function olver_2f1_dlmf_15_8_12(a, b, m, x, tol) result(res)
    !! Returns Olver's hypergeometric function F(a,b;a+b-m;x)
    !! for integral a+b-m<0 ,m > 0

    use rotex__functions, only: logp1
    use rotex__constants, only: macheps => macheps_dp

    implicit none (type, external)

    real(dp), intent(in) :: a, b
    integer,  intent(in) :: m
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: tol
    real(dp) :: res

    integer :: mneg
    real(dp) :: tol_
    real(dp) :: aa, bb, mm

    tol_ = macheps ; if(present(tol)) tol_ = tol

    mm = real(m, kind=dp)
    aa = a - mm
    bb = b - mm

    res = exp( -mm*logp1(-x) ) * olver_2f1_dlmf_15_8_10(aa, bb, m, x, tol_)

  end function olver_2f1_dlmf_15_8_12

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function appell_f2_xy_ac1(a, b1, b2, c1, x) result(res)
    !! Calculate the Appell hypergeometric function for equal arguments x and y, and
    !! equal parameters a and c₂
    !!
    !!   F₂(a,b1,b2;c₁,c₂;x,y) = Σ_{m,n=0}^∞ (a)_{m+n} (b₁)_m (b₂)_n / ( (c₁)_m (c₁)_n m! n!) × x^m y^n
    !!
    !! with
    !!
    !!   a  = -λ + 1 + iξ
    !!   b1 =  l + λ + 1 - iη_i
    !!   b2 =  l + 1 + iη_f
    !!   c1 =  λ + 1 + iξ
    !!   c2 = -λ + 1 + iξ
    !!
    !! The equal arguments collapse and parameters collapse to a single sum over n
    !! modulated by the Gausse hypergeometric function
    !!
    !!   F₂ = Σ_n (b₂)_n / n! * x^n * ₂F₁(a+n, b₁; c₁; x)

    implicit none(type, external)

    complex(dp), intent(in) :: a, b1, b2, c1
    real(dp), intent(in) :: x
    complex(dp) :: res

    integer,  parameter :: NCONSEC = 3          !! consecutive small terms to declare convergence
    integer,  parameter :: NMAX    = 1000       !! max number of loop iterations
    real(qp), parameter :: EPS_REL = epsilon(1.0_dp) / 100.0_dp
    real(qp), parameter :: MAXLOST = 5.0_qp
      !! max number of digits that we will tolerate losing after cancellation

    real(qp) :: sumabs, lost, ratio
    complex(qp) :: asum, term, pochb2, xx
    complex(dp) :: g
    integer :: n, nsmall

    if(abs(x) .ge. 0.5_dp) call die("Appell F2 is only defined for |x|+|y=x|<1 => 2|x| < 1")

    asum   = 0
    sumabs = 0
    pochb2 = 1
    xx   = 1
    nsmall = 0

    do n = 0, NMAX

      g      = F21(a + n, b1, c1, x)
      term   = pochb2 * xx * cmplx(g, kind=qp)
      asum   = asum   + term
      sumabs = sumabs + abs(term)

      if(abs(term) .lt. EPS_REL*max(abs(asum), tiny(1.0_qp))) then
        nsmall = nsmall + 1
        if(nsmall .ge. NCONSEC) exit
      else
        nsmall = 0
      endif

      pochb2 = pochb2 * (cmplx(b2, kind=qp) + n) / cmplx(n+1, kind=qp) ! (b₂)_{n+1}/(n+1)!
      xx = xx * cmplx(x, kind=qp)

      if(n .eq. NMAX) call die("Appel F2 did not converge within NMAX terms. Variable x is&
        & is probably approaching 1/2 (channel threshold)")

    enddo

    res = cmplx(asum, kind=dp)

  end function appell_f2_xy_ac1



! ================================================================================================================================ !
end module rotex__hypergeometric
! ================================================================================================================================ !
