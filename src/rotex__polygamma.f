! ================================================================================================================================ !
module rotex__polygamma
  !! Routines for calculating Γ(z) and ln(Γ(z)).
  use rotex__kinds, only: dp, qp
! #ifdef WITH_STDLIB
!   use stdlib_specialfunctions_gamma, only: gamma, log_gamma
! #endif

  implicit none

  private

  public :: gamma
  public :: log_gamma
  public :: digamma

  real(dp), parameter :: tol_dp = epsilon(1.0_dp)

#ifndef WITH_STDLIB

  interface gamma
    module procedure :: gamma_cdp
  end interface gamma

  interface log_gamma
    module procedure :: l_gamma_rdp
    module procedure :: l_gamma_cdp
  end interface log_gamma

  interface digamma
    module procedure :: digamma_i
    module procedure :: digamma_rdp
    ! module procedure :: digamma_cdp
  end interface digamma

  interface horner
    module procedure :: horner_rdp
    module procedure :: horner_rqp
    module procedure :: horner_cdp
    module procedure :: horner_cqp
  end interface horner

  real(qp), parameter :: DIGAM_ASYMP_LIMIT_DP = 10._qp
    !! The routines to calculate the digamma function will use recursion
    !! to reach values of the argument that are larger than this value.
    !! Then, the expansion over the Bernoulli numbers will be used in the
    !! asymptotic expansion to determine the value of the digamma ψ(x)
    !! at x > digam_asympt_limit

  integer, parameter :: N_DIGAM_XPANSION = 22
    !! The number of terms to include in the asymptotic expansion of the
    !! digamma function over the Bernoulli numbers

  real(qp), parameter :: OEIS_a001067(29) = &
    !! Integer sequence A001067 from the [OEIS](https://oeis.org/) :
    !! numerator of Bernoulli(2n)/(2n)
                                          [ 1._qp                                  &
                                          ,-1._qp                                  &
                                          , 1._qp                                  &
                                          ,-1._qp                                  &
                                          , 1._qp                                  &
                                          ,-691._qp                                &
                                          , 1._qp                                  &
                                          ,-3617._qp                               &
                                          , 43867._qp                              &
                                          ,-174611._qp                             &
                                          , 77683._qp                              &
                                          ,-236364091._qp                          &
                                          , 657931._qp                             &
                                          ,-3392780147._qp                         &
                                          , 1723168255201._qp                      &
                                          ,-7709321041217._qp                      &
                                          , 151628697551._qp                       &
                                          ,-26315271553053477373._qp               &
                                          , 154210205991661._qp                    &
                                          ,-261082718496449122051._qp              &
                                          , 1520097643918070802691._qp             &
                                          ,-2530297234481911294093._qp             &
                                          , 25932657025822267968607._qp            &
                                          ,-5609403368997817686249127547._qp       &
                                          , 19802288209643185928499101._qp         &
                                          ,-61628132164268458257532691681._qp      &
                                          , 29149963634884862421418123812691._qp   &
                                          ,-354198989901889536240773677094747._qp  &
                                          , 2913228046513104891794716413587449._qp &
                                          ]

  real(qp), parameter :: OEIS_a006953(36) = &
    !! Integer sequence A006953 from the [OEIS](https://oeis.org/) : denominator
    !! of Bernoulli(2n)/(2n)
                                          [ 12._qp          &
                                          , 120._qp         &
                                          , 252._qp         &
                                          , 240._qp         &
                                          , 132._qp         &
                                          , 32760._qp       &
                                          , 12._qp          &
                                          , 8160._qp        &
                                          , 14364._qp       &
                                          , 6600._qp        &
                                          , 276._qp         &
                                          , 65520._qp       &
                                          , 12._qp          &
                                          , 3480._qp        &
                                          , 85932._qp       &
                                          , 16320._qp       &
                                          , 12._qp          &
                                          , 69090840._qp    &
                                          , 12._qp          &
                                          , 541200._qp      &
                                          , 75852._qp       &
                                          , 2760._qp        &
                                          , 564._qp         &
                                          , 2227680._qp     &
                                          , 132._qp         &
                                          , 6360._qp        &
                                          , 43092._qp       &
                                          , 6960._qp        &
                                          , 708._qp         &
                                          , 3407203800._qp  &
                                          , 12._qp          &
                                          , 32640._qp       &
                                          , 388332._qp      &
                                          , 120._qp         &
                                          , 9372._qp        &
                                          , 10087262640._qp &
                                          ]


  real(qp), parameter :: digam_xpansion(n_digam_xpansion) = OEIS_a001067(1:n_digam_xpansion) &
                                                          / OEIS_a006953(1:n_digam_xpansion)
    !! The expansion of the digamma function for large values of its argument

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function gamma_cdp(z) result(res)
    !! MIT License
    !!
    !! Copyright (c) 2019-2021 stdlib contributors
    !!
    !! Permission is hereby granted, free of charge, to any person obtaining a copy
    !! of this software and associated documentation files (the "Software"), to deal
    !! in the Software without restriction, including without limitation the rights
    !! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    !! copies of the Software, and to permit persons to whom the Software is
    !! furnished to do so, subject to the following conditions:
    !!
    !! The above copyright notice and this permission notice shall be included in all
    !! copies or substantial portions of the Software.
    !!
    !! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    !! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    !! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    !! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    !! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    !! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    !! SOFTWARE.
    use rotex__utils, only: downcast
    implicit none
    complex(dp), intent(in) :: z
    complex(dp) :: res
    integer :: i

    real(dp), parameter :: zero_k1 = 0.0_dp
    real(qp), parameter :: half = 0.5_qp,             &
                         one = 1.0_qp, pi = acos(- one), sqpi = sqrt(pi)
    complex(qp) :: y, x, sum


    integer, parameter :: n = 24
    real(qp), parameter :: r = 25.617904_qp
    real(qp), parameter :: d(0 : n)=                                         &
                     [1.0087261714899910504854136977047144166e-11_qp,  &
                          1.6339627701280724777912729825256860624_qp,  &
                      -1.4205787702221583745972794018472259342e+1_qp,  &
                       5.6689501646428786119793943350900908698e+1_qp,  &
                      -1.3766376824252176069406853670529834070e+2_qp,  &
                       2.2739972766608392140035874845640820558e+2_qp,  &
                      -2.7058382145757164380300118233258834430e+2_qp,  &
                      2.39614374587263042692333711131832094166e+2_qp,  &
                      -1.6090450559507517723393498276315290189e+2_qp,  &
                      8.27378183187161305711485619113605553100e+1_qp,  &
                      -3.2678977082742592701862249152153110206e+1_qp,  &
                         9.89018079175824824537131521501652931756_qp,  &
                         -2.2762136356329318377213053650799013041_qp,  &
                      3.93265017303573867227590563182750070164e-1_qp,  &
                      -5.0051054352146209116457193223422284239e-2_qp,  &
                      4.57142601898244576789629257292603538238e-3_qp,  &
                      -2.8922592124650765614787233510990416584e-4_qp,  &
                      1.20833375377219592849746118012697473202e-5_qp,  &
                      -3.1220812187551248389268359432609135033e-7_qp,  &
                      4.55117045361638520378367871355819524460e-9_qp,  &
                     -3.2757632817493581828033170342853173968e-11_qp,  &
                     9.49784279240135747819870224486376897253e-14_qp,  &
                     -7.9480594917454410117072562195702526836e-17_qp,  &
                     1.04692819439870077791406760109955648941e-20_qp,  &
                     -5.8990280044857540075384586350723191533e-26_qp]
    ! parameters from above referenced source.


    if(abs(z % im) < tol_dp) then

        res = cmplx(gamma(z % re), kind = dp)
        return

    end if

    if(z % re < zero_k1) then

        x = cmplx(abs(z % re), - z % im, kind = dp)
        y = x - one

    else

        y = z - one

    end if

    sum = cmplx(d(0), kind = qp)

    do i = 1, n

        sum = sum + d(i) / (y + i)

    end do

    y = exp((y + half) * log(y + half + r) - y) * sum

    y = y * 2 / sqpi                         !Re(z) > 0 return

    if(z % re < zero_k1 ) then

        y = - pi / (sin(pi * x) * x * y)     !Re(z) < 0 return

    end if

    ! -- y -> res
    call downcast(y, res)

  end function gamma_cdp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function l_gamma_rdp(x) result(res)
    !! Computes ln(Γ(x)) for real x
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: res
    intrinsic :: log_gamma
    res = log_gamma(x)
  end function l_gamma_rdp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function l_gamma_cdp(z) result (res)
  !
  ! log_gamma function for any complex number, excluding negative whole number
  ! "Computation of special functions", Shanjie Zhang & Jianmin Jin, 1996, p.48
  ! "Computing the principal branch of log-gamma", D.E.G. Hare,
  ! J. of Algorithms, 25(2), 1997 p. 221–236
  !
  ! Fortran 90 program by Jim-215-Fisher
  !
      complex(dp), intent(in) :: z
      complex(dp) :: res, z1, z2
      real(dp) :: d
      integer :: m, i
      complex(qp) :: zr, zr2, sum, s
      real(dp), parameter :: z_limit = 10.0_dp, zero_k1 = 0.0_dp
      integer, parameter :: n = 20
      real(qp), parameter :: zero = 0.0_qp, one = 1.0_qp,              &
                           pi = acos(-one), ln2pi = log(2 * pi)
      real(qp), parameter :: a(n) = [                                          &
                         .8333333333333333333333333333333333333333E-1_qp,&
                        -.2777777777777777777777777777777777777778E-2_qp,&
                         .7936507936507936507936507936507936507937E-3_qp,&
                        -.5952380952380952380952380952380952380952E-3_qp,&
                         .8417508417508417508417508417508417508418E-3_qp,&
                        -.1917526917526917526917526917526917526918E-2_qp,&
                         .6410256410256410256410256410256410256410E-2_qp,&
                        -.2955065359477124183006535947712418300654E-1_qp,&
                         .1796443723688305731649384900158893966944E+0_qp,&
                        -.1392432216905901116427432216905901116427E+1_qp,&
                         .1340286404416839199447895100069013112491E+2_qp,&
                        -.1568482846260020173063651324520889738281E+3_qp,&
                         .2193103333333333333333333333333333333333E+4_qp,&
                        -.3610877125372498935717326521924223073648E+5_qp,&
                         .6914722688513130671083952507756734675533E+6_qp,&
                        -.1523822153940741619228336495888678051866E+8_qp,&
                         .3829007513914141414141414141414141414141E+9_qp,&
                       -.1088226603578439108901514916552510537473E+11_qp,&
                        .3473202837650022522522522522522522522523E+12_qp,&
                       -.1236960214226927445425171034927132488108E+14_qp]
      ! parameters from above reference

      z2 = z

      if(z % re < zero_k1) then

          z2 = cmplx(abs(z % re), - z % im, kind = dp) + 1

      end if

      d = hypot(z2 % re, z2 % im)
      z1 = z2
      m = 0

      if(d <= z_limit) then                       !for small |z|

          m = ceiling(z_limit - d)
          z1 = z2 + m

      end if

      zr = one / z1
      zr2 = zr * zr

      sum = (((a(20) * zr2 + a(19)) * zr2 + a(18)) * zr2 + a(17)) * zr2
      sum = (((sum + a(16)) * zr2 + a(15)) * zr2 + a(14)) * zr2
      sum = (((sum + a(13)) * zr2 + a(12)) * zr2 + a(11)) * zr2
      sum = (((sum + a(10)) * zr2 + a(9)) * zr2 + a(8)) * zr2
      sum = (((sum + a(7)) * zr2 + a(6)) * zr2 + a(5)) * zr2
      sum = (((sum + a(4)) * zr2 + a(3)) * zr2 + a(2)) * zr2
      sum = (sum + a(1)) * zr + ln2pi / 2 - z1 + (z1 - 0.5_qp) * log(z1)

      if(m /= 0) then

          s = cmplx(zero, zero, kind = qp)

          do i = 1, m

              s = s + log(cmplx(z1, kind = qp) - i)

          end do

          sum = sum - s

      end if

      if(z % re < zero_k1) then

          sum = log(pi) - log(sin(pi * z)) - sum
          m = ceiling((2 * z % re - 3) / 4)
          sum % im = sum % im + 2 * pi * m * sign(1.0_dp, z % im)

      end if

      res = cmplx(sum, kind = dp)
  end function l_gamma_cdp

#endif

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function digamma_i(n) result(res)
    !! Returns the digamma function ψ(n) using a truncated Stirling / de Moivre series
    !! for integral n
    implicit none
    integer, intent(in) :: n
    real(dp) :: res
    res = digamma_rdp(real(n, kind=dp))
  end function digamma_i

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function digamma_rdp(x) result(res)
    !! Returns the digamma function ψ(x) using a truncated Stirling / de Moivre series
    !! for real x

    use rotex__utils,     only: downcast, upcast
    use rotex__constants, only: pi => pi_qp, euler_mascheroni
    use rotex__functions, only: cotpi, isinteger, iseven
    use rotex__system,    only: die

    implicit none

    real(dp), intent(in) :: x
    real(dp) :: res
    real(qp) :: resqp

    real(qp), parameter :: twolog2 = 2._qp*log(2._qp)

    integer :: k, m, n
    real(qp) :: x2, xr, xr2

    if(isinteger(x) .AND. nint(x) .le. 0) call die("Digamma ψ(x) not defined for non-positive integers")

    ! -- x = (2n+1)/2
    halfint: if(isinteger(2*x)) then
      m = nint(2*x)
      if(m .lt. 1 .OR. iseven(m)) exit halfint
      n = (m-1)/2
      resqp = -euler_mascheroni - twolog2
      do k=1, n
        resqp = resqp + 2._qp / real(2*k-1, kind=qp)
      enddo
      call downcast(resqp, res)
      return
    endif halfint

    call upcast(x, x2)

    ! -- reflection ψ(1-x) - ψ(x) = π cot(πx) when x < 1/2 to stay away from 0
    resqp = 0._qp
    if(x2 .lt. 0.5_qp) then
      resqp = -pi*cotpi(x2)
      x2 = 1._qp - x2
    endif

    ! -- forward recurrence x -> x+1
    do while(abs(x2) .le. DIGAM_ASYMP_LIMIT_DP)
      resqp = resqp - 1._qp/x2
      x2 = x2 + 1._qp
    enddo

    xr = 1._qp/x2
    xr2 = xr*xr

    ! -- truncated series
    resqp = resqp + log(x2) - xr/2._qp - horner(digam_xpansion, xr2)


    call downcast(resqp, res)

  end function digamma_rdp

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function horner_rdp(coeffs, x) result(res)
    !! Evaluate S(x) = c₁x + c₂x² + ... + cₙxⁿ using Horner's rule :
    !! S(x) = y*(c₁ + y*( c₂ + y*( c₃ + ... ) ))
    implicit none
    real(dp), intent(in) :: coeffs(:)
    real(dp), intent(in) :: x
    real(dp) :: res
    integer :: k, n
    n = size(coeffs, 1)
    res = coeffs(n)
    do k = n-1, 1, -1
      res = coeffs(k) + x*res
    enddo
    res = res * x
  end function horner_rdp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function horner_cdp(coeffs, z) result(res)
    !! Evaluate S(z) = c₁z + c₂z² + ... + cₙzⁿ using Horner's rule :
    !! S(z) = y*(c₁ + y*( c₂ + y*( c₃ + ... ) ))
    implicit none
    real(dp), intent(in) :: coeffs(:)
    complex(dp), intent(in) :: z
    complex(dp) :: res
    integer :: k, n
    n = size(coeffs, 1)
    res = coeffs(n)
    do k = n-1, 1, -1
      res = coeffs(k) + z*res
    enddo
    res = res * z
  end function horner_cdp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function horner_rqp(coeffs, x) result(res)
    !! Evaluate S(x) = c₁x + c₂x² + ... + cₙxⁿ using Horner's rule :
    !! S(x) = y*(c₁ + y*( c₂ + y*( c₃ + ... ) ))
    implicit none
    real(qp), intent(in) :: coeffs(:)
    real(qp), intent(in) :: x
    real(qp) :: res
    integer :: k, n
    n = size(coeffs, 1)
    res = coeffs(n)
    do k = n-1, 1, -1
      res = coeffs(k) + x*res
    enddo
    res = res * x
  end function horner_rqp
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function horner_cqp(coeffs, z) result(res)
    !! Evaluate S(z) = c₁z + c₂z² + ... + cₙzⁿ using Horner's rule :
    !! S(z) = y*(c₁ + y*( c₂ + y*( c₃ + ... ) ))
    implicit none
    real(qp), intent(in) :: coeffs(:)
    complex(qp), intent(in) :: z
    complex(qp) :: res
    integer :: k, n
    n = size(coeffs, 1)
    res = coeffs(n)
    do k = n-1, 1, -1
      res = coeffs(k) + z*res
    enddo
    res = res * z
  end function horner_cqp

! ================================================================================================================================ !
end module rotex__polygamma
! ================================================================================================================================ !
