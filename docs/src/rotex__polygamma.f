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

  real(dp), parameter :: tol_dp = epsilon(1.0_dp)

#ifndef WITH_STDLIB

  interface gamma
    module procedure :: gamma_cdp
  end interface gamma

  interface log_gamma
    module procedure :: l_gamma_cdp
  end interface log_gamma

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

! ================================================================================================================================ !
end module rotex__polygamma
! ================================================================================================================================ !
