! ================================================================================================================================ !
module rotex__wigner
  !! Calculate the Wigner 3j symbols

  implicit none

  private

  public :: wigner3j
  public :: clebsch

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function clebsch(j1, m1, j2, m2, j, m) result(res)
    !! Returns the Clebsch-Gordan coefficient
    !! \(C^{jm}_{j_1m_1,j_2m_2} \equiv \left\langle j_1m_1,j_2_m2|jm \right\rangle\)
    !! by using its relation to the Wigner 3j symbol
    use rotex__types,     only: dp
    use rotex__functions, only: neg
    implicit none
    integer, intent(in) :: j1, m1, j2, m2, j, m
    real(dp) :: res
    if(m .ne. m1+m2) then
      res = 0
      return
    endif
    res = neg(-j1+j2-m) * sqrt(real(2*j+1, kind = dp)) * wigner3j(2*j1, 2*j2, 2*j, 2*m1, 2*m2, -2*m)
  end function clebsch

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function wigner3j(dj1, dj2, dj3, dm1, dm2, dm3) result(res)
    !! returns the Wigner 3j symbol via explicit calculation. These should probably be
    !! precomputed, but this works for now. The [WignerSymbol-f](https://github.com/0382/WignerSymbol-f)
    !! repo seems like a good implementation. The [wigner](https://github.com/ogorton/wigner)
    !! repo by ogorton takes up way too much memory for the high partial waves because
    !! it tries to allocate (2N+1)^6 doubles. NOTE: for the high partial wave, we know that
    !! we'll only need values with m1 = m2 = m3 = 0, so take advantage of this ?

    use rotex__types,     only: dp
    use rotex__system,    only: die
    use rotex__constants, only: two
    use rotex__functions, only: istriangle, iseven, logfac => log_factorial

    implicit none

    integer, intent(in) :: dj1, dj2, dj3
      !! twice the angular momenta j
    integer, intent(in) :: dm1, dm2, dm3
      !! twice the angular momenta m
    real(dp) :: res
      !! the result

    integer :: dk, kmin, kmax
    real(dp) :: phase, summation, logs

    res = 0

    ! -- check basic properties
    if(dj1 .lt. 0) call die("j1 cannot be less than 0 in the Wigner 3j symbol")
    if(dj2 .lt. 0) call die("j2 cannot be less than 0 in the Wigner 3j symbol")
    if(dj3 .lt. 0) call die("j3 cannot be less than 0 in the Wigner 3j symbol")
    if(.not. istriangle(dj1, dj2, dj3)) return ! -- triangle inequality
    if(dm1 + dm2 + dm3 .ne. 0)          return ! -- m1 + m2 + m3 = 0
    if(.not. iseven(dj1 + dj2 + dj3))   return ! -- j1 + j2 + j3 is integer
    if(abs(dm1) .gt. abs(dj1))          return
    if(abs(dm2) .gt. abs(dj2))          return
    if(abs(dm3) .gt. abs(dj3))          return

    phase = (-1.0_dp)**( (dj1 - dj2 - dm3)/two )

    ! -- sqrt of factorials calculated as a logarithm
    logs = (                                                          &
               logfac((dj1+dj2-dj3)/two)  + logfac((dj1-dj2+dj3)/two)     &
             + logfac((-dj1+dj2+dj3)/two) - logfac((dj1+dj2+dj3)/two + 1) &
             + logfac((dj1-dm1)/two)      + logfac((dj1+dm1)/two)         &
             + logfac((dj2-dm2)/two)      + logfac((dj2+dm2)/two)         &
             + logfac((dj3-dm3)/two)      + logfac((dj3+dm3)/two)         &
           ) / two

    ! -- the summation variable z runs over the integers for which
    !    all factorial arguments are non negative
    kmin = max( 0, max( (dj2-dj3-dm1), (dj1-dj3+dm2) ) )
    kmax = min( (dj1+dj2-dj3), min( (dj1-dm1), (dj2+dm2) ) )

    summation = 0
    ! -- this summation could easily be done recursively to avoid recalulating
    !    logfac for each k
    do dk = kmin, kmax, 2
      summation = summation + (-1)**(dk/two) &
                / exp( logfac(dk/two)           + logfac((dj1+dj2-dj3-dk)/two) + logfac((dj1-dm1-dk)/two)       &
                     + logfac((dj2+dm2-dk)/two) + logfac((dj3-dj2+dm1+dk)/two) + logfac((dj3-dj1-dm2+dk)/two) )
    end do

    res = phase * exp(logs) * summation

  end function wigner3j

! ================================================================================================================================ !
end module rotex__wigner
! ================================================================================================================================ !
