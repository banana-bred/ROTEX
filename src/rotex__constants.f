! ================================================================================================================================ !
module rotex__constants
  !! Contains most (but not all) defined constants used throughout the code.

  use rotex__kinds,    only: dp, qp

  implicit none (type, external)

  private

  real(dp), parameter, public :: ahalf = 0.5_dp
  real(dp), parameter, public :: zero  = 0._dp
  real(dp), parameter, public :: one   = 1._dp
  real(dp), parameter, public :: two   = 2._dp
  real(dp), parameter, public :: three = 3._dp
  real(dp), parameter, public :: four  = 4._dp
  real(dp), parameter, public :: pi    = atan(1._dp) * 4._dp
    !! π
  real(dp), parameter, public :: invc = 1._dp/137.035999177_dp
    !! 1/c, where c is the speed of light in atomic units. This is equal to the fine structure constant α.
  real(dp), parameter, public :: macheps_dp = epsilon(one)
  real(dp), parameter, public :: rho = 0.8_dp
  real(dp), parameter, public :: au2ryd = 2._dp
    !! Hartree -> Rydberg
  real(dp), parameter, public :: au2ev = 27.2113834e0_dp
    !! Hartree -> eV
  real(dp), parameter, public :: au2ang = 0.5291772083_dp
    !! Bohr -> Ångstrom
  real(dp), parameter, public :: au2cm = au2ang * 1e-8_dp
    !! Bohr -> cm
  real(dp), parameter, public :: au2invcm = 219474.6313710e0_dp
    !! multiplication factor to convert atomic units of energy (hartree) to wavenumbers (inverse centimeters)
  real(dp), parameter, public :: au2sec = 2.4188843e-17_dp
    !! multiplication factor to convert atomic units of time (hbar / hartree) to seconds
  real(dp), parameter, public :: au2deb = 2.5417462762781615_dp
    !! atomic units (electric dipole) -> Debye

  real(qp), parameter, public :: pi_qp = atan(1._qp) * 4._qp
    !! π
  real(qp), parameter, public :: euler_mascheroni = 0.57721566490153286060651209008240243104215933593992_qp
    !! the Euler Mascheroni constant. The difference between the harmonic series and the natural logarithm

  complex(dp), parameter, public :: im = (zero, one)
    !! the square root of -1

! ================================================================================================================================ !
end module rotex__constants
! ================================================================================================================================ !
