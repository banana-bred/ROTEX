! ================================================================================================================================ !
module rotex__globals
  !! Putting some important global variables in here so we don't have to pass them everywhere
  !! all the time
  use rotex__kinds, only: dp
  use rotex__types, only: cd4_type, cd6_type

  implicit none (type, external)

  private

  save

  ! public :: setglobals
  public :: config_type

  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓ namelist variables declarations ↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  type config_type
    !! Derived type containing data from the namelist variables

    logical :: add_cd4 = .false.
      !! Add centrifugal distortion for fourth order ?
    logical :: add_cd6 = .false.
      !! Add centrifugal distortion for sixth order ?
    logical :: only_einsta
      !! Whether to only calculate the Einstein A coefficients in the Coulomb-Born
      !! cross section routine
    logical :: use_CDMS_einstA
      !! Whether to use Einstein A coefficients obtained from the CDMS or calculate them ourselves
      !! Path for the file containing the CDMS data to be read (if use_CDMS_einstA is .true.)
    logical :: analytic_total_cb(2)
      !! Array of logicals that has the size 2. Choose whether to use the analytic equation
      !! describing the multipole expansions for the dipole (element 1) and the quadrupole (element 2, not yet available)
    logical :: do_xtrap
      !! Do the extrapolation of (de-)excitation cross sections as 1/E to the excitation threshold ?
    logical :: do_dipole
      !! Choose whether to use the dipole term of the potential expansion
    logical :: do_quadrupole
      !! Choose whether to use the quadrupole term of the potential expansion
    logical :: use_kmat
      !! Calculate (de-)excitation cross sections using precomputed K-matrices ?
    logical :: use_cb
      !! Calculate (de-)excitation cross sections using the Coulomb-Born approxiation ?
    logical :: real_spherical_harmonics
      !! Whether the input K-matrices are evaluated in a basis of real spherical harmonics
      !! for the scattering electron. If .true., it will be transformed to a basis of
      !! complex-valued spherical harmonics
    logical :: print_elec_channels
      !! Whether to print electronic channels once they are defined and ordered

    integer :: spin_isomer_kind
      !! Whether and how to enforce ortho/para symmetry for molecules with identical nuclei.
      !!   0: don't
      !!   1: Dsh linear rotor; basically, homonuclear diatomics
      !!   2: C2v rotor (H₂X-like): preserve Ka+Kc parity
      !! Note that this just disables certain transitions from bein calculated
      !! in the CB approx as well as from the S-matrix. This does not affect
      !! the RFT because higher J-blocks of the S-matrix are more affected
      !! by K-mixing (Ka and Kc are not exact quantum numbers)
    integer :: forbidden_states_kind
      !! Some molecules have rotatinonal levels that are forbidden in certain vibrational states, e.g., H₃⁺
      !! In this case, do not construct certain levels:
      !!   0: no restriction
      !!   1: even N, K=0 states cannot exist symmetric top only
    logical :: symtop_reduce_projection
      !! The code internally calculates states for symmetric tops with ±K. If this is
      !! .true., ruce all of these to |K|. For example, the four transitions ±1 -> ±3 will
      !! be reduced to the single transition 1 -> 3
    logical :: edft
      !! Use an energy dependent frame transformation for the K/S-matrices ? This will read all available
      !! energies provided between kmat_Ei and kmat_Ef
    integer :: nE
      !! The number of scattering energies to consider. This does not need to be very high; the CB
      !! cross sections are very smooth and can easily be interpolated.
    integer :: nE_xtrap
      !! Number of extrapolation energies. Excitation cross sections are extrapolated as 1/E to the excitation threshold,
      !! de-excitation cross sections are extrapolated as 1/E to Ei_xtrap. If this is 0, no exptrapolation will be performed.
    integer :: lmax_partial
      !! The maximum value of l to consider in the contribution of the partial CB cross section
      !! from the dipole and the quadrupole. If you're replacing the low-l CB cross sections
      !! with other cross sections, set this to the max l that you have available.
    integer :: lmax_total
      !! The maximum value of l to consider in the contribution of the total CB cross section
      !! in the even that you're not using the analytic expression,
      !! from the dipole  and the quadrupole
    integer :: Nmin
      !! The minimum value of the rotational quantum number (N) to consider
    integer :: Nmax
      !! The maximum value of the rotational quantum number (N) to consider
    integer :: targcharge
      !! The electric charge of the target
    integer :: lmax_kmat
      !! The max partial wave to be included in the K-matrix basis. Cannot exceed the available
      !! basis in the calculation, but can be smaller than the largest available partial wave
    integer :: num_egrid_segs
      !! Number of energy grid segments (evaluation energy for the cross sections)
    integer, allocatable :: num_egrid(:)
      !! Array of number of energies per grid segment (length num_egrid_segs)
    integer, allocatable :: spinmults(:)
      !! Array of spin multiplicities (2S+1) for which the system's (target + e⁻) K-matrices were calculated

    real(dp) :: xs_zero_threshold
      !! Any cross section with value only smaller than this (cm²) will
      !! be ignore and will not be printed
    real(dp) :: eta_thresh
      !! The largest value of η' allowed for evaluating the hypergeometric functions ₂F₁(a,b;c;z)
    real(dp) :: Ef
      !! The last  electron energy for excitation to consider relative to the initial state's energy
    real(dp) :: Ei_xtrap
      !! The lowest electron energy for de-excitation relative to the initial state's energy.
      !! The results will be extrapolated from Ei down to this assuming a 1/E dependence for the
      !! cross section, i.e., constant excitation probablility. If this .le. 0, no extrapolation will not be performed.
      !! Units: (eV)
    real(dp) :: kmat_energy_closest
      !! Input K-matrices are evaluated at a specific energy. If this code is run energy-independently
      !! (most likely the case unless I add energy dependence in the future) The K-matrix that is
      !! selected will be the FIRST ONE whose evaluation energy is CLOSEST to this energy in (eV).
      !! NOTE: UKRMOL+ outputs K-matrix energies in the .kmat files in Rydberg.
      !! If an energy dependent calculation is selected, this is ignored
    real(dp) :: kmat_ei
      !! The lowest K-matrix energy to consider for reading in case of EDFT. (eV)
      !! If 0.0, start at lowest available energy.
    real(dp) :: kmat_ef
      !! The highest K-matrix energy to consider for reading in case of EDFT. (eV)
      !! If 0.0, go to highest available energy
    real(dp) :: abc(3)
      !! Array of reals of length 3
      !! The rotational constants A, B, and C of the target molecule (cm⁻¹).
    real(dp) :: B_rot
      !! Only for linear rotors; the rotational constant B in the expansion
      !! of the rotational energy :
      !!   E(N) = B N(N+1) - D[N(N+1)]² + H[N(N+1)]³ ...
      !! Cannot be 0 for a linear molecule
    real(dp) :: D_rot
      !! Only for linear rotors; the centrifugal distortion coefficient D
      !! in the expansion of the rotational energy :
      !!   E(N) = B N(N+1) - D[N(N+1)]² + H[N(N+1)]³ ...
      !! 0 by default
    real(dp) :: H_rot
      !! Only for linear rotors; the centrifugal distortion coefficient D
      !! in the expansion of the rotational energy :
      !!   E(N) = B N(N+1) - D[N(N+1)]² + H[N(N+1)]³ ...
      !! 0 by default
    real(dp) :: cartesian_dipole_moments(3)
      !! Array of cartesian dipole moments (Debye)
      !! in the order dx, dy, dz
    real(dp) :: cartesian_quadrupole_moments(6)
      !! Array of cartesian quadrupole moments (Debye)
      !! in the order Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
    real(dp), allocatable :: egrid_segs(:)
      !! Array of the bounds (non-degenerate) of the energy grid segments (length num_egrid_segs + 1)
    real(dp) :: post_rft_sincos2s_imag_tol
      !! After the energy dependent frame transformation, the sin and cosine matrices are expected
      !! to have a imaginary values that is no larger than this in magnitude

    character(1) :: rotor_kind
      !! The kind of rotor that describes the targer. Character(1).
      !! Choice of :
      !!  "l"inear
      !!  "a"symmetric top
      !!  "s"ymmetric  top
    character(1) :: rotor_zaxis
      !! The molecular axis (a, b, or c) along which the z-axis is oriented
      !! For asymmetric tops, this should be the main symmetry axis
      !! For symmetric tops, this should be highest symmetry axis
    character(1) :: rotor_c2axis
      !! The C₂ symmetry axis. UKRMol+ calculations for certain point groups require z to be
      !! along the C₂ axis, in which case we'll need to rotate our rotational eigenvectors.
      !! Probably mostly relevant for symmetric tops
    character(1) :: channel_energy_units_override
      !! The units of the channel energies in the file that holds channels. Options are :
      !!  - "r" for Rydberg, "h" for hartree, "e" for eV
      !! By default, this is not set and will allow the code to determine
      !! channel energies on its own based on KMAT_OUTPUT_TYPE, but can
      !! be forcibly overridden with this
    character(1) :: kmat_energy_units_override
      !! The units of the K-matrix evaluation energies in the kmat file. Options are :
      !!  - "r" for Rydberg, "h" for hartree, "e" for eV
      !! By default, this is not set and will allow the code to determine
      !! channel energies on its own based on KMAT_OUTPUT_TYPE, but can
      !! be forcibly overridden with this
    character(3) :: egrid_spacing
      !! The kind of spacing for the energy grid segments. "lin" for linear and "log" for logarithmic
    character(1) :: egrid_xtrap_pre
      !! What kind of extrapolation to perform for the K-matrix evaluation energy grid for scattering energies
      !! that are below the evaluation energy grid that we read.
      !! "l"inear
      !! "c"onstant
    character(1) :: egrid_xtrap_post
      !! What kind of extrapolation to perform for the K-matrix evaluation energy grid for scattering energies
      !! that are above the evaluation energy grid that we read.
      !! "l"inear
      !! "c"onstant
    character(:), allocatable :: point_group
      !! The point group in which the K-matrices were calculated
    character(:), allocatable :: kmat_dir
      !! Path for the file containing the K-matrix to be read. Absolute or relative
    character(:), allocatable :: channels_dir
      !! Path for the file containing the channels for the K-matrix to be read. Absolute or relative
      !! This is only used if kmat_output_type is ukrmol+ because the channel and K-matrix files are separate
    character(:), allocatable :: output_directory
      !! The directory in which to write the output data
      !! This directory must already exist
    character(:), allocatable :: CDMS_file
      !! The file containing CDMS transitions
    character(7) :: kmat_output_type
      !! Determines what kind of K-matrices we're reading. Two possible values:
      !!   'UKRMOL+': default UKRmol+ .kmat file
      !!   'MQDTR2K': a specific format given in the writeup. K-matrices are generated
      !!     directly from the R-matrix, possibly with channel elimination and differently
      !!     normalized Coulomb wavefunctions

    type(cd4_type) :: cd4
      !! Centrifugal distortion parameters (4th order) for the rigid rotor Hamiltonian correction
    type(cd6_type) :: cd6
      !! Centrifugal distortion parameters (6th order) for the rigid rotor Hamiltonian correction

  end type config_type
  ! -- ↑↑↑↑↑↑↑↑↑↑↑↑↑ namelist variable declarations ↑↑↑↑↑↑↑↑↑↑↑↑↑

  type(config_type), public :: G
    !! G for Global. All global variables live in here so that we easily have access to all of them
    !! whenever we want

  ! character(1), parameter :: GLOBAL_ROTOR_KIND_VALS(*)   = ["a", "s", "l"]
  ! character(1), parameter :: GLOBAL_ROTOR_ZAXIS_VALS(*)  = ["a", "b", "c"]
  ! character(1), parameter :: GLOBAL_ROTOR_C2AXIS_VALS(*) = ["a", "b", "c"]
  ! integer,      parameter :: GLOBAL_FORBIDDEN_STATES_KIND_VALS(*) = [0, 1]

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! module subroutine setglobals( rotor_kind       &
  !                               , rotor_zaxis      &
  !                               , c2axis           &
  !                               , spin_isomer_kind &
  !                               , forbidden_states_kind &
  !   )
  !   !! Copies values from CFG into this module's globals
  !   use rotex__characters, only: lower
  !   implicit none (type, external)
  !   character(1), intent(in) :: rotor_kind
  !   character(1), intent(in) :: rotor_zaxis
  !   character(1), intent(in) :: c2axis
  !   integer,      intent(in) :: spin_isomer_kind
  !   integer,      intent(in) :: forbidden_states_kind
  !   GLOBAL_ROTOR_KIND            = lower(rotor_kind)
  !   GLOBAL_ROTOR_ZAXIS           = lower(rotor_zaxis)
  !   GLOBAL_C2AXIS                = lower(c2axis)
  !   GLOBAL_SPIN_ISOMER_KIND      = spin_isomer_kind
  !   GLOBAL_FORBIDDEN_STATES_KIND = forbidden_states_kind
  ! end subroutine setglobals

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine checkglobals
  !   use rotex__system, only: stdout, stderr
  !   implicit none (type, external)
  !   select case(GLOBAL_ROTOR_KIND)
  !   case("a", "s", "l")
  !     write(stdout, '("GLOBAL_ROTOR_KIND: ", A, " ✔️")') GLOBAL_ROTOR_KIND
  !   case default
  !     write(stderr, '("GLOBAL_ROTOR_KIND: ", A, " ❌")') GLOBAL_ROTOR_KIND
  !   end select
  ! end subroutine checkglobals

! ================================================================================================================================ !
end module rotex__globals
! ================================================================================================================================ !
