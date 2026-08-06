! ================================================================================================================================ !
module rotex__globals
  !! Putting some important global variables in here so we don't have to pass them everywhere
  !! all the time
  use rotex__kinds, only: dp
  use rotex__types, only: cd4_type, cd6_type
  use rotex__system, only: stdin, stdout, stderr, die

  implicit none (type, external)

  private

  save

  public :: read_namelists
  public :: is_unset

  ! public :: setglobals
  public :: config_type

  integer, parameter, public :: DEFAULT_INT      = huge(1)
  real(dp), parameter, public :: DEFAULT_REAL = huge(1.0_dp)
  integer, parameter, public :: SJ_COMPUTE       = 0
    !! Value for SJ_mode: compute
  integer, parameter, public :: SJ_COMPUTE_WRITE = 1
    !! Value for SJ_mode: compute + write
  integer, parameter, public :: SJ_READ          = 2
    !! Value for SJ_mode: read

  real(dp), parameter, public :: CB_MINT_IMAG_THRESH = 1e-8_dp
    !! Warn user if CB cross section M integrals have imaginary components
    !! larger than this value (they should be real). The code will try to skip such values.
  real(dp), parameter, public :: HYPGEO_ABC_THRESHOLD = 160._dp
    !! Threshold for warning the user about large values of a,b,c in
    !! ₂F₁(a,b;c;z)

  character(1), parameter, public :: PACKMAT_TRIANGLE  = "L"
  character(*), parameter, public :: CHAR_CR           = achar(13)
  character(*), parameter, public :: DEFAULT_CHAR1     = "x"
  character(*), parameter, public :: UKRMOLX           = 'ukrmol+'
  character(*), parameter, public :: MQDTR2K           = 'mqdtr2k'
  character(*), parameter, public :: SPINMULT_NAMES(5) = &
    [ 'singlet', 'doublet', 'triplet', 'quartet', 'quintet' ]

  interface is_unset
    module procedure :: is_unset_i
    module procedure :: is_unset_r
    module procedure :: is_unset_char
  end interface is_unset

  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓ namelist variables declarations ↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  type config_type
    !! Derived type containing data from the namelist variables

    logical :: ALLOW_EDFT_EGRID_OUT_OF_BOUNDS
      !! If this is true, clamp the total energy grid to the K-matrix evaluation
      !! energy grid. This effectively results in an energy independent treatment past the endpoints.
      !! If this is false, any energies that are outside of the evaluation energy grid will
      !! result in errors.
    logical :: ADD_CD4 = .false.
      !! Add centrifugal distortion for fourth order ?
    logical :: ADD_CD6 = .false.
      !! Add centrifugal distortion for sixth order ?
    logical :: ONLY_EINSTA
      !! Whether to only calculate the Einstein A coefficients in the Coulomb-Born
      !! cross section routine
    logical :: USE_CDMS_EINSTA
      !! Whether to use Einstein A coefficients obtained from the CDMS or calculate them ourselves
      !! Path for the file containing the CDMS data to be read (if use_CDMS_einstA is .true.)
    logical :: ANALYTIC_TOTAL_CB(2)
      !! Array of logicals that has the size 2. Choose whether to use the analytic equation
      !! describing the multipole expansions for the dipole (element 1) and the quadrupole (element 2, not yet available)
    logical :: DO_XTRAP
      !! Do the extrapolation of (de-)excitation cross sections as 1/E to the excitation threshold ?
    logical :: DO_DIPOLE
      !! Choose whether to use the dipole term of the potential expansion
    logical :: DO_QUADRUPOLE
      !! Choose whether to use the quadrupole term of the potential expansion
    logical :: USE_KMAT
      !! Calculate (de-)excitation cross sections using precomputed K-matrices ?
    logical :: USE_CB
      !! Calculate (de-)excitation cross sections using the Coulomb-Born approxiation ?
    logical :: REAL_SPHERICAL_HARMONICS
      !! Whether the input K-matrices are evaluated in a basis of real spherical harmonics
      !! for the scattering electron. If .true., it will be transformed to a basis of
      !! complex-valued spherical harmonics
    logical :: PRINT_ROT_STATES
      !! Whether to print rotational target states once they are defined and ordered
    logical :: PRINT_ELEC_CHANNELS
      !! Whether to print electronic channels once they are defined and ordered
    logical :: PRINT_MEMINFO
      !! Whether to print memory storage information during the EDFT/MQDT step
    logical :: PRINT_CHUNKINFO
      !! Whether to print chunk information during the EDFT/MQDT step
    logical :: ENFORCE_SPIN_ISOMER
      !! Whether do enforce symmetry with respect to nuclear spin. If the calculation is done
      !! in an Abelian group with the highest possible symmetry, e.g., C2v for H3, then this should
      !! be recoverable simply by specifying the symmetry axis, scattering, point group, and target
      !! point group. This is .FALSE. by default
    logical :: PROJECT_ELECMAT_ONTO_TARGET_PG
      !! Sometimes the K-matrix is calculated in a lower symmetry than that of the target.
      !! For example, H3+ would be a C2v calculation, but the target is a  D3h molecule.
      !! This is because C2v is the highest Abelian subgroup of D3h. The K-matrix will be
      !! rotated so that the z-axis is along the C3 symmetry axis. After the rotation, the K-matrix
      !! can be projected onto this higher point group (consists of setting elements to 0 if they are
      !! not in the same symmetry class in the higher point group; symmetry class determined by each
      !! channel's spherical harmonic projection λ).

    integer :: SJ_MODE
      !! The calculation mode for the frame-transformed S-matrix. The general structure is
      !! 0 (SJ_COMPUTE): keep it all in memory (only one S^J subblock is stored at a time)
      !! 1 (SJ_COMPUTE_WRITE) write the transition probabilities to disk for each J
      !! 2 (SJ_READ): read the transition probabilities from disk for each J and DO NOT COMPUTE
    integer :: FORBIDDEN_STATES_KIND
      !! Some molecules have rotatinonal levels that are forbidden in certain vibrational states, e.g., H₃⁺
      !! In this case, do not construct certain levels:
      !!   0: no restriction
      !!   1: even N, K=0 states cannot exist symmetric top only
    logical :: SYMTOP_REDUCE_PROJECTION
      !! The code internally calculates states for symmetric tops with ±K. If this is
      !! .true., ruce all of these to |K|. For example, the four transitions ±1 -> ±3 will
      !! be reduced to the single transition 1 -> 3
    logical :: EDFT
      !! Use an energy dependent frame transformation for the K/S-matrices ? This will read all available
      !! energies provided between kmat_Ei and kmat_Ef
    integer :: NE
      !! The number of scattering energies to consider. This does not need to be very high; the CB
      !! cross sections are very smooth and can easily be interpolated.
    integer :: NE_XTRAP
      !! Number of extrapolation energies. Excitation cross sections are extrapolated as 1/E to the excitation threshold,
      !! de-excitation cross sections are extrapolated as 1/E to Ei_xtrap. If this is 0, no exptrapolation will be performed.
    integer :: LMAX_PARTIAL
      !! The maximum value of l to consider in the contribution of the partial CB cross section
      !! from the dipole and the quadrupole. If you're replacing the low-l CB cross sections
      !! with other cross sections, set this to the max l that you have available.
    integer :: LMAX_TOTAL
      !! The maximum value of l to consider in the contribution of the total CB cross section
      !! in the even that you're not using the analytic expression,
      !! from the dipole  and the quadrupole
    integer :: NMIN
      !! The minimum value of the rotational quantum number (N) to consider
    integer :: NMAX
      !! The maximum value of the rotational quantum number (N) to consider
    integer :: TARGCHARGE
      !! The electric charge of the target
    integer :: LMAX_KMAT
      !! The max partial wave to be included in the K-matrix basis. Cannot exceed the available
      !! basis in the calculation, but can be smaller than the largest available partial wave
    integer :: NUM_EGRID_SEGS
      !! Number of energy grid segments (evaluation energy for the cross sections)
    integer, allocatable :: NUM_EGRID(:)
      !! Array of number of energies per grid segment (length num_egrid_segs)
    integer, allocatable :: SPINMULTS(:)
      !! Array of spin multiplicities (2S+1) for which the system's (target + e⁻) K-matrices were calculated
    integer :: EDFT_CHUNK_TARGET_MB
      !! The approximate memory budget in MB for the temporary transformed rotatinal S-matrix buffer used
      !! fduring the EDFT for a given J block. If this value is non-positive, chunking is disabled and the full
      !! evaluation energy grid is stored in memory.

    real(dp) :: C2PRIME_PHI_DEG
      !! The angle φ between the SYMAXIS z-axis and the C₂' z-axis. For example, this is π/2 in the
      !! case of H₃⁺. This can be left unspecified, but is presented as an override. For calculations
      !! making use of K-matrices where there is a SYMAXIS and a SCATTERING frame, the code will
      !! attempt to determine this automatically. If it is somehow incorrect, the supplied value can
      !! be used instead. In the case of multipole-only calculations (e.g., CB with no Kmat), then
      !! this can be omitted and will default to φ=0. Degrees
    real(dp) :: XS_ZERO_THRESHOLD
      !! Any cross section with value only smaller than this (cm²) will
      !! be ignore and will not be printed
    real(dp) :: ETA_THRESH
      !! The largest value of η' allowed for evaluating the hypergeometric functions ₂F₁(a,b;c;z)
    real(dp) :: EF
      !! The last  electron energy for excitation to consider relative to the initial state's energy
    real(dp) :: EI_XTRAP
      !! The lowest electron energy for de-excitation relative to the initial state's energy.
      !! The results will be extrapolated from Ei down to this assuming a 1/E dependence for the
      !! cross section, i.e., constant excitation probablility. If this .le. 0, no extrapolation will not be performed.
      !! Units: (eV)
    real(dp) :: KMAT_ENERGY_CLOSEST
      !! Input K-matrices are evaluated at a specific energy. If this code is run energy-independently
      !! (most likely the case unless I add energy dependence in the future) The K-matrix that is
      !! selected will be the FIRST ONE whose evaluation energy is CLOSEST to this energy in (eV).
      !! NOTE: UKRMOL+ outputs K-matrix energies in the .kmat files in Rydberg.
      !! If an energy dependent calculation is selected, this is ignored
    real(dp) :: KMAT_EI
      !! The lowest K-matrix energy to consider for reading in case of EDFT. (eV)
      !! If 0.0, start at lowest available energy.
    real(dp) :: KMAT_EF
      !! The highest K-matrix energy to consider for reading in case of EDFT. (eV)
      !! If 0.0, go to highest available energy
    real(dp) :: ABC(3)
      !! Array of reals of length 3
      !! The rotational constants A, B, and C of the target molecule (cm⁻¹).
    real(dp) :: B_ROT
      !! Only for linear rotors; the rotational constant B in the expansion
      !! of the rotational energy :
      !!   E(N) = B N(N+1) - D[N(N+1)]² + H[N(N+1)]³ ...
      !! Cannot be 0 for a linear molecule
    real(dp) :: D_ROT
      !! Only for linear rotors; the centrifugal distortion coefficient D
      !! in the expansion of the rotational energy :
      !!   E(N) = B N(N+1) - D[N(N+1)]² + H[N(N+1)]³ ...
      !! 0 by default
    real(dp) :: H_ROT
      !! Only for linear rotors; the centrifugal distortion coefficient D
      !! in the expansion of the rotational energy :
      !!   E(N) = B N(N+1) - D[N(N+1)]² + H[N(N+1)]³ ...
      !! 0 by default
    real(dp) :: DIPOLE_ABC(3)
      !! Array of cartesian dipole moments (Debye) in the order dA, dB, dC.
    real(dp) :: QUADRUPOLE_XYZ(6)
      !! Array of cartesian quadrupole moments (Debye) in the order Qxx, Qxy, Qxz, Qyy, Qyz, Qzz.
      !! Cannot be given with another QUADRUPOLE_??? array
    ! real(dp) :: QUADRUPOLE_SPH(6)
    !   !! Array of cartesian quadrupole moments (Debye) in the order Q(-2)..Q(2)
    !   !! Cannot be given with another QUADRUPOLE_??? array
    real(dp), allocatable :: EGRID_SEGS(:)
      !! Array of the bounds (non-degenerate) of the energy grid segments (length num_egrid_segs + 1)

    character(1) :: ROTOR_KIND
      !! The kind of rotor that describes the targer. Character(1).
      !! Choice of :
      !!  "l"inear
      !!  "a"symmetric top
      !!  "s"ymmetric  top
    character(1) :: RR_DIAG_AXIS
      !! Rigid rotor diagonalization axis
    character(1) :: SYMAXIS
      !! The axis of highest symmetry: one of "A", "B", "C"
    character(1) :: SCATTERING_ZAXIS
      !! The z-axis of the scattering calculations: one of "A", "B", "C"
    character(1) :: SCATTERING_YAXIS
      !! The y-axis of the scattering calculations: one of "A", "B", "C"
    character(1) :: SCATTERING_XAXIS
      !! The x-axis of the scattering calculations: one of "A", "B", "C"
    character(1) :: CHANNEL_ENERGY_UNITS_OVERRIDE
      !! The units of the channel energies in the file that holds channels. Options are :
      !!  - "r" for Rydberg, "h" for hartree, "e" for eV
      !! By default, this is not set and will allow the code to determine
      !! channel energies on its own based on KMAT_OUTPUT_TYPE, but can
      !! be forcibly overridden with this
    character(1) :: KMAT_ENERGY_UNITS_OVERRIDE
      !! The units of the K-matrix evaluation energies in the kmat file. Options are :
      !!  - "r" for Rydberg, "h" for hartree, "e" for eV
      !! By default, this is not set and will allow the code to determine
      !! channel energies on its own based on KMAT_OUTPUT_TYPE, but can
      !! be forcibly overridden with this
    character(3) :: EGRID_SPACING
      !! The kind of spacing for the energy grid segments. "lin" for linear and "log" for logarithmic
    character(1) :: EGRID_XTRAP_PRE
      !! What kind of extrapolation to perform for the K-matrix evaluation energy grid for scattering energies
      !! that are below the evaluation energy grid that we read.
      !! "l"inear
      !! "c"onstant
    character(1) :: EGRID_XTRAP_POST
      !! What kind of extrapolation to perform for the K-matrix evaluation energy grid for scattering energies
      !! that are above the evaluation energy grid that we read.
      !! "l"inear
      !! "c"onstant
    character(:), allocatable :: SCATTERING_POINT_GROUP
      !! The point group in which the K-matrices were calculated
    character(:), allocatable :: TARGET_POINT_GROUP
      !! The point group of the target molecule
    character(:), allocatable :: KMAT_DIR
      !! Path for the file containing the K-matrix to be read. Absolute or relative
    character(:), allocatable :: CHANNELS_DIR
      !! Path for the file containing the channels for the K-matrix to be read. Absolute or relative
      !! This is only used if kmat_output_type is ukrmol+ because the channel and K-matrix files are separate
    character(:), allocatable :: OUTPUT_DIRECTORY
      !! The directory in which to write the output data
      !! This directory must already exist
    character(:), allocatable :: CDMS_FILE
      !! The file containing CDMS transitions
    character(7) :: KMAT_OUTPUT_TYPE
      !! Determines what kind of K-matrices we're reading. Two possible values:
      !!   'UKRMOL+': default UKRmol+ .kmat file
      !!   'MQDTR2K': a specific format given in the writeup. K-matrices are generated
      !!     directly from the R-matrix, possibly with channel elimination and differently
      !!     normalized Coulomb wavefunctions

    type(cd4_type) :: CD4
      !! Centrifugal distortion parameters (4th order) for the rigid rotor Hamiltonian correction
    type(cd6_type) :: CD6
      !! Centrifugal distortion parameters (6th order) for the rigid rotor Hamiltonian correction

  end type config_type
  ! -- ↑↑↑↑↑↑↑↑↑↑↑↑↑ namelist variable declarations ↑↑↑↑↑↑↑↑↑↑↑↑↑

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(config_type), public :: G
    !! G for Global. All global variables live in here so that we easily have access to all of them
    !! whenever we want
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_namelists
    !! Reads user parameters and puts them into the config derived type
    use rotex__kinds,      only: dp
    use rotex__arrays,     only: append, remove_value
    use rotex__system,     only: ds => directory_separator
    use rotex__constants,  only: au2invcm, au2ev, macheps => macheps_dp, au2cm, au2deb
    use rotex__pointgroups, only: is_supported_pg, is_abelian_pg, is_subgroup, write_pg_table
    use rotex__characters, only: add_trailing, to_lower, lower

    implicit none (type, external)

    logical :: dipole_abc_isgiven
    real(dp),    parameter :: DEFAULT_DIPOLE = DEFAULT_REAL

    ! -- namelist: control
    integer :: Nmin
    integer :: Nmax
    logical :: use_kmat
    logical :: use_CB
    logical :: symtop_reduce_projection = .true.
    logical :: enforce_spin_isomer = .false.
    logical :: project_elecmat_onto_target_pg = .false.
    integer :: forbidden_states_kind = 0
    character(:), allocatable :: output_directory
    character(1) :: rotor_kind = DEFAULT_CHAR1
    character(1) :: symaxis = DEFAULT_CHAR1
    character(:), allocatable :: target_point_group
    real(dp) :: c2prime_phi_deg = DEFAULT_REAL
    real(dp) :: abc(3) = 0.0_dp
    real(dp) :: B_rot = 0.0_dp
    real(dp) :: H_rot = 0.0_dp
    real(dp) :: D_rot = 0.0_dp
    integer :: targcharge = DEFAULT_INT
    logical :: add_cd4 = .false.
    logical :: add_cd6 = .false.
    ! -- cd4
    real(dp) :: dn     = 0.0_dp
    real(dp) :: dnk    = 0.0_dp
    real(dp) :: dk     = 0.0_dp
    real(dp) :: deltan = 0.0_dp
    real(dp) :: deltak = 0.0_dp
    ! -- cd6
    real(dp) :: hn     = 0.0_dp
    real(dp) :: hnk    = 0.0_dp
    real(dp) :: hkn    = 0.0_dp
    real(dp) :: hk     = 0.0_dp
    real(dp) :: etan   = 0.0_dp
    real(dp) :: etank  = 0.0_dp
    real(dp) :: etak   = 0.0_dp

    ! -- namelist: kmat
    logical :: real_spherical_harmonics = .true.
    logical :: allow_edft_egrid_out_of_bounds  = .false.
    logical :: edft = .false.
    integer :: lmax_kmat = DEFAULT_INT
    integer :: num_egrid_segs
    integer :: edft_chunk_target_mb = 2048
    integer, allocatable :: num_egrid(:)
    integer, allocatable :: spinmults(:)
    real(dp), allocatable :: egrid_segs(:)
    real(dp) :: kmat_Ei = 0._dp, kmat_Ef = 0._dp
    character(1) :: channel_energy_units_override = DEFAULT_CHAR1
    character(1) :: kmat_energy_units_override    = DEFAULT_CHAR1
    character(1) :: scattering_zaxis = DEFAULT_CHAR1
    character(1) :: scattering_xaxis = DEFAULT_CHAR1
    character(1) :: scattering_yaxis = DEFAULT_CHAR1
    character(3) :: egrid_spacing
    character(7) :: kmat_output_type = "======="
    character(:), allocatable :: scattering_point_group
    character(:), allocatable :: kmat_dir
    character(:), allocatable :: channels_dir

    ! -- namelist: coulomb
    logical :: use_CDMS_einstA = .false.
    logical :: only_einsta = .false.
    logical :: do_xtrap = .false.
    logical :: do_dipole = .true.
    logical :: do_quadrupole = .false.
    logical :: analytic_total_cb(2)
    integer :: nE = DEFAULT_INT
    integer :: SJ_MODE = SJ_COMPUTE
    integer :: nE_xtrap = DEFAULT_INT
    integer :: lmax_partial = DEFAULT_INT
    integer :: lmax_total = DEFAULT_INT
    real(dp) :: eta_thresh = 0.0_dp
    real(dp) :: Ef = 0.0_dp
    real(dp) :: Ei_xtrap = 0.0_dp
    real(dp) :: dipole_abc(3) = DEFAULT_DIPOLE
    ! real(dp) :: qudarupole_xyz(6)
    ! real(dp) :: qudarupole_sph(6)
    real(dp) :: xs_zero_threshold = 0.0_dp   ! include all cross sections by default
    real(dp) :: kmat_energy_closest = 0.0_dp ! just take the first one
    character(1) :: egrid_xtrap_pre, egrid_xtrap_post
    character(:), allocatable :: CDMS_file

    namelist / control_namelist /           &
      !! Contains parameters and values that are necessary to run the program
        output_directory                    &
      , forbidden_states_kind               &
      , symtop_reduce_projection            &
      , target_point_group                  &
      , nmin                                &
      , nmax                                &
      , use_kmat                            &
      , use_cb                              &
      , symaxis                             &
      , enforce_spin_isomer                 &
      , rotor_kind                          &
      , targcharge                          &
      , abc                                 &
      , c2prime_phi_deg                     &
      , B_rot                               &
      , D_rot                               &
      , H_rot                               &
      , add_cd4                             &
      , add_cd6                             &
      , dn, dnk, dk, deltan, deltak         &
      , hn, hnk, hkn, hk, etan, etank, etak &
      , xs_zero_threshold

    namelist / kmat_namelist /         &
      !! Parameters regarding the K-matrces used for (de-excitation)
        kmat_dir                       &
      , channels_dir                   &
      , lmax_kmat                      &
      , num_egrid_segs                 &
      , num_egrid                      &
      , egrid_xtrap_pre                &
      , egrid_xtrap_post               &
      , scattering_point_group         &
      , egrid_segs                     &
      , project_elecmat_onto_target_pg      &
      , edft                           &
      , kmat_ei                        &
      , kmat_ef                        &
      , scattering_zaxis               &
      , scattering_xaxis               &
      , scattering_yaxis               &
      , egrid_spacing                  &
      , spinmults                      &
      , kmat_output_type               &
      , kmat_energy_closest            &
      , real_spherical_harmonics       &
      , channel_energy_units_override  &
      , SJ_mode                        &
      , allow_edft_egrid_out_of_bounds &
      , edft_chunk_target_mb           &
      , kmat_energy_units_override

    namelist / coulomb_namelist /                 &
      !! Parameters regarding the Coulomb-Born approximation
      !! used for (de-)excitation
        use_CDMS_einstA          &
      , only_einsta              &
      , cdms_file                &
      , eta_thresh               &
      , ef                       &
      , ne                       &
      , ne_xtrap                 &
      , do_xtrap                 &
      , ei_xtrap                 &
      , dipole_abc               &
      ! , quadrupole_xyz           &
      ! , quadrupole_sph           &
      , do_dipole                &
      , do_quadrupole            &
      , analytic_total_cb        &
      , lmax_partial             &
      , lmax_total

    !!!!!!!!!!!!!!!!!!!!!! CONTROL_NAMELIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(character(1000) :: output_directory)
    allocate(character(10)   :: target_point_group)
    ! -- read
    read(stdin, control_namelist)
    ! -- check defaults
    if(Nmin          .eq. DEFAULT_INT)   call die("Must specify NMIN in CONTROL_NAMELIST")
    if(Nmax          .eq. DEFAULT_INT)   call die("Must specify NMAX in CONTROL_NAMELIST")
    if(rotor_kind    .eq. DEFAULT_CHAR1) call die("Must specify ROTOR_KIND in CONTROL_NAMELIST")
    if(targcharge .eq. DEFAULT_INT)   call die("Must specify TARGCHARGE in CONTROL_NAMELIST")
    if(symaxis             .eq. DEFAULT_CHAR1) call die("Must specify SYMAXIS in CONTROL_NAMELIST")
    if(lower(rotor_kind) .eq. "l") then
      if(B_rot .le. 0.0_dp) call die("Must have a positive rotational constant B_rot for a linear molecule")
    else
      if(any(ABC .eq. 0.0_dp)) call die("Must specify nonzero rotational constants ABC in CONTROL_NAMELIST")
    endif
    call to_lower(target_point_group)
    if(allocated(target_point_group)) target_point_group  = trim(target_point_group)
    if(is_supported_pg(target_point_group) .eqv. .false.) then
      call write_pg_table(stderr)
      call die("Supplied TARGET_POINT_GROUP is not valid: " // target_point_group)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! KMAT_NAMELIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(use_kmat .eqv. .true.) then
      ! -- prepare for reading
      rewind(stdin)
      allocate(character(1000) :: kmat_dir)
      allocate(character(1000) :: channels_dir)
      allocate(num_egrid(100))
      allocate(egrid_segs(101))
      allocate(character(10)   :: scattering_point_group)
      allocate(spinmults(10))
      spinmults = DEFAULT_INT
      kmat_dir(1:1)     = DEFAULT_CHAR1
      channels_dir(1:1) = DEFAULT_CHAR1
      kmat_output_type(1:1) = DEFAULT_CHAR1
      ! -- read
      read(stdin, kmat_namelist)
      ! -- trim arrays
      num_egrid  = num_egrid(1:num_egrid_segs)
      egrid_segs = egrid_segs(1:num_egrid_segs+1) / au2ev
      ! -- normalize characters
      call to_lower(egrid_spacing)
      call to_lower(kmat_output_type)
      select case(kmat_output_type)
      case(UKRMOLX, MQDTR2K)
        continue
      case default
        call die("KMAT_OUTPUT_TYPE in KMAT_NAMELIST must be one of " // UKRMOLX // " or " // MQDTR2K)
      end select
      select case(egrid_spacing)
        case("lin", "log") ; continue
        case default ; call die("EGRID_SPACING (" // egrid_spacing // ") must be LIN or LOG in KMAT_NAMELIST")
      end select
      ! -- check defaults
      if(kmat_dir(1:1) .eq. DEFAULT_CHAR1) call die("Must specify KMAT_DIR in KMAT_NAMELIST")
      if(scattering_zaxis    .eq. DEFAULT_CHAR1) call die("Must specify scattering_zaxis in CONTROL_NAMELIST")
      if(lmax_kmat .eq. DEFAULT_INT .OR. lmax_kmat .lt. 0) &
        call die("LMAX_KMAT in KMAT_NAMELIST must be defined and be non-negative")
      call remove_value(spinmults, DEFAULT_INT)
      if(any(spinmults .lt. 1)) call die("SPINMULTS in KMAT_NAMELIST cannot have values that are < 1")
      if(     channels_dir(1:1) .eq. DEFAULT_CHAR1 &
        .AND. kmat_output_type  .eq. UKRMOLX) call die("Must specify CHANNELS_DIR in KMAT_NAMELIST with&
          & KMAT_OUPUT_TYPE = " // UKRMOLX)
      if(is_supported_pg(scattering_point_group) .eqv. .false.) then
        call write_pg_table(stderr)
        call die("Supplied SCATTERING_POINT_GROUP is not valid: " // scattering_point_group)
      elseif(is_subgroup(scattering_point_group, target_point_group) .eqv. .false.) then
        call write_pg_table(stderr)
        call die("SCATTERING_POINT_GROUP ("//scattering_point_group//") is not a valid subgroup&
          & of TARGET_POINT_GROUP ("//target_point_group//")")
      endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! COULOMB_NAMELIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(use_CB .eqv. .true.) then

      dipole_abc_isgiven = .false.

      allocate(character(1000) :: cdms_file)
      ! -- prepare for reading
      cdms_file(1:1) = DEFAULT_CHAR1
      rewind(stdin)
      ! -- read
      read(stdin, coulomb_namelist)
      ! -- check values
      if(use_CDMS_einstA) then
        if(CDMS_file(1:1) .eq. DEFAULT_CHAR1) call die("Must define CDMS_FILE in COULOMB_NAMELIST&
          & when USE_CDMS_EINSTA is .TRUE.")
      endif
      if(do_dipole .eqv. .false. .AND. do_quadrupole) call die("DO_DIPOLE and DO_QUADRUPOLE in COULOMB_NAMELIST should not both&
        & be set to .FALSE.; nothing would be done")
      if(do_quadrupole) call die("DO_QUADRUPOLE in COULOMB_NAMELIST should not be set to true; it is not implemented")
      if(eta_thresh .eq. 0.0_dp) call die("eta_thresh in COULOMB_NAMELIST must be defined and be positive")
      if(Ef .eq. 0.0_dp) call die("EF in COULOMB_NAMELIST must be defined and be positive")
      if(nE .eq. DEFAULT_INT .OR. ne .le. 0) call die("NE in COULOMB_NAMELIST must be defined and positive")
      if(lmax_partial .eq. DEFAULT_INT) call die("LMAX_PARTIAL in COULOMB_NAMELIST must be defined and nonnegative")
      if(lmax_total .eq. DEFAULT_INT .AND. (analytic_total_cb(1) .eqv. .false.)) &
        call die("LMAX_TOTAL in COULOMB_NAMELIST must be defined and nonnegative if ANALYTIC_TOTAL_CB(1) is .false.")
      if(do_xtrap) then
        if(ne_xtrap .eq. DEFAULT_INT .OR. ne_xtrap .le. 0) then
          call die("NE_XTRAP in COULOMB_NAMELIST must be defined and positive if DO_XTRAP is .TRUE.")
        endif
        if(Ei_xtrap .eq. 0.0_dp) then
          call die("EI_XTRAP in COULOMB_NAMELIST must be defined, and nonzero if DO_XTRAP is .TRUE.")
        endif
      endif
      ! -- only exactly one of the following must be true
      dipole_abc_isgiven = any(dipole_abc .ne. DEFAULT_DIPOLE)
      if(dipole_abc_isgiven .eqv. .false.) call die("Please supply DIPOLE_ABC")
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call to_lower(symaxis)
    call to_lower(scattering_zaxis)
    call to_lower(scattering_xaxis)
    call to_lower(scattering_yaxis)

    ! -- convert to atomic units
    Ef                = Ef                / au2ev
    Ei_xtrap          = Ei_xtrap          / au2ev
    ABC(:)            = ABC(:)            / au2invcm
    B_rot             = B_rot             / au2invcm
    D_rot             = D_rot             / au2invcm
    H_rot             = H_rot             / au2invcm
    kmat_ei           = kmat_ei           / au2ev
    kmat_ef           = kmat_ef           / au2ev
    xs_zero_threshold = xs_zero_threshold / (au2cm*au2cm)

    ! -- convert to lower case
    call to_lower(rotor_kind)
    if(use_kmat .eqv. .true.) call to_lower(scattering_point_group)
    call to_lower(kmat_energy_units_override)
    call to_lower(channel_energy_units_override)

    ! -- remove spaces
    if(allocated(scattering_point_group)) scattering_point_group = trim(scattering_point_group)
    if(allocated(kmat_dir))               kmat_dir               = trim(kmat_dir)
    if(allocated(channels_dir))           channels_dir           = trim(channels_dir)
    if(allocated(output_directory))       output_directory       = trim(output_directory)

    ! -- add trailing directory separator to directories if needed, make directories as needed
    call add_trailing(output_directory, ds)
    call add_trailing(kmat_dir,         ds)
    call add_trailing(channels_dir,     ds)

    write(stdout, '(A)') "--------------------------------------------------------------------------------------------------------"
    write(stdout, *)
    write(stdout, control_namelist)
    write(stdout, *)
    if(use_kmat .eqv. .true.) then
      write(stdout, kmat_namelist)
      write(stdout, *)
    endif
    if(use_CB .eqv. .true.) then
      write(stdout, coulomb_namelist)
      write(stdout, *)
    endif

    ! -- checks
    if(Nmin .gt. Nmax) call die("Nmin > Nmax not allowed")
    if(targcharge .eq. DEFAULT_INT) call die("Must set the charge of the target in namelist CONTROL !")

    ! -- namelist: control
    G%NMIN                  = nmin
    G%NMAX                  = nmax
    G%USE_KMAT              = use_kmat
    G%USE_CB                = use_cb
    G%FORBIDDEN_STATES_KIND = forbidden_states_kind
    G%SYMTOP_REDUCE_PROJECTION = symtop_reduce_projection
    G%OUTPUT_DIRECTORY      = output_directory
    G%ROTOR_KIND            = rotor_kind
    G%SYMAXIS               = symaxis
    G%SCATTERING_ZAXIS      = scattering_zaxis
    G%SCATTERING_XAXIS      = scattering_xaxis
    G%SCATTERING_YAXIS      = scattering_yaxis
    G%ABC                   = abc(:)
    G%B_ROT                 = b_rot
    G%D_ROT                 = d_rot
    G%H_ROT                 = h_rot
    G%TARGCHARGE            = targcharge
    G%ADD_CD4               = add_cd4
    G%ADD_CD6               = add_cd6
    G%XS_ZERO_THRESHOLD     = xs_zero_threshold
    G%TARGET_POINT_GROUP    = target_point_group
    G%ENFORCE_SPIN_ISOMER   = enforce_spin_isomer
    G%C2PRIME_PHI_DEG       = c2prime_phi_deg
    if(add_cd4 .eqv. .true.) then
      dn      = dn     / au2invcm
      dnk     = dnk    / au2invcm
      dk      = dk     / au2invcm
      deltan  = deltan / au2invcm
      deltak  = deltak / au2invcm
      G%CD4 = cd4_type(dn = dn, dnk = dnk, dk = dk, deltan = deltan, deltak = deltak)
    endif
    if(add_cd6 .eqv. .true.) then
      if(add_cd4 .eqv. .false.) call die("Don't add the sextic correction while omitting the quartic correction !")
      hn    = hn    / au2invcm
      hnk   = hnk   / au2invcm
      hkn   = hkn   / au2invcm
      hk    = hk    / au2invcm
      etan  = etan  / au2invcm
      etank = etank / au2invcm
      etak  = etak  / au2invcm
      G%CD6 = cd6_type(hn = hn, hnk = hnk, hkn = hkn, hk = hk, etan = etan, etank = etank, etak = etak)
    endif

    ! -- namelist: kmat
    if(use_kmat .eqv. .true.) then
      if(kmat_output_type .eq. "=======") then
        call die("Must speficy KMAT_OUTPUT_TYPE. It should be one of "// UKRMOLX //" or "// MQDTR2K)
      elseif(all(kmat_output_type .ne. [UKRMOLX, MQDTR2K])) then
        call die("Poorly specified KMAT_OUTPUT_TYPE. It should be one of "// UKRMOLX //" or "// MQDTR2K)
      endif
      G%KMAT_DIR                       = kmat_dir
      G%CHANNELS_DIR                   = channels_dir
      G%LMAX_KMAT                      = lmax_kmat
      G%SCATTERING_POINT_GROUP         = scattering_point_group
      G%SPINMULTS                      = spinmults(:)
      G%NUM_EGRID_SEGS                 = num_egrid_segs
      G%NUM_EGRID                      = num_egrid(:)
      G%EGRID_SEGS                     = egrid_segs(:)
      G%EGRID_SPACING                  = egrid_spacing
      G%KMAT_EI                        = kmat_ei
      G%KMAT_Ef                        = kmat_ef
      G%EGRID_XTRAP_PRE                = lower(egrid_xtrap_pre)
      G%EGRID_XTRAP_POST               = lower(egrid_xtrap_post)
      G%REAL_SPHERICAL_HARMONICS       = real_spherical_harmonics
      G%KMAT_ENERGY_CLOSEST            = kmat_energy_closest / au2ev
      G%KMAT_OUTPUT_TYPE               = kmat_output_type
      G%KMAT_ENERGY_UNITS_OVERRIDE     = kmat_energy_units_override
      G%CHANNEL_ENERGY_UNITS_OVERRIDE  = channel_energy_units_override
      G%EDFT                           = edft
      G%ALLOW_EDFT_EGRID_OUT_OF_BOUNDS = ALLOW_EDFT_EGRID_OUT_OF_BOUNDS
      G%EDFT_CHUNK_TARGET_MB           = EDFT_CHUNK_TARGET_MB
      G%PROJECT_ELECMAT_ONTO_TARGET_PG   = project_elecmat_onto_target_pg
    endif

    ! -- namelist: coulomb
    if(use_CB .eqv. .true.) then
      if(do_quadrupole) call die("DO_QUADRUPOLE exists as an option, but I'm yet confident in its&
        & implementation. Remove this call if you want and see what happens.")
#ifndef USE_CDMSREADER
      if(use_cdms_einsta .eqv. .true.) call die("User requested use of CDMS data, but the code is&
        & not compiled with that capability. Build with 'USE_CDMSREADER=1' to change this.")
#endif
      G%USE_CDMS_EINSTA   = use_cdms_einsta
      G%ANALYTIC_TOTAL_CB = analytic_total_cb(:)
      G%ETA_THRESH        = eta_thresh
      G%EF                = ef
      G%NE                = ne
      G%NE_XTRAP          = ne_xtrap
      G%EI_XTRAP          = ei_xtrap
      G%DO_XTRAP          = do_xtrap
      G%DO_DIPOLE         = do_dipole
      G%DO_QUADRUPOLE     = do_quadrupole
      G%LMAX_PARTIAL      = lmax_partial
      G%LMAX_TOTAL        = lmax_total
      G%DIPOLE_ABC        = dipole_abc / au2deb
      G%CDMS_FILE         = cdms_file
      G%ONLY_EINSTA       = only_einsta

    endif

    call read_info_namelist(stdin)

  end subroutine read_namelists

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_info_namelist(funit)
    implicit none (type, external)
    integer, intent(in) :: funit
    ! -- namelist: info
    logical :: print_elec_channels = .true.
    logical :: print_rot_states    = .true.
    logical :: print_meminfo       = .false.
    logical :: print_chunkinfo     = .true.
    namelist / info_namelist / &
      !! Parameters used to control the printing of information
      !! during program execution
        print_elec_channels    &
      , print_rot_states       &
      , print_meminfo          &
      , print_chunkinfo
    read(funit, info_namelist)
    rewind(funit)
    G%PRINT_ROT_STATES      = print_rot_states
    G%PRINT_ELEC_CHANNELS   = print_elec_channels
    G%PRINT_MEMINFO         = print_meminfo
    G%PRINT_CHUNKINFO       = print_chunkinfo
    write(stdout, '(A)') "--------------------------------------------------------------------------------------------------------"
    write(stdout, *)
    write(stdout, info_namelist)
    write(stdout, *)
    write(stdout, '(A)') "--------------------------------------------------------------------------------------------------------"
  end subroutine read_info_namelist

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_unset_i(i) result(res)
    implicit none(type, external)
    integer, intent(in) :: i
    logical :: res
    res = .true.
    if(i .eq. DEFAULT_INT) return
    res = .false.
  end function is_unset_i
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_unset_r(x) result(res)
    implicit none(type, external)
    real(dp), intent(in) :: x
    logical :: res
    res = .true.
    if(x .eq. DEFAULT_REAL) return
    res = .false.
  end function is_unset_r
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_unset_char(c) result(res)
    implicit none(type, external)
    character(*), intent(in) :: c
    logical :: res
    res = .true.
    if(c .eq. DEFAULT_CHAR1) return
    res = .false.
  end function is_unset_char

! ================================================================================================================================ !
end module rotex__globals
! ================================================================================================================================ !
