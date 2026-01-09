! ================================================================================================================================ !
module rotex__types
  !! Contains type definitions and procedures for those types used throughout the program
  use rotex__kinds,     only: dp, xdp, qp
  use rotex__constants, only: IQ_DEFAULT

  implicit none

  private

  ! -- kinds
  public :: dp, xdp, qp

  ! -- types
  public :: eigenH_type
  public :: N_states_type
  public :: channel_type
  public :: elec_channel_type
  public :: asymtop_rot_channel_type
  public :: asymtop_rot_channel_l_type
  public :: asymtop_rot_channel_l_vector_type
  public :: ivector_type
  public :: rvector_type
  public :: rmatrix_type
  public :: cmatrix_type
  public :: asymtop_rot_transition_type
  public :: config_type
  public :: cd4_type
  public :: cd6_type

  ! -- procedures
  public :: sort_channels
  public :: permsort_channels
  public :: sort_channels_by_energy
  public :: operator(.eq.)
  public :: operator(.ne.)
  public :: operator(.isin.)
  public :: assignment(=)
  public :: trim_channel_l
  public :: get_channel_index
  public :: findloc_transitions

  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓ type definitions ↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

  type eigenH_type
    !! Contains the eigenvectors and eigenvalues of a hamiltonian
    real(dp), allocatable :: eigvals(:)
    complex(dp), allocatable :: eigvecs(:,:)
  end type eigenH_type

  type N_states_type
    !! The rotational state of the system described by its eigenvectors, eigenvalues,
    !! and state labels
    integer :: N
    type(eigenH_type) :: eigenH
      !! The decomposed Hamiltonian for this rotational level
    integer, allocatable :: Ka(:)
      !! The projections Ka
    integer, allocatable :: Kc(:)
      !! The projections Kc
    real(dp), allocatable :: EinstA(:)
      !! The Einstein coefficients for transitions to all lower states (0 if none)
  end type N_states_type

  type, abstract :: channel_type
    !! |nelec> (E)
    integer :: nelec
      !! Electronic state
    real(dp) :: E
      !! The channel energy
  end type channel_type

  type, extends(channel_type) :: elec_channel_type
    !! Describes the electronic channel basis of the (optional) K and S-matrices by adding l and its projection ml
    !! |nelec,l,ml> (iq, E)
    integer :: l
      !! Partial wave degree
    integer :: ml
      !! Partial wave order (projection on body-frame ẑ-axis)
    integer :: iq = IQ_DEFAULT
      !! The kind of normalization for the Coulomb f/g functions:
      !!   4: usual normalization
      !!   0: f₀/g₀ normalization
  end type elec_channel_type

  type, extends(channel_type) :: asymtop_rot_channel_type
    !! Describes the rotational and electronic channel basis of the (optional) S-matrix
    !! after the rotational frame transformation by adding the and rotational quantum numbers to the channel type
    !! |nelec,N,Ka,Kc> (E)
    integer :: N
      !! The rotatinal quantum number of the target
    integer :: Ka
      !! The projection Ka of N
    integer :: Kc
      !! The projection Kc of N
    integer :: sym
      !! The nuclear spin symmetry
  end type asymtop_rot_channel_type

  type, extends(asymtop_rot_channel_type) :: asymtop_rot_channel_l_type
    !! Describes the rotational and electronic channel basis of the (optional) S-matrix
    !! after the rotational frame transformation by adding the partial wave degree
    !! to the rotational quantum numbers
    !! |nelec,l,N,Ka,Kc> (iq, E)
    integer :: l
      !! Partial wave degree
    integer :: iq
      !! The kind of normalization for the Coulomb f/g functions:
      !!   4: usual normalization
      !!   0: f₀/g₀ normalization
  end type asymtop_rot_channel_l_type

  type asymtop_rot_channel_l_vector_type
    !! Contains a vector of channels. The idea is that this type is indexed at
    !! each value of the angular momentum quantum number J, for which a different
    !! combination of channels exists than for other Js
    type(asymtop_rot_channel_l_type), allocatable :: channels(:)
  end type asymtop_rot_channel_l_vector_type

  type asymtop_rot_transition_type
    !! The a type containing the indices for a pair of initial and final rotational states
    type(asymtop_rot_channel_type) :: lo
    type(asymtop_rot_channel_type) :: up
  end type asymtop_rot_transition_type

  type ivector_type
    !! The type of an integer vector
    integer, allocatable :: vec(:)
  end type ivector_type

  type rvector_type
    !! The type of a real vector
    real(dp), allocatable :: vec(:)
  end type rvector_type

  type rmatrix_type
    !! The type of a real matrix
    real(dp), allocatable :: mtrx(:,:)
  end type rmatrix_type

  type cmatrix_type
    !! The type of a complex matrix
    complex(dp), allocatable :: mtrx(:,:)
  end type cmatrix_type

  type cd4_type
    !! Centrifugal Distortion parameters for quartric (4) order
    real(dp) :: dn
      !! ΔN  (AKA ΔJ)
    real(dp) :: dnk
      !! ΔNK (AKA ΔJK)
    real(dp) :: dk
      !! ΔK
    real(dp) :: deltan
      !! δn  (AKA δJ)
    real(dp) :: deltak
      !! δK
  end type cd4_type

  type cd6_type
    !! Centrifugal Distortion parameters for sextic (6) order
    real(dp) :: hn
      !! HN   (N²)³
    real(dp) :: hnk
      !! HNK [(N²)² Nz²]
    real(dp) :: hkn
      !! HKN [ N²   Nz⁴]
    real(dp) :: hk
      !! HK         Nz⁶
    real(dp) :: etan
      !! ηN  [N⁴,    (J₊)²+(J₋)²]₊ / 2
    real(dp) :: etank
      !! ηNK [N²Nz², (J₊)²+(J₋)²]₊ / 2
    real(dp) :: etak
      !! ηK  [Nz⁴,   (J₊)²+(J₋)²]₊ / 2
  end type cd6_type

  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓ namelist variables declarations ↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  type config_type
    !! Derived type containing data from the namelist variables

    logical :: add_cd4 = .true.
      !! Add centrifugal distortion for fourth order ?
    logical :: add_cd6 = .true.
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

    integer :: spin_isomer_kind
      !! Whether and how to enforce ortho/para symmetry for molecules with identical nuclei.
      !!   0: don't
      !!   2: C2v rotor (H₂X-like): preserve Ka+Kc parity
      !! Note that this just disables certain transitions from bein calculated
      !! in the CB approx as well as from the S-matrix. This does not affect
      !! the RFT because higher J-blocks of the S-matrix are more affected
      !! by K-mixing (Ka and Kc are not exact quantum numbers)
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
    integer :: target_charge
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
    real(dp) :: abc(3)
      !! Array of reals of length 3
      !! The rotational constants A, B, and C of the target molecule (cm⁻¹).
    real(dp) :: cartesian_dipole_moments(3)
      !! Array of cartesian dipole moments (Debye)
      !! in the order dx, dy, dz
    real(dp) :: cartesian_quadrupole_moments(6)
      !! Array of cartesian quadrupole moments (Debye)
      !! in the order Qxx, Qxy, Qxz, Qyy, Qyz, Qzz
    real(dp), allocatable :: egrid_segs(:)
      !! Array of the bounds (non-degenerate) of the energy grid segments (length num_egrid_segs + 1)

    character(1) :: rotor_kind
      !! The kind of rotor that describes the targer. Character(1).
      !! Choice of :
      !!  "l"inear
      !!  "a"symmetric top
      !!  "s"ymmetric  top
    character(1) :: zaxis
      !! The molecular axis (a, b, or c) along which the z-axis is oriented
      !! **This should also be the symmetry axis**
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

  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓ interfaces ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
  ! -- ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

  interface sort_channels
    module procedure :: sort_elec_channels
  end interface sort_channels

  interface permsort_channels
    module procedure :: permsort_elec_channels
  end interface permsort_channels

  interface operator(.eq.)
    module procedure :: channel_iseq
    module procedure :: transition_iseq
  end interface operator(.eq.)

  interface operator(.ne.)
    module procedure :: channel_isne
  end interface operator(.ne.)

  interface operator(.isin.)
    module procedure channel_isin
    module procedure transition_isin
  end interface operator(.isin.)

  interface assignment(=)
    module procedure :: channel_set_eq
  end interface assignment(=)


! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function get_channel_index(channels, channel, reverse) result(i)
    !! Return the first index i where channel .eq. channels(i) is .true., or
    !! the last index is reverse is .true.
    use rotex__system, only: die
    implicit none
    class(channel_type), intent(in) :: channel, channels(:)
    logical, intent(in), optional :: reverse
    integer :: i, istart, iend, istep
    logical :: reverse_local
    reverse_local = .false. ; if(present(reverse)) reverse_local = reverse
    if(reverse_local .eqv. .true.) then
      istart = lbound(channels, 1)
      iend   = ubound(channels, 1)
      istep  = 1
    else
      istart = ubound(channels, 1)
      iend   = lbound(channels, 1)
      istep  = -1
    endif
    do i = istart, iend, istep
      if(channel .ne. channels(i)) cycle
      return
    enddo
    call die("Failed finding a channel in the channel array")
  end function get_channel_index

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure elemental module subroutine channel_set_eq(channel_out, channel_in)
    !! Sets the channel left equal to the channel right
    use rotex__system, only: die
    implicit none
    class(channel_type), intent(out) :: channel_out
    class(channel_type), intent(in)  :: channel_in

    ! -- start at base class
    channel_out % nelec = channel_in % nelec
    channel_out % E     = channel_in % E

    select type (left => channel_out)

    ! -- return if this is just a base channel
    class default
      return

    ! -- electronic channels
    type is (elec_channel_type)
      select type (right => channel_in)
      type is (elec_channel_type)
        left % l  = right % l
        left % ml = right % ml
        left % iq = right % iq
      class default
        call die("Trying to assign a non-electronic channel to an electronic channel.")
      end select

    ! -- rotational channels + l
    type is (asymtop_rot_channel_l_type)
      select type (right => channel_in)
      type is (asymtop_rot_channel_l_type)
        left % N   = right % N
        left % Ka  = right % Ka
        left % Kc  = right % Kc
        left % l   = right % l
        left % iq  = right % iq
        left % sym = right % sym
      class default
        call die("Trying to assign rotation+l channel to a different kind of channel.")
      end select

    ! -- rotational channels (no l)
    type is (asymtop_rot_channel_type)
      select type (right => channel_in)
      type is (asymtop_rot_channel_type)
        left % N   = right % N
        left % Ka  = right % Ka
        left % Kc  = right % Kc
        left % sym = right % sym
      class default
        call die("Trying to assign rotational channel (no l) to a different kind of channel.")
      end select

    end select
  end subroutine channel_set_eq

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function channel_isin(channel, channels) result(res)
    !! Check if a channel is in the array channels
    implicit none
    class(channel_type), intent(in) :: channel, channels(:)
    logical :: res
    res = .true.
    if(any(channel .eq. channels)) return
    res = .false.
  end function channel_isin
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function transition_isin(transition, transitions) result(res)
    !! Check if a transition is in the array transitions
    implicit none
    type(asymtop_rot_transition_type), intent(in) :: transition, transitions(:)
    logical :: res
    res = .true.
    if(any(transition .eq. transitions)) return
    res = .false.
  end function transition_isin

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure elemental module function trim_channel_l(channel_with_l) result(channel_without_l)
    !! Given a rotational channel with the l quantum number, return the equvalent channel without the l quantum number
    implicit none
    type(asymtop_rot_channel_l_type), intent(in) :: channel_with_l
    type(asymtop_rot_channel_type) :: channel_without_l
    integer :: nelec, N, Ka, Kc, sym
    real(dp) :: E
    nelec = channel_with_l % nelec
    N     = channel_with_l % N
    Ka    = channel_with_l % Ka
    Kc    = channel_with_l % Kc
    sym   = channel_with_l % sym
    E     = channel_with_l % E
    channel_without_l = asymtop_rot_channel_type(nelec = nelec, N = N, Ka = Ka, Kc = Kc, E = E, sym = sym)
  end function trim_channel_l

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function transition_iseq(transition1, transition2) result(res)
    !! Test if two transitions are equal
    implicit none
    type(asymtop_rot_transition_type), intent(in) :: transition1, transition2
    logical :: res
    res = .false.
    if(transition1 % lo .ne. transition2 % lo) return
    if(transition1 % up .ne. transition2 % up) return
    res = .true.
  end function transition_iseq

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function channel_iseq(channel1, channel2) result(res)
    !! Test for channel equality on the basis of their quantum numbers only
    use rotex__system, only: die
    implicit none
    class(channel_type), intent(in) :: channel1, channel2
    logical :: res
    res = .false.
    if(channel1 % nelec .ne. channel2 % nelec)  return
    select type(channel1)
    ! -- electronic channel equality comparison
    type is (elec_channel_type)
      select type(channel2)
      type is (elec_channel_type)
        if(channel1 % l  .ne. channel2 % l)      return
        if(channel1 % ml .ne.  channel2 % ml) return
      class default
        call die("Cannot compare an electronic channel to a different channel")
      end select
    ! -- rotational channel equality comparison (with l)
    type is (asymtop_rot_channel_l_type)
      select type(channel2)
      type is (asymtop_rot_channel_l_type)
        if(channel1 % N  .ne. channel2 % N)  return
        if(channel1 % Ka .ne. channel2 % Ka) return
        if(channel1 % Kc .ne. channel2 % Kc) return
        if(channel1 % l  .ne. channel2 % l)  return
        if(channel1 % sym.ne. channel2 % sym)return
      type is (asymtop_rot_channel_type)
        if(channel1 % N  .ne. channel2 % N)  return
        if(channel1 % Ka .ne. channel2 % Ka) return
        if(channel1 % Kc .ne. channel2 % Kc) return
        if(channel1 % sym.ne. channel2 % sym)return
      class default
        call die("Cannot compare a rotational channel (with l) to a non-rotational channel")
      end select
    ! -- rotational channel equality comparison (no l)
    type is (asymtop_rot_channel_type)
      select type(channel2)
      class is (asymtop_rot_channel_type)
        if(channel1 % N  .ne. channel2 % N)  return
        if(channel1 % Ka .ne. channel2 % Ka) return
        if(channel1 % Kc .ne. channel2 % Kc) return
        if(channel1 % sym.ne. channel2 % sym)return
      class default
        call die("Cannot compare a rotational channel (no l) to a non-rotational channel")
      end select
    end select
    res = .true.
  end function channel_iseq
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function channel_isne(channel1, channel2) result(res)
    implicit none
    class(channel_type), intent(in) :: channel1, channel2
    logical :: res
    res = .not. (channel1 .eq. channel2)
  end function channel_isne

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine permsort_elec_channels(elec_channels, idx)
    !! Sorts the inout array based on the quantum numbers of the channels, and returns
    !! the permutation array that would produce the same output
    implicit none
    type(elec_channel_type), intent(inout) :: elec_channels(:)
    integer, allocatable :: idx(:)
    integer :: ichan, jchan, nchan
    nchan = size(elec_channels, 1)
    idx = [(ichan,ichan=1,nchan)]
    ! -- sort by nelec
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ge. elec_channels(ichan) % nelec) cycle
      call swap_channels(elec_channels, ichan, jchan)
      call swap_ints(idx, ichan, jchan)
    enddo ; enddo
    ! -- sort by l
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ge. elec_channels(ichan) % l) cycle
      call swap_channels(elec_channels, ichan, jchan)
      call swap_ints(idx, ichan, jchan)
    enddo ; enddo
    ! -- sort by λ (ml)
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ne. elec_channels(ichan) % l) cycle
      if(elec_channels(jchan) % ml .ge. elec_channels(ichan) % ml) cycle
      call swap_channels(elec_channels, ichan, jchan)
      call swap_ints(idx, ichan, jchan)
    enddo ; enddo
  end subroutine permsort_elec_channels


  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine sort_elec_channels(elec_channels)
    !! Sorts the inout array based on the quantum numbers of the channels
    implicit none
    type(elec_channel_type), intent(inout) :: elec_channels(:)
    integer :: ichan, jchan, nchan
    nchan = size(elec_channels, 1)
    ! -- sort by nelec
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ge. elec_channels(ichan) % nelec) cycle
      call swap_channels(elec_channels, ichan, jchan)
    enddo ; enddo
    ! -- sort by l
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ge. elec_channels(ichan) % l) cycle
      call swap_channels(elec_channels, ichan, jchan)
    enddo ; enddo
    ! -- sort by λ (ml)
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ne. elec_channels(ichan) % l) cycle
      if(elec_channels(jchan) % ml .ge. elec_channels(ichan) % ml) cycle
      call swap_channels(elec_channels, ichan, jchan)
    enddo ; enddo
  end subroutine sort_elec_channels

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine sort_channels_by_energy(channels)
    !! Bubble sort the array of channels such that the channel energies are in ascending order
    implicit none
    class(channel_type), intent(inout) :: channels(:)
    logical :: swapped
    integer :: nchans
    integer :: j, i
    nchans = size(channels, 1)
    do j=nchans-1, 1, -1
      swapped = .false.
      do i=1,j
        if(channels(i) % E .le. channels(i+1) % E) cycle
        call swap_channels(channels, i, i+1)
        swapped = .true.
      enddo
      if(.not. swapped) exit
    enddo
  end subroutine sort_channels_by_energy

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine swap_ints(arr, ichan, jchan)
    !! Swaps the array elements at arr(ichan) and arr(jchan)
    use rotex__system, only: die
    implicit none
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: ichan, jchan
    integer :: tmp
    if(ichan .lt. lbound(arr, 1) .OR. ichan .gt. ubound(arr, 1)) call die("Trying to swap ints,&
      & but one of the indices exceeds the bounds of the int array")
    if(jchan .lt. lbound(arr, 1) .OR. jchan .gt. ubound(arr, 1)) call die("Trying to swap ints,&
      & but one of the indices exceeds the bounds of the int array")
    tmp = arr(ichan)
    arr(ichan) = arr(jchan)
    arr(jchan) = tmp
  end subroutine swap_ints

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine swap_channels(channels, ichan, jchan)
    !! Swaps the channels at channels(ichan) and channels(jchan)
    use rotex__system, only: die
    implicit none
    class(channel_type), intent(inout) :: channels(:)
    integer, intent(in) :: ichan, jchan
    class(channel_type), allocatable :: tmp
    if(ichan .lt. lbound(channels, 1) .OR. ichan .gt. ubound(channels, 1)) call die("Trying to swap channels,&
      & but one of the indices exceeds the bounds of the channel array")
    if(jchan .lt. lbound(channels, 1) .OR. jchan .gt. ubound(channels, 1)) call die("Trying to swap channels,&
      & but one of the indices exceeds the bounds of the channel array")
    allocate(tmp, source = channels(ichan))
    tmp = channels(ichan)
    channels(ichan) = channels(jchan)
    channels(jchan) = tmp
    deallocate(tmp)
  end subroutine swap_channels

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function findloc_transitions(targs, search) result(idxtarg)
    !! Find the indices for each element in targs that map to the elements in search.
    !! Return 0 if there is no such mapping.
    implicit none
    type(asymtop_rot_transition_type), intent(in) :: targs(:), search(:)
    integer, allocatable :: idxtarg(:)
      !! TARGS -> SEARCH mapping
    integer, parameter :: IDX_NOT_FOUND = 0
    integer :: isearch, itarg, ntargs, nsearch
    integer, allocatable :: idxsearch(:)
    logical, allocatable :: mask(:)
    ntargs  = size(targs, 1)
    nsearch = size(search, 1)
    allocate(idxtarg(ntargs), source = IDX_NOT_FOUND)
    do itarg=1, ntargs
      mask = search .eq. targs(itarg)
      if(all(mask .eqv. .false.)) cycle ! cycle if no match
      idxsearch = pack([(isearch, isearch=1, nsearch)], mask)
      idxtarg(itarg) = idxsearch(1) ! take first match
    enddo
  end function findloc_transitions

  ! ! ------------------------------------------------------------------------------------------------------------------------------- !
  ! pure module subroutine swap_elec_channel_values(elec_channels, ichan, jchan)
  !   !! Swaps the channels at elec_channels(ichan) and elec_channels(jchan)
  !   use rotex__system, only: die
  !   implicit none
  !   type(elec_channel_type), intent(inout) :: elec_channels(:)
  !   integer, intent(in) :: ichan, jchan
  !   type(elec_channel_type) :: tmp
  !   if(ichan .lt. lbound(elec_channels, 1) .OR. ichan .gt. ubound(elec_channels, 1)) call die("Trying to swap channels,&
  !     & but one of the indices exceeds the bounds of the channel array")
  !   if(jchan .lt. lbound(elec_channels, 1) .OR. jchan .gt. ubound(elec_channels, 1)) call die("Trying to swap channels,&
  !     & but one of the indices exceeds the bounds of the channel array")
  !   tmp = elec_channels(ichan)
  !   elec_channels(ichan) = elec_channels(jchan)
  !   elec_channels(jchan) = tmp
  ! end subroutine swap_elec_channel_values
  ! ! ------------------------------------------------------------------------------------------------------------------------------- !
  ! pure module subroutine swap_asymtop_rot_channel_values(asymtop_rot_channels, ichan, jchan)
  !   !! Swaps the channels at asymtop_rot_channels(ichan) and asymtop_rot_channels(jchan)
  !   use rotex__system, only: die
  !   implicit none
  !   type(asymtop_rot_channel_type), intent(inout) :: asymtop_rot_channels(:)
  !   integer, intent(in) :: ichan, jchan
  !   type(asymtop_rot_channel_type) :: tmp
  !   if(ichan .lt. lbound(asymtop_rot_channels, 1) .OR. ichan .gt. ubound(asymtop_rot_channels, 1)) &
  !     call die("Trying to swap channels, but one of the indices exceeds the bounds of the channel array")
  !   if(jchan .lt. lbound(asymtop_rot_channels, 1) .OR. jchan .gt. ubound(asymtop_rot_channels, 1)) &
  !     call die("Trying to swap channels, but one of the indices exceeds the bounds of the channel array")
  !   tmp = asymtop_rot_channels(ichan)
  !   asymtop_rot_channels(ichan) = asymtop_rot_channels(jchan)
  !   asymtop_rot_channels(jchan) = tmp
  ! end subroutine swap_asymtop_rot_channel_values

! ================================================================================================================================ !
end module rotex__types
! ================================================================================================================================ !
