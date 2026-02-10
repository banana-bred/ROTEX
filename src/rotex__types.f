! ================================================================================================================================ !
module rotex__types
  !! Contains type definitions and procedures for those types used throughout the program
  use rotex__kinds,     only: dp, xdp, qp
  use rotex__constants, only: IQ_DEFAULT

  implicit none (type, external)

  private

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
  public :: cd4_type
  public :: cd6_type

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
    type(eigenH_type) :: eigenH
      !! The decomposed Hamiltonian for this rotational level
    real(dp), allocatable :: EinstA(:)
      !! The Einstein coefficients for transitions to all lower states (0 if none)
    integer, allocatable :: Ka(:)
      !! The projections Ka
    integer, allocatable :: Kc(:)
      !! The projections Kc
    integer :: N
      !! Rotational quantum number
  end type N_states_type

  type, abstract :: channel_type
    !! |nelec> (E)
    real(dp) :: E
      !! The channel energy
    integer :: nelec
      !! Electronic state
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

! ================================================================================================================================ !
end module rotex__types
! ================================================================================================================================ !
