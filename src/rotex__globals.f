! ================================================================================================================================ !
module rotex__globals
  !! Putting some important global variables in here so we don't have to pass them everywhere
  !! all the time

  implicit none

  private

  save

  public :: copy2globals

  character(1), public :: GLOBAL_ROTOR_KIND
    !! "a"symmetric top, "s"ymmetric top, or "l"inear rotor
  character(1), public :: GLOBAL_C2AXIS
    !! The C₂ axis: "A", "B", "C"
  character(1), public :: GLOBAL_ROTOR_ZAXIS
    !! The z axis for determining projections of N: "A", "B", "C"
  integer, public :: GLOBAL_SPIN_ISOMER_KIND
    !! The kind of spin isomer to enforce. See symmetry module for details
  integer, public :: GLOBAL_FORBIDDEN_STATES_KIND
    !! Special case exclusion of certain rotational levels based on the
    !! symmetry of the ground-state vibrational wavefunction.
    !!   0: none
    !!   1: even N, K=0 forbidden (e.g., H₃⁺)

  ! character(1), parameter :: GLOBAL_ROTOR_KIND_VALS(*)   = ["a", "s", "l"]
  ! character(1), parameter :: GLOBAL_ROTOR_ZAXIS_VALS(*)  = ["a", "b", "c"]
  ! character(1), parameter :: GLOBAL_ROTOR_C2AXIS_VALS(*) = ["a", "b", "c"]

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine copy2globals( rotor_kind       &
                                , rotor_zaxis      &
                                , c2axis           &
                                , spin_isomer_kind &
                                , forbidden_states_kind &
    )
    !! Copies values from CFG into this module's globals
    use rotex__types,      only: config_type
    use rotex__characters, only: lower
    implicit none
    character(1), intent(in) :: rotor_kind
    character(1), intent(in) :: rotor_zaxis
    character(1), intent(in) :: c2axis
    integer,      intent(in) :: spin_isomer_kind
    integer,      intent(in) :: forbidden_states_kind

    GLOBAL_ROTOR_KIND            = lower(rotor_kind)
    GLOBAL_ROTOR_ZAXIS           = lower(rotor_zaxis)
    GLOBAL_C2AXIS                = lower(c2axis)
    GLOBAL_SPIN_ISOMER_KIND      = spin_isomer_kind
    GLOBAL_FORBIDDEN_STATES_KIND = forbidden_states_kind

  end subroutine copy2globals

! ================================================================================================================================ !
end module rotex__globals
! ================================================================================================================================ !
