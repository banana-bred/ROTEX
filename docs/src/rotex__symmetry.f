! ================================================================================================================================ !
module rotex__symmetry
  !! All things related to symmetry

  implicit none

  private

  public :: irrep_name
  public :: get_group_irreps
  public :: group_size
  public :: spin_symmetry
  public :: possible_spin_symmetries
  public :: is_spin_allowed
  public :: is_spin_forbidden

  interface is_spin_allowed
    module procedure :: is_spin_allowed_chan
    module procedure :: is_spin_allowed_qnums
  end interface is_spin_allowed

  interface is_spin_forbidden
    module procedure :: is_spin_forbidden_chan
    module procedure :: is_spin_forbidden_qnums
  end interface is_spin_forbidden

  integer, allocatable, public, save :: m_parity(:)
    !! Array containing the parity (1/even or -1/odd) of an electronic channel
    !! based on its label m. This array is indexed by m directly. This is only
    !! for calculation in the Cs point group

  character(33), parameter :: abelian_point_groups = "C1, Cs, C2, Ci, C2v, C2h, D2, D2h"

  ! ---------------------------------------- !
  ! The integer labels of the various irreps !
  ! ---------------------------------------- !
  ! -- Cs
  integer, parameter, public :: Ap  = 1
  integer, parameter, public :: App = 2
  ! -- Ci, C2h
  integer, parameter, public :: Ag = 1
  integer, parameter, public :: Au = 2
  integer, parameter, public :: Bg = 3
  integer, parameter, public :: Bu = 4
  ! -- C1, C2, C2v, D2, D2h
  integer, parameter, public :: A  = 1
  integer, parameter, public :: A1 = 1
  integer, parameter, public :: A2 = 4
  integer, parameter, public :: B  = 2
  integer, parameter, public :: B1 = 2
  integer, parameter, public :: B2 = 3
  integer, parameter, public :: B3 = 4
  integer, parameter, public :: B1g = 3
  integer, parameter, public :: B1u = 4
  integer, parameter, public :: B2g = 5
  integer, parameter, public :: B2u = 6
  integer, parameter, public :: B3g = 7
  integer, parameter, public :: B3u = 8

  integer, parameter, public :: even = 1
  integer, parameter, public :: odd  =-1

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function group_size(point_group) result(n)
    !! Return the number of elements in point_group
    use rotex__system,     only: die
    use rotex__characters, only: to_upper
    implicit none
    character(*), intent(in) :: point_group
    integer :: n
    character(:), allocatable :: pg
    pg = trim(point_group)
    call to_upper(pg)
    select case(pg)
      case("C1")
        n = 1
      case("CS", "C2", "CI")
        n = 2
      case("C2V", "C2H", "D2")
        n = 4
      case("D2H")
        n = 8
      case default
        call die("Bad point group (" // pg // ") supplied. Please choose one of " // abelian_point_groups)
    end select
  end function group_size

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine get_group_irreps(point_group, irreps)
    !! Given the point group, output an array containing the names of the irreps in the supplied point_group.
    !! Only Abelian point groups are considered. Irreps in the code will be referred to by their indicies

    use rotex__system,     only: die
    use rotex__characters, only: to_upper

    implicit none

    character(*), intent(in) :: point_group
    character(:), intent(out), allocatable :: irreps(:)

    integer :: nirreps
    character(:), allocatable :: pg

    pg = trim(point_group)
    call to_upper(pg)
    nirreps = group_size(pg)

    select case(pg)

      case("C1")
        allocate(character(1) :: irreps(nirreps))
        irreps(A) = "A"

      case("CS")
        allocate(character(3) :: irreps(nirreps))
        irreps(Ap)  = "Ap"
        irreps(App) = "App"

      case("C2")
        allocate(character(1) :: irreps(nirreps))
        irreps(A) = "A"
        irreps(B) = "B"

      case("CI")
        allocate(character(2) :: irreps(nirreps))
        irreps(Ag) = "Ag"
        irreps(Au) = "Au"

      case("C2V")
        allocate(character(2) :: irreps(nirreps))
        irreps(A1) =  "A1"
        irreps(B1) =  "B1"
        irreps(B2) =  "B2"
        irreps(A2) =  "A2"

      case("C2H")
        allocate(character(2) :: irreps(nirreps))
        irreps(A)  = "Ag"
        irreps(B1) = "Au"
        irreps(B2) = "Bg"
        irreps(B3) = "Bu"

      case("D2")
        allocate(character(2) :: irreps(nirreps))
        irreps(A)  = "A"
        irreps(B1) = "B1"
        irreps(B2) = "B2"
        irreps(B3) = "B3"

      case("D2H")
        allocate(character(3) :: irreps(nirreps))
        irreps(Ag)  = "Ag"
        irreps(Au)  = "Au"
        irreps(B1g) = "B1g"
        irreps(B1u) = "B1u"
        irreps(B2g) = "B2g"
        irreps(B2u) = "B2u"
        irreps(B3g) = "B3g"
        irreps(B3u) = "B3u"

      case default
        call die("Unacceptable point group '" // point_group // "' given. Please choose one of " // abelian_point_groups // ".")

    end select

  end subroutine get_group_irreps

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function irrep_name(irrep, point_group) result(output)
    !! Given an irrep index in point_group, return the name of the corresponding irrep

    use rotex__system,     only: die
    use rotex__characters, only: to_upper

    implicit none

    integer, intent(in) :: irrep
    character(*), intent(in) :: point_group
    character(:), allocatable :: output

    character(:), allocatable :: pg

    character(:), allocatable :: irreps(:)

    pg     = trim(point_group)
    call to_upper(pg)
    call get_group_irreps(pg, irreps)
    ! irreps = group_irreps(pg)
    output = trim(irreps(irrep))

  end function irrep_name

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function possible_spin_symmetries(kind) result(res)
    !! Returns an array of possible spin symmetry values
    use rotex__system, only: die
    implicit none
    integer, intent(in) :: kind
    integer, allocatable :: res(:)
    select case(kind)
    case(0)   ; res = [0]
    ! case(1,2) ; res = [0,1]
    case(2) ; res = [0,1]
    ! case(3)   ; res = [0,1]
    case default
      ! call die("Symmetry kind not supported. Must be one of 0,1,2,3")
      call die("Symmetry kind not supported. Must be one of 0,2")
    end select
  end function possible_spin_symmetries

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function spin_symmetry(n, ka, kc, kind, symaxis) result(res)
    !! Returns the spin symmetry of the current N, Ka, Kc state
    use rotex__system, only: die
    implicit none
    integer,      intent(in) :: n, ka, kc, kind
    character(1), intent(in) :: symaxis
    integer :: ksym
    integer :: res
    select case(kind)
    case(0)
      ! -- no restriction
      res = 0
    case(1)
      ! -- linear: N-parity
      res = iand(n, 1)
    case(2)
      ! -- Ka+Kc parity
      !    Ortho: Ka+Kc odd
      !    Para:  Ka+Kc even
      res = iand(Ka+Kc, 1)
    case(3)
      ! -- must preserve Ks partity (mod 3) (s=a,c)
      !    Ortho: Ks (mod 3) = 0
      !    Para:  Ks (mod 3) ≠ 0
      select case (symaxis)
      case("a", "A")
        ksym = abs(ka)
        call die("Symaxis is B for symmetry rule 3. Are you sure ?")
      case("b", "B")
        call die("Symaxis is B for symmetry rule 3. Are you sure ?")
      case("c", "C")
        ksym = abs(kc)
      case default
        call die("Undetermined SYMAXIS: "//symaxis)
      end select
      res = merge(0, 1, mod(ksym, 3) .eq. 0)
    case default
      ! call die("Illegal symmetry rule. Must be 0,1,2,3.")
      call die("Illegal symmetry rule. Must be 0,2.")
    end select
  end function spin_symmetry

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function is_spin_allowed_qnums(nlo, kalo, kclo, nup, kaup, kcup, kind, symaxis) result(res)
    !! Determine if the transition Nlo,Kalo,Kclo -> Nup,Kaup,Kcup is allowed by nuclear spin symmetry
    !! selection rules
    implicit none
    integer,      intent(in) :: nlo, kalo, kclo, nup, kaup, kcup, kind
    character(1), intent(in) :: symaxis
    logical :: res
    res = spin_symmetry(nlo, kalo, kclo, kind, symaxis) .eq. spin_symmetry(nup, kaup, kcup, kind, symaxis)
  end function is_spin_allowed_qnums
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function is_spin_allowed_chan(channel1, channel2, spin_isomer_kind, symaxis) result(res)
    !! Test if two rotational channels respect ortho/para symmetry
    use rotex__types, only: asymtop_rot_channel_type
    implicit none
    type(asymtop_rot_channel_type), intent(in) :: channel1, channel2
    integer, intent(in) :: spin_isomer_kind
    character(1), intent(in) :: symaxis
    logical :: res
    integer :: n1, ka1, kc1
    integer :: n2, ka2, kc2
    n1  = channel1%n ; ka1 = channel1%ka ; kc1 = channel1%kc
    n2  = channel2%n ; ka2 = channel2%ka ; kc2 = channel2%kc
    res = is_spin_allowed_qnums(n1, ka1, kc1, n2, ka2, kc2, spin_isomer_kind, symaxis)
  end function is_spin_allowed_chan

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function is_spin_forbidden_qnums(nlo, kalo, kclo, nup, kaup, kcup, kind, symaxis) result(res)
    !! Determine if the transition Nlo,Kalo,Kclo -> Nup,Kaup,Kcup is forbidden by nuclear spin symmetry
    !! selection rules
    implicit none
    integer,      intent(in) :: nlo, kalo, kclo, nup, kaup, kcup, kind
    character(1), intent(in) :: symaxis
    logical :: res
    res = .not. is_spin_allowed_qnums(nlo, kalo, kclo, nup, kaup, kcup, kind, symaxis)
  end function is_spin_forbidden_qnums
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function is_spin_forbidden_chan(channel1, channel2, spin_isomer_kind, symaxis) result(res)
    !! Test if two rotational channels respect ortho/para symmetry
    use rotex__types, only: asymtop_rot_channel_type
    implicit none
    type(asymtop_rot_channel_type), intent(in) :: channel1, channel2
    integer, intent(in) :: spin_isomer_kind
    character(1), intent(in) :: symaxis
    logical :: res
    integer :: n1, ka1, kc1
    integer :: n2, ka2, kc2
    n1  = channel1%n ; ka1 = channel1%ka ; kc1 = channel1%kc
    n2  = channel2%n ; ka2 = channel2%ka ; kc2 = channel2%kc
    res = .not. is_spin_allowed_qnums(n1, ka1, kc1, n2, ka2, kc2, spin_isomer_kind, symaxis)
  end function is_spin_forbidden_chan

! ================================================================================================================================ !
end module rotex__symmetry
! ================================================================================================================================ !
