! ================================================================================================================================ !
module rotex__frames
  !! Procedures involving reference frames
  use rotex__kinds,      only: dp
  use rotex__types,      only: xyz_type
  use rotex__system,     only: stderr, stdout, die
  use rotex__characters, only: lower

  private

  public :: get_euler_angles
  public :: xyz_is_valid
  public :: xyz_from_z

  public :: operator(.ne.)
  public :: operator(.eq.)

  interface operator(.eq.)
    module procedure :: xyz_is_eq
  end interface operator(.eq.)

  interface operator(.ne.)
    module procedure :: xyz_is_ne
  end interface operator(.ne.)

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  subroutine get_euler_angles(xyz_from, xyz_to, a, b, g)
    !! Calculates the euler angles to rotate between two ABC frames defined by xyz_from and xyz_to
    use rotex__types,      only: xyz_type
    use rotex__characters, only: lower
    use rotex__constants,  only: pi
    implicit none(type, external)
    type(xyz_type), intent(in) :: xyz_from, xyz_to
      !! Right-handed axes: ABC, BCA, or CAB
    real(dp), intent(out) :: a, b, g
      !! Euler angles α, β, γ
    character(8) :: s

    if(xyz_is_valid(xyz_from) .eqv. .false.) then
      write(stderr, '("XYZ_from: ", 3A)') xyz_from
      call die("Invalid xyz_from")
    endif
    if(xyz_is_valid(xyz_to) .eqv. .false.) then
      write(stderr, '("XYZ_to: ", 3A)') xyz_to
      call die("Invalid xyz_to")
    endif

    a = 0.0_dp ; b = 0.0_dp ; g = 0.0_dp

    if(xyz_from .eq. xyz_to) return

    associate( xfrom => xyz_from % x, yfrom => xyz_from % y, zfrom => xyz_from % z &
             , xto   => xyz_to   % x, yto   => xyz_to   % y, zto   => xyz_to   % z)

      s = lower(xfrom//yfrom//zfrom//"->"//xto//yto//zto)

      ! -- active rotation of ABC inside of xyz
      select case(s)
      case ("abc->cab", "cab->bca", "bca->abc")
        a = 0.0_dp
        b = pi/2.0_dp
        g = pi/2.0_dp
      case ("abc->bca", "bca->cab", "cab->abc")
        a = pi/2.0_dp
        b = pi/2.0_dp
        g = 0.0_dp
      case default
        write(stderr, '("Somehow, the string S is not properly formatted as `XYZ->XYZ`: ", A)') s
        write(stderr, '("XYZ_FROM % X`: ", A)') xyz_from % x
        write(stderr, '("XYZ_FROM % Y`: ", A)') xyz_from % y
        write(stderr, '("XYZ_FROM % Z`: ", A)') xyz_from % z
        write(stderr, '("XYZ_TO   % X`: ", A)') xyz_to % x
        write(stderr, '("XYZ_TO   % Y`: ", A)') xyz_to % y
        write(stderr, '("XYZ_TO   % Z`: ", A)') xyz_to % z
        call die("Couldn't determine reference frame(s)")
      end select

    end associate

  end subroutine get_euler_angles

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental module function xyz_is_valid(xyz) result(res)
    !! Test if xyz is a valid right-handed coordinate system
    implicit none(type, external)
    type(xyz_type), intent(in) :: xyz
    logical :: res
    character(1) :: x,y,z
    character(1) :: abc_arr(3) = ["a", "b", "c"]
    character(3) :: abc

    x = lower(xyz%x)
    y = lower(xyz%y)
    z = lower(xyz%z)

    ! -- non-ABC axes ?
    res = any(x .eq. abc_arr) .AND. any(y .eq. abc_arr) .AND. any(z .eq. abc_arr)
    if(res .eqv. .false.) then
      write(stderr, '("xyz % x: ", A)') x
      write(stderr, '("xyz % y: ", A)') y
      write(stderr, '("xyz % z: ", A)') z
      call die("xyz is not abc: non-ABC axis/axes detected")
    endif

    ! -- redundant axes ?
    if(x.eq.y .OR. x.eq.z .OR. y.eq.z) then
      write(stderr, '("xyz % x: ", A)') x
      write(stderr, '("xyz % y: ", A)') y
      write(stderr, '("xyz % z: ", A)') z
      call die("xyz is not abc: redundant axes")
    endif

    ! -- right handed ?
    abc = x//y//z
    select case(abc)
    case("abc", "bca", "cab")
      res = .true.
    case default
      res = .false.
    end select

  end function xyz_is_valid

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function xyz_from_z(z) result(xyz)
    !! Constructs a valid right-handed xyz frame in terms of interial axes ABC
    !! given only the z-axis
    use rotex__characters, only: lower
    implicit none(type, external)
    character(1), intent(in) :: z
      !! The z-axis
    type(xyz_type) :: xyz
    select case(lower(z))
    case("a")
      xyz = xyz_type(x = "b", y = "c", z = "a")
    case("b")
      xyz = xyz_type(x = "c", y = "a", z = "b")
    case("c")
      xyz = xyz_type(x = "a", y = "b", z = "c")
    case default
      call die("Z-axis '"//z//"' is not one of ABC")
    end select
  end function xyz_from_z

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function xyz_is_eq(xyz1, xyz2) result(res)
    implicit none(type, external)
    type(xyz_type), intent(in) :: xyz1, xyz2
    logical :: res
    res = .false.
    if(xyz1%x .ne. xyz2%x) return
    if(xyz1%y .ne. xyz2%y) return
    if(xyz1%z .ne. xyz2%z) return
    res = .true.
  end function xyz_is_eq

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function xyz_is_ne(xyz1, xyz2) result(res)
    implicit none(type, external)
    type(xyz_type), intent(in) :: xyz1, xyz2
    logical :: res
    res = .not. xyz_is_eq(xyz1, xyz2)
  end function xyz_is_ne

! ================================================================================================================================ !
end module rotex__frames
! ================================================================================================================================ !
