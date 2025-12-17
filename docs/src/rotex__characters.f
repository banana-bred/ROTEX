! =================================================================================================================================!
module rotex__characters
  !! Contains procedures related to characters and character arrays, such as converting integers to characters

  use rotex__kinds,     only: dp
  use rotex__functions, only: iseven

  implicit none

  private

  ! -- procedures
  public :: ndigits
  public :: s2hms
  public :: int2char
  public :: dJ2char
  public :: add_trailing
  public :: upper
  public :: to_upper
  public :: lower
  public :: to_lower
  public :: sup
  public :: sub


  integer, parameter :: big_char = 100

  interface s2hms
    !! convert seconds to hours, minutes, seconds
    module procedure int_s2hms
    module procedure real_s2hms
  end interface

  interface int2char
    module procedure :: scalar_int2char
    module procedure :: vector_int2char
  end interface int2char

! =================================================================================================================================!
contains
! =================================================================================================================================!

  ! ---------------------------------------------------------------------------------------------------------------------------------!
  pure module function int_s2hms(s) result(time)
    !! Given an integer 's' in seconds, convert to the format hh:mm:ss.

    implicit none

    integer, intent(in) :: s
    character(:), allocatable :: time

    character(big_char / 10) :: tmp
    integer :: hh,mm,ss

    select case(s)
      case(:-1)     ; time = 'negative time' ; return
      case(0:59)    ; hh = 0      ; mm = 0              ; ss = s                ; write(tmp,'(I2,"h ",I2,"m ",I2,"s")') hh,mm,ss
      case(60:3599) ; hh = 0      ; mm = s/60           ; ss = s - mm*60        ; write(tmp,'(I2,"h ",I2,"m ",I2,"s")') hh,mm,ss
      case(3600:)   ; hh = s/3600 ; mm = (s-hh*3600)/60 ; ss = s-hh*3600-mm*60  ; write(tmp,'(I0,"h ",I2,"m ",I2,"s")') hh,mm,ss
    end select

    time = trim(tmp)

  end function int_s2hms

  ! ---------------------------------------------------------------------------------------------------------------------------------!
  pure module function real_s2hms(s_re) result(time)
    !! Given an integer in seconds, convert to the format hh:mm:ss. Input is a real, gets converted to int
    implicit none

    real(dp), intent(in) :: s_re
    character(:), allocatable :: time

    character(big_char) :: tmp
    integer :: hh,mm,ss
    integer :: s

    s = int(s_re)

    select case(s)
    case(:-1)     ; time = 'negative time' ; return
    case(0:59)    ; hh = 0      ; mm = 0              ; ss = s
    case(60:3599) ; hh = 0      ; mm = s/60           ; ss = s - mm*60
    case(3600:)   ; hh = s/3600 ; mm = (s-hh*3600)/60 ; ss = s-hh*3600-mm*60
    end select

    write(tmp,'(I0,"h ",I0,"m ",I0,"s")') hh,mm,ss

    time = trim(tmp)

  end function real_s2hms

  ! ---------------------------------------------------------------------------------------------------------------------------------!
  pure elemental function ndigits(n) result(num)
    !! Returns number of characters an integer will occupy
    use rotex__constants, only: one
    implicit none
    integer, intent(in) :: n
    integer :: num
    num = 1
    if(n .eq. 0) return
    num = floor(log10(abs(n) * one)) + 1
    ! -- account for minus sign
    if(n.lt.1) num = num + 1
  end function ndigits

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function scalar_int2char(i) result(res)
    !! Writes the value i to a character as I0
    implicit none
    integer, intent(in) :: i
    character(:), allocatable :: res
    allocate(character(ndigits(i)) :: res)
    write(res, '(I0)') i
  end function scalar_int2char

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function vector_int2char(i) result(res)
    !! Writes the value i to a character as I0
    implicit none
    integer, intent(in) :: i(:)
    character(:), allocatable :: res
    integer :: n
    if(size(i,1) .eq. 0) then
      res = "()"
      return
    endif
    ! -- n-1 commas, and all the digits
    n = (size(i,1) - 1) + sum(ndigits(i))
    allocate(character(n) :: res)
    write(res, '(*(I0,:,","))') i
    res = "(" // res // ")"
  end function vector_int2char

  ! ---------------------------------------------------------------------------------------------------------------------------------!
  pure module subroutine add_trailing(chr, trail)
    !! Add a trailing character `trail` to the character `chr` if it is not already the
    !! last character

    implicit none

    character(:), allocatable,  intent(inout) :: chr
    character(*), intent(in) :: trail

    integer :: n
    integer :: m

    n = len(chr) - len(trail) + 1
    m = len(chr)

    if(n < 1) return

    if(chr(n:m) .eq. trail) return

    chr = chr // trail

  end subroutine add_trailing

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module function dJ2char(dJ) result(res)
    !! Takes an integer dJ and results the character representing half of it.
    !! dJ2char(2) -> "1"
    !! dJ2char(3) -> "3/2"

    implicit none

    integer, intent(in) :: dJ
      !! Twice the angular momentum
    character(:), allocatable :: res
      !! The output character representation

    if(iseven(dJ) .eqv. .true.) then
      res = int2char(dJ/2)
    else
      res = int2char(dJ) // "|2"
    endif

  end function dJ2char

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function lower(chr) result(res)
    !! returns a lower case character
    implicit none
    character(*), intent(in) :: chr
    character(:), allocatable :: res
    integer, parameter :: shift = ichar('a') - ichar("A")
    integer, parameter :: uppercase_a = ichar('A')
    integer, parameter :: uppercase_z = ichar('Z')
    integer :: i, n, ic
    n = len(chr)
    res = chr
    do i = 1, n
      ic = ichar(res(i:i))
      ! -- cycle if the character isn't in [A,Z]
      if(ic .lt. uppercase_a) cycle
      if(ic .gt. uppercase_z) cycle
      res(i:i) = char(ic + shift)
    enddo
  end function lower
  pure elemental subroutine to_lower(chr)
    !! converts a character to lower case
    implicit none
    character(*), intent(inout) :: chr
    integer, parameter :: shift = ichar('a') - ichar("A")
    integer, parameter :: uppercase_a = ichar('A')
    integer, parameter :: uppercase_z = ichar('Z')
    integer :: i, n, ic
    n = len(chr)
    do i = 1, n
      ic = ichar(chr(i:i))
      ! -- cycle if the character isn't in [A,Z]
      if(ic .lt. uppercase_a) cycle
      if(ic .gt. uppercase_z) cycle
      chr(i:i) = char(ic + shift)
    enddo
  end subroutine to_lower
  pure function upper(chr) result(res)
    !! returns an upper case character
    implicit none
    character(*), intent(in) :: chr
    character(:), allocatable :: res
    integer, parameter :: shift = ichar('a') - ichar("A")
    integer, parameter :: lowercase_a = ichar('a')
    integer, parameter :: lowercase_z = ichar('z')
    integer :: i, n, ic
    n = len(chr)
    res = chr
    do i = 1, n
      ic = ichar(res(i:i))
      ! -- cycle if the character isn't in [A,Z]
      if(ic .lt. lowercase_a) cycle
      if(ic .gt. lowercase_z) cycle
      res(i:i) = char(ic - shift)
    enddo
  end function upper
  pure elemental subroutine to_upper(chr)
    !! converts a character to upper case
    implicit none
    character(*), intent(inout) :: chr
    integer, parameter :: shift = ichar('a') - ichar("A")
    integer, parameter :: lowercase_a = ichar('a')
    integer, parameter :: lowercase_z = ichar('z')
    integer :: i, n, ic
    n = len(chr)
    do i = 1, n
      ic = ichar(chr(i:i))
      ! -- cycle if the character isn't in [A,Z]
      if(ic .lt. lowercase_a) cycle
      if(ic .gt. lowercase_z) cycle
      chr(i:i) = char(ic - shift)
    enddo
  end subroutine to_upper

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function sub(x) result(res)
    !! Returns the subscript version of the integer x
    implicit none
    integer, intent(in) :: x
    character(:), allocatable :: res
    integer :: i
    character(:), allocatable :: xchar
    allocate(character(ndigits(x)) :: res)
    allocate(character(ndigits(x)) :: xchar)
    xchar(:) = int2char(x)
    do i=1, len(xchar)
      select case(xchar(i:i))
        case("-") ; res(i:i) = "_"
        case("0") ; res(i:i) = "₀"
        case("1") ; res(i:i) = "₁"
        case("2") ; res(i:i) = "₂"
        case("3") ; res(i:i) = "₃"
        case("4") ; res(i:i) = "₄"
        case("5") ; res(i:i) = "₅"
        case("6") ; res(i:i) = "₆"
        case("7") ; res(i:i) = "₇"
        case("8") ; res(i:i) = "₈"
        case("9") ; res(i:i) = "₉"
      end select
    enddo
  end function sub

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function sup(x) result(res)
    !! Returns the superscript version of the integer x
    implicit none
    integer, intent(in) :: x
    character(:), allocatable :: res
    integer :: i
    character(:), allocatable :: xchar
    allocate(character(ndigits(x)) :: res)
    allocate(character(ndigits(x)) :: xchar)
    xchar(:) = int2char(x)
    do i=1, len(xchar)
      select case(xchar(i:i))
        case("-") ; res(i:i) = "⁻"
        case("0") ; res(i:i) = "⁰"
        case("1") ; res(i:i) = "¹"
        case("2") ; res(i:i) = "²"
        case("3") ; res(i:i) = "³"
        case("4") ; res(i:i) = "⁴"
        case("5") ; res(i:i) = "⁵"
        case("6") ; res(i:i) = "⁶"
        case("7") ; res(i:i) = "⁷"
        case("8") ; res(i:i) = "⁸"
        case("9") ; res(i:i) = "⁹"
      end select
    enddo
  end function sup


! =================================================================================================================================!
end module rotex__characters
! =================================================================================================================================!
