! ================================================================================================================================ !
module rotex__system
  !! Contains the definitions of stdout, stdin, stderr, and procedures to interact with the program/system
  !! such as producing warnings and stopping the execution of the code while producing error messages

  use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, error_unit, iostat_end

  implicit none

  private

  save

  ! -- procedures
  public :: die
  public :: warn
  public :: error
  public :: determine_system_properties
  public :: mkdir


                      public :: OS_NAME
  integer, parameter, public :: OS_ALL     = -1   ! "all" flag for profile support
  integer, parameter, public :: OS_UNKNOWN = 0
  integer, parameter, public :: OS_LINUX   = 1
  integer, parameter, public :: OS_MACOS   = 2
  integer, parameter, public :: OS_WINDOWS = 3
  integer, parameter, public :: OS_CYGWIN  = 4
  integer, parameter, public :: OS_SOLARIS = 5
  integer, parameter, public :: OS_FREEBSD = 6
  integer, parameter, public :: OS_OPENBSD = 7

  logical, public :: OS_is_windows
    !! Is the current operating system Windows ?

  integer, parameter, public :: STDIN  = input_unit
    !! The file unit associated with standard input
  integer, parameter, public :: STDOUT = output_unit
    !! The file unit associated with standard output
  integer, parameter, public :: STDERR = error_unit
    !! The file unit associated with standard error

  integer, parameter, public :: IOSTAT_OK = 0
    !! The expected iostat result from a successful call to read()
  integer, public :: shell_ok
    !! The expected return value for the current environment and shell. Used in system calls.

  character(5), parameter, public :: PROGNAME = "ROTEX"
    !! The program name

  character(1), public :: DIRECTORY_SEPARATOR
    !! The OS-dependent directory separator

  character(:), allocatable:: mkdir_command
    !! The OS-dependent command used to make directories

  interface die
    module procedure :: die_1
    ! module procedure :: die_2
  end interface die

  interface warn
    module procedure :: warn_1
    ! module procedure :: warn_2
  end interface warn

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine error(message)
    !! Print error messages to the screen without the WARNING prompt. This will typically precede a call to DIE
    character(*), intent(in) :: message
    write(stderr,*)
    write(stderr,'("ERROR :: ", A)') message
    write(stderr,*)
  end subroutine error

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine die_1(message)
    !! Stop program execution with a message
    implicit none
    character(*), intent(in), optional :: message
    if(.not.present(message)) error stop ; error stop message
  end subroutine die_1

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine warn_1(message)
    !! Print a warning message, but don't stop the program's execution
    implicit none
    character(*), intent(in) :: message
    write(stderr,*)
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,*)
    write(stderr,'("WARN",X,"::",X,A)') message
    write(stderr,*)
    write(stderr,'("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")')
    write(stderr,*)
  end subroutine warn_1

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine determine_system_properties()
    !! Detects the type of the operating system. As far as system calls and directory structure go,
    !! this basically resolved to Windows or not Windows.

    implicit none

    integer :: OS

    OS = get_os_type()

    select case(OS)
    case(OS_LINUX, OS_MACOS, OS_CYGWIN, OS_SOLARIS, OS_FREEBSD, OS_OPENBSD)

      OS_is_windows = .false.

    case(OS_UNKNOWN)

      OS_is_windows = .false.
      call warn("Operating system unknown. Assuming it is of type unix.")

    case(OS_WINDOWS)

      OS_is_windows = .true.

    case default

      OS_is_windows = .false.
      call warn("Unable to detect the fact that the operating system is unknown. Assuming it is of type unix.")

    end select

    call system("", status = shell_ok)

    write(stdout, '(A)') "Detected operating system type :: " // OS_NAME(OS)
    write(stdout, *)

    if(OS_is_windows) then
      directory_separator = "\"
      mkdir_command = "md "
    else
      directory_separator = "/"
      mkdir_command = "mkdir -p "
    endif

  end subroutine determine_system_properties

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine mkdir(directory)
    !! Makes the directory "directory" and checks that it exists and is writeable

    implicit none

    character(*) :: directory

    integer :: stat
    character(:), allocatable :: cstat

    call execute_command_Line(mkdir_command // directory, exitstat = stat)

    cstat = "         "
    write(cstat, '(I0)') stat
    cstat = trim(cstat)

    if(stat .eq. shell_ok) return

    call die("Trying to make directory '" // directory // "' returned status code " // cstat )

  end subroutine mkdir

  ! ---------------------------------------------------------------------------------------------- !
  !  MIT License
  !
  !  Copyright (c) 2020 fpm contributors
  !
  !  Permission is hereby granted, free of charge, to any person obtaining a copy
  !  of this software and associated documentation files (the "Software"), to deal
  !  in the Software without restriction, including without limitation the rights
  !  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  !  copies of the Software, and to permit persons to whom the Software is
  !  furnished to do so, subject to the following conditions:
  !
  !  The above copyright notice and this permission notice shall be included in all
  !  copies or substantial portions of the Software.
  !
  !  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  !  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  !  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  !  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  !  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  !  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  !  SOFTWARE.
  ! ---------------------------------------------------------------------------------------------- !
  integer function get_os_type() result(r)
    !! Returns one of OS_UNKNOWN, OS_LINUX, OS_MACOS, OS_WINDOWS, OS_CYGWIN,
    !! OS_SOLARIS, OS_FREEBSD, OS_OPENBSD.
    !!
    !! At first, the environment variable `OS` is checked, which is usually
    !! found on Windows. Then, `OSTYPE` is read in and compared with common
    !! names. If this fails too, check the existence of files that can be
    !! found on specific system types only.
    !!
    !! Returns OS_UNKNOWN if the operating system cannot be determined.
    character(len=255) :: val
    integer            :: length, rc
    logical            :: file_exists
    logical, save      :: first_run = .true.
    integer, save      :: ret = OS_UNKNOWN

    if (.not. first_run) then
        r = ret
        return
    end if

    first_run = .false.
    r = OS_UNKNOWN

    ! Check environment variable `OSTYPE`.
    call get_environment_variable('OSTYPE', val, length, rc)

    if (rc == 0 .and. length > 0) then
        ! Linux
        if (index(val, 'linux') > 0) then
            r = OS_LINUX
            ret = r
            return
        end if

        ! macOS
        if (index(val, 'darwin') > 0) then
            r = OS_MACOS
            ret = r
            return
        end if

        ! Windows, MSYS, MinGW, Git Bash
        if (index(val, 'win') > 0 .or. index(val, 'msys') > 0) then
            r = OS_WINDOWS
            ret = r
            return
        end if

        ! Cygwin
        if (index(val, 'cygwin') > 0) then
            r = OS_CYGWIN
            ret = r
            return
        end if

        ! Solaris, OpenIndiana, ...
        if (index(val, 'SunOS') > 0 .or. index(val, 'solaris') > 0) then
            r = OS_SOLARIS
            ret = r
            return
        end if

        ! FreeBSD
        if (index(val, 'FreeBSD') > 0 .or. index(val, 'freebsd') > 0) then
            r = OS_FREEBSD
            ret = r
            return
        end if

        ! OpenBSD
        if (index(val, 'OpenBSD') > 0 .or. index(val, 'openbsd') > 0) then
            r = OS_OPENBSD
            ret = r
            return
        end if
    end if

    ! Check environment variable `OS`.
    call get_environment_variable('OS', val, length, rc)

    if (rc == 0 .and. length > 0 .and. index(val, 'Windows_NT') > 0) then
        r = OS_WINDOWS
        ret = r
        return
    end if

    ! Linux
    inquire (file='/etc/os-release', exist=file_exists)

    if (file_exists) then
        r = OS_LINUX
        ret = r
        return
    end if

    ! macOS
    inquire (file='/usr/bin/sw_vers', exist=file_exists)

    if (file_exists) then
        r = OS_MACOS
        ret = r
        return
    end if

    ! FreeBSD
    inquire (file='/bin/freebsd-version', exist=file_exists)

    if (file_exists) then
        r = OS_FREEBSD
        ret = r
        return
    end if
  end function get_os_type

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function OS_NAME(os)
      integer, intent(in) :: os
      character(len=:), allocatable :: OS_NAME

      select case (os)
          case (OS_LINUX);   OS_NAME =  "Linux"
          case (OS_MACOS);   OS_NAME =  "macOS"
          case (OS_WINDOWS); OS_NAME =  "Windows"
          case (OS_CYGWIN);  OS_NAME =  "Cygwin"
          case (OS_SOLARIS); OS_NAME =  "Solaris"
          case (OS_FREEBSD); OS_NAME =  "FreeBSD"
          case (OS_OPENBSD); OS_NAME =  "OpenBSD"
          case (OS_UNKNOWN); OS_NAME =  "Unknown"
          case (OS_ALL)    ; OS_NAME =  "all"
          case default     ; OS_NAME =  "UNKNOWN"
      end select
  end function OS_NAME

! ================================================================================================================================ !
end module rotex__system
! ================================================================================================================================ !
