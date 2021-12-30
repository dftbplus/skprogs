!> Module that handles command line argument parsing.
module cmdargs

  implicit none
  private

  public :: parse_command_arguments

  character(len=*), parameter :: programName = 'sktwocnt'
  character(len=*), parameter :: programVersion = '0.9'


contains

  !> Parses command line arguments or prints program/version information.
  subroutine parse_command_arguments()

    !> number of command line arguments and length buffer
    integer :: nArgs, argLen

    !> string representation of a single command line argument
    character(len=:), allocatable :: arg

    nArgs = command_argument_count()
    if (nArgs > 0) then
      call get_command_argument(1, length=argLen)
      allocate(character(argLen) :: arg)
      call get_command_argument(1, arg)
      select case (arg)
      case ('--version')
        write(*, '(A,1X,A)') programName, programVersion
        stop
      case default
        write(*, '(A,A,A)') "Invalid command line argument '", arg, "'"
        error stop
      end select
    end if

  end subroutine parse_command_arguments

end module cmdargs
