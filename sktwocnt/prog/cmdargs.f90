module cmdargs
  implicit none

  character(*), parameter :: programName = 'sktwocnt'
  character(*), parameter :: programVersion = '0.9'

contains

  subroutine parse_command_arguments()

    integer :: nArgs, argLen
    character(:), allocatable :: arg

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
