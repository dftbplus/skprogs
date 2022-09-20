!> Module that contains routines for throwing messages (warnings/errors) due to runtime events.
module common_message

  implicit none
  private


  !> Recoverable error warnings.
  interface warning
    module procedure warning_single
    module procedure warning_array
  end interface


  !> Fatal error warnings, terminating the code.
  interface error
    module procedure error_single
    module procedure error_array
  end interface error


  public :: warning, error


contains

  !> Throws a single warning message.
  subroutine warning_single(message)

    !> Warning message to print to standard out
    character (len=*), intent(in) :: message

    write(*, '(1a)') 'WARNING!'
    write(*, '(2a)') '-> ', trim(message)

  end subroutine warning_single


  !> Throws multiple warning messages.
  subroutine warning_array(messages)

    !> Lines of the warning message to print to standard out
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    write(*, '(1a)') 'WARNING!'
    do ii = 1, size(messages)
      write(*, '(2a)') '-> ', trim(messages(ii))
    end do

  end subroutine warning_array


  !> Throws a single error message and stops the code.
  subroutine error_single(message)

    !> Error message to print to standard out
    character (len=*), intent(in) :: message

    write(*, '(1a)') 'ERROR!'
    write(*, '(2a)') '-> ', trim(message)
    error stop

  end subroutine error_single


  !> Throws multiple error messages and stops the code.
  subroutine error_array(messages)

    !> Lines of the error message to print to standard out
    character(len=*), intent(in) :: messages(:)

    integer :: ii

    write(*, '(1a)') 'ERROR!'
    do ii = 1, size(messages)
      write(*, '(2a)') '-> ', trim(messages(ii))
    end do
    error stop

  end subroutine error_array

end module common_message
