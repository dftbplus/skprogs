#:include 'common.fypp'

!> Auxiliary subroutines for the ASSERT command
module common_assert

  implicit none

  private
#:block DEBUG_CODE
  public :: assertError
#:endblock DEBUG_CODE

contains

#:block DEBUG_CODE

  !> Prints assertion error and abort program execution.
  subroutine assertError(fileName, lineNr, message)

    !> Name of the file in which the error occurred.
    character(*), intent(in) :: fileName

    !> Nr. of the line at which the error occurred.
    integer, intent(in) :: lineNr

    !> Additional message for error
    character(*), intent(in), optional :: message

    write(*, '(A)') "!!! UNFULLFILLED ASSERTION"
    write(*, '(A,A)') "!!! FILE:      ", fileName
    write(*, '(A,I0)') "!!! LINE NR.:  ", lineNr
    if (present(message)) then
      write(*, '(A,A,A)') '!!! MESSAGE:  "', trim(message), '"'
    end if
    error stop

  end subroutine assertError

#:endblock DEBUG_CODE

end module common_assert
