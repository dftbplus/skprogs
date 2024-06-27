#:include 'common.fypp'

!> Contains F90 wrapper functions for some commonly used blas calls needed in the code.
!! The interface of all BLAS calls must be defined in the module blas.
module blasroutines

  use common_accuracy, only : rdp
  use blas, only : dger
  implicit none

  private
  public :: ger


  !> Rank 1 update of a matrix A := alpha*x*y' + A
  !! Wrapper for the level 2 blas routine xger to perform the rank 1 update of a general matrix
  interface ger
    module procedure ger_dble
  end interface ger


contains

  !> Double precision rank 1 update of a general matrix
  subroutine ger_dble(a, alpha, x, y)

    !> Contains the matrix for the update
    real(rdp), intent(inout) :: a(:,:)

    !> Scaling value for the update contribution
    real(rdp), intent(in) :: alpha

    !> Vector of values for the update
    real(rdp), intent(in) :: x(:)

    !> Vector of values for the update
    real(rdp), intent(in) :: y(:)

    integer :: n, m

    @:ASSERT(size(a,dim=1) == size(x))
    @:ASSERT(size(a,dim=2) == size(y))

    m = size(x)
    n = size(y)

    call dger(m, n, alpha, x, 1, y, 1, a, m)

  end subroutine ger_dble

end module blasroutines
