!> Interface wrapper for the blas routines.
!!
!! ALL BLAS routines which are called from the main code must be included here.
module blas

  use common_accuracy, only : rdp
  public

  interface


    !> Performs the rank 1 operation
    !! A := alpha*x*y**T + A,
    subroutine dger(mm, nn, alpha, xx, incx, yy, incy, aa, lda)
      import rdp

      !> Matrix sizing
      integer, intent(in) :: mm

      !> Matrix  size
      integer, intent(in) :: nn

      !> Scale factor
      real(rdp), intent(in) :: alpha

      !> Vector
      real(rdp), intent(in) :: xx(*)

      !> Stride
      integer, intent(in) :: incx

      !> Vector
      real(rdp), intent(in) :: yy(*)

      !> Stride
      integer, intent(in) :: incy

      !> Leading matrix dimension
      integer, intent(in) :: lda

      !> Matrix A
      real(rdp), intent(inout) :: aa(lda, *)

    end subroutine dger

  end interface

end module blas
