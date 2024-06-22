!> Interface wrapper for the lapack routines. See the <a href="http://www.netlib.org/lapack/">lapack
!! project documentation</a> for more details
module lapack

  use common_accuracy, only : rsp, rdp
  implicit none
  public


  !> Computes LU factorization of double precision matrix
  interface dgetrf

    !> Computes LU factorization of double precision matrix
    subroutine dgetrf(mm, nn, aa, lda, ipiv, info)
      import rdp

      !> number of rows of the matrix
      integer, intent(in) :: mm

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> pivot array
      integer, intent(out) :: ipiv(*)

      !> state of routine on return
      integer, intent(out) :: info
    end subroutine dgetrf

  end interface dgetrf

end module lapack
