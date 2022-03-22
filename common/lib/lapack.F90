!> Interface wrapper for the lapack routines. See the <a href="http://www.netlib.org/lapack/">lapack
!> project documentation</a> for more details
module common_lapack

  use common_accuracy, only : rsp, rdp

  implicit none
  private

  public :: dsyev, dsygv


  !> Double precision symmetric eigensolver
  interface dsyev

    !> Double precision symmetric eigensolver
    subroutine dsyev(jobz, uplo, nn, aa, lda, ww, work, lwork, info)
      import rdp

      !> job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> workspace
      real(rdp), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine dsyev

  end interface dsyev


  !> Double precision generalised symmetric eigensolver
  interface dsygv

    !> Double precision generalised symmetric eigensolver
    subroutine dsygv(itype, jobz, uplo, nn, aa, lda, bb, ldb, ww, work, lwork, info)
      import rdp

      !> Specifies the problem type to be solved
      integer, intent(in) :: itype

      !> job type
      character, intent(in) :: jobz

      !> Upper 'U' or lower 'L' triangle
      character, intent(in) :: uplo

      !> matrix dimension
      integer, intent(in) :: nn

      !> Leading dimension of A
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(inout) :: aa(lda, *)

      !> leading dimension of B
      integer, intent(in) :: ldb

      !> matrix B
      real(rdp), intent(inout) :: bb(ldb, *)

      !> eigenvalues
      real(rdp), intent(out) :: ww(*)

      !> workspace
      real(rdp), intent(inout) :: work(*)

      !> workspace sizing
      integer, intent(in) :: lwork

      !> state of routine on return
      integer, intent(out) :: info

    end subroutine dsygv

  end interface dsygv

end module common_lapack
