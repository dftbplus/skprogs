!> Contains F90 wrapper functions for some commonly used lapack calls needed in the code.
!! Contains some fixes for lapack 3.0 bugs, if this gets corrected in lapack 4.x they should be
!! removed.
module common_eigensolver

  use common_accuracy, only : rsp, rdp

  implicit none
  private

  public :: heev, hegv


  !> Simple eigensolver for a symmetric/Hermitian matrix
  !> Caveat: the matrix a is overwritten
  interface heev
    module procedure dble_dsyev
  end interface heev


  !> Simple eigensolver for a symmetric/Hermitian generalized matrix problem
  !> caveat: the matrix a is overwritten
  !> caveat: the matrix b is overwritten with Cholesky factorization
  interface hegv
    module procedure dble_dsygv
  end interface hegv

contains


  !> Double precision eigensolver for a symmetric matrix
  subroutine dble_dsyev(a, w, uplo, jobz)

    !> contains the matrix for the solver, returns as eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rdp), intent(inout) :: a(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of the matrix
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    real(rdp), allocatable :: work(:)
    integer n, info
    integer :: int_idealwork
    real(rdp) :: idealwork(1)
    character(len=100) :: error_string

    n = size(a, dim=1)

    call dsyev(jobz, uplo, n, a, n, w, idealwork, -1, info)

    if (info /= 0) then
      write(*,*) "Failure in DSYEV to determine optimum workspace"
      error stop
    end if

    int_idealwork = floor(idealwork(1))
    allocate(work(int_idealwork))

    call dsyev(jobz, uplo, n, a, n, w, work, int_idealwork, info)
    if (info /= 0) then
      if (info < 0) then
99020   format('Failure in diagonalisation routine dsyev, illegal argument at position ', i6)
        write(error_string, 99020) info
        write(*,*) error_string
        error stop
      else
99030   format('Failure in diagonalisation routine dsyev, diagonal element ', i6,&
            & ' did not converge to zero.')
        write(error_string, 99030) info
        write(*,*) error_string
        error stop
      end if
    end if

  end subroutine dble_dsyev


  !> Double precision eigensolver for generalized symmetric matrix problem
  subroutine dble_dsygv(a, b, w, uplo, jobz, itype)

    !> contains the matrix for the solver, returns eigenvectors if requested (matrix always
    !> overwritten on return anyway)
    real(rdp), intent(inout) :: a(:,:)

    !> contains the second matrix for the solver (overwritten by Cholesky factorization)
    real(rdp), intent(inout) :: b(:,:)

    !> eigenvalues
    real(rdp), intent(out) :: w(:)

    !> upper or lower triangle of both matrices
    character, intent(in) :: uplo

    !> compute eigenvalues 'N' or eigenvalues and eigenvectors 'V'
    character, intent(in) :: jobz

    !> specifies the problem type to be solved 1:A*x=(lambda)*B*x, 2:A*B*x=(lambda)*x,
    !> 3:B*A*x=(lambda)*x default is 1
    integer, optional, intent(in) :: itype

    real(rdp), allocatable :: work(:)
    integer n, lda, info, iitype, ldb
    integer :: int_idealwork
    real(rdp) :: idealwork(1)
    character(len=100) :: error_string

    n = size(a, dim=2)
    lda = size(a, dim=1)
    ldb = size(b, dim=1)

    if (present(itype)) then
      iitype = itype
    else
      iitype = 1
    end if

    call dsygv(iitype, jobz, uplo, n, a, lda, b, ldb, w, idealwork, -1, info)
    if (info /= 0) then
      write(*,*) "Failure in dsygv to determine optimum workspace"
      error stop
    end if

    int_idealwork = nint(idealwork(1))
    allocate(work(int_idealwork))

    call dsygv(iitype, jobz, uplo, n, a, lda, b, ldb, w, work, int_idealwork, info)

    if (info /= 0) then
      if (info < 0) then
        write(error_string, "('Failure in diagonalisation routine dsygv, illegal ',&
            & 'argument at position ',i6)") info
        write(*,*) error_string
        error stop
      else if (info <= n) then
        write(error_string, "('Failure in diagonalisation routine dsygv, diagonal ',&
            & 'element ', i6, ' did not converge to zero.')") info
        write(*,*) error_string
        error stop
      else
        write(error_string, "('Failure in diagonalisation routine dsygv,', &
            & ' non-positive definite overlap! Minor ',i6,' responsible.')") info - n
        write(*,*) error_string
        error stop
      end if
    end if

  end subroutine dble_dsygv

end module common_eigensolver
