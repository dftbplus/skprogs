#:include 'common.fypp'

!> Contains F90 wrapper functions for some commonly used lapack calls needed in the code.
!! The interface of all LAPACK calls must be defined in the module lapack.
module lapackroutines

  use common_accuracy, only : dp, rdp
  implicit none

  private
  public :: getrf, getrs


  !> Computes the LU decomposition of a general rectangular matrix using partial pivoting with row
  !! interchanges.
  !! The decomposition has the form: A = P*L*U, where P is a permutation matrix, L is a lower
  !! triangular matrix with unit diagonal elements and U is an upper triangular matrix.
  interface getrf
    module procedure getrf_dble
  end interface getrf


  !> Solves a system of linear equations
  !! A * X = B  or  A**T * X = B
  !! with a general N-by-N matrix A using the LU factorization computed by getrf.
  interface getrs
    module procedure :: getrs_dble
    module procedure :: getrs1_dble
  end interface getrs


contains

  !> Double precision version of getrf.
  subroutine getrf_dble(aa, ipiv, nRow, nColumn, iError)

    !> Matrix to decompose on entry, L and U on exit. Unit diagonal elements of L are not stored.
    real(rdp), intent(inout) :: aa(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(out) :: ipiv(:)

    !> Number of rows of the matrix to decomposea. (Necessary if different from the number of rows
    !> of the passed matrix)
    integer, optional, intent(in) :: nRow

    !> Number of rows of the matrix to decompose. (Necessary if different from the number of columns
    !> of the passed matrix)
    integer, optional, intent(in) :: nColumn

    !> Error flag. Zero on successful exit. If not present, any lapack error causes program
    !> termination. If passed only fatal lapack errors with error flag < 0 cause abort.
    integer, optional, intent(out) :: iError

    integer :: mm, nn, lda, info
    character(len=100) :: error_string

    lda = size(aa, dim=1)
    nn = size(aa, dim=2)
    if (present(nRow)) then
      @:ASSERT(nRow >= 1 .and. nRow <= lda)
      mm = nRow
    else
      mm = lda
    end if
    if (present(nColumn)) then
      @:ASSERT(nColumn >= 1 .and. nColumn <= nn)
      nn = nColumn
    end if
    @:ASSERT(size(ipiv) == min(mm, nn))

    call dgetrf(mm, nn, aa, lda, ipiv, info)

    if (info < 0) then
99060 format('Failure in LU factorisation dgetrf,', ' illegal argument at position ', i10)
      write(error_string, 99060) info
      error stop error_string
    else
      if (present(iError)) then
        iError = info
      elseif (info > 0) then
99070   format('Factor U is exactly zero in dgetrf,', ' info flag is ', i10)
        write(error_string, 99070) info
        error stop error_string
      end if
    end if

  end subroutine getrf_dble


  !> Solves a system of linear equations with multiple right hand sides
  subroutine getrs_dble(amat, ipiv, bmat, trans, iError)

    !> Matrix of the linear system
    real(rdp), intent(in) :: amat(:, :)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(in) :: ipiv(:)

    !> Matrix of the right hand side vectors
    real(rdp), intent(inout) :: bmat(:, :)

    !> Optional transpose (defaults to 'n')
    character(len=1), intent(in), optional :: trans

    !> Error flag, zero on successful exit
    integer, intent(out), optional :: iError

    character(len=1) :: atr
    integer :: info, nn, nrhs, lda, ldb

    @:ASSERT(size(amat, 1) == size(amat, dim=2))
    @:ASSERT(size(amat, 1) == size(bmat, dim=1))

    if (present(trans)) then
      @:ASSERT(any(trans == ['n', 'N', 't', 'T', 'c', 'C']))
      atr = trans
    else
      atr = 'n'
    end if

    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    nn = size(amat, 2)
    nrhs = size(bmat, 2)

    call dgetrs(atr, nn, nrhs, amat, lda, ipiv, bmat, ldb, info)

    if (present(iError)) then
      iError = info
    else
      if (info /= 0) then
        error stop "Failed to solve linear system by diagonal pivoting"
      end if
    end if

  end subroutine getrs_dble


  !> Solves a system of linear equations with one right hand side
  subroutine getrs1_dble(amat, ipiv, bvec, trans, iError)

    !> Matrix of the linear system
    real(rdp), intent(in) :: amat(:,:)

    !> Pivot indices, row i of the matrix was interchanged with row ipiv(i).
    integer, intent(in) :: ipiv(:)

    !> Right hand side vector
    real(rdp), intent(inout), target :: bvec(:)

    !> optional transpose (defaults to 'n')
    character(len=1), intent(in), optional :: trans

    !> Error flag, zero on successful exit
    integer, intent(out), optional :: iError

    real(rdp), pointer :: bptr(:,:)

    bptr(1:size(bvec, 1), 1:1) => bvec(1:size(bvec, 1))
    call getrs(amat, ipiv, bptr, trans, iError)

  end subroutine getrs1_dble

end module lapackroutines
