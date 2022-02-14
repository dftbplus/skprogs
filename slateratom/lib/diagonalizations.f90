!> Module that provides routines for matrix diagonalization.
module diagonalizations

  use common_accuracy, only : dp
  use common_eigensolver, only : heev, hegv

  implicit none
  private

  public :: diagonalize_overlap, diagonalize


contains

  !> Diagonalizes overlap matrix to check for linear dependency of basis set.
  !! Implicitely LAPACK's dsyev is called.
  subroutine diagonalize_overlap(max_l, num_alpha, poly_order, ss)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> overlap supervector
    real(dp), intent(in) :: ss(0:, :,:)

    !! overlap matrices
    real(dp), allocatable :: overlap(:,:)

    !! eigenvalues of overlap matrices
    real(dp), allocatable :: eigenvalues(:)

    !> auxiliary variables
    integer :: ll, diagsize

    do ll = 0, max_l

      diagsize = num_alpha(ll) * poly_order(ll)

      allocate(eigenvalues(diagsize))
      eigenvalues(:) = 0.0_dp

      overlap = ss(ll, :,:)

      call heev(overlap, eigenvalues, 'U', 'N')

      write(*, '(A,I3,A,E16.8)') 'Smallest eigenvalue of overlap for l= ', ll, ' : ', eigenvalues(1)

      if (eigenvalues(1) < 1.0e-10_dp) then
        write(*, '(A)') ' '
        write(*, '(A)') 'Basis set is nearly linear dependent, reduction necessary'
        write(*, '(A)') ' '
        stop
      end if

      deallocate(overlap)
      deallocate(eigenvalues)

    end do
    write(*,*) ' '

  end subroutine diagonalize_overlap


  !> This is a driver for hegv. The idea is that the matrices are allocated in the main program
  !! for the maximum size of the problem but hegv is only fed with a matrix of the current size of
  !! the eigenproblem.
  subroutine diagonalize(max_l, num_alpha, poly_order, ff, ss, cof_new, eigval)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> fock matrix supervector
    real(dp), intent(in) :: ff(:, 0:, :,:)

    !> overlap supervector
    real(dp), intent(in) :: ss(0:, :,:)

    !> new wavefunction coefficients
    real(dp), intent(inout) :: cof_new(:, 0:, :,:)

    !> eigenvalues
    real(dp), intent(inout) :: eigval(:, 0:, :)

    !! fock matrices
    real(dp), allocatable :: fock(:,:)

    !! overlap matrices
    real(dp), allocatable :: overlap(:,:)

    !! eigenvalues of overlap matrices
    real(dp), allocatable :: eigenvalues(:)

    !! auxiliary variables
    integer :: iSpin, ll, diagsize

    do iSpin = 1, 2
      do ll = 0, max_l

        diagsize = num_alpha(ll) * poly_order(ll)

        allocate(eigenvalues(diagsize))
        eigenvalues(:) = 0.0_dp

        fock = ff(iSpin, ll, :,:)
        overlap = ss(ll, :,:)

        call hegv(fock, overlap, eigenvalues, 'U', 'V', itype=1)

        ! fock matrix overwritten by eigenvectors
        cof_new(iSpin, ll, :,:) = fock
        eigval(iSpin, ll, :) = eigenvalues

        deallocate(fock, overlap, eigenvalues)

        ! Debug:
        ! real(dp), allocatable :: res1(:,:), res2(:,:)
        ! integer, allocatable :: loc(:)
        ! integer :: i, j, alpha, beta

        ! allocate(res1(diagsize, diagsize))
        ! allocate(res2(diagsize, diagsize))

        ! res1(:,:) = 0.0_dp
        ! res2(:,:) = 0.0_dp
        ! do alpha = 1, diagsize
        !   do beta = 1, diagsize
        !     do i = 1, diagsize
        !       do j = 1, diagsize
        !         res2(alpha, beta) = res2(alpha, beta)&
        !             & + fock(i, alpha) * ss(ll, i, j) * fock(j, beta)
        !         res1(alpha, beta) = res1(alpha, beta)&
        !             & + temp1(i, alpha) * ss(ll, i, j) * temp1(j, beta)
        !       end do
        !     end do
        !   end do
        ! end do
        ! loc = maxloc(abs(res1 - res2))
        ! print *, loc, res1(loc(1), loc(2)), res2(loc(1), loc(2))
        ! deallocate(res1, res2)

        ! print '(12F10.6)', matmul(reshape(fock(:, 1), [diagsize, 1]),&
        !     & reshape(matmul(ss(ll, :,:), fock(:, 1)), [1, diagsize]))

        ! do alpha = 1, diagsize
        !   print *, maxval((matmul(ff(iSpin, ll, :,:), temp1(:, alpha)) / matmul(ss(ll, :,:),&
        !       & temp1(:, alpha))) - temp2(alpha))
        !   print *, (matmul(ff(iSpin, ll, :,:), fock(:, alpha)) / matmul(ss(ll, :,:), fock(:, alpha)))&
        !       & - eigenvalues(alpha)
        ! end do

        ! print *, matmul(fock(:, 1), transpose(matmul(overlap, fock(:, 1))))

        ! print *, 'eigprec: ', maxval(abs(abs(eigenvalues) - abs(temp2)))
        ! print *, 'vecprec: ', maxval(abs(abs(fock) - abs(temp1)))

      end do
    end do

  end subroutine diagonalize

end module diagonalizations
