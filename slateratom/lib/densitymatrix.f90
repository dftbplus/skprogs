!> Module that sets up the density matrix, based on occupations and wavefunction coefficients.
module densitymatrix

  use common_accuracy, only : dp

  implicit none
  private

  public :: densmatrix


contains

  !> Get density matrix from wavefunction coefficients.
  pure subroutine densmatrix(problemsize, max_l, occ, cof, pp)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> occupation numbers
    real(dp), intent(in) :: occ(:,0:,:)

    !> wavefunction coefficients
    real(dp), intent(in) :: cof(:,0:,:,:)

    !> density matrix supervector
    real(dp), intent(out) :: pp(:,0:,:,:)

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm

    pp(:,:,:,:) = 0.0_dp

    do ii = 1, 2
      do jj = 0, max_l
        do kk = 1, problemsize
          do ll = kk, problemsize
            do mm = 1, problemsize
              pp(ii, jj, kk, ll) = pp(ii, jj, kk, ll)&
                  & + occ(ii, jj, mm) * cof(ii, jj, kk, mm) * cof(ii, jj, ll, mm)
              pp(ii, jj, ll, kk) = pp(ii, jj, kk, ll)
            end do
          end do
        end do
      end do
    end do

  end subroutine densmatrix

end module densitymatrix
