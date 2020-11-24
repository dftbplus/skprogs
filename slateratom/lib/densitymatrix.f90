module densitymatrix
  use accuracy
  use constants
  use utilities
  implicit none
  private

  public :: densmatrix
  
contains

  subroutine densmatrix(problemsize,max_l,occ,cof,p)

    ! Get density matrix from wavefunction coefficients.

    real(dp), intent(in) :: cof(:,0:,:,:),occ(:,0:,:)
    integer, intent(in) :: problemsize,max_l
    real(dp), intent(out) :: p(:,0:,:,:)
    integer :: ii,jj,kk,ll,mm

    p=0.0d0

    do ii=1,2
      do jj=0,max_l
        do kk=1,problemsize
          do ll=kk,problemsize
            do mm=1,problemsize
              p(ii,jj,kk,ll)=p(ii,jj,kk,ll)+occ(ii,jj,mm)*&
                  &cof(ii,jj,kk,mm)*cof(ii,jj,ll,mm)
              p(ii,jj,ll,kk)=p(ii,jj,kk,ll)
            end do
          end do
        end do
      end do
    end do

    !  write(*,*) 'DENSITY MATRIX'
    !  write(*,*) p

  end subroutine densmatrix

end module densitymatrix
