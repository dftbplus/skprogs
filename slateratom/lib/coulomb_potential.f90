module coulomb_potential

! the routines in this module server output purposes only
! during SCF except in the ZORA case, but even then the Coulomb matrix 
! (J supermatrix) elements are calculated directly

  use accuracy
  use utilities
  use integration
  use core_overlap
  implicit none
  private

  public :: cou_pot
  
contains

  subroutine cou_pot(p,max_l,num_alpha,poly_order,alpha,problemsize,&
      &num_points,abcissa,cpot)
    ! calculate coulomb potential on arbitraty set of points
    ! by analytical evaluation of the integrals indicated
    !                _                                         _
    !               |                                           |
    !               | 1  r      2            rmax               |
    !  V(r)= 4*PI * | - int * r' * rho(r') + int  r' * rho (r') |
    !               | r  0                    r                 |
    !               |_                                         _|
    !                          help1                help2

    implicit none

    real(dp), intent(in) :: p(0:,:,:),abcissa(:),alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),num_points
    integer, intent(in) :: problemsize
    real(dp), intent(out) :: cpot(:)
    real(dp), allocatable :: help1(:,:,:,:),help2(:,:,:,:)
    real(dp) :: alpha1
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq

    allocate(help1(num_points,0:max_l,problemsize,problemsize))
    allocate(help2(num_points,0:max_l,problemsize,problemsize))

    help1=0.0d0
    help2=0.0d0
    cpot=0.0d0

    ! get integrals for pairs of basis functions
    do ii=0,max_l
      ll=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ll=ll+1
          oo=0
          nlp=kk+ii
          do mm=1,num_alpha(ii)

            ! exp_int has no notion of implicit "-" of alpha 
            alpha1=-(alpha(ii,jj)+alpha(ii,mm))

            do nn=1,poly_order(ii)
              oo=oo+1
              nlq=nn+ii

              ! integrals as indicated in comment, no normalization
              do pp=1,num_points
                help1(pp,ii,ll,oo)=(exp_int(alpha1,nlp+nlq,abcissa(pp))-&
                    &exp_int(alpha1,nlp+nlq,0.0d0))/abcissa(pp)
                help2(pp,ii,ll,oo)=&
                    &-exp_int(alpha1,nlp+nlq-1,abcissa(pp))
              end do

              ! add normalization of basis functions
              ! watch out for 2**(nlp+nlq+1) needed because variable integration ranges
              help1(:,ii,ll,oo)=help1(:,ii,ll,oo)*float(2**(nlp+nlq+1))/&
                  &sqrt(v(alpha(ii,jj),2*nlp)*v(alpha(ii,mm),2*nlq))
              help2(:,ii,ll,oo)=help2(:,ii,ll,oo)*float(2**(nlp+nlq+1))/&
                  &sqrt(v(alpha(ii,jj),2*nlp)*v(alpha(ii,mm),2*nlq))

            end do
          end do
        end do
      end do
    end do

    ! now actually get potential, multiply with density matrix
    do pp=1,num_points
      do ii=0,max_l
        ll=0
        do jj=1,num_alpha(ii)
          do kk=1,poly_order(ii)
            ll=ll+1
            oo=0
            do mm=1,num_alpha(ii)
              do nn=1,poly_order(ii)
                oo=oo+1
                cpot(pp)=cpot(pp)+p(ii,ll,oo)*&
                    &(help1(pp,ii,ll,oo)+help2(pp,ii,ll,oo))
              end do
            end do
          end do
        end do
      end do
    end do

    !  write(*,*) 'CPOT'
    !  write(*,*) cpot

    deallocate(help1)
    deallocate(help2)

  end subroutine cou_pot

end module coulomb_potential
