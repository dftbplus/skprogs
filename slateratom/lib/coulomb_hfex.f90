module coulomb_hfex

  use common_accuracy, only : dp
  use common_constants
  use utilities
  use core_overlap

  implicit none
  private

  public :: coulomb, hfex
  
  
contains

  subroutine coulomb(j,max_l,num_alpha,alpha,poly_order,u,s)

    ! Coulomb supermatrix, see rmp_32_186_1960.pdf eqn. 6 and eqn. 21

    
    real(dp), intent(out) :: j(0:,:,:,0:,:,:)
    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:4,10)
    real(dp), intent(in) :: u(0:,:,:)
    real(dp), intent(in) :: s(0:,:,:)
    real(dp) :: alpha1,alpha2
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,xx,yy,zz
    integer :: nlpq,nmrs

    j=0.0d0

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              do nn=0,max_l
                uu=0
                do oo=1,num_alpha(nn)
                  do pp=1,poly_order(nn)
                    uu=uu+1
                    vv=0
                    do qq=1,num_alpha(nn)
                      do rr=1,poly_order(nn)
                        vv=vv+1

                        alpha1=(alpha(ii,jj)+alpha(ii,ll))/&
                            &(alpha(nn,oo)+alpha(nn,qq))
                        alpha2=(alpha(nn,oo)+alpha(nn,qq))/&
                            &(alpha(ii,jj)+alpha(ii,ll))
                        nlpq=kk+mm+2*ii
                        nmrs=pp+rr+2*nn

                        j(ii,ss,tt,nn,uu,vv)=&
                            &u(ii,ss,tt)*s(nn,uu,vv)*&
                            &c(nlpq-1,nmrs,alpha1)+&
                            &u(nn,uu,vv)*s(ii,ss,tt)*&
                            &c(nmrs-1,nlpq,alpha2)
                        !   write(*,'(A,F12.8,6I3)') 'j ',j(ii,ss,tt,nn,uu,vv),ii,ss,tt,nn,uu,vv
                        !   write(*,'(A,F12.8,3I3)') 's1',s(ii,ss,tt),ii,ss,tt
                        !   write(*,'(A,F12.8,3I3)') 's2',s(nn,uu,vv),nn,uu,vv
                        !   write(*,'(A,F12.8,3I3)') 'u1',u(ii,ss,tt),ii,ss,tt
                        !   write(*,'(A,F12.8,3I3)') 'u2',u(nn,uu,vv),nn,uu,vv
                        !   write(*,'(A,F12.8,2I3,F12.8)') 'c1',c(kk+mm+2*ii-1,pp+rr+2*nn,alpha1),&
                        !   &kk+mm+2*ii-1,pp+rr+2*nn,alpha1
                        !   write(*,'(A,F12.8,2I3,F12.8)') 'c2',c(pp+rr+2*nn-1,kk+mm+2*ii,alpha2),&
                        !   &pp+rr+2*ii-1,kk+mm+2*nn,alpha2
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    !  do ii=0,max_l
    !    do jj=0,max_l
    !      j(ii,:,:,jj,:,:)=j(ii,:,:,jj,:,:)/&
    !         &((2.0d0*float(ii)+1.0d0)*(2.0d0*float(jj)+1.0d0))
    !    end do
    !  end do

    !  write(*,*) 'COULOMB'
    !  write(*,*) j

  end subroutine coulomb

  subroutine hfex(k,max_l,num_alpha,alpha,poly_order,problemsize)

    ! HF Exchange supermatrix, see rmp_32_186_1960.pdf eqn. 7/8 and eqn. 21

    
    real(dp), intent(out) :: k(0:,:,:,0:,:,:)
    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:4,10)
    real(dp),allocatable :: knu(:,:,:,:,:,:,:)
    real(dp) :: alpha1,alpha2,alpha3,alpha4,beta1,beta2,beta3,beta4
    real(dp) :: pre,t1,t2,t3,t4
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,xx,yy,zz
    integer :: nu,problemsize
    integer :: nlp,nlq,nmr,nms

    allocate(knu(0:max_l,problemsize,problemsize,0:max_l,problemsize,&
        &problemsize,0:2*max_l+2))

    k=0.0d0
    knu=0.0d0

    ! Build knu according to eqn. 8

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              do nn=0,max_l
                uu=0
                do oo=1,num_alpha(nn)
                  do pp=1,poly_order(nn)
                    uu=uu+1
                    vv=0
                    do qq=1,num_alpha(nn)
                      do rr=1,poly_order(nn)
                        vv=vv+1

                        alpha1=0.5d0*(alpha(ii,jj)+alpha(nn,oo))
                        alpha2=0.5d0*(alpha(ii,ll)+alpha(nn,qq))
                        alpha3=0.5d0*(alpha(ii,jj)+alpha(nn,qq))
                        alpha4=0.5d0*(alpha(ii,ll)+alpha(nn,oo))
                        beta1=alpha1/alpha2
                        beta2=alpha2/alpha1
                        beta3=alpha3/alpha4
                        beta4=alpha4/alpha3
                        nlp=kk+ii
                        nlq=mm+ii
                        nmr=pp+nn
                        nms=rr+nn

                        pre=1.0d0/sqrt(v(alpha(ii,jj),2*(kk+ii))*&
                            &       v(alpha(ii,ll),2*(mm+ii))*&
                            &       v(alpha(nn,oo),2*(pp+nn))*&
                            &       v(alpha(nn,qq),2*(rr+nn)))

                        do nu=abs(ii-nn),ii+nn,2

                          t1=v(alpha1,nlp+nmr-nu-1)*v(alpha2,nlq+nms+nu)*&
                              &c(nlp+nmr-nu-1,nlq+nms+nu,beta1)
                          t2=v(alpha2,nlq+nms-nu-1)*v(alpha1,nlp+nmr+nu)*&
                              &c(nlq+nms-nu-1,nlp+nmr+nu,beta2)
                          t3=v(alpha3,nlp+nms-nu-1)*v(alpha4,nlq+nmr+nu)*&
                              &c(nlp+nms-nu-1,nlq+nmr+nu,beta3)
                          t4=v(alpha4,nlq+nmr-nu-1)*v(alpha3,nlp+nms+nu)*&
                              &c(nlq+nmr-nu-1,nlp+nms+nu,beta4)

                          knu(ii,ss,tt,nn,uu,vv,nu)=pre*(t1+t2+t3+t4)

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    ! Build k according to eqn. 7

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              do nn=0,max_l
                uu=0
                do oo=1,num_alpha(nn)
                  do pp=1,poly_order(nn)
                    uu=uu+1
                    vv=0
                    do qq=1,num_alpha(nn)
                      do rr=1,poly_order(nn)
                        vv=vv+1

                        do nu=abs(ii-nn),ii+nn,2

                          k(ii,ss,tt,nn,uu,vv)=k(ii,ss,tt,nn,uu,vv)+&
                              &almn(ii,nn,nu)*knu(ii,ss,tt,nn,uu,vv,nu)

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    !  do ii=0,max_l
    !    do jj=0,max_l
    !      k(ii,:,:,jj,:,:)=k(ii,:,:,jj,:,:)/&
    !         &((2.0d0*float(ii)+1.0d0)*(2.0d0*float(jj)+1.0d0))
    !    end do
    !  end do


    !  write(*,*) 'HF EXCHANGE'
    !  write(*,*) k

  end subroutine hfex

  function c(alpha,beta,t)

    ! Auxilliary function, see rmp_32_186_1960.pdf eqn. 22 and eqn. 23

    integer, intent(in)  :: alpha
    integer, intent(in)  :: beta
    real(dp), intent(in) :: t
    real(dp) :: c,factor
    real(dp), allocatable :: carray(:,:)
    integer :: ii,jj

    ! early return if index smaller than zero

    if (alpha<0) then
      c=0.0d0
      return
    end if

    if (beta<0) then
      c=0.0d0
      return
    end if

    allocate(carray(0:alpha,0:beta))

    factor=1.0d0/(1.0d0+t)

    ! Overall this is naive, the matrix could be reused to some extent ...
    ! OTOH, the matrices are relatively small.

    ! first handle Kronecker delta, three cases
    carray(0,0)=factor
    do ii=1,alpha
      carray(ii,0)=factor*(t*carray(ii-1,0)+1.0d0)
    end do
    do ii=1,beta
      carray(0,ii)=factor*(carray(0,ii-1))
    end do

    ! now build up from 1
    do ii=1,alpha
      do jj=1,beta
        carray(ii,jj)=factor*(t*carray(ii-1,jj)+carray(ii,jj-1))
      end do
    end do

    c=carray(alpha,beta)

    return
  end function c

  function a(rho)

    ! Auxilliary function, see rmp_32_186_1960.pdf eqn. 9

    
    integer, intent(in) :: rho
    real(dp) :: a

    a=fak(rho)/((fak(rho/2))**2)

  end function a

  function almn(lambda,mu,nu)

    ! Auxilliary function, see rmp_32_186_1960.pdf eqn. 9

    
    integer, intent(in) :: lambda,mu,nu
    real(dp) :: almn

    almn=a(lambda+mu-nu)*a(lambda-mu+nu)*a(mu-lambda+nu)/&
        &(float(lambda+mu+nu+1)*a(lambda+mu+nu))

  end function almn

end module coulomb_hfex
