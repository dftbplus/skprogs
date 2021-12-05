module core_overlap

  use common_accuracy, only : dp
  use common_constants
  use utilities
  use integration

  implicit none
  private

  public :: overlap, kinetic, nuclear, moments, v, confinement

contains

  subroutine overlap(s,max_l,num_alpha,alpha,poly_order)

    ! Overlap matrix elements, see rmp_32_186_1960.pdf eqn. 5 and eqn. 19


    ! Definition of the primitive basis functions based on Roothaan:
    ! R_{\lambda p}=1/sqrt((2n_{\lambda p})!)*
    !               (2*\zeta_{\lambda p})**(n_{\lambda p}+0.5)*
    !               r**(n_{\lambda p}-1)*exp(-\zeta_{\lambda p}*r)
    !
    ! For every exponent \zeta_{\lambda p} there are num_power coefficients,
    ! each connected to one r**(n_{\lambda p}-1). The sum over all
    ! coefficients, e.g. implicitely \zeta and r**n, gives the usual DFTB
    ! basis function.
    !
    ! Note: in DFTB one usually has r**(n+l-1) explicitely, here the angular
    ! momentum index l is implicit. Result:
    ! for l=0, e.g. s, n_{\lambda p}=0,1,...,num_power
    ! for l=1, e.g. p, n_{\lambda p}=1,2,...,num_power+1
    ! for l=2, e.g. d, n_{\lambda p}=2,3,...,num_power+2
    !

    real(dp), intent(out) :: s(0:,:,:)
    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq
    real(dp) :: alpha1

    s=0.0d0

    ! These loops define the indizes S_{\lambda p q}
    ! p=alpha1/n=0+l,alpha1/n=1+l,...,alpha2/n=0+l,alpha2/n=1+l...
    !
    !  write(*,*) 'max_l',max_l
    !  write(*,*) 'num_alpha',num_alpha
    !  write(*,*) 'poly_order',poly_order
    !  write(*,'(A)') 'ii jj ll kk mm nn oo'
    do ii=0,max_l
      nn=0
      do jj=1,num_alpha(ii)
        do ll=1,poly_order(ii)
          nn=nn+1
          oo=0
          nlp=ll+ii
          do kk=1,num_alpha(ii)
            alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
            do mm=1,poly_order(ii)
              oo=oo+1
              nlq=mm+ii
              !            write(*,'(I2,I2,I2,I2,I2,I2,I2)') ii,jj,ll,kk,mm,nn,oo
              !
              ! use ll+ii and mm+ii becaue of DFTB basis function definition
              s(ii,nn,oo)=1.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                  &v(alpha(ii,kk),2*nlq))*v(alpha1,nlp+nlq)
            end do
          end do
        end do
      end do
    end do

    !  write(*,*) 'OVERLAP'
    !  write(*,*) s

  end subroutine overlap

  subroutine nuclear(u,max_l,num_alpha,alpha,poly_order)

    ! Nuclear attraction matrix elements, see rmp_32_186_1960.pdf eqn. 5 and eqn.19


    real(dp), intent(out) :: u(0:,:,:)
    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq
    real(dp) :: alpha1

    u=0.0d0

    do ii=0,max_l
      nn=0
      do jj=1,num_alpha(ii)
        do ll=1,poly_order(ii)
          nn=nn+1
          oo=0
          nlp=ll+ii
          do kk=1,num_alpha(ii)
            alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
            do mm=1,poly_order(ii)
              oo=oo+1
              nlq=mm+ii
              u(ii,nn,oo)=2.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                  &v(alpha(ii,kk),2*nlq))*v(alpha1,nlp+nlq-1)
            end do
          end do
        end do
      end do
    end do

    !  write(*,*) 'NUCLEAR'
    !  write(*,*) u

  end subroutine nuclear

  ! WARNING: a finite nucleus is a bad idea with the currently implemented ZORA,
  ! because the integration by parts done there does certainly fail with a finite
  ! nucleus. Second: this routine does not even work without ZORA, unknown bug.
  !
  !  subroutine nuclear_finite(u,nuc,max_l,num_alpha,alpha,poly_order)
  !! simple finite nucleus
  !! v=-Z/(2R_0)*(3-r^2/R_0^2) for r<=R_0
  !! v=-Z/r for r>R_0
  !
  !  implicit none
  !
  !  real(dp), intent(out) :: u(0:,:,:)
  !  integer, intent(in) :: max_l,nuc
  !  integer, intent(in) :: num_alpha(0:)
  !  integer, intent(in) :: poly_order(0:)
  !  real(dp), intent(in) :: alpha(0:,:)
  !  integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq
  !  real(dp) :: alpha1,alpha2,part1,part2,part3,part4,part5,part6,r0,normalization
  !  integer :: iso(109)
  !  DATA iso/&
  !  &1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40,&
  !  &45, 48, 51, 52, 55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79, 84, 85, 88, 89,&
  !  &90, 93, 98, 98, 102, 103, 106, 107, 114, 115, 120, 121, 130, 127, 132, 133,&
  !  &138, 139, 140, 141, 144, 145, 152, 153, 158, 159, 162, 162, 168, 169, 174, &
  !  &175, 180, 181, 184, 187, 192, 193, 195, 197, 202, 205, 208, 209, 209, 210,&
  !  &222, 223, 226, 227, 232, 231, 238, 237, 244, 243, 247, 247, 251, 252, 257,&
  !  &258, 259, 262, 261, 262, 263, 262, 265, 266/
  !
  !  r0=sqrt(5.0d0/3.0d0)*(0.836*(iso(nuc)**(1.0d0/3.0d0))+0.570)*1.0d-5/0.529177d0
  !
  !  write(*,'(A,E)') 'FINITE NUCLEUS MODEL, RADIUS ',r0
  !
  !  u=0.0d0
  !
  !  do ii=0,max_l
  !    nn=0
  !    do jj=1,num_alpha(ii)
  !      do ll=1,poly_order(ii)
  !        nn=nn+1
  !        oo=0
  !        nlp=ll+ii
  !        do kk=1,num_alpha(ii)
  !          alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
  !          alpha2=-(alpha(ii,jj)+alpha(ii,kk))
  !          do mm=1,poly_order(ii)
  !            oo=oo+1
  !            nlq=mm+ii
  !
  !            normalization=real(2**(nlp+nlq+1),dp)/&
  !                sqrt(v(alpha(ii,jj),2*nlp)*v(alpha(ii,kk),2*nlq))
  !
  !            part1=exp_int(alpha2,nlp+nlq-1,r0)-exp_int(alpha2,nlp+nlq-1,0.0d0)
  !            part2=(exp_int(alpha2,nlp+nlq,r0)-&
  !                    &exp_int(alpha2,nlp+nlq,0.0d0))*3.0d0/(2.0d0*r0)
  !            part3=(exp_int(alpha2,nlp+nlq+2,r0)-&
  !                    &exp_int(alpha2,nlp+nlq+2,0.0d0))/(2.0d0*(r0**3))
  !
  !            u(ii,nn,oo)=2.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
  !                  &v(alpha(ii,kk),2*nlq))*v(alpha1,nlp+nlq-1)&
  !                  &-normalization*(+part1-part2+part3)
  !            write(*,*) 'part1',part1
  !            write(*,*) 'part2',part2
  !            write(*,*) 'part3',part3
  !            write(*,*) 'norma',normalization
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !
  !  write(*,*) 'NUCLEAR FINITE'
  !  write(*,*) u
  !
  !  end subroutine nuclear_finite

  subroutine kinetic(t,max_l,num_alpha,alpha,poly_order)

    ! Kinetic matrix elements, see rmp_32_186_1960.pdf eqn. 5 and eqn. 19

    real(dp), intent(out) :: t(0:,:,:)
    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq
    real(dp) :: alpha1

    t=0.0d0

    do ii=0,max_l
      nn=0
      do jj=1,num_alpha(ii)
        do ll=1,poly_order(ii)
          nn=nn+1
          oo=0
          nlp=ll+ii
          do kk=1,num_alpha(ii)
            alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
            do mm=1,poly_order(ii)
              oo=oo+1
              nlq=mm+ii
              t(ii,nn,oo)=0.5d0*alpha(ii,jj)*alpha(ii,kk)/&
                  &sqrt(v(alpha(ii,jj),2*nlp)*v(alpha(ii,kk),2*nlq))*&
                  &(v(alpha1,nlp+nlq)-&
                  &(w(alpha(ii,jj),ii,nlp)+w(alpha(ii,kk),ii,nlq))*&
                  &v(alpha1,nlp+nlq-1)+&
                  &(w(alpha(ii,jj),ii,nlp)*w(alpha(ii,kk),ii,nlq))*&
                  &v(alpha1,nlp+nlq-2)&
                  &)
            end do
          end do
        end do
      end do
    end do

    !  write(*,*) 'KINETIC'
    !  write(*,*) t

  end subroutine kinetic

  subroutine confinement(vconf,max_l,num_alpha,alpha,poly_order,&
      &conf_r0,conf_power)

    ! Analytic matrix elements of confining potential
    ! No checking for power, e.g. power==0 or power<0 etc. !

    real(dp), intent(out) :: vconf(0:,:,:)
    integer, intent(in) :: max_l,conf_power(0:)
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:,:),conf_r0(0:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq
    real(dp) :: alpha1

    vconf=0.0d0

    do ii=0,max_l
      if (conf_power(ii)/=0) then
        nn=0
        do jj=1,num_alpha(ii)
          do ll=1,poly_order(ii)
            nn=nn+1
            oo=0
            nlp=ll+ii
            do kk=1,num_alpha(ii)
              alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
              do mm=1,poly_order(ii)
                oo=oo+1
                nlq=mm+ii
                vconf(ii,nn,oo)=1.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                    &v(alpha(ii,kk),2*nlq))/(conf_r0(ii)*2.0d0)**conf_power(ii)*&
                    &v(alpha1,nlp+nlq+conf_power(ii))
              end do
            end do
          end do
        end do
      end if
    end do

    !  write(*,*) 'CONFINEMENT'
    !  write(*,*) vconf

  end subroutine confinement

  subroutine moments(moment,max_l,num_alpha,alpha,poly_order,problemsize,cof,&
      &power)

    ! Arbitrary moments of electron distribution, e.g. expectation values
    ! of <r>, <r^2> etc.; this is implemented analytically for arbitrary
    ! powers

    real(dp), intent(out) :: moment(:,0:,:)
    integer, intent(in) :: max_l,problemsize
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:),power
    real(dp), intent(in) :: alpha(0:,:),cof(:,0:,:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,nlp,nlq
    real(dp) :: alpha1

    moment=0.0d0

    ! <r^-3> only computed for p-functions and higher
    if (power>-3) then
      do ii=0,max_l
        do pp=1,num_alpha(ii)*poly_order(ii)
          nn=0
          do jj=1,num_alpha(ii)
            do ll=1,poly_order(ii)
              nn=nn+1
              oo=0
              nlp=ll+ii
              do kk=1,num_alpha(ii)
                alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
                do mm=1,poly_order(ii)
                  oo=oo+1
                  nlq=mm+ii

                  moment(1,ii,pp)=moment(1,ii,pp)+1.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                      &v(alpha(ii,kk),2*nlq))/(2.0d0**power)*&
                      &v(alpha1,nlp+nlq+power)*cof(1,ii,nn,pp)*cof(1,ii,oo,pp)

                  moment(2,ii,pp)=moment(2,ii,pp)+1.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                      &v(alpha(ii,kk),2*nlq))/(2.0d0**power)*&
                      &v(alpha1,nlp+nlq+power)*cof(2,ii,nn,pp)*cof(2,ii,oo,pp)

                end do
              end do
            end do
          end do
        end do
      end do
    else if (power==-3) then
      do ii=1,max_l
        do pp=1,num_alpha(ii)*poly_order(ii)
          nn=0
          do jj=1,num_alpha(ii)
            do ll=1,poly_order(ii)
              nn=nn+1
              oo=0
              nlp=ll+ii
              do kk=1,num_alpha(ii)
                alpha1=0.5d0*(alpha(ii,jj)+alpha(ii,kk))
                do mm=1,poly_order(ii)
                  oo=oo+1
                  nlq=mm+ii

                  moment(1,ii,pp)=moment(1,ii,pp)+1.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                      &v(alpha(ii,kk),2*nlq))/(2.0d0**power)*&
                      &v(alpha1,nlp+nlq+power)*cof(1,ii,nn,pp)*cof(1,ii,oo,pp)

                  moment(2,ii,pp)=moment(2,ii,pp)+1.0d0/sqrt(v(alpha(ii,jj),2*nlp)*&
                      &v(alpha(ii,kk),2*nlq))/(2.0d0**power)*&
                      &v(alpha1,nlp+nlq+power)*cof(2,ii,nn,pp)*cof(2,ii,oo,pp)

                end do
              end do
            end do
          end do
        end do
      end do
    end if

    !  write(*,*) 'MOMENT'
    !  write(*,*) moment

  end subroutine moments

  function v(x,i) ! V_{i}(x)

    ! Auxilliary function, see rmp_32_186_1960.pdf eqn. 20

    real(dp), intent(in) :: x
    integer, intent(in) :: i
    real(dp) :: v

    v=fak(i)/(x**(i+1))

    return  
  end function v

  function w(x,i,j) ! W_{ij}(x)

    ! Auxilliary function, see rmp_32_186_1960.pdf eqn. 20

    real(dp), intent(in) :: x
    integer, intent(in) :: i,j
    real(dp) :: w

    w=2.0d0*real((j-i-1),dp)/x

    return
  end function w

end module core_overlap
