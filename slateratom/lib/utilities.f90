module utilities

  use common_accuracy, only : dp
  use common_constants

  implicit none
  private

  public :: check_convergence, check_electron_number, vector_length
  public :: fak, polcart, cartpol, dscalar
  public :: v

  interface v
    module procedure v_int, v_real
  end interface
  
contains

  subroutine check_convergence(pot_old,pot_new,max_l,problemsize,iter,&
      &change_max,final)

    ! check SCF convergence

    real(dp), intent(out) :: change_max
    real(dp), intent(in) :: pot_old(:,0:,:,:),pot_new(:,0:,:,:)
    integer, intent(in) :: max_l,problemsize,iter
    logical, intent(out) :: final
    integer ii,jj,kk,ll

    change_max=0.0d0
    if (iter<3) then
      final=.false.
    end if

    do ii=1,2
      do jj=0,max_l
        do kk=1,problemsize
          do ll=1,problemsize
            change_max=max(change_max,&
                &abs(pot_old(ii,jj,kk,ll)-pot_new(ii,jj,kk,ll)))
          end do
        end do
      end do
    end do

    if (change_max<1.0d-8) then
      final=.true.
    end if

  end subroutine check_convergence

  subroutine check_electron_number(cof,s,occ,max_l,num_alpha,poly_order,&
      &problemsize)

    ! check conservation of electron number during SCF
    ! if this fluctuates you are in deep trouble 

    real(dp) :: cof(:,0:,:,:)
    real(dp), intent(in) :: s(0:,:,:),occ(:,0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),problemsize
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq
    real(dp) :: electron_number
    real(dp) :: scaling

    ! get actual number per shell by multiplying dens matrix and overlap
    do mm=1,2
      do ii=0,max_l
        do qq=1,problemsize
          electron_number=0.0d0
          ll=0
          do jj=1,num_alpha(ii)
            do kk=1,poly_order(ii)
              ll=ll+1
              pp=0
              do nn=1,num_alpha(ii)
                do oo=1,poly_order(ii)
                  pp=pp+1

                  electron_number=electron_number+&
                      &occ(mm,ii,qq)*cof(mm,ii,ll,qq)*cof(mm,ii,pp,qq)*s(ii,ll,pp)

                end do
              end do
            end do
          end do

          if (abs(occ(mm,ii,qq)-electron_number)>1.0d-8) then
            write(*,*) 'Electron number fluctuation',&
                &occ(mm,ii,qq)-electron_number
          end if

        end do
      end do
    end do

  end subroutine check_electron_number

  function vector_length(vector,size)

    real(dp) :: vector_length
    real(dp), intent(in) :: vector(:)
    integer, intent(in) :: size
    integer :: ii

    vector_length=0.0d0

    do ii=1,size
      vector_length=vector_length+vector(ii)*vector(ii)
    end do

    vector_length=sqrt(vector_length)

  end function vector_length

  pure FUNCTION fak(n)
    REAL(dp) :: fak
    INTEGER, intent(in) :: n
    INTEGER :: h
    fak = 1.0_dp
    DO h = 1,n 
      fak = fak*real(h,dp)
    END DO
    RETURN
  END function fak
  !
  SUBROUTINE polcart(r,zeta,phi,vec)
    ! zeta=cos(theta)
    REAL(dp) :: vec(3),r,zeta,phi,s_teta
    s_teta=SQRT(1.0_dp-zeta*zeta)
    vec(3)=r*zeta
    vec(2)=r*s_teta*SIN(phi)
    vec(1)=r*s_teta*COS(phi)
    RETURN
  END subroutine polcart
  !
  SUBROUTINE cartpol(vec1,vec2,vec3,r,zeta,phi)
    REAL(dp) :: eps, tol
    PARAMETER ( eps=1.d-8 )
    PARAMETER ( tol=1.d-8 )
    REAL(dp) :: vec(3), r, zeta, phi,vec1,vec2,vec3
    !c      external dscalar
    vec(1)=vec1
    vec(2)=vec2
    vec(3)=vec3
    r = SQRT(dscalar(vec,vec))
    IF(((ABS(vec(1)).LT.eps).AND.(ABS(vec(2)).LT.eps))) &
        & phi = 0.0_dp
    IF(.NOT.((ABS(vec(1)).LT.eps).AND.(ABS(vec(2)).LT.eps))) &
        & phi = ATAN2(vec(2),vec(1))
    IF((ABS(r).LT.eps)) zeta = 0.0_dp
    IF(.NOT.(ABS(r).LT.eps)) zeta = vec(3)/r
    RETURN
  END subroutine cartpol
  !
  FUNCTION dscalar(r1,r2)
    REAL(dp) :: r1(3), r2(3), dscalar
    dscalar = r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)
    RETURN
  END function dscalar

  !> Auxilliary function, see eqn. 20 of Roothan https://link.aps.org/doi/10.1103/RevModPhys.32.186
  pure function v_int(x,i) ! V_{i}(x)

    real(dp), intent(in) :: x
    integer, intent(in) :: i
    real(dp) :: v_int

    v_int=fak(i)/(x**(i+1))

  end function v_int

  ! Auxilliary function, continuation of factorial function
  pure function v_real(x,i) ! V_{i}(x)

    real(dp), intent(in) :: x, i
    real(dp) :: v_real

    v_real=gamma(i+1)/(x**(i+1))

  end function v_real

end module utilities

