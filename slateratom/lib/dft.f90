module dft

  use, intrinsic :: iso_c_binding, only : c_size_t
  use common_accuracy, only : dp
  use common_constants
  use density
  use integration
  use xc_f90_lib_m

  implicit none
  private

  public :: dft_start_pot, density_grid, dft_exc_energy, dft_vxc_energy
  public :: dft_exc_matrixelement, xalpha, pbe_driver
  public :: check_accuracy
  public :: derive, radial_divergence, derive1_5, derive2_5

contains

  subroutine dft_start_pot(abcissa,num_mesh_points,nuc,vxc)

    ! Total potential to initialize a DFT calculation from Thomas-Fermi
    ! Theory. this does not work as intended in the current code, since
    ! we do not have a numerical Coulomb-Potential.

    ! Generalized Thomas-Fermi atomic potential
    ! as published by R. Latter, Phys. Rev. 99, 510 (1955).
    ! and implemented in Dirk Porezags scfatom

    real(dp), intent(in) :: abcissa(:)
    integer, intent(in) :: nuc,num_mesh_points
    real(dp), intent(out) :: vxc(:,:)
    real(dp) :: b,t,x,rtx
    integer :: ii

    b= (0.69395656d0/float(nuc))**(1.0d0/3.0d0)

    do ii=1,num_mesh_points

      x= abcissa(ii)/b
      rtx= sqrt(x)

      t= float(nuc)/(1.0d0+rtx*(0.02747d0-x*(0.1486d0-0.007298d0*x))&
          &+x*(1.243d0+x*(0.2302d0+0.006944d0*x)));
      if (t < 1.0d0) t= 1.0d0

      vxc(ii,1)= (t/abcissa(ii))/2.0d0
      vxc(ii,2)= (t/abcissa(ii))/2.0d0

    end do

  end subroutine dft_start_pot

  subroutine density_grid(p,max_l,num_alpha,poly_order,alpha,num_mesh_points,&
      &abcissa, dzdr, d2zdr2, dz, xcnr, rho,drho,ddrho,vxc,exc,xalpha_const)

    ! Calculate and store density and density derivatives on radial grid.
    ! Also calculate and store exchange-correlation potential and energy
    ! density on grid.

    real(dp), intent(in) :: p(:,0:,:,:),abcissa(:),alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),num_mesh_points
    real(dp), intent(in) :: dzdr(:), d2zdr2(:)
    real(dp), intent(in) :: dz,xalpha_const
    integer, intent(in) :: xcnr
    real(dp), intent(out) :: rho(:,:),drho(:,:),ddrho(:,:),vxc(:,:),exc(:)
    real(dp) :: rhotot,rhodiff,drhotot,ddrhotot,drhodiff,ddrhodiff
    integer :: ii,jj,kk,ll,mm,oo
    integer(c_size_t) :: nn
    type(xc_f90_func_t) :: xcfunc_x, xcfunc_c
    type(xc_f90_func_info_t) :: xcinfo
    real(dp), allocatable :: tmprho(:,:), ex(:), ec(:), vx(:,:), vc(:,:)
    real(dp), allocatable :: tmpsigma(:,:), vxsigma(:,:), vcsigma(:,:)
    real(dp), allocatable :: tmpv(:), tmpv2(:)
    integer :: ispin, ispin2, isigma
    real(dp), parameter :: rec4pi = 1.0_dp / (4.0_dp * pi)


    if (xcnr==0) return
    if (xcnr == 2) then
      call xc_f90_func_init(xcfunc_x, XC_LDA_X, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
      call xc_f90_func_init(xcfunc_c, XC_LDA_C_PW, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
    elseif (xcnr == 3) then
      call xc_f90_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
      call xc_f90_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
    end if

    do ii=1,num_mesh_points

      rho(ii,1)=density_at_point(p(1,:,:,:),max_l,num_alpha,poly_order,alpha,&
          &abcissa(ii))
      rho(ii,2)=density_at_point(p(2,:,:,:),max_l,num_alpha,poly_order,alpha,&
          &abcissa(ii))

    end do

    rho = max(rho, 0.0_dp)
    !rho(:,:) = sign(max(abs(rho), 1e-14_dp), rho)
    !drho(ii,:) = sign(max(abs(drho(ii,:)), 1e-14_dp), drho(ii,:))
    !ddrho(ii,:) = sign(max(abs(ddrho(ii,:)), 1e-14_dp), ddrho(ii,:))
    !if (abs(rho(ii,1))<1.0d-16) rho(ii,1)=0.0d0
    !if (abs(rho(ii,2))<1.0d-16) rho(ii,2)=0.0d0
    !if (abs(drho(ii,1))<1.0d-16) drho(ii,1)=0.0d0
    !if (abs(drho(ii,2))<1.0d-16) drho(ii,2)=0.0d0
    !if (abs(ddrho(ii,1))<1.0d-16) ddrho(ii,1)=0.0d0
    !if (abs(ddrho(ii,2))<1.0d-16) ddrho(ii,2)=0.0d0

    if (xcnr > 2) then

      !call derive2_5(rho(:,1), dz, ddrho(:,1), dzdr, d2zdr2, drho(:,1))
      !call derive2_5(rho(:,2), dz, ddrho(:,2), dzdr, d2zdr2, drho(:,2))

      do ii = 1, num_mesh_points

        drho(ii,1)=density_at_point_1st(p(1,:,:,:),max_l,num_alpha,poly_order,&
            &alpha,abcissa(ii))
        drho(ii,2)=density_at_point_1st(p(2,:,:,:),max_l,num_alpha,poly_order,&
            &alpha,abcissa(ii))

        ddrho(ii,1)=density_at_point_2nd(p(1,:,:,:),max_l,num_alpha,poly_order,&
            &alpha,abcissa(ii))
        ddrho(ii,2)=density_at_point_2nd(p(2,:,:,:),max_l,num_alpha,poly_order,&
            &alpha,abcissa(ii))
      end do

    end if

    ! divide by 4*pi to catch different normalization of spherical harmonics
    if (xcnr==1) then
      do ii = 1, num_mesh_points
        rhotot = (rho(ii,1) + rho(ii,2)) * rec4pi
        rhodiff = (rho(ii,1) - rho(ii,2)) * rec4pi
        call xalpha(rhotot,rhodiff,vxc(ii,:),exc(ii),xalpha_const)
      end do

    else if (xcnr==2) then
      nn = size(rho, dim=1)
      allocate(tmprho(2, nn))
      allocate(ex(nn))
      allocate(ec(nn))
      allocate(vx(2, nn))
      allocate(vc(2, nn))
      tmprho(:,:) = transpose(rho) * rec4pi
      call xc_f90_lda_exc_vxc(xcfunc_x, nn, tmprho(1,1), ex(1), vx(1,1))
      call xc_f90_lda_exc_vxc(xcfunc_c, nn, tmprho(1,1), ec(1), vc(1,1))
      vxc(:,:) = transpose(vx + vc)
      exc = ec + ex
!!! OLD hand coded XC version
!      do ii = 1, num_mesh_points
!        rhotot = (rho(ii,1) + rho(ii,2)) * rec4pi
!        rhodiff = (rho(ii,1) - rho(ii,2)) * rec4pi
!        call pbe_driver(0,rhotot,0.0d0,0.0d0,&
!            &rhodiff,0.0d0,0.0d0,0.0d0,vxc(ii,:),exc(ii))
!      end do
!!!
    else if (xcnr==3) then
      nn = size(rho, dim=1)
      allocate(tmprho(2, nn))
      allocate(ex(nn))
      allocate(ec(nn))
      allocate(vx(2, nn))
      allocate(vc(2, nn))
      allocate(tmpsigma(3, nn))
      allocate(vxsigma(3, nn))
      allocate(vcsigma(3, nn))
      allocate(tmpv(nn))
      allocate(tmpv2(nn))
      tmprho(:,:) = transpose(rho) * rec4pi
      tmpsigma(1,:) = drho(:,1) * drho(:,1) * rec4pi * rec4pi
      tmpsigma(2,:) = drho(:,1) * drho(:,2) * rec4pi * rec4pi
      tmpsigma(3,:) = drho(:,2) * drho(:,2) * rec4pi * rec4pi
      call xc_f90_gga_exc_vxc(xcfunc_x, nn, tmprho(1,1), tmpsigma(1,1),&
          & ex(1), vx(1,1), vxsigma(1,1))
      call xc_f90_gga_exc_vxc(xcfunc_c, nn, tmprho(1,1), tmpsigma(1,1), ec(1), &
          &vc(1,1), vcsigma(1,1))
      vxc = transpose(vx + vc)
      do ispin = 1, 2
        ispin2 = 3 - ispin           ! the other spin
        isigma = 2 * ispin - 1       ! 1 for spin up, 3 for spin down
        tmpv(:) = (vxsigma(isigma,:) + vcsigma(isigma,:)) &
            & * drho(:,ispin) * rec4pi
        call radial_divergence(tmpv, abcissa, dz, tmpv2, dzdr)
        vxc(:,ispin) = vxc(:,ispin) - 2.0_dp * tmpv2
        tmpv(:) = (vxsigma(2,:) +  vcsigma(2,:)) &
            & * drho(:,ispin2) * rec4pi
        call radial_divergence(tmpv, abcissa, dz, tmpv2, dzdr)
        vxc(:,ispin) = vxc(:,ispin) - tmpv2
      end do
      exc = ex + ec
!!! OLD: hand coded xc-version
!      do ii = 1, num_mesh_points
!        rhotot = (rho(ii,1) + rho(ii,2)) * rec4pi
!        rhodiff = (rho(ii,1) - rho(ii,2)) * rec4pi
!        drhotot=(drho(ii,1)+drho(ii,2))/4.0d0/pi
!        ddrhotot=(ddrho(ii,1)+ddrho(ii,2))/4.0d0/pi
!        drhodiff=(drho(ii,1)-drho(ii,2))/4.0d0/pi
!        ddrhodiff=(ddrho(ii,1)-ddrho(ii,2))/4.0d0/pi
!        call pbe_driver(1,rhotot,drhotot,ddrhotot,&
!            &rhodiff,drhodiff,ddrhodiff,abcissa(ii),vxc(ii,:),exc(ii))
!      end do
!!!

    else

      write(*,'(A,I2,A)') 'XCNR= ',xcnr,' not implemented'
      STOP

    end if


    call xc_f90_func_end(xcfunc_x)
    call xc_f90_func_end(xcfunc_c)

  end subroutine density_grid

  subroutine dft_exc_energy(num_mesh_points,rho,exc,weight,abcissa,&
      &xcnr,exc_energy)

    ! Calculate DFT Exc energy from energy density and electron density on
    ! grid.

    real(dp),intent(out) :: exc_energy
    real(dp), intent(in) :: rho(:,:),weight(:),exc(:),abcissa(:)
    integer, intent(in) :: num_mesh_points,xcnr
    integer :: ii,jj,kk,ll,mm,nn,oo
    real(dp) :: rhotot,rhodiff

    exc_energy=0.0d0

    do ii=1,num_mesh_points

      exc_energy=exc_energy+weight(ii)*exc(ii)*(rho(ii,1)+rho(ii,2))*&
          &abcissa(ii)**2

    end do

    !
    ! For usual DFT functionals E_xc=\int \rho \eps(\rho,\zeta) d^3r
    ! so there is only one exchange-correlation energy density \eps(\rho,\zeta) and
    ! exc_energy could be a scalar without problems.
    !

  end subroutine dft_exc_energy

  subroutine dft_vxc_energy(num_mesh_points,rho,vxc,weight,abcissa,&
      &xcnr,vxc_energy)
    ! vxc contribution for double counting correction

    real(dp),intent(out) :: vxc_energy(2)
    real(dp), intent(in) :: rho(:,:),weight(:),vxc(:,:),abcissa(:)
    integer, intent(in) :: num_mesh_points,xcnr
    integer :: ii,jj,kk,ll,mm,nn,oo
    real(dp) :: rhotot,rhodiff

    vxc_energy=0.0d0

    do ii=1,num_mesh_points

      vxc_energy(1)=vxc_energy(1)+weight(ii)*vxc(ii,1)*(rho(ii,1))*&
          &abcissa(ii)**2
      vxc_energy(2)=vxc_energy(2)+weight(ii)*vxc(ii,2)*(rho(ii,2))*&
          &abcissa(ii)**2

    end do


  end subroutine dft_vxc_energy

  subroutine dft_exc_matrixelement(num_mesh_points,weight,abcissa,rho,vxc,&
      &xcnr,alpha1,poly1,alpha2,poly2,l,exc_matrixelement)

    ! Calculate a single matrix element of the exchange correlation potential.

    real(dp),intent(out) :: exc_matrixelement(2)
    real(dp), intent(in) :: weight(:),abcissa(:),rho(:,:),vxc(:,:)
    real(dp), intent(in) :: alpha1,alpha2
    integer, intent(in) :: num_mesh_points,xcnr
    integer, intent(in) :: poly1,poly2,l
    real(dp) :: basis
    integer :: ii,jj,kk,ll,mm,nn,oo

    exc_matrixelement=0.0d0

    do ii=1,num_mesh_points

      basis=basis_times_basis_times_r2(alpha1,poly1,alpha2,poly2,l,abcissa(ii))

      exc_matrixelement(1)=exc_matrixelement(1)-weight(ii)*vxc(ii,1)*basis

      exc_matrixelement(2)=exc_matrixelement(2)-weight(ii)*vxc(ii,2)*basis

    end do


  end subroutine dft_exc_matrixelement

  subroutine xalpha(rhotot,rhodiff,vxc,exc,alpha)

    ! Xalpha potential and energy density.

    ! alpha=2/3 recovers the Gaspar/Kohn/Sham functional commonly used as
    ! exchange part in most current LSDA and GGA functionals
    ! the original Slater exchange is recoverd with alpha=1

    real(dp), intent(in) :: rhotot,rhodiff,alpha
    real(dp), intent(out) :: exc,vxc(2)
    real(dp) :: third,fourthird,vfac,cx,fzeta,dfzeta,eps0,eps1,spinpart,zeta

    third=1.0d0/3.0d0
    fourthird=4.0d0/3.0d0
    vfac=2.0d0**third
    cx=0.75d0*(3.d0/pi)**third

    if (abs(rhotot)<1.0d-12) then
      exc=0.0d0
      vxc(1)=0.0d0
      vxc(2)=0.0d0
      return
    end if

    zeta=rhodiff/rhotot

    if (abs(zeta)>1.0d12) write(*,*) 'ZETA LARGE IN X-ALPHA'

    fzeta=((1.0d0+zeta)**fourthird+(1.0d0-zeta)**fourthird-2.0d0)/(2.0d0*(vfac-1.0d0))
    dfzeta=fourthird*((1.0d0+zeta)**third-(1.0d0-zeta)**third)/(2.0d0*(vfac-1.0d0))

    eps0=-3.0d0/2.0d0*alpha*cx*rhotot**third
    eps1=vfac*eps0

    exc=eps0+(eps1-eps0)*fzeta

    spinpart=(eps1-eps0)*dfzeta*(1.0d0-zeta)

    vxc(1)=fourthird*exc+spinpart
    vxc(2)=fourthird*exc-spinpart

  end subroutine xalpha

  subroutine pbe_driver(xcnr,rho,drho,ddrho,zeta,dzeta,ddzeta,r,vxc,exc)

    ! Driver for the PBE routines. Note: this does a lot of Voodoo but seems
    ! to work.

    integer, intent(in) :: xcnr
    real(dp), intent(in) :: rho,drho,ddrho,zeta,dzeta,ddzeta,r
    real(dp), intent(out) :: vxc(2),exc
    real(dp) :: z,dz,rs,alfa,gg,t,u,v,w,vc(2),rho1,rho2,drho1,drho2,ec
    real(dp) :: ddrho1,ddrho2,rs1,rs2,s1,s2,u1,u2,eps,t1,t2,ex(2),vx(2)
    integer :: igga,idft

    igga=xcnr

    if(abs(rho).lt.1.d-14)then
      vxc(1)=0.0d0
      vxc(2)=0.0d0
      exc=0.0d0
      return
    endif

    ! FROM BURKEs FORTRAN REFERENCE SOURCE
    !
    ! Now do correlation
    ! zet=(up-dn)/rho
    ! g=phi(zeta)
    ! rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
    ! sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
    ! twoksg=2*Ks*phi
    ! t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
    ! uu=delgrad/(rho^2*twoksg^3)
    ! rholap=Laplacian
    ! vv=Laplacian/(rho*twoksg^2)
    ! ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
    ! ec=lsd correlation energy
    ! vcup=lsd up correlation potential
    ! vcdn=lsd down correlation potential
    ! h=gradient correction to correlation energy
    ! dvcup=gradient correction to up correlation potential
    ! dvcdn=gradient correction to down correlation potential

    alfa=(4.d0/(9.d0*pi))**(1.d0/3.d0)

    eps=1.d0-1.0d-12
    z=zeta/rho
    if(z.ge. eps) z= eps
    if(z.le.-eps) z=-eps
    dz=(dzeta*rho-zeta*drho)/rho**2

    rs=(4.d0*pi*rho/3.d0)**(-1.d0/3.d0)
    gg=((1+z)**(2.d0/3.d0)+(1-z)**(2.d0/3.d0))/2.d0
    t=dabs(drho)/rho*dsqrt(pi/4.d0*alfa*rs)/(2.d0*gg)
    u=dabs(drho)*ddrho/(rho**2)*(dsqrt(pi/4.d0*alfa*rs)/(2.d0*gg))**3
    v=(ddrho+2.d0/r*drho)/rho * (dsqrt(pi/4.d0*alfa*rs)/(2.d0*gg))**2
    w=drho*dz/rho*(dsqrt(pi/4.d0*alfa*rs)/(2.d0*gg))**2
    call correlation(rs,z,t,u,v,w,igga,ec,vc(1),vc(2))

    ! rho1=up electron desnity
    rho1  =(rho+zeta)/2.d0

    ! rho2=down electron desnity
    rho2  =(rho-zeta)/2.d0

    ! derivatives
    drho1 =(drho+dzeta)/2.d0
    drho2 =(drho-dzeta)/2.d0
    ddrho1=(ddrho+ddzeta)/2.d0
    ddrho2=(ddrho-ddzeta)/2.d0

    ! FROM BURKEs FORTRAN REFERENCE SOURCE
    !
    ! PBE exchange
    ! use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
    ! do up exchange
    ! fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3)
    ! s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
    ! u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
    ! v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
    !
    ! Wigner-Seitz Radii of up (rs1) and down (rs2) electrons
    !
    ! actually this should be rs=((4*pi*rho)/3)**(-1/3)
    ! but s is calculated correctly later compared to Burkes comments

    rs1=(8.d0*pi*rho1/3.d0)**(-1.d0/3.d0)
    rs2=(8.d0*pi*rho2/3.d0)**(-1.d0/3.d0)

    ! alfa=(4.d0/(9.d0*pi))**(1.d0/3.d0)

    s1=dabs(drho1)*(alfa*rs1/2.d0)/rho1
    s2=dabs(drho2)*(alfa*rs2/2.d0)/rho2
    u1=dabs(drho1)*ddrho1/(rho1**2)*(alfa*rs1/2.d0)**3
    u2=dabs(drho2)*ddrho2/(rho2**2)*(alfa*rs2/2.d0)**3
    t1=(ddrho1+2.d0/r*drho1)/rho1*(alfa*rs1/2.d0)**2
    t2=(ddrho2+2.d0/r*drho2)/rho2*(alfa*rs2/2.d0)**2

    !
    ! use 2.d0*rho1 and 2.d0*rho2 because of spin scaling, see Burkes comment
    !
    call exchange(2.d0*rho1,s1,u1,t1,igga,ex(1),vx(1))
    call exchange(2.d0*rho2,s2,u2,t2,igga,ex(2),vx(2))

    exc=0.5d0*((1.d0+z)*ex(1)+(1.d0-z)*ex(2))+ec

    vxc(1)=vx(1)+vc(1)
    vxc(2)=vx(2)+vc(2)

  end subroutine pbe_driver

  SUBROUTINE CORRELATION(RS,ZET,T,UU,VV,WW,igga,ec,vc1,vc2)

    !
    ! APART FROM COSMETICS THIS IS IN FACT BURKEs FORTRAN REFERENCE IMPLEMENTATION
    !

    ! This is the PBE and PW-LDA Correlation routine.

    IMPLICIT REAL*8 (A-H,O-Z)
    !----------------------------------------------------------------------
    !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
    !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
    !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
    !       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
    !       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
    !       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
    !       : UU,VV,WW, only needed for PBE potential
    !       : igga=flag to do gga (0=>LSD only)
    !  output: ecl=lsd correlation energy from [a]
    !        : ecn=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
    !        : vcup=lsd up correlation potential
    !        : vcdn=lsd dn correlation potential
    !        : dvcup=nonlocal correction to vcup
    !        : dvcdn=nonlocal correction to vcdn
    !----------------------------------------------------------------------
    ! References:
    ! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof,
    !     {\sl Generalized gradient approximation made simple}, sub.
    !     to Phys. Rev.Lett. May 1996.
    ! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
    !     construction of a generalized gradient approximation:  The PW91
    !     density functional}, submitted to Phys. Rev. B, Feb. 1996.
    ! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
    !----------------------------------------------------------------------
    !     bet=coefficient in gradient expansion for correlation, [a](4).
    integer :: igga
    parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
    parameter(GAM=0.5198420997897463295344212145565d0)
    parameter(thrd4=4.d0*thrd, fzz=8.d0/(9.d0*GAM))
    parameter(gamma=0.03109069086965489503494086371273d0)
    parameter(bet=0.06672455060314922d0,delt=bet/gamma)
    dimension u(6),p(6),s(6)
    data u/ 0.03109070D0, 0.2137000D0, 7.5957000D0,&
        &        3.58760000D0, 1.6382000D0, 0.4929400D0/
    data p/ 0.01554535D0, 0.2054800D0,14.1189000D0,&
        &        6.19770000D0, 3.3662000D0, 0.6251700D0/
    data s/ 0.01688690D0, 0.1112500D0,10.3570000D0,&
        &        3.62310000D0, 0.8802600D0, 0.4967100D0/
    !----------------------------------------------------------------------
    !     find LSD energy contributions, using [c](10) .
    !     EU=unpolarized LSD correlation energy , EURS=dEU/drs
    !     EP=fully polarized LSD correlation energy , EPRS=dEP/drs
    !     ALFM=-spin stiffness, [c](3) , ALFRSM=-dalpha/drs .
    !     F=spin-scaling factor from [c](9).
    !     construct ecl, using [c](8) .
    !

    rtrs=dsqrt(rs)
    Q0 = -2.D0*u(1)*(1.D0+u(2)*rtrs*rtrs)
    Q1 = 2.D0*u(1)*rtrs*(u(3)+rtrs*(u(4)+rtrs*(u(5)+u(6)*rtrs)))
    Q2 = DLOG(1.D0+1.D0/Q1)
    Q3 = u(1)*(u(3)/rtrs+2.D0*u(4)+rtrs*(3.D0*u(5)+4.D0*u(6)*rtrs))
    EU = Q0*Q2
    EURS = -2.D0*u(1)*u(2)*Q2-Q0*Q3/(Q1*(1.d0+Q1))
    Q0 = -2.D0*p(1)*(1.D0+p(2)*rtrs*rtrs)
    Q1 = 2.D0*p(1)*rtrs*(p(3)+rtrs*(p(4)+rtrs*(p(5)+p(6)*rtrs)))
    Q2 = DLOG(1.D0+1.D0/Q1)
    Q3 = p(1)*(p(3)/rtrs+2.D0*p(4)+rtrs*(3.D0*p(5)+4.D0*p(6)*rtrs))
    EP = Q0*Q2
    EPRS = -2.D0*p(1)*p(2)*Q2-Q0*Q3/(Q1*(1.d0+Q1))
    Q0 = -2.D0*s(1)*(1.D0+s(2)*rtrs*rtrs)
    Q1 = 2.D0*s(1)*rtrs*(s(3)+rtrs*(s(4)+rtrs*(s(5)+s(6)*rtrs)))
    Q2 = DLOG(1.D0+1.D0/Q1)
    Q3 = s(1)*(s(3)/rtrs+2.D0*s(4)+rtrs*(3.D0*s(5)+4.D0*s(6)*rtrs))
    ALFM = Q0*Q2
    ALFRSM = -2.D0*s(1)*s(2)*Q2-Q0*Q3/(Q1*(1.d0+Q1))

    Z4 = ZET**4
    F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
    ECL= EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
    !----------------------------------------------------------------------
    !     LSD potential from [c](A1)
    !     ECRS = dEc/drs , ECZET=dEc/dzeta , FZ = dF/dzeta   [c](A2-A4)
    !
    ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
    FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
    ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)&
        &        +FZ*(Z4*EP-Z4*EU-(1.D0-Z4)*ALFM/FZZ)
    COMM = ECL -RS*ECRS/3.D0-ZET*ECZET
    VCUP = COMM + ECZET
    VCDN = COMM - ECZET
    if(igga.eq.0)then
      EC=ECL
      VC1=VCUP
      VC2=VCDN
      return
    endif
    !----------------------------------------------------------------------
    !     PBE correlation energy
    !     G=phi(zeta), given after [a](3)
    !     DELT=bet/gamma , B=A of [a](8)
    !
    G=((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
    G3 = G**3
    PON=-ECL/(G3*gamma)
    B = DELT/(DEXP(PON)-1.D0)
    B2 = B*B
    T2 = T*T
    T4 = T2*T2
    Q4 = 1.D0+B*T2
    Q5 = 1.D0+B*T2+B2*T4
    ECN= G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
    EC = ECL + ECN
    !----------------------------------------------------------------------
    !     ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
    !
    G4 = G3*G
    T6 = T4*T2
    RSTHRD = RS/3.D0
    !      GZ=((1.d0+zet)**thirdm-(1.d0-zet)**thirdm)/3.d0
    ! ckoe: hack thirdm never gets defined, but 1-1 should be zero anyway
    GZ=0.0d0
    FAC = DELT/B+1.D0
    BG = -3.D0*B2*ECL*FAC/(BET*G4)
    BEC = B2*FAC/(BET*G3)
    Q8 = Q5*Q5+DELT*Q4*Q5*T2
    Q9 = 1.D0+2.D0*B*T2
    hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
    hRS = -RSTHRD*hB*BEC*ECRS
    FACT0 = 2.D0*DELT-6.D0*B
    FACT1 = Q5*Q9+Q4*Q9*Q9
    hBT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
    hRST = RSTHRD*T2*hBT*BEC*ECRS
    hZ = 3.D0*GZ*ecn/G + hB*(BG*GZ+BEC*ECZET)
    hT = 2.d0*BET*G3*Q9/Q8
    hZT = 3.D0*GZ*hT/G+hBT*(BG*GZ+BEC*ECZET)
    FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
    FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
    hTT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
    COMM = ECN+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
    PREF = HZ-GZ*T2*HT/G
    FACT5 = GZ*(2.D0*HT+T*HTT)/G
    COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
    DVCUP = COMM + PREF
    DVCDN = COMM - PREF
    VC1 = VCUP + DVCUP
    VC2 = VCDN + DVCDN
    !	print*,'c igga is',dvcup

  END subroutine correlation

  subroutine exchange(rho,s,u,t,igga,EX,VX)

    ! APART FROM COSMETICS THIS IS IN FACT BURKEs FORTRAN REFERENCE IMPLEMENTATION

    ! This is the PBE and PW-LDA Exchange routine.

    implicit integer*4 (i-n)
    implicit real*8 (a-h,o-z)

    parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
    parameter(pi=3.14159265358979323846264338327950d0)
    parameter(ax=-0.738558766382022405884230032680836d0)

    parameter(um=0.21951d0,uk=0.8040d0,ul=um/uk)

    parameter(ap=1.647127d0,bp=0.980118d0,cp=0.017399d0)
    parameter(aq=1.523671d0,bq=0.367229d0,cq=0.011282d0)
    parameter(ah=0.19645d0,bh=7.7956d0)
    parameter(ahp=0.27430d0,bhp=0.15084d0,ahq=0.004d0)
    parameter(a1=0.19645d0,a2=0.27430d0,a3=0.15084d0,a4=100.d0)
    parameter(a=7.79560d0,b1=0.004d0,eps=1.d-15)

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !  GGA EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
    !----------------------------------------------------------------------
    !  INPUT rho : DENSITY
    !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
    !  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
    !  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)  (for U,V, see PW86(24))
    !  input igga:  (=0=>don't put in gradient corrections, just LDA)
    !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (LOCAL: EXL, NONLOCAL: EXN,
    !           TOTAL: EX) AND POTENTIAL (VX)
    !----------------------------------------------------------------------
    ! References:
    ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
    ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
    !     {\bf 40},  3399  (1989) (E).
    !----------------------------------------------------------------------
    ! Formulas: e_x[unif]=ax*rho^(4/3)  [LDA]
    !           ax = -0.75*(3/pi)^(1/3)
    !	    e_x[PBE]=e_x[unif]*FxPBE(s)
    !	    FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
    !           uk, ul defined after [a](13)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !     construct LDA exchange energy density

    exunif = ax*rho**thrd
    if((igga.eq.0).or.(s.lt.eps))then
      EXL=exunif
      EXN=0.d0
      EX=EXL+EXN
      VX= exunif*thrd4
      return
    endif
    !----------------------------------------------------------------------
    !     construct GGA enhancement factor
    !     find first and second derivatives of f and:
    !     fs=(1/s)*df/ds  and  fss=dfs/ds = (d2f/ds2 - (1/s)*df/ds)/s

    !
    ! PBE enhancement factors checked against NRLMOL
    !
    if(igga.eq.1)then
      p0 =1.d0+ul*s**2
      f  =1.d0+uk-uk/p0
      fs =2.d0*uk*ul/p0**2
      fss=-4.d0*ul*s*fs/p0
    endif

    !

    EXL= exunif
    EXN= exunif*(f-1.0d0)
    EX = EXL+EXN
    !----------------------------------------------------------------------
    !     energy done. calculate potential from [b](24)
    !
    VX = exunif*(thrd4*f-(u-thrd4*s**3)*fss-t*fs )
    !	 print*,'e igga is',igga,vx,xunif*thrd4


  END subroutine exchange

  subroutine check_accuracy(weight,abcissa,num_mesh_points,max_l,&
      &num_alpha,alpha,poly_order)

    ! Test integration to check the accuracy of the radial mesh by
    ! integrating the square of a primitive Slater basis function which are
    ! analytically normalized to 1.0d0 !

    real(dp), intent(in) :: weight(:),abcissa(:),alpha(0:,:)
    integer, intent(in) :: num_mesh_points,max_l,num_alpha(0:),poly_order(0:)
    real(dp) :: value
    integer :: ii,jj,kk,ll

    do ii=0,max_l
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          value=0.0d0
          do ll=1,num_mesh_points

            value=value+weight(ll)*abcissa(ll)**2*&
                &basis(alpha(ii,jj),kk,ii,abcissa(ll))**2

          end do
          if (abs(1.0d0-value)>1.0d-12) then
            write(*,'(A,F12.6,I3,E12.3)') 'WARNING: Integration bad for basis &
                &function ',alpha(ii,jj),kk+ii-1,abs(1.0d0-value)
            write(*,'(A)') 'Accuracy is not better than 1.0d-12'
          end if
        end do
      end do
    end do

  end subroutine check_accuracy


  subroutine radial_divergence(ff, rr, dr, rdiv, jacobi)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: rr(:)
    real(dp), intent(in) :: dr
    real(dp), intent(out) :: rdiv(:)
    real(dp), intent(in), optional :: jacobi(:)

    call derive1_5(ff, dr, rdiv, jacobi)
    rdiv = rdiv + 2.0_dp / rr * ff

  end subroutine radial_divergence


  subroutine derive(ff, dx, jacobi)
    real(dp), intent(inout) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(in), optional :: jacobi(:)

    real(dp), allocatable :: tmp1(:)
    integer :: nn

    nn = size(ff)
    allocate(tmp1(nn))
    tmp1(:) = ff
    ff(2:nn-1) = (ff(3:nn) - ff(1:nn-2)) / (2.0 * dx)
    ff(1) = (tmp1(2) - tmp1(1)) / dx
    ff(nn) = (tmp1(nn) - tmp1(nn-1)) / dx
    if (present(jacobi)) then
      ff = ff * jacobi
    end if

  end subroutine derive


  subroutine derive1_5(ff, dx, dfdx, dudx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: dfdx(:)
    real(dp), intent(in), optional :: dudx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) = &
        reshape([ &
        &-25.0_dp,  48.0_dp, -36.0_dp,  16.0_dp,  -3.0_dp, &
        & -3.0_dp, -10.0_dp,  18.0_dp,  -6.0_dp,   1.0_dp, &
        &  1.0_dp, -8.0_dp,   0.0_dp,    8.0_dp,  -1.0_dp, &
        & -1.0_dp,  6.0_dp,  -18.0_dp,  10.0_dp,   3.0_dp, &
        &  3.0_dp, -16.0_dp,  36.0_dp, -48.0_dp,  25.0_dp ], [ np, np ])

    integer :: ngrid
    integer :: ii

    ngrid = size(ff)
    do ii = 1, nleft
      dfdx(ii) = dot_product(coeffs(:,ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      dfdx(ii) = dot_product(coeffs(:,imiddle), ff(ii-nleft:ii+nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      dfdx(ii) = dot_product(coeffs(:,np-(ngrid-ii)), ff(ngrid-np+1:ngrid))
    end do

    if (present(dudx)) then
      dfdx = dfdx * (dudx /  (dxprefac * dx))
    else
      dfdx = dfdx / (dxprefac * dx)
    end if

  end subroutine derive1_5



  subroutine derive2_5(ff, dx, d2fdx2, dudx, d2udx2, dfdx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: d2fdx2(:)
    real(dp), intent(in), optional :: dudx(:), d2udx2(:)
    real(dp), intent(out), target, optional :: dfdx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) = &
        reshape([ &
        &  35.0_dp, -104.0_dp,  114.0_dp,  -56.0_dp,   11.0_dp, &
        &  11.0_dp,  -20.0_dp,    6.0_dp,    4.0_dp,   -1.0_dp, &
        &  -1.0_dp,  16.0_dp,   -30.0_dp,   16.0_dp,   -1.0_dp, &
        &  -1.0_dp,   4.0_dp,     6.0_dp,  -20.0_dp,   11.0_dp, &
        &  11.0_dp,  -56.0_dp,  114.0_dp, -104.0_dp,   35.0_dp ], [ np, np ])

    integer :: ngrid
    integer :: ii
    real(dp), allocatable, target :: dfdxlocal(:)
    real(dp), pointer :: pdfdx(:)

    ngrid = size(ff)
    if (present(dfdx)) then
      pdfdx => dfdx
    elseif (present(d2udx2)) then
      allocate(dfdxlocal(ngrid))
      pdfdx => dfdxlocal
    end if

    do ii = 1, nleft
      d2fdx2(ii) = dot_product(coeffs(:,ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      d2fdx2(ii) = dot_product(coeffs(:,imiddle), ff(ii-nleft:ii+nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      d2fdx2(ii) = dot_product(coeffs(:,np-(ngrid-ii)), ff(ngrid-np+1:ngrid))
    end do

    if (present(dudx)) then
      d2fdx2 = d2fdx2 * (dudx * dudx /  (dxprefac * dx * dx))
    else
      d2fdx2 = d2fdx2 / (dxprefac * dx * dx)
    end if

    if (present(d2udx2) .or. present(dfdx)) then
      call derive1_5(ff, dx, pdfdx)
      if (present(d2udx2)) then
        d2fdx2 = d2fdx2 + pdfdx * d2udx2
      end if
      if (present(dfdx) .and. present(dudx)) then
        dfdx = dfdx * dudx
      end if
    end if

  end subroutine derive2_5


  !  subroutine grad_test(p,max_l,num_alpha,poly_order,alpha,num_mesh_points,&
  !              &abcissa,xcnr,rho,drho,ddrho,vxc,exc)
  !
  !  implicit none
  !
  !  real(dp), intent(in) :: p(:,0:,:,:),abcissa(:),alpha(0:,:)
  !  integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),num_mesh_points
  !  integer, intent(in) :: xcnr
  !  real(dp), intent(out) :: rho(:,:),drho(:,:),ddrho(:,:),vxc(:,:),exc(:)
  !  real(dp) :: rhotot,rhodiff,drhotot,ddrhotot
  !  integer :: ii,jj,kk,ll,mm,nn,oo
  !
  !  do ii=1,500
  !
  !    rho(1,ii)=density_at_point(p(1,:,:,:),max_l,num_alpha,poly_order,alpha,&
  !           &0.01d0*ii)
  !    rho(2,ii)=density_at_point(p(2,:,:,:),max_l,num_alpha,poly_order,alpha,&
  !           &0.01d0*ii)
  !
  !    drho(1,ii)=density_at_point_1st(p(1,:,:,:),max_l,num_alpha,poly_order,&
  !           &alpha,0.01d0*ii)
  !    drho(2,ii)=density_at_point_1st(p(2,:,:,:),max_l,num_alpha,poly_order,&
  !           &alpha,0.01d0*ii)
  !
  !    ddrho(1,ii)=density_at_point_2nd(p(1,:,:,:),max_l,num_alpha,poly_order,&
  !           &alpha,0.01d0*ii)
  !    ddrho(2,ii)=density_at_point_2nd(p(2,:,:,:),max_l,num_alpha,poly_order,&
  !           &alpha,0.01d0*ii)
  !
  !    write(*,'(F12.4,3F20.8)') ii*0.01d0,rho(1,ii),drho(1,ii),ddrho(1,ii)
  !  end do
  !    STOP
  !
  !  end subroutine grad_test

end module dft
