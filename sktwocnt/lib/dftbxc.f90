module dftxc

  use, intrinsic :: ieee_arithmetic
  use common_accuracy, only : dp
  use common_constants

  implicit none
  private

  public :: getxcpot_ldapw91, getxcpot_ggapbe

  real(dp), parameter :: rec4pi = 1.0_dp / (4.0_dp * pi)

contains

  subroutine getxcpot_ldapw91(rho4pi, xcpot)
    real(dp), intent(in) :: rho4pi(:)
    real(dp), intent(out) :: xcpot(:)

    integer :: nn, ii
    real(dp), allocatable :: rho(:), rs(:)
    real(dp) :: vcup, vcdn, ec, vx, ex

    nn = size(rho4pi)
    allocate(rs(nn), rho(nn))

    ! Renorm rho (incoming quantity is 4pi normed)
    rho = rho4pi * rec4pi
    ! Note: rho is normed to 4pi, therefore 4*pi missing in rs
    rs = (3.0_dp / rho4pi)**(1.0_dp / 3.0_dp)
    do ii = 1, nn
      if (rho(ii) < epsilon(1.0_dp)) then
        xcpot(ii) = 0.0_dp
      else
        call correlation_pbe(rs(ii), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            &0, ec, vcup, vcdn)
        call exchange_pbe(rho(ii), 0.0_dp, 0.0_dp, 0.0_dp, 0, ex, vx)
        xcpot(ii) = vcup + vx
      end if
    end do

    deallocate(rs, rho)
    
  end subroutine getxcpot_ldapw91



  subroutine getxcpot_ggapbe(rho4pi, absgr4pi, laplace4pi, gr_grabsgr4pi, xcpot)
    real(dp), intent(in) :: rho4pi(:)
    real(dp), intent(in) :: absgr4pi(:), laplace4pi(:), gr_grabsgr4pi(:)
    real(dp), intent(out) :: xcpot(:)

    real(dp), allocatable :: rho(:), absgr(:), laplace(:), gr_grabsgr(:)
    real(dp), allocatable :: rs(:), fac(:), tt(:), uu(:), vv(:)
    real(dp), allocatable :: ss(:), u2(:), v2(:)
    real(dp) :: alpha, zeta, gg, ww
    real(dp) :: ec, vcup, vcdn, ex, vx
    integer :: nn, ii

    nn = size(rho4pi)
    allocate(rho(nn), absgr(nn), laplace(nn), gr_grabsgr(nn))
    allocate(rs(nn), fac(nn), tt(nn), uu(nn), vv(nn), ss(nn), u2(nn), v2(nn))

    ! Renorm rho and derivatives (incoming quantities are 4pi normed)
    rho = rho4pi * rec4pi
    absgr = absgr4pi / rho4pi
    laplace = laplace4pi / rho4pi
    gr_grabsgr = gr_grabsgr4pi / rho4pi**2

    ! Note: rho is normed to 4pi, therefore 4*pi missing in rs
    rs = (3.0_dp / rho4pi)**(1.0_dp / 3.0_dp)
    zeta = 0.0_dp
    gg = 1.0_dp
    alpha = (4.0_dp/(9.0_dp * pi))**(1.0_dp/3.0_dp)
    ! Factors for the correlation routine
    fac = sqrt(pi / 4.0_dp * alpha * rs) / (2.0_dp * gg)
    tt = absgr * fac
    uu = gr_grabsgr * fac**3
    vv = laplace * fac**2
    ww = 0.0_dp
    ! Factors for the exchange routine
    fac = alpha * rs / 2.0_dp
    ss = absgr * fac
    u2 = gr_grabsgr * fac**3
    v2 = laplace * fac**2

    do ii = 1, nn
      if (rho(ii) < epsilon(1.0_dp)) then
        xcpot(ii) = 0.0_dp
      else
        call correlation_pbe(rs(ii), 0.0_dp, tt(ii), uu(ii), vv(ii), ww, 1, &
            &ec, vcup, vcdn)
        call exchange_pbe(rho(ii), ss(ii), u2(ii), v2(ii), 1, ex, vx)
        if (ieee_is_nan(vcup)) then
          print *, "VCUP NAN", ii, rs(ii), tt(ii), uu(ii), vv(ii)
          print *, ":", absgr(ii), gr_grabsgr(ii), laplace(ii)
          stop
        elseif (ieee_is_nan(vx)) then
          print *, "VX NAN", ii
          stop
        end if
        xcpot(ii) = vcup + vx
      end if
    end do

    deallocate(rho, absgr, laplace, gr_grabsgr)
    deallocate(rs, fac, tt, uu, vv)
     
  end subroutine getxcpot_ggapbe


  
  SUBROUTINE CORRELATION_PBE(RS,ZET,T,UU,VV,WW,igga,ec,vc1,vc2)

    !
    ! APART FROM COSMETICS THIS IS IN FACT BURKEs FORTRAN REFERENCE IMPLEMENTATION
    !

    ! This is the PBE and PW-LDA Correlation routine.

    IMPLICIT REAL(8) (A-H,O-Z)
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

    RETURN
  END subroutine CORRELATION_PBE


  subroutine exchange_pbe(rho,s,u,t,igga,EX,VX)

    ! APART FROM COSMETICS THIS IS IN FACT BURKEs FORTRAN REFERENCE IMPLEMENTATION

    ! This is the PBE and PW-LDA Exchange routine.

    implicit integer(4) (i-n)
    implicit real(8) (a-h,o-z)

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


    RETURN
  END subroutine exchange_pbe

end module dftxc
