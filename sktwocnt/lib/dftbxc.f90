!> Module that provides exchange-correlation DFT routines.
module dftxc

  use, intrinsic :: ieee_arithmetic
  use common_accuracy, only : dp
  use common_constants, only : pi

  !! vanderhe: proposed libxc integration
  ! use, intrinsic :: iso_c_binding, only : c_size_t
  ! use xc_f90_lib_m, only : xc_f90_func_t, xc_f90_func_info_t, xc_f90_func_init,&
  !     & xc_f90_func_get_info, xc_f90_lda_vxc, xc_f90_gga_vxc, xc_f90_func_end, XC_LDA_X,&
  !     & XC_LDA_C_PW, XC_GGA_X_PBE, XC_GGA_C_PBE, XC_UNPOLARIZED

  implicit none
  private

  public :: getxcpot_ldapw91, getxcpot_ggapbe

  !> pre-factor for re-normalization
  real(dp), parameter :: rec4pi = 1.0_dp / (4.0_dp * pi)


contains

  !> Calculates xc-potential based on the LDA-PW91 functional.
  subroutine getxcpot_ldapw91(rho4pi, xcpot)

    !> density times 4pi on grid
    real(dp), intent(in) :: rho4pi(:)

    !> resulting xc-potential
    real(dp), intent(out) :: xcpot(:)

    !> density with libxc compatible normalization
    real(dp), allocatable :: rho(:)

    !> local Seitz radius, needed for functional evaluation
    real(dp), allocatable :: rs(:)

    !> exchange and correlation (up, down) potential of a single grid point
    real(dp) :: vx, vcup, vcdn

    !> exchange and correlation energy of a single grid point
    real(dp) :: ex, ec

    !> number of density grid points
    integer :: nn

    !> auxiliary variable
    integer :: ii

    nn = size(rho4pi)
    allocate(rs(nn), rho(nn))

    ! renorm rho (incoming quantity is 4pi normed)
    rho = rho4pi * rec4pi
    ! note: rho is normed to 4pi, therefore 4*pi missing in rs
    rs = (3.0_dp / rho4pi)**(1.0_dp / 3.0_dp)
    do ii = 1, nn
      if (rho(ii) < epsilon(1.0_dp)) then
        xcpot(ii) = 0.0_dp
      else
        call correlation_pbe(rs(ii), 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0, ec, vcup, vcdn)
        call exchange_pbe(rho(ii), 0.0_dp, 0.0_dp, 0.0_dp, 0, ex, vx)
        xcpot(ii) = vcup + vx
      end if
    end do

    !! vanderhe: proposed libxc integration
    !! --> but Hamiltonian matrix elements differ up to 1e-07 a.u. (something is wrong)!?

    ! !> libxc related objects
    ! type(xc_f90_func_t) :: xcfunc_x, xcfunc_c
    ! type(xc_f90_func_info_t) :: xcinfo

    ! !> density with libxc compatible normalization
    ! real(dp), allocatable :: rho(:)

    ! !> exchange and correlation potential on grid
    ! real(dp), allocatable :: vx(:), vc(:)

    ! !> number of density grid points
    ! integer(c_size_t) :: nn

    ! call xc_f90_func_init(xcfunc_x, XC_LDA_X, XC_UNPOLARIZED)
    ! xcinfo = xc_f90_func_get_info(xcfunc_x)
    ! call xc_f90_func_init(xcfunc_c, XC_LDA_C_PW, XC_UNPOLARIZED)
    ! xcinfo = xc_f90_func_get_info(xcfunc_x)

    ! nn = size(rho4pi)
    ! allocate(vx(nn), vc(nn))

    ! rho = rho4pi * rec4pi

    ! call xc_f90_lda_vxc(xcfunc_x, nn, rho, vx)
    ! call xc_f90_lda_vxc(xcfunc_c, nn, rho, vc)

    ! xcpot(:) = vx + vc

    ! call xc_f90_func_end(xcfunc_x)
    ! call xc_f90_func_end(xcfunc_c)

  end subroutine getxcpot_ldapw91


  !> Calculates xc-potential based on the GGA-PBE functional.
  subroutine getxcpot_ggapbe(rho4pi, absgr4pi, laplace4pi, gr_grabsgr4pi, xcpot)

    !> density times 4pi on grid
    real(dp), intent(in) :: rho4pi(:)

    !> absolute gradient of density times 4pi on grid
    real(dp), intent(in) :: absgr4pi(:)

    !> laplace operator acting on density times 4pi on grid
    real(dp), intent(in) :: laplace4pi(:)

    !> (grad rho4pi) * grad(abs(grad rho4pi))
    real(dp), intent(in) :: gr_grabsgr4pi(:)

    !> resulting xc-potential
    real(dp), intent(out) :: xcpot(:)

    !> density with libxc compatible normalization
    real(dp), allocatable :: rho(:)

    !> absolute gradient of density on grid
    real(dp), allocatable :: absgr(:)

    !> laplace operator acting on density on grid
    real(dp), allocatable :: laplace(:)

    !> (grad rho) * grad(abs(grad rho)) / rho**2
    !! actually calculated based on rho4pi, but 4pi cancels out
    real(dp), allocatable :: gr_grabsgr(:)

    !> number of density grid points
    integer :: nn

    !> auxiliary variables
    real(dp), allocatable :: rs(:), fac(:), tt(:), uu(:), vv(:)
    real(dp), allocatable :: ss(:), u2(:), v2(:)
    real(dp) :: alpha, zeta, gg, ww
    real(dp) :: ec, vcup, vcdn, ex, vx
    integer :: ii

    nn = size(rho4pi)
    allocate(rho(nn), absgr(nn), laplace(nn), gr_grabsgr(nn))
    allocate(rs(nn), fac(nn), tt(nn), uu(nn), vv(nn), ss(nn), u2(nn), v2(nn))

    ! renorm rho and derivatives (incoming quantities are 4pi normed)
    rho = rho4pi * rec4pi
    absgr = absgr4pi / rho4pi
    laplace = laplace4pi / rho4pi
    gr_grabsgr = gr_grabsgr4pi / rho4pi**2

    ! note: rho is normed to 4pi, therefore 4*pi missing in rs
    rs = (3.0_dp / rho4pi)**(1.0_dp / 3.0_dp)
    zeta = 0.0_dp
    gg = 1.0_dp
    alpha = (4.0_dp / (9.0_dp * pi))**(1.0_dp / 3.0_dp)

    ! factors for the correlation routine
    fac = sqrt(pi / 4.0_dp * alpha * rs) / (2.0_dp * gg)
    tt = absgr * fac
    uu = gr_grabsgr * fac**3
    vv = laplace * fac**2
    ww = 0.0_dp

    ! factors for the exchange routine
    fac = alpha * rs / 2.0_dp
    ss = absgr * fac
    u2 = gr_grabsgr * fac**3
    v2 = laplace * fac**2

    do ii = 1, nn
      if (rho(ii) < epsilon(1.0_dp)) then
        xcpot(ii) = 0.0_dp
      else
        call correlation_pbe(rs(ii), zeta, tt(ii), uu(ii), vv(ii), ww, 1, ec, vcup, vcdn)
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

    !! vanderhe: proposed libxc integration
    !! --> but Hamiltonian matrix elements differ up to 1e-02 a.u. (something is wrong)!?

    ! !> libxc related objects
    ! type(xc_f90_func_t) :: xcfunc_x, xcfunc_c
    ! type(xc_f90_func_info_t) :: xcinfo

    ! !> density with libxc compatible normalization
    ! real(dp), allocatable :: rho(:)

    ! !> contracted gradients of the density
    ! real(dp), allocatable :: sigma(:)

    ! !> exchange and correlation potential on grid
    ! real(dp), allocatable :: vx(:), vc(:)

    ! !> first partial derivative of the energy per unit volume in terms of sigma
    ! real(dp), allocatable :: vxsigma(:), vcsigma(:)

    ! !> number of density grid points
    ! integer(c_size_t) :: nn

    ! nn = size(rho4pi)
    ! allocate(vx(nn), vc(nn), vxsigma(nn), vcsigma(nn))

    ! rho = rho4pi * rec4pi
    ! sigma = (absgr4pi * rec4pi)**2

    ! call xc_f90_func_init(xcfunc_x, XC_GGA_X_PBE, XC_UNPOLARIZED)
    ! xcinfo = xc_f90_func_get_info(xcfunc_x)
    ! call xc_f90_func_init(xcfunc_c, XC_GGA_C_PBE, XC_UNPOLARIZED)
    ! xcinfo = xc_f90_func_get_info(xcfunc_x)

    ! call xc_f90_gga_vxc(xcfunc_x, nn, rho, sigma, vx, vxsigma)
    ! call xc_f90_gga_vxc(xcfunc_c, nn, rho, sigma, vc, vcsigma)

    ! xcpot(:) = vx + vc

    ! call xc_f90_func_end(xcfunc_x)
    ! call xc_f90_func_end(xcfunc_c)

  end subroutine getxcpot_ggapbe


  SUBROUTINE CORRELATION_PBE(RS, ZET, T, UU, VV, WW, igga, ec, vc1, vc2)

    !
    ! APART FROM COSMETICS THIS IS IN FACT BURKEs FORTRAN REFERENCE IMPLEMENTATION
    !

    ! This is the PBE and PW-LDA Correlation routine.

    IMPLICIT REAL(8) (A - H, O - Z)
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
    parameter(thrd=1._dp / 3._dp, thrdm=-thrd, thrd2=2._dp * thrd)
    parameter(GAM=0.5198420997897463295344212145565_dp)
    parameter(thrd4=4._dp * thrd, fzz=8._dp / (9._dp * GAM))
    parameter(gamma=0.03109069086965489503494086371273_dp)
    parameter(bet=0.06672455060314922_dp, delt=bet / gamma)
    dimension u(6), p(6), s(6)
    data u/0.03109070_dp, 0.2137000_dp, 7.5957000_dp,&
        &        3.58760000_dp, 1.6382000_dp, 0.4929400_dp/
    data p/0.01554535_dp, 0.2054800_dp, 14.1189000_dp,&
        &        6.19770000_dp, 3.3662000_dp, 0.6251700_dp/
    data s/0.01688690_dp, 0.1112500_dp, 10.3570000_dp,&
        &        3.62310000_dp, 0.8802600_dp, 0.4967100_dp/
    !----------------------------------------------------------------------
    !     find LSD energy contributions, using [c](10) .
    !     EU=unpolarized LSD correlation energy , EURS=dEU/drs
    !     EP=fully polarized LSD correlation energy , EPRS=dEP/drs
    !     ALFM=-spin stiffness, [c](3) , ALFRSM=-dalpha/drs .
    !     F=spin-scaling factor from [c](9).
    !     construct ecl, using [c](8) .
    !

    rtrs = dsqrt(rs)
    Q0 = -2._dp * u(1) * (1._dp + u(2) * rtrs * rtrs)
    Q1 = 2._dp * u(1) * rtrs * (u(3) + rtrs * (u(4) + rtrs * (u(5) + u(6) * rtrs)))
    Q2 = DLOG(1._dp + 1._dp / Q1)
    Q3 = u(1) * (u(3) / rtrs + 2._dp * u(4) + rtrs * (3._dp * u(5) + 4._dp * u(6) * rtrs))
    EU = Q0 * Q2
    EURS = -2._dp * u(1) * u(2) * Q2 - Q0 * Q3 / (Q1 * (1._dp + Q1))
    Q0 = -2._dp * p(1) * (1._dp + p(2) * rtrs * rtrs)
    Q1 = 2._dp * p(1) * rtrs * (p(3) + rtrs * (p(4) + rtrs * (p(5) + p(6) * rtrs)))
    Q2 = DLOG(1._dp + 1._dp / Q1)
    Q3 = p(1) * (p(3) / rtrs + 2._dp * p(4) + rtrs * (3._dp * p(5) + 4._dp * p(6) * rtrs))
    EP = Q0 * Q2
    EPRS = -2._dp * p(1) * p(2) * Q2 - Q0 * Q3 / (Q1 * (1._dp + Q1))
    Q0 = -2._dp * s(1) * (1._dp + s(2) * rtrs * rtrs)
    Q1 = 2._dp * s(1) * rtrs * (s(3) + rtrs * (s(4) + rtrs * (s(5) + s(6) * rtrs)))
    Q2 = DLOG(1._dp + 1._dp / Q1)
    Q3 = s(1) * (s(3) / rtrs + 2._dp * s(4) + rtrs * (3._dp * s(5) + 4._dp * s(6) * rtrs))
    ALFM = Q0 * Q2
    ALFRSM = -2._dp * s(1) * s(2) * Q2 - Q0 * Q3 / (Q1 * (1._dp + Q1))

    Z4 = ZET**4
    F = ((1._dp + ZET)**THRD4 + (1._dp - ZET)**THRD4 - 2._dp) / GAM
    ECL = EU * (1._dp - F * Z4) + EP * F * Z4 - ALFM * F * (1._dp - Z4) / FZZ
    !----------------------------------------------------------------------
    !     LSD potential from [c](A1)
    !     ECRS = dEc/drs , ECZET=dEc/dzeta , FZ = dF/dzeta   [c](A2-A4)
    !
    ECRS = EURS * (1._dp - F * Z4) + EPRS * F * Z4 - ALFRSM * F * (1._dp - Z4) / FZZ
    FZ = THRD4 * ((1._dp + ZET)**THRD - (1._dp - ZET)**THRD) / GAM
    ECZET = 4._dp * (ZET**3) * F * (EP - EU + ALFM / FZZ)&
        &        + FZ * (Z4 * EP - Z4 * EU - (1._dp - Z4) * ALFM / FZZ)
    COMM = ECL - RS * ECRS / 3._dp - ZET * ECZET
    VCUP = COMM + ECZET
    VCDN = COMM - ECZET
    if (igga .eq. 0) then
      EC = ECL
      VC1 = VCUP
      VC2 = VCDN
      return
    end if
    !----------------------------------------------------------------------
    !     PBE correlation energy
    !     G=phi(zeta), given after [a](3)
    !     DELT=bet/gamma , B=A of [a](8)
    !
    G = ((1._dp + ZET)**thrd2 + (1._dp - ZET)**thrd2) / 2._dp
    G3 = G**3
    PON = -ECL / (G3 * gamma)
    B = DELT / (DEXP(PON) - 1._dp)
    B2 = B * B
    T2 = T * T
    T4 = T2 * T2
    Q4 = 1._dp + B * T2
    Q5 = 1._dp + B * T2 + B2 * T4
    ECN = G3 * (BET / DELT) * DLOG(1._dp + DELT * Q4 * T2 / Q5)
    EC = ECL + ECN
    !----------------------------------------------------------------------
    !     ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
    !
    G4 = G3 * G
    T6 = T4 * T2
    RSTHRD = RS / 3._dp
    !      GZ=((1._dp+zet)**thirdm-(1._dp-zet)**thirdm)/3._dp
    ! ckoe: hack thirdm never gets defined, but 1-1 should be zero anyway
    GZ = 0.0_dp
    FAC = DELT / B + 1._dp
    BG = -3._dp * B2 * ECL * FAC / (BET * G4)
    BEC = B2 * FAC / (BET * G3)
    Q8 = Q5 * Q5 + DELT * Q4 * Q5 * T2
    Q9 = 1._dp + 2._dp * B * T2
    hB = -BET * G3 * B * T6 * (2._dp + B * T2) / Q8
    hRS = -RSTHRD * hB * BEC * ECRS
    FACT0 = 2._dp * DELT - 6._dp * B
    FACT1 = Q5 * Q9 + Q4 * Q9 * Q9
    hBT = 2._dp * BET * G3 * T4 * ((Q4 * Q5 * FACT0 - DELT * FACT1) / Q8) / Q8
    hRST = RSTHRD * T2 * hBT * BEC * ECRS
    hZ = 3._dp * GZ * ecn / G + hB * (BG * GZ + BEC * ECZET)
    hT = 2._dp * BET * G3 * Q9 / Q8
    hZT = 3._dp * GZ * hT / G + hBT * (BG * GZ + BEC * ECZET)
    FACT2 = Q4 * Q5 + B * T2 * (Q4 * Q9 + Q5)
    FACT3 = 2._dp * B * Q5 * Q9 + DELT * FACT2
    hTT = 4._dp * BET * G3 * T * (2._dp * B / Q8 - (Q9 * FACT3 / Q8) / Q8)
    COMM = ECN + HRS + HRST + T2 * HT / 6._dp + 7._dp * T2 * T * HTT / 6._dp
    PREF = HZ - GZ * T2 * HT / G
    FACT5 = GZ * (2._dp * HT + T * HTT) / G
    COMM = COMM - PREF * ZET - UU * HTT - VV * HT - WW * (HZT - FACT5)
    DVCUP = COMM + PREF
    DVCDN = COMM - PREF
    VC1 = VCUP + DVCUP
    VC2 = VCDN + DVCDN

    RETURN
  END subroutine CORRELATION_PBE


  subroutine exchange_pbe(rho, s, u, t, igga, EX, VX)

    ! APART FROM COSMETICS THIS IS IN FACT BURKEs FORTRAN REFERENCE IMPLEMENTATION

    ! This is the PBE and PW-LDA Exchange routine.

    implicit integer(4) (i - n)
    implicit real(8) (a - h, o - z)

    parameter(thrd=1._dp / 3._dp, thrd4=4._dp / 3._dp)
    parameter(pi=3.14159265358979323846264338327950_dp)
    parameter(ax=-0.738558766382022405884230032680836_dp)

    parameter(um=0.21951_dp, uk=0.8040_dp, ul=um / uk)

    parameter(ap=1.647127_dp, bp=0.980118_dp, cp=0.017399_dp)
    parameter(aq=1.523671_dp, bq=0.367229_dp, cq=0.011282_dp)
    parameter(ah=0.19645_dp, bh=7.7956_dp)
    parameter(ahp=0.27430_dp, bhp=0.15084_dp, ahq=0.004_dp)
    parameter(a1=0.19645_dp, a2=0.27430_dp, a3=0.15084_dp, a4=100._dp)
    parameter(a=7.79560_dp, b1=0.004_dp, eps=1.d-15)

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
    ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submitted to PRL, May96
    ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
    !     {\bf 40},  3399  (1989) (E).
    !----------------------------------------------------------------------
    ! Formulas: e_x[unif]=ax*rho^(4/3)  [LDA]
    !           ax = -0.75*(3/pi)^(1/3)
    !            e_x[PBE]=e_x[unif]*FxPBE(s)
    !            FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
    !           uk, ul defined after [a](13)
    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !     construct LDA exchange energy density

    exunif = ax * rho**thrd
    if ((igga .eq. 0) .or. (s .lt. eps)) then
      EXL = exunif
      EXN = 0._dp
      EX = EXL + EXN
      VX = exunif * thrd4
      return
    end if
    !----------------------------------------------------------------------
    !     construct GGA enhancement factor
    !     find first and second derivatives of f and:
    !     fs=(1/s)*df/ds  and  fss=dfs/ds = (d2f/ds2 - (1/s)*df/ds)/s

    !
    ! PBE enhancement factors checked against NRLMOL
    !
    if (igga .eq. 1) then
      p0 = 1._dp + ul * s**2
      f = 1._dp + uk - uk / p0
      fs = 2._dp * uk * ul / p0**2
      fss = -4._dp * ul * s * fs / p0
    end if

    !

    EXL = exunif
    EXN = exunif * (f - 1.0_dp)
    EX = EXL + EXN
    !----------------------------------------------------------------------
    !     energy done. calculate potential from [b](24)
    !
    VX = exunif * (thrd4 * f - (u - thrd4 * s**3) * fss - t * fs)

    RETURN
  END subroutine exchange_pbe

end module dftxc
