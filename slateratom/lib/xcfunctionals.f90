!> Module related to supported xc-functionals of the slateratom code.
module xcfunctionals

  use, intrinsic :: iso_c_binding, only : c_size_t
  use common_accuracy, only : dp
  use common_constants, only : rec4pi
  use utilities, only : zeroOutCpotOfEmptyDensitySpinChannels
  use xc_f03_lib_m, only : xc_f03_func_t, xc_f03_func_init, xc_f03_func_end, xc_f03_lda_exc_vxc,&
      & xc_f03_gga_exc_vxc, xc_f03_func_set_ext_params, XC_LDA_X, XC_LDA_X_YUKAWA, XC_LDA_C_PW,&
      & XC_GGA_X_PBE, XC_GGA_X_B88, XC_GGA_X_SFAT_PBE, XC_HYB_GGA_XC_B3LYP,&
      & XC_HYB_GGA_XC_CAMY_B3LYP, XC_GGA_C_PBE, XC_GGA_C_LYP, XC_POLARIZED

  implicit none
  private

  public :: xcFunctional, radial_divergence
  public :: getExcVxc_LDA_PW91
  public :: getExcVxc_GGA_PBE96, getExcVxc_GGA_BLYP
  public :: getExcVxc_LCY_PBE96, getExcVxc_LCY_BNL
  public :: getExcVxc_HYB_B3LYP, getExcVxc_HYB_PBE0
  public :: getExcVxc_CAMY_B3LYP, getExcVxc_CAMY_PBEh


  interface libxcVxcToInternalVxc
    module procedure :: libxcVxcToInternalVxc_joined
    module procedure :: libxcVxcToInternalVxc_separate
  end interface libxcVxcToInternalVxc


  !> Enumerates available xc-functionals.
  type :: TXcFunctionalsEnum

    !> Hartree-Fock exchange (for spherically symmetric problems)
    !! Strictly speaking not an xc-functional but nevertheless listed here.
    integer :: HF_Exchange = 0

    !> X-Alpha
    integer :: X_Alpha = 1

    !> LDA-PW91
    integer :: LDA_PW91 = 2

    !> GGA-PBE96
    integer :: GGA_PBE96 = 3

    !> GGA-BLYP
    integer :: GGA_BLYP = 4

    !> LCY-PBE96
    integer :: LCY_PBE96 = 5

    !> LCY-BNL
    integer :: LCY_BNL = 6

    !> HYB-PBE0
    integer :: HYB_PBE0 = 7

    !> HYB-B3LYP
    integer :: HYB_B3LYP = 8

    !> CAMY-B3LYP
    integer :: CAMY_B3LYP = 9

    !> CAMY-PBEh
    integer :: CAMY_PBEh = 10

  contains

    procedure :: isLDA => TXcFunctionalsEnum_isLDA
    procedure :: isGGA => TXcFunctionalsEnum_isGGA
    procedure :: isGlobalHybrid => TXcFunctionalsEnum_isGlobalHybrid
    procedure :: isLongRangeCorrected => TXcFunctionalsEnum_isLongRangeCorrected
    procedure :: isCAMY => TXcFunctionalsEnum_isCAMY
    procedure :: isNotImplemented => TXcFunctionalsEnum_isNotImplemented

  end type TXcFunctionalsEnum


  !> Container for enumerated xc-functional types.
  type(TXcFunctionalsEnum), parameter :: xcFunctional = TXcFunctionalsEnum()


contains

  pure function TXcFunctionalsEnum_isLDA(this, xcnr) result(isLDA)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to an LDA functional
    logical :: isLDA

    isLDA = .false.

    if (xcnr == this%LDA_PW91) isLDA = .true.

  end function TXcFunctionalsEnum_isLDA


  pure function TXcFunctionalsEnum_isGGA(this, xcnr) result(isGGA)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a GGA functional
    logical :: isGGA

    isGGA = .false.

    if (xcnr == this%GGA_PBE96 .or. xcnr == this%GGA_BLYP) isGGA = .true.

  end function TXcFunctionalsEnum_isGGA


  pure function TXcFunctionalsEnum_isLongRangeCorrected(this, xcnr) result(isLongRangeCorrected)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a long-range corrected functional
    logical :: isLongRangeCorrected

    isLongRangeCorrected = .false.

    if (xcnr == this%LCY_PBE96 .or. xcnr == this%LCY_BNL) then
      isLongRangeCorrected = .true.
    end if

  end function TXcFunctionalsEnum_isLongRangeCorrected


  pure function TXcFunctionalsEnum_isGlobalHybrid(this, xcnr) result(isGlobalHybrid)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a global hybrid functional
    logical :: isGlobalHybrid

    isGlobalHybrid = .false.

    if (xcnr == this%HYB_PBE0 .or. xcnr == this%HYB_B3LYP) then
      isGlobalHybrid = .true.
    end if

  end function TXcFunctionalsEnum_isGlobalHybrid


  pure function TXcFunctionalsEnum_isCAMY(this, xcnr) result(isCamy)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a general CAM functional
    logical :: isCamy

    isCamy = .false.

    if (xcnr == this%CAMY_B3LYP .or. xcnr == this%CAMY_PBEh) then
      isCamy = .true.
    end if

  end function TXcFunctionalsEnum_isCAMY


  pure function TXcFunctionalsEnum_isNotImplemented(this, xcnr) result(isNotImplemented)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index is not in the expected range
    logical :: isNotImplemented

    isNotImplemented = .false.

    if (xcnr < 0 .or. xcnr > 10) then
      isNotImplemented = .true.
    end if

  end function TXcFunctionalsEnum_isNotImplemented


  !> Calculates exc and vxc for the LDA-PW91 xc-functional.
  subroutine getExcVxc_LDA_PW91(rho, exc, vxc)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    rhor = transpose(rho) * rec4pi

    allocate(ex(nn))
    allocate(ec(nn))
    allocate(vx(2, nn))
    allocate(vc(2, nn))

    call xc_f03_func_init(xcfunc_x, XC_LDA_X, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_LDA_C_PW, XC_POLARIZED)

    ! exchange
    call xc_f03_lda_exc_vxc(xcfunc_x, nn, rhor(1, 1), ex(1), vx(1, 1))

    ! correlation
    call xc_f03_lda_exc_vxc(xcfunc_c, nn, rhor(1, 1), ec(1), vc(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LDA_PW91


  !> Calculates exc and vxc for the GGA-PBE96 xc-functional.
  subroutine getExcVxc_GGA_PBE96(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_GGA_PBE96


  !> Calculates exc and vxc for the GGA-BLYP xc-functional.
  subroutine getExcVxc_GGA_BLYP(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_B88, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_LYP, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_GGA_BLYP


  !> Calculates exc and vxc for the LCY-PBE96 xc-functional.
  subroutine getExcVxc_LCY_PBE96(abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_SFAT_PBE, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_x, [omega])
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LCY_PBE96


  !> Calculates exc and vxc for the LCY-BNL xc-functional.
  subroutine getExcVxc_LCY_BNL(abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_LDA_X_YUKAWA, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_x, [omega])
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_lda_exc_vxc(xcfunc_x, nn, rhor(1, 1), ex(1), vx(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    exc(:) = ex + ec
    vxc(:,:) = transpose(vx + vc)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_LCY_BNL


  !> Calculates exc and vxc for the HYB-PBE0 xc-functional.
  subroutine getExcVxc_HYB_PBE0(abcissa, dz, dzdr, rho, drho, sigma, camAlpha, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp
    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp
    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1), vxsigma(1, 1))

    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    ! build PBE0 functional
    vxcsigma(:,:) = (1.0_dp - camAlpha) * vxsigma + vcsigma
    vxc(:,:) = transpose((1.0_dp - camAlpha) * vx + vc)
    exc(:) = (1.0_dp - camAlpha) * ex + ec

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_HYB_PBE0


  !> Calculates exc and vxc for the HYB-B3LYP xc-functional.
  subroutine getExcVxc_HYB_B3LYP(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_xc

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exc energy density on grid
    real(dp), allocatable :: exc_tmp(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vxc_tmp(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(exc_tmp(nn))
    exc_tmp(:) = 0.0_dp

    allocate(vxc_tmp(2, nn))
    vxc_tmp(:,:) = 0.0_dp

    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_B3LYP, XC_POLARIZED)
    ! Standard parametrization of B3LYP taken from
    ! J. Phys. Chem. 1994, 98, 45, 11623-11627; DOI: 10.1021/j100096a001
    call xc_f03_func_set_ext_params(xcfunc_xc, [0.20_dp, 0.72_dp, 0.81_dp])

    ! exchange + correlation
    call xc_f03_gga_exc_vxc(xcfunc_xc, nn, rhor(1, 1), sigma(1, 1), exc_tmp(1), vxc_tmp(1, 1),&
        & vxcsigma(1, 1))
    exc(:) = exc_tmp
    vxc(:,:) = transpose(vxc_tmp)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_xc)

  end subroutine getExcVxc_HYB_B3LYP


  !> Calculates exc and vxc for the CAMY-B3LYP xc-functional.
  subroutine getExcVxc_CAMY_B3LYP(abcissa, dz, dzdr, rho, drho, sigma, omega, camAlpha, camBeta,&
      & exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! exc energy density on grid
    real(dp), allocatable :: exc_tmp(:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_xc

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vxc_tmp(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(exc_tmp(nn))
    exc_tmp(:) = 0.0_dp

    allocate(vxc_tmp(2, nn))
    vxc_tmp(:,:) = 0.0_dp

    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_CAMY_B3LYP, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_xc, [0.81_dp, camAlpha + camBeta, -camBeta, omega])

    ! exchange + correlation
    call xc_f03_gga_exc_vxc(xcfunc_xc, nn, rhor(1, 1), sigma(1, 1), exc_tmp(1), vxc_tmp(1, 1),&
        & vxcsigma(1, 1))
    exc(:) = exc_tmp
    vxc(:,:) = transpose(vxc_tmp)

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_xc)

  end subroutine getExcVxc_CAMY_B3LYP


  !> Calculates exc and vxc for the CAMY-PBEh xc-functional.
  subroutine getExcVxc_CAMY_PBEh(abcissa, dz, dzdr, rho, drho, sigma, omega, camAlpha, camBeta,&
      & exc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> contracted gradients of the density
    real(dp), intent(in), allocatable :: sigma(:,:)

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_x, xcfunc_xsr, xcfunc_c

    !! number of density grid points
    integer(c_size_t) :: nn

    !! exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ex_sr(:), ec(:)

    !! exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vx_sr(:,:), vc(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), allocatable :: vxsigma(:,:), vxsigma_sr(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), allocatable :: vcsigma(:,:)

    !! first partial derivative of the energy per unit volume in terms of sigma (x+c)
    real(dp), allocatable :: vxcsigma(:,:)

    nn = size(rho, dim=1)
    ! divide by 4*pi to catch different normalization of spherical harmonics
    allocate(rhor(2, nn))
    rhor(:,:) = transpose(rho) * rec4pi

    allocate(ex(nn))
    ex(:) = 0.0_dp
    allocate(ex_sr(nn))
    ex_sr(:) = 0.0_dp
    allocate(ec(nn))
    ec(:) = 0.0_dp

    allocate(vx(2, nn))
    vx(:,:) = 0.0_dp
    allocate(vx_sr(2, nn))
    vx_sr(:,:) = 0.0_dp
    allocate(vc(2, nn))
    vc(:,:) = 0.0_dp

    allocate(vxsigma(3, nn))
    vxsigma(:,:) = 0.0_dp
    allocate(vxsigma_sr(3, nn))
    vxsigma_sr(:,:) = 0.0_dp
    allocate(vcsigma(3, nn))
    vcsigma(:,:) = 0.0_dp
    allocate(vxcsigma(3, nn))
    vxcsigma(:,:) = 0.0_dp

    ! short-range exchange
    call xc_f03_func_init(xcfunc_xsr, XC_GGA_X_SFAT_PBE, XC_POLARIZED)
    call xc_f03_func_set_ext_params(xcfunc_xsr, [omega])

    ! full-range exchange
    call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)

    ! correlation
    call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)

    ! short-range exchange
    call xc_f03_gga_exc_vxc(xcfunc_xsr, nn, rhor(1, 1), sigma(1, 1), ex_sr(1), vx_sr(1, 1),&
        & vxsigma_sr(1, 1))
    ! full-range exchange
    call xc_f03_gga_exc_vxc(xcfunc_x, nn, rhor(1, 1), sigma(1, 1), ex(1), vx(1, 1),&
        & vxsigma(1, 1))
    ! correlation
    call xc_f03_gga_exc_vxc(xcfunc_c, nn, rhor(1, 1), sigma(1, 1), ec(1), vc(1, 1), vcsigma(1, 1))
    call zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    ! build CAMY-PBEh functional
    vxcsigma(:,:) = camBeta * vxsigma_sr + (1.0_dp - (camAlpha + camBeta)) * vxsigma + vcsigma
    vxc(:,:) = transpose(camBeta * vx_sr + (1.0_dp - (camAlpha + camBeta)) * vx + vc)
    exc(:) = camBeta * ex_sr + (1.0_dp - (camAlpha + camBeta)) * ex + ec

    call libxcVxcToInternalVxc(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    ! finalize libxc objects
    call xc_f03_func_end(xcfunc_x)
    call xc_f03_func_end(xcfunc_xsr)
    call xc_f03_func_end(xcfunc_c)

  end subroutine getExcVxc_CAMY_PBEh


  !> Converts libXC vxc to our internal representation.
  subroutine libxcVxcToInternalVxc_joined(abcissa, dz, dzdr, drho, vxcsigma, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> first partial derivative of the energy per unit volume in terms of sigma (xc)
    real(dp), intent(in) :: vxcsigma(:,:)

    !> xc potential on grid
    real(dp), intent(inout) :: vxc(:,:)

    !! temporary storage
    real(dp), allocatable :: tmpv1(:), tmpv2(:)

    !! auxiliary variables
    integer :: iSpin, iSpin2, iSigma

    allocate(tmpv1(size(drho, dim=1)))
    allocate(tmpv2(size(drho, dim=1)))

    ! derivative of E vs. grad n
    do iSpin = 1, 2
      ! the other spin
      iSpin2 = 3 - iSpin
      ! 1 for spin up, 3 for spin down
      iSigma = 2 * iSpin - 1
      tmpv1(:) = vxcsigma(iSigma, :) * drho(:, iSpin) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - 2.0_dp * tmpv2
      tmpv1(:) = vxcsigma(2, :) * drho(:, iSpin2) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - tmpv2
    end do

  end subroutine libxcVxcToInternalVxc_joined


  !> Converts libXC vxc to our internal representation.
  subroutine libxcVxcToInternalVxc_separate(abcissa, dz, dzdr, drho, vxsigma, vcsigma, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> 1st deriv. of density on grid
    real(dp), intent(in) :: drho(:,:)

    !> first partial derivative of the energy per unit volume in terms of sigma (exchange)
    real(dp), intent(in) :: vxsigma(:,:)

    !> first partial derivative of the energy per unit volume in terms of sigma (correlation)
    real(dp), intent(in) :: vcsigma(:,:)

    !> xc potential on grid
    real(dp), intent(inout) :: vxc(:,:)

    !! temporary storage
    real(dp), allocatable :: tmpv1(:), tmpv2(:)

    !! auxiliary variables
    integer :: iSpin, iSpin2, iSigma

    allocate(tmpv1(size(drho, dim=1)))
    allocate(tmpv2(size(drho, dim=1)))

    ! derivative of E vs. grad n
    do iSpin = 1, 2
      ! the other spin
      iSpin2 = 3 - iSpin
      ! 1 for spin up, 3 for spin down
      iSigma = 2 * iSpin - 1
      tmpv1(:) = (vxsigma(iSigma, :) + vcsigma(iSigma, :)) * drho(:, iSpin) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - 2.0_dp * tmpv2
      tmpv1(:) = (vxsigma(2, :) +  vcsigma(2, :)) * drho(:, iSpin2) * rec4pi
      call radial_divergence(tmpv1, abcissa, dz, tmpv2, dzdr)
      vxc(:, iSpin) = vxc(:, iSpin) - tmpv2
    end do

  end subroutine libxcVxcToInternalVxc_separate


  !>
  pure subroutine radial_divergence(ff, rr, dr, rdiv, jacobi)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: rr(:)
    real(dp), intent(in) :: dr
    real(dp), intent(out) :: rdiv(:)
    real(dp), intent(in), optional :: jacobi(:)

    call derive1_5(ff, dr, rdiv, jacobi)
    rdiv(:) = rdiv + 2.0_dp / rr * ff

  end subroutine radial_divergence


  !>
  pure subroutine derive(ff, dx, jacobi)
    real(dp), intent(inout) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(in), optional :: jacobi(:)

    real(dp), allocatable :: tmp1(:)
    integer :: nn

    nn = size(ff)
    allocate(tmp1(nn))
    tmp1(:) = ff
    ff(2:nn - 1) = (ff(3:nn) - ff(1:nn - 2)) / (2.0 * dx)
    ff(1) = (tmp1(2) - tmp1(1)) / dx
    ff(nn) = (tmp1(nn) - tmp1(nn - 1)) / dx
    if (present(jacobi)) then
      ff = ff * jacobi
    end if

  end subroutine derive


  !>
  pure subroutine derive1_5(ff, dx, dfdx, dudx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: dfdx(:)
    real(dp), intent(in), optional :: dudx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) =&
        reshape([&
        & -25.0_dp,  48.0_dp, -36.0_dp,  16.0_dp, -3.0_dp,&
        &  -3.0_dp, -10.0_dp,  18.0_dp,  -6.0_dp,  1.0_dp,&
        &   1.0_dp,  -8.0_dp,   0.0_dp,   8.0_dp, -1.0_dp,&
        &  -1.0_dp,   6.0_dp, -18.0_dp,  10.0_dp,  3.0_dp,&
        &   3.0_dp,  -16.0_dp, 36.0_dp, -48.0_dp, 25.0_dp], [np, np])

    integer :: ngrid
    integer :: ii

    ngrid = size(ff)
    do ii = 1, nleft
      dfdx(ii) = dot_product(coeffs(:, ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      dfdx(ii) = dot_product(coeffs(:, imiddle), ff(ii - nleft:ii + nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      dfdx(ii) = dot_product(coeffs(:, np - (ngrid - ii)), ff(ngrid - np + 1:ngrid))
    end do

    if (present(dudx)) then
      dfdx = dfdx * (dudx / (dxprefac * dx))
    else
      dfdx = dfdx / (dxprefac * dx)
    end if

  end subroutine derive1_5


  !>
  pure subroutine derive2_5(ff, dx, d2fdx2, dudx, d2udx2, dfdx)
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
        &  35.0_dp, -104.0_dp, 114.0_dp, -56.0_dp, 11.0_dp, &
        &  11.0_dp, -20.0_dp, 6.0_dp, 4.0_dp, -1.0_dp, &
        &  -1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp, &
        &  -1.0_dp, 4.0_dp, 6.0_dp, -20.0_dp, 11.0_dp, &
        &  11.0_dp, -56.0_dp, 114.0_dp, -104.0_dp, 35.0_dp], [np, np])

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
      d2fdx2(ii) = dot_product(coeffs(:, ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      d2fdx2(ii) = dot_product(coeffs(:, imiddle), ff(ii - nleft:ii + nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      d2fdx2(ii) = dot_product(coeffs(:, np - (ngrid - ii)), ff(ngrid - np + 1:ngrid))
    end do

    if (present(dudx)) then
      d2fdx2 = d2fdx2 * (dudx * dudx / (dxprefac * dx * dx))
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

end module xcfunctionals
