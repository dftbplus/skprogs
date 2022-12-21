!> Module related to supported xc-functionals of the slateratom code.
module xcfunctionals

  implicit none
  private

  public :: xcFunctional


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

end module xcfunctionals
