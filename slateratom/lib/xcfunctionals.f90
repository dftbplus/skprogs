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

    !> LCY-PBE
    integer :: LCY_PBE96 = 5

    !> LCY-BNL
    integer :: LCY_BNL = 6

  contains

    procedure :: isLDA => TXcFunctionalsEnum_isLDA
    procedure :: isGGA => TXcFunctionalsEnum_isGGA
    procedure :: isGlobalHybrid => TXcFunctionalsEnum_isGlobalHybrid
    procedure :: isLongRangeCorrected => TXcFunctionalsEnum_isLongRangeCorrected
    procedure :: isCam => TXcFunctionalsEnum_isCam

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

  end function TXcFunctionalsEnum_isGlobalHybrid


  pure function TXcFunctionalsEnum_isCam(this, xcnr) result(isCam)

    !> Class instance
    class(TXcFunctionalsEnum), intent(in) :: this

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> True, if xc-functional index corresponds to a general CAM functional
    logical :: isCam

    isCam = .false.

  end function TXcFunctionalsEnum_isCam

end module xcfunctionals
