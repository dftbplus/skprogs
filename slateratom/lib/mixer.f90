!> Provides a general mixer which contains the desired actual mixer.
module mixer

  use common_accuracy, only : dp
  use broydenmixer, only : TBroydenMixer, TBroydenMixer_mix, TBroydenMixer_reset
  use simplemixer, only : TSimpleMixer, TSimpleMixer_mix, TSimpleMixer_reset
  use diismixer, only : TDiisMixer, TDiisMixer_init, TDiisMixer_mix, TDiisMixer_reset
  implicit none

  private
  public :: TMixer, TMixer_init, TMixer_getMixerType, TMixer_reset, TMixer_mix, mixerTypes


  !> Interface type for mixers
  type TMixer
    private

    !> Numerical type of mixer 1:3
    integer :: mixerType

    !> Simple mixer instance
    type(TSimpleMixer), allocatable :: pSimpleMixer

    !> Broyden mixer instance
    type(TBroydenMixer), allocatable :: pBroydenMixer

    !> DIIS mixer instance
    type(TDiisMixer), allocatable :: pDiisMixer

  end type TMixer


  !> Initialises a specific mixer
  interface TMixer_init
    module procedure TMixer_initSimple
    module procedure TMixer_initBroyden
    module procedure TMixer_initDiis
  end interface TMixer_init


  !> Mixes the given quantity
  interface TMixer_mix
    module procedure TMixer_mix1D
    module procedure TMixer_mix4D
  end interface TMixer_mix


  type :: TMixerTypesEnum
    integer :: simple = 1
    integer :: broyden = 2
    integer :: diis = 3
  end type TMixerTypesEnum

  !> Contains mixer types
  type(TMixerTypesEnum), parameter :: mixerTypes = TMixerTypesEnum()


contains

  !> Initializes a simple mixer.
  subroutine TMixer_initSimple(this, pSimple)

    !> Mixer instance
    type(TMixer), intent(out) :: this

    !> A valid Simple mixer instance on exit
    type(TSimpleMixer), allocatable, intent(inout) :: pSimple

    this%mixerType = mixerTypes%simple
    call move_alloc(pSimple, this%pSimpleMixer)

  end subroutine TMixer_initSimple


  !> Initializes a Broyden mixer.
  subroutine TMixer_initBroyden(this, pBroyden)

    !> Mixer instance
    type(TMixer), intent(out) :: this

    !> A valid Broyden mixer instance on exit
    type(TBroydenMixer), allocatable, intent(inout) :: pBroyden

    this%mixerType = mixerTypes%broyden
    call move_alloc(pBroyden, this%pBroydenMixer)

  end subroutine TMixer_initBroyden


  !> Initializes a Broyden mixer.
  subroutine TMixer_initDiis(this, pDiis)

    !> Mixer instance
    type(TMixer), intent(out) :: this

    !> A valid Broyden mixer instance on exit
    type(TDiisMixer), allocatable, intent(inout) :: pDiis

    this%mixerType = mixerTypes%diis
    call move_alloc(pDiis, this%pDiisMixer)

  end subroutine TMixer_initDiis


  !> Returns mixer type as an int (1: simple, 2: Broyden)
  pure function TMixer_getMixerType(this) result(res)
    !> Mixer instance
    type(TMixer), intent(in) :: this

    integer :: res

    res = this%mixerType

  end function


  !> Resets the mixer.
  subroutine TMixer_reset(this, nElem)

    !> Mixer instance
    type(TMixer), intent(inout) :: this

    !> Size of the vectors to mix
    integer, intent(in) :: nElem

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixer_reset(this%pSimpleMixer, nElem)
    case(mixerTypes%broyden)
      call TBroydenMixer_reset(this%pBroydenMixer, nElem)
    case(mixerTypes%diis)
      call TDiisMixer_reset(this%pDiisMixer, nElem)
    end select

  end subroutine TMixer_reset


  !> Mixes two vectors.
  subroutine TMixer_mix1D(this, inp, diff, errvec)

    !> Mixer instance
    type(TMixer), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    real(dp), intent(inout) :: inp(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    real(dp), intent(in) :: diff(:)

    !> Error vector: measure of lack of convergence
    real(dp), intent(in) :: errvec(:)

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixer_mix(this%pSimpleMixer, inp, diff)
    case(mixerTypes%broyden)
      call TBroydenMixer_mix(this%pBroydenMixer, inp, diff)
    case(mixerTypes%diis)
      call TDiisMixer_mix(this%pDiisMixer, inp, diff, errvec)
    end select

  end subroutine TMixer_mix1D


  !> Mixes two 4D matrices.
  subroutine TMixer_mix4D(this, inp, diff, errvec)

    !> Mixer instance
    type(TMixer), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    real(dp), intent(inout), contiguous, target :: inp(:,0:,:,:)

    !> Difference between input and output vectors
    real(dp), intent(in), contiguous, target :: diff(:,0:,:,:)

    !> Error vector: measure of lack of convergence
    real(dp), intent(in), contiguous, target :: errvec(:,0:,:,:)

    !! Difference between input and output vectors (1D pointer)
    real(dp), pointer :: pDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    real(dp), pointer :: pInp(:)

    !> Error vector: measure of lack of convergence (1D pointers)
    real(dp), pointer :: pErrVec(:)

    pInp(1:size(inp)) => inp
    pDiff(1:size(diff)) => diff
    pErrVec(1:size(diff)) => errvec

    call TMixer_mix1D(this, pInp, pDiff, pErrVec)

  end subroutine TMixer_mix4D

end module mixer
