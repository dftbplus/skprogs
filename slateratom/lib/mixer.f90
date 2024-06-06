!> Provides a general mixer which contains the desired actual mixer.
module mixer

  use common_accuracy, only : dp
  use broydenmixer, only : TBroydenMixer, TBroydenMixer_mix, TBroydenMixer_reset
  use simplemixer, only : TSimpleMixer, TSimpleMixer_mix, TSimpleMixer_reset
  implicit none

  private
  public :: TMixer, TMixer_init, TMixer_reset, TMixer_mix, mixerTypes


  !> Interface type for mixers
  type TMixer
    private

    !> Numerical type of mixer 1:2
    integer :: mixerType

    !> Simple mixer instance
    type(TSimpleMixer), allocatable :: pSimpleMixer

    !> Broyden mixer instance
    type(TBroydenMixer), allocatable :: pBroydenMixer

  end type TMixer


  !> Initialises a specific mixer
  interface TMixer_init
    module procedure TMixer_initSimple
    module procedure TMixer_initBroyden
  end interface TMixer_init


  !> Mixes the given quantity
  interface TMixer_mix
    module procedure TMixer_mix1D
    module procedure TMixer_mix4D
  end interface TMixer_mix


  type :: TMixerTypesEnum
    integer :: simple = 1
    integer :: broyden = 2
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
    end select

  end subroutine TMixer_reset


  !> Mixes two vectors.
  subroutine TMixer_mix1D(this, inp, diff)

    !> Mixer instance
    type(TMixer), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    real(dp), intent(inout) :: inp(:)

    !> Difference between input and output vectors (measure of lack of convergence)
    real(dp), intent(in) :: diff(:)

    select case (this%mixerType)
    case(mixerTypes%simple)
      call TSimpleMixer_mix(this%pSimpleMixer, inp, diff)
    case(mixerTypes%broyden)
      call TBroydenMixer_mix(this%pBroydenMixer, inp, diff)
    end select

  end subroutine TMixer_mix1D


  !> Mixes two 4D matrices.
  subroutine TMixer_mix4D(this, inp, diff)

    !> Mixer instance
    type(TMixer), intent(inout) :: this

    !> Input vector on entry, result vector on exit
    real(dp), intent(inout), contiguous, target :: inp(:,0:,:,:)

    !> Difference between input and output vectors (measure of lack of convergence)
    real(dp), intent(in), contiguous, target :: diff(:,0:,:,:)

    !! Difference between input and output vectors (1D pointer)
    real(dp), pointer :: pDiff(:)

    !! Input vector on entry, result vector on exit (1D pointer)
    real(dp), pointer :: pInp(:)

    pInp(1:size(inp)) => inp
    pDiff(1:size(diff)) => diff

    call TMixer_mix1D(this, pInp, pDiff)

  end subroutine TMixer_mix4D

end module mixer
