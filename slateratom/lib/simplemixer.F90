#:include 'common.fypp'

!> Simple mixer for mixing charges.
module simplemixer

  use common_accuracy, only : dp
  implicit none

  private
  public :: TSimpleMixer
  public :: TSimpleMixer_init, TSimpleMixer_reset, TSimpleMixer_mix


  !> Contains data for a simple mixer
  type TSimpleMixer
    private

    !> Mixing parameter
    real(dp) :: mixParam

  end type TSimpleMixer


contains

  !> Creates a simple mixer.
  subroutine TSimpleMixer_init(this, mixParam)

    !> Simple mixer instance on exit
    type(TSimpleMixer), intent(out) :: this

    !> Mixing parameter
    real(dp), intent(in) :: mixParam

    this%mixParam = mixParam

  end subroutine TSimpleMixer_init


  !> Resets the mixer.
  subroutine TSimpleMixer_reset(this, nElem)

    !> Simple mixer instance
    type(TSimpleMixer), intent(inout) :: this

    !> Length of the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    continue

  end subroutine TSimpleMixer_reset


  !> Does the actual mixing.
  subroutine TSimpleMixer_mix(this, qInpResult, qDiff)

    !> SimpleMixer instance
    type(TSimpleMixer), intent(inout) :: this

    !> Input charge on entry, mixed charge on exit
    real(dp), intent(inout) :: qInpResult(:)

    !> Charge difference
    real(dp), intent(in) :: qDiff(:)

    @:ASSERT(size(qInpResult) == size(qDiff))

    qInpResult(:) = qInpResult + this%mixParam * qDiff

  end subroutine TSimpleMixer_mix

end module simplemixer
