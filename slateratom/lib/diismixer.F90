#:include 'common.fypp'

!> Contains the DIIS mixer
!! The DIIS mixing is done by building a weighted combination over the previous input charges to
!! minimise the residue of the error.
!! Only a specified number of previous charge vectors are considered.
!! The modification based on from Kovalenko et al. (J. Comput. Chem., 20: 928-936 1999) and Patrick
!! Briddon to add a contribution from the gradient vector as well is also used.
!! In order to use the mixer you have to create and reset it.
!!
!! Code is adapted from DFTB+.
module diismixer
  use common_accuracy, only : dp
  use lapackroutines, only : gesv
  implicit none

  private

  public :: TDiisMixer, TDiisMixer_init, TDiisMixer_mix, TDiisMixer_reset

  !> Contains the necessary data for an DIIS mixer.
  type TDiisMixer
    private

    !> Initial mixing parameter
    real(dp) :: initMixParam

    !> Max. nr. of stored prev. vectors
    integer :: mPrevVector

    !> Nr. of stored previous vectors
    integer :: iPrevVector

    !> Nr. of elements in the vectors
    integer :: nElem

    !> Index for the storage
    integer :: indx

    !> Stored previous input quantities
    real(dp), allocatable :: prevInput(:,:)

    !> Stored differences of previous input quantities
    real(dp), allocatable :: prevIncrement(:,:)

    !> Stored prev. error vectors
    real(dp), allocatable :: prevErrorVec(:,:)

    !> True if DIIS used from iteration 2 as well as mixing
    logical :: tFromStart

    !> Alpha factor to add in new information
    real(dp) :: alpha

    contains
      procedure :: reset => TDiisMixer_reset
      procedure :: mix1D => TDiisMixer_mix
  end type TDiisMixer

contains

  !> Initializes a DIIS mixer instance.
  subroutine TDiisMixer_init(this, iGenerations, initMixParam, tFromStart)

    !> Pointer to an initialized DIIS mixer on exit
    type(TDiisMixer), intent(out) :: this

    !> Number of generations to consider (including current)
    integer, intent(in) :: iGenerations

    !> Damping parameter for the first mixing steps
    real(dp), intent(in) :: initMixParam

    !> Use DIIS from step 2 onwards?
    logical, intent(in) :: tFromStart

    @:ASSERT(iGenerations >= 2)

    this%nElem = 0
    this%mPrevVector = iGenerations

    allocate(this%prevInput(this%nElem, this%mPrevVector))
    allocate(this%prevErrorVec(this%nElem, this%mPrevVector))
    allocate(this%prevIncrement(this%nElem, this%mPrevVector))

    this%initMixParam = initMixParam
    this%tFromStart = tFromStart

  end subroutine TDiisMixer_init


  !> Makes the mixer ready for a new SCC cycle.
  subroutine TDiisMixer_reset(this, nElem)

    !> DIIS mixer instance
    class(TDiisMixer), intent(inout) :: this

    !> Nr. of elements in the vectors to mix
    integer, intent(in) :: nElem

    @:ASSERT(nElem > 0)

    if (nElem /= this%nElem) then
      this%nElem = nElem
      deallocate(this%prevInput)
      deallocate(this%prevErrorVec)
      deallocate(this%prevIncrement)
      allocate(this%prevInput(this%nElem, this%mPrevVector))
      allocate(this%prevErrorVec(this%nElem, this%mPrevVector))
      allocate(this%prevIncrement(this%nElem, this%mPrevVector))
    end if
    this%iPrevVector = 0
    this%indx = 0

  end subroutine TDiisMixer_reset


  !> Mixes quantities according to the DIIS method.
  subroutine TDiisMixer_mix(this, inputResult, inputIncrement, errorVector)

    !> Pointer to the diis mixer
    class(TDiisMixer), intent(inout) :: this

    !> Input quantity on entry, mixed quantity on exit.
    real(dp), intent(inout) :: inputResult(:)

    !> Increment vector between input and output quantities
    real(dp), intent(in) :: inputIncrement(:)

    !> Error metric vector
    real(dp), intent(in) :: errorVector(:)


    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer :: ii, jj

    @:ASSERT(size(inputResult) == this%nElem)
    @:ASSERT(size(errorVector) == this%nElem)

    if (this%iPrevVector < this%mPrevVector) then
      this%iPrevVector = this%iPrevVector + 1
    end if

    call storeVectors(this%prevInput, this%prevIncrement, this%prevErrorVec, this%indx,&
        & inputResult, inputIncrement, errorVector, this%mPrevVector)
    if (this%tFromStart .or. this%iPrevVector == this%mPrevVector) then

      allocate(aa(this%iPrevVector + 1, this%iPrevVector + 1))
      allocate(bb(this%iPrevVector + 1, 1))

      aa(:,:) = 0.0_dp
      bb(:,:) = 0.0_dp

      ! (due to the hermitian property of our density matrices, the dot-product below is real)
      do ii = 1, this%iPrevVector
        do jj = 1, this%iPrevVector
          aa(ii, jj) = dot_product(this%prevErrorVec(:, ii), this%prevErrorVec(:, jj))
        end do
      end do
      aa(this%iPrevVector + 1, 1:this%iPrevVector) = -1.0_dp
      aa(1:this%iPrevVector, this%iPrevVector + 1) = -1.0_dp

      bb(this%iPrevVector + 1, 1) = -1.0_dp

      ! Solve DIIS system of linear equations
      call gesv(aa, bb)

      inputResult(:) = 0.0_dp
      do ii = 1, this%iPrevVector
        inputResult(:) = inputResult + bb(ii, 1) * (this%prevInput(:, ii) + this%prevIncrement(:, ii))
      end do

    end if

    if (this%iPrevVector < this%mPrevVector) then
      ! First few iterations return simple mixed vector
      inputResult(:) = inputResult + this%initMixParam * inputIncrement(:)
    end if

  end subroutine TDiisMixer_mix


  !> Stores a vector pair in a limited storage.
  !! If the stack is full, oldest vector pair is overwritten.
  subroutine storeVectors(prevInp, prevIncrement, prevErrorVector, indx, input, increment,&
        & errorVector, mPrevVector)

    !> Contains previous vectors of the first type
    real(dp), intent(inout) :: prevInp(:,:)

    !> Contains previous vectors of the second type
    real(dp), intent(inout) :: prevErrorVector(:,:)

    !> Contains previous differences of vectors of first type
    real(dp), intent(inout) :: prevIncrement(:,:)

    !> Indexing of data
    integer, intent(inout) :: indx

    !> New first vector
    real(dp), intent(in) :: input(:)

    !> New increment of first vector
    real(dp), intent(in) :: increment(:)

    !> New second vector
    real(dp), intent(in) :: errorVector(:)

    !> Size of the stacks.
    integer, intent(in) :: mPrevVector

    indx = mod(indx, mPrevVector) + 1
    prevInp(:, indx) = input
    prevIncrement(:, indx) = increment
    prevErrorVector(:, indx) = errorVector
    
  end subroutine storeVectors

end module diismixer