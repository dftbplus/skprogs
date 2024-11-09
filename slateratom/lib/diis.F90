#:include 'common.fypp'

module diis
  use common_accuracy, only : dp
  use lapackroutines, only : getrs

  implicit none
    
  private
  public :: TDIISMixer, TDIISMixer_init, TDIISMixer_reset, TDIISMixer_mix

  !> DIIS mixer class
  !! see 10.1002/jcc.540030413 for details
  type, public :: TDIISMixer
    private

    !> Actual iteration
    integer :: iIter

    !> Max. number of iterations
    integer :: mIter

    !> N. of vectors to store for DIIS
    integer :: mVec

    !> Mixing parameter
    real(dp) :: mixParam

    !> Potential matrix storage vector
    real(dp), allocatable :: pot_storage(:,:,:,:,:)

    !> Commutator storage vector
    real(dp), allocatable :: comm_storage(:,:,:,:,:)

    !> B matrix
    real(dp), allocatable :: B(:,:)

    !> Constraint vector
    real(dp), allocatable :: c_vec(:)

    contains

end type TDIISMixer

contains
  
  !> Creates a DIIS mixer instance
  !! TODO: complete description
  subroutine TDIISMixer_init(this, mVec, mIter, mixParam, max_l, num_alpha, poly_order)
    
    !> An initialized Broyden mixer on exit
    type(TDIISMixer), intent(out) :: this
    
    !> Max. nr. of vectors to store
    integer, intent(in) :: mVec

    !> Max. nr. of iterations
    integer, intent(in) :: mIter

    !> Mixing parameter
    real(dp), intent(in) :: mixParam

    !> Max. angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !! auxilliary variables
    integer :: problemsize

    @:ASSERT(mVec > 0)
    @:ASSERT(mixParam > 0.0_dp)
    
    problemsize = maxval(num_alpha) * maxval(poly_order)

    this%iIter = 0
    this%mVec = mVec
    this%mIter = mIter
    this%mixParam = mixParam
    allocate(this%pot_storage(mVec, 2, 0:max_l, problemsize, problemsize))
    allocate(this%comm_storage(mVec, 2, 0:max_l, problemsize, problemsize))
    allocate(this%B(mVec+1, mVec+1))
    allocate(this%c_vec(mVec+1))

    this%pot_storage = 0.0_dp
    this%comm_storage = 0.0_dp

    this%B = 0.0_dp
    this%c_vec = 0.0_dp
    this%B(1,:) = -1.0_dp
    this%B(:,1) = -1.0_dp
    this%B(1, 1) = 0.0_dp
    this%c_vec(1) = -1.0_dp
  
  end subroutine TDIISMixer_init


  !> Mixes 4D matrices according to DIIS scheme
  !! TODO: complete description
  subroutine TDIISMixer_mix(this, inpResult, diff, commutator)

    !> DIIS mixer object
    type(TDIISMixer), intent(inout) :: this

    !> Input Fock matrix on entry, result Fock matrix on exit
    real(dp), intent(inout) :: inpResult(:,0:,:,:)

    !> Difference between Fock matrices on this and previous steps
    real(dp), intent(in) :: diff(:,0:,:,:)

    !> Commutator [F,PS]
    real(dp), intent(in) :: commutator(:,0:,:,:)
  
    
  end subroutine TDIISMixer_mix


  !> Makes the mixer ready for a new SCC cycle.
  subroutine TDIISMixer_reset(this)

    !> DIIS mixer object
    type(TDIISMixer), intent(inout) :: this
    
    this%iIter = 0

    this%pot_storage = 0.0_dp
    this%comm_storage = 0.0_dp

    this%B = 0.0_dp
    this%c_vec = 0.0_dp
    this%B(1,:) = -1.0_dp
    this%B(:,1) = -1.0_dp
    this%B(1, 1) = 0.0_dp
    this%c_vec(1) = -1.0_dp
    
  end subroutine TDIISMixer_reset

end module diis