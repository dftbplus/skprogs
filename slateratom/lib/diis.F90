#:include 'common.fypp'

module diis
  use common_accuracy, only : dp
  use lapackroutines, only : getrf, getrs

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
    real(dp), allocatable :: Bmat(:,:)

    contains

end type TDIISMixer

contains
  
  !> Creates a DIIS mixer instance
  !! TODO: complete description
  subroutine TDIISMixer_init(this, mVec, mIter, mixParam, max_l, problemsize)
    
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

    !> Max. size of diagonalization problem
    integer, intent(in) :: problemsize

    @:ASSERT(mVec > 0)
    @:ASSERT(mixParam > 0.0_dp)
    
    this%iIter = 0
    this%mVec = mVec
    this%mIter = mIter
    this%mixParam = mixParam

    allocate(this%pot_storage(mVec, 2, 0:max_l, problemsize, problemsize))
    allocate(this%comm_storage(mVec, 2, 0:max_l, problemsize, problemsize))
    allocate(this%Bmat(mVec+1, mVec+1))

    this%pot_storage = 0.0_dp
    this%comm_storage = 0.0_dp

    this%Bmat = 0.0_dp
    this%Bmat(1,:) = -1.0_dp
    this%Bmat(:,1) = -1.0_dp
    this%Bmat(1, 1) = 0.0_dp
  
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

    !! auxiliary variables
    integer :: jj, ii, iSpin, ll, diag_ind, problemsize, max_l, vec_lim
    integer, allocatable :: ipiv(:)
    real(dp) :: e_metric
    real(dp), allocatable :: aux_mat(:,:), Bmat_copy(:,:), c_vec(:)

    @:ASSERT(size(inpResult, dim=4) .eq. size(this%pot_storage, dim=5))
    @:ASSERT(size(inpResult, dim=3) .eq. size(this%pot_storage, dim=4))
    @:ASSERT(size(inpResult, dim=2) .eq. size(this%pot_storage, dim=3))
    @:ASSERT(size(inpResult, dim=1) .eq. size(this%pot_storage, dim=2))
    

    this%iIter = this%iIter + 1
    if (this%iIter > this%mIter) then
      error stop "DIIS mixer: Maximal nr. of steps exceeded"
    end if

    problemsize = size(this%pot_storage, dim=5)
    max_l = size(this%pot_storage, dim=3) - 1

    allocate(aux_mat(problemsize, problemsize))
    allocate(ipiv(this%mVec + 1))
    aux_mat = 0.0_dp

    ! Move rows in storages and the B matrix
    jj = this%mVec - 1
    do while (jj > 0)
      this%pot_storage(jj + 1, :, 0:, :, :) = this%pot_storage(jj, :, 0:, :, :)
      this%comm_storage(jj + 1, :, 0:, :, :) = this%comm_storage(jj, :, 0:, :, :)
      this%Bmat(jj + 2, :) = this%Bmat(jj + 1, :)
      jj = jj - 1
    end do
    
    ! Update storages
    this%pot_storage(1, :, 0:, :, :) = inpResult(:,0:,:,:)
    this%comm_storage(1, :, 0:, :, :) = commutator(:,0:,:,:)

    ! Update the row in the B matrix if this%mVec > 5;
    ! compute the whole B matrix if this%mVec = 5
    if (this%iIter >= 5) then
  
      vec_lim = 1
      if (this%iIter == this%mVec) then
        vec_lim = this%mVec
      end if

      do ii = 1, vec_lim
        do jj = 1, this%mVec
          e_metric = 0.0_dp
          do iSpin = 1, 2
            do ll = 0, max_l
              ! Compute scalar product between commutators as
              ! Tr(<A|T(B)>)
              ! TODO: is scalar product being computed correctly?
              aux_mat = 0.0_dp
              aux_mat = matmul(this%comm_storage(ii, iSpin, ll, :, :),&
                & transpose(this%comm_storage(jj, iSpin, ll, :, :)))
              e_metric = e_metric + sum( (/ (aux_mat(diag_ind,diag_ind), diag_ind=1, size(aux_mat, 1)) /) )
            end do
          end do
          this%Bmat(ii + 1, jj + 1) = e_metric
        end do
      end do
    end if

    ! Do a simple mixer step until the full B matrix is built
    if (this%iIter < this%mVec) then

      inpResult = inpResult + this%mixParam * diff
    
    ! Solve DIIS problem to get coefficients
    else
      
      ! TODO: check if solving for coefficients behaves correctly
      ! Copy B-matrix to prevent rewriting it in place
      allocate(Bmat_copy, mold=this%Bmat)
      allocate(c_vec(this%mVec + 1))
      Bmat_copy = transpose(this%Bmat)
      write(*, *) "aaaa:", Bmat_copy
      c_vec = 0.0_dp
      c_vec(1) = -1.0_dp

      ! Compute LU decomposition
      call getrf(Bmat_copy, ipiv)

      ! Solve for coefficients
      ! TODO: is transposing necessary?
      call getrs(Bmat_copy, ipiv, c_vec, trans='t')
      write(*, *) "aaaa:", Bmat_copy

      inpResult = 0.0_dp
      do jj = 2, this%mVec + 1
        inpResult = inpResult + c_vec(jj) * this%pot_storage(jj - 1, :, 0:, :, :)
      end do

      deallocate(Bmat_copy)
      deallocate(c_vec)

    end if

    deallocate(aux_mat)
    deallocate(ipiv)
  
  end subroutine TDIISMixer_mix


  !> Makes the mixer ready for a new SCF cycle.
  subroutine TDIISMixer_reset(this)

    !> DIIS mixer object
    type(TDIISMixer), intent(inout) :: this
    
    this%iIter = 0

    this%pot_storage = 0.0_dp
    this%comm_storage = 0.0_dp

    this%Bmat = 0.0_dp
    this%Bmat(1,:) = -1.0_dp
    this%Bmat(:,1) = -1.0_dp
    this%Bmat(1, 1) = 0.0_dp
    
  end subroutine TDIISMixer_reset

end module diis