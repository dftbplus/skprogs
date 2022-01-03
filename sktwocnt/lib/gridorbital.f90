!> Module that implements a grid-type orbital.
module gridorbital

  use common_accuracy, only : dp
  use common_constants, only : pi
  use bisection, only : bisect
  use interpolation, only : polyinter, poly5zero

  implicit none
  private

  public :: TGridorb1, TGridorb1_init, TGridorb2, TGridorb2_init


  !> Contains the data of a grid function.
  type TGridorb1

    !> number of grid points
    integer :: nGrid

    !> r, f(r) values on grid
    real(dp), allocatable :: rvalues(:), fvalues(:)

  contains

    procedure :: getValue => TGridorb1_getValue
    procedure :: destruct => TGridorb1_destruct

  end type TGridorb1


  !> Contains the data of a grid function.
  type TGridorb2

    !> number of grid points
    integer :: nGrid

    !> r, f(r) values on grid
    real(dp), allocatable :: rvalues(:), fvalues(:)

    !> Gauss-Chebyshev pre-factor
    real(dp) :: delta

    !> cutoff radius at which the values f(r) shall vanish
    real(dp) :: rcut

  contains

    procedure :: getValue => TGridorb2_getValue
    procedure :: rescale => TGridorb2_rescale
    procedure :: destruct => TGridorb2_destruct

  end type TGridorb2


  !> Wraps around TGridorb1 pointer.
  type TGridorb1Wrap
    type(TGridorb1), pointer :: ptr => null()
  end type TGridorb1Wrap


  !> Wraps around TGridorb2 pointer.
  type TGridorb2Wrap
    type(TGridorb2), pointer :: ptr => null()
  end type TGridorb2Wrap


  real(dp), parameter :: distfudge = 1.0_dp
  real(dp), parameter :: deltar = 1e-04_dp

  integer, parameter :: ninter = 8
  integer, parameter :: nrightinter = 4

  integer, parameter :: npoint = 10000
  integer, parameter :: ninter2 = 4
  integer, parameter :: nrightinter2 = 2


contains

  !> Initializes a TGridorb1 grid-orbital.
  subroutine TGridorb1_init(this, rvals, fvals)

    !> initialised grid-orbital instance on exit
    type(TGridorb1), intent(out) :: this

    !> r, f(r) values on grid
    real(dp), intent(in) :: rvals(:), fvals(:)

    ! assert(size(values, dim=1) == 2)
    ! assert(size(values, dim=2) > 0)

    this%nGrid = size(rvals)

    this%rvalues = rvals
    this%fvalues = fvals

  end subroutine TGridorb1_init


  !> Destructs a TGridorb1 grid-orbital.
  subroutine TGridorb1_destruct(this)

    !> initialised grid-orbital instance to destruct
    class(TGridorb1), intent(inout) :: this

    if (allocated(this%rvalues)) deallocate(this%rvalues)
    if (allocated(this%fvalues)) deallocate(this%fvalues)

  end subroutine TGridorb1_destruct


  !> Delivers radial part of the orbital at the given distance.
  elemental function TGridorb1_getValue(this, rr) result(rad)

    !> grid-orbital instance
    class(TGridorb1), intent(in) :: this

    !> radius to calculate the value for
    real(dp), intent(in) :: rr

    !! radial part of the orbital at the given distance
    real(dp) :: rad

    !! auxiliary variables
    integer :: ind, iStart, iEnd
    real(dp) :: rmax, f0, f1, f2, f1p, f1pp

    ! sanity check
    ! if (this%nGrid < ninter + 1) then
    !   write(*,*) "Not enough points in the orbital grid!"
    !   stop
    ! end if

    ! find position of the point
    call bisect(this%rvalues, rr, ind, 1e-10_dp)
    rmax = this%rvalues(this%nGrid) + distfudge
    if (rr >= rmax) then
      ! outside of the region -> 0
      rad = 0.0_dp
    elseif (ind < this%nGrid) then
      ! before last gridpoint
      iEnd = min(this%nGrid, ind + nrightinter)
      iEnd = max(iEnd, ninter)
      iStart = iEnd - ninter + 1

      rad = polyinter(this%rvalues(iStart:iEnd), this%fvalues(iStart:iEnd), rr)
    else
      iEnd = this%nGrid
      iStart = iEnd - ninter + 1

      ! calculate 1st und 2nd derivatives at the end
      f1 = this%fvalues(iEnd)
      f0 = polyinter(this%rvalues(iStart:iEnd), this%fvalues(iStart:iEnd),&
          & this%rvalues(iEnd) - deltar)
      f2 = polyinter(this%rvalues(iStart:iEnd), this%fvalues(iStart:iEnd),&
          & this%rvalues(iEnd) + deltar)

      ! 1st order central finite difference --> 1st derivative
      f1p = (f2 - f0) / (2.0_dp * deltar)
      ! 2nd order central finite difference --> 2nd derivative
      f1pp = (f2 + f0 - 2.0_dp * f1) / deltar**2

      rad = poly5zero(f1, f1p, f1pp, rr - rmax, - 1.0_dp * distfudge)
    end if

  end function TGridorb1_getValue


  !> Initializes a TGridorb2 grid-orbital.
  subroutine TGridorb2_init(this, rvals, fvals)

    !> initialised grid-orbital instance on exit
    type(TGridorb2), intent(out) :: this

    !> r, f(r) values on grid
    real(dp), intent(in) :: rvals(:), fvals(:)

    !! grid-orbital instance
    type(TGridorb1) :: orb

    !! Gauss-Chebyshev abscissas and inverse Becke radii
    real(dp) :: xx, rr

    !! auxiliary variable
    integer :: ii

    ! assert(size(values, dim=1) == 2)
    ! assert(size(values, dim=2) > 0)

    call TGridorb1_init(orb, rvals, fvals)

    this%nGrid = npoint

    allocate(this%rvalues(this%nGrid))
    allocate(this%fvalues(this%nGrid))

    ! Gauss-Chebyshev pre-factor
    this%delta = pi / real(this%nGrid + 1, dp)

    do ii = 1, this%nGrid
      ! Gauss-Chebyshev abscissas
      xx = cos(this%delta * real(ii, dp))

      ! inverse Becke radius?
      rr = (1.0_dp - xx) / (1.0_dp + xx)
      this%rvalues(ii) = rr
      this%fvalues(ii) = orb%getValue(rr)
    end do

    ! cutoff radius at which the values f(r) shall vanish
    this%rcut = this%rvalues(this%nGrid) + distfudge

    call orb%destruct()

  end subroutine TGridorb2_init


  !> Destructs a TGridorb2 grid-orbital.
  subroutine TGridorb2_destruct(this)

    !> initialised grid-orbital instance to destruct
    class(TGridorb2), intent(inout) :: this

    if (allocated(this%fvalues)) deallocate(this%fvalues)

  end subroutine TGridorb2_destruct


  !> Delivers radial part of the orbital at the given distance.
  elemental function TGridorb2_getValue(this, rr) result(rad)

    !> grid-orbital instance
    class(TGridorb2), intent(in) :: this

    !> radius to calculate the value for
    real(dp), intent(in) :: rr

    !! radial part of the orbital at the given distance
    real(dp) :: rad

    !! auxiliary variables
    integer :: ind, iStart, iEnd
    real(dp) :: rmax, f0, f1, f2, f1p, f1pp
    real(dp) :: xx

    if (rr > this%rcut) then
      rad = 0.0_dp
    end if

    ! abscissa
    xx = (1.0_dp - rr) / (1.0_dp + rr)

    ! abscissa index
    ind = floor(acos(xx) / this%delta)

    if (ind < this%nGrid) then
      iEnd = min(this%nGrid, ind + nrightinter2)
      iEnd = max(iEnd, ninter2)
      iStart = iEnd - ninter2 + 1
      rad = polyinter(this%rvalues(iStart:iEnd), this%fvalues(iStart:iEnd), rr)
    else
      iEnd = this%nGrid
      iStart = iEnd - ninter2 + 1

      ! calculate 1st und 2nd derivatives at the end
      f1 = this%fvalues(iEnd)
      f0 = polyinter(this%rvalues(iStart:iEnd), this%fvalues(iStart:iEnd),&
          & this%rvalues(iEnd) - deltar)
      f2 = polyinter(this%rvalues(iStart:iEnd), this%fvalues(iStart:iEnd),&
          & this%rvalues(iEnd) + deltar)

      ! 1st order central finite difference --> 1st derivative
      f1p = (f2 - f0) / (2.0_dp * deltar)
      ! 2nd order central finite difference --> 2nd derivative
      f1pp = (f2 + f0 - 2.0_dp * f1) / deltar**2

      rad = poly5zero(f1, f1p, f1pp, rr - rmax, - 1.0_dp * distfudge)
    end if

  end function TGridorb2_getValue


  !> Rescales stored values f(r) of a grid-orbital instance.
  subroutine TGridorb2_rescale(this, fac)

    !> grid-orbital instance
    class(TGridorb2), intent(inout) :: this

    !> rescaling factor for f(r) values
    real(dp), intent(in) :: fac

    this%fvalues = this%fvalues * fac

  end subroutine TGridorb2_rescale

end module gridorbital
