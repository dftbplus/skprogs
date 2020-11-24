!> Implements a grid-type orbital.
module gridorbital
  use accuracy
  use constants
  use bisection
  use interpolation
  implicit none
  private

  public :: gridorb, gridorb2, init, destruct, getvalue, rescale

  !> Contains the data of a grid function.
  type gridorb
    integer :: ngrid
    real(dp), allocatable :: rvalues(:), fvalues(:)
  end type gridorb

  type gridorb2
    integer :: ngrid
    real(dp), allocatable :: rvalues(:), fvalues(:)
    real(dp) :: delta, rcut
  end type gridorb2

  type gridorb_wrap
    type(gridorb), pointer :: ptr => null()
  end type gridorb_wrap

  type gridorb2_wrap
    type(gridorb2), pointer :: ptr => null()
  end type gridorb2_wrap

  interface init
    module procedure gridorb_init
    module procedure gridorb2_init
  end interface  

  interface destruct
    module procedure gridorb_destruct
    module procedure gridorb2_destruct
  end interface

  interface getvalue
    module procedure gridorb_getvalue
    module procedure gridorb2_getvalue
  end interface

  interface rescale
    module procedure gridorb2_rescale
  end interface

  real(dp), parameter :: distfudge = 1.0_dp
  integer, parameter :: ninter = 8
  integer, parameter :: nrightinter = 4
  real(dp), parameter :: deltar = 1e-4_dp

  integer, parameter :: npoint = 10000
  !real(dp), parameter :: tol = 1e-12_dp
  integer, parameter :: ninter2 = 4
  integer, parameter :: nrightinter2 = 2

contains

  !> Initializes the grid orbital.
  !! \param self  initialised instance on exit.
  !! \param values r,f(r) values for the grid
  subroutine gridorb_init(self, rvals, fvals)
    type(gridorb), intent(inout) :: self
    real(dp), intent(in) :: rvals(:), fvals(:)

    !assert(size(values, dim=1) == 2)
    !assert(size(values, dim=2) > 0)

    self%ngrid = size(rvals)
    allocate(self%rvalues(self%ngrid))
    allocate(self%fvalues(self%ngrid))
    self%rvalues = rvals(:)
    self%fvalues = fvals(:)
    
  end subroutine gridorb_init

  
  !> Destructs the instance.
  !! \param self instance.
  subroutine gridorb_destruct(self)
    type(gridorb), intent(inout) :: self
    
    deallocate(self%rvalues)
    deallocate(self%fvalues)
    
  end subroutine gridorb_destruct

  
  !> Delivers the value of the orbital
  !! \param self instance.
  !! \param rr radius at which to calculate the value.
  !! \return rad radial part of the orbital at the given distance.
  elemental function gridorb_getvalue(self, rr) result(rad)
    type(gridorb), intent(in) :: self
    real(dp), intent(in) :: rr
    real(dp) :: rad

    integer :: ind, istart, iend
    real(dp) :: rmax, f0, f1, f2, f1p, f1pp

    ! sanity check
    !if (self%ngrid < ninter + 1) then
    !  write (*,*) "not enough points in the orbital grid!"
    !  stop
    !end if

    ! Find position of the point
    call bisect(self%rvalues, rr, ind, 1e-10_dp)
    rmax = self%rvalues(self%ngrid) + distfudge
    if (rr >=  rmax) then
      ! outside of the region -> 0
      rad = 0.0_dp
    elseif (ind < self%ngrid) then
      ! before last gridpoint
      iend = min(self%ngrid, ind + nrightinter)
      iend = max(iend, ninter)
      istart = iend - ninter + 1
      rad = polyinter(self%rvalues(istart:iend), self%fvalues(istart:iend), rr)
    else
      iend = self%ngrid
      istart = iend - ninter + 1
      ! calculate 1st und 2nd derivatives at the end
      f1 = self%fvalues(iend)
      f0 = polyinter(self%rvalues(istart:iend), self%fvalues(istart:iend), &
          &self%rvalues(iend) - deltar)
      f2 = polyinter(self%rvalues(istart:iend), self%fvalues(istart:iend), &
          &self%rvalues(iend) + deltar)
      f1p = (f2 - f0) / (2.0_dp * deltar)
      f1pp = (f2 + f0 - 2.0_dp * f1) / deltar**2
      rad = poly5zero(f1, f1p, f1pp, rr - rmax, -1.0_dp * distfudge)
    end if
        
  end function gridorb_getvalue



  !> Initializes the grid orbital.
  !! \param self  initialised instance on exit.
  !! \param values r,f(r) values for the grid
  subroutine gridorb2_init(self, rvals, fvals)
    type(gridorb2), intent(inout) :: self
    real(dp), intent(in) :: rvals(:), fvals(:)

    type(gridorb) :: orb
    real(dp) :: xx, rr
    integer :: ii

    !assert(size(values, dim=1) == 2)
    !assert(size(values, dim=2) > 0)

    call init(orb, rvals, fvals)
    self%ngrid = npoint
    allocate(self%rvalues(self%ngrid))
    allocate(self%fvalues(self%ngrid))
    self%delta = pi / real(self%ngrid + 1, dp)
    do ii = 1, self%ngrid
      xx = cos(self%delta * real(ii, dp))
      rr = (1.0_dp - xx) / (1.0_dp + xx)
      self%rvalues(ii) = rr
      self%fvalues(ii) = getvalue(orb, rr)
    end do
    self%rcut = self%rvalues(self%ngrid) + distfudge
    call destruct(orb)
    
  end subroutine gridorb2_init

  
  !> Destructs the instance.
  !! \param self instance.
  subroutine gridorb2_destruct(self)
    type(gridorb2), intent(inout) :: self
    
    deallocate(self%fvalues)
    
  end subroutine gridorb2_destruct

  
  !> Delivers the value of the orbital
  !! \param self instance.
  !! \param rr radius at which to calculate the value.
  !! \return rad radial part of the orbital at the given distance.
  elemental function gridorb2_getvalue(self, rr) result(rad)
    type(gridorb2), intent(in) :: self
    real(dp), intent(in) :: rr
    real(dp) :: rad

    integer :: ind, istart, iend
    real(dp) :: rmax, f0, f1, f2, f1p, f1pp
    real(dp) :: xx

    if (rr > self%rcut) then
      rad = 0.0_dp
    end if
    xx = (1.0_dp - rr) / (1.0_dp + rr)
    ind = floor(acos(xx) / self%delta)
    if (ind < self%ngrid) then
      iend = min(self%ngrid, ind + nrightinter2)
      iend = max(iend, ninter2)
      istart = iend - ninter2 + 1
      rad = polyinter(self%rvalues(istart:iend), self%fvalues(istart:iend), rr)
    else
      iend = self%ngrid
      istart = iend - ninter2 + 1
      ! calculate 1st und 2nd derivatives at the end
      f1 = self%fvalues(iend)
      f0 = polyinter(self%rvalues(istart:iend), self%fvalues(istart:iend), &
          &self%rvalues(iend) - deltar)
      f2 = polyinter(self%rvalues(istart:iend), self%fvalues(istart:iend), &
          &self%rvalues(iend) + deltar)
      f1p = (f2 - f0) / (2.0_dp * deltar)
      f1pp = (f2 + f0 - 2.0_dp * f1) / deltar**2
      rad = poly5zero(f1, f1p, f1pp, rr - rmax, -1.0_dp * distfudge)
    end if
      
  end function gridorb2_getvalue


  subroutine gridorb2_rescale(self, fac)
    type(gridorb2), intent(inout) :: self
    real(dp), intent(in) :: fac

    self%fvalues = self%fvalues * fac
    
  end subroutine gridorb2_rescale


end module gridorbital
