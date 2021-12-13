!> Contains routines to locate a value in an array using bisection.
module bisection

  use common_accuracy, only: dp

  implicit none
  private

  public :: bisect

  !> Bisection driver.
  interface bisect
    module procedure bisect_real
    module procedure bisect_int
  end interface bisect

contains

  !> Real case for bisection search to to find a point in an array xx(:)
  !! between xx(1) and xx(size(xx)) such that element indexed ind is less than
  !! the value x0 queried.
  !! \param xx Array of values in monotonic order to search through.
  !! \param x0 Value to locate ind for.
  !! \param ind Located element such that xx(ind) < x < xx(ind).
  pure subroutine bisect_real(xx, x0, ind, tol)
    real(dp), intent(in) :: xx(:), x0
    integer, intent(out) :: ind
    real(dp), intent(in), optional :: tol

    integer :: nn
    integer :: ilower, iupper, icurr
    real(dp) :: rTol      ! real tolerance
    logical :: ascending

    nn = size(xx)
    if (nn == 0) then
      ind = 0
      return
    end if

    if (present(tol)) then
      rTol = tol
    else
      rTol = epsilon(0.0_dp)
    end if

    if (x0 < xx(1) - rTol) then
      ind = 0
    else if (abs(x0 - xx(1)) <= rTol) then
      ind = 1
    else if (abs(x0 - xx(nn)) <= rTol) then
      ind = nn - 1
    else if (x0 > xx(nn) + rTol) then
      ind = nn
    else
      ascending = (xx(nn) >= xx(1))
      ilower = 0
      icurr = nn + 1
      do while ((icurr - ilower) > 1)
        iupper = (icurr + ilower) / 2
        if (ascending .eqv. (x0 >= xx(iupper) + rTol)) then
          ilower = iupper
        else
          icurr = iupper
        end if
      end do
      ind = ilower
    end if

  end subroutine bisect_real

  !> Integer case for bisection search to to find a point in an array xx(:)
  !! between xx(1) and xx(size(xx)) such that element indexed ind is less than
  !! the value x0 queried
  !! \param xx Array of values in monotonic order to search through.
  !! \param x0 Value to locate ind for.
  !! \param ind Located element such that xx(ind) < x < xx(ind).
  pure subroutine bisect_int(xx, x0, ind)
    integer, intent(in) :: xx(:), x0
    integer, intent(out) :: ind

    integer :: nn
    integer :: ilower, iupper, icurr

    nn = size(xx)
    if (nn == 0) then
      ind = 0
      return
    end if

    if (x0 < xx(1)) then
      ind = 0
    else if (x0 == xx(1)) then
      ind = 1
    else if (x0 == xx(nn)) then
      ind = nn - 1
    else if (x0 > xx(nn)) then
      ind = nn
    else
      ilower = 0
      icurr = nn + 1
      do while ((icurr - ilower) > 1)
        iupper = (icurr + ilower) / 2
        if ((xx(nn) >= xx(1)) .eqv. (x0 >= xx(iupper))) then
          ilower = iupper
        else
          icurr = iupper
        end if
      end do
      ind = ilower
    end if

  end subroutine bisect_int

end module bisection
