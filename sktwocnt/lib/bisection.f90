!> Module that contains routines to locate a value in an array using bisection.
module bisection

  use common_accuracy, only : dp

  implicit none
  private

  public :: bisect

  !> Bisection driver that interfaces integer- and real-valued array routines.
  interface bisect
    module procedure bisect_real
    module procedure bisect_int
  end interface bisect


contains

  !> Real case for bisection search to to find a point in an array xx(:) between xx(1) and
  !! xx(size(xx)) such that element indexed ind is less than the value x0 queried.
  pure subroutine bisect_real(xx, x0, ind, tol)

    !> array of values in monotonic order to search through
    real(dp), intent(in) :: xx(:)

    !> value to locate ind for
    real(dp), intent(in) :: x0

    !> located element such that xx(ind) < x0 < xx(ind)
    integer, intent(out) :: ind

    !> optional, user-specified tolerance for comparisons
    real(dp), intent(in), optional :: tol

    !! length of array to search
    integer :: nn

    !! lower, upper and current value index
    integer :: iLower, iUpper, iCurr

    !! actual tolerance selected
    real(dp) :: rTol

    !! true, if xx(:) is in ascending ordering
    logical :: tAscending

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
      tAscending = (xx(nn) >= xx(1))
      iLower = 0
      iCurr = nn + 1
      do while ((iCurr - iLower) > 1)
        iUpper = (iCurr + iLower) / 2
        if (tAscending .eqv. (x0 >= xx(iUpper) + rTol)) then
          iLower = iUpper
        else
          iCurr = iUpper
        end if
      end do
      ind = iLower
    end if

  end subroutine bisect_real


  !> Integer case for bisection search to to find a point in an array xx(:) between xx(1) and
  !! xx(size(xx)) such that element indexed ind is less than the value x0 queried.
  pure subroutine bisect_int(xx, x0, ind)

    !> array of values in monotonic order to search through
    integer, intent(in) :: xx(:)

    !> value to locate ind for
    integer, intent(in) :: x0

    !> located element such that xx(ind) < x0 < xx(ind)
    integer, intent(out) :: ind

    !! length of array to search
    integer :: nn

    !! lower, upper and current value index
    integer :: iLower, iUpper, iCurr

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
      iLower = 0
      iCurr = nn + 1
      do while ((iCurr - iLower) > 1)
        iUpper = (iCurr + iLower) / 2
        if ((xx(nn) >= xx(1)) .eqv. (x0 >= xx(iUpper))) then
          iLower = iUpper
        else
          iCurr = iUpper
        end if
      end do
      ind = iLower
    end if

  end subroutine bisect_int

end module bisection
