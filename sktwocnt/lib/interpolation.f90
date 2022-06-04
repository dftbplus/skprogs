!> Module that contains routines for inter- and extrapolation.
module interpolation

  use common_accuracy, only : dp

  implicit none
  private

  public :: poly5zero, spline3_free, polyinter


contains

  !> Returns the value of a polynomial of 5th degree at x.
  !! \details The polynomial is created with the following boundary conditions:
  !!   Its value, its 1st and 2nd derivatives are zero at x = 0 and agree with the provided values
  !!   at x = dx.
  pure function poly5zero(y0, y0p, y0pp, xx, dx) result(yy)

    !> value of the polynom at x = dx
    real(dp), intent(in) :: y0

    !> value of the 1st derivative at x = dx
    real(dp), intent(in) :: y0p

    !> value of the 2nd derivative at x = dx
    real(dp), intent(in) :: y0pp

    !> point where the polynomial should be calculated
    real(dp), intent(in) :: xx

    !> point, where the polynomials value and first two derivatives should take the provided values
    real(dp), intent(in) :: dx

    !! value of the polynomial at xx
    real(dp) :: yy

    real(dp) :: dx1, dx2, cc, bb, aa, xr

    ! f(x) = ax^5 + bx^4 + cx^3 + dx^2 + ex + f
    ! f(0) = 0, f'(0) = 0, f''(0) = 0 --> d = e = f = 0

    dx1 = y0p * dx
    dx2 = y0pp * dx**2

    ! c * (dx)**3
    cc = 10.0_dp * y0 - 4.0_dp * dx1 + 0.5_dp * dx2

    ! b * (dx)**4
    bb = - 15.0_dp * y0 + 7.0_dp * dx1 - 1.0_dp * dx2

    ! a * (dx)**5
    aa = 6.0_dp * y0 - 3.0_dp * dx1 + 0.5_dp * dx2

    xr = xx / dx
    yy = ((aa * xr + bb) * xr + cc) * xr**3

  end function poly5zero


  !! Returns the value of a free spline at a certain point.
  !! \details The spline is created with the following boundary conditions:
  !!   Its value, 1st and 2nd derivatives agree with the provided values at
  !!   x = 0 and its value agrees with the provided value at x = dx.
  !! \note If you want the value for a derivative, you have to query them both.
  pure subroutine spline3_free(y0, y0p, y0pp, dx, ydx, xx, yy, yp, ypp)

    !> function value at x = 0
    real(dp), intent(in) :: y0

    !> first derivative at x = 0
    real(dp), intent(in) :: y0p

    !> second derivative at x = 0
    real(dp), intent(in) :: y0pp

    !> function value at dx
    real(dp), intent(in) :: ydx

    !> second fitting point
    real(dp), intent(in) :: dx

    !> point to interpolate
    real(dp), intent(in) :: xx

    !> value of the 3rd order polynomial at xx
    real(dp), intent(out), optional :: yy

    !> first derivative at xx
    real(dp), intent(out), optional :: yp

    !> second derivative at xx
    real(dp), intent(out), optional :: ypp

    !! spline coefficients
    real(dp) :: aa, bb, cc, dd

    !! reciprocal second fitting point
    real(dp) :: dx1

    ! assert(present(yp) .eqv. present(ypp))

    dx1 = 1.0_dp / dx

    aa = y0
    bb = y0p
    cc = 0.5_dp * y0pp
    dd = (((ydx - y0) * dx1 - y0p) * dx1 - 0.5_dp * y0pp) * dx1

    if (present(yy)) then
      yy = ((dd * xx + cc) * xx + bb) * xx + aa
    end if

    if (present(yp)) then
      yp = (3.0_dp * dd * xx + 2.0_dp * cc) * xx + bb
      ypp = 6.0_dp * dd * xx + 2.0_dp * cc
    end if

  end subroutine spline3_free


  !> Polynomial interpolation through given points.
  !! \note The algorithm is based on the Numerical recipes.
  pure function polyinter(xp, yp, xx) result(yy)

    !> x-coordinates of the fit points
    real(dp), intent(in) :: xp(:)

    !> y-coordinates of the fit points
    real(dp), intent(in) :: yp(:)

    !> point, where the polynomial should be evaluated
    real(dp), intent(in) :: xx

    !! value of the polynomial
    real(dp) :: yy

    !! number of interpolation abscissas
    integer :: nn

    !! auxiliary variables
    integer :: icl, ii, mm
    real(dp) :: cc(size(xp)), dd(size(xp))
    real(dp) :: dx, dxnew, dyy, rtmp

    nn = size(xp)

    ! assert(nn > 1)
    ! assert(size(yp) == nn)

    cc(:) = yp
    dd(:) = yp
    icl = 1
    dx = abs(xx - xp(icl))
    do ii = 2, nn
      dxnew = abs(xx - xp(ii))
      if (dxnew < dx) then
        icl = ii
        dx = dxnew
      end if
    end do
    yy = yp(icl)
    icl = icl - 1
    do mm = 1, nn - 1
      do ii = 1, nn - mm
        rtmp = xp(ii) - xp(ii + mm)
        ! if (abs(rtmp) < epsilon(1.0_dp)) then
        !   write(*,*) "Polint failed"
        !   stop
        ! end if
        rtmp = (cc(ii + 1) - dd(ii)) / rtmp
        cc(ii) = (xp(ii) - xx) * rtmp
        dd(ii) = (xp(ii + mm) - xx) * rtmp
      end do
      if (2 * icl < nn - mm) then
        dyy = cc(icl + 1)
      else
        dyy = dd(icl)
        icl = icl - 1
      end if
      yy = yy + dyy
    end do

  end function polyinter

end module interpolation
