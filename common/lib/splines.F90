!> Splines are fully specified by the interpolation points, except that at the ends, we have the
!! freedom to prescribe the second derivatives. If we know a derivative at an end (exactly), then
!! it is best to impose that. Otherwise, it is better to use the "consistent" end conditions:
!! the second derivative is determined such that it is smooth.
!!
!! High level API: spline3, spline3ders
!! Low level API: the rest of public subroutines
!!
!! Use the high level API to obtain cubic spline fit with consistent boundary conditions and
!! optionally the derivatives. Use the low level API if more fine grained control is needed.
!!
!! This module is based on a code written by John E. Pask, LLNL.
module common_splines

  use common_accuracy, only : dp
  use common_utils, only: stop_error

  implicit none
  private

  public :: spline3pars, spline3valder, iix, iixmin, iixun, iixexp, poly3, dpoly3, d2poly3,&
      & spline3, spline3ders


contains

  !> Takes the function values y on the grid x and returns new values ynew at the given grid xnew
  !! using cubic splines interpolation with such boundary conditions so that the 2nd derivative is
  !! consistent with the interpolating cubic.
  function spline3(xx, yy, xnew) result(ynew)

    !> tabulated abscissas of function
    real(dp), intent(in) :: xx(:)

    !> tabulated values of function
    real(dp), intent(in) :: yy(:)

    !> evaluation points for interpolation
    real(dp), intent(in) :: xnew(:)

    !> interpolated values at x = xnew
    real(dp) :: ynew(size(xnew))

    !! parameters defining spline: c(i,j) = ith parameter of jth spline polynomial,
    !! p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data points
    !! dimensions: c(0:4, 1:n-1)
    real(dp) :: cc(0:4, size(xx) - 1)

    !! iterates over all evaluation points
    integer :: ii

    !!
    integer :: ip

    ! get spline parameters: 2nd derivs at ends determined by cubic interpolation
    call spline3pars(xx, yy, [2, 2], [0.0_dp, 0.0_dp], cc)

    ip = 0
    do ii = 1, size(xnew)
      ip = iixmin(xnew(ii), xx, ip)
      ynew(ii) = poly3(xnew(ii), cc(:, ip))
    end do

  end function spline3


  !> Just like spline3(), but also calculate 1st and 2nd derivatives.
  subroutine spline3ders(xx, yy, xnew, ynew, dynew, d2ynew)

    !> tabulated abscissas of function
    real(dp), intent(in) :: xx(:)

    !> tabulated values of function
    real(dp), intent(in) :: yy(:)

    !> evaluation points for interpolation
    real(dp), intent(in) :: xnew(:)

    !> interpolated values at x = xnew
    real(dp), intent(out), optional :: ynew(:)

    !> 1st derivative at interpolated values at x = xnew
    real(dp), intent(out), optional :: dynew(:)

    !> 2nd derivative at interpolated values at x = xnew
    real(dp), intent(out), optional :: d2ynew(:)

    !! parameters defining spline: c(i,j) = ith parameter of jth spline polynomial,
    !! p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data points
    !! dimensions: c(0:4, 1:n-1)
    real(dp) :: cc(0:4, size(xx) - 1)

    !! iterates over evaluation points
    integer :: ii

    !! index ip of interval [xx(ip), xx(ip+1)] containing xnew(ii)
    integer :: ip

    call spline3pars(xx, yy, [2, 2], [0.0_dp, 0.0_dp], cc)

    ip = 0
    do ii = 1, size(xnew)
      ip = iixmin(xnew(ii), xx, ip)
      if (present(ynew)) ynew(ii) = poly3(xnew(ii), cc(:, ip))
      if (present(dynew)) dynew(ii) = dpoly3(xnew(ii), cc(:, ip))
      if (present(d2ynew)) d2ynew(ii) = d2poly3(xnew(ii), cc(:, ip))
    end do

  end subroutine spline3ders


  !> Returns parameters c defining cubic spline interpolating x-y data xi, yi, with boundary
  !! conditions specified by bcytpe, bcvals.
  subroutine spline3pars(xi, yi, bctype, bcval, cc)

    !> x values of data
    real(dp), intent(in) :: xi(:)

    !> y values of data
    real(dp), intent(in) :: yi(:)

    !> type of boundary condition at each end:
    !! bctype(1) = type at left end; bctype(2) = type at right end
    !! 1 = specified 2nd derivative, 2 = 2nd derivative consistent with interpolating cubic.
    integer, intent(in) :: bctype(2)

    !> boundary condition values at each end:
    !! bcval(1) = value at left end; bcval(2) = value at right end
    real(dp), intent(in) :: bcval(2)

    !> parameters defining spline: c(i,j) = ith parameter of jth spline polynomial,
    !! p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data points
    !! dimensions: c(0:4, 1:n-1)
    real(dp), intent(out) :: cc(0:,:)

    !! spline eq. matrix -- LAPACK band form
    real(dp) :: as(5, 2 * size(cc, 2))

    !! spline eq. rhs vector
    real(dp) :: bs(2 * size(cc, 2))

    !! spline eq. solution vector
    real(dp) :: cs(2 * size(cc, 2))

    !! spline intervals
    real(dp) :: hi(size(cc, 2))

    !! end-cubic eq. matrix
    real(dp) :: ae(4, 4)

    !! end-cubic eq. rhs vector
    real(dp) :: be(4)

    !! end-cubic eq. solution vector
    real(dp) :: ce(4)

    !! x, y values at ends
    real(dp) :: xe(4), ye(4)

    !! 2nd derivatives at ends
    real(dp) :: d2p1, d2pn

    !! expansion center
    real(dp) :: x0

    !! expansion coefficients
    real(dp) :: c1, c2, c3, c4

    !! number of data points
    integer :: nn

    !!
    integer :: ii, jj, i2

    !! lapack variables
    integer :: ipiv(4), ipiv2(2 * size(cc, 2))
    real(dp) :: bemat(4, 1), bmat(2 * size(cc, 2), 1)
    integer :: info

    ! check input parameters
    if (bctype(1) < 1 .or. bctype(1) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
    if (bctype(2) < 1 .or. bctype(2) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
    if (size(cc, 1) /= 5) call stop_error("spline3pars error: size(c,1) /= 5.")
    if (size(cc, 2) /= size(xi) - 1) call stop_error("spline3pars error: size(c,2) /= size(xi)-1.")
    if (size(xi) /= size(yi)) call stop_error("spline3pars error: size(xi) /= size(yi)")

    ! get rid of compiler warnings
    d2p1 = 0
    d2pn = 0

    ! initializations
    nn = size(xi)
    do ii = 1, nn - 1
      hi(ii) = xi(ii + 1) - xi(ii)
    end do

    ! compute interpolating-cubic 2nd derivs at ends, if required left end
    if (bctype(1) == 2) then
      if (nn < 4) call stop_error("spline3pars error: n < 4")
      xe(:) = xi(1:4)
      ye(:) = yi(1:4)
      x0 = xe(1) ! center at end
      do ii = 1, 4
        do jj = 1, 4
          ae(ii,jj) = (xe(ii) - x0)**(jj - 1)
        end do
      end do
      ae(:, 1) = 1 ! set 0^0 = 1
      be(:) = ye
      bemat(:, 1) = be
      call dgesv(4, 1, ae, 4, ipiv, bemat, 4, info)
      if (info /= 0) call stop_error("spline3pars error: dgesv error.")
      ce = bemat(:, 1)
      d2p1 = 2 * ce(3)
    end if
    ! right end
    if (bctype(2) == 2) then
      if (nn < 4) call stop_error("spline3pars error: n < 4")
      xe(:) = xi(nn-3:nn)
      ye(:) = yi(nn-3:nn)
      x0 = xe(4) ! center at end
      do ii = 1, 4
        do jj = 1, 4
          ae(ii, jj) = (xe(ii) - x0)**(jj - 1)
        end do
      end do
      ae(:, 1) = 1 ! set 0^0 = 1
      be(:) = ye
      bemat(:, 1) = be
      call dgesv(4, 1, ae, 4, ipiv, bemat, 4, info)
      if (info /= 0) call stop_error("spline3pars error: dgesv error.")
      ce = bemat(:, 1)
      d2pn = 2 * ce(3)
    end if

    ! set 2nd derivs at ends
    if (bctype(1) == 1) d2p1 = bcval(1)
    if (bctype(2) == 1) d2pn = bcval(2)

    ! construct spline equations -- LAPACK band form
    ! basis: phi1 = -(x-x_i)/h_i, phi2 = (x-x_{i+1})/h_i, phi3 = phi1^3-phi1, phi4 = phi2^3-phi2
    ! on interval [x_i,x_{i+1}] of length h_i = x_{i+1}-x_i
    ! A(:,:) = 0  ! full matrix
    as(:,:) = 0
    ! left end condition
    as(4, 1) = 6 / hi(1)**2
    bs(1) = d2p1
    ! internal knot conditions
    do ii = 2, nn - 1
      i2 = 2 * (ii - 1)
      as(5, i2-1) = 1 / hi(ii - 1)
      as(4, i2) = 2 / hi(ii - 1)
      as(3, i2 + 1) = 2 / hi(ii)
      as(2, i2 + 2) = 1 / hi(ii)
      as(5, i2) = 1 / hi(ii - 1)**2
      as(4, i2 + 1) = -1 / hi(ii)**2
      bs(i2) = (yi(ii + 1) - yi(ii)) / hi(ii) - (yi(ii) - yi(ii - 1)) / hi(ii - 1)
      bs(i2 + 1) = 0
    end do
    ! right end condition
    as(4, 2 * (nn - 1)) = 6 / hi(nn - 1)**2
    bs(2 * (nn - 1)) = d2pn

    ! solve spline equations -- LAPACK band form
    bmat(:, 1) = bs
    call dgbsv(2 * (nn - 1), 1, 2, 1, as, 5, ipiv2, bmat, 2 * (nn - 1), info)
    if (info /= 0) call stop_error("spline3pars error: dgbsv error.")
    cs = bmat(:, 1)

    ! transform to (x-x0)^(i-1) basis and return
    do ii = 1, nn - 1
      ! coefficients in spline basis:
      c1 = yi(ii)
      c2 = yi(ii + 1)
      c3 = cs(2 * ii - 1)
      c4 = cs(2 * ii)
      ! coefficients in (x-x0)^(i-1) basis
      cc(0, ii) = xi(ii)
      cc(1, ii) = c1
      cc(2, ii) = -(c1 - c2 + 2 * c3 + c4) / hi(ii)
      cc(3, ii) = 3 * c3 / hi(ii)**2
      cc(4, ii) = (-c3 + c4) / hi(ii)**3
    end do

  end subroutine spline3pars


  !> Returns value and 1st derivative of spline defined by knots xi and parameters c returned by
  !! spline3pars.
  subroutine spline3valder(xx, xi, cc, val, der)

    !> point at which to evaluate spline
    real(dp), intent(in):: xx

    !> spline knots (x values of data)
    real(dp), intent(in):: xi(:)

    !> parameters defining spline: c(i,j) = ith parameter of jth spline polynomial,
    !! p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data points
    !! dimensions: c(0:4, 1:n-1)
    real(dp), intent(in):: cc(0:,:)

    !> value of spline at x
    real(dp), intent(out):: val

    !> 1st derivative of spline at x
    real(dp), intent(out):: der

    !! index i of interval [xi(i), xi(i+1)] containing x
    integer i1

    ! initialize, check input parameters
    if (size(cc, 1) /= 5) call stop_error("spline3 error: size(c,1) /= 5.")
    if (size(cc, 2) /= size(xi) - 1) call stop_error("spline3 error: size(c,2) /= size(xi)-1.")

    ! find interval containing x
    i1 = iix(xx, xi)

    ! return value and derivative
    val = poly3(xx, cc(:, i1))
    der = dpoly3(xx, cc(:, i1))

  end subroutine spline3valder


  !> Returns index i of interval [xi(i), xi(i+1)] containing x in mesh xi, with intervals indexed by
  !! left-most points.
  !! Note: x outside [x1,xn] are indexed to nearest end. Uses bisection, except if "x" lies in the
  !! first or second elements (which is often the case).
  function iix(xx, xi) result(i1)

    !> target value
    real(dp), intent(in) :: xx

    !> mesh, xi(i) < xi(i+1)
    real(dp), intent(in) :: xi(:)

    !> number of mesh points
    integer nn

    !> index i of interval [xi(i), xi(i+1)] containing x
    integer :: i1

    !!
    integer i2, ic

    nn = size(xi)
    i1 = 1
    if (nn < 2) then
      call stop_error("error in iix: nn < 2")
    elseif (nn == 2) then
      i1 = 1
    elseif (nn == 3) then
      if (xx <= xi(2)) then ! first element
        i1 = 1
      else
        i1 = 2
      end if
    elseif (xx <= xi(1)) then ! left end
      i1 = 1
    elseif (xx <= xi(2)) then ! first element
      i1 = 1
    elseif (xx <= xi(3)) then ! second element
      i1 = 2
    elseif (xx >= xi(nn)) then  ! right end
      i1 = nn - 1
    else
      ! bisection: xi(i1) <= x < xi(i2)
      i1 = 3
      i2 = nn
      do
        if (i2 - i1 == 1) exit
        ic = i1 + (i2 - i1) / 2
        if (xx >= xi(ic)) then
          i1 = ic
        else
          i2 = ic
        end if
      end do
    end if

  end function iix


  !> Just like iix(), but assumes that x >= xi(i_min).
  function iixmin(xx, xi, i_min) result(ip)

    !> value to search mesh index for
    real(dp), intent(in) :: xx

    !> mesh to search
    real(dp), intent(in) :: xi(:)

    !> index of assuption x >= xi(i_min)
    integer, intent(in) :: i_min

    !> i of interval [xi(i), xi(i+1)] containing x
    integer :: ip

    if ((i_min >= 1) .and. (i_min <= size(xi) - 1)) then
      ip = iix(xx, xi(i_min:)) + i_min - 1
    else
      ip = iix(xx, xi)
    end if

  end function iixmin


  !> Returns index i of interval [x(i), x(i+1)] containing x in uniform mesh defined by
  !! x(i) = x1 + (i-1)/(n-1)*(xn-x1), i = 1 .. n, with intervals indexed by left-most points.
  !! Notex: x outside [x1, xn] are indexed to nearest end.
  pure function iixun(xx, nn, x1, xn)

    !> target value
    real(dp), intent(in):: xx

    !> number of mesh points
    integer, intent(in):: nn

    !> initial point of mesh
    real(dp), intent(in):: x1

    !> final point of mesh
    real(dp), intent(in):: xn

    !> index i of interval [x(i), x(i+1)] containing x
    integer iixun

    !! index i, but may be outside 1..n
    integer ii

    ! compute index
    ii = int((xx - x1) / (xn - x1) * (nn - 1)) + 1

    ! reset if ouside 1..n
    if (ii < 1) ii = 1
    if (ii > nn - 1) ii = nn - 1
    iixun = ii

  end function iixun


  !> Returns index i of interval [x(i), x(i+1)] containing x in exponential mesh defined by
  !! x(i) = x1 + alpha [exp(beta(i-1)) - 1], i = 1 .. n,
  !! where alpha = (x(n) - x(1))/[exp(beta(n-1)) - 1], beta = log(r)/(n-2),
  !! r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
  !! and intervals indexed by left-most points.
  !! Note: x outside [x1, xn] are indexed to nearest end.
  pure function iixexp(xx, nn, x1, alpha, beta)

    !> target value
    real(dp), intent(in):: xx

    !> number of mesh points
    integer, intent(in):: nn

    !> initial point of mesh
    real(dp), intent(in):: x1

    !> mesh parameter
    real(dp), intent(in):: alpha

    !> mesh parameter
    real(dp), intent(in):: beta

    !> index i of interval [x(i), x(i+1)] containing x
    integer iixexp

    !! index i, but may be outside 1..n
    integer ii

    ! compute index
    ii = int(log((xx - x1) / alpha + 1) / beta) + 1

    ! reset if outside 1..n
    if (ii < 1) ii = 1
    if (ii > nn - 1) ii = nn - 1
    iixexp = ii

  end function iixexp


  !> Returns value of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1).
  pure function poly3(xx, cc)

    !> point at which to evaluate polynomial
    real(dp), intent(in):: xx

    !> coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
    real(dp), intent(in):: cc(0:)

    !> value of polynomial
    real(dp) poly3

    !! x - c(0)
    real(dp) dx

    dx = xx - cc(0)
    poly3 = cc(1) + cc(2) * dx + cc(3) * dx**2 + cc(4) * dx**3

  end function poly3


  !> Returns 1st derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1).
  pure function dpoly3(xx, cc)

    !> point at which to evaluate polynomial
    real(dp), intent(in):: xx

    !> coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
    real(dp), intent(in):: cc(0:)

    !> 1st derivative of polynomial
    real(dp) dpoly3

    !! x - c(0)
    real(dp) dx

    dx = xx - cc(0)
    dpoly3 = cc(2) + 2.0_dp * cc(3) * dx + 3.0_dp * cc(4) * dx**2

  end function dpoly3


  !> Returns 2nd derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1).
  pure function d2poly3(xx, cc)

    !> point at which to evaluate polynomial
    real(dp), intent(in):: xx

    !> coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
    real(dp), intent(in):: cc(0:)

    !> 2nd derivative of polynomial
    real(dp) d2poly3

    !! x - c(0)
    real(dp) dx

    dx = xx - cc(0)
    d2poly3 = 2.0_dp * cc(3) + 6.0_dp * cc(4) * dx

  end function d2poly3

end module common_splines
