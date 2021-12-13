module quadratures

  use common_accuracy, only: dp
  use common_constants

  implicit none

  type quadrature
    real(dp), allocatable :: xx(:)
    real(dp), allocatable :: ww(:)
  end type quadrature

  real(dp), parameter :: eps = 1e-14_dp

contains

  !> Gauss-Legendre quadrature for integration in the interval [-1,1].
  !! \param nn Number of points for the quadrature
  !! \param quad  Quadrature with abscissas and weights.
  !! \sa Numerical Recipes
  subroutine gauss_legendre_quadrature(nn, quad)
    integer, intent(in) :: nn
    type(quadrature), intent(out) :: quad

    integer :: mm, ii, jj
    real(dp) :: zz, z1, pp, p1, p2, p3, rj

    allocate(quad % xx(nn))
    allocate(quad % ww(nn))
    mm = (nn + 1) / 2
    do ii = 1, mm
      zz = cos(pi * (real(ii, dp) - 0.25_dp) / (real(nn, dp) + 0.5_dp))
      do
        p1 = 1.0_dp
        p2 = 0.0_dp
        do jj = 1, nn
          p3 = p2
          p2 = p1
          rj = real(jj, dp)
          p1 = ((2.0_dp * rj - 1.0_dp) * zz * p2 - (rj - 1.0_dp) * p3) / rj
        end do
        pp = real(nn, dp) * (zz * p1 - p2) / (zz * zz - 1.0_dp)
        z1 = zz
        zz = z1 - (p1 / pp)
        if (abs(zz - z1) <= eps) then
          exit
        end if
      end do
      quad % xx(ii) = -zz
      quad % xx(nn + 1 - ii) = zz
      quad % ww(ii) = 2.0_dp / ((1.0_dp - zz * zz) * pp * pp)
      quad % ww(nn + 1 - ii) = quad % ww(ii)
    end do

  end subroutine gauss_legendre_quadrature

  !> Gauss-Chebishev quadrature for integration in the interval [-1,1].
  !!
  !! Integration of functions with Gauss-Chebishev quadrature of second kind.
  !! The weights already contain 1/sqrt(1-x^2) so that it can be directly
  !! used to integrate a function on [-1,1].
  !! See also: J. M. Pérez-Jordá et al., J. Chem. Phys. 100 6520 (1994).
  !!
  !! \param nn Number of points for the quadrature
  !! \param quad  Quadrature with abscissas and weights.
  subroutine gauss_chebyshev_quadrature(nn, quad)
    integer, intent(in) :: nn
    type(quadrature), intent(out) :: quad

    integer :: ii
    real(dp) :: rtmp

    allocate(quad % xx(nn))
    allocate(quad % ww(nn))
    !do ii = 1, nn
    !  quad%xx(ii) = cos(pi * (real(ii, dp) - 0.5_dp) / real(nn, dp))
    !end do
    !quad%ww = pi / real(nn, dp)
    do ii = 1, nn
      rtmp = real(ii, dp) * pi / real(nn + 1, dp)
      quad % xx(ii) = cos(rtmp)
      quad % ww(ii) = sin(rtmp)
    end do
    quad % ww = quad % ww * pi / real(nn + 1, dp)

  end subroutine gauss_chebyshev_quadrature

  !> Trapezoidal quadrature for integration in the interval [-1,1].
  !! \param nn Number of points for the quadrature
  !! \param quad  Quadrature with abscissas and weights.
  !! \sa Numerical Recipes
  subroutine trapezoidal_quadrature(nn, quad)
    integer, intent(in) :: nn
    type(quadrature), intent(out) :: quad

    integer :: ii
    real(dp) :: fac

    allocate(quad % xx(nn))
    allocate(quad % ww(nn))
    fac = 2.0_dp / real(nn, dp)
    do ii = 1, nn
      quad % xx(ii) = -1.0_dp + fac * real(ii - 1, dp)
    end do
    quad % ww = fac

  end subroutine trapezoidal_quadrature

end module quadratures
