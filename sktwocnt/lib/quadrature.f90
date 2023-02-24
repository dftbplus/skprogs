!> Module that provides several quadrature related functionalities.
module quadratures

  use common_accuracy, only : dp
  use common_constants, only : pi

  implicit none
  private

  public :: TQuadrature
  public :: gauss_legendre_quadrature, gauss_chebyshev_quadrature, trapezoidal_quadrature


  !> Holds abscissas and weights for numerical quadrature.
  type TQuadrature

    !> abscissas
    real(dp), allocatable :: xx(:)

    !> weights
    real(dp), allocatable :: ww(:)

    !>
    real(dp), allocatable :: zz(:)

  end type TQuadrature

  !> relative quadrature precision
  real(dp), parameter :: eps = 1e-14_dp


contains

  !> Gauss-Legendre quadrature for integration in the interval [-1,1],
  !! see Numerical Recipes or J. M. Pérez-Jordá et al., J. Chem. Phys. 100 6520 (1994).
  pure subroutine gauss_legendre_quadrature(nn, quad)

    !> number of points for the quadrature
    integer, intent(in) :: nn

    !> at exit, holds abscissas and weights for numerical quadrature
    type(TQuadrature), intent(out) :: quad

    !! number of roots after symmetry is considered
    integer :: mm

    !! initial approximations to the roots
    real(dp) :: zz

    !! auxiliary variables
    integer :: ii, jj
    real(dp) :: z1, pp, p1, p2, p3, rj

    allocate(quad%xx(nn))
    allocate(quad%ww(nn))

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
        if (abs(zz - z1) <= eps) exit
      end do
      quad%xx(ii) = - zz
      quad%xx(nn + 1 - ii) = zz
      quad%ww(ii) = 2.0_dp / ((1.0_dp - zz**2) * pp**2)
      quad%ww(nn + 1 - ii) = quad%ww(ii)
    end do

  end subroutine gauss_legendre_quadrature


  !> Gauss-Chebishev quadrature for integration in the interval [-1,1].
  !!
  !! Integration of functions with Gauss-Chebishev quadrature of second kind. The weights already
  !! contain 1/sqrt(1-x^2) so that it can be directly used to integrate a function on [-1,1],
  !! see J. M. Pérez-Jordá et al., J. Chem. Phys. 100 6520 (1994).
  pure subroutine gauss_chebyshev_quadrature(nn, quad)

    !> number of points for the quadrature
    integer, intent(in) :: nn

    !> at exit, holds abscissas and weights for numerical quadrature
    type(TQuadrature), intent(out) :: quad

    !! recurring argument of trigonometry functions
    real(dp) :: rtmp

    !! auxiliary variable
    integer :: ii

    allocate(quad%xx(nn))
    allocate(quad%ww(nn))

    ! see J. M. Pérez-Jordá et al., J. Chem. Phys. 100 6520 (1994), eqn. 28/29
    do ii = 1, nn
      rtmp = real(ii, dp) * pi / real(nn + 1, dp)
      quad%xx(ii) = cos(rtmp)
      quad%ww(ii) = sin(rtmp)
    end do
    quad%ww(:) = quad%ww * pi / real(nn + 1, dp)

  end subroutine gauss_chebyshev_quadrature


  !> Trapezoidal quadrature for integration in the interval [-1,1],
  !! see Numerical Recipes.
  pure subroutine trapezoidal_quadrature(nn, quad)

    !> number of points for the quadrature
    integer, intent(in) :: nn

    !> at exit, holds abscissas and weights for numerical quadrature
    type(TQuadrature), intent(out) :: quad

    !! discretization stepwidth of interval [-1,1]
    real(dp) :: fac

    !! auxiliary variable
    integer :: ii

    allocate(quad%xx(nn))
    allocate(quad%ww(nn))

    fac = 2.0_dp / real(nn, dp)
    do ii = 1, nn
      quad%xx(ii) = - 1.0_dp + fac * real(ii - 1, dp)
    end do
    quad%ww(:) = fac

  end subroutine trapezoidal_quadrature

end module quadratures
