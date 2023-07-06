!******************************************************************************************
! the module depends now on Lebedev-Laikov.F
!******************************************************************************************

module common_quadratures

  use common_accuracy, only : dp
  use common_constants, only : pi
  use common_lebedev_laikov, only : ld0006, ld0014, ld0026, ld0038, ld0050, ld0074, ld0086,&
      & ld0110, ld0146, ld0170, ld0194, ld0230, ld0266, ld0302, ld0350

  implicit none
  private

  public :: TQuadrature, TQuadrature2D
  public :: lebedev_laikov_quadrature, gauss_legendre_quadrature, gauss_chebyshev_quadrature,&
      & trapezoidal_quadrature


  !> Holds abscissas and weights for numerical quadrature.
  type TQuadrature

    !> abscissas
    real(dp), allocatable :: xx(:)

    !> weights
    real(dp), allocatable :: ww(:)

    !> scaled weight/abscissa index i / (n + 1)
    real(dp), allocatable :: zz(:)

  end type TQuadrature


  !> Quadrature type for nodes given by coordinate pairs.
  type TQuadrature2D

    !> first coordinate vector
    real(dp), allocatable :: c1(:)

    !> second coordinate vector
    real(dp), allocatable :: c2(:)

    !> weights
    real(dp), allocatable :: ww(:)

  end type TQuadrature2D


  !> relative quadrature precision
  real(dp), parameter :: eps = 1e-14_dp


contains

  !> Gauss type quadrature, which integrates spherical harmonics exactly.
  !! Use for effective integration on the unit sphere.
  !!
  !! Includes the code from Lebedev-Laikov.F
  !! The subroutines give nodes in carthesian coordinates -> we need polar.
  !! => ToDo: do hardcoded version, which returns polar coordinates
  !!
  !! references: see Lebedev-Laikov.F
  subroutine lebedev_laikov_quadrature(nn, quad2d)

    !> number of points for the quadrature
    integer, intent(in) :: nn

    !> at exit, holds quadrature type for nodes given by coordinate pairs
    type(TQuadrature2D), intent(out) :: quad2d

    !> carthesian coordinates of nodes
    real(dp), allocatable :: xx(:), yy(:), zz(:)

    !! actual number of points for the quadrature
    integer :: nPoints

    !! auxiliary variable
    integer :: ii

    select case(nn)
    case(1:7)
       nPoints = 6
    case(8:16)
       nPoints = 14
    case(17:28)
       nPoints = 26
    case(29:40)
       nPoints = 38
    case(41:60)
       nPoints = 50
    case(61:80)
       nPoints = 74
    case(81:100)
       nPoints = 86
    case(101:130)
       nPoints = 110
    case(131:160)
       nPoints = 146
    case(161:180)
       nPoints = 170
    case(181:210)
       nPoints = 194
    case(211:250)
       nPoints = 230
    case(251:290)
       nPoints = 266
    case(291:320)
       nPoints = 302
    case(321:400)
       nPoints = 350
    case default
       nPoints = 6
    end select

    allocate(quad2d%c1(nPoints))
    allocate(quad2d%c2(nPoints))
    allocate(quad2d%ww(nPoints))

    allocate(xx(nPoints))
    allocate(yy(nPoints))
    allocate(zz(nPoints))

    select case(nPoints)
    case(6)
       call ld0006(xx, yy, zz, quad2d%ww, nPoints)
    case(14)
       call ld0014(xx, yy, zz, quad2d%ww, nPoints)
    case(26)
       call ld0026(xx, yy, zz, quad2d%ww, nPoints)
    case(38)
       call ld0038(xx, yy, zz, quad2d%ww, nPoints)
    case(50)
       call ld0050(xx, yy, zz, quad2d%ww, nPoints)
    case(74)
       call ld0074(xx, yy, zz, quad2d%ww, nPoints)
    case(86)
       call ld0086(xx, yy, zz, quad2d%ww, nPoints)
    case(110)
       call ld0110(xx, yy, zz, quad2d%ww, nPoints)
    case(146)
       call ld0146(xx, yy, zz, quad2d%ww, nPoints)
    case(170)
       call ld0170(xx, yy, zz, quad2d%ww, nPoints)
    case(194)
       call ld0194(xx, yy, zz, quad2d%ww, nPoints)
    case(230)
       call ld0230(xx, yy, zz, quad2d%ww, nPoints)
    case(266)
       call ld0266(xx, yy, zz, quad2d%ww, nPoints)
    case(302)
       call ld0302(xx, yy, zz, quad2d%ww, nPoints)
    case(350)
       call ld0350(xx, yy, zz, quad2d%ww, nPoints)
    case default
       call ld0006(xx, yy, zz, quad2d%ww, nPoints)
    end select

    ! necessary, since Lebedev-Quadrature is defined for integral
    ! I = \frac{1}{4\pi}\int F(\Omega) d\Omega
    ! see references in lebedev_laikov_f77.f
    quad2d%ww(:) = quad2d%ww * 4.0_dp * pi

    ! transform the nodes in carthesian coordinates to spherical coordinates
    do ii = 1, nPoints
       quad2d%c1(ii) = acos(zz(ii))

       ! compute angle phi in interval [0, 2*pi)
       if (xx(ii) > 0 .and. yy(ii) >= 0) quad2d%c2(ii) = atan(yy(ii) / xx(ii))
       if (xx(ii) > 0 .and. yy(ii) < 0) quad2d%c2(ii) = atan(yy(ii) / xx(ii)) + 2.0_dp * pi
       if (xx(ii) < 0) quad2d%c2(ii) = atan(yy(ii) / xx(ii)) + pi
       if (xx(ii) == 0 .and. yy(ii) > 0) quad2d%c2(ii) = 1.5707963267948966_dp
       if (xx(ii) == 0 .and. yy(ii) < 0) quad2d%c2(ii) = 3.0_dp * 1.5707963267948966_dp
    end do

  end subroutine lebedev_laikov_quadrature


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
    real(dp) :: fac

    allocate(quad%xx(nn))
    allocate(quad%ww(nn))
    allocate(quad%zz(nn))

    ! see J. M. Pérez-Jordá et al., J. Chem. Phys. 100 6520 (1994), eqn. 28/29
    fac = real(nn + 1, dp)
    do ii = 1, nn
       rtmp = real(ii, dp) / fac
       quad%zz(ii) = rtmp
       rtmp = rtmp * pi
       quad%xx(ii) = cos(rtmp)
       quad%ww(ii) = sin(rtmp)
    end do
    fac = pi / real(nn + 1, dp)
    quad%ww(:) = quad%ww * fac

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

end module common_quadratures
