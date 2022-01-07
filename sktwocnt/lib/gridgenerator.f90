!> Module that provides routines for quadrature grid generation.
module gridgenerator

  use common_accuracy, only : dp
  use quadratures, only : TQuadrature
  use coordtrans, only : coordtransFunc
  use partition, only : partitionFunc

  implicit none
  private

  public :: gengrid1_12, gengrid2_12

contains


  !> Generates a 1D (radial) grid around two centers.
  pure subroutine gengrid1_12(quads, coordtrans, grid, weights)

    !> abscissas and weights for numerical quadrature
    type(TQuadrature), intent(in) :: quads(2)

    !> coordinate transformation procedure
    procedure(coordtransFunc) :: coordtrans

    !> two-dimensional atom grid, whereas r = grid(:, 1) and theta = grid(:, 2)
    real(dp), intent(out), allocatable :: grid(:,:)

    !> integration weights
    real(dp), intent(out), allocatable :: weights(:)

    !! atomic and total number of quadrature abscissas
    integer :: n1, n2, nn

    !! auxiliary variables
    integer :: ind, i1, i2
    real(dp) :: coord(2), coordreal(2), jacobi

    n1 = size(quads(1)%xx)
    n2 = size(quads(2)%xx)

    nn = n1 * n2

    allocate(grid(nn, 2))
    allocate(weights(nn))

    ind = 1
    do i2 = 1, n2
      coord(2) = quads(2)%xx(i2)
      do i1 = 1, n1
        coord(1) = quads(1)%xx(i1)
        call coordtrans(coord, coordreal, jacobi)
        grid(ind, 1) = coordreal(1)
        grid(ind, 2) = coordreal(2)
        weights(ind) = quads(1)%ww(i1) * quads(2)%ww(i2) * jacobi
        ind = ind + 1
      end do
    end do

  end subroutine gengrid1_12


  !> Generates a 2D (radial and azimuthal) grid around two centers.
  pure subroutine gengrid2_12(quads, coordtrans, partition, partparams, dist, grid1, grid2, dots,&
      & weights)

    !> abscissas and weights for numerical quadrature
    type(TQuadrature), intent(in) :: quads(2)

    !> coordinate transformation procedure
    procedure(coordtransFunc) :: coordtrans

    !> partitioning procedure
    procedure(partitionFunc) :: partition

    !> arbitrary dummy real array, unused in this routine
    real(dp), intent(in) :: partparams(:)

    !> distance between centers
    real(dp), intent(in) :: dist

    !> two-dimensional atom grids, whereas r = grid(:, 1) and theta = grid(:, 2)
    real(dp), intent(out), allocatable :: grid1(:,:), grid2(:,:)

    !> ???
    real(dp), intent(out), allocatable :: dots(:)

    !> integration weights
    real(dp), intent(out), allocatable :: weights(:)

    !! atomic and total number of quadrature abscissas
    integer :: n1, n2, nn

    !! auxiliary variables
    integer :: ind, i1, i2
    real(dp) :: coord(2), coordreal(2)
    real(dp) :: r1, theta1, r2a, r2b, theta2a, theta2b, rtmpa, rtmpb, jacobi

    n1 = size(quads(1)%xx)
    n2 = size(quads(2)%xx)

    nn = n1 * n2

    allocate(grid1(2 * nn, 2))
    allocate(grid2(2 * nn, 2))
    allocate(dots(2 * nn))
    allocate(weights(2 * nn))

    ind = 1
    do i2 = 1, n2
      coord(2) = quads(2)%xx(i2)
      do i1 = 1, n1
        coord(1) = quads(1)%xx(i1)
        call coordtrans(coord, coordreal, jacobi)
        r1 = coordreal(1)
        theta1 = coordreal(2)

        rtmpa = dist**2 + r1**2
        rtmpb = 2.0_dp * r1 * dist * cos(theta1)

        r2a = sqrt(rtmpa - rtmpb) ! dist > 0
        r2b = sqrt(rtmpa + rtmpb) ! dist < 0

        rtmpa = -0.5_dp * (dist**2 + r2a**2 - r1**2) / (dist * r2a)
        rtmpb = 0.5_dp * (dist**2 + r2b**2 - r1**2) / (dist * r2b)

        ! make sure, we are not sliding out from [-1,1] range for acos
        rtmpa = min(rtmpa, 1.0_dp)
        rtmpa = max(rtmpa, -1.0_dp)
        rtmpb = min(rtmpb, 1.0_dp)
        rtmpb = max(rtmpb, -1.0_dp)

        theta2a = acos(rtmpa)
        theta2b = acos(rtmpb)

        grid1(ind, 1) = r1
        grid1(ind, 2) = theta1
        grid1(ind + nn, 1) = r2b
        grid1(ind + nn, 2) = theta2b

        grid2(ind, 1) = r2a
        grid2(ind, 2) = theta2a
        grid2(ind + nn, 1) = r1
        grid2(ind + nn, 2) = theta1

        dots(ind) = cos(theta1 - theta2a)
        dots(ind + nn) = cos(theta2b - theta1)

        rtmpa = quads(1)%ww(i1) * quads(2)%ww(i2) * jacobi
        weights(ind) = rtmpa * partition(r1, r2a, dist, partparams)
        weights(ind + nn) = rtmpa * partition(r1, r2b, -dist, partparams)
        ind = ind + 1
      end do
    end do

  end subroutine gengrid2_12

end module gridgenerator
