!> Module that provides routines for quadrature grid generation.
module common_gridgenerator

  use common_accuracy, only : dp
  use common_constants, only : pi
  use common_quadratures, only : TQuadrature, TQuadrature2D
  use common_coordtrans, only : coordtransFunc_0d, coordtransFunc_1d
  use common_partition, only : partitionFunc

  implicit none
  private

  public :: gengrid1_1, gengrid1_2, gengrid1_3, gengrid2_2, gengrid2_3


contains

  !> Generates a one-center grid, for radial integration only.
  subroutine gengrid1_1(radQuad, rm, coordtrans, grid, weights)

    interface
      !> General interface for utilized coordinate transformation.
      subroutine coordtrans(oldc, rm, newc, jacobi)
        use common_accuracy, only : dp

        !> old coordinate
        real(dp), intent(in) :: oldc

        !> midpoint of the integration interval,
        !! allows adjustment of the radial point distribution to a suitable physical scale
        real(dp), intent(in) :: rm

        !> new coordinate after transformation
        real(dp), intent(out) :: newc

        !> Jacobi determinant
        real(dp), intent(out) :: jacobi
      end subroutine coordtrans
    end interface

    !> holds abscissas and weights for radial quadrature
    type(TQuadrature), intent(in) :: radQuad

    !> midpoint of the integration interval,
    !! allows adjustment of the radial point distribution to a suitable physical scale
    real(dp), intent(in) :: rm

    !> integration grid around center
    real(dp), allocatable, intent(out) :: grid(:,:)

    !> integration weights
    real(dp), allocatable, intent(out) :: weights(:)

    !! current abscissa and transformed equivalent
    real(dp) :: coord, coordreal

    !! current Jacobi determinant
    real(dp) :: jacobi

    !! number of radial quadrature points
    integer :: nn

    !! iterates over all grid points
    integer :: i1

    nn = size(radQuad%xx)

    allocate(grid(nn, 1))
    allocate(weights(nn))

    do i1 = 1, nn
      coord = radQuad%xx(i1)
      ! transform the computational coordinate to the real one
      call coordtrans(coord, rm, coordreal, jacobi)
      ! grid is a set of real coordinates
      grid(i1, 1) = coordreal
      weights(i1) = radQuad%ww(i1) * jacobi
    end do

  end subroutine gengrid1_1


  !> Generates a one-center grid for two variables.
  subroutine gengrid1_2(quads, coordtrans, grid, weights)

    !> abscissas and weights for numerical quadrature
    type(TQuadrature), intent(in) :: quads(2)

    !> coordinate transformation procedure
    procedure(coordtransFunc_1d) :: coordtrans

    !> integration grid around center
    real(dp), allocatable, intent(out) :: grid(:,:)

    !> integration weights
    real(dp), allocatable, intent(out) :: weights(:)

    !! number of radial and angular quadrature points
    integer :: n1, n2

    !! total number of quadrature points
    integer :: nn

    !! (real) current radial and angular coordinate
    real(dp) :: coord(2), coordreal(2)

    !! Jacobi determinant
    real(dp) :: jacobi

    !! loop counter
    integer :: ind

    !! iterates over radial and angular quadrature points
    integer :: i1, i2

    n1 = size(quads(1)%xx)
    n2 = size(quads(2)%xx)

    nn = n1 * n2

    allocate(grid(nn, 2))
    allocate(weights(nn))

    ind = 1
    do i2 = 1, n2
      coord(2) = quads(2)%xx(i2)
      do i1 = 1, n1
        ! coord contains computational coordinate in [-1,1]
        coord(1) = quads(1)%xx(i1)
        ! transform the computational coordinate to the real one
        call coordtrans(coord, coordreal, jacobi)
        ! grid is a set of real coordiantes
        grid(ind, 1) = coordreal(1)
        grid(ind, 2) = coordreal(2)
        weights(ind) = quads(1)%ww(i1) * quads(2)%ww(i2) * jacobi
        ind = ind + 1
      end do
    end do

  end subroutine gengrid1_2


  !> Generates a one-center grid for three variables.
  subroutine gengrid1_3(angQuad, radQuad, coordtrans, grid, weights)

    !> abscissas and weights for angular quadrature with coordinate pairs
    type(TQuadrature2D), intent(in) :: angQuad

    !> abscissas and weights for radial quadrature
    type(TQuadrature), intent(in) :: radQuad

    !> coordinate transformation procedure
    procedure(coordtransFunc_0d) :: coordtrans

    !> integration grid around center
    real(dp), allocatable, intent(out) :: grid(:,:)

    !> integration weights
    real(dp), allocatable, intent(out) :: weights(:)

    !! number of radial and angular quadrature points
    integer :: n1, n2

    !! total number of quadrature points
    integer :: nn

    !! (real) current radial and angular coordinate
    real(dp) :: coord, coordreal

    !! Jacobi determinant
    real(dp) :: jacobi

    !! loop counter
    integer :: ind

    !! iterates over radial and angular quadrature points
    integer :: i1, i2

    n1 = size(radQuad%xx)
    n2 = size(angQuad%c1)

    nn = n1 * n2

    allocate(grid(nn, 3))
    allocate(weights(nn))

    ind = 1
    do i1 = 1, n1
      coord = radQuad%xx(i1)
      ! transform the computational coordinate to the real one
      call coordtrans(coord, coordreal, jacobi)
      do i2 = 1, n2
        ! grid is a set of real coordiantes
        grid(ind, 1) = coordreal
        grid(ind, 2) = angQuad%c1(i2)
        grid(ind, 3) = angQuad%c2(i2)
        weights(ind) = radQuad%ww(i1) * angQuad%ww(i2) * jacobi
        ind = ind + 1
      end do
    end do

  end subroutine gengrid1_3


  !> Generates a two-center grid for two spherical variables (r, theta) each.
  pure subroutine gengrid2_2(quads, coordtrans, partition, partparams, dist, grid1, grid2, dots,&
      & weights)

    !> abscissas and weights for numerical quadrature
    type(TQuadrature), intent(in) :: quads(2)

    !> coordinate transformation procedure
    procedure(coordtransFunc_1d) :: coordtrans

    !> partitioning procedure
    procedure(partitionFunc) :: partition

    !> arbitrary dummy real array, unused in this routine
    real(dp), intent(in) :: partparams(:)

    !> distance between centers
    real(dp), intent(in) :: dist

    !> two-dimensional atom grids, whereas r = grid(:, 1) and theta = grid(:, 2)
    real(dp), intent(out), allocatable :: grid1(:,:), grid2(:,:)

    !> scalar products of unit vectors from each perspective
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

        rtmpa = - 0.5_dp * (dist**2 + r2a**2 - r1**2) / (dist * r2a)
        rtmpb = 0.5_dp * (dist**2 + r2b**2 - r1**2) / (dist * r2b)

        ! make sure, we are not sliding out from [-1,1] range for acos
        rtmpa = min(rtmpa, 1.0_dp)
        rtmpa = max(rtmpa, - 1.0_dp)
        rtmpb = min(rtmpb, 1.0_dp)
        rtmpb = max(rtmpb, - 1.0_dp)

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
        weights(ind + nn) = rtmpa * partition(r1, r2b, - dist, partparams)
        ind = ind + 1
      end do
    end do

  end subroutine gengrid2_2


  !> Generates two radial grids at (0,0,0) and (0,0,dist).
  !! Each grid consists of points given by quadratures. Also the mapping between grids is computed,
  !! i.e. for each point in coordiante system 1 (reference frame for grid) the corresponding point
  !! in coordinate system 2 is calculated and vice versa.
  !!
  !! \note Both grids consist of the same set of quadrature points in their local coordinate system.
  !!
  !! The structure of grid object:
  !!  grid1      grid2    weight
  !! [Quad]  -> [Link]   [Quad1]
  !! [Link]  <- [Quad]   [Quad2]
  !!
  !! Quad: quadrature points
  !! Link: corresponding coordinates in other system
  subroutine gengrid2_3(angQuad, radQuad, coordtrans, partition, partparams,&
      & dist, grid1, grid2, weights, part)
    type(TQuadrature2D), intent(in) :: angQuad
    type(TQuadrature), intent(in) :: radQuad

    !> coordinate transformation procedure
    procedure(coordtransFunc_0d) :: coordtrans

    !> partitioning procedure
    procedure(partitionFunc) :: partition

    real(dp), intent(in) :: partparams(:)
    real(dp), intent(in) :: dist
    real(dp), allocatable, intent(out) :: grid1(:,:), grid2(:,:)
    real(dp), allocatable, intent(out) :: weights(:), part(:)

    !! number of radial and angular quadrature points
    integer :: n1, n2

    !! total number of quadrature points
    integer :: nn

    !! (real) current radial coordinate
    real(dp) :: coord, coordreal

    !! Jacobi determinant
    real(dp) :: jacobi

    !! loop counter
    integer :: ind

    !! iterates over radial and angular quadrature points
    integer :: i1, i2

    !! auxiliary variables
    real(dp) :: r1, theta1, r2a, r2b, theta2a, theta2b, rtmpa, rtmpb
    real(dp) :: phi, rm, z2a, z2b

    rm = 1.0_dp

    n1 = size(radQuad%xx)
    n2 = size(angQuad%c1)

    nn = n1 * n2

    allocate(grid1(2 * nn, 4))
    allocate(grid2(2 * nn, 4))
    allocate(weights(2 * nn))
    allocate(part(2 * nn))

    ind = 1
    do i1 = 1, n1
      ! the radial coordinate should be transformed
      coord = radQuad%xx(i1)
      call coordtrans(coord, coordreal, jacobi)
      r1 = coordreal
      do i2 = 1, n2
        ! Lebedev quadrature returns spherical coordinates
        theta1 = angQuad%c1(i2)
        phi = angQuad%c2(i2)

        ! now calculate the link to the other system, use spherical prolate coordinates
        rtmpa = dist * dist + r1 * r1
        rtmpb = 2.0_dp * r1 * dist * cos(theta1)
        r2a = sqrt(rtmpa - rtmpb) ! dist > 0
        r2b = sqrt(rtmpa + rtmpb) ! dist < 0
        rtmpa = -0.5_dp * (dist * dist + r2a * r2a - r1 * r1) / (dist * r2a)
        rtmpb = 0.5_dp * (dist * dist + r2b * r2b - r1 * r1) / (dist * r2b)

        ! here the z coordinate is calculated additionally
        z2a = acos((r2a - rm)/(r2a + rm))/pi
        z2b = acos((r2b - rm)/(r2b + rm))/pi

        !! make sure, we are not sliding out from [-1,1] range for acos
        rtmpa = min(rtmpa, 1.0_dp)
        rtmpa = max(rtmpa, - 1.0_dp)
        rtmpb = min(rtmpb, 1.0_dp)
        rtmpb = max(rtmpb, - 1.0_dp)

        theta2a = acos(rtmpa)
        theta2b = acos(rtmpb)

        ! fill the grid with values, according to description above
        ! note that phi is equal for both systems
        grid1(ind, 1) = r1
        grid1(ind, 2) = theta1
        grid1(ind, 3) = phi

        ! here the z coordinate is calculated additionally
        grid1(ind, 4) = radQuad%zz(i1)

        grid1(ind + nn, 1) = r2b
        grid1(ind + nn, 2) = theta2b
        grid1(ind + nn, 3) = phi

        ! here the z coordinate is calculated additionally
        grid1(ind + nn, 4) = z2b

        grid2(ind, 1) = r2a
        grid2(ind, 2) = theta2a
        grid2(ind, 3) = phi

        ! here the z coordinate is calculated additionally
        grid2(ind, 4) = z2a

        grid2(ind + nn, 1) = r1
        grid2(ind + nn, 2) = theta1
        grid2(ind + nn, 3) = phi

        ! here the z coordinate is calculated additionally
        grid2(ind + nn, 4) = radQuad%zz(i1)

        ! generate weights, including space partition for full 3D integration
        rtmpa = radQuad%ww(i1) * angQuad%ww(i2) * jacobi

        ! additionaly generate the partition function
        part(ind) = partition(r1, r2a, dist, partparams)
        part(ind + nn) = partition(r1, r2b, - dist, partparams)
        weights(ind) =  rtmpa * part(ind)
        weights(ind + nn) = rtmpa * part(ind + nn)
        ind = ind + 1
      end do
    end do

  end subroutine gengrid2_3

end module common_gridgenerator
