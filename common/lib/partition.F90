!> Module that provides (Becke's) space partitioning functions.
module common_partition

  use common_accuracy, only : dp

  implicit none
  private

  public :: partitionFunc
  public :: partition_becke_homo, partition_becke_hetero, beckepar


  abstract interface

    !> General interface of (Becke's) partitioning functions.
    pure function partitionFunc(r1, r2, dist, partparams) result(res)

      use common_accuracy, only : dp

      implicit none

      !> distance from 1st center
      real(dp), intent(in) :: r1

      !> distance from 2nd center
      real(dp), intent(in) :: r2

      !> distance between centers
      real(dp), intent(in) :: dist

      !> holds partitioning parameters, if required
      real(dp), intent(in) :: partparams(:)

      !! resulting value of the partition function, between [0,1]
      real(dp) :: res

    end function partitionFunc

  end interface


contains

  !> Becke partition function for 2 homonuclear centers, Voronoi polyhedra bisect internuclear axes,
  !! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
  pure function partition_becke_homo(r1, r2, dist, partparams) result(res)

    !> distance from 1st center
    real(dp), intent(in) :: r1

    !> distance from 2nd center
    real(dp), intent(in) :: r2

    !> distance between centers
    real(dp), intent(in) :: dist

    !> arbitrary dummy real array, unused in this routine
    real(dp), intent(in) :: partparams(:)

    !! resulting value of the partition function, between [0,1]
    real(dp) :: res

    !! auxiliary variable
    integer :: ii

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 11
    res = (r1 - r2) / abs(dist)

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 19/20, choosing k=3
    do ii = 1, 3
      res = 1.5_dp * res - 0.5 * res**3
    end do

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 21
    res = 0.5_dp * (1.0_dp - res)

  end function partition_becke_homo


  !> Becke partition function for 2 heteronuclear centers, cell boundaries shifted away from
  !! internuclear midpoints, see A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
  pure function partition_becke_hetero(r1, r2, dist, partparams) result(res)

    !> distance from 1st center
    real(dp), intent(in) :: r1

    !> distance from 2nd center
    real(dp), intent(in) :: r2

    !> distance between centers
    real(dp), intent(in) :: dist

    !> real array containing the parameter aij in the Becke partitioning scheme
    real(dp), intent(in) :: partparams(:)

    !! resulting value of the partition function, between [0,1]
    real(dp) :: res

    !! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 11
    real(dp) :: mu

    !! auxiliary variable
    integer :: ii

    ! assert(abs(partparams(1)) <= 0.5_dp)

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 11
    mu = (r1 - r2) / abs(dist)

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. A2
    res = mu + partparams(1) * (1.0_dp - mu**2)

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 19/20, choosing k=3
    do ii = 1, 3
      res = 1.5_dp * res - 0.5 * res**3
    end do

    ! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. 21
    res = 0.5_dp * (1.0_dp - res)

  end function partition_becke_hetero


  !> Delivers parameter aij in the becke partition scheme for given atomic radii,
  !! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
  pure function beckepar(r1, r2) result(res)

    !> Bragg-Slater radius of first and second atom
    real(dp), intent(in) :: r1, r2

    !! parameter a_{ij}, see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. A5
    real(dp) :: res

    !! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. A4
    real(dp) :: chi

    !! see A. D. Becke, J. Chem. Phys. 88, 2547 (1988), eqn. A6
    real(dp) :: uu

    ! chi = r1 / r2
    chi = sqrt(r1 / r2)

    uu = (chi - 1.0_dp) / (chi + 1.0_dp)

    res = uu / (uu**2 - 1.0_dp)

    if (abs(res) > 0.5_dp) then
      res = sign(0.5_dp, res)
    end if

  end function beckepar

end module common_partition
