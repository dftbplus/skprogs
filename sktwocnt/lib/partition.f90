!> Conains space partioning functions.
module partition
  use accuracy
  implicit none
  private

  public :: partition_becke, partition_becke_hetero, beckepar


contains

  !> Becke partition function for 2 centers.
  !! \param r1 Distance from 1st center.
  !! \param r2 Distance from 2nd center.
  !! \param dist Distance between centers.
  !! \param partparams Arbitrary dummy real array.
  !! \return Value of the partition function (between [0,1])
  !! \sa A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
  function partition_becke(r1, r2, dist, partparams) result(res)
    real(dp), intent(in) :: r1, r2, dist, partparams(:)
    real(dp) :: res

    integer :: ii

    res = (r1 - r2) / abs(dist)
    do ii = 1, 3
      res = 1.5_dp * res - 0.5 * res**3
    end do
    res = 0.5_dp * (1.0_dp - res)

  end function partition_becke


  !> Becke partition function for 2 heteronuclear centers.
  !! \param r1 Distance from 1st center.
  !! \param r2 Distance from 2nd center.
  !! \param dist Distance between centers.
  !! \param partparams Real array containing the parameter aij in the
  !!   Becke partitioning scheme.
  !! \return Value of the partition function (between [0,1])
  !! \sa A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
  function partition_becke_hetero(r1, r2, dist, partparams) result(res)
    real(dp), intent(in) :: r1, r2, dist, partparams(:)
    real(dp) :: res

    integer :: ii
    real(dp) :: mu

    mu = (r1 - r2) / abs(dist)
    res = mu + partparams(1) * (1.0_dp - mu**2)
    do ii = 1, 3
      res = 1.5_dp * res - 0.5 * res**3
    end do
    res = 0.5_dp * (1.0_dp - res)

  end function partition_becke_hetero


  !> Delivers parameter aij in the becke partition scheme for given atomic
  !! radii.
  !! \param r1 Radius of the first atom.
  !! \param r2 Radius of the second atom.
  !! \return Value of aij.
  function beckepar(r1, r2) result(res)
    real(dp), intent(in) :: r1, r2
    real(dp) :: res

    real(dp) :: chi, uu

    chi = sqrt(r1 / r2)
    uu = (chi - 1.0_dp) / (chi + 1.0_dp)
    res = uu / (uu**2 - 1.0_dp)
    if (abs(res) > 0.5_dp) then
      res = sign(0.5_dp, res)
    end if

  end function beckepar
    
end module partition
