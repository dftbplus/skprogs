!> Module that provides several routines related to coordinate transformation.
module coordtrans

  use common_accuracy, only : dp
  use common_constants, only : pi

  implicit none
  private

  public :: coordtransFunc, coordtrans_becke, coordtrans_becke_12, coordtrans_becke_23,&
      & coordtrans_ahlrichs1, coordtrans_ahlrichs1_2d, coordtrans_ahlrichs2,&
      & coordtrans_ahlrichs2_2d, coordtrans_identity


  abstract interface

    !> General interface for (Bekce's) coordinate transformations.
    pure subroutine coordtransFunc(oldc, newc, jacobi)

      use common_accuracy, only : dp

      implicit none

      !> old coordinate vector
      real(dp), intent(in) :: oldc(:)

      !> new coordinate vector after transformation
      real(dp), intent(out) :: newc(:)

      !> Jacobi determinant
      real(dp), intent(out) :: jacobi

    end subroutine coordtransFunc

  end interface


contains

  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical coordinates, using
  !! the Becke algorithm, see A. D. Becke, J. Chem. Phys. 88, 2547 (1988)
  !! or  J. Chem. Phys. 100, 6520 (1994).
  pure subroutine coordtrans_becke(c11, spheric, jacobi)

    !> 3d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    !! midpoint of the integration interval,
    !! allows adjustment of the radial point distribution to a suitable physical scale
    real(dp), parameter :: rm = 1.5_dp

    !! recurring factors
    real(dp) :: rtmp1, rtmp2

    ! assert(size(c11) == 3)
    ! assert(size(spheric) == 3)

    rtmp1 = 1.0_dp + c11(1)
    rtmp2 = 1.0_dp - c11(1)
    spheric(1) = rm * (rtmp1 / rtmp2)
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)
    jacobi = 2.0_dp * pi * rm**3 * rtmp1**2 / rtmp2**4

  end subroutine coordtrans_becke


  !> Transforms a 2 dimensional vector with coordinates in [-1,1] onto spherical coordinates
  !! (r, theta), using the Becke algorithm, see A. D. Becke, J. Chem. Phys. 88, 2547 (1988)
  !! or  J. Chem. Phys. 100, 6520 (1994).
  pure subroutine coordtrans_becke_12(c11, spheric, jacobi)

    !> 2d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates (r, theta)
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    !! midpoint of the integration interval,
    !! allows adjustment of the radial point distribution to a suitable physical scale
    real(dp), parameter :: rm = 1.5_dp

    !! recurring factors
    real(dp) :: rtmp1, rtmp2

    ! assert(size(c11) == 2)
    ! assert(size(spheric) == 2)

    rtmp1 = 1.0_dp + c11(1)
    rtmp2 = 1.0_dp - c11(1)
    spheric(1) = rm * (rtmp1 / rtmp2)
    spheric(2) = acos(c11(2))
    jacobi = 2.0_dp * rm**3 * rtmp1**2 / rtmp2**4

  end subroutine coordtrans_becke_12


  !> Transforms a 2 dimensional vector with coordinates in [-1,1] onto spherical coordinates
  !! (theta, phi), using the Becke algorithm, see A. D. Becke, J. Chem. Phys. 88, 2547 (1988)
  !! or  J. Chem. Phys. 100, 6520 (1994).
  pure subroutine coordtrans_becke_23(c11, spheric, jacobi)

    !> 2d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates (theta, phi)
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    ! assert(size(c11) == 2)
    ! assert(size(spheric) == 2)

    spheric(1) = acos(c11(1))
    spheric(2) = pi * (c11(2) + 1.0_dp)
    jacobi = pi

  end subroutine coordtrans_becke_23


  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical coordinates, using
  !! the Ahlrichs algorithm (cf. Ahlrichs paper).
  pure subroutine coordtrans_ahlrichs1(c11, spheric, jacobi)

    !> 3d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.20_dp
    real(dp) :: rr

    ! assert(size(c11) == 3)
    ! assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)
    jacobi = (zeta / log(2.0_dp)) / (1.0_dp - c11(1)) * rr * rr * pi

  end subroutine coordtrans_ahlrichs1


  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical coordinates, using
  !! the Ahlrichs algorithm (cf. Ahlrichs paper).
  pure subroutine coordtrans_ahlrichs1_2d(c11, spheric, jacobi)

    !> 3d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.20_dp
    real(dp) :: rr

    ! assert(size(c11) == 3)
    ! assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    ! spheric(3) = pi * (c11(3) + 1.0_dp)
    jacobi = (zeta / log(2.0_dp)) / (1.0_dp - c11(1)) * rr * rr

  end subroutine coordtrans_ahlrichs1_2d


  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical coordinates, using
  !! the Ahlrichs algorithm (cf. Ahlrichs paper).
  pure subroutine coordtrans_ahlrichs2(c11, spheric, jacobi)

    !> 3d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.1_dp
    real(dp), parameter :: alpha = 0.6_dp
    real(dp) :: rr

    !assert(size(c11) == 3)
    !assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * (1.0_dp + c11(1))**alpha * log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)

    jacobi = (zeta * (1.0_dp + c11(1))**alpha / log(2.0_dp))&
        & * (alpha * log(2.0_dp / (1.0_dp - c11(1))) / (1.0_dp + c11(1)) + 1.0_dp&
        & / (1.0_dp - c11(1))) * rr * rr * pi

  end subroutine coordtrans_ahlrichs2


  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical coordinates, using
  !! the Ahlrichs algorithm (cf. Ahlrichs paper).
  pure subroutine coordtrans_ahlrichs2_2d(c11, spheric, jacobi)

    !> 3d coordinate vector, each coordinate in interval [-1,1]
    real(dp), intent(in) :: c11(:)

    !> corresponding spherical coordinates
    real(dp), intent(out) :: spheric(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.1_dp
    real(dp), parameter :: alpha = 0.6_dp
    real(dp) :: rr

    ! assert(size(c11) == 3)
    ! assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * (1.0_dp + c11(1))**alpha * log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)

    jacobi = (zeta * (1.0_dp + c11(1))**alpha / log(2.0_dp))&
        & * (alpha * log(2.0_dp / (1.0_dp - c11(1))) / (1.0_dp + c11(1))&
        & + 1.0_dp / (1.0_dp - c11(1))) * rr * rr

  end subroutine coordtrans_ahlrichs2_2d


  !> Identity coordinate transformation.
  pure subroutine coordtrans_identity(c11, ctarget, jacobi)

    !> coordinate vector
    real(dp), intent(in) :: c11(:)

    !> target vector
    real(dp), intent(out) :: ctarget(:)

    !> Jacobi determinant
    real(dp), intent(out) :: jacobi

    ctarget(:) = c11
    jacobi = 1.0_dp

  end subroutine coordtrans_identity

end module coordtrans
