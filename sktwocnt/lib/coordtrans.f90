module coordtrans

  use common_accuracy, only: dp
  use common_constants

  implicit none

contains

  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates, using the Becke algorithm.
  !! \param crd11 3d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates.
  !! \param jacobi Jacobi determinant.
  !! \sa Becke paper.
  subroutine coordtrans_becke(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: rm = 1.5_dp; 
    real(dp) :: rtmp1, rtmp2

    !assert(size(c11) == 3)
    !assert(size(spheric) == 3)

    rtmp1 = 1.0_dp + c11(1)
    rtmp2 = 1.0_dp - c11(1)
    spheric(1) = rm * (rtmp1 / rtmp2)
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)
    jacobi = 2.0_dp * rm**3 * rtmp1**2 / rtmp2**4 * pi

  end subroutine coordtrans_becke

  !> Transforms a 2 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates (r, theta), using the Becke algorithm.
  !! \param crd11 2d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates (r, theta)
  !! \param jacobi Jacobi determinant.
  !! \sa Becke paper.
  subroutine coordtrans_becke_12(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: rm = 1.5_dp; 
    real(dp) :: rtmp1, rtmp2

    !assert(size(c11) == 2)
    !assert(size(spheric) == 2)

    rtmp1 = 1.0_dp + c11(1)
    rtmp2 = 1.0_dp - c11(1)
    spheric(1) = rm * (rtmp1 / rtmp2)
    spheric(2) = acos(c11(2))
    jacobi = 2.0_dp * rm**3 * rtmp1**2 / rtmp2**4

  end subroutine coordtrans_becke_12

  !> Transforms a 2 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates (theta, phi), using the Becke algorithm.
  !! \param crd11 2d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates (theta, phi).
  !! \param jacobi Jacobi determinant.
  !! \sa Becke paper.
  subroutine coordtrans_becke_23(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    !assert(size(c11) == 2)
    !assert(size(spheric) == 2)

    spheric(1) = acos(c11(1))
    spheric(2) = pi * (c11(2) + 1.0_dp)
    jacobi = pi

  end subroutine coordtrans_becke_23

  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates, using the Ahlrichs algorithm.
  !! \param c11 3d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates.
  !! \param jacobi Jacobi determinant.
  !! \sa Ahlrichs paper.
  subroutine coordtrans_ahlrichs1(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.20_dp
    real(dp) :: rr

    !assert(size(c11) == 3)
    !assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)
    jacobi = (zeta / log(2.0_dp)) / (1.0_dp - c11(1)) * rr * rr * pi

  end subroutine coordtrans_ahlrichs1

  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates, using the Ahlrichs algorithm.
  !! \param c11 3d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates.
  !! \param jacobi Jacobi determinant.
  !! \sa Ahlrichs paper.
  subroutine coordtrans_ahlrichs1_2d(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.20_dp
    real(dp) :: rr

    !assert(size(c11) == 3)
    !assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    !spheric(3) = pi * (c11(3) + 1.0_dp)
    jacobi = (zeta / log(2.0_dp)) / (1.0_dp - c11(1)) * rr * rr

  end subroutine coordtrans_ahlrichs1_2d

  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates, using the Ahlrichs algorithm.
  !! \param c11 3d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates.
  !! \param jacobi Jacobi determinant.
  !! \sa Ahlrichs paper.
  subroutine coordtrans_ahlrichs2(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.1_dp
    real(dp), parameter :: alpha = 0.6_dp
    real(dp) :: rr

    !assert(size(c11) == 3)
    !assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * (1.0_dp + c11(1))**alpha &
        &* log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)

    jacobi = (zeta * (1.0_dp + c11(1))**alpha / log(2.0_dp)) &
        &* (alpha * log(2.0_dp / (1.0_dp - c11(1))) / (1.0_dp + c11(1)) &
        &+ 1.0_dp / (1.0_dp - c11(1))) * rr * rr * pi

  end subroutine coordtrans_ahlrichs2

  !> Transforms a 3 dimensional vector with coordinates in [-1,1] onto spherical
  !! coordinates, using the Ahlrichs algorithm.
  !! \param c11 3d coordinate vector, each coordinate in [-1,1].
  !! \param spheric Corresponding spherical coordinates.
  !! \param jacobi Jacobi determinant.
  !! \sa Ahlrichs paper.
  subroutine coordtrans_ahlrichs2_2d(c11, spheric, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: spheric(:)
    real(dp), intent(out) :: jacobi

    real(dp), parameter :: zeta = 1.1_dp
    real(dp), parameter :: alpha = 0.6_dp
    real(dp) :: rr

    !assert(size(c11) == 3)
    !assert(size(spheric) == 3)

    rr = (zeta / log(2.0_dp)) * (1.0_dp + c11(1))**alpha &
        &* log(2.0_dp / (1.0_dp - c11(1)))
    spheric(1) = rr
    spheric(2) = acos(c11(2))
    spheric(3) = pi * (c11(3) + 1.0_dp)

    jacobi = (zeta * (1.0_dp + c11(1))**alpha / log(2.0_dp)) &
        &* (alpha * log(2.0_dp / (1.0_dp - c11(1))) / (1.0_dp + c11(1)) &
        &+ 1.0_dp / (1.0_dp - c11(1))) * rr * rr

  end subroutine coordtrans_ahlrichs2_2d

  subroutine coordtrans_identity(c11, ctarget, jacobi)
    real(dp), intent(in) :: c11(:)
    real(dp), intent(out) :: ctarget(:)
    real(dp), intent(out) :: jacobi

    ctarget = c11
    jacobi = 1.0_dp

  end subroutine coordtrans_identity

end module coordtrans
