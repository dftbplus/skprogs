!> Module that provides the functionality for real tesseral spherical harmonics.
module sphericalharmonics

  use common_accuracy, only : dp

  implicit none
  private

  public :: TRealTessY, TRealTessY_init


  !> Real tesseral spherical harmonics.
  type TRealTessY

    !> angular momentum
    integer :: ll

    !> magnetic quantum number
    integer :: mm

  contains

    procedure :: getValue => TRealTessY_getValue
    procedure :: getValue_1d => TRealTessY_getValue_1d
    procedure :: destruct => TRealTessY_destruct

  end type TRealTessY


contains

  !> Initialises a TRealTessY object.
  subroutine TRealTessY_init(this, ll, mm)

    !> real tesseral spherical harmonics instance
    type(TRealTessY), intent(out) :: this

    !> angular momentum (l)
    integer, intent(in) :: ll

    !> magnetic quantum number (m)
    integer, intent(in) :: mm

    this%ll = ll
    this%mm = mm

  end subroutine TRealTessY_init


  !> Destroys an initialised instance.
  subroutine TRealTessY_destruct(this)

    !> real tesseral spherical harmonics instance
    class(TRealTessY), intent(inout) :: this

    continue

  end subroutine TRealTessY_destruct


  !> Returns value of real tesseral spherical harmonic function.
  elemental function TRealTessY_getValue(this, theta, phi) result(ang)

    !> real tesseral spherical harmonics instance
    class(TRealTessY), intent(in) :: this

    !> spherical coordinate theta
    real(dp), intent(in) :: theta

    !> spherical coordinate phi
    real(dp), intent(in) :: phi

    !> value of real tesseral spherical harmonic function
    real(dp) :: ang

    ang = calc_realtessy(this%ll, this%mm, theta, phi)

  end function TRealTessY_getValue


  !> Returns value of real tesseral spherical harmonic function.
  elemental function TRealTessY_getValue_1d(this, theta) result(ang)

    !> real tesseral spherical harmonics instance
    class(TRealTessY), intent(in) :: this

    !> spherical coordinate theta
    real(dp), intent(in) :: theta

    !> value of real tesseral spherical harmonic function
    real(dp) :: ang

    ang = calc_realtessy_1d(this%ll, this%mm, theta)

  end function TRealTessY_getValue_1d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! private functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Real tesseral spherical harmonics up to angular momentum f.
  elemental function calc_realtessy(ll, mm, theta, phi) result(rty)

    !> angular momentum (l)
    integer, intent(in) :: ll

    !> magnetic quantum number (m)
    integer, intent(in) :: mm

    !> spherical coordinate theta
    real(dp), intent(in) :: theta

    !> spherical coordinate phi
    real(dp), intent(in) :: phi

    !> value of real tesseral spherical harmonic function
    real(dp) :: rty

    ! assert(ll >= 0 .and. ll <= 3)
    ! assert(abs(mm) <= ll)

    select case (ll)
    case (0)
      rty = 0.2820947917738782_dp
    case (1)
      select case (mm)
      case (-1)
        rty = 0.4886025119029198_dp * sin(theta) * sin(phi)
      case (0)
        rty = 0.4886025119029198_dp * cos(theta)
      case (1)
        rty = 0.4886025119029198_dp * sin(theta) * cos(phi)
      end select
    case (2)
      select case (mm)
      case (-2)
        rty = 0.5462742152960395_dp * sin(theta)**2 * sin(2.0_dp * phi)
      case (-1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta) * sin(phi)
      case (0)
        rty = 0.9461746957575600_dp * cos(theta)**2 - 0.3153915652525200_dp
      case (1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta) * cos(phi)
      case (2)
        rty = 0.5462742152960395_dp * sin(theta)**2 * cos(2.0_dp * phi)
      end select
    case (3)
      select case (mm)
      case (-3)
        rty = 0.5900435899266435_dp * sin(theta)**3 * sin(3.0_dp * phi)
      case (-2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta) &
             &* sin(2.0_dp * phi)
      case (-1)
        rty = 0.4570457994644658_dp * sin(theta) &
             &* (5.0_dp * cos(theta)**2 - 1.0_dp) * sin(phi)
      case (0)
        rty = 0.3731763325901155_dp * cos(theta) &
             &* (5.0_dp * cos(theta)**2 - 3.0_dp)
      case (1)
        rty = 0.4570457994644658_dp * sin(theta) &
             &* (5.0_dp * cos(theta)**2 - 1.0_dp) * cos(phi)
      case (2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta) &
             &* cos(2.0_dp * phi)
      case (3)
        rty = 0.5900435899266435_dp * sin(theta)**3 * cos(3.0_dp * phi)
      end select
    end select

  end function calc_realtessy


  !> Real tesseral spherical harmonics up to angular momentum f.
  elemental function calc_realtessy_1d(ll, mm, theta) result(rty)

    !> angular momentum (l)
    integer, intent(in) :: ll

    !> magnetic quantum number (m)
    integer, intent(in) :: mm

    !> spherical coordinate theta
    real(dp), intent(in) :: theta

    !> value of real tesseral spherical harmonic function
    real(dp) :: rty

    ! assert(ll >= 0 .and. ll <= 3)
    ! assert(abs(mm) <= ll)

    select case (ll)
    case (0)
      rty = 0.2820947917738782_dp
    case (1)
      select case (mm)
      case (-1)
        rty = 0.4886025119029198_dp * sin(theta)
      case (0)
        rty = 0.4886025119029198_dp * cos(theta)
      case (1)
        rty = 0.4886025119029198_dp * sin(theta)
      end select
    case (2)
      select case (mm)
      case (-2)
        rty = 0.5462742152960395_dp * sin(theta)**2
      case (-1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta)
      case (0)
        rty = 0.9461746957575600_dp * cos(theta)**2 - 0.3153915652525200_dp
      case (1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta)
      case (2)
        rty = 0.5462742152960395_dp * sin(theta)**2
      end select
    case (3)
      select case (mm)
      case (-3)
        rty = 0.5900435899266435_dp * sin(theta)**3
      case (-2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta)
      case (-1)
        rty = 0.4570457994644658_dp * sin(theta) &
             &* (5.0_dp * cos(theta)**2 - 1.0_dp)
      case (0)
        rty = 0.3731763325901155_dp * cos(theta) &
             &* (5.0_dp * cos(theta)**2 - 3.0_dp)
      case (1)
        rty = 0.4570457994644658_dp * sin(theta) &
             &* (5.0_dp * cos(theta)**2 - 1.0_dp)
      case (2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta)
      case (3)
        rty = 0.5900435899266435_dp * sin(theta)**3
      end select
    end select

  end function calc_realtessy_1d

end module sphericalharmonics
