!> Spherical harmonics.
module sphericalharmonics
  use accuracy
  implicit none
  private

  public :: realtess, init, destruct, getvalue, getvalue_1d

  !> Real tessereal shperical.
  type realtess
    private
    integer :: ll, mm
  end type realtess

  interface init
    module procedure realtess_init
  end interface

  interface destruct
    module procedure realtess_destruct
  end interface

  interface getvalue
    module procedure realtess_getvalue
  end interface

  interface getvalue_1d
    module procedure realtess_getvalue_1d
  end interface

contains

  !> Initialises realtess.
  !! \param self instance.
  !! \param ll angulam momentum (l)
  !! \param mm magnetic quantum number (m)
  subroutine realtess_init(self, ll, mm)
    type(realtess), intent(inout) :: self
    integer, intent(in) :: ll, mm

    self%ll = ll
    self%mm = mm
    
  end subroutine realtess_init


  !> Destroys the instance.
  !! \param self instance.
  subroutine realtess_destruct(self)
    type(realtess), intent(inout) :: self

    continue
    
  end subroutine realtess_destruct

   
  !> returns the value of the tessereal function.
  !! \param self instance.
  !! \param theta spherical coordinate theta.
  !! \param phi spherical coordinate phi.
  elemental function realtess_getvalue(self, theta, phi) result(ang)
    type(realtess), intent(in) :: self
    real(dp), intent(in) :: theta, phi
    real(dp) :: ang

    ang = calc_realtess(self%ll, self%mm, theta, phi)
    
  end function realtess_getvalue


  elemental function realtess_getvalue_1d(self, theta) result(ang)
    type(realtess), intent(in) :: self
    real(dp), intent(in) :: theta
    real(dp) :: ang

    ang = calc_realtess_1d(self%ll, self%mm, theta)
    
  end function realtess_getvalue_1d
  
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! private functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Real tessereal spherical harmonics up to f.
  !! \param ll angular momentum (l).
  !! \param mm magnetic moment (m)
  !! \param theta spherical coordinate theta.
  !! \param phi spherical coordinate phi.
  !! \return value of the real tesseral harmonics.
  elemental function calc_realtess(ll, mm, theta, phi) result (rty)
    integer, intent(in) :: ll
    integer, intent(in) :: mm
    real(dp), intent(in) :: theta, phi
    real(dp) :: rty
    
    !assert(ll >= 0 .and. ll <= 3)
    !assert(abs(mm) <= ll)
    
    select case (ll)
    case(0)
      rty = 0.2820947917738782_dp
    case(1)
      select case(mm)
      case(-1)
        rty = 0.4886025119029198_dp * sin(theta) * sin(phi)
      case(0)
        rty = 0.4886025119029198_dp * cos(theta)
      case(1)
        rty = 0.4886025119029198_dp * sin(theta) * cos(phi)
      end select
    case(2)
      select case(mm)
      case(-2)
        rty = 0.5462742152960395_dp * sin(theta)**2 * sin(2.0_dp * phi)
      case(-1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta) * sin(phi)
      case(0)
        rty = 0.9461746957575600_dp * cos(theta)**2 - 0.3153915652525200_dp
      case(1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta) * cos(phi)
      case(2)
        rty = 0.5462742152960395_dp * sin(theta)**2 * cos(2.0_dp * phi)
      end select
    case(3)
      select case (mm)
      case(-3)
        rty = 0.5900435899266435_dp * sin(theta)**3 * sin(3.0_dp * phi)
      case(-2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta) &
            &* sin(2.0_dp * phi)
      case(-1)
        rty = 0.4570457994644658_dp * sin(theta) &
            &* (5.0_dp * cos(theta)**2 - 1.0_dp) * sin(phi)
      case(0)
        rty = 0.3731763325901155_dp * cos(theta) &
            &* (5.0_dp * cos(theta)**2 - 3.0_dp)
      case(1)
        rty = 0.4570457994644658_dp * sin(theta) &
            &* (5.0_dp * cos(theta)**2 - 1.0_dp) * cos(phi)
      case(2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta) &
            &* cos(2.0_dp * phi)
      case(3)
        rty = 0.5900435899266435_dp * sin(theta)**3 * cos(3.0_dp * phi)
      end select
    end select
    
  end function calc_realtess


  !> Real tessereal spherical harmonics up to f.
  !! \param ll angular momentum (l).
  !! \param mm magnetic moment (m)
  !! \param theta spherical coordinate theta.
  !! \param phi spherical coordinate phi.
  !! \return value of the real tesseral harmonics.
  elemental function calc_realtess_1d(ll, mm, theta) result (rty)
    integer, intent(in) :: ll
    integer, intent(in) :: mm
    real(dp), intent(in) :: theta
    real(dp) :: rty
    
    !assert(ll >= 0 .and. ll <= 3)
    !assert(abs(mm) <= ll)
    
    select case (ll)
    case(0)
      rty = 0.2820947917738782_dp
    case(1)
      select case(mm)
      case(-1)
        rty = 0.4886025119029198_dp * sin(theta)
      case(0)
        rty = 0.4886025119029198_dp * cos(theta)
      case(1)
        rty = 0.4886025119029198_dp * sin(theta)
      end select
    case(2)
      select case(mm)
      case(-2)
        rty = 0.5462742152960395_dp * sin(theta)**2
      case(-1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta)
      case(0)
        rty = 0.9461746957575600_dp * cos(theta)**2 - 0.3153915652525200_dp
      case(1)
        rty = 1.092548430592079_dp * sin(theta) * cos(theta)
      case(2)
        rty = 0.5462742152960395_dp * sin(theta)**2
      end select
    case(3)
      select case (mm)
      case(-3)
        rty = 0.5900435899266435_dp * sin(theta)**3
      case(-2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta)
      case(-1)
        rty = 0.4570457994644658_dp * sin(theta) &
            &* (5.0_dp * cos(theta)**2 - 1.0_dp)
      case(0)
        rty = 0.3731763325901155_dp * cos(theta) &
            &* (5.0_dp * cos(theta)**2 - 3.0_dp)
      case(1)
        rty = 0.4570457994644658_dp * sin(theta) &
            &* (5.0_dp * cos(theta)**2 - 1.0_dp)
      case(2)
        rty = 1.445305721320277_dp * sin(theta)**2 * cos(theta)
      case(3)
        rty = 0.5900435899266435_dp * sin(theta)**3
      end select
    end select
    
  end function calc_realtess_1d

end module sphericalharmonics
