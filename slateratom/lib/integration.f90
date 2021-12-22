!> Module that provides numerical integration routines.
module integration

  use common_accuracy, only : dp
  use common_constants, only : pi
  use utilities, only : fak

  implicit none
  private

  public :: gauss_chebyshev_becke_mesh
  public :: get_abcissas, get_abcissas_z_1st, get_abcissas_z_2nd
  public :: reverse_abcissas, reverse_abcissas_1st, reverse_abcissas_2nd
  public :: exp_int


contains

  !> Generate Beckes Gauss-Chebyschev mesh, e.g. radial points and weights.
  !! Becke mapping: r \in [-1, 1] --> [0, inf)
  pure subroutine gauss_chebyshev_becke_mesh(num_mesh_points, nuc, weight, abcissa, dzdr, d2zdr2,&
      & dz)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> numerical integration weight factors of mesh
    real(dp), intent(out) :: weight(:)

    !> numerical integration abcissas (radii), Becke mapping
    real(dp), intent(out) :: abcissa(:)

    !> dz/dr
    real(dp), intent(out) :: dzdr(:)

    !> d2z/dr2
    real(dp), intent(out) :: d2zdr2(:)

    !> step width in linear coordinates
    real(dp), intent(out) :: dz

    !> determinental factor of mapping
    real(dp), allocatable :: fak(:)

    !> gauss-chebyshev abcissas
    real(dp), allocatable :: xx(:)

    !> auxiliary variables
    integer :: ii
    real(dp) :: temp, zz, cosz, cosz2, sinz

    allocate(xx(num_mesh_points))
    allocate(fak(num_mesh_points))

    temp = pi / real(num_mesh_points + 1, dp)
    dz = temp

    do ii = 1, num_mesh_points
      zz = dz * real(ii, dp)
      cosz = cos(zz)
      cosz2 = cosz * cosz
      sinz = sqrt(1.0_dp - cosz2)
      ! NOTE prefactor
      xx(ii) = (-1.0_dp) * cosz ! gauss-chebyshev abcissas
      abcissa(ii) = (1.0_dp + xx(ii)) / (1.0_dp - xx(ii)) * bragg(nuc)
      !dzdr(ii) = (1.0_dp + 2.0_dp * cos(zz) + cos(zz)**2) / (2.0_dp * bragg(nuc) * sin(zz))
      dzdr(ii) = (1.0_dp + cosz)**2 / (2.0_dp * bragg(nuc) * sinz)
      d2zdr2(ii) = ((2.0_dp + cosz - cosz2) * (1.0_dp + cosz)**2) / (4.0_dp * bragg(nuc)**2&
          & * (-1.0_dp + cosz) * sinz)

      ! r**2 times first derivative of x -> r mapping function
      weight(ii) = temp * (sin(real(ii, dp) * temp))
      ! fak(ii)=2.0_dp*abcissa(ii)**2*bragg(nuc)/(1.0_dp-xx(ii))**2
      fak(ii) = 2.0_dp * bragg(nuc) / (1.0_dp - xx(ii))**2

      ! put fak into weight
      weight(ii) = weight(ii) * fak(ii)
    end do

  end subroutine gauss_chebyshev_becke_mesh


  pure subroutine get_abcissas(num_mesh_points, nuc, rr, step)
    ! r(x)=bragg*(1-x)/(1+x)
    ! x(z)=cos(pi*z)
    ! r(x(z))=bragg*(1-cos(pi*z))/(1+cos(pi*z)), z=ii/(num_mesh_points+1)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> numerical integration abcissas (radii), Becke mapping
    real(dp), intent(out) :: rr(:)

    !> generator step size
    integer, intent(out) :: step

    !> gauss-chebyshev abcissas
    real(dp), allocatable :: xx(:)

    !> auxiliary variable
    integer :: ii

    allocate(xx(num_mesh_points))

    step = int(pi / real(num_mesh_points + 1, dp))

    do ii = 1, num_mesh_points
      ! NOTE prefactor
      xx(ii) = (-1.0_dp) * cos(step * real(ii, dp))
      rr(ii) = (1.0_dp + xx(ii)) / (1.0_dp - xx(ii)) * bragg(nuc)
    end do

  end subroutine get_abcissas


  !> 1st derivative of r(x(z)) with respect to z, see grid_differentiation_sign_2.txt
  pure subroutine get_abcissas_z_1st(num_mesh_points, nuc, dr, step)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> 1st dderiv. of abcissas, Becke mapping
    real(dp), intent(out) :: dr(:)

    !> generator step size
    integer, intent(out) :: step

    !> auxiliary variable
    integer :: ii

    step = int(pi / real(num_mesh_points + 1, dp))

    do ii = 1, num_mesh_points
      dr(ii) = 2.0_dp * bragg(nuc) * pi * sin(step * real(ii, dp))&
          & / (1.0_dp + 2.0_dp * cos(step * real(ii, dp)) + cos(step * real(ii, dp))**2)
    end do

  end subroutine get_abcissas_z_1st


  !> 2nd derivative of r(x) with respect to x, see grid_differentiation_sign_2.txt
  pure subroutine get_abcissas_z_2nd(num_mesh_points, nuc, ddr, step)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> 2nd deriv. of abcissas, Becke mapping
    real(dp), intent(out) :: ddr(:)

    !> generator step size
    integer, intent(out) :: step

    !> auxiliary variable
    integer :: ii

    step = int(pi / real(num_mesh_points + 1, dp))

    do ii = 1, num_mesh_points
      ddr(ii) = (- 2.0_dp * bragg(nuc) * pi**2) * (cos(step * real(ii, dp)) - 2.0_dp)&
          & / (1.0_dp + 2.0_dp * cos(step * real(ii, dp)) + cos(step * real(ii, dp))**2)
    end do

  end subroutine get_abcissas_z_2nd


  !> z(x(r)) reverse mapping function, see grid_differentiation_sign_2.txt
  !!
  !! z=1/pi*arccos((a-r)/(a+r))
  pure function reverse_abcissas(nuc, rr)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> numerical integration abcissa (radius), Becke mapping
    real(dp), intent(in) :: rr

    !> reverse mapping
    real(dp) :: reverse_abcissas

    reverse_abcissas = 1.0_dp / pi * acos((bragg(nuc) - rr) / (bragg(nuc) + rr))

  end function reverse_abcissas


  !> 1st derivative of z(x(r)) reverse mapping function w.r.t. r,
  !! see grid_differentiation_sign_2.txt
  !!
  !! be careful: can easily overflow
  pure function reverse_abcissas_1st(nuc, rr)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> numerical integration abcissa (radius), Becke mapping
    real(dp), intent(in) :: rr

    !> 1st derivative of reverse mapping
    real(dp) :: reverse_abcissas_1st

    reverse_abcissas_1st = 1.0_dp / pi * sqrt(bragg(nuc) / rr) / (rr + bragg(nuc))

  end function reverse_abcissas_1st


  !> 2nd derivative of z(x(r)) reverse mapping function w.r.t. r, see
  !! grid_differentiation_sign_2.txt
  !!
  !! be careful: can easily overflow
  pure function reverse_abcissas_2nd(nuc, rr)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> numerical integration abcissa (radius), Becke mapping
    real(dp), intent(in) :: rr

    !> 2nd derivative of reverse mapping
    real(dp) :: reverse_abcissas_2nd

    reverse_abcissas_2nd = - 1.0_dp / (2.0_dp * pi) * sqrt(bragg(nuc) / rr) / rr&
        & * (bragg(nuc) + 3.0_dp * rr) / (bragg(nuc) + rr)**2

  end function reverse_abcissas_2nd


  pure function bragg(nuc)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> prefactor corresponding to the requested element
    real(dp) :: bragg

    !> element-resolved Becke prefactors
    real(dp), parameter :: braggd(110) = [&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp,&
        & 1.0_dp, 1.0_dp]

    bragg = braggd(nuc)

  end function bragg


  !> evaluate \int x**power*exp(alpha*x) dx at point r, for formula see Bronstein
  !! (assumes alpha<0 and power=>0, WATCH OUT FOR SIGN OF alpha)
  function exp_int(alpha, power, rr)

    !> exponential coefficient
    real(dp), intent(in) :: alpha

    !> power coefficient
    integer, intent(in) :: power

    !> evaluation point
    real(dp), intent(in) :: rr

    !> integral value
    real(dp) :: exp_int

    !> auxiliary variable
    integer :: ii

    exp_int = 0.0_dp

    ! catch power<0
    if (power < 0) then
      write(*,*) 'NEGATIVE POWERS NOT IMPLEMENTED !'
      stop
    end if

    ! catch alpha>0
    if (alpha > 0.0_dp) then
      write(*,*) 'POSITIVE ALPHAS NOT IMPLEMENTED !'
      stop
    end if

    ! catch r=0
    if (rr == 0.0_dp) then
      exp_int = fak(power) / (alpha**(power + 1)) * (- 1.0_dp)**power
      return
    end if

    ! catch r=infty and alpha<0 (should always be !)
    if (abs(alpha * rr) > 75.0_dp) then
      exp_int = 0.0_dp
      return
    end if

    exp_int = 1.0_dp / alpha * exp(alpha * rr)

    do ii = 1, power
      exp_int = 1.0_dp / alpha * rr**ii * exp(alpha * rr) - real(ii, dp) / alpha * exp_int
    end do

  end function exp_int

end module integration
