!> Module that provides routines for numerical differentiation.
module numerical_differentiation

  use common_accuracy, only : dp
  use utilities, only : fak
  use integration, only : reverse_abcissas_1st

  implicit none
  private

  public :: numerical_1st_derivative, six_point


contains

  !> Calculates numerical 1st derivative of Beckes Gauss-Chebyschev mesh values.
  pure subroutine numerical_1st_derivative(num_mesh_points, abcissa, nuc, step, input, output)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> differentiation stepwidth
    real(dp), intent(in) :: step

    !> input on a grid
    real(dp), intent(in) :: input(:)

    !> numerical 1st derivative
    real(dp), intent(out) :: output(:)

    !> contains the six function values for 6-point differentiation
    real(dp) :: stencil(6)

    !> auxiliary variable
    integer :: ii

    output(:) = 0.0_dp
    stencil(:) = 0.0_dp

    ! handle lower mesh bound

    do ii = 1, 6
      stencil(ii) = input(ii)
    end do

    output(1) = six_point(stencil, 1, 0, step)
    output(2) = six_point(stencil, 1, 1, step)

    ! handle upper mesh bound

    do ii = 1, 6
      stencil(ii) = input(num_mesh_points - 6 + ii)
    end do

    output(num_mesh_points - 2) = six_point(stencil, 1, 3, step)
    output(num_mesh_points - 1) = six_point(stencil, 1, 4, step)
    output(num_mesh_points) = six_point(stencil, 1, 5, step)

    ! handle rest of mesh

    do ii = 3, num_mesh_points - 3

      stencil(1) = input(ii - 2)
      stencil(2) = input(ii - 1)
      stencil(3) = input(ii)
      stencil(4) = input(ii + 1)
      stencil(5) = input(ii + 2)
      stencil(6) = input(ii + 3)

      output(ii) = six_point(stencil, 1, 2, step)

    end do

    ! now remember: df(x)/dx=df(z)/dz*dz/dx, e.g. x is the abcissa which is
    ! not equally spaced and z is the generating variable of the Becke mesh
    ! which is equally spaced; so multiply by dz/dx

    do ii = 1, num_mesh_points
      output(ii) = output(ii) * reverse_abcissas_1st(nuc, abcissa(ii))
    end do

  end subroutine numerical_1st_derivative


  !> Numerical k-th derivative of tabulated function from six point
  !! formula; Bickley, Math. Gaz. vol. 25 (1941) 19-27
  !!          Abramowitz, Stegun, Handbook of Mathematical functions
  !! The function is assumed to be tabulated on equally spaced abcissas.
  !!
  !! INPUT: points      contains the six function values
  !!        kk          order of derivative, 0<k<=3
  !!        offset      0<=offset<=5
  !!                    determines on which point the point f(x) is where
  !!                    the derivative is taken
  !!                    offset=0: points(1)=f(x)
  !!                              points(2)=f(x+h)
  !!                              points(3)=f(x+2h)
  !!                              points(4)=f(x+3h)
  !!                              points(5)=f(x+4h)
  !!                              points(6)=f(x+5h)
  !!                    offset=1: points(1)=f(x-h)
  !!                              points(2)=f(x)
  !!                              points(3)=f(x+h)
  !!                              points(4)=f(x+2h)
  !!                              points(5)=f(x+3h)
  !!                              points(6)=f(x+4h)
  !!                    offset=2: points(1)=f(x-2h)
  !!                              points(2)=f(x-h)
  !!                              points(3)=f(x)
  !!                              points(4)=f(x+h)
  !!                              points(5)=f(x+2h)
  !!                              points(6)=f(x+3h)
  !!                    offset=3: points(1)=f(x-3h)
  !!                              points(2)=f(x-2h)
  !!                              points(3)=f(x-h)
  !!                              points(4)=f(x)
  !!                              points(5)=f(x+h)
  !!                              points(6)=f(x+2h)
  !!                    offset=4: points(1)=f(x-4h)
  !!                              points(2)=f(x-3h)
  !!                              points(3)=f(x-2h)
  !!                              points(4)=f(x-h)
  !!                              points(5)=f(x)
  !!                              points(6)=f(x+h)
  !!                    offset=5: points(1)=f(x-5h)
  !!                              points(2)=f(x-4h)
  !!                              points(3)=f(x-3h)
  !!                              points(4)=f(x-2h)
  !!                              points(5)=f(x-h)
  !!                              points(6)=f(x)
  !!        hh          stepwidth
  !! OUTPUT: six_point  numerical k-th derivative of tabulated function from six point formula
  pure function six_point(points, kk, offset, hh)

    !> contains the six function values for 6-point differentiation
    real(dp), intent(in) :: points(6)

    !> order of derivative, 0<k<=3
    integer, intent(in) :: kk

    !> determines on which point the point f(x) is where the derivative is taken, 0<=offset<=5
    integer, intent(in) :: offset

    !> stepwidth
    real(dp), intent(in) :: hh

    !> coefficients of differentiation formula
    !! first 36 are first derivative, second 36 are second derivative and
    !! third 36 are third derivative
    real(dp), parameter :: coeff(108) = [&
        & -274.0_dp,  600.0_dp, -600.0_dp,  400.0_dp, -150.0_dp,  24.0_dp,&
        &  -24.0_dp, -130.0_dp,  240.0_dp, -120.0_dp,   40.0_dp,  -6.0_dp,&
        &    6.0_dp,  -60.0_dp,  -40.0_dp,  120.0_dp,  -30.0_dp,   4.0_dp,&
        &   -4.0_dp,   30.0_dp, -120.0_dp,   40.0_dp,   60.0_dp,  -6.0_dp,&
        &    6.0_dp,  -40.0_dp,  120.0_dp, -240.0_dp,  130.0_dp,  24.0_dp,&
        &  -24.0_dp,  150.0_dp, -400.0_dp,  600.0_dp, -600.0_dp, 274.0_dp,&
        &  225.0_dp, -770.0_dp, 1070.0_dp, -780.0_dp,  305.0_dp, -50.0_dp,&
        &   50.0_dp,  -75.0_dp,  -20.0_dp,   70.0_dp,  -30.0_dp,   5.0_dp,&
        &   -5.0_dp,   80.0_dp, -150.0_dp,   80.0_dp,   -5.0_dp,   0.0_dp,&
        &    0.0_dp,   -5.0_dp,   80.0_dp, -150.0_dp,   80.0_dp,  -5.0_dp,&
        &    5.0_dp,  -30.0_dp,   70.0_dp,  -20.0_dp,  -75.0_dp,  50.0_dp,&
        &  -50.0_dp,  305.0_dp, -780.0_dp, 1070.0_dp, -770.0_dp, 225.0_dp,&
        &  -85.0_dp,  355.0_dp, -590.0_dp,  490.0_dp, -205.0_dp,  35.0_dp,&
        &  -35.0_dp,  125.0_dp, -170.0_dp,  110.0_dp,  -35.0_dp,   5.0_dp,&
        &   -5.0_dp,   -5.0_dp,   50.0_dp,  -70.0_dp,   35.0_dp,  -5.0_dp,&
        &    5.0_dp,  -35.0_dp,   70.0_dp,  -50.0_dp,    5.0_dp,   5.0_dp,&
        &   -5.0_dp,   35.0_dp, -110.0_dp,  170.0_dp, -125.0_dp,  35.0_dp,&
        &  -35.0_dp,  205.0_dp, -490.0_dp,  590.0_dp, -355.0_dp,  85.0_dp]

    !> pre-factor
    real(dp) :: pp

    !> numerical k-th derivative of tabulated function from six point formula
    real(dp) :: six_point

    ! get prefactor of sum formula k!/(m!*h^k), m=5 for six-point formula
    pp = fak(kk) / (real(120, dp) * (hh**kk))

    six_point = pp * (coeff((kk - 1) * 36 + offset * 6 + 1) * points(1)&
        &           + coeff((kk - 1) * 36 + offset * 6 + 2) * points(2)&
        &           + coeff((kk - 1) * 36 + offset * 6 + 3) * points(3)&
        &           + coeff((kk - 1) * 36 + offset * 6 + 4) * points(4)&
        &           + coeff((kk - 1) * 36 + offset * 6 + 5) * points(5)&
        &           + coeff((kk - 1) * 36 + offset * 6 + 6) * points(6))

  end function six_point

end module numerical_differentiation
