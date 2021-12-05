module numerical_differentiation

  use common_accuracy, only : dp
  use common_constants
  use utilities
  use integration

  implicit none
  private

  public :: numerical_1st_derivative, six_point
  
  
contains

  subroutine numerical_1st_derivative(num_mesh_points,abcissa,nuc,step,&
      &input,output)

    integer, intent(in) :: num_mesh_points,nuc
    real(dp), intent(in) :: input(:),step,abcissa(:)
    real(dp), intent(out) :: output(:)
    real(dp) :: stencil(6)
    integer :: ii

    output=0.0d0
    stencil=0.0d0

    ! handle lower mesh bound

    do ii=1,6
      stencil(ii)=input(ii)
    end do

    output(1)=six_point(stencil,1,0,step)
    output(2)=six_point(stencil,1,1,step)

    ! handle upper mesh bound

    do ii=1,6
      stencil(ii)=input(num_mesh_points-6+ii)
    end do

    output(num_mesh_points-2)=six_point(stencil,1,3,step)
    output(num_mesh_points-1)=six_point(stencil,1,4,step)
    output(num_mesh_points)=six_point(stencil,1,5,step)

    ! handle rest of mesh

    !$OMP PARALLEL DO PRIVATE(ii,stencil)
    do ii=3,num_mesh_points-3

      stencil(1)=input(ii-2)
      stencil(2)=input(ii-1)
      stencil(3)=input(ii)
      stencil(4)=input(ii+1)
      stencil(5)=input(ii+2)
      stencil(6)=input(ii+3)

      output(ii)=six_point(stencil,1,2,step)

    end do
    !$OMP END PARALLEL DO

    ! now remember: df(x)/dx=df(z)/dz*dz/dx, e.g. x is the abcissa which is
    ! not equally spaced and z is the generating variable of the Becke mesh
    ! which is equally spaced; so multiply by dz/dx

    !$OMP PARALLEL DO PRIVATE(ii,stencil)
    do ii=1,num_mesh_points

      output(ii)=output(ii)*reverse_abcissas_1st(nuc,abcissa(ii))

    end do
    !$OMP END PARALLEL DO

  end subroutine numerical_1st_derivative

  function six_point(points,k,offset,h)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    ! Numerical k-th derivative of tabulated function from six point
    ! formula; Bickley, Math. Gaz. vol. 25 (1941) 19-27
    !          Abramowitz, Stegun, Handbook of Mathematical functions
    ! The function is assumed to be tabulated on equally spaced abcissas
    !
    ! INPUT: points      contains the six function values
    !        k           order of derivative, 0<k<=3
    !        offset      0<=offset<=5
    !                    determines on which point the point f(x) is where 
    !                    the derivative is taken
    !                    offset=0: points(1)=f(x)
    !                              points(2)=f(x+h)
    !                              points(3)=f(x+2h)
    !                              points(4)=f(x+3h)
    !                              points(5)=f(x+4h)
    !                              points(6)=f(x+5h)
    !                    offset=1: points(1)=f(x-h)
    !                              points(2)=f(x)
    !                              points(3)=f(x+h)
    !                              points(4)=f(x+2h)
    !                              points(5)=f(x+3h)
    !                              points(6)=f(x+4h)
    !                    offset=2: points(1)=f(x-2h)
    !                              points(2)=f(x-h)
    !                              points(3)=f(x)
    !                              points(4)=f(x+h)
    !                              points(5)=f(x+2h)
    !                              points(6)=f(x+3h)
    !                    offset=3: points(1)=f(x-3h)
    !                              points(2)=f(x-2h)
    !                              points(3)=f(x-h)
    !                              points(4)=f(x)
    !                              points(5)=f(x+h)
    !                              points(6)=f(x+2h)
    !                    offset=4: points(1)=f(x-4h)
    !                              points(2)=f(x-3h)
    !                              points(3)=f(x-2h)
    !                              points(4)=f(x-h)
    !                              points(5)=f(x)
    !                              points(6)=f(x+h)
    !                    offset=5: points(1)=f(x-5h)
    !                              points(2)=f(x-4h)
    !                              points(3)=f(x-3h)
    !                              points(4)=f(x-2h)
    !                              points(5)=f(x-h)
    !                              points(6)=f(x)
    !        h           stepwidth
    ! OUTPUT:
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    real(dp) :: points(6),h,six_point,coeff(108),p
    integer :: k,offset
    ! coefficients of differentiation formula
    ! first 36 are first derivative, second 36 are second derivative and
    ! third 36 are third derivative
    DATA coeff/-274.0d0, 600.0d0,-600.0d0, 400.0d0,-150.0d0,  24.0d0,&
        &            -24.0d0,-130.0d0, 240.0d0,-120.0d0,  40.0d0,  -6.0d0,&
        &              6.0d0, -60.0d0, -40.0d0, 120.0d0, -30.0d0,   4.0d0,&
        &             -4.0d0,  30.0d0,-120.0d0,  40.0d0,  60.0d0,  -6.0d0,&
        &              6.0d0, -40.0d0, 120.0d0,-240.0d0, 130.0d0,  24.0d0,&
        &            -24.0d0, 150.0d0,-400.0d0, 600.0d0,-600.0d0, 274.0d0,&
        &            225.0d0,-770.0d0,1070.0d0,-780.0d0, 305.0d0, -50.0d0,&
        &             50.0d0, -75.0d0, -20.0d0,  70.0d0, -30.0d0,   5.0d0,&
        &             -5.0d0,  80.0d0,-150.0d0,  80.0d0,  -5.0d0,   0.0d0,&
        &              0.0d0,  -5.0d0,  80.0d0,-150.0d0,  80.0d0,  -5.0d0,&
        &              5.0d0, -30.0d0,  70.0d0, -20.0d0, -75.0d0,  50.0d0,&
        &            -50.0d0, 305.0d0,-780.0d0,1070.0d0,-770.0d0, 225.0d0,&
        &            -85.0d0, 355.0d0,-590.0d0, 490.0d0,-205.0d0,  35.0d0,&
        &            -35.0d0, 125.0d0,-170.0d0, 110.0d0, -35.0d0,   5.0d0,&
        &             -5.0d0,  -5.0d0,  50.0d0, -70.0d0,  35.0d0,  -5.0d0,&
        &              5.0d0, -35.0d0,  70.0d0, -50.0d0,   5.0d0,   5.0d0,&
        &             -5.0d0,  35.0d0,-110.0d0, 170.0d0,-125.0d0,  35.0d0,&
        &            -35.0d0, 205.0d0,-490.0d0, 590.0d0,-355.0d0,  85.0d0/

    ! get prefactor of sum formula k!/(m!*h^k), m=5 for six-point formula
    p=fak(k)/(120*(h**k))

    six_point=p*(coeff((k-1)*36+offset*6+1)*points(1)+&
        &             coeff((k-1)*36+offset*6+2)*points(2)+&
        &             coeff((k-1)*36+offset*6+3)*points(3)+&
        &             coeff((k-1)*36+offset*6+4)*points(4)+&
        &             coeff((k-1)*36+offset*6+5)*points(5)+&
        &             coeff((k-1)*36+offset*6+6)*points(6))

  end function six_point

end module numerical_differentiation
