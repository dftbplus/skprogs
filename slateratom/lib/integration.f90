module integration
  use accuracy
  use constants
  use utilities
  implicit none
  private

  public :: gauss_chebyshev_becke_mesh
  public :: get_abcissas, get_abcissas_z_1st, get_abcissas_z_2nd
  public :: reverse_abcissas, reverse_abcissas_1st, reverse_abcissas_2nd
  public :: exp_int
 
contains
  
  subroutine gauss_chebyshev_becke_mesh(N,nuc,w,r, dzdr, d2zdr2, dz)

    ! Generate Beckes Gauss-Chebyschev mesh, e.g. radial points and weights.

    integer, intent(in) :: N ! number of mesh points
    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(out) :: w(:) ! weight factors of mesh
    real(dp), intent(out) :: r(:) ! radii of abcissas, Becke mapping !
    real(dp), intent(out) :: dzdr(:) ! dz/dr
    real(dp), intent(out) :: d2zdr2(:) ! d^2 z / dr^2
    real(dp), intent(out) :: dz

    real(dp), allocatable :: fak(:) ! determinental factor of mapping
    real(dp), allocatable :: x(:)
    real(dp) :: temp
    integer :: ii
    real(dp) :: zz, cosz, cosz2, sinz
    !
    allocate(x(N)) 
    allocate(fak(N)) 
    !
    temp=pi/float(N+1)
    dz = temp
    !    
    do ii=1,N
      zz = dz * real(ii, dp)
      cosz = cos(zz)
      cosz2 = cosz * cosz
      sinz = sqrt(1.0_dp - cosz2)
      ! NOTE prefactor 
      x(ii)=(-1.0_dp) * cosz ! gauss-chebyshev abcissas
      r(ii)= (1.0_dp + x(ii)) / (1.0_dp - x(ii)) * bragg(nuc)
      !dzdr(ii) = (1.0_dp + 2.0_dp * cos(zz) + cos(zz)**2) &
      !    &/ (2.0_dp * bragg(nuc) * sin(zz))
      dzdr(ii) = (1.0_dp + cosz)**2 / (2.0_dp * bragg(nuc) * sinz)
      d2zdr2(ii) = ((2.0_dp + cosz - cosz2) * (1.0_dp + cosz)**2) &
          &/ (4.0_dp * bragg(nuc)**2 * (-1.0_dp + cosz) * sinz)

      ! r**2 times first derivative of x -> r mapping function
      w(ii)=temp*(sin(float(ii)*temp))
      !      fak(ii)=2.0_dp*r(ii)**2*bragg(nuc)/(1.0_dp-x(ii))**2
      fak(ii)=2.0_dp*bragg(nuc)/(1.0_dp-x(ii))**2

      ! put fak into weight
      w(ii)=w(ii)*fak(ii)
    end do

    deallocate(x) 
    deallocate(fak) 

  end subroutine gauss_chebyshev_becke_mesh

  subroutine get_abcissas(N,nuc,r,step)
    ! r(x)=bragg*(1-x)/(1+x)
    ! x(z)=cos(pi*z)
    ! r(x(z))=bragg*(1-cos(pi*z))/(1+cos(pi*z)), z=ii/(N+1)

    integer, intent(in) :: N ! number of mesh points
    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(out) :: r(:) ! radii of abcissas, Becke mapping !
    integer, intent(out) :: step ! generator step size
    real(dp), allocatable :: x(:)
    integer :: ii

    allocate(x(N)) 

    step=pi/float(N+1)

    do ii=1,N

      ! NOTE prefactor 
      x(ii)=(-1.0_dp)*cos(step*float(ii)) ! gauss-chebyshev abcissas
      r(ii)=(1.0_dp+x(ii))/(1.0_dp-x(ii))*bragg(nuc)

    end do

    deallocate(x) 

  end subroutine get_abcissas

  subroutine get_abcissas_z_1st(N,nuc,dr,step)
    ! 1st derivative of r(x(z)) with respect to z, see
    ! grid_differentiation_sign_2.txt

    integer, intent(in) :: N ! number of mesh points
    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(out) :: dr(:) ! 1st dderiv. of abcissas, Becke mapping !
    integer, intent(out) :: step ! generator step size
    integer :: ii

    step=pi/float(N+1)

    do ii=1,N

      dr(ii)=2.0d0*bragg(nuc)*pi*sin(step*float(ii))/&
          &(1.0d0+2.0d0*cos(step*float(ii))+cos(step*float(ii))**2)

    end do

  end subroutine get_abcissas_z_1st

  subroutine get_abcissas_z_2nd(N,nuc,ddr,step)
    ! 2nd derivative of r(x) with respect to x, see
    ! grid_differentiation_sign_2.txt

    integer, intent(in) :: N ! number of mesh points
    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(out) :: ddr(:) ! 2nd deriv. of abcissas, Becke mapping !
    integer, intent(out) :: step ! generator step size
    integer :: ii

    step=pi/float(N+1)

    do ii=1,N

      ddr(ii)=(-2.0d0*bragg(nuc)*pi**2)*(cos(step*float(ii))-2.0d0)/&
          &(1.0d0+2.0d0*cos(step*float(ii))+cos(step*float(ii))**2)

    end do

  end subroutine get_abcissas_z_2nd

  function reverse_abcissas(nuc,r)
    ! z(x(r)) reverse mapping function, see
    ! grid_differentiation_sign_2.txt
    !
    ! z=1/pi*arccos((a-r)/(a+r))

    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(in) :: r ! radii of abcissas, Becke mapping !
    real(dp) :: reverse_abcissas

    reverse_abcissas=1.0d0/pi*acos((bragg(nuc)-r)/(bragg(nuc)+r))

  end function reverse_abcissas

  function reverse_abcissas_1st(nuc,r)
    ! 1st derivative of z(x(r)) reverse mapping function with resp. to r, see
    ! grid_differentiation_sign_2.txt
    !
    ! be careful: can easily overflow

    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(in) :: r ! radii of abcissas, Becke mapping !
    real(dp) :: reverse_abcissas_1st

    reverse_abcissas_1st=1.0d0/pi*sqrt(bragg(nuc)/r)/(r+bragg(nuc))

  end function reverse_abcissas_1st

  function reverse_abcissas_2nd(nuc,r)
    ! 2nd derivative of z(x(r)) reverse mapping function with resp. to r, see
    ! grid_differentiation_sign_2.txt
    !
    ! be careful: can easily overflow

    integer, intent(in) :: nuc ! nuclear charge 
    real(dp), intent(in) :: r ! radii of abcissas, Becke mapping !
    real(dp) :: reverse_abcissas_2nd

    reverse_abcissas_2nd=-1.0d0/(2.0d0*pi)*sqrt(bragg(nuc)/r)/r*&
        &(bragg(nuc)+3.0d0*r)/(bragg(nuc)+r)**2

  end function reverse_abcissas_2nd

  FUNCTION bragg(nuc)

    INTEGER :: nuc
    REAL(dp) :: bragg,braggd(110)
    DATA braggd/&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,&
        &1.0_dp,1.0_dp/
    bragg=braggd(nuc)
    RETURN
  END FUNCTION bragg

  function exp_int(alpha,power,r)
    ! evaluate \int x**power*exp(alpha*x) dx at point r
    ! for formula see Bronstein
    ! assumes alpha<0 and power=>0 !
    ! WATCH OUT FOR SIGN OF alpha !

    real(dp), intent(in) :: alpha,r
    integer, intent(in) :: power
    real(dp) :: exp_int
    integer :: ii

    exp_int=0.0d0

    ! catch power<0
    if (power<0) then
      write(*,*) 'NEGATIVE POWERS NOT IMPLEMENTED !'
      STOP
    end if

    ! catch alpha>0
    if (alpha>0.0d0) then
      write(*,*) 'POSITIVE ALPHAS NOT IMPLEMENTED !'
      STOP
    end if

    ! catch r=0
    if (r==0.0d0) then
      exp_int=fak(power)/(alpha**(power+1))*(-1.0d0)**(power)
      return
    end if

    ! catch r=infty and alpha<0 (should always be !)
    if (abs(alpha*r)>75.0d0) then
      exp_int=0.0d0
      return
    end if

    exp_int=1.0d0/alpha*exp(alpha*r)

    do ii=1,power
      exp_int=1.0d0/alpha*r**ii*exp(alpha*r)-float(ii)/alpha*exp_int
    end do

  end function exp_int

end module integration
