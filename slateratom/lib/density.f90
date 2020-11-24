module density
  use accuracy
  use utilities
  implicit none
  private

  public :: density_at_point, density_at_point_1st, density_at_point_2nd
  public :: wavefunction, wavefunction_1st, wavefunction_2nd
  public :: basis, basis_1st, basis_2nd
  public :: basis_times_basis, basis_times_basis_1st, basis_times_basis_2nd
  public :: basis_1st_times_basis_1st, basis_2nd_times_basis_2nd
  public :: basis_times_basis_times_r2, basis_times_basis_1st_times_r2, &
      &basis_times_basis_2nd_times_r2, basis_times_basis_1st_times_r, &
      &basis_1st_times_basis_1st_times_r2

contains

  function density_at_point(p,max_l,num_alpha,poly_order,alpha,r)

    ! Calculate electron density at a radial point in space

    real(dp), intent(in) :: p(0:,:,:),r,alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:)
    real(dp) :: density_at_point
    integer :: ii,jj,kk,ll,mm,nn,oo,start

    density_at_point=0.0d0

    do ii=0,max_l
      ll=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ll=ll+1

          ! set global index correctly
          oo=ll-1
          do mm=jj,num_alpha(ii)

            ! catch start index for polynomials, different depending on alpha block
            start=1
            if (mm==jj) start=kk

            do nn=start,poly_order(ii)
              oo=oo+1

              if (ll==oo) then
                density_at_point=density_at_point+p(ii,ll,oo)*&
                    &basis_times_basis(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)
              end if

              if (ll/=oo) then
                density_at_point=density_at_point+2.0d0*p(ii,ll,oo)*&
                    &basis_times_basis(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)
              end if

            end do
          end do
        end do
      end do
    end do

  end function density_at_point

  function density_at_point_1st(p,max_l,num_alpha,poly_order,alpha,r)

    ! Calculate 1st derivative at a radial point in space analytically


    real(dp), intent(in) :: p(0:,:,:),r,alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:)
    real(dp) :: density_at_point_1st
    integer :: ii,jj,kk,ll,mm,nn,oo,start

    density_at_point_1st=0.0d0
    !
    !  do ii=0,max_l
    !    ll=0
    !    do jj=1,num_alpha(ii)
    !      do kk=1,poly_order(ii)
    !        ll=ll+1
    !        oo=0
    !        do mm=1,num_alpha(ii)
    !          do nn=1,poly_order(ii)
    !            oo=oo+1
    !            density_at_point_1st=density_at_point_1st+p(ii,ll,oo)*(&
    !!              &basis(alpha(ii,jj),kk,ii,r)*basis_1st(alpha(ii,mm),nn,ii,r)&
    !!              &+basis_1st(alpha(ii,jj),kk,ii,r)*basis(alpha(ii,mm),nn,ii,r))
    !            &basis_times_basis_1st(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
    !            &basis_times_basis_1st(alpha(ii,mm),nn,alpha(ii,jj),kk,ii,r))
    !          end do
    !        end do
    !      end do
    !    end do
    !  end do

    do ii=0,max_l
      ll=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ll=ll+1

          ! set global index correctly
          oo=ll-1
          do mm=jj,num_alpha(ii)

            ! catch start index for polynomials, different depending on alpha block
            start=1
            if (mm==jj) start=kk

            do nn=start,poly_order(ii)
              oo=oo+1

              if (ll==oo) then
                density_at_point_1st=density_at_point_1st+p(ii,ll,oo)*(&
                    &basis_times_basis_1st(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
                    &basis_times_basis_1st(alpha(ii,mm),nn,alpha(ii,jj),kk,ii,r))
              end if

              if (ll/=oo) then
                density_at_point_1st=density_at_point_1st+2.0d0*p(ii,ll,oo)*(&
                    &basis_times_basis_1st(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
                    &basis_times_basis_1st(alpha(ii,mm),nn,alpha(ii,jj),kk,ii,r))
              end if

            end do
          end do
        end do
      end do
    end do

  end function density_at_point_1st

  function density_at_point_2nd(p,max_l,num_alpha,poly_order,alpha,r)

    ! Calculate 2nd derivative at a radial point in space analytically


    real(dp), intent(in) :: p(0:,:,:),r,alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:)
    real(dp) :: density_at_point_2nd
    integer :: ii,jj,kk,ll,mm,nn,oo,start

    density_at_point_2nd=0.0d0
    !
    !  do ii=0,max_l
    !    ll=0
    !    do jj=1,num_alpha(ii)
    !      do kk=1,poly_order(ii)
    !        ll=ll+1
    !        oo=0
    !        do mm=1,num_alpha(ii)
    !          do nn=1,poly_order(ii)
    !            oo=oo+1
    !            density_at_point_2nd=density_at_point_2nd+p(ii,ll,oo)*(&
    !            &basis_times_basis_2nd(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
    !            &+2.0d0*basis_1st_times_basis_1st(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
    !            &basis_times_basis_2nd(alpha(ii,mm),nn,alpha(ii,jj),kk,ii,r))
    !          end do
    !        end do
    !      end do
    !    end do
    !  end do

    do ii=0,max_l
      ll=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ll=ll+1

          ! set global index correctly
          oo=ll-1
          do mm=jj,num_alpha(ii)

            ! catch start index for polynomials, different depending on alpha block
            start=1
            if (mm==jj) start=kk

            do nn=start,poly_order(ii)
              oo=oo+1

              if (ll==oo) then
                density_at_point_2nd=density_at_point_2nd+p(ii,ll,oo)*(&
                    &basis_times_basis_2nd(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
                    &+2.0d0*basis_1st_times_basis_1st(alpha(ii,jj),kk,alpha(ii,mm),&
                    &nn,ii,r)+&
                    &basis_times_basis_2nd(alpha(ii,mm),nn,alpha(ii,jj),kk,ii,r))
              end if

              if (ll/=oo) then
                density_at_point_2nd=density_at_point_2nd+2.0d0*p(ii,ll,oo)*(&
                    &basis_times_basis_2nd(alpha(ii,jj),kk,alpha(ii,mm),nn,ii,r)+&
                    &+2.0d0*basis_1st_times_basis_1st(alpha(ii,jj),kk,alpha(ii,mm),&
                    &nn,ii,r)+&
                    &basis_times_basis_2nd(alpha(ii,mm),nn,alpha(ii,jj),kk,ii,r))
              end if

            end do
          end do
        end do
      end do
    end do

  end function density_at_point_2nd

  function wavefunction(cof,alpha,num_alpha,poly_order,ang,r)

    ! Calculate value of wavefunction at a radial point in space


    integer, intent(in) :: num_alpha(0:),poly_order(0:)
    integer, intent(in) :: ang
    real(dp), intent(in) :: cof(:),alpha(0:,:),r
    real(dp) :: wavefunction
    integer :: ii,jj,kk

    wavefunction=0.0d0
    kk=0

    do ii=1,num_alpha(ang)
      do jj=1,poly_order(ang)
        kk=kk+1
        !      write(*,'(3I3,F12.6,I3,F12.6)') ang,ii,jj,alpha(ang,ii),jj+ang,cof(kk)
        wavefunction=wavefunction+cof(kk)*basis(alpha(ang,ii),jj,ang,r)
      end do
    end do

  end function wavefunction

  function wavefunction_1st(cof,alpha,num_alpha,poly_order,ang,r)

    ! Calculate value of 1st derivative of wavefunction at a radial point in
    ! space analytically


    integer, intent(in) :: num_alpha(0:),poly_order(0:)
    integer, intent(in) :: ang
    real(dp), intent(in) :: cof(:),alpha(0:,:),r
    real(dp) :: wavefunction_1st
    integer :: ii,jj,kk

    wavefunction_1st=0.0d0
    kk=0

    do ii=1,num_alpha(ang)
      do jj=1,poly_order(ang)
        kk=kk+1
        wavefunction_1st=wavefunction_1st+cof(kk)*basis_1st(alpha(ang,ii),jj,ang,r)
      end do
    end do

  end function wavefunction_1st

  function wavefunction_2nd(cof,alpha,num_alpha,poly_order,ang,r)

    ! Calculate value of 2nd derivative of wavefunction at a radial point in
    ! space analytically


    integer, intent(in) :: num_alpha(0:),poly_order(0:)
    integer, intent(in) :: ang
    real(dp), intent(in) :: cof(:),alpha(0:,:),r
    real(dp) :: wavefunction_2nd
    integer :: ii,jj,kk

    wavefunction_2nd=0.0d0
    kk=0

    do ii=1,num_alpha(ang)
      do jj=1,poly_order(ang)
        kk=kk+1
        wavefunction_2nd=wavefunction_2nd+cof(kk)*basis_2nd(alpha(ang,ii),jj,ang,r)
      end do
    end do

  end function wavefunction_2nd

  function basis(alpha,poly_order,l,r)

    ! Value of a primitive Slater basis function at a radial point in space
    ! See rmp_32_186_1960.pdf eqn. 3


    integer, intent(in) :: l,poly_order
    real(dp), intent(in) :: alpha,r
    real(dp) :: basis,normalization

    normalization=(2.0d0*alpha)**(poly_order+l)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*(poly_order+l)))

    ! catch 0^0
    if ((r==0.0d0).and.((poly_order+l-1)==0)) then
      basis=normalization*exp(-alpha*r)
    else
      basis=normalization*r**(poly_order+l-1)*exp(-alpha*r)
    end if

  end function basis

  function basis_1st(alpha,poly_order,l,r)

    ! Value of 1st derivative of a primitive Slater basis function at a radial 
    ! point in space


    integer, intent(in) :: l,poly_order
    real(dp), intent(in) :: alpha,r
    real(dp) :: basis_1st,normalization

    normalization=(2.0d0*alpha)**(poly_order+l)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*(poly_order+l)))

    ! catch 0^0, setting 0^0=1 and 0^-1=0.0
    if ((r==0.0d0).and.((poly_order+l-1)==0)) then
      basis_1st=normalization*(-alpha*exp(-alpha*r))
    else if ((r==0.0d0).and.((poly_order+l-2)==0)) then
      basis_1st=normalization*(float(poly_order+l-1)*&
          &exp(-alpha*r)-alpha*r**(poly_order+l-1)*exp(-alpha*r))
    else
      basis_1st=normalization*(float(poly_order+l-1)*r**(poly_order+l-2)*&
          &exp(-alpha*r)-alpha*r**(poly_order+l-1)*exp(-alpha*r))
    end if

  end function basis_1st

  function basis_2nd(alpha,poly_order,l,r)

    ! Value of 2nd derivative of a primitive Slater basis function at a radial 
    ! point in space


    integer, intent(in) :: l,poly_order
    real(dp), intent(in) :: alpha,r
    real(dp) :: basis_2nd,normalization

    normalization=(2.0d0*alpha)**(poly_order+l)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*(poly_order+l)))

    ! catch 0^0
    if ((r==0.0d0).and.((poly_order+l-3)==0)) then
      basis_2nd=normalization*(float(poly_order+l-1)*float(poly_order+l-2)*&
          &exp(-alpha*r))
    else if ((r==0.0d0).and.((poly_order+l-2)==0)) then
      basis_2nd=normalization*(-2.0d0*alpha*float(poly_order+l-1)*&
          &exp(-alpha*r))
    else if ((r==0.0d0).and.((poly_order+l-1)==0)) then
      basis_2nd=normalization*(alpha**2*exp(-alpha*r))
    else
      basis_2nd=normalization*(float(poly_order+l-1)*float(poly_order+l-2)*&
          &r**(poly_order+l-3)*exp(-alpha*r)-2.0d0*alpha*float(poly_order+l-1)*&
          &r**(poly_order+l-2)*exp(-alpha*r)+alpha**2*r**(poly_order+l-1)*&
          &exp(-alpha*r))
    end if

  end function basis_2nd

  function basis_times_basis(alpha,poly1,beta,poly2,l,r)
    ! Value of a product of two primitive Slater basis functions at a radial 
    ! point in space
    ! r^(m-1)*e^(-alpha*r)*r^(n-1)*exp(-beta*r)=r^(m+n-2)*exp(-(alpha+beta)*r)


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis,normalization1,normalization2
    real(dp) :: ab
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))


    ! catch 0^0
    if ((r==0.0d0).and.((m+n-2)==0)) then
      basis_times_basis=normalization1*normalization2*exp(ab*r)
    else
      basis_times_basis=normalization1*normalization2*&
          &r**(m+n-2)*exp(ab*r)
    end if

    if (abs(basis_times_basis)<1.0d-20) basis_times_basis=0.0d0

  end function basis_times_basis

  function basis_times_basis_1st(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a basis function with first the derivative of another
    ! basis function
    ! beta and poly2 are the arguments of the derivative


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis_1st,normalization1,normalization2
    real(dp) :: ab
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !

    ! catch 0^0, setting 0^0=1 and 0^1=0
    if ((r==0.0d0).and.((m+n-2)==0)) then
      basis_times_basis_1st=normalization1*normalization2*&
          &(-beta)*exp(ab*r)
    else if ((r==0.0d0).and.((m+n-3)==0)) then
      basis_times_basis_1st=normalization1*normalization2*&
          &(float(n-1))*exp(ab*r)
    else
      basis_times_basis_1st=normalization1*normalization2*&
          &(float(n-1)*r**(m+n-3)-beta*r**(n+m-2))*exp(ab*r)
    end if

    if (abs(basis_times_basis_1st)<1.0d-20) basis_times_basis_1st=0.0d0

  end function basis_times_basis_1st

  function basis_times_basis_2nd(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a basis function with the second derivative of 
    ! another basis function
    ! beta and poly2 are the arguments of the derivative


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis_2nd,normalization1,normalization2
    real(dp) :: ab,positive,negative
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !
    positive=float((n-1)*(n-2))*r**(m+n-4)+beta**2*r**(m+n-2)
    negative=float(2*(n-1))*beta*r**(n+m-3)

    basis_times_basis_2nd=normalization1*normalization2*&
        &(positive-negative)*exp(ab*r)

    if (abs(basis_times_basis_2nd)<1.0d-20) basis_times_basis_2nd=0.0d0

  end function basis_times_basis_2nd

  function basis_1st_times_basis_1st(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a first derivatives of basis functions


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_1st_times_basis_1st,normalization1,normalization2
    real(dp) :: ab,positive,negative
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !

    ! catch 0^0
    if ((r==0.0d0).and.((m+n-2)==0)) then
      positive=alpha*beta
    else if ((r==0.0d0).and.((m+n-4)==0)) then
      positive=float((m-1)*(n-1))
    else
      positive=float((m-1)*(n-1))*r**(m+n-4)+&
          &alpha*beta*r**(m+n-2)
    end if

    if ((r==0.0d0).and.((m+n-3)==0)) then
      negative=(alpha*float(n-1)+beta*float(m-1))
    else
      negative=(alpha*float(n-1)+beta*float(m-1))*r**(m+n-3)
    end if

    basis_1st_times_basis_1st=normalization1*normalization2*&
        &(positive-negative)*exp(ab*r)

    if (abs(basis_1st_times_basis_1st)<1.0d-20) basis_1st_times_basis_1st=0.0d0

  end function basis_1st_times_basis_1st

  function basis_2nd_times_basis_2nd(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a first derivatives of basis functions


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_2nd_times_basis_2nd,normalization1,normalization2
    real(dp) :: ab,positive,negative
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !
    positive=float((m-1)*(m-2)*(n-1)*(n-2))*r**(n+m-6)+&
        &r**(m+n-4)*(beta**2*float((m-1)*(m-2))+alpha**2*float((n-1)*(n-2))+&
        &alpha*beta*float(4*(m-1)*(n-1)))+&
        &alpha**2*beta**2*r**(m+n-2)

    negative=r**(m+n-5)*(beta*float(2*(n-1)*(m-1)*(m-2))+&
        &alpha*float(2*(m-1)*(n-1)*(n-2)))+&
        &r**(m+n-3)*(alpha*beta**2*float(2*(m-1))+&
        &beta*alpha**2*float(2*(n-1)))

    basis_2nd_times_basis_2nd=normalization1*normalization2*&
        &(positive-negative)*exp(ab*r)

    if (abs(basis_2nd_times_basis_2nd)<1.0d-20) basis_2nd_times_basis_2nd=0.0d0

  end function basis_2nd_times_basis_2nd

  function basis_times_basis_times_r2(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of two basis functions and r^2 in one go
    ! r^(m-1)*e^(-alpha*r)*r^(n-1)*exp(-beta*r) *r^2=r^(m+n)*exp(-(alpha+beta)*r)


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis_times_r2,normalization1,normalization2
    real(dp) :: ab
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    basis_times_basis_times_r2=normalization1*normalization2*&
        &r**(m+n)*exp(ab*r)

    if (abs(basis_times_basis_times_r2)<1.0d-20) basis_times_basis_times_r2=0.0d0

  end function basis_times_basis_times_r2

  function basis_times_basis_1st_times_r2(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a basis function with first the derivative of another
    ! basis function and r^2
    ! beta and poly2 are the arguments of the derivative


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis_1st_times_r2,normalization1,normalization2
    real(dp) :: ab
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !
    basis_times_basis_1st_times_r2=normalization1*normalization2*&
        &(float(n-1)*r**(m+n-1)-beta*r**(n+m))*exp(ab*r)

    if (abs(basis_times_basis_1st_times_r2)<1.0d-20) &
        &basis_times_basis_1st_times_r2=0.0d0

  end function basis_times_basis_1st_times_r2

  function basis_times_basis_2nd_times_r2(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a basis function with the second derivative of 
    ! another basis function and r^2
    ! beta and poly2 are the arguments of the derivative


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis_2nd_times_r2,normalization1,normalization2
    real(dp) :: ab,positive,negative
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !
    positive=float((n-1)*(n-2))*r**(m+n-2)+beta**2*r**(m+n)
    negative=float(2*(n-1))*beta*r**(n+m-1)

    basis_times_basis_2nd_times_r2=normalization1*normalization2*&
        &(positive-negative)*exp(ab*r)

    if (abs(basis_times_basis_2nd_times_r2)<1.0d-20) &
        &basis_times_basis_2nd_times_r2=0.0d0

  end function basis_times_basis_2nd_times_r2

  function basis_times_basis_1st_times_r(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a basis function with first the derivative of another
    ! basis function and r
    ! beta and poly2 are the arguments of the derivative


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_times_basis_1st_times_r,normalization1,normalization2
    real(dp) :: ab
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !
    basis_times_basis_1st_times_r=normalization1*normalization2*&
        &(float(n-1)*r**(m+n-2)-beta*r**(n+m-1))*exp(ab*r)

    if (abs(basis_times_basis_1st_times_r)<1.0d-20) &
        &basis_times_basis_1st_times_r=0.0d0

  end function basis_times_basis_1st_times_r

  function basis_1st_times_basis_1st_times_r2(alpha,poly1,beta,poly2,l,r)
    ! evaluation of product of a first derivatives of basis functions and r^2


    integer, intent(in) :: l,poly1,poly2
    real(dp), intent(in) :: alpha,beta,r
    real(dp) :: basis_1st_times_basis_1st_times_r2,normalization1,normalization2
    real(dp) :: ab,positive,negative
    integer :: m,n

    m=poly1+l
    n=poly2+l
    ab=-(alpha+beta)

    normalization1=(2.0d0*alpha)**(m)*sqrt(2.0d0*alpha)/&
        &sqrt(fak(2*m))
    normalization2=(2.0d0*beta)**(n)*sqrt(2.0d0*beta)/&
        &sqrt(fak(2*n))

    ! WARNING: without summing negative and positive contributions independently
    ! zora becomes completely unstable !
    positive=float((m-1)*(n-1))*r**(m+n-2)+&
        &alpha*beta*r**(m+n)
    negative=(alpha*float(n-1)+beta*float(m-1))*r**(m+n-1)

    basis_1st_times_basis_1st_times_r2=normalization1*normalization2*&
        &(positive-negative)*exp(ab*r)

    if (abs(basis_1st_times_basis_1st_times_r2)<1.0d-20) &
        &basis_1st_times_basis_1st_times_r2=0.0d0

  end function basis_1st_times_basis_1st_times_r2

end module density
