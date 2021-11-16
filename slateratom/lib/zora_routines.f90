module zora_routines

  use common_accuracy, only : dp
  use common_constants
  use coulomb_potential
  use density

  implicit none
  private

  public :: zora_t_correction,scaled_zora

contains

  subroutine zora_t_correction(mode,t,max_l,num_alpha,alpha,poly_order,&
      &num_mesh_points,weight,abcissa,vxc,rho,nuc,p,problemsize)

    ! ZORA relativistic correction to kinetic energy matrix elements
    ! mode=1: correction to kinetic energy matrix elements
    ! mode=2: additional terms for scaling matrix elements

    real(dp), intent(out) :: t(:,0:,:,:)
    integer, intent(in) :: max_l,num_mesh_points,mode
    integer, intent(in) :: num_alpha(0:),nuc,problemsize
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:,:),weight(:),abcissa(:),vxc(:,:),rho(:,:)
    real(dp), intent(in) :: p(:,0:,:,:)
    real(dp), allocatable :: kappa(:,:),kappa2(:,:),vtot(:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,start

    allocate(kappa(2,num_mesh_points))
    allocate(kappa2(2,num_mesh_points))
    allocate(vtot(2,num_mesh_points))

    t=0.0d0

    call potential_to_mesh(num_mesh_points,abcissa,&
        &vxc,rho,nuc,p,max_l,num_alpha,poly_order,alpha,problemsize,vtot)

    call kappa_to_mesh(num_mesh_points,vtot,kappa,kappa2)

    do ii=0,max_l
      nn=0
      do jj=1,num_alpha(ii)
        do ll=1,poly_order(ii)
          nn=nn+1

          oo=nn-1
          do kk=jj,num_alpha(ii)

            start=1
            if (kk==jj) start=ll

            do mm=start,poly_order(ii)
              oo=oo+1

              ! kinetic energy correction depends on spin via potential 

              if (mode==1) then

                t(1,ii,nn,oo)=kinetic_part_1(num_mesh_points,weight,abcissa,&
                    &kappa(1,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)+kinetic_part_2(num_mesh_points,weight,abcissa,&
                    &kappa(1,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)*dfloat(ii*(ii+1))

                t(2,ii,nn,oo)=kinetic_part_1(num_mesh_points,weight,abcissa,&
                    &kappa(2,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)+kinetic_part_2(num_mesh_points,weight,abcissa,&
                    &kappa(2,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)*dfloat(ii*(ii+1))

              end if

              if (mode==2) then

              ! calculate matrix elements needed for scaled ZORA 
              ! prefactor 1/2 is included as the same subroutines as for t are
              ! used

                t(1,ii,nn,oo)=kinetic_part_1(num_mesh_points,weight,abcissa,&
                    &kappa2(1,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)+kinetic_part_2(num_mesh_points,weight,abcissa,&
                    &kappa2(1,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)*dfloat(ii*(ii+1))

                t(2,ii,nn,oo)=kinetic_part_1(num_mesh_points,weight,abcissa,&
                    &kappa2(2,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)+kinetic_part_2(num_mesh_points,weight,abcissa,&
                    &kappa2(2,:),alpha(ii,jj),ll,alpha(ii,kk),&
                    &mm,ii)*dfloat(ii*(ii+1))

              end if

              t(1,ii,oo,nn)=t(1,ii,nn,oo)
              t(2,ii,oo,nn)=t(2,ii,nn,oo)

            end do
          end do
        end do
      end do
    end do

!    write(*,'(A)') 'SR-ZORA KINETIC ENERGY CORRECTION'

    deallocate(kappa)
    deallocate(kappa2)
    deallocate(vtot)

  end subroutine zora_t_correction

  subroutine scaled_zora(eigval,max_l,num_alpha,alpha,&
             &poly_order,problemsize,num_mesh_points,weight,abcissa,&
             &vxc,rho,nuc,p,t,cof,occ,eigval_scaled,zora_ekin)

    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),problemsize
    real(dp), intent(in) :: eigval(:,0:,:),alpha(0:,:),cof(:,0:,:,:)
    real(dp), intent(in) :: occ(:,0:,:),t(0:,:,:)
    integer, intent(in) :: num_mesh_points,nuc
    real(dp), intent(in) :: weight(:),abcissa(:),vxc(:,:),rho(:,:),p(:,0:,:,:)
    real(dp), intent(out) :: eigval_scaled(:,0:,:),zora_ekin
    real(dp), allocatable :: zscale(:,:,:,:),zscale2(:,:,:,:)
    real(dp) :: dummy1,dummy2,tsol2,zora_ekin1,zora_ekin2
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww

    allocate(zscale(2,0:max_l,problemsize,problemsize))
    allocate(zscale2(2,0:max_l,problemsize,problemsize))
    zscale=0.0d0
    zscale2=0.0d0
    eigval_scaled=0.0d0
    zora_ekin=0.0d0
    zora_ekin1=0.0d0
    zora_ekin2=0.0d0
    tsol2=1.0_dp/cc**2

    call zora_t_correction(1,zscale,max_l,num_alpha,alpha,poly_order,&
      &num_mesh_points,weight,abcissa,vxc,rho,nuc,p,problemsize)
    call zora_t_correction(2,zscale2,max_l,num_alpha,alpha,poly_order,&
      &num_mesh_points,weight,abcissa,vxc,rho,nuc,p,problemsize)

! First get scaled eigenvalues

! Sum over all angular momenta
    do ii=0,max_l
! Sum over all eigenvectors
    do jj=1,num_alpha(ii)*poly_order(ii)
      oo=0
      dummy1=0.0d0
      dummy2=0.0d0
! sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
      do kk=1,num_alpha(ii)
        do ll=1,poly_order(ii)
          oo=oo+1
          pp=0
! other sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
          do mm=1,num_alpha(ii)
            do nn=1,poly_order(ii)
              pp=pp+1
! occupation numbers do not enter here
              dummy1=dummy1+cof(1,ii,pp,jj)*cof(1,ii,oo,jj)*&
                     &tsol2*(zscale(1,ii,oo,pp)+&
                     &0.5d0*(zscale2(1,ii,oo,pp)+t(ii,oo,pp)))
              dummy2=dummy2+cof(2,ii,pp,jj)*cof(2,ii,oo,jj)*&
                     &tsol2*(zscale(2,ii,oo,pp)+&
                     &0.5d0*(zscale2(2,ii,oo,pp)+t(ii,oo,pp)))
            end do
          end do
        end do
      end do



      eigval_scaled(1,ii,jj)=eigval(1,ii,jj)/(1.0d0+dummy1)
      eigval_scaled(2,ii,jj)=eigval(2,ii,jj)/(1.0d0+dummy2)
    end do  
    end do

! Now ZORA kinetic energy

      dummy1=0.0d0
      dummy2=0.0d0
! Sum over all angular momenta
    do ii=0,max_l
! Sum over all eigenvectors
    do jj=1,num_alpha(ii)*poly_order(ii)
      oo=0
! sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
      do kk=1,num_alpha(ii)
        do ll=1,poly_order(ii)
          oo=oo+1
          pp=0
! other sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
          do mm=1,num_alpha(ii)
            do nn=1,poly_order(ii)
              pp=pp+1
! dummy contains the non-relativistic kinetic energy operator applied
! to the relativistic ZORA wavefunction, debug only
!              dummy1=dummy1+occ(1,ii,jj)*cof(1,ii,pp,jj)*cof(1,ii,oo,jj)*t(ii,oo,pp)
!              dummy2=dummy2+occ(2,ii,jj)*cof(2,ii,pp,jj)*cof(2,ii,oo,jj)*t(ii,oo,pp)
              zora_ekin1=zora_ekin1+&
                &occ(1,ii,jj)*cof(1,ii,pp,jj)*cof(1,ii,oo,jj)*&
                &(t(ii,oo,pp)+zscale(1,ii,oo,pp)-&
                &eigval_scaled(1,ii,jj)*tsol2*(0.5d0*(&
                &zscale2(1,ii,oo,pp)+t(ii,oo,pp))+zscale(1,ii,oo,pp)))
              zora_ekin2=zora_ekin2+&
                &occ(2,ii,jj)*cof(2,ii,pp,jj)*cof(2,ii,oo,jj)*&
                &(t(ii,oo,pp)+zscale(2,ii,oo,pp)-&
                &eigval_scaled(2,ii,jj)*tsol2*(0.5d0*(&
                &zscale2(2,ii,oo,pp)+t(ii,oo,pp))+zscale(2,ii,oo,pp)))
            end do
          end do
        end do
      end do
    end do  
    end do
!    write(*,*) 'SCAL2 ',dummy1,dummy2,zora_ekin1,zora_ekin2

    zora_ekin=zora_ekin1+zora_ekin2

    deallocate(zscale)
    deallocate(zscale2)

  end subroutine scaled_zora

  function kinetic_part_1(num_mesh_points,weight,abcissa,kappa,&
      &alpha1,poly1,alpha2,poly2,l)

    ! get 0.5*\int_0^\inf r^2 kappa (d/dr R_A) (d/dr R_B) dr
    ! pass either up or down total potential as kappa

    real(dp), intent(in) :: weight(:),abcissa(:),kappa(:)
    real(dp), intent(in) :: alpha1,alpha2
    integer, intent(in) :: num_mesh_points
    integer, intent(in) :: poly1,poly2,l
    integer :: ii,jj,kk,ll,mm,nn,oo
    real(dp) :: kinetic_part_1

    kinetic_part_1=0.0d0

    do ii=1,num_mesh_points

      kinetic_part_1=kinetic_part_1+weight(ii)*kappa(ii)*&
          &basis_1st_times_basis_1st_times_r2(alpha1,poly1,alpha2,poly2,l,abcissa(ii))

    end do

    kinetic_part_1=kinetic_part_1*0.5d0

  end function kinetic_part_1

  function kinetic_part_2(num_mesh_points,weight,abcissa,kappa,alpha1,&
      &poly1,alpha2,poly2,l)

    ! get \int_0^\inf R_B R_A kappa dr; multiply by l(l+1) in calling routine
    ! pass either up or down total potential as kappa

    real(dp), intent(in) :: weight(:),abcissa(:),kappa(:)
    real(dp), intent(in) :: alpha1,alpha2
    integer, intent(in) :: num_mesh_points
    integer, intent(in) :: poly1,poly2,l
    integer :: ii,jj,kk,ll,mm,nn,oo
    real(dp) :: kinetic_part_2

    kinetic_part_2=0.0d0

    do ii=1,num_mesh_points

      kinetic_part_2=kinetic_part_2+weight(ii)*kappa(ii)*&
          &basis_times_basis(alpha1,poly1,alpha2,poly2,l,abcissa(ii))

    end do

    kinetic_part_2=kinetic_part_2*0.5d0

  end function kinetic_part_2

  subroutine kappa_to_mesh(num_mesh_points,vtot,kappa,kappa2)

    ! kappa=V/(2*c^2-V), V total potential, c speed of light
    ! kappa2=kappa^2, i.e. square of kappa
    
    integer, intent(in) :: num_mesh_points
    real(dp), intent(in) :: vtot(:,:)
    real(dp), intent(out) :: kappa(:,:),kappa2(:,:)
    integer :: ii

    real(dp), parameter :: tsol2 =2.0_dp*cc**2

    do ii=1,num_mesh_points

      kappa(1,ii)=vtot(1,ii)/(tsol2-vtot(1,ii))
      kappa(2,ii)=vtot(2,ii)/(tsol2-vtot(2,ii))

      kappa2(1,ii)=kappa(1,ii)**2
      kappa2(2,ii)=kappa(2,ii)**2

    end do

  end subroutine kappa_to_mesh

  subroutine potential_to_mesh(num_mesh_points,abcissa,&
      &vxc,rho,nuc,p,max_l,num_alpha,poly_order,alpha,problemsize,vtot)

    ! get total potential on mesh, spinpolarized

    real(dp), intent(in) :: abcissa(:),vxc(:,:),p(:,0:,:,:),alpha(0:,:)
    real(dp), intent(in) :: rho(:,:)
    integer, intent(in) :: num_mesh_points,nuc,max_l,num_alpha(0:)
    integer, intent(in) :: poly_order(0:),problemsize
    real(dp), intent(out) :: vtot(:,:)
    real(dp), allocatable :: cpot(:),ptot(:,:,:)
    integer :: ii

    allocate(cpot(num_mesh_points))
    allocate(ptot(0:max_l,problemsize,problemsize))

    cpot=0.0d0
    ptot=0.0d0
    vtot=0.0d0

    ptot(:,:,:)=p(1,:,:,:)+p(2,:,:,:)

    call cou_pot(ptot(:,:,:),max_l,num_alpha,poly_order,alpha,problemsize,&
        &num_mesh_points,abcissa,cpot)

    do ii=1,num_mesh_points

      vtot(1,ii)=-float(nuc)/abcissa(ii)+cpot(ii)+vxc(ii,1)
      vtot(2,ii)=-float(nuc)/abcissa(ii)+cpot(ii)+vxc(ii,2)

    end do

    deallocate(cpot)
    deallocate(ptot)

  end subroutine potential_to_mesh
!
end module zora_routines
