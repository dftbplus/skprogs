module totalenergy
  use accuracy
  use constants
  use dft
  implicit none
  private

  public :: total_energy, zora_total_energy
  
contains

  subroutine total_energy(t,u,nuc,vconf,j,k,p,max_l,num_alpha,poly_order,&
      &problemsize,xcnr,num_mesh_points,weight,abcissa,rho,exc,&
      &kinetic,nuclear,coulomb,exchange,confinement,etot)

    ! Calculate total energy for non-ZORA calculations

    real(dp), intent(in) :: t(0:,:,:),u(0:,:,:),j(0:,:,:,0:,:,:),k(0:,:,:,0:,:,:)
    real(dp), intent(in) :: vconf(0:,:,:),abcissa(:)
    real(dp), intent(in) :: p(:,0:,:,:),weight(:),rho(:,:),exc(:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),problemsize,nuc,xcnr
    integer, intent(in) :: num_mesh_points
    real(dp), intent(out) :: etot,kinetic,nuclear,coulomb,exchange,confinement
    real(dp) :: dummy1,dummy2,dummy3
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww
    real(dp), allocatable :: p_total(:,:,:)

    allocate(p_total(0:max_l,problemsize,problemsize))
    p_total=0.0d0

    etot=0.0d0
    kinetic=0.0d0
    nuclear=0.0d0
    confinement=0.0d0
    coulomb=0.0d0
    exchange=0.0d0
    dummy1=0.0d0
    dummy2=0.0d0

    ! Build total density matrix
    do ii=0,max_l
      do jj=1,problemsize
        do kk=1,problemsize
          p_total(ii,jj,kk)=p(1,ii,jj,kk)+p(2,ii,jj,kk)
        end do
      end do
    end do

    ! get total energy

    call core_hamiltonian_energies(t,u,vconf,p_total,max_l,num_alpha,&
        &poly_order,nuc,kinetic,nuclear,confinement)

    dummy1=nuclear+kinetic+confinement

    call coulomb_hf_ex_energy(j,k,p_total,max_l,num_alpha,poly_order,xcnr,&
        &coulomb,exchange)

    if (xcnr>0) then

      exchange=0.0d0
      call dft_exc_energy(num_mesh_points,rho,exc,weight,abcissa,&
          &xcnr,exchange)

    end if

    ! make sure total energy breakdown agrees with total energy

    if (xcnr==0) then
      etot=dummy1+0.5d0*coulomb+0.5d0*exchange
    else
      etot=dummy1+0.5d0*coulomb+exchange
    end if

    !  write(*,*) 'TOTAL ENERGY',hf_total_energy

  end subroutine total_energy

  subroutine zora_total_energy(t,u,nuc,vconf,j,k,p,max_l,num_alpha,poly_order,&
      &problemsize,xcnr,num_mesh_points,weight,abcissa,rho,exc,vxc,&
      &eigval_scaled,occ,kinetic,nuclear,coulomb,exchange,confinement,etot)

    ! Calculate total energy for ZORA, note that total energy is not well defined
    ! here ...

    real(dp), intent(in) :: t(0:,:,:),u(0:,:,:),j(0:,:,:,0:,:,:),k(0:,:,:,0:,:,:)
    real(dp), intent(in) :: vconf(0:,:,:),abcissa(:),eigval_scaled(:,0:,:)
    real(dp), intent(in) :: occ(:,0:,:)
    real(dp), intent(in) :: p(:,0:,:,:),weight(:),rho(:,:),exc(:),vxc(:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),problemsize,nuc,xcnr
    integer, intent(in) :: num_mesh_points
    real(dp), intent(out) :: etot,kinetic,nuclear,coulomb,exchange,confinement
    real(dp) :: dummy1,dummy2,dummy3(2),eigsum
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww
    real(dp), allocatable :: p_total(:,:,:)

    allocate(p_total(0:max_l,problemsize,problemsize))
    p_total=0.0d0

    etot=0.0d0
    kinetic=0.0d0
    nuclear=0.0d0
    confinement=0.0d0
    coulomb=0.0d0
    exchange=0.0d0
    dummy1=0.0d0
    dummy2=0.0d0
    eigsum=0.0d0

    ! Build total density matrix
    do ii=0,max_l
      do jj=1,problemsize
        do kk=1,problemsize
          p_total(ii,jj,kk)=p(1,ii,jj,kk)+p(2,ii,jj,kk)
        end do
      end do
    end do

    ! get total energy

    call core_hamiltonian_energies(t,u,vconf,p_total,max_l,num_alpha,&
        &poly_order,nuc,kinetic,nuclear,confinement)

    ! sum of occupied eigenvalues
    do ii=1,2
      do jj=0,max_l
        do kk=1,problemsize
          eigsum=eigsum+eigval_scaled(ii,jj,kk)*occ(ii,jj,kk)
        end do
      end do
    end do

    kinetic=eigsum

    call coulomb_hf_ex_energy(j,k,p_total,max_l,num_alpha,poly_order,xcnr,&
        &coulomb,exchange)

    exchange=0.0d0
    call dft_exc_energy(num_mesh_points,rho,exc,weight,abcissa,&
        &xcnr,exchange)

    call dft_vxc_energy(num_mesh_points,rho,vxc,weight,abcissa,&
        &xcnr,dummy3)

    dummy2=dummy3(1)+dummy3(2)

    etot=eigsum-0.5d0*coulomb+exchange-dummy2

    !  write(*,*) 'ZORA TOTAL ENERGY'

  end subroutine zora_total_energy

  subroutine coulomb_hf_ex_energy(j,k,p_total,max_l,num_alpha,poly_order,xcnr,&
      &coulomb,exchange)

    ! get Hartee-Fock exchange and Coulomb contributions to total energy
    ! by multiplying j and k supermatrixes with the density matrix supervector
    ! twice

    real(dp), intent(in) :: j(0:,:,:,0:,:,:),k(0:,:,:,0:,:,:)
    real(dp), intent(in) :: p_total(0:,:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),xcnr
    real(dp), intent(out) :: coulomb,exchange
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww


    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              do nn=0,max_l
                uu=0
                do oo=1,num_alpha(nn)
                  do pp=1,poly_order(nn)
                    uu=uu+1
                    vv=0
                    do qq=1,num_alpha(nn)
                      do rr=1,poly_order(nn)
                        vv=vv+1

                        coulomb=coulomb+p_total(ii,ss,tt)*j(ii,ss,tt,nn,uu,vv)*&
                            &p_total(nn,uu,vv)

                        if (xcnr==0) then
                          exchange=exchange-0.5d0*p_total(ii,ss,tt)*&
                              &k(ii,ss,tt,nn,uu,vv)*p_total(nn,uu,vv)
                        end if

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine coulomb_hf_ex_energy

  subroutine core_hamiltonian_energies(t,u,vconf,p_total,max_l,num_alpha,&
      &poly_order,nuc,kinetic,nuclear,confinement)

    ! Core Hamiltonian contributions to total energy by multiplying the
    ! t,u,vconf supervectors with the density matrix supervector once

    real(dp), intent(in) :: t(0:,:,:),u(0:,:,:)
    real(dp), intent(in) :: vconf(0:,:,:)
    real(dp), intent(in) :: p_total(0:,:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),nuc
    real(dp), intent(out) :: kinetic,nuclear,confinement
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              kinetic=kinetic+t(ii,ss,tt)*p_total(ii,ss,tt)
              nuclear=nuclear-float(nuc)*u(ii,ss,tt)*p_total(ii,ss,tt)
              confinement=confinement+vconf(ii,ss,tt)*p_total(ii,ss,tt)
            end do
          end do
        end do
      end do
    end do

  end subroutine core_hamiltonian_energies

end module totalenergy
