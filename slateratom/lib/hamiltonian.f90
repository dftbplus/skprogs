module hamiltonian

  use common_accuracy, only : dp
  use common_constants
  use dft
  use broyden
  use utilities
  use zora_routines

  implicit none
  private

  public :: build_fock, build_coulomb_matrix
  public :: build_hf_ex_matrix, build_dft_exc_matrix

contains

  subroutine build_fock(iter,t,u,nuc,vconf,j,k,p,max_l,num_alpha,poly_order,&
      &problemsize,xcnr,num_mesh_points,weight,abcissa,rho,vxc,alpha,pot_old,&
      &pot_new,zora,broyden,mixing_factor,f)

    ! Main driver routine for Fock matrix build-up. Calls also mixer with
    ! potential matrix.

    real(dp), intent(in) :: t(0:,:,:),u(0:,:,:),j(0:,:,:,0:,:,:),k(0:,:,:,0:,:,:)
    real(dp), intent(in) :: vconf(0:,:,:)
    real(dp), intent(in) :: p(:,0:,:,:),weight(:),abcissa(:),alpha(0:,:),rho(:,:)
    real(dp), intent(in) :: pot_old(:,0:,:,:),vxc(:,:),mixing_factor
    real(dp), intent(out) :: f(:,0:,:,:),pot_new(:,0:,:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),problemsize,nuc,xcnr
    integer, intent(in) :: num_mesh_points,iter
    logical, intent(in) :: zora,broyden
    real(dp), allocatable :: j_matrix(:,:,:),k_matrix(:,:,:,:),p_total(:,:,:)
    real(dp), allocatable :: t_zora(:,:,:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,biter

    f=0.0d0

    allocate(j_matrix(0:max_l,problemsize,problemsize))
    allocate(k_matrix(2,0:max_l,problemsize,problemsize))
    allocate(p_total(0:max_l,problemsize,problemsize))
    allocate(t_zora(2,0:max_l,problemsize,problemsize))
    p_total=0.0d0
    t_zora=0.0d0

    ! form total densitymatrix supervector
    do ii=0,max_l
      do jj=1,problemsize
        do kk=1,problemsize
          p_total(ii,jj,kk)=p(1,ii,jj,kk)+p(2,ii,jj,kk)
        end do
      end do
    end do

    ! build coulomb and exchange potential matrices

    call build_coulomb_matrix(j,p_total,max_l,num_alpha,poly_order,alpha,j_matrix)

    if (xcnr==0) then
      call build_hf_ex_matrix(k,p,max_l,num_alpha,poly_order,alpha,k_matrix)
    else
      call build_dft_exc_matrix(max_l,num_alpha,poly_order,alpha,&
          &num_mesh_points,abcissa,weight,rho,vxc,xcnr,k_matrix)
    end if

    ! build mixer input 

      pot_new(1,:,:,:)=-real(nuc,dp)*u(:,:,:)+j_matrix(:,:,:)-k_matrix(1,:,:,:)
      pot_new(2,:,:,:)=-real(nuc,dp)*u(:,:,:)+j_matrix(:,:,:)-k_matrix(2,:,:,:)


    ! mixer
    biter=int((iter)/40)
    call mixing_driver(pot_old,pot_new,max_l,num_alpha,&
        &poly_order,problemsize,iter-biter*40,broyden,mixing_factor)

! Not sure: before or after mixer .... ? Potential .ne. Matrix elements
! Should be irrelevant once self-consistency is reached
    if (zora) then

      call zora_t_correction(1,t_zora,max_l,num_alpha,alpha,poly_order,&
          &num_mesh_points,weight,abcissa,vxc,rho,nuc,p,problemsize)

    end if


    ! finally build Fock matrix
    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1

              f(1,ii,ss,tt)=t(ii,ss,tt)+pot_new(1,ii,ss,tt)+vconf(ii,ss,tt)
              f(2,ii,ss,tt)=t(ii,ss,tt)+pot_new(2,ii,ss,tt)+vconf(ii,ss,tt)

              if (zora) then
                f(1,ii,ss,tt)=f(1,ii,ss,tt)+t_zora(1,ii,ss,tt)
                f(2,ii,ss,tt)=f(2,ii,ss,tt)+t_zora(2,ii,ss,tt)
              end if

            end do
          end do
        end do
      end do
    end do

    !  write(*,*) 'FOCK MATRIX'
    !  write(*,*) f

    deallocate(j_matrix)
    deallocate(k_matrix)
    deallocate(p_total)
    deallocate(t_zora)

  end subroutine build_fock

  subroutine build_coulomb_matrix(j,p,max_l,num_alpha,poly_order,alpha,j_matrix)

    ! Build Coulomb matrix to be added to the Fock matrix from Coulomb Supermatrix
    ! by multiplying with density matrix supervector

    real(dp), intent(in) :: j(0:,:,:,0:,:,:),p(0:,:,:)
    real(dp), intent(in) :: alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:)
    real(dp), intent(out) :: j_matrix(0:,:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww

    j_matrix=0.0d0

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

                        ! multiply coulomb supermatrix with total densitymatrix supervector
                        j_matrix(ii,ss,tt)=j_matrix(ii,ss,tt)+&
                            &j(ii,ss,tt,nn,uu,vv)*p(nn,uu,vv)

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

  end subroutine build_coulomb_matrix

  subroutine build_hf_ex_matrix(k,p,max_l,num_alpha,poly_order,alpha,k_matrix)

    ! Build Hartree-Fock exchange matrix to be added to the Fock matrix from 
    ! supermatrix by multiplying with density matrix supervector

    real(dp), intent(in) :: k(0:,:,:,0:,:,:),p(:,0:,:,:)
    real(dp), intent(in) :: alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:)
    real(dp), intent(out) :: k_matrix(:,0:,:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww

    k_matrix=0.0d0

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

                        ! multiply hf exchange supermatrix with densitymatrix supervector per spin
                        k_matrix(1,ii,ss,tt)=k_matrix(1,ii,ss,tt)+&
                            &k(ii,ss,tt,nn,uu,vv)*p(1,nn,uu,vv)
                        k_matrix(2,ii,ss,tt)=k_matrix(2,ii,ss,tt)+&
                            &k(ii,ss,tt,nn,uu,vv)*p(2,nn,uu,vv)

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

  end subroutine build_hf_ex_matrix

  subroutine build_dft_exc_matrix(max_l,num_alpha,poly_order,alpha,&
      &num_mesh_points,abcissa,weight,rho,vxc,xcnr,k_matrix)

    ! Build DFT exchange matrix to be added to the Fock matrix by calculating
    ! the single matrix elements and putting them together

    real(dp), intent(in) :: alpha(0:,:)
    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),xcnr,num_mesh_points
    real(dp), intent(in) :: weight(:),abcissa(:),rho(:,:),vxc(:,:)
    real(dp), intent(out) :: k_matrix(:,0:,:,:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,start
    real(dp) :: exc_matrixelement(2)

    k_matrix=0.0d0
    exc_matrixelement=0.0d0

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1

          tt=ss-1
          do ll=jj,num_alpha(ii)

            start=1
            if (ll==jj) start=kk

            do mm=start,poly_order(ii)
              tt=tt+1

              call dft_exc_matrixelement(num_mesh_points,weight,abcissa,rho,&
                  &vxc,xcnr,alpha(ii,jj),kk,&
                  &alpha(ii,ll),mm,ii,exc_matrixelement)

              k_matrix(1,ii,ss,tt)=exc_matrixelement(1)
              k_matrix(2,ii,ss,tt)=exc_matrixelement(2)
              k_matrix(1,ii,tt,ss)=exc_matrixelement(1)
              k_matrix(2,ii,tt,ss)=exc_matrixelement(2)

            end do
          end do
        end do
      end do
    end do


  end subroutine build_dft_exc_matrix

end module hamiltonian
