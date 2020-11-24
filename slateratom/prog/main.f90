program HFAtom
  use accuracy
  use globals
  use integration
  use input
  use core_overlap
  use coulomb_hfex
  use densitymatrix
  use hamiltonian
  use diagonalizations
  use output
  use totalenergy
  use density
  use dft
  use utilities
  use zora_routines
  use cmdargs
  implicit none

  integer :: iter
  integer, allocatable :: qnvalorbs(:,:)

  call parse_command_arguments()
  call read_input_1(nuc,max_l,occ_shells,maxiter,poly_order,&
       &min_alpha,max_alpha,num_alpha,generate_alpha,alpha,&
       &conf_r0,conf_power,num_occ,num_power,num_alphas,xcnr,&
       &eigprint,zora,broyden,mixing_factor,xalpha_const)

  problemsize=num_power*num_alphas

! first index reserved for spin
  allocate(occ(2,0:max_l,problemsize))
  allocate(qnvalorbs(2, 0:max_l))

  call read_input_2(occ,max_l,occ_shells, qnvalorbs)

! fix number of mesh points depending on nuclear charge
  num_mesh_points=500
  if (nuc>10) num_mesh_points=750
  if (nuc>18) num_mesh_points=1000
  if (nuc>36) num_mesh_points=1250
  if (nuc>54) num_mesh_points=1500

  call echo_input(nuc,max_l,occ_shells,maxiter,poly_order,num_alpha,alpha,&
       &conf_r0,conf_power,occ,num_occ,num_power,num_alphas,xcnr,zora,&
       &num_mesh_points,xalpha_const)

! allocate global stuff and zero out
  call allocate_globals

! generate radial integration mesh
  call gauss_chebyshev_becke_mesh(num_mesh_points,nuc,weight,abcissa, dzdr, &
      &d2zdr2, dz)

! check mesh accuracy
  call check_accuracy(weight,abcissa,num_mesh_points,max_l,&
      &num_alpha,alpha,poly_order)

  if (xcnr >= 2) then
    write (*, "(A,/)") "LDA/PBE ROUTINES: LIBXC IMPLEMENTATION"
  end if
!!! OLD hand-coded xc implementation
!  if (xcnr == 2) then
!    write (*, "(A,/)") "LDA ROUTINES: BURKE IMPLEMENTATION"
!  end if
!  if (xcnr == 3) then
!    write (*, "(A,/)") "PBE ROUTINES: BURKE IMPLEMENTATION"
!  end if
!  if (xcnr > 3) then
!    write (*, "(A,/)") "STOP: Only xcnr <=3 supported without libxc"
!  end if
!!!
  

! Build supervectors

  write(*,'(A)') 'Startup: Building Supervectors'
  call overlap(s,max_l,num_alpha,alpha,poly_order)
  call nuclear(u,max_l,num_alpha,alpha,poly_order)
  call kinetic(t,max_l,num_alpha,alpha,poly_order)
  call confinement(vconf,max_l,num_alpha,alpha,poly_order,conf_r0,conf_power)

! test for linear dependency
  call diagonalize_overlap(max_l,num_alpha,poly_order,s)

! Build supermatrices

  write(*,'(A)') 'Startup: Building Supermatrices'
  call coulomb(j,max_l,num_alpha,alpha,poly_order,u,s)
  if (xcnr==0) call hfex(k,max_l,num_alpha,alpha,poly_order,problemsize)

! convergence flag
  final=.false.
    
! dft start potential
  if (xcnr>0) call dft_start_pot(abcissa,num_mesh_points,nuc,vxc)

! build initial fock matrix, core hamiltonian only
  write(*,'(A)') 'Startup: Building Initial Fock Matrix'
  write(*,'(A)') ' '

! do not confuse mixer
  pot_old=0.0d0

! kinetic energy, nuclear-electron, and confinement matrix elements
! which are constant during SCF
  call build_fock(0,t,u,nuc,vconf,j,k,p,max_l,num_alpha,poly_order,problemsize,&
       &xcnr,num_mesh_points,weight,abcissa,rho,vxc,alpha,pot_old,pot_new,&
       &zora,broyden,mixing_factor,f)

  do iter=1,maxiter

    write(*,'(A,I5)') 'Iteration :',iter

    pot_old=pot_new

! diagonalize 
    call diagonalize(max_l,num_alpha,poly_order,f,s,cof,eigval)


! build density matrix  
    call densmatrix(problemsize,max_l,occ,cof,p)

! get electron density and derivatives and exc related potentials and
! energy densities

    call density_grid(p,max_l,num_alpha,poly_order,alpha,num_mesh_points,&
        &abcissa, dzdr, d2zdr2, dz, xcnr,rho,drho,ddrho,vxc,exc,xalpha_const)

! Build Fock matrix and get total energy during SCF
    call build_fock(iter,t,u,nuc,vconf,j,k,p,max_l,num_alpha,poly_order,&
         &problemsize,xcnr,num_mesh_points,weight,abcissa,rho,vxc,alpha,&
         &pot_old,pot_new,zora,broyden,mixing_factor,f)
    
    call total_energy(t,u,nuc,vconf,j,k,p,max_l,num_alpha,poly_order,&
           &problemsize,xcnr,num_mesh_points,weight,abcissa,rho,exc,&
           &kinetic_energy,nuclear_energy,coulomb_energy,exchange_energy,&
           &conf_energy,total_ene)

    if (.not.zora) then
! non-rel. total energy during SCF meaningless for ZORA
! but energy contributions needed once SCF converged, so surpress output
       write(*,'(A,F18.6,A)') 'TOTAL ENERGY',total_ene,' Hartree'
    end if

    call check_convergence(pot_old,pot_new,max_l,problemsize,iter,&
         &change_max,final)

    write(*,'(A,E20.12)') 'CHANGE in potential matrix', change_max

! converged, nuke
    if (final) exit

! check conservation of number of electrons during SCF
    call check_electron_number(cof,s,occ,max_l,num_alpha,&
         &poly_order,problemsize)

    write(*,*) ' '

  end do

! output

  if (eigprint) then 
    call write_eigvec(max_l,num_alpha,alpha,poly_order,&
         &eigval,cof)
    call write_moments(max_l,num_alpha,alpha,poly_order,problemsize,cof)
    call cusp_values(max_l,occ,cof,p,alpha,num_alpha,poly_order,nuc)
  end if


  call write_eigval(max_l,num_alpha,poly_order,eigval)
  call write_energies(kinetic_energy,nuclear_energy,coulomb_energy,&
       &exchange_energy,conf_energy,total_ene,.false.)

  if (zora) then
    call scaled_zora(eigval,max_l,num_alpha,alpha,&
         &poly_order,problemsize,num_mesh_points,weight,abcissa,&
         &vxc,rho,nuc,p,t,cof,occ,eigval_scaled,zora_ekin)

    write(*,'(A)') 'Scaled Scalar-Relativistic ZORA EIGENVALUES and ENERGY'
    write(*,'(A)') '------------------------------------------------------'
    call write_eigval(max_l,num_alpha,poly_order,eigval_scaled)
    call write_energies(zora_ekin,nuclear_energy,coulomb_energy,&
       &exchange_energy,conf_energy,total_ene,.true.)
  end if
!
  write(*,'(A,E20.12)') 'Potential Matrix Elements converged to ', change_max
  write(*,'(A)') ' '

  if (zora) then
    call write_energies_tagged(zora_ekin,nuclear_energy,coulomb_energy,&
        &exchange_energy,conf_energy,0.0d0,zora, eigval_scaled, occ)
  else
    call write_energies_tagged(kinetic_energy,nuclear_energy,coulomb_energy,&
        &exchange_energy,conf_energy,total_ene,zora, eigval, occ)
  end if

  call write_potentials_file_standard(num_mesh_points,abcissa,weight,&
       &vxc,rho,nuc,p,max_l,num_alpha,poly_order,alpha,problemsize)

  call write_densities_file_standard(num_mesh_points,abcissa,weight,&
             &rho,drho,ddrho)

  ! Write wave functions and eventually invert to have positive starting
  ! gradient
  call write_waves_file_standard(num_mesh_points, abcissa, weight,&
      &alpha, num_alpha, poly_order,max_l, problemsize, occ, qnvalorbs, cof)

  call write_wave_coeffs_file(max_l, num_alpha, poly_order, cof, alpha, &
      &occ, qnvalorbs)

end program HFAtom
