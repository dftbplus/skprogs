program HFAtom

  use common_accuracy, only : dp
  use common_message, only : error
  use integration, only : gauss_chebyshev_becke_mesh
  use input, only : read_input_1, read_input_2, echo_input
  use core_overlap, only : overlap, nuclear, kinetic, confinement
  use coulomb_hfex, only : coulomb, hfex, hfex_lr
  use densitymatrix, only : densmatrix
  use hamiltonian, only : build_hamiltonian
  use diagonalizations, only : diagonalize, diagonalize_overlap
  use output, only : write_eigvec, write_eigval, write_moments, write_energies,&
      & write_energies_tagged, write_potentials_file_standard, write_densities_file_standard,&
      & write_waves_file_standard, write_wave_coeffs_file, cusp_values
  use totalenergy, only : getTotalEnergy, getTotalEnergyZora
  use dft, only : check_accuracy, dft_start_pot, density_grid
  use utilities, only : check_electron_number, check_convergence
  use zora_routines, only : scaled_zora
  use cmdargs, only : parse_command_arguments
  use common_poisson, only : TBeckeGridParams
  use xcfunctionals, only : xcFunctional
  use globals

  implicit none


  !! current SCF step
  integer :: iScf

  !! quantum numbers of wavefunctions to be written
  integer, allocatable :: qnvalorbs(:,:)

  !!
  real(dp) :: x_en_2

  !! range-separation parameter
  real(dp) :: omega

  !! CAM alpha parameter
  real(dp) :: camAlpha

  !! CAM beta parameter
  real(dp) :: camBeta

  !! holds parameters, defining a Becke integration grid
  type(TBeckeGridParams) :: grid_params

  call parse_command_arguments()
  call read_input_1(nuc, max_l, occ_shells, maxiter, poly_order, min_alpha, max_alpha, num_alpha,&
      & tAutoAlphas, alpha, conf_r0, conf_power, num_occ, num_power, num_alphas, xcnr,&
      & tPrintEigvecs, tZora, tBroyden, mixing_factor, xalpha_const, omega, camAlpha, camBeta,&
      & grid_params)

  problemsize = num_power * num_alphas

  ! first index reserved for spin
  allocate(occ(2, 0:max_l, problemsize))
  allocate(qnvalorbs(2, 0:max_l))

  call read_input_2(occ, max_l, occ_shells, qnvalorbs)

  ! fix number of mesh points depending on nuclear charge
  num_mesh_points = 500
  if (nuc > 10) num_mesh_points = 750
  if (nuc > 18) num_mesh_points = 1000
  if (nuc > 36) num_mesh_points = 1250
  if (nuc > 54) num_mesh_points = 1500

  call echo_input(nuc, max_l, occ_shells, maxiter, poly_order, num_alpha, alpha, conf_r0,&
      & conf_power, occ, num_occ, num_power, num_alphas, xcnr, tZora, num_mesh_points, xalpha_const)

  ! allocate global stuff and zero out
  call allocate_globals()

  ! generate radial integration mesh
  call gauss_chebyshev_becke_mesh(num_mesh_points, nuc, weight, abcissa, dzdr, d2zdr2, dz)

  ! check mesh accuracy
  call check_accuracy(weight, abcissa, num_mesh_points, max_l, num_alpha, alpha, poly_order)

  ! build supervectors
  write(*, '(A)') 'Startup: Building Supervectors'
  call overlap(ss, max_l, num_alpha, alpha, poly_order)
  call nuclear(uu, max_l, num_alpha, alpha, poly_order)
  call kinetic(tt, max_l, num_alpha, alpha, poly_order)
  call confinement(vconf, max_l, num_alpha, alpha, poly_order, conf_r0, conf_power)

  ! test for linear dependency
  call diagonalize_overlap(max_l, num_alpha, poly_order, ss)

  ! build supermatrices
  write(*, '(A)') 'Startup: Building Supermatrices'
  call coulomb(jj, max_l, num_alpha, alpha, poly_order, uu, ss)
  if (xcnr == xcFunctional%HF_Exchange) then
    call hfex(kk, max_l, num_alpha, alpha, poly_order, problemsize)
  elseif (xcFunctional%isLongRangeCorrected(xcnr)) then
    call hfex_lr(kk_lr, max_l, num_alpha, alpha, poly_order, problemsize, omega, grid_params)
  elseif (xcFunctional%isGlobalHybrid(xcnr)) then
    call hfex(kk, max_l, num_alpha, alpha, poly_order, problemsize)
  elseif (xcFunctional%isCAMY(xcnr)) then
    call hfex(kk, max_l, num_alpha, alpha, poly_order, problemsize)
    call hfex_lr(kk_lr, max_l, num_alpha, alpha, poly_order, problemsize, omega, grid_params)
  end if

  ! convergence flag
  tConverged = .false.

  ! DFT start potential
  if (.not. (xcnr == xcFunctional%HF_Exchange)) then
    call dft_start_pot(abcissa, num_mesh_points, nuc, vxc)
  end if

  ! build initial fock matrix, core hamiltonian only
  write(*, '(A)') 'Startup: Building Initial Fock Matrix'
  write(*, '(A)') ' '

  ! do not confuse mixer
  pot_old(:,:,:,:) = 0.0_dp

  ! kinetic energy, nuclear-electron, and confinement matrix elements which are constant during SCF
  call build_hamiltonian(0, tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha, poly_order,&
      & problemsize, xcnr, num_mesh_points, weight, abcissa, vxc, alpha, pot_old, pot_new, tZora,&
      & tBroyden, mixing_factor, ff, camAlpha, camBeta)

  ! self-consistency cycles
  write(*,*) 'Energies in Hartree'
  write(*,*)
  write(*,*) ' Iter |   Total energy  |   HF-X energy  |   XC energy   |   Change in pot'
  write(*,*) '--------------------------------------------------------------------------'
  lpScf: do iScf = 1, maxiter

    pot_old(:,:,:,:) = pot_new

    ! diagonalize
    call diagonalize(max_l, num_alpha, poly_order, ff, ss, cof, eigval)

    ! build density matrix
    call densmatrix(problemsize, max_l, occ, cof, pp)

    ! get electron density, derivatives, exc related potentials and energy densities
    call density_grid(pp, max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa, dzdr,&
        & dz, xcnr, omega, camAlpha, camBeta, rho, drho, ddrho, vxc, exc, xalpha_const)

    ! build Fock matrix and get total energy during SCF
    call build_hamiltonian(iScf, tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha,&
        & poly_order, problemsize, xcnr, num_mesh_points, weight, abcissa, vxc, alpha, pot_old,&
        & pot_new, tZora, tBroyden, mixing_factor, ff, camAlpha, camBeta)

    if (tZora) then
      call getTotalEnergyZora(tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha, poly_order,&
          & problemsize, xcnr, num_mesh_points, weight, abcissa, rho, exc, vxc, eigval_scaled, occ,&
          & camAlpha, camBeta, kinetic_energy, nuclear_energy, coulomb_energy, exchange_energy,&
          & x_en_2, conf_energy, total_ene)
    else
      call getTotalEnergy(tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha, poly_order,&
          & xcnr, num_mesh_points, weight, abcissa, rho, exc, camAlpha, camBeta, kinetic_energy,&
          & nuclear_energy, coulomb_energy, exchange_energy, x_en_2, conf_energy, total_ene)
    end if

    call check_convergence(pot_old, pot_new, max_l, problemsize, iScf, change_max, tConverged)

    write(*, '(I4,2X,3(1X,F16.9),3X,E16.9)') iScf, total_ene, exchange_energy, x_en_2, change_max

    ! if self-consistency is reached, exit loop
    if (tConverged) exit lpScf

    ! check conservation of number of electrons during SCF
    call check_electron_number(cof, ss, occ, max_l, num_alpha, poly_order, problemsize)

    write(*,*) ' '

  end do lpScf

  ! handle non-converged calculations
  if (.not. tConverged) then
    call error('SCC is NOT converged, maximal SCC iterations exceeded.')
  end if

  ! handle output of requested data

  if (tPrintEigvecs) then
    call write_eigvec(max_l, num_alpha, alpha, poly_order, eigval, cof)
    call write_moments(max_l, num_alpha, alpha, poly_order, problemsize, cof)
    call cusp_values(max_l, cof, pp, alpha, num_alpha, poly_order)
  end if

  call write_eigval(max_l, num_alpha, poly_order, eigval)
  call write_energies(kinetic_energy, nuclear_energy, coulomb_energy, exchange_energy, x_en_2,&
      & conf_energy, total_ene, tZora)

  if (tZora) then
    call scaled_zora(eigval, max_l, num_alpha, alpha, poly_order, problemsize, num_mesh_points,&
        & weight, abcissa, vxc, nuc, pp, tt, cof, occ, eigval_scaled, zora_ekin)

    write(*, '(A)') 'Scaled Scalar-Relativistic ZORA EIGENVALUES and ENERGY'
    write(*, '(A)') '------------------------------------------------------'
    call write_eigval(max_l, num_alpha, poly_order, eigval_scaled)
  end if

  if (tZora) then
    call getTotalEnergyZora(tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha, poly_order,&
        & problemsize, xcnr, num_mesh_points, weight, abcissa, rho, exc, vxc, eigval_scaled, occ,&
        & camAlpha, camBeta, zora_ekin, nuclear_energy, coulomb_energy, exchange_energy, x_en_2,&
        & conf_energy, total_ene)
  end if

  write(*, '(A,E20.12)') 'Potential Matrix Elements converged to ', change_max
  write(*, '(A)') ' '

  if (tZora) then
    call write_energies_tagged(zora_ekin, nuclear_energy, coulomb_energy, exchange_energy,&
        & conf_energy, 0.0_dp, tZora, eigval_scaled, occ)
  else
    call write_energies_tagged(kinetic_energy, nuclear_energy, coulomb_energy, x_en_2, conf_energy,&
        & total_ene, tZora, eigval, occ)
  end if

  call write_potentials_file_standard(num_mesh_points, abcissa, weight, vxc, rho, nuc, pp, max_l,&
      & num_alpha, poly_order, alpha, problemsize)

  call write_densities_file_standard(num_mesh_points, abcissa, weight, rho, drho, ddrho)

  ! write wave functions and eventually invert to have positive starting gradient
  call write_waves_file_standard(num_mesh_points, abcissa, weight, alpha, num_alpha, poly_order,&
      & max_l, problemsize, occ, qnvalorbs, cof)

  call write_wave_coeffs_file(max_l, num_alpha, poly_order, cof, alpha, occ, qnvalorbs)

end program HFAtom
