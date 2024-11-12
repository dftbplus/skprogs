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
      & write_waves_file_standard, write_wave_coeffs_file, cusp_values, writeAveragePotential
  use sap, only : sap_start_pot
  use totalenergy, only : getTotalEnergy, getTotalEnergyZora
  use dft, only : check_accuracy, density_grid
  use utilities, only : check_electron_number, check_convergence_energy,&
      & check_convergence_commutator, check_convergence_eigenspectrum,&
      & compute_commutator, check_convergence_pot
  use zora_routines, only : scaled_zora
  use cmdargs, only : parse_command_arguments
  use common_poisson, only : TBeckeGridParams
  use xcfunctionals, only : xcFunctional
  use average_potential, only : getAveragePotential
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

  !! Kinetic energy reference for average potential calculation
  real(dp), allocatable :: kinetic_energy_ref

  !! holds parameters, defining a Becke integration grid
  type(TBeckeGridParams) :: grid_params

  ! deactivate average potential calculation for now
  isAvgPotNeeded = .false.

  call parse_command_arguments()
  call read_input_1(nuc, max_l, occ_shells, maxiter, scftol, poly_order, min_alpha, max_alpha,&
      & num_alpha, tAutoAlphas, alpha, conf_r0, conf_power, num_occ, num_power, num_alphas, xcnr,&
      & tPrintEigvecs, tZora, mixnr, mixing_factor, xalpha_const, omega, camAlpha, camBeta,&
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

  ! meta-GGA functionals sometimes exhibit extraordinary slow
  ! convergence w.r.t the number of radial grid points
  ! see 10.1063/5.0121187
  ! WARNING: too high number of grid points somehow
  ! manages to break the Broyden mixer!
  if (xcFunctional%isMGGA(xcnr)) then
    ! num_mesh_points = num_mesh_points + 2000
    num_mesh_points = num_mesh_points + 500
  end if

  call echo_input(nuc, max_l, occ_shells, maxiter, scftol, poly_order, num_alpha, alpha, conf_r0,&
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

  ! convergence flags
  tCommutatorConverged = .false.
  tEnergyConverged = .false.
  tEigenspectrumConverged = .false.

  ! Generate guess for DFT;
  ! Thomas-Fermi guess potential is currently disabled
  if (.not. (xcnr == xcFunctional%HF_Exchange)) then
    ! SAP potential
    call sap_start_pot(abcissa, num_mesh_points, nuc, vxc)
  end if

  ! build initial fock matrix, core hamiltonian only
  write(*, '(A)') 'Startup: Building Initial Fock Matrix'
  write(*, '(A)') ' '

  ! do not confuse mixer
  pot_old(:,:,:,:) = 0.0_dp

  ! kinetic energy, nuclear-electron, and confinement matrix elements which are constant during SCF
  call build_hamiltonian(pMixer, 0, tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha,&
      & poly_order, problemsize, xcnr, num_mesh_points, weight, abcissa, vxc, vtau, alpha, pot_old,&
      & pot_new, tZora, ff, commutator, camAlpha, camBeta)

  ! self-consistency cycles
  write(*,*) 'Energies in Hartree'
  write(*,*)
  write(*,*) ' Iter |   Total energy  |   HF-X energy  |   XC energy   |   max(abs([F,PS]))   &
      & | Delta (Total energy) | Delta (spectrum)'
  write(*,*) '-------------------------------------------------------------------------------------&
      & ------------------------------------------'
  lpScf: do iScf = 1, maxiter

    pot_old(:,:,:,:) = pot_new
    total_ene_old = total_ene
    eigval_old(:,:,:) = eigval

    ! diagonalize
    call diagonalize(max_l, num_alpha, poly_order, ff, ss, cof, eigval)

    ! build density matrix
    call densmatrix(problemsize, max_l, occ, cof, pp)

    ! get electron density, derivatives, exc related potentials and energy densities
    call density_grid(pp, max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa, dzdr,&
        & dz, xcnr, omega, camAlpha, camBeta, rho, drho, ddrho, tau, vxc, vtau, exc, xalpha_const)

    ! build Fock matrix and get total energy during SCF
    call build_hamiltonian(pMixer, iScf, tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha,&
        & poly_order, problemsize, xcnr, num_mesh_points, weight, abcissa, vxc, vtau, alpha,&
        & pot_old, pot_new, tZora, ff, commutator, camAlpha, camBeta)

    ! compute [F,PS]
    call compute_commutator(max_l, num_alpha, poly_order, ff, pp, ss, commutator)

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

    call check_convergence_commutator(commutator, scftol, iScf, commutator_max,&
        & tCommutatorConverged)
    ! For debugging:
    ! call check_convergence_pot(pot_old, pot_new, max_l, problemsize, scftol, iScf,&
        ! & commutator_max, tCommutatorConverged)
    call check_convergence_energy(total_ene_old, total_ene, scftol, iScf, total_ene_diff,&
        & tEnergyConverged)
    call check_convergence_eigenspectrum(max_l, num_alpha, poly_order, eigval, eigval_old, occ,&
        & scftol, iScf, eigval_diff, tEigenspectrumConverged)

    ! Print SCF loop information
    write(*, '(I4,2X,3(1X,F16.9),7X,E16.9,8X,E16.9,8X,E16.9)') iScf, total_ene, exchange_energy, x_en_2,&
        & commutator_max, total_ene_diff, eigval_diff

    ! if self-consistency is reached, exit loop
    if (tCommutatorConverged .and. tEnergyConverged .and. tEigenspectrumConverged) exit lpScf

    ! check conservation of number of electrons during SCF
    call check_electron_number(cof, ss, occ, max_l, num_alpha, poly_order, problemsize)

    write(*,*) ' '

  end do lpScf

  ! handle non-converged calculations
  if (.not. (tEnergyConverged .and. tCommutatorConverged .and. tEigenspectrumConverged)) then
    call error('SCF is NOT converged, maximal SCF iterations exceeded.')
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

  write(*, '(A,E20.12)') '[F,PS] converged to', commutator_max
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

  call write_densities_file_standard(xcnr, num_mesh_points, abcissa, weight, rho, drho, ddrho, tau)

  ! write wave functions and eventually invert to have positive starting gradient
  call write_waves_file_standard(num_mesh_points, abcissa, weight, alpha, num_alpha, poly_order,&
      & max_l, problemsize, occ, qnvalorbs, cof)

  call write_wave_coeffs_file(max_l, num_alpha, poly_order, cof, alpha, occ, qnvalorbs)

  if (isAvgPotNeeded) then
    if (.not. tZora) kinetic_energy_ref = kinetic_energy
    call getAveragePotential(cof, eigval, occ, abcissa, weight, max_l, num_alpha, alpha,&
        & poly_order, problemsize, scftol, maxiter, avgPot, kinetic_energy_ref=kinetic_energy_ref)
    call writeAveragePotential(abcissa, avgPot)
  end if

end program HFAtom
