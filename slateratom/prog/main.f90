program HFAtom

  use common_accuracy, only : dp
  use integration, only : gauss_chebyshev_becke_mesh
  use input, only : read_input_1, read_input_2, echo_input
  use core_overlap, only : overlap, nuclear, kinetic, confinement
  use coulomb_hfex, only : coulomb, hfex
  use densitymatrix, only : densmatrix
  use hamiltonian, only : build_fock
  use diagonalizations, only : diagonalize, diagonalize_overlap
  use output, only : write_eigvec, write_eigval, write_moments, write_energies,&
      & write_energies_tagged, write_potentials_file_standard, write_densities_file_standard,&
      & write_waves_file_standard, write_wave_coeffs_file, cusp_values
  use totalenergy, only : total_energy
  use dft, only : check_accuracy, dft_start_pot, density_grid
  use utilities, only : check_electron_number, check_convergence
  use zora_routines, only : scaled_zora
  use cmdargs, only : parse_command_arguments
  use globals

  implicit none


  !> current SCF step
  integer :: iScf

  !> quantum numbers of wavefunctions to be written
  integer, allocatable :: qnvalorbs(:,:)


  call parse_command_arguments()
  call read_input_1(nuc, max_l, occ_shells, maxiter, poly_order, min_alpha, max_alpha, num_alpha,&
      & tAutoAlphas, alpha, conf_r0, conf_power, num_occ, num_power, num_alphas, xcnr,&
      & tPrintEigvecs, tZora, tBroyden, mixing_factor, xalpha_const)

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

  call echo_input(nuc, max_l, occ_shells, poly_order, num_alpha, alpha, conf_r0, conf_power, occ,&
      & num_occ, num_power, num_alphas, xcnr, tZora, num_mesh_points, xalpha_const)

  ! allocate global stuff and zero out
  call allocate_globals()

  ! generate radial integration mesh
  call gauss_chebyshev_becke_mesh(num_mesh_points, nuc, weight, abcissa, dzdr, d2zdr2, dz)

  ! check mesh accuracy
  call check_accuracy(weight, abcissa, num_mesh_points, max_l, num_alpha, alpha, poly_order)

  if (xcnr >= 2) then
    write(*, "(A,/)") "LDA/PBE ROUTINES: LIBXC IMPLEMENTATION"
  end if

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
  if (xcnr == 0) call hfex(kk, max_l, num_alpha, alpha, poly_order, problemsize)

  ! convergence flag
  tConverged = .false.

  ! dft start potential
  if (xcnr > 0) call dft_start_pot(abcissa, num_mesh_points, nuc, vxc)

  ! build initial fock matrix, core hamiltonian only
  write(*, '(A)') 'Startup: Building Initial Fock Matrix'
  write(*, '(A)') ' '

  ! do not confuse mixer
  pot_old(:,:,:,:) = 0.0_dp

  ! kinetic energy, nuclear-electron, and confinement matrix elements which are constant during SCF
  call build_fock(0, tt, uu, nuc, vconf, jj, kk, pp, max_l, num_alpha, poly_order, problemsize,&
      & xcnr, num_mesh_points, weight, abcissa, vxc, alpha, pot_old, pot_new, tZora, tBroyden,&
      & mixing_factor, ff)

  lpScf: do iScf = 1, maxiter

    write(*, '(A,I5)') 'Iteration :', iScf

    pot_old = pot_new

    ! diagonalize
    call diagonalize(max_l, num_alpha, poly_order, ff, ss, cof, eigval)

    ! build density matrix
    call densmatrix(problemsize, max_l, occ, cof, pp)

    ! get electron density, derivatives, exc related potentials and energy densities
    call density_grid(pp, max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa, dzdr,&
        & dz, xcnr, rho, drho, ddrho, vxc, exc, xalpha_const)

    ! build Fock matrix and get total energy during SCF
    call build_fock(iScf, tt, uu, nuc, vconf, jj, kk, pp, max_l, num_alpha, poly_order,&
        & problemsize, xcnr, num_mesh_points, weight, abcissa, vxc, alpha, pot_old, pot_new, tZora,&
        & tBroyden, mixing_factor, ff)

    call total_energy(tt, uu, nuc, vconf, jj, kk, pp, max_l, num_alpha, poly_order, problemsize,&
        & xcnr, num_mesh_points, weight, abcissa, rho, exc, kinetic_energy, nuclear_energy,&
        & coulomb_energy, exchange_energy, conf_energy, total_ene)

    if (.not. tZora) then
      ! non-rel. total energy during SCF meaningless for ZORA
      ! but energy contributions needed once SCF converged, so surpress output
      write(*, '(A,F18.6,A)') 'TOTAL ENERGY', total_ene, ' Hartree'
    end if

    call check_convergence(pot_old, pot_new, max_l, problemsize, iScf, change_max, tConverged)

    write(*, '(A,E20.12)') 'CHANGE in potential matrix', change_max

    ! if self-consistency is reached, exit loop
    if (tConverged) exit

    ! check conservation of number of electrons during SCF
    call check_electron_number(cof, ss, occ, max_l, num_alpha, poly_order, problemsize)

    write(*, '(A)') ' '

  end do lpScf

  ! handle output of requested data

  if (tPrintEigvecs) then
    call write_eigvec(max_l, num_alpha, alpha, poly_order, eigval, cof)
    call write_moments(max_l, num_alpha, alpha, poly_order, problemsize, cof)
    call cusp_values(max_l, cof, pp, alpha, num_alpha, poly_order)
  end if

  call write_eigval(max_l, num_alpha, poly_order, eigval)
  call write_energies(kinetic_energy, nuclear_energy, coulomb_energy, exchange_energy, conf_energy,&
      & total_ene, .false.)

  if (tZora) then
    call scaled_zora(eigval, max_l, num_alpha, alpha, poly_order, problemsize, num_mesh_points,&
        & weight, abcissa, vxc, nuc, pp, tt, cof, occ, eigval_scaled, zora_ekin)

    write(*, '(A)') 'Scaled Scalar-Relativistic ZORA EIGENVALUES and ENERGY'
    write(*, '(A)') '------------------------------------------------------'
    call write_eigval(max_l, num_alpha, poly_order, eigval_scaled)
    call write_energies(zora_ekin, nuclear_energy, coulomb_energy, exchange_energy, conf_energy,&
        & total_ene, .true.)
  end if

  write(*, '(A,E20.12)') 'Potential Matrix Elements converged to ', change_max
  write(*, '(A)') ' '

  if (tZora) then
    call write_energies_tagged(zora_ekin, nuclear_energy, coulomb_energy, exchange_energy,&
        & conf_energy, 0.0_dp, tZora, eigval_scaled, occ)
  else
    call write_energies_tagged(kinetic_energy, nuclear_energy, coulomb_energy, exchange_energy,&
        & conf_energy, total_ene, tZora, eigval, occ)
  end if

  call write_potentials_file_standard(num_mesh_points, abcissa, weight, vxc, rho, nuc, pp, max_l,&
      & num_alpha, poly_order, alpha, problemsize)

  call write_densities_file_standard(num_mesh_points, abcissa, weight, rho, drho, ddrho)

  ! write wave functions and eventually invert to have positive starting gradient
  call write_waves_file_standard(num_mesh_points, abcissa, weight, alpha, num_alpha, poly_order,&
      & max_l, problemsize, occ, qnvalorbs, cof)

  call write_wave_coeffs_file(max_l, num_alpha, poly_order, cof, alpha, occ, qnvalorbs)

end program HFAtom
