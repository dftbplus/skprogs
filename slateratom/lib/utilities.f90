!> Module that provides common utilities used by the slateratom program.
module utilities

  use common_accuracy, only : dp

  implicit none
  private

  public :: check_electron_number, check_convergence_energy, check_convergence_commutator,&
      & check_convergence_eigenspectrum, compute_commutator, check_convergence_pot,&
      & check_convergence_density
  public :: vector_length, fak, zeroOutCpotOfEmptyDensitySpinChannels


contains

  pure subroutine zeroOutCpotOfEmptyDensitySpinChannels(rho, vc)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> correlation potential on grid, shape: (nSpin, nGridPoints)
    real(dp), intent(inout) :: vc(:,:)

    !! maximum density in up/down channels
    real(dp) :: maxSpinUp, maxSpinDn

    maxSpinUp = maxval(abs(rho(:, 1)))
    maxSpinDn = maxval(abs(rho(:, 2)))

    if (maxSpinUp < 1e-16_dp) vc(1, :) = 0.0_dp
    if (maxSpinDn < 1e-16_dp) vc(2, :) = 0.0_dp

  end subroutine zeroOutCpotOfEmptyDensitySpinChannels


  !> Checks SCF convergence by comparing new and old potential.
  pure subroutine check_convergence_pot(pot_old, pot_new, max_l, problemsize, scftol, iScf,&
      & change_max, tConverged)

    !> old and new potential to compare
    real(dp), intent(in) :: pot_old(:,0:,:,:), pot_new(:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> size of the problem at hand
    integer, intent(in) :: problemsize

    !> scf tolerance, i.e. convergence criteria
    real(dp), intent(in) :: scftol

    !> current SCF step
    integer, intent(in) :: iScf

    !> obtained maximum change in potential
    real(dp), intent(out) :: change_max

    !> true, if SCF converged
    logical, intent(out) :: tConverged

    !> auxiliary variables
    integer ii, jj, kk, ll

    change_max = 0.0_dp

    do ii = 1, 2
      do jj = 0, max_l
        do kk = 1, problemsize
          do ll = 1, problemsize
            change_max = max(change_max, abs(pot_old(ii, jj, kk, ll) - pot_new(ii, jj, kk, ll)))
          end do
        end do
      end do
    end do

    if (change_max < scftol) then
      tConverged = .true.
    else
      tConverged = .false.
    end if

    if (iScf < 3) then
      tConverged = .false.
    end if

  end subroutine check_convergence_pot


  !> Check convergence by finding the maximum absolute value of the orbital commutator [F,PS]
  pure subroutine check_convergence_commutator(commutator, scftol, iScf, commutator_max,&
      & tConverged)

    !> commutator [F, PS]
    real(dp), intent(in) :: commutator(:, 0:, :, :)

    !> scf tolerance, i.e. convergence criteria
    real(dp), intent(in) :: scftol

    !> current SCF step
    integer, intent(in) :: iScf

    !> orbital gradient norm value
    real(dp), intent(out) :: commutator_max

    !> true, if SCF converged
    logical, intent(out) :: tConverged

    commutator_max = maxval(abs(commutator))

    tConverged = commutator_max < scftol

    if (iScf < 3) then
      tConverged = .false.
    end if

  end subroutine check_convergence_commutator


  !> Compute the commutator [F,PS] in MO basis  
  subroutine compute_commutator(max_l, num_alpha, poly_order, ff, pp, ss, commutator)
    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> Fock matrix
    real(dp), intent(in) :: ff(:, 0:, :, :)

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

    !> overlap supervector
    real(dp), intent(in) :: ss(0:,:,:)

    !> commutator [F,PS]
    real(dp), intent(out) :: commutator(:, 0:, :, :)

    !! auxilliary variables
    integer :: iSpin, ll, diagsize
    real(dp), allocatable :: fps(:,:), spf(:,:)

    commutator = 0.0_dp

    diagsize = size(ff, dim=4)
    allocate(fps(diagsize, diagsize))
    allocate(spf(diagsize, diagsize))

    do iSpin = 1, 2
      do ll = 0, max_l
        ! Compute [F,PS]= FPS - SPF
        fps = matmul(ff(iSpin, ll, :, :), matmul(pp(iSpin, ll, :, :), ss(ll, :, :)))
        spf = matmul(matmul(ss(ll, :, :), pp(iSpin, ll, :, :)), ff(iSpin, ll, :, :))
        commutator(iSpin, ll, :, :) = fps - spf
      end do
    end do

    deallocate(fps)
    deallocate(spf)

  end subroutine 


  !> Checks convergence by computing max. density change from the last SCF iteration.
  pure subroutine check_convergence_density(pp_old, pp_new, scftol, iScf, res, tConverged)

    !> old and new energy to compare
    real(dp), intent(in) :: pp_old(:,:,:,:), pp_new(:,:,:,:)

    !> scf tolerance, i.e. convergence criteria
    real(dp), intent(in) :: scftol

    !> current SCF step
    integer, intent(in) :: iScf

    !> obtained maximum change in energy
    real(dp), intent(out) :: res

    !> true, if SCF converged
    logical, intent(out) :: tConverged

    res = 0.0_dp

    res = maxval(abs(pp_new - pp_old))
    tConverged = res < scftol

    if (iScf < 3) then
      tConverged = .false.
    end if

  end subroutine check_convergence_density


    !> Checks convergence by computing energy change from the last SCF iteration.
  pure subroutine check_convergence_energy(energy_old, energy_new, scftol, iScf, change,&
      & tConverged)

    !> old and new energy to compare
    real(dp), intent(in) :: energy_new, energy_old

    !> scf tolerance, i.e. convergence criteria
    real(dp), intent(in) :: scftol

    !> current SCF step
    integer, intent(in) :: iScf

    !> obtained maximum change in energy
    real(dp), intent(out) :: change

    !> true, if SCF converged
    logical, intent(out) :: tConverged

    change = 0.0_dp

    change = energy_new - energy_old
    tConverged = abs(change) < scftol

    if (iScf < 3) then
      tConverged = .false.
    end if

  end subroutine check_convergence_energy


  !> Checks convergence by evaluating change in the occupied part of the eigenspectrum
  pure subroutine check_convergence_eigenspectrum(max_l, num_alpha, poly_order, eigval_new,&
      & eigval_old, occ, scftol, iScf, res, tConverged)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)
  
    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> old and new eigenspectra to compare
    real(dp), intent(in) :: eigval_new(:,0:,:), eigval_old(:,0:,:)

    !> occupations
    real(dp), intent(in) :: occ(:,0:,:)

    !> scf tolerance, i.e. convergence criteria
    real(dp), intent(in) :: scftol

    !> current SCF step
    integer, intent(in) :: iScf

    !> obtained change
    real(dp), intent(out) :: res

    !> true, if SCF converged
    logical, intent(out) :: tConverged

    ! Loop indices
    integer :: iSpin, ll, diagsize, ii

    ! Occupation
    real(dp) :: occ_ii

    ! Difference
    real(dp) :: e_diff

    ! Max. difference, only occupied orbitals contribute
    res = 0.0_dp
    do iSpin = 1, 2
      do ll = 0, max_l
        diagsize = num_alpha(ll) * poly_order(ll)
        do ii = 1, diagsize
          occ_ii = abs(occ(iSpin, ll, ii))
          if (occ_ii > 1e-16) then
            ! change = change + (eigval_new(iSpin, ll, ii) - eigval_old(iSpin, ll, ii))**2 
            e_diff = abs(eigval_new(iSpin, ll, ii) - eigval_old(iSpin, ll, ii)) 
            if (e_diff > res) then
              res = e_diff
            end if
          end if
        end do
      end do
    end do

    tConverged = res < scftol

    if (iScf < 3) then
      tConverged = .false.
    end if

  end subroutine check_convergence_eigenspectrum


  !> Checks conservation of electron number during SCF. If this fluctuates you are in deep trouble.
  subroutine check_electron_number(cof, ss, occ, max_l, num_alpha, poly_order, problemsize)

    !> wavefunction coefficients
    real(dp), intent(in) :: cof(:,0:,:,:)

    !> overlap supervector
    real(dp), intent(in) :: ss(0:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occ(:,0:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> size of the problem at hand
    integer, intent(in) :: problemsize

    !> actual number of electrons during SCF
    real(dp) :: electron_number

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, pp, qq

    ! get actual number per shell by multiplying dens matrix and overlap
    do mm = 1, 2
      do ii = 0, max_l
        do qq = 1, problemsize
          electron_number = 0.0_dp
          ll = 0
          do jj = 1, num_alpha(ii)
            do kk = 1, poly_order(ii)
              ll = ll + 1
              pp = 0
              do nn = 1, num_alpha(ii)
                do oo = 1, poly_order(ii)
                  pp = pp + 1

                  electron_number = electron_number + occ(mm, ii, qq)&
                      & * cof(mm, ii, ll, qq) * cof(mm, ii, pp, qq) * ss(ii, ll, pp)

                end do
              end do
            end do
          end do

          if (abs(occ(mm, ii, qq) - electron_number) > 1.0e-08_dp) then
            write(*,*) 'Electron number fluctuation', occ(mm, ii, qq) - electron_number
          end if

        end do
      end do
    end do

  end subroutine check_electron_number


  !> Calculates the Euclidean vector norm.
  pure function vector_length(vector) result(norm)

    !> vector to calculate the Euclidean norm for
    real(dp), intent(in) :: vector(:)

    !> norm of vector
    real(dp) :: norm

    norm = norm2(vector)

  end function vector_length


  !> Calculates the factorial of nn.
  pure function fak(nn)

    !> integer argument to calculate factorial for
    integer, intent(in) :: nn

    !> resulting factorial
    real(dp) :: fak

    !> auxiliary variable
    integer :: ii

    fak = 1.0_dp

    do ii = 1, nn
      fak = fak * real(ii, dp)
    end do

  end function fak

end module utilities
