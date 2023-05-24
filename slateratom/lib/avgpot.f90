!> Module to calculate the approximate average potential and effective orbital energies proposed in:
!! Baer R, Livshits E, Salzner U. Tuned range-separated hybrids in density functional theory.
!! Annual review of physical chemistry. 2010 May 5;61:85-109.
module average_potential

  use common_accuracy, only : dp, mc
  use common_constants, only : pi
  use common_message, only : error

  use utilities, only : index_heap_sort
  use density, only : wavefunction, wavefunction_1st, wavefunction_2nd

  implicit none
  private

  public :: getAveragePotential


contains


  !> Tries to infer the energy of the highest occupied atomic orbital (HOAO) and the associated
  !! principal and angular quantum number from the eigenvalues and occupations handed over.
  !!
  !! The effort made here would be an argument for moving this functionality to the sktools
  !! sometime in the future.
  subroutine getHoaoOrLowestNl(eigval, occ, max_l, num_alpha, poly_order, hoaoN, hoaoL, eHoao)

    !> Eigenvalues of selected spin-channel
    real(dp), intent(in) :: eigval(0:, :)

    !> Occupation numbers of selected spin-channel
    real(dp), intent(in) :: occ(0:, :)

    !> Maximum angular momentum
    integer, intent(in) :: max_l

    !> Number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> Highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> Principal and angular quantum number of HOAO eigenvalue
    integer, intent(out) :: hoaoN, hoaoL

    !> Energy of highest occupied atomic orbital
    real(dp), intent(out) :: eHoao

    !! Flattened, one-dimensional versions of occupations and eigenvalues
    real(dp), allocatable :: occFlat(:), eigvalFlat(:)

    !! Index array that sorts the flattened eigenvalues
    integer, allocatable :: eigvalFlatSortInds(:)

    !! True, if all occupations of current spin channel are below threshold
    logical :: isUnoccupied

    !! Auxiliary variables
    integer :: ll, nn, iDummy1, iDummy2, iHoao, ii, ind

    eigvalFlat = reshape(eigval, [size(eigval)])
    occFlat = reshape(occ, [size(occ)])

    allocate(eigvalFlatSortInds(size(eigvalFlat)))
    call index_heap_sort(eigvalFlatSortInds, eigvalFlat)

    ! reverse order, since we would like to iterate over the eigenvalues from unoccupied to occupied
    eigvalFlatSortInds(:) = eigvalFlatSortInds(size(eigvalFlatSortInds):1:-1)

    isUnoccupied = all(occFlat < 1.0e-08_dp)

    if (isUnoccupied) then
      iHoao = eigvalFlatSortInds(1)
    else
      lpIndx: do ii = 1, size(eigvalFlatSortInds)
        ind = eigvalFlatSortInds(ii)
        if (occFlat(ind) >= 1.0e-08_dp) then
          iHoao = ind
          exit lpIndx
        end if
      end do lpIndx
    end if

    eHoao = eigvalFlat(iHoao)

    lpAng: do ll = 0, max_l
      nn = 0
      do iDummy1 = 1, num_alpha(ll)
        do iDummy2 = 1, poly_order(ll)
          nn = nn + 1

          if (abs(eigval(ll, nn) - eHoao) < 1e-10_dp) then
            hoaoN = nn
            hoaoL = ll
            exit lpAng
          end if

        end do
      end do
    end do lpAng

  end subroutine getHoaoOrLowestNl


  !> Determines an approximate, average, local potential and associated effective orbital energies
  !! by minimizing the deviance from local Schrödinger equations.
  !!
  !! See the following reference for full details:
  !!
  !! Roi Baer, Ester Livshits, and Ulrike Salzner.
  !! “Tuned Range-Separated Hybrids in Density Functional Theory”.
  !! In: Annu. Rev. Phys. Chem. 61.1 (2010), pp. 85–109.
  !! DOI: 10.1146/annurev.physchem.012809.103321
  subroutine getAveragePotential(cof, eigval, occ, abcissa, weights, max_l, num_alpha, alpha,&
      & poly_order, problemsize, scftol, maxiter, avgPot, kinetic_energy_ref)

    !> wavefunction coefficients
    real(dp), intent(in) :: cof(:, 0:, :,:)

    !> eigenvalues
    real(dp), intent(in) :: eigval(:, 0:, :)

    !> occupation numbers
    real(dp), intent(in) :: occ(:, 0:, :)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> numerical integration weights
    real(dp), intent(in) :: weights(:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> Tolerance for self-consistency
    real(dp), intent(in) :: scftol

    !> Maximum number of self-consistency iterations
    integer, intent(in) :: maxiter

    !> Average local, effective potential
    real(dp), intent(out) :: avgPot(:,:)

    !> Reference kinetic energy
    real(dp), intent(in), optional :: kinetic_energy_ref

    !! HOAO shift required to resolve the ambiguity caused by the common constant
    real(dp) :: hoaoShift(2)

    !! Maximum change in effective orbital energies
    real(dp) :: change_max

    !! Wave function and 2nd derivative
    real(dp), allocatable :: rad(:,:,:,:), radp(:,:,:,:), radpp(:,:,:,:)

    !! Effective orbital energies of current and previous SC iteration
    real(dp), allocatable :: eps(:,:,:), epsLast(:,:,:)

    !! Electron density
    real(dp), allocatable :: rho(:,:)

    !! Self-consitent iteration
    integer :: iSC

    !! Iterates over radial grid points
    integer :: iRad

    !! Principal quantum number
    integer :: nn

    !! Angular quantum number
    integer :: ll

    !! Dummy loop indices
    integer :: iDummy1, iDummy2

    !! Number of radial integration points
    integer :: nRadial

    !! Principal and angular quantum number of HOAO
    integer :: hoaoN(2), hoaoL(2)

    !! HOAO eigenvalues
    real(dp) :: eHoao(2)

    !! Spin index and number of spin channels
    integer :: iSpin, nSpin

    !! True, if SC loop reached convergency
    logical :: tConverged

    !! Number of electrons by integrating up density
    real(dp) :: electron_number

    !! Kinetic energy is evaluated with quantities from this module
    real(dp) :: kinetic_energy

    tConverged = .false.

    nRadial = size(abcissa)
    nSpin = size(occ, dim=1)

    do iSpin = 1, nSpin
      call getHoaoOrLowestNl(eigval(iSpin, 0:, :), occ(iSpin, 0:, :), max_l, num_alpha,&
          & poly_order, hoaoN(iSpin), hoaoL(iSpin), eHoao(iSpin))
    end do

    allocate(rad(nRadial, problemsize, 0:max_l, nSpin), source=0.0_dp)
    allocate(radp(nRadial, problemsize, 0:max_l, nSpin), source=0.0_dp)
    allocate(radpp(nRadial, problemsize, 0:max_l, nSpin), source=0.0_dp)

    allocate(eps, mold=eigval)
    eps(:, 0:, :) = 0.0_dp
    allocate(epsLast, mold=eigval)
    epsLast(:, 0:, :) = 0.0_dp

    allocate(rho(nRadial, nSpin), source=0.0_dp)

    ! Build radial orbitals on Becke's Gauss-Chebyschev mesh
    do iSpin = 1, nSpin
      do ll = 0, max_l
        nn = 0
        do iDummy1 = 1, num_alpha(ll)
          do iDummy2 = 1, poly_order(ll)
            nn = nn + 1

            do iRad = 1, nRadial
              rad(iRad, nn, ll, iSpin) = wavefunction(cof(iSpin, ll, :, nn), alpha, num_alpha,&
                  & poly_order, ll, abcissa(iRad))
              radp(iRad, nn, ll, iSpin) = wavefunction_1st(cof(iSpin, ll, :, nn), alpha, num_alpha,&
                  & poly_order, ll, abcissa(iRad))
              radpp(iRad, nn, ll, iSpin) = wavefunction_2nd(cof(iSpin, ll, :, nn), alpha,&
                  & num_alpha, poly_order, ll, abcissa(iRad))
            end do

            ! Electron density as occupation-weighted absolute squares of the orbitals
            rho(:, iSpin) = rho(:, iSpin) + occ(iSpin, ll, nn) * rad(:, nn, ll, iSpin)**2

          end do
        end do
      end do
    end do

    ! Check if established density integrates up to the correct number of electrons
    electron_number = sum(weights * (rho(:, 1) + rho(:, 2)) * abcissa**2)
    if (abs(sum(occ) - electron_number) > 1.0e-08_dp) then
      call error("Average-potential: Mismatch in number of electrons.")
    end if

    ! Optionally check if the original kinetic energy is reproduced by the local quantities
    if (present(kinetic_energy_ref)) then
      kinetic_energy = getKineticEnergy(occ, abcissa, weights, max_l, num_alpha, poly_order, rad,&
          & radp, radpp)
      if (abs(kinetic_energy - kinetic_energy_ref) > 1.0e-06_dp) then
        call error("Average-potential: Mismatch in kinetic energy.")
      end if
    end if

    ! As an initial guess for the effective orbital energies we may just use the eigenvalues from
    ! the converged GKS calculation:
    eps(:, 0:, :) = eigval(:, 0:, :)
    epsLast(:, 0:, :) = eps(:, 0:, :)

    ! Start self-consistency iterations
    lpSC: do iSC = 1, maxiter

      ! Build average potential
      avgPot(:,:) = 0.0_dp
      do iSpin = 1, nSpin
        do ll = 0, max_l
          nn = 0
          do iDummy1 = 1, num_alpha(ll)
            do iDummy2 = 1, poly_order(ll)
              nn = nn + 1

              ! uses +laplace / 2
              avgPot(:, iSpin) = avgPot(:, iSpin) + occ(iSpin, ll, nn)&
                  & * (epsLast(iSpin, ll, nn) * rad(:, nn, ll, iSpin)**2&
                  & + rad(:, nn, ll, iSpin) * (0.5_dp * radpp(:, nn, ll, iSpin)&
                  & + radp(:, nn, ll, iSpin) / abcissa - 0.5_dp * ll * (ll + 1) / abcissa**2&
                  & * rad(:, nn, ll, iSpin)))
            end do
          end do
        end do
      end do

      ! Divide by the density, but try to avoid singular terms
      do iSpin = 1, nSpin
        where (rho(:, iSpin) > 0.0_dp)
          avgPot(:, iSpin) = avgPot(:, iSpin) / rho(:, iSpin)
        end where
      end do

      ! Build effective orbital energies
      do iSpin = 1, nSpin
        do ll = 0, max_l
          nn = 0
          do iDummy1 = 1, num_alpha(ll)
            do iDummy2 = 1, poly_order(ll)
              nn = nn + 1
              ! uses -laplace / 2
              eps(iSpin, ll, nn)&
                  & = sum(weights * abcissa**2 * (rad(:, nn, ll, iSpin) * (-0.5_dp&
                  & * radpp(:, nn, ll, iSpin) - radp(:, nn, ll, iSpin) / abcissa&
                  & + 0.5_dp * ll * (ll + 1) / abcissa**2 * rad(:, nn, ll, iSpin)&
                  & + avgPot(:, iSpin) * rad(:, nn, ll, iSpin))))
            end do
          end do
        end do
      end do

      ! Shift current effective orbital energies so that HOAO agrees with GKS calculation
      do iSpin = 1, nSpin
        hoaoShift(iSpin) = eHoao(iSpin) - eps(iSpin, hoaoL(iSpin), hoaoN(iSpin))
        do ll = 0, max_l
          nn = 0
          do iDummy1 = 1, num_alpha(ll)
            do iDummy2 = 1, poly_order(ll)
              nn = nn + 1
              eps(iSpin, ll, nn) = eps(iSpin, ll, nn) + hoaoShift(iSpin)
            end do
          end do
        end do
      end do

      ! Probe convergence
      change_max = 0.0_dp
      do iSpin = 1, nSpin
        do ll = 0, max_l
          nn = 0
          do iDummy1 = 1, num_alpha(ll)
            do iDummy2 = 1, poly_order(ll)
              nn = nn + 1
              change_max = max(change_max, abs(eps(iSpin, ll, nn) - epsLast(iSpin, ll, nn)))
            end do
          end do
        end do
      end do

      tConverged = change_max <= scftol

      ! If self-consistency is reached, exit loop
      if (tConverged) exit lpSC

      do iSpin = 1, nSpin
        do ll = 0, max_l
          nn = 0
          do iDummy1 = 1, num_alpha(ll)
            do iDummy2 = 1, poly_order(ll)
              nn = nn + 1
              epsLast(iSpin, ll, nn) = eps(iSpin, ll, nn)
            end do
          end do
        end do
      end do

    end do lpSC

    ! Handle non-converged calculations
    if (.not. tConverged) then
      call error('Average potential NOT converged, maximal SC iterations exceeded.')
    end if

    ! Clean up and deal with numerical instability
    do iSpin = 1, nSpin
      where (abs(rho(:, iSpin)) < 1.0e-20_dp)
        avgPot(:, iSpin) = 0.0_dp
      end where
    end do

  end subroutine getAveragePotential


  !> Calculates the kinetic energy contribution based on the properties present in the average
  !! potential module. The result should always match the one obtained from the kinetic supervector.
  !! While this calculation is not strictly required, it is an important sanity-check.
  pure function getKineticEnergy(occ, abcissa, weights, max_l, num_alpha, poly_order, rad, radp,&
      & radpp) result(kinetic_energy)

    !> Occupation numbers
    real(dp), intent(in) :: occ(:, 0:, :)

    !> Numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> Numerical integration weights
    real(dp), intent(in) :: weights(:)

    !> Maximum angular momentum
    integer, intent(in) :: max_l

    !> Number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> Highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !! Wave function, 1st and 2nd derivative
    real(dp), intent(in) :: rad(:,:, 0:, :), radp(:,:, 0:, :), radpp(:,:, 0: ,:)

    !> Kinetic energy
    real(dp) :: kinetic_energy

    !! Spin index and number of spin channels
    integer :: iSpin, nSpin

    !! Principal quantum number
    integer :: nn

    !! Angular quantum number
    integer :: ll

    !! Dummy loop indices
    integer :: iDummy1, iDummy2

    nSpin = size(occ, dim=1)

    kinetic_energy = 0.0_dp

    do iSpin = 1, nSpin
      do ll = 0, max_l
        nn = 0
        do iDummy1 = 1, num_alpha(ll)
          do iDummy2 = 1, poly_order(ll)
            nn = nn + 1
            kinetic_energy = kinetic_energy&
                & + occ(iSpin, ll, nn) * sum(weights * (rad(:, nn, ll, iSpin) * (-0.5_dp&
                & * radpp(:, nn, ll, iSpin) - radp(:, nn, ll, iSpin) / abcissa&
                & + 0.5_dp * ll * (ll + 1) / abcissa**2 * rad(:, nn, ll, iSpin))) * abcissa**2)
          end do
        end do
      end do
    end do

  end function getKineticEnergy

end module average_potential
