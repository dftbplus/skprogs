!> Module that provides common utilities used by the slateratom program.
module utilities

  use common_accuracy, only : dp

  implicit none
  private

  public :: check_convergence, check_electron_number
  public :: vector_length, fak, zeroOutCpotOfEmptyDensitySpinChannels
  public :: index_heap_sort


  !> Heap sort returning an index array.
  interface index_heap_sort
    module procedure index_heap_sort_real
  end interface


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
  pure subroutine check_convergence(pot_old, pot_new, max_l, problemsize, scftol, iScf, change_max,&
      & tConverged)

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

    if (iScf < 3) then
      tConverged = .false.
    end if

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

  end subroutine check_convergence


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


  !> Real case heap sort returning an index array.
  !> based on Numerical Recipes Software 1986-92
  subroutine index_heap_sort_real(indx, array, tolerance)

    !> Indexing array on return
    integer, intent(out) :: indx(:)

    !> Array of values to be sorted
    real(dp), intent(in) :: array(:)

    !> Tolerance for equality of two elements
    real(dp), intent(in), optional :: tolerance

    integer :: nn, ir, ij, il, ii, ik
    integer :: indxTmp
    real(dp) :: arrayTmp, tol

    if (present(tolerance)) then
      tol = tolerance
    else
      tol = epsilon(0.0_dp)
    end if

    do ii = 1, size(indx)
      indx(ii) = ii
    end do

    nn = size(array)

    if (nn <= 1) return
    il = nn / 2 + 1
    ir = nn
    ik = 1
    do while (ik == 1)
      if (il .gt. 1) then
        il = il - 1
        indxTmp = indx(il)
        arrayTmp = array(indxTmp)
      else
        indxTmp = indx(ir)
        arrayTmp = array(indxTmp)
        indx(ir) = indx(1)
        ir = ir - 1
        if (ir .lt. 1) then
          indx(1) = indxTmp
          return
        end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
        if (ij < ir) then
          if (array(indx(ij)) < array(indx(ij+1)) - tol) then
            ij = ij + 1
          end if
        end if
        if(arrayTmp < array(indx(ij)) - tol) then
          indx(ii) = indx(ij)
          ii = ij
          ij = 2 * ij
        else
          ij = ir + 1
        end if
      end do
      indx(ii) = indxTmp
    end do

  end subroutine index_heap_sort_real

end module utilities
