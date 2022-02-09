!> Module that serves output purposes during SCF, except in the ZORA case,
!! but even then the Coulomb matrix (J supermatrix) elements are calculated directly.
module coulomb_potential

  use common_accuracy, only : dp
  use integration, only : exp_int
  use core_overlap, only : v

  implicit none
  private

  public :: cou_pot


contains

  !> Calculates Coulomb potential on arbitraty set of points,
  !! by an analytical evaluation of the integrals indicated.
  !!                _                                         _
  !!               |                                           |
  !!               | 1  r      2            rmax               |
  !!  V(r)= 4*PI * | - int * r' * rho(r') + int  r' * rho (r') |
  !!               | r  0                    r                 |
  !!               |_                                         _|
  !!                          help1                help2
  subroutine cou_pot(ptot, max_l, num_alpha, poly_order, alpha, problemsize, num_mesh_points,&
      & abcissa, cpot)

    !> total density matrix supervector (spins summed up)
    real(dp), intent(in) :: ptot(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> coulomb potential on mesh
    real(dp), intent(out) :: cpot(:)

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, pp, nlp, nlq
    real(dp) :: alpha1
    real(dp), allocatable :: help1(:,:,:,:), help2(:,:,:,:)

    allocate(help1(num_mesh_points, 0:max_l, problemsize, problemsize))
    allocate(help2(num_mesh_points, 0:max_l, problemsize, problemsize))

    help1(:,:,:,:) = 0.0_dp
    help2(:,:,:,:) = 0.0_dp
    cpot(:) = 0.0_dp

    ! get integrals for pairs of basis functions
    do ii = 0, max_l
      ll = 0
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          ll = ll + 1
          oo = 0
          nlp = kk + ii
          do mm = 1, num_alpha(ii)

            ! exp_int has no notion of implicit "-" of alpha
            alpha1 = -(alpha(ii, jj) + alpha(ii, mm))

            do nn = 1, poly_order(ii)
              oo = oo + 1
              nlq = nn + ii

              ! integrals as indicated in comment, no normalization
              do pp = 1, num_mesh_points
                help1(pp, ii, ll, oo) = (exp_int(alpha1, nlp + nlq, abcissa(pp))&
                    & - exp_int(alpha1, nlp + nlq, 0.0_dp)) / abcissa(pp)
                help2(pp, ii, ll, oo) = - exp_int(alpha1, nlp + nlq - 1, abcissa(pp))
              end do

              ! add normalization of basis functions
              ! watch out for 2**(nlp+nlq+1) needed because variable integration ranges
              help1(:, ii, ll, oo) = help1(:, ii, ll, oo) * real(2**(nlp + nlq + 1), dp)&
                  & / sqrt(v(alpha(ii, jj), 2 * nlp) * v(alpha(ii, mm), 2 * nlq))
              help2(:, ii, ll, oo) = help2(:, ii, ll, oo) * real(2**(nlp + nlq + 1), dp)&
                  & / sqrt(v(alpha(ii, jj), 2 * nlp) * v(alpha(ii, mm), 2 * nlq))

            end do
          end do
        end do
      end do
    end do

    ! now actually get potential, multiply with density matrix
    do pp = 1, num_mesh_points
      do ii = 0, max_l
        ll = 0
        do jj = 1, num_alpha(ii)
          do kk = 1, poly_order(ii)
            ll = ll + 1
            oo = 0
            do mm = 1, num_alpha(ii)
              do nn = 1, poly_order(ii)
                oo = oo + 1
                cpot(pp) = cpot(pp) + ptot(ii, ll, oo)&
                    & * (help1(pp, ii, ll, oo) + help2(pp, ii, ll, oo))
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine cou_pot

end module coulomb_potential
