!> Module for zero-order regular approximation for relativistic effects (ZORA) related routines.
module zora_routines

  use common_accuracy, only : dp
  use common_constants, only : cc
  use coulomb_potential, only : cou_pot
  use density, only : basis_times_basis, basis_1st_times_basis_1st_times_r2

  implicit none
  private

  public :: zora_t_correction, scaled_zora


contains

  !> Evaluates ZORA relativistic correction to kinetic energy matrix elements.
  !! mode=1: correction to kinetic energy matrix elements
  !! mode=2: additional terms for scaling matrix elements
  subroutine zora_t_correction(mode, tt, max_l, num_alpha, alpha, poly_order, num_mesh_points,&
      & weight, abcissa, vxc, nuc, pp, problemsize)

    !> determines correction mode, see above
    integer, intent(in) :: mode

    !> kinetic supervector
    real(dp), intent(out) :: tt(:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, start
    real(dp), allocatable :: kappa(:,:), kappa2(:,:), vtot(:,:)

    allocate(kappa(2, num_mesh_points))
    allocate(kappa2(2, num_mesh_points))
    allocate(vtot(2, num_mesh_points))

    tt(:,:,:,:) = 0.0_dp

    call potential_to_mesh(num_mesh_points, abcissa, vxc, nuc, pp, max_l, num_alpha, poly_order,&
        & alpha, problemsize, vtot)

    call kappa_to_mesh(num_mesh_points, vtot, kappa, kappa2)

    do ii = 0, max_l
      nn = 0
      do jj = 1, num_alpha(ii)
        do ll = 1, poly_order(ii)
          nn = nn + 1

          oo = nn - 1
          do kk = jj, num_alpha(ii)

            start = 1
            if (kk == jj) start = ll

            do mm = start, poly_order(ii)
              oo = oo + 1

              ! kinetic energy correction depends on spin via potential

              if (mode == 1) then

                tt(1, ii, nn, oo) = kinetic_part_1(num_mesh_points, weight, abcissa,&
                    & kappa(1, :), alpha(ii, jj), ll, alpha(ii, kk),&
                    & mm, ii) + kinetic_part_2(num_mesh_points, weight, abcissa,&
                    & kappa(1, :), alpha(ii, jj), ll, alpha(ii, kk),&
                    & mm, ii) * real(ii * (ii + 1), dp)

                tt(2, ii, nn, oo) = kinetic_part_1(num_mesh_points, weight, abcissa, kappa(2, :),&
                    & alpha(ii, jj), ll, alpha(ii, kk), mm, ii)&
                    & + kinetic_part_2(num_mesh_points, weight, abcissa, kappa(2, :),&
                    & alpha(ii, jj), ll, alpha(ii, kk), mm, ii) * real(ii * (ii + 1), dp)

              end if

              if (mode == 2) then

                ! calculate matrix elements needed for scaled ZORA
                ! prefactor 1/2 is included as the same subroutines as for t are used

                tt(1, ii, nn, oo) = kinetic_part_1(num_mesh_points, weight, abcissa, kappa2(1, :),&
                    & alpha(ii, jj), ll, alpha(ii, kk), mm, ii)&
                    & + kinetic_part_2(num_mesh_points, weight, abcissa, kappa2(1, :),&
                    & alpha(ii, jj), ll, alpha(ii, kk), mm, ii) * real(ii * (ii + 1), dp)

                tt(2, ii, nn, oo) = kinetic_part_1(num_mesh_points, weight, abcissa, kappa2(2, :),&
                    & alpha(ii, jj), ll, alpha(ii, kk), mm, ii)&
                    & + kinetic_part_2(num_mesh_points, weight, abcissa, kappa2(2, :),&
                    & alpha(ii, jj), ll, alpha(ii, kk), mm, ii) * real(ii * (ii + 1), dp)

              end if

              tt(1, ii, oo, nn) = tt(1, ii, nn, oo)
              tt(2, ii, oo, nn) = tt(2, ii, nn, oo)

            end do
          end do
        end do
      end do
    end do

  end subroutine zora_t_correction


  !>
  subroutine scaled_zora(eigval, max_l, num_alpha, alpha, poly_order, problemsize, num_mesh_points,&
      & weight, abcissa, vxc, nuc, pp, tt, cof, occ, eigval_scaled, zora_ekin)

    !> eigenvalues
    real(dp), intent(in) :: eigval(:,0:,:)

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

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

    !> kinetic supervector
    real(dp), intent(in) :: tt(0:,:,:)

    !> wavefunction coefficients
    real(dp), intent(in) :: cof(:,0:,:,:)

    !> occupation numbers
    real(dp), intent(in) :: occ(:,0:,:)

    !> zora scaled eigenvalues
    real(dp), intent(out) :: eigval_scaled(:,0:,:)

    !> zora kinetic energy contribution to total energy
    real(dp), intent(out) :: zora_ekin

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, ppp
    real(dp), allocatable :: zscale(:,:,:,:), zscale2(:,:,:,:)
    real(dp) :: dummy1, dummy2, tsol2, zora_ekin1, zora_ekin2

    allocate(zscale(2, 0:max_l, problemsize, problemsize))
    allocate(zscale2(2, 0:max_l, problemsize, problemsize))
    zscale(:,:,:,:) = 0.0_dp
    zscale2(:,:,:,:) = 0.0_dp
    eigval_scaled(:,:,:) = 0.0_dp
    zora_ekin = 0.0_dp
    zora_ekin1 = 0.0_dp
    zora_ekin2 = 0.0_dp
    tsol2 = 1.0_dp / cc**2

    call zora_t_correction(1, zscale, max_l, num_alpha, alpha, poly_order, num_mesh_points, weight,&
        & abcissa, vxc, nuc, pp, problemsize)
    call zora_t_correction(2, zscale2, max_l, num_alpha, alpha, poly_order, num_mesh_points,&
        & weight, abcissa, vxc, nuc, pp, problemsize)

    ! first get scaled eigenvalues

    ! sum over all angular momenta
    do ii = 0, max_l
      ! sum over all eigenvectors
      do jj = 1, num_alpha(ii) * poly_order(ii)
        oo = 0
        dummy1 = 0.0_dp
        dummy2 = 0.0_dp
        ! sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
        do kk = 1, num_alpha(ii)
          do ll = 1, poly_order(ii)
            oo = oo + 1
            ppp = 0
            ! other sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
            do mm = 1, num_alpha(ii)
              do nn = 1, poly_order(ii)
                ppp = ppp + 1
                ! occupation numbers do not enter here
                dummy1 = dummy1 + cof(1, ii, ppp, jj) * cof(1, ii, oo, jj) * tsol2&
                    & * (zscale(1, ii, oo, ppp) + 0.5_dp * (zscale2(1, ii, oo, ppp)&
                    & + tt(ii, oo, ppp)))
                dummy2 = dummy2 + cof(2, ii, ppp, jj) * cof(2, ii, oo, jj) * tsol2&
                    & * (zscale(2, ii, oo, ppp) + 0.5_dp * (zscale2(2, ii, oo, ppp)&
                    & + tt(ii, oo, ppp)))
              end do
            end do
          end do
        end do

        eigval_scaled(1, ii, jj) = eigval(1, ii, jj) / (1.0_dp + dummy1)
        eigval_scaled(2, ii, jj) = eigval(2, ii, jj) / (1.0_dp + dummy2)
      end do
    end do

    ! now, ZORA kinetic energy

    dummy1 = 0.0_dp
    dummy2 = 0.0_dp
    ! sum over all angular momenta
    do ii = 0, max_l
      ! sum over all eigenvectors
      do jj = 1, num_alpha(ii) * poly_order(ii)
        oo = 0
        ! sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
        do kk = 1, num_alpha(ii)
          do ll = 1, poly_order(ii)
            oo = oo + 1
            ppp = 0
            ! other sum over all basis functions in alpha and polynomial, i.e. prim. Slaters
            do mm = 1, num_alpha(ii)
              do nn = 1, poly_order(ii)
                ppp = ppp + 1
                ! dummy contains the non-relativistic kinetic energy operator applied to the
                ! relativistic ZORA wavefunction, debug only
                ! dummy1 = dummy1 + occ(1, ii, jj) * cof(1, ii, ppp, jj) * cof(1, ii, oo, jj)&
                !     & * tt(ii, oo, ppp)
                ! dummy2 = dummy2 + occ(2, ii, jj) * cof(2, ii, ppp, jj) * cof(2, ii, oo, jj)&
                !     & * tt(ii, oo, ppp)
                zora_ekin1 = zora_ekin1&
                  & + occ(1, ii, jj) * cof(1, ii, ppp, jj) * cof(1, ii, oo, jj)&
                  & * (tt(ii, oo, ppp) + zscale(1, ii, oo, ppp)&
                  & - eigval_scaled(1, ii, jj) * tsol2 * (0.5_dp * (&
                  & zscale2(1, ii, oo, ppp) + tt(ii, oo, ppp)) + zscale(1, ii, oo, ppp)))
                zora_ekin2 = zora_ekin2&
                  & + occ(2, ii, jj) * cof(2, ii, ppp, jj) * cof(2, ii, oo, jj)&
                  & * (tt(ii, oo, ppp) + zscale(2, ii, oo, ppp)&
                  & - eigval_scaled(2, ii, jj) * tsol2 * (0.5_dp * (&
                  & zscale2(2, ii, oo, ppp) + tt(ii, oo, ppp)) + zscale(2, ii, oo, ppp)))
              end do
            end do
          end do
        end do
      end do
    end do

    zora_ekin = zora_ekin1 + zora_ekin2

  end subroutine scaled_zora


  !> Calculates 0.5*\int_0^\inf r^2 kappa (d/dr R_A) (d/dr R_B) dr
  !! Pass either up or down total potential as kappa.
  pure function kinetic_part_1(num_mesh_points, weight, abcissa, kappa, alpha1, poly1, alpha2,&
      & poly2, ll)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weight factors of mesh
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas (radii), Becke mapping
    real(dp), intent(in) :: abcissa(:)

    !> either up or down total potential
    real(dp), intent(in) :: kappa(:)

    !> !> basis exponent of 1st basis derivative
    real(dp), intent(in) :: alpha1

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis derivative
    real(dp), intent(in) :: alpha2

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> result
    real(dp) :: kinetic_part_1

    !> auxiliary variable
    integer :: ii

    kinetic_part_1 = 0.0_dp

    do ii = 1, num_mesh_points

      kinetic_part_1 = kinetic_part_1 + weight(ii) * kappa(ii)&
          & * basis_1st_times_basis_1st_times_r2(alpha1, poly1, alpha2, poly2, ll, abcissa(ii))

    end do

    kinetic_part_1 = kinetic_part_1 * 0.5_dp

  end function kinetic_part_1


  !> Calculates \int_0^\inf R_B R_A kappa dr; multiply by l(l+1) in calling routine.
  !! Pass either up or down total potential as kappa.
  pure function kinetic_part_2(num_mesh_points, weight, abcissa, kappa, alpha1, poly1, alpha2,&
      & poly2, ll)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weight factors of mesh
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas (radii), Becke mapping
    real(dp), intent(in) :: abcissa(:)

    !> either up or down total potential
    real(dp), intent(in) :: kappa(:)

    !> basis exponent of 1st basis
    real(dp), intent(in) :: alpha1

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis
    real(dp), intent(in) :: alpha2

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> result
    real(dp) :: kinetic_part_2

    !> auxiliary variable
    integer :: ii

    kinetic_part_2 = 0.0_dp

    do ii = 1, num_mesh_points

      kinetic_part_2 = kinetic_part_2 + weight(ii) * kappa(ii)&
          & * basis_times_basis(alpha1, poly1, alpha2, poly2, ll, abcissa(ii))

    end do

    kinetic_part_2 = kinetic_part_2 * 0.5_dp

  end function kinetic_part_2


  !> kappa=V/(2*c^2-V), V total potential, c speed of light
  !! kappa2=kappa^2, i.e. square of kappa
  pure subroutine kappa_to_mesh(num_mesh_points, vtot, kappa, kappa2)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> total, spinpolarized potential on mesh
    real(dp), intent(in) :: vtot(:,:)

    !> kappa and kappa^2
    real(dp), intent(out) :: kappa(:,:), kappa2(:,:)

    !> auxiliary variables
    integer :: ii
    real(dp), parameter :: tsol2 = 2.0_dp * cc**2

    do ii = 1, num_mesh_points

      kappa(1, ii) = vtot(1, ii) / (tsol2 - vtot(1, ii))
      kappa(2, ii) = vtot(2, ii) / (tsol2 - vtot(2, ii))

      kappa2(1, ii) = kappa(1, ii)**2
      kappa2(2, ii) = kappa(2, ii)**2

    end do

  end subroutine kappa_to_mesh


  !> Calculates total, spinpolarized potential on mesh.
  subroutine potential_to_mesh(num_mesh_points, abcissa, vxc, nuc, pp, max_l, num_alpha,&
      & poly_order, alpha, problemsize, vtot)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas (radii), Becke mapping
    real(dp), intent(in) :: abcissa(:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

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

    !> total, spinpolarized potential on mesh
    real(dp), intent(out) :: vtot(:,:)

    !> total density matrix supervector (spins summed up)
    real(dp), allocatable :: ptot(:,:,:)

    !> coulomb potential on mesh
    real(dp), allocatable :: cpot(:)

    !> auxiliary variable
    integer :: ii

    allocate(cpot(num_mesh_points))

    cpot(:) = 0.0_dp
    vtot(:,:) = 0.0_dp

    ptot(:,:,:) = pp(1, :,:,:) + pp(2, :,:,:)

    call cou_pot(ptot, max_l, num_alpha, poly_order, alpha, problemsize, num_mesh_points, abcissa,&
        & cpot)

    do ii = 1, num_mesh_points

      vtot(1, ii) = - real(nuc, dp) / abcissa(ii) + cpot(ii) + vxc(ii, 1)
      vtot(2, ii) = - real(nuc, dp) / abcissa(ii) + cpot(ii) + vxc(ii, 2)

    end do

  end subroutine potential_to_mesh

end module zora_routines
