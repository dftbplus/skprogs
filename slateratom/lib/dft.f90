!> Module that provides various functionality related to density functional theory.
module dft

  use, intrinsic :: iso_c_binding, only : c_size_t
  use common_accuracy, only : dp
  use common_constants, only : pi, rec4pi
  use xcfunctionals, only : xcFunctional, getExcVxc_LDA_PW91,&
      & getExcVxc_GGA_PBE96, getExcVxc_GGA_BLYP, getExcVxc_LCY_PBE96, getExcVxc_LCY_BNL,&
      & getExcVxc_HYB_B3LYP, getExcVxc_HYB_PBE0, getExcVxc_CAMY_B3LYP, getExcVxc_CAMY_PBEh, &
      & getExcVxc_MGGA_TPSS, getExcVxc_MGGA_SCAN, getExcVxc_MGGA_r2SCAN, getExcVxc_MGGA_r4SCAN,&
      & getExcVxc_MGGA_TASK, getExcVxc_MGGA_TASK_CC
  use density, only : basis, basis_times_basis_times_r2, density_at_point, density_at_point_1st,&
      & density_at_point_2nd, tau_at_point, basis_1st_times_basis_1st_times_r2,&
      & basis_1st_times_basis_1st

  implicit none
  private

  public :: thomas_fermi_start_pot, density_grid, dft_exc_energy, dft_vxc_energy
  public :: dft_exc_matrixelement, xalpha
  public :: check_accuracy


contains

  !> Total potential to initialize a DFT calculation from Thomas-Fermi theory. This does not work as
  !! intended in the current code, since we do not have a numerical Coulomb-Potential.
  !!
  !! Generalized Thomas-Fermi atomic potential as published by R. Latter,
  !! Phys. Rev. 99, 510 (1955) eqn. 5/9 and implemented in Dirk Porezags scfatom.
  pure subroutine thomas_fermi_start_pot(abcissa, num_mesh_points, nuc, vxc)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> xc potential on grid for two spin channels vxc(:, 1) and vxc(:, 2)
    real(dp), intent(out) :: vxc(:,:)

    !> auxiliary variables
    integer :: ii
    real(dp) :: bb, tt, xx, rtx

    bb = (0.69395656_dp / real(nuc, dp))**(1.0_dp / 3.0_dp)

    do ii = 1, num_mesh_points

      xx = abcissa(ii) / bb
      rtx = sqrt(xx)

      tt = real(nuc, dp) / (1.0_dp + rtx * (0.02747_dp - xx * (0.1486_dp - 0.007298_dp * xx))&
          & + xx * (1.243_dp + xx * (0.2302_dp + 0.006944_dp * xx)))
      if (tt < 1.0_dp) tt = 1.0_dp

      vxc(ii, 1) = (tt / abcissa(ii)) / 2.0_dp
      vxc(ii, 2) = (tt / abcissa(ii)) / 2.0_dp

    end do

  end subroutine thomas_fermi_start_pot


  !> Calculate and store density and density derivatives on radial grid.
  !! Further calculates and stores exchange-correlation potential and energy density on grid.
  subroutine density_grid(pp, max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa, dzdr,&
      & dz, xcnr, omega, camAlpha, camBeta, rho, drho, ddrho, tau, vxc, vtau, exc, xalpha_const)

    !> density matrix supervector
    real(dp), intent(in) :: pp(:, 0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> dz/dr
    real(dp), intent(in) :: dzdr(:)

    !> step width in linear coordinates
    real(dp), intent(in) :: dz

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> range-separation parameter
    real(dp), intent(in) :: omega

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> density on grid
    real(dp), intent(out) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(out) :: drho(:,:)

    !> 2nd deriv. of density on grid
    real(dp), intent(out) :: ddrho(:,:)

    !> kinetic energy density on grid
    real(dp), intent(out) :: tau(:,:)

    !> orbital-dependent tau potential on grid
    real(dp), intent(out) :: vtau(:,:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> exchange parameter for X-Alpha exchange
    real(dp), intent(in) :: xalpha_const

    !! total density on grid (spins summed up)
    real(dp) :: rhotot

    !! difference between spin densities
    real(dp) :: rhodiff

    !! number of density grid points
    integer(c_size_t) :: nn

    !! density in libxc compatible format, i.e. rho/(4pi)
    real(dp), allocatable :: rhor(:,:)

    !! temporary storage
    real(dp), allocatable :: sigma(:,:)

    !! auxiliary variable
    integer :: ii

    nn = size(rho, dim=1)

    rho(:,:) = 0.0_dp
    drho(:,:) = 0.0_dp
    ddrho(:,:) = 0.0_dp
    exc(:) = 0.0_dp
    vxc(:,:) = 0.0_dp
    vtau(:,:) = 0.0_dp
    tau(:,:) = 0.0_dp

    ! get density on grid
    do ii = 1, num_mesh_points
      rho(ii, 1) = density_at_point(pp(1, :,:,:), max_l, num_alpha, poly_order, alpha, abcissa(ii))
      rho(ii, 2) = density_at_point(pp(2, :,:,:), max_l, num_alpha, poly_order, alpha, abcissa(ii))
    end do
    rho = max(rho, 0.0_dp)

    ! get density derivatives on grid
    if (xcFunctional%isGGA(xcnr) .or. xcFunctional%isMGGA(xcnr)&
        & .or. xcFunctional%isLongRangeCorrected(xcnr) .or. xcFunctional%isGlobalHybrid(xcnr)&
        & .or. xcFunctional%isCAMY(xcnr)) then
      do ii = 1, num_mesh_points

        drho(ii, 1) = density_at_point_1st(pp(1, :,:,:), max_l, num_alpha, poly_order, alpha,&
            & abcissa(ii))
        drho(ii, 2) = density_at_point_1st(pp(2, :,:,:), max_l, num_alpha, poly_order, alpha,&
            & abcissa(ii))

        ddrho(ii, 1) = density_at_point_2nd(pp(1, :,:,:), max_l, num_alpha, poly_order, alpha,&
            & abcissa(ii))
        ddrho(ii, 2) = density_at_point_2nd(pp(2, :,:,:), max_l, num_alpha, poly_order, alpha,&
            & abcissa(ii))
      end do
    end if

    ! divide by 4*pi to catch different normalization of spherical harmonics
    rhor = transpose(rho) * rec4pi

    ! get contracted gradients of the density
    if (xcFunctional%isGGA(xcnr) .or. xcFunctional%isMGGA(xcnr)&
        & .or. xcFunctional%isLongRangeCorrected(xcnr) .or. xcFunctional%isGlobalHybrid(xcnr)&
        & .or. xcFunctional%isCAMY(xcnr)) then
      allocate(sigma(3, nn))
      sigma(1, :) = drho(:, 1) * drho(:, 1) * rec4pi**2
      sigma(2, :) = drho(:, 1) * drho(:, 2) * rec4pi**2
      sigma(3, :) = drho(:, 2) * drho(:, 2) * rec4pi**2
    end if

    if (xcFunctional%isMGGA(xcnr)) then
      do ii = 1, num_mesh_points
        tau(ii, 1) = tau_at_point(pp(1, :,:,:), max_l, num_alpha, poly_order, alpha, abcissa(ii))
        tau(ii, 2) = tau_at_point(pp(2, :,:,:), max_l, num_alpha, poly_order, alpha, abcissa(ii))
      end do
    end if

    select case (xcnr)
    case(xcFunctional%HF_Exchange)
      return
    case(xcFunctional%X_Alpha)
      do ii = 1, num_mesh_points
        rhotot = (rho(ii, 1) + rho(ii, 2)) * rec4pi
        rhodiff = (rho(ii, 1) - rho(ii, 2)) * rec4pi
        ! Xalpha handled by our in-house routine
        call xalpha(rhotot, rhodiff, vxc(ii, :), exc(ii), xalpha_const)
      end do
      return
    case(xcFunctional%LDA_PW91)
      call getExcVxc_LDA_PW91(rho, exc, vxc)
    case(xcFunctional%GGA_PBE96)
      call getExcVxc_GGA_PBE96(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)
    case(xcFunctional%GGA_BLYP)
      call getExcVxc_GGA_BLYP(abcissa, dz, dzdr, rho, drho, sigma, exc, vxc)
    case(xcFunctional%LCY_PBE96)
      call getExcVxc_LCY_PBE96(abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)
    case(xcFunctional%LCY_BNL)
      call getExcVxc_LCY_BNL(abcissa, dz, dzdr, rho, drho, sigma, omega, exc, vxc)
    case(xcFunctional%HYB_PBE0)
      call getExcVxc_HYB_PBE0(abcissa, dz, dzdr, rho, drho, sigma, camAlpha, exc, vxc)
    case(xcFunctional%HYB_B3LYP)
      call getExcVxc_HYB_B3LYP(abcissa, dz, dzdr, rho, drho, sigma, camAlpha, exc, vxc)
    case(xcFunctional%CAMY_B3LYP)
      call getExcVxc_CAMY_B3LYP(abcissa, dz, dzdr, rho, drho, sigma, omega, camAlpha, camBeta, exc,&
          & vxc)
    case(xcFunctional%CAMY_PBEh)
      call getExcVxc_CAMY_PBEh(abcissa, dz, dzdr, rho, drho, sigma, omega, camAlpha, camBeta, exc,&
          & vxc)
    case(xcFunctional%MGGA_TPSS)
      call getExcVxc_MGGA_TPSS(abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)
    case(xcFunctional%MGGA_SCAN)
      call getExcVxc_MGGA_SCAN(abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)
    case(xcFunctional%MGGA_r2SCAN)
      call getExcVxc_MGGA_r2SCAN(abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)
    case(xcFunctional%MGGA_r4SCAN)
      call getExcVxc_MGGA_r4SCAN(abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)
    case(xcFunctional%MGGA_TASK)
      call getExcVxc_MGGA_TASK(abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)
    case(xcFunctional%MGGA_TASK_CC)
      call getExcVxc_MGGA_TASK_CC(abcissa, dz, dzdr, rho, drho, sigma, tau, exc, vxc, vtau)
    case default
      write(*, '(A,I2,A)') 'XCNR=', xcnr, ' not implemented!'
      stop
    end select

 end subroutine density_grid


  !> Calculates DFT xc-energy from energy density and electron density on grid.
  pure subroutine dft_exc_energy(num_mesh_points, rho, exc, weight, abcissa, exc_energy)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> exc energy density on grid
    real(dp), intent(in) :: exc(:)

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> xc-energy obtained from the energy density and electron density on the grid
    real(dp), intent(out) :: exc_energy

    !> auxiliary variable
    integer :: ii

    exc_energy = 0.0_dp

    do ii = 1, num_mesh_points
      exc_energy = exc_energy + weight(ii) * exc(ii) * (rho(ii, 1) + rho(ii, 2)) * abcissa(ii)**2
    end do

    !
    ! For usual DFT functionals E_xc=\int \rho \eps(\rho,\zeta) d^3r
    ! so there is only one exchange-correlation energy density \eps(\rho,\zeta) and
    ! exc_energy could be a scalar without problems.
    !

  end subroutine dft_exc_energy


  !> Calculates vxc contribution for double counting correction.
  pure subroutine dft_vxc_energy(num_mesh_points, rho, vxc, weight, abcissa, vxc_energy)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> vxc contribution for double counting correction
    real(dp), intent(out) :: vxc_energy(2)

    !> auxiliary variable
    integer :: ii

    vxc_energy(:) = 0.0_dp

    do ii = 1, num_mesh_points
      vxc_energy(1) = vxc_energy(1) + weight(ii) * vxc(ii, 1) * (rho(ii, 1)) * abcissa(ii)**2
      vxc_energy(2) = vxc_energy(2) + weight(ii) * vxc(ii, 2) * (rho(ii, 2)) * abcissa(ii)**2
    end do

  end subroutine dft_vxc_energy


  !> Calculates a single matrix element of the exchange correlation potential.
  pure subroutine dft_exc_matrixelement(xcnr, num_mesh_points, weight, abcissa, vxc, vtau, alpha1,&
      & poly1, alpha2, poly2, ll, exc_matrixelement)

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> orbital-dependent tau potential on grid
    real(dp), intent(in) :: vtau(:,:)

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

    !> single matrix element of the exchange correlation potential
    real(dp), intent(out) :: exc_matrixelement(2)

    !! stores product of two basis functions and r^2
    real(dp) :: basis

    !! auxiliary variable
    integer :: ii

    exc_matrixelement(:) = 0.0_dp

    do ii = 1, num_mesh_points
      basis = basis_times_basis_times_r2(alpha1, poly1, alpha2, poly2, ll, abcissa(ii))

      exc_matrixelement(1) = exc_matrixelement(1) - weight(ii) * vxc(ii, 1) * basis
      exc_matrixelement(2) = exc_matrixelement(2) - weight(ii) * vxc(ii, 2) * basis

      if (xcFunctional%isMGGA(xcnr)) then
        ! «basis» for MGGA is a product of gradients
        basis = basis_1st_times_basis_1st_times_r2(alpha1, poly1, alpha2, poly2, ll, abcissa(ii))

        exc_matrixelement(1) = exc_matrixelement(1) - weight(ii) * vtau(ii, 1) * basis * 0.5_dp
        exc_matrixelement(2) = exc_matrixelement(2) - weight(ii) * vtau(ii, 2) * basis * 0.5_dp
      end if
    end do

  end subroutine dft_exc_matrixelement


  !> Calculates Xalpha potential and energy density.
  !!
  !! alpha=2/3 recovers the Gaspar/Kohn/Sham functional commonly used as exchange part in most
  !! current LSDA and GGA functionals. The original Slater exchange is recoverd with alpha=1.
  subroutine xalpha(rhotot, rhodiff, vxc, exc, xalpha_const)

    !> total density on grid (spins summed up)
    real(dp), intent(in) :: rhotot

    !> difference between spin densities
    real(dp), intent(in) :: rhodiff

    !> Xalpha energy density and potential
    real(dp), intent(out) :: exc, vxc(2)

    !> exchange parameter for X-Alpha exchange
    real(dp), intent(in) :: xalpha_const

    !> auxiliary variables
    real(dp) :: third, fourthird, vfac, cx, fzeta, dfzeta, eps0, eps1, spinpart, zeta

    third = 1.0_dp / 3.0_dp
    fourthird = 4.0_dp / 3.0_dp
    vfac = 2.0_dp**third
    cx = 0.75_dp * (3.0_dp / pi)**third

    if (abs(rhotot) < 1.0d-12) then
      exc = 0.0_dp
      vxc(:) = 0.0_dp
      return
    end if

    zeta = rhodiff / rhotot

    if (abs(zeta) > 1.0d12) write(*,*) 'ZETA LARGE IN X-ALPHA'

    fzeta = ((1.0_dp + zeta)**fourthird + (1.0_dp - zeta)**fourthird - 2.0_dp)&
        & / (2.0_dp * (vfac - 1.0_dp))
    dfzeta = fourthird * ((1.0_dp + zeta)**third - (1.0_dp - zeta)**third)&
        & / (2.0_dp * (vfac - 1.0_dp))

    eps0 = -3.0_dp / 2.0_dp * xalpha_const * cx * rhotot**third
    eps1 = vfac * eps0

    exc = eps0 + (eps1 - eps0) * fzeta

    spinpart = (eps1 - eps0) * dfzeta * (1.0_dp - zeta)

    vxc(1) = fourthird * exc + spinpart
    vxc(2) = fourthird * exc - spinpart

  end subroutine xalpha


  !> Tests integration to check the accuracy of the radial mesh by integrating the square of a
  !! primitive Slater basis function which are analytically normalized to 1.0_dp.
  subroutine check_accuracy(weight, abcissa, num_mesh_points, max_l, num_alpha, alpha, poly_order)

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> integrated primitive Slater basis function
    real(dp) :: value

    !> auxiliary variables
    integer :: ii, jj, kk, ll

    do ii = 0, max_l
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          value = 0.0_dp
          do ll = 1, num_mesh_points

            value = value + weight(ll) * abcissa(ll)**2&
                & * basis(alpha(ii, jj), kk, ii, abcissa(ll))**2

          end do
          if (abs(1.0_dp - value) > 1.0d-12) then
            write(*, '(A,F12.6,I3,E12.3)') 'WARNING: Integration bad for basis function ',&
                & alpha(ii, jj), kk + ii - 1, abs(1.0_dp - value)
            write(*, '(A)') 'Accuracy is not better than 1.0d-12'
          end if
        end do
      end do
    end do

  end subroutine check_accuracy

end module dft
