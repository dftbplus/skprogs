!> Module that provides various functionality related to density functional theory.
module dft

  use, intrinsic :: iso_c_binding, only : c_size_t
  use common_accuracy, only : dp
  use common_constants, only : pi
  use density, only : basis, basis_times_basis_times_r2, density_at_point, density_at_point_1st,&
      & density_at_point_2nd
  use xc_f90_lib_m, only : xc_f90_func_t, xc_f90_func_info_t, xc_f90_func_init,&
      & xc_f90_func_get_info, xc_f90_lda_exc_vxc, xc_f90_gga_exc_vxc, xc_f90_func_end, XC_LDA_X,&
      & XC_LDA_C_PW, XC_GGA_X_PBE, XC_GGA_C_PBE, XC_POLARIZED

  implicit none
  private

  public :: dft_start_pot, density_grid, dft_exc_energy, dft_vxc_energy
  public :: dft_exc_matrixelement, xalpha
  public :: check_accuracy
  public :: derive, radial_divergence, derive1_5, derive2_5


contains

  !> Total potential to initialize a DFT calculation from Thomas-Fermi theory. This does not work as
  !! intended in the current code, since we do not have a numerical Coulomb-Potential.
  !!
  !! Generalized Thomas-Fermi atomic potential as published by R. Latter,
  !! Phys. Rev. 99, 510 (1955) eqn. 5/9 and implemented in Dirk Porezags scfatom.
  pure subroutine dft_start_pot(abcissa, num_mesh_points, nuc, vxc)

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

  end subroutine dft_start_pot


  !> Calculate and store density and density derivatives on radial grid.
  !! Further calculates and stores exchange-correlation potential and energy density on grid.
  subroutine density_grid(pp, max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa, dzdr,&
      & dz, xcnr, rho, drho, ddrho, vxc, exc, xalpha_const)

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

    !> density on grid
    real(dp), intent(out) :: rho(:,:)

    !> 1st deriv. of density on grid
    real(dp), intent(out) :: drho(:,:)

    !> 2nd deriv. of density on grid
    real(dp), intent(out) :: ddrho(:,:)

    !> xc potential on grid
    real(dp), intent(out) :: vxc(:,:)

    !> exc energy density on grid
    real(dp), intent(out) :: exc(:)

    !> exchange parameter for X-Alpha exchange
    real(dp), intent(in) :: xalpha_const

    !> total density on grid (spins summed up)
    real(dp) :: rhotot

    !> difference between spin densities
    real(dp) :: rhodiff

    !> libxc related objects
    type(xc_f90_func_t) :: xcfunc_x, xcfunc_c
    type(xc_f90_func_info_t) :: xcinfo

    !> number of density grid points
    integer(c_size_t) :: nn

    !> exchange and correlation energy on grid
    real(dp), allocatable :: ex(:), ec(:)

    !> exchange and correlation potential on grid
    real(dp), allocatable :: vx(:,:), vc(:,:)

    !> density in libxc compatible shape
    real(dp), allocatable :: tmprho(:,:)

    !> temporary storage
    real(dp), allocatable :: tmpv(:), tmpv2(:)
    real(dp), allocatable :: tmpsigma(:,:), vxsigma(:,:), vcsigma(:,:)

    !> auxiliary variables
    integer :: ii, ispin, ispin2, isigma

    !> pre-factor to catch different normalization of spherical harmonics
    real(dp), parameter :: rec4pi = 1.0_dp / (4.0_dp * pi)

    if (xcnr == 0) return
    if (xcnr == 2) then
      call xc_f90_func_init(xcfunc_x, XC_LDA_X, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
      call xc_f90_func_init(xcfunc_c, XC_LDA_C_PW, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
    elseif (xcnr == 3) then
      call xc_f90_func_init(xcfunc_x, XC_GGA_X_PBE, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
      call xc_f90_func_init(xcfunc_c, XC_GGA_C_PBE, XC_POLARIZED)
      xcinfo = xc_f90_func_get_info(xcfunc_x)
    end if

    do ii = 1, num_mesh_points
      rho(ii, 1) = density_at_point(pp(1, :,:,:), max_l, num_alpha, poly_order, alpha, abcissa(ii))
      rho(ii, 2) = density_at_point(pp(2, :,:,:), max_l, num_alpha, poly_order, alpha, abcissa(ii))
    end do

    rho = max(rho, 0.0_dp)

    if (xcnr > 2) then

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
    if (xcnr == 1) then
      do ii = 1, num_mesh_points
        rhotot = (rho(ii, 1) + rho(ii, 2)) * rec4pi
        rhodiff = (rho(ii, 1) - rho(ii, 2)) * rec4pi
        call xalpha(rhotot, rhodiff, vxc(ii, :), exc(ii), xalpha_const)
      end do
    else if (xcnr == 2) then
      nn = size(rho, dim=1)
      allocate(tmprho(2, nn))
      allocate(ex(nn))
      allocate(ec(nn))
      allocate(vx(2, nn))
      allocate(vc(2, nn))
      tmprho(:,:) = transpose(rho) * rec4pi
      call xc_f90_lda_exc_vxc(xcfunc_x, nn, tmprho(1, 1), ex(1), vx(1, 1))
      call xc_f90_lda_exc_vxc(xcfunc_c, nn, tmprho(1, 1), ec(1), vc(1, 1))
      vxc(:,:) = transpose(vx + vc)
      exc = ec + ex
    else if (xcnr == 3) then
      nn = size(rho, dim=1)
      allocate(tmprho(2, nn))
      allocate(ex(nn))
      allocate(ec(nn))
      allocate(vx(2, nn))
      allocate(vc(2, nn))
      allocate(tmpsigma(3, nn))
      allocate(vxsigma(3, nn))
      allocate(vcsigma(3, nn))
      allocate(tmpv(nn))
      allocate(tmpv2(nn))
      tmprho(:,:) = transpose(rho) * rec4pi
      tmpsigma(1, :) = drho(:, 1) * drho(:, 1) * rec4pi * rec4pi
      tmpsigma(2, :) = drho(:, 1) * drho(:, 2) * rec4pi * rec4pi
      tmpsigma(3, :) = drho(:, 2) * drho(:, 2) * rec4pi * rec4pi
      call xc_f90_gga_exc_vxc(xcfunc_x, nn, tmprho(1, 1), tmpsigma(1, 1), ex(1), vx(1, 1),&
          & vxsigma(1, 1))
      call xc_f90_gga_exc_vxc(xcfunc_c, nn, tmprho(1, 1), tmpsigma(1, 1), ec(1), vc(1, 1),&
          & vcsigma(1, 1))
      vxc = transpose(vx + vc)
      do ispin = 1, 2
        ! the other spin
        ispin2 = 3 - ispin
        ! 1 for spin up, 3 for spin down
        isigma = 2 * ispin - 1
        tmpv(:) = (vxsigma(isigma, :) + vcsigma(isigma, :)) * drho(:, ispin) * rec4pi
        call radial_divergence(tmpv, abcissa, dz, tmpv2, dzdr)
        vxc(:, ispin) = vxc(:, ispin) - 2.0_dp * tmpv2
        tmpv(:) = (vxsigma(2, :) + vcsigma(2, :)) * drho(:, ispin2) * rec4pi
        call radial_divergence(tmpv, abcissa, dz, tmpv2, dzdr)
        vxc(:, ispin) = vxc(:, ispin) - tmpv2
      end do
      exc = ex + ec
    else
      write(*, '(A,I2,A)') 'XCNR=', xcnr, ' not implemented'
      stop
    end if

    call xc_f90_func_end(xcfunc_x)
    call xc_f90_func_end(xcfunc_c)

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
  pure subroutine dft_exc_matrixelement(num_mesh_points, weight, abcissa, vxc, alpha1, poly1,&
      & alpha2, poly2, ll, exc_matrixelement)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

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

    !> stores product of two basis functions and r^2
    real(dp) :: basis

    !> auxiliary variable
    integer :: ii

    exc_matrixelement(:) = 0.0_dp

    do ii = 1, num_mesh_points
      basis = basis_times_basis_times_r2(alpha1, poly1, alpha2, poly2, ll, abcissa(ii))

      exc_matrixelement(1) = exc_matrixelement(1) - weight(ii) * vxc(ii, 1) * basis
      exc_matrixelement(2) = exc_matrixelement(2) - weight(ii) * vxc(ii, 2) * basis
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


  !>
  pure subroutine radial_divergence(ff, rr, dr, rdiv, jacobi)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: rr(:)
    real(dp), intent(in) :: dr
    real(dp), intent(out) :: rdiv(:)
    real(dp), intent(in), optional :: jacobi(:)

    call derive1_5(ff, dr, rdiv, jacobi)
    rdiv = rdiv + 2.0_dp / rr * ff

  end subroutine radial_divergence


  !>
  pure subroutine derive(ff, dx, jacobi)
    real(dp), intent(inout) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(in), optional :: jacobi(:)

    real(dp), allocatable :: tmp1(:)
    integer :: nn

    nn = size(ff)
    allocate(tmp1(nn))
    tmp1(:) = ff
    ff(2:nn - 1) = (ff(3:nn) - ff(1:nn - 2)) / (2.0 * dx)
    ff(1) = (tmp1(2) - tmp1(1)) / dx
    ff(nn) = (tmp1(nn) - tmp1(nn - 1)) / dx
    if (present(jacobi)) then
      ff = ff * jacobi
    end if

  end subroutine derive


  !>
  pure subroutine derive1_5(ff, dx, dfdx, dudx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: dfdx(:)
    real(dp), intent(in), optional :: dudx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) =&
        reshape([&
        & -25.0_dp,  48.0_dp, -36.0_dp,  16.0_dp, -3.0_dp,&
        &  -3.0_dp, -10.0_dp,  18.0_dp,  -6.0_dp,  1.0_dp,&
        &   1.0_dp,  -8.0_dp,   0.0_dp,   8.0_dp, -1.0_dp,&
        &  -1.0_dp,   6.0_dp, -18.0_dp,  10.0_dp,  3.0_dp,&
        &   3.0_dp,  -16.0_dp, 36.0_dp, -48.0_dp, 25.0_dp], [np, np])

    integer :: ngrid
    integer :: ii

    ngrid = size(ff)
    do ii = 1, nleft
      dfdx(ii) = dot_product(coeffs(:, ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      dfdx(ii) = dot_product(coeffs(:, imiddle), ff(ii - nleft:ii + nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      dfdx(ii) = dot_product(coeffs(:, np - (ngrid - ii)), ff(ngrid - np + 1:ngrid))
    end do

    if (present(dudx)) then
      dfdx = dfdx * (dudx / (dxprefac * dx))
    else
      dfdx = dfdx / (dxprefac * dx)
    end if

  end subroutine derive1_5


  !>
  pure subroutine derive2_5(ff, dx, d2fdx2, dudx, d2udx2, dfdx)
    real(dp), intent(in) :: ff(:)
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: d2fdx2(:)
    real(dp), intent(in), optional :: dudx(:), d2udx2(:)
    real(dp), intent(out), target, optional :: dfdx(:)

    integer, parameter :: np = 5
    integer, parameter :: nleft = np / 2
    integer, parameter :: nright = nleft
    integer, parameter :: imiddle = nleft + 1
    real(dp), parameter :: dxprefac = 12.0_dp
    real(dp), parameter :: coeffs(np, np) = &
        reshape([ &
        &  35.0_dp, -104.0_dp, 114.0_dp, -56.0_dp, 11.0_dp, &
        &  11.0_dp, -20.0_dp, 6.0_dp, 4.0_dp, -1.0_dp, &
        &  -1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp, &
        &  -1.0_dp, 4.0_dp, 6.0_dp, -20.0_dp, 11.0_dp, &
        &  11.0_dp, -56.0_dp, 114.0_dp, -104.0_dp, 35.0_dp], [np, np])

    integer :: ngrid
    integer :: ii
    real(dp), allocatable, target :: dfdxlocal(:)
    real(dp), pointer :: pdfdx(:)

    ngrid = size(ff)
    if (present(dfdx)) then
      pdfdx => dfdx
    elseif (present(d2udx2)) then
      allocate(dfdxlocal(ngrid))
      pdfdx => dfdxlocal
    end if

    do ii = 1, nleft
      d2fdx2(ii) = dot_product(coeffs(:, ii), ff(1:np))
    end do
    do ii = nleft + 1, ngrid - nright
      d2fdx2(ii) = dot_product(coeffs(:, imiddle), ff(ii - nleft:ii + nright))
    end do
    do ii = ngrid - nright + 1, ngrid
      d2fdx2(ii) = dot_product(coeffs(:, np - (ngrid - ii)), ff(ngrid - np + 1:ngrid))
    end do

    if (present(dudx)) then
      d2fdx2 = d2fdx2 * (dudx * dudx / (dxprefac * dx * dx))
    else
      d2fdx2 = d2fdx2 / (dxprefac * dx * dx)
    end if

    if (present(d2udx2) .or. present(dfdx)) then
      call derive1_5(ff, dx, pdfdx)
      if (present(d2udx2)) then
        d2fdx2 = d2fdx2 + pdfdx * d2udx2
      end if
      if (present(dfdx) .and. present(dudx)) then
        dfdx = dfdx * dudx
      end if
    end if

  end subroutine derive2_5

end module dft
