!> Module that contains the two-center integrator routines for tabulating Hamiltonian and overlap.
module twocnt

  use common_accuracy, only : dp
  use common_constants, only : pi, rec4pi
  use common_anglib, only : initGaunt, freeGaunt, realGaunt
  use common_coordtrans, only : coordtrans_becke_12, coordtrans_radial_becke2
  use common_sphericalharmonics, only : TRealTessY, TRealTessY_init
  use common_quadratures, only : TQuadrature, gauss_legendre_quadrature, gauss_chebyshev_quadrature
  use common_gridgenerator, only : gengrid1_1, gengrid2_2
  use common_partition, only : partition_becke_homo
  use common_splines, only : spline3ders
  use common_fifo, only : TFiFoReal2
  use common_interpolation, only : get2ndNaturalSplineDerivs, get_cubic_spline
  use common_poisson, only : solvePoisson, TBeckeGridParams, TBeckeIntegrator,&
      & TBeckeIntegrator_init, TBeckeIntegrator_setKernelParam, TBeckeIntegrator_precompFdMatrix,&
      & TBeckeIntegrator_buildLU, TBeckeIntegrator_getCoords, TBeckeIntegrator_solveHelmholz

  use gridorbital, only : TGridorb2

  use, intrinsic :: iso_c_binding, only : c_size_t

  use xc_f03_lib_m, only : xc_f03_func_t, xc_f03_func_info_t, xc_f03_func_init, xc_f03_func_end,&
      & xc_f03_func_get_info, xc_f03_lda_vxc, xc_f03_gga_vxc, XC_LDA_X, XC_LDA_X_YUKAWA,&
      & XC_LDA_C_PW, XC_GGA_X_PBE, XC_GGA_C_PBE, XC_GGA_X_B88, XC_GGA_C_LYP, XC_GGA_X_SFAT_PBE,&
      & XC_HYB_GGA_XC_PBEH, XC_HYB_GGA_XC_B3LYP, XC_HYB_GGA_XC_CAMY_B3LYP,&
      & XC_HYB_GGA_XC_CAMY_PBEH, XC_UNPOLARIZED, xc_f03_func_set_ext_params

  implicit none
  private

  public :: TTwocntInp, TAtomdata, TIntegMap
  public :: get_twocenter_integrals


  ! Holds properties associated with a single atom.
  type TAtomdata

    !> number of basis functions
    integer :: nBasis

    !> angular momenta
    integer, allocatable :: angmoms(:)

    !> number of core orbitals
    integer :: nCore

    !> angular momenta of core orbitals
    integer, allocatable :: coreAngmoms(:)

    !> occupation of core orbitals
    real(dp), allocatable :: coreOcc(:)

    !> radial grid-orbital portion of core orbitals
    type(TGridorb2), allocatable  :: coreRad(:)

    !> radial grid-orbital portion and 1st/2nd derivative
    type(TGridorb2), allocatable :: rad(:), drad(:), ddrad(:)

    !> atomic potential on grid
    type(TGridorb2) :: pot

    !> atomic density and 1st/2nd derivative on grid
    type(TGridorb2) :: rho, drho, ddrho

  end type TAtomdata


  !> Holds parsed input for twocnt.
  type TTwocntInp

    !> true, if heteronuclear dimer is present
    logical :: tHetero

    !> true, if density superposition is requested, otherwise potential superposition is applied
    logical :: tDensitySuperpos

    !> start grid distance
    real(dp) :: r0

    !> grid separation, i.e. stepwidth
    real(dp) :: dr

    !> convergence criteria for Hamiltonian and overlap matrix elements
    real(dp) :: epsilon

    !> maximum grid distance
    real(dp) :: maxdist

    !> number of integration points
    integer :: ninteg1, ninteg2

    !> range-separation parameter
    real(dp) :: omega

    !> CAM alpha parameter
    real(dp) :: camAlpha

    !> CAM beta parameter
    real(dp) :: camBeta

    !> number of radial and angular Becke integration points
    integer :: nRadial, nAngular

    !> maximum angular momentum for Becke integration
    integer :: ll_max

    !> scaling factor of Becke transformation
    real(dp) :: rm

    !> xc-functional type
    !! (1: LDA-PW91, 2: GGA-PBE, 3: GGA-BLYP, 4: LCY-PBE, 5: LCY-BNL, 6: PBE0, 7: B3LYP,
    !! 8: CAMY-B3LYP, 9: CAMY-PBEh)
    integer :: iXC

    !! true, if a global hybrid functional is requested
    logical :: tGlobalHybrid

    !! true, if a (long-range corrected) range-separated hybrid functional is requested
    logical :: tLC

    !! true, if a CAM functional is requested
    logical :: tCam

    !> atomic properties of slateratom code, in the homonuclear case only atom1 is read
    type(TAtomdata) :: atom1, atom2

  end type TTwocntInp


  !> Type for mapping integrals.
  type TIntegMap

    !> number of all nonzero two-center integrals between orbitals of two atoms
    integer :: ninteg

    !> Indicates for every integral the integrands:
    !!
    !! o type(1, ii): index of orbital on first atom for integral ii
    !! o type(2, ii): index of orbital on second atom for integral ii
    !! o type(3, ii): interaction type for integral ii: (0 - sigma, 1 - pi, ...)
    integer, allocatable :: type(:,:)

    !> Indicates which integral corresponds to a given (i1, i2, mm) combination,
    !! where i1 and i2 are the orbital indices on the two atoms and mm the
    !! interaction type. If the integral vanishes, the corresponding element is 0.
    integer, allocatable :: index(:,:,:)

  end type TIntegMap


contains

  !> Calculates Hamiltonian and overlap matrix elements for different dimer distances.
  subroutine get_twocenter_integrals(inp, imap, skham, skover)

    !> parsed twocnt input instance
    type(TTwocntInp), intent(in), target :: inp

    !> integral mapping instance
    type(TIntegMap), intent(out) :: imap

    !> resulting Hamiltonian and overlap matrices
    real(dp), intent(out), allocatable :: skham(:,:), skover(:,:)

    !! abscissas and weight instances for numerical quadrature
    type(TQuadrature) :: quads(2)

    !! pointer to atomic properties of dimer atoms
    type(TAtomdata), pointer :: atom1, atom2

    !! database that holds Hamiltonian and overlap matrices
    type(TFiFoReal2) :: hamfifo, overfifo

    !! radial full-range Hartree-Fock quadrature information
    type(TQuadrature) :: radialHFQuadrature

    !! integration grids of dimer atoms, holding spherical coordinates (r, theta)
    real(dp), allocatable, target :: grid1(:,:), grid2(:,:), rr3(:,:)

    !! dot product of unit distance vectors and integration weights
    real(dp), allocatable :: dots(:), weights(:), dummyWeights(:)

    !! relative density integration error for all dimer distances of a batch
    real(dp), allocatable :: denserr(:)

    !! buffer holding Hamiltonian and overlap of current distance batch
    real(dp), allocatable :: skhambuffer(:,:), skoverbuffer(:,:)

    !! arbitrary dummy real array, unused for homonuclear Becke partitioning
    real(dp) :: beckepars(1)

    !! maximal density integration error
    real(dp) :: denserrmax

    !! current dimer distance
    real(dp) :: dist

    !! maximum absolute Hamiltonian or overlap matrix element
    real(dp) :: maxabs

    !! maximum dimer distance
    real(dp) :: maxdist

    !! iterates through a batch of dimer distances
    integer :: ir

    !! number of batches for which SK-integrals got calculated
    integer :: nBatch

    !! number of dimer distances in a single batch
    integer :: nBatchline

    !! true, if dimer distances are shall dynamically be extended if convergency isn't reached
    logical :: tDynlen

    !! true, if maximum absolute Hamiltonian or overlap matrix element is below given tolerance
    logical :: tConverged

    !! libxc related objects
    type(xc_f03_func_t) :: xcfunc_xc, xcfunc_x, xcfunc_c
    type(xc_f03_func_info_t) :: xcinfo

    !! Becke integrator instances
    type(TBeckeIntegrator) :: beckeInt

    !! grid characteristics
    type(TBeckeGridParams) :: beckeGridParams

    !! number of radial and angular integration abscissas
    integer :: nRad, nAng

    if (inp%iXC == 1) then
      call xc_f03_func_init(xcfunc_x, XC_LDA_X, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_x)
      call xc_f03_func_init(xcfunc_c, XC_LDA_C_PW, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_c)
    elseif (inp%iXC == 2) then
      call xc_f03_func_init(xcfunc_x, XC_GGA_X_PBE, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_x)
      call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_c)
    elseif (inp%iXC == 3) then
      call xc_f03_func_init(xcfunc_x, XC_GGA_X_B88, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_x)
      call xc_f03_func_init(xcfunc_c, XC_GGA_C_LYP, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_c)
    elseif (inp%iXC == 4) then
      call xc_f03_func_init(xcfunc_x, XC_GGA_X_SFAT_PBE, XC_UNPOLARIZED)
      call xc_f03_func_set_ext_params(xcfunc_x, [inp%omega])
      xcinfo = xc_f03_func_get_info(xcfunc_x)
      call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_c)
    elseif (inp%iXC == 5) then
      call xc_f03_func_init(xcfunc_x, XC_LDA_X_YUKAWA, XC_UNPOLARIZED)
      call xc_f03_func_set_ext_params(xcfunc_x, [inp%omega])
      xcinfo = xc_f03_func_get_info(xcfunc_x)
      call xc_f03_func_init(xcfunc_c, XC_GGA_C_PBE, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_c)
    elseif (inp%iXC == 6) then
      call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_PBEH, XC_UNPOLARIZED)
      xcinfo = xc_f03_func_get_info(xcfunc_xc)
    elseif (inp%iXC == 7) then
      call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_B3LYP, XC_UNPOLARIZED)
      call xc_f03_func_set_ext_params(xcfunc_xc, [0.20_dp, 0.72_dp, 0.81_dp])
      xcinfo = xc_f03_func_get_info(xcfunc_xc)
    elseif (inp%iXC == 8) then
      call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_CAMY_B3LYP, XC_UNPOLARIZED)
      call xc_f03_func_set_ext_params(xcfunc_xc, [inp%camAlpha + inp%camBeta, -inp%camBeta,&
          & inp%omega, 0.81_dp])
      xcinfo = xc_f03_func_get_info(xcfunc_xc)
    elseif (inp%iXC == 9) then
      call xc_f03_func_init(xcfunc_xc, XC_HYB_GGA_XC_CAMY_PBEH, XC_UNPOLARIZED)
      call xc_f03_func_set_ext_params(xcfunc_xc, [inp%camAlpha + inp%camBeta, -inp%camBeta,&
          & inp%omega])
      xcinfo = xc_f03_func_get_info(xcfunc_xc)
    end if

    if (inp%tLC .or. inp%tCam) then
      beckeGridParams%nRadial = inp%nRadial
      beckeGridParams%nAngular = inp%nAngular
      beckeGridParams%ll_max = inp%ll_max
      beckeGridParams%rm = inp%rm

      ! inititalize the Becke integrator
      call TBeckeIntegrator_init(beckeInt, beckeGridParams)

      call TBeckeIntegrator_setKernelParam(beckeInt, 0.1e-16_dp)
      call TBeckeIntegrator_precompFdMatrix(beckeInt)
      call TBeckeIntegrator_buildLU(beckeInt)

      !== Note: this is a workaround! ==
      ! we generate the LU decomposition for \omega \approx 0 and copy the LU decomposition into H4
      ! and permutation matrix into ipiv4
      beckeInt%fdmat%H4 = beckeInt%fdmat%H3
      beckeInt%fdmat%ipiv4 = beckeInt%fdmat%ipiv2
      call TBeckeIntegrator_setKernelParam(beckeInt, inp%omega)
      call TBeckeIntegrator_precompFdMatrix(beckeInt)
      call TBeckeIntegrator_buildLU(beckeInt)

    end if

    if (inp%tGlobalHybrid .or. inp%tCam) then
      call gauss_chebyshev_quadrature(inp%nRadial, radialHFQuadrature)
      ! generate spherical coordinate (r) for full-range Hartree-Fock contribution to H0
      call gengrid1_1(radialHFQuadrature, inp%rm, coordtrans_radial_becke2, rr3, dummyWeights)
    end if

    if (inp%tGlobalHybrid) then
      call initGaunt(inp%ll_max)
    end if

    call gauss_legendre_quadrature(inp%ninteg1, quads(1))
    call gauss_legendre_quadrature(inp%ninteg2, quads(2))

    atom1 => inp%atom1
    if (inp%tHetero) then
      atom2 => inp%atom2
    else
      atom2 => inp%atom1
    end if

    call TIntegMap_init(imap, atom1, atom2)

    ! calculate lines for 1 Bohr in one batch.
    dist = 0.0_dp
    tDynlen = (inp%maxdist > 0.0_dp)
    if (tDynlen) then
      nBatchline = ceiling(1.0_dp / inp%dr)
      maxdist = inp%maxdist + real(nBatchline, dp) * inp%dr
    else
      maxdist = abs(inp%maxdist)
      nBatchline = ceiling((maxdist - inp%r0) / inp%dr)
    end if
    nBatch = 0
    denserrmax = 0.0_dp
    allocate(denserr(nBatchline))
    lpBatch: do
      allocate(skhambuffer(imap%ninteg, nBatchline))
      allocate(skoverbuffer(imap%ninteg, nBatchline))
      write(*, "(A,I0,A,F6.3,A,F6.3)") "Calculating ", nBatchline, " lines: r0 = ",&
          & inp%r0 + inp%dr * real(nBatch * nBatchline, dp), " dr = ", inp%dr
      lpDist: do ir = 1, nBatchline
        dist = inp%r0 + inp%dr * real(nBatch * nBatchline + ir - 1, dp)
        write(*, "(A,F6.2,A)") 'Calculating dimer distance: ', dist, ' Bohr'
        call gengrid2_2(quads, coordtrans_becke_12, partition_becke_homo, beckepars, dist, grid1,&
            & grid2, dots, weights)
        nRad = size(quads(1)%xx)
        nAng = size(quads(2)%xx)
        call getskintegrals(beckeInt, radialHFQuadrature, nRad, nAng, atom1, atom2, grid1, grid2,&
            & rr3(:, 1), dots, weights, inp%tDensitySuperpos, inp%ll_max, inp%iXC, inp%camAlpha,&
            & inp%camBeta, inp%tGlobalHybrid, inp%tLC, inp%tCam, imap, xcfunc_xc, xcfunc_x,&
            & xcfunc_c, skhambuffer(:, ir), skoverbuffer(:, ir), denserr(ir))
      end do lpDist
      denserrmax = max(denserrmax, maxval(denserr))
      maxabs = max(maxval(abs(skhambuffer)), maxval(abs(skoverbuffer)))
      if (tDynlen) then
        tConverged = (maxabs < inp%epsilon)
        ! if new batch gave no contributions above tolerance: omit it and exit
        if (tConverged .or. dist > maxdist) exit
        nBatch = nBatch + 1
        call hamfifo%push_alloc(skhambuffer)
        call overfifo%push_alloc(skoverbuffer)
      else
        tConverged = .true.
        call hamfifo%push_alloc(skhambuffer)
        call overfifo%push_alloc(skoverbuffer)
        exit
      end if
    end do lpBatch

    if (.not. tConverged) then
      write(*, "(A,F6.2,A,ES10.3)") "Warning, maximal distance ", inp%maxdist,&
          & " reached! Max integral value:", maxabs
    end if
    write(*, "(A,ES10.3)") "Maximal integration error: ", denserrmax

    ! hand over Hamiltonian and overlap
    call hamfifo%popall_concat(skham)
    call overfifo%popall_concat(skover)

    ! finalize libxc objects
    if (inp%tGlobalHybrid .or. inp%tCam) then
      call xc_f03_func_end(xcfunc_xc)
    else
      call xc_f03_func_end(xcfunc_x)
      call xc_f03_func_end(xcfunc_c)
    end if

    ! finalize anglib (global Hybrids only)
    if (inp%tGlobalHybrid) then
      call freeGaunt()
    end if

  end subroutine get_twocenter_integrals


  !> Calculates SK-integrals.
  subroutine getskintegrals(beckeInt, radialHFQuadrature, nRad, nAng, atom1, atom2, grid1, grid2,&
      & rr3, dots, weights, tDensitySuperpos, ll_max, iXC, camAlpha, camBeta, tGlobalHybrid, tLC,&
      & tCam, imap, xcfunc_xc, xcfunc_x, xcfunc_c, skham, skover, denserr)

    !> Becke integrator instances
    type(TBeckeIntegrator), intent(inout) :: beckeInt

    !! radial full-range Hartree-Fock quadrature information
    type(TQuadrature), intent(in) :: radialHFQuadrature

    !> number of radial and angular integration abscissas
    integer, intent(in) :: nRad, nAng

    !> atomic property instances of dimer atoms
    type(TAtomdata), intent(in), pointer :: atom1, atom2

    !> integration grids of dimer atoms, holding spherical coordinates (r, theta)
    real(dp), intent(in), target :: grid1(:,:), grid2(:,:)

    !> evaluation points of radial wavefunction (full-range Hartree-Fock)
    real(dp), intent(in) :: rr3(:)

    !> dot product of unit distance vectors
    real(dp), intent(in) :: dots(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !> true, if density superposition is requested, otherwise potential superposition is applied
    logical, intent(in) :: tDensitySuperpos

    !> maximum angular momentum for Becke integration
    integer, intent(in) :: ll_max

    !> xc-functional type
    !! (1: LDA-PW91, 2: GGA-PBE, 3: GGA-BLYP, 4: LCY-PBE, 5: LCY-BNL, 6: PBE0, 7: B3LYP,
    !! 8: CAMY-B3LYP, 9: CAMY-PBEh)
    integer, intent(in) :: iXC

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> true, if a global hybrid functional is requested
    logical, intent(in) :: tGlobalHybrid

    !> true, if a (long-range corrected) range-separated hybrid functional is requested
    logical, intent(in) :: tLC

    !> true, if a CAM functional is requested
    logical, intent(in) :: tCam

    !> two-center integration mapping instance
    type(TIntegMap), intent(in) :: imap

    !! libxc related objects for exchange and correlation
    type(xc_f03_func_t), intent(in) :: xcfunc_xc, xcfunc_x, xcfunc_c

    !> resulting Hamiltonian and overlap matrix
    real(dp), intent(out) :: skham(:), skover(:)

    !> relative density integration error
    real(dp), intent(out) :: denserr

    !! instance of real tesseral spherical harmonics
    type(TRealTessY) :: tes1, tes2

    !! spherical coordinates (r, theta) of atom 1 and atom 2 on grid
    real(dp), pointer :: r1(:), r2(:), theta1(:), theta2(:)

    !! radial grid-orbital portion for all basis functions of atom 1
    real(dp), allocatable :: radval1(:,:)

    !! radial grid-orbital portion and 1st/2nd derivative for all basis functions of atom 2
    real(dp), allocatable :: radval2(:,:), radval2p(:,:), radval2pp(:,:)

    !! total potential and electron density of two atoms
    real(dp), allocatable :: potval(:), densval(:)

    !! atomic 1st density derivatives of atom 1
    real(dp), allocatable :: densval1p(:)

    !! atomic 1st density derivatives of atom 2
    real(dp), allocatable :: densval2p(:)

    !! real tesseral spherical harmonic for spherical coordinate (theta) of atom 1 and atom 2
    real(dp), allocatable :: spherval1(:), spherval2(:)

    !! temporary storage for Hamiltonian, overlap, density and pre-factors
    real(dp) :: integ1, integ2, dens, prefac

    !! full-/long-range Hartree-Fock exchange contribution
    real(dp) :: frx, lrx

    !! number of integration points
    integer :: nGrid

    !> number of density grid points, compatible with libxc signature
    integer(c_size_t) :: nGridLibxc

    !!  orbital indices/angular momenta on the two atoms and interaction type
    integer :: i1, i2, l1, l2, mm

    !! integral index
    integer :: ii

    !! libxc related objects
    real(dp), allocatable :: vxc(:), vx(:), vc(:)
    real(dp), allocatable :: rhor(:), sigma(:), vxcsigma(:), vxsigma(:), vcsigma(:)
    real(dp), allocatable :: divvxc(:), divvx(:), divvc(:)

    r1 => grid1(:, 1)
    theta1 => grid1(:, 2)
    r2 => grid2(:, 1)
    theta2 => grid2(:, 2)
    nGrid = size(r1)
    nGridLibxc = nGrid

    allocate(radval1(nGrid, atom1%nbasis))
    allocate(radval2(nGrid, atom2%nbasis))
    allocate(radval2p(nGrid, atom2%nbasis))
    allocate(radval2pp(nGrid, atom2%nbasis))
    allocate(spherval1(nGrid))
    allocate(spherval2(nGrid))

    ! get radial portions of all basis functions of atom 1
    do ii = 1, size(radval1, dim=2)
      radval1(:, ii) = atom1%rad(ii)%getValue(r1)
    end do

    ! get radial portions (and derivatives) of all basis functions of atom 2
    do ii = 1, size(radval2, dim=2)
      radval2(:, ii) = atom2%rad(ii)%getValue(r2)
      radval2p(:, ii) = atom2%drad(ii)%getValue(r2)
      radval2pp(:, ii) = atom2%ddrad(ii)%getValue(r2)
    end do

    ifPotSup: if (.not. tDensitySuperpos) then
      potval = atom1%pot%getValue(r1) + atom2%pot%getValue(r2)
    else
      allocate(densval(nGrid))
      densval(:) = atom1%rho%getValue(r1) + atom2%rho%getValue(r2)

      ! prepare xc-functional specific arrays
      ! care about correct 4pi normalization of density
      rhor = getLibxcRho(densval)
      allocate(vx(nGrid))
      allocate(vc(nGrid))
      if (iXC /= 1) then
        allocate(vxsigma(nGrid))
        allocate(vcsigma(nGrid))
        densval1p = atom1%drho%getValue(r1)
        densval2p = atom2%drho%getValue(r2)
        ! care about correct 4pi normalization of density and compute sigma
        sigma = getLibxcSigma(densval1p, densval2p, dots)
      end if
      if (tGlobalHybrid .or. tCam) then
        allocate(vxc(nGrid))
        allocate(vxcsigma(nGrid))
      end if

      select case (iXC)
      ! 1: LDA-PW91
      case(1)
        call xc_f03_lda_vxc(xcfunc_x, nGridLibxc, rhor(1), vx(1))
        call xc_f03_lda_vxc(xcfunc_c, nGridLibxc, rhor(1), vc(1))
        potval = vx + vc
      ! 2: GGA-PBE, 3: GGA-BLYP, 4: LCY-PBE
      case(2:4)
        call xc_f03_gga_vxc(xcfunc_x, nGridLibxc, rhor(1), sigma(1), vx(1), vxsigma(1))
        call xc_f03_gga_vxc(xcfunc_c, nGridLibxc, rhor(1), sigma(1), vc(1), vcsigma(1))
        call getDivergence(nRad, nAng, densval1p, densval2p, r1, r2, theta1, theta2, vxsigma, divvx)
        call getDivergence(nRad, nAng, densval1p, densval2p, r1, r2, theta1, theta2, vcsigma, divvc)
        potval = vx + vc + divvx + divvc
      ! 5: LCY-BNL
      case(5)
        call xc_f03_lda_vxc(xcfunc_x, nGridLibxc, rhor(1), vx(1))
        call xc_f03_gga_vxc(xcfunc_c, nGridLibxc, rhor(1), sigma(1), vc(1), vcsigma(1))
        call getDivergence(nRad, nAng, densval1p, densval2p, r1, r2, theta1, theta2, vcsigma, divvc)
        potval = vx + vc + divvc
      ! 6: PBE0, 7: B3LYP, 8: CAMY-B3LYP, 9: CAMY-PBEh
      case(6:9)
        call xc_f03_gga_vxc(xcfunc_xc, nGridLibxc, rhor(1), sigma(1), vxc(1), vxcsigma(1))
        call getDivergence(nRad, nAng, densval1p, densval2p, r1, r2, theta1, theta2, vxcsigma,&
            & divvxc)
        potval = vxc + divvxc
      end select
      ! add nuclear and coulomb potential to obtain the effective potential
      potval(:) = potval + atom1%pot%getValue(r1) + atom2%pot%getValue(r2)
    end if ifPotSup

    denserr = 0.0_dp
    do ii = 1, imap%ninteg
      i1 = imap%type(1, ii)
      l1 = atom1%angmoms(i1)
      i2 = imap%type(2, ii)
      l2 = atom2%angmoms(i2)
      mm = imap%type(3, ii) - 1
      call TRealTessY_init(tes1, l1, mm)
      call TRealTessY_init(tes2, l2, mm)

      ! Y_{l1 mm}(\theta_1, \phi = 0)
      spherval1(:) = tes1%getValue_1d(theta1)
      ! Y_{l2 mm}(\theta_2, \phi = 0)
      spherval2(:) = tes2%getValue_1d(theta2)

      ! calculate SK-quantities
      ! Hamiltonian
      integ1 = getHamiltonian(radval1(:, i1), radval2(:, i2), radval2p(:, i2), radval2pp(:, i2),&
          & r2, l2, spherval1, spherval2, potval, weights)
      ! overlap integral: \sum_{r,\Omega} R_1(r) Y_1(\Omega) R_2(r) Y_2(\Omega) weight
      integ2 = getOverlap(radval1(:, i1), radval2(:, i2), spherval1, spherval2, weights)
      ! total density: \int (|\phi_1|^2 + |\phi_2|^2)
      dens = getDensity(radval1(:, i1), radval2(:, i2), spherval1, spherval2, weights)

      if (tGlobalHybrid) then
        ! full-range Hartree-Fock exchange contribution
        frx = 0.5_dp * getFullRangeHFContribution(radialHFQuadrature%xx, rr3, ll_max, atom1, atom2,&
            & imap, ii, r1, theta1, r2, theta2, weights)
        ! add up full-range exchange to the Hamiltonian
        integ1 = integ1 - frx
      elseif (tLC) then
        ! long-range Hartree-Fock exchange contribution
        lrx = 0.5_dp * getLongRangeHFContribution(beckeInt, atom1, atom2, imap, ii, r1, theta1, r2,&
            & theta2, weights)
        ! add up long-range exchange to the Hamiltonian
        integ1 = integ1 - lrx
      elseif (tCam) then
        ! full-range Hartree-Fock exchange contribution
        frx = 0.5_dp * camAlpha * getFullRangeHFContribution(radialHFQuadrature%xx, rr3, ll_max,&
            & atom1, atom2, imap, ii, r1, theta1, r2, theta2, weights)
        ! long-range Hartree-Fock exchange contribution
        lrx = 0.5_dp * camBeta * getLongRangeHFContribution(beckeInt, atom1, atom2, imap, ii, r1,&
            & theta1, r2, theta2, weights)
        ! add up full-/long-range exchange to the Hamiltonian
        integ1 = integ1 - frx - lrx
      end if

      if (mm == 0) then
        prefac = 2.0_dp * pi
      else
        prefac = pi
      end if
      skham(ii) = prefac * integ1
      skover(ii) = prefac * integ2
      dens = prefac * dens
      denserr = max(denserr, abs(dens - 2.0_dp) / 2.0_dp)
    end do

  end subroutine getskintegrals


  !> Calculates libXC renormalized density superposition of dimer.
  pure function getLibxcRho(rho) result(rhor)

    !> superposition of atomic densities of atom 1 and atom 2
    real(dp), intent(in) :: rho(:)

    !> renormalized density
    real(dp), allocatable :: rhor(:)

    ! renorm rho (incoming quantities are 4pi normed)
    rhor = rho * rec4pi

  end function getLibxcRho


  !> Calculates libXC sigma of dimer.
  pure function getLibxcSigma(drho1, drho2, dots) result(sigma)

    !> 1st derivative of atomic densities
    real(dp), intent(in) :: drho1(:), drho2(:)

    !> dot product of unit distance vectors
    real(dp), intent(in) :: dots(:)

    !! libXC's contracted gradients of the density
    real(dp), allocatable :: sigma(:)

    !! number of tabulated grid points
    integer :: nn

    !! recurring factors
    real(dp), allocatable :: f1(:), f2(:)

    nn = size(drho1)

    f1 = drho1 + dots * drho2
    f2 = drho2 + dots * drho1

    ! get dot product of density gradients
    sigma = (drho1 * f1 + drho2 * f2) * rec4pi**2

  end function getLibxcSigma


  !> Computes contribution div(v) to the xc-potential due to vsigma = deps/dsigma returned by libxc.
  !! div(v) = -2 div(vsigma grad(n)), evaluated in spherical coordinates r1, theta1 as:
  !! div(v) = -2 [(1/r1^2) d(r1^2 vsigma grad(n)_r1)/dr1 +
  !!             (1/r1*sin(theta1)) d(sin(theta1 vsigma grad(n)_theta1)/dtheta1]
  !! where grad(n) = (drho1(r1) + cos(theta2-theta1) drho2(r2)) \vec{e}_r1 +
  !!                  sin(theta2-theta1) drho2(r2) \vec{e}_theta1
  subroutine getDivergence(nRad, nAng, drho1, drho2, r1, r2, theta1, theta2, vsigma, divv)

    !> # radial points of grid => r1
    integer, intent(in)  :: nRad

    !> # angular points of grid => theta1
    integer, intent(in)  :: nAng

    !> radial derivative density atom1 on grid (nRad, nAng)
    real(dp), intent(in) :: drho1(:)

    !> radial derivative density atom2 on grid (nRad, nAng)
    real(dp), intent(in) :: drho2(:)

    !> values of r1 on grid points (for index <= nn), values of r2b (for index > nn)
    real(dp), intent(in) :: r1(:)

    !> values of r2 on grid points (for index <= nn), values of r1 (for index > nn)
    real(dp), intent(in) :: r2(:)

    !> values of theta1 on grid points (for index <= nn), values of theta2b (for index > nn)
    real(dp), intent(in) :: theta1(:)

    !> values of theta2 on grid points (for index <= nn), values of theta1 (index > nn)
    real(dp), intent(in) :: theta2(:)

    !> deps/dsigma returned by libxc on grid
    real(dp), intent(in) :: vsigma(:)

    !> -2 div(vsigma grad(n)) on grid
    real(dp), intent(out), allocatable :: divv(:)

    !!
    integer :: nn, ia, ir

    !! radii and theta values of the grid
    real(dp), allocatable :: rval(:), tval(:)

    !!
    real(dp), allocatable :: aa(:,:), dar(:), dat(:), bb(:,:)

    allocate(divv(size(drho1)))
    divv(:) = 0.0_dp
    nn = size(drho1) / 2

    allocate(dar(nRad), dat(nAng), rval(nRad), bb(nRad, nAng))

    ! rval holds radii of the grid
    rval = r1(1:nRad)
    aa = reshape(theta1(1:nn), [nRad, nAng])

    !! tval holds theta values of the grid
    tval = aa(1, :)

    if ((rval(2) < rval(1)) .or. (tval(2) > tval(1))) then
      write(*,*) 'getDivergence: Expected ascending order in radii and descending order in theta!'
      stop
    end if

    ! div n = (drho1(r1) + cos(theta2-theta1) drho2(r2)) \vec{e}_r1
    !         + sin(theta2-theta1) drho2(r2) \vec{e}_theta1
    ! elements <= nn refer to: atom1(r1) -- atom2(r2a)
    aa(:,:) = reshape(rec4pi * vsigma(1:nn)&
        & * (drho1(1:nn) + cos(theta2(1:nn) - theta1(1:nn)) * drho2(1:nn)) * r1(1:nn)**2,&
        & [nRad, nAng])

    ! take numerical derivative w.r.t. r1
    do ia = 1, nAng
      call spline3ders(rval, aa(:, ia), rval, dynew=dar)
      bb(:, ia) = dar
    end do
    divv(1:nn) = reshape(bb, [nRad * nAng]) / r1(1:nn)**2

    !! elements > nn refer to: atom1(r2b) -- atom2(r1)
    aa(:,:) = reshape(rec4pi * vsigma(nn+1:2*nn)&
        & * (drho2(nn+1:2*nn) + cos(theta1(nn+1:2*nn) - theta2(nn+1:2*nn)) * drho1(nn+1:2*nn))&
        & * r2(nn+1:2*nn)**2, [nRad, nAng])

    do ia = 1, nAng
      call spline3ders(rval, aa(:, ia), rval, dynew=dar)
      bb(:, ia) = dar
    end do
    divv(nn+1:2*nn) = reshape(bb, [nRad * nAng]) / r2(nn+1:2*nn)**2

    ! take numerical derivative w.r.t. theta1
    aa(:,:) = reshape(rec4pi * vsigma(1:nn)&
        & * (sin(theta2(1:nn) - theta1(1:nn)) * drho2(1:nn)) * sin(theta1(1:nn)),&
        & [nRad, nAng])

    !! spline3der requires data in ascending order
    do ir = 1, nRad
      call spline3ders(tval(nAng:1:-1), aa(ir,nAng:1:-1), tval(nAng:1:-1), dynew=dat)
      bb(ir, :) = dat(nAng:1:-1)
    end do
    divv(1:nn) = divv(1:nn) + reshape(bb, [nRad * nAng]) / (r1(1:nn) * sin(theta1(1:nn)))

    aa(:,:) = reshape(rec4pi * vsigma(nn+1:2*nn)&
        & * (sin(theta1(nn+1:2*nn) - theta2(nn+1:2*nn)) * drho1(nn+1:2*nn)) * sin(theta1(1:nn)),&
        & [nRad, nAng])

    do ir = 1, nRad
      call spline3ders(tval(nAng:1:-1), aa(ir,nAng:1:-1), tval(nAng:1:-1), dynew=dat)
      bb(ir, :) = dat(nAng:1:-1)
    end do
    divv(nn+1:2*nn) = divv(nn+1:2*nn) + reshape(bb, [nRad * nAng]) / (r1(1:nn) * sin(theta1(1:nn)))

    ! pre-factor
    divv(:) = -2.0_dp * divv

  end subroutine getDivergence


  !> Calculates overlap for a fixed orbital and interaction configuration.
  pure function getOverlap(rad1, rad2, spher1, spher2, weights) result(res)

    !> radial grid-orbital portion of atom 1 and atom 2
    real(dp), intent(in) :: rad1(:), rad2(:)

    !> real tesseral spherical harmonic for spherical coordinate (theta) of atom 1 and atom 2
    real(dp), intent(in) :: spher1(:), spher2(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !! resulting orbital overlap
    real(dp) :: res

    res = sum(rad1 * rad2 * spher1 * spher2 * weights)

  end function getOverlap


  !> Calculates density for a fixed orbital and interaction configuration.
  pure function getDensity(rad1, rad2, spher1, spher2, weights) result(res)

    !> radial grid-orbital portion of atom 1 and atom 2
    real(dp), intent(in) :: rad1(:), rad2(:)

    !> real tesseral spherical harmonic for spherical coordinate (theta) of atom 1 and atom 2
    real(dp), intent(in) :: spher1(:), spher2(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !! resulting electron density
    real(dp) :: res

    res = sum(((rad1 * spher1)**2 + (rad2 * spher2)**2) * weights)

  end function getdensity


  !> Calculates Hamiltonian for a fixed orbital and interaction configuration.
  pure function getHamiltonian(rad1, rad2, rad2p, rad2pp, r2, l2, spher1, spher2, pot, weights)&
      & result(res)

    !> radial grid-orbital portion of atom 1 and atom 2
    real(dp), intent(in) :: rad1(:), rad2(:)

    !> radial grid-orbital portion's 1st and 2nd derivative of atom 2
    real(dp), intent(in) :: rad2p(:), rad2pp(:)

    !> radial spherical coordinates of atom 2 on grid
    real(dp), intent(in) :: r2(:)

    !> angular momentum corresponding to current orbital index of atom 2
    integer, intent(in) :: l2

    !> real tesseral spherical harmonic for spherical coordinate (theta) of atom 1 and atom 2
    real(dp), intent(in) :: spher1(:), spher2(:)

    !> effective potential on grid
    real(dp), intent(in) :: pot(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !! resulting Hamiltonian matrix element
    real(dp) :: res

    res = sum((rad1 * spher1)&
        & * (- 0.5_dp * rad2pp&
        & - rad2p / r2&
        & + 0.5_dp * l2 * (l2 + 1) * rad2 / r2**2&
        & + pot * rad2)&
        & * spher2 * weights)

  end function getHamiltonian


  !> Calculates full-range HF contribution to H0, i.e. two-center integral over Coulomb kernel.
  function getFullRangeHFContribution(radialQuadratureXx, rr3, ll_max, atom1, atom2, imap, iInt,&
      & rr1, theta1, rr2, theta2, weights) result(frContribution)

    !> radial quadrature abscissae
    real(dp), intent(in) :: radialQuadratureXx(:)

    !> evaluation points of radial wavefunction
    real(dp), intent(in) :: rr3(:)

    !> maximum angular momentum for Becke integration
    integer, intent(in) :: ll_max

    !> atomic property instances of dimer atoms
    type(TAtomdata), intent(in) :: atom1, atom2

    !> integral mapping instance
    type(TIntegMap), intent(in) :: imap

    !> current integral index
    integer, intent(in) :: iInt

    !! spherical coordinates (r, theta) of atom 1 and atom 2 on grid
    real(dp), intent(in) :: rr1(:), rr2(:), theta1(:), theta2(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !> FR contribution to the Hamiltonian H0
    real(dp) :: frContribution

    !! number of radial evaluation and grid points
    integer :: nGrid, nRadial

    !! interpolated (cubic splines) rho_lm at back-transformed coordinates
    real(dp) :: yy2

    !! T-symbol, see Vitalij's thesis
    real(dp) :: tsymbol

    !! back-transformed, radial spherical coordinates of atom 1 and atom 2
    real(dp), allocatable :: tmp11(:), tmp22(:)

    !! radial (core) parts of the inner integrand
    real(dp), allocatable :: rrin(:), rrin2(:)

    !! solution (and temporary storage) of inner integrands
    real(dp), allocatable :: V(:), V_sum(:)

    !! scaled weight/abscissa index i / (n + 1)
    real(dp), allocatable :: zi(:)

    !! l-resolved integral solution on fdiff nodes (+ 2nd derivatives) and lm-resolved density
    real(dp), allocatable :: V_l(:,:,:), rho_lm(:)

    !! 2nd derivatives of natural splines
    real(dp), allocatable :: secondDerivs(:)

    !! angular and magnetic momenta
    integer :: ll, ll_a, ll_nu, mm_nu, ll_mu, mm_mu

    !! instance of real tesseral spherical harmonics
    type(TRealTessY) :: tess

    !! grid point and core iterator
    integer :: iGridPt, iCore

    ! \int phi^B_a(r2) * phi^A_b(r2) \int phi^A_b(r1) phi^A_c(r1) / |r2-r1|

    nGrid = size(rr2)
    nRadial = size(rr3)

    allocate(V(nGrid))
    allocate(rrin(nRadial))
    allocate(rrin2(nRadial))
    allocate(rho_lm(nRadial))
    allocate(zi(nRadial))
    allocate(V_l(nRadial, ll_max, 2))

    allocate(tmp22(nGrid))
    allocate(tmp11(nGrid))
    tmp22(:) = acos((rr2 - 1.0_dp) / (rr2 + 1.0_dp)) / pi
    tmp11(:) = acos((rr1 - 1.0_dp) / (rr1 + 1.0_dp)) / pi

    ! nu-orbital
    ll_nu = atom1%angmoms(imap%type(1, iInt))
    mm_nu = imap%type(3, iInt) - 1
    rrin2(:) = atom1%rad(imap%type(1, iInt))%getValue(rr3)

    zi(:) = acos(radialQuadratureXx) / pi

    allocate(V_sum(nGrid))
    V_sum(:) = 0.0_dp

    ! start the sigma loop for all core electrons of atom1
    do iCore = 1, atom1%nCore

      ! evaluate the radial part of the inner integrand
      rrin(:) = rrin2 * atom1%corerad(iCore)%getValue(rr3)
      ll_a = atom1%coreAngmoms(iCore)

      ! solve the inner integral
      V(:) = 0.0_dp
      do ll = 0, ll_max - 1
        tsymbol = getTsymbol(ll_a, ll, ll_nu, mm_nu) * atom1%coreOcc(iCore)
        if (tsymbol >= 1.0e-16_dp) then
          rho_lm(:) = rrin
          ! solve the equation for ll
          call solvePoisson(ll, rho_lm, zi, 0.0_dp)
          call get2ndNaturalSplineDerivs(zi, rho_lm, secondDerivs)
          V_l(:, ll+1, 1) = rho_lm
          V_l(:, ll+1, 2) = secondDerivs
          do iGridPt = 1, nGrid
            call get_cubic_spline(zi, V_l(:, ll+1, 1), V_l(:, ll+1, 2), tmp11(iGridPt), yy2)
            V(iGridPt) = V(iGridPt) + yy2 * tsymbol
          end do
        end if
      end do
      ! end of evaluation of the inner integral

      ! radial part of sigma
      V(:) = V * atom1%corerad(iCore)%getValue(rr1)
      V_sum(:) = V_sum + V
    end do

    ! end sigma loop
    V(:) = V_sum

    ! angular part Y(ll_nu, mm_mu)
    ll_mu = atom2%angmoms(imap%type(2, iInt))
    mm_mu = imap%type(3, iInt) - 1
    call TRealTessY_init(tess, ll_nu, mm_mu)
    V(:) = V / rr1 * tess%getValue_1d(theta1)

    ! mu-orbital
    call TRealTessY_init(tess, ll_mu, mm_mu)
    V(:) = V * atom2%rad(imap%type(2, iInt))%getValue(rr2) * tess%getValue_1d(theta2)

    frContribution = sum(V * weights)

    ! nu-orbital -> mu-orbital
    ll_nu = atom2%angmoms(imap%type(2, iInt))
    mm_nu = imap%type(3, iInt) - 1
    rrin2(:) = atom2%rad(imap%type(2, iInt))%getValue(rr3)

    V_sum(:) = 0.0_dp

    ! start the sigma loop for all core electrons of atom2
    do iCore = 1, atom2%nCore

      ! evaluate the radial part of the inner integrand
      rrin(:) = rrin2 * atom2%corerad(iCore)%getValue(rr3)
      ll_a = atom2%coreAngmoms(iCore)
      zi(:) = acos(radialQuadratureXx) / pi

      ! solve the inner integral
      V(:) = 0.0_dp
      do ll = 0, ll_max - 1
        tsymbol = getTsymbol(ll_a, ll, ll_nu, mm_nu) * atom2%coreOcc(iCore)
        if (tsymbol >= 1.0e-16_dp) then
          rho_lm(:) = rrin
          ! solve the equation for ll
          call solvePoisson(ll, rho_lm, zi, 0.0_dp)
          call get2ndNaturalSplineDerivs(zi, rho_lm, secondDerivs)
          V_l(:, ll+1, 1) = rho_lm
          V_l(:, ll+1, 2) = secondDerivs
          do iGridPt = 1, nGrid
            call get_cubic_spline(zi, V_l(:, ll+1, 1), V_l(:, ll+1, 2), tmp22(iGridPt), yy2)
            V(iGridPt) = V(iGridPt) + yy2 * tsymbol
          end do
        end if
      end do
      ! end of evaluation of the inner integral

      ! radial part of sigma
      V(:) = V * atom2%corerad(iCore)%getValue(rr2)
      V_sum(:) = V_sum + V
    end do

    ! end sigma loop
    V(:) = V_sum

    ! angular part Y(ll_nu, mm_mu)
    ll_mu = atom1%angmoms(imap%type(1, iInt))
    mm_mu = imap%type(3, iInt) - 1
    call TRealTessY_init(tess, ll_nu, mm_mu)
    V(:) = V / rr2 * tess%getValue_1d(theta2)

    ! mu-orbital
    call TRealTessY_init(tess, ll_mu, mm_mu)
    V(:) = V * atom1%rad(imap%type(1, iInt))%getValue(rr1) * tess%getValue_1d(theta1)

    frContribution = frContribution + sum(V * weights)

  end function getFullRangeHFContribution


  !> Calculates long-range HF contribution to H0, i.e. two-center integral over descreened Yukawa.
  function getLongRangeHFContribution(beckeInt, atom1, atom2, imap, iInt, rr1, theta1, rr2, theta2,&
      & weights) result(lcContribution)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: beckeInt

    !> atomic property instances of dimer atoms
    type(TAtomdata), intent(in) :: atom1, atom2

    !> integral mapping instance
    type(TIntegMap), intent(in) :: imap

    !> current integral index
    integer, intent(in) :: iInt

    !! spherical coordinates (r, theta) of atom 1 and atom 2 on grid
    real(dp), intent(in) :: rr1(:), rr2(:), theta1(:), theta2(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !> LC contribution to the Hamiltonian H0
    real(dp) :: lcContribution

    !!
    real(dp), pointer :: rr3(:)

    !!
    integer :: nGrid, nRadial, iGridPt, iCore

    !!
    real(dp) :: yy2, tsymbol

    !!
    real(dp), allocatable :: V(:), tmp11(:), tmp22(:), V_sum(:), AA2(:), rrin(:), rrin2(:)

    !!
    real(dp), allocatable :: CCC(:)

    !!
    real(dp), allocatable :: zi(:), V_l(:,:,:), rho_lm(:), res(:), rho_lm2(:), rho_lm3(:)

    !!
    integer :: ll, ll_a, ll_nu, mm_nu, ll_mu, mm_mu

    !!
    real(dp) :: charge

    !!
    type(TRealTessY) :: tess

    ! \int phi^B_a(r2) * phi^A_b(r2) \int phi^A_b(r1) phi^A_c(r1) / |r2-r1|

    ! get the coordinates for the inner one-center integral
    call TBeckeIntegrator_getCoords(beckeInt, [3, 1, 1], rr3)

    nGrid = size(rr2)
    nRadial = size(rr3)

    allocate(V(nGrid))
    allocate(V_l(nRadial, beckeInt%beckeGridParams%ll_max, 2))

    allocate(tmp22(nGrid))
    allocate(tmp11(nGrid))
    tmp22(:) = acos((rr2 - 1.0_dp) / (rr2 + 1.0_dp)) / pi
    tmp11(:) = acos((rr1 - 1.0_dp) / (rr1 + 1.0_dp)) / pi

    allocate(rrin(nRadial))
    allocate(rrin2(nRadial))
    allocate(rho_lm(nRadial))
    allocate(rho_lm2(nRadial))
    allocate(rho_lm3(nRadial))
    allocate(zi(nRadial))

    ! nu-orbital
    ll_nu = atom1%angmoms(imap%type(1, iInt))
    mm_nu = imap%type(3, iInt) - 1
    rrin2(:) = atom1%rad(imap%type(1, iInt))%getValue(rr3)

    zi(:) = acos(beckeInt%radialQuadrature%xx) / pi

    allocate(V_sum(nGrid))
    allocate(AA2(nGrid))
    allocate(CCC(nRadial))
    V_sum(:) = 0.0_dp
    AA2(:) = 0.0_dp
    CCC(:) = 0.0_dp

    ! start the sigma loop for all core electrons of atom1
    do iCore = 1, atom1%nCore

      ! evaluate the radial part of the inner integrand
      rrin(:) = rrin2 * atom1%corerad(iCore)%getValue(rr3)
      ll_a = atom1%coreAngmoms(iCore)

      ! evaluate the correction to the poisson solver, necessary only for small omega
      rho_lm3(:) = rrin
      ll = 0
      if (ll_a == ll_nu) then
        rho_lm3(:) = rrin
        charge = sum(rho_lm3 * beckeInt%beckeGrid(3)%weight)
        call solvePoisson(ll, rho_lm3, zi, charge)
        CCC(:) = rho_lm3
        rho_lm3(:) = rrin
        charge = 0.0_dp
        call solvePoisson(ll, rho_lm3, zi, charge)
        CCC(:) = (CCC - rho_lm3) / sqrt(4.0_dp * pi)

        call get2ndNaturalSplineDerivs(zi, CCC, res)
        V_l(:, ll+1, 1) = CCC
        V_l(:, ll+1, 2) = res
        do iGridPt = 1, nGrid
          call get_cubic_spline(zi, V_l(:, ll+1, 1), V_l(:, ll+1, 2), tmp11(iGridPt), yy2)
          AA2(iGridPt) = AA2(iGridPt) + yy2 / rr1(iGridPt)
        end do
      end if

      ! solve the inner integral
      V(:) = 0.0_dp
      do ll = 0, beckeInt%beckeGridParams%ll_max - 1
        tsymbol = getTsymbol(ll_a, ll, ll_nu, mm_nu) * atom1%coreOcc(iCore)
        if (tsymbol >= 1.0e-16_dp) then
          rho_lm(:) = rrin
          rho_lm2(:) = rrin
          ! solve the equation for ll
          charge = 0.0_dp
          call solvePoisson(ll, rho_lm2, zi, charge)
          call TBeckeIntegrator_solveHelmholz(beckeInt, ll, rho_lm)
          ! interpolate
          rho_lm(:) = rho_lm2 - rho_lm
          call get2ndNaturalSplineDerivs(zi, rho_lm, res)
          V_l(:, ll+1, 1) = rho_lm
          V_l(:, ll+1, 2) = res
          do iGridPt = 1, nGrid
            call get_cubic_spline(zi, V_l(:, ll+1, 1), V_l(:, ll+1, 2), tmp11(iGridPt), yy2)
            V(iGridPt) = V(iGridPt) + yy2 * tsymbol
          end do
        end if
      end do
      ! end of evaluation of the inner integral

      ! radial part of sigma
      V(:) = V * atom1%corerad(iCore)%getValue(rr1)
      V_sum(:) = V_sum + V
    end do

    ! end sigma loop
    V(:) = V_sum

    ! angular part Y(ll_nu, mm_mu)
    ll_mu = atom2%angmoms(imap%type(2, iInt))
    mm_mu = imap%type(3, iInt) - 1
    call TRealTessY_init(tess, ll_nu, mm_mu)
    V(:) = V / rr1 * tess%getValue_1d(theta1)

    ! evaluate the correction term
    AA2(:) = AA2 * atom1%rad(imap%type(1, iInt))%getValue(rr1) * tess%getValue_1d(theta1)

    ! mu-orbital
    call TRealTessY_init(tess, ll_mu, mm_mu)
    V(:) = V * atom2%rad(imap%type(2, iInt))%getValue(rr2) * tess%getValue_1d(theta2)

    ! evaluate the correction term
    AA2(:) = AA2 * atom2%rad(imap%type(2, iInt))%getValue(rr2) * tess%getValue_1d(theta2)

    lcContribution = sum(V * weights) + sum(AA2 * weights)

    ! nu-orbital -> mu_orbital
    ll_nu = atom2%angmoms(imap%type(2, iInt))
    mm_nu = imap%type(3, iInt) - 1
    rrin2(:) = atom2%rad(imap%type(2, iInt))%getValue(rr3)

    AA2(:) = 0.0_dp
    V_sum(:) = 0.0_dp

    ! start the sigma loop for all core electrons of atom2
    do iCore = 1, atom2%nCore
      ! evaluate the radial part of the inner integrand
      rrin(:) = rrin2 * atom2%corerad(iCore)%getValue(rr3)
      ll_a = atom2%coreAngmoms(iCore)
      zi(:) = acos(beckeInt%radialQuadrature%xx) / pi

      ! evaluate the correction
      rho_lm3(:) = rrin
      ll = 0
      if (ll_a == ll_nu) then
        rho_lm3(:) = rrin
        charge = sum(rho_lm3 * beckeInt%beckeGrid(3)%weight)
        call solvePoisson(ll, rho_lm3, zi, charge)
        CCC(:) = rho_lm3
        rho_lm3(:) = rrin
        charge = 0.0_dp
        call solvePoisson(ll, rho_lm3, zi, charge)
        CCC(:) = (CCC - rho_lm3) / sqrt(4.0_dp * pi)

        call get2ndNaturalSplineDerivs(zi, CCC, res)
        V_l(:, ll+1, 1) = CCC
        V_l(:, ll+1, 2) = res
        do iGridPt = 1, nGrid
          call get_cubic_spline(zi, V_l(:, ll+1, 1), V_l(:, ll+1, 2), tmp22(iGridPt), yy2)
          AA2(iGridPt) = AA2(iGridPt) + yy2 / rr2(iGridPt)
        end do
      end if

      ! solve the inner integral
      V = 0.0_dp
      do ll = 0, beckeInt%beckeGridParams%ll_max - 1
        tsymbol = getTsymbol(ll_a, ll, ll_nu, mm_nu) * atom2%coreOcc(iCore)
        if (tsymbol >= 1.0e-16_dp) then
          rho_lm(:) = rrin
          rho_lm2(:) = rrin

          ! solve the equation for ll
          charge = 0.0_dp
          call solvePoisson(ll, rho_lm2, zi, charge)
          call TBeckeIntegrator_solveHelmholz(beckeInt, ll, rho_lm)

          ! interpolate
          rho_lm(:) = rho_lm2 - rho_lm

          call get2ndNaturalSplineDerivs(zi, rho_lm, res)
          V_l(:, ll+1, 1) = rho_lm
          V_l(:, ll+1, 2) = res
          do iGridPt = 1, nGrid
            call get_cubic_spline(zi, V_l(:, ll+1, 1), V_l(:, ll+1, 2), tmp22(iGridPt), yy2)
            V(iGridPt) = V(iGridPt) + yy2 * tsymbol
          end do
        end if
      end do
      ! end of evaluation of the inner integral

      ! radial part of sigma
      V(:) = V * atom2%corerad(iCore)%getValue(rr2)
      V_sum(:) = V_sum + V
    end do

    ! end sigma loop
    V(:) = V_sum

    ! angular part Y(ll_nu, mm_mu)
    ll_mu = atom1%angmoms(imap%type(1, iInt))
    mm_mu = imap%type(3, iInt) - 1
    call TRealTessY_init(tess, ll_nu, mm_mu)
    V(:) = V / rr2 * tess%getValue_1d(theta2)

    ! evaluate the correction term
    AA2(:) = AA2 * atom2%rad(imap%type(2,iInt))%getValue(rr2) * tess%getValue_1d(theta2)

    ! mu-orbital
    call TRealTessY_init(tess, ll_mu, mm_mu)
    V(:) = V * atom1%rad(imap%type(1, iInt))%getValue(rr1) * tess%getValue_1d(theta1)

    ! evaluate the correction term
    AA2(:) = AA2 * atom1%rad(imap%type(1, iInt))%getValue(rr1) * tess%getValue_1d(theta1)

    lcContribution = lcContribution + sum(V * weights) + sum(AA2 * weights)

  end function getLongRangeHFContribution


  !> Auxiliary function to evaluate the t-symbol.
  function getTsymbol(ll1, ll2, ll3, mm3) result(res)

    !>
    integer, intent(in) :: ll1, ll2, ll3, mm3

    !>
    real(dp) :: res

    !!
    integer :: mm1, mm2

    res = 0.0_dp
    do mm1 = -ll1, ll1
      do mm2 = -ll2, ll2
        res = res + realGaunt(ll1, mm1, ll2, mm2, ll3, mm3)**2
      end do
    end do

    if (ll1 == 1) then
      res = res * (1.0_dp / 3.0_dp)
    else
      res = 1.0_dp * res
    end if

  end function getTsymbol


  !> Initializes the two-center integration map based on the basis on two atoms.
  subroutine TIntegMap_init(this, atom1, atom2)

    !> two-center integration mapping instance
    type(TIntegMap), intent(out) :: this

    !> atomic property instances of dimer atoms
    type(TAtomdata), intent(in) :: atom1, atom2

    !! number of all nonzero two-center integrals between orbitals of two atoms
    integer :: ninteg

    !! maximum mutual angular momentum
    integer :: mmax

    !!  orbital indices/angular momenta on the two atoms and interaction type
    integer :: i1, i2, l1, l2, mm

    !! auxiliary variable
    integer :: ind

    mmax = min(maxval(atom1%angmoms), maxval(atom2%angmoms))
    allocate(this%index(atom1%nbasis, atom2%nbasis, mmax + 1))
    this%index = 0
    ninteg = 0
    do i1 = 1, atom1%nbasis
      l1 = atom1%angmoms(i1)
      do i2 = 1, atom2%nbasis
        l2 = atom2%angmoms(i2)
        do mm = 0, min(l1, l2)
          ninteg = ninteg + 1
          this%index(i1, i2, mm + 1) = ninteg
        end do
      end do
    end do
    this%ninteg = ninteg
    allocate(this%type(3, ninteg))
    ind = 0
    do i1 = 1, atom1%nbasis
      l1 = atom1%angmoms(i1)
      do i2 = 1, atom2%nbasis
        l2 = atom2%angmoms(i2)
        do mm = 1, min(l1, l2) + 1
          ind = ind + 1
          this%type(:, ind) = [i1, i2, mm]
        end do
      end do
    end do

  end subroutine TIntegMap_init

end module twocnt
