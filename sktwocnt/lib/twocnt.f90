!> Module that contains the two-center integrator routines for tabulating Hamiltonian and overlap.
module twocnt

  use common_accuracy, only : dp
  use common_constants, only : pi
  use coordtrans, only : coordtrans_becke_12
  use gridorbital, only : TGridorb2
  use sphericalharmonics, only : TRealTessY, TRealTessY_init
  use quadratures, only : TQuadrature, gauss_legendre_quadrature
  use gridgenerator, only : gengrid2_12
  use partition, only : partition_becke_homo
  use dftxc, only : getxcpot_ldapw91, getxcpot_ggapbe
  use common_fifo, only : TFiFoReal2

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

    !> xc-functional type (0: potential superposition, 1: LDA-PW91, 2: GGA-PBE)
    integer :: ixc

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

    !! integration grids of dimer atoms, holding spherical coordinates (r, theta)
    real(dp), allocatable :: grid1(:,:), grid2(:,:)

    !! ??? and integration weights
    real(dp), allocatable :: dots(:), weights(:)

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
    do
      allocate(skhambuffer(imap%ninteg, nBatchline))
      allocate(skoverbuffer(imap%ninteg, nBatchline))
      write(*, "(A,I0,A,F6.3,A,F6.3)") "Calculating ", nBatchline, " lines: r0 = ",&
          & inp%r0 + inp%dr * real(nBatch * nBatchline, dp), " dr = ", inp%dr
      do ir = 1, nBatchline
        dist = inp%r0 + inp%dr * real(nBatch * nBatchline + ir - 1, dp)
        call gengrid2_12(quads, coordtrans_becke_12, partition_becke_homo, beckepars, dist, grid1,&
            & grid2, dots, weights)
        call getskintegrals(atom1, atom2, grid1, grid2, dots, weights, inp%tDensitySuperpos,&
            & inp%ixc, imap, skhambuffer(:, ir), skoverbuffer(:, ir), denserr(ir))
      end do
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
    end do
    if (.not. tConverged) then
      write(*, "(A,F6.2,A,ES10.3)") "Warning, maximal distance ", inp%maxdist,&
          & " reached! Max integral value:", maxabs
    end if
    write(*, "(A,ES10.3)") "Maximal integration error: ", denserrmax

    ! hand over Hamiltonian and overlap
    call hamfifo%popall_concat(skham)
    call overfifo%popall_concat(skover)

  end subroutine get_twocenter_integrals


  !> Calculates SK-integrals.
  subroutine getskintegrals(atom1, atom2, grid1, grid2, dots, weights, tDensitySuperpos, ixc, imap,&
      & skham, skover, denserr)

    !> atomic property instances of dimer atoms
    type(TAtomdata), intent(in) :: atom1, atom2

    !> integration grids of dimer atoms, holding spherical coordinates (r, theta)
    real(dp), intent(in), target :: grid1(:,:), grid2(:,:)

    !> ???
    real(dp), intent(in) :: dots(:)

    !> integration weights
    real(dp), intent(in) :: weights(:)

    !> true, if density superposition is requested, otherwise potential superposition is applied
    logical, intent(in) :: tDensitySuperpos

    !> xc-functional type (0: potential superposition, 1: LDA-PW91, 2: GGA-PBE)
    integer, intent(in) :: ixc

    !> two-center integration mapping instance
    type(TIntegMap), intent(in) :: imap

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

    !! atomic 1st and 2nd density derivatives of atom 1
    real(dp), allocatable :: densval1p(:), densval1pp(:)

    !! atomic 1st and 2nd density derivatives of atom 2
    real(dp), allocatable :: densval2p(:), densval2pp(:)

    !! real tesseral spherical harmonic for spherical coordinate (theta) of atom 1 and atom 2
    real(dp), allocatable :: spherval1(:), spherval2(:)

    !! higher-level density expressions
    real(dp), allocatable :: absgr(:), laplace(:), gr_grabsgr(:)

    !! temporary storage for Hamiltonian, overlap, density and pre-factors
    real(dp) :: integ1, integ2, dens, prefac

    !! number of integration points
    integer :: nGrid

    !!  orbital indices/angular momenta on the two atoms and interaction type
    integer :: i1, i2, l1, l2, mm

    !! auxiliary variable
    integer :: ii

    r1 => grid1(:, 1)
    theta1 => grid1(:, 2)
    r2 => grid2(:, 1)
    theta2 => grid2(:, 2)
    nGrid = size(r1)

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

    allocate(potval(nGrid))
    ifPotSup: if (.not. tDensitySuperpos) then
      potval(:) = atom1%pot%getValue(r1) + atom2%pot%getValue(r2)
    else
      allocate(densval(nGrid))
      densval(:) = atom1%rho%getValue(r1) + atom2%rho%getValue(r2)
      select case (ixc)
      case (1)
        ! LDA-PW91 xc-functional
        call getxcpot_ldapw91(densval, potval)
      case (2)
        ! GGA-PBE xc-functional
        allocate(densval1p(nGrid))
        allocate(densval1pp(nGrid))
        allocate(densval2p(nGrid))
        allocate(densval2pp(nGrid))
        densval1p(:) = atom1%drho%getValue(r1)
        densval1pp(:) = atom1%ddrho%getValue(r1)
        densval2p(:) = atom2%drho%getValue(r2)
        densval2pp(:) = atom2%ddrho%getValue(r2)
        allocate(absgr(nGrid))
        allocate(laplace(nGrid))
        allocate(gr_grabsgr(nGrid))
        ! calculate derivatives for combined density
        call getDerivs(densval1p, densval1pp, densval2p, densval2pp, r1, r2, dots, absgr, laplace,&
            & gr_grabsgr)
        ! get xc-potential
        call getxcpot_ggapbe(densval, absgr, laplace, gr_grabsgr, potval)
      case default
        write(*,*) "Unknown functional type!"
        stop
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
      spherval1(:) = tes1%getValue_1d(theta1)
      spherval2(:) = tes2%getValue_1d(theta2)
      integ1 = getHamiltonian(radval1(:, i1), radval2(:, i2), radval2p(:, i2), radval2pp(:, i2),&
          & r2, l2, spherval1, spherval2, potval, weights)
      integ2 = getOverlap(radval1(:, i1), radval2(:, i2), spherval1, spherval2, weights)
      dens = getDensity(radval1(:, i1), radval2(:, i2), spherval1, spherval2, weights)
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


  !> Calculates higher-level expressions based on the density's 1st and 2nd derivatives.
  pure subroutine getDerivs(drho1, d2rho1, drho2, d2rho2, r1, r2, dots, absgr, laplace, gr_grabsgr)

    !> 1st and 2nd atomic density derivatives on grid
    real(dp), intent(in) :: drho1(:), d2rho1(:), drho2(:), d2rho2(:)

    !> radial spherical coordinates of atom 1 and atom 2 on grid
    real(dp), intent(in) :: r1(:), r2(:)

    !> ???
    real(dp), intent(in) :: dots(:)

    !> absolute total density gradient
    real(dp), intent(out) :: absgr(:)

    !> laplace operator acting on total density
    real(dp), intent(out) :: laplace(:)

    !> (grad rho4pi) * grad(abs(grad rho4pi))
    real(dp), intent(out) :: gr_grabsgr(:)

    !! temporary storage
    real(dp), allocatable :: f1(:), f2(:)

    !! number of grid points
    integer :: nn

    nn = size(drho1)
    allocate(f1(nn), f2(nn))

    f1(:) = drho1 + dots * drho2
    f2(:) = drho2 + dots * drho1

    absgr(:) = sqrt(drho1 * f1 + drho2 * f2)
    laplace(:) = d2rho1 + d2rho2 + 2.0_dp * (drho1 / r1 + drho2 / r2)
    where (absgr > epsilon(1.0_dp))
      gr_grabsgr = (d2rho1 * f1 * f1 + d2rho2 * f2 * f2&
          & + (1.0_dp - dots**2) * drho1 * drho2 * (drho2 / r1 + drho1 / r2))&
          & / absgr
    elsewhere
      gr_grabsgr = 0.0_dp
    end where

  end subroutine getDerivs


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
