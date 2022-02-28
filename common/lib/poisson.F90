!> Module that solves the Poisson equation.
!!    ref: Becke, J. Chem. Phys., 89, 2993(1988)
!!    implementation: V.L. 2011
!!    Related Project: Range-Separated Hybrids (RSH)
!!
!!    Purpose: The program should solve the poisson equation with density, given
!!             at non-equidistant points (e.g. Chebyshev radial quadrature points),
!!             by transforming the equation to the equidistant mesh and solving it
!!             by finite difference method.
module common_poisson

  use common_accuracy, only : dp
  use common_constants, only : pi, pi_hlf

  use common_quadratures, only : TQuadrature, TQuadrature2D, gauss_chebyshev_quadrature
  use common_quadratures, only : lebedev_laikov_quadrature
  use common_gridgenerator, only : gengrid1_1, gengrid1_3, gengrid2_3
  use common_coordtrans, only : coordtrans_radial_becke1, coordtrans_radial_becke2
  use common_partition, only : partition_becke_homo
  use common_anglib, only : initGaunt
  use common_finitedifferences, only : makeFDMatrix7P, P_BFDM7P

  implicit none
  private

  public :: solvePoisson, solveHelmholz
  public :: TBeckeIntegrator_solvePoisson, TBeckeIntegrator_solveHelmholz
  public :: TBeckeIntegrator_init, TBeckeIntegrator, TBeckeGridParams, TBeckeIntegrator_buildLU
  public :: TBeckeIntegrator_getCoords, TBeckeIntegrator_setKernelParam
  public :: TBeckeIntegrator_precompFdMatrix


  !> Stores finite differences matrix.
  type TFdiffMatrix

    real(dp), allocatable :: d(:)
    real(dp), allocatable :: g(:)
    real(dp), allocatable :: zi(:)
    real(dp), allocatable :: b(:)
    real(dp), allocatable :: H1(:,:), H2(:,:)

    ! LU decomposition supermatrix
    real(dp), allocatable :: H3(:,:,:)
    real(dp), allocatable :: H4(:,:,:)

  end type TFdiffMatrix


  !> Defines a Becke integration grid.
  type TBeckeGridParams

    !> number of radial integration points
    integer :: nRadial

    !> number of angular integration points
    integer :: nAngular

    !> maximum angular momentum
    integer :: ll_max

    !> midpoint of the integration interval
    real(dp) :: rm

  end type TBeckeGridParams


  !> Wraps around Becke data array.
  type TBeckeSubgrid

    !> data on Becke grid
    real(dp), allocatable :: data(:,:)

  end type TBeckeSubgrid


  !> Wraps around multiple Becke integration grids.
  type TBeckeGrid

    !> Becke subgrids
    type(TBeckeSubgrid), allocatable :: subgrid(:)

    !> integration weights
    real(dp), allocatable :: weight(:)

    !> space partition for full 3D integration
    real(dp), allocatable :: partition(:)

  end type TBeckeGrid


  !> Contains information, needed for setting up and performing Becke integration.
  type TBeckeIntegrator

    !> Becke grid parameters
    type(TBeckeGridParams) :: beckeGridParams

    !> radial quadrature information
    type(TQuadrature) :: radialQuadrature

    !> angular quadrature information
    type(TQuadrature2D) :: angularQuadrature

    !> kernel screening parameter
    real(dp) :: kernelParameter

    !> Becke integration grids
    type(TBeckeGrid), allocatable :: beckeGrid(:)

    !> finite differences matrix
    type(TFdiffMatrix) :: fdmat

    !> midpoint of the integration interval
    real(dp) :: rm

  end type TBeckeIntegrator


contains

  !> Precomputes the finite differences matrix.
  subroutine TBeckeIntegrator_precompFdMatrix(this)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !! number of radial grid points
    integer :: nRadial

    !! screening parameter
    real(dp) :: kappa

    !! midpoint of the integration interval
    real(dp) :: rm

    !! auxiliary variables
    real(dp) :: step, step_2, tmp1, llp1_pi_2_rm_4, f0, zz,pi_rm_4_3
    real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf, a2c

    !! iterates over radial points
    integer :: ii

    kappa = this%kernelParameter
    nRadial=this%beckeGridParams%nRadial
    rm = this%rm

    this%fdmat%g(:) = 0.0_dp
    this%fdmat%H1(:,:) = 0.0_dp
    this%fdmat%H2(:,:) = 0.0_dp
    this%fdmat%zi(:) = 0.0_dp

    tmp1 = (1.0_dp / real(nRadial + 1, dp))**2 * 180.0_dp

    this%fdmat%zi(:) = acos(this%radialQuadrature%xx) / pi

    step = 1.0_dp / real(nRadial + 1, dp)
    step_2 = step**2

    f0 = 0.0_dp
    a2c = kappa**2 * rm**2 * pi**2 * 0.25_dp
    llp1_pi_2_rm_4 = 4.0_dp * pi**2
    pi_rm_4_3 = pi**3 * rm**3 * 4.0_dp * tmp1

    zz = this%fdmat%zi(1)
    sin_pi = sin(pi * zz)
    sin_pi_hlf = sin(pi_hlf * zz)
    cos_pi_hlf = cos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2

    sin_pi_hlf = sin_pi_hlf**2 ! ^4
    sin_pi_hlf = sin_pi_hlf**2 ! ^8
    gama = -(a2c * sin_pi / sin_pi_hlf)
    this%fdmat%d(1) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
    this%fdmat%g(1) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

    this%fdmat%H1(7, 1) = -300.0_dp - beta * 150.0_dp + gama * step_2 * 180.0_dp
    this%fdmat%H1(6, 2) = 90.0_dp + beta * 270.0_dp
    this%fdmat%H1(5, 3) = 60.0_dp - beta * 90.0_dp
    this%fdmat%H1(4, 4) = -15.0_dp + beta * 15.0_dp

    zz = this%fdmat%zi(2)

    sin_pi = sin(pi * zz)
    sin_pi_hlf = sin(pi_hlf * zz)
    cos_pi_hlf = cos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2

    sin_pi_hlf = sin_pi_hlf**2 ! ^4
    sin_pi_hlf = sin_pi_hlf**2 ! ^8
    gama = -(a2c * sin_pi / sin_pi_hlf)
    this%fdmat%d(2) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
    this%fdmat%g(2) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

    this%fdmat%H1(7, 2) = -450.0_dp - beta * 60.0_dp
    this%fdmat%H1(6, 3) = 240.0_dp + beta * 180.0_dp
    this%fdmat%H1(5, 4) = -15.0_dp - beta * 45.0_dp
    this%fdmat%H1(4, 5) = beta * 6.0_dp
    this%fdmat%H1(8, 1) = 240.0_dp - beta * 90.0_dp + gama * step_2 * 180.0_dp

    zz = this%fdmat%zi(3)

    sin_pi = sin(pi * zz)
    sin_pi_hlf = sin(pi_hlf * zz)
    cos_pi_hlf = cos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2

    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
    gama = -(a2c * sin_pi / sin_pi_hlf)
    this%fdmat%d(3) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
    this%fdmat%g(3) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

    this%fdmat%H1(7, 3) = -490.0_dp + gama * step_2 * 180.0_dp
    this%fdmat%H1(6, 4) = 270.0_dp + beta * 135.0_dp
    this%fdmat%H1(5, 5) = -27.0_dp - beta * 27.0_dp
    this%fdmat%H1(4, 6) = 2.0_dp + beta * 3.0_dp
    this%fdmat%H1(8, 2) = 270.0_dp - beta * 135.0_dp
    this%fdmat%H1(9, 1) = -27.0_dp + beta * 27.0_dp

   do ii = 4, nRadial - 3
     zz = this%fdmat%zi(ii)

     sin_pi = sin(pi * zz)
     sin_pi_hlf = sin(pi_hlf * zz)
     cos_pi_hlf = cos(pi_hlf * zz)
     sin_pi_hlf = sin_pi_hlf**2 ! ^2
     cos_pi_hlf = cos_pi_hlf**4 ! ^4
     beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
     sin_pi = sin_pi**2 ! ^2

     sin_pi_hlf = sin_pi_hlf**2 ! ^4
     sin_pi_hlf = sin_pi_hlf**2 ! ^8
     gama = -(a2c * sin_pi / sin_pi_hlf)
     this%fdmat%d(ii) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
     this%fdmat%g(ii) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

     this%fdmat%H1(7, ii) = -490.0_dp + gama * step_2 * 180.0_dp
     this%fdmat%H1(6, ii + 1) = 270.0_dp + beta * 135.0_dp
     this%fdmat%H1(5, ii + 2) = -27.0_dp - beta * 27.0_dp
     this%fdmat%H1(4, ii + 3) = 2.0_dp + beta * 3.0_dp
     this%fdmat%H1(8, ii - 1) = 270.0_dp - beta * 135.0_dp
     this%fdmat%H1(9, ii - 2) = -27.0_dp + beta * 27.0_dp
     this%fdmat%H1(10, ii - 3) = 2.0_dp - beta * 3.0_dp
   end do

   zz = this%fdmat%zi(nRadial - 2)

   sin_pi = sin(pi * zz)
   sin_pi_hlf = sin(pi_hlf * zz)
   cos_pi_hlf = cos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf**2 ! ^2
   cos_pi_hlf = cos_pi_hlf**4 ! ^4
   beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
   sin_pi = sin_pi**2 ! ^2

   sin_pi_hlf = sin_pi_hlf**2 ! ^4
   sin_pi_hlf = sin_pi_hlf**2 ! ^8
   gama = -(a2c * sin_pi / sin_pi_hlf)
   this%fdmat%d(nRadial - 2) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
   this%fdmat%g(nRadial - 2) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

   this%fdmat%H1(7, nRadial - 2) = -490.0_dp + gama * step_2 * 180.0_dp
   this%fdmat%H1(6, nRadial - 1) = 270.0_dp + beta * 135.0_dp
   this%fdmat%H1(5, nRadial) = -27.0_dp - beta * 27.0_dp
   this%fdmat%H1(8, nRadial - 3) = 270.0_dp - beta * 135.0_dp
   this%fdmat%H1(9, nRadial - 4) = -27.0_dp + beta * 27.0_dp
   this%fdmat%H1(10, nRadial - 5) = 2.0_dp - beta * 3.0_dp

   zz = this%fdmat%zi(nRadial - 1)

   sin_pi = sin(pi * zz)
   sin_pi_hlf = sin(pi_hlf * zz)
   cos_pi_hlf = cos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf**2 ! ^2
   cos_pi_hlf = cos_pi_hlf**4 ! ^4
   beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
   sin_pi = sin_pi**2 ! ^2

   sin_pi_hlf = sin_pi_hlf**2 ! ^4
   sin_pi_hlf = sin_pi_hlf**2 ! ^8
   gama = -(a2c * sin_pi / sin_pi_hlf)
   this%fdmat%d(nRadial - 1) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
   this%fdmat%g(nRadial - 1) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

   this%fdmat%H1(7, nRadial - 1) = -450.0_dp + beta * 60.0_dp
   this%fdmat%H1(6, nRadial) = 240.0_dp + gama * step_2 * 180.0_dp + beta * 90.0_dp
   this%fdmat%H1(8, nRadial - 2) = 240.0_dp - beta * 180.0_dp
   this%fdmat%H1(9, nRadial - 3) = -15.0_dp + beta * 45.0_dp
   this%fdmat%H1(10, nRadial - 4) = -beta * 6.0_dp

   zz = this%fdmat%zi(nRadial)

   sin_pi = sin(pi * zz)
   sin_pi_hlf = sin(pi_hlf * zz)
   cos_pi_hlf = cos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf**2 ! ^2
   cos_pi_hlf = cos_pi_hlf**4 ! ^4
   beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
   sin_pi = sin_pi**2 ! ^2

   sin_pi_hlf = sin_pi_hlf**2 ! ^4
   sin_pi_hlf = sin_pi_hlf**2 ! ^8

   gama = -(a2c * sin_pi / sin_pi_hlf)
   this%fdmat%d(nRadial) = -pi_rm_4_3 * cos_pi_hlf / sin_pi_hlf
   this%fdmat%g(nRadial) = -llp1_pi_2_rm_4 / sin_pi * step_2 * 180.0_dp

   this%fdmat%H1(7, nRadial) = -300.0_dp + gama * step_2 * 180.0_dp + beta * 150.0_dp
   this%fdmat%H1(8, nRadial - 1) = 90.0_dp - beta * 270.0_dp
   this%fdmat%H1(9, nRadial - 2) = 60.0_dp + beta * 90.0_dp
   this%fdmat%H1(10, nRadial - 3) = -15.0_dp - beta * 15.0_dp

  end subroutine TBeckeIntegrator_precompFdMatrix


  !>
  subroutine TBeckeIntegrator_buildFdMatrix(this, ll)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !> angular momentum
    integer, intent(in) :: ll

    !! number of radial points
    integer :: nRadial

    !! iterates over radial points
    integer :: ii

    !! stores: ll * (ll + 1)
    real(dp) :: lll

    this%fdmat%H2 = this%fdmat%H1
    lll = real(ll * (ll + 1), dp)
    nRadial = this%beckeGridParams%nRadial

    this%fdmat%H2(7, 1) = this%fdmat%H2(7, 1) + this%fdmat%g(1) * lll
    this%fdmat%H2(8, 1) = this%fdmat%H2(8, 1) + this%fdmat%g(2) * lll
    this%fdmat%H2(7, 3) = this%fdmat%H2(7, 3) +  this%fdmat%g(3) * lll
    do ii = 4, nRadial - 3
       this%fdmat%H2(7, ii) = this%fdmat%H2(7, ii) + this%fdmat%g(ii) * lll
    end do
    this%fdmat%H2(7, nRadial - 2) = this%fdmat%H2(7, nRadial - 2) + this%fdmat%g(nRadial - 2) * lll
    this%fdmat%H2(6, nRadial) = this%fdmat%H2(6, nRadial) + this%fdmat%g(nRadial - 1) * lll
    this%fdmat%H2(7, nRadial) = this%fdmat%H2(7, nRadial) + this%fdmat%g(nRadial) * lll

  end subroutine TBeckeIntegrator_buildFdMatrix


  !> Solves the modified Helmholz equation for angular momentum ll and density rho_lm,
  !! using LU decomposed finite differences matrix fdmat.
  subroutine TBeckeIntegrator_solveHelmholz(this, ll, rho_lm)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !> angular momentum
    integer, intent(in) :: ll

    !> lm-resolved density
    real(dp), intent(inout), allocatable :: rho_lm(:)

    !! pivot indices: for 1 <= i <= N, row i of the matrix was interchanged with row ipiv(i)
    integer, allocatable :: ipiv(:)

    !! number of radial points
    integer :: nRadial

    !! error status
    integer :: info

    info = 0

    allocate(ipiv(this%beckeGridParams%nRadial))
    ipiv(:) = 0

    nRadial = this%beckeGridParams%nRadial
    rho_lm = rho_lm * this%fdmat%d

    call DGBTRS('No transpose', nRadial, 3, 3, 1, this%fdmat%H3(:,:, ll + 1), 10, ipiv, rho_lm,&
        & nRadial, info)

  end subroutine TBeckeIntegrator_solveHelmholz


  !> Solves the Poisson equation for angular momentum ll and density rho_lm,
  !! using LU decomposed finite differences matrix fdmat.
  subroutine TBeckeIntegrator_solvePoisson(this, ll, rho_lm)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !> angular momentum
    integer, intent(in) :: ll

    !> lm-resolved density
    real(dp), intent(inout), allocatable :: rho_lm(:)

    !! pivot indices: for 1 <= i <= N, row i of the matrix was interchanged with row ipiv(i)
    integer, allocatable :: ipiv(:)

    !! number of radial points
    integer :: nRadial

    !! error status
    integer :: info

    info = 0

    allocate(ipiv(this%beckeGridParams%nRadial))
    ipiv(:) = 0

    nRadial = this%beckeGridParams%nRadial
    rho_lm = rho_lm * this%fdmat%d

    call DGBTRS('No transpose', nRadial, 3, 3, 1, this%fdmat%H4(:,:, ll + 1), 10, ipiv, rho_lm,&
        & nRadial, info)

  end subroutine TBeckeIntegrator_solvePoisson


  !> Initializes the integration module.
  subroutine TBeckeIntegrator_init(this, beckeGridParams)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !> Becke grid parameters
    type(TBeckeGridParams), intent(in) :: beckeGridParams

    this%beckeGridParams = beckeGridParams
    this%kernelParameter = 0.0_dp
    this%rm = beckeGridParams%rm

    call gauss_chebyshev_quadrature(beckeGridParams%nRadial, this%radialQuadrature)
    call lebedev_laikov_quadrature(beckeGridParams%nAngular, this%angularQuadrature)

    allocate(this%beckeGrid(3))
    call TBeckeGrid_init(this%beckeGrid(1), 1, this%radialQuadrature, this%angularQuadrature,&
        & 0.5_dp, this%rm)
    call TBeckeGrid_init(this%beckeGrid(2), 2, this%radialQuadrature, this%angularQuadrature,&
        & 0.5_dp, this%rm)
    call TBeckeGrid_init(this%beckeGrid(3), 11, this%radialQuadrature, this%angularQuadrature,&
        & 0.5_dp, this%rm)

    call initGaunt(beckeGridParams%ll_max)

    ! allocate the finite differences matrix
    allocate(this%fdmat%g(beckeGridParams%nRadial))
    allocate(this%fdmat%d(beckeGridParams%nRadial))
    allocate(this%fdmat%b(beckeGridParams%nRadial))
    allocate(this%fdmat%zi(beckeGridParams%nRadial))
    allocate(this%fdmat%H1(10, beckeGridParams%nRadial))
    allocate(this%fdmat%H2(10, beckeGridParams%nRadial))
    allocate(this%fdmat%H3(10, beckeGridParams%nRadial, beckeGridParams%ll_max))
    allocate(this%fdmat%H4(10, beckeGridParams%nRadial, beckeGridParams%ll_max))

  end subroutine TBeckeIntegrator_init


  !> Precomputes the LU decompositions.
  subroutine TBeckeIntegrator_buildLU(this)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !! pivot indices: for 1 <= i <= N, row i of the matrix was interchanged with row ipiv(i)
    integer, allocatable :: ipiv(:)

    !! number of radial points
    integer :: nRadial

    !! angular momentum
    integer :: ll

    !! error status
    integer :: info

    nRadial = this%beckeGridParams%nRadial

    allocate(ipiv(nRadial))
    ipiv(:) = 0

    do ll = 1, this%beckeGridParams%ll_max
       call TBeckeIntegrator_buildFdMatrix(this, ll - 1)
       this%fdmat%H3(:,:, ll) = this%fdmat%H2
       call DGBTRF(nRadial, nRadial, 3, 3, this%fdmat%H3(:,:, ll), 10, ipiv, info)
       if (info /= 0) write(*, '(A)') "ERROR: LU decomposition failed!"
    end do

  end subroutine TBeckeIntegrator_buildLU


  !> Returns pointer to selected grid coordinates.
  !! grid_no(grid_nr, subgrid_nr, coordinate_nr)
  !! grid_nr: 1 for 1_3 grid
  !!          2 for 2_3 grid
  !! subgrid_nr: 1 for 1_3
  !!             2 for 2_3
  !! coordiante_nr: 1 for rr
  !!                2 for theta
  !!                3 for phi
  subroutine TBeckeIntegrator_getCoords(this, grid_no, pCoords)

    !> Becke integrator instance
    type(TBeckeIntegrator), target, intent(in) :: this

    !> selection indices, i.e. [grid_nr, subgrid_nr, coordinate_nr]
    integer, intent(in) :: grid_no(3)

    !> pointer to coordinates
    real(dp), intent(out), pointer :: pCoords(:)

    pCoords => this%beckeGrid(grid_no(1))%subgrid(grid_no(2))%data(:, grid_no(3))

  end subroutine TBeckeIntegrator_getCoords


  !> Sets the screening parameter in the integration kernel.
  subroutine TBeckeIntegrator_setKernelParam(this, omega)

    !> Becke integrator instance
    type(TBeckeIntegrator), intent(inout) :: this

    !> screening parameter
    real(dp), intent(in) :: omega

    this%kernelParameter = omega

  end subroutine TBeckeIntegrator_setKernelParam


  !> Initializes a Becke grid instance with a given number of subgrids.
  subroutine TBeckeGrid_init(this, nSubgrids, radquad, angquad, dist, rm)

    !> Becke grid instance
    type(TBeckeGrid), intent(inout) :: this

    !> number of subgrids in instance
    integer, intent(in) :: nSubgrids

    !> abscissas and weight instances for radial quadrature
    type(TQuadrature), intent(in) :: radquad

    !> abscissas and weight instances for angular quadrature
    type(TQuadrature2D), intent(in) :: angquad

    !> distance between centers
    real(dp), intent(in) :: dist

    !> midpoint of the integration interval
    real(dp), intent(in) :: rm

    !! arbitrary dummy real array, unused for homonuclear Becke partitioning
    real(dp) :: beckepars(1)

    ! TBeckeGrid contains only one subgrid
    if (nSubgrids == 11) then
      allocate(this%subgrid(1))
      call gengrid1_1(radQuad, rm, coordtrans_radial_becke2, this%subgrid(1)%data, this%weight)
    ! TBeckeGrid contains only one subgrid
    elseif (nSubgrids == 1) then
      allocate(this%subgrid(1))
      call gengrid1_3(angquad, radquad, coordtrans_radial_becke1, this%subgrid(1)%data, this%weight)
    ! TBeckeGrid contains two subgrids
    elseif (nSubgrids == 2) then
      allocate(this%subgrid(2))
      call gengrid2_3(angquad, radquad, coordtrans_radial_becke1, partition_becke_homo,beckepars,&
          & dist, this%subgrid(1)%data, this%subgrid(2)%data, this%weight, this%partition)
    end if

  end subroutine TBeckeGrid_init


  !> Kills a Becke grid instance by deallocating arrays.
  subroutine TBeckeGrid_kill(this)

    !> Becke grid instance
    type(TBeckeGrid), intent(inout) :: this

    !! iterates over all subgrids
    integer :: ii

    do ii = 1, size(this%subgrid)
      deallocate(this%subgrid(ii)%data)
    end do

    deallocate(this%subgrid)
    deallocate(this%weight)

  end subroutine TBeckeGrid_kill


  !> Solves the poisson equation.
  subroutine solvePoisson(ll, rho, zi, charge)

    !> angular momentum
    integer, intent(in) :: ll

    !> lm-component of density (contains solution at exit)
    real(dp), intent(inout), allocatable :: rho(:)

    !>
    real(dp), intent(in), allocatable :: zi(:)

    !> integral over density (boundary condition)
    real(dp), intent(in) :: charge

    !! iterates over grid points
    integer :: ii

    !! error status
    integer :: info

    !! pivot indices: for 1 <= i <= N, row i of the matrix was interchanged with row ipiv(i)
    integer, allocatable :: ipiv(:)

    !! number of tabulated density grid points
    integer :: nGridPts

    !!
    real(dp), allocatable :: alpha(:), beta(:), gama(:), delta(:)
    real(dp), allocatable :: H2(:,:), B(:)
    real(dp), allocatable :: solution(:), work(:), dummy(:)
    real, allocatable :: swork(:)
    real(dp) :: f0, f1, tmp1
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf

    !! midpoint of the integration interval
    real(dp), parameter :: rm = 1.0_dp

    nGridPts = size(rho)

    !********************************
    ! matrix/vector allocation
    !********************************
    allocate(alpha(nGridPts))
    allocate(beta(nGridPts))
    allocate(gama(nGridPts))
    allocate(delta(nGridPts))
    allocate(B(nGridPts))
    allocate(ipiv(nGridPts))
    allocate(work(nGridPts))
    allocate(solution(nGridPts))
    allocate(swork(nGridPts * (nGridPts + 1)))
    allocate(dummy(nGridPts))
    allocate(H2(10, nGridPts))
    B(:) = 0.0_dp
    H2(:,:) = 0.0_dp

    !***********************************
    ! boundary conditions
    !***********************************
    if (ll == 0) then
      ! r --> oo
      f0 = charge * sqrt(4.0_dp * pi)
    else
      ! r --> oo
      f0 = 0.0_dp
    end if
    ! r --> 0
    f1 = 0.0_dp

    pi_rm_4_3 = pi**3 * rm**3 * 4.0_dp
    llp1_pi_2_rm_4 = 4.0_dp * pi**2 * rm * real(ll * (ll + 1), dp)

    do ii = 1, nGridPts

       sin_pi = sin(pi * zi(ii))
       sin_pi_hlf = sin(pi_hlf * zi(ii))
       cos_pi_hlf = cos(pi_hlf * zi(ii))

       sin_pi_hlf = sin_pi_hlf**2 ! ^2
       cos_pi_hlf = cos_pi_hlf**4 ! ^4

       alpha(ii) = 1.0_dp
       beta(ii) = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf)

       sin_pi = sin_pi**2 ! ^2
       gama(ii) = -llp1_pi_2_rm_4 / sin_pi

       sin_pi_hlf = sin_pi_hlf**4 ! ^8
       delta(ii) = -pi_rm_4_3 * rho(ii) * cos_pi_hlf / sin_pi_hlf

    end do

    call P_BFDM7P(H2, B, nGridPts, ll, zi, rm, charge)

    tmp1 = (1.0_dp / real(nGridPts + 1, dp))**2 * 180.0_dp
    B(:) = B + tmp1 * delta

    !***********************************************************
    ! solve the band diagonal matrix
    !***********************************************************
    call DGBSV(nGridPts, 3, 3, 1, H2, 10, ipiv, B, nGridPts, info)

    rho(:) = B

  end subroutine solvePoisson


  !> Solves the Helmholz equation.
  subroutine solveHelmholz(ll, rho, zi, omega)

    !> angular momentum
    integer, intent(in) :: ll

    !>
    real(dp), intent(inout) :: rho(:)

    !>
    real(dp), intent(in) :: zi(:)

    !> screening parameter
    real(dp), intent(in) :: omega

    integer :: ii, info, iter

    !! pivot indices: for 1 <= i <= N, row i of the matrix was interchanged with row ipiv(i)
    integer, allocatable :: ipiv(:)

    !! number of tabulated density grid points
    integer :: nGridPts

    !!
    real(dp), allocatable :: H(:,:), alpha(:), beta(:), gama(:), delta(:)

    !!
    real(dp), allocatable :: solution(:), work(:), dummy(:)

    !!
    real, allocatable :: swork(:)

    !! boundary conditions
    real(dp) :: f0, f1

    !! midpoint of the integration interval
    real(dp), parameter :: rm = 1.0_dp

    !! auxiliary variables
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf, a2c

    nGridPts = size(rho)

    !********************************
    ! boundary conditions
    !********************************
    f0 = 0.0_dp
    f1 = 0.0_dp

    !********************************
    ! matrix/vector allocation
    !********************************
    allocate(H(nGridPts, nGridPts))
    allocate(alpha(nGridPts))
    allocate(beta(nGridPts))
    allocate(gama(nGridPts))
    allocate(delta(nGridPts))
    allocate(ipiv(nGridPts))
    allocate(work(nGridPts))
    allocate(solution(nGridPts))
    allocate(swork(nGridPts * (nGridPts + 1)))
    allocate(dummy(nGridPts))
    H(:,:) = 0.0_dp

    !****************************************
    ! set up the equation:
    !****************************************
    pi_rm_4_3 = pi**3 * rm**3 * 4.0_dp
    llp1_pi_2_rm_4 = 4.0_dp * pi**2 * rm * real(ll * (ll + 1), dp)
    a2c = omega**2 * rm**2 * pi**2 * 0.25_dp
    do ii = 1, nGridPts
      sin_pi = sin(pi * zi(ii))
      sin_pi_hlf = sin(pi_hlf * zi(ii))
      cos_pi_hlf = cos(pi_hlf * zi(ii))
      sin_pi_hlf = sin_pi_hlf**2 ! ^2
      cos_pi_hlf = cos_pi_hlf**4 ! ^4
      alpha(ii) = 1.0_dp
      beta(ii) = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf)
      sin_pi = sin_pi**2 ! ^2
      sin_pi_hlf = sin_pi_hlf**4 ! ^8
      gama(ii) = -(llp1_pi_2_rm_4 / sin_pi + a2c * sin_pi / sin_pi_hlf)
      delta(ii) = -pi_rm_4_3 * rho(ii) * cos_pi_hlf / sin_pi_hlf
    end do

    !******************************
    ! generate the FD Matrix
    !******************************
    ! 7-point FD scheme, ref: Bickley
    call makeFDMatrix7P(H, delta, nGridPts, alpha, beta, gama, dummy, f0, f1)

    !******************************
    ! call LAPACK (d)sgesv routine
    ! to solve the linear eqn Hx=b
    !******************************
    call DSGESV(nGridPts, 1, H, nGridPts, ipiv, delta, nGridPts, solution, nGridPts, work, swork,&
        & iter, info)

    rho(:) = solution

  end subroutine solveHelmholz

end module common_poisson
