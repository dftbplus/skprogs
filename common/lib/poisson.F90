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

  use common_sphericalharmonics, only : TRealTessY, TRealTessY_init
  use common_quadratures, only : TQuadrature, TQuadrature2D, gauss_chebyshev_quadrature
  use common_quadratures, only : lebedev_laikov_quadrature
  use common_gridgenerator, only : gengrid1_1, gengrid1_3, gengrid2_3
  use common_interpolation, only : ipl_tst, get_cubic_spline
  use common_coordtrans, only : coordtrans_radial_becke1, coordtrans_radial_becke2
  use common_partition, only : partition_becke_homo
  use common_anglib, only : initGaunt, realGaunt
  use common_finitedifferences, only : makeFDMatrix7P, H_BFDM7P, P_BFDM7P

  implicit none
  ! private

  public :: solve_poisson, solve_helmholz
  public :: integrator_solve_poisson, integrator_solve_helmholz
  public :: integrator_init, becke_integrator, becke_grid_params, integrator_build_LU
  public :: integrator_get_coords, integrator_set_kernel_param, integrator_precomp_fdmat


  !>
  type separable_integrand

    real(dp), allocatable :: radial(:)
    integer :: ll1
    integer :: mm1
    integer :: ll2
    integer :: mm2

  end type separable_integrand


  !>
  type fdiff_matrix

    real(dp), allocatable :: d(:)
    real(dp), allocatable :: g(:)
    real(dp), allocatable :: zi(:)
    real(dp), allocatable :: b(:)
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: H1(:,:)
    real(dp), allocatable :: H2(:,:)

    ! LU decomposition supermatrix
    real(dp), allocatable :: H3(:,:,:)
    integer, allocatable :: ipiv2(:,:)
    real(dp), allocatable :: H4(:,:,:)
    integer, allocatable :: ipiv4(:,:)

  end type fdiff_matrix


  !>
  type becke_grid_params

    integer :: N_radial
    integer :: N_angular
    integer :: ll_max
    real(dp) :: rm

  end type becke_grid_params


  !>
  type becke_subgrid

    !>
    real(dp), allocatable :: data(:,:)

  end type becke_subgrid


  !>
  type becke_grid

    integer :: type

    !>
    type(becke_subgrid), allocatable :: subgrid(:)
    real(dp), allocatable :: weight(:), part(:)

  end type becke_grid


  !> Contains information, needed for setting up and performing the integration.
  type becke_integrator

    type(becke_grid_params) :: grid_params
    type(TQuadrature) :: radial_quadrature
    type(TQuadrature2D) :: angular_quadrature

    !>
    integer :: kernel_type ! Coulomb(0) or Yukawa(1)
    real(dp) :: kernel_parameter ! irrelevant for kernel_type=0

    !>
    type(becke_grid), allocatable :: integration_grid(:)

    !>
    type(fdiff_matrix) :: fdmat

    !>
    real(dp) :: rm

    !> amount of messages (0 for no messages)
    integer :: verbosity

  end type becke_integrator


contains

  !> precomputes the finite differences matrix.
  subroutine integrator_precomp_fdmat(self)
    type(becke_integrator), intent(inout) :: self

    integer :: N
    real(dp) :: kappa,rm

    integer :: ii
    real(dp) :: step, step_2, tmp1, llp1_pi_2_rm_4, f0, zz,pi_rm_4_3
    real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf,a2c

    kappa = self%kernel_parameter

    if (self%verbosity > 0) then
       write(*,'(a,F12.7,a)',advance="no") "Precompute the FD_Matrix for kappa="&
           &,kappa,"..."
    end if

    N=self%grid_params%N_radial
    rm=self%rm

    self%fdmat%g = 0.0_dp
    self%fdmat%H1 = 0.0_dp
    self%fdmat%H2 = 0.0_dp
    self%fdmat%ipiv = 0.0_dp
    self%fdmat%zi = 0.0_dp

    tmp1 = 1.0_dp/real(N+1, dp)
    tmp1 = tmp1*tmp1*180.0_dp

    self%fdmat%zi = dacos(self%radial_quadrature%xx)/pi

    step = 1.0_dp/real(N+1, dp)
    step_2 = step*step

    f0=0.0_dp
    a2c = kappa*kappa*rm*rm*pi*pi*0.25_dp
    llp1_pi_2_rm_4 = 4.0_dp*pi*pi!*rm
    pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp*tmp1

    zz = self%fdmat%zi(1)
    sin_pi = dsin(pi * zz)
    sin_pi_hlf = dsin(pi_hlf * zz)
    cos_pi_hlf = dcos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
    beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
    sin_pi = sin_pi*sin_pi ! ^2

    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
    gama = -(a2c*sin_pi/sin_pi_hlf)
    self%fdmat%d(1) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
    self%fdmat%g(1) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

    self%fdmat%H1(7, 1) = (-300.0_dp)     + (beta*(-150.0_dp)) + gama*step_2*180.0_dp
    self%fdmat%H1(6, 2) = (90.0_dp)       + beta*270.0_dp
    self%fdmat%H1(5, 3) = (60.0_dp)       + beta*(-90.0_dp)
    self%fdmat%H1(4, 4) = (-15.0_dp)      + beta*15.0_dp

    zz = self%fdmat%zi(2)

    sin_pi = dsin(pi * zz)
    sin_pi_hlf = dsin(pi_hlf * zz)
    cos_pi_hlf = dcos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
    beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
    sin_pi = sin_pi*sin_pi ! ^2

    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
    gama = -(a2c*sin_pi/sin_pi_hlf)
    self%fdmat%d(2) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
    self%fdmat%g(2) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

    self%fdmat%H1(7, 2) = (-450.0_dp)     + beta*(-60.0_dp)
    self%fdmat%H1(6, 3) = ( 240.0_dp)     + beta*(180.0_dp)
    self%fdmat%H1(5, 4) = ( -15.0_dp)     + beta*(-45.0_dp)
    self%fdmat%H1(4, 5) =                 + beta*(6.0_dp)
    self%fdmat%H1(8, 1) = (240.0_dp)    + (beta*(-90.0_dp))  + gama*step_2*180.0_dp

    zz = self%fdmat%zi(3)

    sin_pi = dsin(pi * zz)
    sin_pi_hlf = dsin(pi_hlf * zz)
    cos_pi_hlf = dcos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
    beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
    sin_pi = sin_pi*sin_pi ! ^2

   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(3) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(3) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7, 3) = (-490.0_dp)   + gama*step_2*180.0_dp
   self%fdmat%H1(6, 4) = (270.0_dp)    + beta*(135.0_dp)
   self%fdmat%H1(5, 5) = (-27.0_dp)    + beta*(-27.0_dp)
   self%fdmat%H1(4, 6) = (2.0_dp)      + beta*3.0_dp
   self%fdmat%H1(8, 2) = (270.0_dp)    + beta*(-135.0_dp)
   self%fdmat%H1(9, 1) = (-27.0_dp)    + beta*(27.0_dp)

   do ii=4, N-3
      zz = self%fdmat%zi(ii)

      sin_pi = dsin(pi * zz)
      sin_pi_hlf = dsin(pi_hlf * zz)
      cos_pi_hlf = dcos(pi_hlf * zz)
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
      beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
      sin_pi = sin_pi*sin_pi ! ^2

      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
      gama = -(a2c*sin_pi/sin_pi_hlf)
      self%fdmat%d(ii) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
      self%fdmat%g(ii) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

      self%fdmat%H1(7,ii) = (-490.0_dp)    + gama*step_2*180.0_dp !H(ii,ii)
      self%fdmat%H1(6, ii+1) = (270.0_dp)     + beta*(135.0_dp)!H(ii,ii+1)
      self%fdmat%H1(5, ii+2) = (-27.0_dp)     + beta*(-27.0_dp)!H(ii,ii+2)
      self%fdmat%H1(4, ii+3) = (2.0_dp)       + beta*(3.0_dp)!H(ii,ii+3)
      self%fdmat%H1(8, ii-1) = (270.0_dp)     + beta*(-135.0_dp)!H(ii,ii-1)
      self%fdmat%H1(9, ii-2) = (-27.0_dp)     + beta*(27.0_dp)!H(ii,ii-2)
      self%fdmat%H1(10, ii-3) =  (2.0_dp)     + beta*(-3.0_dp)!H(ii,ii-3)
   end do

   zz = self%fdmat%zi(N-2)
   sin_pi = dsin(pi * zz)
   sin_pi_hlf = dsin(pi_hlf * zz)
   cos_pi_hlf = dcos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2

   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(N-2) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(N-2) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7, N-2) = (-490.0_dp)  + gama*step_2*180.0_dp ! H(N-2,N-2)
   self%fdmat%H1(6, N-1) = (270.0_dp)   + beta*(135.0_dp)! H(N-2,N-1)
   self%fdmat%H1(5, N) =  (-27.0_dp)    + beta*(-27.0_dp)! H(N-2, N)
   self%fdmat%H1(8, N-3) =(270.0_dp)    + beta*(-135.0_dp)! H(N-2,N-3)
   self%fdmat%H1(9, N-4) =(-27.0_dp)    + beta*(27.0_dp)! H(N-2,N-4)
   self%fdmat%H1(10, N-5) =(2.0_dp)     + beta*(-3.0_dp)! H(N-2,N-5)

   zz = self%fdmat%zi(N-1)
   sin_pi = dsin(pi * zz)
   sin_pi_hlf = dsin(pi_hlf * zz)
   cos_pi_hlf = dcos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2

   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(N-1) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(N-1) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7,N-1) = (-450.0_dp)   + beta*(60.0_dp)!H(N-1,N-1)
   self%fdmat%H1(6, N) =  (240.0_dp)    + gama*step_2*180.0_dp + beta*(90.0_dp)!H(N-1,N)
   self%fdmat%H1(8, N-2) =(240.0_dp)    + beta*(-180.0_dp)! H(N-1,N-2)
   self%fdmat%H1(9, N-3) =(-15.0_dp)    + beta*(45.0_dp)! H(N-1,N-3)
   self%fdmat%H1(10, N-4) =             + beta*(-6.0_dp)! H(N-1,N-4)

   zz = self%fdmat%zi(N)
   sin_pi = dsin(pi * zz)
   sin_pi_hlf = dsin(pi_hlf * zz)
   cos_pi_hlf = dcos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2

   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8

   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(N) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(N) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7,N) =   (-300.0_dp)    + gama*step_2*180.0_dp + beta*150.0_dp!H(N,N)
   self%fdmat%H1(8, N-1) =(90.0_dp)      + beta*(-270.0_dp)!H(N,N-1)
   self%fdmat%H1(9, N-2) =(60.0_dp)      + beta*(90.0_dp)!H(N,N-2)
   self%fdmat%H1(10, N-3) = (-15.0_dp)   + beta*(-15.0_dp)!H(N,N-3)

   if (self%verbosity > 0) then
     write(*,'(a)') "OK"
   end if

  end subroutine integrator_precomp_fdmat


  !>
  subroutine integrator_build_fdmat(self,ll)
    type(becke_integrator), intent(inout) :: self
    integer, intent(in) :: ll

    integer :: ii,N
    real(dp) :: lll

    self%fdmat%H2=self%fdmat%H1
    lll=real(ll*(ll+1),dp)
    N = self%grid_params%N_radial
    self%fdmat%H2(7,1) = self%fdmat%H2(7,1) + self%fdmat%g(1)*lll
    self%fdmat%H2(8,1) = self%fdmat%H2(8,1) + self%fdmat%g(2)*lll
    self%fdmat%H2(7,3) = self%fdmat%H2(7,3) +  self%fdmat%g(3)*lll
    do ii=4, N-3
       self%fdmat%H2(7,ii) = self%fdmat%H2(7,ii) + self%fdmat%g(ii)*lll
    end do
    self%fdmat%H2(7,N-2) = self%fdmat%H2(7,N-2) + self%fdmat%g(N-2)*lll
    self%fdmat%H2(6,N) = self%fdmat%H2(6,N) + self%fdmat%g(N-1)*lll
    self%fdmat%H2(7,N) = self%fdmat%H2(7,N) + self%fdmat%g(N)*lll

  end subroutine integrator_build_fdmat


  ! solves the modified Helmholz equation for
  ! angular momentum ll and density rho_lm
  ! using LU decomposed finite differences matrix fdmat
  subroutine integrator_solve_helmholz(self,ll,rho_lm)
    type(becke_integrator), intent(inout) :: self
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho_lm(:)

    integer :: N,info

    N=self%grid_params%N_radial
    rho_lm = rho_lm*self%fdmat%d
    info=0
    call DGBTRS( 'No transpose', N, 3, 3, 1, self%fdmat%H3(:,:,ll+1), 10, self&
        &%fdmat%ipiv2(:,ll+1), rho_lm, N, info)

  end subroutine integrator_solve_helmholz


  ! solves the modified Helmholz equation for
  ! angular momentum ll and density rho_lm
  ! using LU decomposed finite differences matrix fdmat
  subroutine integrator_solve_poisson(self,ll,rho_lm)
    type(becke_integrator), intent(inout) :: self
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho_lm(:)

    integer :: N,info

    N=self%grid_params%N_radial
    rho_lm = rho_lm*self%fdmat%d
    info=0
    call DGBTRS( 'No transpose', N, 3, 3, 1, self%fdmat%H4(:,:,ll+1), 10, self&
        &%fdmat%ipiv4(:,ll+1), rho_lm, N, info)
  end subroutine integrator_solve_poisson


  !> Initializes the integration module.
  subroutine integrator_init(self, grid_params)
    type(becke_integrator), intent(inout) :: self
    type(becke_grid_params), intent(in) :: grid_params

    if (self%verbosity > 0) then
       write(*, '(a)') "=================================================="
       write(*, '(a)') "Initializing the integrator module"
       write(*, '(a, 3I4)') "Becke quadrature parameters: ", &
           &grid_params%N_radial, grid_params%N_angular, grid_params%ll_max
    end if

    self%grid_params=grid_params
    self%kernel_parameter = 0.0_dp
    self%rm=grid_params%rm
    self%verbosity = 0

    if(self%verbosity > 0) then
       write(*,'(2(a,F12.7))') "rm=", self%rm, " kappa=", self%kernel_parameter
    end if

    if(self%verbosity > 0) then
       write(*, '(a)',advance="no") "generating quadrature..."
    end if
    call gauss_chebyshev_quadrature(grid_params%N_radial, self%radial_quadrature)
    call lebedev_laikov_quadrature(grid_params%N_angular, self&
        &%angular_quadrature)
    if(self%verbosity > 0) then
       write(*,'(a)') "done."
       write(*, '(a)') "generating grids..."
    end if
    allocate(self%integration_grid(3)) ! allocate two grids
    call becke_grid_init(self%integration_grid(1), 1, &
        &self%radial_quadrature, self%angular_quadrature, 0.5_dp,self%rm)
    call becke_grid_init(self%integration_grid(2), 2, &
        &self%radial_quadrature, self%angular_quadrature, 0.5_dp,self%rm)
    call becke_grid_init(self%integration_grid(3), 11, &
        &self%radial_quadrature, self%angular_quadrature, 0.5_dp,self%rm)
    if(self%verbosity > 0) then
       write(*,'(a)') "done."
    end if

    call initGaunt(grid_params%ll_max)

    ! allocate the fd_mat
    allocate(self%fdmat%g(grid_params%N_radial))
    allocate(self%fdmat%d(grid_params%N_radial))
    allocate(self%fdmat%b(grid_params%N_radial))
    allocate(self%fdmat%ipiv(grid_params%N_radial))
    allocate(self%fdmat%ipiv2(grid_params%N_radial,grid_params%ll_max))
    allocate(self%fdmat%ipiv4(grid_params%N_radial,grid_params%ll_max))
    allocate(self%fdmat%zi(grid_params%N_radial))
    allocate(self%fdmat%H1(10,grid_params%N_radial))
    allocate(self%fdmat%H2(10,grid_params%N_radial))
    allocate(self%fdmat%H3(10,grid_params%N_radial,grid_params%ll_max))
    allocate(self%fdmat%H4(10,grid_params%N_radial,grid_params%ll_max))

    if (self%verbosity > 0) then
      write(*,'(a)') "integrator init done."
    end if

  end subroutine integrator_init


  !> Precomputes the LU decompositions.
  subroutine integrator_build_LU(self)
    type(becke_integrator), intent(inout) :: self

    integer :: N,ll,info

    N = self%grid_params%N_radial

    do ll=1, self%grid_params%ll_max
       call integrator_build_fdmat(self,ll-1)
       self%fdmat%H3(:,:,ll) = self%fdmat%H2(:,:)
       call DGBTRF( N, N, 3, 3, self%fdmat%H3(:,:,ll), 10, self%fdmat%ipiv2(:,ll), info )
       if (info /= 0) write(*,'(a)') "ERROR: in LU decomposition!!!"
    end do

  end subroutine integrator_build_LU


  !> grid_no(grid_nr, subgrid_nr, coordinate_nr)
  !! grid_nr: 1 for 1_3 grid
  !!          2 for 2_3 grid
  !! subgrid_nr: 1 for 1_3
  !!             2 for 2_3
  !! coordiante_nr: 1 for rr
  !!                2 for theta
  !!                3 for phi
  subroutine integrator_get_coords(self, grid_no, coords)
    type(becke_integrator), target, intent(in) :: self
    integer, intent(in) :: grid_no(3)
    real(dp), pointer, intent(out) :: coords(:)
    ! build in the error checking
    coords => self%integration_grid(grid_no(1))%subgrid(grid_no(2))%data(:,grid_no(3))
  end subroutine integrator_get_coords


  !>
  subroutine integrator_set_kernel_param(self, kappa)
    type(becke_integrator), intent(inout) :: self
    real(dp), intent(in) :: kappa
    self%kernel_parameter = kappa
  end subroutine integrator_set_kernel_param


  !>
  subroutine becke_grid_init(self, n_subgrids, radquad, angquad, dist, rm)
    type(becke_grid), intent(inout) :: self
    integer, intent(in) :: n_subgrids
    type(TQuadrature), intent(in) :: radquad
    type(TQuadrature2D), intent(in) :: angquad
    real(dp), intent(in) :: dist
    real(dp), intent(in) :: rm

    real(dp) :: beckepars(1)

    if(n_subgrids==11) then
       self%type = 11
       allocate( self%subgrid(1)) ! becke_grid contains only one subgrid
       call gengrid1_1(radQuad, rm, coordtrans_radial_becke2, self%subgrid(1)%data, self%weight)
    end if

    if(n_subgrids==1) then
       self%type = 1
       allocate( self%subgrid(1)) ! becke_grid contains only one subgrid
       call gengrid1_3(angquad, radquad,&
           & coordtrans_radial_becke1, self%subgrid(1)%data, self%weight)
    end if
    if(n_subgrids==2) then
       self%type = 2
       allocate( self%subgrid(2)) ! becke_grid contains two subgrids
       call gengrid2_3(angquad, radquad, &
           &  coordtrans_radial_becke1, partition_becke_homo,beckepars,&
           & dist, self%subgrid(1)%data, self%subgrid(2)%data, &
           & self%weight, self%part)
     end if

  end subroutine becke_grid_init


  !>
  subroutine becke_grid_kill(self)
    type(becke_grid), intent(inout) :: self

    integer :: ii

    do ii=1, size(self%subgrid)
       deallocate(self%subgrid(ii)%data)
    end do
    deallocate(self%subgrid)
    deallocate(self%weight)

  end subroutine becke_grid_kill


  !> Solves the poisson equation.
  !! \param ll: angular momentum
  !! \param rho: lm-component of density
  !! \param u: solution
  !! \param charge: integral over density (boundary condition)
  subroutine solve_poisson(ll, rho, zi, charge)
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho(:)
    real(dp), allocatable, intent(in) :: zi(:)
    real(dp), intent(in) :: charge

    integer :: ii, info
    integer, allocatable :: ipiv(:)
    integer :: N
    real(dp), allocatable :: alpha(:), beta(:), gama(:), delta(:)!,
    real(dp), allocatable :: H2(:,:), B(:)
    real(dp), allocatable :: solution(:), work(:), dummy(:)
    real, allocatable :: swork(:)
    real(dp) :: f0, f1, tmp1
    real(dp), parameter :: rm = 1.0_dp
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf

    ! number of grid points
    N=size(rho)

    !***********************************
    ! boundary conditions
    !***********************************
    if (ll == 0) then
       f0 = charge*sqrt(4.0_dp*pi) ! r->oo
    else
       f0 = 0.0_dp ! r->oo
    end if
    f1 = 0.0_dp ! r->0

    !********************************
    ! matrix/vector allocation
    !********************************
    allocate(alpha(N))
    allocate(beta(N))
    allocate(gama(N))
    allocate(delta(N))
    allocate(B(N))
    allocate(ipiv(N))
    allocate(work(N))
    allocate(solution(N))
    allocate(swork(N*(N+1)))
    allocate(dummy(N))
    allocate(H2(10,N))
    H2 = 0.0_dp

    pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
    llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)
    do ii=1, N

       sin_pi = sin(pi * zi(ii))
       sin_pi_hlf = sin(pi_hlf * zi(ii))
       cos_pi_hlf = cos(pi_hlf * zi(ii))

       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
       cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
       cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4

       alpha(ii) = 1.0_dp
       beta(ii) = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)

       sin_pi = sin_pi*sin_pi ! ^2
       gama(ii) = -llp1_pi_2_rm_4/sin_pi

       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
       delta(ii) = -pi_rm_4_3*rho(ii)*cos_pi_hlf/sin_pi_hlf

    end do

    B = 0.0_dp
    call P_BFDM7P(H2, B, N, ll, zi, rm, charge)

    tmp1 = 1.0_dp/real(N+1, dp)
    tmp1 = tmp1*tmp1*180.0_dp
    B = B + tmp1*delta

    !***********************************************************
    ! solve the band diagonal matrix
    !***********************************************************
    call DGBSV(N, 3, 3, 1, H2, 10, ipiv, B, N, info)

    do ii=1, N
       rho(ii) =  B(ii) ! delta(ii)
    end do

  end subroutine solve_poisson


  !>
  subroutine solve_helmholz(ll, rho, zi, kappa)
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho(:)
    real(dp), allocatable, intent(in) :: zi(:)
    real(dp), intent(in) :: kappa

    integer :: ii, info, iter
    integer, allocatable :: ipiv(:)
    integer :: N
    real(dp), allocatable :: H(:,:), alpha(:), beta(:), gama(:), delta(:)!, zi(:), ri(:)
    real(dp), allocatable :: solution(:), work(:), dummy(:)
    real, allocatable :: swork(:)
    real(dp) :: f0, f1
    real(dp), parameter :: rm = 1.0_dp
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf, a2c

    ! number of grid points
    N=size(rho)

    !********************************
    ! boundary conditions
    !********************************
    f1 = 0.0_dp
    f0 = 0.0_dp

    !********************************
    ! matrix/vector allocation
    !********************************
    allocate(H(N,N))
    allocate(alpha(N))
    allocate(beta(N))
    allocate(gama(N))
    allocate(delta(N))
    allocate(ipiv(N))
    allocate(work(N))
    allocate(solution(N))
    allocate(swork(N*(N+1)))
    allocate(dummy(N))
    H = 0.0_dp

    !****************************************
    ! set up the equation:
    !****************************************
    pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
    llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)
    a2c = kappa*kappa*rm*rm*pi*pi*0.25_dp
    do ii=1, N
       sin_pi = sin(pi * zi(ii))
       sin_pi_hlf = sin(pi_hlf * zi(ii))
       cos_pi_hlf = cos(pi_hlf * zi(ii))
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
       cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
       cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
       alpha(ii) = 1.0_dp
       beta(ii) = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)
       sin_pi = sin_pi*sin_pi ! ^2
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
       gama(ii) = -(llp1_pi_2_rm_4/sin_pi + a2c*sin_pi/sin_pi_hlf)
       delta(ii) = -pi_rm_4_3*rho(ii)*cos_pi_hlf/sin_pi_hlf
    end do

    !******************************
    ! generate the FD Matrix
    !******************************
    ! 7-point FD scheme, ref: Bickley
    call makeFDMatrix7P(H, delta, N, alpha, beta, gama, dummy, f0, f1)

    !******************************
    ! call LAPACK (d)sgesv routine
    ! to solve the linear eqn Hx=b
    !******************************
    call DSGESV(N, 1, H, N, ipiv, delta, N, solution, N, work, swork, iter, info)

    do ii = 1, N
      rho(ii) = solution(ii)
    end do

  end subroutine solve_helmholz

end module common_poisson
