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

  public :: poisson_solve_2cnt, poisson_solve_1cnt, helmholtz_solve_1cnt, helmholtz_solve_2cnt
  public :: solve_poisson, solve_helmholz, set_tess
  public :: integrator_init, becke_integrator, becke_grid_params
  public :: integrator_get_coords, integrator_process
  public :: integrator_set_kernel_param, integrator_set_dist, integrator_precomp_fdmat
  public :: integrator_get_v1cnt

  !
  type separable_integrand
    real(dp), allocatable :: radial(:)
    integer :: ll1
    integer :: mm1
    integer :: ll2
    integer :: mm2
  end type separable_integrand

  !
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

  !
  type becke_grid_params
    integer :: N_radial
    integer :: N_angular
    integer :: ll_max
    real(dp) :: rm
  end type becke_grid_params

  !
  type becke_subgrid
    !
    real(dp), allocatable :: data(:,:)
  end type becke_subgrid

  !
  type becke_grid
    integer :: type
    !
    type(becke_subgrid), allocatable :: subgrid(:)
    real(dp), allocatable :: weight(:), part(:)
  end type becke_grid

  ! structure, containing the information
  ! needed for setup and perform the integration
  type becke_integrator
    type(becke_grid_params) :: grid_params
    type(TQuadrature) :: radial_quadrature
    type(TQuadrature2D) :: angular_quadrature
    !
    integer :: kernel_type ! Coulomb(0) or Yukawa(1)
    real(dp) :: kernel_parameter ! irrelevant for kernel_type=0
    !
    type(becke_grid), allocatable :: integration_grid(:)
    !
    type(fdiff_matrix) :: fdmat
    !
    real(dp) :: rm
    ! amount of messages (0 for no messages)
    integer :: verbosity
  end type becke_integrator

  interface integrator_init
    module procedure integrator_init
  end interface integrator_init

  interface integrator_get_coords
    module procedure integrator_get_coords
  end interface integrator_get_coords

  interface integrator_process
    module procedure integrator_process_1cnt
    module procedure integrator_process_2cnt
    module procedure integrator_process_11
  end interface integrator_process

  interface integrator_set_kernel_param
    module procedure integrator_set_kernel_param
  end interface integrator_set_kernel_param

  interface integrator_set_dist
    module procedure integrator_set_dist
  end interface integrator_set_dist


contains

  subroutine integrator_get_v1cnt(self,radial,ll1,mm1,ll2,mm2,V)
    type(becke_integrator), intent(inout) :: self
    real(dp), allocatable, intent(in) :: radial(:)
    integer, intent(in) :: ll1,mm1,ll2,mm2
    real(dp), allocatable, intent(out) :: V(:)
    !
    integer :: ngrid, N_radial
    real(dp), pointer :: rr1(:), theta1(:), phi1(:), rr3(:)
    real(dp), pointer :: rr2(:), theta2(:), phi2(:)
    type(TRealTessY) :: tess

    real(dp), allocatable, target :: zi(:)
    integer :: ii, ll, mm, kk, ind
    real(dp) :: dist, sum, tmp11, rr, yy1, yy2, yy3, tmp
    integer :: ll_max
    real(dp) :: kappa, charge

    real(dp) :: gauntcof
    real(dp), allocatable :: rho_lm(:), V_l(:,:,:), res(:)
    real(dp), allocatable :: rho2_lm(:),V2_l(:,:,:)
    real(dp), allocatable :: rrin(:), AA(:), tmp22(:), ylm(:,:,:)

    call integrator_get_coords(self,[3,1,1],rr3)

    !
    call integrator_get_coords(self, [2,2,1], rr2)
    call integrator_get_coords(self, [2,2,2], theta2)
    call integrator_get_coords(self, [2,2,3], phi2)

    !
    ll_max = self%grid_params%ll_max
    N_radial = size(rr3)
    if(N_radial .ne. size(radial)) then
       write(*,*) "ERROR: size mismatch of arrays in integrator_get_v1cnt!"
    end if
    !
    allocate(rho_lm(N_radial))
    allocate(rho2_lm(N_radial))
    allocate(zi(N_radial))
    allocate(rrin(N_radial))
    allocate(V_l(N_radial,ll_max,2))
    allocate(V2_l(N_radial,ll_max,2))
    !
    zi = acos(self%radial_quadrature%xx)/pi

    ! ==> evaluate radial part of inner integrand
    rrin = radial

    ! solve the inner integral
    do ll=0, ll_max-1
       rho_lm = rrin
 !      rho2_lm = rrin
       ! solve the equation for ll
!       call integrator_solve_helmholz(self, ll, rho_lm)
!                     call integrator_solve_poisson(self, ll, rho2_lm)
       charge = 0.0_dp
       if( ll == 0) then
          if(ll1==ll2 .and. mm1==mm2) then
             charge = sum(rho_lm*self%integration_grid(3)%weight)
          end if
!          charge = sum(rho2_lm*self%integration_grid(3)%weight)*realGaunt(ll1,mm1,ll2,mm2,ll,0)
!          charge = sum(rho2_lm*self%integration_grid(3)%weight)*realGaunt(ll1,0,ll2,0,ll,0)
!!          write (*,*) "charge:", charge, realGaunt(ll1,mm1,ll2,mm2,ll,0)
       end if
       call solve_poisson(ll, rho_lm, zi, charge)
       ! interpolate
!       rho_lm = rho2_lm - rho_lm
       call ipl_tst(zi, rho_lm, res)
       V_l(:,ll+1,1) = rho_lm
       V_l(:,ll+1,2) = res
!       call ipl_tst(zi, rho2_lm, res)
!       V2_l(:,ll+1,1) = rho2_lm
!       V2_l(:,ll+1,2) = res
    end do
    !
    ngrid = size(rr2)
!    write(*,*) ngrid
    allocate(V(ngrid))
    allocate(AA(ngrid))
    allocate(tmp22(ngrid))
    tmp22 = acos( (rr2-1.0_dp)/(rr2+1.0_dp) )/pi
    ! assemble the inner integrand
    V=0.0_dp
    do ll=0, ll_max-1
       AA=0.0_dp
       do mm=-ll,ll
          call TRealTessY_init(tess, ll, mm)
          gauntcof = realGaunt(ll1,mm1,ll2,mm2,ll,mm)
          if(abs(gauntcof) .ge. 1.0e-16_dp) then
             do ii=1, ngrid
                AA(ii) = AA(ii) + gauntcof * tess%getValue(theta2(ii), phi2(ii))
             end do
          end if
       end do
       !
       do ii=1, ngrid
          tmp = tmp22(ii)
          call get_cubic_spline(zi, V_l(:,ll+1,1), V_l(:,ll+1,2), tmp, yy2)
!          call get_cubic_spline(zi, V2_l(:,ll+1,1), V2_l(:,ll+1,2), tmp, yy3)
!          V(ii) = V(ii) + (yy3-yy2)*AA(ii)
          V(ii) = V(ii) + (yy2)*AA(ii)
       end do
    end do
  end subroutine integrator_get_v1cnt

  !! precompute the finite differences matrix
  !!
  subroutine integrator_precomp_fdmat(self)
    type(becke_integrator), intent(inout) :: self
    !
    integer :: N
    real(dp) :: kappa,rm

    integer :: ii
    real(dp) :: step, step_2, tmp1, llp1_pi_2_rm_4, f0, zz,pi_rm_4_3
    real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf,a2c

    kappa=self%kernel_parameter

    if(self%verbosity > 0) then
       write(*,'(a,F12.7,a)',advance="no") "Precompute the FD_Matrix for kappa="&
           &,kappa,"..."
    end if
    !
    N=self%grid_params%N_radial
!!=> rm !!!
    rm=self%rm!1.0_dp

    self%fdmat%g = 0.0_dp
    self%fdmat%H1 = 0.0_dp
    self%fdmat%H2 = 0.0_dp
    self%fdmat%ipiv = 0.0_dp
    self%fdmat%zi = 0.0_dp
    !
    tmp1 = 1.0_dp/real(N+1, dp)
    tmp1 = tmp1*tmp1*180.0_dp
    !
    self%fdmat%zi = dacos(self%radial_quadrature%xx)/pi
    !
    step = 1.0_dp/real(N+1, dp)
    step_2 = step*step
    !
    f0=0.0_dp
    a2c = kappa*kappa*rm*rm*pi*pi*0.25_dp
    !     llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)
    llp1_pi_2_rm_4 = 4.0_dp*pi*pi!*rm
    pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp*tmp1
    !
    zz = self%fdmat%zi(1)
    !  ii=1
    sin_pi = dsin(pi * zz)
    sin_pi_hlf = dsin(pi_hlf * zz)
    cos_pi_hlf = dcos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
    beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
    sin_pi = sin_pi*sin_pi ! ^2
    !
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
    !    gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
    gama = -(a2c*sin_pi/sin_pi_hlf)
    self%fdmat%d(1) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
    self%fdmat%g(1) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp
    !
    self%fdmat%H1(7, 1) = (-300.0_dp)     + (beta*(-150.0_dp)) + gama*step_2*180.0_dp
    self%fdmat%H1(6, 2) = (90.0_dp)       + beta*270.0_dp
    self%fdmat%H1(5, 3) = (60.0_dp)       + beta*(-90.0_dp)
    self%fdmat%H1(4, 4) = (-15.0_dp)      + beta*15.0_dp
    !
    zz = self%fdmat%zi(2)
    !  ii=2
    sin_pi = dsin(pi * zz)
    sin_pi_hlf = dsin(pi_hlf * zz)
    cos_pi_hlf = dcos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
    beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
    sin_pi = sin_pi*sin_pi ! ^2
    !
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
    !   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
    gama = -(a2c*sin_pi/sin_pi_hlf)
    self%fdmat%d(2) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
    self%fdmat%g(2) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp
    !
    self%fdmat%H1(7, 2) = (-450.0_dp)     + beta*(-60.0_dp)
    self%fdmat%H1(6, 3) = ( 240.0_dp)     + beta*(180.0_dp)
    self%fdmat%H1(5, 4) = ( -15.0_dp)     + beta*(-45.0_dp)
    self%fdmat%H1(4, 5) =                 + beta*(6.0_dp)
    self%fdmat%H1(8, 1) = (240.0_dp)    + (beta*(-90.0_dp))  + gama*step_2*180.0_dp

    zz = self%fdmat%zi(3)
   !   ii=3
    sin_pi = dsin(pi * zz)
    sin_pi_hlf = dsin(pi_hlf * zz)
    cos_pi_hlf = dcos(pi_hlf * zz)
    sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
    cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
    beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
    sin_pi = sin_pi*sin_pi ! ^2
   !
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
!   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(3) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(3) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp
!
   self%fdmat%H1(7, 3) = (-490.0_dp)   + gama*step_2*180.0_dp
   self%fdmat%H1(6, 4) = (270.0_dp)    + beta*(135.0_dp)
   self%fdmat%H1(5, 5) = (-27.0_dp)    + beta*(-27.0_dp)
   self%fdmat%H1(4, 6) = (2.0_dp)      + beta*3.0_dp
   self%fdmat%H1(8, 2) = (270.0_dp)    + beta*(-135.0_dp)
   self%fdmat%H1(9, 1) = (-27.0_dp)    + beta*(27.0_dp)

   !
   do ii=4, N-3
      zz = self%fdmat%zi(ii)
      !
      sin_pi = dsin(pi * zz)
      sin_pi_hlf = dsin(pi_hlf * zz)
      cos_pi_hlf = dcos(pi_hlf * zz)
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
      beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
      sin_pi = sin_pi*sin_pi ! ^2
      !
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
!      gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
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

   !    ii=N-2
   zz = self%fdmat%zi(N-2)
   sin_pi = dsin(pi * zz)
   sin_pi_hlf = dsin(pi_hlf * zz)
   cos_pi_hlf = dcos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   !
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
!   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(N-2) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(N-2) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7, N-2) = (-490.0_dp)  + gama*step_2*180.0_dp ! H(N-2,N-2)
   self%fdmat%H1(6, N-1) = (270.0_dp)   + beta*(135.0_dp)! H(N-2,N-1)
   self%fdmat%H1(5, N) =  (-27.0_dp)    + beta*(-27.0_dp)! H(N-2, N)
   self%fdmat%H1(8, N-3) =(270.0_dp)    + beta*(-135.0_dp)! H(N-2,N-3)
   self%fdmat%H1(9, N-4) =(-27.0_dp)    + beta*(27.0_dp)! H(N-2,N-4)
   self%fdmat%H1(10, N-5) =(2.0_dp)     + beta*(-3.0_dp)! H(N-2,N-5)

   !    ii=N-1
   zz = self%fdmat%zi(N-1)
   sin_pi = dsin(pi * zz)
   sin_pi_hlf = dsin(pi_hlf * zz)
   cos_pi_hlf = dcos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   !
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
!   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(N-1) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(N-1) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7,N-1) = (-450.0_dp)   + beta*(60.0_dp)!H(N-1,N-1)
   self%fdmat%H1(6, N) =  (240.0_dp)    + gama*step_2*180.0_dp + beta*(90.0_dp)!H(N-1,N)
   self%fdmat%H1(8, N-2) =(240.0_dp)    + beta*(-180.0_dp)! H(N-1,N-2)
   self%fdmat%H1(9, N-3) =(-15.0_dp)    + beta*(45.0_dp)! H(N-1,N-3)
   self%fdmat%H1(10, N-4) =             + beta*(-6.0_dp)! H(N-1,N-4)

   !   ii=N
   zz = self%fdmat%zi(N)
   sin_pi = dsin(pi * zz)
   sin_pi_hlf = dsin(pi_hlf * zz)
   cos_pi_hlf = dcos(pi_hlf * zz)
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   !
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
!   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)
   gama = -(a2c*sin_pi/sin_pi_hlf)
   self%fdmat%d(N) = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf
   self%fdmat%g(N) = -llp1_pi_2_rm_4/sin_pi*step_2*180.0_dp

   self%fdmat%H1(7,N) =   (-300.0_dp)    + gama*step_2*180.0_dp + beta*150.0_dp!H(N,N)
   self%fdmat%H1(8, N-1) =(90.0_dp)      + beta*(-270.0_dp)!H(N,N-1)
   self%fdmat%H1(9, N-2) =(60.0_dp)      + beta*(90.0_dp)!H(N,N-2)
   self%fdmat%H1(10, N-3) = (-15.0_dp)   + beta*(-15.0_dp)!H(N,N-3)

   if(self%verbosity > 0) then
      write(*,'(a)') "OK"
   end if
  end subroutine integrator_precomp_fdmat

  !
  subroutine integrator_build_fdmat(self,ll)
    type(becke_integrator), intent(inout) :: self
    integer, intent(in) :: ll
    !
    integer :: ii,N
    real(dp) :: lll
    !
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
    !
    integer :: N,info
    !
    N=self%grid_params%N_radial
    rho_lm = rho_lm*self%fdmat%d
    !    call integrator_build_fdmat(self,ll)
    !    call DGBSV(N, 3, 3, 1, self%fdmat%H2, 10, self%fdmat%ipiv, rho_lm, N, info)
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
    !
    integer :: N,info
    !
    N=self%grid_params%N_radial
    rho_lm = rho_lm*self%fdmat%d
    !    call integrator_build_fdmat(self,ll)
    !    call DGBSV(N, 3, 3, 1, self%fdmat%H2, 10, self%fdmat%ipiv, rho_lm, N, info)
    info=0
    call DGBTRS( 'No transpose', N, 3, 3, 1, self%fdmat%H4(:,:,ll+1), 10, self&
        &%fdmat%ipiv4(:,ll+1), rho_lm, N, info)
  end subroutine integrator_solve_poisson

  ! initialize the integration module
  !
  subroutine integrator_init(self, grid_params)
    type(becke_integrator), intent(inout) :: self
    type(becke_grid_params), intent(in) :: grid_params
    !
    if(self%verbosity > 0) then
       write(*, '(a)') "=================================================="
       write(*, '(a)') "Initializing the integrator module"
       write(*, '(a, 3I4)') "Becke quadrature parameters: ", &
           &grid_params%N_radial, grid_params%N_angular, grid_params%ll_max
    end if
    !
    self%grid_params=grid_params
    self%kernel_parameter = 0.0_dp
    self%rm=grid_params%rm
    self%verbosity = 0
    if(self%verbosity > 0) then
       write(*,'(2(a,F12.7))') "rm=", self%rm, " kappa=", self%kernel_parameter
    end if
    !
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
    ! init anglib
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
    !
    !    call integrator_precomp_fdmat(self)
    !    call integrator_build_LU(self)
    !
    if(self%verbosity > 0) then
       write(*,'(a)') "integrator init done."
    end if
  end subroutine integrator_init

  ! precompute the LU decompositions
  subroutine integrator_build_LU(self)
    type(becke_integrator), intent(inout) :: self
    !
    integer :: N,ll,info

    N = self%grid_params%N_radial
    !
    do ll=1, self%grid_params%ll_max
       call integrator_build_fdmat(self,ll-1)
       self%fdmat%H3(:,:,ll) = self%fdmat%H2(:,:)
       call DGBTRF( N, N, 3, 3, self%fdmat%H3(:,:,ll), 10, self%fdmat%ipiv2(:,ll), info )
       if(info .ne. 0) write(*,'(a)') "ERROR: in LU decomposition!!!"
    end do
  end subroutine integrator_build_LU

  !
  subroutine integrator_set_dist(self, dist)
    type(becke_integrator), intent(inout) :: self
    real(dp), intent(in)  :: dist
    !
    call becke_grid_kill(self%integration_grid(2))
    call becke_grid_init(self%integration_grid(2), 2, &
        &self%radial_quadrature, self%angular_quadrature, dist, self%rm)
  end subroutine integrator_set_dist

  ! grid_no(grid_nr, subgrid_nr, coordinate_nr)
  ! grid_nr: 1 for 1_3 grid
  !          2 for 2_3 grid
  ! subgrid_nr: 1 for 1_3
  !             2 for 2_3
  ! coordiante_nr: 1 for rr
  !                2 for theta
  !                3 for phi
  subroutine integrator_get_coords(self, grid_no, coords)
    type(becke_integrator), target, intent(in) :: self
    integer, intent(in) :: grid_no(3)
    real(dp), pointer, intent(out) :: coords(:)
    ! build in the error checking
    coords=>self%integration_grid(grid_no(1))%subgrid(grid_no(2))%data(:,grid_no(3))
  end subroutine integrator_get_coords

  subroutine integrator_set_kernel_param(self, kappa)
    type(becke_integrator), intent(inout) :: self
    real(dp), intent(in) :: kappa
    self%kernel_parameter = kappa
  end subroutine integrator_set_kernel_param

  function integrator_process_1cnt(self, cnt1, type) result(res)
    type(becke_integrator), intent(in) :: self
    real(dp), intent(in) :: cnt1(:)
    integer, intent(in) :: type
    real(dp) :: res
    !
    select case (type)
    case (1)
       res = sum(self%integration_grid(1)%weight*cnt1)
    end select
  end function integrator_process_1cnt

  !
  function integrator_process_11(self, inner, outer) result(res)
    type(becke_integrator), intent(inout) :: self
    type(separable_integrand), intent(in) :: inner
    type(separable_integrand), intent(in) :: outer
    real(dp) :: res
    !
    integer :: ll,mm,N_radial,ll_max,ii,ind, count
    real(dp), pointer :: rr2(:)
    real(dp), allocatable :: rho1_lm(:)
    real(dp), allocatable :: V1_lm(:,:,:)
    real(dp), allocatable :: zi(:)
    real(dp), allocatable :: deriv(:)
    real(dp) :: gaunt_lm, kappa,yy2, gaunt2_lm, tmp, gaunt,charge
    integer :: alpha, beta
    real(dp), allocatable, target :: V(:)
    !
    logical :: first
    !
    res=0.0_dp

    N_radial=self%grid_params%N_radial
    ll_max=self%grid_params%ll_max
    !    ll_max=max(inner%ll1+inner%ll2,outer%ll1+outer%ll2)
    allocate(rho1_lm(N_radial))
    !    allocate(V1_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )
    !    allocate(deriv(N_radial))
    allocate(zi(N_radial))

    kappa=self%kernel_parameter
    zi = dacos(self%radial_quadrature%xx)/pi
    call integrator_get_coords(self,[3,1,1],rr2)

    allocate(V(N_radial))
    V = 0.0_dp

   ! charge=0.0_dp
   ! if(inner%ll1 == inner%ll2) then
   !    charge = sum(inner%radial*self%integration_grid(3)%weight)!/sqrt(4.0_dp*pi)
   ! end if

!       ll=0
 !      write(*,*) "L=", ll

!       rho1_lm = inner%radial
!       call solve_poisson(ll, rho1_lm, zi, charge)
!       call integrator_solve_helmholz(self,ll,rho1_lm)

!    res = sum(rho1_lm*outer%radial/rr2*self%integration_grid(3)%weight)/(4.0_dp*pi)

!       stop
!    do ll=0, ll_max-1
    do ll=abs(inner%ll1-inner%ll2), (inner%ll1+inner%ll2),2
       gaunt_lm=0.0_dp
       gaunt2_lm=0.0_dp
       rho1_lm = inner%radial
       !       call solve_helmholz_opt(ll, rho1_lm, zi, kappa)
       call integrator_solve_helmholz(self,ll,rho1_lm)
!                   call solve_poisson(ll, rho1_lm, zi, charge)
       do alpha=-inner%ll1, inner%ll1
          do beta=-inner%ll2, inner%ll2
             do mm=-ll,ll
                gaunt_lm = realGaunt(inner%ll1,alpha,inner%ll2,beta,ll,mm)
                gaunt_lm=gaunt_lm*gaunt_lm
                gaunt2_lm=gaunt2_lm + gaunt_lm
             end do
          end do
       end do
       V=V+rho1_lm * gaunt2_lm
    end do
    res = sum(V*outer%radial/rr2*self%integration_grid(3)%weight)


!              if(inner%ll1 == inner%ll2 .and. inner%mm1 == inner%mm2) then
!                 charge = sum(inner%radial*self%integration_grid(3)%weight)
!              end if


!     do ll=0, ll_max
!        first = .true.
!        gaunt = 0.0_dp
!        do mm=-ll, ll
!           gaunt_lm = realGaunt(inner%ll1,inner%mm1,inner%ll2,inner%mm2,ll,mm)
!           gaunt2_lm = realGaunt(outer%ll1,outer%mm1,outer%ll2,outer%mm2,ll,mm)
!           if(abs(gaunt_lm) < 1.0e-15_dp .and. abs(gaunt2_lm) < 1.0e-15_dp ) cycle
!           if(first) then
!              rho1_lm = inner%radial
! !                         call solve_helmholz_opt(ll, rho1_lm, zi, kappa)
!     !         call integrator_solve_helmholz(self,ll,rho1_lm)
!              charge=0.0_dp
!              call solve_poisson(ll, rho1_lm, zi, charge)
!              first = .false.
!           end if
!           gaunt = gaunt + gaunt2_lm * gaunt_lm
!        end do
!        V = V + rho1_lm * gaunt
!     end do
!     res = sum(V*outer%radial/rr2*self%integration_grid(3)%weight)

!     allocate(V(N_radial))
!     V = 0.0_dp
!     do ll=0, ll_max
!        do mm=-ll, ll
!           gaunt_lm = realGaunt(inner%ll1,inner%mm1,inner%ll2,inner%mm2,ll,mm)
!           gaunt2_lm = realGaunt(outer%ll1,outer%mm1,outer%ll2,outer%mm2,ll,mm)
!           if(abs(gaunt_lm) < 1.0e-15_dp .and. abs(gaunt2_lm) < 1.0e-15_dp ) cycle
!           rho1_lm = inner%radial * gaunt_lm
!           call solve_helmholz(ll, rho1_lm, zi, kappa)
! !          call ipl_tst(zi, rho1_lm, deriv)
! !          V1_lm(:,(ll*(ll+1) + mm + 1),1) = rho1_lm
! !          V1_lm(:,(ll*(ll+1) + mm + 1),2) = deriv
! !          ind = ll*(ll+1) + mm +1
! !          do ii=1, N_radial
! !             tmp = acos( (rr2(ii)-1.0_dp)/(rr2(ii)+1.0_dp) )/pi
! !             call get_cubic_spline(zi, V1_lm(:,ind,1), V1_lm(:,ind,2), tmp, yy2)
! !             V(ii) = V(ii) + yy2 * gaunt2_lm
! !          end do
!           V = V + rho1_lm* gaunt2_lm
!        end do
!     end do
  end function integrator_process_11

  !
  function integrator_process_2cnt(self, cnt1_1, cnt1_2, cnt2, type) result(res)
    type(becke_integrator), intent(in) :: self
    real(dp), intent(in) :: cnt1_1(:)
    real(dp), intent(in) :: cnt1_2(:)
    real(dp), intent(in) :: cnt2(:)
    integer, intent(in) :: type
    real(dp) :: res
    !
    real(dp), allocatable :: V1_lm(:,:,:)
    real(dp), allocatable, target :: zi(:), V(:)
    integer :: ll, mm, ngrid, ind, ii
    type(TRealTessY) :: tess
    real(dp) :: tmp, yy2
    !
    real(dp), pointer :: rr1(:), theta1(:), phi1(:)
    real(dp), pointer :: rr2(:), theta2(:), phi2(:)
    select case (type)
    case (1)
       call helmholtz_solve_1cnt_p2(self%radial_quadrature, &
           &self%angular_quadrature, self%grid_params%ll_max, &
           & self%integration_grid(1)%subgrid(1)%data,&
           & self%kernel_parameter, cnt1_1, V1_lm, zi)
       call integrator_get_coords(self, [2,2,1], rr2)
       call integrator_get_coords(self, [2,2,2], theta2)
       call integrator_get_coords(self, [2,2,3], phi2)
       ngrid=size(rr2)
       allocate(V(ngrid))
       V = 0.0_dp
       do ll=0, self%grid_params%ll_max
          do mm=-ll,ll
             ind = ll*(ll+1) + mm +1
             call TRealTessY_init(tess, ll, mm)
             do ii=1, ngrid
                tmp = acos( (rr2(ii)-1.0_dp)/(rr2(ii)+1.0_dp) )/pi
                call get_cubic_spline(zi, V1_lm(:,ind,1), V1_lm(:,ind,2), tmp, yy2)
                V(ii) = V(ii) + yy2/rr2(ii) * tess%getValue(theta2(ii), phi2(ii))
             end do
          end do
       end do
       res = sum(V*cnt2*cnt1_2*self%integration_grid(2)%weight)
 !      res = sum(V*cnt2*self%integration_grid(2)%weight)
       deallocate(V)
    end select
  end function integrator_process_2cnt

  subroutine becke_grid_init(self, n_subgrids, radquad, angquad, dist, rm)
    type(becke_grid), intent(inout) :: self
    integer, intent(in) :: n_subgrids
    type(TQuadrature), intent(in) :: radquad
    type(TQuadrature2D), intent(in) :: angquad
    real(dp), intent(in) :: dist
    real(dp), intent(in) :: rm
    !
    real(dp) :: beckepars(1)

    if(n_subgrids==11) then
       self%type = 11
!       if(self%verbosity > 0) then
!          write(*,'(a,F12.7)') "generating the one-center radial grid, rm=", rm
!       end if
       allocate( self%subgrid(1)) ! becke_grid contains only one subgrid
       call gengrid1_1(radQuad, rm, coordtrans_radial_becke2, self%subgrid(1)%data, self%weight)
    end if

    if(n_subgrids==1) then
       self%type = 1
 !      if(self%verbosity > 0) then
 !         write(*,'(a)') "generating the one-center 3d grid"
 !      end if
       allocate( self%subgrid(1)) ! becke_grid contains only one subgrid
       call gengrid1_3(angquad, radquad,&
           & coordtrans_radial_becke1, self%subgrid(1)%data, self%weight)
    end if
    if(n_subgrids==2) then
       self%type = 2
  !     if(self%verbosity > 0) then
  !        write(*,'(a, F12.7)') "generating the two-center 3d grid, dist=",&
  !            & dist
  !     end if
       allocate( self%subgrid(2)) ! becke_grid contains two subgrids
       call gengrid2_3(angquad, radquad, &
           &  coordtrans_radial_becke1, partition_becke_homo,beckepars,&
           & dist, self%subgrid(1)%data, self%subgrid(2)%data, &
           & self%weight, self%part)
    end if
  end subroutine becke_grid_init

  subroutine becke_grid_kill(self)
    type(becke_grid), intent(inout) :: self
    !
    integer :: ii

    do ii=1, size(self%subgrid)
       deallocate(self%subgrid(ii)%data)
    end do
    deallocate(self%subgrid)
    deallocate(self%weight)
    !deallocate(self%part)
  end subroutine becke_grid_kill

  !> solves the poisson equation
  !! \param ll: angular momentum
  !! \param rho: lm-component of density
  !! \param u: solution
  !! \param charge: integral over density(boundary condition)
  subroutine solve_poisson(ll, rho, zi, charge)
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho(:)
    !    real(dp), allocatable, intent(out) :: u(:)
    real(dp), allocatable, intent(in) :: zi(:)
    real(dp), intent(in) :: charge

    !**********************************************************************
    ! local variables
    !**********************************************************************
    integer :: ii, jj, info, ios
    integer,allocatable :: ipiv(:)
    integer :: N
    real(dp), allocatable :: H(:,:), alpha(:), beta(:), gama(:), delta(:)!,
    ! zi(:), ri(:)
    real(dp), allocatable :: H2(:,:), B(:)
    real(dp), allocatable :: solution(:), work(:), dummy(:)
    real, allocatable :: swork(:)
    real(dp) :: step, rr, rr_sqr, f0, f1, one_min_cos, one_pls_cos, sin_1, sin_2
    real(dp), parameter :: rm = 1.0_dp
    real(dp) :: tmp1, tmp2
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf
    integer :: stat, iter

    ! real(dp), allocatable :: AFB(:,:), SOL(:,:)
    ! real(dp), allocatable :: R(:), C(:), ERR_BNDS_NORM(:,:)
    ! real(dp), allocatable :: ERR_BNDS_COMP(:,:), PARAMS(:)! WORK(:)
    ! real(dp) :: RCOND, RPVGRW, BERR
    ! integer, allocatable :: IWORK(:)
    ! integer :: N_ERR_BNDS, NPARAMS
    ! character(1) :: EQUED

    ! allocate(AFB(10,N))
    ! allocate(R(N))
    ! allocate(C(N))
    ! allocate(SOL(N,1))
    ! allocate(WORK(4*N))
    ! allocate(IWORK(N))

    N=size(rho) ! number of grid points
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
    !    allocate(zi(N))
    !    allocate(ri(N))
    allocate(work(N))
    allocate(solution(N))
    allocate(swork(N*(N+1)))
    allocate(dummy(N))
    !    allocate(u(N))
    !******************************
    !
    !    allocate(H(N,N))
    allocate(H2(10,N))
    !    allocate(H(10,N))
    !    H = 0.0_dp
    H2 = 0.0_dp
    !*************************************************
    ! generate real and computational grids
    ! -> ri contains the real coordinates
    ! -> zi is the grid in computational domain
    !*************************************************
    !09.12.11 nicht noetig.
    !  do ii=1, N
    !     zi(ii) = real(ii, dp)/real(N+1, dp)
    !     ri(ii) = rm*(1.0_dp + cos(pi * zi(ii)))/(1.0_dp - cos(pi * zi(ii)))
    ! !    write(*,*) "# ", zi(ii), ri(ii)
    !  end do

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

       !**************************************************
!!! old version !!!
       !**************************************************
       ! one_pls_cos = 1.0_dp + cos(pi*zi(ii))
       ! one_min_cos = 2.0_dp - one_pls_cos
       ! sin_1 = sin(pi*zi(ii))
       ! sin_2 = sin_1*sin_1

       ! alpha(ii) = 1.0_dp

       ! beta(ii) = pi*(one_min_cos + sin_2)/(sin_1*one_min_cos)
       ! gama(ii) =  -rm*real(ll*(ll+1),dp)*sin_2/(one_min_cos*one_min_cos*one_pls_cos*one_pls_cos)*4.0_dp*pi*pi
       ! delta(ii) = -16.0_dp*pi*rm*rm*rm*pi*pi*sin_2*one_pls_cos*rho(ii)&
       !     &/(one_min_cos*one_min_cos*one_min_cos*one_min_cos*one_min_cos)
       !**************************************************
!!! old version !!!
       !**************************************************
    end do

    !******************************
    ! generate the FD Matrix
    !******************************
    !    call makeFDMatrix7P(H, delta, N, alpha, beta, gama, dummy, f0, f1) ! 7-point FD scheme, ref: Bickley
    ! call makeFDMatrix(H, delta, N, alpha, beta, gama, dummy, f0, f1)  ! usual 3-point FD scheme
    !    call makeBFDM7P(H2, delta, N, alpha, beta, gama, dummy, f0, f1)

    B = 0.0_dp
    call P_BFDM7P(H2, B, N, ll, zi, rm, charge)
    !    call TST_P_BFDM7P(H2, B, N, ll, zi, rm, charge)

    tmp1 = 1.0_dp/real(N+1, dp)
    tmp1 = tmp1*tmp1*180.0_dp
    B = B + tmp1*delta


    !******************************
    ! call LAPACK (d)sgesv routine
    ! to solve the linear eqn Hx=b
    !******************************
    !    call DSGESV(N, 1, H, N, ipiv, delta, N, solution, N, work, swork, iter, info)
    !call sgesv(N, 1, H, N, ipiv, delta, N, info)

    !***********************************************************
    ! solve the band diagonal matrix
    !***********************************************************
!!!==>
    call DGBSV(N, 3, 3, 1, H2, 10, ipiv, B, N, info)

    !N_ERR_BNDS = 1
    !    NPARAMS = 0
    !    PARAMS = 1

    ! call DGBSVXX ('E', 'N', N, 3, 3, 1, H2, 7, AFB, 10, IPIV,&
    !     & EQUED, R, C, B, N, SOL, N, RCOND, RPVGRW, BERR, N_ERR_BNDS,&
    !     & ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO)

    !    if(info == 0) write(*,*) "OK"

    do ii=1, N
       !             rho(ii) = solution(ii)
       rho(ii) =  B(ii) !delta(ii)
       !09.12.11: nicht noetig.   u(ii) = solution(ii)/ri(ii)
    end do

    deallocate(H2, alpha, beta, gama, ipiv, delta, work, swork, B)
  end subroutine solve_poisson



    !> solves the poisson equation
  !! \param ll: angular momentum
  !! \param rho: lm-component of density
  !! \param u: solution
  !! \param charge: integral over density(boundary condition)
  subroutine solve_helmholz_opt(ll, rho, zi, kappa)
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho(:)
    real(dp), allocatable, intent(in) :: zi(:)
    real(dp), intent(in) :: kappa

    !**********************************************************************
    ! local variables
    !**********************************************************************
    integer :: ii, jj, info, ios
    integer,allocatable :: ipiv(:)
    integer :: N
    real(dp), allocatable :: H(:,:), alpha(:), beta(:), gama(:), delta(:)!,
    ! zi(:), ri(:)
    real(dp), allocatable :: H2(:,:), B(:)
    real(dp), allocatable :: solution(:), work(:), dummy(:)
    real, allocatable :: swork(:)
    real(dp) :: step, rr, rr_sqr, f0, f1, one_min_cos, one_pls_cos, sin_1, sin_2
    real(dp), parameter :: rm = 1.0_dp
    real(dp) :: tmp1, tmp2
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf, a2c
    integer :: stat, iter

    N=size(rho) ! number of grid points
    !***********************************
    ! boundary conditions
    !***********************************
    f0 = 0.0_dp ! r->oo
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
    !
    allocate(H2(10,N))
    H2 = 0.0_dp
    !*************************************************
    ! generate real and computational grids
    ! -> ri contains the real coordinates
    ! -> zi is the grid in computational domain
    !*************************************************
    pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
 !   llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)
 !   a2c = kappa*kappa*rm*rm*pi*pi*0.25_dp
    do ii=1, N
!       sin_pi = sin(pi * zi(ii))
       sin_pi_hlf = sin(pi_hlf * zi(ii))
       cos_pi_hlf = cos(pi_hlf * zi(ii))
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
       cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
       cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
!       alpha(ii) = 1.0_dp
!       beta(ii) = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)
!       sin_pi = sin_pi*sin_pi ! ^2
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
       sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
!       gama(ii) = -(llp1_pi_2_rm_4/sin_pi + a2c*sin_pi/sin_pi_hlf)
       delta(ii) = -pi_rm_4_3*rho(ii)*cos_pi_hlf/sin_pi_hlf
    end do

    !******************************
    ! generate the FD Matrix
    !******************************
    B = 0.0_dp
    call H_BFDM7P(H2, B, N, ll, kappa, zi, rm)
    tmp1 = 1.0_dp/real(N+1, dp)
    tmp1 = tmp1*tmp1*180.0_dp
    B = B + tmp1*delta
    !***********************************************************
    ! solve the band diagonal matrix
    !***********************************************************
    call DGBSV(N, 3, 3, 1, H2, 10, ipiv, B, N, info)
 !   do ii=1, N
 !      rho(ii) =  B(ii) !delta(ii)
    !   end do
    rho = B
    deallocate(H2, alpha, beta, gama, ipiv, delta, work, swork, B)
  end subroutine solve_helmholz_opt

  !
  subroutine solve_helmholz(ll, rho, zi, kappa)
    integer, intent(in) :: ll
    real(dp), allocatable, intent(inout) :: rho(:)
    real(dp), allocatable, intent(in) :: zi(:)
    real(dp), intent(in) :: kappa

    integer :: ii, jj, info, ios
    integer, allocatable :: ipiv(:)
    integer :: N
    real(dp), allocatable :: H(:,:), alpha(:), beta(:), gama(:), delta(:)!, zi(:), ri(:)
    real(dp), allocatable :: solution(:), work(:), dummy(:)
    real, allocatable :: swork(:)
    real(dp) :: step, rr, rr_sqr, f0, f1, one_min_cos, one_pls_cos, sin_1, sin_2
    real(dp), parameter :: rm = 1.0_dp
    real(dp) :: tmp1, tmp2
    real(dp) :: pi_rm_4_3, llp1_pi_2_rm_4, sin_pi, sin_pi_hlf, cos_pi_hlf, a2c
    integer :: stat, iter
    N=size(rho)    ! number of grid points

    !********************************
    ! boundary conditions
    !********************************
    !  if(ll == 0) then
    !     f0 = sqrt(4.0_dp*pi)
    !  else
    !     f0=0.0_dp
    !  end if
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
    !    allocate(zi(N))
    !    allocate(ri(N))
    allocate(work(N))
    allocate(solution(N))
    allocate(swork(N*(N+1)))
    allocate(dummy(N))
    !    allocate(u(N))
    H = 0.0_dp

    !*************************************************
    ! generate real and computational grids
    ! -> ri contains the real coordinates
    ! -> zi is the grid in computational domain
    !*************************************************
    ! do ii=1, N
    !    zi(ii) = real(ii, dp)/real(N+1, dp)
    !    ri(ii) = rm*(1.0_dp + cos(pi * zi(ii)))/(1.0_dp - cos(pi * zi(ii)))
    ! end do

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
    !
    !******************************
    ! generate the FD Matrix
    !******************************
    call makeFDMatrix7P(H, delta, N, alpha, beta, gama, dummy, f0, f1) ! 7-point FD scheme, ref: Bickley
    !******************************
    ! call LAPACK (d)sgesv routine
    ! to solve the linear eqn Hx=b
    !******************************
    call DSGESV(N, 1, H, N, ipiv, delta, N, solution, N, work, swork, iter, info)
    !
!    call DGBSV(N, 3, 3, 1, H, 10, ipiv, B, N, info)
    !
    do ii=1, N
       rho(ii) = solution(ii)
       !       u(ii) = solution(ii)/ri(ii)
    end do
    deallocate(H, alpha, beta, gama, ipiv, delta, work, swork)
  end subroutine solve_helmholz



  !> Solves the one center poisson equation with Becke method
  !! \param N_rad: size of radial quadrature
  !! \param N_ang: size of angular quadrature
  !! \param ll_max: sperhical harmonics truncation
  !! \param dist: distance between centers
  !! \param func_cnt1: subroutine, which gives values of the 1. center
  !! \param func_cnt2: subroutine, which gives values of the 2. center
  !! \param V1_lm: array, which contains the expansion coefficients of
  !!  potential 1
  !! \param V2_lm: array, which contains the expansion coefficients of
  !!  potential 2
  !! \param zi: quadrature values in z-space(needed for interpolation)
  !! \ref beck paper
  !!
  subroutine helmholtz_solve_1cnt(N_radial, N_angular, ll_max, kappa, set_func_cnt1, V1_lm, zi)
    integer, intent(in) :: N_radial, N_angular, ll_max
    real(dp), intent(in) :: kappa
    interface
      subroutine set_func_cnt1(rr, theta, phi, gridval)
        use common_accuracy, only : dp
        real(dp), intent(in) :: rr(:), theta(:), phi(:)
        real(dp), intent(out) :: gridval(:)
      end subroutine set_func_cnt1
    end interface
    real(dp), allocatable, intent(out) :: V1_lm(:,:,:)
    real(dp), allocatable, intent(out) :: zi(:)

    type(TQuadrature2D) :: angQuad
    type(TQuadrature) :: radQuad
    real(dp), allocatable, target :: grid1(:,:)
    real(dp), allocatable, target :: weights(:)
    real(dp) :: charge1
    integer :: n1, n2, ngrid
    real(dp), pointer :: rr1(:), theta1(:), phi1(:)
    real(dp), pointer :: lcnt1(:), lcnt2(:), lprt(:), lwght(:)
    real(dp), allocatable, target :: center1(:), center2(:)
    real(dp), allocatable, target :: tess1(:), tess2(:)
    real(dp), allocatable :: rho1_lm(:), u1_lm(:)
    real(dp), allocatable :: rho2_lm(:), u2_lm(:)
    real(dp) :: beckepars(1)
    real(dp), allocatable :: res(:) !, zi(:)
    integer :: ii, lbound, hbound
    integer :: NN
    real(dp) :: sum, tmp11, rr, yy1, yy2, tmp
    integer :: ll, mm

    ! write(*,*) "#********************************************************"
    ! write(*,*) "# Solution of the two center poisson equation with "
    ! write(*,*) "# Becke's method. "
    ! write(*,*) "#********************************************************"

    NN = N_radial*N_angular

    ! write(*,*) "# Radial points: ", N_radial
    ! write(*,*) "# Angular points: ", N_angular
    ! write(*,*) "# l_max: ", ll_max
    ! write(*,*) "# Distance: ", dist

    ! write(*,*) "#========================================================"
    ! write(*,*) "# calculating charge..."

    !**********************************************************************
    ! calculate charge of one-center density
    !**********************************************************************
    call gauss_chebyshev_quadrature(N_radial, radQuad)
    call lebedev_laikov_quadrature(N_angular, angQuad)
    call gengrid1_3(angQuad, radQuad, coordtrans_radial_becke1, grid1, weights)

    rr1 => grid1(:,1)
    theta1 => grid1(:,2)
    phi1 => grid1(:,3)

    ngrid = size(rr1)
    allocate(center1(ngrid))
    call set_func_cnt1(rr1, theta1, phi1, center1)

    !    charge1 = sum(center1*weights)

    !     write (*,*) "# 1 -> ", charge1
    !**********************************************************************
    ! end calculate charge
    !**********************************************************************

    allocate(rho1_lm(N_radial))
    allocate(u1_lm(N_radial))
    allocate(V1_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )

    allocate(tess1(N_angular))

    allocate(res(N_radial))
    allocate(zi(N_radial))

    zi = acos(radQuad%xx)/pi

    do ll=0, ll_max
       do mm=-ll, ll
          do ii=1, N_radial
             lbound = (1+N_angular*(ii-1))
             hbound = (ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             theta1 => grid1( lbound:hbound, 2)
             phi1   => grid1( lbound:hbound, 3)
             call set_tess(theta1, phi1, ll, mm, tess1)
             rho1_lm(ii) = sum( (lcnt1)*tess1*angQuad%ww )
             !           write(*,*) grid1(lbound-NN,1), rho1_lm(ii), rho2_lm(ii)
          end do
          !         call solve_poisson(ll, rho1_lm, zi, charge1)
          call solve_helmholz(ll, rho1_lm, zi, kappa)

          ! write(*,*)
          ! do ii=1, N_radial
          !    lbound = (1+N_angular*(ii-1))
          !    write(*,*) grid1(lbound,1), rho1_lm(ii), u1_lm(ii)*grid1(lbound,1)
          ! end do

          ! calculate the second derivatives for the cubic spline interpolating
          ! function.
          !          call set_cubic_spline(zi, rho1_lm, 0.0_dp, 0.0_dp, res, .true.)
          call ipl_tst(zi, rho1_lm, res)

          V1_lm(:,(ll*(ll+1) + mm + 1),1) = rho1_lm
          V1_lm(:,(ll*(ll+1) + mm + 1),2) = res

       end do
    end do

    deallocate(center1, rho1_lm, tess1)
    deallocate(u1_lm)

  end subroutine helmholtz_solve_1cnt




  !> Solves the one center poisson equation with Becke method ver 2
  !! \param N_rad: size of radial quadrature
  !! \param N_ang: size of angular quadrature
  !! \param ll_max: sperhical harmonics truncation
  !! \param dist: distance between centers
  !! \param func_cnt1: subroutine, which gives values of the 1. center
  !! \param func_cnt2: subroutine, which gives values of the 2. center
  !! \param V1_lm: array, which contains the expansion coefficients of
  !!  potential 1
  !! \param V2_lm: array, which contains the expansion coefficients of
  !!  potential 2
  !! \param zi: quadrature values in z-space(needed for interpolation)
  !! \ref beck paper
  !!
  subroutine helmholtz_solve_1cnt_p2(radQuad, angQuad, ll_max, grid1, kappa, center1, V1_lm, zi)
    type(TQuadrature), intent(in) :: radQuad
    type(TQuadrature2D), intent(in) :: angQuad
    integer, intent(in) :: ll_max
    real(dp), target, intent(in) :: grid1(:,:)
    real(dp), intent(in) :: kappa
    ! interface
    !   subroutine set_func_cnt1(rr, theta, phi, gridval)
    !     use common_accuracy, only : dp
    !     real(dp), intent(in) :: rr(:), theta(:), phi(:)
    !     real(dp), intent(out) :: gridval(:)
    !   end subroutine set_func_cnt1
    ! end interface
    real(dp), target, intent(in) :: center1(:)
    real(dp), allocatable, intent(out) :: V1_lm(:,:,:)
    real(dp), allocatable, intent(out) :: zi(:)

    !    type(TQuadrature2D) :: angQuad
    !    type(TQuadrature) :: radQuad
    !    real(dp), allocatable, target :: grid1(:,:)
    integer :: N_radial, N_angular
    real(dp), allocatable, target :: weights(:)
    real(dp) :: charge1
    integer :: n1, n2, ngrid
    real(dp), pointer :: rr1(:), theta1(:), phi1(:)
    real(dp), pointer :: lcnt1(:), lcnt2(:), lprt(:), lwght(:)
    !    real(dp), allocatable, target :: center1(:), center2(:)
    real(dp), allocatable, target :: tess1(:), tess2(:)
    real(dp), allocatable :: rho1_lm(:), u1_lm(:)
    real(dp), allocatable :: rho2_lm(:), u2_lm(:)
    real(dp) :: beckepars(1)
    real(dp), allocatable :: res(:) !, zi(:)
    integer :: ii, lbound, hbound
    integer :: NN
    real(dp) :: sum, tmp11, rr, yy1, yy2, tmp
    integer :: ll, mm

    N_radial=size(radQuad%xx)
    N_angular=size(angQuad%c1)
    NN = N_radial*N_angular

    ! call gauss_chebyshev_quadrature(N_radial, radQuad)
    ! call lebedev_laikov_quadrature(N_angular, angQuad)
    ! call gengrid1_3(angQuad, radQuad, coordtrans_radial_becke1, grid1, weights)

    ! rr1 => grid1(:,1)
    ! theta1 => grid1(:,2)
    ! phi1 => grid1(:,3)

    ! ngrid = size(rr1)
    ! allocate(center1(ngrid))
    ! call set_func_cnt1(rr1, theta1, phi1, center1)
    !
    allocate(rho1_lm(N_radial))
    allocate(u1_lm(N_radial))
    allocate(V1_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )
    allocate(tess1(N_angular))
    allocate(res(N_radial))
    allocate(zi(N_radial))

    zi = acos(radQuad%xx)/pi

    do ll=0, ll_max
       do mm=-ll, ll
          do ii=1, N_radial
             lbound = (1+N_angular*(ii-1))
             hbound = (ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             theta1 => grid1( lbound:hbound, 2)
             phi1   => grid1( lbound:hbound, 3)
             call set_tess(theta1, phi1, ll, mm, tess1)
             rho1_lm(ii) = sum( (lcnt1)*tess1*angQuad%ww )
          end do
          call solve_helmholz(ll, rho1_lm, zi, kappa)

          ! calculate the second derivatives for the cubic spline interpolating
          ! function.
          call ipl_tst(zi, rho1_lm, res)

          V1_lm(:,(ll*(ll+1) + mm + 1),1) = rho1_lm
          V1_lm(:,(ll*(ll+1) + mm + 1),2) = res
       end do
    end do

    deallocate(rho1_lm, tess1)
    deallocate(u1_lm)

  end subroutine helmholtz_solve_1cnt_p2




  !> Solves the two center poisson equation with Becke method
  !! \param N_rad: size of radial quadrature
  !! \param N_ang: size of angular quadrature
  !! \param ll_max: sperhical harmonics truncation
  !! \param dist: distance between centers
  !! \param func_cnt1: subroutine, which gives values of the 1. center
  !! \param func_cnt2: subroutine, which gives values of the 2. center
  !! \param V1_lm: array, which contains the expansion coefficients of
  !!  potential 1
  !! \param V2_lm: array, which contains the expansion coefficients of
  !!  potential 2
  !! \param zi: quadrature values in z-space(needed for interpolation)
  !! \ref becke paper
  !!
  subroutine helmholtz_solve_2cnt(radQuad, angQuad, kappa, ll_max, dist, set_func_cnt1,&
      & set_func_cnt2, V1_lm, V2_lm, zi)
    type(TQuadrature), intent(in) :: radQuad
    type(TQuadrature2D), intent(in) :: angQuad
    real(dp), intent(in) :: kappa
    !    integer, intent(in) :: N_radial, N_angular
    integer, intent(in) :: ll_max
    real(dp), intent(in) :: dist
    interface
      subroutine set_func_cnt1(rr, theta, phi, gridval)
        use common_accuracy, only : dp
        real(dp), intent(in) :: rr(:), theta(:), phi(:)
        real(dp), intent(out) :: gridval(:)
      end subroutine set_func_cnt1
      subroutine set_func_cnt2(rr, theta, phi, gridval)
        use common_accuracy, only : dp
        real(dp), intent(in) :: rr(:), theta(:), phi(:)
        real(dp), intent(out) :: gridval(:)
      end subroutine set_func_cnt2
    end interface
    real(dp), allocatable, intent(out) :: V1_lm(:,:,:), V2_lm(:,:,:)
    real(dp), allocatable, intent(out) :: zi(:)

    !    type(TQuadrature2D) :: angQuad
    !    type(TQuadrature) :: radQuad
    integer  :: N_radial, N_angular
    real(dp), allocatable, target :: grid1(:,:), grid2(:,:)
    real(dp), allocatable, target :: weights(:), part(:)
    integer, allocatable :: ipiv(:)
    real(dp) :: charge1, charge2
    integer :: n1, n2, ngrid, info
    real(dp), pointer :: rr1(:), theta1(:), phi1(:)
    real(dp), pointer :: rr2(:), theta2(:), phi2(:)
    real(dp), pointer :: lcnt1(:), lcnt2(:), lprt(:), lwght(:)
    real(dp), allocatable, target :: center1(:), center2(:)
    real(dp), allocatable, target :: tess1(:), tess2(:)
    real(dp), allocatable :: rho1_lm(:)!, u1_lm(:)
    real(dp), allocatable :: rho2_lm(:)!, u2_lm(:)
    real(dp) :: beckepars(1)
    real(dp), allocatable :: res(:), res2(:) !, zi(:)
    integer :: ii, lbound, hbound
    integer :: NN
    real(dp) :: sum, tmp11, rr, yy1, yy2, tmp, tmp1, tmp2
    real(dp) :: sin_pi_hlf, cos_pi_hlf, pi_rm_4_3, rm
    integer :: ll, mm
    real(dp), allocatable :: B1(:), B2(:), H1(:,:), H2(:,:)

    rm = 1.0_dp

    N_radial = size(radQuad%xx)
    N_angular = size(angQuad%c1)
    NN = N_radial*N_angular

    !**********************************************************************
    ! calculate charge of two-center density
    !**********************************************************************
    !call gauss_chebyshev_quadrature(N_radial, radQuad)
    !call lebedev_laikov_quadrature(N_angular, angQuad)

    call gengrid2_3(angQuad, radQuad, coordtrans_radial_becke1, partition_becke_homo,&
        &beckepars, dist, grid1, grid2, weights, part)

    rr1 => grid1(:,1)
    theta1 => grid1(:,2)
    phi1 => grid1(:,3)

    rr2 => grid2(:,1)
    theta2 => grid2(:,2)
    phi2 => grid2(:,3)

    ngrid = size(rr1)
    allocate(center1(ngrid))
    allocate(center2(ngrid))
    ! note: there is no need to separate the integrands.
    !       The density in either coordinate system
    !       would also do it.

    ! call set_integrand(rr1, theta1, phi1, center1)
    ! call set_integrand(rr2, theta2, phi2, center2)
    call set_func_cnt1(rr1, theta1, phi1, center1)
    call set_func_cnt2(rr2, theta2, phi2, center2)

!!!===> Ab hier die eigentliche poisson funktionalitaet!!!
!!! Parameter: center1, center2, weights, NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lbound = 1
    hbound = NN
    lcnt1 => center1( lbound:hbound )
    lcnt2 => center2( lbound:hbound )
    lwght => weights( lbound:hbound )

!!!!=> replace the sum of integrands to product
    charge1 = sum( (lcnt1 * lcnt2) * lwght)

    lbound = NN+1
    hbound = 2*NN
    lcnt1 => center1( lbound:hbound )
    lcnt2 => center2( lbound:hbound )
    lwght => weights( lbound:hbound )

!!!!=> replace the sum of integrands to product
    charge2 = sum( (lcnt1 * lcnt2) * lwght)

    ! write (*,*) "# 1 -> ", charge1
    ! write (*,*) "# 2 -> ", charge2
    !**********************************************************************
    ! end calculate charge
    !**********************************************************************

    ! now two radial grids are set up
    ! grid(:,1) contains radial quadrature points according to the quadrature
    ! used. Now we should calculate the projections of density on the spherical
    ! harmonics.

    ! ngrid = size(theta)
    ! allocate(center1(ngrid))
    allocate(rho1_lm(N_radial))
    allocate(rho2_lm(N_radial))
    !    allocate(u1_lm(N_radial))
    !    allocate(u2_lm(N_radial))
    allocate(V1_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )
    allocate(V2_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )

    allocate(tess1(N_angular))
    allocate(tess2(N_angular))

    allocate(res(N_radial))
    allocate(res2(N_radial))
    allocate(zi(N_radial))

    allocate(B1(N_radial))
    allocate(B2(N_radial))
    allocate(H1(10,N_radial))
    allocate(H2(10,N_radial))

    allocate(ipiv(N_radial))

    !    zi = acos(radQuad%xx)/pi ! =grid1(:,4)
    zi = radQuad%zz

    !==>
    ! tmp1 = 1.0_dp/real(N_radial+1, dp)
    ! tmp1 = tmp1*tmp1*180.0_dp
    ! pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
    !==>

    do ll=0, ll_max
       !==>
       ! B1 = 0.0_dp
       ! B2 = 0.0_dp
       ! call P_BFDM7P(H1, B1, N_radial, ll, zi, rm, charge1)
       ! call P_BFDM7P(H2, B2, N_radial, ll, zi, rm, charge2)
       ! call DGBTRF(N_radial,N_radial,3,3,H1,10,ipiv,info)
       !       if(info /= 0) write(*,*) "DGTRF, ERR"
       ! call DGBTRF(N_radial,N_radial,3,3,H2,10,ipiv,info)
       !       if(info /= 0) write(*,*) "DGTRF, ERR"
       !==>
       do mm=-ll, ll
          do ii=1, N_radial
             ! take the values of the center1 and center2, corresponding to the
             ! rr, multiply by Y_lm and sum them up with Lebedev weights and
             ! partition function

             ! lcnt points to the portion of one center integrand
             ! which correponds to particular radial value
             lbound = (1+N_angular*(ii-1))
             hbound = (ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             lcnt2 => center2( lbound:hbound )
             theta1 => grid1( lbound:hbound, 2)
             phi1   => grid1( lbound:hbound, 3)
             lprt => part( lbound:hbound )
             call set_tess(theta1, phi1, ll, mm, tess1)
!!!!=> replace the sum of integrands to product
             rho1_lm(ii) = sum( (lcnt1 * lcnt2)*tess1*angQuad%ww*lprt )
             !==>
             !              sin_pi_hlf = sin(pi_hlf * zi(ii))
             !              cos_pi_hlf = cos(pi_hlf * zi(ii))
             !              sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
             !              cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
             !              cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
             !              sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
             !              sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
             !              tmp2 = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf*rho1_lm(ii)
             ! !             B1(ii) = B1(ii) + tmp1*rho1_lm(ii)
             !              rho1_lm(ii) = B1(ii) + tmp2
             !==>

             lbound = (NN+1+N_angular*(ii-1))
             hbound = (NN+ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             lcnt2 => center2( lbound:hbound )
             theta2 => grid2( lbound:hbound, 2)
             phi2   => grid2( lbound:hbound, 3)
             lprt => part( lbound:hbound )
             call set_tess(theta2, phi2, ll, mm, tess2)
!!!!=> replace the sum of integrands to product
             rho2_lm(ii) = sum( (lcnt1 * lcnt2)*tess2*angQuad%ww*lprt )
             !             tmp2 = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf*rho2_lm(ii)
             !             rho2_lm(ii) = B2(ii) + tmp2
             !           write(*,*) grid1(lbound-NN,1), rho1_lm(ii), rho2_lm(ii)

          end do

          ! rho1_lm, rho2_lm contain projections for both centers
          ! now solve the poisson equation
          ! ==>
          call solve_helmholz(ll, rho1_lm, zi, kappa)
          call solve_helmholz(ll, rho2_lm, zi, kappa)
          !          call DGBTRS('No transpose', N_radial, 3,3,1, H1, 10, ipiv, rho1_lm, N_radial, info)
          !          if(info /= 0) write(*,*) "DGTRS, ERR"
          !          call DGBTRS('No transpose', N_radial, 3,3,1, H2, 10, ipiv, rho2_lm, N_radial, info)
          !          if(info /= 0) write(*,*) "DGTRS, ERR"


          ! write(*,*)
          ! do ii=1, N_radial
          !    lbound = (1+N_angular*(ii-1))
          !    write(*,*) grid1(lbound,1), rho1_lm(ii), u1_lm(ii)*grid1(lbound,1)
          ! end do

          ! calculate the second derivatives for the cubic spline interpolating
          ! function.
          !          call set_cubic_spline(zi, rho1_lm, 0.0_dp, 0.0_dp, res, .true.)
          call ipl_tst(zi, rho1_lm, res)
          !          call  interpolation_scseq(zi, rho1_lm, res)
          !         res = res+res2
          V1_lm(:,(ll*(ll+1) + mm + 1),1) = rho1_lm
          V1_lm(:,(ll*(ll+1) + mm + 1),2) = res

          !          call set_cubic_spline(zi, rho2_lm, 0.0_dp, 0.0_dp, res, .true.)
          call ipl_tst(zi, rho2_lm, res)
          !          call  interpolation_scseq(zi, rho2_lm, res)
          !          res = res+res2

          V2_lm(:,(ll*(ll+1) + mm + 1),1) = rho2_lm
          V2_lm(:,(ll*(ll+1) + mm + 1),2) = res

       end do
    end do

    deallocate(center1, center2, rho1_lm, rho2_lm, tess1, tess2)
    !    deallocate(u1_lm, u2_lm, tess2)

  end subroutine helmholtz_solve_2cnt

!!!***************************************************************




  !> Solves the two center poisson equation with Becke method
  !! \param N_rad: size of radial quadrature
  !! \param N_ang: size of angular quadrature
  !! \param ll_max: sperhical harmonics truncation
  !! \param dist: distance between centers
  !! \param func_cnt1: subroutine, which gives values of the 1. center
  !! \param func_cnt2: subroutine, which gives values of the 2. center
  !! \param V1_lm: array, which contains the expansion coefficients of
  !!  potential 1
  !! \param V2_lm: array, which contains the expansion coefficients of
  !!  potential 2
  !! \param zi: quadrature values in z-space(needed for interpolation)
  !! \ref beck paper
  !!
  subroutine poisson_solve_2cnt(radQuad, angQuad, ll_max, dist, set_func_cnt1,&
      & set_func_cnt2, V1_lm, V2_lm, zi)
    type(TQuadrature), intent(in) :: radQuad
    type(TQuadrature2D), intent(in) :: angQuad
    !    integer, intent(in) :: N_radial, N_angular
    integer, intent(in) :: ll_max
    real(dp), intent(in) :: dist
    interface
      subroutine set_func_cnt1(rr, theta, phi, gridval)
        use common_accuracy, only : dp
        real(dp), intent(in) :: rr(:), theta(:), phi(:)
        real(dp), intent(out) :: gridval(:)
      end subroutine set_func_cnt1
      subroutine set_func_cnt2(rr, theta, phi, gridval)
        use common_accuracy, only : dp
        real(dp), intent(in) :: rr(:), theta(:), phi(:)
        real(dp), intent(out) :: gridval(:)
      end subroutine set_func_cnt2
    end interface
    real(dp), allocatable, intent(out) :: V1_lm(:,:,:), V2_lm(:,:,:)
    real(dp), allocatable, intent(out) :: zi(:)

    !    type(TQuadrature2D) :: angQuad
    !    type(TQuadrature) :: radQuad
    integer  :: N_radial, N_angular
    real(dp), allocatable, target :: grid1(:,:), grid2(:,:)
    real(dp), allocatable, target :: weights(:), part(:)
    integer, allocatable :: ipiv(:)
    real(dp) :: charge1, charge2
    integer :: n1, n2, ngrid, info
    real(dp), pointer :: rr1(:), theta1(:), phi1(:)
    real(dp), pointer :: rr2(:), theta2(:), phi2(:)
    real(dp), pointer :: lcnt1(:), lcnt2(:), lprt(:), lwght(:)
    real(dp), allocatable, target :: center1(:), center2(:)
    real(dp), allocatable, target :: tess1(:), tess2(:)
    real(dp), allocatable :: rho1_lm(:)!, u1_lm(:)
    real(dp), allocatable :: rho2_lm(:)!, u2_lm(:)
    real(dp) :: beckepars(1)
    real(dp), allocatable :: res(:), res2(:) !, zi(:)
    integer :: ii, lbound, hbound
    integer :: NN
    real(dp) :: sum, tmp11, rr, yy1, yy2, tmp, tmp1, tmp2
    real(dp) :: sin_pi_hlf, cos_pi_hlf, pi_rm_4_3, rm
    integer :: ll, mm
    real(dp), allocatable :: B1(:), B2(:), H1(:,:), H2(:,:)

    rm = 1.0_dp

    N_radial = size(radQuad%xx)
    N_angular = size(angQuad%c1)
    NN = N_radial*N_angular

    !**********************************************************************
    ! calculate charge of two-center density
    !**********************************************************************
    !call gauss_chebyshev_quadrature(N_radial, radQuad)
    !call lebedev_laikov_quadrature(N_angular, angQuad)

    call gengrid2_3(angQuad, radQuad, coordtrans_radial_becke1, partition_becke_homo,&
        &beckepars, dist, grid1, grid2, weights, part)

    rr1 => grid1(:,1)
    theta1 => grid1(:,2)
    phi1 => grid1(:,3)

    rr2 => grid2(:,1)
    theta2 => grid2(:,2)
    phi2 => grid2(:,3)

    ngrid = size(rr1)
    allocate(center1(ngrid))
    allocate(center2(ngrid))
    ! note: there is no need to separate the integrands.
    !       The density in either coordinate system
    !       would also do it.

    ! call set_integrand(rr1, theta1, phi1, center1)
    ! call set_integrand(rr2, theta2, phi2, center2)
    call set_func_cnt1(rr1, theta1, phi1, center1)
    call set_func_cnt2(rr2, theta2, phi2, center2)

!!!===> Ab hier die eigentliche poisson funktionalitaet!!!
!!! Parameter: center1, center2, weights, NN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lbound = 1
    hbound = NN
    lcnt1 => center1( lbound:hbound )
    lcnt2 => center2( lbound:hbound )
    lwght => weights( lbound:hbound )

!!!!=> replace the sum of integrands to product
    charge1 = sum( (lcnt1 * lcnt2) * lwght)

    lbound = NN+1
    hbound = 2*NN
    lcnt1 => center1( lbound:hbound )
    lcnt2 => center2( lbound:hbound )
    lwght => weights( lbound:hbound )

!!!!=> replace the sum of integrands to product
    charge2 = sum( (lcnt1 * lcnt2) * lwght)

    ! write (*,*) "# 1 -> ", charge1
    ! write (*,*) "# 2 -> ", charge2
    !**********************************************************************
    ! end calculate charge
    !**********************************************************************

    ! now two radial grids are set up
    ! grid(:,1) contains radial quadrature points according to the quadrature
    ! used. Now we should calculate the projections of density on the spherical
    ! harmonics.

    ! ngrid = size(theta)
    ! allocate(center1(ngrid))
    allocate(rho1_lm(N_radial))
    allocate(rho2_lm(N_radial))
    !    allocate(u1_lm(N_radial))
    !    allocate(u2_lm(N_radial))
    allocate(V1_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )
    allocate(V2_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )

    allocate(tess1(N_angular))
    allocate(tess2(N_angular))

    allocate(res(N_radial))
    allocate(res2(N_radial))
    allocate(zi(N_radial))

    allocate(B1(N_radial))
    allocate(B2(N_radial))
    allocate(H1(10,N_radial))
    allocate(H2(10,N_radial))

    allocate(ipiv(N_radial))

    !    zi = acos(radQuad%xx)/pi ! =grid1(:,4)
    zi = radQuad%zz

    !==>
    ! tmp1 = 1.0_dp/real(N_radial+1, dp)
    ! tmp1 = tmp1*tmp1*180.0_dp
    ! pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
    !==>

    do ll=0, ll_max
       !==>
       ! B1 = 0.0_dp
       ! B2 = 0.0_dp
       ! call P_BFDM7P(H1, B1, N_radial, ll, zi, rm, charge1)
       ! call P_BFDM7P(H2, B2, N_radial, ll, zi, rm, charge2)
       ! call DGBTRF(N_radial,N_radial,3,3,H1,10,ipiv,info)
       !       if(info /= 0) write(*,*) "DGTRF, ERR"
       ! call DGBTRF(N_radial,N_radial,3,3,H2,10,ipiv,info)
       !       if(info /= 0) write(*,*) "DGTRF, ERR"
       !==>
       do mm=-ll, ll
          do ii=1, N_radial
             ! take the values of the center1 and center2, corresponding to the
             ! rr, multiply by Y_lm and sum them up with Lebedev weights and
             ! partition function

             ! lcnt points to the portion of one center integrand
             ! which correponds to particular radial value
             lbound = (1+N_angular*(ii-1))
             hbound = (ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             lcnt2 => center2( lbound:hbound )
             theta1 => grid1( lbound:hbound, 2)
             phi1   => grid1( lbound:hbound, 3)
             lprt => part( lbound:hbound )
             call set_tess(theta1, phi1, ll, mm, tess1)
!!!!=> replace the sum of integrands to product
             rho1_lm(ii) = sum( (lcnt1 * lcnt2)*tess1*angQuad%ww*lprt )
             !==>
             !              sin_pi_hlf = sin(pi_hlf * zi(ii))
             !              cos_pi_hlf = cos(pi_hlf * zi(ii))
             !              sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
             !              cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
             !              cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
             !              sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
             !              sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
             !              tmp2 = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf*rho1_lm(ii)
             ! !             B1(ii) = B1(ii) + tmp1*rho1_lm(ii)
             !              rho1_lm(ii) = B1(ii) + tmp2
             !==>

             lbound = (NN+1+N_angular*(ii-1))
             hbound = (NN+ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             lcnt2 => center2( lbound:hbound )
             theta2 => grid2( lbound:hbound, 2)
             phi2   => grid2( lbound:hbound, 3)
             lprt => part( lbound:hbound )
             call set_tess(theta2, phi2, ll, mm, tess2)
!!!!=> replace the sum of integrands to product
             rho2_lm(ii) = sum( (lcnt1 * lcnt2)*tess2*angQuad%ww*lprt )
             !             tmp2 = -pi_rm_4_3*cos_pi_hlf/sin_pi_hlf*rho2_lm(ii)
             !             rho2_lm(ii) = B2(ii) + tmp2
             !           write(*,*) grid1(lbound-NN,1), rho1_lm(ii), rho2_lm(ii)

          end do

          ! rho1_lm, rho2_lm contain projections for both centers
          ! now solve the poisson equation
          ! ==>
          call solve_poisson(ll, rho1_lm, zi, charge1)
          call solve_poisson(ll, rho2_lm, zi, charge2)
          !          call DGBTRS('No transpose', N_radial, 3,3,1, H1, 10, ipiv, rho1_lm, N_radial, info)
          !          if(info /= 0) write(*,*) "DGTRS, ERR"
          !          call DGBTRS('No transpose', N_radial, 3,3,1, H2, 10, ipiv, rho2_lm, N_radial, info)
          !          if(info /= 0) write(*,*) "DGTRS, ERR"


          ! write(*,*)
          ! do ii=1, N_radial
          !    lbound = (1+N_angular*(ii-1))
          !    write(*,*) grid1(lbound,1), rho1_lm(ii), u1_lm(ii)*grid1(lbound,1)
          ! end do

          ! calculate the second derivatives for the cubic spline interpolating
          ! function.
          !          call set_cubic_spline(zi, rho1_lm, 0.0_dp, 0.0_dp, res, .true.)
          call ipl_tst(zi, rho1_lm, res)
          !          call  interpolation_scseq(zi, rho1_lm, res)
          !         res = res+res2
          V1_lm(:,(ll*(ll+1) + mm + 1),1) = rho1_lm
          V1_lm(:,(ll*(ll+1) + mm + 1),2) = res

          !          call set_cubic_spline(zi, rho2_lm, 0.0_dp, 0.0_dp, res, .true.)
          call ipl_tst(zi, rho2_lm, res)
          !          call  interpolation_scseq(zi, rho2_lm, res)
          !          res = res+res2

          V2_lm(:,(ll*(ll+1) + mm + 1),1) = rho2_lm
          V2_lm(:,(ll*(ll+1) + mm + 1),2) = res

       end do
    end do

    deallocate(center1, center2, rho1_lm, rho2_lm, tess1, tess2)
    !    deallocate(u1_lm, u2_lm, tess2)

  end subroutine poisson_solve_2cnt


  !> Solves the one center poisson equation with Becke method
  !! \param N_rad: size of radial quadrature
  !! \param N_ang: size of angular quadrature
  !! \param ll_max: sperhical harmonics truncation
  !! \param dist: distance between centers
  !! \param func_cnt1: subroutine, which gives values of the 1. center
  !! \param func_cnt2: subroutine, which gives values of the 2. center
  !! \param V1_lm: array, which contains the expansion coefficients of
  !!  potential 1
  !! \param V2_lm: array, which contains the expansion coefficients of
  !!  potential 2
  !! \param zi: quadrature values in z-space(needed for interpolation)
  !! \ref beck paper
  !!
  subroutine poisson_solve_1cnt(N_radial, N_angular, ll_max, set_func_cnt1, V1_lm, zi)
    integer, intent(in) :: N_radial, N_angular, LL_max
    interface
      subroutine set_func_cnt1(rr, theta, phi, gridval)
        use common_accuracy, only : dp
        real(dp), intent(in) :: rr(:), theta(:), phi(:)
        real(dp), intent(out) :: gridval(:)
      end subroutine set_func_cnt1
    end interface
    real(dp), allocatable, intent(out) :: V1_lm(:,:,:)
    real(dp), allocatable, intent(out) :: zi(:)

    type(TQuadrature2D) :: angQuad
    type(TQuadrature) :: radQuad
    real(dp), allocatable, target :: grid1(:,:)
    real(dp), allocatable, target :: weights(:)
    real(dp) :: charge1
    integer :: n1, n2, ngrid
    real(dp), pointer :: rr1(:), theta1(:), phi1(:)
    real(dp), pointer :: lcnt1(:), lcnt2(:), lprt(:), lwght(:)
    real(dp), allocatable, target :: center1(:), center2(:)
    real(dp), allocatable, target :: tess1(:), tess2(:)
    real(dp), allocatable :: rho1_lm(:), u1_lm(:)
    real(dp), allocatable :: rho2_lm(:), u2_lm(:)
    real(dp) :: beckepars(1)
    real(dp), allocatable :: res(:) !, zi(:)
    integer :: ii, lbound, hbound
    integer :: NN
    real(dp) :: sum, tmp11, rr, yy1, yy2, tmp
    integer :: ll, mm

    ! write(*,*) "#********************************************************"
    ! write(*,*) "# Solution of the two center poisson equation with "
    ! write(*,*) "# Becke's method. "
    ! write(*,*) "#********************************************************"

    NN = N_radial*N_angular

    ! write(*,*) "# Radial points: ", N_radial
    ! write(*,*) "# Angular points: ", N_angular
    ! write(*,*) "# l_max: ", ll_max
    ! write(*,*) "# Distance: ", dist

    ! write(*,*) "#========================================================"
    ! write(*,*) "# calculating charge..."

    !**********************************************************************
    ! calculate charge of one-center density
    !**********************************************************************
    call gauss_chebyshev_quadrature(N_radial, radQuad)
    call lebedev_laikov_quadrature(N_angular, angQuad)
    call gengrid1_3(angQuad, radQuad, coordtrans_radial_becke1, grid1, weights)

    rr1 => grid1(:,1)
    theta1 => grid1(:,2)
    phi1 => grid1(:,3)

    ngrid = size(rr1)
    allocate(center1(ngrid))
    call set_func_cnt1(rr1, theta1, phi1, center1)

    charge1 = sum(center1*weights)

    !     write (*,*) "# 1 -> ", charge1
    !**********************************************************************
    ! end calculate charge
    !**********************************************************************

    allocate(rho1_lm(N_radial))
    allocate(u1_lm(N_radial))
    allocate(V1_lm(N_radial, (ll_max*(ll_max+2))+1, 2) )

    allocate(tess1(N_angular))

    allocate(res(N_radial))
    allocate(zi(N_radial))

    zi = acos(radQuad%xx)/pi

    do ll=0, ll_max
       do mm=-ll, ll
          do ii=1, N_radial
             lbound = (1+N_angular*(ii-1))
             hbound = (ii*N_angular)
             lcnt1 => center1( lbound:hbound )
             theta1 => grid1( lbound:hbound, 2)
             phi1   => grid1( lbound:hbound, 3)
             call set_tess(theta1, phi1, ll, mm, tess1)
             rho1_lm(ii) = sum( (lcnt1)*tess1*angQuad%ww )
             !           write(*,*) grid1(lbound-NN,1), rho1_lm(ii), rho2_lm(ii)
          end do
          call solve_poisson(ll, rho1_lm, zi, charge1)

          ! write(*,*)
          ! do ii=1, N_radial
          !    lbound = (1+N_angular*(ii-1))
          !    write(*,*) grid1(lbound,1), rho1_lm(ii), u1_lm(ii)*grid1(lbound,1)
          ! end do

          ! calculate the second derivatives for the cubic spline interpolating
          ! function.
          !          call set_cubic_spline(zi, rho1_lm, 0.0_dp, 0.0_dp, res, .true.)

          call ipl_tst(zi, rho1_lm, res)
          V1_lm(:,(ll*(ll+1) + mm + 1),1) = rho1_lm
          V1_lm(:,(ll*(ll+1) + mm + 1),2) = res

       end do
    end do

    deallocate(center1, rho1_lm, tess1)
    deallocate(u1_lm)

  end subroutine poisson_solve_1cnt

  subroutine set_tess(theta, phi, ll, mm, gridval)
    real(dp), intent(in) :: theta(:), phi(:)
    integer, intent(in) :: ll, mm
    real(dp), intent(out) :: gridval(:)

    type(TRealTessY) :: yy
    call TRealTessY_init(yy, ll, mm)
    gridval = yy%getValue(theta, phi)
  end subroutine set_tess

end module common_poisson
