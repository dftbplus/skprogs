module globals

  use common_accuracy, only : dp

  implicit none

  real(dp) :: conf_r0(0:4) ! confinement radius
  integer  :: conf_power(0:4) ! power of confinement
  real(dp) :: alpha(0:4,10) ! exponents
  integer  :: occ_shells(0:4) ! number of occupied shells
  integer  :: num_alpha(0:4) ! number of exponents in each shell
  integer  :: poly_order(0:4) ! highest polynomial order + l in each shell
  integer  :: nuc ! nuclear charge
  integer  :: max_l ! maximum angular momentum
  integer  :: maxiter ! maximum number of SCF calculations
  logical  :: generate_alpha ! generate alphas automatically
  logical  :: eigprint ! print eigenvectors to stdout
  real(dp) :: min_alpha ! smallest exponent if generate_alpha
  real(dp) :: max_alpha ! largest exponent if generate_alpha
  integer :: num_occ ! maximal occupied shell
  integer :: num_power ! maximum number of coefficients
  integer :: num_alphas ! maximum number of exponents
  real(dp), allocatable :: occ(:,:,:) ! occupation numbers

  real(dp), allocatable :: s(:,:,:) ! overlap supervector
  real(dp), allocatable :: u(:,:,:) ! nucleus-electron supervector
  real(dp), allocatable :: t(:,:,:) ! kinetic supervector
  real(dp), allocatable :: vconf(:,:,:) ! confinement supervector

  real(dp), allocatable :: j(:,:,:,:,:,:) ! coulomb supermatrix
  real(dp), allocatable :: k(:,:,:,:,:,:) ! (hf) exchange supermatrix

  real(dp), allocatable :: cof(:,:,:,:) ! wavefunction coefficients
  real(dp) :: change_max ! relative changes during scf
  real(dp), allocatable :: p(:,:,:,:)   ! density matrix supervector

  real(dp), allocatable :: f(:,:,:,:) ! fock matrix supervector
  real(dp), allocatable :: pot_new(:,:,:,:) ! potential matrix supervector
  real(dp), allocatable :: pot_old(:,:,:,:) ! potential matrix supervector

  real(dp), allocatable :: eigval(:,:,:) ! eigenvalues
  real(dp), allocatable :: eigval_scaled(:,:,:) ! zora scaled eigenvalues

  real(dp) :: total_ene,kinetic_energy,nuclear_energy,conf_energy
  real(dp) :: coulomb_energy,exchange_energy

  integer :: xcnr ! switch exchange-correlation
  real(dp) :: xalpha_const ! exchange parameter for X-Alpha exchange

  integer :: num_mesh_points ! number of numerical integration points
  real(dp), allocatable :: weight(:) ! numerical integration weights
  real(dp), allocatable :: abcissa(:) ! numerical integration abcissas
  real(dp), allocatable :: dzdr(:)  ! dz/dr
  real(dp), allocatable :: d2zdr2(:)  ! d2z/dr2
  real(dp) :: dz ! step width in linear coordinates
  real(dp), allocatable :: rho(:,:) ! density on grid
  real(dp), allocatable :: drho(:,:) ! 1st deriv. of density on grid
  real(dp), allocatable :: ddrho(:,:) ! 2nd deriv. of density on grid
  real(dp), allocatable :: vxc(:,:) ! xc potential on grid
  real(dp), allocatable :: exc(:) ! exc energy density on grid

  logical :: zora,final

  logical :: broyden ! switch broyden/simplemix
  real(dp) :: mixing_factor ! mixing factor
  real(dp) :: zora_ekin ! zora kinetic energy contribution to total energy

  integer :: problemsize

contains

  subroutine allocate_globals

    ! Allocate all the variables in the globals module

    allocate(weight(num_mesh_points))
    allocate(abcissa(num_mesh_points))
    allocate(dzdr(num_mesh_points))
    allocate(d2zdr2(num_mesh_points))
    allocate(rho(num_mesh_points,2))
    allocate(drho(num_mesh_points,2))
    allocate(ddrho(num_mesh_points,2))
    allocate(exc(num_mesh_points))
    allocate(vxc(num_mesh_points,2))

    allocate(s(0:max_l,problemsize,problemsize))
    allocate(u(0:max_l,problemsize,problemsize))
    allocate(t(0:max_l,problemsize,problemsize))
    allocate(vconf(0:max_l,problemsize,problemsize))
    allocate(f(2,0:max_l,problemsize,problemsize))
    allocate(pot_old(2,0:max_l,problemsize,problemsize))
    allocate(pot_new(2,0:max_l,problemsize,problemsize))
    write(*,'(A,I0,A)') 'Size of one Supervectors is ',size(s),' &
        &double precision elements'

    allocate(eigval(2,0:max_l,problemsize))
    allocate(eigval_scaled(2,0:max_l,problemsize))

    allocate(j(0:max_l,problemsize,problemsize,0:max_l,problemsize,problemsize))
    allocate(k(0:max_l,problemsize,problemsize,0:max_l,problemsize,problemsize))
    write(*,'(A,I0,A)') 'Size of one Supermatrix is ',size(j),' &
        &double precision elements'

    write(*,'(A,I3)') 'MAXIMUM SIZE OF EIGENPROBLEM IS ',problemsize
    write(*,'(A)') ' '

    ! first index reserved for spin
    ! fourth index of cof is the eigenvalue index, see densmatrix build
    allocate(cof(2,0:max_l,problemsize,problemsize))
    allocate(p(2,0:max_l,problemsize,problemsize))

    weight=0.0d0
    abcissa=0.0d0
    rho=0.0d0
    drho=0.0d0
    ddrho=0.0d0

    eigval=0.0d0
    eigval_scaled=0.0d0

    cof=0.0d0
    p=0.0d0

  end subroutine allocate_globals

end module globals
