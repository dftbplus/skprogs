!> Module that defines all global variables of the project.
module globals

  use common_accuracy, only : dp

  implicit none

  !> confinement radii
  real(dp) :: conf_r0(0:4)

  !> power of confinement
  integer :: conf_power(0:4)

  !> basis exponents
  real(dp) :: alpha(0:4, 10)

  !> number of occupied shells
  integer :: occ_shells(0:4)

  !> number of exponents in each shell
  integer :: num_alpha(0:4)

  !> highest polynomial order + l in each shell
  integer :: poly_order(0:4)

  !> nuclear charge, i.e. atomic number
  integer :: nuc

  !> maximum angular momentum
  integer :: max_l

  !> maximum number of SCF calculations
  integer :: maxiter

  !> smallest exponent if generate_alpha
  real(dp) :: min_alpha

  !> largest exponent if generate_alpha
  real(dp) :: max_alpha

  !> maximal occupied shell
  integer :: num_occ

  !> maximum number of coefficients
  integer :: num_power

  !> maximum number of exponents
  integer :: num_alphas

  !> occupation numbers
  real(dp), allocatable :: occ(:,:,:)

  !> overlap supervector
  real(dp), allocatable :: ss(:,:,:)

  !> nucleus-electron supervector
  real(dp), allocatable :: uu(:,:,:)

  !> kinetic supervector
  real(dp), allocatable :: tt(:,:,:)

  !> confinement supervector
  real(dp), allocatable :: vconf(:,:,:)

  !> coulomb supermatrix
  real(dp), allocatable :: jj(:,:,:,:,:,:)

  !> (hf) exchange supermatrix
  real(dp), allocatable :: kk(:,:,:,:,:,:)

  !> wavefunction coefficients
  real(dp), allocatable :: cof(:,:,:,:)

  !> relative changes during scf
  real(dp) :: change_max

  !> density matrix supervector
  real(dp), allocatable :: pp(:,:,:,:)

  !> fock matrix supervector
  real(dp), allocatable :: ff(:,:,:,:)

  !> potential matrix supervectors
  real(dp), allocatable :: pot_new(:,:,:,:), pot_old(:,:,:,:)

  !> eigenvalues
  real(dp), allocatable :: eigval(:,:,:)

  !> zora scaled eigenvalues
  real(dp), allocatable :: eigval_scaled(:,:,:)

  !> total energy
  real(dp) :: total_ene

  !> kinetic energy
  real(dp) :: kinetic_energy

  !> nuclear energy
  real(dp) :: nuclear_energy

  !> confinement energy
  real(dp) :: conf_energy

  !> coulomb energy
  real(dp) :: coulomb_energy

  !> exchange energy
  real(dp) :: exchange_energy

  !> identifier of exchange-correlation type
  integer :: xcnr

  !> exchange parameter for X-Alpha exchange
  real(dp) :: xalpha_const

  !> number of numerical integration points
  integer :: num_mesh_points

  !> numerical integration weights
  real(dp), allocatable :: weight(:)

  !> numerical integration abcissas
  real(dp), allocatable :: abcissa(:)

  !> dz/dr
  real(dp), allocatable :: dzdr(:)

  !> d2z/dr2
  real(dp), allocatable :: d2zdr2(:)

  !> step width in linear coordinates
  real(dp) :: dz

  !> density on grid
  real(dp), allocatable :: rho(:,:)

  !> 1st deriv. of density on grid
  real(dp), allocatable :: drho(:,:)

  !> 2nd deriv. of density on grid
  real(dp), allocatable :: ddrho(:,:)

  !> xc potential on grid
  real(dp), allocatable :: vxc(:,:)

  !> exc energy density on grid
  real(dp), allocatable :: exc(:)

  !> generate alphas automatically
  logical :: tAutoAlphas

  !> print eigenvectors to stdout
  logical :: tPrintEigvecs

  !> true, if zero-order regular approximation for relativistic effects is desired
  logical :: tZora

  !> true, if SCF cycle reached convergency
  logical :: tConverged

  !> true, if Broyden mixing is desired, otherwise simple mixing is applied
  logical :: tBroyden

  !> mixing factor
  real(dp) :: mixing_factor

  !> zora kinetic energy contribution to total energy
  real(dp) :: zora_ekin

  !> maximum size of the eigenproblem
  integer :: problemsize


contains

  !> Allocates all the variables in the globals module.
  subroutine allocate_globals()

    allocate(weight(num_mesh_points))
    allocate(abcissa(num_mesh_points))
    allocate(dzdr(num_mesh_points))
    allocate(d2zdr2(num_mesh_points))
    allocate(rho(num_mesh_points, 2))
    allocate(drho(num_mesh_points, 2))
    allocate(ddrho(num_mesh_points, 2))
    allocate(exc(num_mesh_points))
    allocate(vxc(num_mesh_points, 2))

    allocate(ss(0:max_l, problemsize, problemsize))
    write(*, '(A,I0,A)') 'Size of one Supervectors is ', size(ss), ' double precision elements'

    allocate(uu(0:max_l, problemsize, problemsize))
    allocate(tt(0:max_l, problemsize, problemsize))
    allocate(vconf(0:max_l, problemsize, problemsize))
    allocate(ff(2, 0:max_l, problemsize, problemsize))
    allocate(pot_old(2, 0:max_l, problemsize, problemsize))
    allocate(pot_new(2, 0:max_l, problemsize, problemsize))

    allocate(eigval(2, 0:max_l, problemsize))
    allocate(eigval_scaled(2, 0:max_l, problemsize))

    allocate(jj(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize))
    write(*, '(A,I0,A)') 'Size of one Supermatrix is ', size(jj), ' double precision elements'

    allocate(kk(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize))

    write(*, '(A,I3)') 'MAXIMUM SIZE OF EIGENPROBLEM IS ', problemsize
    write(*, '(A)') ' '

    ! first index reserved for spin
    ! fourth index of cof is the eigenvalue index, see densmatrix build
    allocate(cof(2, 0:max_l, problemsize, problemsize))
    allocate(pp(2, 0:max_l, problemsize, problemsize))

    weight(:) = 0.0_dp
    abcissa(:) = 0.0_dp
    rho(:,:) = 0.0_dp
    drho(:,:) = 0.0_dp
    ddrho(:,:) = 0.0_dp

    eigval(:,:,:) = 0.0_dp
    eigval_scaled(:,:,:) = 0.0_dp

    cof(:,:,:,:) = 0.0_dp
    pp(:,:,:,:) = 0.0_dp

  end subroutine allocate_globals

end module globals
