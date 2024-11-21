!> Module that defines all global variables of the project.
module globals

  use common_accuracy, only : dp
  use mixer, only : TMixer, TMixer_init, TMixer_reset, mixerTypes
  use broydenmixer, only : TBroydenMixer, TBroydenMixer_init
  use simplemixer, only : TSimpleMixer, TSimpleMixer_init
  use diismixer, only : TDiisMixer, TDiisMixer_init

  implicit none

  !> confinement radii
  real(dp) :: conf_r0(0:4)

  !> power of confinement
  real(dp) :: conf_power(0:4)

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

  !> scf tolerance, i.e. convergence criteria
  real(dp) :: scftol

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

  !> (hf) exchange supermatrix (long-range, range-separated version)
  real(dp), allocatable :: kk_lr(:,:,:,:,:,:)

  !> wavefunction coefficients
  real(dp), allocatable :: cof(:,:,:,:)

  !> maximum (absolute) value in [F,PS] 
  real(dp) :: commutator_max

  !> density matrix supervector
  real(dp), allocatable :: pp(:,:,:,:), pp_old(:,:,:,:)

  !> density matrix max difference
  real(dp) :: pp_diff

  !> fock matrix supervector
  real(dp), allocatable :: ff(:,:,:,:)

  !> commutator [F, PS]
  real(dp), allocatable :: commutator(:,:,:,:)

  !> potential matrix supervectors
  real(dp), allocatable :: pot_new(:,:,:,:), pot_old(:,:,:,:)

  !> eigenvalues
  real(dp), allocatable :: eigval(:,:,:), eigval_old(:,:,:)

  !> eigenvalue difference vector norm
  real(dp) :: eigval_diff

  !> zora scaled eigenvalues
  real(dp), allocatable :: eigval_scaled(:,:,:)

  !> total energy
  real(dp) :: total_ene, total_ene_old, total_ene_diff

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

  !> kinetic energy density on grid
  real(dp), allocatable :: tau(:,:)

  !> orbital-dependent tau potential on grid
  real(dp), allocatable :: vtau(:,:)

  !> xc potential on grid
  real(dp), allocatable :: vxc(:,:)

  !> exc energy density on grid
  real(dp), allocatable :: exc(:)

  !> average local effective potential equivalent to non-local GKS potential (if present)
  real(dp), allocatable :: avgPot(:,:)

  !> generate alphas automatically
  logical :: tAutoAlphas

  !> print eigenvectors to stdout
  logical :: tPrintEigvecs

  !> true, if zero-order regular approximation for relativistic effects is desired
  logical :: tZora

  !> true, if SCF cycle reached convergency on a given quantity
  logical :: tCommutatorConverged, tEnergyConverged, tEigenspectrumConverged, tDensityConverged

  !> identifier of mixer
  integer :: mixnr

  !> DIIS subspace size
  integer, parameter :: n_vec_diis = 10

  !> mixer instance
  type(TMixer), allocatable :: pMixer

  !> simple mixer (if used)
  type(TSimpleMixer), allocatable :: pSimpleMixer

  !> broyden mixer (if used)
  type(TBroydenMixer), allocatable :: pBroydenMixer

  !> DIIS mixer (if used)
  type(TDiisMixer), allocatable :: pDiisMixer

  !> mixing factor
  real(dp) :: mixing_factor

  !> true, if average local, effective potential should be calculated
  logical :: isAvgPotNeeded

  !> zora kinetic energy contribution to total energy
  real(dp) :: zora_ekin

  !> maximum size of the eigenproblem
  integer :: problemsize


contains

  !> Allocates all the variables in the globals module.
  subroutine allocate_globals()

    !! auxiliary variables to count the elements to mix
    integer :: ind, idx1, idx2, idx3, idx4

    allocate(weight(num_mesh_points))
    allocate(abcissa(num_mesh_points))
    allocate(dzdr(num_mesh_points))
    allocate(d2zdr2(num_mesh_points))
    allocate(rho(num_mesh_points, 2))
    allocate(drho(num_mesh_points, 2))
    allocate(ddrho(num_mesh_points, 2))
    allocate(tau(num_mesh_points, 2))
    allocate(vtau(num_mesh_points, 2))
    allocate(exc(num_mesh_points))
    allocate(vxc(num_mesh_points, 2))

    allocate(ss(0:max_l, problemsize, problemsize))
    write(*, '(A,I0,A)') 'Size of one Supervectors is ', size(ss), ' double precision elements'

    allocate(uu(0:max_l, problemsize, problemsize))
    allocate(tt(0:max_l, problemsize, problemsize))
    allocate(vconf(0:max_l, problemsize, problemsize))
    allocate(ff(2, 0:max_l, problemsize, problemsize))
    allocate(commutator(2, 0:max_l, problemsize, problemsize))
    allocate(pot_old(2, 0:max_l, problemsize, problemsize))
    allocate(pot_new(2, 0:max_l, problemsize, problemsize))

    allocate(eigval(2, 0:max_l, problemsize))
    allocate(eigval_old(2, 0:max_l, problemsize))
    allocate(eigval_scaled(2, 0:max_l, problemsize))

    allocate(jj(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize))
    write(*, '(A,I0,A)') 'Size of one Supermatrix is ', size(jj), ' double precision elements'

    allocate(kk(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize))
    allocate(kk_lr(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize))

    write(*, '(A,I3)') 'MAXIMUM SIZE OF EIGENPROBLEM IS ', problemsize
    write(*, '(A)') ' '

    ! first index reserved for spin
    ! fourth index of cof is the eigenvalue index, see densmatrix build
    allocate(cof(2, 0:max_l, problemsize, problemsize))
    allocate(pp(2, 0:max_l, problemsize, problemsize))
    allocate(pp_old(2, 0:max_l, problemsize, problemsize))

    weight(:) = 0.0_dp
    abcissa(:) = 0.0_dp
    rho(:,:) = 0.0_dp
    drho(:,:) = 0.0_dp
    ddrho(:,:) = 0.0_dp
    tau(:,:) = 0.0_dp
    vtau(:,:) = 0.0_dp

    eigval(:,:,:) = 0.0_dp
    eigval_scaled(:,:,:) = 0.0_dp

    cof(:,:,:,:) = 0.0_dp
    pp(:,:,:,:) = 0.0_dp

    if (isAvgPotNeeded) allocate(avgPot(num_mesh_points, 2), source=0.0_dp)

    ! initialize mixer
    allocate(pMixer)
    select case(mixnr)
    case(mixerTypes%simple)
      allocate(pSimplemixer)
      call TSimpleMixer_init(pSimpleMixer, mixing_factor)
      call TMixer_init(pMixer, pSimpleMixer)
    case(mixerTypes%broyden)
      allocate(pBroydenMixer)
      ! defaults taken from DFTB+
      call TBroydenMixer_init(pBroydenMixer, maxiter, mixing_factor, 0.01_dp, 1.0_dp, 1.0e5_dp,&
          & 1.0e-2_dp)
      call TMixer_init(pMixer, pBroydenMixer)
    case(mixerTypes%diis)
      allocate(pDiisMixer)
      call TDiisMixer_init(pDiisMixer, n_vec_diis, mixing_factor, .false.)
      call TMixer_init(pMixer, pDiisMixer)
    case default
      error stop "Unknown mixer type."
    end select

    ! count elements to mix
    ind = 0
    do idx1 = 1, 2
      do idx2 = 0, max_l
        do idx3 = 1, num_alpha(idx2) * poly_order(idx2)
          do idx4 = 1, problemsize
            ind = ind + 1
          end do
        end do
      end do
    end do

    call TMixer_reset(pMixer, ind)

  end subroutine allocate_globals

end module globals
