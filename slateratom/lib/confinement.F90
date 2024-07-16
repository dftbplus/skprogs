#:include 'common.fypp'

!> Module that builds up various supervectors.
module confinement

  use common_accuracy, only : dp
  use utilities, only : fak
  use core_overlap, only : v
  use density, only : basis_times_basis_times_r2

  implicit none
  private

  public :: TConf, TConfInp, confType
  public :: TPowerConf, TPowerConf_init
  public :: TWsConf, TWsConf_init


  !> Generic wrapper for confinement potentials.
  type, abstract :: TConf
  contains

    procedure(getPotOnGrid), deferred :: getPotOnGrid
    procedure(getSupervec), deferred :: getSupervec

  end type TConf


  abstract interface

    !> Tabulates the (shell-resolved) confinement potential on a given grid.
    subroutine getPotOnGrid(this, max_l, num_mesh_points, abcissa, vconf)
      import :: TConf, dp
      implicit none

      !> instance
      class(TConf), intent(in) :: this

      !> maximum angular momentum
      integer, intent(in) :: max_l

      !> number of numerical integration points
      integer, intent(in) :: num_mesh_points

      !> numerical integration abcissas
      real(dp), intent(in) :: abcissa(:)

      !> confinement potential on grid
      real(dp), intent(out) :: vconf(:, 0:)

    end subroutine getPotOnGrid


    !> Calculates supervector matrix elements of the confining potential.
    subroutine getSupervec(this, max_l, num_mesh_points, abcissa, weight, num_alpha, alpha,&
        & poly_order, vconf, vconf_matrix)
      import :: TConf, dp
      implicit none

      !> instance
      class(TConf), intent(in) :: this

      !> maximum angular momentum
      integer, intent(in) :: max_l

      !> number of numerical integration points
      integer, intent(in) :: num_mesh_points

      !> numerical integration abcissas
      real(dp), intent(in) :: abcissa(:)

      !> numerical integration weights
      real(dp), intent(in) :: weight(:)

      !> number of exponents in each shell
      integer, intent(in) :: num_alpha(0:)

      !> basis exponents
      real(dp), intent(in) :: alpha(0:,:)

      !> highest polynomial order + l in each shell
      integer, intent(in) :: poly_order(0:)

      !> Woods-Saxon potential on grid
      real(dp), intent(in) :: vconf(:, 0:)

      !> Woods-Saxon confinement supervector
      real(dp), intent(out) :: vconf_matrix(0:,:,:)

    end subroutine getSupervec

  end interface


  !> Input for the Power confinement.
  type :: TPowerConfInp

    !> confinement radii (power compression)
    real(dp) :: r0(0:4)

    !> power of confinement (power compression)
    real(dp) :: power(0:4)

  end type TPowerConfInp


  !> Input for the Woods-Saxon confinement.
  type :: TWsConfInp

    !> potential heights (W)
    real(dp) :: ww(0:4)

    !> potential slopes (a)
    real(dp) :: aa(0:4)

    !> half-height radii (r0)
    real(dp) :: r0(0:4)

  end type TWsConfInp


  !> Input for the Power confinement.
  type :: TConfInp

    !> Power confinement input structure
    type(TPowerConfInp) :: power

    !> Woods-Saxon confinement input structure
    type(TWsConfInp) :: ws

  end type TConfInp


  !> Power confinement.
  type, extends(TConf) :: TPowerConf
    private

    !> confinement radii
    real(dp) :: r0(0:4)

    !> power of confinement
    real(dp) :: power(0:4)

  contains

    procedure :: getPotOnGrid => TPowerConf_getPotOnGrid
    procedure :: getSupervec => TPowerConf_getSupervec

  end type TPowerConf


  !> Woods-Saxon confinement.
  type, extends(TConf) :: TWsConf
    private

    !> potential heights (W)
    real(dp) :: ww(0:4)

    !> potential slopes (a)
    real(dp) :: aa(0:4)

    !> half-height radii (r0)
    real(dp) :: r0(0:4)

  contains

    procedure :: getPotOnGrid => TWsConf_getPotOnGrid
    procedure :: getSupervec => TWsConf_getSupervec

  end type TWsConf


  !> Enumerator for type of confinement potential.
  type :: TConfEnum

    !> no compression
    integer :: none = 0

    !> power compression
    integer :: power = 1

    !> Woods-Saxon compression
    integer :: ws = 2

  end type TConfEnum

  !> Container for enumerated types of confinement potentials.
  type(TConfEnum), parameter :: confType = TConfEnum()


contains

  !> Initializes a TPowerConf instance.
  subroutine TPowerConf_init(this, input)

    !> instance
    type(TPowerConf), intent(out) :: this

    !> input data
    type(TPowerConfInp), intent(inout) :: input

    this%r0(0:) = input%r0
    this%power(0:) = input%power

  end subroutine TPowerConf_init


  !> Initializes a TWsConf instance.
  subroutine TWsConf_init(this, input)

    !> instance
    type(TWsConf), intent(out) :: this

    !> input data
    type(TWsConfInp), intent(inout) :: input

    this%ww(0:) = input%ww
    this%aa(0:) = input%aa
    this%r0(0:) = input%r0

  end subroutine TWsConf_init


  !> Tabulates the (shell-resolved) Power confinement potential on the grid.
  subroutine TPowerConf_getPotOnGrid(this, max_l, num_mesh_points, abcissa, vconf)

    !> instance
    class(TPowerConf), intent(in) :: this

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> Woods-Saxon potential on grid
    real(dp), intent(out) :: vconf(:, 0:)

    !! angular momentum
    integer :: ll

    @:ASSERT(size(vconf, dim=1) == num_mesh_points)
    @:ASSERT(size(vconf, dim=2) == max_l + 1)

    do ll = 0, max_l
      vconf(:, ll) = getPowerPotential(abcissa, this%r0(ll), this%power(ll))
    end do

  end subroutine TPowerConf_getPotOnGrid


  !> Tabulates the (shell-resolved) Woods-Saxon confinement potential on the grid.
  subroutine TWsConf_getPotOnGrid(this, max_l, num_mesh_points, abcissa, vconf)

    !> instance
    class(TWsConf), intent(in) :: this

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> Woods-Saxon potential on grid
    real(dp), intent(out) :: vconf(:, 0:)

    !! angular momentum
    integer :: ll

    @:ASSERT(size(vconf, dim=1) == num_mesh_points)
    @:ASSERT(size(vconf, dim=2) == max_l + 1)

    do ll = 0, max_l
      vconf(:, ll) = getWSPotential(abcissa, this%ww(ll), this%aa(ll), this%r0(ll))
    end do

  end subroutine TWsConf_getPotOnGrid


  !> Tabulates the (shell-resolved) Power confinement potential on the grid.
  subroutine TPowerConf_getSupervec(this, max_l, num_mesh_points, abcissa, weight, num_alpha,&
      & alpha, poly_order, vconf, vconf_matrix)

    !> instance
    class(TPowerConf), intent(in) :: this

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> Power potential on grid
    real(dp), intent(in) :: vconf(:, 0:)

    !> Power confinement supervector
    real(dp), intent(out) :: vconf_matrix(0:,:,:)

    @:ASSERT(size(vconf, dim=1) == num_mesh_points)
    @:ASSERT(size(vconf, dim=2) == max_l + 1)

    @:ASSERT(size(vconf_matrix, dim=1) == max_l + 1)

    call build_dft_conf_matrix(max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa,&
        & weight, vconf, vconf_matrix)

  end subroutine TPowerConf_getSupervec


  !> Tabulates the (shell-resolved) Woods-Saxon confinement potential on the grid.
  subroutine TWsConf_getSupervec(this, max_l, num_mesh_points, abcissa, weight, num_alpha, alpha,&
      & poly_order, vconf, vconf_matrix)

    !> instance
    class(TWsConf), intent(in) :: this

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> Woods-Saxon potential on grid
    real(dp), intent(in) :: vconf(:, 0:)

    !> Woods-Saxon confinement supervector
    real(dp), intent(out) :: vconf_matrix(0:,:,:)

    @:ASSERT(size(vconf, dim=1) == num_mesh_points)
    @:ASSERT(size(vconf, dim=2) == max_l + 1)

    @:ASSERT(size(vconf_matrix, dim=1) == max_l + 1)

    call build_dft_conf_matrix(max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa,&
        & weight, vconf, vconf_matrix)

  end subroutine TWsConf_getSupervec


  !> Builds DFT confinement matrix to be added to the Fock matrix by calculating the single matrix
  !! elements and putting them together.
  subroutine build_dft_conf_matrix(max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa,&
      & weight, vconf, conf_matrix)

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

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> confinement potential on grid
    real(dp), intent(in) :: vconf(:, 0:)

    !> DFT confinement matrix
    real(dp), intent(out) :: conf_matrix(0:,:,:)

    !> single matrix element of the confinement potential
    real(dp) :: conf_matrixelement

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, ss, tt, start

    conf_matrix(:,:,:) = 0.0_dp
    conf_matrixelement = 0.0_dp

    do ii = 0, max_l
      ss = 0
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          ss = ss + 1

          tt = ss - 1
          do ll = jj, num_alpha(ii)

            start = 1
            if (ll == jj) start = kk

            do mm = start, poly_order(ii)
              tt = tt + 1

              call dft_conf_matrixelement(num_mesh_points, weight, abcissa, vconf(:, ii),&
                  & alpha(ii, jj), kk, alpha(ii, ll), mm, ii, conf_matrixelement)

              conf_matrix(ii, ss, tt) = conf_matrixelement
              conf_matrix(ii, tt, ss) = conf_matrixelement

            end do
          end do
        end do
      end do
    end do

  end subroutine build_dft_conf_matrix


  !> Calculates a single matrix element of the exchange correlation potential.
  pure subroutine dft_conf_matrixelement(num_mesh_points, weight, abcissa, vconf, alpha1, poly1,&
      & alpha2, poly2, ll, conf_matrixelement)

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> confinement potential on grid
    real(dp), intent(in) :: vconf(:)

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
    real(dp), intent(out) :: conf_matrixelement

    !> stores product of two basis functions and r^2
    real(dp) :: basis

    !> auxiliary variable
    integer :: ii

    conf_matrixelement = 0.0_dp

    do ii = 1, num_mesh_points
      basis = basis_times_basis_times_r2(alpha1, poly1, alpha2, poly2, ll, abcissa(ii))
      conf_matrixelement = conf_matrixelement + weight(ii) * vconf(ii) * basis
    end do

  end subroutine dft_conf_matrixelement


  !> Returns the Power potential for a given radial distance.
  elemental impure function getPowerPotential(rr, r0, power) result(pot)

    !> radial distance
    real(dp), intent(in) :: rr

    !> confinement radius
    real(dp), intent(in) :: r0

    !> confinement power
    real(dp), intent(in) :: power

    !> Power potential at evaluation point
    real(dp) :: pot

    @:ASSERT(rr >= 0.0_dp)

    pot = (rr / r0)**power

  end function getPowerPotential


  !> Returns the Woods-Saxon potential for a given radial distance.
  !! see J. Chem. Theory Comput. 12, 1, 53-64 (2016) eqn. 4.
  elemental impure function getWSPotential(rr, ww, aa, r0) result(pot)

    !> radial distance
    real(dp), intent(in) :: rr

    !> potential height (W)
    real(dp), intent(in) :: ww

    !> potential slope (a)
    real(dp), intent(in) :: aa

    !> half-height radius (r0)
    real(dp), intent(in) :: r0

    !> Woods-Saxon potential at evaluation point
    real(dp) :: pot

    @:ASSERT(rr >= 0.0_dp)

    pot = ww / (1.0 + exp(-aa * (rr - r0)))

  end function getWSPotential


  ! !> Calculates analytic matrix elements of confining potential.
  ! !! No checking for power, e.g. power==0 or power<0 etc. !
  ! subroutine getConfPower_analytical(this, max_l, num_alpha, alpha, poly_order, vconf_matrix)

  !   !> instance
  !   class(TPowerConf), intent(in) :: this

  !   !> maximum angular momentum
  !   integer, intent(in) :: max_l

  !   !> number of exponents in each shell
  !   integer, intent(in) :: num_alpha(0:)

  !   !> basis exponents
  !   real(dp), intent(in) :: alpha(0:,:)

  !   !> highest polynomial order + l in each shell
  !   integer, intent(in) :: poly_order(0:)

  !   !> confinement supervector
  !   real(dp), intent(out) :: vconf_matrix(0:,:,:)

  !   !! temporary storage
  !   real(dp) :: alpha1

  !   !! auxiliary variables
  !   integer :: ii, jj, kk, ll, mm, nn, oo, nlp, nlq

  !   vconf_matrix(:,:,:) = 0.0_dp

  !   do ii = 0, max_l
  !     if (this%power(ii) > 1.0e-06_dp) then
  !       nn = 0
  !       do jj = 1, num_alpha(ii)
  !         do ll = 1, poly_order(ii)
  !           nn = nn + 1
  !           oo = 0
  !           nlp = ll + ii
  !           do kk = 1, num_alpha(ii)
  !             alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
  !             do mm = 1, poly_order(ii)
  !               oo = oo + 1
  !               nlq = mm + ii
  !               vconf_matrix(ii, nn, oo) = 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
  !                   & * v(alpha(ii, kk), 2 * nlq)) / (this%r0(ii) * 2.0_dp)**this%power(ii)&
  !                   & * v(alpha1, nlp + nlq + this%power(ii))
  !             end do
  !           end do
  !         end do
  !       end do
  !     end if
  !   end do

  ! end subroutine TConf_getConfPower_analytical

end module confinement
