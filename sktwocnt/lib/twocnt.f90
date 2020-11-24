!> Contains the twocenter integrator routines.
module twocnt
  use omp_lib
  use accuracy
  use constants
  use quadratures
  use coordtrans
  use gridorbital
  use sphericalharmonics
  use gridgenerator
  use partition
  use dftxc
  use fifo_module
  implicit none
  private

  public :: twocnt_in, atomdata, integmap
  public :: get_twocenter_integrals

  ! Data associated with atoms
  type atomdata
    integer :: nbasis
    integer, allocatable :: angmoms(:)
    type(gridorb2), allocatable  :: rad(:), drad(:), ddrad(:)
    type(gridorb2) :: pot, rho, drho, ddrho
  end type atomdata

  !> Parsed input for twocnt.
  type twocnt_in
    logical :: hetero
    logical :: density
    integer :: ixc
    real(dp) :: r0, dr, epsilon, maxdist
    integer :: ninteg1, ninteg2
    type(atomdata) :: atom1, atom2
  end type twocnt_in


  !> Type for mapping integrals.
  type integmap
    !> Nr. of all nonzero twocenter integrals between orbitals of two atoms.
    integer :: ninteg

    !> Indicates for every integral the integrands:
    !! 
    !! o type(1,ii): index of orbital on first atom for integral ii.
    !! o type(2,ii): index of orbital on second atom for integral ii
    !! o type(3,ii): interaction type for integral ii: (0 - sigma, 1 - pi, ...)
    integer, allocatable :: type(:,:)

    !> Indicates which integral corresponds to a given (i1, i2, mm) combination,
    !! where i1 and i2 are the orbital indices on the two atoms and mm the
    !! interaction type. If the integral vanishes, the corresponding elemet is 0.
    integer, allocatable :: index(:,:,:)
  contains
    procedure :: init => integmap_init
  end type integmap

  
contains

  subroutine get_twocenter_integrals(inp, imap, skham, skover)
    type(twocnt_in), target, intent(in) :: inp
    type(integmap), intent(out) :: imap
    real(dp), allocatable, intent(out) :: skham(:,:), skover(:,:)
    
    type(quadrature) :: quads(2)

    type(atomdata), pointer :: atom1, atom2
    type(fifo_real2) :: hamfifo, overfifo
    real(dp), allocatable :: grid1(:,:), grid2(:,:)
    real(dp), allocatable :: dots(:), weights(:)
    real(dp), allocatable :: denserr(:)
    real(dp), allocatable :: skhambuffer(:,:), skoverbuffer(:,:)
    real(dp) :: beckepars(1)
    real(dp) :: dist, maxdist, denserrmax, maxabs
    integer :: ir, nbatch, nbatchline
    logical :: converged, dynlen

    call gauss_legendre_quadrature(inp%ninteg1, quads(1))
    call gauss_legendre_quadrature(inp%ninteg2, quads(2))

    atom1 => inp%atom1
    if (inp%hetero) then
      atom2 => inp%atom2
    else
      atom2 => inp%atom1
    end if
    call imap%init(atom1, atom2)

    ! Calculate lines for 1 Bohr in one batch.
    dist = 0.0_dp
    dynlen = (inp%maxdist > 0.0_dp)
    if (dynlen) then
      nbatchline = ceiling(1.0_dp / inp%dr)
      maxdist = inp%maxdist + real(nbatchline, dp) * inp%dr
    else
      maxdist = abs(inp%maxdist)
      nbatchline = ceiling((maxdist - inp%r0) / inp%dr)
    end if
    nbatch = 0
    denserrmax = 0.0_dp
    allocate(denserr(nbatchline))
    do
      allocate(skhambuffer(imap%ninteg, nbatchline))
      allocate(skoverbuffer(imap%ninteg, nbatchline))
      write(*, "(A,I0,A,F6.3,A,F6.3)") "Calculating ", nbatchline,&
          & " lines: r0 = ", inp%r0 + inp%dr * real(nbatch * nbatchline, dp),&
          & " dr = ", inp%dr
      do ir = 1, nbatchline
        dist = inp%r0 + inp%dr * real(nbatch * nbatchline + ir - 1, dp)
        call gengrid2_12(quads, coordtrans_becke_12, partition_becke,&
            & beckepars, dist, grid1, grid2, dots, weights)
        call getskintegrals(atom1, atom2, grid1, grid2, dots, weights,&
            &inp%density, inp%ixc, imap, skhambuffer(:,ir), skoverbuffer(:,ir),&
            & denserr(ir))
      end do
      denserrmax = max(denserrmax, maxval(denserr))
      maxabs = max(maxval(abs(skhambuffer)), maxval(abs(skoverbuffer)))
      if (dynlen) then
        converged = (maxabs < inp%epsilon)
        ! If new batch gave no contributions above tolerance: omit it and exit
        if (converged .or. dist > maxdist) then
          exit
        end if
        nbatch = nbatch + 1
        call hamfifo%push_alloc(skhambuffer)
        call overfifo%push_alloc(skoverbuffer)
      else
        converged = .true.
        call hamfifo%push_alloc(skhambuffer)
        call overfifo%push_alloc(skoverbuffer)
        exit
      end if
    end do
    if (.not. converged) then
      write(*, "(A,F6.2,A,ES10.3)") "Warning, maximal distance ", inp%maxdist,&
          & " reached! Max integral value:", maxabs
    end if
    write(*, "(A,ES10.3)") "Maximal integration error:", denserrmax

    call hamfifo%popall_concat(skham)
    call overfifo%popall_concat(skover)

  end subroutine get_twocenter_integrals


  !> Calculate SK-integrals.
  subroutine getskintegrals(atom1, atom2, grid1, grid2, dots, weights,&
      & densitysuper, ixc, imap, skham, skover, denserr)
    type(atomdata), intent(in) :: atom1, atom2
    real(dp), intent(in), target :: grid1(:,:), grid2(:,:), dots(:), weights(:)
    logical, intent(in) :: densitysuper
    integer, intent(in) :: ixc
    type(integmap), intent(in) :: imap
    real(dp), intent(out) :: skham(:), skover(:), denserr

    type(realtess) :: tes1, tes2
    real(dp), pointer :: r1(:), r2(:), theta1(:), theta2(:)
    real(dp), allocatable :: radval1(:,:)
    real(dp), allocatable :: radval2(:,:), radval2p(:,:), radval2pp(:,:)
    real(dp), allocatable :: potval(:), densval(:)
    real(dp), allocatable :: densval1p(:), densval1pp(:)
    real(dp), allocatable :: densval2p(:), densval2pp(:)
    real(dp), allocatable :: spherval1(:), spherval2(:)
    real(dp), allocatable :: absgr(:), laplace(:), gr_grabsgr(:)
    
    real(dp) :: integ1, integ2, dens, prefac
    integer :: ngrid
    integer :: ii, i1, i2, l1, l2, mm

    r1 => grid1(:,1)
    theta1 => grid1(:,2)
    r2 => grid2(:,1)
    theta2 => grid2(:,2)
    ngrid = size(r1)

    allocate(radval1(ngrid, atom1%nbasis))
    allocate(radval2(ngrid, atom2%nbasis))
    allocate(radval2p(ngrid, atom2%nbasis))
    allocate(radval2pp(ngrid, atom2%nbasis))
    allocate(spherval1(ngrid))
    allocate(spherval2(ngrid))
    do ii = 1, size(radval1, dim=2)
      radval1(:,ii) = getvalue(atom1%rad(ii), r1)
    end do
    do ii = 1, size(radval2, dim=2)
      radval2(:,ii) = getvalue(atom2%rad(ii), r2)
      radval2p(:,ii) = getvalue(atom2%drad(ii), r2)
      radval2pp(:,ii) = getvalue(atom2%ddrad(ii), r2)
    end do

    allocate(potval(ngrid))
    ifPotSup: if (.not. densitysuper) then
      potval = getvalue(atom1%pot, r1) + getvalue(atom2%pot, r2)
    else
      allocate(densval(ngrid))
      densval = getvalue(atom1%rho, r1) + getvalue(atom2%rho, r2)
      select case(ixc)
      case(1)
        call getxcpot_ldapw91(densval, potval)
      case(2)
        allocate(densval1p(ngrid))
        allocate(densval1pp(ngrid))
        allocate(densval2p(ngrid))
        allocate(densval2pp(ngrid))
        densval1p = getvalue(atom1%drho, r1)
        densval1pp = getvalue(atom1%ddrho, r1)
        densval2p = getvalue(atom2%drho, r2)
        densval2pp = getvalue(atom2%ddrho, r2)
        allocate(absgr(ngrid))
        allocate(laplace(ngrid))
        allocate(gr_grabsgr(ngrid))
        ! Calculate derivatives for combined density
        call getderivs(densval1p, densval1pp, densval2p, densval2pp, r1, r2,&
            &dots, absgr, laplace, gr_grabsgr)
        ! Get XC potential
        call getxcpot_ggapbe(densval, absgr, laplace, gr_grabsgr, potval)
      case default
        write(*,*) "Unknown functional type"
        stop
      end select
      ! Add nuclear and coulomb potential
      potval = potval + getvalue(atom1%pot, r1) + getvalue(atom2%pot, r2)
    end if ifPotSup

    denserr = 0.0_dp
    do ii = 1, imap%ninteg
      i1 = imap%type(1, ii)
      l1 = atom1%angmoms(i1)
      i2 = imap%type(2, ii)
      l2 = atom2%angmoms(i2)
      mm = imap%type(3, ii) - 1
      call init(tes1, l1, mm)
      call init(tes2, l2, mm)
      spherval1 = getvalue_1d(tes1, theta1)
      spherval2 = getvalue_1d(tes2, theta2)
      integ1 = gethamiltonian(radval1(:,i1), radval2(:,i2), &
          &radval2p(:,i2), radval2pp(:,i2), r2, l2, spherval1, &
          &spherval2, potval, weights)
      integ2 = getoverlap(radval1(:,i1), radval2(:,i2), spherval1, &
          &spherval2, weights)
      dens = getdensity(radval1(:,i1), radval2(:,i2), spherval1, &
          &spherval2, weights)
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






  function getoverlap(rad1, rad2, spher1, spher2, weights) result(res)
    real(dp), intent(in) :: rad1(:), rad2(:), spher1(:), spher2(:), weights(:)
    real(dp) :: res

    res = sum(rad1 * rad2 * spher1 * spher2 * weights)
    
  end function getoverlap


  function getdensity(rad1, rad2, spher1, spher2, weights) result(res)
    real(dp), intent(in) :: rad1(:), rad2(:), spher1(:), spher2(:), weights(:)
    real(dp) :: res

    res = sum(((rad1 * spher1)**2 + (rad2 * spher2)**2) * weights)
    
  end function getdensity


  function gethamiltonian(rad1, rad2, rad2p, rad2pp, r2, l2, spher1, spher2, &
      &pot, weights) result(res)
    real(dp), intent(in) :: rad1(:), rad2(:), rad2p(:), rad2pp(:), r2(:)
    integer, intent(in) :: l2
    real(dp), intent(in) :: spher1(:), spher2(:), pot(:), weights(:)
    real(dp) :: res
    
    res = sum((rad1 * spher1) &
        &* (-0.5_dp * rad2pp &
        &- rad2p / r2 &
        &+ 0.5_dp * l2 * (l2 + 1) * rad2 / r2**2&
        &+ pot * rad2) &
        &* spher2 * weights)
    
  end function gethamiltonian



  subroutine getderivs(drho1, d2rho1, drho2, d2rho2, r1, r2, dots, &
      &absgr, laplace, gr_grabsgr)
    real(dp), intent(in) :: drho1(:), d2rho1(:), drho2(:), d2rho2(:)
    real(dp), intent(in) :: r1(:), r2(:), dots(:)
    real(dp), intent(out) :: absgr(:), laplace(:), gr_grabsgr(:)

    integer :: nn
    real(dp), allocatable :: f1(:), f2(:)

    nn = size(drho1)
    allocate(f1(nn), f2(nn))
    f1 = drho1 + dots * drho2
    f2 = drho2 + dots * drho1
    absgr = sqrt(drho1 * f1 + drho2 * f2)
    laplace = d2rho1 + d2rho2 + 2.0_dp * (drho1 / r1 + drho2 / r2)
    where (absgr > epsilon(1.0_dp))
      gr_grabsgr = (d2rho1 * f1 * f1 +  d2rho2 * f2 * f2 &
          &+(1.0_dp - dots**2) * drho1 * drho2 * (drho2 / r1 + drho1 / r2)) &
          &/ absgr
    elsewhere
      gr_grabsgr = 0.0_dp
    end where
        
  end subroutine getderivs

  
  !> Initializes the twocenter integration map based on the basis on two atoms.
  !! \param self  Instance.
  !! \param atom1  Properties of atom1.
  !! \param atom2  Properties of atom2.
  subroutine integmap_init(self, atom1, atom2)
    class(integmap), intent(out) :: self
    type(atomdata), intent(in) :: atom1, atom2

    integer :: mmax, ninteg, ind, i1, l1, i2, l2, mm
  
    mmax = min(maxval(atom1%angmoms), maxval(atom2%angmoms))
    allocate(self%index(atom1%nbasis, atom2%nbasis, mmax+1))
    self%index = 0
    ninteg = 0
    do i1 = 1, atom1%nbasis
      l1 = atom1%angmoms(i1)
      do i2 = 1, atom2%nbasis
        l2 = atom2%angmoms(i2)
        do mm = 0, min(l1, l2)
          print *, l1, l2, mm
          ninteg = ninteg + 1
          self%index(i1, i2, mm+1) = ninteg
        end do
      end do
    end do
    self%ninteg = ninteg
    allocate(self%type(3, ninteg))
    ind = 0
    do i1 = 1, atom1%nbasis
      l1 = atom1%angmoms(i1)
      do i2 = 1, atom2%nbasis
        l2 = atom2%angmoms(i2)
        do mm = 1, min(l1, l2) + 1
          ind = ind + 1
          self%type(:, ind) = [ i1, i2, mm ]
        end do
      end do
    end do

  end subroutine integmap_init


end module twocnt
    
