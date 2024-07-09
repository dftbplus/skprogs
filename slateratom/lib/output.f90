!> Module that provides basic routines to write out various results.
module output

  use common_accuracy, only : dp, mc
  use common_constants, only : pi
  use core_overlap, only : moments
  use density, only : wavefunction, wavefunction_1st, wavefunction_2nd, density_at_point,&
      & density_at_point_1st
  use coulomb_potential, only : cou_pot
  use common_taggedout, only : TTaggedwriter, TTaggedwriter_init, writetag
  use utilities, only : fak

  implicit none
  private

  public :: write_energies, write_eigval, write_eigvec, write_moments
  public :: write_potentials_file_standard, write_densities_file_standard
  public :: write_waves_file_standard, cusp_values, write_energies_tagged
  public :: write_wave_coeffs_file
  public :: writeAveragePotential

  character(len=1), parameter :: orbnames(0:4) = ["s", "p", "d", "f", "g"]


contains

  subroutine write_energies(kinetic_energy, nuclear_energy, coulomb_energy, exchange_energy,&
      & xc_energy, conf_energy, total_ene, tZora)

    real(dp), intent(in) :: kinetic_energy, nuclear_energy, coulomb_energy
    real(dp), intent(in) :: exchange_energy, xc_energy, total_ene, conf_energy

    !> true, if zero-order regular approximation for relativistic effects is desired
    logical, intent(in) :: tZora

    write(*, '(A)') 'FINAL ENERGIES :                   '
    write(*, '(A)') '-----------------------------------'
    write(*, '(A)') ' '
    write(*, '(A,F22.6,A)') 'KINETIC ENERGY     ', kinetic_energy, ' Hartree'
    write(*, '(A,F22.6,A)') 'NUC. ATTR. ENERGY  ', nuclear_energy, ' Hartree'
    write(*, '(A,F22.6,A)') 'COULOMB ENERGY     ', 0.5_dp * coulomb_energy, ' Hartree'
    write(*, '(A,F22.6,A)') 'HF X ENERGY        ', 0.5_dp * exchange_energy, ' Hartree'
    write(*, '(A,F22.6,A)') 'DFT XC ENERGY      ', xc_energy, ' Hartree'
    write(*, '(A,F22.6,A)') 'CONFINEMENT ENERGY ', conf_energy, ' Hartree'
    write(*,*) ' '

    if (.not. tZora) then
      write(*, '(A,F22.6,A)') 'TOTAL ENERGY       ', total_ene, ' Hartree'
      write(*, '(A)') ' '
      write(*, '(A,F22.6)') 'DFT VIRIAL ', (nuclear_energy + 0.5_dp * coulomb_energy&
          & + 0.5_dp * exchange_energy + xc_energy + conf_energy) / kinetic_energy
    else
      write(*, '(A,F22.6,A)') 'TOTAL ENERGY       ',&
        & kinetic_energy + nuclear_energy + 0.5_dp * coulomb_energy + exchange_energy&
        & + conf_energy, ' Hartree'
    end if
    write(*, '(A)') ' '

    !
    ! NOTE: For DFT
    !
    ! (Sum of Eigenvalues)-0.5*(Coulomb Energy)+(XC Energy)-(\int v_xc^up \rho^up +
    ! \int v_xc^dwn+rho^dwn)=Total Energy
    !
    ! The integrals (\int v_xc^up \rho^up) and (\int v_xc^dwn+rho^dwn) are
    ! computed and printed later ("V_xc integrals from pot.dat, Up/Dwn:")
    !

  end subroutine write_energies


  subroutine write_eigval(max_l, num_alpha, poly_order, eigval)

    implicit none

    integer, intent(in) :: max_l, num_alpha(0:), poly_order(0:)
    real(dp), intent(in) :: eigval(:,0:,:)
    integer :: ii, jj, full, rest

    write(*, '(A)') ' '
    write(*, '(A)') 'EIGENVALUES UP/ Hartree'
    do ii = 0, max_l

      full = num_alpha(ii) * poly_order(ii) / 5
      rest = num_alpha(ii) * poly_order(ii) - full * 5

      write(*, '(A,I3)') 'l=', ii
      do jj = 1, full
        write(*, '(5F15.6)') eigval(1, ii, jj * 5 - 4), eigval(1, ii, jj * 5 - 3),&
            &eigval(1, ii, jj * 5 - 2), eigval(1, ii, jj * 5 - 1), eigval(1, ii, jj * 5)
        write(*, '(A)') ' '
      end do
      if (rest == 4) then
        write(*, '(4F15.6)') eigval(1, ii, 5 * full + 1), eigval(1, ii, 5 * full + 2),&
            &eigval(1, ii, 5 * full + 3), eigval(1, ii, 5 * full + 4)
      else if (rest == 3) then
        write(*, '(3F15.6)') eigval(1, ii, 5 * full + 1), eigval(1, ii, 5 * full + 2),&
            &eigval(1, ii, 5 * full + 3)
      else if (rest == 2) then
        write(*, '(2F15.6)') eigval(1, ii, 5 * full + 1), eigval(1, ii, 5 * full + 2)
      else if (rest == 1) then
        write(*, '(1F15.6)') eigval(1, ii, 5 * full + 1)
      end if
      write(*, '(A)') ' '
    end do
    write(*, '(A)') ' '
    write(*, '(A)') 'EIGENVALUES DWN/ Hartree'
    do ii = 0, max_l

      full = num_alpha(ii) * poly_order(ii) / 5
      rest = num_alpha(ii) * poly_order(ii) - full * 5

      write(*, '(A,I3)') 'l=', ii
      do jj = 1, full
        write(*, '(5F15.6)') eigval(2, ii, jj * 5 - 4), eigval(2, ii, jj * 5 - 3),&
            &eigval(2, ii, jj * 5 - 2), eigval(2, ii, jj * 5 - 1), eigval(2, ii, jj * 5)
        write(*, '(A)') ' '
      end do
      if (rest == 4) then
        write(*, '(4F15.6)') eigval(2, ii, 5 * full + 1), eigval(2, ii, 5 * full + 2),&
            &eigval(2, ii, 5 * full + 3), eigval(2, ii, 5 * full + 4)
      else if (rest == 3) then
        write(*, '(3F15.6)') eigval(2, ii, 5 * full + 1), eigval(2, ii, 5 * full + 2),&
            &eigval(2, ii, 5 * full + 3)
      else if (rest == 2) then
        write(*, '(2F15.6)') eigval(2, ii, 5 * full + 1), eigval(2, ii, 5 * full + 2)
      else if (rest == 1) then
        write(*, '(1F15.6)') eigval(2, ii, 5 * full + 1)
      end if
      write(*, '(A)') ' '
      write(*, '(A)') ' '
    end do
    write(*, '(A)') ' '

  end subroutine write_eigval


  subroutine write_eigvec(max_l, num_alpha, alpha, poly_order, eigval, cof)

    integer, intent(in) :: max_l, num_alpha(0:)
    real(dp), intent(in) :: alpha(0:,:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: eigval(:,0:,:), cof(:,0:,:,:)

    integer :: ispin, ii, jj, kk, ll, nn
    character(4), parameter :: spinname(2) = ["UP  ", "DOWN"]

    write(*,*)
    write(*, "(A)") "EIGENVALUES:"
    write(*,*)
    do ispin = 1, 2
      do jj = 0, max_l
        do ii = 1, num_alpha(jj) * poly_order(jj)
          write(*, "(A,A,A,I1,A,I4,A,F20.8)") "Spin ", trim(spinname(ispin)),&
              & ", l = ", jj, &
              & ", eigenvalue", ii, ":", eigval(ispin, jj, ii)
          write(*, '(A20,1X,A3,1X,A20)') "Alpha", "Pwr", "Coefficient"
          nn = 0
          do kk = 1, num_alpha(jj)
            do ll = 1, poly_order(jj)
              nn = nn + 1
              write(*, "(F20.8,1X,I3,1X,F20.8)") alpha(jj, kk), ll + jj - 1, &
                  &cof(ispin, jj, nn, ii)
            end do
          end do
          write(*,*)
        end do
      end do
    end do

  end subroutine write_eigvec


  subroutine write_moments(max_l, num_alpha, alpha, poly_order, problemsize, cof)

    integer, intent(in) :: max_l, problemsize
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:,:), cof(:,0:,:,:)
    real(dp), allocatable :: moment(:,:,:,:)
    integer, parameter :: exponents(5) = [-3, -1, 0, 1, 2]
    integer :: ii, jj, kk

    allocate(moment(2, 0:max_l, 5, problemsize))
    moment(:,:,:,:) = 0.0_dp

    ! get expectation values <r^-3>, <r^-1>, <1>, <r>, <r^2>
    call moments(moment(:,:, 1, :), max_l, num_alpha, alpha, poly_order, cof, -3)
    call moments(moment(:,:, 2, :), max_l, num_alpha, alpha, poly_order, cof, -1)
    !  call moments(moment(:,: ,3, :), max_l, num_alpha, alpha, poly_order, cof, 0)
    call moments(moment(:,:, 4, :), max_l, num_alpha, alpha, poly_order, cof, 1)
    call moments(moment(:,:, 5, :), max_l, num_alpha, alpha, poly_order, cof, 2)

    write(*, '(A)') 'WAVEFUNCTION EXPECTATION VALUES'
    write(*, '(A)') '-------------------------------'
    write(*, '(A)') ' '

    write(*,*) 'UP Electrons'
    write(*, '(A)') ' '
    do ii = 1, 5
      if (ii /= 3) then
        write(*, '(A,I2,A)') '<r^', exponents(ii), '>, '
        do jj = 0, max_l
          write(*, '(A,I3)') 'l= ', jj
          do kk = 1, num_alpha(jj) * poly_order(jj)
            write(*, '(F12.4)') moment(1, jj, ii, kk)
          end do
          write(*, '(A)') ' '
          write(*, '(A)') ' '
        end do
        write(*, '(A)') ' '
      end if
    end do

    write(*, '(A)') ' '
    write(*, '(A)') ' '

    write(*,*) 'DOWN Electrons'
    write(*, '(A)') ' '
    do ii = 1, 5
      if (ii /= 3) then
        write(*, '(A,I2,A)') '<r^', exponents(ii), '>, '
        do jj = 0, max_l
          write(*, '(A,I3)') 'l= ', jj
          do kk = 1, num_alpha(jj) * poly_order(jj)
            write(*, '(F12.4)') moment(2, jj, ii, kk)
          end do
          write(*, '(A)') ' '
          write(*, '(A)') ' '
        end do
        write(*, '(A)') ' '
      end if
    end do

  end subroutine write_moments

  subroutine write_potentials_file_standard(num_mesh_points, abcissa, weight,&
      &vxc, rho, nuc, p, max_l, num_alpha, poly_order, alpha, problemsize)
    ! write potentials and mesh info to file on standard (internal) integration mesh
    ! in principle one could read in the points from another file to have
    ! other meshes !

    real(dp), intent(in) :: abcissa(:), weight(:), vxc(:,:), p(:,0:,:,:), alpha(0:,:)
    real(dp), intent(in) :: rho(:,:)
    integer, intent(in) :: num_mesh_points, nuc, max_l, num_alpha(0:)
    integer, intent(in) :: poly_order(0:), problemsize
    real(dp), allocatable :: cpot(:), ptot(:,:,:), rhotot(:)
    real(dp) :: ecou, enuc, vxcint(2)
    integer :: ii

    allocate(cpot(num_mesh_points))
    allocate(ptot(0:max_l, problemsize, problemsize))
    allocate(rhotot(num_mesh_points))

    cpot(:) = 0.0_dp
    ptot(:,:,:) = 0.0_dp
    rhotot(:) = 0.0_dp
    ecou = 0.0_dp
    enuc = 0.0_dp
    vxcint(:) = 0.0_dp

    ptot(:,:,:) = p(1, :,:,:) + p(2, :,:,:)
    rhotot(:) = rho(:, 1) + rho(:, 2)

    call cou_pot(ptot, max_l, num_alpha, poly_order, alpha, problemsize, num_mesh_points, abcissa,&
        & cpot)

    open(95, FILE='pot.dat', FORM='formatted', STATUS='unknown')
    write(95, '(A)') '# 1st line: number of mesh points'
    write(95, '(A)') '# abcissa weight nuclear coulomb dft-vxc_up dft-vxc_down'
    write(95, '(I0)') num_mesh_points

    do ii = 1, num_mesh_points
      write(95, '(6ES21.12E3)') abcissa(ii), weight(ii), real(- nuc, dp) / abcissa(ii), cpot(ii),&
          & vxc(ii, 1), vxc(ii, 2)
    end do
    close(95)

    do ii = 1, num_mesh_points
      ecou = ecou + weight(ii) * rhotot(ii) * cpot(ii) * abcissa(ii)**2
      enuc = enuc - weight(ii) * rhotot(ii) * real(nuc, dp) * abcissa(ii)
      vxcint(1) = vxcint(1) + weight(ii) * rho(ii, 1) * vxc(ii, 1) * abcissa(ii)**2
      vxcint(2) = vxcint(2) + weight(ii) * rho(ii, 2) * vxc(ii, 2) * abcissa(ii)**2
    end do

    write(*, '(A,F18.6)') 'Nuc. attr. energy from potential in pot.dat: ', enuc
    write(*, '(A,F18.6)') 'Coulomb    energy from potential in pot.dat: ', 0.5_dp * ecou
    write(*, '(A,2F18.6)') 'V_xc integrals from pot.dat, Up/Dwn: ', vxcint(1), vxcint(2)

  end subroutine write_potentials_file_standard


  !> Writes potentials and mesh info to file on standard (internal) integration mesh;
  !! in principle one could read in the points from another file to have other meshes!
  subroutine write_densities_file_standard(num_mesh_points, abcissa, weight, rho, drho, ddrho)

    real(dp), intent(in) :: abcissa(:), weight(:)
    real(dp), intent(in) :: rho(:,:), drho(:,:), ddrho(:,:)
    integer, intent(in) :: num_mesh_points
    real(dp) :: enumber, zeta, r_seitz
    integer :: ii

    open(95, file='dens.dat', form='formatted', status='unknown')
    write(95, '(A)') '# 1st line: number of mesh points'
    write(95, '(A)') '# rho and r_seitz are calculated from total density'
    write(95, '(A)') '# zeta and r_seitz only correct of rho > 1d-12'
    write(95, '(A)') ''
    write(95, '(A)') '# abcissa weight rho drho ddrho zeta r_seitz'
    write(95, '(I0)') num_mesh_points

    enumber = 0.0_dp

    ! note division of total density by 4*pi in calculation of r_seitz
    ! commonly r_seitz=((4*pi*rho)/3)**(-1/3) but our rho is from the
    ! radial part only and the angular part must be taken into account
    ! explicitely; during integration this happens implicitely, see enumber

    do ii = 1, num_mesh_points

      if ((rho(ii, 1) + rho(ii, 2)) > 1.0d-12) then
        zeta = (rho(ii, 1) - rho(ii, 2)) / (rho(ii, 1) + rho(ii, 2))
        r_seitz = (4.0_dp * pi / 3.0_dp&
            & * ((rho(ii, 1) + rho(ii, 2)) / 4.0_dp / pi))**(-1.0_dp / 3.0_dp)
      else
        zeta = 0.0_dp
        r_seitz = 0.0_dp
      end if

      write(95, '(7ES21.12E3)') abcissa(ii), weight(ii), rho(ii, 1) + rho(ii, 2),&
          & drho(ii, 1) + drho(ii, 2), ddrho(ii, 1) + ddrho(ii, 2), zeta, r_seitz
      enumber = enumber + weight(ii) * (rho(ii, 1) + rho(ii, 2)) * abcissa(ii)**2
    end do

    close(95)

    write(*, '(A,F18.6)') 'Total number of electrons from dens.dat : ', enumber

  end subroutine write_densities_file_standard


  !> Writes potentials and mesh info to file on standard (internal) integration mesh;
  !! in principle one could read in the points from another file to have other meshes!
  subroutine write_waves_file_standard(num_mesh_points, abcissa, weight, alpha, num_alpha,&
      & poly_order, max_l, problemsize, occ, qnvalorbs, cof)

    real(dp), intent(in) :: abcissa(:), weight(:), alpha(0:,:)
    real(dp), intent(in) :: occ(:,0:,:)
    integer, intent(in) :: num_mesh_points, num_alpha(0:), poly_order(0:), max_l
    integer, intent(in) :: problemsize
    integer, intent(in) :: qnvalorbs(:,0:)
    real(dp), intent(inout) :: cof(:,0:,:,:)

    real(dp), allocatable :: coftot(:)
    real(dp) :: xx
    integer :: ii, jj, kk, ll, mm, ispin, imax
    character(len=20) :: fname
    real(dp), allocatable :: wavedata(:)

    allocate(wavedata(num_mesh_points))
    allocate(coftot(problemsize))

    do jj = 0, max_l
      mm = 0
      do kk = 1, num_alpha(jj)
        do ll = 1, poly_order(jj)
          mm = mm + 1
          if (mm < qnvalorbs(1, jj) .or. mm > qnvalorbs(2, jj)) then
            cycle
          end if
          do ispin = 1, 2
            if (ispin == 1) then
              write(fname, "(A,I2.2,A,A)") "wave_", mm + jj, orbnames(jj), "_up.dat"
            else
              write(fname, "(A,I2.2,A,A)") "wave_", mm + jj, orbnames(jj), "_dn.dat"
            end if
            open(95, file=fname, form='formatted', status='unknown')
            write(95, '(A)') '# 1st line: number of mesh points'
            write(95, '(A)') '# abcissa weight wavefunction wavefunction_1st wavefunction_2nd'
            write(95, '(I0)') num_mesh_points
            write(95, '(A,I3,A,I3,A,F8.4)') '# Principal QN= ', mm, ' , l= ', jj,&
                & ' , Occupation= ', occ(1, jj, mm) + occ(2, jj, mm)

            coftot(:) = cof(ispin, jj, :, mm)

            do ii = 1, num_mesh_points
              xx = abcissa(ii)
              wavedata(ii) = wavefunction(coftot, alpha, num_alpha, poly_order, jj, xx)
            end do
            imax = maxloc(abs(abcissa * wavedata), dim=1)
            if (wavedata(imax) < 0.0_dp) then
              coftot = -coftot
              cof(1, jj, :, mm) = coftot
              write(*, "(A,I3,A,I3)") "Changing wavefunction sign: n =", mm + jj, ", l =", jj
            end if

            do ii = 1, num_mesh_points
              xx = abcissa(ii)
              write(95, '(5ES21.12E3)') xx, weight(ii),&
                  & wavefunction(coftot, alpha, num_alpha, poly_order, jj, xx),&
                  & wavefunction_1st(coftot, alpha, num_alpha, poly_order, jj, xx),&
                  & wavefunction_2nd(coftot, alpha, num_alpha, poly_order, jj, xx)
            end do
            close(95)
          end do
        end do
      end do
    end do

  end subroutine write_waves_file_standard


  !> Writes average local potential to disk.
  subroutine writeAveragePotential(abcissa, avgPot)

    !> Numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !! Average local, effective potential
    real(dp), intent(in) :: avgPot(:,:)

    !! Number of numerical integration points
    integer :: num_mesh_points

    !! Iterates over radial grid points
    integer :: iRad

    !! File name and identifier
    character(mc) :: fname
    integer :: fp

    fname = "avgpot.dat"
    num_mesh_points = size(abcissa)

    open(newunit=fp, file=fname, status="replace", action="write")

    write(fp, "(A)") "# 1st line: number of mesh points"
    write(fp, "(A)") "# abcissa avgpot_up avgpot_down"
    write(fp, "(I0)") num_mesh_points

    do iRad = 1, size(avgPot, dim=1)
      write(fp, "(3ES21.12E3)") abcissa(iRad), avgPot(iRad, 1), avgPot(iRad, 2)
    end do

    close(fp)

  end subroutine writeAveragePotential


  subroutine cusp_values(max_l, cof, p, alpha, num_alpha, poly_order)

    integer, intent(in) :: max_l, num_alpha(0:), poly_order(0:)
    real(dp), intent(in) :: cof(:,0:,:,:), alpha(0:,:), p(:,0:,:,:)
    integer :: ii

    write(*, '(A)') 'Cusp Values '
    write(*, '(A)') '------------'

    ii = 0

    write(*, '(A,F14.6)') '1s, UP  ',&
        & - wavefunction_1st(cof(1, ii, :, 1), alpha, num_alpha, poly_order, ii, 0.0_dp)&
        & / wavefunction(cof(1, ii, :, 1), alpha, num_alpha, poly_order, ii, 0.0_dp)
    write(*, '(A,F14.6)') '1s, DWN ',&
        & - wavefunction_1st(cof(2, ii, :, 1), alpha, num_alpha, poly_order, ii, 0.0_dp)&
        & / wavefunction(cof(2, ii, :, 1), alpha, num_alpha, poly_order, ii, 0.0_dp)

    write(*, '(A,F14.6)') 'Total density UP ',&
        & - density_at_point_1st(p(1, :, :, :), max_l, num_alpha, poly_order, alpha, 0.0_dp)&
        & / density_at_point(p(1, :, :, :), max_l, num_alpha, poly_order, alpha, 0.0_dp) / 2.0_dp
    write(*, '(A,F14.6)') 'Total density DWN ',&
        & - density_at_point_1st(p(2, :, :, :), max_l, num_alpha, poly_order, alpha, 0.0_dp)&
        & / density_at_point(p(2, :, :, :), max_l, num_alpha, poly_order, alpha, 0.0_dp) / 2.0_dp

    write(*, '(A)') ' '

  end subroutine cusp_values


  subroutine write_energies_tagged(ekin, enuc, ecoul, exc, econf, etot, zora, eigvals, occ)

    real(dp), intent(in) :: ekin, enuc, ecoul, exc, etot, econf
    logical, intent(in) :: zora
    real(dp), intent(in) :: eigvals(:,0:,:), occ(:,0:,:)

    integer :: fp
    type(TTaggedwriter) :: twriter

    call TTaggedwriter_init(twriter)

    open(newunit=fp, file="energies.tag", status="replace", action="write")
    call writetag(twriter, fp, "zora", zora)
    call writetag(twriter, fp, "kinetic_energy", ekin)
    call writetag(twriter, fp, "nuclear_energy", enuc)
    call writetag(twriter, fp, "coulomb_energy", 0.5_dp * ecoul)
    call writetag(twriter, fp, "xc_energy", exc)
    call writetag(twriter, fp, "confinement_energy", econf)
    call writetag(twriter, fp, "total_energy", etot)
    !! Transposing eigenvalues to appear in a more convinient order
    call writetag(twriter, fp, "eigenlevels_up", transpose(eigvals(1, :, :)))
    call writetag(twriter, fp, "eigenlevels_dn", transpose(eigvals(2, :, :)))
    call writetag(twriter, fp, "occupations_up", transpose(occ(1, :, :)))
    call writetag(twriter, fp, "occupations_dn", transpose(occ(2, :, :)))
    close(fp)

  end subroutine write_energies_tagged


  subroutine write_wave_coeffs_file(max_l, num_alpha, poly_order, cof, alpha, occ, qnvalorbs)

    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:), poly_order(0:)
    real(dp), intent(in) :: cof(:,0:,:,:), alpha(0:,:), occ(:,0:,:)
    integer, intent(in) :: qnvalorbs(:,0:)

    integer :: fp, ii, ll
    type(TTaggedwriter) :: twriter
    character(len=20) :: fname
    real(dp), allocatable :: coeffs(:,:)

    call TTaggedwriter_init(twriter)

    do ll = 0, max_l
      allocate(coeffs(poly_order(ll), num_alpha(ll)))
      do ii = 1, num_alpha(ll) * poly_order(ll)
        if (ii < qnvalorbs(1, ll) .or. ii > qnvalorbs(2, ll)) then
          cycle
        end if
        write(fname, "(A,I2.2,A,A)") "coeffs_", ii + ll, orbnames(ll), ".tag"
        open(newunit=fp, file=fname, status="replace", action="write")
        call writetag(twriter, fp, "exponents", alpha(ll, :num_alpha(ll)))
        call convcoeffs(cof(1, ll, :, ii), alpha(ll, :num_alpha(ll)), ll, coeffs)
        call writetag(twriter, fp, "coefficients", coeffs)
        call writetag(twriter, fp, "occupation", sum(occ(:, ll, ii)))
        close(fp)
      end do
      deallocate(coeffs)
    end do

  contains

    subroutine convcoeffs(cof, alpha, angmom, normcoeffs)

      real(dp), intent(in) :: cof(:), alpha(:)
      integer, intent(in) :: angmom
      real(dp), intent(out) :: normcoeffs(:,:)

      integer :: npow, nalpha, ialpha, ipow
      real(dp) :: aa, normfac

      npow = size(normcoeffs, dim=1)
      nalpha = size(normcoeffs, dim=2)
      normcoeffs = reshape(cof, [npow, nalpha])
      do ialpha = 1, nalpha
        aa = alpha(ialpha)
        do ipow = 1, npow
          normfac = (2.0_dp * aa)**(ipow + angmom) * sqrt(2.0_dp * aa)&
              & / sqrt(fak(2 * (ipow + angmom)))
          normcoeffs(ipow, ialpha) = normfac * normcoeffs(ipow, ialpha)
        end do
      end do

    end subroutine convcoeffs

  end subroutine write_wave_coeffs_file

end module output
