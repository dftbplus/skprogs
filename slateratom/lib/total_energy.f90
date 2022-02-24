!> Module that provides routines for calculating the total energy of a system.
module totalenergy

  use common_accuracy, only : dp
  use dft, only : dft_exc_energy, dft_vxc_energy

  implicit none
  private

  public :: total_energy, zora_total_energy


contains

  !> Calculates total energy for non-ZORA calculations.
  subroutine total_energy(tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha, poly_order,&
      & problemsize, xcnr, num_mesh_points, weight, abcissa, rho, exc, camAlpha, camBeta, kinetic,&
      & nuclear, coulomb, exchange, dft_xc_energy, confinement, etot)

    !> kinetic supervector
    real(dp), intent(in) :: tt(0:,:,:)

    !> nucleus-electron supervector
    real(dp), intent(in) :: uu(0:,:,:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> confinement supervector
    real(dp), intent(in) :: vconf(0:,:,:)

    !> coulomb supermatrix
    real(dp), intent(in) :: jj(0:,:,:,0:,:,:)

    !> (hf) exchange supermatrix
    real(dp), intent(in) :: kk(0:,:,:,0:,:,:)

    !> (hf) exchange supermatrix (long-range, range-separated version)
    real(dp), intent(in) :: kk_lr(0:,:,:,0:,:,:)

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> exc energy density on grid
    real(dp), intent(in) :: exc(:)

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> resulting energy portions
    real(dp), intent(out) :: kinetic, nuclear, coulomb, exchange, dft_xc_energy, confinement, etot

    !! exchange contribution of long-range Hartree-Fock exchange matrix
    real(dp) :: exchange_lr

    !! density matrix supervector (spins summed up)
    real(dp), allocatable :: p_total(:,:,:)

    !! true, if a (semi-)local functional is requested
    logical :: tGgaLda

    !! true, if a (long-range corrected) range-separated hybrid functional is requested
    logical :: tLC

    !! true, if a CAM functional is requested
    logical :: tCam

    !! true, if a global hybrid functional is requested
    logical :: tGlobalHybrid

    !! auxiliary variable, i.e. (nuclear + kinetic + confinement)
    real(dp) :: dummy1

    tGgaLda = ((xcnr > 0) .and. (xcnr <= 4))
    tLC = ((xcnr == 5) .or. (xcnr == 6))
    tCam = ((xcnr == 9) .or. (xcnr == 10))
    tGlobalHybrid = ((xcnr == 7) .or. (xcnr == 8))

    etot = 0.0_dp
    kinetic = 0.0_dp
    nuclear = 0.0_dp
    dft_xc_energy = 0.0_dp
    confinement = 0.0_dp
    coulomb = 0.0_dp
    exchange = 0.0_dp
    exchange_lr = 0.0_dp

    ! build total density matrix
    p_total = pp(1, :,:,:) + pp(2, :,:,:)

    ! get total energy
    call core_hamiltonian_energies(tt, uu, vconf, p_total, max_l, num_alpha, poly_order, nuc,&
        & kinetic, nuclear, confinement)

    dummy1 = nuclear + kinetic + confinement

    ! get Coulomb and HF exchange contributions
    coulomb = coulomb_energy(jj, p_total, max_l, num_alpha, poly_order)
    if (xcnr == 0) then
      exchange = hf_ex_energy(kk, pp, max_l, num_alpha, poly_order)
    elseif (tLC) then
      exchange = hf_ex_energy(kk_lr, pp, max_l, num_alpha, poly_order)
    elseif (tGlobalHybrid) then
      exchange = hf_ex_energy(kk, pp, max_l, num_alpha, poly_order)
      if (xcnr == 7) then
        ! PBE0 requires 1/4 Hartree-Fock exchange
        exchange = 1.0_dp / 4.0_dp * exchange
      elseif (xcnr == 8) then
        ! B3LYP requires 0.20 * HF exchange
        exchange = 0.20_dp * exchange
      end if
    elseif (tCam) then
      exchange = hf_ex_energy(kk, pp, max_l, num_alpha, poly_order)
      exchange_lr = hf_ex_energy(kk_lr, pp, max_l, num_alpha, poly_order)
      if (xcnr == 9) then
        ! CAMY-B3LYP parameters a=0.20, b=0.72, c=0.81 (libXC defaults)
        ! exchange = camAlpha * 0.20_dp * exchange
        exchange = camAlpha * exchange
        exchange_lr = camBeta * exchange_lr
        exchange = exchange + exchange_lr
      elseif (xcnr == 10) then
        ! CAMY-PBE0
        exchange = camAlpha * exchange
        exchange_lr = camBeta * exchange_lr
        exchange = exchange + exchange_lr
      end if
    end if

    ! pure HF:
    ! coulomb = P^Tot J P^Tot. The Coulomb energy is half of that.
    ! exchange = -P^up K P^up - -P^dn K P^dn. The exchange energy is half of that.

    ! DFT exchange-correlation energy
    if (xcnr > 0) then
      call dft_exc_energy(num_mesh_points, rho, exc, weight, abcissa, dft_xc_energy)
    end if

    ! pure HF
    if (xcnr == 0) then
      etot = dummy1 + 0.5_dp * coulomb + 0.5_dp * exchange
    ! (semi-)local functionals
    elseif ((xcnr > 0) .and. (xcnr <= 4)) then
      etot = dummy1 + 0.5_dp * coulomb + dft_xc_energy
    ! global hybrids, LC, CAM functionals
    else
      etot = dummy1 + 0.5_dp * coulomb + dft_xc_energy + 0.5_dp * exchange
    end if

  end subroutine total_energy


  !> Calculates total energy for ZORA, note that the total energy is not well defined here.
  subroutine zora_total_energy(tt, uu, nuc, vconf, jj, kk, kk_lr, pp, max_l, num_alpha,&
      & poly_order, problemsize, xcnr, num_mesh_points, weight, abcissa, rho, exc, vxc,&
      & eigval_scaled, occ, camAlpha, camBeta, kinetic, nuclear, coulomb, exchange, xc_energy,&
      & confinement, etot)

    !> kinetic supervector
    real(dp), intent(in) :: tt(0:,:,:)

    !> nucleus-electron supervector
    real(dp), intent(in) :: uu(0:,:,:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> confinement supervector
    real(dp), intent(in) :: vconf(0:,:,:)

    !> coulomb supermatrix
    real(dp), intent(in) :: jj(0:,:,:,0:,:,:)

    !> (hf) exchange supermatrix
    real(dp), intent(in) :: kk(0:,:,:,0:,:,:)

    !> (hf) exchange supermatrix (long-range, range-separated version)
    real(dp), intent(in) :: kk_lr(0:,:,:,0:,:,:)

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> number of numerical integration points
    integer, intent(in) :: num_mesh_points

    !> numerical integration weights
    real(dp), intent(in) :: weight(:)

    !> numerical integration abcissas
    real(dp), intent(in) :: abcissa(:)

    !> density on grid
    real(dp), intent(in) :: rho(:,:)

    !> exc energy density on grid
    real(dp), intent(in) :: exc(:)

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> zora scaled eigenvalues
    real(dp), intent(in) :: eigval_scaled(:,0:,:)

    !> occupation numbers
    real(dp), intent(in) :: occ(:,0:,:)

    !> CAM alpha parameter
    real(dp), intent(in) :: camAlpha

    !> CAM beta parameter
    real(dp), intent(in) :: camBeta

    !> resulting energy portions
    real(dp), intent(out) :: kinetic, nuclear, coulomb, exchange, xc_energy, confinement, etot

    !! exchange contribution of long-range Hartree-Fock exchange matrix
    real(dp) :: exchange_lr

    !! density matrix supervector (spins summed up)
    real(dp), allocatable :: p_total(:,:,:)

    !! true, if a (semi-)local functional is requested
    logical :: tGgaLda

    !! true, if a (long-range corrected) range-separated hybrid functional is requested
    logical :: tLC

    !! true, if a CAM functional is requested
    logical :: tCam

    !! true, if a global hybrid functional is requested
    logical :: tGlobalHybrid

    !> auxiliary variables
    integer :: mm, nn, oo
    real(dp) :: xc_pot, dummy(2), eigsum

    tGgaLda = ((xcnr > 0) .and. (xcnr <= 4))
    tLC = ((xcnr == 5) .or. (xcnr == 6))
    tCam = (xcnr == 9)
    tGlobalHybrid = ((xcnr == 7) .or. (xcnr == 8))

    etot = 0.0_dp
    kinetic = 0.0_dp
    nuclear = 0.0_dp
    confinement = 0.0_dp
    coulomb = 0.0_dp
    exchange = 0.0_dp
    exchange_lr = 0.0_dp
    xc_energy = 0.0_dp
    eigsum = 0.0_dp

    ! build total density matrix
    p_total = pp(1, :,:,:) + pp(2, :,:,:)

    ! get total energy
    call core_hamiltonian_energies(tt, uu, vconf, p_total, max_l, num_alpha, poly_order, nuc,&
        & kinetic, nuclear, confinement)

    ! sum of occupied eigenvalues
    do mm = 1, 2
      do nn = 0, max_l
        do oo = 1, problemsize
          eigsum = eigsum + eigval_scaled(mm, nn, oo) * occ(mm, nn, oo)
        end do
      end do
    end do

    kinetic = eigsum

    ! get Coulomb and HF exchange contributions
    coulomb = coulomb_energy(jj, p_total, max_l, num_alpha, poly_order)
    if (tLC) then
      exchange = hf_ex_energy(kk, pp, max_l, num_alpha, poly_order)
    end if

    call dft_exc_energy(num_mesh_points, rho, exc, weight, abcissa, xc_energy)
    call dft_vxc_energy(num_mesh_points, rho, vxc, weight, abcissa, dummy)

    xc_pot = dummy(1) + dummy(2)

    ! pure Hartree-Fock
    if (xcnr == 0) then
      etot = eigsum - 0.5_dp * coulomb - 0.5_dp * exchange
    ! (semi-)local functionals
    elseif (tGgaLda) then
      etot = eigsum - 0.5_dp * coulomb + xc_energy - xc_pot
    ! range-separated (long-range corrected) hybrid functionals
    elseif (tLC) then
      etot = eigsum - 0.5_dp * coulomb - 0.5_dp * exchange + xc_energy - xc_pot
    elseif (tGlobalHybrid) then
      etot = eigsum - 0.5_dp * coulomb - 0.5_dp * exchange + xc_energy - xc_pot
    elseif (tCam) then
      etot = eigsum - 0.5_dp * coulomb - 0.5_dp * exchange + xc_energy - xc_pot
    end if

  end subroutine zora_total_energy


  !> Calculates Coulomb contribution to total energy by multiplying the jj supermatrix with the
  !! density matrix supervector twice.
  pure function coulomb_energy(jj, p_total, max_l, num_alpha, poly_order) result(coulomb)

    !> coulomb supermatrix
    real(dp), intent(in) :: jj(0:,:,:,0:,:,:)

    !> density matrix supervector (spins summed up)
    real(dp), intent(in) :: p_total(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> resulting Coulomb energy contribution
    real(dp) :: coulomb

    !! auxiliary variables
    integer :: ii, jjj, kkk, ll, mm, nn, oo, ppp, qq, rr, ss, tt, uu, vv

    coulomb = 0.0_dp
    do ii = 0, max_l
      ss = 0
      do jjj = 1, num_alpha(ii)
        do kkk = 1, poly_order(ii)
          ss = ss + 1
          tt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              tt = tt + 1
              do nn = 0, max_l
                uu = 0
                do oo = 1, num_alpha(nn)
                  do ppp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        coulomb = coulomb + p_total(ii, ss, tt) * jj(ii, ss, tt, nn, uu, vv)&
                            & * p_total(nn, uu, vv)

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end function coulomb_energy


  !> Calculates Hartee-Fock exchange contribution to total energy by the kk supermatrices with the
  !! density matrix supervector twice.
  pure function hf_ex_energy(kk, pp, max_l, num_alpha, poly_order) result(exchange)

    !> (hf) exchange supermatrix
    real(dp), intent(in) :: kk(0:,:,:,0:,:,:)

    !> density matrix supervector
    real(dp), intent(in) :: pp(:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> resulting exchange energy contribution
    real(dp) :: exchange

    !! auxiliary variables
    integer :: ii, jj, kkk, ll, mm, nn, oo, ppp, qq, rr, ss, tt, uu, vv

    exchange = 0.0_dp

    do ii = 0, max_l
      ss = 0
      do jj = 1, num_alpha(ii)
        do kkk = 1, poly_order(ii)
          ss = ss + 1
          tt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              tt = tt + 1
              do nn = 0, max_l
                uu = 0
                do oo = 1, num_alpha(nn)
                  do ppp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        exchange = exchange - kk(ii, ss, tt, nn, uu, vv)&
                            & * (pp(1, ii, ss, tt) * pp(1, nn, uu, vv)&
                            & + pp(2, ii, ss, tt) * pp(2, nn, uu, vv))

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end function hf_ex_energy


  !> Core Hamiltonian contributions to total energy by multiplying the tt, uu, vconf supervectors
  !! with the density matrix supervector once.
  pure subroutine core_hamiltonian_energies(tt, uu, vconf, p_total, max_l, num_alpha, poly_order,&
      & nuc, kinetic, nuclear, confinement)

    !> kinetic supervector
    real(dp), intent(in) :: tt(0:,:,:)

    !> nucleus-electron supervector
    real(dp), intent(in) :: uu(0:,:,:)

    !> confinement supervector
    real(dp), intent(in) :: vconf(0:,:,:)

    !> density matrix supervector (spins summed up)
    real(dp), intent(in) :: p_total(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> nuclear charge, i.e. atomic number
    integer, intent(in) :: nuc

    !> resulting energy portions
    real(dp), intent(out) :: kinetic, nuclear, confinement

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, ss, ttt

    kinetic = 0.0_dp
    nuclear = 0.0_dp
    confinement = 0.0_dp
    do ii = 0, max_l
      ss = 0
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          ss = ss + 1
          ttt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              ttt = ttt + 1
              kinetic = kinetic + tt(ii, ss, ttt) * p_total(ii, ss, ttt)
              nuclear = nuclear - real(nuc, dp) * uu(ii, ss, ttt) * p_total(ii, ss, ttt)
              confinement = confinement + vconf(ii, ss, ttt) * p_total(ii, ss, ttt)
            end do
          end do
        end do
      end do
    end do

  end subroutine core_hamiltonian_energies

end module totalenergy
