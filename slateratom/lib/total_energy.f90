!> Module that provides routines for calculating the total energy of a system.
module totalenergy

  use common_accuracy, only : dp
  use dft, only : dft_exc_energy, dft_vxc_energy

  implicit none
  private

  public :: total_energy, zora_total_energy


contains

  !> Calculates total energy for non-ZORA calculations.
  pure subroutine total_energy(tt, uu, nuc, vconf, jj, kk, pp, max_l, num_alpha, poly_order,&
      & problemsize, xcnr, num_mesh_points, weight, abcissa, rho, exc, kinetic, nuclear, coulomb,&
      & exchange, confinement, etot)

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

    !> resulting energy portions
    real(dp), intent(out) :: kinetic, nuclear, coulomb, exchange, confinement, etot

    !> density matrix supervector (spins summed up)
    real(dp), allocatable :: p_total(:,:,:)

    !> auxiliary variables
    real(dp) :: dummy1, dummy2

    etot = 0.0_dp
    kinetic = 0.0_dp
    nuclear = 0.0_dp
    confinement = 0.0_dp
    coulomb = 0.0_dp
    exchange = 0.0_dp
    dummy1 = 0.0_dp
    dummy2 = 0.0_dp

    ! Build total density matrix
    p_total = pp(1, :,:,:) + pp(2, :,:,:)

    ! get total energy

    call core_hamiltonian_energies(tt, uu, vconf, p_total, max_l, num_alpha, poly_order, nuc,&
        & kinetic, nuclear, confinement)

    dummy1 = nuclear + kinetic + confinement

    call coulomb_hf_ex_energy(jj, kk, p_total, max_l, num_alpha, poly_order, xcnr, coulomb,&
        & exchange)

    if (xcnr > 0) then
      exchange = 0.0_dp
      call dft_exc_energy(num_mesh_points, rho, exc, weight, abcissa, exchange)
    end if

    ! make sure total energy breakdown agrees with total energy

    if (xcnr == 0) then
      etot = dummy1 + 0.5_dp * coulomb + 0.5_dp * exchange
    else
      etot = dummy1 + 0.5_dp * coulomb + exchange
    end if

  end subroutine total_energy


  !> Calculates total energy for ZORA, note that the total energy is not well defined here.
  pure subroutine zora_total_energy(tt, uu, nuc, vconf, jj, kk, pp, max_l, num_alpha, poly_order,&
      & problemsize, xcnr, num_mesh_points, weight, abcissa, rho, exc, vxc, eigval_scaled, occ,&
      & kinetic, nuclear, coulomb, exchange, confinement, etot)

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

    !> resulting energy portions
    real(dp), intent(out) :: kinetic, nuclear, coulomb, exchange, confinement, etot

    !> density matrix supervector (spins summed up)
    real(dp), allocatable :: p_total(:,:,:)

    !> auxiliary variables
    integer :: mm, nn, oo
    real(dp) :: dummy1, dummy2, dummy3(2), eigsum

    etot = 0.0_dp
    kinetic = 0.0_dp
    nuclear = 0.0_dp
    confinement = 0.0_dp
    coulomb = 0.0_dp
    exchange = 0.0_dp
    dummy1 = 0.0_dp
    dummy2 = 0.0_dp
    eigsum = 0.0_dp

    ! Build total density matrix
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

    call coulomb_hf_ex_energy(jj, kk, p_total, max_l, num_alpha, poly_order, xcnr, coulomb,&
        & exchange)

    exchange = 0.0_dp
    call dft_exc_energy(num_mesh_points, rho, exc, weight, abcissa, exchange)

    call dft_vxc_energy(num_mesh_points, rho, vxc, weight, abcissa, dummy3)

    dummy2 = dummy3(1) + dummy3(2)

    etot = eigsum - 0.5_dp * coulomb + exchange - dummy2

  end subroutine zora_total_energy


  !> Calculates Hartee-Fock exchange and Coulomb contributions to total energy by multiplying jj and
  !! kk supermatrixes with the density matrix supervector twice.
  pure subroutine coulomb_hf_ex_energy(jj, kk, p_total, max_l, num_alpha, poly_order, xcnr,&
      & coulomb, exchange)

    !> coulomb supermatrix
    real(dp), intent(in) :: jj(0:,:,:,0:,:,:)

    !> (hf) exchange supermatrix
    real(dp), intent(in) :: kk(0:,:,:,0:,:,:)

    !> density matrix supervector (spins summed up)
    real(dp), intent(in) :: p_total(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> identifier of exchange-correlation type
    integer, intent(in) :: xcnr

    !> resulting energy portions
    real(dp), intent(out) :: coulomb, exchange

    !> auxiliary variables
    integer :: ii, jjj, kkk, ll, mm, nn, oo, pp, qq, rr, ss, tt, uu, vv

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
                  do pp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        coulomb = coulomb + p_total(ii, ss, tt) * jj(ii, ss, tt, nn, uu, vv)&
                            & * p_total(nn, uu, vv)

                        if (xcnr == 0) then
                          exchange = exchange - 0.5_dp * p_total(ii, ss, tt)&
                              & * kk(ii, ss, tt, nn, uu, vv) * p_total(nn, uu, vv)
                        end if

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

  end subroutine coulomb_hf_ex_energy


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
