!> Module that provides the functionality to build-up the HF-Hamiltonian.
module hamiltonian

  use common_accuracy, only : dp
  use dft, only : dft_exc_matrixelement
  use broyden, only : mixing_driver
  use zora_routines, only : zora_t_correction

  implicit none
  private

  public :: build_fock, build_coulomb_matrix
  public :: build_hf_ex_matrix, build_dft_exc_matrix


contains

  !> Main driver routine for Fock matrix build-up. Also calls mixer with potential matrix.
  subroutine build_fock(iScf, tt, uu, nuc, vconf, jj, kk, pp, max_l, num_alpha, poly_order,&
      & problemsize, xcnr, num_mesh_points, weight, abcissa, vxc, alpha, pot_old, pot_new,&
      & tZora, tBroyden, mixing_factor, ff)

    !> current SCF step
    integer, intent(in) :: iScf

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

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> old potential
    real(dp), intent(in) :: pot_old(:,0:,:,:)

    !> new potential
    real(dp), intent(out) :: pot_new(:,0:,:,:)

    !> true, if zero-order regular approximation for relativistic effects is desired
    logical, intent(in) :: tZora

    !> true, if Broyden mixing is desired, otherwise simple mixing is applied
    logical, intent(in) :: tBroyden

    !> mixing factor
    real(dp), intent(in) :: mixing_factor

    !> fock matrix supervector
    real(dp), intent(out) :: ff(:,0:,:,:)

    !> auxiliary matrices
    real(dp), allocatable :: j_matrix(:,:,:), k_matrix(:,:,:,:), p_total(:,:,:), t_zora(:,:,:,:)

    !> auxiliary matrix for hybrids
    real(dp), allocatable :: k_matrix2(:,:,:,:)

    !> auxiliary variables
    integer :: ii, jjj, kkk, ll, mm, ss, ttt, iMix

    ff(:,:,:,:) = 0.0_dp

    allocate(j_matrix(0:max_l, problemsize, problemsize))
    allocate(k_matrix(2, 0:max_l, problemsize, problemsize))

    ! additional k_matrix for hybrids
    allocate(k_matrix2(2, 0:max_l, problemsize, problemsize))

    allocate(t_zora(2, 0:max_l, problemsize, problemsize))
    t_zora(:,:,:,:) = 0.0_dp

    ! form total densitymatrix supervector from spin channels
    allocate(p_total(0:max_l, problemsize, problemsize))
    p_total(:,:,:) = pp(1, :,:,:) + pp(2, :,:,:)

    ! build coulomb matrices
    call build_coulomb_matrix(jj, p_total, max_l, num_alpha, poly_order, j_matrix)

    ! build exchange(-correlation) potential matrices:

    ! pure Hartree-Fock
    if (xcnr == 0) then
      call build_hf_ex_matrix(kk, pp, max_l, num_alpha, poly_order, k_matrix)
    end if

    ! pure DFT
    if ((xcnr > 0) .and. (xcnr <= 4)) then
      call build_dft_exc_matrix(max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa,&
          & weight, vxc, k_matrix)
    end if

    ! HF - DFT hybrid
    if (xcnr >= 5) then
       call build_hf_ex_matrix(kk, pp, max_l, num_alpha, poly_order, k_matrix)
       call build_dft_exc_matrix(max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa,&
           & weight, vxc, k_matrix2)
       k_matrix(1, :,:,:) = k_matrix(1, :,:,:) + k_matrix2(1, :,:,:)
       k_matrix(2, :,:,:) = k_matrix(2, :,:,:) + k_matrix2(2, :,:,:)
    end if

    ! build mixer input
    pot_new(1, :,:,:) = - real(nuc, dp) * uu + j_matrix - k_matrix(1, :,:,:)
    pot_new(2, :,:,:) = - real(nuc, dp) * uu + j_matrix - k_matrix(2, :,:,:)

    ! mixer
    iMix = int(iScf / 40)
    call mixing_driver(pot_old, pot_new, max_l, num_alpha, poly_order, problemsize,&
        & iScf - iMix * 40, tBroyden, mixing_factor)

    ! Not sure: before or after mixer .... ? Potential .ne. Matrix elements
    ! Should be irrelevant once self-consistency is reached
    if (tZora .and. (iScf /= 0)) then
      call zora_t_correction(1, t_zora, max_l, num_alpha, alpha, poly_order, num_mesh_points,&
          & weight, abcissa, vxc, nuc, pp, problemsize)
    end if

    ! finally build Fock matrix
    do ii = 0, max_l
      ss = 0
      do jjj = 1, num_alpha(ii)
        do kkk = 1, poly_order(ii)
          ss = ss + 1
          ttt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              ttt = ttt + 1

              ff(1, ii, ss, ttt) = tt(ii, ss, ttt) + pot_new(1, ii, ss, ttt) + vconf(ii, ss, ttt)
              ff(2, ii, ss, ttt) = tt(ii, ss, ttt) + pot_new(2, ii, ss, ttt) + vconf(ii, ss, ttt)

              if (tZora) then
                ff(1, ii, ss, ttt) = ff(1, ii, ss, ttt) + t_zora(1, ii, ss, ttt)
                ff(2, ii, ss, ttt) = ff(2, ii, ss, ttt) + t_zora(2, ii, ss, ttt)
              end if

            end do
          end do
        end do
      end do
    end do

  end subroutine build_fock


  !> Builds Coulomb matrix to be added to the Fock matrix from Coulomb supermatrix by multiplying
  !! with density matrix supervector.
  subroutine build_coulomb_matrix(jj, pp, max_l, num_alpha, poly_order, j_matrix)

    !> coulomb supermatrix
    real(dp), intent(in) :: jj(0:,:,:,0:,:,:)

    !> density matrix supervector
    real(dp), intent(in) :: pp(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> coulomb matrix
    real(dp), intent(out) :: j_matrix(0:,:,:)

    !> auxiliary variables
    integer :: ii, jjj, kk, ll, mm, nn, oo, ppp, qq, rr, ss, tt, uu, vv

    j_matrix(:,:,:) = 0.0_dp

    do ii = 0, max_l
      ss = 0
      do jjj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
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

                        ! multiply Coulomb supermatrix with total density matrix supervector
                        j_matrix(ii, ss, tt) = j_matrix(ii, ss, tt)&
                            & + jj(ii, ss, tt, nn, uu, vv) * pp(nn, uu, vv)

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

  end subroutine build_coulomb_matrix


  !> Builds Hartree-Fock exchange matrix to be added to the Fock matrix from supermatrix by
  !! multiplying with density matrix supervector.
  subroutine build_hf_ex_matrix(kk, pp, max_l, num_alpha, poly_order, k_matrix)

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

    !> Hartree-Fock exchange matrix
    real(dp), intent(out) :: k_matrix(:,0:,:,:)

    !> auxiliary variables
    integer :: ii, jj, kkk, ll, mm, nn, oo, ppp, qq, rr, ss, tt, uu, vv

    k_matrix(:,:,:,:) = 0.0_dp

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

                        ! multiply hf exchange supermatrix with density matrix supervector per spin
                        k_matrix(1, ii, ss, tt) = k_matrix(1, ii, ss, tt)&
                            & + kk(ii, ss, tt, nn, uu, vv) * pp(1, nn, uu, vv)
                        k_matrix(2, ii, ss, tt) = k_matrix(2, ii, ss, tt)&
                            & + kk(ii, ss, tt, nn, uu, vv) * pp(2, nn, uu, vv)

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

  end subroutine build_hf_ex_matrix


  !> Builds DFT exchange matrix to be added to the Fock matrix by calculating the single matrix
  !! elements and putting them together.
  subroutine build_dft_exc_matrix(max_l, num_alpha, poly_order, alpha, num_mesh_points, abcissa,&
      & weight, vxc, k_matrix)

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

    !> xc potential on grid
    real(dp), intent(in) :: vxc(:,:)

    !> DFT exchange matrix
    real(dp), intent(out) :: k_matrix(:,0:,:,:)

    !> single matrix element of the exchange correlation potential
    real(dp) :: exc_matrixelement(2)

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, ss, tt, start

    k_matrix(:,:,:,:) = 0.0_dp
    exc_matrixelement(:) = 0.0_dp

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

              call dft_exc_matrixelement(num_mesh_points, weight, abcissa, vxc, alpha(ii, jj), kk,&
                  & alpha(ii, ll), mm, ii, exc_matrixelement)

              k_matrix(1, ii, ss, tt) = exc_matrixelement(1)
              k_matrix(2, ii, ss, tt) = exc_matrixelement(2)
              k_matrix(1, ii, tt, ss) = exc_matrixelement(1)
              k_matrix(2, ii, tt, ss) = exc_matrixelement(2)

            end do
          end do
        end do
      end do
    end do

  end subroutine build_dft_exc_matrix

end module hamiltonian
