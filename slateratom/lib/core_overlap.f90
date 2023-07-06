!> Module that builds up various supervectors.
module core_overlap

  use common_accuracy, only : dp
  use utilities, only : fak

  implicit none
  private

  public :: overlap, kinetic, nuclear, moments, v, confinement


  interface v
    module procedure v_int, v_real
  end interface


contains

  !> Calculates overlap matrix elements,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 5 and eqn. 19.
  !!
  !! Definition of the primitive basis functions based on Roothaan:
  !! R_{\lambda p}=1/sqrt((2n_{\lambda p})!)*
  !!               (2*\zeta_{\lambda p})**(n_{\lambda p}+0.5)*
  !!               r**(n_{\lambda p}-1)*exp(-\zeta_{\lambda p}*r)
  !!
  !! For every exponent \zeta_{\lambda p} there are num_power coefficients,
  !! each connected to one r**(n_{\lambda p}-1). The sum over all
  !! coefficients, e.g. implicitely \zeta and r**n, gives the usual DFTB
  !! basis function.
  !!
  !! Note: in DFTB one usually has r**(n+l-1) explicitely, here the angular
  !! momentum index l is implicit. Result:
  !! for l=0, e.g. s, n_{\lambda p}=0,1,..., num_power
  !! for l=1, e.g. p, n_{\lambda p}=1,2,..., num_power + 1
  !! for l=2, e.g. d, n_{\lambda p}=2,3,..., num_power + 2
  pure subroutine overlap(ss, max_l, num_alpha, alpha, poly_order)

    !> overlap supervector
    real(dp), intent(out) :: ss(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> temporary storage
    real(dp) :: alpha1

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, nlp, nlq

    ss(:,:,:) = 0.0_dp

    ! these loops define the indizes S_{\lambda p q}
    ! p=alpha1/n=0+l,alpha1/n=1+l,...,alpha2/n=0+l,alpha2/n=1+l...

    do ii = 0, max_l
      nn = 0
      do jj = 1, num_alpha(ii)
        do ll = 1, poly_order(ii)
          nn = nn + 1
          oo = 0
          nlp = ll + ii
          do kk = 1, num_alpha(ii)
            alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
            do mm = 1, poly_order(ii)
              oo = oo + 1
              nlq = mm + ii
              ! use ll+ii and mm+ii because of DFTB basis function definition
              ss(ii, nn, oo) = 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                  & * v(alpha(ii, kk), 2 * nlq)) * v(alpha1, nlp + nlq)
            end do
          end do
        end do
      end do
    end do

  end subroutine overlap


  !> Calculates nuclear attraction matrix elements,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 5 and eqn. 19.
  pure subroutine nuclear(uu, max_l, num_alpha, alpha, poly_order)

    !> nucleus-electron supervector
    real(dp), intent(out) :: uu(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> temporary storage
    real(dp) :: alpha1

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, nlp, nlq

    uu = 0.0_dp

    do ii = 0, max_l
      nn = 0
      do jj = 1, num_alpha(ii)
        do ll = 1, poly_order(ii)
          nn = nn + 1
          oo = 0
          nlp = ll + ii
          do kk = 1, num_alpha(ii)
            alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
            do mm = 1, poly_order(ii)
              oo = oo + 1
              nlq = mm + ii
              uu(ii, nn, oo) = 2.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp) *&
                  &v(alpha(ii, kk), 2 * nlq)) * v(alpha1, nlp + nlq - 1)
            end do
          end do
        end do
      end do
    end do

  end subroutine nuclear


  !> Calculates the kinetic matrix elements,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 5 and eqn. 19.
  pure subroutine kinetic(tt, max_l, num_alpha, alpha, poly_order)

    !> kinetic supervector
    real(dp), intent(out) :: tt(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> temporary storage
    real(dp) :: alpha1

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, nlp, nlq

    tt(:,:,:) = 0.0_dp

    do ii = 0, max_l
      nn = 0
      do jj = 1, num_alpha(ii)
        do ll = 1, poly_order(ii)
          nn = nn + 1
          oo = 0
          nlp = ll + ii
          do kk = 1, num_alpha(ii)
            alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
            do mm = 1, poly_order(ii)
              oo = oo + 1
              nlq = mm + ii
              tt(ii, nn, oo) = 0.5_dp * alpha(ii, jj) * alpha(ii, kk) /&
                  & sqrt(v(alpha(ii, jj), 2 * nlp) * v(alpha(ii, kk), 2 * nlq)) *&
                  & (v(alpha1, nlp + nlq) -&
                  & (w(alpha(ii, jj), ii, nlp) + w(alpha(ii, kk), ii, nlq)) *&
                  & v(alpha1, nlp + nlq - 1) +&
                  & (w(alpha(ii, jj), ii, nlp) * w(alpha(ii, kk), ii, nlq)) *&
                  & v(alpha1, nlp + nlq - 2))
            end do
          end do
        end do
      end do
    end do

  end subroutine kinetic


  !> Calculates analytic matrix elements of confining potential.
  !! No checking for power, e.g. power==0 or power<0 etc. !
  pure subroutine confinement(vconf, max_l, num_alpha, alpha, poly_order, conf_r0, conf_power)

    !> confinement supervector
    real(dp), intent(out) :: vconf(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

     !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> confinement radii
    real(dp), intent(in) :: conf_r0(0:)

    !> power of confinement
    real(dp), intent(in) :: conf_power(0:)

    !> temporary storage
    real(dp) :: alpha1

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, nlp, nlq

    vconf(:,:,:) = 0.0_dp

    do ii = 0, max_l
      if (conf_power(ii) > 1.0e-06_dp) then
        nn = 0
        do jj = 1, num_alpha(ii)
          do ll = 1, poly_order(ii)
            nn = nn + 1
            oo = 0
            nlp = ll + ii
            do kk = 1, num_alpha(ii)
              alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
              do mm = 1, poly_order(ii)
                oo = oo + 1
                nlq = mm + ii
                vconf(ii, nn, oo) = 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                    & * v(alpha(ii, kk), 2 * nlq)) / (conf_r0(ii) * 2.0_dp)**conf_power(ii)&
                    & * v(alpha1, nlp + nlq + conf_power(ii))
              end do
            end do
          end do
        end do
      end if
    end do

  end subroutine confinement


  !> Calculates arbitrary moments of electron distribution, e.g. expectation values of <r>, <r^2>
  !! etc.; this is implemented analytically for arbitrary powers.
  pure subroutine moments(moment, max_l, num_alpha, alpha, poly_order, cof, power)

    !> moment of electron distribution
    real(dp), intent(out) :: moment(:,0:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> wavefunction coefficients
    real(dp), intent(in) :: cof(:,0:,:,:)

    !> power of moment
    integer, intent(in) :: power

    !> temporary storage
    real(dp) :: alpha1

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, pp, nlp, nlq

    moment(:,:,:) = 0.0_dp

    ! <r^-3> only computed for p-functions and higher
    if (power > -3) then
      do ii = 0, max_l
        do pp = 1, num_alpha(ii) * poly_order(ii)
          nn = 0
          do jj = 1, num_alpha(ii)
            do ll = 1, poly_order(ii)
              nn = nn + 1
              oo = 0
              nlp = ll + ii
              do kk = 1, num_alpha(ii)
                alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
                do mm = 1, poly_order(ii)
                  oo = oo + 1
                  nlq = mm + ii

                  moment(1, ii, pp) = moment(1, ii, pp) + 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                      & * v(alpha(ii, kk), 2 * nlq)) / (2.0_dp**power)&
                      & * v(alpha1, nlp + nlq + power) * cof(1, ii, nn, pp) * cof(1, ii, oo, pp)

                  moment(2, ii, pp) = moment(2, ii, pp) + 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                      & * v(alpha(ii, kk), 2 * nlq)) / (2.0_dp**power)&
                      & * v(alpha1, nlp + nlq + power) * cof(2, ii, nn, pp) * cof(2, ii, oo, pp)

                end do
              end do
            end do
          end do
        end do
      end do
    else if (power == -3) then
      do ii = 1, max_l
        do pp = 1, num_alpha(ii) * poly_order(ii)
          nn = 0
          do jj = 1, num_alpha(ii)
            do ll = 1, poly_order(ii)
              nn = nn + 1
              oo = 0
              nlp = ll + ii
              do kk = 1, num_alpha(ii)
                alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(ii, kk))
                do mm = 1, poly_order(ii)
                  oo = oo + 1
                  nlq = mm + ii

                  moment(1, ii, pp) = moment(1, ii, pp) + 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                      & * v(alpha(ii, kk), 2 * nlq)) / (2.0_dp**power)&
                      & * v(alpha1, nlp + nlq + power) * cof(1, ii, nn, pp) * cof(1, ii, oo, pp)

                  moment(2, ii, pp) = moment(2, ii, pp) + 1.0_dp / sqrt(v(alpha(ii, jj), 2 * nlp)&
                      & * v(alpha(ii, kk), 2 * nlq)) / (2.0_dp**power)&
                      & * v(alpha1, nlp + nlq + power) * cof(2, ii, nn, pp) * cof(2, ii, oo, pp)

                end do
              end do
            end do
          end do
        end do
      end do
    end if

  end subroutine moments


  !> Auxiliary function V_{i}(x), see Rev. Mod. Phys. 32, 186 (1960) eqn. 20.
  pure function v_int(xx, ii) result(res)

    !> argument
    real(dp), intent(in) :: xx

    !> index
    integer, intent(in) :: ii

    !> result
    real(dp) :: res

    res = fak(ii) / (xx**(ii + 1))

  end function v_int


  !> Auxiliary function V_{i}(x), see Rev. Mod. Phys. 32, 186 (1960) eqn. 20.
  pure function v_real(xx, ii) result(res)

    !> argument
    real(dp), intent(in) :: xx

    !> index
    real(dp), intent(in) :: ii

    !> result
    real(dp) :: res

    res = gamma(ii + 1.0_dp) / xx**(ii + 1.0_dp)

  end function v_real


  !> Auxiliary function W_{ij}(x), see Rev. Mod. Phys. 32, 186 (1960) eqn. 20.
  pure function w(xx, ii, jj) result(res)

    !> argument
    real(dp), intent(in) :: xx

    !> indices
    integer, intent(in) :: ii, jj

    !> result
    real(dp) :: res

    res = 2.0_dp * real(jj - ii - 1, dp) / xx

  end function w

end module core_overlap
