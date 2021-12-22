!> Module that builds up Coulomb and Hartree-Fock exchange supermatrix elements.
module coulomb_hfex

  use common_accuracy, only : dp
  use utilities, only : fak
  use core_overlap, only : v

  implicit none
  private

  public :: coulomb, hfex


contains

  !> Calculates Coulomb supermatrix elements,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 6 and eqn. 21
  pure subroutine coulomb(jj, max_l, num_alpha, alpha, poly_order, uu, ss)

    !> coulomb supermatrix
    real(dp), intent(out) :: jj(0:,:,:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> nucleus-electron supervector
    real(dp), intent(out) :: uu(0:,:,:)

    !> overlap supervector
    real(dp), intent(out) :: ss(0:,:,:)

    !> temporary storage
    real(dp) :: alpha1, alpha2

    !> auxiliary variables
    integer :: ii, jjj, kk, ll, mm, nn, oo, pp, qq, rr, sss, tt, uuu, vv
    integer :: nlpq, nmrs

    jj(:,:,:,:,:,:) = 0.0_dp

    do ii = 0, max_l
      sss = 0
      do jjj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          sss = sss + 1
          tt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              tt = tt + 1
              do nn = 0, max_l
                uuu = 0
                do oo = 1, num_alpha(nn)
                  do pp = 1, poly_order(nn)
                    uuu = uuu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        alpha1 = (alpha(ii, jjj) + alpha(ii, ll)) /&
                            &(alpha(nn, oo) + alpha(nn, qq))
                        alpha2 = (alpha(nn, oo) + alpha(nn, qq)) /&
                            &(alpha(ii, jjj) + alpha(ii, ll))
                        nlpq = kk + mm + 2 * ii
                        nmrs = pp + rr + 2 * nn

                        jj(ii, sss, tt, nn, uuu, vv) =&
                            & uu(ii, sss, tt) * ss(nn, uuu, vv) *&
                            & cc(nlpq - 1, nmrs, alpha1) +&
                            & uu(nn, uuu, vv) * ss(ii, sss, tt) *&
                            & cc(nmrs - 1, nlpq, alpha2)

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

  end subroutine coulomb


  !> Builds HF exchange supermatrix,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 7/8 and eqn. 21
  pure subroutine hfex(kk, max_l, num_alpha, alpha, poly_order, problemsize)

    !> Hartree-Fock exchange supermatrix
    real(dp), intent(out) :: kk(0:,:,:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> Rev. Mod. Phys. 32, 186 (1960) eqn. 21
    real(dp), allocatable :: knu(:,:,:,:,:,:,:)

    !> temporary storage
    real(dp) :: alpha1, alpha2, alpha3, alpha4, beta1, beta2, beta3, beta4

    !> four terms of Rev. Mod. Phys. 32, 186 (1960) eqn. 21
    real(dp) :: pre, t1, t2, t3, t4

    !> auxiliary variables
    integer :: ii, jj, kkk, ll, mm, nn, oo, pp, qq, rr, ss, tt, uu, vv
    integer :: nu, nlp, nlq, nmr, nms

    allocate(knu(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize,&
        & 0:2 * max_l + 2))

    kk(:,:,:,:,:,:) = 0.0_dp
    knu(:,:,:,:,:,:,:) = 0.0_dp

    ! build knu according to Rev. Mod. Phys. 32, 186 (1960) eqn. 8

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
                  do pp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(nn, oo))
                        alpha2 = 0.5_dp * (alpha(ii, ll) + alpha(nn, qq))
                        alpha3 = 0.5_dp * (alpha(ii, jj) + alpha(nn, qq))
                        alpha4 = 0.5_dp * (alpha(ii, ll) + alpha(nn, oo))
                        beta1 = alpha1 / alpha2
                        beta2 = alpha2 / alpha1
                        beta3 = alpha3 / alpha4
                        beta4 = alpha4 / alpha3
                        nlp = kkk + ii
                        nlq = mm + ii
                        nmr = pp + nn
                        nms = rr + nn

                        pre = 1.0_dp / sqrt(&
                            & v(alpha(ii, jj), 2 * (kkk + ii))&
                            & * v(alpha(ii, ll), 2 * (mm + ii))&
                            & * v(alpha(nn, oo), 2 * (pp + nn))&
                            & * v(alpha(nn, qq), 2 * (rr + nn)))

                        do nu = abs(ii - nn), ii + nn, 2

                          t1 = v(alpha1, nlp + nmr - nu - 1) * v(alpha2, nlq + nms + nu) *&
                              & cc(nlp + nmr - nu - 1, nlq + nms + nu, beta1)
                          t2 = v(alpha2, nlq + nms - nu - 1) * v(alpha1, nlp + nmr + nu) *&
                              & cc(nlq + nms - nu - 1, nlp + nmr + nu, beta2)
                          t3 = v(alpha3, nlp + nms - nu - 1) * v(alpha4, nlq + nmr + nu) *&
                              & cc(nlp + nms - nu - 1, nlq + nmr + nu, beta3)
                          t4 = v(alpha4, nlq + nmr - nu - 1) * v(alpha3, nlp + nms + nu) *&
                              & cc(nlq + nmr - nu - 1, nlp + nms + nu, beta4)

                          knu(ii, ss, tt, nn, uu, vv, nu) = pre * (t1 + t2 + t3 + t4)

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
    end do

    ! build kk according to Rev. Mod. Phys. 32, 186 (1960) eqn. 7

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
                  do pp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        do nu = abs(ii - nn), ii + nn, 2

                          kk(ii, ss, tt, nn, uu, vv) = kk(ii, ss, tt, nn, uu, vv) +&
                              & almn(ii, nn, nu) * knu(ii, ss, tt, nn, uu, vv, nu)

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
    end do

  end subroutine hfex


  !> Auxiliary function,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 22 and eqn. 23
  pure function cc(alpha, beta, tt)

    !> indices of Rev. Mod. Phys. 32, 186 (1960) eqn. 22 and eqn. 23
    integer, intent(in) :: alpha, beta

    !> argument of Rev. Mod. Phys. 32, 186 (1960) eqn. 22 and eqn. 23
    real(dp), intent(in) :: tt

    !> result
    real(dp) :: cc

    !> prefactor of Rev. Mod. Phys. 32, 186 (1960) eqn. 23
    real(dp) :: factor

    !> array that holds all C_{\alpha, \beta}(t)
    real(dp), allocatable :: carray(:,:)

    !> alpha and beta indices to carry out the recursion
    integer :: aa, bb

    ! early return if index smaller than zero

    if (alpha < 0) then
      cc = 0.0_dp
      return
    end if

    if (beta < 0) then
      cc = 0.0_dp
      return
    end if

    allocate(carray(0:alpha, 0:beta))

    factor = 1.0_dp / (1.0_dp + tt)

    ! Overall this is naive, the matrix could be reused to some extent ...
    ! OTOH, the matrices are relatively small.

    ! first handle Kronecker delta, three cases
    carray(0, 0) = factor
    do aa = 1, alpha
      carray(aa, 0) = factor * (tt * carray(aa - 1, 0) + 1.0_dp)
    end do
    do bb = 1, beta
      carray(0, bb) = factor * (carray(0, bb - 1))
    end do

    ! now build up from alpha, beta = 1
    do aa = 1, alpha
      do bb = 1, beta
        carray(aa, bb) = factor * (tt * carray(aa - 1, bb) + carray(aa, bb - 1))
      end do
    end do

    cc = carray(alpha, beta)

  end function cc


  !> Auxiliary function,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 9
  pure function aa(rho)

    !> even integer, according to Rev. Mod. Phys. 32, 186 (1960) eqn. 9
    integer, intent(in) :: rho

    !> result
    real(dp) :: aa

    aa = fak(rho) / ((fak(rho / 2))**2)

  end function aa


  !> Auxiliary function,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 9
  pure function almn(lambda, mu, nu)

    !> indices of Rev. Mod. Phys. 32, 186 (1960) eqn. 9
    integer, intent(in) :: lambda, mu, nu

    !> result
    real(dp) :: almn

    almn = aa(lambda + mu - nu) * aa(lambda - mu + nu) * aa(mu - lambda + nu)&
        & / (real(lambda + mu + nu + 1, dp) * aa(lambda + mu + nu))

  end function almn

end module coulomb_hfex
