!> Module that contains subroutines for the calculation of Gaunt coefficients for real and complex
!! spherical harmonics.
module common_anglib

  use common_accuracy, only : dp

  implicit none
  private

  public :: initGaunt, freeGaunt, realGaunt

  !> true, if Gaunt is initialized
  logical :: tGauntInit_ = .false.

  !> maximum 1st dimension of storeGaunt_
  integer :: iGauMax1_

  !> maximum 2nd dimension of storeGaunt_
  integer :: iGauMax2_

  !> array that stores Gaunts
  real(dp), allocatable :: storeGaunt_(:,:)

  !> transformation complx -> real Gaunt
  complex(dp), allocatable :: uMat_(:,:)


contains

  !> Initializes the calculation of Gaunt coefficients, precalculates non-zero values and stores
  !! these in array storeGaunt_. The transformation matrix from complex to real spherical harmonics
  !! is set up.
  !!
  !! @desc implements the index function approach of Pinchon et al., IJQC 107, pp. 2186 (2007)
  !! @note allows computation of Gaunts with angular momenta up to lMax, lMax, 2lMax in any order
  subroutine initGaunt(lMax)

    !> maximum of angular momenta of the two STOs to be fourier transformed
    integer, intent(in) :: lMax

    !! 1st and 2nd index of storeGaunt_
    integer :: i1, i2

    !! angular momentum quantum numbers
    integer :: j1, j2, j3

    !! magnetic quantum numbers
    integer :: m2, m3

    !! upper boundary for m2
    integer :: bMax

    !! index of transformation matrix uMat_
    integer :: jj

    !! 1/sqrt(2)
    real(dp), parameter :: rs2 = 0.7071067811865475_dp

    if (tGauntInit_) then
      print *, 'Gaunt already initialized'
      stop
    else
      tGauntInit_ = .true.
      write(*, '(A,I4)') "==> Initializing anglib: lMax = ", lMax
    end if

    call indexGaunt(2 * lMax, lMax, lMax, lMax, lMax, iGauMax1_, iGauMax2_)

    allocate(storeGaunt_(iGauMax1_,iGauMax2_))
    storeGaunt_(:,:) = 0.0_dp

    do j3 = 0, lMax
      do j2 = j3, lMax
        do j1 = j2, 2 * lMax
          if (mod(j1 + j2 + j3, 2) /= 0 .or. j1 > j2 + j3 .or. j2 > j1 + j3 .or. j3 > j1 + j2) cycle
          do m3 = 0, j3
            if (j1 == j2) then
              bMax = - (m3 + 1) / 2
            elseif (m3 == 0) then
              bMax = 0
            elseif (j2 == j3) then
              bMax = min(j1 - m3, m3)
            else
              bMax = min(j1 - m3, j2)
            end if
            do m2 = - j2, bMax
              call indexGaunt(j1, j2, j3, m2, m3, i1, i2)
              !! phase definition of Pinchon differs from usual one
              !! g_Pinchon(j1, m1, j2, m2, j3, m3) = - 1^m3 g(j1, m1, j2, m2, j3, - m3)
              if (mod(m3, 2) == 0) then
                storeGaunt_(i1, i2) = gaunt(j1, - m2 - m3, j2, m2, j3, - m3)
              else
                storeGaunt_(i1, i2) = - gaunt(j1, - m2 - m3, j2, m2, j3, - m3)
              end if
            end do
          end do
        end do
      end do
    end do

    !! Setup of matrix for the transformation of gaunt coefficients based on complex spherical
    !! harmonics to real harmonics, see realGaunt for details.
    allocate(uMat_(-2*lMax:2*lMax, -2*lMax:2*lMax))
    uMat_(:,:) = (0.0_dp, 0.0_dp)
    uMat_(0, 0) = (1.0_dp, 0.0_dp)

    do jj = 1, 2 * lMax
      uMat_( jj, jj) = (1.0_dp, 0.0_dp) * rs2
      uMat_( jj, -jj) = (1.0_dp, 0.0_dp) * rs2 * real((-1)**jj, dp)
      uMat_(-jj, jj) = (0.0_dp, 1.0_dp) * rs2 * (- 1.0_dp)
      uMat_(-jj, -jj) = (0.0_dp, 1.0_dp) * rs2 * real((-1)**jj, dp)
    end do

  end subroutine initGaunt


  !> Frees an initialized Gaunt.
  subroutine freeGaunt()

    if (tGauntInit_) then
      deallocate(storeGaunt_)
      deallocate(uMat_)
      tGauntInit_ = .false.
    else
      write(*,*) 'freeGaunt: Gaunt not yet initialized!'
      stop
    end if

  end subroutine freeGaunt


  !> Implements the index function of Pinchon et al., IJQC 107, pp. 2186 (2007).
  !!
  !! @note Restrictions exist for the allowed parameter ranges, see Pinchon (e.g. j1 >= j2 >= j3).
  !!       These are not checked by the code!
  subroutine indexGaunt(j1, j2, j3, m2, m3, ind1, ind2)

    !> angular momentum quantum numbers
    integer, intent(in) :: j1, j2, j3

    !> magnetic quantum numbers
    integer, intent(in) :: m2, m3

    !> corresponding first index in storeGaunt_(:,:)
    integer, intent(out) :: ind1

    !> corresponding second index in storeGaunt_(:,:)
    integer, intent(out) :: ind2

    if (mod(j1, 2) == 0) then
      ind1 = (j1 * (3 * j1**3 + (32 * j1**2 - (12 * j1 + 32)))&
          & - 48 * (j2 * (j1 + 1) * (j1 - j2 - 2) - j3 * j3)) / 192 + m3 + 1
    else
      ind1 = (j1 * (3 * j1**3 + (32 * j1**2 + (18 * j1 + 64)))&
          & - 48 * (j1 * j2 * (j1 - j2 - 1) - j3 * j3) + 219) / 192 + m3
    end if

    ind2 = j2 + m2 + 1

  end subroutine indexGaunt


  !> Retrieves Gaunt coefficients \int Y_l1m1 Y_l2m2 Y_l3m3 dOmega for complex spherical
  !! harmonics Y_lm from precomputed array storeGaunt_ using the index function indexGaunt.
  !!
  !! @note Requires initialization with previous call to initGaunt(lMax).
  function getGaunt(j1, m1, j2, m2, j3, m3) result(res)

    !> angular momentum quantum numbers
    integer, intent(in) :: j1, j2, j3

    !> magnetic quantum numbers
    integer, intent(in) :: m1, m2, m3

    !! angular momentum quantum numbers in standard order
    integer :: l1, l2, l3

    !! magnetic quantum numbers in standard order
    integer :: n1, n2, n3

    !! 1st and 2nd index in storeGaunt_(:,:)
    integer :: ind1, ind2

    !! resulting Gaunt coefficient
    real(dp) :: res

    if (.not. tGauntInit_) then
      print *, 'getGaunt: Gaunt not yet initialized!'
      stop
    end if

    if ((mod(j1 + j2 + j3, 2) /= 0) .or. (m1+m2-m3 /= 0) .or. (j1 > j2+j3) .or. (j2 > j1+j3)&
        & .or. (j3 > j1+j2) .or. (abs(m1) > j1) .or. (abs(m2) > j2) .or. (abs(m3) > j3)) then
      res = 0.0_dp
      return
    end if

    l1 = j1
    l2 = j2
    l3 = j3
    n1 = m1
    n2 = m2
    n3 = - m3

    !! bring quantum numbers in standard order, with m3 > 0
    if (l2 == max(l1, l2, l3)) then
      call iSwap(l1, l2)
      call iSwap(n1, n2)
    elseif (l3 > l1) then
      call iSwap(l1, l3)
      call iSwap(n1, n3)
    end if
    if (l2 == min(l1, l2, l3)) then
      call iSwap(l2, l3)
      call iSwap(n2, n3)
    end if
    if (n3 < 0) then
      n1 = - n1
      n2 = - n2
      n3 = - n3
    end if

    call indexGaunt(l1, l2, l3, n2, n3, ind1, ind2)
    if (l1 == l2) then
      if (n2 > - (n3 + 1) / 2) then
        call indexGaunt(l1, l2, l3, - n2 - n3, n3, ind1, ind2)
      end if
    elseif (n3 == 0) then
      if (n2 > 0) then
        call indexGaunt(l1, l2, l3, - n2, 0, ind1, ind2)
      end if
    elseif (l2 == l3) then
      if (n2 > min(l1 - n3, n3)) then
        call indexGaunt(l1, l2, l3, n3, n2, ind1, ind2)
      end if
    end if
    if ((ind1 > iGauMax1_) .or. (ind2 > iGauMax2_)) then
      print *, 'getGaunt: Index out of range.'
      stop
    else
      !! Pinchons definition of Gaunt differs from usual one
      if (mod(m3, 2) == 0) then
        res = storeGaunt_(ind1, ind2)
      else
        res = - storeGaunt_(ind1, ind2)
      end if
    end if

  end function getGaunt


  !> Calculates a real Gaunt coefficient <j1 m1|j2 m2|j3 m3> = \int YR_l1m1 YR_l2m2 YR_l3m3 dOmega
  !! for real spherical harmonics YR_lm. Arguments are integer and NOT twice the true value!
  !!
  !! @desc  Uses transformation highlighted in Homeier et al., J. Mol. Struct. 368, pp. 31 (1996).
  !! @note  Requires initialization with previous call to initGaunt(lMax).
  function realGaunt(j1, m1, j2, m2, j3, m3) result(res)

    !> angular momentum quantum numbers
    integer, intent(in) :: j1, j2, j3

    !> magnetic quantum numbers
    integer, intent(in) :: m1, m2, m3

    !! resulting real Gaunt coefficient
    real(dp) :: res

    if (.not. tGauntInit_) then
      print *, 'realGaunt: Gaunt not yet initialized!'
      stop
    end if

    res = 0.0_dp

    if ((m1 /= abs(m2 + m3)) .and. (m1 /= abs(m2 - m3)) .and. (m2 /= abs(m3 + m1))&
        & .and. (m2 /= abs(m3 - m1)) .and. (m3 /= abs(m1 + m2)) .and. (m3 /= abs(m1 - m2))) return
    if (mod(j1 + j2 + j3, 2) /= 0) return

    if ((m1 /= 0) .and. (m2 /= 0) .and. (m3 /= 0)) then
      if (abs(m2 + m3) <= j1) then
        res = 2.0_dp * getGaunt(j2, m2, j3, m3, j1, m2 + m3)&
            & * real(conjg(uMat_(m1, m2 + m3)) * uMat_(m2, m2) * uMat_(m3, m3), dp)
      end if
      if (abs(m2 - m3) <= j1) then
        res = res + 2.0_dp * getGaunt(j2, m2, j3, -m3, j1, m2 - m3)&
            & * real(conjg(uMat_(m1, m2 - m3)) * uMat_(m2, m2) * uMat_(m3, - m3), dp)
      end if
    elseif (m1 == 0 .and. m2 /= 0 .and. m3 /= 0) then
      if (abs(m3) <= j2) then
        res = 2.0_dp * getGaunt(j3, m3, j1, 0, j2, m3)&
            & * real(conjg(uMat_(m2, m3)) * uMat_(m3, m3), dp)
      end if
    elseif (m2 == 0 .and. m3 /= 0 .and. m1 /= 0) then
      if (abs(m1) <= j3) then
        res = 2.0_dp * getGaunt(j1, m1, j2, 0, j3, m1)&
            & * real(conjg(uMat_(m3, m1)) * uMat_(m1, m1), dp)
      end if
    elseif ((m3 == 0) .and. (m1 /= 0) .and. (m2 /= 0)) then
      if (abs(m2) <= j1) then
        res = 2.0_dp * getGaunt(j2, m2, j3, 0, j1, m2)&
            & * real(conjg(uMat_(m1, m2)) * uMat_(m2, m2), dp)
      end if
    elseif ((m1 == 0) .and. (m2 == 0) .and. (m3 == 0)) then
      res = getGaunt(j2, 0, j3, 0, j1, 0)
    else
      res = 0.0_dp
    end if

  end function realGaunt


  !> Computes Gaunt coefficients int Y_l1m1 Y_l2m2 Y_l3m3 dOmega for complex spherical harmonics
  !! Y_lm from scratch.
  !!
  !! @desc Uses Clebsch-Gordan coefficients.
  pure function gaunt(j1, m1, j2, m2, j3, m3) result(res)

    !> angular momentum quantum numbers
    integer, intent(in) :: j1, j2, j3

    !> magnetic quantum numbers
    integer, intent(in) :: m1, m2, m3

    !! floating point angular/magnetic momenta
    real(dp) :: rj1, rm1, rj2, rm2, rj3, rm3

    !! real Clebsch-Gordan coefficients
    real(dp) :: rcg1, rcg2

    !! resulting Gaunt coefficient
    real(dp) :: res

    !! 1 / sqrt(4pi)
    real(dp), parameter :: rs4pi = 0.2820947917738782_dp

    rj1 = real(j1, dp)
    rm1 = real(m1, dp)
    rj2 = real(j2, dp)
    rm2 = real(m2, dp)
    rj3 = real(j3, dp)
    rm3 = real(m3, dp)

    call ned(rj1, rj2, rj3, 0.0_dp, 0.0_dp, 0.0_dp, rcg1)
    call ned(rj1, rj2, rj3, rm1, rm2, rm3, rcg2)

    res = rs4pi * sqrt((2.0_dp * rj1 + 1.0_dp) * (2.0_dp * rj2 + 1.0_dp) / (2.0_dp * rj3 + 1.0_dp))&
        & * rcg1 * rcg2

  end function gaunt


  !> ARTURO QUIRANTES SIERRA
  !! Department of Applied Physics, Faculty of Sciences
  !! University of Granada, 18071 Granada (SPAIN)
  !! http://www.ugr.es/local/aquiran/codigos.htm
  !! aquiran@ugr.es
  !!
  !! Last update: 20 May 2.003
  !! Ported to F90 (niehaus@bccms.uni-bremen.de)
  !!
  !! Subroutine NED
  !! to calculate Clebsch-Gordan coefficients
  !!
  !! You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine.
  !! Input:
  !! 	AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices)
  !! Output:
  !!    CG=C-G(AJ,BJ,CJ,AM,BM,CM)
  !!
  !! @param cg       Clebsch-Gordan coefficient
  !! @note Parameters are real not integer!
  pure subroutine ned(aj, bj, cj, am, bm, cm, cg)

    !> angular momentum quantum numbers
    real(dp), intent(in) :: aj, bj, cj

    !> magnetic quantum numbers
    real(dp), intent(in) :: am, bm, cm

    !> Clebsch-Gordan coefficient
    real(dp), intent(out) :: cg

    !!
    real(dp) :: qq(100, 100)

    !!
    real(dp) :: dd, fn, xx, pp

    !!
    integer :: ii, zz, ld

    !!
    integer :: ja, ma, jb, mb, jc, mc, la, lb, lc, lt

    !!
    integer :: ja2, jb2, jc2, i2, k0, k1, kk, ip

    zz = max(2 * aj + 1, 2 * bj + 1, 2 * cj + 1, aj + bj + cj, aj + am, bj + bm, cj + cm) + 2

    do ii = 1, zz
      qq(ii, 1) = 1.0_dp
      qq(ii, ii) = 1.0_dp
    end do

    do ii = 2, zz - 1
      do kk = 2, ii
        qq(ii + 1, kk) = qq(ii, kk - 1) + qq(ii, kk)
      end do
    end do

    cg = 0.0_dp
    ja = aj + am + 1.01_dp
    ma = aj - am + 1.01_dp
    jb = bj + bm + 1.01_dp
    mb = bj - bm + 1.01_dp
    jc = cj + cm + 1.01_dp
    mc = cj - cm + 1.01_dp
    la = bj + cj - aj + 1.01_dp
    lb = cj + aj - bj + 1.01_dp
    lc = aj + bj - cj + 1.01_dp
    lt = aj + bj + cj + 1.01_dp
    dd = abs(am + bm - cm) - 0.01_dp

    if (dd <= 0.0_dp) then
      ld = min(ja, jb, jc, ma, mb, mc, la, lb, lc)
    else
      return
    end if

    if (ld > 0.0_dp) then
      ja2 = aj + aj + am + am
    else
      return
    end if

    jb2 = bj + bj + bm + bm
    jc2 = cj + cj - cm - cm
    i2 = ja2 + jb2 + jc2 - ja2 / 2 * 2 - jb2 / 2 * 2 - jc2 / 2 * 2

    if (i2 == 0) then
      fn = qq(ja + ma - 1, lc) / qq(lt, jc + mc - 1)
    else
      return
    end if

    fn = fn * qq(jb + mb - 1, lc) / qq(lt + 1, 2)
    fn = fn / qq(ja + ma - 1, ja)
    fn = fn / qq(jb + mb - 1, jb)
    fn = fn / qq(jc + mc - 1, jc)

    k0 = max(0, lc - ja, lc - mb) + 1
    k1 = min(lc, ma, jb)

    xx = 0.0_dp

    do kk = k0, k1
      xx = - xx - qq(lc, kk) * qq(lb, ma - kk + 1) * qq(la, jb - kk + 1)
    end do

    ip = k1 + lb + jc
    pp = 1 - 2 * (ip - ip / 2 * 2)
    cg = pp * xx * dsqrt(fn)

    !! What we've calculated is a Wigner 3-j coefficient.
    !! Next, we'll turn it into a Clebsch-Gordan coefficient.
    cg = cg * sqrt(2 * cj + 1) * (-1)**nint(aj - bj - cm)

  end subroutine ned


  !> Swaps two indices.
  pure subroutine iSwap(ii, jj)

    !> indices to swap
    integer, intent(inout) :: ii, jj

    !! temporary storage for ii
    integer :: iTmp

    iTmp = ii
    ii = jj
    jj = iTmp

  end subroutine iSwap

end module common_anglib
