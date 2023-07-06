!> Module that provides functionality to calculate the electron density at a radial point in space.
module density

  use common_accuracy, only : dp
  use utilities, only : fak

  implicit none
  private

  public :: density_at_point, density_at_point_1st, density_at_point_2nd
  public :: wavefunction, wavefunction_1st, wavefunction_2nd
  public :: basis, basis_1st, basis_2nd
  public :: basis_times_basis, basis_times_basis_1st, basis_times_basis_2nd
  public :: basis_1st_times_basis_1st, basis_2nd_times_basis_2nd
  public :: basis_times_basis_times_r2, basis_times_basis_1st_times_r2,&
      & basis_times_basis_2nd_times_r2, basis_times_basis_1st_times_r,&
      & basis_1st_times_basis_1st_times_r2


contains

  !> Calculates electron density at a radial point in space.
  pure function density_at_point(pp, max_l, num_alpha, poly_order, alpha, rr)

    !> density matrix supervector
    real(dp), intent(in) :: pp(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> resulting electron density at a radial point in space
    real(dp) :: density_at_point

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, start

    density_at_point = 0.0_dp

    do ii = 0, max_l
      ll = 0
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          ll = ll + 1

          ! set global index correctly
          oo = ll - 1
          do mm = jj, num_alpha(ii)

            ! catch start index for polynomials, different depending on alpha block
            start = 1
            if (mm == jj) start = kk

            do nn = start, poly_order(ii)
              oo = oo + 1

              if (ll == oo) then
                density_at_point = density_at_point + pp(ii, ll, oo)&
                    & * basis_times_basis(alpha(ii, jj), kk, alpha(ii, mm), nn, ii, rr)
              end if

              if (ll /= oo) then
                density_at_point = density_at_point + 2.0_dp * pp(ii, ll, oo)&
                    & * basis_times_basis(alpha(ii, jj), kk, alpha(ii, mm), nn, ii, rr)
              end if

            end do
          end do
        end do
      end do
    end do

  end function density_at_point


  !> Calculates 1st derivative at a radial point in space analytically.
  pure function density_at_point_1st(pp, max_l, num_alpha, poly_order, alpha, rr)

    !> density matrix supervector
    real(dp), intent(in) :: pp(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> resulting analytical 1st derivative at a radial point in space
    real(dp) :: density_at_point_1st

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, start

    density_at_point_1st = 0.0_dp

    do ii = 0, max_l
      ll = 0
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          ll = ll + 1

          ! set global index correctly
          oo = ll - 1
          do mm = jj, num_alpha(ii)

            ! catch start index for polynomials, different depending on alpha block
            start = 1
            if (mm == jj) start = kk

            do nn = start, poly_order(ii)
              oo = oo + 1

              if (ll == oo) then
                density_at_point_1st = density_at_point_1st + pp(ii, ll, oo) * (&
                    &basis_times_basis_1st(alpha(ii, jj), kk, alpha(ii, mm), nn, ii, rr) +&
                    &basis_times_basis_1st(alpha(ii, mm), nn, alpha(ii, jj), kk, ii, rr))
              end if

              if (ll /= oo) then
                density_at_point_1st = density_at_point_1st + 2.0_dp * pp(ii, ll, oo) * (&
                    &basis_times_basis_1st(alpha(ii, jj), kk, alpha(ii, mm), nn, ii, rr) +&
                    &basis_times_basis_1st(alpha(ii, mm), nn, alpha(ii, jj), kk, ii, rr))
              end if

            end do
          end do
        end do
      end do
    end do

  end function density_at_point_1st


  !> Calculates 2nd derivative at a radial point in space analytically.
  pure function density_at_point_2nd(pp, max_l, num_alpha, poly_order, alpha, rr)

    !> density matrix supervector
    real(dp), intent(in) :: pp(0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> resulting analytical 2nd derivative at a radial point in space
    real(dp) :: density_at_point_2nd

    !> auxiliary variables
    integer :: ii, jj, kk, ll, mm, nn, oo, start

    density_at_point_2nd = 0.0_dp

    do ii = 0, max_l
      ll = 0
      do jj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          ll = ll + 1

          ! set global index correctly
          oo = ll - 1
          do mm = jj, num_alpha(ii)

            ! catch start index for polynomials, different depending on alpha block
            start = 1
            if (mm == jj) start = kk

            do nn = start, poly_order(ii)
              oo = oo + 1

              if (ll == oo) then
                density_at_point_2nd = density_at_point_2nd + pp(ii, ll, oo) * (&
                    & basis_times_basis_2nd(alpha(ii, jj), kk, alpha(ii, mm), nn, ii, rr)&
                    & + 2.0_dp * basis_1st_times_basis_1st(alpha(ii, jj), kk, alpha(ii, mm),&
                    & nn, ii, rr)&
                    & + basis_times_basis_2nd(alpha(ii, mm), nn, alpha(ii, jj), kk, ii, rr))
              end if

              if (ll /= oo) then
                density_at_point_2nd = density_at_point_2nd + 2.0_dp * pp(ii, ll, oo) * (&
                    & basis_times_basis_2nd(alpha(ii, jj), kk, alpha(ii, mm), nn, ii, rr)&
                    & + 2.0_dp * basis_1st_times_basis_1st(alpha(ii, jj), kk, alpha(ii, mm),&
                    & nn, ii, rr) +&
                    & basis_times_basis_2nd(alpha(ii, mm), nn, alpha(ii, jj), kk, ii, rr))
              end if

            end do
          end do
        end do
      end do
    end do

  end function density_at_point_2nd


  !> Calculates wavefunction at a radial point in space analytically.
  pure function wavefunction(cof, alpha, num_alpha, poly_order, ang, rr)

    !> expansion coefficients
    real(dp), intent(in) :: cof(:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> angular momentum
    integer, intent(in) :: ang

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> wavefunction at a radial point in space
    real(dp) :: wavefunction

    !> auxiliary variables
    integer :: ii, jj, kk

    wavefunction = 0.0_dp
    kk = 0

    do ii = 1, num_alpha(ang)
      do jj = 1, poly_order(ang)
        kk = kk + 1
        wavefunction = wavefunction + cof(kk) * basis(alpha(ang, ii), jj, ang, rr)
      end do
    end do

  end function wavefunction


  !> Calculates 1st derivative of wavefunction at a radial point in space analytically.
  pure function wavefunction_1st(cof, alpha, num_alpha, poly_order, ang, rr)

    !> expansion coefficients
    real(dp), intent(in) :: cof(:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> angular momentum
    integer, intent(in) :: ang

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> 1st derivative of wavefunction at a radial point in space
    real(dp) :: wavefunction_1st

    !> auxiliary variables
    integer :: ii, jj, kk

    wavefunction_1st = 0.0_dp
    kk = 0

    do ii = 1, num_alpha(ang)
      do jj = 1, poly_order(ang)
        kk = kk + 1
        wavefunction_1st = wavefunction_1st + cof(kk) * basis_1st(alpha(ang, ii), jj, ang, rr)
      end do
    end do

  end function wavefunction_1st


  !> Calculates 2nd derivative of wavefunction at a radial point in space analytically.
  pure function wavefunction_2nd(cof, alpha, num_alpha, poly_order, ang, rr)

    !> expansion coefficients
    real(dp), intent(in) :: cof(:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> angular momentum
    integer, intent(in) :: ang

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> 2nd derivative of wavefunction at a radial point in space
    real(dp) :: wavefunction_2nd

    !> auxiliary variables
    integer :: ii, jj, kk

    wavefunction_2nd = 0.0_dp
    kk = 0

    do ii = 1, num_alpha(ang)
      do jj = 1, poly_order(ang)
        kk = kk + 1
        wavefunction_2nd = wavefunction_2nd + cof(kk) * basis_2nd(alpha(ang, ii), jj, ang, rr)
      end do
    end do

  end function wavefunction_2nd


  !> Evaluates a primitive Slater basis function at a radial point in space,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 3.
  pure function basis(alpha, poly_order, ll, rr)

    !> basis exponent
    real(dp), intent(in) :: alpha

    !> highest polynomial order in shell
    integer, intent(in) :: poly_order

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> value of a primitive Slater basis function at a radial point in space
    real(dp) :: basis

    !> normalization pre-factor
    real(dp) :: normalization

    normalization = (2.0_dp * alpha)**(poly_order + ll) * sqrt(2.0_dp * alpha)&
        & / sqrt(fak(2 * (poly_order + ll)))

    ! catch 0^0
    if ((rr == 0.0_dp) .and. ((poly_order + ll - 1) == 0)) then
      basis = normalization * exp(- alpha * rr)
    else
      basis = normalization * rr**(poly_order + ll - 1) * exp(- alpha * rr)
    end if

  end function basis


  !> Evaluates 1st derivative of a primitive Slater basis function at a radial point in space.
  pure function basis_1st(alpha, poly_order, ll, rr)

    !> basis exponent
    real(dp), intent(in) :: alpha

    !> highest polynomial order in shell
    integer, intent(in) :: poly_order

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> 1st derivative of a primitive Slater basis function at a radial point in space
    real(dp) :: basis_1st

    !> normalization pre-factor
    real(dp) :: normalization

    normalization = (2.0_dp * alpha)**(poly_order + ll) * sqrt(2.0_dp * alpha)&
        & / sqrt(fak(2 * (poly_order + ll)))

    ! catch 0^0, setting 0^0=1 and 0^-1=0.0
    if ((rr == 0.0_dp) .and. ((poly_order + ll - 1) == 0)) then
      basis_1st = normalization * (- alpha * exp(- alpha * rr))
    else if ((rr == 0.0_dp) .and. ((poly_order + ll - 2) == 0)) then
      basis_1st = normalization * (real(poly_order + ll - 1, dp)&
          & * exp(- alpha * rr) - alpha * rr**(poly_order + ll - 1) * exp(- alpha * rr))
    else
      basis_1st = normalization * (real(poly_order + ll - 1, dp) * rr**(poly_order + ll - 2)&
          & * exp(- alpha * rr) - alpha * rr**(poly_order + ll - 1) * exp(- alpha * rr))
    end if

  end function basis_1st


  !> Evaluates 2nd derivative of a primitive Slater basis function at a radial point in space.
  pure function basis_2nd(alpha, poly_order, ll, rr)

    !> basis exponent
    real(dp), intent(in) :: alpha

    !> highest polynomial order in shell
    integer, intent(in) :: poly_order

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> 2nd derivative of a primitive Slater basis function at a radial point in space
    real(dp) :: basis_2nd

    !> normalization pre-factor
    real(dp) :: normalization

    normalization = (2.0_dp * alpha)**(poly_order + ll) * sqrt(2.0_dp * alpha)&
        & / sqrt(fak(2 * (poly_order + ll)))

    ! catch 0^0
    if ((rr == 0.0_dp) .and. ((poly_order + ll - 3) == 0)) then
      basis_2nd = normalization * (real(poly_order + ll - 1, dp) * real(poly_order + ll - 2, dp)&
          & * exp(- alpha * rr))
    else if ((rr == 0.0_dp) .and. ((poly_order + ll - 2) == 0)) then
      basis_2nd = normalization * (- 2.0_dp * alpha * real(poly_order + ll - 1, dp)&
          & * exp(- alpha * rr))
    else if ((rr == 0.0_dp) .and. ((poly_order + ll - 1) == 0)) then
      basis_2nd = normalization * (alpha**2 * exp(- alpha * rr))
    else
      basis_2nd = normalization * (real(poly_order + ll - 1, dp) * real(poly_order + ll - 2, dp)&
          & * rr**(poly_order + ll - 3) * exp(- alpha * rr) - 2.0_dp * alpha&
          & * real(poly_order + ll - 1, dp) * rr**(poly_order + ll - 2) * exp(- alpha * rr)&
          & + alpha**2 * rr**(poly_order + ll - 1) * exp(- alpha * rr))
    end if

  end function basis_2nd


  !> Evaluates product of two primitive Slater basis functions at a radial point in space.
  !! r^(m-1)*e^(-alpha*r)*r^(n-1)*exp(-beta*r)=r^(m+n-2)*exp(-(alpha+beta)*r)
  pure function basis_times_basis(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of 1st basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis
    real(dp), intent(in) :: beta

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of two primitive Slater basis functions at a radial point in space
    real(dp) :: basis_times_basis

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! catch 0^0
    if ((rr == 0.0_dp) .and. ((mm + nn - 2) == 0)) then
      basis_times_basis = normalization1 * normalization2 * exp(ab * rr)
    else
      basis_times_basis = normalization1 * normalization2 * rr**(mm + nn - 2) * exp(ab * rr)
    end if

    if (abs(basis_times_basis) < 1.0d-20) basis_times_basis = 0.0_dp

  end function basis_times_basis


  !> Evaluates product of a basis function with 1st derivative of another basis function.
  !! beta and poly2 are the arguments of the derivative.
  pure function basis_times_basis_1st(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in basis shell
    integer, intent(in) :: poly1

    !> basis exponent of basis 1st derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in basis 1st derivative shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of a basis function with the 1st derivative of another basis function
    real(dp) :: basis_times_basis_1st

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !

    ! catch 0^0, setting 0^0=1 and 0^1=0
    if ((rr == 0.0_dp) .and. ((mm + nn - 2) == 0)) then
      basis_times_basis_1st = normalization1 * normalization2 * (- beta) * exp(ab * rr)
    elseif ((rr == 0.0_dp) .and. ((mm + nn - 3) == 0)) then
      basis_times_basis_1st = normalization1 * normalization2 * real(nn - 1, dp) * exp(ab * rr)
    else
      basis_times_basis_1st = normalization1 * normalization2&
          & * (real(nn - 1, dp) * rr**(mm + nn - 3) - beta * rr**(nn + mm - 2)) * exp(ab * rr)
    end if

    if (abs(basis_times_basis_1st) < 1.0d-20) basis_times_basis_1st = 0.0_dp

  end function basis_times_basis_1st


  !> Evaluates product of a basis function with 2nd derivative of another basis function.
  !! beta and poly2 are the arguments of the derivative.
  pure function basis_times_basis_2nd(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in basis shell
    integer, intent(in) :: poly1

    !> basis exponent of basis 2nd derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in basis 2nd derivative shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of a basis function with the 2nd derivative of another basis function
    real(dp) :: basis_times_basis_2nd

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab, positive, negative

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !
    positive = real((nn - 1) * (nn - 2), dp) * rr**(mm + nn - 4) + beta**2 * rr**(mm + nn - 2)
    negative = real(2 * (nn - 1), dp) * beta * rr**(nn + mm - 3)

    basis_times_basis_2nd = normalization1 * normalization2 * (positive - negative) * exp(ab * rr)

    if (abs(basis_times_basis_2nd) < 1.0d-20) basis_times_basis_2nd = 0.0_dp

  end function basis_times_basis_2nd


  !> Evaluates product of 1st derivatives of basis functions.
  !! beta and poly2 are the arguments of the 2nd basis.
  pure function basis_1st_times_basis_1st(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of 1st basis derivative
    real(dp), intent(in) :: alpha

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of 1st derivatives of basis functions
    real(dp) :: basis_1st_times_basis_1st

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab, positive, negative

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable!

    ! catch 0^0
    if ((rr == 0.0_dp) .and. ((mm + nn - 2) == 0)) then
      positive = alpha * beta
    elseif ((rr == 0.0_dp) .and. ((mm + nn - 4) == 0)) then
      positive = real((mm - 1) * (nn - 1), dp)
    else
      positive = real((mm - 1) * (nn - 1), dp) * rr**(mm + nn - 4)&
          & + alpha * beta * rr**(mm + nn - 2)
    end if

    if ((rr == 0.0_dp) .and. ((mm + nn - 3) == 0)) then
      negative = (alpha * real(nn - 1, dp) + beta * real(mm - 1, dp))
    else
      negative = (alpha * real(nn - 1, dp) + beta * real(mm - 1, dp)) * rr**(mm + nn - 3)
    end if

    basis_1st_times_basis_1st = normalization1 * normalization2&
        & * (positive - negative) * exp(ab * rr)

    if (abs(basis_1st_times_basis_1st) < 1.0d-20) basis_1st_times_basis_1st = 0.0_dp

  end function basis_1st_times_basis_1st


  !> Evaluates product of two 2nd derivatives of basis functions.
  !! beta and poly2 are the arguments of the 2nd basis.
  pure function basis_2nd_times_basis_2nd(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of 1st basis derivative
    real(dp), intent(in) :: alpha

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of 2nd derivatives of basis functions
    real(dp) :: basis_2nd_times_basis_2nd

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab, positive, negative

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !
    positive = real((mm - 1) * (mm - 2) * (nn - 1) * (nn - 2), dp) * rr**(nn + mm - 6)&
        & + rr**(mm + nn - 4) * (beta**2 * real((mm - 1) * (mm - 2), dp) + alpha**2&
        & * real((nn - 1) * (nn - 2), dp) + alpha * beta * real(4 * (mm - 1) * (nn - 1), dp))&
        & + alpha**2 * beta**2 * rr**(mm + nn - 2)

    negative = rr**(mm + nn - 5) * (beta * real(2 * (nn - 1) * (mm - 1) * (mm - 2), dp)&
        & + alpha * real(2 * (mm - 1) * (nn - 1) * (nn - 2), dp)) + rr**(mm + nn - 3)&
        & * (alpha * beta**2 * real(2 * (mm - 1), dp) + beta * alpha**2 * real(2 * (nn - 1), dp))

    basis_2nd_times_basis_2nd = normalization1 * normalization2&
        & * (positive - negative) * exp(ab * rr)

    if (abs(basis_2nd_times_basis_2nd) < 1.0d-20) basis_2nd_times_basis_2nd = 0.0_dp

  end function basis_2nd_times_basis_2nd


  !> Evaluates product of two basis functions and r^2.
  !! r^(m-1)*e^(-alpha*r)*r^(n-1)*exp(-beta*r)*r^2=r^(m+n)*exp(-(alpha+beta)*r)
  pure function basis_times_basis_times_r2(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of 1st basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis
    real(dp), intent(in) :: beta

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of two basis functions and r^2
    real(dp) :: basis_times_basis_times_r2

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    basis_times_basis_times_r2 = normalization1 * normalization2 * rr**(mm + nn) * exp(ab * rr)

    if (abs(basis_times_basis_times_r2) < 1.0d-20) basis_times_basis_times_r2 = 0.0_dp

  end function basis_times_basis_times_r2


  !> Evaluates product of a basis function with 1st derivative of another basis function and r^2.
  !! beta and poly2 are the arguments of the derivative.
  pure function basis_times_basis_1st_times_r2(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in basis shell
    integer, intent(in) :: poly1

    !> basis exponent of basis derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in basis derivative shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of a basis function with 1st derivative of another basis function and r^2
    real(dp) :: basis_times_basis_1st_times_r2

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !
    basis_times_basis_1st_times_r2 = normalization1 * normalization2&
        & * (real(nn - 1, dp) * rr**(mm + nn - 1) - beta * rr**(nn + mm)) * exp(ab * rr)

    if (abs(basis_times_basis_1st_times_r2) < 1.0d-20) basis_times_basis_1st_times_r2 = 0.0_dp

  end function basis_times_basis_1st_times_r2


  !> Evaluates product of a basis function with 2nd derivative of another basis function and r^2.
  !! beta and poly2 are the arguments of the derivative.
  pure function basis_times_basis_2nd_times_r2(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in basis shell
    integer, intent(in) :: poly1

    !> basis exponent of basis 2nd derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in basis 2nd derivative shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of a basis function with 2nd derivative of another basis function and r^2
    real(dp) :: basis_times_basis_2nd_times_r2

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab, positive, negative

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !
    positive = real((nn - 1) * (nn - 2), dp) * rr**(mm + nn - 2) + beta**2 * rr**(mm + nn)
    negative = real(2 * (nn - 1), dp) * beta * rr**(nn + mm - 1)

    basis_times_basis_2nd_times_r2 = normalization1 * normalization2&
        & * (positive - negative) * exp(ab * rr)

    if (abs(basis_times_basis_2nd_times_r2) < 1.0d-20) basis_times_basis_2nd_times_r2 = 0.0_dp

  end function basis_times_basis_2nd_times_r2


  !> Evaluates product of a basis function with 1st derivative of another and r.
  !! beta and poly2 are the arguments of the derivative.
  pure function basis_times_basis_1st_times_r(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of basis
    real(dp), intent(in) :: alpha

    !> highest polynomial order in basis shell
    integer, intent(in) :: poly1

    !> basis exponent of basis derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in basis derivative shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of a basis function with 1st derivative of another and r
    real(dp) :: basis_times_basis_1st_times_r

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !
    basis_times_basis_1st_times_r = normalization1 * normalization2&
        & * (real(nn - 1, dp) * rr**(mm + nn - 2) - beta * rr**(nn + mm - 1)) * exp(ab * rr)

    if (abs(basis_times_basis_1st_times_r) < 1.0d-20) basis_times_basis_1st_times_r = 0.0_dp

  end function basis_times_basis_1st_times_r


  !> Evaluates product of 1st derivatives of basis functions and r^2.
  !! beta and poly2 are the arguments of the 2nd basis.
  pure function basis_1st_times_basis_1st_times_r2(alpha, poly1, beta, poly2, ll, rr)

    !> basis exponent of 1st basis derivative
    real(dp), intent(in) :: alpha

    !> highest polynomial order in 1st basis shell
    integer, intent(in) :: poly1

    !> basis exponent of 2nd basis derivative
    real(dp), intent(in) :: beta

    !> highest polynomial order in 2nd basis shell
    integer, intent(in) :: poly2

    !> angular momentum
    integer, intent(in) :: ll

    !> radial point in space, i.e. abcissa
    real(dp), intent(in) :: rr

    !> product of 1st derivatives of basis functions and r^2
    real(dp) :: basis_1st_times_basis_1st_times_r2

    !> normalization pre-factors
    real(dp) :: normalization1, normalization2

    !> auxiliary variables
    integer :: mm, nn
    real(dp) :: ab, positive, negative

    mm = poly1 + ll
    nn = poly2 + ll
    ab = - (alpha + beta)

    normalization1 = (2.0_dp * alpha)**(mm) * sqrt(2.0_dp * alpha) / sqrt(fak(2 * mm))
    normalization2 = (2.0_dp * beta)**(nn) * sqrt(2.0_dp * beta) / sqrt(fak(2 * nn))

    ! WARNING: without summing negative and positive contributions independently,
    ! zora becomes completely unstable !
    positive = real((mm - 1) * (nn - 1), dp) * rr**(mm + nn - 2) + alpha * beta * rr**(mm + nn)
    negative = (alpha * real(nn - 1, dp) + beta * real(mm - 1, dp)) * rr**(mm + nn - 1)

    basis_1st_times_basis_1st_times_r2 = normalization1 * normalization2&
        & * (positive - negative) * exp(ab * rr)

    if (abs(basis_1st_times_basis_1st_times_r2) < 1.0d-20)&
        & basis_1st_times_basis_1st_times_r2 = 0.0_dp

  end function basis_1st_times_basis_1st_times_r2

end module density
