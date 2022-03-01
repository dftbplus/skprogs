!> Module that generates finite differences matrix and rhs: Mx=B
!! alpha(z)f'' + beta(z)f' + gamma(z)f = B(rho)
!! M: FD matrix
!! B: rhs
!! Mx = B
module common_finitedifferences

  use common_accuracy, only : dp
  use common_constants, only : pi, pi_hlf

  implicit none
  private

  public :: makeHelmholzFDMatrix7P, makePoissonFDMatrix7P


contains

  !> Seven-point finite differences scheme (ref: Bickley).
  !! On boundaries non-centered five and six point formulae are used.
  subroutine makeHelmholzFDMatrix7P(M, B, alpha, beta, gama, f0, f1)

    !> finite differences matrix, 7-point scheme
    real(dp), intent(out), allocatable :: M(:,:)

    !> right-hand side of equation
    real(dp), intent(inout) :: B(:)

    !>
    real(dp), intent(in) :: alpha(:)

    !>
    real(dp), intent(in) :: beta(:)

    !>
    real(dp), intent(in) :: gama(:)

    !> boundary condition, r --> oo
    real(dp), intent(in) :: f0

    !> boundary condition, r --> 0
    real(dp), intent(in) :: f1

    !! number of grid points
    integer :: N

    !! iterates over grid points
    integer :: ii

    !! stepwidth and stepwidth**2
    real(dp) :: step, step_2

    N = size(B)

    allocate(M(N, N))
    M(:,:) = 0.0_dp

    step = 1.0_dp / real(N + 1, dp)
    step_2 = step**2

    !========================================================
    M(1,4) =  alpha(1)*(-15.0_dp)      + beta(1)*step*15.0_dp
    M(1,3) =  alpha(1)*(60.0_dp)       + beta(1)*step*(-90.0_dp)
    M(1,2) =  alpha(1)*(90.0_dp)       + beta(1)*step*270.0_dp
    M(1,1) =  alpha(1)*(-300.0_dp)     + gama(1)*step_2*180.0_dp + (beta(1)*step*(-150.0_dp))
    !================================
    M(N,N-1) = alpha(N)*(90.0_dp)      + beta(N)*(step)*(-270.0_dp)
    M(N,N-2) = alpha(N)*(60.0_dp)      + beta(N)*(step)*(90.0_dp)
    M(N,N-3) = alpha(N)*(-15.0_dp)     + beta(N)*(step)*(-15.0_dp)
    M(N,N) =   alpha(N)*(-300.0_dp)    + gama(N)*step_2*180.0_dp + beta(N)*step*150.0_dp
    !================================
    M(2,5) =                          + beta(2)*step*(6.0_dp)
    M(2,4) = alpha(2)*( -15.0_dp)     + beta(2)*step*(-45.0_dp)
    M(2,3) = alpha(2)*( 240.0_dp)     + beta(2)*step*(180.0_dp)
    M(2,2) = alpha(2)*(-450.0_dp)     + beta(2)*step*(-60.0_dp)
    M(2,1) = alpha(2)*( 240.0_dp)     + gama(2)*step_2*180.0_dp + (beta(2)*step*(-90.0_dp))
    !================================
    M(N-1,N-1) = alpha(N-1)*(-450.0_dp)   + beta(N-1)*(step)*(60.0_dp)
    M(N-1,N-2) = alpha(N-1)*(240.0_dp)    + beta(N-1)*(step)*(-180.0_dp)
    M(N-1,N-3) = alpha(N-1)*(-15.0_dp)    + beta(N-1)*(step)*(45.0_dp)
    M(N-1,N-4) =                          + beta(N-1)*(step)*(-6.0_dp)
    M(N-1,N) =   alpha(N-1)*(240.0_dp)    + gama(N-1)*step_2*180.0_dp + beta(N-1)*step*(90.0_dp)
    !================================

    M(3,3) = alpha(3)*(-490.0_dp)   + gama(3)*step_2*180.0_dp
    M(3,2) = alpha(3)*(270.0_dp)    + beta(3)*(step)*(-135.0_dp)
    M(3,1) = alpha(3)*(-27.0_dp)    + beta(3)*(step)*(27.0_dp)
    M(3,4) = alpha(3)*(270.0_dp)    + beta(3)*(step)*(135.0_dp)
    M(3,5) = alpha(3)*(-27.0_dp)    + beta(3)*(step)*(-27.0_dp)
    M(3,6) = alpha(3)*(2.0_dp)      + beta(3)*(step)*3.0_dp

    M(N-2,N-2) = alpha(N-2)*(-490.0_dp)   + gama(N-2)*step_2*180.0_dp
    M(N-2,N-3) = alpha(N-2)*(270.0_dp)    + beta(N-2)*(step)*(-135.0_dp)
    M(N-2,N-4) = alpha(N-2)*(-27.0_dp)    + beta(N-2)*(step)*(27.0_dp)
    M(N-2,N-5) = alpha(N-2)*(2.0_dp)      + beta(N-2)*(step)*(-3.0_dp)
    M(N-2,N-1) = alpha(N-2)*(270.0_dp)    + beta(N-2)*(step)*(135.0_dp)
    M(N-2,N) =   alpha(N-2)*(-27.0_dp)    + beta(N-2)*(step)*(-27.0_dp)

    do ii=4, N-3
      M(ii,ii) =   alpha(ii)*(-490.0_dp)    + gama(ii)*step_2*180.0_dp
      M(ii,ii-1) = alpha(ii)*(270.0_dp)     + beta(ii)*(step)*(-135.0_dp)
      M(ii,ii-2) = alpha(ii)*(-27.0_dp)     + beta(ii)*(step)*(27.0_dp)
      M(ii,ii-3) = alpha(ii)*(2.0_dp)       + beta(ii)*(step)*(-3.0_dp)
      M(ii,ii+1) = alpha(ii)*(270.0_dp)     + beta(ii)*(step)*(135.0_dp)
      M(ii,ii+2) = alpha(ii)*(-27.0_dp)     + beta(ii)*(step)*(-27.0_dp)
      M(ii,ii+3) = alpha(ii)*(2.0_dp)       + beta(ii)*(step)*(3.0_dp)
    end do

    B(:) = B * step_2 * 180.0_dp

    !==========================================================================================
    ! Boundary conditions
    !==========================================================================================
    B(N)=B(N)-f1*(45.0_dp*beta(N)*step + 165.0_dp*alpha(N))
    B(N-1)=B(N-1)-f1*(- 9.0_dp*beta(N-1)*step - 15.0_dp*alpha(N-1))
    B(N-2)=B(N-2)-f1*(  3.0_dp*beta(N-2)*step + 2.0_dp*alpha(N-2))

    B(1)=B(1)-f0*(-45.0_dp*beta(1)*step + 165.0_dp*alpha(1))
    B(2)=B(2)-f0*(9.0_dp*beta(2)*step  - 15.0_dp*alpha(2))
    B(3)=B(3)-f0*(-3.0_dp*beta(3)*step + 2.0_dp*alpha(3))
    !==========================================================================================

  end subroutine makeHelmholzFDMatrix7P


  !>
  subroutine makePoissonFDMatrix7P(H2, B, N, ll, zi, rm, charge)

    real(dp), intent(inout) :: H2(:,:)
    real(dp), intent(inout) :: B(:)
    integer, intent(in) :: N
    integer, intent(in) :: ll
    real(dp), intent(in) :: zi(:)
    real(dp), intent(in) :: rm
    real(dp), intent(in) :: charge

    integer :: ii
    real(dp) :: step, step_2, llp1_pi_2_rm_4, f0
    real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf

    step = 1.0_dp / real(N + 1, dp)
    step_2 = step**2

    if (ll == 0) then
      ! r --> oo
      f0 = charge * sqrt(4.0_dp * pi)
    else
      ! r->oo
      f0 = 0.0_dp
    end if

    llp1_pi_2_rm_4 = 4.0_dp * pi**2 * rm * real(ll * (ll + 1), dp)

    ! ii = 1
    sin_pi = sin(pi * zi(1))
    sin_pi_hlf = sin(pi_hlf * zi(1))
    cos_pi_hlf = cos(pi_hlf * zi(1))
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2
    gama = -llp1_pi_2_rm_4 / sin_pi

    H2(7, 1) = -300.0_dp + gama * step_2 * 180.0_dp - beta * 150.0_dp
    H2(6, 2) = 90.0_dp + beta * 270.0_dp
    H2(5, 3) = 60.0_dp - beta * 90.0_dp
    H2(4, 4) = -15.0_dp + beta * 15.0_dp
    B(1) = f0 * (45.0_dp * beta - 165.0_dp)

    ! ii = 2
    sin_pi = sin(pi * zi(2))
    sin_pi_hlf = sin(pi_hlf * zi(2))
    cos_pi_hlf = cos(pi_hlf * zi(2))
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2
    gama = -llp1_pi_2_rm_4 / sin_pi

    H2(7, 2) = -450.0_dp - beta * 60.0_dp
    H2(6, 3) = 240.0_dp + beta * 180.0_dp
    H2(5, 4) = -15.0_dp - beta * 45.0_dp
    H2(4, 5) = beta * 6.0_dp
    H2(8, 1) = 240.0_dp + gama * step_2 * 180.0_dp - beta * 90.0_dp
    B(2) = -f0 * (9.0_dp * beta - 15.0_dp)

    ! ii = 3
    sin_pi = sin(pi * zi(3))
    sin_pi_hlf = sin(pi_hlf * zi(3))
    cos_pi_hlf = cos(pi_hlf * zi(3))
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2
    gama = -llp1_pi_2_rm_4 / sin_pi

    H2(7, 3) = -490.0_dp + gama * step_2 * 180.0_dp
    H2(6, 4) = 270.0_dp + beta * 135.0_dp
    H2(5, 5) = -27.0_dp - beta * 27.0_dp
    H2(4, 6) = 2.0_dp + beta * 3.0_dp
    H2(8, 2) = 270.0_dp - beta * 135.0_dp
    H2(9, 1) = -27.0_dp + beta * 27.0_dp
    B(3) = f0 * (3.0_dp*beta - 2.0_dp)

    do ii = 4, N - 3
      sin_pi = sin(pi * zi(ii))
      sin_pi_hlf = sin(pi_hlf * zi(ii))
      cos_pi_hlf = cos(pi_hlf * zi(ii))
      sin_pi_hlf = sin_pi_hlf**2 ! ^2
      cos_pi_hlf = cos_pi_hlf**4 ! ^4
      beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
      sin_pi = sin_pi**2 ! ^2
      gama = -llp1_pi_2_rm_4 / sin_pi

      H2(7, ii) = -490.0_dp + gama * step_2 * 180.0_dp
      H2(6, ii + 1) = 270.0_dp + beta * 135.0_dp
      H2(5, ii + 2) = -27.0_dp - beta * 27.0_dp
      H2(4, ii + 3) = 2.0_dp + beta * 3.0_dp
      H2(8, ii - 1) = 270.0_dp - beta * 135.0_dp
      H2(9, ii - 2) = -27.0_dp + beta * 27.0_dp
      H2(10, ii - 3) = 2.0_dp - beta * 3.0_dp
    end do

    ! ii = N - 2
    sin_pi = sin(pi * zi(N-2))
    sin_pi_hlf = sin(pi_hlf * zi(N-2))
    cos_pi_hlf = cos(pi_hlf * zi(N-2))
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2
    gama = -llp1_pi_2_rm_4 / sin_pi

    H2(7, N - 2) = -490.0_dp + gama * step_2 * 180.0_dp
    H2(6, N - 1) = 270.0_dp + beta * 135.0_dp
    H2(5, N) = -27.0_dp - beta * 27.0_dp
    H2(8, N - 3) = 270.0_dp - beta * 135.0_dp
    H2(9, N - 4) = -27.0_dp + beta * 27.0_dp
    H2(10, N - 5) = 2.0_dp - beta * 3.0_dp

    ! ii = N - 1
    sin_pi = sin(pi * zi(N-1))
    sin_pi_hlf = sin(pi_hlf * zi(N-1))
    cos_pi_hlf = cos(pi_hlf * zi(N-1))
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2
    gama = -llp1_pi_2_rm_4 / sin_pi

    H2(7, N - 1) = -450.0_dp + beta * 60.0_dp
    H2(6, N) = 240.0_dp + gama * step_2 * 180.0_dp + beta * 90.0_dp
    H2(8, N - 2) = 240.0_dp - beta * 180.0_dp
    H2(9, N - 3) = -15.0_dp + beta * 45.0_dp
    H2(10, N - 4) = -beta * 6.0_dp

    ! ii = N
    sin_pi = sin(pi * zi(N))
    sin_pi_hlf = sin(pi_hlf * zi(N))
    cos_pi_hlf = cos(pi_hlf * zi(N))
    sin_pi_hlf = sin_pi_hlf**2 ! ^2
    cos_pi_hlf = cos_pi_hlf**4 ! ^4
    beta = (pi / sin_pi + pi_hlf * sin_pi / sin_pi_hlf) * step
    sin_pi = sin_pi**2 ! ^2
    gama = -llp1_pi_2_rm_4 / sin_pi

    H2(7, N) = -300.0_dp + gama * step_2 * 180.0_dp + beta * 150.0_dp
    H2(8, N - 1) = 90.0_dp - beta * 270.0_dp
    H2(9, N - 2) = 60.0_dp + beta * 90.0_dp
    H2(10, N - 3) = -15.0_dp - beta * 15.0_dp

  end subroutine makePoissonFDMatrix7P

end module common_finitedifferences
