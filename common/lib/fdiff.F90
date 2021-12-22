!******************************************************************************************
! general convention:
!         module m_name
!         all methods should have prefix: m_name_method
!******************************************************************************************

module common_finitedifferences

  use common_accuracy, only : dp
  use common_constants, only : pi, pi_hlf

  implicit none
  private

  public :: makeFDMatrix, makeFDMatrix7P, H_BFDM7P, P_BFDM7P


contains

!*********************************************************************
! Generates finite differences matrix and rhs: Mx=B
! alpha(z)f'' + beta(z)f' + gamma(z)f = B(rho)
! M: FD matrix
! B: rhs
! Mx=B
!*********************************************************************

subroutine makeFDMatrix(M, B, N, alpha, beta, gama, rho, f0, fNp1)
implicit none
real(dp), intent(inout) :: M(:,:)
real(dp), intent(inout) :: B(:)
integer, intent(in) :: N
real(dp), intent(in) :: alpha(:)
real(dp), intent(in) :: beta(:)
real(dp), intent(in) :: gama(:)
real(dp), intent(in) :: rho(:)
real(dp), intent(in) :: f0
real(dp), intent(in) :: fNp1

integer :: ii
real(dp) :: step, step_2, zi, zi_q

step = 1.0_dp/real(N+1, dp)
step_2 = step*step

!write (*,*) "# 3-point FD, step: ", step, " N=", N


 M(1,2) = alpha(1) + beta(1)*0.5_dp*step
 M(1,1) = alpha(1)*(-2.0_dp) + gama(1)*step_2

 M(N,N-1) = alpha(N)  - beta(N)*0.5_dp*(step)
 M(N,N) = alpha(N)*(-2.0_dp) + gama(N)*step_2

 do ii=2, N-1

  M(ii,ii) =   alpha(ii)*(-2.0_dp) + gama(ii)*step_2
  M(ii,ii-1) = alpha(ii)  - beta(ii)*0.5_dp*(step)
  M(ii,ii+1) = alpha(ii)  + beta(ii)*0.5_dp*(step)

 end do

!*****************************************************
! right hand side
!*****************************************************
 do ii=1, N
  B(ii) = B(ii)*step_2
 end do

!*****************************************************
! boundary conditions
!*****************************************************
B(N)=B(N)-fNp1*(alpha(N) + 0.5_dp*step*beta(N))
B(1)=B(1)-f0*(alpha(1) - 0.5_dp*step*beta(1))

end subroutine makeFDMatrix

!*****************************************************
! seven point finite differences scheme
! ref: Bickley
! on boundaries non-centered five and six point
! formulae are used
!*****************************************************

 subroutine makeFDMatrix7P(M, B, N, alpha, beta, gama, rho, f0, fNp1)
 implicit none
 real(dp), intent(inout) :: M(:,:)
 real(dp), intent(inout) :: B(:)
 integer, intent(in) :: N
 real(dp), intent(in) :: alpha(:)
 real(dp), intent(in) :: beta(:)
 real(dp), intent(in) :: gama(:)
 real(dp), intent(in) :: rho(:)
 real(dp), intent(in) :: f0
 real(dp), intent(in) :: fNp1

 integer :: ii
 real(dp) :: step, step_2

 step = 1.0_dp/real(N+1, dp)
 step_2 = step*step

!write (*,*) "# 7-point FD, step: ", step

 !========================================================
 ! ToDo: do an automatic routine for filling the FD
 !       coefficients instead of hard coding (optional)
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

 do ii=1, N
 B(ii) = B(ii)*step_2*180.0_dp
 end do

 !==========================================================================================
 ! Boundary conditions
 !==========================================================================================
 B(N)=B(N)-fNp1*(45.0_dp*beta(N)*step + 165.0_dp*alpha(N))
 B(N-1)=B(N-1)-fNp1*(- 9.0_dp*beta(N-1)*step - 15.0_dp*alpha(N-1))
 B(N-2)=B(N-2)-fNp1*(  3.0_dp*beta(N-2)*step + 2.0_dp*alpha(N-2))

 B(1)=B(1)-f0*(-45.0_dp*beta(1)*step + 165.0_dp*alpha(1))
 B(2)=B(2)-f0*(9.0_dp*beta(2)*step  - 15.0_dp*alpha(2))
 B(3)=B(3)-f0*(-3.0_dp*beta(3)*step + 2.0_dp*alpha(3))

 !==========================================================================================
 end subroutine makeFDMatrix7P

 !
 subroutine makeBFDM7P(H2, B, N, alpha, beta, gama, rho, f0, fNp1)
   implicit none
   real(dp), intent(inout) :: H2(:,:)
   real(dp), intent(inout) :: B(:)
   integer, intent(in) :: N
   real(dp), intent(in) :: alpha(:)
   real(dp), intent(in) :: beta(:)
   real(dp), intent(in) :: gama(:)
   real(dp), intent(in) :: rho(:)
   real(dp), intent(in) :: f0
   real(dp), intent(in) :: fNp1

   integer :: ii, jj
   real(dp) :: step, step_2, tmp1

   step = 1.0_dp/real(N+1, dp)
   step_2 = step*step

   !write (*,*) "# 7-point FD, step: ", step


 ! !========================================================
 ! ! ToDo: do an automatic routine for filling the FD
 ! !       coefficients instead of hard coding (optional)
 ! !========================================================

 ! do ii=4, N-3
 ! M(ii,ii) =   alpha(ii)*(-490.0_dp)    + gama(ii)*step_2*180.0_dp
 ! M(ii,ii-1) = alpha(ii)*(270.0_dp)     + beta(ii)*(step)*(-135.0_dp)
 ! M(ii,ii-2) = alpha(ii)*(-27.0_dp)     + beta(ii)*(step)*(27.0_dp)
 ! M(ii,ii-3) = alpha(ii)*(2.0_dp)       + beta(ii)*(step)*(-3.0_dp)
 ! M(ii,ii+1) = alpha(ii)*(270.0_dp)     + beta(ii)*(step)*(135.0_dp)
 ! M(ii,ii+2) = alpha(ii)*(-27.0_dp)     + beta(ii)*(step)*(-27.0_dp)
 ! M(ii,ii+3) = alpha(ii)*(2.0_dp)       + beta(ii)*(step)*(3.0_dp)
 ! end do

!========
  !  ii=1
    H2(7, 1) = alpha(1)*(-300.0_dp)     + gama(1)*step_2*180.0_dp + (beta(1)*step*(-150.0_dp))
    H2(6, 2) = alpha(1)*(90.0_dp)       + beta(1)*step*270.0_dp
    H2(5, 3) = alpha(1)*(60.0_dp)       + beta(1)*step*(-90.0_dp)
    H2(4, 4) = alpha(1)*(-15.0_dp)      + beta(1)*step*15.0_dp

  !  ii=2
    H2(7, 2) = alpha(2)*(-450.0_dp)     + beta(2)*step*(-60.0_dp)
    H2(6, 3) = alpha(2)*( 240.0_dp)     + beta(2)*step*(180.0_dp)
    H2(5, 4) = alpha(2)*( -15.0_dp)     + beta(2)*step*(-45.0_dp)
    H2(4, 5) =                          + beta(2)*step*(6.0_dp)
    H2(8, 1) = alpha(2)*( 240.0_dp)     + gama(2)*step_2*180.0_dp + (beta(2)*step*(-90.0_dp))

    ii=3
    H2(7, 3) = alpha(3)*(-490.0_dp)   + gama(3)*step_2*180.0_dp
    H2(6, 4) = alpha(3)*(270.0_dp)    + beta(3)*(step)*(135.0_dp)
    H2(5, 5) = alpha(3)*(-27.0_dp)    + beta(3)*(step)*(-27.0_dp)
    H2(4, 6) = alpha(3)*(2.0_dp)      + beta(3)*(step)*3.0_dp
    H2(8, 2) = alpha(3)*(270.0_dp)    + beta(3)*(step)*(-135.0_dp)
    H2(9, 1) = alpha(3)*(-27.0_dp)    + beta(3)*(step)*(27.0_dp)

    tmp1=step_2*180.0_dp
    B(1) = B(1)*tmp1
    B(2) = B(2)*tmp1
    B(3) = B(3)*tmp1
    B(N-2) = B(N-2)*tmp1
    B(N-1) = B(N-1)*tmp1
    B(N) = B(N)*tmp1

    do ii=4, N-3
       H2(7,ii) = alpha(ii)*(-490.0_dp)    + gama(ii)*step_2*180.0_dp !H(ii,ii)
       H2(6, ii+1) = alpha(ii)*(270.0_dp)     + beta(ii)*(step)*(135.0_dp)!H(ii,ii+1)
       H2(5, ii+2) = alpha(ii)*(-27.0_dp)     + beta(ii)*(step)*(-27.0_dp)!H(ii,ii+2)
       H2(4, ii+3) = alpha(ii)*(2.0_dp)       + beta(ii)*(step)*(3.0_dp)!H(ii,ii+3)
       H2(8, ii-1) = alpha(ii)*(270.0_dp)     + beta(ii)*(step)*(-135.0_dp)!H(ii,ii-1)
       H2(9, ii-2) = alpha(ii)*(-27.0_dp)     + beta(ii)*(step)*(27.0_dp)!H(ii,ii-2)
       H2(10, ii-3) =  alpha(ii)*(2.0_dp)       + beta(ii)*(step)*(-3.0_dp)!H(ii,ii-3)
       B(ii) = B(ii)*tmp1
    end do


!    ii=N-2
    H2(7, N-2) = alpha(N-2)*(-490.0_dp)   + gama(N-2)*step_2*180.0_dp ! H(N-2,N-2)
    H2(6, N-1) = alpha(N-2)*(270.0_dp)    + beta(N-2)*(step)*(135.0_dp)! H(N-2,N-1)
    H2(5, N) =  alpha(N-2)*(-27.0_dp)    + beta(N-2)*(step)*(-27.0_dp)! H(N-2, N)
    H2(8, N-3) =alpha(N-2)*(270.0_dp)    + beta(N-2)*(step)*(-135.0_dp)! H(N-2,N-3)
    H2(9, N-4) =alpha(N-2)*(-27.0_dp)    + beta(N-2)*(step)*(27.0_dp)! H(N-2,N-4)
    H2(10, N-5) =alpha(N-2)*(2.0_dp)      + beta(N-2)*(step)*(-3.0_dp)! H(N-2,N-5)

!    ii=N-1
    H2(7,N-1) = alpha(N-1)*(-450.0_dp)   + beta(N-1)*(step)*(60.0_dp)!H(N-1,N-1)
    H2(6, N) =  alpha(N-1)*(240.0_dp)    + gama(N-1)*step_2*180.0_dp + beta(N-1)*step*(90.0_dp)!H(N-1,N)
    H2(8, N-2) = alpha(N-1)*(240.0_dp)    + beta(N-1)*(step)*(-180.0_dp)! H(N-1,N-2)
    H2(9, N-3) =alpha(N-1)*(-15.0_dp)    + beta(N-1)*(step)*(45.0_dp)! H(N-1,N-3)
    H2(10, N-4) =                         + beta(N-1)*(step)*(-6.0_dp)! H(N-1,N-4)

 !   ii=N
    H2(7,N) =   alpha(N)*(-300.0_dp)    + gama(N)*step_2*180.0_dp + beta(N)*step*150.0_dp!H(N,N)
    H2(8, N-1) = alpha(N)*(90.0_dp)      + beta(N)*(step)*(-270.0_dp)!H(N,N-1)
    H2(9, N-2) = alpha(N)*(60.0_dp)      + beta(N)*(step)*(90.0_dp)!H(N,N-2)
    H2(10, N-3) = alpha(N)*(-15.0_dp)     + beta(N)*(step)*(-15.0_dp)!H(N,N-3)


   ! do ii=1, N
   !    B(ii) = B(ii)*step_2*180.0_dp
   ! end do

   !==========================================================================================
   ! Boundary conditions
   !==========================================================================================
   B(N)=B(N)-fNp1*(45.0_dp*beta(N)*step + 165.0_dp*alpha(N))
   B(N-1)=B(N-1)-fNp1*(- 9.0_dp*beta(N-1)*step - 15.0_dp*alpha(N-1))
   B(N-2)=B(N-2)-fNp1*(  3.0_dp*beta(N-2)*step + 2.0_dp*alpha(N-2))

   B(1)=B(1)-f0*(-45.0_dp*beta(1)*step + 165.0_dp*alpha(1))
   B(2)=B(2)-f0*(9.0_dp*beta(2)*step  - 15.0_dp*alpha(2))
   B(3)=B(3)-f0*(-3.0_dp*beta(3)*step + 2.0_dp*alpha(3))

   !==========================================================================================

 end subroutine makeBFDM7P


 subroutine P_BFDM7P(H2, B, N, ll, zi, rm, charge)
   implicit none
   real(dp), intent(inout) :: H2(:,:)
   real(dp), intent(inout) :: B(:)
   integer, intent(in) :: N
   integer, intent(in) :: ll
   real(dp), intent(in) :: zi(:)
   real(dp), intent(in) :: rm
   real(dp), intent(in) :: charge


   integer :: ii
   real(dp) :: step, step_2, tmp1, llp1_pi_2_rm_4, f0
   real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf

   step = 1.0_dp/real(N+1, dp)
   step_2 = step*step

   if (ll == 0) then
      f0 = charge*sqrt(4.0_dp*pi) ! r->oo
   else
      f0 = 0.0_dp ! r->oo
   end if

   !   charge = step_2*180.0_dp

   !pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
   llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)

   !  ii=1
   sin_pi = sin(pi * zi(1))
   sin_pi_hlf = sin(pi_hlf * zi(1))
   cos_pi_hlf = cos(pi_hlf * zi(1))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(7, 1) = (-300.0_dp)     + gama*step_2*180.0_dp + (beta*(-150.0_dp))
   H2(6, 2) = (90.0_dp)       + beta*270.0_dp
   H2(5, 3) = (60.0_dp)       + beta*(-90.0_dp)
   H2(4, 4) = (-15.0_dp)      + beta*15.0_dp
   B(1)=f0*(45.0_dp*beta - 165.0_dp)

   !  ii=2
   sin_pi = sin(pi * zi(2))
   sin_pi_hlf = sin(pi_hlf * zi(2))
   cos_pi_hlf = cos(pi_hlf * zi(2))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(7, 2) = (-450.0_dp)     + beta*(-60.0_dp)
   H2(6, 3) = ( 240.0_dp)     + beta*(180.0_dp)
   H2(5, 4) = ( -15.0_dp)     + beta*(-45.0_dp)
   H2(4, 5) =                 + beta*(6.0_dp)
   H2(8, 1) = (240.0_dp)     + gama*step_2*180.0_dp + (beta*(-90.0_dp))
   B(2)=-f0*(9.0_dp*beta - 15.0_dp)

   !   ii=3
   sin_pi = sin(pi * zi(3))
   sin_pi_hlf = sin(pi_hlf * zi(3))
   cos_pi_hlf = cos(pi_hlf * zi(3))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(7, 3) = (-490.0_dp)   + gama*step_2*180.0_dp
   H2(6, 4) = (270.0_dp)    + beta*(135.0_dp)
   H2(5, 5) = (-27.0_dp)    + beta*(-27.0_dp)
   H2(4, 6) = (2.0_dp)      + beta*3.0_dp
   H2(8, 2) = (270.0_dp)    + beta*(-135.0_dp)
   H2(9, 1) = (-27.0_dp)    + beta*(27.0_dp)
   B(3)=f0*(3.0_dp*beta - 2.0_dp)

   do ii=4, N-3
      sin_pi = sin(pi * zi(ii))
      sin_pi_hlf = sin(pi_hlf * zi(ii))
      cos_pi_hlf = cos(pi_hlf * zi(ii))
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
      beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
      sin_pi = sin_pi*sin_pi ! ^2
      gama = -llp1_pi_2_rm_4/sin_pi

      H2(7,ii) = (-490.0_dp)    + gama*step_2*180.0_dp !H(ii,ii)
      H2(6, ii+1) = (270.0_dp)     + beta*(135.0_dp)!H(ii,ii+1)
      H2(5, ii+2) = (-27.0_dp)     + beta*(-27.0_dp)!H(ii,ii+2)
      H2(4, ii+3) = (2.0_dp)       + beta*(3.0_dp)!H(ii,ii+3)
      H2(8, ii-1) = (270.0_dp)     + beta*(-135.0_dp)!H(ii,ii-1)
      H2(9, ii-2) = (-27.0_dp)     + beta*(27.0_dp)!H(ii,ii-2)
      H2(10, ii-3) =  (2.0_dp)     + beta*(-3.0_dp)!H(ii,ii-3)
   end do

   !    ii=N-2
   sin_pi = sin(pi * zi(N-2))
   sin_pi_hlf = sin(pi_hlf * zi(N-2))
   cos_pi_hlf = cos(pi_hlf * zi(N-2))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(7, N-2) = (-490.0_dp)  + gama*step_2*180.0_dp ! H(N-2,N-2)
   H2(6, N-1) = (270.0_dp)   + beta*(135.0_dp)! H(N-2,N-1)
   H2(5, N) =  (-27.0_dp)    + beta*(-27.0_dp)! H(N-2, N)
   H2(8, N-3) =(270.0_dp)    + beta*(-135.0_dp)! H(N-2,N-3)
   H2(9, N-4) =(-27.0_dp)    + beta*(27.0_dp)! H(N-2,N-4)
   H2(10, N-5) =(2.0_dp)     + beta*(-3.0_dp)! H(N-2,N-5)

   !    ii=N-1
   sin_pi = sin(pi * zi(N-1))
   sin_pi_hlf = sin(pi_hlf * zi(N-1))
   cos_pi_hlf = cos(pi_hlf * zi(N-1))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(7,N-1) = (-450.0_dp)   + beta*(60.0_dp)!H(N-1,N-1)
   H2(6, N) =  (240.0_dp)    + gama*step_2*180.0_dp + beta*(90.0_dp)!H(N-1,N)
   H2(8, N-2) =(240.0_dp)    + beta*(-180.0_dp)! H(N-1,N-2)
   H2(9, N-3) =(-15.0_dp)    + beta*(45.0_dp)! H(N-1,N-3)
   H2(10, N-4) =             + beta*(-6.0_dp)! H(N-1,N-4)

   !   ii=N
   sin_pi = sin(pi * zi(N))
   sin_pi_hlf = sin(pi_hlf * zi(N))
   cos_pi_hlf = cos(pi_hlf * zi(N))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(7,N) =   (-300.0_dp)    + gama*step_2*180.0_dp + beta*150.0_dp!H(N,N)
   H2(8, N-1) =(90.0_dp)      + beta*(-270.0_dp)!H(N,N-1)
   H2(9, N-2) =(60.0_dp)      + beta*(90.0_dp)!H(N,N-2)
   H2(10, N-3) = (-15.0_dp)   + beta*(-15.0_dp)!H(N,N-3)
 end subroutine P_BFDM7P





 !
 subroutine H_BFDM7P(H2, B, N, ll, kappa, zi, rm)
   implicit none
   real(dp), intent(inout) :: H2(:,:)
   real(dp), intent(inout) :: B(:)
   integer, intent(in) :: N
   integer, intent(in) :: ll
   real(dp), intent(in) :: kappa
   real(dp), intent(in) :: zi(:)
   real(dp), intent(in) :: rm
   !
   integer :: ii
   real(dp) :: step, step_2, tmp1, llp1_pi_2_rm_4, f0
   real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf,a2c

   step = 1.0_dp/real(N+1, dp)
   step_2 = step*step

   f0=0.0_dp
   llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)
   a2c = kappa*kappa*rm*rm*pi*pi*0.25_dp

   !  ii=1
   sin_pi = sin(pi * zi(1))
   sin_pi_hlf = sin(pi_hlf * zi(1))
   cos_pi_hlf = cos(pi_hlf * zi(1))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   !
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

   H2(7, 1) = (-300.0_dp)     + gama*step_2*180.0_dp + (beta*(-150.0_dp))
   H2(6, 2) = (90.0_dp)       + beta*270.0_dp
   H2(5, 3) = (60.0_dp)       + beta*(-90.0_dp)
   H2(4, 4) = (-15.0_dp)      + beta*15.0_dp
!   B(1)=f0*(45.0_dp*beta - 165.0_dp)

   !  ii=2
   sin_pi = sin(pi * zi(2))
   sin_pi_hlf = sin(pi_hlf * zi(2))
   cos_pi_hlf = cos(pi_hlf * zi(2))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
!   gama = -llp1_pi_2_rm_4/sin_pi
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

   H2(7, 2) = (-450.0_dp)     + beta*(-60.0_dp)
   H2(6, 3) = ( 240.0_dp)     + beta*(180.0_dp)
   H2(5, 4) = ( -15.0_dp)     + beta*(-45.0_dp)
   H2(4, 5) =                 + beta*(6.0_dp)
   H2(8, 1) = (240.0_dp)     + gama*step_2*180.0_dp + (beta*(-90.0_dp))
!   B(2)=-f0*(9.0_dp*beta - 15.0_dp)

   !   ii=3
   sin_pi = sin(pi * zi(3))
   sin_pi_hlf = sin(pi_hlf * zi(3))
   cos_pi_hlf = cos(pi_hlf * zi(3))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
!   gama = -llp1_pi_2_rm_4/sin_pi
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

   H2(7, 3) = (-490.0_dp)   + gama*step_2*180.0_dp
   H2(6, 4) = (270.0_dp)    + beta*(135.0_dp)
   H2(5, 5) = (-27.0_dp)    + beta*(-27.0_dp)
   H2(4, 6) = (2.0_dp)      + beta*3.0_dp
   H2(8, 2) = (270.0_dp)    + beta*(-135.0_dp)
   H2(9, 1) = (-27.0_dp)    + beta*(27.0_dp)
!   B(3)=f0*(3.0_dp*beta - 2.0_dp)

   do ii=4, N-3
      sin_pi = sin(pi * zi(ii))
      sin_pi_hlf = sin(pi_hlf * zi(ii))
      cos_pi_hlf = cos(pi_hlf * zi(ii))
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
      beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
      sin_pi = sin_pi*sin_pi ! ^2
 !     gama = -llp1_pi_2_rm_4/sin_pi
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

      H2(7,ii) = (-490.0_dp)    + gama*step_2*180.0_dp !H(ii,ii)
      H2(6, ii+1) = (270.0_dp)     + beta*(135.0_dp)!H(ii,ii+1)
      H2(5, ii+2) = (-27.0_dp)     + beta*(-27.0_dp)!H(ii,ii+2)
      H2(4, ii+3) = (2.0_dp)       + beta*(3.0_dp)!H(ii,ii+3)
      H2(8, ii-1) = (270.0_dp)     + beta*(-135.0_dp)!H(ii,ii-1)
      H2(9, ii-2) = (-27.0_dp)     + beta*(27.0_dp)!H(ii,ii-2)
      H2(10, ii-3) =  (2.0_dp)     + beta*(-3.0_dp)!H(ii,ii-3)
   end do

   !    ii=N-2
   sin_pi = sin(pi * zi(N-2))
   sin_pi_hlf = sin(pi_hlf * zi(N-2))
   cos_pi_hlf = cos(pi_hlf * zi(N-2))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
 !  gama = -llp1_pi_2_rm_4/sin_pi
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

   H2(7, N-2) = (-490.0_dp)  + gama*step_2*180.0_dp ! H(N-2,N-2)
   H2(6, N-1) = (270.0_dp)   + beta*(135.0_dp)! H(N-2,N-1)
   H2(5, N) =  (-27.0_dp)    + beta*(-27.0_dp)! H(N-2, N)
   H2(8, N-3) =(270.0_dp)    + beta*(-135.0_dp)! H(N-2,N-3)
   H2(9, N-4) =(-27.0_dp)    + beta*(27.0_dp)! H(N-2,N-4)
   H2(10, N-5) =(2.0_dp)     + beta*(-3.0_dp)! H(N-2,N-5)

   !    ii=N-1
   sin_pi = sin(pi * zi(N-1))
   sin_pi_hlf = sin(pi_hlf * zi(N-1))
   cos_pi_hlf = cos(pi_hlf * zi(N-1))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
 !  gama = -llp1_pi_2_rm_4/sin_pi
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

   H2(7,N-1) = (-450.0_dp)   + beta*(60.0_dp)!H(N-1,N-1)
   H2(6, N) =  (240.0_dp)    + gama*step_2*180.0_dp + beta*(90.0_dp)!H(N-1,N)
   H2(8, N-2) =(240.0_dp)    + beta*(-180.0_dp)! H(N-1,N-2)
   H2(9, N-3) =(-15.0_dp)    + beta*(45.0_dp)! H(N-1,N-3)
   H2(10, N-4) =             + beta*(-6.0_dp)! H(N-1,N-4)

   !   ii=N
   sin_pi = sin(pi * zi(N))
   sin_pi_hlf = sin(pi_hlf * zi(N))
   cos_pi_hlf = cos(pi_hlf * zi(N))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
 !  gama = -llp1_pi_2_rm_4/sin_pi
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^4
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^8
   gama = -(llp1_pi_2_rm_4/sin_pi+ a2c*sin_pi/sin_pi_hlf)

   H2(7,N) =   (-300.0_dp)    + gama*step_2*180.0_dp + beta*150.0_dp!H(N,N)
   H2(8, N-1) =(90.0_dp)      + beta*(-270.0_dp)!H(N,N-1)
   H2(9, N-2) =(60.0_dp)      + beta*(90.0_dp)!H(N,N-2)
   H2(10, N-3) = (-15.0_dp)   + beta*(-15.0_dp)!H(N,N-3)
 end subroutine H_BFDM7P

 !
 subroutine TST_P_BFDM7P(H2, B, N, ll, zi, rm, charge)
   implicit none
   real(dp), intent(inout) :: H2(:,:)
   real(dp), intent(inout) :: B(:)
   integer, intent(in) :: N
   integer, intent(in) :: ll
!   real(dp), intent(in) :: kappa
   real(dp), intent(in) :: zi(:)
   real(dp), intent(in) :: rm
   real(dp), intent(in) :: charge


   integer :: ii
   real(dp) :: step, step_2, tmp1, llp1_pi_2_rm_4, f0
   real(dp) :: beta, gama, sin_pi, sin_pi_hlf, cos_pi_hlf

   step = 1.0_dp/real(N+1, dp)
   step_2 = step*step

   if (ll == 0) then
      f0 = charge*sqrt(4.0_dp*pi) ! r->oo
   else
      f0 = 0.0_dp ! r->oo
   end if

   !   charge = step_2*180.0_dp

   !pi_rm_4_3 = pi*pi*pi*rm*rm*rm*4.0_dp
   llp1_pi_2_rm_4 = 4.0_dp*pi*pi*rm*real(ll*(ll+1), dp)

   !  ii=1
   sin_pi = sin(pi * zi(1))
   sin_pi_hlf = sin(pi_hlf * zi(1))
   cos_pi_hlf = cos(pi_hlf * zi(1))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(4, 1) = (-300.0_dp)     + gama*step_2*180.0_dp + (beta*(-150.0_dp))
   H2(3, 2) = (90.0_dp)       + beta*270.0_dp
   H2(2, 3) = (60.0_dp)       + beta*(-90.0_dp)
   H2(1, 4) = (-15.0_dp)      + beta*15.0_dp
   B(1)=f0*(45.0_dp*beta - 165.0_dp)

   !  ii=2
   sin_pi = sin(pi * zi(2))
   sin_pi_hlf = sin(pi_hlf * zi(2))
   cos_pi_hlf = cos(pi_hlf * zi(2))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(4, 2) = (-450.0_dp)     + beta*(-60.0_dp)
   H2(3, 3) = ( 240.0_dp)     + beta*(180.0_dp)
   H2(2, 4) = ( -15.0_dp)     + beta*(-45.0_dp)
   H2(1, 5) =                 + beta*(6.0_dp)
   H2(5, 1) = (240.0_dp)     + gama*step_2*180.0_dp + (beta*(-90.0_dp))
   B(2)=-f0*(9.0_dp*beta - 15.0_dp)

   !   ii=3
   sin_pi = sin(pi * zi(3))
   sin_pi_hlf = sin(pi_hlf * zi(3))
   cos_pi_hlf = cos(pi_hlf * zi(3))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(4, 3) = (-490.0_dp)   + gama*step_2*180.0_dp
   H2(3, 4) = (270.0_dp)    + beta*(135.0_dp)
   H2(2, 5) = (-27.0_dp)    + beta*(-27.0_dp)
   H2(1, 6) = (2.0_dp)      + beta*3.0_dp
   H2(5, 2) = (270.0_dp)    + beta*(-135.0_dp)
   H2(6, 1) = (-27.0_dp)    + beta*(27.0_dp)
   B(3)=f0*(3.0_dp*beta - 2.0_dp)



   do ii=4, N-3
      sin_pi = sin(pi * zi(ii))
      sin_pi_hlf = sin(pi_hlf * zi(ii))
      cos_pi_hlf = cos(pi_hlf * zi(ii))
      sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
      cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
      beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
      sin_pi = sin_pi*sin_pi ! ^2
      gama = -llp1_pi_2_rm_4/sin_pi

      H2(4,ii) = (-490.0_dp)    + gama*step_2*180.0_dp !H(ii,ii)
      H2(3, ii+1) = (270.0_dp)     + beta*(135.0_dp)!H(ii,ii+1)
      H2(2, ii+2) = (-27.0_dp)     + beta*(-27.0_dp)!H(ii,ii+2)
      H2(1, ii+3) = (2.0_dp)       + beta*(3.0_dp)!H(ii,ii+3)
      H2(5, ii-1) = (270.0_dp)     + beta*(-135.0_dp)!H(ii,ii-1)
      H2(6, ii-2) = (-27.0_dp)     + beta*(27.0_dp)!H(ii,ii-2)
      H2(7, ii-3) =  (2.0_dp)     + beta*(-3.0_dp)!H(ii,ii-3)
   end do

   !    ii=N-2
   sin_pi = sin(pi * zi(N-2))
   sin_pi_hlf = sin(pi_hlf * zi(N-2))
   cos_pi_hlf = cos(pi_hlf * zi(N-2))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(4, N-2) = (-490.0_dp)  + gama*step_2*180.0_dp ! H(N-2,N-2)
   H2(3, N-1) = (270.0_dp)   + beta*(135.0_dp)! H(N-2,N-1)
   H2(2, N) =  (-27.0_dp)    + beta*(-27.0_dp)! H(N-2, N)
   H2(5, N-3) =(270.0_dp)    + beta*(-135.0_dp)! H(N-2,N-3)
   H2(6, N-4) =(-27.0_dp)    + beta*(27.0_dp)! H(N-2,N-4)
   H2(7, N-5) =(2.0_dp)     + beta*(-3.0_dp)! H(N-2,N-5)

   !    ii=N-1
   sin_pi = sin(pi * zi(N-1))
   sin_pi_hlf = sin(pi_hlf * zi(N-1))
   cos_pi_hlf = cos(pi_hlf * zi(N-1))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(4, N-1) = (-450.0_dp)   + beta*(60.0_dp)!H(N-1,N-1)
   H2(3, N) =  (240.0_dp)    + gama*step_2*180.0_dp + beta*(90.0_dp)!H(N-1,N)
   H2(5, N-2) =(240.0_dp)    + beta*(-180.0_dp)! H(N-1,N-2)
   H2(6, N-3) =(-15.0_dp)    + beta*(45.0_dp)! H(N-1,N-3)
   H2(7, N-4) =             + beta*(-6.0_dp)! H(N-1,N-4)

   !   ii=N
   sin_pi = sin(pi * zi(N))
   sin_pi_hlf = sin(pi_hlf * zi(N))
   cos_pi_hlf = cos(pi_hlf * zi(N))
   sin_pi_hlf = sin_pi_hlf*sin_pi_hlf ! ^2
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf !
   cos_pi_hlf = cos_pi_hlf*cos_pi_hlf ! ^4
   beta = (pi/sin_pi + pi_hlf*sin_pi/sin_pi_hlf)*step
   sin_pi = sin_pi*sin_pi ! ^2
   gama = -llp1_pi_2_rm_4/sin_pi

   H2(4,N) =   (-300.0_dp)    + gama*step_2*180.0_dp + beta*150.0_dp!H(N,N)
   H2(5, N-1) =(90.0_dp)      + beta*(-270.0_dp)!H(N,N-1)
   H2(6, N-2) =(60.0_dp)      + beta*(90.0_dp)!H(N,N-2)
   H2(7, N-3) = (-15.0_dp)   + beta*(-15.0_dp)!H(N,N-3)


 end subroutine TST_P_BFDM7P

end module common_finitedifferences
