!> Module that provides routines for matrix diagonalization.
module diagonalizations

  use common_accuracy, only : dp

  implicit none
  private

  public :: diagonalize_overlap, diagonalize


contains

  !> Diagonalizes overlap matrix to check for linear dependency of basis set.
  !! Implicitely ewevge is called, but with a unit matrix instead of a real overlap.
  subroutine diagonalize_overlap(max_l, num_alpha, poly_order, ss)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> overlap supervector
    real(dp), intent(in) :: ss(0:, :,:)

    !> temporary storage container
    real(dp), allocatable :: temp1(:,:), temp2(:), dummy2(:,:), dummy1(:)

    !> auxiliary variables
    integer :: ii, ll, diagsize, iErr

    do ll = 0, max_l

      diagsize = num_alpha(ll) * poly_order(ll)

      allocate(temp1(diagsize, diagsize))
      allocate(temp2(diagsize))
      allocate(dummy1(diagsize))
      allocate(dummy2(diagsize, diagsize))

      temp1(:,:) = 0.0_dp
      temp2(:) = 0.0_dp
      dummy1(:) = 0.0_dp
      dummy2(:,:) = 0.0_dp

      do ii = 1, diagsize
        dummy2(ii, ii) = 1.0_dp
      end do

      temp1(:,:) = ss(ll, :,:)

      call ewevge(diagsize, diagsize, diagsize, temp1, dummy2, temp2, dummy1, 1, -1, iErr)

      if (iErr /= 0) then
        write(*,*) 'Error in Diagonalization', iErr
        stop
      end if

      write(*, '(A,I3,A,E16.8)') 'Smallest eigenvalue of overlap for l= ', ll, ' : ', temp2(1)

      if (temp2(1) < 1.0e-10_dp) then
        write(*, '(A)') ' '
        write(*, '(A)') 'Basis set is nearly linear dependent, reduction necessary'
        write(*, '(A)') ' '
        stop
      end if

      deallocate(temp1)
      deallocate(temp2)
      deallocate(dummy2)
      deallocate(dummy1)

    end do
    write(*,*) ' '

  end subroutine diagonalize_overlap


  !> This is a driver for ewevge. The idea is that the matrices are allocated in the main program
  !! for the maximum size of the problem but ewevge is only fed with a matrix of the current size of
  !! the eigenproblem.
  subroutine diagonalize(max_l, num_alpha, poly_order, ff, ss, cof_new, eigval)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> fock matrix supervector
    real(dp), intent(in) :: ff(:, 0:, :,:)

    !> overlap supervector
    real(dp), intent(in) :: ss(0:, :,:)

    !> new wavefunction coefficients
    real(dp) :: cof_new(:, 0:, :,:)

    !> eigenvalues
    real(dp) :: eigval(:, 0:, :)

    !> temporary storage container
    real(dp), allocatable :: temp1(:,:), temp2(:), dummy2(:,:), dummy1(:)

    !> auxiliary variables
    integer :: ii, jj, diagsize, iErr

    do ii = 1, 2
      do jj = 0, max_l

        diagsize = num_alpha(jj) * poly_order(jj)

        allocate(temp1(diagsize, diagsize))
        allocate(temp2(diagsize))
        allocate(dummy2(diagsize, diagsize))
        allocate(dummy1(4 * diagsize))
        temp1(:,:) = 0.0_dp
        temp2(:) = 0.0_dp
        dummy1(:) = 0.0_dp
        dummy2(:,:) = 0.0_dp

        temp1(:,:) = ff(ii, jj, :,:)
        dummy2(:,:) = ss(jj, :,:)

        call ewevge(diagsize, diagsize, diagsize, temp1, dummy2, temp2, dummy1, 1, -1, iErr)

        if (iErr /= 0) then
          write(*,*) 'Error in Diagonalization', iErr
          stop
        end if

        cof_new(ii, jj, :,:) = temp1
        eigval(ii, jj, :) = temp2

        deallocate(temp1)
        deallocate(temp2)
        deallocate(dummy2)
        deallocate(dummy1)

      end do
    end do

  end subroutine diagonalize


!
! **********************************************************************
!
!  This is a collection of subroutines designated to solve the real*8
!  general symmetric eigenvalue problem with or without eigenvectors.
!  The routines have been taken from different freeware FORTRAN
!  libraries and optimized by hand (or eye ?! ;-)). Most of the
!  optimizations have been done with respect to stride minimization
!  for the innermost loops of the subroutines. Problems with
!  bugs, roaches and other lifestock please report to
!
!  Dirk Porezag   porezag@physik.tu-chemnitz.de
!
!  or to your nearest pest control agency (I doubt they will help).
!  Have fun !!
!
!  Copyright for this file by Dirk Porezag
!  Washington, DC, Janurary 8th, 1995
!
! Modifications with some fortran90 features by ckoe
!
! **********************************************************************
!
!     SUBROUTINE EWEVGE
!     =================
!
! **********************************************************************
!
!  Evevge calculates eigenvalues and eigenvectors of the general
!  symmetric eigenvalue problem.
!
!  Method:  *  A*C = E*S*C
!           *  Choleski decomposition  S = R'*R
!           *  A*C = E*R'*R*C  ->  INV(R')*A*C = E*R*C
!           *  Transformation Y = R*C  ->  C = INV(R)*Y
!           *  Solve INV(R')*A*INV(R)*Y = E*Y  (Householder + IQL)
!           *  Back transformation C = INV(R)*Y
!           *  Sorting of eigenvalues and eigenvectors
!
!     Parameters:
!
!       NA      (I) :  Dimension of A
!       NB      (I) :  Dimension of B
!       N       (I) :  Dimension of Problem
!       A       (I) :  Matrix A (lower triangle)
!               (O) :  Eigenvector matrix
!       B       (I) :  Matrix B (lower triangle)
!               (O) :  R where B = R'*R (upper triangle)
!       EW      (O) :  Eigenvalues
!       H       (-) :  Auxiliary vector
!       IEV     (I) :  0: No eigenvectors
!       IORD    (I) :  1: Descending order of eigenvalues
!                     -1: Ascending order of eigenvalues
!                      otherwise: no sorting
!       IER     (O) :  Error indication
!                      0: No error
!                      K: (K <= N)  B is not positive definite
!                      K: (K > N) Convergence failure for eigenvalue
!                                 (K-N), (K-N-1) eigenvalues are correct
!
! **********************************************************************
!
      SUBROUTINE EWEVGE (NA,NB,N,A,B,EW,H,IEV,IORD,IER)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer, intent(in) :: NA,NB,N
        integer, intent(in) :: iev,iord
        integer :: IER,ii,i,j
        real(dp) :: a,b,ew,h,eps
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  A(NA,N),B(NB,N),EW(N),H(N)
!
!        do i=1,n
!          do j=1,n
!            write(*,*) 'we',i,j,a(i,j),b(i,j)
!          end do
!        end do
        IER = 0
        EPS = 0.0_dp
        CALL CHOLES(N,B,NB,IER)
        IF (IER .NE. 0) RETURN
        CALL MATRAF(N,A,NA,B,NB,H)
        CALL TRIDIA(NA,N,EW,H,A,IEV)
        CALL IQLDIA(NA,N,EW,H,A,IEV,IER)
        IF (IER .GT. 0) IER = IER+N
        IF (IER .NE. 0)  RETURN
        IF (IEV .NE. 0) CALL BACKTR(N,N,B,NB,A,NA,A,NA,H)
        II = 0
        IF (IEV .NE. 0) II = 1
        CALL SORTVC(NA,N,N,EW,A,IORD,II,H)
        RETURN
      END SUBROUTINE EWEVGE
!
! ******************************************************************
!
!     SUBROUTINE CHOLES
!     =================
!
! ******************************************************************
!
!  Choles calculates the Choleski decomposition B = R' * R of B
!  into an upper triangle matrix R for the symmetric positive
!  definite Matrix B. The elements of the main diagonal are
!  stored inverted.
!
!     Parameters:
!
!       N       (I) :  Dimension of problem
!       B       (I) :  Matrix B (lower triangle)
!               (O) :  Matrix R (upper triangle), inverted main diagonal
!       NB      (I) :  Dimension of B
!       ICHO    (I) :  ICHO - 1 is the dimension of the submatrix that
!                      is available as Choleski decomposition ( < 1 = 1)
!               (O) :  Row number where decomposition failed (0 if success)
!
! ******************************************************************
!
      SUBROUTINE CHOLES (N,B,NB,ICHO)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer :: N,NB,ICHO,i,ii,j,K,i1
        real(dp) :: B,d,s
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  B(NB,N)
!
        IF (ICHO .GT. N)  GOTO 200
        IF (ICHO .LT. 1)  ICHO = 1
        DO I = ICHO,N
          I1 = I - 1
          DO J = I,N
            S = B(J,I)
            DO K = 1,I1
              S = S - B(K,I) * B(K,J)
            END DO
            IF (I .NE. J) GOTO 40
            IF (S .LE. 0.0_dp) GOTO 100
            S = 1.0_dp / SQRT(S)
            D = S
            GOTO 60
   40       S = S * D
   60       B(I,J) = S
          END DO
        END DO
        ICHO = 0
        GOTO 200
  100   ICHO = I
  200   RETURN
      END SUBROUTINE CHOLES
!
! ******************************************************************
!
!     SUBROUTINE MATRAF
!     =================
!
! ******************************************************************
!
!  Matraf calculates out of the symmetric matrix A and the
!  upper triangular matrix R the product INV(R') * A * INV(R),
!  where the main diagonal of R is given inverted.
!
!     Parameters:
!
!       N       (I) :  Dimension of problem
!       A       (I) :  Matrix A (lower triangle)
!               (O) :  Transformed matrix (lower triangle)
!       NA      (I) :  Dimension of A
!       B       (I) :  Matrix R (upper triangle), inverted main diagonal
!       NB      (I) :  Dimension of B
!
! *********************************************************************
!
      SUBROUTINE MATRAF (N,A,NA,B,NB,H)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer :: N,NA,NB,i,j,ii,k,i1
        real(dp) :: A,B,H,s,d
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  A(NA,N),B(NB,N),H(N)
!
!  FILL MATRIX
!
          DO I = 1,N
            DO J = I+1,N
               A(I,J) = A(J,I)
            END DO
          END DO
!
!  CALCULATION OF A = INV(R') * A
!
          DO I = 1,N
            I1 = I-1
            D = B(I,I)
            DO J = 1,N
              S = A(I,J)
              DO K = 1,I1
                S = S - B(K,I) * A(K,J)
              END DO
              A(I,J) = S * D
            END DO
          END DO
!
!  CALCULATION OF A = A * INV(R) (USE BUFFER FOR STRIDE OPTIMIZATION)
!
          DO I = 1,N
            I1 = I-1
            D = B(I,I)
            DO J = I,N
              H(J) = A(J,I)
            END DO
            DO K = 1,I1
              S = B(K,I)
              DO J = I,N
                H(J) = H(J) - S * A(J,K)
              END DO
            END DO
            DO J = I,N
              A(J,I) = H(J) * D
            END DO
          END DO
          RETURN
        END SUBROUTINE MATRAF
!
! ******************************************************************
!
!     SUBROUTINE TRIDIA
!     =================
!
! ******************************************************************
!
!  Tridiagonalization of a given symmetric matrix A using Householder
!
!     Parameters:
!
!       NM      (I) :  Dimension of A
!       N       (I) :  Dimension of problem
!       D       (O) :  Diagonal of tridiagonal matrix
!       E       (O) :  Subdiagonal of tridiagonal matrix (E(1) = 0.0)
!       A       (I) :  Matrix A (lower triangle)
!               (O) :  Transformation Matrix
!       IEV     (I) :  0: No eigenvectors
!
! ******************************************************************
!
      SUBROUTINE TRIDIA (NM,N,D,E,A,IEV)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer :: NM,N,iev,i,j,ii,K,JP1,L
        real(dp) :: A,D,E,H,HH,G,F,scale
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  A(NM,N),D(N),E(N)
!
        DO I = 1,N
          D(I) = A(N,I)
        END DO
        IF (N .EQ. 1) GOTO 510
!
!  FOR I = N STEP -1 UNTIL 2 DO
!
        DO II = 2,N
          I = N + 2 - II
          L = I - 1
          H = 0.0_dp
          SCALE = 0.0_dp
          IF (L .LT. 2) GOTO 130
!
!  SCALE ROW
!
          DO K = 1,L
            SCALE = SCALE + ABS(D(K))
          END DO
!
          IF (SCALE .NE. 0.0_dp) GOTO 140
  130     E(I) = D(L)
          DO J = 1,L
            D(J) = A(L,J)
            A(I,J) = 0.0_dp
            A(J,I) = 0.0_dp
          END DO
          GOTO 290
!
  140     DO K = 1,L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
          END DO
          F = D(L)
          G = -SIGN(SQRT(H),F)
          E(I) = SCALE * G
          H = H - F * G
          D(L) = F - G
!
!  FORM A * U
!
          DO J = 1,L
            E(J) = 0.0_dp
          END DO
          DO J = 1,L
            F = D(J)
            A(J,I) = F
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            DO K = JP1,L
              G = G + A(K,J) * D(K)
              E(K) = E(K) + A(K,J) * F
            END DO
            E(J) = G
          END DO
!
!  FORM P
!
          F = 0.0_dp
          DO J = 1,L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
          END DO
          HH = F / (H + H)
!
!  FORM Q
!
          DO J = 1,L
            E(J) = E(J) - HH * D(J)
          END DO
!
!  FORM REDUCED A
!
          DO J = 1,L
            F = D(J)
            G = E(J)
            DO K = J,L
              A(K,J) = A(K,J) - F * E(K) - G * D(K)
            END DO
            D(J) = A(L,J)
            A(I,J) = 0.0_dp
          END DO
!
!  DONE WITH THIS TRANSFORMATION
!
  290     D(I) = H
        END DO
!
!  ACCUMULATION OF TRANSFORMATION MATRICES
!
        IF (IEV .EQ. 0) GOTO 600
        DO I = 2,N
          L = I - 1
          A(N,L) = A(L,L)
          A(L,L) = 1.0_dp
          H = D(I)
          IF (H .EQ. 0.0_dp) GOTO 380
          DO K = 1,L
            D(K) = A(K,I) / H
          END DO
          DO J = 1,L
            G = 0.0_dp
            DO K = 1,L
              G = G + A(K,I) * A(K,J)
            END DO
            DO K = 1,L
              A(K,J) = A(K,J) - G * D(K)
            END DO
          END DO
!
  380     DO K = 1,L
            A(K,I) = 0.0_dp
          END DO
        END DO
  510   DO I = 1,N
         D(I) = A(N,I)
         A(N,I) = 0.0_dp
        END DO
        GOTO 700
!
!  DEAL WITH EIGENVALUES ONLY
!
  600   DO I = 1,N
          D(I) = A(I,I)
        END DO
!
  700   A(N,N) = 1.0_dp
        E(1) = 0.0_dp
        RETURN
      END SUBROUTINE TRIDIA
!
! ******************************************************************
!
!     SUBROUTINE IQLDIA
!     =================
!
! ******************************************************************
!
!  Iqldia calculates eigenvalues and eigenvectors of a tridiagonal
!  matrix using the QL algorithm with implicit shifting.
!
!     Parameters:
!
!       NM      (I) :  Dimension of Z
!       N       (I) :  Dimension of the problem
!       D       (I) :  Diagonal of tridiagonal matrix
!               (O) :  Eigenvalues
!       E       (I) :  Subdiagonal of tridiagonal matrix
!       Z       (I) :  Transformation matrix
!               (O) :  Eigenvectors according to Z
!       IEV     (I) :  0: No eigenvectors
!       IER     (O) :  Error indication
!                      0: no error
!                      K: Convergence failure for the eigenvalue
!                         number k, k-1 eigenvalues are correct
!
! **********************************************************************
!
      SUBROUTINE IQLDIA (NM,N,D,E,Z,IEV,IER)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer :: NM,N,iev,ier,i,j,ii,k,M,L,MM1,KK,MML
        real(dp) :: E,Z,D,DD,P,G,R,S,T,PSI,PSJ,F,B,C,anorm
        real(dp) :: big,eps4,eps,epss
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  D(N),E(N),Z(NM,N)
!
        IER = 0
        IF (N .EQ. 1) RETURN
!
!  GET MACHINE EPSILON AND BIG
!
        EPS = 1.0e-2_dp
   10   IF ((1.0_dp + EPS) .EQ. 1.0_dp) GOTO 20
        EPS = 0.5_dp * EPS
        GOTO 10
   20   EPS = 2.0_dp * EPS
        EPSS = SQRT(EPS)
        EPS4 = EPS * 1.0e-4_dp
        BIG = 1.0_dp/EPS4
!
        ANORM = 0.0_dp
        R = 0.0_dp
        DO I = 2, N
          S = E(I)
          E(I-1) = S
          S = ABS(S)
          P = ABS(D(I-1)) + R + S
          IF (P .GT. ANORM) ANORM = P
          R = S
        END DO
        P = ABS(D(N)) + R
        IF (P .GT. ANORM) ANORM = P
        E(N) = 0.0_dp
        DO 250 L = 1, N
          J = 0
!
!  LOOK FOR SMALL SUBDIAGONAL ELEMENT
!
   50     DO M = L, N-1
            DD = ABS(D(M)) + ABS(D(M+1))
            IF (ABS(E(M)) .LE. (EPS * DD)) GOTO 70
            IF (ABS(E(M)) .LE. (EPS4 * ANORM)) GOTO 70
          END DO
          M = N
   70     P = D(L)
          MM1 = M - 1
          IF (M .EQ. L) GOTO 250
          IF (J .EQ. 30) GOTO 900
          J = J + 1
!
!  FORM SHIFT. THIS IS A SLIGHTLY ADVANCED FORM OF SHIFTING MAKING
!  THE ROUTINE ABOUT 20 PERCENT FASTER THAN THE USUAL STUFF.
!
          G = (D(L+1) - P) / (2.0_dp * E(L))
          R = SQRT (G * G + 1.0_dp)
          S = P - E(L) / (G + SIGN (R, G))
          IF (M .EQ. L+1) GOTO 120
          T = S
          R = MAX(ABS(S),(ANORM / N))
          DO I = 1, 6
            PSI = D(M) - T
            PSJ = -1.0_dp
            DO 90 KK = L, MM1
              K = L + MM1 - KK
              IF (ABS(PSI) .GE. (EPS * ABS(E(K)))) GOTO 80
              PSI = BIG
              PSJ = BIG * BIG
              GOTO 90
   80         P = E(K) / PSI
              PSI = D(K) - T - P * E(K)
              PSJ = P * P * PSJ - 1.0_dp
   90       CONTINUE
            IF (ABS(PSJ) .LE. EPS4) GOTO 120
            P = PSI / PSJ
            C = P
            IF (ABS(P) .GT. (0.5_dp * R)) C = SIGN(R,P)
            T = T - C
            IF (ABS(P) .LE. (EPSS * R)) GOTO 110
          END DO
          GOTO 120
  110     S = T
  120     G = D(M) - S
          S = 1.0_dp
          C = 1.0_dp
          P = 0.0_dp
          MML = M - L
!
!  FOR I = M - 1 STEP -1 UNTIL L DO
!
          DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
!
!  SAFE CALCULATION OF SQRT(G * G + F * F) AND SIMILAR STUFF
!
            IF (ABS(F) .LT. ABS(G)) GOTO 150
            C = G / F
            R = SQRT(1.0_dp + C * C)
            E(I+1) = F * R
            S = 1.0_dp / R
            C = C * S
            GOTO 160
  150       S = F / G
            R = SQRT (1.0_dp + S * S)
            E(I+1) = G * R
            C = 1.0_dp / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0_dp * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
            IF (IEV .EQ. 0) GOTO 200
!
!  FORM VECTOR
!
            DO K = 1,N
              F = Z(K,I+1)
              B = Z(K,I)
              Z(K,I+1) = S * B + C * F
              Z(K,I)   = C * B - S * F
          END DO
  200     CONTINUE
          D(L) = D(L) - P
          E(L) = G
          E(M) = 0.0_dp
          GOTO 50
  250   CONTINUE
        RETURN
  900   IER = L
        RETURN
      END SUBROUTINE IQLDIA
!
! ******************************************************************
!
!  This is another version of Iqldia using a less sophisticated
!  shifting algorithm. It is much simpler but 20 percent slower.
!
! ******************************************************************
!
!     SUBROUTINE IQLDIA (NM,N,D,E,Z,IEV,IER)
!       IMPLICIT REAL*8 (A-H,O-Z)
!       DIMENSION  D(N),E(N),Z(NM,N)
!
!       IER = 0
!       IF (N .EQ. 1) RETURN
!       DO 10 I = 2, N
!         E(I-1) = E(I)
!  10   CONTINUE
!       E(N) = 0.0d0
!       DO 250 L = 1, N
!         ITER = 0
!
!  LOOK FOR SMALL SUBDIAGONAL ELEMENT
!
! 100     DO 110 M = L, N-1
!           DD = ABS(D(M)) + ABS(D(M+1))
!           IF ((ABS(E(M)) + DD) .EQ. DD) GOTO 120
! 110     CONTINUE
!         M = N
! 120     IF (M .EQ. L) GOTO 250
!         IF (ITER .EQ. 30) GOTO 900
!         ITER = ITER + 1
!
!  FORM SHIFT
!
!         G = (D(L+1) - D(L)) / (2.0 * E(L))
!         R = SQRT (G * G + 1.0)
!         G = D(M) - D(L) + E(L) / (G + SIGN(R,G))
!         S = 1.0
!         C = 1.0
!         P = 0.0
!
!  FOR I = M - 1 STEP -1 UNTIL L DO
!
!         DO 200 II = 1, M-L
!           I = M - II
!           F = S * E(I)
!           B = C * E(I)
!
!  SAFE CALCULATION OF SQRT(G * G + F * F) AND SIMILAR STUFF
!
!           IF (ABS(F) .LT. ABS(G)) GOTO 150
!           C = G / F
!           R = SQRT(1.0 + C * C)
!           E(I+1) = F * R
!           S = 1.0 / R
!           C = C * S
!           GOTO 160
! 150       S = F / G
!           R = SQRT (1.0d0 + S * S)
!           E(I+1) = G * R
!           C = 1.0d0 / R
!           S = S * C
! 160       G = D(I+1) - P
!           R = (D(I) - G) * S + 2.0d0 * C * B
!           P = S * R
!           D(I+1) = G + P
!           G = C * R - B
!           IF (IEV .EQ. 0) GOTO 200
!
!  FORM VECTOR
!
!           DO 180 K = 1, N
!             F = Z(K,I+1)
!             Z(K,I+1) = S * Z(K,I) + C * F
!             Z(K,I) =   C * Z(K,I) - S * F
! 180       CONTINUE
! 200     CONTINUE
!         D(L) = D(L) - P
!         E(L) = G
!         E(M) = 0.0d0
!         GOTO 100
! 250   CONTINUE
!       RETURN
! 900   IER = L
!       RETURN
!     END
!
! ******************************************************************
!
!     SUBROUTINE BACKTR
!     =================
!
! ******************************************************************
!
!  Backtr solves the system R * X = Y (R upper triangular matrix),
!  where the main diagonal of R is given inverted.
!
!     Parameters:
!       N       (I) :  Dimension of problem
!       M       (I) :  Number of columns in X and Y
!       R       (I) :  Matrix R (upper triangle)
!       NR      (I) :  Dimension of R
!       X       (O) :  Matrix X (solution of system)
!       NX      (I) :  Dimension of X
!       Y       (I) :  Matrix Y (right side)
!       NY      (I) :  Dimension of Y
!       H       (I) :  Auxiliary vector
!
! **********************************************************************
!
      SUBROUTINE BACKTR(N,M,R,NR,X,NX,Y,NY,H)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer :: N,M,NR,NX,NY,i,j,ii,I1,K
        real(dp) :: R,X,Y,H,D,S
!        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION  R(NR,N),X(NX,M),Y(NY,M),H(N)
!
!  CALCULATION OF X = INV(R) * Y
!
        DO II = 1,N
          I = N + 1 - II
          I1 = I + 1
          D = R(I,I)
          DO J= I,N
            H(J)= R(I,J)
          END DO
          DO J = 1,M
            S = Y(I,J)
            DO K = I1,N
              S = S - H(K) * X(K,J)
            END DO
            X(I,J) = S * D
          END DO
        END DO
        RETURN
      END SUBROUTINE BACKTR
!
! ******************************************************************
!
!     SUBROUTINE SORTVC
!     =================
!
! ******************************************************************
!
!  Sortvc sorts D and (if required) E and the columns of Q.
!
!     Prameters:
!
!       NM      (I) :  Dimension of Q
!       N       (I) :  Dimension of problem (size of one vector in Q)
!       NQ      (I) :  Number of elements in D (or columns in Q)
!       D       (I) :  Vector to sort
!               (O) :  Sorted vector
!       Q       (I) :  Matrix to sort (vectors in columns)
!               (O) :  Sorted matrix (vectors in columns)
!       M       (I) :  1: Descending order in D
!                     -1: Ascending order in D
!                      otherwise: no sorting
!       IEV     (I) :  0: No sorting of Q and E
!                      1: Sorting of Q, no sorting of E
!                      2: Sorting of Q and E
!       E       (I) :  Additional Vector to sort
!               (O) :  Sorted additional vector
!
! **********************************************************************
!
      SUBROUTINE SORTVC (NM,N,NQ,D,Q,M,IEV,E)
        use common_accuracy, only : dp
        IMPLICIT NONE
        integer :: NM,M,NQ,IEV,i,j,ii,KK,K,N
        real(dp) :: D,Q,E,H,S
!        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL    LMIN,LMAX
        DIMENSION  D(NQ),E(NQ),Q(NM,NQ)
!
        IF (NQ .LT. 2) RETURN
        LMAX = (M .EQ.  1)
        LMIN = (M .EQ. -1)
        IF (.NOT. (LMAX .OR. LMIN)) RETURN
        DO 40 KK = 2,NQ
          K = KK - 1
          J = K
          H = D(K)
!
!  FIND EXTREMUM
!
          DO 10 I = KK,NQ
            S = D(I)
            IF (LMIN .AND. (S .GE. H)) GOTO 10
            IF (LMAX .AND. (S .LE. H)) GOTO 10
            J = I
            H = S
   10     CONTINUE
          IF (J .EQ. K) GOTO 40
!
!  SORT D
!
          D(J) = D(K)
          D(K) = H
          IF (IEV .EQ. 0) GOTO 40
!
!  SORT Q
!
          DO I = 1,N
            H = Q(I,K)
            Q(I,K) = Q(I,J)
            Q(I,J) = H
          END DO
          IF (IEV .LT. 2) GOTO 40
!
!  SORT E
!
          H    = E(K)
          E(K) = E(J)
          E(J) = H
   40   CONTINUE
        RETURN
      END SUBROUTINE SORTVC

end module diagonalizations
