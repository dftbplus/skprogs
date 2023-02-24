!> Module that provides the density mixing functionality (simple, Broyden).
module broyden

  use common_accuracy, only : dp

  implicit none
  private

  public :: mixing_driver


contains

  !> This is the main driver for simple and broyden mixers, both mix one big one-dimensional array.
  subroutine mixing_driver(pot_old, pot_new, max_l, num_alpha, poly_order, problemsize, iter,&
      & tBroyden, mixing_factor)

    !> old potential
    real(dp), intent(in) :: pot_old(:,0:,:,:)

    !> contains current potential on entry and mixed potential on exit
    real(dp), intent(inout) :: pot_new(:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> current SCF iteration
    integer, intent(in) :: iter

    !> true, if Broyden mixing is desired, otherwise simple mixing is applied
    logical, intent(in) :: tBroyden

    !> mixing factor
    real(dp), intent(in) :: mixing_factor

    !> serialized potentials
    real(dp), allocatable :: vecin(:), vecout(:)

    !> equals current SCF iteration or next one if (iter == 0)
    integer :: titer

    !> auxiliary variables
    integer :: ii, jj, kk, ll, pp

    allocate(vecout(size(pot_old)))
    allocate(vecin(size(pot_old)))

    ! serialize potentials
    pp = 0
    do ii = 1, 2
      do jj = 0, max_l
        do kk = 1, num_alpha(jj) * poly_order(jj)
          do ll = 1, problemsize
            pp = pp + 1
            vecin(pp) = pot_old(ii, jj, kk, ll)
            vecout(pp) = pot_new(ii, jj, kk, ll)
          end do
        end do
      end do
    end do

    ! this check is still necessary, since max. 10000 entries is hardcoded in broyden_mixer
    if (pp > 10000) then
      write(*,*) 'Static dimensions in broyden_mixer too small: ', pp
      stop
    end if

    titer = iter
    ! broyden returns if (iter == 0)
    if (iter == 0) titer = 1

    if (tBroyden) then
      call broyden_mixer(titer, mixing_factor, size(vecin), vecin, vecout)
    else
      call simple_mix(vecin, vecout, mixing_factor)
    end if

    ! deserialize obtained potential
    pp = 0
    do ii = 1, 2
      do jj = 0, max_l
        do kk = 1, num_alpha(jj) * poly_order(jj)
          do ll = 1, problemsize
            pp = pp + 1
            pot_new(ii, jj, kk, ll) = vecin(pp)
          end do
        end do
      end do
    end do

  end subroutine mixing_driver


  !> Simple mixer for last and current density.
  subroutine simple_mix(last, cur, factor)

    !> old vector, holds mixed values at exit
    real(dp), intent(inout) :: last(:)

    !> new vector
    real(dp), intent(in) :: cur(:)

    !> mixing factor
    real(dp), intent(in) :: factor

    last(:) = factor * cur + (1.0_dp - factor) * last

  end subroutine simple_mix


!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BROYDEN_mixer(NITER,ALPHA,JTOP,VECIN,VECOUT)

! This is the Broyden routine as also implemented in the old DFTB code.

       IMPLICIT real(dp) (A-H,O-Z)
       IMPLICIT INTEGER (I-N)
!
!************************************************************
!*  THE VECTORS UI(MAXSIZ) AND VTI(MAXSIZ) ARE JOHNSON'S    *
!*  U(OF I ) AND DF(TRANSPOSE), RESPECTIVELY. THESE ARE     *
!*  CONTINUALLY UPDATED. ALL ITERATIONS ARE S7ORED ON TAPE  *
!*  32 . THIS IS DONE TO PREVENT THE PROHIBITIVE STORAGE    *
!*  COSTS ASSOCIATED WITH HOLDING ONTO THE ENTIRE JACOBIAN. *
!*  VECTOR TL IS THE VT OF EARLIER ITERATIONS. VECTOR F IS: *
!*  VECTOR(OUTPUT) - VECTOR(IN). VECTOR DF IS:  F(M+1)-F(M) *
!*  FINALLY,VECTOR DUMVI(MAXSIZ) IS THE PREDICTED VECTOR.   *
!*  ON RETURN, VECIN CONTAINS THE NEW TRIAL VECTOR.         *
!************************************************************
!*  FOR THE CRAY2-CIVIC ENVIRONMENT , FILES 32 AND 31       *
!*  SHOULD BE INTRODUCED IN THE LINK STATEMENT.             *
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,IMATSZ=40,maxsiz=10000)
!                         formerly    IMATSZ=90
!
! ADDED PARAMETER MAXITER. POREZAG, MAY 1995
!
      PARAMETER(MAXITER=15)
!
! replaced writing to disk by storing values in
! arrays UNIT31, UNIT32  hajnal@scientist.com 2000-10-4
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!      CHARACTER*7 NAMES
!
! SCRATCH COMMON BLOCK FOR LOCAL VARIABLES
!
      DIMENSION VECIN(*),VECOUT(*)
      DIMENSION F(MAXSIZ),UI(MAXSIZ),VTI(MAXSIZ),T1(MAXSIZ),&
     &          VECTOR(MAXSIZ,2),DUMVI(MAXSIZ),DF(MAXSIZ)
!      DIMENSION NAMES(3)
      DIMENSION A(IMATSZ,IMATSZ),B(IMATSZ,IMATSZ),CM(IMATSZ)
      DIMENSION D(IMATSZ,IMATSZ),W(IMATSZ)
      DIMENSION UNIT31(MAXSIZ,2),UNIT32(MAXSIZ,2,MAXITER+15)
!      DATA NAMES/'BROYD01','BROYD02','BROYD03'/
      real(dp) UAMIX,WTMP
      INTEGER ILASTIT
      common /broyd/ uamix, w, WTMP, unit31, unit32, ilastit
      save
!
!     PRINT *,'IN MIXING, WHERE ARE YOU?'
!
! NEW LINES JULY 1996
!
      IF (JTOP .GT. MAXSIZ) THEN
       PRINT *,'BROYDEN: JTOP > MAXSIZ'
       STOP
      END IF
!
! NEW LINES POREZAG, MAY 1995
!
      ITER=NITER
      IF(NITER.GT.MAXITER)ITER=MOD(ITER,MAXITER)+1
      IF(ITER.EQ.0)RETURN
!
! END NEW LINES
!
!      OPEN(66,FILE=NAMES(1),STATUS='UNKNOWN',FORM='FORMATTED')
!      REWIND(66)
!      OPEN(31,FILE=NAMES(2),STATUS='UNKNOWN',FORM='UNFORMATTED')
!      OPEN(32,FILE=NAMES(3),STATUS='UNKNOWN',FORM='UNFORMATTED')
!      REWIND(31)
!      REWIND(32)

!      IF(ITER.EQ.1)THEN
!       ENDFILE 31
!       ENDFILE 32
!      END IF
!
!
!++++++ SET UP THE VECTOR OF THE CURRENT ITERATION FOR MIXING ++++++
!
!  FOR THIS METHOD WE HAVE ONLY SAVED INPUT/OUTPUT CHG. DENSITIES,
      DO K=1,JTOP
      VECTOR(K,1)= VECIN(K)
      VECTOR(K,2)= VECOUT(K)
      END DO
!++++++ END OF PROGRAM SPECIFIC LOADING OF VECTOR FROM MAIN ++++++++
!
!  IVSIZ IS THE LENGTH OF THE VECTOR
      IVSIZ=JTOP
!     IF(ITER.LT.3)WRITE( 6,1001)IVSIZ
      IF(IVSIZ.GT.MAXSIZ)THEN
       PRINT *,'MIXING: EXCEEDED MAXIMAL VECTOR LENGTH'
       STOP
      END IF
!
!
!*******************  BEGIN BROYDEN'S METHOD  **********************
!
!  WEIGHTING FACTOR FOR THE ZEROTH ITERATION
      W0=0.01D0
!
!      F:  THE DIFFERENCE OF PREVIOUS OUTPUT AND INPUT VECTORS
!  DUMVI:  A DUMMY VECTOR, HERE IT IS THE PREVIOUS INPUT VECTOR
!      REWIND(31)
!      READ(31,END=119,ERR=119)AMIX,LASTIT
      IF (ITER .EQ. 1) THEN
        GOTO 119
      ELSE
        AMIX=UAMIX
        LASTIT=ILASTIT
      END IF
!      READ(31)(F(K),K=1,IVSIZ)
      DO k=1,IVSIZ
        F(k)=UNIT31(K,1)
      END DO
!      READ(31)(DUMVI(K),K=1,IVSIZ)
      DO k=1,IVSIZ
        DUMVI(k)=UNIT31(K,2)
      END DO
!      IF(ITER.EQ.1 .AND. LASTIT.GT.1)THEN
!      READ(31)LTMP,((A(I,J),I=1,LTMP),J=1,LTMP)
!      READ(31)(W(I),I=1,LTMP)
!      ENDIF
!
!  ALPHA(OR AMIX)IS SIMPLE MIXING PARAMETERS
!      WRITE(66,1002)AMIX,ITER+1
!
      DO K=1,IVSIZ
      DUMVI(K)=VECTOR(K,1)-DUMVI(K)
      DF(K)=VECTOR(K,2)-VECTOR(K,1)-F(K)
      END DO
      DO K=1,IVSIZ
        F(K)=VECTOR(K,2)-VECTOR(K,1)
      END DO
!
!  FOR I-TH ITER.,DFNORM IS ( F(I) MINUS F(I-1) ), USED FOR NORMALIZATION
!
      DFNORM=ZERO
      FNORM=ZERO
      DO K=1,IVSIZ
        DFNORM=DFNORM + DF(K)*DF(K)
        FNORM=FNORM + F(K)*F(K)
      END DO
      DFNORM=SQRT(DFNORM)
      FNORM=SQRT(FNORM)
!      WRITE(66,'(''  DFNORM '',E12.6,'' FNORM '',E12.6)')DFNORM,FNORM
!
      FAC2=ONE/DFNORM
      FAC1=AMIX*FAC2
!
      DO K=1,IVSIZ
        UI(K) = FAC1*DF(K) + FAC2*DUMVI(K)
        VTI(K)= FAC2*DF(K)
      END DO
!
!*********** CALCULATION OF COEFFICIENT MATRICES *************
!***********    AND THE SUM FOR CORRECTIONS      *************
!
! RECALL: A(I,J) IS A SYMMETRIC MATRIX
!       : B(I,J) IS THE INVERSE OF [ W0**2 I + A ]
!
         LASTIT=LASTIT+1
         LASTM1=LASTIT-1
         LASTM2=LASTIT-2
!
! DUMVI IS THE U(OF I) AND T1 IS THE VT(OF I)
! FROM THE PREVIOUS ITERATIONS
!      REWIND(32)
!      WRITE(66,1003)LASTIT,LASTM1
      IF(LASTIT.GT.2)THEN
      DO J=1,LASTM2
!      READ(32)(DUMVI(K),K=1,IVSIZ)
      DO k=1,IVSIZ
       DUMVI(k)=UNIT32(k,1,J)
      END DO
!      READ(32)(T1(K),K=1,IVSIZ)
      DO k=1,IVSIZ
       T1(k)=UNIT32(k,2,J)
      END DO
!
      AIJ=ZERO
      CMJ=ZERO
      DO  K=1,IVSIZ
        CMJ=CMJ + T1(K)*F(K)
        AIJ=AIJ + T1(K)*VTI(K)
      END DO
      A(LASTM1,J)=AIJ
      A(J,LASTM1)=AIJ
            CM(J)=CMJ
      END DO
      ENDIF
!
      AIJ=ZERO
      CMJ=ZERO
      DO K=1,IVSIZ
        CMJ= CMJ + VTI(K)*F(K)
        AIJ= AIJ + VTI(K)*VTI(K)
      END DO
      A(LASTM1,LASTM1)=AIJ
            CM(LASTM1)=CMJ
!
!      WRITE(32)(UI(K),K=1,IVSIZ)
      DO k=1,IVSIZ
       UNIT32(k,1,LASTM1)=UI(k)
      END DO
!      WRITE(32)(VTI(K),K=1,IVSIZ)
      DO k=1,IVSIZ
       UNIT32(k,2,LASTM1)=VTI(k)
      END DO
!      REWIND(32)
!
! THE WEIGHTING FACTORS FOR EACH ITERATION HAVE BEEN CHOSEN
! EQUAL TO ONE OVER THE R.M.S. ERROR. THIS NEED NOT BE THE CASE.
       IF(FNORM .GT. 1.0D-7)THEN
       WTMP=0.010D0/FNORM
       ELSE
       WTMP=1.0D5
       END IF
       IF(WTMP.LT. 1.00D0) then
         WTMP=1.00D0
       end if
!       print *,wtmp,lastm1,w(lastm1)
       W(LASTM1)=WTMP
!       WRITE(66,'(''  WEIGHTING SET =  '',E12.6)')WTMP
!
!
! WITH THE CURRENT ITERATIONS F AND VECTOR CALCULATED,
! WRITE THEM TO UNIT 31 FOR USE LATER.
!      REWIND(31)
!      WRITE(31)AMIX,LASTIT
      UAMIX=AMIX
      ILASTIT=LASTIT
!      WRITE(31)(F(K),K=1,IVSIZ)
      DO k=1,IVSIZ
        UNIT31(K,1)=F(k)
      END DO
!      WRITE(31)(VECTOR(K,1),K=1,IVSIZ)
      DO k=1,IVSIZ
        UNIT31(K,2)=VECTOR(K,1)
      END DO
!      WRITE(31)LASTM1,((A(I,J),I=1,LASTM1),J=1,LASTM1)
!      WRITE(31)(W(I),I=1,LASTM1)
!
! SET UP AND CALCULATE BETA MATRIX
      DO LM=1,LASTM1
        DO LN=1,LASTM1
          D(LN,LM)= A(LN,LM)*W(LN)*W(LM)
          B(LN,LM)= ZERO
        END DO
        B(LM,LM)= ONE
        D(LM,LM)= W0**2 + A(LM,LM)*W(LM)*W(LM)
      END DO
!
      CALL INVERSE(D,B,LASTM1)
!
!  CALCULATE THE VECTOR FOR THE NEW ITERATION
      DO K=1,IVSIZ
        DUMVI(K)= VECTOR(K,1) + AMIX*F(K)
      END DO
!
      DO I=1,LASTM1
!      READ(32)(UI(K),K=1,IVSIZ)
        DO k=1,IVSIZ
         UI(k)=UNIT32(k,1,I)
        END DO
!       READ(32)(VTI(K),K=1,IVSIZ)
        DO k=1,IVSIZ
         VTI(k)=UNIT32(k,2,I)
        END DO
        GMI=ZERO
        DO IP=1,LASTM1
          GMI=GMI + CM(IP)*B(IP,I)*W(IP)
        END DO
        DO K=1,IVSIZ
          DUMVI(K)=DUMVI(K)-GMI*UI(K)*W(I)
        END DO
      END DO
!  END OF THE CALCULATION OF DUMVI, THE NEW VECTOR
!
!      REWIND(31)
!      REWIND(32)
!
      GOTO 120
! IF THIS IS THE FIRST ITERATION, THEN LOAD
!    F=VECTOR(OUT)-VECTOR(IN) AND VECTOR(IN)
  119 CONTINUE
!     PRINT*,'SIMPLE MIXING THIS ITERATION'
!      REWIND(31)
      LASTIT=1
      AMIX=ALPHA
!      WRITE(31)AMIX,LASTIT
      UAMIX=AMIX
      ILASTIT=LASTIT
      DO K=1,IVSIZ
        F(K)=VECTOR(K,2)-VECTOR(K,1)
      END DO
!      WRITE(31)(F(K),K=1,IVSIZ)
      DO k=1,IVSIZ
        UNIT31(K,1)=F(k)
      END DO
!      WRITE(31)(VECTOR(K,1),K=1,IVSIZ)
      DO k=1,IVSIZ
        UNIT31(K,2)=VECTOR(K,1)
      END DO
!
! SINCE WE ARE ON THE FIRST ITERATION, SIMPLE MIX THE VECTOR.
      DO K=1,IVSIZ
       DUMVI(K)= VECTOR(K,1) + AMIX*F(K)
      END DO
!     WRITE( 6,1000)
  120 CONTINUE
!
!      CLOSE(31,STATUS='KEEP')
!      CLOSE(32,STATUS='KEEP')
!
!*************  THE END OF THE BROYDEN METHOD **************
!
!+++++++ PROGRAM SPECIFIC CODE OF RELOADING ARRAYS +++++++++
!
! NEED TO UNLOAD THE NEW VECTOR INTO THE APPROPRIATE ARRAYS.
      DO K=1,JTOP
        VECIN(K)=DUMVI(K)
      END DO
!
!+++++++++ END OF PROGRAM SPECIFIC RELOADING OF ARRAYS +++++++++
!
!      WRITE(66,1004)ITER+1
!      CLOSE(66)
      RETURN
!
 1000 FORMAT(' ---->  STRAIGHT MIXING ON THIS ITERATION')
 1001 FORMAT(' IN MIXING:   IVSIZ =',I7,/)
 1002 FORMAT(' IN MIXING: SIMPLE MIX PARAMETER',1(F10.6,',')&
     &      ,'  FOR ITER=',I5)
 1003 FORMAT(' CURRENT ITER= ',I5,' INCLUDES VALUES FROM ITER=',I5)
 1004 FORMAT(10X,'DENSITY FOR ITERATION',I4,' PREPARED')
      END subroutine broyden_mixer
!
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INVERSE(A,B,M)
      IMPLICIT real(dp) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

!     =============================================================
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PARAMETER (IMATSZ=40)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      DIMENSION A(IMATSZ,IMATSZ),B(IMATSZ,IMATSZ)
      DIMENSION TD(IMATSZ),AD(IMATSZ),BD(IMATSZ)
      SAVE
!
! SUBROUTINE TO PREFORM GAUSSIAN ELIMINATION
!            NO ZEROS ALONG THE DIAGONAL
!
      N=M
      IF(N.GT.IMATSZ)THEN
       PRINT *,'INVERT: MATRIX A TOO LARGE'
       STOP
      END IF
!
      DO I=1,N
      ATMP=A(I,I)
      IF(ABS(ATMP) .LT. 1.0D-08)THEN
!        WRITE(66,'('' INVERT: MATRIX HAS ZERO DIAGONAL'',
!     &            '' ELEMENT IN THE '',I4,'' ROW'')')I
        STOP
      ENDIF
      END DO
!
      IF(N.EQ.1) GO TO 605
!
      DO I=1,N
!
        DO J=1,N
         TD(J)=A(J,I)/A(I,I)
        END DO
!
!       TD(I)=(0.0E+00,0.0E+00)
        TD(I)=0.0D0
!
        DO K=1,N
           BD(K)=B(I,K)
           AD(K)=A(I,K)
        END DO
!
        DO K=1,N
          DO J=1,N
            B(J,K)=B(J,K)-(TD(J)*BD(K))
            A(J,K)=A(J,K)-(TD(J)*AD(K))
          END DO
        END DO
!
      END DO
!
      DO I=1,N
        DO J=1,N
          B(J,I)=B(J,I)/A(J,J)
        END DO
      END DO
!
      RETURN
!
 605  B(1,1)=1.0D0/A(1,1)
      RETURN
      END subroutine inverse
!

end module broyden
