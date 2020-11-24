module broyden
  use accuracy
  implicit none
  private

  public :: mixing_driver

contains

  ! This is the main driver for simple and broyden mixers, both mix one
  ! big one-dimensional array.
  subroutine mixing_driver(pot_old,pot_new,max_l,num_alpha,&
      &poly_order,problemsize,iter,broyden,mixing_factor)

    integer, intent(in) :: max_l,num_alpha(0:),poly_order(0:),problemsize,iter
    logical, intent(in) :: broyden
    real(dp), intent(in) :: mixing_factor

    integer :: actualsize,titer
    real(dp) :: pot_old(:,0:,:,:),pot_new(:,0:,:,:)
    real(dp), allocatable :: vecin(:),vecout(:)
    integer :: ii,jj,kk,ll,mm,nn,oo,pp

    allocate(vecout(10000))
    allocate(vecin(10000))
    vecout=0.0d0
    vecin=0.0d0

    pp=0
    do ii=1,2
      do jj=0,max_l
        do kk=1,num_alpha(jj)*poly_order(jj)
          do ll=1,problemsize
            pp=pp+1
            vecin(pp)=pot_old(ii,jj,kk,ll)
            vecout(pp)=pot_new(ii,jj,kk,ll)
          end do
        end do
      end do
    end do

    if (pp>10000) then
      write(*,*) 'Static dimensions in broyden_mixer too small',pp
      STOP
    end if

    titer=iter
    ! broyden returns if iter==0
    if (iter==0) titer=1

    if (broyden) then
      call broyden_mixer(titer,mixing_factor,10000,vecin,vecout)
    else
      call simple_mix(vecin,vecout,mixing_factor)
    end if

    pp=0
    do ii=1,2
      do jj=0,max_l
        !      do kk=1,problemsize
        do kk=1,num_alpha(jj)*poly_order(jj)
          do ll=1,problemsize
            pp=pp+1
            !          cof_alt(ii,jj,kk,ll)=vecin(pp)
            !          cof_neu(ii,jj,kk,ll)=vecout(pp)
            pot_new(ii,jj,kk,ll)=vecin(pp)
          end do
        end do
      end do
    end do

    deallocate(vecout)
    deallocate(vecin)

  end subroutine mixing_driver

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BROYDEN_mixer(NITER,ALPHA,JTOP,VECIN,VECOUT)

! This is the Broyden routine as also implemented in the old DFTB code.

       IMPLICIT REAL*8 (A-H,O-Z)
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
      REAL*8 UAMIX,WTMP
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
      DO 38  K=1,JTOP
      VECTOR(K,1)= VECIN(K)
   38 VECTOR(K,2)= VECOUT(K)
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
      DO 104 K=1,IVSIZ
      DUMVI(K)=VECTOR(K,1)-DUMVI(K)
  104 DF(K)=VECTOR(K,2)-VECTOR(K,1)-F(K)
      DO 114 K=1,IVSIZ
  114 F(K)=VECTOR(K,2)-VECTOR(K,1)
!
!  FOR I-TH ITER.,DFNORM IS ( F(I) MINUS F(I-1) ), USED FOR NORMALIZATION
!
      DFNORM=ZERO
      FNORM=ZERO
      DO 113 K=1,IVSIZ
      DFNORM=DFNORM + DF(K)*DF(K)
  113 FNORM=FNORM + F(K)*F(K)
      DFNORM=SQRT(DFNORM)
      FNORM=SQRT(FNORM)
!      WRITE(66,'(''  DFNORM '',E12.6,'' FNORM '',E12.6)')DFNORM,FNORM
!
      FAC2=ONE/DFNORM
      FAC1=AMIX*FAC2
!
      DO 105 K=1,IVSIZ
      UI(K) = FAC1*DF(K) + FAC2*DUMVI(K)
 105  VTI(K)= FAC2*DF(K)
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
      DO 500 J=1,LASTM2
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
      DO 501 K=1,IVSIZ
      CMJ=CMJ + T1(K)*F(K)
  501 AIJ=AIJ + T1(K)*VTI(K)
      A(LASTM1,J)=AIJ
      A(J,LASTM1)=AIJ
            CM(J)=CMJ
  500 CONTINUE
      ENDIF
!
      AIJ=ZERO
      CMJ=ZERO
      DO 106 K=1,IVSIZ
      CMJ= CMJ + VTI(K)*F(K)
  106 AIJ= AIJ + VTI(K)*VTI(K)
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
      DO 506 LM=1,LASTM1
      DO 507 LN=1,LASTM1
         D(LN,LM)= A(LN,LM)*W(LN)*W(LM)
 507     B(LN,LM)= ZERO
         B(LM,LM)= ONE
 506     D(LM,LM)= W0**2 + A(LM,LM)*W(LM)*W(LM)
!
      CALL INVERSE(D,B,LASTM1)
!
!  CALCULATE THE VECTOR FOR THE NEW ITERATION
      DO 505 K=1,IVSIZ
  505 DUMVI(K)= VECTOR(K,1) + AMIX*F(K)
!
      DO 504 I=1,LASTM1
!      READ(32)(UI(K),K=1,IVSIZ)
      DO k=1,IVSIZ
       UI(k)=UNIT32(k,1,I)
      END DO
!      READ(32)(VTI(K),K=1,IVSIZ)
      DO k=1,IVSIZ
       VTI(k)=UNIT32(k,2,I)
      END DO
      GMI=ZERO
      DO 503 IP=1,LASTM1
  503 GMI=GMI + CM(IP)*B(IP,I)*W(IP)
      DO 504 K=1,IVSIZ
  504 DUMVI(K)=DUMVI(K)-GMI*UI(K)*W(I)
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
      DO 101 K=1,IVSIZ
  101 F(K)=VECTOR(K,2)-VECTOR(K,1)
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
      DO 102 K=1,IVSIZ
  102 DUMVI(K)= VECTOR(K,1) + AMIX*F(K)
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
      DO 606 K=1,JTOP
      VECIN(K)=DUMVI(K)
 606  CONTINUE
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
      IMPLICIT REAL*8 (A-H,O-Z)
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
      DO 14 I=1,N
      ATMP=A(I,I)
      IF(ABS(ATMP) .LT. 1.0D-08)THEN
!        WRITE(66,'('' INVERT: MATRIX HAS ZERO DIAGONAL'',
!     &            '' ELEMENT IN THE '',I4,'' ROW'')')I
        STOP
      ENDIF
  14  CONTINUE
!
      IF(N.EQ.1) GO TO 605
!
      DO 23 I=1,N
!
      DO 35 J=1,N
 35      TD(J)=A(J,I)/A(I,I)
!
!     TD(I)=(0.0E+00,0.0E+00)
      TD(I)=0.0D0
!
      DO 71 K=1,N
         BD(K)=B(I,K)
 71      AD(K)=A(I,K)
!
      DO 601 K=1,N
      DO 601 J=1,N
         B(J,K)=B(J,K)-(TD(J)*BD(K))
 601     A(J,K)=A(J,K)-(TD(J)*AD(K))
!
 23   CONTINUE
!
      DO 603 I=1,N
      DO 603 J=1,N
 603     B(J,I)=B(J,I)/A(J,J)
!
      RETURN
!
 605  B(1,1)=1.0D0/A(1,1)
      RETURN
      END subroutine inverse
!

  ! Simple mix, nothing else.
  subroutine simple_mix(alt,neu,factor)
    real(dp), intent(inout) :: alt(:)
    real(dp), intent(in) :: neu(:), factor


! simple mix
    alt=factor*neu+(1.0d0-factor)*alt

  end subroutine simple_mix

end module broyden
