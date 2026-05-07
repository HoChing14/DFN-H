C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
C      I T E R A T I V E    S O L V E R S
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DBCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
     &     MSOLVE, MTSOLV, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     &     R, Z, P, RR, ZZ, PP, DZ, RWORK, IWORK)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      integer*8 ITER, IERR, IUNIT, IWORK(*)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N)
      REAL*8 RR(N), ZZ(N), PP(N), DZ(N), RWORK(*)
      REAL*8 DMACH(5)
c      DATA DMACH(3) / 1.1101827117665d-16 /
      DATA DMACH(3) / 1.1101827117665d-30 /
      EXTERNAL MATVEC, MTTVEC, MSOLVE, MTSOLV
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DBCG
C
C
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      FUZZ   = DMACH(3)
      TOLMIN = 500.d0*FUZZ
      FUZZ   = FUZZ*FUZZ
      IF( TOL.LT.TOLMIN ) THEN
         TOL  = TOLMIN
         IERR = 4
      ENDIF
C         
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I)  = B(I) - R(I)
         RR(I) = R(I)
 10   CONTINUE
 
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      CALL MTSOLV(N, RR, ZZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C                  
      ISDBCG = 0
      ITOL   = 1
      IF(ITER .EQ. 0) THEN
         BNRM = 0.d0
         DO 11 I5=1,N
            BNRM = BNRM+B(I5)*B(I5)
   11    CONTINUE
         BNRM = SQRT(BNRM)
      ENDIF
       ERR = 0.0d0
       DO 12 I5=1,N
         ERR = ERR+R(I5)*R(I5)
   12 CONTINUE
      ERR = SQRT(ERR)/BNRM
C         
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF
      IF(ERR .LE. TOL) ISDBCG = 1
C         
      IF(ISDBCG.NE.0) GO TO 200
C         
C         
      DO 100 K=1,ITMAX
         ITER = K
C         
         DDOT = 0.D0
         DO 15 I = 1,N
           DDOT = DDOT + Z(I)*RR(I)
   15    CONTINUE
         BKNUM = DDOT
         IF( ABS(BKNUM).LE.FUZZ ) THEN
            IERR = 6
            RETURN
         ENDIF
         IF(ITER .EQ. 1) THEN
            DO 18 I = 1,N
              P(I) = Z(I)
              PP(I) = ZZ(I)
   18       CONTINUE         
         ELSE
            BK = BKNUM/BKDEN
            DO 20 I = 1, N
               P(I)  = Z(I) + BK*P(I)
               PP(I) = ZZ(I) + BK*PP(I)
 20         CONTINUE
         ENDIF
         BKDEN = BKNUM
C         
         CALL MATVEC(N, P, Z, NELT, IA, JA, A, ISYM)
         DDOT = 0.D0
         DO 25 I = 1,N
           DDOT = DDOT + PP(I)*Z(I)
   25    CONTINUE
         AKDEN = DDOT
         AK = BKNUM/AKDEN
         IF( ABS(AKDEN).LE.FUZZ ) THEN
            IERR = 6
            RETURN
         ENDIF
C         
          DO 26 I = 1,N
             X(I) = X(I) + AK*P(I)
   26     CONTINUE
C
          DO 27 I = 1,N
             R(I) = R(I) - AK*Z(I)
   27     CONTINUE
C         
         CALL MTTVEC(N, PP, ZZ, NELT, IA, JA, A, ISYM)
C         
          DO 28 I = 1,N
             RR(I) = RR(I) - AK*ZZ(I)
   28     CONTINUE
C         
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL MTSOLV(N, RR, ZZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
         ISDBCG = 0
         ITOL   = 1
         ERR    = 0.0d0
         DO 33 I5=1,N
            ERR = ERR+R(I5)*R(I5)
   33    CONTINUE
         ERR = SQRT(ERR)/BNRM
C         
         IF(IUNIT .NE. 0) THEN
            WRITE(IUNIT,1010) ITER, ERR, AK, BK
         ENDIF
         IF(ERR .LE. TOL) ISDBCG = 1
C         
         IF(ISDBCG.NE.0) GO TO 200
C         
 100  CONTINUE
 1000 FORMAT(' Preconditioned BiConjugate Gradient for N, ITOL = ',
     $     I5,I5,/' ITER','   Error Estimate','            Alpha',
     $     '             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
C         
C         
      ITER = ITMAX + 1
      IERR = 2
C         
 200  RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSLUBC(N,B,X,NELT,IA,JA,A,TOL,ITMAX,ITER,ERR,
     &                  IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      integer*8 LENIW
      integer*8 IERR, IUNIT, LENW, IWORK(LENIW)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSMTV, DSLUI, DSLUTI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSLUBC
C
C
      ISYM = 0
      ITOL = 1
C
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM, IUNIT )
C
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCR = LOCU + NU
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCRR = LOCP + N
      LOCZZ = LOCRR + N
      LOCPP = LOCZZ + N
      LOCDZ = LOCPP + N
      LOCW = LOCDZ + N
C
      CALL DCHKW( 'DSLUBC',LOCIW,LENIW,LOCW,LENW,IERR,ITER,ERR)
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C         
      CALL DBCG(N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSMTV,
     $     DSLUI, DSLUTI, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $     RWORK(LOCR), RWORK(LOCZ), RWORK(LOCP),
     $     RWORK(LOCRR), RWORK(LOCZZ), RWORK(LOCPP),
     $     RWORK(LOCDZ), RWORK, IWORK )
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC,
     &     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     &     R, R0, P, Q, U, V1, V2, RWORK, IWORK)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      integer*8 ITER, IERR, IUNIT, IWORK(*)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
      REAL*8 Q(N), U(N), V1(N), V2(N), RWORK(*)
      REAL*8 DMACH(5)
c      DATA DMACH(3) / 1.1101827117665d-16 /
      DATA DMACH(3) / 1.1101827117665d-30 /
      EXTERNAL MATVEC, MSOLVE
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DCGS
C
C
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500.d0*DMACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
C         
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         V1(I)  = R(I) - B(I)
 10   CONTINUE
      CALL MSOLVE(N, V1, R, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         
      ISDCGS = 0
      ITOL   = 2
      IF(ITER.EQ. 0) THEN
         CALL MSOLVE(N, B, V2, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         BNRM = 0.0d0
         DO 11 I5=1,N
            BNRM = BNRM+V2(I5)*V2(I5)
   11    CONTINUE
         BNRM = SQRT(BNRM)
      ENDIF
      ERR = 0.0d0
      DO 12 I5=1,N
         ERR = ERR+R(I5)*R(I5)
   12 CONTINUE
      ERR = SQRT(ERR)/BNRM
C         
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF
      IF(ERR .LE. TOL) ISDCGS = 1
C         
      IF(ISDCGS.NE.0) GO TO 200
C
C
      FUZZ = DMACH(3)**2
      DO 20 I = 1, N
         R0(I) = R(I)
 20   CONTINUE
      RHONM1 = 1.0d0
C         
C         
      DO 100 K=1,ITMAX
         ITER = K
C
         DDOT = 0.D0
         DO 15 I = 1,N
           DDOT = DDOT + R0(I)*R(I)
   15    CONTINUE
         RHON = DDOT
C
         IF( ABS(RHONM1).LT.FUZZ ) GOTO 998
         BK = RHON/RHONM1
         IF( ITER.EQ.1 ) THEN
            DO 30 I = 1, N
               U(I) = R(I)
               P(I) = R(I)
 30         CONTINUE
         ELSE
            DO 40 I = 1, N
               U(I) = R(I) + BK*Q(I)
               V1(I) = Q(I) + BK*P(I)
 40         CONTINUE
            DO 50 I = 1, N
               P(I) = U(I) + BK*V1(I)
 50         CONTINUE
         ENDIF
C         
         CALL MATVEC(N, P, V2, NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
         DDOT = 0.D0
         DO 25 I = 1,N
           DDOT = DDOT + R0(I)*V1(I)
   25    CONTINUE
         SIGMA = DDOT
C
         IF( ABS(SIGMA).LT.FUZZ ) GOTO 999
         AK = RHON/SIGMA
         AKM = -AK
         DO 60 I = 1, N
            Q(I) = U(I) + AKM*V1(I)
 60      CONTINUE
 
         DO 70 I = 1, N
            V1(I) = U(I) + Q(I)
 70      CONTINUE
C
          DO 72 I = 1,N
             X(I) = X(I) + AKM*V1(I)
   72     CONTINUE
C                     -1
         CALL MATVEC(N, V1, V2, NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
         DO 73 I = 1,N
            R(I) = R(I) + AKM*V1(I)
   73    CONTINUE
C         
         ISDCGS = 0
         ITOL   = 2
         ERR    = 0.0d0
         DO 78 I5=1,N
            ERR = ERR+R(I5)*R(I5)
   78    CONTINUE
         ERR = SQRT(ERR)/BNRM
C         
         IF(IUNIT .NE. 0) THEN
            WRITE(IUNIT,1010) ITER, ERR, AK, BK
         ENDIF
         IF(ERR .LE. TOL) ISDCGS = 1
C         
         IF(ISDCGS.NE.0) GO TO 200
C
         RHONM1 = RHON
 100  CONTINUE
 1000 FORMAT(' Preconditioned BiConjugate Gradient Squared for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Error Estimate','            Alpha',
     $     '             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
C         
      ITER = ITMAX + 1
      IERR = 2
 200  RETURN
C
 998  IERR = 5
      RETURN
C
 999  IERR = 6
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSLUCS(N,B,X,NELT,IA,JA,A,TOL,ITMAX,ITER,ERR,
     &                  IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      integer*8 LENIW
      integer*8 IERR, IUNIT, LENW, IWORK(LENIW)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSLUI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSLUCS
C
C
      ISYM = 0
      ITOL = 2
C
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM, IUNIT )
C
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL   = LOCRB
      LOCDIN = LOCL + NL
      LOCUU  = LOCDIN + N
      LOCR   = LOCUU + NU
      LOCR0  = LOCR + N
      LOCP   = LOCR0 + N
      LOCQ   = LOCP + N
      LOCU   = LOCQ + N
      LOCV1  = LOCU + N
      LOCV2  = LOCV1 + N
      LOCW   = LOCV2 + N
C
      CALL DCHKW( 'DSLUCS', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCUU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCUU), IWORK(LOCNR), IWORK(LOCNC) )
C         
      CALL DCGS(N, B, X, NELT, IA, JA, A, ISYM, DSMV,
     $     DSLUI, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $     RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP),
     $     RWORK(LOCQ), RWORK(LOCU), RWORK(LOCV1),
     $     RWORK(LOCV2), RWORK, IWORK )
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DCGSTB(N, B, X, NELT, IA, JA, A, ISYM, MATVEC,
     &     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     &     R, R0, P, Q, U, V1, V2, RWORK, IWORK)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      integer*8 ITER, IERR, IUNIT, IWORK(*)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
      REAL*8 Q(N), U(N), V1(N), V2(N), RWORK(*)
      REAL*8 DMACH(5)
c      DATA DMACH(3) / 1.1101827117665d-16 /
      DATA DMACH(3) / 1.1101827117665d-30 /
      EXTERNAL MATVEC, MSOLVE
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DCGSTB
C
C
      ITER = 0
      IERR = 0
C
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500.d0*DMACH(3)
      IF(TOL.LT.TOLMIN) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
C         
      DO 2 I = 1, N
         R(I) = B(I)
    2 CONTINUE
C
      ERR1 = 0.0d0
      DO 4 i=1,n
         err1 = err1+x(i)*x(i)
    4 CONTINUE
C
C
C
      IF(err1.NE.0.0d0) THEN
         CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
C
         ERR2 = 0.0d0
         DO 5 i=1,n
            R(I)  = B(I)-R(I)
            ERR2  = ERR2+R(I)*R(I)
    5    CONTINUE
         ERR2 = SQRT(ERR2)
cels11/6/01 added initialization of IBCGST
         IBCGST = 0
         IF(ERR2.LE.TOL) IBCGST = 1
         IF(IBCGST.NE.0) GO TO 200
      END IF
C
C
C
      RHONM1 = 1.0d0
      ALPHA  = 1.0d0
      OMEGA  = 1.0d0
      BETA   = 0.0d0
C
C
C         
      IBCGST = 0
      ITOL   = 2
      BNRM   = 0.0d0
      DO 11 I5=1,N
         BNRM = BNRM+B(I5)*B(I5)
   11 CONTINUE
      BNRM = SQRT(BNRM)
C
      ERR = 0.0d0
      DO 12 I5=1,N
         ERR = ERR+R(I5)*R(I5)
   12 CONTINUE
      ERR = SQRT(ERR)/BNRM
C         
      IF(IUNIT.NE.0) THEN
         IF(ITER.EQ.0) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, Alpha, BETA, OMEGA
      ENDIF
      IF(ERR.LE.TOL) IBCGST = 1
C         
      IF(IBCGST.NE.0) GO TO 200
C
C
      FUZZ = DMACH(3)**2
      DO 14 I = 1, N
         R0(I) = R(I)
   14 CONTINUE
C         
C
C         
      DO 100 K=1,ITMAX
         ITER = K
C
         DDOT = 0.D0
         DO 15 I = 1,N
           DDOT = DDOT + R0(I)*R(I)
   15    CONTINUE
         RHON = DDOT
C
         IF(ABS(RHONM1).LT.FUZZ) GOTO 998
         BETA = (RHON/RHONM1)*(ALPHA/OMEGA)
C
         IF( ITER.EQ.1 ) THEN
            DO 20 I = 1, N
               P(I)  = R(I)
   20       CONTINUE
         ELSE
            DO 25 I = 1, N
               U(I) = P(I)-OMEGA*V1(I)
   25       CONTINUE
            DO 26 I = 1, N
               P(I) = R(I)+BETA*U(I)
   26       CONTINUE
         ENDIF
C         
         CALL MSOLVE(N, P, V2, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL MATVEC(N, V2, V1, NELT, IA, JA, A, ISYM)
C
         DDOT = 0.0d0
         DO 30 I = 1,N
           DDOT = DDOT + R0(I)*V1(I)
   30    CONTINUE
         SIGMA = DDOT
C
         IF( ABS(SIGMA).LT.FUZZ ) GOTO 999
         ALPHA = RHON/SIGMA
C
         DO 35 I=1,N
            Q(I) = R(I)-ALPHA*V1(I)
   35    CONTINUE         
C
         ERR = 0.0
         DO 40 I5=1,N
            ERR = ERR+Q(I5)*Q(I5)
   40    CONTINUE
         ERR = SQRT(ERR)/BNRM
C         
         IF(ERR.LE.TOL) IBCGST = 1
C
         DO 42 I5=1,N
            X(I5) = X(I5)+ALPHA*V2(I5)
   42    CONTINUE
C         
         IF(IBCGST.NE.0) THEN 
            IF(IUNIT.NE.0) THEN
               WRITE(IUNIT,1010) ITER, ERR, Alpha, BETA, OMEGA
            ENDIF
            GO TO 200
         END IF       
C                      
         CALL MSOLVE(N, Q, V2, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL MATVEC(N, V2, U, NELT, IA, JA, A, ISYM)
C
         SUM1=0.0d0 
         SUM2=0.0d0 
         DO 50 I = 1, N
            SUM1 = SUM1+U(I)*Q(I)
            SUM2 = SUM2+U(I)*U(I)
   50    CONTINUE
         OMEGA = SUM1/SUM2
C
          DO 60 I = 1,N
             X(I) = X(I)+OMEGA*V2(I)
   60     CONTINUE
C
         DO 65 I = 1,N
            R(I) = Q(I) - OMEGA*U(I)
   65    CONTINUE
C         
         IBCGST = 0
         ITOL   = 2
         ERR    = 0.0d0
         DO 78 I5=1,N
            ERR = ERR+R(I5)*R(I5)
   78    CONTINUE
         ERR = SQRT(ERR)/BNRM
C         
         IF(IUNIT .NE. 0) THEN
            WRITE(IUNIT,1010) ITER, ERR, Alpha, BETA, OMEGA
         ENDIF
         IF(ERR .LE. TOL) IBCGST = 1
C         
         IF(IBCGST.NE.0) GO TO 200
C
         RHONM1 = RHON
 100  CONTINUE
 1000 FORMAT(' Preconditioned BiConjugate Gradient Stabilized for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Error Estimate','            Alpha',
     $     '             Beta','             Omega')
 1010 FORMAT(1X,I4,1X,1pE16.7,1X,1pE16.7,1X,1pE16.7,1X,1pE16.7)
C         
      ITER = ITMAX + 1
      IERR = 2
 200  RETURN
C
 998  IERR = 5
      RETURN
C
 999  IERR = 6
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DLUSTB(N,B,X,NELT,IA,JA,A,TOL,ITMAX,ITER,ERR,
     &                  IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      integer*8 LENIW
      integer*8 IERR, IUNIT, LENW, IWORK(LENIW)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSLUI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DLUSTB
C
C
      ISYM = 0
      ITER = 2
C
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM, IUNIT )
C
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
cels9/29/09 just added line below
cels10/15/09 doesn't seem to make a difference     IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL   = LOCRB
      LOCDIN = LOCL + NL
      LOCUU  = LOCDIN + N
      LOCR   = LOCUU + NU
      LOCR0  = LOCR + N
      LOCP   = LOCR0 + N
      LOCQ   = LOCP + N
      LOCU   = LOCQ + N
      LOCV1  = LOCU + N
      LOCV2  = LOCV1 + N
      LOCW   = LOCV2 + N
C
      CALL DCHKW('DLUSTB',LOCIW,LENIW,LOCW,LENW,IERR,ITER,ERR)
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCUU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCUU), IWORK(LOCNR), IWORK(LOCNC) )
C         
      CALL DCGSTB(N, B, X, NELT, IA, JA, A, ISYM, DSMV,
     $     DSLUI, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $     RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP),
     $     RWORK(LOCQ), RWORK(LOCU), RWORK(LOCV1),
     $     RWORK(LOCV2), RWORK, IWORK )
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
     &     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, 
     &     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK )
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      integer*8  LIGW
      integer*8  IERR, IUNIT, LRGW, IGWK(LIGW)
      integer*8  IWORK(*)
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N)
      REAL*8 RGWK(LRGW), RWORK(*)
      EXTERNAL MATVEC, MSOLVE
      integer*8 JPRE, KMP, MAXL, NMS, MAXLP1, NMSL, NRSTS, NRMAX
      integer*8 I, IFLAG, LR, LDL, LHES, LGMR, LQ, LV, LW
      REAL*8 BNRM, RHOL, SUM
      REAL*8 DMACH(5)
c      DATA DMACH(3) / 1.1101827117665d-16 /
      DATA DMACH(3) / 1.1101827117665d-30 /
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGMRES
C
C
      IERR = 0
      MAXL = IGWK(1)
      IF (MAXL .EQ. 0) MAXL = 10
      IF (MAXL .GT. N) MAXL = N
      KMP = IGWK(2)
      IF (KMP .EQ. 0) KMP = MAXL
      IF (KMP .GT. MAXL) KMP = MAXL
      JSCAL = IGWK(3)
      JPRE = IGWK(4)
C
      IF( ITOL.EQ.1 .AND. JPRE.LT.0 ) GOTO 650
      IF( ITOL.EQ.2 .AND. JPRE.GE.0 ) GOTO 650
      NRMAX = IGWK(5)
C      
      IF( NRMAX.EQ.0 ) NRMAX = 10
      IF( NRMAX.EQ.-1 ) NRMAX = 0
      IF( TOL.EQ.0.0d0 ) TOL = 500.d0*DMACH(3)
C
      ITER  = 0
      NMS   = 0
      NRSTS = 0
C
C
C
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + N*MAXLP1
      LHES = LR + N + 1
      LQ = LHES + MAXL*MAXLP1
      LDL = LQ + 2*MAXL
      LW = LDL + N
      LXL = LW + N
      LZ = LXL + N
C
      IGWK(6) = LZ + N - 1
      IF( LZ+N-1.GT.LRGW ) GOTO 640
      
      IF (JPRE .LT. 0) THEN
         CALL MSOLVE(N,B,RGWK(LR),NELT,IA,JA,A,ISYM,RWORK,IWORK)
         NMS = NMS + 1
      ELSE
         DO  7 I = 1,N
           RGWK(LR+I-1) = B(I)
    7    CONTINUE               
      ENDIF
      IF( JSCAL.EQ.2 .OR. JSCAL.EQ.3 ) THEN
         BNRM = 0.D0
         DO 10 I = 1,N
            BNRM = BNRM +(RGWK(LR-1+I)*SB(I))
     &                  *(RGWK(LR-1+I)*SB(I))
 10      CONTINUE
         BNRM = SQRT(BNRM)
      ELSE
         BNRM = 0.D0
         DO 12 I = 1,N
            BNRM = BNRM + RGWK(LR-1+I)*RGWK(LR-1+I)
 12      CONTINUE
         BNRM = SQRT(BNRM)
      ENDIF
      
      CALL MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)
      
      DO 50 I = 1,N
         RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)
 50   CONTINUE
C
 100  CONTINUE
      IF( NRSTS.GT.NRMAX ) GOTO 610      
      IF( NRSTS.GT.0 ) THEN
         DO 17 I = 1,N
           RGWK(LR+I-1) = RGWK(LDL+I-1)
   17    CONTINUE               
      ENDIF
C
      CALL DPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,
     $       NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV),
     $       RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW),
     $       RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL,
     $       TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
      ITER = ITER + LGMR
      NMS = NMS + NMSL
C
C
      LZM1 = LZ - 1
      DO 110 I = 1,N
         X(I) = X(I) + RGWK(LZM1+I)
 110  CONTINUE
      IF( IFLAG.EQ.0 ) GOTO 600
      IF( IFLAG.EQ.1 ) THEN
         NRSTS = NRSTS + 1
         GOTO 100
      ENDIF
      IF( IFLAG.EQ.2 ) GOTO 620
C
 600  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 0
      RETURN
C
 610  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 1
      RETURN
C
 620  CONTINUE
      IGWK(7) = NMS
      RGWK(1) = RHOL
      IERR = 2
      RETURN
C      
 640  CONTINUE
      ERR = TOL
      IERR = -1
      RETURN
C      
 650  CONTINUE
      ERR = TOL
      IERR = -2
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSLUGM(N,B,X,NELT,IA,JA,A,NSAVE,TOL,ITMAX,ITER,
     &                  ERR,IERR,IUNIT,RWORK,LENW,IWORK,LENIW)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL
      integer*8  ITMAX, ITER, IERR, IUNIT, LENW, LENIW, IWORK(LENIW)
      integer*8 locrgw
      REAL*8 B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV, DSLUI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSLUGM
C
C
      ISYM = 0
      ITOL = 0
C
      IERR = 0
      ERR  = 0.d0
      IF( NSAVE.LE.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM, IUNIT )
C
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIGW = LOCIB
      LOCIL = LOCIGW + 20
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCRGW = LOCU + NU
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)
C
      CALL DCHKW( 'DSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C         
      IWORK(LOCIGW  ) = NSAVE
      IWORK(LOCIGW+1) = NSAVE
      IWORK(LOCIGW+2) = 0
      IWORK(LOCIGW+3) = -1
      IWORK(LOCIGW+4) = ITMAX/NSAVE
      MYITOL = 0
C      
      lenlocdf = LENW-LOCRGW
      i20 = 20
      CALL DGMRES( N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI,
     $     MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,
     $     RWORK(LOCRGW), lenlocdf, IWORK(LOCIGW), i20,
     $     RWORK, IWORK )
C
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DHELS(A, LDA, N, Q, B)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 LDA, N
      REAL*8 A(LDA,*), B(*), Q(*)
C
      integer*8 IQ, K, KB, KP1
      REAL*8 C, S, T, T1, T2
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DHELS
C
C
      DO 20 K = 1, N
         KP1 = K + 1
         IQ     = 2*(K-1) + 1
         C      = Q(IQ)
         S      = Q(IQ+1)
         T1     = B(K)
         T2     = B(KP1)
         B(K)   = C*T1 - S*T2
         B(KP1) = S*T1 + C*T2
 20   CONTINUE
C
C
      DO 40 KB = 1, N
         K    = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T    =-B(K)
         DO 38 I = 1,K-1
            B(I) = B(I) + T*A(I,K)
   38    CONTINUE
 40   CONTINUE
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DHEQR(A, LDA, N, Q, INFO, IJOB)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 LDA, N, INFO, IJOB
      REAL*8 A(LDA,*), Q(*)
C
C
      integer*8 I, IQ, J, K, KM1, KP1, NM1
      REAL*8 C, S, T, T1, T2
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DHEQR
C
C
      IF (IJOB .GT. 1) GO TO 70
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
         IF (KM1 .LT. 1) GO TO 20
         DO 10 J = 1, KM1
            I = 2*(J-1) + 1
            T1 = A(J,K)
            T2 = A(J+1,K)
            C = Q(I)
            S = Q(I+1)
            A(J,K) = C*T1 - S*T2
            A(J+1,K) = S*T1 + C*T2
 10      CONTINUE
C
C
 20      CONTINUE
         IQ = 2*KM1 + 1
         T1 = A(K,K)
         T2 = A(KP1,K)
         IF( T2.EQ.0.0d0 ) THEN
            C = 1.0d0
            S = 0.0d0
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
            T = T1/T2
            S = -1.0d0/SQRT(1.0d0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0d0/SQRT(1.0d0+T*T)
            S = -C*T
         ENDIF
         Q(IQ) = C
         Q(IQ+1) = S
         A(K,K) = C*T1 - S*T2
         IF( A(K,K).EQ.0.0d0 ) INFO = K
 60   CONTINUE
      RETURN
C      
 70   CONTINUE
      NM1 = N - 1
C      
      DO 100 K = 1,NM1
         I = 2*(K-1) + 1
         T1 = A(K,N)
         T2 = A(K+1,N)
         C = Q(I)
         S = Q(I+1)
         A(K,N) = C*T1 - S*T2
         A(K+1,N) = S*T1 + C*T2
 100  CONTINUE
C
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF ( T2.EQ.0.0d0 ) THEN
         C = 1.0d0
         S = 0.0d0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T = T1/T2
         S = -1.0d0/SQRT(1.0d0+T*T)
         C = -S*T
      ELSE
         T = T2/T1
         C = 1.0d0/SQRT(1.0d0+T*T)
         S = -C*T
      ENDIF
C      
      IQ      = 2*N - 1
      Q(IQ)   = C
      Q(IQ+1) = S
      A(N,N)  = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0d0) INFO = N
C
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, LL, LDHES, KMP
      REAL*8 VNEW, V, HES, SNORMW
      DIMENSION VNEW(*), V(N,*), HES(LDHES,*)
C
      integer*8 I, I0
      REAL*8 ARG, SUMDSQ, TEM, VNRM
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DORTH
C
C
       VNRM = 0.D0
       DO 2 I = 1,N
          VNRM = VNRM + VNEW(I)*VNEW(I)
  2    CONTINUE
       VNRM = SQRT(VNRM)
c---------- The following function gives integer type mismatch warning
c      I0 = MAX0(1,LL-KMP+1)
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
         DDOT = 0.D0
         DO 5 J = 1,N
            DDOT = DDOT + V(J,I)*VNEW(J)
    5    CONTINUE
         HES(I,LL) = DDOT
C
         TEM = -HES(I,LL)
         DO 8 J = 1,N
            VNEW(J) = VNEW(J) + TEM*V(J,I)
    8    CONTINUE
 10   CONTINUE
C 
      SNORMW = 0.D0
      DO 9 I = 1,N
         SNORMW = SNORMW + VNEW(I)*VNEW(I)
  9   CONTINUE
      SNORMW = SQRT(SNORMW)
C         
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0d0
      DO 30 I = I0,LL
         DDOT = 0.D0
         DO 25 J = 1,N
            DDOT = DDOT + V(J,I)*VNEW(J)
   25    CONTINUE
         TEM = -DDOT
C
         IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
         HES(I,LL) = HES(I,LL) - TEM
C         
         DO 28 J = 1,N
            VNEW(J) = VNEW(J) + TEM*V(J,I)
   28    CONTINUE
         SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
C 
      IF (SUMDSQ .EQ. 0.0d0) RETURN
      ARG = MAX(0.0d0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
C
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, 
     &     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
     &     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
     &     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      EXTERNAL MATVEC, MSOLVE
      integer*8 N,MAXL,MAXLP1,KMP,JPRE,NMSL,LGMR,IFLAG,JSCAL,NRSTS
      integer*8 NRMAX,ITOL,NELT,ISYM
      REAL*8 RHOL, BNRM, TOL
      REAL*8 R0(*), SR(*), SZ(*), Z(*), V(N,*)
      REAL*8 HES(MAXLP1,*), Q(*), RPAR(*), WK(*), DL(*)
      REAL*8 A(NELT), B(*), X(*), XL(*)
      integer*8 IPAR(*), IA(NELT), JA(NELT)
C
      integer*8 I, INFO, IP1, I2, J, K, LL, LLP1
      REAL*8 R0NRM,C,DLNRM,PROD,RHO,S,SNORMW,TEM
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of LINEQ
C
C
      DO 5 I = 1,N
         Z(I) = 0.0d0
 5    CONTINUE
C
      IFLAG = 0
      LGMR  = 0
      NMSL  = 0
      ITMAX =(NRMAX+1)*MAXL
      
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN
      
         DO 8 I = 1,N
            WK(I) = R0(I)
    8    CONTINUE         
         
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN
         DO 10 I = 1,N
            V(I,1) = R0(I)*SR(I)
 10      CONTINUE
      ELSE
         DO 20 I = 1,N
            V(I,1) = R0(I)
 20      CONTINUE
      ENDIF
      
      R0NRM = 0.D0
      DO 22 I = 1,N
         R0NRM = R0NRM + V(I,1)*V(I,1)
 22   CONTINUE
      R0NRM = SQRT(R0NRM)
      ITER = NRSTS*MAXL
C
C         
      ISDGMR = 0
      ITOL   = 0
      ERR    = R0NRM/BNRM
C         
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, ERR
      ENDIF
      IF ( ERR.LE.TOL ) ISDGMR = 1
C         
      IF(ISDGMR.NE.0) RETURN
C
      TEM = 1.0d0/R0NRM
C      
      DO 33 I = 1,N
        V(I,1) = TEM*V(I,1)
   33 CONTINUE
C
C
      DO 50 J = 1,MAXL
         DO 40 I = 1,MAXLP1
            HES(I,J) = 0.0d0
 40      CONTINUE
 50   CONTINUE
C
      PROD = 1.0d0
      DO 90 LL = 1,MAXL
         LGMR = LL         
         IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
            DO 60 I = 1,N
               WK(I) = V(I,LL)/SZ(I)
 60         CONTINUE
         ELSE
            DO 63 I = 1,N
               WK(I) = V(I,LL)
   63       CONTINUE         
         ENDIF
C        
        IF (JPRE .GT. 0) THEN
           CALL MSOLVE(N,WK,Z,NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
           CALL MATVEC(N,Z,V(1,LL+1),NELT,IA,JA,A,ISYM)
        ELSE
           CALL MATVEC(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM)
        ENDIF
C        
        IF (JPRE .LT. 0) THEN
C        
           DO 64 I = 1,N
              WK(I) = V(I,LL+1)
   64      CONTINUE         
C           
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
        ENDIF
C        
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN
           DO 65 I = 1,N
              V(I,LL+1) = V(I,LL+1)*SR(I)
 65        CONTINUE
        ENDIF
C
        CALL DORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
C
        HES(LL+1,LL) = SNORMW
        CALL DHEQR(HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C        
        PROD = PROD*Q(2*LL)
        RHO  = ABS(PROD*R0NRM)
C        
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              DO 68 I = 1,N
                 DL(I) = V(I,1)
   68         CONTINUE         
C
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,N
                    DL(K) = S*DL(K) + C*V(K,IP1)
 70              CONTINUE
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,N
              DL(K) = S*DL(K) + C*V(K,LLP1)
 80        CONTINUE
C
           DLNRM = 0.D0
           DO 82 I = 1,N
              DLNRM = DLNRM + DL(I)*DL(I)
 82        CONTINUE
           DLNRM = SQRT(DLNRM)
           RHO   = RHO*DLNRM
        ENDIF
        RHOL = RHO
C        
        ITER = NRSTS*MAXL + LGMR
C
C         
        ISDGMR = 0
        ITOL   = 0
        ERR    = RHOL/BNRM      
C         
        IF( IUNIT.NE.0 ) THEN
           IF( ITER.EQ.0 ) THEN
              WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
           ENDIF
           WRITE(IUNIT,1010) ITER, ERR, ERR
        ENDIF
        IF ( ERR.LE.TOL ) ISDGMR = 1
C         
        IF(ISDGMR.NE.0) GO TO 200
C
        IF (LL .EQ. MAXL) GO TO 100
        TEM = 1.0d0/SNORMW
        
      DO 88 I = 1,N
        V(I,LL+1) = TEM*V(I,LL+1)
   88 CONTINUE
C        
 90   CONTINUE
 100  CONTINUE
C
      IF (RHO .LT. R0NRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
C
C
      DO 130 I = 1,N
         Z(I) = 0.D0
 130  CONTINUE
      RETURN
 150  IFLAG = 1
C
C
      IF (NRMAX .GT. 0) THEN
C
         IF (KMP .EQ. MAXL) THEN
            DO 158 I = 1,N
               DL(I) = V(I,1)
  158       CONTINUE         
C
            LLM1 = MAXL - 1
            DO 165 I = 1,LLM1
               IP1 = I + 1
               I2  = I*2
               S   = Q(I2)
               C   = Q(I2-1)
               DO 160 K = 1,N
                  DL(K) = S*DL(K) + C*V(K,IP1)
160            CONTINUE
165         CONTINUE
            S    = Q(2*MAXL)
            C    = Q(2*MAXL-1)/SNORMW
            LLP1 = MAXL + 1
            DO 170 K = 1,N
               DL(K) = S*DL(K) + C*V(K,LLP1)
170         CONTINUE
         ENDIF
C
C
         TEM = R0NRM*PROD
         DO 175 I = 1,N
           DL(I) = TEM*DL(I)
  175    CONTINUE
C
      ENDIF
C
C
C
 200  CONTINUE
C
      LL   = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
         R0(K) = 0.0d0
 210  CONTINUE
C
      R0(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, R0)
      DO 220 K = 1,N
         Z(K) = 0.0d0
 220  CONTINUE
C
      DO 230 I = 1,LL
          DO 228 J = 1,N
             Z(J) = Z(J) + R0(I)*V(J,I)
  228     CONTINUE
 230  CONTINUE
C
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 240 I = 1,N
            Z(I) = Z(I)/SZ(I)
 240     CONTINUE
      ENDIF
C
      IF (JPRE .GT. 0) THEN
         DO 250 I = 1,N
            WK(I) = Z(I)
  250    CONTINUE         
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
C
 1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Natral Err Est','   Error Estimate')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7)
 1020 FORMAT(1X,' ITER = ',I5, ' IELMAX = ',I5,
     $     ' |R(IELMAX)/X(IELMAX)| = ',E12.5)
C
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSMV( N, X, Y, NELT, IA, JA, A, ISYM )
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM
      REAL*8 A(NELT), X(N), Y(N)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSMV
C
C
      DO 10 I = 1, N
         Y(I) = 0.0d0
 10   CONTINUE
C
      DO 30 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 20 I = IBGN, IEND
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20      CONTINUE
 30   CONTINUE
C
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSMTV( N, X, Y, NELT, IA, JA, A, ISYM )
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM
      REAL*8 X(N), Y(N), A(NELT)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSMTV
C
C
      DO 10 I = 1, N
         Y(I) = 0.0d0
 10   CONTINUE
C
      DO 30 IROW = 1, N
         IBGN = JA(IROW)
         IEND = JA(IROW+1)-1
         DO 20 I = IBGN, IEND
            Y(IROW) = Y(IROW) + A(I)*X(IA(I))
 20      CONTINUE
 30   CONTINUE
C
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSLUI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(*)
      REAL*8 B(N), X(N), A(NELT), RWORK(*)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSLUI
C
C
      LOCIL  = IWORK(1)
      LOCJL  = IWORK(2)
      LOCIU  = IWORK(3)
      LOCJU  = IWORK(4)
      LOCL   = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU   = IWORK(7)
C         
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE
C
      DO 30 IROW = 2, N
         JBGN = IWORK(LOCIL+IROW-1)
         JEND = IWORK(LOCIL+IROW)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               X(IROW) = X(IROW)-RWORK(LOCL+J-1)*X(IWORK(LOCJL+J-1))
 20         CONTINUE
         ENDIF
 30   CONTINUE
C         
      DO 40 I=1,N
         X(I) = X(I)*RWORK(LOCDIN+I-1)
 40   CONTINUE
C         
      DO 60 ICOL = N, 2, -1
         JBGN = IWORK(LOCJU+ICOL-1)
         JEND = IWORK(LOCJU+ICOL)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               JJUU    = IWORK(LOCIU+J-1)
               X(JJUU) = X(JJUU)-RWORK(LOCU+J-1)*X(ICOL)
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSLUTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(*)
      REAL*8 B(N), X(N), A(N), RWORK(*)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSLUTI
C
C
      LOCIL  = IWORK(1)
      LOCJL  = IWORK(2)
      LOCIU  = IWORK(3)
      LOCJU  = IWORK(4)
      LOCL   = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU   = IWORK(7)
C
C         
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
C         
      DO 80 IROW = 2, N
         JBGN = IWORK(LOCJU+IROW-1)
         JEND = IWORK(LOCJU+IROW)-1
         IF( JBGN.LE.JEND ) THEN
            DO 70 J = JBGN, JEND
               X(IROW) = X(IROW)-RWORK(LOCU+J-1)*X(IWORK(LOCIU+J-1))
 70         CONTINUE
         ENDIF
 80   CONTINUE
C         
      DO 90 I = 1, N
         X(I) = X(I)*RWORK(LOCDIN+I-1)
 90   CONTINUE
C         
      DO 110 ICOL = N, 2, -1
         JBGN = IWORK(LOCIL+ICOL-1)
         JEND = IWORK(LOCIL+ICOL)-1
         IF( JBGN.LE.JEND ) THEN
            DO 100 J = JBGN, JEND
               JJUU    = IWORK(LOCJL+J-1)
               X(JJUU) = X(JJUU)-RWORK(LOCL+J-1)*X(ICOL)
 100        CONTINUE
         ENDIF
 110  CONTINUE
C
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSILUS(N,NELT,IA,JA,A,ISYM,NL,IL,JL,
     &                  L,DINV,NU,IU,JU,U,NROW,NCOL)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
c changes below for single grid block or separate column problems
cels6/13/08 orig     integer*8 N, NELT,IA(NELT),JA(NELT),ISYM,NL,IL(NL),JL(NL)
cels6/13/08 els, stefans below  integer*8 N, NELT,IA(NELT),JA(NELT),ISYM,NL,IL(N+1),JL(NL)
cns      integer*8 N, NELT,IA(NELT),JA(NELT),ISYM,NL,IL(NL+2),JL(NL+2)
      integer*8 N, NELT,IA(NELT),JA(NELT),ISYM,NL,IL(N+1),JL(NL+2)
cels6/13/08 orig      integer*8 NU, IU(NU), JU(NU), NROW(N), NCOL(N)
cels6/13/08 els, stefans below   integer*8 NU, IU(NU), JU(N+1), NROW(N), NCOL(N)
cns      integer*8 NU, IU(NU+2), JU(NU+2), NROW(N), NCOL(N)
      integer*8 NU, IU(NU+2), JU(N+1), NROW(N), NCOL(N)
      REAL*8 A(NELT), L(NL), DINV(N), U(NU)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSILUS
C
C
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
C         
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
C
      DO 300 I=2,N
C         
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
C         
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
C         
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
C         
 300  CONTINUE
C         
      DO 430 I=1,N
         DINV(I) = 1.d0/DINV(I)
 430  CONTINUE
C         
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DCHKW( NAME, LOCIW, LENIW, LOCW, LENW,
     &     IERR, ITER, ERR )
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      CHARACTER*(*) NAME
      CHARACTER*72 MESG
      integer*8 LOCIW, LENIW, LOCW, LENW, IERR, ITER
      REAL*8 ERR
C
      REAL*8 DMACH(5)
      DATA DMACH(2) / 1.79769313486231D+308 /
c      DATA DMACH(2) / 1.79769313486231D+38 /
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DCHKW
C
C
      IERR = 0
      IF(LOCIW.GT.LENIW) THEN
         IERR = 1
         ITER = 0
         ERR = DMACH(2)
         MESG = NAME // ': integer work array too short. '//
     $        ' IWORK needs i1: have allocated i2.'
      ENDIF
C
      IF(LOCW.GT.LENW) THEN
         IERR = 1
         ITER = 0
         ERR = DMACH(2)
         MESG = NAME // ': REAL work array too short. '//
     $        ' RWORK needs i1: have allocated i2.'
      ENDIF
C      
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE QS2I1D( IA, JA, A, N, KFLAG, IUNIT )
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 IL(21),IU(21)
      integer*8 kflag,n,iunit,kk
      integer*8  IA(N),JA(N),IT,IIT,JT,JJT
      REAL*8 A(N), TA, TTA
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of QS2I1D
C
C
      NN = N
      IF (NN.LT.1) THEN
         WRITE(IUNIT,6100) 
 6100 FORMAT(/,'QS2I1D- the number of values to',  
     &         ' be sorted was NOT POSITIVE.')
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         WRITE(IUNIT,6101) 
 6101 FORMAT(/,'QS2I1D- the sort control parameter, k, ',
     &         'was not 1 OR -1.')
         RETURN
      ENDIF
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
      M = 1
      I = 1
      J = NN
      R = 3.75D-1
 210  IF( R.LE.5.898437D-1 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-2.1875D-1
      ENDIF
 225  K = I
C
C
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
C
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
C                           
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
C
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
C
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
C
C
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
C
C                                  
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
C
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DS2Y(N, NELT, IA, JA, A, ISYM,IUNIT )
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 N, NELT, IA(NELT), JA(NELT), ISYM, iunit
      REAL*8 A(NELT)
      parameter (ione = 1)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DS2Y
C
C
      IF( JA(N+1).EQ.NELT+1 ) RETURN
      CALL QS2I1D( JA, IA, A, NELT, ione,IUNIT )
      JA(1) = 1
      DO 20 ICOL = 1, N-1
         DO 10 J = JA(ICOL)+1, NELT
            IF( JA(J).NE.ICOL ) THEN
               JA(ICOL+1) = J
               GOTO 20
            ENDIF
 10      CONTINUE
 20   CONTINUE
      JA(N+1) = NELT+1
C         
      JA(N+2) = 0
C
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
      RETURN
      END
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
C      D I R E C T    S O L V E R S
C
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
C
      SUBROUTINE DGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 IPIV(*)
      REAL*8 AB(LDAB,*)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
      PARAMETER (NBMAX = 32,LDWORK = NBMAX+1 )
      REAL*8 WORK13(LDWORK,NBMAX),WORK31(LDWORK,NBMAX)
c     added parameter for passing same precision constant
      parameter (ione=1)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGBTRF
C
C
      KV = KU + KL
      IF(KU.LE.64) THEN
         NB = 1
      ELSE
         NB = 32
      END IF
C
      NB = MIN(NB,NBMAX)
      IF(NB.LE.1.OR.NB.GT.KL) THEN
         CALL DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
      ELSE
         DO 20 J = 1,NB
            DO 10 I = 1,J-1
               WORK13(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
C
         DO 40 J = 1,NB
            DO 30 I = J+1,NB
               WORK31(I,J) = ZERO
   30       CONTINUE
   40    CONTINUE
C
C
C
         JU = 1
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
            DO 80 JJ = J, J + JB - 1
               KM = MIN( KL, M-JJ )
C               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               dmax  = abs(AB(KV+1,JJ))
               idamx = 1
               do 31 i8 = 2,KM+2
                  if(i8.gt.km+1) go to 32
                  if(abs(AB(KV+i8,JJ)).le.dmax) go to 31
                  idamx = i8
                  dmax  = abs(AB(KV+i8,JJ))
   31          continue
   32          jp = idamx
               IPIV( JJ ) = JP + JJ - J
               IF( AB(KV+JP,JJ).NE.ZERO ) THEN
                  JU = MAX(JU,MIN(JJ+KU+JP-1,N))
                  IF(JP.NE.1) THEN
                     IF( JP+JJ-1.LT.J+KL ) THEN
                        CALL DSWAP(JB,AB(KV+1+JJ-J,J),LDAB-1,
     $                             AB(KV+JP+JJ-J,J),LDAB-1)
                     ELSE
                        CALL DSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,
     $                             WORK31(JP+JJ-J-KL,1),LDWORK)
                        CALL DSWAP(J+JB-JJ,AB(KV+1,JJ),LDAB-1,
     $                             AB(KV+JP,JJ),LDAB-1)
                     END IF
                  END IF
C
                  u4 = ONE/AB(KV+1,JJ)
                  do 78 i4=1,km
                     ab(KV+1+i4,JJ)=u4*ab(KV+1+i4,JJ)
   78             continue
C
C                  CALL DSCAL(KM,ONE/AB(KV+1,JJ),AB(KV+2,JJ),1)
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
cels4/24/06 modified for passing same precision
cels4/24/06     $               CALL DGER(KM,JM-JJ,-ONE,AB(KV+2,JJ),1,
     $               CALL DGER(KM,JM-JJ,-ONE,AB(KV+2,JJ),ione,
     $                         AB(KV,JJ+1),LDAB-1,
     $                         AB(KV+1,JJ+1),LDAB-1)
               ELSE
                  IF(INFO.EQ.0) INFO = JJ
               END IF
C
               NW = MIN( JJ-J+1, I3 )
               IF(NW.GT.0) then
                  do 79 i6=1,nw
                     work31(i6,jj-j+1) = ab(kv+kl-jj+j+i6,jj)
   79             continue
               end if
   80       CONTINUE
C
C
C
            IF( J+JB.LE.N ) THEN
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
cels4/24/06 modified to pass same precision constant
cels4/24/06               CALL DLASWP(J2,AB(KV+1-JB,J+JB),LDAB-1,1,JB,IPIV(J),1)
         CALL DLASWP(J2,AB(KV+1-JB,J+JB),LDAB-1,ione,JB,IPIV(J),ione)
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
C
C
C
               IF( J2.GT.0 ) THEN
                  CALL DTRSM( JB, J2, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
*
                  IF( I2.GT.0 ) THEN
                     CALL DGEMM( I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, 
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
                     CALL DGEMM( I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, 
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
C
C
C
               IF( J3.GT.0 ) THEN
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
                  CALL DTRSM( JB, J3, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
*
                  IF( I2.GT.0 ) THEN
                     CALL DGEMM( I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
                     CALL DGEMM( I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
C
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
C
C
C
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
                  IF( JP+JJ-1.LT.J+KL ) THEN
                     CALL DSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,
     $                          AB(KV+JP+JJ-J,J),LDAB-1)
                  ELSE
                     CALL DSWAP(JJ-J,AB(KV+1+JJ-J,J),LDAB-1,
     $                          WORK31(JP+JJ-J-KL,1),LDWORK)
                  END IF
               END IF
               NW = MIN( I3, JJ-J+1 )
               IF(NW.GT.0) then
                  do 169 i6=1,nw
                     ab(kv+kl-jj+j+i6,jj) = work31(i6,jj-j+1)
  169             continue
               end if
  170       CONTINUE
  180    CONTINUE
      END IF
C
C
C
      RETURN
      END
C
C
C
C***********************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***********************************************************************
C
C
C
      SUBROUTINE DGBTF2(M,N,KL,KU,AB,LDAB,IPIV,INFO)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 IPIV(*)
      REAL*8 AB(LDAB,*)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
cels4/24/06 added for same precision constant
      parameter (ione=1)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGBTF2
C
C
      KV = KU + KL
C
      JU = 1
      DO 100 J = 1, MIN( M, N )
         KM = MIN(KL,M-J)
         dmax  = abs(AB(KV+1,J))
         idamx = 1
         do 10 i8 = 2,KM+2
            if(i8.gt.km+1) go to 11
            if(abs(AB(KV+i8,J)).le.dmax) go to 10
            idamx = i8
            dmax  = abs(AB(KV+i8,J))
   10    continue
   11    jp = idamx
C
         IPIV( J ) = JP+J-1
         IF(AB(KV+JP,J).NE.ZERO) THEN
C
            JU = MAX(JU,MIN(J+KU+JP-1,N))
            IF(JP.NE.1) then
               do 20 i6=1,(JU-J+1)*(LDAB-1),LDAB-1
                  ddttmm           = ab(KV+JP-1+i6,J)
                  ab(KV+JP-1+i6,J) = ab(KV+i6,J)
                  ab(KV+i6,J)      = ddttmm
   20          continue
            end if
C
            IF( KM.GT.0 ) THEN
               u4 = ONE/AB(KV+1,J)
               do 40 i4=1,km
                  ab(KV+1+i4,J)=u4*ab(KV+1+i4,J)
   40          continue
               IF( JU.GT.J )
cels4/24/06 added for same precision constant
cels4/24/06     $            CALL DGER(KM,JU-J,-ONE,AB(KV+2,J),1,AB(KV,J+1),
     $            CALL DGER(KM,JU-J,-ONE,AB(KV+2,J),ione,AB(KV,J+1),
     $                      LDAB-1,AB(KV+1,J+1),LDAB-1)
            END IF
C
         ELSE
            IF( INFO.EQ.0 ) INFO = J
         END IF
  100 CONTINUE
C
      RETURN
      END
C
C
C
C***********************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***********************************************************************
C
C
C
      SUBROUTINE DGBTRS(N,KL,KU,AB,LDAB,IPIV,B,LDB)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 IPIV(*)
      REAL*8   AB(LDAB,*), B(LDB,1)
      PARAMETER (ONE = 1.0d0)
cels4/24/06 added for passing same precision
      parameter (ione = 1)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGBTRS
C
C
      KD = KU+KL+1
      DO 10 J = 1,N-1
         LM = MIN(KL,N-J)
         L = IPIV(J)
cels4/24/06 modified for passing same precision
cels4/24/06         IF(L.NE.J) CALL DSWAP(1,B(L,1),LDB,B(J,1),LDB)
cels4/24/06         CALL DGER(LM,1,-ONE,AB(KD+1,J),1,B(J,1),LDB,B(J+1,1),LDB)
        IF(L.NE.J) CALL DSWAP(ione,B(L,1),LDB,B(J,1),LDB)
        CALL DGER(LM,ione,-ONE,AB(KD+1,J),ione,B(J,1),LDB,B(J+1,1),LDB)
   10 CONTINUE
C
      CALL DTBSV(N,KL+KU,AB,LDAB,B(1,1))
C
C
      RETURN
      END
C
C
C
C***********************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***********************************************************************
C
C
C
      SUBROUTINE DLASWP(N,A,LDA,K1,K2,IPIV,INCX)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 IPIV(*)
      REAL*8   A(LDA,*)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DLASWP
C
C
      DO 10 I = K1, K2
         IP = IPIV(I)
         IF(IP.NE.I) CALL DSWAP(N,A(I,1),LDA,A(IP,1),LDA)
   10 CONTINUE
C
C
C
      RETURN
      END
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 m,n,incx,incy,lda
      double precision alpha
      REAL*8 A(LDA,*),X(*),Y(*)
      PARAMETER ( ZERO = 0.0d0 )
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGER
C
C
      INFO = 0
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
C
C
      RETURN
      END
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DTBSV (N,K,A,LDA,X)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 n,k,lda
      REAL*8 A(LDA,*),X(*)
      PARAMETER (ZERO = 0.0d0)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DTBSV
C
C
      INFO = 0
      KPLUS1 = K + 1
      DO 20, J = N, 1, -1
         IF(X(J).NE.ZERO)THEN
            L = KPLUS1 - J
            X(J) = X(J)/A(KPLUS1,J)
            TEMP = X(J)
            DO 10, I = J-1,MAX(1,J-K),-1
               X(I) = X(I) - TEMP*A(L+I,J)
   10       CONTINUE
         END IF
   20 CONTINUE
C
C
C
      RETURN
      END
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DGEMM(M,N,K,ALPHA,A,LDA,B,LDB,C,LDC)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 m,n,k,lda,ldb,ldc
      double precision alpha
      REAL*8 A(LDA,*), B(LDB,*),C(LDC,*)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DGEMM
C
C
      DO 100 J=1,N
         DO 80 L=1,K
            IF(B(L,J).NE.ZERO)THEN
               TEMP = ALPHA*B(L,J)
               DO 60 I=1,M
                  C(I,J) = C(I,J)+TEMP*A(I,L)
   60          CONTINUE
            END IF
   80    CONTINUE
  100 CONTINUE
C
C
C
      RETURN
      END
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DTRSM (M,N,A,LDA,B,LDB)
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 m,n,lda,ldb
      REAL*8 A(LDA,*),B(LDB,*)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DTRSM
C
C
      DO 100 J=1,N
         DO 90 K=1,M
            IF(B(K,J).NE.ZERO)THEN
               DO 80 I=K+1,M
                  B(I,J) = B(I,J)-B(K,J)*A(I,K)
   80          CONTINUE
            END IF
   90    CONTINUE
  100 CONTINUE
C
C
C
      RETURN
      END
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C
C
      SUBROUTINE DSWAP(n,dx,incx,dy,incy)
C
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      integer*8 n,incx,incy
      REAL*8 dx(1),dy(1),dtemp
C
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of DSWAP
C
C
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1) go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
      ix = 1
      iy = 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
C
   20 do 50 i = 1,n
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   50 continue
C
C
C
      return
      end
C
C
C***************************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***************************************************************************
C