C
      SUBROUTINE IO  ! inquire and read data from disk files
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
c
      COMMON/SVZ/NOITE,MOP(24)
      LOGICAL EX
C
      SAVE ICALL
      DATA ICALL/0/
      WRITE(34,892)
  892 FORMAT(1X,60('*'),'call subroutine IO in t2cg22_v2.f',60('*'))
      ICALL=ICALL+1  ! number of times this subroutine is called
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' IO 1.0, 1991.4.15: READ DISK FILES:'/)
C
      OPEN(11,FILE='TOFT',STATUS='UNKNOWN')

  176 INQUIRE(FILE='MESH',EXIST=EX)
      IF(EX) GOTO 2
      WRITE (34,3)
    3 FORMAT(' FILE *MESH* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(4,FILE='MESH',STATUS='NEW')
      GOTO 10
C
    2 WRITE (34,4)
    4 FORMAT(' FILE *MESH* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(4,FILE='MESH',STATUS='OLD')
C
   10 INQUIRE(FILE='INCON',EXIST=EX)
      IF(EX) GOTO 12
      WRITE (34,13)
   13 FORMAT(' FILE *INCON* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(1,FILE='INCON',STATUS='NEW')
      ENDFILE 1
      GOTO 20
C
   12 WRITE (34,14)
   14 FORMAT(' FILE *INCON* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(1,FILE='INCON',STATUS='OLD')
C
   20 INQUIRE(FILE='GENER',EXIST=EX)
      IF(EX) GOTO 22
      WRITE (34,23)
   23 FORMAT(' FILE *GENER* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(3,FILE='GENER',STATUS='NEW')
      ENDFILE 3
      GOTO 30
C
   22 WRITE (34,24)
   24 FORMAT(' FILE *GENER* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(3,FILE='GENER',STATUS='OLD')
C
   30 INQUIRE(FILE='SAVE',EXIST=EX)
      IF(EX) GOTO 32
      WRITE (34,33)
   33 FORMAT(' FILE *SAVE* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(7,FILE='SAVE',STATUS='NEW')
      GOTO 40
C
   32 WRITE (34,34)
   34 FORMAT(' FILE *SAVE* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(7,FILE='SAVE',STATUS='OLD')
C
   40 INQUIRE(FILE='LINEQ',EXIST=EX)
      IF(EX) GOTO 42
      WRITE (34,43)
   43 FORMAT(' FILE *LINEQ* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(15,FILE='LINEQ',STATUS='NEW')
      GOTO 50
C
   42 WRITE (34,44)
   44 FORMAT(' FILE *LINEQ* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(15,FILE='LINEQ',STATUS='OLD')
      REWIND 15
C
   50 CONTINUE
      IF(MOP(15).EQ.0) GOTO 60
      INQUIRE(FILE='TABLE',EXIST=EX)
      IF(EX) GOTO 52
      WRITE (34,53)
   53 FORMAT(' FILE *TABLE* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(8,FILE='TABLE',STATUS='NEW')
      ENDFILE 8
      GOTO 51
C
   52 WRITE (34,54)
   54 FORMAT(' FILE *TABLE* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(8,FILE='TABLE',STATUS='OLD')

   51 INQUIRE(FILE='STIME',EXIST=EX)
      IF(EX) GOTO 56
      WRITE (34,55)
   55 FORMAT(' FILE *STIME* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(111,FILE='STIME',STATUS='NEW')
      ENDFILE 111
      GOTO 60
C
   56 WRITE (34,57)
   57 FORMAT(' FILE *STIME* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(111,FILE='STIME',STATUS='OLD')
C
   60 RETURN
      END
C
      SUBROUTINE FLOP
C
C-----CALCULATE NUMBER OF SIGNIFICANT DIGITS FOR FLOATING POINT
C     PROCESSING. ASSIGN DEFAULT FOR DFAC, AND PRINT APPROPRIATE
C     WARNING WHEN MACHINE ACCURACY IS INSUFFICIENT.
c     
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
c
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      WRITE(34,892)
  892 FORMAT(/1X,60('*'),'call subroutine flop in t2cg22_v2.f',60('*'))
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' FLOP 1.0, 1991.4.11: CALCULATE NUMBER OF SIGNIFICANT'
     X' DIGITS FOR FLOATING POINT ARITHMETIC')
C
      A=SQRT(0.99D0)
      B=A
      DO 1 N=1,260
      B=B/2.D0
      C=A+B
      D=C-A
      IF(D.EQ.0.D0) GOTO 2
    1 CONTINUE
C
    2 B2=B*2.0D0
      N10=-INT(LOG10(B2))
      DF=SQRT(B2)
C
      WRITE (34,3) N10
    3 FORMAT(' FLOATING POINT PROCESSOR HAS APPROXIMATELY',I3,
     X' SIGNIFICANT DIGITS')
C
      IF(DFAC.EQ.0.D0) THEN
         DFAC=DF
         REWIND 1
         WRITE(34,16) DFAC
   16    FORMAT(' USE DEFAULT DFAC=', E10.4)
      ELSE
         WRITE(34,10) DFAC
   10    FORMAT(' USE DFAC=', E10.4)
      END IF
      IF(N10.LE.12.AND.N10.GT.8) WRITE (34,4)
    4 FORMAT(' WARNING: NUMBER OF SIGNIFICANT DIGITS IS MARGINAL;',
     X' EXPECT DETERIORATED CONVERGENCE BEHAVIOR'/)
      IF(N10.LE.8) WRITE (34,5)
    5 FORMAT(' WARNING: NUMBER OF SIGNIFICANT DIGITS IS INSUFFICIENT',
     X'; CONVERGENCE WILL BE POOR OR FAIL'/' CODE SHOULD BE RUN IN',
     X' DOUBLE PRECISION!'/)
C
      RETURN
      END
C
      subroutine sinsub
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/SVZ/NOITE,MOP(24) 
      COMMON/BC/NELA
      common/ff/h1
      character*1 h1
      COMMON/SOLVR1/matslv,nmaxit,nnvvcc,iiuunn,iissoo,nactdi
      COMMON/SOLVR2/ritmax,closur
      COMMON/SOLVR3/ordrng,oprocs,zprocs,coord
      CHARACTER*2 ordrng,oprocs,zprocs
      CHARACTER*5 coord
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' SIN 1.00, 1999.10.1: initialize parameters for the'
     x' solver package, and generate informative printout')
C
c.....imported from INPUT of t2fa.f (GM)
C
C     SETTING SOLVER TYPE AND CORRESPONDING PARAMETERS          
C
      ordrng = 'ST'
      IF(iissoo.eq.0) THEN
         matslv = mop(21)
         zprocs = 'Z1'
            if(matslv.eq.6) zprocs='Z0'
         oprocs = 'O0'
         ritmax = 1.0d-1
         closur = 1.0d-6
      END IF
C
      IF(matslv.gt.6) THEN
         matslv = 3
      ENDIF
C
      IF(matslv.eq.6) GOTO 3085
         IF(ritmax.EQ.0.0d0.or.ritmax.GT.1.0d0) ritmax = 1.0d-1
         IF((abs(closur)).gt.1.0d-6)  closur = 1.0d-6
         IF((abs(closur)).lt.1.0d-12) closur = 1.0d-12
      riter  = NELA*neq*ritmax
      nmaxit = MAX(20,INT(riter))
 3085 continue
c
      if(zprocs.ne.'Z0') then
         if(matslv.eq.6) then
            write (34,6101) matslv,zprocs
            zprocs='Z0'
         endif
      endif
c
      if(oprocs.ne.'O0') then
         if(matslv.eq.6) then
            write (34,6036) matslv,oprocs
            oprocs='O0'
         endif
      endif
c
         IF(zprocs.NE.'Z0'.AND.zprocs.NE.'Z1'.AND.
     &      zprocs.NE.'Z2'.AND.zprocs.NE.'Z3'.AND.
     &      zprocs.NE.'Z4') THEN
               write (34,6102) zprocs
               zprocs = 'Z1'
         END IF
c
         IF(oprocs.NE.'O0'.AND.oprocs.NE.'O1'.AND.
     &      oprocs.NE.'O2'.AND.oprocs.NE.'O3'.AND.
     &      oprocs.NE.'O4') THEN
               write (34,6037) oprocs
               oprocs = 'O0'
         END IF
c
      if(neq.eq.1) then
c.....no preprocessing for NEQ=1
         write (34,6200)
         zprocs='Z0'
         oprocs='O0'
      endif
c
 6101 FORMAT(' WARNING: For MATSLV=',I1,', no Z-preprocessing can be'
     &' used, reset ZPROCS = Z0 and continue execution')
 6036 FORMAT(' WARNING: For MATSLV=',I1,', no O-preprocessing can be'
     &' used, reset OPROCS = O0 and continue execution')
 6102 FORMAT(' WARNING: Unknown matrix preprocessing option ZPROCS = '
     X,a2,', reset ZPROCS = Z1 and continue execution')
 6037 FORMAT(' WARNING: Unknown matrix preprocessing option OPROCS = '
     x,a2,', reset OPROCS = O0 and continue execution')
 6200 format(' NEQ=1; do not perform any matrix preprocessing'/)
c
      return
      end
C
      SUBROUTINE LINEQ
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
C
C-----THIS SUBROUTINE CALLS THE LINEAR EQUATION SOLVER T2CG2.
C     IT HAS LOGIC TO HANDLE FAILURES IN LINEAR EQUATION SOLUTION.
C
C     AFTER SOLUTION, LATEST UPDATED ITERATES ARE OBTAINED FOR
C     ALL PRIMARY DEPENDENT VARIABLES.
C
C$$$$$$$$$ COMMON BLOCKS FOR LINEAR EQUATIONS $$$$$$$$$$$$$$$$$$$$$$$$$
C
      PARAMETER(seed = 1.0d-25)
      PARAMETER(iunit = 0, nvectr = 30)
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/L4/WKAREA(MNEQ*MNEL+10)
      COMMON/L7/JVECT(niwork)
C
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
C
      common/soll/lenw,leniw
C
C     arrays used by luband only.
      COMMON/lub1/AB(nrwork)
      COMMON/lub3/NSUPDI,NSUBDI,mnzp1,mnetp1,mnelp1,nnnbig
      COMMON/lub4/matord,nsubdg,nsupdg,ntotd
C
      COMMON/E1/ELEM(MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P4/R(MNEQ*MNEL+1)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/DG/WUP,WNR
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/BC/NELA
C
      CHARACTER*5 ELEM
      CHARACTER*2 ordrng,oprocs,zprocs
      CHARACTER*5 coord
C
      COMMON/SOLVR1/matslv,nmaxit,nnvvcc,iiuunn,iissoo,nactdi
      COMMON/SOLVR2/ritmax,closur
      COMMON/SOLVR3/ordrng,oprocs,zprocs,coord
C
      SAVE ICALL,N,iteruc,iprpro,izero0
      DATA ICALL,iteruc,iprpro,izero0/0,0,0,0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' LINEQ 2.00, 1999.10.4: Interface for linear equation'
     &' solvers T2CG2, Can call a direct solver or a package of'
     &' conjugate gradient solvers')
C
C =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of LINEQ
      MA = MA+1
      IF(MA.EQ.1) N=NEQ*NELA
C
C*********************************************************************
C*                                                                   *
C*             ACCOUNTING FOR ZEROs ON THE MAIN DIAGONAL             *
C*                                                                   *
C*********************************************************************
C
      izerod = 0
      iprpro = 0
      IF(MATSLV.NE.6) THEN
         NNF = N*NEQ
         DO 20 I=1,NNF
         IF(IRN(I).EQ.ICN(I).AND.ABS(CO(I)).EQ.0.0d0) izerod = izerod+1
   20    CONTINUE
      END IF
C
      if(izerod.ne.izero0) then
         WRITE(15,6003) kcyc,iter,izerod,zprocs,oprocs
         izero0 = izerod
      endif
C
cels03/7/2      zertio=1.0d2*izerod/(neq*nela)
      zertio=1.0d2*dble(izerod/(neq*nela))
      IF(izerod.GT.0) THEN
         IF(zprocs.EQ.'Z1') iprpro = 1
         IF(zprocs.EQ.'Z2') iprpro = 2
         IF(zprocs.EQ.'Z3') iprpro = 3
         IF(zprocs.EQ.'Z4') iprpro = 4
         IF(oprocs.NE.'O0') THEN
            IF(zprocs.EQ.'Z0'.OR.zprocs.EQ.'Z1') THEN
               IF(zertio.GT.2.0d+1) THEN
cels02/5/3                  PRINT 6004, zertio,oprocs
                  write (34,6004) zertio,oprocs
                  oprocs='O0'
               END IF
            END IF
         END IF
      END IF
C
 6003 FORMAT(' At KCYC=',i5,' and ITER=',i5,', IZEROD=',I5,', ZPROCS ='
     &,a2,' and OPROCS = ',a2)
 6004 FORMAT(/'WARNING-',F5.2,'% of the elements on the main diagonal'
     &' of the Jacobian matrix are zeros, the matrix preprocessor'
     &' OPROCS = ',A2,' cannot be used, reset OPROCS = O0 and',
     &' continue execution.')
C
cels10/5/09      INUM  = 0
      IGOOD = 0
C
C-----MA COUNTS CALLS TO LINEQ.
C
      IF(MOP(6).NE.0) write (34,6005) kcyc,iter,icall,N,NZ,izerod,zertio
 6005 FORMAT(' LINEQ at [KCYC, ITER] = [',I5,',',I3,']',' ICALL =',I5,
     &' N =',I4,' NZ =',I8,' IZEROD =',I7,' or',F6.2,' % zeros')
C
      IF(MOP(6).GE.7) THEN
         write (34,6015)
         write (34,6010) (IRN(NN),ICN(NN),CO(NN),NN=1,NZ)
      END IF
C
 6010 FORMAT(5(1X,2I5,E14.6))
 6015 FORMAT(/' MATRIX OF COEFFICIENTS'/)
C
C*********************************************************************
C*                                                                   *
C*      DETERMINATION OF THE BANDWIDTHS OF THE U AND L MATRICES      *
C*      OR PLACEMENT OF THE ELEMENTS INTO THE CG SOLUTION ARRAY      *
C*                                                                   *
C*********************************************************************
C
      IF(icall.EQ.1) THEN
         co(mnzp1) = 0.0d0
         r(mnetp1) = 0.0d0
      END IF
C
      IF(MATSLV.EQ.6.AND.icall.EQ.1) THEN
         iddfup = 0
         iddfdn = 0
         DO 400 j=1,nz
            iddffu = icn(j)-irn(j)
            iddffd = irn(j)-icn(j)
            IF(iddffu.GE.iddfup) iddfup=iddffu
            IF(iddffd.GE.iddfdn) iddfdn=iddffd
  400    CONTINUE
C
         matord = neq*nela
         nsupdg = iddfup
         nsubdg = iddfdn
         ntotd  = 2*nsubdg+nsupdg+1
         navdia = lenw/matord
C
         IF(ntotd.GT.navdia) THEN
            write (34,6020) ntotd,navdia,matord*ntotd
            STOP
         END IF
      END IF
C
C
C
      IF(icall.EQ.1.AND.NEQ.GT.1.AND.MATSLV.NE.6) THEN
         IF(oprocs.NE.'O0'.OR.(izerod.NE.0.AND.iprpro.GT.1)) THEN
            CALL ELINDX
            CALL REASSN
         END IF
      END IF
C
C***********************************************************************
C*                                                                     *
C*                          MATRIX SOLUTION                            *
C*                                                                     *
C***********************************************************************
C
      IF(MATSLV.EQ.6) THEN
         CALL LUBAND(matord,nz,nsubdg,nsupdg,ntotd,ab,
     &               matord,r,JVECT,info)
C
      ELSE
C
         IF(izerod.NE.0.AND.NEQ.GT.1.AND.zprocs.NE.'Z0') THEN
            CALL MTRXIN(iprpro)
         END IF
c
         IF(oprocs.NE.'O0'.AND.NEQ.EQ.1) GO TO 415
         IF(oprocs.eq.'O0') THEN
            GO TO 415
         ELSE IF(oprocs.eq.'O1') THEN
cels4/20/06 don't pass constants
cels4/20/06            CALL MTRXPR(1)
            ione = 1
            CALL MTRXPR(ione)
         ELSE IF(oprocs.eq.'O2') THEN
cels4/20/06             CALL MTRXPR(2)
            itwo = 2
            CALL MTRXPR(itwo)
         ELSE IF(oprocs.eq.'O3') THEN
cels4/20/06            CALL MTRXPR(3)
            ithree = 3
            CALL MTRXPR(ithree)
         ELSE IF(oprocs.eq.'O4') THEN
cels4/20/06            CALL MTRXIN(4)
            ifour = 4
            CALL MTRXIN(ifour)
         END IF
C
  415    CONTINUE
C
         DO 440 i=1,n
            wkarea(i) = 0.0d0
  440    CONTINUE
C
         IF(MOP(6).NE.0) CALL THYME(0,TS,TT)
C
         IF(MATSLV.EQ.2) THEN
            CALL DSLUBC(N,r,wkarea,NZ,irn,icn,co,
     &                  CLOSUR,NMAXIT,ITERU,ERR,IERR,IUNIT,
     &                  AB,LENW,jvect,LENIW)
         ELSE IF(MATSLV.EQ.3.or.matslv.eq.1) THEN
            CALL DSLUCS(N,r,wkarea,NZ,irn,icn,co,
     &                  CLOSUR,NMAXIT,ITERU,ERR,IERR,IUNIT,
     &                  AB,LENW,jvect,LENIW)
         ELSE IF(MATSLV.EQ.4) THEN
            CALL DSLUGM(N,r,wkarea,NZ,irn,icn,co,NVECTR,
     &                  CLOSUR,NMAXIT,ITERU,ERR,
     &                  IERR,IUNIT,AB,LENW,jvect,LENIW)
         ELSE IF(MATSLV.EQ.5) THEN
            CALL DLUSTB(N,r,wkarea,NZ,irn,icn,co,
     &                  CLOSUR,NMAXIT,ITERU,ERR,IERR,IUNIT,
     &                  AB,LENW,jvect,LENIW)
         END IF
c
         DO 450 I=1,N
            R(I) = wkarea(I)
  450    CONTINUE
C
         IF(MOP(6).NE.0) THEN
            CALL THYME(1,TSS,TT)
            write (34,6040) TSS
         END IF
C
         iteruc = iteruc+iteru
C
cels5/4/06         WRITE(15,6045) kcyc,iter,deltex,ierr,err,iteru,iteruc
          IF(MOP(6).NE.0)WRITE(15,6045) kcyc,iter,deltex,ierr,
     +        err,iteru,iteruc
      END IF
C
 6020 FORMAT(/'ERROR-SIMULATION ABORTED,DECLARED BANDWIDTH SMALLER THAN'
     &' NEEDED, PLEASE CORRECT AND TRY AGAIN'/'THE NUMBER OF NEEDED'
     &' DIAGONALS IS',I4,' WHILE THE AVAILABLE NUMBER IS',I4,'THE'
     &' PARAMETER nrwork MUST BE INCREASED TO AT LEAST',I10)
 6040 FORMAT(' SOLUTION TIME = ',1PE12.4,' SECONDS')
 6045 FORMAT(' At [',I4,',',I3,']',' DELT=',E12.6,' IERR=',I1,
     &' ERR=',1PE12.6,' IT=',I5,' ITC=',I10)
C
C*********************************************************************
C*                                                                   *
C*                          UPDATE CHANGES                           *
C*                                                                   *
C*********************************************************************
C
  455 IF(MOP(6).GE.5) write (34,6050)
 6050 FORMAT(' ===== INCREMENTS ==== in order of primary variables')
C
      DO 500 NN=1,NELA
         NLOC  = (NN-1)*NEQ
         NLOCP = (NN-1)*NK1
C
         IF(MOP(6).GE.5) write(34,6055) ELEM(NN),(R(NLOC+K),K=1,NEQ)
C
         DO 480 K=1,NEQ
            DX(NLOCP+K) = DX(NLOCP+K)+WNR*R(NLOC+K)
  480    CONTINUE
  500 CONTINUE
C
 6055 FORMAT('       AT ELEMENT *',A5,'*   ',8(1X,E12.6))
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of LINEQ
C
      RETURN
      END
C
      SUBROUTINE LUBAND(N,NZ,KL,KU,LDAB,AB,LDB,B,IPIV8,INFO)
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
C
C***********************************************************************
C*                                                                     *
C*!               C O M M O N    D E C L A R A T I O N S               *
C*                                                                     *
C***********************************************************************
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
C
      COMMON/SVZ/NOITE,MOP(24)
C
      CHARACTER*2 ordrng,oprocs,zprocs
      CHARACTER*5 coord
      COMMON/SOLVR3/ordrng,oprocs,zprocs,coord
C
      integer*8 IPIV8(*)
c     REAL*8 AB(LDAB,*),B(LDB,*)
      Dimension AB(LDAB,*),B(LDB,1)
C
      SAVE ICALL
      DATA ICALL/0/
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of LUBAND
C
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(' LUBAND 1.0, 1997.1.12: Direct banded matrix solver using'
     &' LU decomposition')
C
C
C*********************************************************************
C*                                                                   *
C*                   I N I T I A L I Z A T I O N                     *
C*                                                                   *
C*********************************************************************
C
      info = 0
      DO 502 j=1,n
         DO 502 i=1,LDAB
            AB(i,j) = 0.0D0
  502 CONTINUE
C
C*********************************************************************
C*                                                                   *
C*           PLACEMENT OF MATRIX ELEMENTS INTO THE AB ARRAY          *
C*                                                                   *
C*********************************************************************
C
      DO 503 i=1,NZ
         nfrst = KL+KU+1+irn(i)-icn(i)
         AB(nfrst,icn(i)) = co(i)
  503 CONTINUE
C
      IF(MOP(6).NE.0) CALL THYME(0,TS,TT)
C
C*********************************************************************
C*                                                                   *
C*                         LU DECOMPOSITION                          *
C*                                                                   *
C*********************************************************************
C
      CALL DGBTRF(N,N,KL,KU,AB,LDAB,IPIV8,INFO)
C
      IFLAG=0         !!!! forced
C
      IF(MOP(6).NE.0) THEN
         CALL THYME(1,TD,TT)
         WRITE (34,6001) TD,IFLAG
      END IF
C
 6001 FORMAT(' LU DECOMPOSITION TIME = ',1PE12.4,
     &       '  SECONDS','   **********  INFO =',I4)
C
C*********************************************************************
C*                                                                   *
C*      SOLUTION USING THE LU FACTORIZATION COMPUTED BY DGBTRF       *
C*                                                                   *
C*********************************************************************
C
      IF(info.eq.0) THEN
C
         IF(MOP(6).NE.0) CALL THYME(0,TS,TT)
C
         CALL DGBTRS(N,KL,KU,AB,LDAB,IPIV8,B,LDB)
C
         IF(MOP(6).NE.0) THEN
            CALL THYME(1,TSS,TT)
            WRITE (34,6002) TSS
         END IF
C
      ELSE
         WRITE (34,6003) info
         STOP
      END IF
C
 6002 FORMAT('      SOLUTION TIME = ',1PE12.4,'  SECONDS')
 6003 FORMAT(//,20('ERROR-'),//,T33,
     &             '       S I M U L A T I O N   A B O R T E D',
     &       /,T24,'THE UPPER TRIANGULAR MATRIX ELEMENT U(I,I) ',
     &             'IS EXACTLY ZERO AT I = ',I5,
     &       /,T38,'THE FACTORIZATION IS COMPLETED BUT DIVISION',
     &       /,T38,'BY ZERO WILL OCCUR IF SOLUTION IS ATTEMPTED',
     &       //,20('ERROR-'))         
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of LUBAND
C
      RETURN
      END
C
      SUBROUTINE MTRXPR(ilevel)
C
C
C***********************************************************************
C*                                                                     *
C*                    THE FILE 'flowpar_v2.inc' IS INCLUDED                        *
C*                                                                     *
C***********************************************************************
C
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
c
      INCLUDE 'flowpar_v2.inc'
C
C***********************************************************************
C*                                                                     *
C*!               C O M M O N    D E C L A R A T I O N S               *
C*                                                                     *
C***********************************************************************
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/P4/R(MNEQ*MNEL+1)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
C
      CHARACTER*2 ordrng,oprocs,zprocs
      CHARACTER*5 coord
      COMMON/SOLVR3/ordrng,oprocs,zprocs,coord
      COMMON/SOLVR1/matslv,nmaxit,nnvvcc,iiuunn,iissoo,nactdi
C
      SAVE ICALL
      DATA ICALL/0/
C
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of MTRXPR
C
C
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(6X,'MTRXPR   1.0      19 January   1997',6X,
     &          'Routine for O-preprocessing ',
     &          'of the Jacobian')
C
      MXMYMZ = NELA
      IPHASE = NEQ
 9801 format(t5,'=>=>=> MTRXPR Flag # ',i2)
C
C
C
      IF(IPHASE.EQ.1) GO TO 50
C
C***********************************************************************
C*                                                                     *
C*               ELIMINATION OF LEFT-OFF-M.BAND ELEMENTS               *
C*                                                                     *
C***********************************************************************
C
      DO 45 IJK=1,MXMYMZ
         DO 15 LB=1,IPHASE-1
            LB1=LB+1
            DO 12 LA=LB1,IPHASE
               nb1 = NB0(LB,LA,IJK)
               nb2 = NB0(LB,LB,IJK)
c
               QUOT           = CO(nb1)/CO(nb2)
               R(NX0(LA,IJK)) = R(NX0(LA,IJK))-QUOT*R(NX0(LB,IJK))
C
               DO 4 L=LB1,IPHASE
                  CO(NB0(L,LA,IJK)) = CO(NB0(L,LA,IJK))
     &                               -QUOT*CO(NB0(L,LB,IJK))
    4          CONTINUE
C
               DO 10 L=1,IPHASE
                  IF(IJKMM(ijk).LT.1) GO TO 6
                  DO 5 NDM=1,IJKMM(ijk)
                    NA0LOC               = NDM+ISUMMM(ijk-1)
                    CO(NA0(NA0LOC,L,LA)) = CO(NA0(NA0LOC,L,LA))
     &                                    -CO(NA0(NA0LOC,L,LB))*QUOT
    5             CONTINUE
c
    6             IF(IJKPP(ijk).LT.1) GO TO 10
                  DO 8 NDM=1,IJKPP(ijk)
                    NC0LOC               = NDM+ISUMPP(ijk-1)
                    CO(NC0(NC0LOC,L,LA)) = CO(NC0(NC0LOC,L,LA))
     &                                    -CO(NC0(NC0LOC,L,LB))*QUOT
    8             CONTINUE
   10          CONTINUE
   12       CONTINUE
C
            DO 14 L=LB1,IPHASE
               CO(NB0(LB,L,IJK)) = 0.0d0
   14       CONTINUE
   15    CONTINUE
         IF(ilevel.EQ.1) GO TO 45
C
C***********************************************************************
C*                                                                     *
C*              ELIMINATION OF RIGHT-OFF-M.BAND ELEMENTS               *
C*                                                                     *
C***********************************************************************
C
         DO 25 LB=IPHASE,2,-1
            LB1=LB-1
            DO 23 LA=LB1,1,-1
               QUOT           = CO(NB0(LB,LA,IJK))/CO(NB0(LB,LB,IJK))
               R(NX0(LA,IJK)) = R(NX0(LA,IJK))-QUOT*R(NX0(LB,IJK))
C
               DO 22 L=1,IPHASE
                 IF(IJKMM(ijk).LT.1) GO TO 20
                 DO 18 NDM=1,IJKMM(ijk)
                    NA0LOC               = NDM+ISUMMM(ijk-1)
                    CO(NA0(NA0LOC,L,LA)) = CO(NA0(NA0LOC,L,LA))
     &                                    -CO(NA0(NA0LOC,L,LB))*QUOT
   18            CONTINUE
c
   20            IF(IJKPP(ijk).LT.1) GO TO 22
                 DO 21 NDM=1,IJKPP(ijk)
                    NC0LOC               = NDM+ISUMPP(ijk-1)
                    CO(NC0(NC0LOC,L,LA)) = CO(NC0(NC0LOC,L,LA))
     &                                    -CO(NC0(NC0LOC,L,LB))*QUOT
   21            CONTINUE
   22          CONTINUE
   23       CONTINUE
C
            DO 24 L=LB1,1,-1
               CO(NB0(LB,L,IJK)) = 0.0d0
   24       CONTINUE
   25    CONTINUE
         IF(ilevel.EQ.2) GO TO 45
C
C***********************************************************************
C*                                                                     *
C*                     N O R M A L I Z A T I O N                       *
C*                                                                     *
C***********************************************************************
C
         DO 30 LB=1,IPHASE
            QQQ                =-1.0d0/CO(NB0(LB,LB,IJK))
            CO(NB0(LB,LB,IJK)) =-1.0d0
            R(NX0(LB,IJK))     = R(NX0(LB,IJK))*QQQ
C
            DO 29 L=1,IPHASE
               IF(IJKMM(ijk).LT.1) GO TO 27
               DO 26 NDM=1,IJKMM(ijk)
                 NA0LOC               = NDM+ISUMMM(ijk-1)
                 CO(NA0(NA0LOC,L,LB)) = CO(NA0(NA0LOC,L,LB))*QQQ
   26          CONTINUE
c
   27          IF(IJKPP(ijk).LT.1) GO TO 29
               DO 28 NDM=1,IJKPP(ijk)
                 NC0LOC               = NDM+ISUMPP(ijk-1)
                 CO(NC0(NC0LOC,L,LB)) = CO(NC0(NC0LOC,L,LB))*QQQ
   28          CONTINUE
   29       CONTINUE
   30    CONTINUE
   45 CONTINUE
      GO TO 999
C
C***********************************************************************
C*                                                                     *
C*                            THE NEQ=1 CASE                           *
C*                                                                     *
C***********************************************************************
C
   50 DO 100 IJK=1,MXMYMZ
         QQQ =-1.0d0/CO(NB0(1,1,IJK))
c
         CO(NB0(1,1,IJK)) =-1.0d0
         R(NX0(1,IJK))    = R(NX0(1,IJK))*QQQ
C
         IF(IJKMM(ijk).LT.1) GO TO 56
         DO 55 NDM=1,IJKMM(ijk)
            NA0LOC              = NDM+ISUMMM(ijk-1)
            CO(NA0(NA0LOC,1,1)) = CO(NA0(NA0LOC,1,1))*QQQ
   55    CONTINUE
c
   56    IF(IJKPP(ijk).LT.1) GO TO 100
         DO 60 NDM=1,IJKPP(ijk)
            NC0LOC              = NDM+ISUMPP(ijk-1)
            CO(NC0(NC0LOC,1,1)) = CO(NC0(NC0LOC,1,1))*QQQ
   60    CONTINUE
c
  100 CONTINUE
C
  999 CONTINUE
C
C =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of MTRXPR
C
      RETURN
      END
C
      SUBROUTINE REASSN
C
C***********************************************************************
C*                                                                     *
C*                   THE FILE 'flowpar_v2.inc' IS INCLUDED                         *
C*                                                                     *
C***********************************************************************
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
C***********************************************************************
C*                                                                     *
C*!               C O M M O N    D E C L A R A T I O N S               *
C*                                                                     *
C***********************************************************************
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/L7/JVECT(niwork)
      COMMON/P4/R(MNEQ*MNEL+1)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
C
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
      COMMON/lub3/NSUPDI,NSUBDI,mnzp1,mnetp1,mnelp1,nnnbig
      common/soll/lenw,leniw
C
      SAVE ICALL
      DATA ICALL/0/
C
C =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of REASSN
C
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
!  899 FORMAT(6X,'REASSN   1.0      17 January   1997',6X,
  899 FORMAT(6X,'REASSN   1.1      09 September 2000',6X,
     &          'Establish the sparsity pattern ',
     &          'of the Jacobian for the Z preprocessors')
C
C***********************************************************************
C*                                                                     *
C*            NUMBERING TO MAP CO TO THE A0,C0,B0,X0 SYSTEM            *
C*                                                                     *
C***********************************************************************
C
      DO 100 iee=1,neq*nela
         IJKR  = iee/neq
         IMODR = MOD(iee,neq)
         IF(IMODR.EQ.0) THEN
            IJKMEL = IJKR
            IPH1   = neq
         ELSE
            IJKMEL = IJKR+1
            IPH1   = IMODR
         END IF
         NX0(IPH1,IJKMEL) = iee
  100 CONTINUE
C
C***********************************************************************
C*                                                                     *
C*                           INITIALIZATION                            *
C*                                                                     *
C***********************************************************************
C
      ISUMPP(0)  = 0
      ISUMMM(0)  = 0
      DO 136 ijk=1,MNEL+1
         IJKPP(ijk)  = 0
         ISUMPP(ijk) = 0
c
         IJKMM(ijk)  = 0
         ISUMMM(ijk) = 0
         DO 136 i2=1,neq
            DO 136 i1=1,neq
               NB0(i1,i2,ijk) = mnzp1
  136 CONTINUE
C
C
C
      DO 138 i2=1,neq
         DO 138 i1=1,neq
            DO 138 ic8=1,NCONUP
                  NA0(ic8,i1,i2) = mnzp1
  138 CONTINUE
C
      DO 140 i2=1,neq
         DO 140 i1=1,neq
            DO 140 ic8=1,NCONDN
                  NC0(ic8,i1,i2) = mnzp1
  140 CONTINUE
C
C***********************************************************************
C*                                                                     *
C*      DETERMINING NEIGBORING ELEMENT SPECIFICS ALONG SAME ROW        *
C*                                                                     *
C***********************************************************************
C
      DO 200 n=1,nela
c
         IF(no(iijjkk(n)+1).GT.n) THEN
            IJKPP(n)      = iijjkk(n+1)-iijjkk(n)-1
            jvect(n)      = iijjkk(n)+1
            jvect(nela+n) = iijjkk(n+1)-1
c
            IJKMM(n)        = 0
            jvect(2*nela+n) = 0
            jvect(3*nela+n) = 0
            GO TO 200
         END IF
c
         IF(no(iijjkk(n+1)-1).LT.n) THEN
            IJKPP(n)      = 0
            jvect(n)      = 0
            jvect(nela+n) = 0
c
            IJKMM(n)        = iijjkk(n+1)-iijjkk(n)-1
            jvect(2*nela+n) = iijjkk(n)+1
            jvect(3*nela+n) = iijjkk(n+1)-1
            GO TO 200
         END IF
c
         DO 150 index = iijjkk(n)+1, iijjkk(n+1)-1
            IF(n.GT.no(index).AND.n.LT.no(index+1)) THEN
               IJKPP(n)      = iijjkk(n+1)-1-index
               jvect(n)      = index+1
               jvect(nela+n) = iijjkk(n+1)-1
c
               IJKMM(n)        = index-iijjkk(n)
               jvect(2*nela+n) = iijjkk(n)+1
               jvect(3*nela+n) = index
               GO TO 200
            END IF
  150    CONTINUE
  200 CONTINUE
C
      DO 220 n=1,nela
         ISUMPP(n)  = ISUMPP(n-1)+IJKPP(n)
         ISUMMM(n)  = ISUMMM(n-1)+IJKMM(n)
  220 CONTINUE
C
      IF(ISUMPP(nela).GT.NCONUP) THEN
cels02/5/3         PRINT 6001, ISUMPP(nela)
         write (34,6001) ISUMPP(nela)
         STOP
      END IF
C
      IF(ISUMMM(nela).GT.NCONDN) THEN
cels02/5/3         PRINT 6002, ISUMMM(nela)
         write (34,6002) ISUMMM(nela)
         STOP
      END IF
C
C***********************************************************************
C*                                                                     *
C* MAPPING THE GLOBAL MATRIX ELEMENT NUMBER TO CELL-SPECIFIC POINTERS  *
C*                                                                     *
C***********************************************************************
C
      DO 400 izz=1,nz
         IJKR  = irn(izz)/neq
         IMODR = MOD(irn(izz),neq)
         IF(IMODR.EQ.0) THEN
            IJKMEL = IJKR
            IPH2   = neq
         ELSE
            IJKMEL = IJKR+1
            IPH2   = IMODR
         END IF
         i8 = ijkmel
C
         IJKC  = icn(izz)/neq
         IMODC = MOD(icn(izz),neq)
         IF(IMODC.EQ.0) THEN
            IJKCON = IJKC
            IPH1   = neq
         ELSE
            IJKCON = IJKC+1
            IPH1   = IMODC
         END IF
C
         IF(IJKMEL.EQ.IJKCON) THEN
            NB0(IPH1,IPH2,i8) = izz
         ELSE
c
            IF(i8.LT.IJKCON) THEN
C
               IF(jvect(i8).EQ.0) GO TO 400
               nabov = 0
               DO 340 index = jvect(i8),jvect(nela+i8)
                  id    = no(index)
                  nabov = nabov+1
                  IF(id.EQ.ijkcon) THEN
                     jjj = nabov
                     GO TO 342
                  END IF
  340          CONTINUE
  342          NC0(jjj+ISUMPP(i8-1),IPH1,IPH2) = izz
            ELSE
C
               IF(jvect(2*nela+i8).EQ.0) GO TO 400
               nbelo = 0
               DO 380 index = jvect(2*nela+i8),jvect(3*nela+i8)
                  id    = no(index)
                  nbelo = nbelo+1
                  IF(id.EQ.ijkcon) THEN
                     jjj = nbelo
                     GO TO 382
                  END IF
  380          CONTINUE
  382          NA0(jjj+ISUMMM(i8-1),IPH1,IPH2) = izz
            END IF
         END IF
C
  400 CONTINUE
C
      DO 500 ii = 1,leniw
         jvect(ii) = 0
  500 CONTINUE
C
 6001 FORMAT(//,20('ERROR-'),//,T33,
     &             '       S I M U L A T I O N   A B O R T E D',
     &       /,T24,'THE FIRST DIMENSION OF THE isumpp ARRAY ',
     &             'IS INSUFFICIENT FOR THIS PROBLEM',
     &       /,T12,'THE PARAMETER nconup MUST BE AT LEAST ',I10,
     &       //,T33,'               CORRECT AND TRY AGAIN',
     &       //,20('ERROR-'))
 6002 FORMAT(//,20('ERROR-'),//,T33,
     &             '       S I M U L A T I O N   A B O R T E D',
     &       /,T24,'THE FIRST DIMENSION OF THE isummm ARRAY ',
     &             'IS INSUFFICIENT FOR THIS PROBLEM',
     &       /,T12,'THE PARAMETER ncondn MUST BE AT LEAST ',I10,
     &       //,T33,'               CORRECT AND TRY AGAIN',
     &       //,20('ERROR-'))
C
C =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of REASSN
C
  999 RETURN
      END
C
      SUBROUTINE MTRXIN(IOPT)
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
      PARAMETER (SEED = 1.0d-25)
C
C***********************************************************************
C*                                                                     *
C*!               C O M M O N    D E C L A R A T I O N S               *
C*                                                                     *
C***********************************************************************
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/P4/R(MNEQ*MNEL+1)
C
      DIMENSION B0INV(5,5),TA0(5,5),TC0(5,5),TR0(5)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
C
      SAVE ICALL
      DATA ICALL/0/
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of MTRXIN
C
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(6X,'MTRXIN   1.1      10 September 2000',6X,
     &          'Routine for Z-preprocessing ',
     &          'of the Jacobian')
C
      MXMYMZ = NELA
      IPHASE = NEQ
C
C***********************************************************************
C*                                                                     *
C*            REPLACEMENT OF 0 BY SEED ON THE MAIN DIAGONAL            *
C*                                                                     *
C***********************************************************************
C
      IF(IOPT.GT.1) GO TO 61
      DO 20 I=1,NZ
         IF(IRN(I).EQ.ICN(I).AND.ABS(CO(I)).EQ.0.0d0) THEN
            CO(I) = SEED
         END IF
   20 CONTINUE
      GO TO 999
C
C***********************************************************************
C*                                                                     *
C*     ADDITION OF EQUATIONS TO PREVENT ZEROS ON THE MAIN DIAGONAL     *
C*                                                                     *
C***********************************************************************
C
   61 IF(IOPT.GT.2) GO TO 101
      DO 100 IJK=1,MXMYMZ
        DO 90 LA=1,IPHASE
C
          IF(CO(NB0(LA,LA,IJK)).EQ.0.0d0) THEN
             DO 80 LB=1,IPHASE
               IF(LB.EQ.LA) GO TO 80
C
               IF(CO(NB0(LA,LB,IJK)).NE.0.0d0) THEN
C
                  R(NX0(LA,IJK)) = R(NX0(LA,IJK))+1.75d0*R(NX0(LB,IJK))
                  DO 70 L1=1,IPHASE
                    CO(NB0(L1,LA,IJK))
     &              =        CO(NB0(L1,LA,IJK))
     &               +1.75d0*CO(NB0(L1,LB,IJK))
   70             CONTINUE
                  DO 75 L1=1,IPHASE
                    IF(IJKMM(ijk).LT.1) GO TO 73
                    DO 72 NDM=1,IJKMM(ijk)
                      NA0LOC = NDM+ISUMMM(ijk-1)
                      CO(NA0(NA0LOC,L1,LA))
     &                =        CO(NA0(NA0LOC,L1,LA))
     &                 +1.75d0*CO(NA0(NA0LOC,L1,LB))
   72               CONTINUE
   73               IF(IJKPP(ijk).LT.1) GO TO 75
                    DO 74 NDM=1,IJKPP(ijk)
                      NC0LOC = NDM+ISUMPP(ijk-1)
                      CO(NC0(NC0LOC,L1,LA))
     &                =        CO(NC0(NC0LOC,L1,LA))
     &                 +1.75d0*CO(NC0(NC0LOC,L1,LB))
   74               CONTINUE
   75             CONTINUE
                  GO TO 90
                END IF
   80         CONTINUE
C
            END IF
C
   90   CONTINUE
  100 CONTINUE
      GO TO 999
C
C***********************************************************************
C*                                                                     *
C*                NORMALIZATION + ADDITION OF EQUATIONS                *
C*                                                                     *
C***********************************************************************
C
  101 IF(IOPT.GT.4) GO TO 301
      DO 200 IJK=1,MXMYMZ
C
        ROWMAX = 0.0d0
        DO 150 LB=1,IPHASE
          DO 115 LA=1,IPHASE
            ABB0 = ABS(CO(NB0(LA,LB,IJK)))
            IF(ABB0.GT.ROWMAX) THEN
              ROWMAX = ABB0
              ROMX   = SIGN(ABB0,CO(NB0(LA,LB,IJK)))
            END IF
            IF(iopt.EQ.3) GO TO 115
c
            IF(IJKMM(ijk).LT.1) GO TO 113
            DO 112 NDM=1,IJKMM(ijk)
              NA0LOC = NDM+ISUMMM(ijk-1)
              ABA0   = ABS(CO(NA0(NA0LOC,LA,LB)))
              IF(ABA0.GT.ROWMAX) THEN
                ROWMAX = ABA0
                ROMX   = SIGN(ABA0,CO(NA0(NA0LOC,LA,LB)))
              END IF
  112       CONTINUE
c
  113       IF(IJKPP(ijk).LT.1) GO TO 115
            DO 114 NDM=1,IJKPP(ijk)
              NC0LOC = NDM+ISUMPP(ijk-1)
              ABC0   = ABS(CO(NC0(NC0LOC,LA,LB)))
              IF(ABC0.GT.ROWMAX) THEN
                ROWMAX = ABC0
                ROMX   = SIGN(ABC0,CO(NC0(NC0LOC,LA,LB)))
              END IF
  114       CONTINUE
c
  115     CONTINUE
C
          QUOT           = 1.0d0/ROMX
C
          R(NX0(LB,IJK)) = R(NX0(LB,IJK))*QUOT
          DO 120 LA=1,IPHASE
             CO(NB0(LA,LB,IJK)) = CO(NB0(LA,LB,IJK))*QUOT
            IF(IJKMM(ijk).LT.1) GO TO 117
            DO 116 NDM=1,IJKMM(ijk)
              NA0LOC                = NDM+ISUMMM(ijk-1)
              CO(NA0(NA0LOC,LA,LB)) = CO(NA0(NA0LOC,LA,LB))*QUOT
  116       CONTINUE
c
  117       IF(IJKPP(ijk).LT.1) GO TO 120
            DO 118 NDM=1,IJKPP(ijk)
              NC0LOC                = NDM+ISUMPP(ijk-1)
              CO(NC0(NC0LOC,LA,LB)) = CO(NC0(NC0LOC,LA,LB))*QUOT
  118       CONTINUE
c
  120     CONTINUE
  150   CONTINUE
  200 CONTINUE
C
C
C
  255 DO 300 IJK=1,MXMYMZ
        DO 290 LA=1,IPHASE
C
          IF(CO(NB0(LA,LA,IJK)).EQ.0.0d0) THEN
             DO 280 LB=1,IPHASE
               IF(LB.EQ.LA) GO TO 280
C
               IF(CO(NB0(LA,LB,IJK)).NE.0.0d0) THEN
C
                  R(NX0(LA,IJK)) = R(NX0(LA,IJK))+1.75d0*R(NX0(LB,IJK))
                  DO 270 L1=1,IPHASE
                    CO(NB0(L1,LA,IJK)) = CO(NB0(L1,LA,IJK))
     &                                  +1.75d0*CO(NB0(L1,LB,IJK))
  270             CONTINUE
                  DO 275 L1=1,IPHASE
c
                    IF(IJKMM(ijk).LT.1) GO TO 273
                    DO 272 NDM=1,IJKMM(ijk)
                      NA0LOC = NDM+ISUMMM(ijk-1)
                      CO(NA0(NA0LOC,L1,LA))
     &                =        CO(NA0(NA0LOC,L1,LA))
     &                 +1.75d0*CO(NA0(NA0LOC,L1,LB))
  272               CONTINUE
c
  273               IF(IJKPP(ijk).LT.1) GO TO 275
                    DO 274 NDM=1,IJKPP(ijk)
                      NC0LOC = NDM+ISUMPP(ijk-1)
                      CO(NC0(NC0LOC,L1,LA))
     &                =        CO(NC0(NC0LOC,L1,LA))
     &                 +1.75d0*CO(NC0(NC0LOC,L1,LB))
  274               CONTINUE
c
  275             CONTINUE
                  GO TO 290
                END IF
  280         CONTINUE
C
            END IF
C
  290   CONTINUE
  300 CONTINUE
      GO TO 999
C
C
C***********************************************************************
C*                                                                     *
C*    MULTIPLICATION BY THE INVERSE OF THE MAIN-DIAGONAL SUBMATRIX     *
C*                                                                     *
C***********************************************************************
C
C
  301 IF(NEQ.EQ.3) GO TO 401
C
C
C
      DO 400 IJK=1,MXMYMZ
C
        DETR = CO(NB0(1,1,IJK))*CO(NB0(2,2,IJK))-
     &         CO(NB0(1,2,IJK))*CO(NB0(2,1,IJK))
        IF(ABS(DETR).LE.SEED) THEN
           write(34,6001)
           STOP
        END IF
        DETI = 1.0d0/DETR
C
        B0INV(1,1) = CO(NB0(2,2,IJK))*DETI
        B0INV(1,2) =-CO(NB0(1,2,IJK))*DETI
        B0INV(2,1) =-CO(NB0(2,1,IJK))*DETI
        B0INV(2,2) = CO(NB0(1,1,IJK))*DETI
C
        CO(NB0(1,1,IJK)) = 1.0d0
        CO(NB0(2,2,IJK)) = 1.0d0
        CO(NB0(1,2,IJK)) = 0.0d0
        CO(NB0(2,1,IJK)) = 0.0d0
C
        DO 390 LB=1,IPHASE
          TR0(LB) = 0.0d0
          DO 355 L1=1,IPHASE
            TR0(LB) = TR0(LB)+B0INV(L1,LB)*R(NX0(L1,IJK))
  355     CONTINUE
          R(NX0(LB,IJK)) = TR0(LB)
C
          DO 380 LA=1,IPHASE
            TA0(LA,LB) = 0.0d0
            TC0(LA,LB) = 0.0d0
c
            IF(IJKMM(ijk).LT.1) GO TO 371
            DO 370 NDM=1,IJKMM(ijk)
              NA0LOC = NDM+ISUMMM(ijk-1)
              DO 360 L1=1,IPHASE
                 TA0(LA,LB) = TA0(LA,LB)
     &                       +CO(NA0(NA0LOC,L1,LB))*B0INV(LA,L1)
  360         CONTINUE
              CO(NA0(NA0LOC,LA,LB)) = TA0(LA,LB)
  370       CONTINUE
c
  371       IF(IJKPP(ijk).LT.1) GO TO 380
            DO 375 NDM=1,IJKPP(ijk)
              NC0LOC = NDM+ISUMPP(ijk-1)
              DO 372 L1=1,IPHASE
                 TC0(LA,LB) = TC0(LA,LB)
     &                       +CO(NC0(NC0LOC,L1,LB))*B0INV(LA,L1)
  372         CONTINUE
              CO(NC0(NC0LOC,LA,LB)) = TC0(LA,LB)
  375       CONTINUE
c
  380     CONTINUE
  390   CONTINUE
  400 CONTINUE
      GO TO 999
C
C
C
  401 DO 500 IJK=1,MXMYMZ
C
        DET1 = CO(NB0(2,2,IJK))*CO(NB0(3,3,IJK))
     &        -CO(NB0(2,3,IJK))*CO(NB0(3,2,IJK))
        DET2 = CO(NB0(2,1,IJK))*CO(NB0(3,3,IJK))
     &        -CO(NB0(2,3,IJK))*CO(NB0(3,1,IJK))
        DET3 = CO(NB0(2,1,IJK))*CO(NB0(3,2,IJK))
     &        -CO(NB0(2,2,IJK))*CO(NB0(3,1,IJK))
        DETR = CO(NB0(1,1,IJK))*DET1
     &        +CO(NB0(1,2,IJK))*DET2
     &        +CO(NB0(1,3,IJK))*DET3
C
        IF(ABS(DETR).LE.SEED) THEN
           write(34,6001)
           STOP
        END IF
        DETI = 1.0d0/DETR
C
        B0INV(1,1) = DET1*DETI
        B0INV(1,2) =(CO(NB0(1,3,IJK))*CO(NB0(3,2,IJK))
     &              -CO(NB0(1,2,IJK))*CO(NB0(3,3,IJK)))*DETI
        B0INV(1,3) =(CO(NB0(1,2,IJK))*CO(NB0(2,3,IJK))
     &              -CO(NB0(1,3,IJK))*CO(NB0(2,2,IJK)))*DETI
C
        B0INV(2,1) =-DET2*DETI
        B0INV(2,2) =(CO(NB0(1,1,IJK))*CO(NB0(3,3,IJK))
     &              -CO(NB0(1,3,IJK))*CO(NB0(3,1,IJK)))*DETI
        B0INV(2,3) =(CO(NB0(1,3,IJK))*CO(NB0(2,1,IJK))
     &              -CO(NB0(1,1,IJK))*CO(NB0(2,3,IJK)))*DETI
C
        B0INV(3,1) = DET3*DETI
        B0INV(3,2) =(CO(NB0(1,2,IJK))*CO(NB0(3,1,IJK))
     &              -CO(NB0(1,1,IJK))*CO(NB0(3,2,IJK)))*DETI
        B0INV(3,3) =(CO(NB0(1,1,IJK))*CO(NB0(2,2,IJK))
     &              -CO(NB0(1,2,IJK))*CO(NB0(2,1,IJK)))*DETI
C
        CO(NB0(1,1,IJK)) =-1.0d0
        CO(NB0(2,2,IJK)) =-1.0d0
        CO(NB0(3,3,IJK)) =-1.0d0
        CO(NB0(1,2,IJK)) = 0.0d0
        CO(NB0(1,3,IJK)) = 0.0d0
        CO(NB0(2,1,IJK)) = 0.0d0
        CO(NB0(2,3,IJK)) = 0.0d0
        CO(NB0(3,1,IJK)) = 0.0d0
        CO(NB0(3,2,IJK)) = 0.0d0
C
        DO 440 LB=1,IPHASE
          TR0(LB) = 0.0d0
          DO 405 L1=1,IPHASE
            TR0(LB) = TR0(LB)+B0INV(L1,LB)*R(NX0(L1,IJK))
  405     CONTINUE
          R(NX0(LB,IJK)) =-TR0(LB)
C
          DO 430 LA=1,IPHASE
            TA0(LA,LB) = 0.0d0
            TC0(LA,LB) = 0.0d0
c
            IF(IJKMM(ijk).LT.1) GO TO 421
            DO 420 NDM=1,IJKMM(ijk)
              NA0LOC = NDM+ISUMMM(ijk-1)
              DO 410 L1=1,IPHASE
                 TA0(LA,LB) = TA0(LA,LB)
     &                       +CO(NA0(NA0LOC,L1,LB))*B0INV(LA,L1)
  410         CONTINUE
              CO(NA0(NA0LOC,LA,LB)) =-TA0(LA,LB)
  420       CONTINUE
c
  421       IF(IJKPP(ijk).LT.1) GO TO 430
            DO 425 NDM=1,IJKPP(ijk)
              NC0LOC = NDM+ISUMPP(ijk-1)
              DO 422 L1=1,IPHASE
                 TC0(LA,LB) = TC0(LA,LB)
     &                       +CO(NC0(NC0LOC,L1,LB))*B0INV(LA,L1)
  422         CONTINUE
              CO(NC0(NC0LOC,LA,LB)) =-TC0(LA,LB)
  425       CONTINUE
c
  430     CONTINUE
  440   CONTINUE
  500 CONTINUE
C
 6001 FORMAT(//,20('ERROR-'),//,T40,
     &             '      S I M U L A T I O N   H A L T E D',
     &       /,T35,'THE DETERMINANT OF THE MAIN DIAGONAL SUBMATRIX ',
     &             'IS ZERO',/,
     &         T31,'Preprocessors ZPROCS = Z4 and OPROCS = O4 cannot ',
     &             'be used! ',/,
     &         T40,'        PLEASE CORRECT AND TRY AGAIN',/,
     &       //,20('ERROR-'))
C
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of MTRXIN
C
C
  999 RETURN
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
      SUBROUTINE ELINDX
C
C*********************************************************************
C*                                                                   *
C*!           DETERMINE THE POSITIONING ARRAYS ia AND ja             *     
C*!         WHICH IDENTIFY THE LOCATIONS OF THE NEIGHBORS            *     
C*!                 Version 1.00, January 16, 1998                   *     
C*                                                                   *
C*********************************************************************
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
C*********************************************************************
C*                                                                   *
C*!              C O M M O N    D E C L A R A T I O N S              *     
C*                                                                   *
C*********************************************************************
C
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
C
      COMMON/L7/JVECT(niwork)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
      COMMON/soll/lenw,leniw
C
      SAVE ICALL
      DATA ICALL/0/
C
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of ELINDX
C
C
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(/6X,'ELINDX   1.00     16 January   1998',6X,
     &           'Routine for neighbor element indexing')
C
C*********************************************************************
C*                                                                   *
C*!                 I N I T I A L I Z A T I O N S                    *
C*                                                                   *
C*********************************************************************
C
      DO 100 ii = 1,nel
         iijjkk(ii) = 1
  100 CONTINUE
C
      DO 110 ii = 1,mncon+mnel+1
         no(ii) = 0
  110 CONTINUE
C
      DO 120 ii = 1,leniw
         jvect(ii) = 0
  120 CONTINUE
C
C*********************************************************************
C*                                                                   *
C*!     DETERMINE NUMBER OF NEIGHBORING ELEMENTS AND STORE IN iw     *
C*                                                                   *
C*********************************************************************
C
      DO 200 n = 1,ncon
         n1 = nex1(n)
         n2 = nex2(n)
c
         IF(n1.EQ.0.OR.n2.EQ.0) GO TO 200
c
         IF(n1.LE.nela) iijjkk(n1) = iijjkk(n1)+1         
         IF(n2.LE.nela) iijjkk(n2) = iijjkk(n2)+1         
c
  200 CONTINUE
C
C*********************************************************************
C*                                                                   *
C*!   STORE TEMPORARILY THE # OF NEIGHBORS PER ELEMENT ii IN jvect   *
C*!                  (FIRST nela ELEMENTS OF jvect)                  *
C*                                                                   *
C*********************************************************************
C
      DO 250 ii = 1,nel
         jvect(ii) = iijjkk(ii)
  250 CONTINUE
C
C*********************************************************************
C*                                                                   *
C*!      DETERMINE THE CUMULATIVE NUMBER OF NEIGBORING ELEMENTS      *
C*!     STORE IN jvect - FROM jvect(mshift) TO jvect(mshift+nel)     *
C*                                                                   *
C*********************************************************************
C
      mshift     = nel+1
      no(mshift) = 0
C 
      isum = 0
      DO 300 ii = 1,nel
         isum             = isum+iijjkk(ii)
         jvect(mshift+ii) = isum
  300 CONTINUE
C
C*********************************************************************
C*                                                                   *
C*!   DETERMINE THE iw AND ja ARRAYS WITH THE NEIGHBOR INFORMATION   *
C*                                                                   *
C*********************************************************************
C
      DO 320 ii = 1,nel
         iijjkk(ii)               = 1 
         no(jvect(mshift+ii-1)+1) = ii
  320 CONTINUE
C
      DO 350 n = 1,ncon
         n1 = nex1(n)         
         n2 = nex2(n)         
c    
         IF(n1.GT.nela.AND.n2.GT.nela) GO TO 350
         IF(n1.EQ.0.OR.n2.EQ.0)        GO TO 350
c    
         iijjkk(n1) = iijjkk(n1)+1         
         inx1       = jvect(mshift+n1-1)+iijjkk(n1)
         no(inx1)   = n2
c    
         iijjkk(n2) = iijjkk(n2)+1
         inx2       = jvect(mshift+n2-1)+iijjkk(n2)
         no(inx2)   = n1
c    
  350 CONTINUE
C
C!    
C
      iijjkk(1) = 1
      no(1)     = 1
      DO 400 ii = 1,nel
c    
         itemp            = 0
         iijjkk(ii+1)     = iijjkk(1)+jvect(mshift+ii)
         no(iijjkk(ii+1)) = ii+1
c  
c ------------
c ...... Sort the subarrays ja(ii), ii=iijjkk(ii),...,iijjkk(ii+1)-1
c ------------
c  
         istart = iijjkk(ii)
         numsrt = iijjkk(ii+1)-iijjkk(ii)
c  
         CALL SORT(no(istart),numsrt)
c  
c ------------
c ...... Put the diagonal element first,  
c ...... swap with element in sorted subarray 
c ------------
c  
         DO 380 i5 = iijjkk(ii),iijjkk(ii+1)-1
            IF(no(i5).EQ.ii) THEN
               itemp = i5
               GO TO 390
            END IF
  380    CONTINUE
c  
  390    no(itemp)      = no(iijjkk(ii))
         no(iijjkk(ii)) = ii
c  
  400 CONTINUE
C
      DO 420 ii = 1,leniw
         jvect(ii) = 0       
  420 CONTINUE
C
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of ELINDX
C
C
 9999 RETURN
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
      SUBROUTINE SORT(iar,n)
C
C*********************************************************************
C*                                                                   *
C*!           SORTING THE VECTOR iar IN ASCENDING ORDER              *     
C*!                 Version 1.00, January 14, 1998                   *     
C*                                                                   *
C*********************************************************************
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
C
C*********************************************************************
C*                                                                   *
C*!                    L O C A L    A R R A Y S                      *
C*                                                                   *
C*********************************************************************
C
      DIMENSION iar(n)
C
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of SORT
C
C
      m = n
c 
      DO 1000 i5 = 1,n
         m = m/2
         IF(m.EQ.0) GO TO 9999
c 
         max = n-m
         DO 100 j = 1,max
c 
            DO 50 k = j,1,-m
               iert = iar(k+m)
               IF(iert.GE.iar(k)) go to 100
c 
               itemp    = iar(k+m)
               iar(k+m) = iar(k)
               iar(k)   = itemp
   50       CONTINUE
c 
  100    CONTINUE
c 
 1000 CONTINUE
C
C
C! =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   End of SORT
C
C
 9999 RETURN
       END
C
C
C
C***********************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C***********************************************************************
C
