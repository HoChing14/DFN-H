C
      SUBROUTINE INPUT
C
C-----READ ALL DATA PROVIDED THROUGH THE INPUT FILE.
C
C
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
      COMMON/C1/NEX1(MNCON)               ! 连接内第一个网格的编号
      COMMON/C2/NEX2(MNCON)               ! 连接内第二个网格的编号
      COMMON/C3/DEL1(MNCON)               ! 第一个网格中心点到连接面的垂直距离
      COMMON/C4/DEL2(MNCON)               ! 第二个网格中心点到连接面的垂直距离
      COMMON/C5/AREA(MNCON)               ! 连接面的面积
      COMMON/C6/BETA(MNCON)               ! 两个网格中心点连线与垂直方向夹角的余弦值
      COMMON/C7/ISOX(MNCON)               ! 表征连接面各向同性的指标
      COMMON/C8/GLO(MNCON)                ! 热流速率
      COMMON/C9/ELEM1(MNCON)              ! 连接内第一个网格的名称
      COMMON/C10/ELEM2(MNCON)             ! 连接内第二个我网格的名称
      COMMON/C11/FVD(MNCON)               ! FVD是每个网格的气相流量
      common/c12/sig(MNCON)               ! sig(MNCON)是保存辐射热传递的数组
      common/c13/ALPHA(3),IPMAT
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/nnn/nkin,nkin1
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
C
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G1/F1(MGTAB)                 ! 源汇项时间的列表 
      COMMON/G2/F2(MGTAB)                 ! 源汇项流速的列表
      COMMON/G3/F3(MGTAB)                 ! 注入焓值的列表
      COMMON/G9/NEXG(MNOGN)               ! 这TM又是啥 ?
c
c.... Additional local variables for permeability (els10/23/00)
      double precision perx,pery,perz     ! 三个方向的渗透率
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SVZ/NOITE,MOP(24)
!
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +              CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +            GK(MAXMAT)
      COMMON/SOLII/XKD3(MAXMAT),XKD4(MAXMAT),SII3(MAXMAT)
!
c.....Add tortuosity exponent (ptort) and critical porosity (phicrit)
      common/torpar/ptort(maxmat),phicrit(maxmat)
c
c     for T2VOC
      COMMON/ADSORP/OCK,FOX,FOC(MAXMAT),ALAM
      COMMON/SOIN/YIN(10,MAXMAT)
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
      COMMON/STEP/ELST
      
c.....Added commons for writing out element specific data (TOUGH2 v1.6)
      COMMON/STEPrk1/eplist(200)
      COMMON/STEPrk2/nstrick(200),nelist
      character*5  eplist
c
      COMMON/STE1/NST
      COMMON/DEFINI/DEP(10)
      COMMON/DG/WUP,WNR
      COMMON/DFM/TIMAX,REDLT
      COMMON/DLT/NDLT,DLT(100)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
C      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO,TIMP1
      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO
      COMMON/DOP/ENTH,KDATA,QUAL
      COMMON/DX/K,NI,SCALE
      COMMON/POV6/TSTART
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/PATCH/SING
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
C-----NOW FOR A COMMON BLOCK WITH EOS-SELECTION PARAMETERS.
      COMMON/EOSEL/IE(16),FE(512)
      COMMON/RPCAP/IRP(MAXMAT),RP(7,MAXMAT),ICP(MAXMAT),CP(7,MAXMAT),
     A      IRPD,RPD(7),ICPD,CPD(7)
      COMMON/BIND/DIFF0,TEXP,BE
      COMMON/V/IS
c
c.....Added for correct variable size
c     real*4 tzero
      real*4 rtzero
C------------------------------------------------- for using EOS9
      common/ref_tc/ Tcref     ! reference temperature
c----------------------------------------------------------------    
c.....for T2VOC
      COMMON/CRITP/TCRIT,PCRIT,ZCRIT,OMEGA,DIPOLM
      COMMON/VPOIL/TBOIL,VPA,VPB,VPC,VPD
      COMMON/HCAPL/AMO,CPA,CPB,CPC,CPDD
      COMMON/DENOIL/RHOREF,TDENRF,DIFV0,TDIFRF,TEXPO
      COMMON/VISOIL/VLOA,VLOB,VLOC,VLOD,VOLCRT
      COMMON/HCOIL/SOLA,SOLB,SOLC,SOLD
c....
C
      COMMON/SOLVR1/matslv,nmaxit,nnvvcc,iiuunn,iissoo,nactdi
      COMMON/SOLVR2/ritmax,closur
      COMMON/SOLVR3/ordrng,oprocs,zprocs,coord
C
      common/fgt1/ioft,iofu,igoft,igofu,noft(100),ngoft(100)
      common/fgt2/eoft(100),egoft(100),ecoft(100)
      common/fgt3/icoft,icofu,ncoft(100)
      common/ff/h1
      common/dipa/fddiag(3,5)
      common/dipa1/iddiag
c
      character*1 h1
      character*5 eoft,egoft,desc,fdf,jat
      character*10 ecoft
C
      CHARACTER*5 ELEM1,ELEM2,MAT,ELST,VER,WORD,MA12,NAM
      CHARACTER*3 EL,MA1,EL1,EL2,SL
      CHARACTER MA2*2,TYPE*4,ITAB*1
      CHARACTER*2 ordrng,oprocs,zprocs
      CHARACTER*5 coord
      character w75*75,comm*80
C
c.... Commons for effective thermal conductivity
      common/efkth/timkth(mgtab),fackth(mgtab)
      common/kthtable/ktftb(mnogn)
      
C----------------------------------------- for coupling with reactive transport
C
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
C
      DIMENSION  VER(22),XIN(10,maxmat),NAM(maxmat),comm(50),xx(10)
C--------------------------------------------------------------------------
      SAVE ICALL,icom
      DATA ICALL,icom/0,0/
      DATA VER /'PARAM','GENER','INCON','START','RESTA','ENDCY','ENDFI'
     A,'ROCKS','MULTI','TIMES','SELEC','RPCAP','MESHM','INDOM','FOFT ',
     X'GOFT ','DIFFU','COFT ','CHEMP','SOLVR','REACT','REFPT'/
      WRITE(34,892)
  892 FORMAT(/1X,60('*'),'call subroutine INPUT in t2f_v2.f',60('*'))
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' INPUT 1.50, 2000.5.25: READ DATA PROVIDED THROUGH flow.'
     x'inp'/)

!.....for coupling with reactive transport
     
      MOPR(1)=2            ! default value for only flow
      MOPR(2)=0
      MOPR(3)=0
      MOPR(4)=1
      DO I=5,20
         MOPR(I)=0
      END DO
!
!-------------------------------------------------------------------------
c.....initialize parameters for time slices
      ioft=0
      igoft=0
      icoft=0
c
      NI=0  ! SET NI=1 FOR START (ALLOW *INCON* EXISTS ARBITRARILY)
C
      DELTMX=0.d0
c------------
      ITI=0
C-----FLAG FOR THE *TIMES* DATA BLOCK. IF THIS BLOCK IS PRESENT,
C     ITI WILL BE SET EQUAL TO THE NUMBER OF TIME DATA PROVIDED.
c------------      
      KIN=0
C-----FLAG FOR THE *INDOM* DATA BLOCK. SET EQUAL TO THE NUMBER OF
c     DOMAINS DEFINED IN IT.

c------------      
      IS=0  ! IS=1 WHEN KEYWORD "ENDFI" IS ENCOUNTERED IN THE INPUT FILE
C     AND FLOW SIMULATION WILL BE BYPASSED. OFTEN USED IN MESHMAKER.
C
C     ASSIGN DEFAULTS FOR RELATIVE PERMEABILITY AND CAPILLARY PRESSURE.
C     DEFAULTS ARE KREL = 1, PCAP = 0, APPROPRIATE FOR SINGLE PHASE FLOW
      
      IRPD=5
      ICPD=1
      CPD(1)=0.d0
      CPD(2)=0.d0
      CPD(3)=1.d0
C
      iissoo = 0
c
c.....initialize arrays for primary variables
      do 1910 i=1,10
      dep(i)=0.d0
 1910 xx(i)=0.d0
C
c.....initialize component and equation counters
      NEQ1=NEQ+1                  ! number of equations in iterations
      NK1=NK+1                    ! NK components and energy balance equations
      NFLUX=2*NEQ+1
      NBK=NB+NK                   ! number of first set secondary parameters
      NSEC=NPH*NBK+2              ! total number of secondary parameters
c.....set component counter for intializing thermodynamic conditions; 
c     default values correspond to what is desired for the flow simulation
      nkin=nk
      nkin1=nk1
c.....flag for 'DIFFU'-data (diffusion coefficients)
      iddiag=0
 5019 READ (33,5020) WORD,w75  ! read keywords
 5020 FORMAT(A5,A75)
C
C      DO900 K=1,23
      DO 900 K=1,24
C--------------------------------------
  900 IF(WORD.EQ.VER(K)) GOTO 920
c.....come here for unknown keyword
      icom=icom+1
      if(icom.le.50) comm(icom)=word//w75
c     WRITE (34,901) WORD
c 901 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL "',A5,
c    X'" --- IGNORE THIS, AND CONTINUE READING INPUT DATA')
      GOTO 5019
C
  920 GOTO(1000,1300,1900,1450,5019,1500,1600,1700,1800,2000,
     A2100,2200,2300,2400,80,180,190,280,2600,3000,4000,1452),K
C-------------------------------------------------
C
C*****READ ROCK PROPERTIES.*********************************************
C                                                                      *
c     keyword 'ROCKS'
 1700 NM=1
    2 READ (33,1) MAT(NM),NAD,DM(NM),POR(NM),(PER(I,NM),I=1,3),
     ACWET(NM),SH(NM)  ! record ROCKS-1
    1 FORMAT(A5,I5,7E10.4)  

C------------------------------------------------ For using EOS9
      if (mat(nm).eq.'REFCO') then
         if(por(nm).gt.0.d0) Tcref=por(nm)
      end if
c---------------------------------------------------------------     
      IF(MAT(NM).EQ.'     ') GOTO 3 ! blank record for close this block
C
      YIN(1,NM)=0.d0              
      COM(NM)=0.d0               
      EXPAN(NM)=0.d0              
      CDRY(NM)=0.d0               
      TORT(NM)=0.d0              
      GK(NM)=0.d0                
      IRP(NM)=0                  
      ICP(NM)=0                  
      
cels6/8/09 Add tortuosity exponent (ptort) and critical porosity (phicrit)
      ptort(nm) = 0.d0
      phicrit(nm) = 0.d0
      IF(NAD.GE.1) READ (33,7) COM(NM),EXPAN(NM),CDRY(NM),TORT(NM),
     x    GK(NM),XKD3(NM),XKD4(NM),ptort(nm),phicrit(nm) ! record rocks-2
c.....5-25-00: assign fraction organic carbon for T2VOC
      foc(nm)=xkd3(nm)    
    7 FORMAT(9E10.4)
      IF(CDRY(NM).EQ.0.d0) CDRY(NM)=CWET(NM)
      IF(NAD.LE.1) GOTO 1701
C-----READ PARAMETERS FOR RELATIVE PERMEABILITY AND CAPILLARY PRESSURE (rocks-3)
      READ (33,1702) IRP(NM),(RP(I,NM),I=1,7)
      READ (33,1702) ICP(NM),(CP(I,NM),I=1,7)
 1702 FORMAT(I5,5X,7E10.4)
 1701 CONTINUE
C
      NAD=0
C
      WRITE (34,5) NM,MAT(NM)
    5 FORMAT(' DOMAIN NO.',I3,': ',A5)
      NM=NM+1
      GOTO 2
    3 CONTINUE
      NM=NM-1
      GOTO 5019
C
C-----END OF ROCK PROPERTIES.-------------------------------------------
C
c     keyword 'RPCAP'
 2200 CONTINUE
C-----COME HERE TO READ DEFAULT ASSIGNMENTS FOR RELATIVE PERMEABILITY
C     AND CAPILLARY PRESSURE PARAMETERS.
      READ (33,1702) IRPD,(RPD(I),I=1,7)
      READ (33,1702) ICPD,(CPD(I),I=1,7)
      GOTO 5019
C
c     keyword 'MESHM': INTERNAL MESH GENERATION MODULE.
 2300 CONTINUE 
      CALL MESHM
      GOTO 5019
C
c     keyword 'REFPT': GET REFERENCE POINTS' WEIGHT INFORMATION
 1452 CONTINUE
      READ(33,517) IPMAT,(ALPHA(I),I=1,3)  ! MATRIX THERMAL PROPERTIES AND REFERENCE TEMPERATURE
  517 FORMAT(I2,3X,3E10.4)
      GOTO 5019

c     keyword 'START'
 1450 NI=1
      GOTO 5019
C*****READ COMPUTATION PARAMETERS.**************************************
c     keyword 'PARAM'
 1000 READ (33,200) NOITE,KDATA,MCYC,MSEC,MCYPR,(MOP(I),I=1,24),
     XDIFF0,TEXP,BE
  200 FORMAT(2I2,3I4,24I1,4E10.4)
C
      READ (33,8) TSTART,TIMAX,DELTEN,DELTMX,ELST,GF,REDLT,SCALE
    8 FORMAT(4E10.4,A5,5X,3E10.4)
c
c           from tough2 v1.6
      nelist=0
      if(elst.eq.'wdata') then
          read(33,*)   nelist
          read(33,2455)  (eplist(i),i=1,nelist)
        OPEN(UNIT=67,FILE='GASOBS.DAT',FORM='FORMATTED'
     1      ,STATUS='UNKNOWN')
      endif
 2455 format(a5)
c
      NDLT=0
      IF(DELTEN.GE.0.d0) GOTO 4               
      NDLT=-INT(DELTEN)                       
      DO 6 N=1,NDLT
    6 READ (33,7) (DLT(J+8*(N-1)),J=1,8)
      DELTEN=DLT(1)
    4 CONTINUE
C
      READ (33,7) RE1,RE2,U,WUP,WNR,DFAC,FOR
C
C-----ASSIGN DEFAULT VALUES.
C
      IF(NOITE.EQ.0) NOITE=8
      IF(KDATA.EQ.0) KDATA=1
      IF(MCYPR.EQ.0) MCYPR=1
c.....3-11-98: next statement commented out (re/EOS7R)
      IF(REDLT.EQ.0.d0) REDLT=4.d0
      IF(SCALE.EQ.0.d0) SCALE=1.d0
      IF(RE1.EQ.0.d0) RE1=1.d-5
      IF(RE2.EQ.0.d0) RE2=1.d0
      IF(FOR.EQ.0.d0) FOR=1.d0
      IF(WUP.LE.0.d0) WUP=1.d0
      IF(WNR.EQ.0.d0) WNR=1.d0
C     IF(DFAC.EQ.0.) DFAC=1.E-8
      IF(U.EQ.0.d0) U=0.1d0
C
C-----INITIALIZE SOME PARAMETERS.
C
      KCYC=0
      ITER=0
      ITERC=0
      TIMIN=TSTART
C
      READ (33,1908) (DEP(i),i=1,nkin1)     
C
C-----END OF COMPUTATION PARAMETERS.------------------------------------
C
      GOTO 5019
C
c     keyword 'GENER'
 1300 CONTINUE
      WRITE (34,*) '    '
      REWIND 3
      write (34,*) 'WRITE FILE *GENER* FROM INPUT DATA'
      WRITE(3,"(5HGENER)")
      NOGN=0
c     Modifications for time-dependent factor for thermal conductivity

   35 READ(33,30)EL,NE,SL,NS,NSEQ,NADD,NADS,LTAB,TYPE,ITAB,GX,
     + EX,HX,ktab
c
   30 FORMAT(A3,I2,A3,I2,4I5,5X,A4,A1,3E10.4,i2)
      IF(EL.EQ.'   '.AND.NE.EQ.0) GOTO 60
      IF(EL.EQ.'+++') GOTO 61
      NSEQ1=NSEQ+1
      IF(LTAB.EQ.0) LTAB=1
C
      DO 33 I=1,NSEQ1
        NOGN=NOGN+1
        N1=NE+(I-1)*NADD
        N2=NS+(I-1)*NADS
        WRITE(3,34) EL,N1,SL,N2,LTAB,TYPE,ITAB,GX,EX,HX,ktab
   34   FORMAT(A3,I2,A3,I2,15X,I5,5X,A4,A1,3G10.4,i2)
        LTABA=ABS(LTAB)
c Times and factors for effective thermal conductivity
        if(ktab.gt.0)then
           READ (33,36)(timkth(L),L=1,ktab)
           READ (33,36)(fackth(L),L=1,ktab)
           write(3,36)(timkth(L),L=1,ktab)
           write(3,36)(fackth(L),L=1,ktab)
        endif
        IF (LTABA.LE.1.OR.TYPE.EQ.'DELV') GO TO 33
        IF (I.GE.2) GO TO 133
C-----COME HERE TO READ TABLES OF TIMES AND RATES.
        READ (33,36) (F1(L),L=1,LTABA)
        READ (33,36) (F2(L),L=1,LTABA)
        IF (ITAB.NE.' ') READ (33,36) (F3(L),L=1,LTABA)
   36   FORMAT(4E14.7)
  133   CONTINUE
        WRITE(3,36) (F1(L),L=1,LTABA)
        WRITE(3,36) (F2(L),L=1,LTABA)
        IF (ITAB.NE.' ') WRITE(3,36) (F3(L),L=1,LTABA)
   33   CONTINUE
      GOTO 35
C
C-----END OF SINK/SOURCE DATA.-----------------------------------------
C
   60 WRITE(3,*) '     '
      GOTO 5019
   61 WRITE(3,*) '+++  '
      READ (33,1505) (NEXG(N),N=1,NOGN)
 1505 FORMAT(16I5)
      WRITE(3,1505) (NEXG(N),N=1,NOGN)
      GOTO5019
C
C*****READ MULTICOMPONENT PARAMETERS************************************
C
c     keyword 'MULTI'
 1800 CONTINUE
      READ (33,1801) NK,NEQ,NPH,NB,nkin
      NEQ1=NEQ+1
      NK1=NK+1
      NFLUX=2*NEQ+1
      NBK=NB+NK
      NSEC=NPH*NBK+2
      if(nkin.eq.0) nkin=nk
      nkin1=nkin+1
 1801 FORMAT(16I5)
C
      GOTO 5019
C
C-----END OF MULTICOMPONENT PARAMETERS----------------------------------
C
C***** READ MULTI-COMPONENT INITIAL CONDITIONS *************************
C
c     keyword 'INCON'
 1900 CONTINUE
      REWIND 1
      WRITE (34,*) 'WRITE FILE *INCON* FROM INPUT DATA'
      WRITE(1,1901)
 1901 FORMAT(5HINCON)
C
C------------------------------------------------------------------------
c... Read in permeability components (perx,pery,perz)
 1906 READ (33,1902) EL,NE,NSEQ,NADD,PORX,perx,pery,perz
 1902 FORMAT(A3,I2,2I5,4E15.8,1G15.4)
C------------------------------------------------------------------------
C
      IF(EL.EQ.'   '.AND.NE.EQ.0) GOTO 1903
      IF(EL.EQ.'+++') GOTO 1904
C
      READ (33,1908) (xx(i),i=1,nkin1)
 1908 FORMAT(4E20.13)
C
      NSEQ1=NSEQ+1
      DO 1905 I=1,NSEQ1
      N1=NE+(I-1)*NADD
C
C-------------------------------------------------------------------------
c... Write permeability components (perx,pery,perz)
      WRITE(1,1907) EL,N1,PORX,perx,pery,perz
 1907 FORMAT(A3,I2,10X,4E15.8,1E15.4)
C-------------------------------------------------------------------------
C
      WRITE(1,1908) (xx(j),j=1,nkin1)
C
 1905 CONTINUE
C
      GOTO 1906
 1904 WRITE(1,*) '+++  '
C
      READ (33,74) KCYCX,ITERCX,NMX,TSTX,TIMINX
   74 FORMAT(2I10,I5,2E15.8)
      WRITE(1,74) KCYCX,ITERCX,NMX,TSTX,TIMINX
C
      GOTO 5019
C
 1903 WRITE(1,1403)
 1403 FORMAT(5H     )
      GOTO 5019
C
C-----END OF MULTICOMPONENT INITIAL CONDITIONS--------------------------
C
C***** READ DOMAIN-SPECIFIC INITIAL CONDITIONS ***********************
c     keyword 'INDOM'
 2400 READ (33,1) NAM(KIN+1)
      IF(NAM(KIN+1).EQ.'     ') GOTO 5019
      KIN=KIN+1
      READ (33,1908) (XIN(I,KIN),I=1,nkin1)
      GOTO 2400
C
C-----END OF DOMAIN-SPECIFIC INITIAL CONDITIONS------------------------
C
c     keyword 'NOVER'
 2500 CONTINUE
C-----SET FLAG FOR SUPPRESSING VERSION-PRINTOUT AT END OF RUN.
      IV=0
      GOTO 5019
C
C***********************************************************************
C*                                                                     *
C*            READ SOLVER TYPE AND CORRESPONDING PARAMETERS            *
C*                                                                     *
C***********************************************************************
c     keyword 'SOLVR'
 3000 iissoo = 1
      READ (33,3050) matslv,zprocs,oprocs,ritmax,closur
C
 3050 FORMAT(i1,2x,a2,3x,a2,2(e10.4))
      GOTO 5019
C
C-----END OF SOLVER DATA.-------------------------------------------
C
C***************Addition by Tianfu Xu for coupling with reactive transport
c     keyword 'REACT'
c
 4000 CONTINUE
      READ (33,4050) (MOPR(IR),IR=1,20)
 4050 FORMAT(20i1)
C
      GOTO 5019
C
C-----END OF REACT DATA.-------------------------------------------
C************************************************************************
C
c     化学反应相关的参数
c     keyword 'CHEMP'
 2600 CONTINUE
C-----ASSIGN VOC CHEMICAL PARAMETERS ('CHEMP')
      READ (33,2002) TCRIT,PCRIT,ZCRIT,OMEGA,DIPOLM
      READ (33,2002) TBOIL,VPA,VPB,VPC,VPD
      READ (33,2002) AMO,CPA,CPB,CPC,CPDD
      READ (33,2002) RHOREF,TDENRF,DIFV0,TDIFRF,TEXPO
      READ (33,2002) VLOA,VLOB,VLOC,VLOD,VOLCRT
      READ (33,2002) SOLA,SOLB,SOLC,SOLD
      READ (33,2002) OCK,FOX,ALAM
      GOTO 5019
C
C*****READ DATA BLOCK WITH TIMES.
      
c     keyword 'TIMES'
 2000 CONTINUE
      READ (33,2001) ITI,ITE,DELAF,TINTER
C
C     ITI IS THE NUMBER OF TIMES PROVIDED.
C     ITE IS THE TOTAL NUMBER OF TIMES DESIRED.
C     DELAF IS THE MAXIMUM TIME STEP TO BE TAKEN AFTER ANY OF THE
C           PRESCRIBED TIMES HAVE BEEN REACHED.
C     TINTER IS THE TIME INCREMENT TO BE APPLIED FOR THE TIME VALUES
C            WITH INDEX ITI+1,ITI+2, ... ,ITE
C
 2001 FORMAT(2I5,3E10.4)
      READ (33,2002) (TIS(I),I=1,ITI)
 2002 FORMAT(8E10.4)
      IF(ITE.LE.ITI.OR.TINTER.LE.0.d0) GOTO 5019
C
C-----COME HERE TO ASSIGN INCREMENTED TIMES.
      ITI1=ITI+1
      DO2003 I=ITI1,ITE
 2003 TIS(I)=TIS(I-1)+TINTER
C
      ITI=ITE
C
      GOTO 5019
C
C-----END OF TIMES BLOCK------------------------------------------------
C
C*****READ DATA BLOCK WITH SELECTION PARAMETERS*************************
c keyword 'SELEC'
C used with certain EOS-modules to supply thermophyscial property data
 2100 CONTINUE
      READ (33,1801) (IE(I),I=1,16)
      IE1=1
      IF(IE(1).GT.1) IE1=IE(1)
      IF(IE1.GT.64) WRITE (34,2101) IE1
 2101 FORMAT(' WARNING: IE(1) =',I5,' IN BLOCK *SELEC* > 64, CANNOT '
     X'READ ALL RECORDS'/)
      IE18=MIN(IE1*8,512)
      READ (33,2002) (FE(I),I=1,IE18)
      GOTO 5019
C
C-----END OF SELECTION PARAMETERS---------------------------------------
C
c*****come here for data block *FOFT* to read element names for which
c     time slices are desired.
c
c     keyword 'FOFT '
   80 continue
      do81 i=1,100                    ! FOFT的上限值在这里 !!!
      read (33,5020) eoft(i)
      if(eoft(i).eq.'     ') goto82
   81 continue
   82 ioft=i-1                        ! ioft是FOFT文件中的检测单元点个数
      iofu=ioft                       ! iofu也是FOFT文件中的检测点个数
      if(iofu.gt.0) open(12,file='FOFT',status='unknown')
      rewind 12
      goto5019
c
c*****come here to open file for time-dependent generation data.
c     keyword 'GOFT '
  180 continue
      open(13,file='GOFT',status='unknown')
      rewind 13
      do 181 i=1,100
      read (33,5020) egoft(i)
      if(egoft(i).eq.'     ') goto 182
  181 continue
  182 igoft=i-1                       ! igoft存储了goft检测点的个数
c.....when no element data are given in block GOFT, set igoft
c     to a flag that will cause all elements with generation
c     data to be tabulated
      if(i.eq.1) igoft=-1
      igofu=igoft
c     WRITE (34,183) igoft,(egoft(i),i=1,igoft)
  183 format(' INPUT !!!  igoft =',I4,' EGOFT ='/(20(1X,A5)))
      goto5019
c
c*****come here for data block *COFT* to read connections for which
c     time slices are desired.
c
c     keyword 'COFT '
  280 continue
      do 281 i=1,100
      read (33,5021) ecoft(i)
 5021 format(A10)
      if(ecoft(i).eq.'          ') goto 282
  281 continue
  282 icoft=i-1
      icofu=icoft             ! icofu和icoft都存储了COFT的检测点个数
      if(icofu.gt.0) open(14,file='COFT',status='unknown')
      rewind 14
      goto5019
c
c     keyword 'DIFFU'
  190 continue
c.....come here to read diffusion coefficients
      iddiag=iddiag+1
      do191 k=1,nk
  191 read (33,7) (fddiag(i,k),i=1,nph)
      goto 5019
c
c     keyword 'ENDFI'
 1600 WRITE (34,1601)
 1601 FORMAT(' HAVE READ KEYWORD *ENDFI* IN FILE INPUT --- BYPASS'
     X' FLOW SIMULATION'/' OPTIONAL PRINTOUT OF INPUT DATA (FOR'
     X' MOP(7).NE.0) IS AVAILABLE')
      IS=IS+1
c
c     keyword 'ENDCY'
 1500 CONTINUE
C
      IF(KIN.EQ.0) GOTO 2410
C***** COME HERE TO PROCESS DOMAIN-SPECIFIC INITIAL CONDITIONS *******
      DO 2401 J=1,KIN
      DO 2402 I=1,NM
      IF(MAT(I).EQ.NAM(J)) GOTO 2403
 2402 CONTINUE
      GOTO 2401
 2403 DO 2404 L=1,nkin1
 2404 YIN(L,I)=XIN(L,J)
 2401 CONTINUE
C
 2410 CONTINUE
C-----END OF DOMAIN-SPECIFIC INITIAL CONDITIONS ----------------------
      RETURN
      END
c
      SUBROUTINE RFILE
C
C-----READ INPUT DATA FROM DISK FILES, WHICH ARE EITHER PROVIDED
C     BY THE USER, OR ARE INTERNALLY GENERATED IN SUBROUTINE
C     INPUT FROM DATA GIVEN IN THE JOB DECK.
C
c.....2-10-98: modified to read multiple data sets with flowing
c              wellbore pressure corrections.
C
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
C
      COMMON/C1/NEX1(MNCON)                   
      COMMON/C2/NEX2(MNCON)                   
      COMMON/C3/DEL1(MNCON)                   
      COMMON/C4/DEL2(MNCON)                   
      COMMON/C5/AREA(MNCON)                   
      COMMON/C6/BETA(MNCON)                   
      COMMON/C7/ISOX(MNCON)                   
      COMMON/C8/GLO(MNCON)                    
      COMMON/C9/ELEM1(MNCON)                  
      COMMON/C10/ELEM2(MNCON)                 
      COMMON/C11/FVD(MNCON)                   
      common/c12/sig(MNCON)
      common/c13/ALPHA(3),IPMAT                  
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G1/F1(MGTAB)                     ! 抽出点的时间列表
      COMMON/G2/F2(MGTAB)                     ! 抽出点的流速列表
      COMMON/G3/F3(MGTAB)                     ! 注入点的焓值列表
      common/g3a/pw(mgtab)                    ! ？？？
      COMMON/G4/ELEG(MNOGN)                   !  
      COMMON/G5/SOURCE(MNOGN)                 !
      COMMON/G6/LTABG(MNOGN)                  !
      COMMON/G7/G(MNOGN)                      !
      COMMON/G8/EG(MNOGN)                     !
      COMMON/G9/NEXG(MNOGN)                   !
      COMMON/G10/ITABG(MNOGN)                 !
      COMMON/G11/NGIND(MNOGN)                 !
      COMMON/G12/LCOM(MNOGN)                  !
      COMMON/G13/PI(MNOGN)                    !
      COMMON/G14/PWB(MNOGN)                   !
      COMMON/G15/HG(MNOGN)                    !
      COMMON/G16/GPO(MNOGN)                   !
      COMMON/G17/SDENS(MNOGN)                 !
      COMMON/G18/SSAT(MNOGN)                  !
      COMMON/G19/GVOL(MNOGN)                  !
      COMMON/G20/HL(MNOGN)                    !
      COMMON/G21/HS(MNOGN)                    !
      COMMON/G22/QVGC(MNOGN)                  !
      COMMON/G23/QVWC(MNOGN)                  !
      COMMON/G24/QVOC(MNOGN)                  !
      COMMON/G25/GRAD(MNOGN)                  !
      COMMON/G26/FF(MNPH*MNOGN)               !
      common/g27/fnam(mnogn)                  !
      common/g28/nftab(mnogn)                 !
      common/g29/iftit(mnogn)                 !
      common/g30/jftit(mnogn)                 !
      common/g31/ijf(mnogn)                   !
C
c... Commons for effective thermal conductivity
      common/efkth/timkth(mgtab),fackth(mgtab)
      common/kthtable/ktftb(mnogn)
c
c... Additional local variables for permeability
      double precision perx,pery,perz
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
      COMMON/SVZ/NOITE,MOP(24)                
      COMMON/P1/X((MNK+1)*MNEL)  ! converged primary variables
      COMMON/E1/ELEM(MNEL) ! name of elements
      COMMON/E2/MATX(MNEL) ! name or index of rock domain
      COMMON/E3/EVOL(MNEL) ! 网格体积
      COMMON/E4/PHI(MNEL) ! porosity
      COMMON/E5/P(MNEL) ! pressure converged
      COMMON/E6/T(MNEL) ! temperature converged
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL) ! heat exchange area
c.....7-3-95: include array to be used for permeability modifiers
      common/e7/pm(mnel)
c.....7-20-93: define coordinate arrays
      common/xyz1/x1(mnel)
      common/xyz2/x2(mnel)
      common/xyz3/x3(mnel)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)         ! ???
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +            GK(MAXMAT)          ! ???
      COMMON/SOIN/YIN(10,MAXMAT)      ! ???
      COMMON/STEP/ELST                ! ???
c.....Added commons for writing out element specific data (TOUGH2 v1.6)
      COMMON/STEPrk1/eplist(200)
      COMMON/STEPrk2/nstrick(200),nelist
      character*5  eplist
c
      COMMON/STE1/NST
      COMMON/DEFINI/DEP(10)
      COMMON/DFM/TIMAX,REDLT          
      COMMON/DLT/NDLT,DLT(100)        
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO            
      COMMON/DOP/ENTH,KDATA,QUAL      
c added for correct variable size
      real*4 rtzero
c
      COMMON/DX/K,NI,SCALE
      COMMON/POV6/TSTART
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/nnn/nkin,nkin1
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/BC/NELA
      COMMON/MMC/MM,MC
      COMMON/V/IS
      COMMON/SOLVR1/matslv,nmaxit,nnvvcc,iiuunn,iissoo,nactdi
      common/fgt1/ioft,iofu,igoft,igofu,noft(100),ngoft(100)
      common/fgt2/eoft(100),egoft(100),ecoft(100)
      common/fgt3/icoft,icofu,ncoft(100)
      common/ran2/iran
c.....Missing declaration for type(n)
c      character*5 type(mnogn)
      common/source_type/type5(mnogn)
      character*5 type5

      character*5 eoft,egoft,ffi,typen,MA12,fnam
      character*10 ecoft
c
C----------------------------------------------------- For using EOS9
      common/ref_tc/ Tcref           ! reference temperature
      common/TEM_EOS9/Tc_EOS9(MNEL)  ! initial temperature (oC)
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
c--------------------------------------------------------------------
c
      CHARACTER*5 ELEM1,ELEM2,ELEG,SOURCE,ELEM,MAT,ELST,DENT,EL1,EL2,EL
      CHARACTER ITABG*1,typ*4,MA1*3,MA2*2,ty1*1,titl*80

      double precision DEPU(10)
      SAVE ICALL
      DATA ICALL/0/
      WRITE(34,892)
  892 FORMAT(/1X,60('*'),'call subroutine rfile in t2f_v2.f',60('*'))
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' RFILE 1.5, 1999.6.18: initializes permeability '
     x'modifiers and coordinate arrays, and optionally reads tables'
     x' with flowing wellbore pressures from disk files')
      SCALEA=SCALE*SCALE
      SCALEV=SCALEA*SCALE
C
      WRITE (34,901)
  901 FORMAT(1H )
c      OPEN(UNIT=17,FILE='TQLOS',STATUS='UNKNOWN')
C
C*****READ GEOMETRY DATA FROM FILE MESH OR MINC.***********************
c.....2-4-93: initialize flag for "dummy" connections w/ unknown blocks
      ndum=0
      im=10 ! read file *minc*
      IF(MC.EQ.0) THEN    
      IM=4
      ENDIF
C
      REWIND IM
      READ(IM,5020) DENT
      N=0
      N0=0
 1492 CONTINUE
      N=N+1
      READ(IM,1499) dent,MA1,MA2,EVOLx,AHTx,pmx,xx,yy,zz
c      READ(17,7572) amttc
 7572 FORMAT(E20.13)
 1499 FORMAT(A5,10X,A3,A2,6E10.4)
      IF(dent.EQ.'     ') GOTO 1502
      IF(N.LE.MNEL) GOTO 40
      WRITE (34,41) MNEL
   41 FORMAT(' NUMBER OF ELEMENTS SPECIFIED IN DATA BLOCK "ELEME"',
     X' EXCEEDS ALLOWABLE MAXIMUM OF PARAMETER *MNEL*=',I8,
     x' PLEASE RECOMPILE'/'*****SKIP FLOW SIMULATION*****')
      IS=IS+1
      RETURN

   40 CONTINUE
      elem(n)=dent                
      evol(n)=evolx               
      aht(n)=ahtx                 
      pm(n)=pmx                   
      x1(n)=xx
      x2(n)=yy
      x3(n)=zz
      amtt(n)=amttc
      IF(MA1.EQ.'   ') GOTO 15
C-----FIND MATERIAL INDEX
      MA12=MA1//MA2
      DO16 M=1,NM
      IF(MA12.EQ.MAT(M)) GOTO 17
   16 CONTINUE
      if(mop(7).ne.0) WRITE (34,18) MA12,ELEM(N)
   18 FORMAT(' REFERENCE TO UNKNOWN MATERIAL ',A5,' AT ELEMENT *',A5,
     X'* --- WILL IGNORE ELEMENT')
      n=n-1
      GOTO 1492
C
   17 MATX(N)=M
      GOTO 19
   15 CONTINUE
C-----COME HERE FOR ELEMENTS WITH DOMAIN NUMBER
      READ(MA2,'(I2)') MATX(N)
   19 CONTINUE
      IF(MATX(N).NE.0) GOTO 1498
C     COME HERE FOR ELEMENTS WITHOUT DOMAIN ASSIGNMENT, AND ASSIGN THEM TO DOMAIN # 1.
      MATX(N)=1
      WRITE (34,1497) ELEM(N),MAT(1)
 1497 FORMAT(' WARNING: ELEMENT ',A5,' HAS NO DOMAIN; ASSIGN TO DOMAIN',
     X' #1, *',A5,'*')
 1498 CONTINUE
C ========================================================================
      IF(EVOL(N).LE.0..AND.N0.EQ.0) N0=N-1
C     N0 IS THE INDEX PRECEDING THE FIRST ELEMENT WITH V = 0.
      IF(SCALE.EQ.1.) GOTO 1490
      EVOL(N)=EVOL(N)*SCALEV
c.....10-26-95: add proper scaling of heat transfer areas
      aht(n)=aht(n)*scalea
 1490 CONTINUE
      GOTO 1492
 1502 NEL=N-1
      NELA=NEL
      IF(N0.NE.0) NELA=N0
C     NELA IS THE NUMBER OF "ACTIVE" ELEMENTS; IT IS ASSUMED THAT
C     ALL ELEMENTS WITH INDEX N > NELA ARE BOUNDARY ELEMENTS.
C
      NST=0
      IF(ELST.EQ.'     ') GOTO 80
C
C from tough2 v1.6 (ysw & rick 7-13-95)
C
      elst='     '
c
c this is much better
      do i=1,nelist
        do n=1,nel
          IF(ELEM(N).EQ.eplist(i)) nstrick(i)=n
        enddo
      enddo
C
      DO 81 N=1,NEL
      IF(ELEM(N).EQ.ELST) GOTO 82
   81 CONTINUE
      GOTO 80
C
   82 NST=N
   80 CONTINUE
C
c.....FOFT-handling
      if(ioft.gt.0) then
c     come here to assign index numbers of elements for which time
c     plots are desired
      do 50 i=1,ioft                       
      do 51 n=1,nel
      if(elem(n).ne.eoft(i)) goto51
      noft(i)=n                           
      goto 50
   51 continue
c
      noft(i)=0
      iofu=iofu-1
      WRITE (34,52) eoft(i)
   52 format(' have encountered unknown element *',A5,'* in data block'
     x' "FOFT" - will ignore this element')
c
   50 continue
      endif
c
      READ(IM,5020) DENT
      NCON=0
C! ---GJM - Begin: determine active dimensions
      isumd1 = 0
      isumd2 = 0
      isumd3 = 0
C! ---GJM - End: determine active dimensions
 1493 CONTINUE
     
      IF(NCON.LT.MNCON) GOTO 42
      WRITE (34,43) MNCON
   43 FORMAT('NUMBER OF CONNECTIONS SPECIFIED IN DATA BLOCK "CONNE"',
     X' EXCEEDS ALLOWABLE MAXIMUM OF PARAMETER *MNCON*=',I8,'PLEASE'
     X' RECOMPILE'/'*****SKIP FLOW SIMULATION*****')
      IS=IS+1
      RETURN
   42 CONTINUE
C
      N=NCON+1
      READ(IM,1503) ELEM1(N),ELEM2(N),ISOX(N),DEL1(N),DEL2(N),AREA(N),
     BBETA(N),sig(n)
      
 1503 FORMAT(2A5,15X,I5,5E10.4)
      IF(ELEM1(N).EQ.'+++  ') GOTO 1504
      IF(ELEM1(N).EQ.'     ') GOTO 1507
      NCON=NCON+1
c
c ----------------
c! ..... GJM - Begin: Determining dimensions
c ----------------
c
         IF(isox(n).eq.1) isumd1 = isumd1+1
         IF(isox(n).eq.2) isumd2 = isumd2+1
         IF(isox(n).eq.3) isumd3 = isumd3+1
c
c ----------------
c! ..... GJM - End: Determining dimensions
c ----------------
c
      IF(SCALE.EQ.1.) GOTO 1491
      DEL1(N)=DEL1(N)*SCALE           
      DEL2(N)=DEL2(N)*SCALE
      AREA(N)=AREA(N)*SCALEA
 1491 CONTINUE
      GOTO 1493
C
 1504 READ(IM,1505) (NEX1(N),NEX2(N),N=1,NCON)        
      do33 j=1,ncon
      if(nex1(j).ne.0.and.nex2(j).ne.0) goto 33
      ndum=ndum+1
 1505 FORMAT(16I5)
   33 continue
      GOTO 1479
C
 1507 BACKSPACE IM
      WRITE(IM,1508)
 1508 FORMAT('+++  ')
c
c ----------------
c! ..... GJM - Begin: Determining dimensions
c ----------------
c
      IF(isumd1.NE.0.AND.isumd2.NE.0.AND.isumd3.NE.0) nactdi = 3
      IF(isumd2.EQ.0.OR.isumd3.EQ.0)  nactdi = 2
      IF(isumd2.EQ.0.AND.isumd3.EQ.0) nactdi = 1
c
      WRITE (34,6081) isumd1,isumd2,isumd3
 6081 FORMAT(' *The numbers of connections in the X, Y and Z directions'
     &' are (',i6,',',i6,',',i6,')*'/)
c
c ----------------
c! ..... GJM - End: Determining dimensions
c ----------------
c
      DO 26 J=1,NCON
      EL1=ELEM1(J)
      EL2=ELEM2(J)
      DO 27 N=1,NEL
      IF(ELEM(N).NE.EL1) GOTO 27
      NEX1(J)=N
      GOTO 28
   27 CONTINUE
      ndum=ndum+1
      NEX1(J)=0
      if(mop(7).ne.0) WRITE (34,127) EL1,J
  127 FORMAT(' REFERENCE TO UNKNOWN ELEMENT',A5,' AT CONNECTION',I7,
     A' --- WILL IGNORE THIS CONNECTION')
      goto 26
C
   28 DO 29 N=1,NEL
      IF(ELEM(N).NE.EL2) GOTO 29
      NEX2(J)=N
      GOTO 300
   29 CONTINUE
      ndum=ndum+1
      NEX2(J)=0
      if(mop(7).ne.0) WRITE (34,127) EL2,J
  300 CONTINUE
   26 CONTINUE
C
   22 WRITE(IM,1505) (NEX1(N),NEX2(N),N=1,NCON)
      ENDFILE IM
 1479 CONTINUE
c
c.....handling of "dummy" connections
      if(ndum.ne.0) then
c        come here when dummy connections are present
         n=0
      do 30 j=1,ncon
        if(nex1(j).ne.0.and.nex2(j).ne.0) then
          n=n+1
          elem1(n)=elem1(j)
          elem2(n)=elem2(j)
          nex1(n)=nex1(j)
          nex2(n)=nex2(j)
          isox(n)=isox(j)
          del1(n)=del1(j)
          del2(n)=del2(j)
          area(n)=area(j)
          beta(n)=beta(j)
          sig(n)=sig(j)
       endif
   30    continue
       WRITE (34,31) ncon,ndum,n
   31  format(I7,' connections read from file *MESH*, NDUM =',I7,
     x' connections make reference to unknown elements and will be'
     x' ignored'/' initialize a total of NCON =',I7,' connections to'
     x' the data arrays'/)
       ncon=n
      else
      WRITE (34,32) ncon
   32 format(' all NCON =',I7,' connections read from file *MESH*'
     x' reference known elements, and have been initialized to the'
     x' data arrays'/)
      endif
C
C-----END OF FILE MESH.------------------------------------------------
c
c.....COFT-handling: tabulate time-dependent connection data
      if(icoft.gt.0) then
c     come here to assign index numbers of connections for which time
c     plots are desired
      do 250 i=1,icoft
      do 251 n=1,ncon
      if(elem1(n)//elem2(n).ne.ecoft(i)) goto251
      ncoft(i)=n
      goto 250
  251 continue
c
      ncoft(i)=0
      icofu=icofu-1
      WRITE (34,252) ecoft(i)(1:5),ecoft(i)(6:10)
  252 format(' have encountered unknown connection (',A5,')-(',A5,')',
     x' in data block "COFT" - will ignore this connection')
  250 continue
      endif
c
c.....9-26-97: optional block-by-block permeability modification
      call pmin
      WRITE(34,*) '+++end of call subroutine pmin+++'
C     INITIALIZE DIFFUSIVE VAPOR FLUX AS ZERO.
      DO 25 J=1,NCON
      FVD(J)=0.d0                 
   25 CONTINUE
C
C*****READ SINKS/SOURCES FROM FILE GENER.******************************
C
      REWIND 3                    
      NOGN=0                      
      READ(3,5020,END=1415) DENT  
C5020 FORMAT(A5,A75)
      NGL=0
      kf=0
      kij=0
      ngkt=0  
C
 1483 CONTINUE
    
      IF(NOGN.LT.MNOGN) GOTO 44
      WRITE (34,45) MNOGN
   45 FORMAT(' NUMBER OF SINKS/SOURCES SPECIFIED IN DATA BLOCK "GENER"',
     X' EXCEEDS ALLOWABLE MAXIMUM, SKIP FLOW SIMULATION')
      IS=IS+1
      RETURN
   44 CONTINUE
C 
      N=NOGN+1
      READ(3,1481) ELEG(N),SOURCE(N),LTABG(N),type5(N)
     A,G(N),EG(N),HG(N),ktftb(n)
 1481 FORMAT(2A5,15X,I5,5X,A5,3G10.4,i2)
      IF(ELEG(N).EQ.'+++  ') GOTO 1482
      IF(ELEG(N).EQ.'     ') GOTO 1406
      NOGN=NOGN+1
C
      NGIND(N)=NGL
      LCOM(N)=0
      typen=type5(n)
      typ=typen(1:4)
      ty1=typen(1:1)
      itabg(n)=typen(5:5)
c
      if(ty1.eq.char(70).or.ty1.eq.char(102)) then
c.....come here when first character of a TYPE designation
c     is an upper or lower-case F; signifying a well on
c     deliverability with flowing wellbore pressure correction
c
      do 140 i=1,5
         if(typen(i:i).eq.' ') then
            ffi=typen(1:i-1)
            goto 141
         endif
  140 continue
         ffi=typen
  141 continue
c
      open(9,file=ffi,status='old',err=142)
      goto 149
c.....if file specified in the F-type GENER item is not
c     available; ignore the item
  142 continue
      WRITE (34,138) ffi
  138 format(' Invalid F-type GENER item;  file *',A5,'* is',
     x' not available')
      goto 139
c
  149 continue
c.....come here for valid F-type item
         LCOM(n)=nk1+2
         pi(n)=g(n)
c
c.....check whether the current wellbore pressure file has
c     already been read.
         if(kf.ge.1) then
            do 143 k=1,kf
            if(ffi.eq.fnam(k)) goto 144
  143       continue
            goto 145
  144       continue
c.....assign index of previously read wellbore pressure
c     file to current GENER item
            nftab(n)=k
            goto 1483
c
         endif
c
  145 continue
c.....read new wellbore pressure file
         kf=kf+1
         fnam(kf)=ffi
         nftab(n)=kf
         read(9,146) titl
  146    format(A80)
         read(9,147) ng,nh
         iftit(kf)=ng
         jftit(kf)=nh
  147    format(2I5)
         read(9,148) (pw(kij+i),i=1,ng)
  148    format(8E10.4)
         read(9,148) (pw(kij+ng+j),j=1,nh)
      do136 i=1,ng
  136 read(9,148) (pw(kij+ng+nh+(i-1)*nh+j),j=1,nh)
         ijf(kf)=kij
         kij=kij+ng+nh+ng*nh
      endif
c
c.....end of initializing F-type GENER items
c
      IF(typ.EQ.'WATR') typ='COM2'
      IF(typ.EQ.'NACL') typ='COM1'
      IF(typ.EQ.'HEAT') LCOM(N)=NK1
      IF(typ.EQ.'MASS') LCOM(N)=1
      IF(typ.EQ.'WATE') LCOM(N)=1
      IF(typ.EQ.'AIR ') LCOM(N)=2
      IF(typ.EQ.'COM1') LCOM(N)=1
      IF(typ.EQ.'COM2') LCOM(N)=2
      IF(typ.EQ.'COM3') LCOM(N)=3
      IF(typ.EQ.'COM4') LCOM(N)=4
      IF(typ.EQ.'COM5') LCOM(N)=5
      IF(typ.EQ.'COM6') LCOM(N)=6
c
      IF(typ.EQ.'DELV') then
C-----COME HERE FOR TYPE = DELV, IN WHICH CASE A PRODUCTIVITY INDEX
C     (PI) AND A FLOWING WELLBORE PRESSURE (PWB) ARE SPECIFIED.
         LCOM(N)=NK1+1
         PI(N)=G(N)
         PWB(N)=EG(N)
      endif
c
      IF(typ.EQ.'VOL.') then
C-----COME HERE IF VOLUMETRIC PRODUCTION RATE IS SPECIFIED.
         LCOM(N)=NK1+3
         GVOL(N)=G(N)
      endif
C
      IF (LCOM(N).NE.0) GO TO 13
  139 WRITE (34,14) typen,ELEG(N),SOURCE(N)
   14 FORMAT(' IGNORE UNKNOWN GENERATION OPTION *',A5,'* AT ELEMENT *',
     1  A5,'* SOURCE *',A5,'*')
      NOGN=NOGN-1
      GO TO 1483
   13 CONTINUE
C
C
      LTABA=ABS(LTABG(N))
c... new counter for k-thermal factors (els9/8/99)
        if(ktftb(n).gt.0)then
          ktaba = ktftb(n)
        endif
c
c Times and factors for effective thermal conductivity
        if(ktftb(n).gt.0)then
           read(3,36)(timkth(ngkt+L),L=1,ktaba)
           read(3,36)(fackth(ngkt+L),L=1,ktaba)
c           ngkt = ngkt + 1
           ngkt = ngkt + ktaba
        endif
c
      IF(LTABA.LE.1.OR.typ.EQ.'DELV') GOTO 1483
      IF(NGL+LTABA.LE.MGTAB) GOTO 46
      WRITE (34,47) MGTAB
   47 FORMAT(' NUMBER OF TABULAR GENERATION DATA SPECIFIED IN DATA',
     X' BLOCK "GENER" EXCEEDS ALLOWABLE MAXIMUM OF ',I7,' INCREASE'
     X' PARAMETER *MGTAB* IN MAIN PROGRAM, SKIP FLOW SIMULATION')
      IS=IS+1
      RETURN
   46 CONTINUE
C
      READ(3,36) (F1(NGL+L),L=1,LTABA)
      READ(3,36) (F2(NGL+L),L=1,LTABA)
      IF(itabg(n).NE.' ') READ(3,36) (F3(NGL+L),L=1,LTABA)
      NGL=NGL+LTABA
      GOTO 1483
 1482 CONTINUE
      IF(NOGN.EQ.0) GOTO 1415
C
      READ(3,1505) (NEXG(N),N=1,NOGN)
      GOTO 1415
C
 1406 if(nogn.eq.0) goto 1415
      BACKSPACE 3
      WRITE(3,1508)
      DO 37 N=1,NOGN
      EL=ELEG(N)
      DO 38 J=1,NEL
      IF(ELEM(J).NE.EL) GOTO 38
      NEXG(N)=J
      GOTO 37
   38 CONTINUE
      NEXG(N)=0
      WRITE (34,128) EL,N
  128 FORMAT(' REFERENCE TO UNKNOWN ELEMENT *',A5,'* AT SOURCE #',I5,
     A' --- WILL IGNORE')
   37 CONTINUE
C
   39 WRITE(3,1505) (NEXG(N),N=1,NOGN)
      ENDFILE 3
 1415 CONTINUE
c
c.....GOFT-handling (tabulating generation data)
      if(igoft.gt.0) then
         do 150 i=1,igoft
         do 151 n=1,nogn
c     WRITE (34,183) i,n,egoft(i),eleg(n)
  183 format(' RFILE I =',I4,' N =',I4,' ELEG = ',A5,' EGOFT = ',A5)
         if(eleg(n).ne.egoft(i).or.nexg(n).eq.0) goto 151
         ngoft(i)=n
         goto 150
  151    continue
         ngoft(i)=0
         igofu=igofu-1
         WRITE (34,152) egoft(i)
  152 format(' element *',A5,'* in data block "GOFT" does not',
     x' have a valid "GENER" item - ignore this and proceed')
  150    continue
      endif
C
C-----END OF FILE GENER.-----------------------------------------------
c     WRITE (34,137) (pw(ij),ij=1,kij)
  137 format(' PW-array ======================================='/
     x(10(3X,E10.4)))
C
C*****READ INITIAL CONDITIONS FROM FILE INCON.*************************
C
 1300 REWIND 1            ! #1 <INCON>        
      IF(NI.EQ.0) GOTO 1412
C
C-----COME HERE FOR START.
      DO 1409 N=1,NEL
c Added definition of nmat for speedup
        nmat = matx(n)
      PHI(N)=POR(nmat)
c First set permeabilities to those in rocks block
      perm(1,n)=per(1, nmat)
      perm(2,n)=per(2, nmat)
      perm(3,n)=per(3, nmat)
C------------------------------------------------------------------------
C
      NLOC=(N-1)*nkin1
C
      IF(YIN(1,MATX(N)).NE.0.d0) GOTO 1410
      DO 2000 I=1,nkin1
 2000 X(NLOC+I)=DEP(I)
c
C----------------------------------------------------For using EOS9
      TC_DUM=DEP(2)
      IF (DEP(2).EQ.0.0D0) TC_DUM=Tcref
      Tc_EOS9(N)=TC_DUM    ! initial temperature (oC)
c------------------------------------------------------------------
c
      GOTO 1409
C
 1410 CONTINUE
      DO 1411 I=1,nkin1
 1411 X(NLOC+I)=YIN(I,MATX(N))
C
 1409 CONTINUE
      READ(1,5020,END=1703) DENT
C
C-----COME HERE FOR ASSIGNING INITIAL CONDITIONS.
C
      iun=0
 2002 CONTINUE
C
C--------------------------------------------------------------------------
      READ(1,1403,END=1700) EL,PORX,perx,pery,perz
 1403 FORMAT(A5,10X,4E15.8,1E15.4)
C 
C--------------------------------------------------------------------------
C
      IF(EL.EQ.'     ') then
         if(iun.ne.0) WRITE (34,1407) iun
 1407    format(/' WARNING: INCON at ',I5,' elements have been ignored')
         RETURN
      endif
      IF(EL.EQ.'+++  ') GOTO 2004
C
      READ(1,1402,END=1700) (DEPU(I),I=1,nkin1)
 1402 FORMAT(4E20.13)
C
      DO 2005 J=1,NEL
      IF(ELEM(J).NE.EL) GOTO 2005
      JLOC=(J-1)*nkin1
!
!.....TX, 02/14/2013 modification:
!.....Take porosity and permeabilities from INCON file for default case: MOPR(17)=0
!.....Otherwise, skip
!
      if (mopr(17) .eq. 0) then
         IF(PORX.NE.0.d0) PHI(J)=PORX
         IF(perx.NE.0.D0) perm(1,j)=perx
         IF(pery.NE.0.D0) perm(2,j)=pery
         IF(perz.NE.0.D0) perm(3,j)=perz
      end if
!
C----------------------------------------------------------------------------
C
      DO 2006 I=1,nkin1
 2006 X(JLOC+I)=DEPU(I)
c
C--------------------------------------------------------------For using EOS9
      TC_DUM=DEPU(2)
      IF (DEPU(2).EQ.0.0D0) TC_DUM=Tcref
       Tc_EOS9(J)=TC_DUM    ! initial temperature (oC)
c----------------------------------------------------------------------------
c
      GOTO 2002
 2005 CONTINUE
C
      iun=iun+1
      if(iun.le.50) WRITE (34,1408) EL
 1408 FORMAT(50H WILL IGNORE INITIAL CONDITION AT UNKNOWN ELEMENT ,A7)
C
      GOTO 2002
C
 2004 continue
         if(iun.ne.0) WRITE (34,1407) iun
      READ(1,74) KCYC,ITERC,nmx,TSTART,TIMIN
c
c.....1-16-92: interrogate number of domains fed through *INCON*
c              always use NM from input data
      if(nmx.ne.nm) then
      WRITE (34,2008) nmx,nm
 2008 FORMAT(' WARNING: number of domains read from *INCON* is',I3,
     x', different from NM =',I3,'provided through input data, the'
     x' latter value will be used'/)
      endif
      IF(TIMAX.NE.0.d0.AND.TIMIN.GE.TIMAX) TIMIN=0.d0
      IF(MOP(19).EQ.0) RETURN
      WRITE (34,2007) MOP(19)
 2007 FORMAT(' SUBROUTINE RFILE: HAVE ENCOUNTERED "+++" IN *INCON* ---'
     x' RESET MOP(19) =',I2,' TO DEFAULT VALUE 0.'/)
      MOP(19)=0
      RETURN
C
 1412 CONTINUE
C-----COME HERE TO READ RESTART INFORMATION.
      READ(1,5020,END=1700) DENT
      IF(DENT.NE.'INCON') GOTO 1700
      DO 1702 N=1,NEL
      NLOC=(N-1)*nkin1
C
C-----------------------------------------------------------------------------
      READ(1,1403,END=1700) EL,PHI(N),perm(1,n),perm(2,n),perm(3,n)
      READ(1,1402,END=1700) (X(NLOC+I),I=1,nkin1)

C---------------------------------------------------------------For using EOS9
        Tc_EOS9(N)=X(NLOC+2)    ! initial temperature (oC)
c-----------------------------------------------------------------------------
c
 1702 CONTINUE
C
      READ(1,5020) DENT
      IF(DENT.EQ.'     ') RETURN
      READ(1,74) KCYC,ITERC,nmx,TSTART,TIMIN
      if(nmx.ne.nm) then
      WRITE (34,2008) nmx,nm,nm
      endif
   74 FORMAT(2I10,I5,2E15.8)
      IF(TIMAX.NE.0..AND.TIMIN.GE.TIMAX) TIMIN=0.d0
      IF(MOP(19).EQ.0) RETURN
      WRITE (34,2007) MOP(19)
      MOP(19)=0
   75 FORMAT(5E15.8)
C
C-----END OF FILE INCON.-----------------------------------------------
C
 1703 RETURN
 5020 FORMAT(A5)
   36 FORMAT(4E14.7)
C
 1700 WRITE (34,1701) DENT
 1701 FORMAT(' IMPROPER FORMAT ON FILE *INCON* --- IDENTIFIER "',A5,'"')
      IS=IS+1
      RETURN
      END
c
      subroutine pmin
c.....initialize block-by-block permeability modifiers, and
c     generate informative printout.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      common/E7/pm(mnel)
C
C$$$$$$$$$ COMMON BLOCKS FOR ROCK PROPERTIES $$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)     ! ??????
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +             GK(MAXMAT)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/SVZ/NOITE,MOP(24)
      common/ran2/iran                ! ???
      common/ff/h1                    ! ???
      character*1 h1                  ! ???
c
      character mat*5
      save icall
      data icall/0/
      ICALL=ICALL+1
      WRITE (34,*) '+++subroutine pmin in t2f_v2.f is called+++'
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' PMIN 1.0, 1997.9.26: initialize block-by-block'
     x' permeability modifiers')
c
c.....4-18-94: use array PM as random permeability modifiers.
      iran=0
      do 28 n=1,nm
      if(mat(n).eq.'SEED ') goto 222
   28 continue
      WRITE (34,47)
   47 format(' domain = "SEED " is not present, no permeability'
     x' modification will be made')
      goto 27
c
  222 continue
      WRITE (34,42)
   42 format(' domain = "SEED " is present, permeability modification'
     x' will be made'/)
      if(dm(n).ne.0.d0) then
c.....use "linear" permeability modifiers
      iran=2
      s=dm(n)
      elseif(por(n).ne.0.d0) then
c.....use "logarithmic" permeability modifiers
      iran=3
      s=por(n)
      else
c.....use externally supplied PM-data as permeability
c     modifiers
      iran=1
      WRITE (34,38)
   38 format(' Option 1: Externally supplied'/)
      goto 37
      endif
c
c.....optional scale factor
      sran=1.d0
      if(per(1,n).ne.0.d0) sran=per(1,n)
      WRITE (34,223) nel,s,sran
  223 format(' Generate random permeability modifiers for',I5,' grid'
     x' blocks with seed S = ',E12.6,', scale factor SRAN = ',E12.6)
      if(iran.eq.2) WRITE (34,29)
      if(iran.eq.3) WRITE (34,230)
   29 format(' Option 2: Linear modification'/)
  230 format(' Option 3: Logarithmic modification'/)
c
      do 24 i=1,nel
c Conflicts with intrinsic      call rand(s)
      if(i.eq.1) write(34,1777) nel
 1777 format(' subroutine t2rand is called',I5,'times in pmin')
      call t2rand(s)
      if(iran.eq.2) pm(i)=sran*s
      if(iran.eq.3) pm(i)=exp(-sran*s)
   24 continue
c
   37 continue
c.....5-17-94
      if (per(2,n).ne.0.d0) then
         im=0
         do 40 i=1,nel
            pm(i)=pm(i)-per(2,n)
            if(pm(i).lt.0.d0) then
               im=im+1
               pm(i)=0.d0
            endif
   40    continue
c)
         pim=100.d0*dble(im)/dble(nel)
         WRITE (34,41) nel,im,pim
   41    format(' of NEL = ',I5,' grid blocks, a total',
     x   ' of IM = ',I5,' or ',F5.2,' % is impermeable'/)
      endif
   27 return
      end
c
c  Conflicts with intrinsic  subroutine rand(x)
      subroutine t2rand(x)
c.....11-29-93: from Meissner, Organick, FORTRAN77 (1984), p. 335.
c Added implicits
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      save icall
      data k,j,m,rm/5701,3612,566927,566927.0d0/
      data icall/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***RAND 1.0, 1993.11.29: random numbers (cf. Meissner &'
     x' Organick, FORTRAN77, Addison-Wesley, 1984, p.335***')
      ix=int(x*rm)
      irand=mod(j*ix+k,m)
c      x=(real(irand)+0.5d0)/rm
      x=(dble(irand)+0.5d0)/rm
      return
      end
c
      SUBROUTINE INDATA
C
C-----THIS ROUTINE PROVIDES A PRINTOUT OF MOST OF THE INPUT DATA.
C     IT IS CALLED ONLY IF MOP(7).NE.0.
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)
C
C$$$$$$$$$ COMMON BLOCKS FOR PRIMARY VARIABLES $$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NEQ
C
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P3/DELX((MNK+1)*MNEL)
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/P5/DOLD(MNEQ*MNEL)
C
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)
      COMMON/C4/DEL2(MNCON)
      COMMON/C5/AREA(MNCON)
      COMMON/C6/BETA(MNCON)
      COMMON/C7/ISOX(MNCON)
      COMMON/C8/GLO(MNCON)
      COMMON/C9/ELEM1(MNCON)
      COMMON/C10/ELEM2(MNCON)
      common/c12/sig(mncon)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/TITLE/TITLE
      CHARACTER*80 TITLE
C
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
C
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G1/F1(MGTAB)
      COMMON/G2/F2(MGTAB)
      COMMON/G3/F3(MGTAB)
      COMMON/G4/ELEG(MNOGN)
      COMMON/G5/SOURCE(MNOGN)
      COMMON/G6/LTABG(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G10/ITABG(MNOGN)
      COMMON/G11/NGIND(MNOGN)
      COMMON/G12/LCOM(MNOGN)
      COMMON/G13/PI(MNOGN)
      COMMON/G14/PWB(MNOGN)
      COMMON/G15/HG(MNOGN)
      COMMON/G16/GPO(MNOGN)
      COMMON/G17/SDENS(MNOGN)
      COMMON/G18/SSAT(MNOGN)
      COMMON/G19/GVOL(MNOGN)
      COMMON/G20/HL(MNOGN)
      COMMON/G21/HS(MNOGN)
      COMMON/G22/QVGC(MNOGN)
      COMMON/G23/QVWC(MNOGN)
      COMMON/G24/QVOC(MNOGN)
      COMMON/G25/GRAD(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
C
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +              SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +             GK(MAXMAT)
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
      COMMON/STEP/ELST
      COMMON/STE1/NST
      COMMON/DEFINI/DEP(10)
      COMMON/DG/WUP,WNR
      COMMON/DFM/TIMAX,REDLT
      COMMON/DLT/NDLT,DLT(100)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
C      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO,TIMP1
      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO
      COMMON/DOP/ENTH,KDATA,QUAL
      COMMON/DX/K,NI,SCALE
      COMMON/POV6/TSTART
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/PATCH/SING
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BIND/DIFF0,TEXP,BE
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      common/ff/h1
      character*1 h1
      CHARACTER*5 ELEM1,ELEM2,ELEM,ELEG,SOURCE,MAT,ELST
      CHARACTER ITABG*1
c     Added for correct variable size
      real*4 rtzero
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      WRITE(34,892)
  892 FORMAT(/1X,60('*'),'call subroutine indata in t2f_v2.f',60('*'))
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' INDATA 1.0, 1991.3.5: READ DATA PROVIDED THROUGH DISK '
     x'FILES')
C
      h1=char(12)
      WRITE (34,100) h1
      WRITE (34,150) TITLE
      WRITE (34,200)
      WRITE (34,250)
      WRITE (34,300) NOITE,KDATA,MCYC,MSEC,MCYPR,(MOP(I),I=1,24),
     XDIFF0,TEXP,BE
      WRITE (34,350)
      WRITE (34,401) TSTART,TIMAX,DELTEN,DELTMX,ELST,GF,REDLT,SCALE
  401 FORMAT(4(5X,E10.4),10X,A5,3(5X,E10.4))
      IF(NDLT.EQ.0) GOTO 10
      WRITE (34,450)
      WRITE (34,451) ((DLT(J+8*(N-1)),J=1,8),N=1,NDLT)
  451 FORMAT(5X,8(5X,E10.4))
      GO TO 20
   10 WRITE (34,500)
   20 WRITE (34,550)
      WRITE (34,400) RE1,RE2,U,WUP,WNR,DFAC,FOR
      WRITE (34,600)
      WRITE (34,650) DEP(1),DEP(2),DEP(3),DEP(4)
      WRITE (34,700)
      DO 1 I=1,NM
      WRITE (34,750)
      WRITE (34,800)
     AI,MAT(I),DM(I),POR(I),CWET(I),SH(I),COM(I),EXPAN(I)
      WRITE (34,810)
      WRITE (34,820) PER(1,I),PER(2,I),PER(3,I)
    1 CONTINUE
      WRITE (34,830)
      WRITE (34,840)
      DO 2 I=1,NEL
      WRITE (34,850) ELEM(I),MATX(I),EVOL(I)
    2 CONTINUE
      WRITE (34,860)
      WRITE (34,870)
      DO 3 I=1,NCON
      WRITE (34,880) ELEM1(I),ELEM2(I),ISOX(I),DEL1(I),DEL2(I),AREA(I),
     ABETA(I)
    3 CONTINUE
      IF(NOGN .LE. 0 ) GOTO 30
      WRITE (34,890)
      IDELV=0
      DO 4 I=1,NOGN
      IF(LCOM(I).EQ.NK1+1) GOTO 6
      WRITE (34,895)
      WRITE (34,900) ELEG(I),SOURCE(I),G(I),EG(I)
      IF(LTABG(I) .GT. 1)GO TO 40
      WRITE (34,896)
      GO TO 50
   40 WRITE (34,897)
   50 CONTINUE
      GOTO 4
    6 CONTINUE
      IDELV=IDELV+1
      IF(IDELV.EQ.1) LAY=LTABG(I)
      IF(IDELV.EQ.1) WRITE (34,885) LAY
  885 FORMAT(5X,43HWELL ON DELIVERABILITY   $$$$$$$$   OPEN IN,I3,19H  L
     AAYERS   $$$$$$$$//5X,'   ELEMENT         SOURCE           PI',
     B'             PWB            DEL(Z)'/)
      WRITE (34,900) ELEG(I),SOURCE(I),PI(I),PWB(I),HG(I)
    4 CONTINUE
   30 WRITE (34,910)
      WRITE (34,920)
      DO 5 I=1,NEL
      WRITE (34,930) ELEM(I),PHI(I),(X((I-1)*NK1+J),J=1,NK1)
    5 CONTINUE
      WRITE (34,940)
      RETURN
  100 FORMAT(A1/' TOUGH2 INPUT DATA')
  150 FORMAT(' PROBLEM TITLE:' ,A80/)
  200 FORMAT(/' PROBLEM SPECIFICATIONS:'/)
  250 FORMAT('     NOITE     KDATA      MCYC      MSEC     MCYPR
     A      MOP                   DIFF0          TEXP            BE')
  300 FORMAT(5(4X,I5,1X),6X,24I1,3(5X,E10.4),/)
  350 FORMAT(117H       TSTART         TIMAX          DELTEN         DEL
     ATMX            ELST        GF             REDLT         SCALE )
  400 FORMAT(5X,8(5X,E10.4))
  450 FORMAT(/,5X,35H VARIABLE TIME STEPS ARE PRESCRIBED)
  500 FORMAT(/,5X,45H A CONSTANT TIME STEP OF DELTEN IS PRESCRIBED,/)
  550 FORMAT('            RE1            RE2            U              W
     AUP            WNR            DFAC           FOR')
  600 FORMAT(/79H             DEP(1)              DEP(2)              DE
     AP(3)              DEP(4))
  650 FORMAT(5X,4(5X,E15.8))
  700 FORMAT(/,5X,15HROCK PROPERTIES,/)
  750 FORMAT(5X,102HDOMAIN     MAT        DENSITY        POROSITY     CO
     ANDUCTIVITY     HEAT CAP       COMPR          EXPAN)
  800 FORMAT(5X,I4,6X,A5,6(5X,E10.4),/)
  810 FORMAT(2X,45H          PERM1          PERM2          PERM3)
  820 FORMAT(5X,3(5X,E10.4),/)
  830 FORMAT(/,5X,8HELEMENTS,/)
  840 FORMAT(5X,45H        ELEMENT       MATERIAL         VOLUME)
  850 FORMAT(14X,A5,8X,I5,10X,E10.4)
  860 FORMAT(/,5X,11HCONNECTIONS,/)
  870 FORMAT(5X,'          ELEM1          ELEM2           ISOT    ',
     A'       DEL1           DEL2           AREA           BETA'/)
  880 FORMAT(5X,2(10X,A5),9X,I5,6X,4(5X,E10.4))
  890 FORMAT(/,5X,15HGENERATION DATA)
  895 FORMAT('        ELEMENT         SOURCE           RATE',
     A'           ENTHALPY'/)
  896 FORMAT(5X,26HCONSTANT RATE IS SPECIFIED,/)
  897 FORMAT(5X,26HVARIABLE RATE IS SPECIFIED,/)
  900 FORMAT(5X,5X,A5,10X,A5,4X,4(5X,E10.4))
  910 FORMAT(/,5X,18HINITIAL CONDITIONS,/)
  920 FORMAT(5X,84H        ELEMENT       POROSITY       X1             X
     A2             X3             X4,/)
  930 FORMAT(5X,10X,A5,5(5X,E10.4))
  940 FORMAT(/,5X,17HEND OF INPUT DATA,/)
      END
C
      SUBROUTINE QU
C
C-----THIS SUBROUTINE COMPUTES ALL TERMS ARISING FROM SINKS AND SOURCES.
C
C          *****************************************************
C          *              modular version of 1998              *
C          *****************************************************
c
c     generation TYPE is processed in RFILE for checking and
c     assignment of major generation options through array LCOM
c
c     TYPE is conveyed to QU for sub-choices, especially the
c     assignment of file names to tables of flowing bottomhole
c     pressures, for geothermal production at specified wellhead
c     conditions
c
c     a new subroutine TTAB was written to perform interpolation
c     on time-dependent tables
c
c     new subroutines PHAS and PHASD were written to determine
c     composition and flowing enthalpy of multiphase fluids
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/L3/CO(mnz+1)     ! 就这个CO干啥用的我到现在都不知道...
      COMMON/DM/DELTEN,DELTEX,FOR,FORD    ! 时间步长的长度,时间步长的上限,后面的俩就不知道了
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL) ! 主要变量
      COMMON/SVZ/NOITE,MOP(24)        ! 每个时间步内有多少个迭代步(默认是8),24个MOP
      COMMON/ BCIJT /    IHD,    MAX,    MID,    SET
      COMMON/KONIT/KON,DELT,IGOOD     ! KON是收敛标志,DELT是时间步长(初始值等于DELTEN)
C
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G1/F1(MGTAB)             ! 抽出点的时间列表 
      COMMON/G2/F2(MGTAB)             ! 抽出点的流速列表
      COMMON/G3/F3(MGTAB)             ! 注入点的焓值列表
      COMMON/G4/ELEG(MNOGN)           ! 网格名
      COMMON/G5/SOURCE(MNOGN)         ! 源汇项的名
      COMMON/G6/LTABG(MNOGN)          ! 时间列表个数
      COMMON/G7/G(MNOGN)              ! 流量
      COMMON/G8/EG(MNOGN)             ! 焓值流量
      COMMON/G9/NEXG(MNOGN)           ! ???  
      COMMON/G10/ITABG(MNOGN)         ! itabg := ITAB 如果不为空则读取焓值列表
      COMMON/G11/NGIND(MNOGN)         ! ???
      COMMON/G12/LCOM(MNOGN)          ! 产流的TYPE在模块RFILE中被处理以通过数组LCOM检测和声明产流项的类别
      COMMON/G13/PI(MNOGN)            ! ??
      COMMON/G14/PWB(MNOGN)           ! 注入的层厚(应该是开孔段厚度之类的)
      COMMON/G15/HG(MNOGN)
      COMMON/G16/GPO(MNOGN)
      COMMON/G17/SDENS(MNOGN)
      COMMON/G18/SSAT(MNOGN)
      COMMON/G19/GVOL(MNOGN)
      COMMON/G20/HL(MNOGN)
      COMMON/G21/HS(MNOGN)
      COMMON/G22/QVGC(MNOGN)
      COMMON/G23/QVWC(MNOGN)
      COMMON/G24/QVOC(MNOGN)
      COMMON/G25/GRAD(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/E3/EVOL(MNEL)                ! 网格体积
      COMMON/P1/X((MNK+1)*MNEL)           ! X是迭代步之前的变量值
      COMMON/P2/DX((MNK+1)*MNEL)          ! 最新的主要变量的增量
      COMMON/P3/DELX((MNK+1)*MNEL)        ! 数值微分过程中的主要变量的增量小量
      COMMON/P4/R(MNEQ*MNEL+1)            ! 残差
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/DFM/TIMAX,REDLT
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      COMMON/BC/NELA
      common/dkm/d(11,12)
      common/source_type/type5(mnogn)
      character*5 type5
      CHARACTER ELEG*5,SOURCE*5,ITABG*1,ITABA*1,typen*5
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **QU 1.1, 1998.1.23: ASSEMBLE SOURCE AND SINK TERMS,'
     x'"rigorous" step rate capability for MOP(12) = 2, and capability'
     x' for flowing wellbore pressure corrections******')
C
      IF(MOP(4).GE.1) WRITE (34,101) KCYC,ITER
  101 FORMAT(/' SUBROUTINE QU --- [KCYC,ITER] = [',I4,',',I3,']')
C
      DO 3901 N=1,NOGN
C
      J=NEXG(N)
      IF(J.EQ.0.OR.J.GT.NELA) GOTO 3900
      FAC=FORD/EVOL(J)
C
      JLOC=(J-1)*NEQ
      JLOCP=(J-1)*NK1
C
      MN=LCOM(N)
C
      JMN=JLOC+MN
      JNK1=JLOC+NK1
      J2LOC=(J-1)*NSEC*NEQ1
C
      LTABA=ABS(LTABG(N))
      typen=type5(n)
      ITABA=typen(5:5)
c
c----- now jump to proper sub-section
c
      if(mn.le.nk) goto 1
      if(mn.eq.nk1) goto 2
      lnk=mn-nk1
c
c***** choose generation option
      goto(3,4,5,6),lnk
c
    1 continue
C     ************************************************************
C     *  come here when production or injection of a mass        *
C     *  component at prescribed rate is specified (LTABA = 1)   *
C     ************************************************************
c
      if(ltaba.le.1) then
         gn=g(n)
         gpo(n)=g(n)
         if(gn.ge.0.d0.or.(gn.lt.0.d0.and.itaba.ne.' ')) then
c...........mass injection (GN > 0), or production (GN < 0)
c           with prescribed enthalpy
            r(jmn)=r(jmn)-fac*gn
            if(neq.eq.nk1) then
              if(eg(n).le.100.d0) then
                DO 1051 MM=1,NEQ1
                NGJNK1=J2LOC+(MM-1)*NSEC
                DO 1050 K=1,NK1
 1050           D(K,MM)=0.D0
                D(NK1,MM)=D(NK1,MM)-FAC*GN*(eg(n)-
     A          PAR(NGJNK1+NSEC-1))*4.2d3
 1051           CONTINUE
                GOTO 3898
              else
                R(JNK1)=R(JNK1)-FAC*G(N)*EG(N)
              end if
            end if
            goto 3900
         else
c...........mass production, with enthalpy to be determined
c           from conditions in producing block
            call phas(n,j2loc,gn,fac)
            goto 3898
C
         endif
      else
C     ************************************************************
C     *  come here when a table of time-dependent generation     *
C     *  data is specified (LTABA > 1)                           *
C     ************************************************************
c
      call ttab(n,ltaba,itaba,gn,egn)
      if(igood.ne.0) return
c
      g(n)=gn
      gpo(n)=gn
c Changed to allow non-tabular constant enthalpies
c      if(itaba.ne.' ') then
c          eg(n)=egn
      if(itaba.ne.' ') eg(n)=egn
      if(gn.ge.0.d0.or.(gn.lt.0.d0.and.itaba.ne.' ')) then
c...........mass injection (GN > 0), or production (GN < 0)
c           with prescribed enthalpy
            r(jmn)=r(jmn)-fac*gn
            if(neq.eq.nk1) R(JNK1)=R(JNK1)-FAC*G(N)*EG(N)
            goto 3900
      else
         call phas(n,j2loc,gn,fac)
         goto 3898
      endif
C
C
      endif
c
    2 continue
C     ************************************************************
C     *  come here when withdrawal or injection of heat          *
C     *  at prescribed rate is specified                         *
C     ************************************************************
c
c
C-----IGNORE HEAT SINKS/SOURCES WHEN NOT SOLVING ENERGY EQUATION.
      if(neq.lt.nk1) goto 3900
c
      if(ltaba.gt.1) then
C*******************************************************
C* come here when a table of time-dependent generation *
C* data is specified (LTABA > 1) *
C************************************************************
c
         call ttab(n,ltaba,itaba,gn,egn)
         if(igood.ne.0) return
         g(n)=gn
       endif
c
      gpo(n)=g(n)
      r(jmn)=r(jmn)-fac*g(n)
c
      goto 3900
c
    3 continue
C     ************************************************************
C     *  deliverability option: production well with specified   *
C     *  downhole pressure; a simple gravity correction for      *
C     *  flowing wellbore pressure may be performed*             *
C     ************************************************************
c
c.....optionally perform simple gravity correction for flowing
c     downhole pressure of multi-feedzone wells.
      if(iter.eq.1.and.ltaba.gt.1) call gcor(n,ltaba)
c
c.....obtain source rate, phase composition, flowing enthalpy
      call phasd(n,jlocp,j2loc,fac)
c
      gpo(n)=g(n)
      goto 3898
c
    4 continue
C     ************************************************************
C     *  production well with specified wellhead pressure, and   *
C     *  flowing wellbore pressure correction                    *
C     ************************************************************
c
c.....on first call, initialize starting guess for well rate
      if(icall.eq.1) sdens(n)=10.d0
c
c.....obtain source rate, phase composition, flowing enthalpy
      call phasf(n,jlocp,j2loc,fac)
c
      gpo(n)=g(n)
      goto 3898
c
    5 continue
C-----COME HERE IF VOLUMETRIC PRODUCTION RATE IS SPECIFIED.
c
    6 continue
C
C
 3898 continue
C
C-----ASSIGN TERMS FOR LINEAR EQUATIONS.
C
      DO 105 K=1,NEQ
      R(JLOC+K)=R(JLOC+K)+D(K,1)
C
      DO 106 L=1,NEQ
      JKL=(J-1)*NEQ*NEQ+(K-1)*NEQ+L
      CO(JKL)=CO(JKL)-(D(K,L+1)-D(K,1))/DELX(JLOCP+L)
  106 CONTINUE
C
  105 CONTINUE
C
 3900 CONTINUE
C
      IF(MOP(4).GE.2) WRITE(34,102) ELEG(N),SOURCE(N),G(N),EG(N),
     A(FF((N-1)*NPH+NP),NP=1,NPH)
  102 FORMAT(9H ELEMENT ,A5,10H   SOURCE ,A5,21H   ---   FLOW RATE = ,E1
     A2.6,23H   SPECIFIC ENTHALPY = ,E12.6/8H FNP1 = ,E12.6,10H   FNP2 =
     B ,E12.6,10H   FNP3 = ,E12.6,10H   FNP4 = ,E12.6)
C
 3901 CONTINUE
      RETURN
      END
c
c
c
      subroutine phas(n,j2loc,gn,fac)
c
c-----this subroutine calculates phase composition and flowing
c     enthalpy of multiphase source fluid.
c
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/G8/EG(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/dkm/d(11,12)
!
!.....Modify for H2 generation problems
!
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
!............................................
c
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *PHAS 1.0, 1999.11.5: calculate composition and enthalpy'
     x' of source fluid**********')
c
      ic=mop(9)+1
C-----LOOP OVER NEQ1 SETS OF PARAMETERS.
      DO 101 M=1,NEQ1
      JLM2=J2LOC+(M-1)*NSEC
C
C-----ZERO OUT BLOCK (J,J).
      DO 1011 K=1,NK1
 1011 D(K,M)=0.d0
C
      FFS=0.d0
C-----FFS IS THE SUM OVER MOBILITY*DENSITY IN ALL PHASES.
C-----NOW COMPUTE FRACTIONAL FLOWS IN THE VARIOUS PHASES.
C
      DO 30 NP=1,NPH
      J2LNP=JLM2+(NP-1)*NBK
c
      FF((N-1)*NPH+NP)=0.d0
c
      IF(IC.EQ.1.AND.PAR(J2LNP+3).NE.0.d0)
     AFF((N-1)*NPH+NP)=PAR(J2LNP+2)*PAR(J2LNP+4)/PAR(J2LNP+3)
c
      IF(IC.EQ.2) FF((N-1)*NPH+NP)=PAR(J2LNP+1)*PAR(J2LNP+4)
      FFS=FFS+FF((N-1)*NPH+NP)
c
      IF(MOP(4).GE.5) WRITE (34,3) M,NP,PAR(J2LNP+2),PAR(J2LNP+3),
c     APAR(J2LNP+4),FF((N-1)*NPH+NP),VOLR,FFS,DELP
c    3 FORMAT(3H M=,I2,4H NP=,I1,3H K=,E12.6,5H VIS=,E12.6,3H D=,E12.6,4H
c     A FF=,E12.6,6H VOLR=,E12.6,5H FFS=,E12.6,6H DELP=,E12.6)
     APAR(J2LNP+4),FF((N-1)*NPH+NP),FFS
    3 FORMAT(3H M=,I2,4H NP=,I1,3H K=,E12.6,5H VIS=,E12.6,3H D=,E12.6,4H
     A FF=,E12.6,5H FFS=,E12.6)
C
   30 CONTINUE
C
      EG(N)=0.d0
C
C-----RENORMALIZE FRACTIONAL FLOWS TO 1, AND COMPUTE CONTRIBUTIONS
C     OF COMPONENTS IN PHASES.
C
      DO 31 NP=1,NPH
      J2LNP=JLM2+(NP-1)*NBK
      IF(FFS.NE.0.d0)
     AFF((N-1)*NPH+NP)=FF((N-1)*NPH+NP)/FFS
      EG(N)=EG(N)+FF((N-1)*NPH+NP)*PAR(J2LNP+5)
!
!------------
!.....Modify for NAGRA H2 generation and H2O consumption 
!.....due to iron corrsion (EOS5)
!------------
!
Crsb  Specify H2O only extraction  
!
cc      if (ieos .eq. 5 .and. source(n).eq.'H2O 1')   then
cc        if (np.eq.1) FF((N-1)*NPH+NP)=0.D0
cc        if (np.eq.2) FF((N-1)*NPH+NP)=1.D0
cc      end if
!
Crse
!
!-------------
!
      DO 32 K=1,NK
      D(K,M)=D(K,M)-FAC*GN*FF((N-1)*NPH+NP)*PAR(J2LNP+NB+K)
   32 CONTINUE
C
   31 CONTINUE
C
      D(NK1,M)=D(NK1,M)-FAC*GN*EG(N)
C
  101 CONTINUE
      return
      end
c
      subroutine phasd(n,jlocp,j2loc,fac)
c
c-----this subroutine calculates phase composition, flow rate,
c     and flowing enthalpy for wells on deliverability
c
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P3/DELX((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G13/PI(MNOGN)
      COMMON/G14/PWB(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/dkm/d(11,12)
      common/ech/eosn(20)
      character*10 eosn
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***PHASD 1.0, 2001.4.23: calculate composition, rate and'
     x' enthalpy for wells on deliverability*********')
c
      ic=mop(9)+1
      pin=pi(n)
C-----LOOP OVER NEQ1 SETS OF PARAMETERS.
      DO 101 M=1,NEQ1
      JLM2=J2LOC+(M-1)*NSEC
      if(eosn(1).eq.'EWASG     ') then
c.....   renormalize productivity index for permeability reduction
         pin=pi(n)*par(jlm2+2*nbk+3)
c        pin=pi(n)*par(jlm2+2*nbk+6)
      endif
C
C-----ZERO OUT BLOCK (J,J).
      DO 1011 K=1,NK1
 1011 D(K,M)=0.d0
C
      FFS=0.0D0
C-----FFS IS THE SUM OVER MOBILITY*DENSITY IN ALL PHASES.
C-----NOW COMPUTE FRACTIONAL FLOWS IN THE VARIOUS PHASES.
      g(n)=0.d0
C
      DO 30 NP=1,NPH
      J2LNP=JLM2+(NP-1)*NBK
c
      FF((N-1)*NPH+NP)=0.d0
c
      IF(IC.EQ.1.AND.PAR(J2LNP+3).NE.0.d0)
     AFF((N-1)*NPH+NP)=PAR(J2LNP+2)*PAR(J2LNP+4)/PAR(J2LNP+3)
c
      IF(IC.EQ.2) FF((N-1)*NPH+NP)=PAR(J2LNP+1)*PAR(J2LNP+4)
c
C
C=====COMPUTE FLOW RATES FOR EACH PHASE AND TOTAL RATE FOR WELL ON
C     DELIVERABILITY.===================================================
C
      DPRES=0.d0
      IF(M.EQ.2) DPRES=DELX(JLOCP+1)
      DELP=X(JLOCP+1)+DX(JLOCP+1)+DPRES+PAR(J2LNP+6)-PWB(N)
      IF(DELP.LE.0.d0) DELP=0.d0
C-----COMPUTE MASS FLOW RATE.
      FF((N-1)*NPH+NP)=FF((N-1)*NPH+NP)*pin*DELP
      G(N)=G(N)-FF((N-1)*NPH+NP)
      FFS=FFS+FF((N-1)*NPH+NP)
C
C=====END OF DELIVERABILITY SECTION.====================================
c
      IF(MOP(4).GE.5) WRITE (34,3) M,NP,PAR(J2LNP+2),PAR(J2LNP+3),
     APAR(J2LNP+4),delp,FF((N-1)*NPH+NP),g(n)
    3 FORMAT(3H M=,I2,4H NP=,I1,3H K=,E12.6,5H VIS=,E12.6,3H D=,E12.6,
     A' DELP=',E12.6,' FF=',E12.6,' GN=',E12.6)
C
   30 CONTINUE
C
      gn=g(n)
      EG(N)=0.d0
C
C-----RENORMALIZE FRACTIONAL FLOWS TO 1, AND COMPUTE CONTRIBUTIONS
C     OF COMPONENTS IN PHASES.
C
      DO 31 NP=1,NPH
      J2LNP=JLM2+(NP-1)*NBK
      IF(FFS.NE.0.d0)
     AFF((N-1)*NPH+NP)=FF((N-1)*NPH+NP)/FFS
      EG(N)=EG(N)+FF((N-1)*NPH+NP)*PAR(J2LNP+5)
C
      DO 32 K=1,NK
      D(K,M)=D(K,M)-FAC*GN*FF((N-1)*NPH+NP)*PAR(J2LNP+NB+K)
   32 CONTINUE
C
   31 CONTINUE
C
      D(NK1,M)=D(NK1,M)-FAC*GN*EG(N)
C
  101 CONTINUE
c
      return
      end
c
c
c
      subroutine phasf(n,jlocp,j2loc,fac)
c
c-----this subroutine calculates flowing bottomhole pressure,
c     phase composition, flow rate, and flowing enthalpy for
c     wells on deliverability with specified wellhead pressure
c
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P3/DELX((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/G4/ELEG(MNOGN)
      COMMON/G5/SOURCE(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G13/PI(MNOGN)
      COMMON/G14/PWB(MNOGN)
      COMMON/G17/SDENS(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/dkm/d(11,12)
      common/ech/eosn(20)
      character*10 eosn
c     common below was missing
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      character eleg*5,source*5
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ****PHASF 1.0, 2001.4.23: calculate bottomhole pressure,'
     x' fluid composition, rate and enthalpy for production with'
     x' specified wellhead pressure********')
c
c.....initialize well rate
      gw0=sdens(n)
c
      ic=mop(9)+1
      pin=pi(n)
C-----LOOP OVER NEQ1 SETS OF PARAMETERS.
      DO 101 M=1,NEQ1
      JLM2=J2LOC+(M-1)*NSEC
      if(eosn(1).eq.'EWASG     ') then
c.....   renormalize productivity index for permeability reduction
         pin=pi(n)*par(jlm2+2*nbk+3)
c        pin=pi(n)*par(jlm2+2*nbk+6)
      endif
C
C-----ZERO OUT BLOCK (J,J).
      DO 1011 K=1,NK1
 1011 D(K,M)=0.d0
C
      FFS=0.0D0
C-----FFS IS THE SUM OVER MOBILITY*DENSITY IN ALL PHASES.
C-----NOW COMPUTE FRACTIONAL FLOWS IN THE VARIOUS PHASES.
      eg(n)=0.d0
C
      DO 30 NP=1,NPH
      J2LNP=JLM2+(NP-1)*NBK
c
      FF((N-1)*NPH+NP)=0.d0
c
      IF(IC.EQ.1.AND.PAR(J2LNP+3).NE.0.d0)
     AFF((N-1)*NPH+NP)=PAR(J2LNP+2)*PAR(J2LNP+4)/PAR(J2LNP+3)
c
      IF(IC.EQ.2) FF((N-1)*NPH+NP)=PAR(J2LNP+1)*PAR(J2LNP+4)
c
      FFS=FFS+FF((N-1)*NPH+NP)
      EG(N)=EG(N)+FF((N-1)*NPH+NP)*PAR(J2LNP+5)
c
      IF(MOP(4).GE.5) WRITE (34,3) M,NP,PAR(J2LNP+2),PAR(J2LNP+3),
     APAR(J2LNP+4),FF((N-1)*NPH+NP),FFS,eg(n)
    3 FORMAT(' M=',I2,' NP=',I1,' kr=',E12.6,' VIS=',E12.6,' D=',E12.6,
     A' FF=',E12.6,' FFS=',E12.6,' egffs=',E12.6)
C
   30 CONTINUE
C
      if(ffs.ne.0.d0) EG(N)=eg(n)/ffs
c
c.....6-10-97: now iterate on flow rate to find bottomhole pressure
c
      pj=x(jlocp+1)+dx(jlocp+1)
      if(m.eq.2) pj=pj+delx(jlocp+1)
c
c     initialize well flow rate (positive for production)
      gw=gw0
c
c     Newton/Raphson iteration
      ig=0
   34 continue
      ig=ig+1
c
      call wf(n,gw,eg(n),pwbn,dpdg)
      if(igood.ne.0) goto 4002
c     reservoir rate (negative for production)
      gres=-pin*ffs*(pj-pwbn)
      rg=gres+gw
      if(mop(4).ge.5) WRITE (34,35) ig,pwbn,gw,gres,dpdg,rg
   35 format(' IG=',I2,'   Pwb=',E12.6,'   Gwell=',E12.6,
     x'   Gres=',E12.6,'   dPwb/dG =',E12.6,'   RG =',E12.6)
      if(abs(rg/gw).le.1.d-10) goto 33
c
      if(ig.gt.20) goto 4004
      drgdg=pin*ffs*dpdg+1.0d0
      gw=gw-rg/drgdg
      goto 34
c
   33 continue
         gn=gres
      if(m.eq.1) then
         pwb(n)=pwbn
         g(n)=gres
c.....store current well rate, to be used for next initialization
         sdens(n)=-gres
      endif
C
C-----RENORMALIZE FRACTIONAL FLOWS TO 1, AND COMPUTE CONTRIBUTIONS
C     OF COMPONENTS IN PHASES.
C
      DO 31 NP=1,NPH
      J2LNP=JLM2+(NP-1)*NBK
      IF(FFS.NE.0.d0)
     AFF((N-1)*NPH+NP)=FF((N-1)*NPH+NP)/FFS
C
      DO 32 K=1,NK
      D(K,M)=D(K,M)-FAC*GN*FF((N-1)*NPH+NP)*PAR(J2LNP+NB+K)
   32 CONTINUE
C
   31 CONTINUE
C
      D(NK1,M)=D(NK1,M)-FAC*GN*EG(N)
C
  101 CONTINUE
c
      return
c
 4002 WRITE (34,4003) kcyc,eleg(n),source(n)
 4003 format(' at KCYC = ',I4,' exceed *WFLO* table data at element ',
     xA5,' (source ',A5,') -- will reduce DELTEX')
      return
c
 4004 igood=3
      WRITE (34,4005) kcyc,eleg(n),source(n),gres,gn,pwbn
 4005 format(' at KCYC = ',I4,', element ',A5,' (source ',A5,') no',
     x' convergence for Pwb = ',E12.6,'   Gres = ',E12.6,'   Gwel = ',
     xE12.6)
      return
c
      end
c
c
c
      SUBROUTINE WF(n,g,h,p,dpdg)
C
C-----THIS ROUTINE PERFORMS BIVARIATE INTERPOLATION TO FIND
C     bottomhole pressure as function of flow rate and flowing
c     enthalpy.
c     (after routine KTP of EOSYYY)
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      common/g3a/pw(mgtab)
      common/g28/nftab(mnogn)
      common/g29/iftit(mnogn)
      common/g30/jftit(mnogn)
      common/g31/ijf(mnogn)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/KC/KC
      SAVE ICALL
c
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) then
      WRITE(11,899)
c 899 FORMAT(6X,'WF       1.0       3 July      1997',6X,
  899 FORMAT(6X,'WF       1.1      13 February  1998',6X,
     x'perform bivariate interpolation from flowing wellbore',
     x' pressure table')
      endif
c
      kf=nftab(n)
      kij=ijf(kf)
      ng=iftit(kf)
      nh=jftit(kf)
c
      IF(g.LT.pw(kij+1).OR.g.GT.pw(kij+ng)) GOTO 100
      IF(h.LT.pw(kij+ng+1).OR.h.GT.pw(kij+ng+nh)) GOTO 200
C
      KL=1
      KR=ng
    2 IF(KR-KL.EQ.1) GOTO 10
      KM=KL+(KR-KL)/2
      IF(g.LE.pw(kij+KM)) GOTO 1
      KL=KM
      GOTO 2
    1 KR=KM
      GOTO 2
   10 CONTINUE
C
C
      LL=1
      LR=nh
   12 IF(LR-LL.EQ.1) GOTO 20
      LM=LL+(LR-LL)/2
      IF(h.LE.pw(kij+ng+LM)) GOTO 11
      LL=LM
      GOTO 12
   11 LR=LM
      GOTO 12
   20 CONTINUE
C
C
      Q=(g-pw(kij+KL))/(pw(kij+KR)-pw(kij+KL))
      R=(h-pw(kij+ng+LL))/(pw(kij+ng+LR)-pw(kij+ng+LL))
C
      p=(1.d0-Q)*(1.d0-R)*pw(kij+ng+nh+(KL-1)*nh+LL)
     A+Q*(1.d0-R)*pw(kij+ng+nh+(KR-1)*nh+LL)
     B+(1.d0-Q)*R*pw(kij+ng+nh+(KL-1)*nh+LR)
     C+R*Q*pw(kij+ng+nh+(KR-1)*nh+LR)
C
C-----NOW COMPUTE DERIVATIVE OF pressure WITH RESPECT TO flow rate.
      DQDg=1.d0/(pw(kij+KR)-pw(kij+KL))
      DpDQ=(1.d0-R)*(pw(kij+ng+nh+(KR-1)*nh+LL)
     x-pw(kij+ng+nh+(KL-1)*nh+LL))
     A+R*(pw(kij+ng+nh+(KR-1)*nh+LR)-pw(kij+ng+nh+(KL-1)*nh+LR))
      DpDg=DpDQ*DQDg
C
c     WRITE (34,3) g,h,kl,ll,q,r,p,pw(kij+kl),pw(kij+ng+ll),
c    xpw(kij+ng+nh+(kl-1)*nh+ll)
    3 format(' g = ',E12.6,'   h = ',E12.6,'   kl =',I3,'   ll =',I3,
     x'   q = ',E12.6,'   r = ',E12.6,'   p = ',E12.6/
     x' pw(kij+kl) = ',E12.6,'   pw(kij+ng+ll) = ',E12.6,
     x'   pw(kij+ng+nh+(kl-1)*nh+ll) = ',E12.6)
c     WRITE (34,4) pw(kij+ng+nh+(KL-1)*nh+LL),pw(kij+ng+nh+(KR-1)*nh+LL),
c    xpw(kij+ng+nh+(KL-1)*nh+LR),pw(kij+ng+nh+(KR-1)*nh+LR)
    4 format(' pw(kl,ll) = ',E12.6,'   pw(kr,ll) = ',E12.6,
     x'   pw(kl,lr) = ',E12.6,'   pw(kr,lr) = ',E12.6)
      RETURN
C
  100 CONTINUE
      WRITE(6,101)g,pw(kij+1),pw(kij+ng)
  101 FORMAT(' &&& WF &&&   g = ',E12.6,' IS OUTSIDE THE RANGE (',
     xE12.6,',',E12.6,')   &&&&&&&&&&')
      IF(KC.EQ.0) STOP
      IGOOD=3
      RETURN
C
C
  200 CONTINUE
      WRITE(6,201)h,pw(kij+ng+1),pw(kij+ng+nh)
  201 FORMAT(' &&& WF &&&   h = ',E12.6,' IS OUTSIDE THE RANGE (',
     xE12.6,',',E12.6,')   &&&&&&&&&&')
      IF(KC.EQ.0) STOP
      IGOOD=3
      RETURN
C
      END
c
c
c
      subroutine ttab(n,ltaba,itaba,gx,egx)
c
c-----This subroutine performs interpolation of sink/source
c     rates and enthalpies from time-dependent tables.
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/ BCIJT /    IHD,    MAX,    MID,    SET
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/G1/F1(MGTAB)
      COMMON/G2/F2(MGTAB)
      COMMON/G3/F3(MGTAB)
      COMMON/G4/ELEG(MNOGN)
      COMMON/G5/SOURCE(MNOGN)
      COMMON/G11/NGIND(MNOGN)
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
c
      character itaba*1,eleg*5,source*5
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **TTAB 1.0, 1997.11.5: Interpolate sink/source rates and'
     x' enthalpies from tables**********')
c
      ihd=16
      MAX=LTABA
      seti=sumtim
      SET=SETI
      setf=sumtim+deltex
C
      NBG=NGIND(N)
      CALL FINDL(NBG)
      IF(IHD.EQ.0) GOTO4000
      if(mop(12).eq.2) then
      mid1=mid
      else
C
      IF(MOP(12).EQ.1) Q1=F2(NBG+MID)
      IF(MOP(12).EQ.0) CALL QINTER(Q1,NBG)
C-----Q1 IS THE GENERATION RATE AT BEGINNING OF CURRENT TIME STEP.
C
      IF(ITABA.NE.' '.AND.MOP(12).EQ.1) R1=F3(NBG+MID)
      IF(ITABA.NE.' '.AND.MOP(12).EQ.0) CALL HINTER(R1,NBG)
      endif
C
      SET=SETF
      MAX=LTABA
      CALL FINDL(NBG)
      IF(IHD.EQ.0) GOTO 4000
C
      if(mop(12).eq.2) then
      mid2=mid
      else
      IF(MOP(12).EQ.1) Q2=F2(NBG+MID)
      IF(MOP(12).EQ.0) CALL QINTER(Q2,NBG)
C-----Q2 IS THE GENERATION RATE AT END OF CURRENT TIME STEP.
C
      IF(ITABA.NE.' '.AND.MOP(12).EQ.1) R2=F3(NBG+MID)
      IF(ITABA.NE.' '.AND.MOP(12).EQ.0) CALL HINTER(R2,NBG)
C
      gx=(Q1+Q2)/2.d0
C-----GENERATION RATE IS ASSUMED TO BE THE AVERAGE OF THE VALUES AT
C     BEGINNING AND END OF TIME STEP, RESPECTIVELY.
C
      IF(ITABA.NE.' ') egx=(R1+R2)/2.d0
C-----FOR (ITABA.NE.' ') A VALUE FOR (PRODUCED OR INJECTED) ENTHALPY
C     IS ASSIGNED.
      endif
      if(mop(12).eq.2) then
         if(mid1.eq.mid2) then
         gx=f2(nbg+mid1)
         if(itaba.ne.' ') egx=f3(nbg+mid1)
         else
         qdt=f2(nbg+mid1)*(f1(nbg+mid1+1)-seti)
     x   +f2(nbg+mid2)*(setf-f1(nbg+mid2))
         if(itaba.ne.' ')
     x   hqdt=f3(nbg+mid1)*f2(nbg+mid1)*(f1(nbg+mid1+1)-seti)
     x   +f3(nbg+mid2)*f2(nbg+mid2)*(setf-f1(nbg+mid2))
            if(mid2.gt.mid1+1) then
            mid11=mid1+1
            mid21=mid2-1
            do3902 m=mid11,mid21
            qdt=qdt+f2(nbg+m)*(f1(nbg+m+1)-f1(nbg+m))
            if(itaba.ne.' ')
     x      hqdt=hqdt+f3(nbg+m)*f2(nbg+m)*(f1(nbg+m+1)-f1(nbg+m))
 3902       continue
            endif
         gx=qdt/deltex
         egx=0.d0
         if(itaba.ne.' '.and.qdt.ne.0.d0) egx=hqdt/qdt
         endif
      endif
c
      return
C
 4000 IGOOD=3
      WRITE (34,4001) KCYC,ELEG(N),SOURCE(N)
 4001 FORMAT(11H AT KCYC = ,I4,41H EXCEED GENERATION TIME TABLE AT ELEME
     ANT ,A5,9H (SOURCE ,A5,22H) -- WILL REDUCE DELTE)
      RETURN
C
      end
c
c
c
      subroutine gcor(n,ltaba)
c
c-----this subroutine calculates a very simple gravity correction
c     to flowing bottomhole pressure of multi-feedzone wells
c     on deliverability.
c
C=====COME HERE FOR WELL ON DELIVERABILITY, AND COMPUTE GRAVITY
C     CORRECTION TO FLOWING BOTTOMHOLE PRESSURE.========================
C     CODING CAN HANDLE UP TO THREE PHASES FLOWING SIMULTANEOUSLY
C
C     COMPUTE THIS CORRECTION ONCE AND FOR ALL, AT BEGINNING OF EACH
C     TIME STEP, USING NON-INCREMENTED VARIABLES.
C
C     THE GRAVITY-CORRECTIONS FOR FLOWING WELLBORE PRESSURE ARE
C     COMPUTED FOR ALL LAYERS IN WHICH A WELL ON DELIVERABILITY
C     IS OPEN. THE COMPUTATION IS MADE WHEN THE SOURCE CORRESPONDING
C     TO THE BOTTOM LAYER IS ENCOUNTERED. THIS IS SINGLED OUT BY
C     HAVING "TYPE = DELV" AND "LTAB = NUMBER OF LAYERS".
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
c
      include 'flowpar_v2.inc'
c
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G13/PI(MNOGN)
      COMMON/G14/PWB(MNOGN)
      COMMON/G15/HG(MNOGN)
      COMMON/G22/QVGC(MNOGN)
      COMMON/G23/QVWC(MNOGN)
      COMMON/G24/QVOC(MNOGN)
      COMMON/G25/GRAD(MNOGN)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/ech/eosn(20)
      character*10 eosn
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' GCOR 1.0, 2001.4.23: perform simple gravity correction'
     x' for flowing bottomhole pressure')
c
      NLTAB1=N+LTABA-1
C
C
C-----LOOP OVER ALL LAYERS FOR THIS WELL.
      DO 50 ND=N,NLTAB1
C     ND IS THE GENERATION INDEX.
      JD=NEXG(ND)
c not used      JDLOC=(JD-1)*NEQ
      JDLOCP=(JD-1)*NK1
      JDLOC2=(JD-1)*NSEC*NEQ1
      pind=pi(nd)
      if(eosn(1).eq.'EWASG     ') then
c.....   renormalize productivity index for permeability reduction
         pind=pi(nd)*par(jdloc2+2*nbk+3)
      endif
C     COMPUTE PHASE MOBILITIES.
      AMG=0.d0
      IF(PAR(JDLOC2+3).NE.0.d0) AMG=PAR(JDLOC2+2)/PAR(JDLOC2+3)
      AMW=0.d0
      IF(NPH.GE.2.AND.
     XPAR(JDLOC2+3+NBK).NE.0.d0) AMW=PAR(JDLOC2+2+NBK)/PAR(JDLOC2+3+NBK)
      AMO=0.d0
      IF(NPH.GE.3.AND.
     XPAR(JDLOC2+3+2*NBK).NE.0.d0) AMO=PAR(JDLOC2+2+2*NBK)/PAR(JDLOC2
     A+3+2*NBK)
C
C*****COMPUTE VOLUMETRIC PRODUCTION RATES FOR INDIVIDUAL PHASES.
      QVG=pind  *AMG*X(JDLOCP+1)
      QVW=pind  *AMW*X(JDLOCP+1)
      QVO=pind  *AMO*X(JDLOCP+1)
C*****ARRAYS FOR CUMULATIVE PHASE RATES.
      IF(ND.EQ.N) GOTO 52
      QVGC(ND)=QVGC(ND-1)+QVG
      QVWC(ND)=QVWC(ND-1)+QVW
      QVOC(ND)=QVOC(ND-1)+QVO
      GOTO 53
   52 QVGC(ND)=QVG
      QVWC(ND)=QVW
      QVOC(ND)=QVO
   53 CONTINUE
C*****ARRAY OF GRADIENTS.
      GRAD(ND)=(QVGC(ND)*PAR(JDLOC2+4)+QVWC(ND)*PAR(JDLOC2+NBK+4)
     A+QVOC(ND)*PAR(JDLOC2+2*NBK+4))*9.80665d0
     B/(QVGC(ND)+QVWC(ND)+QVOC(ND))
   50 CONTINUE
C
C-----NOW PERFORM LOOP GOING FROM TOP LAYER DOWN.
C
      N1=N+1
      DO 51 ND=N1,NLTAB1
      NDD=NLTAB1+N-ND
C     NDD IS THE RUNNING SOURCE INDEX, GOING FROM THE LAYER JUST BELOW
C     THE TOP (INDEX N+LTABA-2) DOWN TO THE BOTTOM LAYER (INDEX N).
C
      PWB(NDD)=PWB(NDD+1)+(HG(NDD+1)*GRAD(NDD+1)+HG(NDD)*GRAD(NDD))/2.d0
   51 CONTINUE
      return
      end
c
      SUBROUTINE CONVER
C
C-----THIS SUBROUTINE IS CALLED AFTER SUCCESFULL COMPLETION OF
C     A TIME STEP.
C     IT UPDATES PRIMARY VARIABLES, AND DEFINES THE NEXT TIME STEP.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
c
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
C
      COMMON/G7/G(MNOGN)
      COMMON/G12/LCOM(MNOGN)
      COMMON/G15/HG(MNOGN)
C
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +             GK(MAXMAT)
      COMMON/KC/KC
      COMMON/DFM/TIMAX,REDLT
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DLT/NDLT,DLT(100)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      COMMON/BC/NELA
C
      DIMENSION DXM(10)
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *****CONVER 1.0, 1991.3.4: UPDATE PRIMARY VARIABLES'
     X' AFTER CONVERGENCE IS ACHIEVED*****')
C
      DO 10 K=1,NK1
   10 DXM(K)=0.d0
C
      DO 3 N=1,NELA
      NLOC=(N-1)*NK1
      NLOC2=(N-1)*NSEC*NEQ1
C-----COMPUTE CHANGES IN POROSITY.
      PHIN=PHI(N)
      NMAT=MATX(N)
      DPHI=PHIN*(COM(NMAT)*DX(NLOC+1)+EXPAN(NMAT)*(PAR(NLOC2+NSEC-1)
     A-T(N)))
      PHI(N)=PHIN+DPHI
C
C-----UPDATE ELEMENT PRESSURES AND TEMPERATURES.
      P(N)=X(NLOC+1)+DX(NLOC+1)
      T(N)=PAR(NLOC2+NSEC-1)
C
C-----INCREMENT PRIMARY VARIABLES.
      DO3 M=1,NEQ
      NLM=NLOC+M
      X(NLM)=X(NLM)+DX(NLM)
      DXM(M)=MAX(DXM(M),ABS(DX(NLM)))
    3 CONTINUE
C
C-----FOR PERCENTAGE INJECTION, ASSIGN INJECTION RATE FOR NEXT
C     TIME STEP.
C%    DO 30 N=1,NOGN
C%    IF(LCOM(N).NE.NEQ+4) GOTO 30
C%    G(N)=-HG(N)*G(N-1)
C% 30 CONTINUE
C
      SUMTIM=SUMTIM+DELTEX
      IF(TIMAX.NE.0.d0.AND.TIMAX.EQ.SUMTIM) NOWTIM=1
C-----AFTER CONVERGENCE UPDATE TOTAL TIME AND ASSIGN NEW TIME STEP.
      IF(NDLT.EQ.0) GOTO20
c Add line below so that DLT(KC+1) does not exceed its bounds
        if(kc.ge.99)go to 20
      IF(KC+1.GT.8*NDLT.OR.DLT(KC+1).EQ.0.d0) GOTO20
C-----IF NO FURTHER TIME STEP INSTRUCTIONS ARE PROVIDED, KEEP
C     GOING WITH LAST TIME STEP.
C-----COME HERE FOR NEW TIME STEP ASSIGNMENT.
      DELT=DLT(KC+1)
      GOTO 22
   20 DELT=DELTEX
      IF(ITER.LE.MOP(16)) DELT=2.d0*DELTEX
   22 IF(TIMAX.NE.0.d0) DELT=MIN(DELT,TIMAX-SUMTIM)
      IF(DELTMX.NE.0.d0) DELT=MIN(DELT,DELTMX)
cels9/18/06 add lower limit
      delt = max(delt,0.01d0)
C
      RETURN
      END
c
      SUBROUTINE WRIFI
C
C-----THIS SUBROUTINE IS CALLED AT THE COMPLETION OF A TOUGH2-RUN.
C     IT GENERATES A FILE *SAVE* TO BE USED FOR RESTARTING.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
c
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +             SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/SVZ/NOITE,MOP(24)
c      COMMON/BALA/TFM0(27),TVM0(27),TFE0(27),TVE0(27),
c     +               TSE0(27)
      COMMON/POV6/TSTART
c
C--------------------------------------------------------For using EOS9
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
      COMMON/TEM_EOS9/Tc_EOS9(MNEL)  ! initial temperature (oC)
c---------------------------------------------------------------------
c
      CHARACTER*5 ELEM,MAT
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***WRIFI 1.0, 1990.1.22: WRITE PRIMARY VARIABLES ON FILE'
     X' *SAVE* AT THE COMPLETION OF A TOUGH2 RUN***')
C
      kcycloc = kcyc
C
      OPEN(UNIT=111,FILE='STIME',STATUS='UNKNOWN')
      REWIND 7
      REWIND 111
      if (ntsave.gt.1) WRITE (34,5) kcycloc,SUMTIM
    5 FORMAT(/' WRITE FILE *SAVE* AFTER',I7,' TIME STEPS'/)
      WRITE(7,1) NEL,SUMTIM
    1 FORMAT(31HINCON -- INITIAL CONDITIONS FOR,I5,17H ELEMENTS AT TIME,
     AE14.6)
      DO 10 N=1,NEL
      NLOC=(N-1)*NK1
      IF(MOP(15).GE.1) WRITE(111,1461) STIME(N),MLAGNR(N)
 1461 FORMAT(2X,E12.6,2X,I1)
C
C-----------------------------------------------------------For using EOS9
      IF (IEOS.EQ.9) THEN
         X(NLOC+2)=Tc_EOS9(N) ! initial temperature distribution (oC)
      END IF
c
C-------------------------------------------------------------------------
c            Modified to write out permeabilities
c               WRITE(7,2) ELEM(N),PHI(N),(X(NLOC+I),I=1,NK1)
c               2 FORMAT(A5,10X,E15.8/(4E20.13))
      IF(X(NLOC+1).LE.-2.7315D2) THEN
         WRITE(*,*) 'parameter failure'
         STOP
      END IF
      WRITE(7,2) ELEM(N),PHI(N),perm(1,n),perm(2,n),perm(3,n),
     +   (X(NLOC+I),I=1,NK1)
    2 FORMAT(A5,10X,4E15.8/(4E20.13))
C-------------------------------------------------------------------------
C
   10 CONTINUE
C-----write restart information
      WRITE(7,3)
    3 FORMAT(5H+++  )
      WRITE(7,4) KCYCLOC,ITERC,NM,TSTART,SUMTIM
    4 FORMAT(2I10,I5,2E15.8)
C
      ENDFILE 7
      CLOSE(111)

      RETURN
      END
C
      SUBROUTINE FINDER(VECT,N)
C
C-----FINDER LOCATES GIVEN VALUE AMONG ENTRIES OF GIVEN VECTOR.         FINDER.4
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
C                                                                       FINDER.5
      COMMON/ BCIJT /    IHD,    MAX,    MID,    SET
      COMMON/  CDT7 /     HC
C                                                                       FINDER.9
      DIMENSION VECT(N)
C                                                                       FINDER.1
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ******FINDER 1.0, 1990.1.22: INTERPOLATE FROM A TABLE OF'
     x' TIME-DEPENDENT DATA******')
C
      MIN = 1
      IF(SET    .LT. VECT(MIN) .OR. SET .GT. VECT(MAX)) GO TO 5
 1    MID = (MIN + MAX)/2
      IF(SET - VECT(MID)) 2, 4, 3
 2    MAX = MID
      IF(MAX - 2) 4, 1, 1
 3    MIN = MID
      IF(MAX - MIN .GE. 2) GO TO 1
 4    RETURN
    5 CONTINUE
      WRITE (34,5000) SET,VECT(MIN),VECT(MAX)
      IHD    = 0
      RETURN
C                                                                       FINDER.2
 5000 FORMAT(1X,E12.6,24H IS OUTSIDE THE RANGE  (,E11.4,1H,,E11.4,1H))
C                                                                       FINDER.3
      END
C
C
C
      SUBROUTINE FINDL(NB)
C                                                                       FINDER.3
C-----FINDER LOCATES GIVEN VALUE AMONG ENTRIES OF GIVEN VECTOR.         FINDER.4
C                                                                       FINDER.5
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      COMMON/G1/F1(MGTAB)
      COMMON/ BCIJT /    IHD,    MAX,    MID,    SET
      COMMON/  CDT7 /     HC
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(6X,'FINDL 1.0, 1990.1.22: INTERPOLATE FROM A TABLE OF'
     x' TIME-DEPENDENT DATA')
C                                                                       FINDER.1
      MIN = 1
      IF(SET    .LT. F1(NB+MIN) .OR. SET .GT. F1(NB+MAX)) GO TO 5
 1    MID = (MIN + MAX)/2
      IF(SET - F1(NB+MID)) 2, 4, 3
 2    MAX = MID
      IF(MAX - 2) 4, 1, 1
 3    MIN = MID
      IF(MAX - MIN .GE. 2) GO TO 1
 4    RETURN
    5 CONTINUE
      WRITE (34,5000) SET,F1(NB+MIN),F1(NB+MAX)
 5000 FORMAT(1X,E12.6,24H IS OUTSIDE THE RANGE  (,E11.4,1H,,E11.4,1H))
      IHD    = 0
      RETURN
      END
C
      SUBROUTINE QINTER(Q,NBG)
C
C*****THIS ROUTINE USES LINEAR INTERPOLATION TO COMPUTE GENERATION
C     RATE AT TIME T=SET.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      COMMON/G1/F1(MGTAB)
      COMMON/G2/F2(MGTAB)
      COMMON/BCIJT/IHD,MAX,MID,SET
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' QINTER 1.0, 1990.1.22: PERFORM LINEAR INTERPOLATION')
      Q=F2(NBG+MID)+(SET-F1(NBG+MID))*(F2(NBG+MID+1)-F2(NBG+MID))
     A                               /(F1(NBG+MID+1)-F1(NBG+MID))
C
      RETURN
      END
c
      SUBROUTINE HINTER(Q,NBG)
C
C*****THIS ROUTINE USES LINEAR INTERPOLATION TO COMPUTE FLOWING
C     ENTHALPY AT TIME T=SET.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      COMMON/G1/F1(MGTAB)
      COMMON/G3/F3(MGTAB)
      COMMON/BCIJT/IHD,MAX,MID,SET
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' HINTER 1.0, 1990.1.22: PERFORM LINEAR INTERPOLATION')
      Q=F3(NBG+MID)+(SET-F1(NBG+MID))*(F3(NBG+MID+1)-F3(NBG+MID))
     A                               /(F1(NBG+MID+1)-F1(NBG+MID))
      RETURN
      END
c
      SUBROUTINE THYME(IFLAG,DT,TOTAL)
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      real*4 tsec,t1
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' THYME 1.0, 1990.1.22: CALCULATE ELAPSED CPU TIME')
      t1 = 0.0
      call cpu_time(tsec)
      IF(IFLAG.NE.0) GOTO 1
      T1=Tsec
      RETURN
    1 DT = dble(Tsec-T1)
      TOTAL = dble(Tsec)
      RETURN
      END
c
      SUBROUTINE TSTEPT
C
C-----THIS ROUTINE MODIFIES TIME STEPS TO COINCIDE WITH USER-
C     INPUT VALUES AT WHICH PRINTOUT IS DESIRED.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/KONIT/KON,DELT,IGOOD
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' TSTEP 1.0, 1991.3.4: ADJUST TIME STEPS TO COINCIDE WITH'
     x' USER-DEFINED TARGET TIMES')
C
      IF(TIS(ITI).LE.SUMTIM) GOTO 100
C
C-----FIND TIME FOLLOWING SUMTIM.
      DO 1 I=1,ITI
c after ysw      IF(TIS(I).EQ.SUMTIM) GOTO 4
      IF(abs(TIS(I)-SUMTIM).le.0.d0) GOTO 4
c...
      IF(TIS(I).LT.SUMTIM) GOTO 1
      GOTO 2
    1 CONTINUE
C
    2 IF(SUMTIM+DELT.LT.TIS(I)) GOTO 10
C
C-----COME HERE TO ADJUST DELT.
      DELT=TIS(I)-SUMTIM
c     add lower limit
      delt = max(delt,0.01d0)
      NOWTIM=1
      GOTO 10
C
C-----COME HERE AFTER NEXT TIME HAS BEEN REACHED.
    4 ITCO=I
C
      IF(DELAF.GT.0.d0)
     ADELT=MIN(DELT,DELAF)
C
      IF(TIS(I+1).GT.TIS(I)+DELT) RETURN
      DELT=MIN(DELT,TIS(I+1)-TIS(I))
c     add lower limit
      delt = max(delt,0.01d0)
      NOWTIM=1
C
      RETURN
   10 ITCO=I-1
      RETURN
  100 ITCO=ITI
      RETURN
      END
C
      SUBROUTINE RELP(SG,REPL,REPG,NMAT,K,NLOC,SG0)
C
C************ Modified for active fracture model: Tianfu Xu, 11/15/2001
C
C-----THIS ROUTINE COMPUTES RELATIVE PERMEABILITIES FOR LIQUID
C     AND GASEOUS PHASES.
C
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/P3/DELX((MNK+1)*MNEL)
      COMMON/RPCAP/IRP(MAXMAT),RP(7,MAXMAT),ICP(MAXMAT),CP(7,MAXMAT),
     XIRPD,RPD(7),ICPD,CPD(7)
c       Added definitions below
      double precision sg,repl,repg,sg0
      integer*8 nmat,k,nloc
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' RELP 1.1, 2009.3.20: LIQUID AND GAS PHASE RELATIVE'
     x' PERMEABILITIES AS FUNCTIONS OF SATURATION'/' for IRP=7, use'
     x' Corey-krg when RP(4).ne.0, with Sgr = RP(4)')
C
      SL=1.d0-SG
      GOTO(10,11,12,12,13,14,15,16,17,18,19),IRP(NMAT)
   10 CONTINUE
C-----LINEAR FUNCTIONS.
C
C     CHECK IF INCREMENT NEEDS TO BE ADJUSTED AT LOWER LIQUID CUTOFF.
      IF(K.NE.3) GOTO 20
      IF((SL-RP(1,NMAT))*(1.d0-SG0-RP(1,NMAT)).GE.0.d0) GOTO 20
C     ADJUST INCREMENT.
      DELX(NLOC+2)=-DELX(NLOC+2)
      SG=SG0+DELX(NLOC+2)
      SL=1.d0-SG
   20 CONTINUE
C
      REPL=(SL-RP(1,NMAT))/(RP(3,NMAT)-RP(1,NMAT))
      IF(SL.GE.RP(3,NMAT)) REPL=1.d0
      IF(SL.LE.RP(1,NMAT)) REPL=0.d0
      REPG=(SG-RP(2,NMAT))/(RP(4,NMAT)-RP(2,NMAT))
      IF(SG.GE.RP(4,NMAT)) REPG=1.d0
      IF(SG.LE.RP(2,NMAT)) REPG=0.d0
C
      RETURN
C
   11 CONTINUE
C-----RELATIVE PERMEABILITY OF PICKENS ET AL.
C
      REPG=1.d0
      REPL=(1.d0-SG)**RP(1,NMAT)
C
      RETURN
C
   12 CONTINUE
C-----COREY@S OR GRANT@S CURVES.
C
      SSTAR=(SL-RP(1,NMAT))/(1.d0-RP(1,NMAT)-RP(2,NMAT))
      REPL=SSTAR**4
      REPG=(1.d0-SSTAR**2)*(1.d0-SSTAR)**2
      IF(SG.GE.RP(2,NMAT)) GOTO 50
      REPG=0.d0
      REPL=1.d0
      GOTO 102
   50 IF(SG.LT.(1.d0-RP(1,NMAT))) GOTO 102
      REPL=0.d0
      REPG=1.d0
  102 CONTINUE
      IF(IRP(NMAT).EQ.4) REPG=1.d0-REPL
      RETURN
C
   13 CONTINUE
C-----BOTH PHASES ARE PERFECTLY MOBILE.
C
      REPL=1.d0
      REPG=1.d0
C
      RETURN
   14 CONTINUE
C-----RELATIVE PERMEABILITIES OF FATT AND KLIKOFF (1959), AS REPORTED
C     BY K. UDELL (BERKELEY, 1982).
C
      SS=0.d0
      IF(SL.GT.RP(1,NMAT)) SS=(SL-RP(1,NMAT))/(1.d0-RP(1,NMAT))
      REPL=SS**3
      REPG=(1.d0-SS)**3
      RETURN
C
   15 CONTINUE
C-----RELATIVE PERMEABILITY OF VAN GENUCHTEN, SOIL SCI. SOC. AM. J. 44,
C     PP. 892-898, 1980.
C
      IF(SL.GE.RP(3,NMAT)) GOTO 150
      SS=(SL-RP(2,NMAT))/(RP(3,NMAT)-RP(2,NMAT))
      REPL=0.d0
C
c----------------------------Active fracture model (Liu et al., 1998)
c----------------------------added by Tianfu Xu on 16-August-1999
c      IF(SS.GT.0.D0)
c     XREPL=SQRT(SS)*(1.D0-(1.D0-SS**(1.D0/RP(1,NMAT)))**RP(1,NMAT))**2
      IF(SS.GT.0.D0) THEN
c     Added from yu-shu eos3v1.4
        IF(cp(6,nmat).le.0.0d0) then
           REPL=SQRT(SS)*(1.d0-(1.d0-SS**(1.d0/RP(1,NMAT)))**
     +        RP(1,NMAT))**2
        else
c This is written differently than in eos3v1.4, but is equivalent
           gm=cp(6,nmat)  !Gam factor in active fracture model
           REPL=(SS**((1.0d0+gm)*0.5D0))*
     +      (1.d0-(1.d0-SS**((1.d0-gm)/RP(1,NMAT)))**RP(1,NMAT))**2
        endif
      END IF
c-------------------------------------------------------------------
c
c     11-23-94: for RP(4).ne.0, take Sgr=RP(4) and use Corey krg.
      if(rp(4,nmat).le.0.d0) then
           REPG=1.d0-REPL
         else
c.....7-26-95
              if(1.d0-sl.le.rp(4,nmat)) then
                 repg=0.d0
              else
                 SSTAR=(SL-RP(2,NMAT))/(1.d0-RP(2,NMAT)-RP(4,NMAT))
                 sstar=max(0.d0,sstar)
                 sstar=min(1.d0,sstar)
                 REPG=(1.d0-SSTAR**2)*(1.d0-SSTAR)**2
              endif
      endif
c-------------------------------------------------------------------
c Yu-shu's corey function
      if(rp(4,nmat).le.0.0d0.and.rp(5,nmat).le.0.0d0) then
         REPG=1.d0-REPL
      elseif(rp(4,nmat).gt.0.0d0.and.rp(5,nmat).le.0.0d0) then
        swbar=(sl-rp(2,nmat))/(1.d0-rp(2,nmat)-rp(4,nmat))
        repg=0.d0
        if(sg.gt.rp(4,nmat)) then
          swbar=min(1.d0,swbar)
          swbar=max(0.d0,swbar)
          ss1=1.d0-swbar
          repg=ss1*ss1*(1.d0-swbar**2)
        endif
      elseif(rp(4,nmat).le.0.0d0.and.rp(5,nmat).gt.0.0d0) then
c
C USING modified COREY-TYPE FUNCTION
        IF (SL.LE.RP(2,NMAT)) THEN
          REPG=1.0D0
        ELSE IF (SL.GE.RP(3,NMAT)) THEN
          REPG=0.0D0
        ELSE
          XM=RP(1,NMAT)/(1.0D0-RP(1,NMAT))
          REPG=(1.0D0-SS)*(1.0D0-SS)*(1.0D0-SS**((2.0D0+XM)/XM))
        ENDIF
c
      else
         stop 'check rp(4,nmat) and rp(5,nmat) values '
      endif
      RETURN
C
  150 REPL=1.d0
      REPG=0.d0
      RETURN
C
   16 CONTINUE
C     RELATIVE PERMEABILITIES AS MEASURED BY VERMA ET AL. IN
C     LABORATORY FLOW EXPERIMENTS FOR STEAM-WATER MIXTURES
C
      SS=(SL-RP(1,NMAT))/(RP(2,NMAT)-RP(1,NMAT))
      IF(SS.GT.1.d0) SS=1.d0
      IF(SS.LT.0.d0) SS=0.d0
      REPL=SS**3
      REPG=RP(3,NMAT)+RP(4,NMAT)*SS+RP(5,NMAT)*SS*SS
      IF(REPG.GT.1.d0) REPG=1.d0
      IF(REPG.LT.0.d0) REPG=0.d0
      RETURN
C
c     Options 9 and 10 not used now
   17 CONTINUE
      RETURN
c
   18 CONTINUE
      RETURN
c
c Directly from S. Finsterle (implemented in ITOUGH)
c             IRP = 11
c Note!!! MUST BE USED IN CONJUNCTION WITH ICP=11 !!!!
c
   19 CONTINUE
C-----Van Genuchten/Mualem (only together with ICP=11)
C     RP(1)=Residual liquid saturation
C     RP(2)=Residual gas saturation
C     RP(3).gt.0 krg=1.-krl
C     RP(4)=eta (exponent of Se in rel. liq. perm.)
C     RP(5)=Epsilon, linear interpolation between krl(1-eps) and 1
C     RP(7)=zeta (exponent of Se in rel. gas. perm.)
C     CP(1)=n
C     CP(4)=m (if 0, m=1-1/n)
C     CP(6)=gamma (Active Fracture Model)
C
C      
      IF (CP(4,NMAT).GT.1.0D-20) THEN
         XM=CP(4,NMAT)
      ELSE
         XM=1.0D0-1.0D0/CP(1,NMAT)
      ENDIF
      IF (RP(4,NMAT).EQ.0.d0) THEN
         ETA=0.5D0
      ELSE
         ETA=RP(4,NMAT)
      ENDIF
      SLRL=ABS(RP(1,NMAT))
      SGRL=MAX(RP(2,NMAT),0.0D0)
      SEL=(SL-SLRL)/(1.0D0-SLRL-SGRL)
      IF (SEL.GE.1.0D0) THEN
         REPL=1.0D0
      ELSE IF (SEL.LE.0.0D0) THEN
         REPL=0.0D0
      ELSE IF (SEL.LE.1.0D0-RP(5,NMAT)) THEN
         SSS=SEL**CP(6,NMAT)
         SS=SEL**(1.0D0-CP(6,NMAT))
c         REPL=SSS*SS**ETA*(1.0D0-(1.0D0-SS**(1.0D0/XM))**XM)**2.0D0
         REPL=SSS*SS**ETA*(1.0D0-(1.0D0-SS**(1.0D0/XM))**XM)**2
      ELSE
         SSS=(1.0D0-RP(5,NMAT))**CP(6,NMAT)
         C1=(1.0D0-RP(5,NMAT))**(1.0D0-CP(6,NMAT))
c         C2=SSS*C1**ETA*(1.0D0-(1.0D0-C1**(1.0D0/XM))**XM)**2.0D0
         C2=SSS*C1**ETA*(1.0D0-(1.0D0-C1**(1.0D0/XM))**XM)**2
         REPL=C2+(SEL-C1)*(1.0D0-C2)/RP(5,NMAT)
      ENDIF
      IF (RP(3,NMAT).GT.1.0D-10) THEN
         REPG=1.0D0-REPL
      ELSE
         IF (RP(7,NMAT).EQ.0.d0) THEN
            GAMMA=1.0D0/3.0D0
         ELSE
            GAMMA=RP(7,NMAT)
         ENDIF
         SLRG=MAX(RP(1,NMAT),0.0D0)
         SGRG=ABS(RP(2,NMAT))
         SEG=(SL-SLRG)/(1.0D0-SLRG-SGRG)
         IF (SEG.GT.1.0D0) THEN
            REPG=0.0D0
         ELSE IF (SEG.LT.0.0D0) THEN
            REPG=1.0D0
         ELSE
            REPG=(1.0D0-SEG)**GAMMA*
     &           (1.0D0-SEG**(1.0D0/XM))**(2.0D0*XM)
         ENDIF
      ENDIF
c  not used    SEFFL=SEL
c  not used      SEFFG=1.0D0-SEG
      RETURN
c
      END
c
c
c
      SUBROUTINE PCAP(SL,T,PC,NMAT)
C
c     Added options from yu-shu wu and stefan finsterle (itough2)
c
C************ Modified for active fracture model
C
C
C-----THIS ROUTINE COMPUTES CAPILLARY PRESSURE AS FUNCTION OF LIQUID
C     SATURATION SL AND TEMPERATURE T.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/RPCAP/IRP(MAXMAT),RP(7,MAXMAT),ICP(MAXMAT),CP(7,MAXMAT),
     AIRPD,RPD(7),ICPD,CPD(7)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +              CWET(MAXMAT),SH(MAXMAT)
C
      double precision sl,t,pc
      integer*8 nmat
c
      double precision pc_e(maxmat),pc_slop(maxmat)
      save pc_e,pc_slop
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' PCAP 1.1, 2009.4.14: CAPILLARY PRESSURE AS FUNCTION OF'
     X' SATURATION')
C
c--------------------------------------------
      IF(ICALL.EQ.1) then
c
c for linear extrapolation of Pc, ysw 2/12/97
c
         do i=1,nm
           if(icp(i).eq.10) then
             if(cp(6,i).le.0.0d0) then
                epsl=cp(4,i)
                if(epsl.le.0.d0) epsl=1.0d-5
                sbar=epsl/(cp(5,i)-cp(2,i))
                pc_e(i)=1.d0/cp(3,i)*(sbar
     B               **(-1.d0/cp(1,i))-1.d0)**(1.d0-cp(1,i))
                envg=1.0d0/(1.0d0-cp(1,i))
                pc_slop(i)=-1.0d0/(cp(3,i)*cp(1,i)*envg)
     A                   /(cp(5,i)-cp(2,i))*(sbar**(-1.d0/cp(1,i))
     B                   -1.0d0)**((1.0d0-envg)/envg)*sbar
     C                   **(-(1.0d0+cp(1,i))/cp(1,i))
             else
c  active fracture model
                epsl=cp(4,i)
                gama=cp(6,i)
                if(epsl.le.0.d0) epsl=1.0d-5
                sbar=epsl/(cp(5,i)-cp(2,i))
                pc_e(i)=1.d0/cp(3,i)*(sbar
     B               **((gama-1.0d0)/cp(1,i))-1.d0)**(1.d0-cp(1,i))
                pc_slop(i)=1.0d0/(cp(3,i)*cp(1,i))*(1.0d0-cp(1,i))
     1               *(gama-1.0d0)
     A               /(cp(5,i)-cp(2,i))*(sbar**((gama-1.d0)/cp(1,i))
     B               -1.0d0)**(-cp(1,i))*sbar
     C               **((gama-1.0d0-cp(1,i))/cp(1,i))
             endif
           endif
         enddo
      endif
c
c--------------------------------------------
c
c      GOTO(10,11,12,13,14,15,16,17),ICP(NMAT)
      GOTO(10,11,12,13,14,15,16,17,18,19,20),ICP(NMAT)
C
   10 CONTINUE
C-----LINEAR FUNCTION.
      PC=-CP(1,NMAT)*(CP(3,NMAT)-SL)/(CP(3,NMAT)-CP(2,NMAT))
      IF(SL.GE.CP(3,NMAT)) PC=0.d0
      IF(SL.LE.CP(2,NMAT)) PC=-CP(1,NMAT)
      RETURN
   11 CONTINUE
C-----CAPILLARY PRESSURE FUNCTION OF PICKENS ET AL, AS GIVEN IN
C     J. HYDROLOGY 40, 243-264, 1979.
C
      SLX=MAX(SL,1.001d0*CP(2,NMAT))
      IF(SLX.GT..999d0*CP(3,NMAT)) SLX=.999d0*CP(3,NMAT)
      A=(1.d0+SLX/CP(3,NMAT))*(CP(3,NMAT)-CP(2,NMAT))/
     A(CP(3,NMAT)+CP(2,NMAT))
      B=(1.d0-SLX/CP(3,NMAT))
      PC=-CP(1,NMAT)*DLOG(A*(1.d0+DSQRT(1.d0-B*B/(A*A)))/B)**
     A(1.d0/CP(4,NMAT))
      IF(SL.GT..999d0*CP(3,NMAT)) PC=PC*(1.d0-SL)/.001d0
      RETURN
C
C
   12 CONTINUE
C-----CAPILLARY PRESSURE FUNCTION AS USED IN THE TRUST-PROGRAM, WHICH
C     WAS DEVELOPED BY T.N. NARASIMHAN AT LAWRENCE BERKELEY LABORATORY.
C
      IF(SL.NE.1.d0) GOTO 120
      PC=0.d0
      RETURN
C
  120 SLX=SL
      IF(CP(5,NMAT).EQ.0.d0)SLX=MAX(SL,1.001d0*CP(2,NMAT))
      PC=-ABS(CP(5,NMAT))
      IF(SLX.GT.CP(2,NMAT))
     APC=-CP(4,NMAT)-CP(1,NMAT)*((1.d0-SLX)/(SLX-CP(2,NMAT)))
     B**(1.d0/CP(3,NMAT))
      IF(CP(5,NMAT).NE.0.d0)PC=MAX(PC,-ABS(CP(5,NMAT)))
      IF(SL.GT..999d0) PC=PC*(1.d0-SL)/.001d0
      RETURN
C
   13 CONTINUE
C-----CAPILLARY PRESSURE OF YOLO CLAY AFTER CHRIS MILLY,
C     WATER RES. RES., VOL. 18 NO.3 (JUNE 1982), PP. 489-498.
C
      IF(SL-CP(1,NMAT).GE..371d0) GOTO 130
      SLX=MAX(SL,1.001d0*CP(1,NMAT))
      EX=(0.371d0/(SLX-CP(1,NMAT))-1.d0)**.25d0
      EX=2.26d0*EX-2.d0
      PC=-9.7783d3*10.d0**EX
      RETURN
C
  130 PC=-97.783d0
      RETURN
   14 CONTINUE
   15 CONTINUE
C-----LEVERETT@S J-FUNCTION.
      SS=0.d0
      IF(SL.GT.CP(2,NMAT)) SS=(SL-CP(2,NMAT))/(1.d0-CP(2,NMAT))
      OSS=1.d0-SS
      F=1.417d0*OSS-2.120d0*OSS**2+1.263d0*OSS**3
      CALL SIGMA(T,ST)
      PC=-CP(1,NMAT)*ST*F
      RETURN
   16 CONTINUE
C-----CAPILLARY FUNCTION OF VAN GENUCHTEN, SOIL SCI. SOC. AM. J. 44,
C     PP.892-898, 1980.
C
      IF(SL.lt.1.d0)GO TO 160
      PC=0.d0
      RETURN
C
  160 SLX=SL
      IF(SLX.GE.CP(5,NMAT)) GOTO 161
      IF(CP(4,NMAT).EQ.0.d0)SLX=MAX(SL,1.001d0*CP(2,NMAT))
      PC=-ABS(CP(4,NMAT))
c
c----------------------------Active fracture model (Liu et al., 1998)
c
      IF(SLX.GT.CP(2,NMAT))  THEN
         gm=cp(6,nmat)  !Gam factor in active fracture model
         PC=-1.D0/ABS(CP(3,NMAT))*(((SL-CP(2,NMAT))/
     +      (CP(5,NMAT)-CP(2,NMAT)))**
     +      ((gm-1.d0)/CP(1,NMAT))-1.D0)**(1.D0-CP(1,NMAT))
      END IF
c------------------------------------------------------------------
c
      IF(CP(4,NMAT).NE.0.d0) PC=MAX(PC,-ABS(CP(4,NMAT)))
      IF(SL.GT..999d0) PC=PC*(1.d0-SL)/.001d0
      RETURN
  161 PC=0.d0
      RETURN
c
c... ICP = 8: PC = 0.
   17 continue
      pc=0.d0
      return
C
c... ICP = 9: not used
   18 continue
      return
c
c.... ICP = 10: Interpolated saturation
   19 continue
c
      SLX=SL
c
c come here for linear extrapolation at Sw <= Swir
c     icp=10
c     cp(1,nmat)=m, v.G parameter
c     cp(2,nmat)=Swir
c     cp(3,nmat)=alpha, v.G. parameter
c     cp(4,nmat)=e, small number (defaul=1.e-5),
c                for Sw <= Swir+e using linear extrapolation
c     cp(5,nmat)=Sls
c
      IF(SLX.GE.CP(5,NMAT)) GOTO 1161
      IF(SLX.GT.(CP(2,NMAT)+cp(4,nmat))) then
        CPPP=1.0D0-CP(6,NMAT)
        PC=-1.0D0/ABS(CP(3,NMAT))*(((SL-CP(2,NMAT))/
     &   (CP(5,NMAT)-CP(2,NMAT)))
     &   **(-CPPP/CP(1,NMAT))-1.0D0)**(1.0D0-CP(1,NMAT))
      else
c
c  for linear extrapolation with sw <= swir
c
         pc=pc_e(nmat)+pc_slop(nmat)*(slx-cp(2,nmat)-cp(4,nmat))
         pc=-pc
      endif
         pc_cmax=-pc_e(nmat)+pc_slop(nmat)*(cp(2,nmat)+cp(4,nmat))
         pc=max(pc,-abs(pc_cmax))
      RETURN
 1161 PC=0.d0
      RETURN
c
c ICP=11: from S. Finsterle for new V.G. functions
   20 CONTINUE
C     Van Genuchten (only together with IRP=3,4,6,8,10,11);
C     RP(1)=Residual liquid saturation (if CP(6).eq.0)
C     CP(1)=n
C     CP(2)=1/alpha 
C           can also be provided through USERX(2,N), 
C           if CP(2) is negative and USERX(1,N).ne.0,
C           Leverett's scaling rule is applied
C     CP(3)=Epsilon for linear extrapolation for Sl < Slr+Epsilon
C           If greater than 1, CP(3) is maximum pressure
C           If less than zero, use log-linear extrapolation  for Sl < Slr+|Epsilon|
C     CP(4)=m (default: m=1-1/n)
C     CP(5)=if neg. then reference temperature for temperature dependency
C     CP(6)=gamma (Active Fracture Model)
C     CP(7)=Residual liquid saturation (if zero use RP(1))
c     added option for gas res. sat. if cp(7) is negative
c     disabled permeability modifiers from ITOUGH (can add later using
c    equivalent arrays from TOUGH2 V2 (pm)
C     USERX(1,N)=permeability or permeability modifier
C     USERX(2,N)=1/alpha or 1/alpha modifier
C
c
         AE=ABS(CP(2,NMAT))
c........initialize slr
         slr = 0.d0
      IF (CP(7,NMAT).le.0.0D0) THEN
         SLR=abs(RP(1,NMAT))
c changed for option of gas res. sat.    ELSE
      ELSEif(CP(7,NMAT).gt.0.d0)then
         SLR=CP(7,NMAT)
      ENDIF
c Option for gas res. sat.
      if(CP(7,NMAT).lt.0.d0)then
         sgrcp = cp(7,nmat)
      else
         sgrcp = 0.d0
      endif
      SLX=SL
      SE=(SLX-SLR+sgrcp)/(1.0D0-SLR+sgrcp)
      IF (CP(3,NMAT).LE.0.0D0) THEN
         EPSL=CP(3,NMAT)
         SCUT=SLR+ABS(EPSL)
         PCMAX=1.0D50
      ELSE IF (CP(3,NMAT).GE.1.d0) THEN
         EPSL=-1.0D0
         SCUT=SLR+EPSL
         PCMAX=CP(3,NMAT)
      ELSE
         EPSL=CP(3,NMAT)
         SCUT=SLR+EPSL
         PCMAX=1.0D50
      ENDIF
      IF (CP(4,NMAT).LE.0.0D0) THEN
         XN=CP(1,NMAT)
         XM=1.0D0-1.0D0/XN
      ELSE
         XM=CP(4,NMAT)
         XN=1.0D0/(1.0D0-XM)
      ENDIF
c
      IF (SE.GE.1.0D0) THEN
         PC=0.0D0
      ELSE IF (SLX.GT.SCUT) THEN
         IF (SE.LE.0.0D0) THEN
            PC=-PCMAX
         ELSE
            CPPP=1.0D0-CP(6,NMAT)
            PC=-AE*(SE**(-CPPP/XM)-1.0D0)**(1.0D0/XN)
         ENDIF
      ELSE IF (EPSL.GE.0.0D0) THEN
C --- linear extension
         SBAR=EPSL/(1.0D0-SLR)
         CPPP=1.0D0-CP(6,NMAT)
         PCE=AE*(SBAR**(-CPPP/XM)-1.0D0)**(1.0D0/XN)
         PCSLOPE=-CPPP*AE/(XM*XN*(1.0D0-SLR))
     &           *(SBAR**(-CPPP/XM)-1.0D0)**(1.0D0/XN-1.0D0)
     &           *SBAR**((-CPPP-XM)/XM)
         PC=-PCE-PCSLOPE*(SLX-SLR-EPSL)
      ELSE IF (EPSL.LT.0.D0) THEN
C --- log-linear extension 
         SLSTAR=SLR+ABS(EPSL)
c         SE=(SLSTAR-SLR)/(1.0D0-SLR)
         SE=(SLSTAR-SLR+sgrcp)/(1.0D0-SLR+sgrcp)
         CPPP=1.0D0-CP(6,NMAT)
         PCE=AE*(SE**(-CPPP/XM)-1.0D0)**(1.0D0/XN)
         PCSLOPE=-LOG10(EXP(1.0D0))*(-CPPP)/ABS(EPSL)*(1.0D0-XM)/XM/
     &           (SE**(CPPP/XM)-1.0D0)
         PC=-PCE*1.0D1**(PCSLOPE*(SLX-SLSTAR))
      ENDIF
      IF (CP(5,NMAT).LT.-1.0D0) THEN
         ENVG=-CP(5,NMAT)
         CALL SIGMA(ENVG,SBAR)
         CALL SIGMA(T,ST)
         PC=PC*(ST/SBAR-0.0017D0*(T-ENVG))
      ENDIF
      PC=MAX(PC,-PCMAX)
      RETURN
c
      END
c
c
c
      SUBROUTINE SAT(T,P)         
C--------- Fast SAT M.J.O'Sullivan - 17 SEPT 1990 ---------
C
C     20 September 1990.  VAX needs double precision, CRAY does not.
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      COMMON/KONIT/KON,DELT,IGOOD
c
      double precision t, p
C
      save icall
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9/
     1-7.691234564d0,-2.608023696d1,-1.681706546d2,6.423285504d1,
     2-1.189646225d2,4.167117320d0,2.097506760d1,1.d9,6.d0/
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' +++++acquire Psat with given TX in t2f_v2.f+++++'/' SAT'
     x' 1.0S, 1990.9.17: STEAM TABLE EQUATION: SATURATION PRESSURE AS'
     x' FUNCTION OF TEMPERATURE (M. OS.)')
      IF(T.LT.1.d0.OR.T.GT.500.d0) GOTO 10
      TC=(T+273.15d0)/647.3d0
      X1=1.d0-TC
      X2=X1*X1
      SC=A5*X1+A4
      SC=SC*X1+A3
      SC=SC*X1+A2
      SC=SC*X1+A1
      SC=SC*X1
      PC=EXP(SC/(TC*(1.d0+A6*X1+A7*X2))-X1/(A8*X2+A9))
      P=PC*2.212d7
      RETURN
   10 IGOOD=2
      WRITE(34,1) T
    1 FORMAT(' T = ',E12.6,' +++DEGREES CELSIUS, OUT OF RANGE IN SAT, '
     X 'PARAMETER ACQUISITION FAILED, EXIT SUBROUTINE+++')
      RETURN
      END
C
C
C
      SUBROUTINE COWAT(TF,PP,D,U)
C--------- Fast COWAT M.J.O'Sullivan - 17 SEPT 1990 ---------
C     20 September 1990.  VAX needs double precision, CRAY does not.
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      COMMON/KONIT/KON,DELT,IGOOD
      double precision tf,pp,d,u
      save icall
      DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12 /
     16.824687741d3,-5.422063673d2,-2.096666205d4,3.941286787d4,
     2-13.466555478d4,29.707143084d4,-4.375647096d5,42.954208335d4,
     3-27.067012452d4,9.926972482d4,-16.138168904d3,7.982692717d0/
      DATA A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23 /
     4-2.616571843d-2,1.522411790d-3,2.284279054d-2,2.421647003d2,
     51.269716088d-10,2.074838328d-7,2.174020350d-8,1.105710498d-9,
     61.293441934d1,1.308119072d-5,6.047626338d-14/
      DATA SA1,SA2,SA3,SA4,SA5,SA6,SA7,SA8,SA9,SA10,SA11,SA12 /
     18.438375405d-1,5.362162162d-4,1.720000000d0,7.342278489d-2,
     24.975858870d-2,6.537154300d-1,1.150d-6,1.51080d-5,
     31.41880d-1,7.002753165d0,2.995284926d-4,2.040d-1   /
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(6X,'COWAT 1.0S, 1990.9.17: LIQUID WATER DENSITY AND INT.'
     X' ENERGY VERSUS TEMPERATURE AND PRESSURE (M. OS.)')
C
      TKR=(TF+273.15d0)/647.3d0
      TKR2=TKR*TKR
      TKR3=TKR*TKR2
      TKR4=TKR2*TKR2
c not needed because TKR9 not used  TKR5=TKR2*TKR3
      TKR6=TKR4*TKR2
      TKR7=TKR4*TKR3
      TKR8=TKR4*TKR4
c not used      TKR9=TKR4*TKR5
      TKR10=TKR4*TKR6
      TKR11=TKR*TKR10
      TKR19=TKR8*TKR11
      TKR18=TKR8*TKR10
      TKR20=TKR10*TKR10
      PNMR=PP/2.212d7
      PNMR2=PNMR*PNMR
      PNMR3=PNMR*PNMR2
      PNMR4=PNMR*PNMR3
      Y=1.d0-SA1*TKR2-SA2/TKR6
      ZP=SA3*Y*Y-2.d0*SA4*TKR+2.d0*SA5*PNMR
      IF(ZP.LT.0.d0) GOTO 1
      Z=Y+SQRT(ZP)
      CZ=Z**(5.d0/17.d0)
      PAR1=A12*SA5/CZ
      CC1=SA6-TKR
      CC2=CC1*CC1
      CC4=CC2*CC2
      CC8=CC4*CC4
      CC10=CC2*CC8
      AA1=SA7+TKR19
      PAR2=A13+A14*TKR+A15*TKR2+A16*CC10+A17/AA1
      PAR3=(A18+2.d0*A19*PNMR+3.d0*A20*PNMR2)/(SA8+TKR11)
      DD1=SA10+PNMR
      DD2=DD1*DD1
      DD4=DD2*DD2
      PAR4=A21*TKR18*(SA9+TKR2)*(-3.d0/DD4+SA11)
      PAR5=3.d0*A22*(SA12-TKR)*PNMR2+4.d0*A23/TKR20*PNMR3
      VMKR=PAR1+PAR2-PAR3-PAR4+PAR5
      V=VMKR*3.17d-3
      D=1.d0/V
      YD=-2.d0*SA1*TKR+6.d0*SA2/TKR7
      SNUM= A10+A11*TKR
      SNUM=SNUM*TKR + A9
      SNUM=SNUM*TKR + A8
      SNUM=SNUM*TKR + A7
      SNUM=SNUM*TKR + A6
      SNUM=SNUM*TKR + A5
      SNUM=SNUM*TKR + A4
      SNUM=SNUM*TKR2 - A2
      PRT1=A12*(Z*(17.d0*(Z/29.d0-Y/12.d0)+5.d0*TKR*YD/12.d0)+SA4*TKR-
     1(SA3-1.d0)*TKR*Y*YD)/CZ
      PRT2=PNMR*(A13-A15*TKR2+A16*(9.d0*TKR+SA6)*CC8*CC1
     2+A17*(19.d0*TKR19+AA1)/(AA1*AA1))
c     2+A17*(19.*TKR19+AA1)/(AA1*AA1))
      BB1=SA8+TKR11
      BB2=BB1*BB1
      PRT3=(11.d0*TKR11+BB1)/BB2*(A18*PNMR+A19*PNMR2+A20*PNMR3)
      EE1=SA10+PNMR
      EE3=EE1*EE1*EE1
      PRT4=A21*TKR18*(17.d0*SA9+19.d0*TKR2)*(1.d0/EE3+SA11*PNMR)
      PRT5=A22*SA12*PNMR3+21.d0*A23/TKR20*PNMR4
      ENTR= A1*TKR - SNUM +PRT1+PRT2-PRT3+PRT4+PRT5
      H=ENTR*70120.4d0
      U=H-PP*V
      RETURN
    1 IGOOD=2
      WRITE(34,2)TF
C      WRITE(7,2)TF
    2 FORMAT(' T = ',E12.6,' +++DEGREES CELSIUS, OUT OF RANGE IN COWAT'
     X', PARAMETER ACQUISITION FAILED, EXIT SUBROUTINE+++')
  100 FORMAT(1H ,5X,A6,2X,E20.10)
  102 FORMAT(1H ,5X,A6,5X,I2,2X,E20.10)
      RETURN
      END
c
c
c
      SUBROUTINE SUPST(T,P,D,U)
C--------- Fast SUPST M.J.O'Sullivan - 17 SEPT 1990 ---------
C     20 September 1990.  VAX needs double precision, CRAY does not.
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      double precision I1
      double precision t,p,d,u
      save icall
      DATA B0,B01,B02,B03,B04,B05/
     116.83599274d0,28.56067796d0,-54.38923329d0,0.4330662834d0,
     2-0.6547711697d0,8.565182058d-2/
      DATA B11,B12,B21,B22,B23,B31,B32,B41,B42/
     16.670375918d-2,1.388983801d0,8.390104328d-2,2.614670893d-2,
     2-3.373439453d-2,4.520918904d-1,1.069036614d-1,
     3-5.975336707d-1,-8.847535804d-2/
      DATA B51,B52,B53,B61,B62,B71,B72,B81,B82/
     15.958051609d-1,-5.159303373d-1,2.075021122d-1,1.190610271d-1,
     2-9.867174132d-2,1.683998803d-1,-5.809438001d-2,
     36.552390126d-3,5.710218649d-4/
      DATA B90,B91,B92,B93,B94,B95,B96/
     11.936587558d2,-1.388522425d3,4.126607219d3,-6.508211677d3,
     25.745984054d3,-2.693088365d3,5.235718623d2/
      DATA SB,SB61,SB71,SB81,SB82/
     17.633333333d-1,4.006073948d-1,8.636081627d-2,-8.532322921d-1,
     23.460208861d-1/
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' SUPST 1.0S, 1991.2.1: VAPOR DENSITY AND INT. ENERGY AS'
     X' FUNCTION OF TEMPERATURE AND PRESSURE (M. OS.)')
C
      THETA=(T+273.15d0)/647.3d0
      BETA=P/2.212d7
      I1=4.260321148d0
      X=EXP(SB*(1.d0-THETA))
C
      X2=X*X
      X3=X2*X
      X4=X3*X
      X5=X4*X
      X6=X5*X
      X8=X6*X2
      X10=X6*X4
      X11=X10*X
      X14=X10*X4
      X18=X14*X4
      X19=X18*X
      X24=X18*X6
      X27=X24*X3
C
      THETA2=THETA*THETA
      THETA3=THETA2*THETA
      THETA4=THETA3*THETA
C
      BETA2=BETA*BETA
      BETA3=BETA2*BETA
      BETA4=BETA3*BETA
      BETA5=BETA4*BETA
      BETA6=BETA5*BETA
      BETA7=BETA6*BETA
C
      BETAL=15.74373327d0-34.17061978d0*THETA+19.31380707d0*THETA2
      DBETAL=-34.17061978d0+38.62761414d0*THETA
      R=BETA/BETAL
      R2=R*R
      R4=R2*R2
      R6=R4*R2
      R10=R6*R4
C
      CHI2=I1*THETA/BETA
      SC=(B11*X10+B12)*X3
      CHI2=CHI2-SC
      SC=B21*X18+B22*X2+B23*X
      CHI2=CHI2-2.d0*BETA*SC
      SC=(B31*X8+B32)*X10
      CHI2=CHI2-3.d0*BETA2*SC
      SC=(B41*X11+B42)*X14
      CHI2=CHI2-4.d0*BETA3*SC
      SC=(B51*X8+B52*X4+B53)*X24
      CHI2=CHI2-5.d0*BETA4*SC
C
      SD1=1.d0/BETA4+SB61*X14
      SD2=1.d0/BETA5+SB71*X19
      SD3=1.d0/BETA6+(SB81*X27+SB82)*X27
C
      SN=(B61*X+B62)*X11
      chi2=chi2-(sn/sd1*4.d0/beta5)/sd1
      SN=(B71*X6+B72)*X18
      chi2=chi2-(sn/sd2*5.d0/beta6)/sd2
      SN=(B81*X10+B82)*X14
      chi2=chi2-(sn/sd3*6.d0/beta7)/sd3
      SC=B96
      SC=SC*X+B95
      SC=SC*X+B94
      SC=SC*X+B93
      SC=SC*X+B92
      SC=SC*X+B91
      SC=SC*X+B90
      CHI2=CHI2+11.d0*R10*SC
      V=CHI2*0.00317d0
      D=1.d0/V
C
      OS1=SB*THETA
c.....add some constants to speed this up
      os1x24 = 1.d0 + os1*24.d0
      EPS2=0.0d0+B0*THETA-(-B01+B03*THETA2+2*B04*THETA3+3*B05*THETA4)
      SC=(B11*(1.d0+13.d0*OS1)*X10+B12*(1.d0+3.d0*OS1))*X3
      EPS2=EPS2-BETA*SC
      SC=B21*(1.d0+18.d0*OS1)*X18+B22*(1.d0+2.d0*OS1)*X2+B23*
     + (1.d0+OS1)*X
      EPS2=EPS2-BETA2*SC
      SC=(B31*(1.d0+18.d0*OS1)*X8+B32*(1.d0+10.d0*OS1))*X10
      EPS2=EPS2-BETA3*SC
      SC=(B41*(1.d0+25.d0*OS1)*X11+B42*(1.d0+14.d0*OS1))*X14
      EPS2=EPS2-BETA4*SC
      SC=(B51*(1.d0+32.d0*OS1)*X8+B52*(1.d0+28.d0*OS1)*X4+
     1 B53*os1x24)*X24
c     1 B53*(1.d0+24.d0*OS1))*X24
      EPS2=EPS2-BETA5*SC
C
      SN6=14.d0*SB61*X14
      SN7=19.d0*SB71*X19
      SN8=(54.d0*SB81*X27+27.d0*SB82)*X27
      OS5= 1.d0+11.d0*OS1-OS1*SN6/SD1
      SC=(B61*X*(OS1+OS5)+B62*OS5)*(X11/SD1)
      EPS2=EPS2-SC
c      OS6= 1.d0+24.d0*OS1-OS1*SN7/SD2
      OS6= os1x24-OS1*SN7/SD2
      SC=(B71*X6*OS6+B72*(OS6-6.d0*OS1))*(X18/SD2)
      EPS2=EPS2-SC
c      OS7= 1.d0+24.d0*OS1-OS1*SN8/SD3
      OS7= os1x24-OS1*SN8/SD3
      SC=(B81*X10*OS7+B82*(OS7-10.d0*OS1))*(X14/SD3)
      EPS2=EPS2-SC
      OS2=1.d0+THETA*10.d0*DBETAL/BETAL
      SC= (OS2+6.d0*OS1)*B96
c NS3/06 add .d0 to integers below
      SC=SC*X + (OS2+5.d0*OS1)*B95
      SC=SC*X + (OS2+4.d0*OS1)*B94
      SC=SC*X + (OS2+3.d0*OS1)*B93
      SC=SC*X + (OS2+2.d0*OS1)*B92
      SC=SC*X + (OS2+OS1)*B91
      SC=SC*X + OS2*B90
      EPS2=EPS2+BETA*R10*SC
      H=EPS2*70120.4d0
      U=H-P*V
      RETURN
      END
!
!
!
      SUBROUTINE TSAT(PX,TX00,TS)
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      COMMON/KONIT/KON,DELT,IGOOD
C
      double precision px,tx00,ts,ps,tsd,psd
c
C-----FIND SATURATION TEMPERATURE TS AT PRESSURE PX.
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' TSAT 1.0, 1991.3.14: SATURATION TEMPERATURE AS FUNCTION'
     X' OF PRESSURE')
C
      TX0=TX00 ! STARTING TEMPERATURE FOR ITERATION
      IF(TX0.NE.0.d0) GOTO 2
C
C-----COME HERE TO OBTAIN ROUGH STARTING VALUE FOR ITERATION.
      TX0=4606.d0/(24.02d0-LOG(PX)) - 273.15d0
      TX0=MAX(TX0,5.d0)
C
    2 CONTINUE
      TS=TX0
      DT=TS*1.d-8
      TSD=TS+DT
C
    1 CONTINUE
C
      CALL SAT(TS,PS)
      IF(IGOOD.NE.0) RETURN
      IF(ABS((PX-PS)/PX).LE.1.d-10) RETURN
C
      TSD=TS+DT
      CALL SAT(TSD,PSD)
      TS=TS+(PX-PS)*DT/(PSD-PS)
      GOTO 1
C
      END
c
c
c
      SUBROUTINE SIGMA(T,ST)
C
C-----COMPUTE SURFACE TENSION OF WATER, USING THE
C     "INTERNATIONAL REPRESENTATION OF THE SURFACE TENSION OF
C                                               WATER SUBSTANCE" (1975).
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      double precision t,st
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *******SIGMA 1.0, 1990.10.19:SURFACE TENSION OF WATER AS'
     x' FUNCTION OF TEMPERATURE*******')
C
      IF(T.GE.374.15d0) GOTO 1
      ST=1.d0-0.625d0*(374.15d0-T)/647.3d0
      ST=ST*.2358d0*((374.15d0-T)/647.3d0)**1.256d0
      RETURN
C
    1 CONTINUE
      ST=0.d0
      RETURN
      END
c
      SUBROUTINE VIS(T,P,D,VW,VS,PS)
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      double precision T,P,D,VW,VS,PS
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' VIS 1.0, 1990.1.22: VISCOSITY OF LIQUID WATER AND VAPOR'
     X' AS FUNCTION OF TEMPERATURE AND PRESSURE*****')
      EX=247.8d0/(T+133.15d0)
      PHI=1.0467d0*(T-31.85d0)
      AM=1.d0+PHI*(P-PS)*1.d-11
      VW=1.d-7*AM*241.4d0*10.d0**EX
      V1=.407d0*T+80.4d0
      IF(T.LE.350.d0) VS=1.d-7*(V1-D*(1858.d0-5.9d0*T)*1.d-3)
      IF(T.GT.350.d0) VS=1.d-7*(V1+.353d0*D+676.5d-6*D**2+102.1d-9*D**3)
      RETURN
      END
C
      SUBROUTINE VISW(T,P,PS,VW)
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      double precision T,P,VW,PS
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' VISW 1.0, 1990.1.22: VISCOSITY OF LIQUID WATER AS'
     X' FUNCTION OF TEMPERATURE AND PRESSURE')
      EX=247.8d0/(T+133.15d0)
      PHI=1.0467d0*(T-31.85d0)
      AM=1.d0+PHI*(P-PS)*1.d-11
      VW=1.d-7*AM*241.4d0*10.d0**EX
      RETURN
      END
c
      SUBROUTINE VISS(T,P,D,VS)
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      double precision T,P,D,VS
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **VISS 1.0, 1990.1.22: VISCOSITY OF VAPOR AS FUNCTION OF'
     X' TEMPERATURE AND PRESSURE**********')
      V1=.407d0*T+80.4d0
      IF(T.LE.350.d0) VS=1.d-7*(V1-D*(1858.d0-5.9d0*T)*1.d-3)
      IF(T.GT.350.d0) VS=1.d-7*(V1+.353d0*D+676.5d-6*D**2+
     +    102.1d-9*D**3)
      RETURN
      END
c
      SUBROUTINE THERC(T,P,D,CONW,CONS,PS)
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      SAVE ICALL,A0,A1,A2,A3,A4,B0,B1,B2,B3,C0,C1,C2,C3,T0
      DATA A0,A1,A2,A3,A4/-922.47d0,2839.5d0,-1800.7d0,
     +   525.77d0,-73.440d0/
      DATA B0,B1,B2,B3/-.94730d0,2.5186d0,-2.0012d0,.51536d0/
      DATA C0,C1,C2,C3/1.6563d-3,-3.8929d-3,2.9323d-3,-7.1693d-4/
      DATA T0/273.15d0/
C
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **THERC 1.0, 1991.3.4: THERMAL CONDUCTIVITY OF WATER AND'
     X' VAPOR AS FUNCTION OF TEMPERATURE AND PRESSURE********')
C
      T1=(T+273.15d0)/T0
      T2=T1*T1
      T3=T2*T1
      T4=T3*T1
C
c speed this up
      pmps = p - ps
      tsq = t*t
c
C     IF(P-PS.LT.0.d0) GOTO1
      CON1=A0+A1*T1+A2*T2+A3*T3+A4*T4
c      CON2=(P-PS)*(B0+B1*T1+B2*T2+B3*T3)*1.d-5
      CON2=(pmps)*(B0+B1*T1+B2*T2+B3*T3)*1.d-5
c      CON3=(P-PS)*(P-PS)*(C0+C1*T1+C2*T2+C3*T3)*1.d-10
      CON3=pmps*pmps*(C0+C1*T1+C2*T2+C3*T3)*1.d-10
      CONW=(CON1+CON2+CON3)*1.d-3
      CON1=17.6d0+5.87d-2*T
c      CON2=1.04d-4*T*T
      CON2=1.04d-4*tsq
c      CON3=4.51d-8*T**3
      CON3=4.51d-8*tsq*t
      CONS1=1.d-3*(CON1+CON2-CON3)
c      CONS=CONS1+1.d-6*(103.51d0+.4198d0*T-2.771d-5*T*T)*D
      CONS=CONS1+1.d-6*(103.51d0+.4198d0*T-2.771d-5*tsq)*D
     A+1.d-9*D*D*2.1482d14/T**4.2d0
C
C     WRITE (34,10) T,P,PS,CON
   10 FORMAT(5H T = ,E12.6,5H P = ,E12.6,6H PS = ,E12.6,7H CON = ,E12.6)
      RETURN
    1 CONTINUE
    2 FORMAT(8H AT T = ,E12.6,5H P = ,E12.6,19H IS LESS THAN PS = ,E12.6
     A)
      RETURN
      END
c
      SUBROUTINE OUTDF
C
C-----THIS SUBROUTINE GENERATES A PRINTOUT OF DIFFUSIVE FLUXES.
C
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C9/ELEM1(MNCON)
      COMMON/C10/ELEM2(MNCON)
      COMMON/FMOLDIF/FDIF(MNCON*MNPH*MNK)
      COMMON/TITLE/ TITLE
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      CHARACTER*1 HB,H0
      CHARACTER*80 TITLE
      CHARACTER PHAC*12
      CHARACTER*5 ELEM1,ELEM2
      COMMON/FF/H1
      CHARACTER*1 H1
      CHARACTER*3 IJ(5),ijnph
C
      SAVE ICALL,HB,H0,PHAC,IJ
      DATA HB,H0/' ',' '/
      DATA ICALL/0/
      DATA PHAC/'  PHASE COMP'/
      DATA IJ/'-1-','-2-','-3-','-4-','-5-'/
      data ijnph/'all'/
C
C     DEFINE ARITHMETIC FUNCTION II FOR CONVENIENCE.
C
      II(IKK,INP,IN)=(IN-1)*NPH*NK+(INP-1)*NK+IKK
C
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *****OUTDF 1.0, 2001.1.7: PRINT OUT INTERBLOCK DIFFUSIVE'
     X' FLOW RATES**********')
C
C-----PRINT A SHORT HEADER.
C
      WRITE (34,5000) H1,TITLE,KCYC,ITER,SUMTIM
c
      if(mop(24).eq.1) then
      WRITE (34,5060) H0,((PHAC,I=1,NK),J=1,NPH)
      WRITE (34,5061) ((IJ(J),IJ(I),I=1,NK),J=1,NPH)
C
      DO 3030 N=1,NCON
        IF (NEX1(N).EQ.0.OR.NEX2(N).EQ.0) GO TO 3030
C
        IF (MOD(N,57).EQ.53)WRITE(34,5062) H1,((PHAC,I=1,NK),J=1,NPH)
        IF (MOD(N,57).EQ.53)WRITE(34,5061)
     &                              ((IJ(J),IJ(I),I=1,NK),J=1,NPH)
C
        WRITE(34,5071)ELEM1(N),ELEM2(N),((FDIF(II(KK,NP,N)),KK=1,NK),
     &    NP=1,NPH)
C
C
 3030   CONTINUE
C
      endif
C
      if(mop(24).eq.0) then
      WRITE (34,5060) H0,(PHAC,I=1,NK)
      WRITE (34,5061)(IJnph,IJ(I),I=1,NK)
C
      DO 3031 N=1,NCON
        IF (NEX1(N).EQ.0.OR.NEX2(N).EQ.0) GO TO 3031
C
        IF (MOD(N,57).EQ.53) WRITE (34,5062)H1,(PHAC,I=1,NK)
        IF (MOD(N,57).EQ.53) WRITE (34,5061)(IJnph,IJ(I),I=1,NK)
C
c       PRINT 5071,ELEM1(N),ELEM2(N),(FDIF(II(KK,NPH,N)),KK=1,NK)
c        WRITE (34,5071)ELEM1(N),ELEM2(N),(FDIF(II(KK,2,N)),KK=1,NK)
           npheq2 = 2
        WRITE (34,5071)ELEM1(N),ELEM2(N),(FDIF(II(KK,npheq2,N)),KK=1,NK)
C
C
 3031   CONTINUE
C
      endif
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8
 5000 FORMAT(A1/10X,A80/80X,'KCYC =',I7,'  - ITER =',I5,' - TIME ='
     1,E12.6,' MASS FLOW RATES (KG/S) FROM DIFFUSION'/)
 5060 FORMAT(A1,'ELEM1 ELEM2',10A12)
 5061 FORMAT(12X,10('   ',A3,2X,A3,' ')/)
 5062 FORMAT(A1/' ELEM1 ELEM2',10A12)
 5071 FORMAT(1X,A5,1X,A5,(10(1X,E11.5)))
      RETURN
      END
c
c
c
c-----------------------------------------------------------------------
c
      SUBROUTINE INTP1(XT,YT,X,Y,NX)
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer*8 (i-n)
      COMMON/KC/KC
      COMMON/KONIT/KON,DELT,IGOOD
      DIMENSION X(NX),Y(NX)
C
c     IGOOD=0
      IF(NX.EQ.1) GOTO10
      IF(XT.GE.X(NX)) GOTO11
      IF(XT.LE.X(1)) GOTO12
C
      KL=1
      KR=NX
    2 KM=(KL+KR)/2
      IF((KR-KL).EQ.1) GOTO3
      IF(XT.LT.X(KM)) KR=KM
      IF(XT.GE.X(KM)) KL=KM
      GOTO2
C
   11 KL=NX-1
      GOTO3
   12 KL=1
C
    3 F=(XT-X(KL))/(X(KL+1)-X(KL))
      YT=F*(Y(KL+1)-Y(KL))+Y(KL)
C
      RETURN
C
   10 YT=Y(1)
      RETURN
      END
c
c
c
c-----------------------------------------------------------------------
c
c
c
      subroutine kthermfac(ngn,ngkt,timetmp,tceff)
c
c... Finds effective thermal conductivity by linear weighting over time
c      ELS written 8/23/99 to account for combined radiative-conductive
c        heat transfer
c      Uses a linear interpolation between given times and factors that
c        are read from the GENER file
c
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
c... Commons for effective thermal conductivity factors
      common/efkth/timkth(mgtab),fackth(mgtab)
      common/kthtable/ktftb(mnogn)
      double precision tceff,timetmp
      integer*8 ngn,ngkt,nn
c
        nn = 0
      do k = 1, ktftb(ngn) - 1
        nn = ngkt + k
        if( (timetmp.ge.timkth(nn).and.timetmp.lt.timkth(nn+1) )
     +   .or. (timetmp.eq.timkth(nn+1)).and.k.eq.ktftb(ngn) - 1 )then
          factint = (timetmp - timkth(nn))/(timkth(nn+1)-timkth(nn))
          fackthi = factint*(fackth(nn+1)-fackth(nn)) + fackth(nn)
          tceff = tceff*fackthi
        endif
      enddo
c
      return
      end
!
      subroutine fgtab
!
      implicit real*8  (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
!
!.....Other EOSs except for ECO2 and ECO2N
!
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G14/PWB(MNOGN)
      COMMON/G16/GPO(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
      COMMON/C8/GLO(MNCON)
      COMMON/COMPO/FLO(MNPH*MNCON)
      common/fgt1/ioft,iofu,igoft,igofu,noft(100),ngoft(100)
      common/fgt3/icoft,icofu,ncoft(100)
      dimension xd(500),ind(100)
c
      SAVE ICALL,nogu
      DATA ICALL/0/
c
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' FGTAB 1.00, 1998.5.28: Tabulate element, connection, and'
     X' generation data vs. time for plotting**********')
      if(icall.eq.1) nogu=0
c
      if(iofu.gt.0) then
c
C     *****************************************************
C     *              tabulate ELEME data                  *
C     *****************************************************
c
c-----come here to write time data for grid blocks selected in *FOFT*
            pco2=0.d0
            i=0
            do82 n=1,ioft
               if(noft(n).gt.0) then
                  i=i+1
                  nloc=(noft(n)-1)*nk1
                  nloc2=(noft(n)-1)*nsec*neq1
                  if(neq.eq.3) pco2=x(nloc+3)
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8
                  ind(i)=noft(n)
                  xd((i-1)*4+1)=x(nloc+1)
                  xd((i-1)*4+2)=par(nloc2+nsec-1)
                  xd((i-1)*4+3)=par(nloc2+1)
                  xd((i-1)*4+4)=pco2
               endif
   82       continue
c
         endif
         if(iofu.gt.0) write(12,81) kcyc,sumtim,(ind(i),(xd((i-1)*4+k),
     x   k=1,4),i=1,iofu)
   81    format(I5,' , ',E12.6,100(' , ',I5,4(' , ',E14.8)))
c
c
         if(icofu.gt.0) then
c
C     *****************************************************
C     *              tabulate CONNE data                  *
C     *****************************************************
c
c-----come here to write time data for connections selected in *COFT*
            i=0
            do282 j=1,icoft
                  n=ncoft(j)
               if(n.gt.0) then
                  i=i+1
                  ind(i)=n
                  nnp=(n-1)*nph
                  xd((i-1)*3+1)=flo(nnp+1)
                  xd((i-1)*3+2)=flo(nnp+2)
                  xd((i-1)*3+3)=glo(n)
               endif
  282       continue
c
         endif
      if(icofu.gt.0) write(14,281) kcyc,sumtim,(ind(i),(xd((i-1)*3+k),
     xk=1,3),i=1,icofu)
  281    format(I5,' , ',E12.6,100(' , ',I5,3(' , ',E10.4)))
c
c
c     print 183,icall,nogn,igoft,igofu,(ngoft(n),n=1,igoft)
  183 format(' FGTAB !!!  icall =',I3,'   NOGN =',I3,'   IGOFT =',
     xI3,'   IGOFU =',I3,' NGOFT ='/(25(1X,I4)))
c
         if(nogn.gt.0.and.igoft.eq.-1) then
c
c
C     *****************************************************
C     *              tabulate GENER data                  *
C     *****************************************************
c
c-----come here when generation data vs. time for all sinks/sources
c     are to be written on unit 13
            i=0
            do85 n=1,nogn
            j=nexg(n)
            if(j.gt.0) then
               i=i+1
               if(icall.eq.1) nogu=nogu+1
               jloc2=(j-1)*nsec*neq1
c..............compute flowing CO2 mass fraction
               cgn=0.d0
               fft=ff((n-1)*nph+1)+ff((n-1)*nph+2)
               if(nk.eq.2.and.fft.ne.0.d0) cgn=ff((n-1)*nph+1)
     x         *par(jloc2+8)+ff((n-1)*nph+2)*par(jloc2+nbk+8)/fft
               ind(i)=nexg(n)
               xd((i-1)*5+1)=gpo(n)
               xd((i-1)*5+2)=eg(n)
               xd((i-1)*5+3)=cgn
               xd((i-1)*5+4)=ff((n-1)*nph+1)
               xd((i-1)*5+5)=pwb(n)
            endif
   85       continue
         endif
c-----come here when generation data vs. time for only the sinks/
c     sources listed in data block GOFT are to be written on unit 13
         if(nogn.gt.0.and.igofu.gt.0) then
         i=0
         do185 l=1,igoft
c.....assign index of generation item
            n=ngoft(l)
c.....check whether the element listed in GOFT is valid
c     (n.gt.0)
            if(n.gt.0) then
c.....assign element index
            j=nexg(n)
            jloc2=(j-1)*nsec*neq1
            i=i+1
               if(icall.eq.1) nogu=nogu+1
c..............compute flowing CO2 mass fraction
               cgn=0.d0
               fft=ff((n-1)*nph+1)+ff((n-1)*nph+2)
               if(nk.eq.2.and.fft.ne.0.d0) cgn=ff((n-1)*nph+1)
     x         *par(jloc2+8)+ff((n-1)*nph+2)*par(jloc2+nbk+8)/fft
               ind(i)=nexg(n)
               xd((i-1)*5+1)=gpo(n)
               xd((i-1)*5+2)=eg(n)
               xd((i-1)*5+3)=cgn
               xd((i-1)*5+4)=ff((n-1)*nph+1)
               xd((i-1)*5+5)=pwb(n)
            endif
  185       continue
         endif
c
         if(nogu.gt.0) write(13,86) kcyc,sumtim,(ind(n),
     x   (xd((n-1)*5+k),k=1,5),n=1,nogu)
   86    format(I5,' , ',E12.6,100(' , ',I5,5(' , ',E10.4)))
c
      return
      end
!
!
!
      subroutine fgtab_ECO2
!
!.....For ECO2 and ECO2N
!
!
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
c
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G14/PWB(MNOGN)
      COMMON/G16/GPO(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
      COMMON/C8/GLO(MNCON)
      COMMON/COMPO/FLO(MNPH*MNCON)
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      common/fgt1/ioft,iofu,igoft,igofu,noft(100),ngoft(100)
      common/fgt3/icoft,icofu,ncoft(100)
      dimension xd(500),ind(100)
c
      SAVE ICALL,nogu
      DATA ICALL/0/
c
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' FGTAB 1.00, 2005.6.19: Tabulate element, connection, and'
     X' generation data vs. time for plotting, special for ECO2N******')
      if(icall.eq.1) nogu=0
c
      if(iofu.gt.0) then
c
C     *****************************************************
C     *              tabulate ELEME data                  *
C     *****************************************************
c
c-----come here to write time data for grid blocks selected in *FOFT*
            i=0
            do82 n=1,ioft
               if(noft(n).gt.0) then
                  i=i+1
                  nloc=(noft(n)-1)*nk1
                  nloc2=(noft(n)-1)*nsec*neq1
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8
                  ind(i)=noft(n)
c.... 1: pressure, 2: XCO2liq, 3: Sgas, 4: XSliq, 5: SS
                  xd((i-1)*5+1)=x(nloc+1)
                  xd((i-1)*5+2)=par(nloc2+nbk+nb+3)
                  xd((i-1)*5+3)=par(nloc2+1)
                  xd((i-1)*5+4)=par(nloc2+nbk+nb+2)
                  xd((i-1)*5+5)=par(nloc2+2*nbk+1)
               endif
   82       continue
c
         endif
         if(iofu.gt.0) write(12,81) kcyc,sumtim,(ind(i),(xd((i-1)*5+k),
     x   k=1,5),i=1,iofu)
   81    format(I5,' , ',E12.6,100(' , ',I5,5(' , ',E10.4)))
c
c
         if(icofu.gt.0) then
c
C     *****************************************************
C     *              tabulate CONNE data                  *
C     *****************************************************
c
c-----come here to write time data for connections selected in *COFT*
            i=0
            do282 j=1,icoft
                  n=ncoft(j)
               if(n.gt.0) then
c
                  n1=nex1(n)
                  n2=nex2(n)
                  n1loc2=(n1-1)*nsec*neq1
                  n2loc2=(n2-1)*nsec*neq1
                  Xgas1=par(n1loc2+nb+3)
                  Xgas2=par(n2loc2+nb+3)
                  Xliq1=par(n1loc2+nbk+nb+3)
                  Xliq2=par(n2loc2+nbk+nb+3)
c
                  i=i+1
                  ind(i)=n
                  nnp=(n-1)*nph
c.... 1: gas flow, 2: liquid flow, 3: CO2 flow
                  xd((i-1)*3+1)=flo(nnp+1)
                  xd((i-1)*3+2)=flo(nnp+2)
c                 xd((i-1)*3+3)=glo(n)
         if(flo(nnp+1).ge.0.d0) then
            Xfgas=flo(nnp+1)*Xgas2
         else
            Xfgas=flo(nnp+1)*Xgas1
         endif
c
         if(flo(nnp+2).ge.0.d0) then
            Xfliq=flo(nnp+2)*Xliq2
         else
            Xfliq=flo(nnp+2)*Xliq1
         endif
                  Xf=Xfgas+Xfliq
                  xd((i-1)*3+3)=Xf
               endif
  282       continue
c
         endif
      if(icofu.gt.0) write(14,281) kcyc,sumtim,(ind(i),(xd((i-1)*3+k),
     xk=1,3),i=1,icofu)
  281    format(I5,' , ',E12.6,100(' , ',I5,3(' , ',E10.4)))
c
c
c     WRITE (34, 183,icall,nogn,igoft,igofu,(ngoft(n),n=1,igoft)
  183 format(' FGTAB !!!  icall =',I3,'   NOGN =',I3,'   IGOFT =',
     xI3,'   IGOFU =',I3,' NGOFT ='/(25(1X,I4)))
c
         if(nogn.gt.0.and.igoft.eq.-1) then
c
c
C     *****************************************************
C     *              tabulate GENER data                  *
C     *****************************************************
c
c-----come here when generation data vs. time for all sinks/sources
c     are to be written on unit 13
            i=0
            do85 n=1,nogn
            nnp=(n-1)*nph
            j=nexg(n)
            if(j.gt.0) then
               i=i+1
               if(icall.eq.1) nogu=nogu+1
               jloc2=(j-1)*nsec*neq1
c..............compute flowing CO2 mass fraction
               cgn=0.0d0
               fft=ff(nnp+1)+ff(nnp+2)
               if(fft.ne.0.0d0) cgn=(ff(nnp+1)
     x         *par(jloc2+nb+3)+ff(nnp+2)*par(jloc2+nbk+nb+3))/fft
               ind(i)=nexg(n)
c.... 1: mass flow rate, 2: flowing enthalpy, 3: flowing CO2 mass
c        fraction, 4: fractional flow in gas phase, 5: flowing
c        wellbore pressure
               xd((i-1)*5+1)=gpo(n)
               xd((i-1)*5+2)=eg(n)
               xd((i-1)*5+3)=cgn
               xd((i-1)*5+4)=ff((n-1)*nph+1)
               xd((i-1)*5+5)=pwb(n)
            endif
   85       continue
         endif
c-----come here when generation data vs. time for only the sinks/
c     sources listed in data block GOFT are to be written on unit 13
         if(nogn.gt.0.and.igofu.gt.0) then
         i=0
         do185 l=1,igoft
c.....assign index of generation item
            n=ngoft(l)
c.....check whether the element listed in GOFT is valid
c     (n.gt.0)
            if(n.gt.0) then
               nnp=(n-1)*nph
c.....assign element index
               j=nexg(n)
               jloc2=(j-1)*nsec*neq1
               i=i+1
               if(icall.eq.1) nogu=nogu+1
c..............compute flowing CO2 mass fraction
               cgn=0.0d0
               fft=ff(nnp+1)+ff(nnp+2)
               if(fft.ne.0.0d0) cgn=(ff(nnp+1)
     x         *par(jloc2+nb+3)+ff(nnp+2)*par(jloc2+nbk+nb+3))/fft
               ind(i)=nexg(n)
               xd((i-1)*5+1)=gpo(n)
               xd((i-1)*5+2)=eg(n)
               xd((i-1)*5+3)=cgn
               xd((i-1)*5+4)=ff((n-1)*nph+1)
               xd((i-1)*5+5)=pwb(n)
            endif
  185       continue
         endif
c
         if(nogu.gt.0) write(13,86) kcyc,sumtim,(ind(n),
     x   (xd((n-1)*5+k),k=1,5),n=1,nogu)
   86    format(I5,' , ',E12.6,100(' , ',I5,5(' , ',E10.4)))
c
      return
      end
!
