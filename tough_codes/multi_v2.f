
      SUBROUTINE MULTI
c
c------------------------For ECO2 and also shared by EOS3, EOS4, and EOS2
c
c
c     *het/mulh.f*, from MULTI of dvlp/t2f.f
C
C-----SET UP THE COUPLED LINEAR EQUATIONS ARISING AT EACH ITERATION
C     STEP.
C
C***** N O T A T I O N ********************
C
C     NEQ IS THE NUMBER OF EQUATIONS, AND THE NUMBER OF PRIMARY
C         DEPENDENT VARIABLES (PER ELEMENT).
C
C     NPH IS THE NUMBER OF PHASES WHICH CAN BE PRESENT.
C
C     NK  IS THE NUMBER OF COMPONENTS (SPECIES).
C
C     NB  IS THE NUMBER OF PHASE-DEPENDENT SECONDARY VARIABLES OTHER
C         THAN COMPONENT MASS FRACTIONS.
C
C     NBK = NB+NK IS THE TOTAL NUMBER OF PHASE-DEPENDENT SECONDARY
C         VARIABLES.
C
C     NSEC = NPH*NBK+2 IS THE TOTAL NUMBER OF SECONDARY VARIABLES.
C
C******************************************
C
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-n)
C
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
C
      COMMON/E1/ELEM(MNEL)            ! 网格名称
      COMMON/E2/MATX(MNEL)            ! 网格岩性
      COMMON/E3/EVOL(MNEL)            ! 网格体积
      COMMON/E4/PHI(MNEL)             ! 孔隙度
      COMMON/E5/P(MNEL)               ! 压力
      COMMON/E6/T(MNEL)               ! 温度
      common/E7/pm(MNEL)              ! 渗透率
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR PRIMARY VARIABLES $$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NK1
C
      COMMON/P1/X((MNK+1)*MNEL)       ! 主要变量数组
      COMMON/P2/DX((MNK+1)*MNEL)      ! 主要变化的变量数组
      COMMON/P3/DELX((MNK+1)*MNEL)    ! 主要变量的数值微分线性主部系数
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR RESIDUALS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NEQ
C
      COMMON/P4/R(MNEQ*MNEL+1)        ! 残差
      COMMON/P5/DOLD(MNEQ*MNEL)       ! 时间步开始时的累计量
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)           ! 连接中第一个网格的编号
      COMMON/C2/NEX2(MNCON)           ! 连接中第二个网格的编号
      COMMON/C3/DEL1(MNCON)           ! 第一个网格到连接面的距离
      COMMON/C4/DEL2(MNCON)           ! 第二个网格到连接面的距离
      COMMON/C5/AREA(MNCON)           ! 连接面的面积
      COMMON/C6/BETA(MNCON)           ! 连接面法线方向与水平方向家教的余弦值
      COMMON/C7/ISOX(MNCON)           ! 各向异性的参数
      COMMON/C8/GLO(MNCON)            ! 热流速率
      COMMON/C9/ELEM1(MNCON)          ! 连接中第一个网格的名称 
      COMMON/C10/ELEM2(MNCON)         ! 连接中第二个网格的名称
      COMMON/C11/FVD(MNCON)           ! FVD是每个网格的气相流量
      common/c12/sig(MNCON)           ! sig(MNCON)是保存辐射热传递的数组
      COMMON/FMOLDIF/FDIF(MNCON*MNPH*MNK)     ! ???
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)     ! 次级变量的数组
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR FLOW RATES AND VELOCITIES $$$$$$$$$$$$$$$$$
C
      COMMON/COMPO/FLO(MNPH*MNCON)    ! 气相和液相在界面上的流速
      COMMON/PORVEL/VEL(MNPH*MNCON)   ! 气相和液相的孔隙流速
      COMMON/DARVELM/DVELM((2*MNEQ+1)*MNPH*MNCON)     ! 应该是空隙流速的增量
C
CRNT  BEGIN ADDITION FOR RNT
C
      COMMON/EOSEL/IE(16),FE(512)
      COMMON/VOLBC/VBC(MNEL)
C
CRNT  END ADDITION FOR RNT
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR ROCK PROPERTIES $$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +            GK(MAXMAT)
      COMMON/SOLII/XKD3(MAXMAT),XKD4(MAXMAT),SII3(MAXMAT)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR LINEAR EQUATIONS $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/L4/WKAREA(MNEQ*MNEL+10)
      COMMON/L5/IKEEP(5*MNEQ*MNEL)
      COMMON/L6/IW(8*MNEQ*MNEL)
      COMMON/L7/JVECT(niwork)
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)    ! darcy velocity
C
C-----------------------------end   addition for solute transport
C------------------------------------Interface area reduction factor
      common/afactor/a_fm2(mncon) ! advection area reduction (flow from F to M)
      common/afactord/a_fmd(mncon) ! diffusion area reduction (Both sides)
c... Save modified active fracture area for reaction --
c    limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
      common/constraints/sl1min,stimax,dlstmx
c
c------------------------------ For active fracture model
      COMMON/RPCAP/IRP(maxmat),RP(7,maxmat),ICP(maxmat),
     +   CP(7,maxmat),IRPD,RPD(7),ICPD,CPD(7)
c-------------------------------------------------------------------
      COMMON/ICO2/ICO2H2O   ! CO2 and H2O reaction sources considered in the flow
c
c-----------------------------------------Indicators from EOS module
C
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
C
C--------------------------------------------------------------------
C
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
C
C--------------------------------------------------------------------
      COMMON/PARNP/NPL,NPG          ! Phase index
!
      common/control2/ispia,inibound,kcpl
!
c---------------------------------------------------------------------
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/DG/WUP,WNR
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/PATCH/SING
      COMMON/BIND/DIFF0,TEXP,BE
      COMMON/GASLAW/RGAS,AMS,AMA,CVGAS
      COMMON/DFM/TIMAX,REDLT
      COMMON/BC/NELA
      common/ran2/iran
      common/dipa/fddiag(3,5)
      common/ech/eosn(20)
      character*10 eosn
      CHARACTER*5 ELEM,ELEM1,ELEM2,MAT,ELEG
C
c--------------------------- For consistency with ymp version
c... Commons for effective thermal conductivity
      common/efkth/timkth(mgtab),fackth(mgtab)
      common/kthtable/ktftb(mnogn)
c
c... Needed for effective thermal conductivity
      COMMON/G4/ELEG(MNOGN)
c
c... Definitions for passed conductivities and nodes
      double precision con1,con2,timetmp
      integer*8 ngk,ngkt
c
C######### LOCAL ARRAYS ###############################################
C
      DIMENSION D(11,12),F(11,37)
CRNT  BEGIN ADDITION FOR RNT.
C
      DIMENSION XLAM(6),XKD(5,MAXMAT)
      SAVE XLAM,XKD,XMW3,XMW4,XMW43,deg3zero,deg4zero
C
CRNT  END ADDITION FOR RNT.
C
C######################################################################
C
      SAVE ICALL,M11,steb
c     steb is the Stefan-Boltzmann constant for radiative heat
c     emission, in units of J/(m^2 K^4 s)
      DATA ICALL,steb/0,5.6687d-8/
C
c---------- eps is for active fracture model
      DATA eps/1.d-100/
      DATA AA4,AA3,AA2,AA1,AA0/3.62536437E-7,-5.38186884E-4,0.309866854,
     x-80.4072879,12010.1471/
C
      ICALL=ICALL+1
c
c===== ICALL = 1 section
      IF(ICALL.EQ.1) THEN
        WRITE(34,899)
  899   FORMAT(' **MULTI 1.2, 2009.3.20: ASSEMBLE ALL ACCUMULATION AND'
     X' FLOW TERMS, includes capabilities for radiative heat transfer'
     x' and diffusion in all phases with local equilibrium phase'
     X' partitioning between gas and liquid, allows block-by-block'
     x' permeability modification***************')
C
C--------------- For reactive surface area (useful for no-connection blocks)
c------------------ Active fracture model factor for reactions
        DO N=1,NEL
           a_fmr(n) = 1.d0
        END DO
C----------------------------------------------------------------------
        M11=MOD(MOP(11),2)
c
C -----------------EOS7R----------------------------------        
         
        if(eosn(1).eq.'EOS7R     ') then
c+++++ section for EOS7R
           IF(FE(33).EQ.0.d0.OR.FE(41).EQ.0.d0) then
              WRITE (34,8891)
 8891 FORMAT(/'EEEEEEE  INPUT ERROR.  RADIONUCLIDE TRANSPORT',
     &' REQUIRES NON-ZERO ',
     &'VALUES FOR THE HALF-LIVES (FE(33) AND FE(41)).   EEEEEE'/)
              STOP
           endif
           XLAM(1)=0.d0
           XLAM(2)=0.d0
           XLAM(3)=0.693147d0/FE(33)
           XLAM(4)=0.693147d0/FE(41)
           XLAM(5)=0.d0
           XLAM(6)=0.d0
           DO N=1,NM
               XKD(1,N)=0.d0
               XKD(2,N)=0.d0
               XKD(3,N)=XKD3(N)
               XKD(4,N)=XKD4(N)
               XKD(5,N)=0.d0
           END DO
           XMW3=FE(34)
           XMW4=FE(42)
           XMW43=XMW4/XMW3
C
c     if gridblock volume is greater than 1.e50 m^3, assume the
c     radionuclide component does not decay with time in the
c     gridblock.
c
           DO N=1,NELA            ! NELA=NEL
               VBC(N)=1.d0
               IF(EVOL(N).GE.1.d50) THEN
                  VBC(N)=0.d0
               ENDIF
           END DO
      end if
c+++++ end of EOS7R section

C -----------------EOS7R----------------------------------

c     default NER (if converged right away, NER is not defined)
      ner=1

      endif
c
C$$$$$FOR IAB=0 COLUMN INDICES WILL BE STORED IN ICN, OTHERWISE IN JVECT
C-----INITIALIZE COUNTER FOR MATRIX ELEMENTS.
      NZ=0
C
C***** LOOP OVER ELEMENTS ******************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND ONLY UPON VARIABLES PERTAINING
C     TO ONE VOLUME ELEMENT.
C
      IF(MOP(3).GE.2) WRITE (34,201) KCYC,ITER
  201 FORMAT(/' SUBROUTINE MULTI --- [KCYC,ITER] = [',I4,I3,']'/)
      IF(MOP(3).GE.3) WRITE (34,204)
  204 FORMAT(' ===== ACCUMULATION TERMS ===== MASS BALANCES FIRST,',
     X' ENERGY BALANCE LAST'/)
      DO 100 N=1,NELA
C     IDENTIFY MATRIX BLOCK.
      NLOC=(N-1)*NEQ      
      NLOCP=(N-1)*NK1
C     IDENTIFY THE START OF SECONDARY VARIABLES FOR ELEMENT N.
      NLOC2=(N-1)*NSEC*NEQ1
C
C-----ZERO OUT RESIDUALS FOR ELEMENT N.
      DO 20 K=1,NEQ
   20 R(NLOC+K)=0.d0              ! 残差
C
C     ASSIGN QUANTITIES WHICH DEPEND ONLY UPON ELEMENT INDEX N.
      PHIN=PHI(N)
      NMAT=MATX(N)
C-------------------------------------------------------------------
      CD=SH(NMAT)*DM(NMAT)*(1.D0-PHIN)
C-------------------------------------------------------------------
C
C     LOOP THROUGH INCREMENTS FOR DERIVATIVES.

      DO 101 M=1,NEQ1
C     IDENTIFY BEGINNING OF SECONDARY VARIABLES CORRESPONDING
C     TO INCREMENTING PRIMARY VARIABLE (M-1).
      NLM2=NLOC2+(M-1)*NSEC
C
C-----ZERO OUT BLOCK (N,N).
      DO 1011 K=1,NK1
c          1011 D(K,M)=0.D0
        D(K,M)=0.D0
 1011 continue
C
C-----COMPUTE CHANGE IN POROSITY.
      DPRES=0.d0
      IF(M.EQ.2) DPRES=DELX(NLOCP+1)
      PRES=X(NLOCP+1)+DX(NLOCP+1)+DPRES
      DPHI=PHIN*(COM(NMAT)*(PRES-P(N))+EXPAN(NMAT)*(PAR(NLM2+NSEC-1)
     A-T(N)))
      PHINN=PHIN+DPHI
C
      DO 102 NP=1,NPH
C     IDENTIFY BEGINNING OF SECONDARY VARIABLES FOR PHASE NP.
      NL2NP=NLM2+(NP-1)*NBK  ! 每个网格的次级变量的开始位置保存在NL2NP内
C     SATURATION.     饱和度
      SNP=PAR(NL2NP+1)                
      IF(SNP.EQ.0.d0) GOTO 102        ! 如果饱和度等于0就不读取后面的数据并进行计算
C     DENSITY.        密度
      RHONP=PAR(NL2NP+4)  
      PHISRO=PHINN*SNP*RHONP      ! 计算某一相的质量
C
C-----SUM OVER COMPONENTS IN EACH PHASE.
      DO 103 K=1,NK
C     MASS FRACTION OF COMPONENT K IN PHASE NP.
      XNPKM=PAR(NL2NP+NB+K)
      D(K,M)=D(K,M)+XNPKM*PHISRO
      if(eosn(1).eq.'EOS7R     ') then
C
C     INCLUDE ADSORPTION TERM FOR LIQUID PHASE.
C     液相的吸附条件
C
      IF(NP.EQ.2) THEN
          D(K,M)=D(K,M)+(1.d0-PHINN)*DM(NMAT)*
     &           RHONP*XNPKM*XKD(K,NMAT)
      ENDIF
C
C     CALCULATE CURRENT DECAY/DEGRADATION TERM.
c     vbc constant is zero for gridblocks with volume greater
c     than 1.e50 m^3 making decay not occur in that gridblock
C
          DEG=DABS(XLAM(K))*DELT*D(K,M)*VBC(N)
          IF(K.EQ.3) DEG3=DEG
          IF(K.EQ.4) DEG4=DEG-DEG3*XMW43
C
      endif
  103 continue
C
C     INTERNAL ENERGY IN PHASE NP.
      ENNP=PHISRO*PAR(NL2NP+5)-PHINN*SNP*PRES     ! 相NP的内能
      D(NK1,M)=D(NK1,M)+ENNP
C
  102 CONTINUE
C
C     ADD ROCK ENERGY.
      ENR=CD*PAR(NLM2+NSEC-1)     ! 岩石的内能
      D(NK1,M)=D(NK1,M)+ENR
C
C     ASSIGN ACCUMULATION TERM AT BEGINNING OF TIME STEP.
C
      IF(ITER.EQ.1.AND.M.EQ.1) THEN
!
!        IF (MOPR(1) .GT. 0)     THEN            ! not include reactive transport    !!!!!
!     set  kcpl = 0     in CYCIT
!     copy common block of kcpl in MULTI 
!
        IF (MOPR(1) .GT. 0  .or. 
     &      MOPR(1) .eq. 0 .and. KCPL .NE. 1 )     THEN    
!
          DO 105 K=1,NEQ
             DOLD(NLOC+K)=D(K,1)
 105      CONTINUE
        END IF
C     DOLD(NLOC+K) HOLDS THE QUANTITY OF COMPONENT K IN ELEMENT N
C     AT THE END OF THE LAST COMPLETED TIME STEP.
C
      if(eosn(1).eq.'EOS7R     ') then
CRNT  BEGIN ADDITION FOR RNT.
C
C     STORE DECAY/DEGRADATION AT BEGINNING OF TIME STEP.
C
           DEG3ZERO=DEG3
           DEG4ZERO=DEG4
C
      endif
      ENDIF
C
      if(eosn(1).eq.'EOS7R     ') then
      IF(M.EQ.1.AND.MOP(3).GE.3)
     & WRITE (34,205) ELEM(N),(D(KK,1),KK=1,NEQ),DEG
  205 FORMAT(' AT ELEMENT *',A5,'* ',5(1X,E20.14))
C
C     ADD DECAY/DEGRADATION TERM TO ACCOUNT FOR RADIOACTIVE DECAY.
c     if xlam(3) or xlam(4) is input as a negative number, use the
c     midpoint average decay scheme.  if positive, use the default
c     implicit scheme for radioactive decay.
C
      if(xlam(3).lt.0.d0.or.xlam(4).lt.0.d0) then
          D(3,M)=D(3,M) + 0.5d0*(DEG3 + DEG3ZERO)
          D(4,M)=D(4,M) + 0.5d0*(DEG4 + DEG4ZERO)
      else
          D(3,M)=D(3,M) + DEG3
          D(4,M)=D(4,M) + DEG4
      endif
C
      endif
CRNT  END ADDITION FOR RNT.
C
  101 CONTINUE
C
      IF(MOP(3).GE.3.and.eosn(1).ne.'EOS7R     ')
     xWRITE (34,200) ELEM(N),(D(K,1),K=1,NEQ)
  200 FORMAT('       AT ELEMENT *',A5,'*   ',8(1X,E12.6))
     
C+++++ASSIGN ALL ONE-ELEMENT TERMS++++++++++++++++++++++++++++++
    
      DO 106 K=1,NEQ
C     K IS THE ROW INDEX IN MATRIX BLOCK (N,N).
C-----INITIALIZE RESIDUALS.
      R(NLOC+K)=R(NLOC+K)+D(K,1)-DOLD(NLOC+K)
C
      DO 106 L=1,NEQ
C     L IS THE COLUMN INDEX IN MATRIX BLOCK (N,N).
C-----COMPUTE MATRIX ELEMENT ARISING FROM DEPENDENCE OF COMPONENT K
C     UPON PRIMARY VARIABLE L.
      IRN(NZ+1)=NLOC+K
      IF(IAB.EQ.0) ICN(NZ+1)=NLOC+L
      IF(IAB.NE.0) JVECT(NZ+1)=NLOC+L
      CO(NZ+1)=-(D(K,L+1)-D(K,1))/DELX(NLOCP+L)
      IF(PHINN.EQ.0.d0.AND.K.EQ.L.AND.K.NE.NK1) CO(NZ+1)=1.d0
  106 NZ=NZ+1
C
C+++++END OF ASSIGNMENT OF ONE-ELEMENT TERMS++++++++++++++++++++
C
  100 CONTINUE
C
      IF(NOGN.NE.0) CALL QU
!
!-----------------------------------For using EOS2 module (IEOS=2)
!-------------------Modify residule terms due to CO2 and H2O reaction sources
!
      ICO2M = 0
      IF (IEOS.EQ.2 .OR. IEOS.EQ.13 .OR. IEOS.EQ.14) ICO2M=1    ! for EOS2 and ECO2
      IF (MOPR(1).EQ.0 .AND. ICO2M.EQ.1 .AND. ICO2H2O.GT.0) THEN
         CALL QLOSS_Rco2
      END IF
!
!
!-------------------Modify residule terms due to (only) H2O reaction sources
!
!
      IH2OM = 0
      IF (IEOS.EQ.1 .OR. IEOS.EQ.3 .OR. IEOS.EQ.4 .OR. 
     &                   IEOS.EQ.7)        IH2OM=1  
      IF (MOPR(1).EQ.0 .AND. IH2OM.EQ.1 .AND. ICO2H2O.GT.0) THEN
         CALL QLOSS_Rh2o
      END IF
!
!-----------------------------------For using EOS5 module (IEOS=5)
!-------------------Modify residule terms due to H2 reaction sources
!
      IF (MOPR(1).EQ.0 .AND. IEOS.EQ.5 .AND. ICO2H2O.GT.0) THEN
       CALL QLOSS_Rco2
!                   Rco2 store H2 reaction sources here
      END IF
c
c--------------------------------------------------------------------
C***** SINKS AND SOURCES ***********************************************
C
      IF(IGOOD.EQ.3) RETURN
      IF(MOP(15).GE.1) CALL QLOSS
C
C+++++ END OF SINKS/SOURCES ++++++++++++++++++++++++++++++++++++++++++++
C
c.....for modified permeability, scale Pcap
      if(iran.ne.0) then
         if(iter.eq.1) call eos
         n1=nel
         if(icall.gt.1) n1=nela
         do107 n=1,n1
         nloc2=(n-1)*nsec*neq1
         do107 m=1,neq1
         nlm2=nloc2+(m-1)*nsec
         do107 np=1,nph
         nlm2p=nlm2+(np-1)*nbk
  107    if(pm(n).ne.0.d0) par(nlm2p+6)=par(nlm2p+6)/dsqrt(pm(n))
         else
         endif
c---------- Do leverett scaling here for porosity/permeability changes
c      if(mopr(6).eq.1.and.mopr(5).eq.0.and.iran.eq.0)then
       if(mopr(6).eq.1.and.iran.eq.0)then
         if(iter.eq.1) call eos
         n1=nel
         if(icall.gt.1) n1=nela
         do207 n=1,n1
         nloc2=(n-1)*nsec*neq1
         do207 m=1,neq1
         nlm2=nloc2+(m-1)*nsec
         do207 np=1,nph
         nlm2p=nlm2+(np-1)*nbk
c  207    par(nlm2p+6)=par(nlm2p+6)*pcfact(n)
  207    par(nlm2p+6)=par(nlm2p+6)
       endif
c--------- end addition
C
C***** LOOP OVER CONNECTIONS ***********************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND UPON VARIABLES FOR TWO VOLUME
C     ELEMENTS ("INTERFACE QUANTITIES").
C
C     AT EACH INTERFACE HAVE ONE FLUX FOR EACH COMPONENT, WHICH IS A SUM
C     OF CONTRIBUTIONS FROM ALL PHASES. FOR PURPOSES OF COMPUTING
C     DERIVATIVES, NEED ALSO THE FLUXES OBTAINED BY INCREMENTING ANY
C     OF THE NEQ VARIABLES FOR THE "FIRST" ELEMENT, AND ANY OF THE
C     NEQ VARIABLES FOR THE "SECOND" ELEMENT.
C
C-----INITIALIZE TOTAL NUMBER OF FLUX TERMS.
C
      DO 1 N=1,NCON
C
      N1=NEX1(N)      ! N1是连接内第一个网格的网格编号
      N2=NEX2(N)      ! N2是连接内第二个网格的网格编号
!
!.....Addition for active fracture model
!
      a_fm   = 1.0d0      ! Factor for flow and advective transport
      a_fmdd = 1.0d0      ! Factor for diffusive transport
!
      a_fm2(n) = a_fm
      a_fmd(n) = a_fmdd 
!	      
!.....Active fracture model factor for reactions
cc      a_fmr(n1) = 1.d0
cc      a_fmr(n2) = 1.d0
!...............................................
!
      IF(N1.EQ.0.OR.N2.EQ.0) GOTO 1
      IF(N1.GT.NELA.AND.N2.GT.NELA) GOTO 1    ! NELA=NEL
C-----IDENTIFY LOCATION OF MATRIX BLOCKS.
      N1LOC=(N1-1)*NEQ
      N2LOC=(N2-1)*NEQ
      N1LOCP=(N1-1)*NK1
      N2LOCP=(N2-1)*NK1
C
C-----IDENTIFY START OF SECONDARY VARIABLES FOR ELEMENTS N1 AND N2.
      N1LOC2=(N1-1)*NSEC*NEQ1
      N2LOC2=(N2-1)*NSEC*NEQ1
C
C-----OBTAIN SOME QUANTITIES PERTAINING TO CONNECTION.
C
      D1=DEL1(N)
      D2=DEL2(N)
C     SPATIAL INTERPOLATION FACTORS.
      WT1=D2/(D1+D2)
      WT2=1.d0-WT1
c.....3-23-98: assign WM1, WM2 to "something" in case advective
c              section is skipped over (needed in DO 4 loop)
      wm1=0.5d0
      wm2=0.5d0
C
      NMAT1=MATX(N1)
      NMAT2=MATX(N2)
      phi1=phi(n1)
      phi2=phi(n2)
      ISO=ISOX(N)
      GX=BETA(N)*GF
      AX=AREA(N)
c
c.....Fix, in case iso is not between 1 and 3
c       per1 = perm(iso,n1)
c       per2 = perm(iso,n2)
c
      IF (ISO.LE.0 .OR. ISO.GE.4) THEN
         per1 = perm(3,n1)
         per2 = perm(3,n2)
      else
         per1 = perm(iso,n1)
         per2 = perm(iso,n2)
      END IF
c.....
      if(iran.ne.0) then
        per1=per1*pm(n1)
        per2=per2*pm(n2)
      endif
      DPERI=WT1*PER1+WT2*PER2
      PERI=0.d0
      IF(DPERI.NE.0.d0) PERI=PER1*PER2/DPERI
C
C--------------------------------------------------------------
c   not consistent with T2V2   IF(D1.EQ.0.D0) PERI= min(PER2,per1)
c   not consistent with T2V2   IF(D2.EQ.0.D0) PERI= min(PER1,per2)
      IF(D1.EQ.0.D0) PERI = PER2
      IF(D2.EQ.0.D0) PERI = PER1
C---------------------------------------------------------------
C
C-----ASSIGN FACTORS FOR FLUX TERMS.
      IF(N1.LE.NELA) FAC1=FORD/EVOL(N1)      ! FORD是跟时间步长有关的一个东西
      IF(N2.LE.NELA) FAC2=FORD/EVOL(N2)
C
      DPX0=X(N2LOCP+1)+DX(N2LOCP+1)-X(N1LOCP+1)-DX(N1LOCP+1)
      DPX0=DPX0/(D1+D2)
C
      IF(MOP(3).GE.4) WRITE (34,199) N,ELEM(N1),ELEM(N2),FAC1,FAC2
  199 FORMAT(' ***** CONNECTION #',I5,'   ELEMENTS (',A5,',',A5,')',
     X'  ***  (DELT/V1) = ',E12.6,'   (DELT/V2) = ',E12.6)
C
C-----NOW GET ALL FLUX TERMS.
C     M=1 IS THE FLUX FOR LATEST UPDATED VARIABLES.
C     M=2, ... ,NEQ+1 ARE THE NEQ FLUXES CORRESPONDING TO INCREMENTING
C     THE VARIABLES FOR ELEMENT N1.
C  
C     M=NEQ+2, ... ,2*NEQ+1 ARE THE NEQ FLUXES CORRESPONDING TO
C     INCREMENTING THA VARIABLES FOR ELEMENT N2.
C
      DO 2 M=1,NFLUX
C
      DPRES1=0.d0
      DPRES2=0.d0
      IF(M.EQ.2) DPRES1=DELX(N1LOCP+1)
      IF(M.EQ.NEQ+2) DPRES2=DELX(N2LOCP+1)
      PRES1=X(N1LOCP+1)+DX(N1LOCP+1)+DPRES1
      PRES2=X(N2LOCP+1)+DX(N2LOCP+1)+DPRES2
      DPX=(PRES2-PRES1)/(D1+D2)
C
C-----IDENTIFY LOCATION OF SECONDARY VARIABLES.
C
C     DEFINE AUXILIARY QUANTITY M2 TO BE 0 FOR M=1, TO BE 1 FOR
C     SUBSEQUENT GROUP OF NEQ FLUXES INCREMENTED AT N1, TO BE 2 FOR
C     SUBSEQUENT GROUP OF NEQ FLUXES INCREMENTED AT N2.
C
      M2=2
      IF(M.EQ.1) M2=0
      IF(2.LE.M.AND.M.LE.NEQ1) M2=1
      N1LM2=N1LOC2+MOD(M2,2)*(M-1)*NSEC
      N2LM2=N2LOC2+(M2/2)*(M-NEQ1)*NSEC
C
      DO 21 K=1,NK1
   21 F(K,M)=0.d0
C
c.....6-8-95: conductive heat flux only when radiative heat transfer
c             coefficient is non-negative.
      if(sig(n).ge.0.d0) then
C
C-----HEAT CONDUCTIVITY.
          S1=PAR(N1LM2+NBK+1)
          S2=PAR(N2LM2+NBK+1)
      IF(MOP(10).NE.0) GOTO 22
      S1X=MAX(S1,0.d0)
      S2X=MAX(S2,0.d0)
          CON1=CDRY(NMAT1)+DSQRT(S1X)*(CWET(NMAT1)-CDRY(NMAT1))
          CON2=CDRY(NMAT2)+DSQRT(S2X)*(CWET(NMAT2)-CDRY(NMAT2))
      GOTO23
C
   22 CONTINUE
          CON1=CDRY(NMAT1)+S1*(CWET(NMAT1)-CDRY(NMAT1))
          CON2=CDRY(NMAT2)+S2*(CWET(NMAT2)-CDRY(NMAT2))
   23 CONTINUE
c---------- Factor for modification of maximum thermal conductivity
c.... to create an effective thermal conductivity
c.... (radiative + conductive)
          ngkt = 0
          timetmp = sumtim + deltex
        do ngk = 1, nogn
          if(eleg(ngk).eq.elem(n1).and.ktftb(ngk).gt.0)then
            call kthermfac(ngk,ngkt,timetmp,con1)
            ngkt = ktftb(ngk) + ngkt
          endif
        enddo
          ngkt = 0
        do ngk = 1, nogn
          if(eleg(ngk).eq.elem(n2).and.ktftb(ngk).gt.0)then
            call kthermfac(ngk,ngkt,timetmp,con2)
            ngkt = ktftb(ngk) + ngkt
          endif
        enddo
c---------------------------------------------------------------
      DCONI=WT1*CON1+WT2*CON2
      CONI=0.d0
      IF(DCONI.NE.0.d0) CONI=CON1*CON2/DCONI
      IF(D1.EQ.0.d0) CONI=CON2
      IF(D2.EQ.0.d0) CONI=CON1
C
C-----TEMPERATURE GRADIENT.
      T1=PAR(N1LM2+NSEC-1)
      T2=PAR(N2LM2+NSEC-1)
      DTX=(T2-T1)/(D1+D2)
c specific heat of water=12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4
c =AA0+AA1*T+AA2*T*T+AA3*T*T*T+AA4*T*T*T*T
C      T11=T1+273.15D0
C      T22=T2+273.15D0
C      CP1=AA4*T11+AA3
C      CP1=CP1*T11+AA2
C      CP1=CP1*T11+AA1
C      CP1=CP1*T11+AA0
C      CP2=AA4*T22+AA3
C      CP2=CP2*T22+AA2
C      CP2=CP2*T22+AA1
C      CP2=CP2*T22+AA0
C      DCP=CP1*WT1+CP2*WT2
C      CP=0.d0
C      IF(DCP.NE.0.d0) CP=CP1*CP2/DCP
C      IF(D1.EQ.0.d0) CP=CP2
C      IF(D2.EQ.0.d0) CP=CP1
C
C-----CONDUCTIVE HEAT FLOW.
      F(NK1,M)=CONI*DTX*AX
      endif
c.....6-8-95: add radiative heat flux
      if(sig(n).ne.0.d0) then
      t1k=PAR(N1LM2+NSEC-1)+273.15d0
      t2k=PAR(N2LM2+NSEC-1)+273.15d0
      f(nk1,m)=f(nk1,m)+dabs(SIG(n))*steb*ax
     x*(t2k*t2k*t2k*t2k-t1k*t1k*t1k*t1k)
      endif
C
      IF((MOP(3).GE.6.AND.M.EQ.1).OR.MOP(3).GE.8)
     XWRITE (34,198) M,N1LM2,N2LM2,PRES1,PRES2,DPX
  198 FORMAT(/' *** FLUX NO.',I3,' SECONDARY INDICES (',I5,',',I5,
     X') P1 = ',E12.6,' P2 = ',E12.6,' DPX = ',E12.6)
      IF((MOP(3).GE.6.AND.M.EQ.1).OR.MOP(3).GE.8)
     XWRITE (34,197) CONI,PAR(N1LM2+NSEC-1),
     APAR(N2LM2+NSEC-1),DTX,F(NK1,M)
  197 FORMAT(' CONI = ',E12.6,' T1 = ',E12.6,' T2 = ',E12.6,
     A' DTX = ',E12.6,' FNK1 = ',E12.6)
C
C-----OBTAIN FLUX FOR EACH PHASE.
      DO 3 NP=1,NPH
C
C-----------------------------------------------------------------
c.....Fix, in case iso is not between 1 and 3
c       per1 = perm(iso,n1)
c       per2 = perm(iso,n2)
c
      IF (ISO.LE.0 .OR. ISO.GE.4) THEN
         per1 = perm(3,n1)
         per2 = perm(3,n2)
      else
         per1 = perm(iso,n1)
         per2 = perm(iso,n2)
      END IF
c-----------------------------------------------------------------
c.....
      if(iran.ne.0) then
        per1=per1*pm(n1)
        per2=per2*pm(n2)
      endif
      IF(NP.NE.1) GOTO 30
C-----COME HERE TO ASSIGN KLINKENBERG FACTORS
      FK1=GK(NMAT1)/PRES1
      FK2=GK(NMAT2)/PRES2
      PER1=PER1*(1.d0+FK1)
      PER2=PER2*(1.d0+FK2)
   30 CONTINUE
C--------------------------------------------------------------------------
cc      if(eosn(1).eq.'EWASG     ') then
      if(eosn(1).eq.'EWASG     '.or.eosn(1).eq.'ECO2N     ') then
C........APPLY PERMEABILITY REDUCTION DIRECTLY TO ABSOLUTE PERMEABILITY
c         PER1=PER1*PAR(N1LM2+2*NBK+6)
c         PER2=PER2*PAR(N2LM2+2*NBK+6)
         PER1=PER1*PAR(N1LM2+2*NBK+3)
         PER2=PER2*PAR(N2LM2+2*NBK+3)
      endif
C
      DPERI=WT1*PER1+WT2*PER2
      PERI=0.d0
      IF(DPERI.NE.0.d0) PERI=PER1*PER2/DPERI
C
C-----IDENTIFY BEGINNING OF SECONDARY VARIABLES FOR PHASE NP.
C
      N1L2NP=N1LM2+(NP-1)*NBK
      N2L2NP=N2LM2+(NP-1)*NBK
C
      IF(M.EQ.1) then
         VEL((N-1)*NPH+NP)=0.0D0
C
C**************************************Begin addition for solute transport
C
         VELDAR((N-1)*NPH+NP)=0.0D0  ! darcy velocity
C
C**************************************End   addition for solute transport
C
         DVELM((M-1)*NCON*NPH+(N-1)*NPH+NP)=0.0D0
      end if
C
      REL1=PAR(N1L2NP+2)
      REL2=PAR(N2L2NP+2)
C
      FNPM=0.d0
C
      S1=PAR(N1L2NP+1)
      VIS1=PAR(N1L2NP+3)
      RHO1=PAR(N1L2NP+4)
      RHO10=PAR(N1LOC2+(NP-1)*NBK+4)
      H1=PAR(N1L2NP+5)
      PCAP1=PAR(N1L2NP+6)
      PCAP10=PAR(N1LOC2+(NP-1)*NBK+6)
C
      S2=PAR(N2L2NP+1)
      VIS2=PAR(N2L2NP+3)
      RHO2=PAR(N2L2NP+4)
      RHO20=PAR(N2LOC2+(NP-1)*NBK+4)
      H2=PAR(N2L2NP+5)
      PCAP2=PAR(N2L2NP+6)
      PCAP20=PAR(N2LOC2+(NP-1)*NBK+6)
C
c.....2-16-96: the next two statements ('if ... goto 31') were
c              moved from in front of the PAR-assignments to
c              behind them; this way thermophysical properties
c              are still assigned even when advective fluxes are
c              skipped over, because diffusive fluxes may be
c              present.
C-----CHECK WHETHER PHASE NP IS MOBILE.
      IF(REL1.EQ.0.d0.AND.REL2.EQ.0.d0) GOTO31
c.....5-7-93: no advective flux when either permeability is zero.
      if(per1.eq.0.d0.or.per2.eq.0.d0) goto31
c
C-----OBTAIN WEIGHTED INTERFACE DENSITY.
      W1=0.5d0
      IF(RHO1.EQ.0.d0) W1=0.d0
      IF(RHO2.EQ.0.d0) W1=1.d0
      W2=1.d0-W1
      RHOX=W1*RHO1+W2*RHO2
      RHOX0=W1*RHO10+W2*RHO20
C
C-----EFFECTIVE PRESSURE GRADIENT.
      DR=DPX+(PCAP2-PCAP1)/(D1+D2)-RHOX*GX
      DR0=DPX0+(PCAP20-PCAP10)/(D1+D2)-RHOX0*GX
C
C
C-----PERFORM APPROPRIATE UPSTREAM WEIGHTING FOR MOBILITIES.
C
      IF(DR0.GT.0.d0.AND.S2.EQ.0.d0) GOTO 31
      IF(DR0.LT.0.d0.AND.S1.EQ.0.d0) GOTO 31
      IF(M11.GE.1) WM1=W1
c      IF(M11.EQ.0.AND.DR0.GT.0.d0) WM1=1.-WUP
      IF(M11.EQ.0.AND.DR0.GT.0.d0) WM1=1.d0-WUP
      IF(M11.EQ.0.AND.DR0.LE.0.d0) WM1=WUP
C
      IF(RHO1.EQ.0.d0) WM1=0.d0
      IF(RHO2.EQ.0.d0) WM1=1.d0
C
      WM2=1.d0-WM1
C
C*******************************Addition for active fracture model
c.......TX 02/21/2010 Moved up according Eric's
      a_fm   = 1.0d0      ! for flow and advetive transport
      a_fmdd = 1.0d0      ! for diffusive transport
c----------- Active fracture model factor for reactions
cc      a_fmr(n1) = 1.d0
cc      a_fmr(n2) = 1.d0
c
      if (np.eq.npl) then      ! only for liquid phase
      if(isox(n).lt.0) then
c
c Constant reduction factor
c
        if (mop(8).eq.1.or.(isox(n).le.-7.and.isox(n).ge.-9)) then
          if(dr0.le.0.d0) then
             a_fm=rp(6,nmat1)
          else
             a_fm=rp(6,nmat2)
          endif
          if (a_fm.eq.0.0D0) a_fm=1.0d0
          goto 7489
        endif
c
C..... Active Fracture Concept
c
c---------- Residual saturation must be defined for ICP
        if(icp(nmat2).eq.11)then
           slres2 = dabs(rp(1,nmat2))
        elseif(icp(nmat2).eq.7.or.icp(nmat2).eq.10)then
           slres2 = cp(2,nmat2)
        endif
        if(icp(nmat1).eq.11)then
           slres1 = dabs(rp(1,nmat1))
        elseif(icp(nmat1).eq.7.or.icp(nmat1).eq.10)then
           slres1 = cp(2,nmat1)
        endif
c--------------------------------end addition
c
        if(isox(n).le.-10.and.isox(n).ge.-12)then
c------------------------For flow and advective transport
        if(dr0.ge.0.d0.and.d2.le.d1.and.(icp(nmat2).eq.7.
     +     or.icp(nmat2).eq.10.or.icp(nmat2).eq.11))then
            shh=(s2-slres2)/(1.d0-slres2)
              if(shh.gt.0.0d0) then
                a_fm=shh**(1.0d0+cp(6,nmat2))
              else
                a_fm=0.0d0
              endif
          else if(dr0.le.0.d0.and.d2.ge.d1.and.
     +      (icp(nmat1).eq.7.or.icp(nmat1).eq.10.or.icp(nmat1).eq.11))
     +         then
            shh=(s1-slres1)/(1.d0-slres1)
              if(shh.gt.0.0d0) then
                a_fm=shh**(1.0d0+cp(6,nmat1))
              else
                a_fm=0.0d0
              endif
          endif
c----------For diffusive transport (both sides of F and M)
c
           if(d2.le.d1.and.(icp(nmat2).eq.7.or.
     +        icp(nmat2).eq.10.or.icp(nmat2).eq.11))then
            shh=(s2-slres2)/(1.d0-slres2)
            shr=(s2-sl1min)/(1.d0-sl1min)
              if(shh.gt.0.0d0) then
                a_fmdd=shh**(1.0d0+cp(6,nmat2))
              else
                a_fmdd=0.0d0
              endif
              if(shr.gt.0.0d0) then
                a_fmr(n2)=shr**(1.0d0+cp(6,nmat2))
              else
                a_fmr(n2)=0.0d0
              endif
          else if(d2.ge.d1.and.(icp(nmat1).eq.7.or.
     +        icp(nmat1).eq.10.or.icp(nmat1).eq.11))then
            shh=(s2-slres1)/(1.d0-slres1)
            shr=(s1-sl1min)/(1.d0-sl1min)
              if(shh.gt.0.0d0) then
                a_fmdd=shh**(1.0d0+cp(6,nmat1))
              else
                a_fmdd=0.0d0
              endif
              if(shr.gt.0.0d0) then
                a_fmr(n1)=shr**(1.0d0+cp(6,nmat1))
              else
                a_fmr(n1)=0.0d0
              endif
          endif
c
          goto 7489
c
        endif
c
c.....for special weighting of absolute k of vitric/zeolitic connection
c
        if(isox(n).le.-13.and.isox(n).ge.-15) then
          if(per1.gt.per2) then
             per1=per2*rel2/(rel1+eps)
             IF(dabs(D1).le.0.d0.AND.dabs(REL1).le.0.d0) REL1=REL2
          else
             per2=per1*rel1/(rel2+eps)
             IF(dabs(D2).le.0.d0.AND.dabs(REL2).le.0.d0) REL2=REL1
          endif
          goto 7489
        endif
c
c.......upstream saturation or relative permeability
c
        if(dr0.le.0.d0) then
           if (isox(n).ge.-3) then
              a_fm=s1
           else if (isox(n).ge.-6) then
              a_fm=par(n1l2np+2)
           endif
        else
           if (isox(n).ge.-3) then
              a_fm=s2
           else if (isox(n).ge.-6) then
              a_fm=par(n2l2np+2)
           endif
        endif
        if (mop(8).eq.2) then
           if(dr0.le.0.d0) then
              a_fm=a_fm*rp(7,nmat1)
           else
              a_fm=a_fm*rp(7,nmat2)
           endif
        endif
      endif
 7489 continue
      a_fm2(n)=a_fm          ! flow and advection
      a_fmd(n)=a_fmdd        ! diffusion
      end if
C
C
C*******************************************************************************
C
C-----IF A NODAL POINT FALLS RIGHT ON THE INTERFACE (I.E., NODAL
C     DISTANCE = 0), USE RELATIVE PERMEABILITY OF THE OTHER BLOCK.
c.....2-21-99: next two lines commented out; coding revised
c     IF(D1.EQ.0..AND.REL1.EQ.0.) REL1=REL2
c     IF(D2.EQ.0..AND.REL2.EQ.0.) REL2=REL1
         rel10=par(n1loc2+(np-1)*nbk+2)
         rel20=par(n2loc2+(np-1)*nbk+2)
c      IF(D1.EQ.0.d0.AND.rel20.ne.0.d0) REL1=REL2
c      IF(D2.EQ.0.d0.AND.rel10.ne.0.d0) REL2=REL1
C
C-----USE UPSTREAM WEIGHTING FOR ABSOLUTE PERMEABILITY FOR
C     MOP(11) = 0 OR 1; KEEP HARMONIC WEIGHTING FOR MOP(11) .GT. 1.
      IF(MOP(11).LE.1.AND.DR0.GT.0.d0) PERI=PER2
      IF(MOP(11).LE.1.AND.DR0.LE.0.d0) PERI=PER1
C     EXCEPTION: WHEN A NODAL DISTANCE IS ZERO, USE ABSOLUTE
C     PERMEABILITY OF THE OTHER BLOCK.
c
c      IF(D1.EQ.0.D0) PERI= min(PER2,per1)
c      IF(D2.EQ.0.D0) PERI= min(PER1,per2)
      IF(D1.EQ.0.D0) PERI = PER2
      IF(D2.EQ.0.D0) PERI = PER1
C
C-----INTERFACE MOBILITY.
      IF(MOP(11).EQ.4) THEN
      XM1=PER1*REL1/VIS1
      XM2=PER2*REL2/VIS2
      DEN=WT1*XM1+WT2*XM2
      DMOBI=0.d0
      IF(DEN.NE.0.d0) DMOBI=XM1*XM2/DEN
C
      ELSE
      DMOBI=(WM1*REL1/VIS1+WM2*REL2/VIS2)*PERI
c.....2-21-99: commented out
c.....2-20-99
c     if(d1.eq.0.) dmobi=WM2*REL2/VIS2*PERI
c     if(d2.eq.0.) dmobi=WM1*REL1/VIS1*PERI
      ENDIF
C
C-----COME HERE TO COMPUTE PORE VELOCITIES FOR GASEOUS AND LIQUID
C     PHASES, RESPECTIVELY.
C
      IF(M.EQ.1) THEN
        PHIS=WM1*PHI(N1)*S1+WM2*PHI(N2)*S2
        IF(PHIS.NE.0.d0)
     A  VEL((N-1)*NPH+NP)=DMOBI*DR/PHIS
C
C**************************************Begin addition for solute transport
C
        IF(PHIS.NE.0.D0) VELDAR((N-1)*NPH+NP)=DMOBI*DR  ! darcy velocity
C**************************************End   addition for solute transport
C
      ENDIF
C
C***  BEGIN ADDITION FOR T2DM
C
      DVELM((M-1)*NCON*NPH+(N-1)*NPH+NP)=-DMOBI*DR
C.....STORE DARCY VELOCITIES FOR USE IN DISF.
C     SIGN CONVENTION IS OPPOSITE TO COMMON TOUGH2 USAGE.
C     IN TOUGH2, ALL VECTORIAL QUANTITIES ARE REFERRED TO LOCAL
C     COORDINATES, IN WHICH BY CONVENTION THE COORDINATE AXIS
C     POINTS FROM THE SECOND TO THE FIRST GRID BLOCK IN EACH
C     FLOW CONNECTION.
C     IN MESHMAKER Y-Z GRIDS, THE FIRST GRID BLOCK IN A CONNECTION
C     IS ALWAYS "LEFT" OR "ABOVE," SO THAT LOCAL COORDINATES POINT
C     LEFT OR UPWARD, AND FLUXES AND VELOCITY COMPONENTS ARE
C     POSITIVE WHEN THEY POINT "LEFT" OR "UP."
C     THE GLOBAL COORDINATES USED IN DISF POINT IN DIRECTION OF
C     INCREASING Y- AND Z-INDEX, I.E., TO THE RIGHT AND DOWNWARD.
C     THEREFORE, VECTOR COMPONENTS IN GLOBAL DISF COORDINATES
C     HAVE THE OPPOSITE SIGN AS IN LOCAL COORDINATES.
C     8-10-92 (KP)
C
C***  END OF ADDITION FOR T2DM
C
c    7 CONTINUE
C
C
C-----REDEFINE INTERFACE DENSITY WITH UPSTREAM WEIGHTING.
      IF(MOP(18).EQ.0.OR.S1.EQ.0.d0.OR.S2.EQ.0.d0)
     XRHOX=WM1*RHO1+WM2*RHO2
C
C-----FLUX IN PHASE NP.
C**********************************************************************
C      FNPM=DMOBI*RHOX*DR*AX
      FNPM=DMOBI*RHOX*DR*AX*A_FM                  !!!!
C**********************************************************************
C
   31 CONTINUE
C
C-----OBTAIN FLUX FOR COMPONENT K IN PHASE NP.
      DO 4 K=1,NK
C
      XNPMK=WM1*PAR(N1L2NP+NB+K)+WM2*PAR(N2L2NP+NB+K)
C
      F(K,M)=F(K,M)+XNPMK*FNPM
      if(nk.eq.1.or.nb.ne.8) goto 4
C
c-----diffusion occurs only when
c        - there are at least 2 components (nk > 1)
c        - parameter NB = 8 (ingredients of diffusivity stored
c                           in PAR(...+7) and PAR(...+8)
c
c      if(mop(20).eq.0) then
      if(mop(24).eq.0) then
c-----come here to compute binary diffusion, using harmonic weighting
c     for full multiphase diffusive strength.
c     coding is for a 2-phase system with phase 1 - gas, phase 2-
c     liquid, and assuming equilibrium phase partitioning.
c
      if(np.ne.2) goto 4
c     assemble terms only when "second phase" (aqueous) is being
c     worked through.
      if(fddiag(np,k).eq.0.d0.and.fddiag(np-1,k).eq.0.d0) goto 4
c     skip section when diffusivities are zero
C
      X1l=PAR(N1L2NP+NB+K)
      X2l=PAR(N2L2NP+NB+K)
      X1g=PAR(N1L2NP+NB+K-nbk)
      X2g=PAR(N2L2NP+NB+K-nbk)
      if((x2l-x1l).ne.0.d0) then
         hgl=(x2g-x1g)/(x2l-x1l)
      else
         hgl=0.d0
      endif
      DXR=(X2l-X1l)/(D1+D2)
C
c.....calculate entire coefficient to multiply liquid mass
c     fraction gradient; for both grid blocks
c.....liquid phase
         fdd=fddiag(np,k)
         if(fdd.lt.0.d0) then
            dico1l=-phi1*s1*rho1*fdd
            dico2l=-phi2*s2*rho2*fdd
         else
            dico1l=par(n1l2np+7)*par(n1l2np+8)*fdd
            dico2l=par(n2l2np+7)*par(n2l2np+8)*fdd
         endif
c.....gas phase
         fdd=fddiag(np-1,k)
            s1g=par(n1l2np+1-nbk)
            s2g=par(n2l2np+1-nbk)
         if(fdd.lt.0.) then
            rho1g=par(n1l2np+4-nbk)
            rho2g=par(n2l2np+4-nbk)
            dico1g=-phi1*s1g*rho1g*fdd
            dico2g=-phi2*s2g*rho2g*fdd
         else
            dico1g=par(n1l2np+7-nbk)*par(n1l2np+8-nbk)*fdd
            dico2g=par(n2l2np+7-nbk)*par(n2l2np+8-nbk)*fdd
            if(be.ne.0.d0.and.k.eq.1) then
c.....come here to assign enhancement factor for vapor diffusion.
               dico1g=be*par(n1l2np+4-nbk)*par(n1l2np+8-nbk)*fdd
               dico2g=be*par(n2l2np+4-nbk)*par(n2l2np+8-nbk)*fdd
c     vapor enthalpies
               hv1=par(n1lm2+nsec)
               hv2=par(n2lm2+nsec)
            endif
         endif
      dico1=dico1l+hgl*dico1g
      dico2=dico2l+hgl*dico2g
C     HARMONICALLY WEIGHT THE INTERFACE DIFFUSIVITY.
      DEN=dico1*D2 + dico2*D1
c      IF(DEN.EQ.0) THEN
      IF(DEN.EQ.0.d0) THEN
          dico=0.d0
      ELSE
          dico=dico1*dico2*(d1+d2)/den
      ENDIF
C
      fdiff=ax*dico*dxr
      F(K,M)=F(K,M) + FDIFF
c
         if(be.ne.0.d0.and.k.eq.1) then
c...........come here for heat transfer due to enhanced
c           vapor diffusion
            if(fdiff.gt.0.d0.and.s2g.ne.0.d0) gdiff=fdiff*hv2
            if(fdiff.lt.0.d0.and.s1g.ne.0.d0) gdiff=fdiff*hv1
            f(nk1,m)=f(nk1,m)+gdiff
         endif
C
c     other heat transfer due to mass diffusion is ignored.
            if(m.eq.1) then
C.....STORE AWAY DIFFUSIVE FLUXES "AT STATE POINT" FOR PRINTING.
                II=(N-1)*NPH*NK+(NP-1)*NK+K
                FDIF(II)=FDIFF
                if(k.eq.1) fvd(n)=fdiff
            endif
c
      if(m.eq.1.and.mop(3).ge.7)
     xWRITE (34,189) k,dico1l,dico2l,dico1g,dico2g,hgl,dico1,dico2,
     xdico,dxr
  189 format(/' k=',I1,'  dico1l     dico2l       dico1g       dico2g',
     x'        hgl         dico1        dico2         ',
     x'dico         dxr'/(10(1X,E12.6)))
c
c-----end of MOP(24)=0 section for harmonic weighting of coupled multiphase
c     diffusion strength
      endif
c
      if(mop(24).eq.1) then
c-----now a section for computing binary diffusion,
c     with separate harmonic weighting of phase fluxes
c
      if(fddiag(np,k).eq.0.d0) goto 4
C
      X1=PAR(N1L2NP+NB+K)
      X2=PAR(N2L2NP+NB+K)
      DXR=(X2-X1)/(D1+D2)
C
c.....calculate entire coefficient to multiply mass
c     fraction gradient; for both grid blocks
         fdd=fddiag(np,k)
         if(fdd.lt.0.d0) then
            dico1=-phi1*s1*rho1*fdd
            dico2=-phi2*s2*rho2*fdd
         else
            dico1=par(n1l2np+7)*par(n1l2np+8)*fdd
            dico2=par(n2l2np+7)*par(n2l2np+8)*fdd
            if(be.ne.0.d0.and.k.eq.1.and.np.eq.1) then
c.....come here to assign enhancement factor for vapor diffusion.
               dico1=be*par(n1l2np+4)*par(n1l2np+8)*fdd
               dico2=be*par(n2l2np+4)*par(n2l2np+8)*fdd
c     vapor enthalpies
               hv1=par(n1lm2+nsec)
               hv2=par(n2lm2+nsec)
            endif
         endif
C     HARMONICALLY WEIGHT THE INTERFACE DIFFUSIVITY.
      DEN=dico1*D2 + dico2*D1
c      IF(DEN.EQ.0) THEN
      IF(DEN.EQ.0.d0) THEN
          dico=0.d0
      ELSE
          dico=dico1*dico2*(d1+d2)/den
      ENDIF
C
      fdiff=ax*dico*dxr
      F(K,M)=F(K,M) + FDIFF
c
         if(be.ne.0.d0.and.k.eq.1.and.np.eq.1) then
c...........come here for heat transfer due to enhanced
c           vapor diffusion
            if(fdiff.gt.0.d0.and.s2.ne.0.d0) gdiff=fdiff*hv2
            if(fdiff.lt.0.d0.and.s1.ne.0.d0) gdiff=fdiff*hv1
            f(nk1,m)=f(nk1,m)+gdiff
         endif
C
c     other heat transfer due to mass diffusion is ignored.
            if(m.eq.1) then
C.....STORE AWAY DIFFUSIVE FLUXES "AT STATE POINT" FOR PRINTING.
                II=(N-1)*NPH*NK+(NP-1)*NK+K
                FDIF(II)=FDIFF
                if(k.eq.1) fvd(n)=fdiff
            endif
c
      if(m.eq.1.and.mop(3).ge.7) WRITE (34,188) np,k,dico1,dico2,
     x      dico,dxr
  188 format(/' phase #',I2,'   component #',I2,'   dico1=',E12.6,
     x'   dico2=',E12.6,'   dico=',E12.6,'   dxr=',E12.6)
c
c-----end of MOP(20)=1 section with separate harmonic weighting of
c     phase fluxes
      endif
c
    4 CONTINUE
C
      HNPM=WM1*H1+WM2*H2
C
C-----HEAT FLUX IN PHASE NP.

      F(NK1,M)=F(NK1,M)+HNPM*FNPM
c      F(NK1,M)=F(NK1,M)+4.2D3*(T1-T2)*FNPM
C
C
      IF(MOP(3).LT.7) GOTO 190
      IF(MOP(3).LT.9.AND.M.GT.1) GOTO 190
      WRITE (34,196) NP,N1L2NP,N2L2NP,PER1,PER2,PERI
  196 FORMAT(9X,'NP =',I2,' SECONDARY INDICEsasS (',I5,1H,,I5,1H),
     X15X,' PER1 = ',E12.6,' PER2= ',E12.6,'   PERI = ',E12.6)
      WRITE (34,195) S1,REL1,VIS1,RHO1,H1,PCAP1
  195 FORMAT(' S1   = ',E12.6,' REL1 = ',E12.6,' VIS1 = ',E12.6,
     A' RHO1 = ',E12.6,' H1 = ',E12.6,' PCAP1 = ',E12.6)
      WRITE (34,194) S2,REL2,VIS2,RHO2,H2,PCAP2
  194 FORMAT(' S2   = ',E12.6,' REL2 = ',E12.6,' VIS2 = ',E12.6,
     A' RHO2 = ',E12.6,' H2 = ',E12.6,' PCAP2 = ',E12.6)
      WRITE (34,193) DMOBI,RHOX,DR,FNPM,HNPM
  193 FORMAT(' DMOBI= ',E12.6,' RHOX = ',E12.6,' DR = ',E12.6,
     X' FNPM = ',E12.6,' HNPM= ',E12.6)
      WRITE (34,192) (F(K,M),K=1,NK1)
  192 FORMAT(' * FLOW TERMS',8(1X,E12.6))
  190 CONTINUE
C
      IF(M.NE.1) GOTO 3
C-----STORE FLUXES IN EACH PHASE AND HEAT FLUX.
      FLO((N-1)*NPH+NP)=FNPM
C
    3 CONTINUE
    2 CONTINUE
C
      GLO(N)=F(NK1,1)
C
C+++++ASSIGN ALL INTERFACE TERMS+++++++++++++++++++++++++++++++++
C
      DO 5 K=1,NEQ
C     K IS THE ROW INDEX WITHIN A BLOCK PERTAINING TO ELEMENT N1 OR N2.
C
C-----COMPUTE FLUX CONTRIBUTIONS TO RESIDUALS.
      IF(N1.LE.NELA)
     XR(N1LOC+K)=R(N1LOC+K)-FAC1*F(K,1)
      IF(N2.LE.NELA)
     XR(N2LOC+K)=R(N2LOC+K)+FAC2*F(K,1)
C
      DO 6 L=1,NEQ
C     L IS THE COLUMN INDEX WITHIN A BLOCK.
      IF(N1.GT.NELA .OR. N2.GT.NELA) GOTO 61
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N2.
      IRN(NZ+1)=N1LOC+K
      IF(IAB.EQ.0) ICN(NZ+1)=N2LOC+L
      IF(IAB.NE.0) JVECT(NZ+1)=N2LOC+L
      CO(NZ+1)=FAC1*(F(K,L+1+NEQ)-F(K,1))/DELX(N2LOCP+L)
      NZ=NZ+1
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N1.
      IRN(NZ+1)=N2LOC+K
      IF(IAB.EQ.0) ICN(NZ+1)=N1LOC+L
      IF(IAB.NE.0) JVECT(NZ+1)=N1LOC+L
      CO(NZ+1)=-FAC2*(F(K,L+1)-F(K,1))/DELX(N1LOCP+L)
      NZ=NZ+1
   61 CONTINUE
C
C-----DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N1.
C
C     NOTE THAT EACH CONNECTION INVOLVING N1 WILL GENERATE A TERM
C     AT THE SAME MATRIX LOCATION.
      N1KL=(N1-1)*NEQ*NEQ+(K-1)*NEQ+L
      IF(N1.LE.NELA)
     XCO(N1KL)=CO(N1KL)+FAC1*(F(K,L+1)-F(K,1))/DELX(N1LOCP+L)
C
C-----DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N2.
      N2KL=(N2-1)*NEQ*NEQ+(K-1)*NEQ+L
      IF(N2.LE.NELA)
     XCO(N2KL)=CO(N2KL)-FAC2*(F(K,L+1+NEQ)-F(K,1))/DELX(N2LOCP+L)
C
    6 CONTINUE
    5 CONTINUE
C
C+++++END OF ASSIGNMENT OF INTERFACE TERMS+++++++++++++++++++++++
C
    1 CONTINUE
C
C-----TEST FOR CONVERGENCE----------------------------------------------
C
      RERM=0.d0
      IF(MOP(3).GE.2) WRITE (34,203)
  203 FORMAT(/' ===== RESIDUALS ===== MASS BALANCES FIRST,',
     X' ENERGY BALANCE LAST'/)
      DO 10 N=1,NELA
      NLOC=(N-1)*NEQ
C
      IF(MOP(3).GE.2) WRITE (34,202) ELEM(N),(R(NLOC+K),K=1,NEQ)
  202 FORMAT(' AT ELEMENT *',A5,'* ',8(1X,E12.6))
      DO 10 K=1,NEQ
      NLM=NLOC+K
      DOA=DABS(DOLD(NLM))
      IF(DOA.LT.RE2) RER=R(NLM)/RE2
      IF(DOA.GE.RE2) RER=R(NLM)/DOLD(NLM)
      IF(DABS(RER).LE.RERM) GOTO 10
      RERM=DABS(RER)
      NER=N
      KER=K
   10 CONTINUE
C-----------------------------------------------------------------------
      RETURN
      END
C
c-------------------------------------------------------------------------------
c
      SUBROUTINE MULTI_EOS9
c
c     *het/mul9.f*, from MULTI of dvlp/eos9.f
C
C-----SET UP THE COUPLED LINEAR EQUATIONS ARISING AT EACH ITERATION
C     STEP.
C
C***** N O T A T I O N ********************
C
C     NEQ IS THE NUMBER OF EQUATIONS, AND THE NUMBER OF PRIMARY
C         DEPENDENT VARIABLES (PER ELEMENT).
C
C     NPH IS THE NUMBER OF PHASES WHICH CAN BE PRESENT.
C
C     NK  IS THE NUMBER OF COMPONENTS (SPECIES).
C
C     NB  IS THE NUMBER OF PHASE-DEPENDENT SECONDARY VARIABLES OTHER
C         THAN COMPONENT MASS FRACTIONS.
C
C     NBK = NB+NK IS THE TOTAL NUMBER OF PHASE-DEPENDENT SECONDARY
C         VARIABLES.
C
C     NSEC = NPH*NBK+2 IS THE TOTAL NUMBER OF SECONDARY VARIABLES.
C
C******************************************
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      implicit double precision (a-h,o-z)
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
      common/E7/pm(MNEL)
C
      COMMON/MOP_REACT/MOPR(20)  ! control parameters for reactive transport
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR PRIMARY VARIABLES $$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NK1
C
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P3/DELX((MNK+1)*MNEL)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR RESIDUALS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NEQ
C
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/P5/DOLD(MNEQ*MNEL)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
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
      COMMON/C11/FVD(MNCON)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR FLOW RATES AND VELOCITIES $$$$$$$$$$$$$$$$$
C
      COMMON/COMPO/FLO(MNPH*MNCON)
      COMMON/PORVEL/VEL(MNPH*MNCON)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR ROCK PROPERTIES $$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),TORT(MAXMAT),
     +             GK(MAXMAT)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR LINEAR EQUATIONS $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/L4/WKAREA(MNEQ*MNEL+10)
      COMMON/L5/IKEEP(5*MNEQ*MNEL)
      COMMON/L6/IW(8*MNEQ*MNEL)
      COMMON/L7/JVECT(niwork)
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
C
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)    ! darcy velocity
C
C-----------------------------end   addition for solute transport
C--------------------Interface area reduction factor (TX, 12-Aug-1999)
      common/afactor/a_fm2(mncon) ! advection area reduction (flow from F to M)
      common/afactord/a_fmd(mncon) ! diffusion area reduction (Both sides)
c.... Save modified active fracture area for reaction --
c     limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
      common/constraints/sl1min,stimax,dlstmx
c
c------------------------- Added for active fracture model
      COMMON/RPCAP/IRP(maxmat),RP(7,maxmat),ICP(maxmat),
     +   CP(7,maxmat),IRPD,RPD(7),ICPD,CPD(7)
!
      common/control2/ispia,inibound,kcpl
!
c-------------------------------------------------------------------
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/DG/WUP,WNR
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/PATCH/SING
      COMMON/BIND/DIFF0,TEXP,BE
      COMMON/GASL/HC
      COMMON/GASLAW/RGAS,AMS,AMA,CVGAS
      COMMON/DFM/TIMAX,REDLT
      COMMON/BC/NELA
      common/ech/eosn(20)
      character*10 eosn
      common/ran2/iran
      CHARACTER*5 ELEM,ELEM1,ELEM2,MAT
C
C######### LOCAL ARRAYS ###############################################
C
      DIMENSION D(11,12),F(11,23)
C
C######################################################################
C
c---------- eps is for active fracture model
      DATA eps/1.d-100/
C
      SAVE ICALL,M11
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *MULTI 1.0 PM, 1997.3.4: special version for EOS9: sat/'
     X'unsat. flow; random k, gravity flow only (no Pcap) at'
     x' connections where rock densities are .le.0********')

cns10/09 default NER (if converged right away, NER is not defined)
      IF(ICALL.EQ.1) ner=1
    
C
C
C-------------- For reactive surface area (useful for no-connection blocks)
c---------------- Active fracture model factor for reactions
         DO N=1,NEL
            a_fmr(n) = 1.d0
         END DO
C----------------------------------------------------------------------
C
      IF(ICALL.EQ.1) M11=MOD(MOP(11),2)
c
c
C$$$$$FOR IAB=0 COLUMN INDICES WILL BE STORED IN ICN, OTHERWISE IN JVECT
C-----INITIALIZE COUNTER FOR MATRIX ELEMENTS.
      NZ=0
C
C***** LOOP OVER ELEMENTS ******************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND ONLY UPON VARIABLES PERTAINING
C     TO ONE VOLUME ELEMENT.
C
      IF(MOP(3).GE.2) WRITE (34,201) KCYC,ITER
  201 FORMAT(/39H $$$$$$$$$$ SUBROUTINE MULTI $$$$$$$$$$,24H   ---   [KC
     AYC,ITER] = [,I4,1H,,I3,1H]/)
C
      IF(MOP(3).GE.3) WRITE (34,204)
  204 FORMAT(' ===== ACCUMULATION TERMS =====   MASS BALANCES FIRST,',
     X' ENERGY BALANCE LAST'/)
      DO 100 N=1,NELA
C     IDENTIFY MATRIX BLOCK.
      NLOC=(N-1)*NEQ
      NLOCP=(N-1)*NK1
C     IDENTIFY THE START OF SECONDARY VARIABLES FOR ELEMENT N.
      NLOC2=(N-1)*NSEC*NEQ1
C
C-----ZERO OUT RESIDUALS FOR ELEMENT N.
      DO 20 K=1,NEQ
   20 R(NLOC+K)=0.d0
C
C     ASSIGN QUANTITIES WHICH DEPEND ONLY UPON ELEMENT INDEX N.
      PHIN=PHI(N)
      NMAT=MATX(N)
c      CD=SH(NMAT)*DM(NMAT)*(1.-POR(NMAT))
cels10/5/09 not used in eos9  CD=sh(NMAT)*DM(NMAT)*(1.D0-PHIN)
C
C     LOOP THROUGH INCREMENTS FOR DERIVATIVES.
      DO 101 M=1,NEQ1
C     IDENTIFY BEGINNING OF SECONDARY VARIABLES CORRESPONDING
C     TO INCREMENTING PRIMARY VARIABLE (M-1).
      NLM2=NLOC2+(M-1)*NSEC
C
C-----ZERO OUT BLOCK (N,N).
      DO 1011 K=1,NK1
 1011 D(K,M)=0.d0
C
      PHINN=PHIN
C
      DO 102 NP=1,NPH
C     IDENTIFY BEGINNING OF SECONDARY VARIABLES FOR PHASE NP.
      NL2NP=NLM2+(NP-1)*NBK
C     SATURATION.
      SNP=PAR(NL2NP+1)
      IF(SNP.EQ.0.d0) GOTO 102
C     DENSITY.
      RHONP=PAR(NL2NP+4)
      PHISRO=PHINN*SNP*RHONP
C
C-----SUM OVER COMPONENTS IN EACH PHASE.
      DO 103 K=1,NK
C     MASS FRACTION OF COMPONENT K IN PHASE NP.
      XNPKM=PAR(NL2NP+NB+K)
  103 D(K,M)=D(K,M)+XNPKM*PHISRO
  102 CONTINUE
  101 CONTINUE
C
      IF(MOP(3).GE.3) WRITE (34,200) ELEM(N),(D(K,1),K=1,NEQ)
  200 FORMAT(' AT ELEMENT *',A5,'* ',8(1X,E12.6))
C
C+++++ASSIGN ALL ONE-ELEMENT TERMS++++++++++++++++++++++++++++++
C
      DO 105 K=1,NEQ
C     K IS THE ROW INDEX IN MATRIX BLOCK (N,N).
C-----INITIALIZE RESIDUALS.
C     DOLD(NLOC+K) HOLDS THE QUANTITY OF COMPONENT K IN ELEMENT N AT THE
C     END OF LAST COMPLETED TIME STEP.
C
C-----AT BEGINNING OF TIME STEP, ACCUMULATION TERM CONTRIBUTION TO
C     RESIDUAL IS ZERO.
      IF(ITER.NE.1)
     AR(NLOC+K)=R(NLOC+K)+D(K,1)-DOLD(NLOC+K)
C
C-----SET DOLD(NLOC+K) AT BEGINNING OF TIME STEP.
!
      IF(ITER.EQ.1.AND.M.EQ.1) THEN
!
!        IF (MOPR(1) .GT. 0)     THEN            ! not include reactive transport    !!!!!
!     set  kcpl = 0     in CYCIT
!     copy common block of kcpl in MULTI 
!
        IF (MOPR(1) .GT. 0  .or. 
     &      MOPR(1) .eq. 0 .and. KCPL .NE. 1 )     THEN    
!
          DOLD(NLOC+K)=D(K,1)
         END IF
      END IF
!
      DO 106 L=1,NEQ
C     L IS THE COLUMN INDEX IN MATRIX BLOCK (N,N).
C-----COMPUTE MATRIX ELEMENT ARISING FROM DEPENDENCE OF COMPONENT K
C     UPON PRIMARY VARIABLE L.
      IRN(NZ+1)=NLOC+K
      IF(IAB.EQ.0) ICN(NZ+1)=NLOC+L
      IF(IAB.NE.0) JVECT(NZ+1)=NLOC+L
      CO(NZ+1)=-(D(K,L+1)-D(K,1))/DELX(NLOCP+L)
      IF(PHINN.EQ.0.d0.AND.K.EQ.L.AND.K.NE.NK1) CO(NZ+1)=1.d0
  106 NZ=NZ+1
  105 CONTINUE
C
C+++++END OF ASSIGNMENT OF ONE-ELEMENT TERMS++++++++++++++++++++
C
  100 CONTINUE
C
C***** SINKS AND SOURCES ***********************************************
C
      IF(NOGN.NE.0) CALL QU
      IF(IGOOD.EQ.3) RETURN
C
C+++++ END OF SINKS/SOURCES ++++++++++++++++++++++++++++++++++++++++++++
C
c.....for modified permeability, scale Pcap
      if(iran.ne.0) then
         if(iter.eq.1) call eos
         n1=nel
         if(icall.gt.1) n1=nela
         do107 n=1,n1
         nloc2=(n-1)*nsec*neq1
         do107 m=1,neq1
         nlm2=nloc2+(m-1)*nsec
         do107 np=1,nph
         nlm2p=nlm2+(np-1)*nbk
  107    if(pm(n).ne.0.d0) par(nlm2p+6)=par(nlm2p+6)/sqrt(pm(n))
         else
         endif
c---------- do leverett scaling here for porosity/permeability changes
c       if(mopr(6).eq.1.and.mopr(5).eq.0.and.iran.eq.0)then
       if(mopr(6).eq.1.and.iran.eq.0)then
         if(iter.eq.1) call eos
         n1=nel
         if(icall.gt.1) n1=nela
         do207 n=1,n1
         nloc2=(n-1)*nsec*neq1
         do207 m=1,neq1
         nlm2=nloc2+(m-1)*nsec
         do207 np=1,nph
         nlm2p=nlm2+(np-1)*nbk
c  207    par(nlm2p+6)=par(nlm2p+6)*pcfact(n)
  207    par(nlm2p+6)=par(nlm2p+6)
       endif
c----------------------------------------------------- end addition
C
C***** LOOP OVER CONNECTIONS ***********************************       
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND UPON VARIABLES FOR TWO VOLUME
C     ELEMENTS ("INTERFACE QUANTITIES").
C
C     AT EACH INTERFACE HAVE ONE FLUX FOR EACH COMPONENT, WHICH IS A SUM
C     OF CONTRIBUTIONS FROM ALL PHASES. FOR PURPOSES OF COMPUTING
C     DERIVATIVES, NEED ALSO THE FLUXES OBTAINED BY INCREMENTING ANY
C     OF THE NEQ VARIABLES FOR THE "FIRST" ELEMENT, AND ANY OF THE
C     NEQ VARIABLES FOR THE "SECOND" ELEMENT.
C
C-----INITIALIZE TOTAL NUMBER OF FLUX TERMS.
C
C
      DO1 N=1,NCON
C
      N1=NEX1(N)
      N2=NEX2(N)
!
!.....Addition for active fracture model
!
      a_fm   = 1.0d0      ! Factor for flow and advective transport
      a_fmdd = 1.0d0      ! Factor for diffusive transport
!
      a_fm2(n) = a_fm
      a_fmd(n) = a_fmdd 
!.......................................	      
!
      IF(N1.EQ.0.OR.N2.EQ.0) GOTO 1
      IF(N1.GT.NELA.AND.N2.GT.NELA) GOTO 1    ! NELA=NEL
C-----IDENTIFY LOCATION OF MATRIX BLOCKS.
      N1LOC=(N1-1)*NEQ
      N2LOC=(N2-1)*NEQ
      N1LOCP=(N1-1)*NK1
      N2LOCP=(N2-1)*NK1
C
C-----IDENTIFY START OF SECONDARY VARIABLES FOR ELEMENTS N1 AND N2.
      N1LOC2=(N1-1)*NSEC*NEQ1
      N2LOC2=(N2-1)*NSEC*NEQ1
C
C-----OBTAIN SOME QUANTITIES PERTAINING TO CONNECTION.
C
      D1=DEL1(N)
      D2=DEL2(N)
C     SPATIAL INTERPOLATION FACTORS.
      WT1=D2/(D1+D2)
      WT2=1.d0-WT1
C
      NMAT1=MATX(N1)
      NMAT2=MATX(N2)
c      PO1=POR(NMAT1)
c      PO2=POR(NMAT2)
      ISO=ISOX(N)
      GX=BETA(N)*GF
      AX=AREA(N)
c-------------------------------------------------------------
c.....Fix, in case iso is not between 1 and 3
c       per1 = perm(iso,n1)
c       per2 = perm(iso,n2)
c
      IF (ISO.LE.0 .OR. ISO.GE.4) THEN
         per1 = perm(3,n1)
         per2 = perm(3,n2)
      else
         per1 = perm(iso,n1)
         per2 = perm(iso,n2)
      END IF
c-------------------------------------------------------------
c.....
      if(iran.ne.0) then
      per1=per1*pm(n1)
      per2=per2*pm(n2)
      else
      endif
      DPERI=WT1*PER1+WT2*PER2
      PERI=0.d0
      IF(DPERI.NE.0.d0) PERI=PER1*PER2/DPERI
c---------- Modified so that flow is limited by block of lower permeability
      IF(D1.EQ.0.D0) PERI = PER2
      IF(D2.EQ.0.D0) PERI = PER1
C
C-----ASSIGN FACTORS FOR FLUX TERMS.
      IF(N1.LE.NELA) FAC1=FORD/EVOL(N1)
      IF(N2.LE.NELA) FAC2=FORD/EVOL(N2)
C
C
      IF(MOP(3).GE.4) WRITE (34,199) N,ELEM(N1),ELEM(N2),FAC1,FAC2
  199 FORMAT(/' ***** CONNECTION #',I5,' ELEMENTS (',A5,',',A5,')',
     X' *** (DELT/V1) = ',E12.6,' (DELT/V2) = ',E12.6)
C
C-----NOW GET ALL FLUX TERMS.
C     M=1 IS THE FLUX FOR LATEST UPDATED VARIABLES.
C     M=2, ... ,NEQ+1 ARE THE NEQ FLUXES CORRESPONDING TO INCREMENTING
C     THE VARIABLES FOR ELEMENT N1.
C
C     M=NEQ+2, ... ,2*NEQ+1 ARE THE NEQ FLUXES CORRESPONDING TO
C     INCREMENTING THA VARIABLES FOR ELEMENT N2.
C
      DO 2 M=1,NFLUX
C
C-----IDENTIFY LOCATION OF SECONDARY VARIABLES.
C
C     DEFINE AUXILIARY QUANTITY M2 TO BE 0 FOR M=1, TO BE 1 FOR
C     SUBSEQUENT GROUP OF NEQ FLUXES INCREMENTED AT N1, TO BE 2 FOR
C     SUBSEQUENT GROUP OF NEQ FLUXES INCREMENTED AT N2.
C
      M2=2
      IF(M.EQ.1) M2=0
      IF(2.LE.M.AND.M.LE.NEQ1) M2=1
      N1LM2=N1LOC2+MOD(M2,2)*(M-1)*NSEC
      N2LM2=N2LOC2+(M2/2)*(M-NEQ1)*NSEC
C
C
      DO 21 K=1,NK1
   21 F(K,M)=0.d0
C
C
      IF((MOP(3).GE.6.AND.M.EQ.1).OR.MOP(3).GE.8)
     X WRITE (34,198) M,N1LM2,N2LM2
  198  FORMAT(/' *** FLUX NO.',I3,' SECONDARY INDICES (',I5,',',I5,')')
C
C-----OBTAIN FLUX FOR EACH PHASE.
      DO 3 NP=1,NPH
C
c      PER1=PER(ISO,NMAT1)
c      PER2=PER(ISO,NMAT2)
c.....Fix, in case iso is not between 1 and 3
c       per1 = perm(iso,n1)
c       per2 = perm(iso,n2)
c
      IF (ISO.LE.0 .OR. ISO.GE.4) THEN
         per1 = perm(3,n1)
         per2 = perm(3,n2)
      else
         per1 = perm(iso,n1)
         per2 = perm(iso,n2)
      END IF
c------------------------------------------------------------------
c.....
      if(iran.ne.0) then
      per1=per1*pm(n1)
      per2=per2*pm(n2)
      else
      endif
C
      DPERI=WT1*PER1+WT2*PER2
      PERI=0.d0
      IF(DPERI.NE.0.d0) PERI=PER1*PER2/DPERI
C
C-----IDENTIFY BEGINNING OF SECONDARY VARIABLES FOR PHASE NP.
C
      N1L2NP=N1LM2+(NP-1)*NBK
      N2L2NP=N2LM2+(NP-1)*NBK
C
      IF(M.EQ.1) then
         VEL((N-1)*NPH+NP)=0.0D0
C
C**************************************Begin addition for solute transport
C
         VELDAR((N-1)*NPH+NP)=0.0D0  ! darcy velocity
C
C**************************************End   addition for solute transport
C
      END IF
C
      REL1=PAR(N1L2NP+2)
      REL2=PAR(N2L2NP+2)
C
      FNPM=0.d0
C-----CHECK WHETHER PHASE NP IS MOBILE.
      IF(REL1.EQ.0.d0.AND.REL2.EQ.0.d0) GOTO31
c.....5-7-93: no advective flux when either permeability is zero.
      if(per1.eq.0.d0.or.per2.eq.0.d0) goto31
C
      S1=PAR(N1L2NP+1)
      VIS1=PAR(N1L2NP+3)
      RHO1=PAR(N1L2NP+4)
      RHO10=PAR(N1LOC2+(NP-1)*NBK+4)
      H1=PAR(N1L2NP+5)
      PCAP1=PAR(N1L2NP+6)
      PCAP10=PAR(N1LOC2+(NP-1)*NBK+6)
C
      S2=PAR(N2L2NP+1)
      VIS2=PAR(N2L2NP+3)
      RHO2=PAR(N2L2NP+4)
      RHO20=PAR(N2LOC2+(NP-1)*NBK+4)
      H2=PAR(N2L2NP+5)
      PCAP2=PAR(N2L2NP+6)
      PCAP20=PAR(N2LOC2+(NP-1)*NBK+6)
C
C-----OBTAIN WEIGHTED INTERFACE DENSITY.
      W1=0.5d0
      IF(RHO1.EQ.0.d0) W1=0.d0
      IF(RHO2.EQ.0.d0) W1=1.d0
      W2=1.d0-W1
      RHOX=W1*RHO1+W2*RHO2
      RHOX0=W1*RHO10+W2*RHO20
C
C-----EFFECTIVE PRESSURE GRADIENT.
c.....11-2-95: revise to allow for possibility of ignoring the
c              capillary part of eff. pressure gradient.
c              Ignore capillary gradient at connections where
c              one or both grid blocks have rock density.le.0.
c              This will allow to model gravity-driven outflow.
      if(dm(nmat1).gt.0.d0.and.dm(nmat2).ge.0.d0) then
         DR=(PCAP2-PCAP1)/(D1+D2)-RHOX*GX
         DR0=(PCAP20-PCAP10)/(D1+D2)-RHOX0*GX
      else
         DR=-RHOX*GX
         DR0=-RHOX0*GX
      endif
C
C-----PERFORM APPROPRIATE UPSTREAM WEIGHTING FOR MOBILITIES.
C
      IF(DR0.GT.0.d0.AND.S2.EQ.0.d0) GOTO 31
      IF(DR0.LT.0.d0.AND.S1.EQ.0.d0) GOTO 31
      IF(M11.GE.1) WM1=W1
      IF(M11.EQ.0.AND.DR0.GT.0.d0) WM1=1.d0-WUP
      IF(M11.EQ.0.AND.DR0.LE.0.d0) WM1=WUP
C
      IF(RHO1.EQ.0.d0) WM1=0.d0
      IF(RHO2.EQ.0.d0) WM1=1.d0
C
      WM2=1.d0-WM1
C
C***************************************************************************
C-----IF A NODAL POINT FALLS RIGHT ON THE INTERFACE (I.E., NODAL **********
C     DISTANCE = 0), USE RELATIVE PERMEABILITY OF THE OTHER BLOCK.
C
C***************************************Addition for active fracture model
      a_fm   = 1.0d0      ! for flow and advetive transport
      a_fmdd = 1.0d0      ! for diffusive transport
      if(isox(n).lt.0) then
c
c Constant reduction factor
c
        if (mop(8).eq.1.or.(isox(n).le.-7.and.isox(n).ge.-9)) then
          if(dr0.le.0.d0) then
             a_fm=rp(6,nmat1)
          else
             a_fm=rp(6,nmat2)
          endif
          if (a_fm.eq.0.0D0) a_fm=1.0d0
          goto 7489
        endif
c
C HHLIUB Active Fracture Concept
c---------------- Modified from ysw and corrected a bug
c
c---------------- Residual saturation must be defined for ICP
        if(icp(nmat2).eq.11)then
           slres2 = abs(rp(1,nmat2))
        elseif(icp(nmat2).eq.7.or.icp(nmat2).eq.10)then
           slres2 = cp(2,nmat2)
        endif
        if(icp(nmat1).eq.11)then
           slres1 = abs(rp(1,nmat1))
        elseif(icp(nmat1).eq.7.or.icp(nmat1).eq.10)then
           slres1 = cp(2,nmat1)
        endif
c------------------------------- end addition
c
        if(isox(n).le.-10.and.isox(n).ge.-12)then
c------------------------For flow and advective transport
        if(dr0.ge.0.d0.and.d2.le.d1.and.(icp(nmat2).eq.7.
     +     or.icp(nmat2).eq.10.or.icp(nmat2).eq.11))then
            shh=(s2-slres2)/(1.d0-slres2)
              if(shh.gt.0.0d0) then
                a_fm=shh**(1.0d0+cp(6,nmat2))
              else
                a_fm=0.0d0
              endif
          else if(dr0.le.0.d0.and.d2.ge.d1.and.
     +      (icp(nmat1).eq.7.or.icp(nmat1).eq.10.or.icp(nmat1).eq.11))
     +         then
            shh=(s1-slres1)/(1.d0-slres1)
              if(shh.gt.0.0d0) then
                a_fm=shh**(1.0d0+cp(6,nmat1))
              else
                a_fm=0.0d0
              endif
          endif
c----------For diffusive transport (both sides of F and M)
           if(d2.le.d1.and.(icp(nmat2).eq.7.or.
     +        icp(nmat2).eq.10.or.icp(nmat2).eq.11))then
            shh=(s2-slres2)/(1.d0-slres2)
            shr=(s2-sl1min)/(1.d0-sl1min)
              if(shh.gt.0.0d0) then
                a_fmdd=shh**(1.0d0+cp(6,nmat2))
              else
                a_fmdd=0.0d0
              endif
              if(shr.gt.0.0d0) then
                a_fmr(n2)=shr**(1.0d0+cp(6,nmat2))
              else
                a_fmr(n2)=0.0d0
              endif
          else if(d2.ge.d1.and.(icp(nmat1).eq.7.or.
     +        icp(nmat1).eq.10.or.icp(nmat1).eq.11))then
            shh=(s2-slres1)/(1.d0-slres1)
            shr=(s1-sl1min)/(1.d0-sl1min)
              if(shh.gt.0.0d0) then
                a_fmdd=shh**(1.0d0+cp(6,nmat1))
              else
                a_fmdd=0.0d0
              endif
              if(shr.gt.0.0d0) then
                a_fmr(n1)=shr**(1.0d0+cp(6,nmat1))
              else
                a_fmr(n1)=0.0d0
              endif
          endif
c
          goto 7489
c
        endif
C HHLIUE
c
c for special weighting of absolute k of vitric/zeolitic connection
c
        if(isox(n).le.-13.and.isox(n).ge.-15) then
          if(per1.gt.per2) then
             per1=per2*rel2/(rel1+eps)
             IF(abs(D1).le.0.d0.AND.abs(REL1).le.0.d0) REL1=REL2
          else
             per2=per1*rel1/(rel2+eps)
             IF(abs(D2).le.0.d0.AND.abs(REL2).le.0.d0) REL2=REL1
          endif
          goto 7489
        endif
c
c Upstream saturation or relative permeability
c
        if(dr0.le.0.d0) then
           if (isox(n).ge.-3) then
              a_fm=s1
           else if (isox(n).ge.-6) then
              a_fm=par(n1l2np+2)
           endif
        else
           if (isox(n).ge.-3) then
              a_fm=s2
           else if (isox(n).ge.-6) then
              a_fm=par(n2l2np+2)
           endif
        endif
        if (mop(8).eq.2) then
           if(dr0.le.0.d0) then
              a_fm=a_fm*rp(7,nmat1)
           else
              a_fm=a_fm*rp(7,nmat2)
           endif
        endif
      endif
 7489 continue
      a_fm2(n)=a_fm          ! flow and advection
      a_fmd(n)=a_fmdd        ! diffusion
C YSWE
C***********************************************************************
c
C-----IF A NODAL POINT FALLS RIGHT ON THE INTERFACE (I.E., NODAL
C     DISTANCE = 0), USE RELATIVE PERMEABILITY OF THE OTHER BLOCK.
c      IF(D1.EQ.0.d0.AND.REL1.EQ.0.d0) REL1=REL2
c      IF(D2.EQ.0.d0.AND.REL2.EQ.0.d0) REL2=REL1
C
C-----USE UPSTREAM WEIGHTING FOR ABSOLUTE PERMEABILITY FOR
C     MOP(11) = 0 OR 1; KEEP HARMONIC WEIGHTING FOR MOP(11) .GT. 1.
      IF(MOP(11).LE.1.AND.DR0.GT.0.d0) PERI=PER2
      IF(MOP(11).LE.1.AND.DR0.LE.0.d0) PERI=PER1
C     EXCEPTION: WHEN A NODAL DISTANCE IS ZERO, USE ABSOLUTE
C     PERMEABILITY OF THE OTHER BLOCK.
c--------- Modified so that flow is limited by block of lower permeability
c      IF(D1.EQ.0.D0) PERI=PER2
c      IF(D2.EQ.0.D0) PERI=PER1
      IF(D1.EQ.0.D0) PERI= min(PER2,per1)
      IF(D2.EQ.0.D0) PERI= min(PER1,per2)
C
C-----INTERFACE MOBILITY.
      IF(MOP(11).EQ.4) THEN
      XM1=PER1*REL1/VIS1
      XM2=PER2*REL2/VIS2
      DEN=WT1*XM1+WT2*XM2
      DMOBI=0.d0
      IF(DEN.NE.0.d0) DMOBI=XM1*XM2/DEN
C
      ELSE
      DMOBI=(WM1*REL1/VIS1+WM2*REL2/VIS2)*PERI
      ENDIF
C
      IF(M.NE.1) GOTO 7
C-----COME HERE TO COMPUTE PORE VELOCITIES FOR GASEOUS AND LIQUID
C     PHASES, RESPECTIVELY.
C
      PHIS=WM1*PHI(N1)*S1+WM2*PHI(N2)*S2
      IF(PHIS.NE.0.d0)   THEN
        VEL((N-1)*NPH+NP)=DMOBI*DR/PHIS
C
C**************************************Begin addition for solute transport
C
c        VELDAR((N-1)*NPH+NP)=DMOBI*DR  ! darcy velocity
      IF(PHIS.NE.0.D0)VELDAR((N-1)*NPH+NP)=DMOBI*DR
C**************************************End   addition for solute transport
C
      END IF
c
    7 CONTINUE
C
C-----REDEFINE INTERFACE DENSITY WITH UPSTREAM WEIGHTING.
      IF(MOP(18).EQ.0.OR.S1.EQ.0.d0.OR.S2.EQ.0.d0)
     XRHOX=WM1*RHO1+WM2*RHO2
C-----FLUX IN PHASE NP.
c      FNPM=DMOBI*RHOX*DR*AX
      FNPM=DMOBI*RHOX*DR*AX*a_fm
C
C
C-----OBTAIN FLUX FOR COMPONENT K IN PHASE NP.
      DO 4 K=1,NK
C
      XNPMK=WM1*PAR(N1L2NP+NB+K)+WM2*PAR(N2L2NP+NB+K)
C
    4 F(K,M)=F(K,M)+XNPMK*FNPM
C
      goto 32
c
   31 CONTINUE
C
C-----NOW A SECTION FOR COMPUTING BINARY DIFFUSION, DRIVEN BY VAPOR
C     OR GAS DENSITY GRADIENTS.
C
C-----END OF BINARY DIFFUSION SECTION.
C
   32 CONTINUE
C
      IF(MOP(3).LT.7) GOTO 190
      IF(MOP(3).LT.9.AND.M.GT.1) GOTO 190
      WRITE (34,196) NP,N1L2NP,N2L2NP,PER1,PER2,PERI
  196 FORMAT(9X,'NP =',I2,'     SECONDARY INDICES (',I5,1H,,I5,1H),
     X15X,'    PER1 = ',E12.6,' PER2= ',E12.6,'   PERI = ',E12.6)
      WRITE (34,195) S1,REL1,VIS1,RHO1,H1,PCAP1
  195 FORMAT('       S1   = ',E12.6,'  REL1 = ',E12.6,'  VIS1 = ',E12.6,
     A'  RHO1 = ',E12.6,'  H1 = ',E12.6,'  PCAP1 = ',E12.6)
      WRITE (34,194) S2,REL2,VIS2,RHO2,H2,PCAP2
  194 FORMAT('       S2   = ',E12.6,'  REL2 = ',E12.6,'  VIS2 = ',E12.6,
     A'  RHO2 = ',E12.6,'  H2 = ',E12.6,'  PCAP2 = ',E12.6)
      WRITE (34,193) DMOBI,RHOX,DR,FNPM
  193 FORMAT('       DMOBI= ',E12.6,'  RHOX = ',E12.6,'    DR = ',E12.6,
     X'  FNPM = ',E12.6)
C
      WRITE (34,192) (F(K,M),K=1,NK1)
  192 FORMAT('     * FLOW TERMS',8(1X,E12.6))
  190 CONTINUE
C
      IF(M.NE.1) GOTO 3
C-----STORE FLUXES IN EACH PHASE AND HEAT FLUX.
      FLO((N-1)*NPH+NP)=FNPM
C
    3 CONTINUE
    2 CONTINUE
C
C+++++ASSIGN ALL INTERFACE TERMS+++++++++++++++++++++++++++++++++
C
      DO 5 K=1,NEQ
C     K IS THE ROW INDEX WITHIN A BLOCK PERTAINING TO ELEMENT N1 OR N2.
C
C-----COMPUTE FLUX CONTRIBUTIONS TO RESIDUALS.
      IF(N1.LE.NELA)
     XR(N1LOC+K)=R(N1LOC+K)-FAC1*F(K,1)
      IF(N2.LE.NELA)
     XR(N2LOC+K)=R(N2LOC+K)+FAC2*F(K,1)
C
      DO 6 L=1,NEQ
C     L IS THE COLUMN INDEX WITHIN A BLOCK.
      IF(N1.GT.NELA .OR. N2.GT.NELA) GOTO 61
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N2.
      IRN(NZ+1)=N1LOC+K
      IF(IAB.EQ.0) ICN(NZ+1)=N2LOC+L
      IF(IAB.NE.0) JVECT(NZ+1)=N2LOC+L
      CO(NZ+1)=FAC1*(F(K,L+1+NEQ)-F(K,1))/DELX(N2LOCP+L)
      NZ=NZ+1
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N1.
      IRN(NZ+1)=N2LOC+K
      IF(IAB.EQ.0) ICN(NZ+1)=N1LOC+L
      IF(IAB.NE.0) JVECT(NZ+1)=N1LOC+L
      CO(NZ+1)=-FAC2*(F(K,L+1)-F(K,1))/DELX(N1LOCP+L)
      NZ=NZ+1
   61 CONTINUE
C
C-----DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N1.
C
C     NOTE THAT EACH CONNECTION INVOLVING N1 WILL GENERATE A TERM
C     AT THE SAME MATRIX LOCATION.
      N1KL=(N1-1)*NEQ*NEQ+(K-1)*NEQ+L
      IF(N1.LE.NELA)
     XCO(N1KL)=CO(N1KL)+FAC1*(F(K,L+1)-F(K,1))/DELX(N1LOCP+L)
C
C-----DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FROM
C     DEPENDENCE OF COMPONENT K-FLUX UPON VARIABLE L IN ELEMENT N2.
      N2KL=(N2-1)*NEQ*NEQ+(K-1)*NEQ+L
      IF(N2.LE.NELA)
     XCO(N2KL)=CO(N2KL)-FAC2*(F(K,L+1+NEQ)-F(K,1))/DELX(N2LOCP+L)
C
    6 CONTINUE
    5 CONTINUE
C
C+++++END OF ASSIGNMENT OF INTERFACE TERMS+++++++++++++++++++++++
C
    1 CONTINUE
C
C-----TEST FOR CONVERGENCE----------------------------------------------
C
      RERM=0.d0
      IF(MOP(3).GE.2) WRITE (34,203)
  203 FORMAT(/' ===== RESIDUALS ===== MASS BALANCES FIRST,',
     X' ENERGY BALANCE LAST'/)
      DO 10 N=1,NELA
      NLOC=(N-1)*NEQ
C
      IF(MOP(3).GE.2) WRITE (34,202) ELEM(N),(R(NLOC+K),K=1,NEQ)
  202 FORMAT('       AT ELEMENT *',A5,'*   ',8(1X,E12.6))
C
      DO10 K=1,NEQ
      NLM=NLOC+K
      DOA=ABS(DOLD(NLM))
      IF(DOA.LT.RE2) RER=R(NLM)/RE2
      IF(DOA.GE.RE2) RER=R(NLM)/DOLD(NLM)
      IF(ABS(RER).LE.RERM) GOTO 10
      RERM=ABS(RER)
      NER=N
      KER=K
   10 CONTINUE
      RETURN
      END
c
