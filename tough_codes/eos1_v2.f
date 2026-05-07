c-------- E. Sonnenthal 3/20/09
c      - Fixed leverett scaling so that minimum does not go below that calculated at SL=0
c      - Removed all YMP-specific changes and calls
c
C
      BLOCK DATA EQOS
C          ......................................................
C          .                                                    .
C          .  TOUGH2, MODULE EOS1, VERSION 1.0, SEPTEMBER 1990  .
C          ......................................................
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT integer*8 (I-N)
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      common/ech/eosn(20)
      character*10 eosn
      data eosn/'EOS1      ','WATER     ','WATER(2)  ','HEAT      ',
     x'          ','          ','          ','          ',
     x'          ','          ',
     x'          ','          ','          ','          ','          ',
     x'          ','          ','          ','          ','          '/
c
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8
      DATA NK,NEQ,NPH,NB/1,2,2,6/
C
C----- NK IS THE NUMBER OF COMPONENTS.
C----- NEQ IS THE NUMBER OF EQUATIONS PER GRID BLOCK.
C      USUALLY WE HAVE NEQ = NK+1, FOR NK MASS- AND ONE ENERGY-BALANCE.
C----- NPH IS THE NUMBER OF PHASES WHICH CAN BE PRESENT.
C----- NB IS THE NUMBER OF SECONDARY PARAMETERS OTHER THAN MASS
C      FRACTIONS.
C      THE TOTAL NUMBER OF SECONDARY PARAMETERS IS NBK = NB+NK.
C
      END



      SUBROUTINE EOS
C
C*****SINGLE OR TWO COMPONENT WATER/STEAM***********************
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT integer*8 (I-N)
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
C
      COMMON/P1/X((MNK+1)*MNEL)               ! 每个网格最新收敛的主要变量值
      COMMON/P2/DX((MNK+1)*MNEL)              ! 最新的主要变量的增量
      COMMON/P3/DELX((MNK+1)*MNEL)            ! 数值微分过程中的主要变量的增量小量
      COMMON/E1/ELEM(MNEL)                    ! 网格名称
      COMMON/E2/MATX(MNEL)                    ! 每个网格的岩性名称
      COMMON/E4/PHI(MNEL)                     ! 每个网格的孔隙度
      COMMON/E5/P(MNEL)                       ! 每个网格的压力
      COMMON/E6/T(MNEL)                       ! 每个网格的温度
      common/E7/pm(mnel)                      ! 渗透率
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/KC/KC
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/TITLE/TITLE
      CHARACTER*80 TITLE
      COMMON/FAIL/IHALVE
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +                CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/SOLI/COM(MAXMAT),EXPAN(MAXMAT),CDRY(MAXMAT),
     +            TORT(MAXMAT),GK(MAXMAT)
      COMMON/RPCAP/IRP(MAXMAT),RP(7,MAXMAT),ICP(MAXMAT),
     +             CP(7,MAXMAT),IRPD,RPD(7),ICPD,CPD(7)
      COMMON/BIND/DIFF0,TEXP,BE
      COMMON/VPL/NOVPL,VPLF
      COMMON/GASLAW/R,AMS,AMA,CVA
      COMMON/BC/NELA
      common/ff/h1
      character*1 h1
      CHARACTER*5 ELEM,MAT
C
C-------------------------------- for coupling with reactive transport
C
      COMMON/MOP_REACT/MOPR(20)  ! controlling parameters for reactive transport
      double precision pcfct,onemsx,pczero      
c
C----------------------------------Begin addition for solute transport
      COMMON/PARNP/NPL,NPG          ! specify in EOS module
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
C       1: EOS1    
C       2: EOS2    
C       3: EOS3       
C       4: EOS4       
C       7: EOS7       
C       9: EOS9    
C      12: EWASG    
C      13: ECO2    
C--------------------------------------------------------------------
C
      DIMENSION XX(10),DP(22)
C
      SAVE ICALL,ZERO,DAMB
      DATA ZERO,DAMB/1.D-6,998.203D0/
      DATA ICALL/0/
C     
C*******************************************Addition for reactive transport
C------------------Phase code
c
      NPL=2                    ! for solute transport
      NPG=1                    !
c
c------------------EOS module indicator   EOS模块类别指示器
c
      IEOS=1        ! EOS1        
c     
C     1: EOS1    
C     2: EOS2    
C     3: EOS3       
C     4: EOS4       
C     7: EOS7       
C     9: EOS9    
C    12: EWASG    
C    13: ECO2    
C    14: ECO2N
C    16: TMGAS   
C
C*************************************************************************
C
      ICALL=ICALL+1
c
C*************************************************************************
C   
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' EOS 1.1, 2009.3.20: *EOS1* ... THERMOPHYSICAL'
     x' PROPERTIES MODULE FOR WATER')
C
      IF(MOP(5).GE.4) WRITE (34,32) KCYC,ITER
   32 FORMAT(/5H EOS ,16H [KCYC,ITER] = [,I4,1H,,I3,1H]/)
C
      XX(3)=0.D0
C
      IF(KC.NE.0) GOTO 40
      h1=char(12)
C     INITIALIZE SOME DATA
      NOVPL=1
      VPLF=1.D0
      R=8314.56D0
      AMS=18.016D0
      AMA=18.016D0
      CVA=0.D0
C
C-----COME HERE TO ASSIGN DEFAULT RELATIVE PERMEABILITY AND CAPILLARY
C     PRESSURE PARAMETERS TO DOMAINS WITHOUT SPECIAL ASSIGNMENTS.
      DO 33 N=1,NM
      IF(IRP(N).NE.0) GOTO 38
      IRP(N)=IRPD
      ICP(N)=ICPD
      DO 34 M=1,7
      RP(M,N)=RPD(M)
   34 CP(M,N)=CPD(M)
   38 CONTINUE
C
C.....THE FOLLOWING LINES UP TO "33 CONTINUE" INCLUDE SOME CODING
C     FOR HANDLING VAPOR PRESSURE LOWERING; HOWEVER, THIS IS NOT
C     PRESENTLY IMPLEMENTED IN THE REST OF THE CODE, AND ONLY THE
C     OPTION NOVPL = 1 (NO VPL) IS USEABLE.
      IF(NOVPL.EQ.0) GOTO 138
C-----COME HERE WHEN NO VAPOR PRESSURE LOWERING IS DESIRED.
      GOTO 33
C
  138 CONTINUE
      TX=20.D0
      TK=293.15D0
      SL=ZERO
C
      CALL PCAP(SL,TX,PC,N)
C     
      FV=EXP(PC*AMS/(R*TK*DAMB))
      WRITE (34,39) N,MAT(N),FV
   39 FORMAT(5X,I2,1X,A5,10X,E12.5)
   33 CONTINUE
C
C-----GENERATE SOME PRINTOUT CONCERNING THE EQUATION-OF-STATE PACKAGE.
      WRITE (34,300) NK,NEQ,NPH,NB
  300 FORMAT(' EOS1: EQUATION OF STATE FOR WATER (OPTIONAL: TWO-WATER'
     X' MIXTURES)'//' OPTIONS SELECTED ARE: (NK,NEQ,NPH,NB) = (',
     X3(I1,','),I1,')'/' WHERE NK - NUMBER OF FLUID COMPONENTS, NEQ -'
     x' NUMBER OF EQUATIONS PER GRID BLOCK, NPH - NUMBER OF PHASES'
     x' THAT CAN BE PRESENT, NB - NUMBER OF SECONDARY PARAMETERS'
     x' (OTHER THAN COMPONENT MASS FRACTIONS)'//' disable diffusion'
     x' when NK = 1 or NB = 6; enable when NK = 2 and NB = 8')
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7----*----8
      WRITE (34,303)
  303 FORMAT(/' AVAILABLE OPTIONS ARE (NK,NEQ,NPH,NB) = '/' (1,2,2,6'
     X' or 8) - ONE WATER COMPONENT, NON-ISOTHERMAL (DEFAULT:1,2,2,6)'/
     X' (1,1,2,6 or 8) - SINGLE PHASE(LIQUID OR VAPOR), ISOTHERMAL'/
     X' (2,3,2,6 or 8) - TWO WATERS (NON-ISOTHERMAL ONLY)'/)
      WRITE (34,304)
  304 FORMAT(' POSSIBLE PRIMARY VARIABLES ARE:'/' P - PRESSURE, T -'
     X' TEMPERATURE, S - GAS PHASE SATURATION'/' PRIMARY VARIABLES'
     X' COMBINITIONS IN DIFFERENT PHASE CONDITIONS: SINGLE-PHASE: P,T;'
     X' TWO-PHASE: P,S'/)
      IF(NK.EQ.2.AND.NEQ.EQ.3) WRITE(34,301)
  301 FORMAT('FOR TWO-WATER CONDITIONS (NK = 2), X2 - MASS FRACTION OF'
     X' WATER "2" SHOULD ALSO BE SERVED AS A PRIMARY VARIABLE')
C
      IF(NOVPL.EQ.0) WRITE (34,22)
   22 FORMAT(55H MAXIMUM VAPOR PRESSURE LOWERING AT AMBIENT TEMPERATURE,
     A12H OF 20 DEG C//42H       DOMAIN        VAPOR PRESSURE FACTOR)
C
C-----COME HERE TO CONVERT INITIAL CONDITIONS FROM T,S TO P,S.
C
      DO 41 N=1,NEL        
      NLOC=(N-1)*NK1                  
      NLOC2=(N-1)*NEQ1*NSEC           
      IF(X(NLOC+1).GT.374.15D0) GOTO 43
      TX=X(NLOC+1) ! THE FIRST PRIMARY VARIABLE IS TEMPERATURE
      CALL SAT(TX,PS)
      X(NLOC+1)=PS
C
   43 CONTINUE ! NOW THE FIRST PRIMARY VARIABLE IS PRESSURE 
      PX=X(NLOC+1)                    
      IF(X(NLOC+2).GT.1.5D0) GOTO 44
      
C-----COME HERE FOR TWO-PHASE POINTS, AND INITIALIZE GAS SATURATION
C     FOR PHASE CHOICE.
      PAR(NLOC2+1)=X(NLOC+2)
      GOTO 41
C
   44 CONTINUE
C-----COME HERE FOR SINGLE PHASE POINTS.
      TX=X(NLOC+2)
      CALL SAT(TX,PS)
      IF(IGOOD.NE.0) GOTO 200
      PAR(NLOC2+1)=0.D0 ! LIQUID ONLY
      IF(PX.GT.PS) GOTO 41
      PAR(NLOC2+1)=1.D0 ! GAS ONLY
C      
   41 CONTINUE
   40 CONTINUE
C
      IF(MOP(5).GE.7) WRITE (34,31)
   31 FORMAT(/' PRIMARY VARIABLES')
C
      N1=NEL
      IF(KC.GT.0) N1=NELA
      
      DO 1 N=1,N1
C
      NMAT=MATX(N)
c
        phi1=phi(n)
cels6/8/09 to keep consistent w/ TOUGH2 V2   tort0=tort(nmat)
        tort0=dabs(tort(nmat))
c
cels10/14/09
c Set factor for capillary pressure scaling

        if(mopr(6).eq.1)then
           pcfct = pcfact(n)
        else
           pcfct = 1.d0
        endif
      NLOC=(N-1)*NK1
      NLOC2=(N-1)*NSEC*NEQ1
C
C-----SET PRIMARY VARIABLES.
      DO 2 M=1,NK1
      XINCR=0.D0                  ! small increment of primary variables
      IF(ITER.NE.0.AND.KON.NE.2) XINCR=DX(NLOC+M)
      XX(M)=X(NLOC+M)+XINCR
    2 CONTINUE
C
      IF(MOP(5).GE.7) WRITE (34,35) ELEM(N),(XX(M),M=1,NK1)
   35 FORMAT(/9H ELEMENT ,A5/(10(1X,E12.5)))
C
C-----CHECK WHETHER PHASE TRANSITIONS OCCUR FOR LATEST UPDATED
C     VARIABLES.
C
      PX=XX(1)
      if(px.le.0.d0) goto 200
      XTS=XX(2)
C
C..... MAKE CHOICE OF PHASE CONDITIONS. ................................
C
      IF(IHALVE.EQ.0) GOTO 136
C-----COME HERE FOR PHASE CHOICE AFTER TIME STEP REDUCTION FOR
C     LACK OF CONVERGENCE.
      IF(XX(2).LT.1.D0) GOTO 10       
      
C-----COME HERE FOR SINGLE-PHASE POINTS.
c.....4-25-95: rename PSAT --> PS
      CALL SAT(XX(2),ps)
      IF(IGOOD.NE.0) GOTO 200 ! TEMPERATURE OUT OF RANGE OF SATUATION CURVE
      IF(XX(1).GT.ps) GOTO 8
      GOTO 9
C
  136 CONTINUE
C   
      SG=PAR(NLOC2+1)                 
      IF(SG.EQ.0.D0.OR.SG.EQ.1.D0) GOTO 6
C
C-----COME HERE FOR TWO-PHASE POINTS.
      IF(0.D0.LE.XTS.AND.XTS.LE.1.D0) GOTO 10
C
C-----COME HERE FOR PHASE TRANSITION TO SINGLE PHASE.
      IF(ITER.EQ.0) GOTO 20
C
C-----TAKE TEMPERATURE FROM LATEST ITERATION.
      TX=PAR(NLOC2+NSEC-1)
      DX(NLOC+2)=TX-X(NLOC+2) ! TEMPERATURE INCREMENT
      XX(2)=TX
C
C-----TAKE PRESSURE TO BE SOMEWHAT BELOW OR ABOVE SATURATION PRESSURE.
      CALL SAT(TX,PS)
      IF(IGOOD.NE.0) GOTO 200
C
      IF(MOP(5).GE.3) WRITE (34,133) ELEM(N),TX,PX,PS
  133 FORMAT(53H .....PHASE TRANSITION TO SINGLE PHASE.. AT ELEMENT *,A5
     A,8H*   T = ,E12.5,7H   P = ,E12.5,8H   PS = ,E12.5)
C
      IF(XTS.LT.0.D0) GOTO 18
C
C-----COME HERE FOR PHASE TRANSITION TO STEAM.
      DX(NLOC+1)=.999999D0*PS-X(NLOC+1)
      XX(1)=.999999D0*PS
      GOTO 9
C
   18 CONTINUE
C-----COME HERE FOR PHASE TRANSITION TO LIQUID.
      DX(NLOC+1)=1.000001D0*PS-X(NLOC+1)
      XX(1)=1.000001D0*PS
      GOTO 8
C
C
    6 CONTINUE
C
C-----COME HERE FOR SINGLE-PHASE POINTS.
C
      CALL SAT(XTS,PS)
      IF(IGOOD.NE.0) GOTO 200
C
      IF(ITER.NE.0) GOTO 7
C-----COME HERE TO INITIALIZE SINGLE-PHASE POINTS.
      IF(PX.GT.PS) GOTO 8
      GOTO 9
    7 CONTINUE
C
      SG=PAR(NLOC2+1)
      IF(SG.EQ.0.D0.AND.PX.GT.PS) GOTO 8
      IF(SG.EQ.1.D0.AND.PX.LT.PS) GOTO 9
C
C-----COME HERE FOR PHASE TRANSITION TO TWO-PHASE.
C
C-----TAKE PRESSURE TO BE EQUAL TO SATURATION PRESSURE AT LATEST
C     TEMPERATURE.
      DX(NLOC+1)=PS-X(NLOC+1)
      XX(1)=PS
C
      IF(MOP(5).GE.3) WRITE (34,134) ELEM(N),XTS,PX,PS
  134 FORMAT(53H .....PHASE TRANSITION TO TWO-PHASE..... AT ELEMENT *,A5
     A,8H*   T = ,E12.5,7H   P = ,E12.5,8H   PS = ,E12.5)
C
      IF(SG.EQ.0.D0) GOTO 19
C
C-----COME HERE FOR PHASE TRANSITION FROM STEAM TO TWO-PHASE.
      tmpomz = 1.D0-ZERO
      DX(NLOC+2)=tmpomz-X(NLOC+2)
      XX(2)=tmpomz
      GOTO 10
C
   19 CONTINUE
C
C-----COME HERE FOR PHASE TRANSITION FROM LIQUID TO TWO-PHASE.
      DX(NLOC+2)=ZERO-X(NLOC+2)
      XX(2)=ZERO
      GOTO 10
C
   20 CONTINUE
C-----COME HERE IF PHASE TRANSITION IS FOUND AT BEGINNING OF TIME
C     STEP (ITER=0) -- THIS SHOULD NEVER HAPPEN +++++++++++++++++
      WRITE (34,221) KCYC,ITER,ELEM(N),PX,XTS
  221 FORMAT(28H PHONY PHASE TRANSITION AT [,I4,1H,,I3,10H]  -----  ,9H
     AELEMENT ,A5,6H   P =,E12.5,7H  XTS =,E12.5)
      IGOOD=2
      RETURN
C
C*****NOW COMPUTE ALL SECONDARY VARIABLES.
C
    8 CONTINUE
C-----LIQUID-----LIQUID-----LIQUID-----LIQUID-----LIQUID
C
      DO 100 K=1,NEQ1
      NLK2=NLOC2+(K-1)*NSEC
      NLK2L=NLK2+NBK
C
      PX=XX(1)
      TX=XX(2)
      X2M=XX(3)
C
      IF(K.EQ.1) GOTO 101
      DELX(NLOC+K-1)=DFAC*XX(K-1)+1.D-10
      IF(K.EQ.2) PX=XX(1)+DELX(NLOC+1)
      IF(K.EQ.3) TX=XX(2)+DELX(NLOC+2)
      IF(K.EQ.4) DELX(NLOC+3)=DFAC
      IF(K.EQ.4) X2M=XX(3)+DELX(NLOC+3)
C
  101 CONTINUE
      X1M=1.D0-X2M
C
      CALL VISW(TX,PX,PS,VW) ! VISCOSITY
      CALL COWAT(TX,PX,D,U) ! LIQUID WATER DENSITY AND INTERNAL ENERGY
      IF(IGOOD.NE.0) GOTO 200
      CALL SUPST(TX,PS,DS,US) ! VAPOR DENSITY AND INT. ENERGY
C
      SL=1.D0
C
      CALL PCAP(SL,TX,PC,NMAT)
C
      PAR(NLK2L+1)=1.D0
      PAR(NLK2L+2)=1.D0
      PAR(NLK2L+3)=VW
      PAR(NLK2L+4)=D
      PAR(NLK2L+5)=U+PX/D
      PAR(NLK2L+6)=PC
      if(nk.ne.1.and.nb.ge.8) then
c         effective diffusivities.
             if(tort0.eq.0.d0) then
c.....Millington-Quirk
                toto=phi1**(1.d0/3.d0)
             else
                toto=tort0
             endif
          par(nlk2l+7)=phi1*toto*d
          par(nlk2l+8)=1.d0
      endif
c
      PAR(NLK2L+nb+1)=X1M
C
      PAR(NLK2+1)=0.D0
      PAR(NLK2+2)=0.D0
      PAR(NLK2+3)=1.D0
      PAR(NLK2+4)=DS
      PAR(NLK2+5)=0.D0
      PAR(NLK2+6)=0.D0
      if(nk.ne.1.and.nb.ge.8) then
c         effective diffusivities.
          par(nlk2+7)=0.d0
          par(nlk2+8)=1.d0
      endif
c
      PAR(NLK2+nb+1)=X1M
C
      PAR(NLK2+NSEC-1)=TX
      PAR(NLK2+NSEC)=0.D0
C
      IF(KC.EQ.0.AND.K.EQ.1) GOTO 118
      GOTO 102
C-----INITIALIZE TEMPERATURE AND PRESSURE.
  118 T(N)=TX
      P(N)=PX
  102 CONTINUE
      IF(NK.EQ.1) GOTO 100
      PAR(NLK2L+nb+2)=X2M
      PAR(NLK2+nb+2)=X2M
C
  100 CONTINUE
C
      GOTO 1
C

      
      
    9 CONTINUE

C-----STEAM-----STEAM-----STEAM-----STEAM-----STEAM-----STEAM
C
      DO 110 K=1,NEQ1
      NLK2=NLOC2+(K-1)*NSEC
      NLK2L=NLK2+NBK
C
      PX=XX(1)
      TX=XX(2)
      X2M=XX(3)
C
      IF(K.EQ.1) GOTO 111
      DELX(NLOC+K-1)=DFAC*XX(K-1)+1.D-10
      IF(K.EQ.2) PX=XX(1)+DELX(NLOC+1)
      IF(K.EQ.3) TX=XX(2)+DELX(NLOC+2)
      IF(K.EQ.4) DELX(NLOC+3)=DFAC
      IF(K.EQ.4) X2M=XX(3)+DELX(NLOC+3)
C
  111 CONTINUE
      X1M=1.D0-X2M
C
      CALL SUPST(TX,PX,D,U)
      CALL COWAT(TX,PS,DW,UW)
      IF(IGOOD.NE.0) DW=0.d0
      IGOOD=0
      CALL VISS(TX,PX,D,VS)           ! 用气体密度,温度和压力计算气态水的粘度
C
      SL=0.D0
       CALL PCAP(SL,TX,PC,NMAT)       ! 用水的饱和度和网格的温度计算了网格内的毛细压力
C
      SG=1.D0
      KX=K
      CALL RELP(SG,REPW,REPS,NMAT,KX,NLOC,SG) ! 计算了气体和液体的相对渗透率
C
      PAR(NLK2+1)=1.D0                ! 次级变量第一位：气相饱和度
      PAR(NLK2+2)=REPS                ! 次级变量第二位：气相的相对渗透率 
      PAR(NLK2+3)=VS                  ! 次级变量第三位：气相的粘度
      PAR(NLK2+4)=D                   ! 次级变量第四位：气相的比焓
      PAR(NLK2+5)=U+PX/D              ! 次级变量第五位：气相的毛细压力
      PAR(NLK2+6)=0.D0                ! 次级变量第六位：气相的扩散系数
      
      if(nk.ne.1.and.nb.ge.8) then
c         effective diffusivities for gas are a function of
c         pressure and temperature.
          peffect=1.0d5/px
          teffect=1.0d0
          if(texp.ne.0.d0) teffect=((tx+273.15d0)/273.15d0)**texp
             if(tort0.eq.0.d0) then
c.....Millington-Quirk
                toto=phi1**(1.d0/3.d0)
             else
                toto=tort0
             endif
          par(nlk2+7)=phi1*toto*d
          par(nlk2+8)=peffect*teffect
      end if
             
      PAR(NLK2+nb+1)=X1M              ! 液相的饱和度
      PAR(NLK2L+1)=0.D0               ! 液相的相对渗透率
      PAR(NLK2L+2)=0.D0               ! 液相的粘度
      PAR(NLK2L+3)=1.D0               ! 液相的密度
      PAR(NLK2L+4)=DW                 ! 气相的比焓
      PAR(NLK2L+5)=0.D0               ! 气相的毛细压力
      PAR(NLK2L+6)=PC                 ! 气相的扩散系数
     
      if(nk.ne.1.and.nb.ge.8) then
c         effective diffusivities.
          par(nlk2l+7)=0.d0
          par(nlk2l+8)=1.d0
      endif
C-----如果次级变量数等于8则激活扩散的赋值
      

      PAR(NLK2L+nb+1)=X1M             ! 气态水的质量分数
C
      PAR(NLK2+NSEC-1)=TX             ! 温度
      PAR(NLK2+NSEC)=U+PX/D           ! 其他,这里存的是气相的毛细压力
C
      IF(KC.EQ.0.AND.K.EQ.1) GOTO 117 ! 要搞清楚:KC等于各种值的时候都是什么条件,K等于各种值的时候是什么条件
      GOTO 112
C-----INITIALIZE TEMPERATURE AND PRESSURE.
  117 T(N)=TX
      P(N)=PX
  112 CONTINUE
      IF(NK.EQ.1) GOTO 110            ! 如果相数不等于1(就是等于2),那么就需要保存两相的质量分数
      PAR(NLK2L+nb+2)=X2M             ! X2M=XX(3),另一种组分的质量分数
      PAR(NLK2+nb+2)=X2M
C
  110 CONTINUE
C
      GOTO 1
C
C
   10 CONTINUE
C*****TWO-PHASE*****TWO-PHASE*****TWO-PHASE*****TWO-PHASE
C
      DO 120 K=1,NEQ1
      NLK2=NLOC2+(K-1)*NSEC
      NLK2L=NLK2+NBK
C
      PX=XX(1)
      SX=XX(2)
      X2M=XX(3)
      IF(KC.EQ.0) TX0=0.D0
      IF(KC.NE.0) TX0=T(N)
C
      IF(K.EQ.1) GOTO 121
      DELX(NLOC+K-1)=DFAC*XX(K-1)+1.D-10
      IF(K.EQ.2) PX=XX(1)+DELX(NLOC+1)
      IF(K.EQ.3) DELX(NLOC+2)=DFAC
      IF(K.EQ.3) SX=XX(2)+DELX(NLOC+2)
      IF(K.EQ.4) DELX(NLOC+3)=DFAC
      IF(K.EQ.4) X2M=XX(3)+DELX(NLOC+3)
C
  121 CONTINUE
      X1M=1.D0-X2M
C
      KX=K
      CALL RELP(SX,REPW,REPS,NMAT,KX,NLOC,XX(2))
      CALL TSAT(PX,TX0,TX)
      IF(IGOOD.NE.0) GOTO 200
C
      SL=1.D0-SX
C
      CALL PCAP(SL,TX,PC,NMAT)
c
c... Modify pcap only if permeability-porosity changes are required
           pc = pc*pcfct
c
C----------------------------------------------------------
C
      CALL COWAT(TX,PX,DW,UW)
      IF(IGOOD.NE.0) GOTO 200
      CALL SUPST(TX,PX,DS,US)
      CALL VIS(TX,PX,DS,VW,VS,PX)
C
      PAR(NLK2+1)=SX
      PAR(NLK2+2)=REPS
      PAR(NLK2+3)=VS
      PAR(NLK2+4)=DS
      PAR(NLK2+5)=US+PX/DS
      PAR(NLK2+6)=0.D0
      if(nk.ne.1.and.nb.ge.8) then
c         effective diffusivities for gas are a function of
c         pressure and temperature.
          peffect=1.0d5/px
          teffect=1.0d0
          if(texp.ne.0.d0) teffect=((tx+273.15d0)/273.15d0)**texp
             if(tort0.eq.0.d0) then
c.....Millington-Quirk
                toto=(phi1*sx)**(1.d0/3.d0)*sx*sx*sx
             else
                toto=tort0*reps
             endif
          par(nlk2+7)=phi1*toto*ds
          par(nlk2+8)=peffect*teffect
      endif
c
      PAR(NLK2+nb+1)=X1M
C
      PAR(NLK2L+1)=SL
      PAR(NLK2L+2)=REPW
      PAR(NLK2L+3)=VW
      PAR(NLK2L+4)=DW
      PAR(NLK2L+5)=UW+PX/DW
      PAR(NLK2L+6)=PC
      if(nk.ne.1.and.nb.ge.8) then
c         effective diffusivities.
             if(tort0.eq.0.d0) then
c.....Millington-Quirk
                toto=(phi1*sl)**(1.d0/3.d0)*sl*sl*sl
             else
                toto=tort0*repw
             endif
          par(nlk2l+7)=phi1*toto*dw
          par(nlk2l+8)=1.d0
      endif
c
      PAR(NLK2L+nb+1)=X1M
C
      PAR(NLK2+NSEC-1)=TX
      PAR(NLK2+NSEC)=US+PX/DS
C
      IF(KC.EQ.0.AND.K.EQ.1) GOTO 119
      GOTO 122
C-----INITIALIZE TEMPERATURE AND PRESSURE.
  119 T(N)=TX
      P(N)=PX
  122 CONTINUE
      IF(NK.EQ.1) GOTO 120
      PAR(NLK2L+nb+2)=X2M
      PAR(NLK2+nb+2)=X2M
C
  120 CONTINUE
C
    1 CONTINUE
C
C
C
C
      IF(MOP(5).LT.8) GOTO 29
      WRITE (34,27)
   27 FORMAT(/21H SECONDARY PARAMETERS)
      DO 28 N=1,NEL
      NLOC2=(N-1)*NSEC*NEQ1
C-----PRINT PARAMETERS AT STATE POINT.
      WRITE (34,30) ELEM(N),(PAR(NLOC2+M),M=1,NSEC)
C
      IF(MOP(5).EQ.9) GOTO 26
C-----COME HERE TO PRINT INCREMENTED PARAMETERS.
      DO 25 K=2,NEQ1
   25 WRITE (34,37) (PAR(NLOC2+(K-1)*NSEC+M),M=1,NSEC)
   30 FORMAT(/9H ELEMENT ,A5/(10(1X,E12.5)))
   37 FORMAT(/(10(1X,E12.5)))
      GOTO 28
C
C-----COME HERE TO PRINT DERIVATIVES.
   26 CONTINUE
      NLOC=(N-1)*NK1
      DO 24 K=2,NEQ1
      DO 23 M=1,NSEC
   23 DP(M)=(PAR(NLOC2+(K-1)*NSEC+M)-PAR(NLOC2+M))/DELX(NLOC+K-1)
      WRITE (34,37) (DP(M),M=1,NSEC)
   24 CONTINUE
   28 CONTINUE
   29 CONTINUE
C
      RETURN
C
  200 WRITE (34,201) ELEM(N),(XX(M),M=1,NK1)
  201 FORMAT(' EOS CANNOT FIND PARAMETERS AT ELEMENT *',A5,'* XX(M) =' 
     A,4(E12.5,1X))
      RETURN
      END
c
c
      SUBROUTINE OUT
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT integer*8 (I-N)
C
C-----THIS SUBROUTINE GENERATES PRINTOUT.
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
C
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)
      COMMON/XYZ11/XXX(mnel)     !!!! mesh coordinates for TECPLOT
      COMMON/XYZ22/YYY(mnel)
      COMMON/XYZ33/ZZZ(mnel)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR PRIMARY VARIABLES $$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NEQ = 4*NEL
C
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P3/DELX((MNK+1)*MNEL)
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
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
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
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      COMMON/COMPO/FLO(MNPH*MNCON)
      COMMON/CHANGE/DEMAX,DDMAX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/ TITLE /   TITLE
      CHARACTER*80 TITLE
      COMMON/  DOP  /ENTH,  KDATA,QUAL
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/  POV6 / TSTART
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +                CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      common/ff/h1
      character*1 h1
      DIMENSION DXM(10)
      CHARACTER*1 HB,H0,ITABG
cels5/1/08      CHARACTER*5 ELEM,ELEM1,ELEM2,ELEG,SOURCE,MAT,type,typen
      CHARACTER*5 ELEM,ELEM1,ELEM2,ELEG,SOURCE,MAT,typen
      common/source_type/type5(mnogn)
      character*5 type5
C
      SAVE ICALL,HB,H0,i57,i58
      DATA HB,H0/' ',' '/
      DATA ICALL,i57,i58/0,57,58/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *******OUT 1.01, 1999.9.2: PRINT RESULTS FOR ELEMENTS,'
     X' CONNECTIONS, AND SINKS/SOURCES')
C
C-----COMPUTE MAXIMUM CHANGES.
C
      DO 1 I=1,NEQ1
    1 DXM(I)=0.D0
C
      DO 2 N=1,NEL
      NLOC=(N-1)*NK1
      DO 3 I=1,NEQ
      ADX=ABS(DX(NLOC+I))
      IF(ADX.GT.DXM(I)) DXM(I)=ADX
    3 CONTINUE
    2 CONTINUE
C
C-----PRINT HEADER INFORMATION.
C
      DAY=SUMTIM/8.64D4
      WRITE (34,5010) H1,TITLE,KCYC,ITER,KON,DAY
      WRITE (34,9010) SUMTIM,KCYC,ITER,ITERC,KON,(DXM(I),I=1,3),
     ARERM,NER,KER,DELTEX
      WRITE (34,1000) H0
C
      WRITE(11,5012) DAY
 5012 FORMAT(1X,'X,Y,Z,P,T'/1X,'AT TIME=',E12.6,'DAYS')
      DO 3021 N=1,NEL
      NLOC=(N-1)*NK1
      NLOC2=(N-1)*NEQ1*NSEC
      NLOC2L=NLOC2+NBK
      X1=X(NLOC+1)
      X2=0.D0
      IF(NK.EQ.2) X2=PAR(NLOC2+nb+2)
      IF(KON.NE.2) X1=X1+DX(NLOC+1)
      IF(MOD(N,i57).EQ.46) WRITE (34,1002) H1
      WRITE(11,5354) XXX(N),YYY(N),ZZZ(N),X(NLOC+1),PAR(NLOC2+NSEC-1)
 5354 FORMAT(3(E10.4,1X),E14.6,1X,F14.6,1X)
      WRITE (34,5040) ELEM(N),N,X(NLOC+1),PAR(NLOC2+NSEC-1),
     APAR(NLOC2+1),PAR(NLOC2L+1),
     BPAR(NLOC2+nb+1),X22,PAR(NLOC2L+6),
     CPAR(NLOC2+4),PAR(NLOC2L+4)
 3021 CONTINUE
C
      i10=10
      IF(MOD(KDATA,i10).LT.2) GOTO 3045
C
C-----NOW PRINT FLOW TERMS.
C
      WRITE (34,5000) h1,TITLE,KCYC,ITER,SUMTIM
      WRITE (34,5060) H0,HB
C
      DO 3030 N=1,NCON
      IF(NEX1(N).EQ.0.OR.NEX2(N).EQ.0) GOTO 3030
C
      N1=NEX1(N)
      N2=NEX2(N)
      NNP=(N-1)*NPH
      FLOF=0.D0
C
      DO 3031 NP=1,NPH
 3031 FLOF=FLOF+FLO(NNP+NP)
      H=0.D0
      IF(FLOF.NE.0.D0) H=GLO(N)/FLOF
C
      IF(MOD(N,i57).EQ.54) WRITE (34,5062) H1,HB
C
      F2=0.D0
      IF(NK.EQ.1) GOTO 3032
      N1LOC2=(N1-1)*NSEC*NEQ1
      N2LOC2=(N2-1)*NSEC*NEQ1
      F2V=FLO(NNP+1)*PAR(N2LOC2+nb+2)
      IF(FLO(NNP+1).LT.0.D0) F2V=FLO(NNP+1)*PAR(N1LOC2+nb+2)
      F2L=FLO(NNP+2)*PAR(N2LOC2+NBK+nb+2)
      IF(FLO(NNP+2).LT.0.D0) F2L=FLO(NNP+2)*PAR(N1LOC2+NBK+nb+2)
      F2=F2V+F2L
 3032 CONTINUE
C
      WRITE (34,5070) ELEM1(N),ELEM2(N),N,GLO(N),H,FLOF,(FLO(NNP+NP)
     A,NP=1,NPH),F2
C
 3030 CONTINUE
C
C*******************************
      if(nk.ne.1.and.nb.ge.8) CALL OUTDF
C*******************************
C
      i10=10
      IF(MOD(KDATA,i10).LT.3) GOTO 3045
C
      WRITE (34,5000) h1,TITLE,KCYC,ITER,SUMTIM
      IF(NK.EQ.2) WRITE (34,5030) H0
      IF(NK.EQ.1) WRITE (34,5031) H0
C
      DO 3020 N = 1, NEL
      NLOC=(N-1)*NK1
      NLOC2=(N-1)*NEQ1*NSEC
      NLOC2L=NLOC2+NBK
      IF(MOD(N,i58).EQ.54.and.nk.eq.2)     WRITE (34,5032) H1
      IF(MOD(N,i58).EQ.54.and.nk.eq.1)     WRITE (34,5033) H1
      PRES=X(NLOC+1)
      IF(KON.NE.2) PRES=PRES+DX(NLOC+1)
      WRITE (34,5050) ELEM(N),N,(X(NLOC+I),I=1,NK1),
     A (DX(NLOC+I),I=1,NK1),PAR(NLOC2+2),PAR(NLOC2L+2)
C
 3020 CONTINUE
C
 3045 CONTINUE
C
      IF(NOGN.EQ.0) RETURN
      WRITE (34,5000) h1,TITLE,KCYC,ITER,SUMTIM
      WRITE (34,5120) H0
C
      GC=0.D0
      HGC=0.D0
C
      DO 3050 N=1,NOGN
cels10/8/09      typen=type(n)
      typen=type5(n)
cels5/1/08      itabg=typen(5:5)
      itabg(n)=typen(5:5)
      J=NEXG(N)
      JLOC2=(J-1)*NSEC*NEQ1
      IF(J.EQ.0) GOTO 3050
      IF(MOD(N,i58).EQ.54) WRITE (34,5122) H1
C
      X1=0.D0
cels5/1/08      IF((LCOM(N).lt.NK1.AND.ITABG.EQ.' ').or.lcom(n).gt.nk1)
      IF((LCOM(N).lt.NK1.AND.ITABG(n).EQ.' ').or.lcom(n).gt.nk1)
     A   X1=FF((N-1)*NPH+1)*PAR(JLOC2+nb+1)+FF((N-1)*NPH+2)
     x   *PAR(JLOC2+NBK+nb+1)
      X2=0.D0
      IF(NK.EQ.2) X2=1.D0-X1
C
      IF(GPO(N).GT.0.D0) then
         x1=1.d0
         if(lcom(n).eq.2) x1=0.d0
         x2=1.d0-x1
      endif
c
      if(lcom(n).eq.nk1) then
         WRITE (34,5130) ELEG(N),SOURCE(N),N,GPO(N)
      else
         if(gpo(n).ge.0.d0) then
            WRITE (34,5130) ELEG(N),SOURCE(N),N,GPO(N),EG(N),X1,X2
         else
               if(lcom(n).lt.nk1) then
      WRITE (34,5130) ELEG(N),SOURCE(N),N,GPO(N),EG(N),X1,X2,
     B(FF((N-1)*NPH+NP),NP=1,NPH)
               else
      WRITE (34,5130) ELEG(N),SOURCE(N),N,G(N),EG(N),X1,X2,
     B(FF((N-1)*NPH+NP),NP=1,NPH),PWB(N)
               endif
         endif
      endif
C
      IF(LCOM(N).NE.NEQ1) GOTO 102
C
C-----NOW FOR A LITTLE SECTION WHICH COMPUTES TOTAL FLOWRATE
C     AND FLOWING ENTHALPY FOR WELLS ON DELIVERABLITY WITH COMPLE-
C     TIONS IN DIFFERENT LAYERS (OR, GENERALLY, ELEMENTS)
C
C     ALL OPEN INTERVALS OF A WELL MUST BE KNOWN BY THE SAME SOURCE
C     NAME, AND MUST BE GIVEN IN UNINTERRUPTED SEQUENCE.
C
      IF(NOGN.EQ.1) GOTO 102
      IF(N.EQ.1.AND.SOURCE(N+1).NE.SOURCE(N)) GOTO 102
      IF(N.EQ.1) GOTO 100
      IF(SOURCE(N-1).EQ.SOURCE(N).OR.SOURCE(N).EQ.SOURCE(N+1)) GOTO 100
C
      GOTO 103
C
C-----COME HERE FOR SOURCE IN A CHAIN AND ACCUMULATE TERMS.
  100 GC=GC+GPO(N)
      HGC=HGC+GPO(N)*EG(N)
C
C-----FIND OUT WHETHER CHAIN TERMINATES.
      IF(N.EQ.NOGN) GOTO 101
      IF(SOURCE(N+1).EQ.SOURCE(N)) GOTO 102
C
C-----COME HERE FOR END OF CHAIN.
  101 IF(GC.NE.0.D0) HGC=HGC/GC
      WRITE (34,110) SOURCE(N),GC,HGC
  110 FORMAT(16H ***** SOURCE  $,A5,12H$     RATE =,E12.5,23H     FLOWIN
     AG ENTHALPY =,E12.5,10H     *****/)
C
C-----COME HERE FOR SOURCES OUTSIDE CHAINS, AND AT END OF CHAIN.
  103 GC=0.D0
      HGC=0.D0
C
  102 CONTINUE
C
C-----END OF SECTION FOR MULTI-LAYER WELLS ON DELIVERABILITY
C
 3050 CONTINUE
      RETURN
C
 5000 FORMAT(A1/A80/,' KCYC =',I7,' ITER =',I5,' TIME =',E12.5/)
 5010 FORMAT(A1/A80/' OUTPUT DATA AFTER (',I4,I3,')',I1,'TIME STEPS'
     A' TIME =',E12.5,'DAYS')
 9010 FORMAT(3X,51HTOTAL TIME   KCYC   ITER  ITERC    KON    DX1M     ,6
     A1H    DX2M         DX3M                      RERM        NER   ,14
     BH KER    DELTEX/1X,E12.5,4I7,3E13.6,13X,E13.6,2I6,E13.6)
 5030 FORMAT(A1,43HELEM.  INDEX    X1           X2          X3,9X,3HDX1,
     A9X,3HDX2,9X,3HDX3,6X,9HKREL(VAP),3X,9HKREL(LIQ)/)
 5031 FORMAT(A1,'ELEM.  INDEX    X1           X2          DX1',
     A'         DX2     KREL(VAP)   KREL(LIQ)'/)
 5032 FORMAT(A1/44H ELEM.  INDEX    X1           X2          X3,9X,
     A3HDX1,9X,3HDX2,9X,3HDX3,6X,9HKREL(VAP),3X,9HKREL(LIQ)/)
 5033 FORMAT(A1/' ELEM.  INDEX    X1           X2          DX1',
     A'         DX2     KREL(VAP)   KREL(LIQ)'/)
 1000 FORMAT(A1,42HELEM.  INDEX     P           T          SG,10X,2HSW,1
     A0X,2HX1,10X,2HX2,7X,4HPCAP,23X,2HDG,10X,2HDW/
     X17X,'(PA)',6X,'(DEG-C)',52X,'(PA)',20X,'(KG/M**3)',3X,
     X'(KG/M**3)'/)
 1002 FORMAT(A1/43H ELEM.  INDEX     P           T          SG,10X,2HSW,
     A10X,2HX1,10X,2HX2,7X,4HPCAP,23X,2HDG,10X,2HDW/
     X17X,'(PA)',6X,'(DEG-C)',52X,'(PA)',20X,'(KG/M**3)',3X,
     X'(KG/M**3)'/)
 5040 FORMAT(1X,A5,I6,7E12.5,12X,2E12.5)
 5050 FORMAT(1X,A5,I6,10E12.5)
 5060 FORMAT( A1,26HELEM1 ELEM2  INDEX    FLOH,A6,9HFLOH/FLOF,7X,4HFLOF
     A,7X,8HFLO(GAS),5X,8HFLO(AQ.),5X,9HFLO(WTR2)/
     X24X,'(W)',8X,'(J/KG)',8X,'(KG/S)',7X,'(KG/S)',7X,'(KG/S)',
     X7X,'(KG/S)'/)
 5062 FORMAT(A1/27H ELEM1 ELEM2  INDEX    FLOH,A6,9HFLOH/FLOF,7X,4HFLOF
     A,7X,8HFLO(GAS),5X,8HFLO(AQ.),5X,9HFLO(WTR2)/
     X24X,'(W)',8X,'(J/KG)',8X,'(KG/S)',7X,'(KG/S)',7X,'(KG/S)',
     X7X,'(KG/S)'/)
 5070 FORMAT(1X,2A6,I6,8E13.6)
 5120 FORMAT(A1,62HELEMENT SOURCE INDEX      GENERATION RATE     ENTHALP
     AY      X1,11X,46HX2          FF(GAS)      FF(AQ.)         P(WB)/
     X29X,'(KG/S) OR (W)',5X,'(J/KG)',62X,'(PA)'/)
 5122 FORMAT(A1/63H ELEMENT SOURCE INDEX      GENERATION RATE     ENTHAL
     APY      X1,11X,46HX2          FF(GAS)      FF(AQ.)         P(WB)/
     X29X,'(KG/S) OR (W)',5X,'(J/KG)',62X,'(PA)'/)
 5130 FORMAT(2X,2A8,I2,10X,E12.5,3X,6(E12.5,1X))
      END
c
      SUBROUTINE BALLA
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT integer*8 (I-N)
C
C-----THIS SUBROUTINE PERFORMS VOLUME- AND MASS-BALANCES.
C
C
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +                CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
      common/ff/h1
      character*1 h1
      CHARACTER*5 MAT
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ******BALLA 1.0, 1990.1.22: PERFORM SUMMARY BALANCES FOR'
     x' VOLUME, MASS AND ENERGY**********')
      WRITE (34,1) KCYC,ITER,SUMTIM
    1 FORMAT(/' VOLUME- AND MASS-BALANCES: IN [KCYC,ITER,TIME'
     X'(SECONDS)] = [',I4,',',I3,',',E12.5,']')
C
      VOLG=0.D0
      VOLW=0.D0
      AMG=0.D0
      AMW=0.D0
      AMWS=0.D0
      AMWTR2=0.D0
      URGW=0.D0
      DMWTR2=0.D0
C
      DO 10 N=1,NELA
C
      NMAT=MATX(N)
      IF(SH(NMAT).GE.1.E4) GOTO 10
C
      NLOC=(N-1)*NK1
      NLOC2=(N-1)*NEQ1*NSEC
C
      PHIV=PHI(N)*EVOL(N)
      DVOLG=PHIV*PAR(NLOC2+1)
      DVOLW=PHIV*PAR(NLOC2+NBK+1)
C
      VOLG=VOLG+DVOLG
      VOLW=VOLW+DVOLW
      DUR=(1.D0-PHI(N))*EVOL(N)*DM(NMAT)*SH(NMAT)*PAR(NLOC2+NSEC-1)
C
      DAMG=DVOLG*PAR(NLOC2+4)         ! 变化的气体的质量
      DUG=DAMG*(PAR(NLOC2+5)-X(NLOC+1)/PAR(NLOC2+4))  ! 变化的气体的比焓
      AMG=AMG+DAMG                    ! 变化后总的气体的焓值
      DAMW=DVOLW*PAR(NLOC2+NBK+4)     ! 变化的液体的质量
      DUW=DAMW*(PAR(NLOC2+NBK+5)-X(NLOC+1)/PAR(NLOC2+NBK+4))          ! 变化的液体的比焓
      URGW=URGW+DUR+DUG+DUW           ! internal energy
      AMW=AMW+DAMW                    ! total aqueous mass
      IF(NK.EQ.2)
     ADMWTR2=DAMG*PAR(NLOC2+NB+2)+DAMW*PAR(NLOC2+NBK+NB+2)
      AMWTR2=AMWTR2+DMWTR2
      AMWS=AMWS+DAMW+DAMG
C
   10 CONTINUE
      WRITE (34,3) VOLG,VOLW
    3 FORMAT(/' PHASE VOLUMES (M^3): GAS ',E11.5,'; AQUEOUS ',E11.5)
      WRITE (34,4) AMG,AMW,AMWS,AMWTR2,URGW
    4 FORMAT(/' MASS (KG): GAS ',E11.5,'; AQUEOUS ',E11.5,' TOTAL H2O '
     X,E11.5,' KG;  H2O(2) ',E11.5,' KG; INT. ENERGY ',E11.5,' J')
      RETURN
      END
c
      SUBROUTINE QLOSS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT integer*8 (I-N)
C
C-----THIS ROUTINE COMPUTES A SEMI-ANALYTICAL APPROXIMATION FOR HEAT
C     EXCHANGE WITH CONFINING LAYERS.
C     IT USES THE METHOD OF VINSOME AND WESTERVELD, J. OF CANADIAN
C     PETROLEUM TECHNOLOGY, JULY-SEPTEMBER 1980, PP. 87-90.
C
C     THE CONFINING LAYERS ARE ASSUMED TO HAVE UNIFORM INITIAL
C     TEMPERATURE, TAKEN TO BE THE TEMPERATURE OF THE VERY LAST
C     ELEMENT APPEARING IN DATA BLOCK *ELEME*. THERMAL CONDUCTI-
C     VITY, HEAT CAPACITY, AND DENSITY OF THE CONFINING LAYERS
C     ARE TAKEN FROM THE DOMAIN TO WHICH THIS LAST ELEMENT BELONGS.
C
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      common/c13/ALPHA(3),IPMAT
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E6/T(MNEL)
      COMMON/VINWES/AI(MNEL)
C-----ARRAY *AI* HOLDS ONE PARAMETER PER ELEMENT WHICH CHARACTERIZES THE
C     TEMPERATURE DISTRIBUTION IN THE CONFINING LAYER. ON INITIALIZATION
C     ALL AI=0. AT THE END OF A RUN THE PARAMETERS *AI* ARE WRITTEN ONTO
C     A DISK FILE CALLED *TABLE*, WHICH NEEDS TO BE PROVIDED FOR A
C     RESTART.
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL)
C-----ARRAY *AHTRAN* HOLDS VALUES OF HEAT TRANSFER AREA FOR ALL
C     ELEMENTS WHICH ARE ADJACENT TO THE CONFINING LAYERS.
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/P3/DELX((MNK+1)*MNEL)
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/L3/CO(mnz+1)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +                CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/QQ/Q(200),QC(200),HTL,HTLC
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/BC/NELA
      CHARACTER*5 ELEM,MAT
C
      SAVE ICALL,T00,DIF
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.NE.1) GOTO 1
      WRITE(34,899)
  899 FORMAT(' ********QLOSS 1.01, 1995.1.10: PERFORM SEMI-ANALYTICAL'
     X' CALCULATION FOR HEAT EXCHANGE WITH CONFINING BEDS**********')
C
      IF(NEQ.EQ.NK1) GOTO 3
      WRITE (34,4)
    4 FORMAT(' ERRONEOUS CALL TO SUBROUTINE QLOSS SEMI-ANALYTICAL HEAT'
     x' EXCHANGE (MOP(15) > 0) IS ONLY AVAILABLE FOR NEQ = NK1, RESET'
     x' MOP(15) TO ZERO AND CONTINUE WITHOUT SEMIANALYTICAL HEAT LOSS')
      MOP(15)=0
      RETURN
    3 CONTINUE
C

      HTL=0.D0
      HTLC=0.D0
C
      DO 7 N=1,NELA
      AI(N)=0.D0
c      AI2(N)=0.D0
c      AI3(N)=0.D0
    7 CONTINUE

      REWIND 8
      READ(8,9,END=2) (AI(N),N=1,NELA)
    9 FORMAT(4E20.13)
    2 CONTINUE

c      REWIND 111
c      READ(111,9,END=13) (AI2(N),N=1,NELA)
c   13 CONTINUE

c      REWIND 112
c      READ(112,9,END=14) (AI3(N),N=1,NELA)
c   14 CONTINUE

C
C
C-----ASSIGN INITIAL TEMPERATURE AND DIFFUSIVITY OF CONFINING LAYERS.
c      DO 16 J=1,3
      T00=ALPHA(2)
C-----REPLACE WITH THE REFERENCE TEMPERATURE AND MATRIX MATERIALS
C-----OF THE CORRESPONDING DIRECTION
c      T00=T(NEL)
c      NELMAT=MATX(NEL)
      DIF=CWET(IPMAT)/(DM(IPMAT)*SH(IPMAT))
      WRITE (34,30) T00,CWET(IPMAT),DM(IPMAT),SH(IPMAT),DIF
   30 FORMAT(/' PERFORM SEMI-ANALYTICAL HEAT EXCHANGE CALCULATION'/
     X' THERMAL PARAMETERS ARE:'/' TEMPERATURE = ',E12.5,' HEAT'
     X' CONDUCTIVITY = ',E12.5,' DENSITY = ',E12.5,' SPECIFIC HEAT = ',
     xE12.5,' DIFFUSIVITY = ',E12.5//)
    1 CONTINUE
C
      IF(MOP(15).GE.2) WRITE (34,11) KCYC,ITER
   11 FORMAT(' SUBROUTINE QLOSS ----- (KCYC,ITER) = (',I4,I3,')')
C
      HTL=0.D0

c      OPEN(UNIT=298,FILE='DS.DAT',STATUS='UNKNOWN')
c      WRITE(298,157) SUMTIM,DELTEX
c  157 FORMAT('TIME=',E12.6,'TIME STEP=',E12.6)
C
      DO 5 N=1,NELA
C
c      IF(AHT(N).EQ.0.D0.OR.T(N).GE.AMTT(N)) GOTO 5
      IF(AHT(N).EQ.0.D0) GOTO 5

      NLOC=(N-1)*NEQ
      NLK1=NLOC+NK1
      NLOC2=(N-1)*NSEC*NEQ1
C
      TCUR=PAR(NLOC2+NSEC-1)
C-----ASSIGN CONDUCTION DISTANCE
      D=SQRT(DIF*(SUMTIM+DELTEX-STIME(N)))/2.D0
      DIFDT=DIF*DELTEX

      THETA=TCUR-T00
      THETAK=T(N)-T00

      PNUM=DIFDT*THETA/D-(THETA-THETAK)*D**3.D0/DIFDT
      IF(KCYC.GE.2) PNUM=PNUM+AI(N)
      PDEN=3.D0*D*D+DIFDT
      PP=PNUM/PDEN
      QNK1=ALPHA(1)*CWET(IPMAT)*(THETA/D-PP)*DELTEX
      QQ=((THETA-THETAK)/DIFDT-THETA/(D*D)+2.D0*PP/D)/2.D0
      IF(KON.EQ.2) AI(N)=THETA*D+PP*D*D+2.D0*QQ*D**3.0D0
      FLOH=-QNK1*AHT(N)/DELTEX  ! in W
      IF(KON.EQ.2) GOTO 6
C
      R(NLK1)=R(NLK1)+QNK1*AHT(N)/EVOL(N)
C
C-----IDENTIFY BEGINNING OF MATRIX ELEMENTS FOR ENERGY EQUATION
C     IN VOLUME ELEMENT N.
      N1KL=(N-1)*NEQ*NEQ+NEQ*NK
C
C-----COMPUTE DERIVATIVE OF LATERAL HEAT FLUX WITH RESPECT TO
C     TEMPERATURE.
      DPPDT=(DIFDT/D-D**3.D0/DIFDT)/(3.D0*D*D+DIFDT)
      DQK1DT=ALPHA(1)*CWET(IPMAT)*DELTEX*(1.D0/D-DPPDT)
      IF(X(NLOC+2)+dx(nloc+2).GT.1.5D0) GOTO 20
C
C-----COME HERE FOR TWO-PHASE CASE.
      DTDP=(PAR(NLOC2+2*NSEC-1)-TCUR)/DELX(NLOC+1)
      CO(N1KL+1)=CO(N1KL+1)-DTDP*DQK1DT*AHT(N)/EVOL(N)
      GOTO 21
C
   20 CONTINUE
C-----COME HERE FOR SINGLE PHASE CASE.
      CO(N1KL+2)=CO(N1KL+2)-DQK1DT*AHT(N)/EVOL(N)
   21 CONTINUE
C
    6 CONTINUE
C
      IF(MOP(15).GE.2)
     XWRITE (34,10) ELEM(N),FLOH,QNK1,AI(N),TCUR,DQK1DT
   10 FORMAT(10H ELEMENT *,A5,16H*   ---   FLOH= ,E12.5,9H   QNK1= ,E12.
     A6,10H   AI(N)= ,E12.5,9H   TCUR= ,E12.5,9H   DQDT= ,E12.5)

c      WRITE(298,155) ELEM(N),D,STIME(N)
c  155 FORMAT('ELEMENT #',A5,6X,'D=',E12.6,6X,'RTIME=',E12.6)
C
      IF(MOP(15).GE.3) WRITE (34,12) PP,QQ,DPPDT
   12 FORMAT(27X,'PP= ',E12.5,'     QQ= ',E12.5,22X,
     X'  DPPDT= ',E12.5)
C
      HTL=HTL-FLOH
C
    5 CONTINUE
c   16 CONTINUE
      CLOSE(298)
      RETURN
      END
C
C
