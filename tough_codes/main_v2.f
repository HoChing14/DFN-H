
      PROGRAM TOUGHREACT_V2
c
C          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C          @                                                   @
C          @  TOUGH2, MODULE T2CG2, VERSION 2.0, June 1999     @
C          @                                                   @
C          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C-----PROGRAM MAIN IS A HIGH-LEVEL EXECUTIVE ROUTINE. THAT CALLS SEVERAL
C     LOWER-LEVEL EXECUTIVE ROUTINES, WHICH EXECUTE ACTUAL CALCULATIONS.
C
C     ALL LARGE ARRAYS ARE DIMENSIONED IN PARAMETER STATEMENTS IN AN
C     INCLUDE FILE 'flowpar_v2.inc'. IF MODIFICATIONS IN ARRAY DIMENSIONS ARE
C     DESIRED, THESE NEED TO BE DONE ONLY IN FILE 'flowpar_v2.inc'.
C=======================================================================
C         related read and write files:
C	      #11#           <VERS>		    
C	      #3#            <GENER>		    
C	      #15#           <LINEQ>		    
C	      #8#            <TABLE>		    
C	      #4#            <MESH>		    
C	      #1#            <INCON>		    
C	      #7#            <SAVE>		    
C             #33#           <flow.inp>		
C             #34#           <flow.out>
C             #10#           <MINC>
C             #15#           <LINEQ>
C             #12#           <FOFT>
C             #14#           <COFT>
C             #13#           <GOFT>
C==================================================================
C$$$$$$$$$ ASSIGN PARAMETERS FOR FLEXIBLE DIMENSIONING OF LARGE ARRAYS $
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
C
      INCLUDE 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'perm_v2.inc'
C
C=======================================================================
C
C
C $$$$$$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)    ! darcy velocity
      COMMON/SOLUTE6/SLOLD(MNEL)           ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)           ! old gas saturation
      COMMON/SOLUTE8/SL1(MNEL)             ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)             ! new gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)         ! porosity at previous time step
      COMMON/XYZ11/XXX(MNEL)               ! for TECPLOT
      COMMON/XYZ22/YYY(MNEL)               ! for TECPLOT
      COMMON/XYZ33/ZZZ(MNEL)               ! for TECPLOT
c
c
      COMMON/Kchange/FK(MNEL) !fraction of original permeability due to dis./pre.
      COMMON/PcChange/ Fc(MNEL) ! fraction of capillary pressure
      COMMON/Rdphi/RPHI(MNEL)         ! Porosity change rate (unit time)
C
C--------------------Interface area reduction factor (TX, 12-Aug-1999)
      common/afactor/a_fm2(mncon)  ! advection area reduction (flow from F to M)
      common/afactord/a_fmd(mncon) ! diffusion area reduction (Both sides)
c.... Save modified active fracture area for reaction --
c     limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
C
c------------------ Permeability-porosity laws
      common/ipplaw/ikplaw(mnod)
      common/dpplaw/aparpp(mnod),bparpp(mnod)
C
c-------Added gaseous species properties (mol wt and mol diam)
      common/gasprop/dmwgas(mgas),diamol(mgas)
C
C----------------CO2 partial pressure   For use Steve's CO2 module
      COMMON/PCO2_ALL/PCO2A(MNEL)
C
      PARAMETER (mnz1=(mnel+2*mncon)*2)
      COMMON/LL1/IRNO(mnz1)
      COMMON/LL2/ICNO(mnz1)
      COMMON/LL3/COO(mnz1)    ! solving solute transport
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      COMMON/E1/ELEM(MNEL)  ! ELEMENT CODE NAMES (IDENTIFIERS)
      COMMON/E2/MATX(MNEL)  ! MATERIAL DOMAIN IDENTIFIER
      COMMON/E3/EVOL(MNEL)  ! VOLUME
      COMMON/E4/PHI(MNEL)   ! POROSITY
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)     ! PRESSURE AND TEMPERATURE
C----------------------------------------------------- For using EOS9
      common/TEM_EOS9/Tc_EOS9(MNEL)  ! initial temperature (oC)
c--------------------------------------------------------------------
c
      common/e7/pm(mnel)
      COMMON/VINWES/AI(MNEL)  ! HEAT EXCHANGE WITH CONFINING LAYERS
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL) ! HEAT TRANSFER AREAS
c.....7-20-93: define coordinate arrays
      common/xyz1/x1(mnel)
      common/xyz2/x2(mnel)
      common/xyz3/x3(mnel)
C
C$$$$$$$$$ COMMON BLOCKS FOR LATEST PRIMARY VARIABLES $$$$$$$$$$$$$$$$$$
C
C     NK IS THE NUMBER OF MASS BALANCE EQUATIONS PER GRID BLOCK.
C     NEQ IS THE NUMBER OF BALANCE EQUATIONS PER GRID BLOCK.
C
      COMMON/P1/X((MNK+1)*MNEL)     ! PRESENT PRIMARY VARIABLE VALUES
      COMMON/P2/DX((MNK+1)*MNEL)    ! INCREMENTS 1
      COMMON/P3/DELX((MNK+1)*MNEL)  ! INCREMENTS 2 FOR DERIVATIVES
      COMMON/P4/R(MNEQ*MNEL+1)      ! RESIDUALS
      COMMON/P5/DOLD(MNEQ*MNEL)     ! ACCUMULATION TERMS AT BEGINNING
      COMMON/P6/ROW(MNEQ*MNEL)
      COMMON/P7/COL(MNEQ*MNEL)  ! SCALING FACTORS FOR MATRIX ROWS/COLUMNS
C
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)  ! INDICES OF TWO ELEMENTS
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)  ! NODAL DISTANCES TO INTERFACE
      COMMON/C4/DEL2(MNCON)
      COMMON/C5/AREA(MNCON)  ! INTERFACE AREAS
      COMMON/C6/BETA(MNCON)  ! COSINE OF ANGLE BETWEEN NODAL LINE AND
C                              GRAVITY
      COMMON/C7/ISOX(MNCON)  ! ISOTROPY INDICES
      COMMON/C8/GLO(MNCON)   ! HEAT FLOW RATES
      COMMON/C9/ELEM1(MNCON)
      COMMON/C10/ELEM2(MNCON)  ! CODE NAMES OF ELEMENTS
      COMMON/C11/FVD(MNCON)
c.....6-8-95: new array SIG(MNCON) to hold coefficients for
c             radiative heat transfer.
      common/c12/sig(mncon)
      common/c13/ALPHA(3),IPMAT
c.....6-4-93: append two common blocks used in the T2VOC version
c             of TOUGH2.
      common/flovp1/flovg(mnph*mncon)
      common/flovp2/flovw(mnph*mncon)
c.....diffusive fluxes
      COMMON/FMOLDIF/FDIF(MNCON*MNPH*MNK)
      COMMON/VOLBC/VBC(MNEL)
c.....6-09-99: COMMON blocks for T2DM
C***  BEGIN ADDITION FOR T2DM
C
      COMMON/NUMGBXYZ/NGB(3)
      COMMON/NDUPPAR/NZMULTI,NZDISF,NZDUP
      COMMON/DARVELM/DVELM((2*MNEQ+1)*MNPH*MNCON)
      COMMON/FMECDIS/FDIS(MNCON*MNPH*MNK)
C
      PARAMETER (NDIMM = 2)
      PARAMETER (NNZERO = (4*NDIMM+1)*MNEQ*NREDM)
C
      COMMON/SD1/N2Z(NNZERO)
      COMMON/SD2/IDLI(NNZERO)
      COMMON/SD3/IDLII(NNZERO)
      COMMON/SD4/IDLIII(NNZERO)
      COMMON/SD5/IFLAGG(NNZERO)
      COMMON/SD6/NFLAGG(NNZERO)
      COMMON/SD7/N3ZMRK,ID
C
C***  END ADDITION FOR T2DM
C
C$$$$$$$$$ COMMON BLOCKS FOR LINEAR EQUATIONS $$$$$$$$$$$$$$$$$$$$$$$$$
C
c     arrays used by the conjugate gradient package. 
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)  ! ROW/COLUMN INDICES OF NON-ZERO MATRIX ELEMENTS
      COMMON/L3/CO(mnz+1)
C     WORKSPACE *WKAREA* HAS A LENGTH OF NEQ*NEL+10.
      COMMON/L4/WKAREA(MNEQ*MNEL+10) 
C     COLUMN INDICES OF NON-ZERO MATRIX ELEMENTS, HAS A LENGTH OF NZ
      COMMON/L7/JVECT(niwork)
c
c     array used by conjugate gradient solvers only.
c     COMMON/CGARA6/RWORK(NRWORK)
c     COMMON/SOLVER/MATSLV,NMAXIT,ICLOSR,CLOSUR,ISYM,IUNIT,NVECTR,seed
      common/soll/lenw,leniw
C
C     arrays used by luband only.
      COMMON/lub1/AB(nrwork)
      COMMON/lub3/NSUPDI,NSUBDI,mnzp1,mnetp1,mnelp1,nnnbig
      COMMON/lub4/matord,nsubdg,nsupdg,ntotd
C
C********* END modification.
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G1/F1(MGTAB)
      COMMON/G2/F2(MGTAB)
      COMMON/G3/F3(MGTAB)
      common/g3a/pw(mgtab)
C     ARRAYS F1, F2, AND F3 HOLD, RESPECTIVELY, TABLES OF TIMES, RATES,
C     AND ENTHALPIES OF TIME-DEPENDENT SINKS OR SOURCES. THEIR DIMENSIONS
C     MUST NOT BE SMALLER THAN THE TOTAL NUMBER OF SUCH DATA.
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
C     THE DIMENSION OF ALL ARRAYS IN COMMON G4 THROUGH G25 MUST BE
C     EQUAL TO OR LARGER THAN THE TOTAL NUMBER OF SINKS AND SOURCES.
C
      COMMON/G26/FF(MNPH*MNOGN)
C     THE DIMENSION OF ARRAY FF MUST BE EQUAL TO OR LARGER THAN THE
C     PRODUCT OF NUMBER OF PHASES AND TOTAL NUMBER OF SINKS AND SOURCES.
c
      common/g27/fnam(mnogn)
      common/g28/nftab(mnogn)
      common/g29/iftit(mnogn)
      common/g30/jftit(mnogn)
      common/g31/ijf(mnogn)
C
      character*1   ITABG
      character*5   fnam
c.....Added Commons for effective thermal conductivity
      common/efkth/timkth(mgtab),fackth(mgtab)
      common/kthtable/ktftb(mnogn)
C
C     *PAR* SECONDARY (THERMOPHYSICAL) PARAMETERS FOR ALL ELEMENTS
C     ARRAY PAR HAS A LENGTH OF (NEQ+1)*NSEC*NEL.
C     NSEC = NPH*(NB+NK)+2 IS THE NUMBER OF SECONDARY PARAMETERS.
C     NPH IS THE NUMBER OF PHASES.
C     NB IS THE NUMBER OF THERMOPHYSICAL PARAMETERS (USUALLY 6).
C     NK IS THE NUMBER OF COMPONENTS (SPECIES).
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/COMPO/FLO(MNPH*MNCON)  ! RATES OF FLOW ACROSS INTERFACES
      COMMON/PORVEL/VEL(MNPH*MNCON)  ! PORE VELOCITIES OF FLOW
      COMMON/TITLE/TITLE
      CHARACTER*80 TITLE
      CHARACTER*5 ELEM,ELEM1,ELEM2,ELEG,SOURCE
      LOGICAL EX
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/LDIM/LICN,LIRN
      COMMON/BC/NELA
      COMMON/V/IS
      common/ran2/iran
      COMMON/SOLVR1/matslv,nmaxit,nnvvcc,iiuunn,iissoo,nactdi
      COMMON/SOLVR2/ritmax,closur
      COMMON/SOLVR3/ordrng,oprocs,zprocs,coord
      CHARACTER*2 ordrng,oprocs,zprocs
      CHARACTER*5 coord
C
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX    ! TOUGH2 V2
c
      double precision elt
      real*4 rtzero,relt
      LOGICAL ECG

      INQUIRE(FILE='flow.inp',EXIST=ECG)
      IF(ECG) GOTO 1
      WRITE (*,*) ' FILE *flow.inp* DOES NOT EXIST --- EXIT PROGRAM'
      STOP
    1 OPEN (UNIT=33,FILE='flow.inp',STATUS='UNKNOWN')
      OPEN (UNIT=34,FILE='flow.out',STATUS='UNKNOWN')   

      write (*,*) '     -------- TOUGHREACT Version 2.0 -------'
      write (*,*)
      WRITE (*,*) '   --> reading multiphase flow input data'
!
      mnelp1 = mnel+1
      mnzp1  = mnz+1
      mnetp1 = nredm+1
      lenw   = nrwork
      leniw  = niwork
C
      WRITE (34,4)
    4 FORMAT(' TOUGH2 IS A PROGRAM FOR NON-ISOTHERMAL MULTIPHASE'
     X' MULTICOMPONENT FLOW IN PERMEABLE MEDIA. IT IS A MEMBER OF THE'
     X' MULKOM FAMILY OF CODES DEVELOPED AT LBNL.'//' PROGRAM'
     X' TOUGHREACT 2.0, 2009.3.30'/' Special version for conjugate'
     x' gradient package T2CG2, includes definition of coordinate'
     x' arrays and radiative heat transfer capability.')
C
c     revise length assignment of MA28 arrays
      LIRN=MNZ
      LICN=MNZ
      MPRIM=(MNK+1)*MNEL
      MSEC=(MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL
      WRITE (34,6) MNEQ,MNK,MNPH,MNB
    6 FORMAT(/' PARAMETERS FOR FLEXIBLE DIMENSIONING OF MAJOR ARRAYS '
     X'(MAIN PROGRAM): MNEQ =',I3,'   MNK =',I3,'   MNPH =',I3,
     X'   MNB =',I3)
      WRITE (34,7) MNEL,MNCON,MPRIM,MNOGN
    7 FORMAT(/' MAXIMUM NUMBER OF VOLUME ELEMENTS (GRID BLOCKS) AND',
     X' CONNECTIONS (INTERFACES):',I8,1X,I8/' MAXIMUM LENGTH OF'
     X' PRIMARY VARIABLE ARRAYS (MPRIM):',I8/' MAXIMUM NUMBER OF'
     X' GENERATION (SINKS/SOURCES) ITEMS (MNOGN):',I8)
      WRITE (34,9) MGTAB,MSEC,MNZ,LIRN,LICN
    9 FORMAT(' MAXIMUM NUMBER OF TABULAR (TIME-DEPENDENT) GENERATION'
     X' DATA',I8/' LENGTH OF SECONDARY PARAMETER ARRAY (MSEC):',I8/
     X' MAXIMUM NUMBER OF JACOBIAN MATRIX ELEMENTS (MNZ):',I8/
     X/' LARGE LINEAR EQUATION ARRAYS: LENGTH OF IRN =',I8,' LENGTH OF' 
     X' ICN AND CO = ',I8//' array dimensioning is made according to'
     X' the needs of the conjugate gradient solvers when using LUBAND,'
     X' only a smaller-size problem can be accommodated'/' restriction'
     X' with MA28 is: NEL+2*NCON<{MNEL+2*MNCON}/4'/)
c
      CALL IO
      call cpu_time(rtzero)
      tzero = dble (rtzero)
      READ(33,'(A80)') TITLE
      WRITE(34,3) TITLE
    3 FORMAT(' TITLE:',A80)
      CALL INPUT
      IF(IS.NE.0) GOTO 100
      CALL FLOP
      CALL RFILE
      WRITE (34,10) NEL,NELA,NCON,NOGN
   10 FORMAT(/' MESH HAS',I7,' ELEMENTS (',I7,' ACTIVE) AND',I7,
     A' CONNECTIONS (INTERFACES) BETWEEN THEM, GENER HAS',I6,
     A' SINKS/SOURCES')
C
c.....set solver type and make informative printout
      write(34,*) ' ****call subroutine sinsub in t2cg22_v2.f****'
      call sinsub
      write(34,*) ' ****end of call subroutine sinsub****'
      IF(MOP(7).NE.0) CALL INDATA
      IF (IS.EQ.0) CALL CYCIT ! SOLVE MASS- AND ENERGY-TRANSPORT EQUATIONS
  100 CONTINUE
      call cpu_time(relt)
      elt = dble (relt)
      ELT=ELT-TZERO
      WRITE (34,12) ELT
   12 FORMAT(/' END OF TOUGH2 SIMULATION RUN --- ELAPSED TIME = ',
     AF10.3,'S')
      STOP
      END
!
      SUBROUTINE CYCIT
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
      include 'common_v2.inc'
      include 'perm_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)   ! darcy velocity
      COMMON/Kchange/FK(MNOD)             ! fraction of original permeability due to dis./pre.
      COMMON/SOLUTE6/SLOLD(MNEL)
      COMMON/SOLUTE7/SGOLD(MNEL)          ! old liquid/gas saturation
      COMMON/SOLUTE8/SL1(MNEL)
      COMMON/SOLUTE9/SG1(MNEL)            ! new liquid/gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)        ! old (initial) porosity
      COMMON/PORVEL/VEL(MNPH*MNCON)
      COMMON/PARNP/NPL,NPG                      ! specify in EOS module
      COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)     ! Water density (kg/dm**3)
      COMMON/PRINTC/NOW                         ! print control
      COMMON/AMMISC/IABC,ISOLVC
C
      COMMON/SOLID/NM,DROK(maxmat),POR(maxmat),PER(3,maxmat),
     +   CWET(maxmat),SH(maxmat)     ! for perm initialization
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
c
c     Add tortuosity exponent (ptort) and critical porosity (phicrit)
      common/torpar/ptort(maxmat),phicrit(maxmat)
c     Add multiphase tortuosity at each grid block  
      common/tortmp/tortliq(mnel),tortgas(mnel)
c
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL)
      COMMON/E5/P(MNEL)
c--------------------------------------------------------------------
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      common/c13/ALPHA(3),IPMAT
      COMMON/G7/G(MNOGN)
C---------------------------------------------------------- For chemical QSS
      COMMON/STEADY/IFLOWSS,JSTEADY
      COMMON/TOL_STEADY/TOLDC,TOLDR    ! Concentration and dis/pre changes
      COMMON/BMNO/NBLOCK,NMINERAL !Block and mineral numbers where mineral exhausted
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
C---------------------------------------------------------------------------
      COMMON/GASLAW/R,AMS,AMA,CVAIR
      DATA AMS /18.016d0/, AMA /28.96D0/
C
c----------------------------------------Indicators from EOS module
C
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
C----------------------------------- For coupling with reactive transport
C
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
      common/TDS_REACT1/ iTDS_REACT
!
      COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
C----------------------------------------------------For using EOS2 Flow module
C----------------------------------And calculating co2 gas fugacity coefficient
!
      common/co2_gene1/ nco2g
      common/co2_gene2/ ico2gt0        !=1: initial Pco2>0
      COMMON/REACTh2o/Rh2o(NMNOD)      ! H2O REACTION SOURCES
      COMMON/REACTco2/Rco2(NMNOD)      ! CO2 REACTION SOURCES
      COMMON/SOLIDco2/SMco2(NMNOD)     ! CO2 TRAPPED in solid phase
!
c-------------------For H2 generation by mineral phase using EOS5 module
!
      common/h2_gene1/ nh2g
      common/h2_gene2/ ih2gt0     !=1: initial Ph2>0
!
C----------------CO2 partial pressure    For use Steve's CO2 module
      COMMON/PCO2_ALL/PCO2A(MNEL)      ! calculated from ECO2
!
c---------------------------------------------------------------------
      COMMON/IKC/IK        ! indicator of change in permeability K
      COMMON/PcChange/ Fc(MNEL)
      COMMON/Rdphi/RPHI(MNEL)           ! Porosity change rate (unit time)
c---------------------------------------------------------------------
c
      COMMON/PLOT_FM/IFM      ! separate printout for fracture and matrix
      COMMON/XYZ11/XXX(mnel)     !!!! mesh coordinates for TECPLOT
      COMMON/XYZ22/YYY(mnel)
      COMMON/XYZ33/ZZZ(mnel)
C
c.....Define coordinate arrays
      common/xyz1/x1(mnel)
      common/xyz2/x2(mnel)
      common/xyz3/x3(mnel)
C
C-----THIS ROUTINE IS THE EXECUTIVE ROUTINE FOR MARCHING IN TIME.
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/E1/ELEM(MNEL)
      COMMON/E6/T(MNEL)
      COMMON/VINWES/AI(MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/STEP/ELST
c     Added commons for writing out element specific data (TOUGH2 v1.6)
      COMMON/STEPrk1/eplist(200)
      COMMON/STEPrk2/nstrick(200),nelist
      character*5  eplist
c
      COMMON/STE1/NST
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
      COMMON/KC/KC
      COMMON/DFM/TIMAX,REDLT
      COMMON/DOP/ENTH,KDATA,QUAL
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/DMN/INUM,IPRINT,MCYC,MCYPR,MSEC,TZERO
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      COMMON/FAIL/IHLV
      COMMON/BC/NELA
      common/pseduo/ps_ads,con
      common/fgt1/ioft,iofu,igoft,igofu,noft(100),ngoft(100)
      common/fgt3/icoft,icofu,ncoft(100)
      CHARACTER*5 ELEM,ELST
c
      COMMON/MODELT/IDELT   ! IDELT=1 if DELT modified by TSTEPT1
c
      real*4 rtcur
      double precision tcur
C---------------------------------------------------------------------------
c... Commons for effective thermal conductivity
      common/efkth/timkth(mgtab),fackth(mgtab)
      common/kthtable/ktftb(mnogn)
c
c...  Added all kinetic common blocks below
      common/minkin3/ cr0(mpri)
      common/isarea/imflg2(mmin),imflag(mnel,mmin)
c... .Dissolution kinetics
      common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
      common/iprkin/ideprec(mmin)
c...  Added common block for rate law designations
      common/irtlaw/nplaw(mmin)
c...  next_tstep not defined previously
      integer*8 next_tstep
c...  Save modified active fracture area for reaction --
c     limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
c
      COMMON/TRANGAS8/dcfgas(mgas,mnel),DIFUNG     ! gaseous species diffusion coefficient
c
c            porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
C-----------------------Common blocks for dryout grid blocks
      COMMON/DRYOUT/IDRY(MNOD),ADRY(MNOD,MPRI)
      COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)    ! residual in precipitates
      common/dry_salt/nsalt,isalt(mmin)
      common/water_activity/aw(mnod)    !water activity
c     Time step-limiting info
      common/dtlim/max_chem_it,delt_conne,id_chem
      character*16 delt_conne
      character*5 id_chem
C
      common/kcycsav/kcyc_last
!
!.....Extract rock density for geochemical calculations such as exchange and sorption
!
      common/rock_density1/denss(mnel)
      common/i_restart1/ i_restart
c
!.....Save Henry's constant for air solubility in water
      SAVE HC
      DATA HC/1.d-10/
      WRITE(34,891)
  891 FORMAT(/1X,60('-'),'call CYCIT for simulation',60('-')/' CYCIT'
     X' 1.01, 1998.5.28: EXECUTIVE ROUTINE OF TOUGH2'/)
c
c-----Write header of elemement data in GASOBS.DAT
      if(icall.eq.1.and.nelist.gt.0)  write(67, 232)
c
      KC=0
      SUMTIM=TIMIN
      TIMOUT=SUMTIM
      KON=1
      DELT=DELTEN
      MA=0
      IAB=0
      IPIV=0
      IGOOD=0
      NSTL=(NST-1)*NK1
      NSTL2=(NST-1)*NEQ1*NSEC
!
      ITER=0
      IHLV=0
      KIT=0
      IQIT=0
!F
!.....Needed to be initialized
      next_tstep = 0
      kcyc_last = 0
!
!.....Modify permeability and porosity according to "SEED" in flow.inp
      if (mopr(15) .ge. 1)   then
         CALL Modify_PoroPerm
      end if
!
!.....Saves grid block connection details for time step-limiting option
! 
      delt_conne = 'Flo_Delta_t     '
      id_chem = '     '
c
C
C*************************************************************************************
c
c     Add these initialization for each grid node.
c     First we assign permeability at each grid node to
c     array perm.  Need this here because perm is now used
c     directly in eos (routine MULTI).
c     per(iso=1 to 3,nmat) is permeability for each material (input)
c     perm(iso=1 to 3,n) is permeability assigned  at each grid node
c      IF(MOP(15).GE.1.AND.SUMTIM.GT.0.D0) OPEN(UNIT=112,FILE='RTIME',
c     &STATUS='UNKNOWN')
      do n=1,nel
c         IF(MOP(15).GE.1.AND.SUMTIM.GT.0.D0) READ(112,1461) STIME(N),
c     &MLAGNR(N)
c 1461    FORMAT(2X,E12.6,2X,I1)
         do i=1,3
c----------save permeability from incon
           perm0(i,n) = perm(i,n)
           permm(i,n) = perm(i,n)
         end do
c----------save porosity from incon
           phi0(n) = phi(n)
           phim(n) = phi(n)
c----------move initialization of phiold up here
           phiold(n) = phi(n)
c
c     pcfact is a factor to multiply the calculated
c     capillary pressure (in routine PCAP) to account for the
c     effect of changing porosity and permeability on capillary
c     pressure.  It must be initialized to 1 and remains 1 if
c     porosity is not coupled to permeability (i.e. if kcpl=0)
         pcfact(n)=1.d0
      end do
c
      nnod=nel
c.....Added Leverett scaling even if perm doesn't change (in rctprop_v2.f)
      if (mopr(6).eq.1)   call levscale
C
C*************************************************************************************
C
      CALL EOS
      CALL BALLA  ! balance calculation

C.....Give initial valus of FK (fraction of original permeability  due to dis./pre)
C
      DO N=1,NEL
         FK(N)   = 1.0D0        ! permeability factor     
         Fc(N)   = 1.0D0        ! capillary pressure factor
         Rphi(N) = 0.0D0        ! porosity change rate (unit time)
         Rh2o(N) = 0.0D0        ! H2O reaction sources
      END DO
c
c+++++++++++++++++++++++++++++++++Begin addition for solute transport
c
      IF (MOPR(1).ne.2) THEN          ! include reactive transport
         DO N=1,NEL
            XXX(N)=X1(N)
            YYY(N)=X2(N)
            ZZZ(N)=X3(N)
         END DO
c
         IABC=0                     
c
c.....keep initial values
c
         DO N=1,NEL
            NLOC2=(N-1)*NSEC*NEQ1     ! start of sec. variables for N
            NP=NPL
            NL2NP=NLOC2+(NP-1)*NBK
            SLOLD(N)=PAR(NL2NP+1)     ! initial liquid saturation
            SGOLD(N)=1.0D0-SLOLD(N)   ! initial gas saturation
c...........to avoid zero values in plot files at t=0
            sl1(n)=slold(n)
            sg1(n)=sgold(n)
c...........Add temperatures
            tc(n) = t(n)
            tkelv(n) = tc(n) + 273.15d0
         END DO
!
! -------------
!........Initial solid rock density
! -------------
!
         do n=1,nel
!
            nmat = matx(n)
            denss(n) = drok(nmat)/1.0d3   ! In kg/dm**3
!
!...........Sometime use a inifinitive rock density to act as constant 
!...........temperature condition, in this case do the following:
!
            if (denss(n) .le. 0.0d0 .or. denss(n) .ge. 20.0d0)  then
               denss(n) = 2.65d0   ! In kg/dm**3
            end if    
         end do
!
         kcpl = 0     ! Initialize variable 
C
C+++++++++++++++++++++++++END addition for solute transport
C+++++++++++++++++++++++++++++++++Begin addition for coupling
C
         NNOD=NEL
         WRITE(34,*) '*******READ TRANSPORT AND CHEMISTRY*******'
         write(34,*) '*****read solute.inp in readsolu_v2.f*****'
         CALL READSOLU
         write(34,*) '****read chemistry input data in inichm_v2.f****'
         CALL INIT ! read initial chemistry data (in inichm_v2.f)
         write(34,*) '***convert indices to names in readsolu_v2.f***'
         call get_prt_dat
         
c--------Set water activity zero if dry
         DO N=1,NEL     !gxzh   V3.3Q  Piz
          if (sl1(n).lt.sl1min*1.0d-03) aw(n)=0.0d0            
         ENDDO
C--------Read restart data   for reactive transport
         i_restart = 0
         IF(KCYC.GE.1)   THEN
            OPEN (UNIT=21,FILE='inchem',STATUS='UNKNOWN')
            CALL READ_RESTART ! read inchem file (in geochem_v2.f)
            i_restart = 1
cns3/2010  need this call on restart for surface complexes (in inichm_v2.f)
              if(nads.ne.0) call surfequil
            CLOSE (UNIT=21)
c           See other use for this flag further down.
            next_tstep=2  ! skip some printing later.
         END IF  
C---------------------------------------------------------------------
C
         ITERFL=ITER
c
         if(kcyc.eq.0) CALL WRITE_ITER(kcyc,sumtim)
c
      END IF
c     skip flow and/or chemistry in main and couple
      if(mopr(1).eq.1) then
        write(32,"(
     &  /5x,'!! Flow: Skip transport and chemistry,',
     &  ' MOPR(1) =', i3)") MOPR(1)
      else if(mopr(1).eq.3) then
        write(32,"(
     &  /5x,'!! Flow & Transport: Skip chemistry,',
     &  ' MOPR(1) =', i3)") MOPR(1)
      else if(mopr(1).eq.4) then
        write(32,"(
     &  /5x,'!! Transport & Chemistry: Skip flow,',
     &  ' MOPR(1) =', i3)") MOPR(1)
      else if(mopr(1).eq.5) then
        write(32,"(
     &  /5x,'!! Transport: Skip flow and chemistry,',
     &  ' MOPR(1) =', i3)") MOPR(1)
      endif

C===============================================================================================
!
!.....Write on computer screen
!
      write(*,*)
      write(*,*) '   --> performing simulation'
      write(*,*)
      WRITE(34,*) '***********MAIN ITERATIVE SIMULATION***************'
!.....Moved below three lines here (and the the beginning of CYCLT) from COUPLE
!.....in order to print out maximum between printouts 
!
      MAXITCH  = 0          ! Maximum iterations of solving whole chemistry
      AVERITCH = 0.d0       ! Average 
      COUNTCH  = 0.d0       ! Iteration counter
!
!+++++++++++++++++++++++++++++++++++End   addition for coupling
!
   30 KCYC=KCYC+1
!
! ----------
!.....Evaluate accumulation terms when changes in porosity due to reaction
!.....Subroutine Evaluate_DOLD called only including reactive transport
! ----------
     
      if (mopr(1).eq.0.or.mopr(1).gt.2)   then        
         CALL Evaluate_DOLD
      end if

! ----------
!.....Computes total dissolved solid from -REACT (for ECO2N)
! ----------
!
!     iTDS_REACT = 0        ! Now read from solute.inp
      if (mopr(1) .eq. 2.or.mopr(1).eq.1) iTDS_REACT = 0
!
!.....iTDS_REACT = 1 :  Pass TDS_REACT to TOUGH primary variable X
!.....iTDS_REACT = 2 :  Pass TDS_REACT to density subroutine in ECO2N
!
      if (iTDS_REACT .gt. 0 .and. ieos .eq. 14 )        then
         CALL TDS_From_REACT ! (in rctprop_v2.f)                   
      end if
C===============================================================================================  
C-----COME HERE FOR NEW TIME STEP.
C
      KC=KC+1
!
      IF(MOP(2).NE.0) WRITE (34,29) KC,KCYC,DELT,SUMTIM
   29 FORMAT(/23H CYCIT ..........  KC =,I6,8H  KCYC =,I7,9H  DELT = ,E1
     A2.6,11H  SUMTIM = ,E12.6/)
      IF(IQIT.NE.0) GOTO 10
c
c... changed according to ysw
      IF(abs(TIMAX).gt.0.d0.AND.abs(SUMTIM-TIMAX).le.0.d0) GOTO 10    
      IF(KC.GT.MCYC .AND. MCYC .LT. 9999) GOTO 10
      call cpu_time(rtcur)
      tcur = dble(rtcur)
      IF(MSEC.NE.0.AND.TCUR-TZERO.GT.MSEC) GOTO 10
C
      IHALVE=0
C
   52 ITER=0
      NOW=0
C
      iterfl=0
c
C-----PRINTOUT OCCURS FOR NOW=1
C
      NOWTIM=0
c
C***********************************************Modification for reactive transport
      IF(ITI.NE.0) CALL TSTEPT
C**********************************************************************************
c     
      DELTEX=DELT
      KON=1
!
c.....Add next_tstep flag for iteration between flow, transport,
c     and chemistry.  This flag is reset in subroutine couple (in treact_v2.f)
c     next_tstep=1 we are at the next time step
c     next_tstep=0 we are still at the current time step
c
      if(next_tstep.ne.2) next_tstep=1
!
    3 CONTINUE
C
C-----COME HERE FOR NEXT ITERATION.
C
      IF(ITER.NE.0) GOTO 19
      DO 20 N=1,NEL
          NLOC=(N-1)*NK1
          DO20 K=1,NK1
              DX(NLOC+K)=0.0D0
   20 CONTINUE
C
   19 CONTINUE
C
      ITER=ITER+1
!
!.....Total number of cumulative iterations for solving the flow
!
      if (ITERC .lt. 9999999)   then
         ITERC =ITERC + 1
      end if
!.....
!
!.....iterfl (number of flow interation for printout)
!     was originally set equal to iter.  iter is reset to zero
!     after each flow convergence.  However, if we iterate between
!     ch,tr and flow, we want iterfl to reflect the total number
!     of flow iterations, including sequential iterations
!
      iterfl=iterfl+1
!
      FORD=FOR*DELTEX
C
C-----COMPUTE ACCUMULATION-, SOURCE-, AND FLOW-TERMS (in multi_v2.f).
C
      IF (IEOS .EQ. 12)  THEN          ! For EWASG module
c         CALL MULTI_EWASG
      ELSE IF (IEOS .EQ. 9)   THEN     ! For EOS9 module
         CALL MULTI_EOS9
      ELSE
         CALL MULTI                 ! For other modules
      END IF
C**********************************************************************
      
      IF(IGOOD.EQ.3) GOTO 50
C
C-----CHECK FOR CONVERGENCE.
C
c.....Added this option to force at least one iteration
c.....modified from Yu Shu's version of t2fecm.f
c.....c'out below and replace with following block.
c.....If sequential iteration, we do not want to force at least
c.....one iteration.
c
      if(rerm.le.re1) then           !converge within given relative error
         if(mopr(4).eq.1) goto 51    !no flow calculation before chemistry
         if(next_tstep.eq.0) goto 51 !case for seq. iterations
         if(iter.gt.1) goto 51
      end if
c........
      IF(MOP(1).NE.0) WRITE (34,22) KCYC,ITER,DELTEX,RERM,ELEM(NER),KER
   22 FORMAT(' ...ITERATING...  AT [',I7,',',I3,'] --- DELTEX = ',E12.6,
     A' MAX. RES. = ',E12.6,'  AT ELEMENT ',A5,'  EQUATION ',I3)
      IF(MOP(1).LT.2) GOTO 23
      NIT=NER
      IF(NST.NE.0) NIT=NST
      NIT1=(NIT-1)*NK1
      NIT2=(NIT-1)*NEQ1*NSEC
      PNIT=X(NIT1+1)+DX(NIT1+1)
      WRITE (34,24) ELEM(NIT),DX(NIT1+1),DX(NIT1+2),
     XPAR(NIT2+NSEC-1),PNIT,PAR(NIT2+1)
   24 FORMAT(31X,' AT *',A5,'* ...   DX1= ',E12.6,
     X6H DX2= ,E12.6,5H T = ,F7.3,5H P = ,F9.0,5H S = ,E12.6)
   23 CONTINUE
C
      IF(ITER.GT.NOITE) GOTO 50  ! iter>8 diverge criterion
C
      CALL LINEQ ! LINEAR EQUATION SOLVER T2CG2 (in t2cg22_v2.f)
      IF(IGOOD.NE.0) GOTO 50
      CALL EOS ! recompute present state, update secondary parameters (in eos1_v2.f)
      IF(IGOOD.EQ.0) GOTO 54  ! parameters are reasonable
!
   50 DELT = DELTEX/REDLT  ! iteration diverged, REDUCE TIME STEP
      IGOOD=0
C     CHECK ON NUMBER OF PREVIOUS TIME STEPS THAT CONVERGED ON ITER=1 (no progress calculation).
      IF(KIT.LT.2.OR.MOP(16).EQ.0) GOTO 25
      WRITE (34,26) KCYC,DELTEX
   26 FORMAT(/' CONVERGENCE FAILURE ON TIME STEP #',I4,'WITH DT = '
     X,E12.6,'SECONDS, FOLLOWING TWO STEPS THAT CONVERGED ON ITER = 1'/
     X' STOP EXECUTION AFTER NEXT TIME STEP')
      IQIT=1
C
   25 CONTINUE
      WRITE (34,53) DELT,ITER,KCYC
      write (*, 53) DELT,ITER,KCYC
   53 FORMAT(' REDUCE TIME STEP TO' ,E12.6,' SECONDS AT (',I4,I3,')')
C
      IHALVE=IHALVE+1
      IF (IHALVE.LE.25.OR.DELTEX.GT.0.1D0) GOTO 11
      WRITE (34,12) DELTEX
   12 FORMAT(//' FAILURE IN CONVERGENCE'/
     X' LAST TIME STEP = ',E12.6,' SECONDS '//' STOP EXECUTION ')
c.....Write output files before stopping (eos1_v2.f and t2f_v2.f).
      if(nel.le.1000) call out
      CALL WRIFI
      IF (MOPR(1).ne.2.and.mopr(1).ne.1) THEN
         OPEN (UNIT=21,FILE='savechem',STATUS='UNKNOWN')
         CALL WRITE_RESTART ! write chemistry output (in geochem_v2.f)
         CLOSE (UNIT=21)
         CALL WRITE_ITER (KCYC,sumtim)  ! write chemistry iteration (in geochem_v2.f)
      END IF
C-------------------------------------------------------------------
c
      IF(MOP(15).NE.0) then
      REWIND 8
      WRITE(8,14) (AI(N),N=1,NELA)
c      REWIND 111
c      WRITE(111,14) (AI2(N),N=1,NELA)
c      REWIND 112
c      WRITE(112,14) (AI3(N),N=1,NELA)
      endif
      STOP
C
   11 ITER=0  ! 11: after reducing time step for divergence
      IHLV=1
      CALL EOS
      IHLV=0
      GOTO 52  ! recompute in new time step size in current time step
C
C-----COME HERE FOR CONVERGENCE.
   51 KON=2
C
      IF(MOP(15).GE.1) CALL QLOSS
c
c      Note: CONVER updates the flow primary variables and time step (in t2f_v2.f)
C********************************************************** for reactive transport

      IF (MOPR(1).EQ.2.or.mopr(1).eq.1) then
         CALL CONVER        ! only flow
      ELSE
         CALL CONVER2       ! both flow and transport
      END IF
c
c*********************************************************************************
c
C     COUNT NUMBER OF CONSECUTIVE TIME STEPS THAT CONVERGE ON ITER = 1.
      IF(ITER.GT.1) KIT=0
      IF(ITER.EQ.1) KIT=KIT+1
c
C********************************************************** for reactive transport
c
      NKIT=10000000
      IF (MOPR(1).EQ.2.or.mopr(1).eq.1)  THEN
         NKIT=100              ! Only flow or original TOUGH2
      END IF
!
      IF(KIT.LT.NKIT) GOTO 55            
!
!*********************************************************************************
! very close to steady-state finishing criteria
      WRITE (34,56) NKIT
   56 FORMAT(' FOR',I9,' CONSECUTIVE TIME STEPS HAVE CONVERGENCE ON'
     X' ITER = 1'/' WRITE OUT CURRENT DATA, THEN STOP EXECUTION')
      IQIT=1
   55 CONTINUE
!
!.....Add flag to limit printout and speed up code
      write(*,238) kcyc,sumtim,deltex,iter
  238 format('step=',I6,1x,'time=',E14.6,1x,'dt=',E14.6,1x,'iter=',I1)
      if(mopr(5).eq.0)then
        IF(KON.EQ.2.AND.NST.NE.0) WRITE (34,49) ELEM(NST),KCYC,ITER,
     A    SUMTIM,DELTEX,DX(NSTL+1),DX(NSTL+2),PAR(NSTL2+NSEC-1),
     B    X(NSTL+1),PAR(NSTL2+1)
        IF(KON.EQ.2.AND.NST.EQ.0) WRITE (34,49) ELEM(NER),KCYC,ITER,
     A    SUMTIM,DELTEX,DX((NER-1)*NK1+1),DX((NER-1)*NK1+2),
     X    PAR((NER-1)*NEQ1*NSEC+NSEC-1),
     B    X((NER-1)*NK1+1),PAR((NER-1)*NEQ1*NSEC+1)
        
   49   FORMAT(1X,A5,1H(,I6,1H,,I3,1H),6H ST = ,E12.6,6H DT = ,E12.6,
     A   ' DX1= ',E12.6,6H DX2= ,E12.6,5H T = ,F7.3,5H P = ,F9.0,
     b   ' S = ',E12.6)
      endif
c      DO 107 NR=1,NELA
c      NRL2=(NR-1)*NEQ1*NSEC
C      IF(KC.EQ.1) WRITE(*,*) ELEM(NR),STIME(NR),MLAGNR(NR)
c      IF(PAR(NRL2+NSEC-1).LT.AMTT(NR).AND.MLAGNR(NR).EQ.0) THEN
c      STIME(NR)=SUMTIM
c      MLAGNR(NR)=1
c      END IF
c  107 CONTINUE
c
c.....Added from TOUGH2 v1.6
c     ysw & rick 7/13/95
c
      do i = 1, nelist
cc
cc    Write additional output into GASOBS.DAT
cc
c.......Reduce calcs by adding temp var nstrm1
        nstrm1 = NSTrick(i)-1
        nloc=nstrm1*NK1
        nloc2=nstrm1*NEQ1*NSEC
        nloc2l=nloc2+nbk
        pres=x(nloc+1)
        tx=par(nloc2+nsec-1)
        sgwr=par(nloc2+1)
        slwr=1.0d0-sgwr
c........Change RH calcs such that they work with eos3 as well as eos4
c Relative humidity (RH) calculations (for output only)
c
         call sat(tx,psat)  !pure water saturation P (in t2f_v2.f)
         if(slwr.eq.1.d0) then   !fully liquid saturated
              pv=psat
         else if(slwr.eq.0.d0) then  !fully gas saturated (eos3 or 4)
              pv=par(nloc2l+5)
         else  !two phase
              if(ieos.eq.4) then
              pv = pres-x(nloc+3)
              else
                pv=psat
              end if
         end if
         RH=pv/psat
c
c        Added to write out capillary pressure
c
         do m=1,neq1
           nlm2=nloc2+(m-1)*nsec
            do np=1,nph
              nlm2p=nlm2+(np-1)*nbk
              pcwr = par(nlm2p+6)
            end do
         end do
cc    Write time in days,pressure,gas saturation, temperature,
cc    rel. humidity,liquid saturation
        write(67, 233) eplist(i),(sumtim/3.15569d7),pres,
     +      sgwr,tx,RH,pcwr, PV
      end do

 232  format(1x,'ELEM',3x,' TIME(YRS)','    GAS PRES.',
     1      '     GAS SAT.', '    TEMPERATURE','   REL. HUM.',
     2      '    CAP. PRES.','    H2O PRES.')
 233  format(a5,1x,e14.6,8(1x,e12.5))
!
      if (kon.eq.2.and.(iofu.gt.0.or.igofu.ne.0.or.icofu.ne.0))  then
         if (ieos .eq. 13 .or. ieos .eq. 14)   then
            call fgtab_ECO2
                                               else
            call fgtab
         end if
      end if
!
!.....DETERMINE WHETHER PRINTOUT IS REQUIRED (in eos1_v2.f).
      IF(IQIT.NE.0) NOW=1  ! no progress or steady-state criteria end of calculation
      IF(NOWTIM.EQ.1.AND.KON.EQ.2) NOW=1  ! sumtime equals to points defined in "TIMES"
      IF(KON.EQ.2.AND.MOD(KCYC,MCYPR).EQ.0.and.mcypr.ne.9999) NOW=1  ! coverge at printout time step
   54 IF(KDATA.GE.10.AND.MOD(KCYC,MCYPR).EQ.0.and.mcypr.ne.9999) NOW=1
      IF(SUMTIM.EQ.TIMAX .OR. KCYC.EQ.MCYC.and.mcypr.ne.9999) NOW=1  ! printout for end of calculation
!----------------------------------------------------
      IF(NOW.EQ.0) GOTO 33
      CALL OUT
      IF(KON.EQ.2) CALL BALLA
   33 CONTINUE
c
      IF(KON.EQ.1)  GOTO 3 ! Begin addition for reactive transport when kon=2 (flow converged)
C
c---------------------------------------------------------------------
      
      IF (MOPR(1).ne.2.and.mopr(1).ne.1) THEN     ! coupled with reactive transport
!
       IFLOWSS=0               ! steady state flow, no transport and chemistry
!
       if (MOPR(1).eq.4.or.mopr(1).eq.5) then   
          IFLOWSS=1   !skips flow calcs
          call out    !to print flow results once
          call wrifi  !to print save file once
          write(34,"(/' MOPR(1)=',i2,' Flow calcs are skipped'/)")
     &    MOPR(1)
       end if
499    CONTINUE
       JSTEADY=0          ! flag for chemical steady-state
C
          IF (IFLOWSS.EQ.1) THEN
            KCYC=KCYC+1
C
C-----COME HERE FOR NEW TIME STEP.
C
            KC=KC+1
C
            IF(TIMAX.NE.0.D0.AND.SUMTIM.EQ.TIMAX)  GOTO 699
            IF(TIMAX.NE.0.D0.AND.TIMETOT.GE.TIMAX) GOTO 699
            IF(KC.GT.MCYC.AND.MCYC.LT.9999) GOTO 699
            call cpu_time(rtcur)
            tcur = dble(rtcur)
            IF(MSEC.NE.0.AND.TCUR-TZERO.GT.MSEC) GOTO 699
C
            ITER=0          ! because no flow is solved
            ITERFL=0
            NOW=0
C-----PRINTOUT OCCURS FOR NOW=1
C
            NOWTIM=0
            IDELT=0
            IF(ITI.NE.0) CALL TSTEPT1
            DELTEX=DELT
c
            deltex = max(deltex,0.01d0)
            CALL CONVER1                                ! define next time step
 509         CONTINUE
            IF(NOWTIM.EQ.1) NOW=1  ! sumtim reach to defined in "TIMES" module
            IF(MOD(KCYC,MCYPR).EQ.0) NOW=1  ! time step printout point
            IF(SUMTIM.EQ.TIMAX .OR. KCYC.EQ.MCYC) NOW=1  ! end of simulation
         END IF
C
C********************************************************************************
c      write(*,*) ctot(1032,3),'a1'
      call couple(next_tstep)
c      write(*,*) ctot(1032,3),'a2'
c      cgold=ctot(1032,3)
c      write(*,*) ctot(1032,3),ctot(1032,3)/cgold,'a3'
c      do i=1,nel
c         con=1-aht(i)*ps_ads*deltex
c         ctot(i,3)=ctot(i,3)*con ! pseudo adsorption
c         ctot(i,4)=ctot(i,4)*con ! to keep charge balance
c      end do
           if(next_tstep.eq.0) then
               iter=0
              goto 3
           end if
!
           CALL CONVER3
           if (dabs(rcour).gt.0.0d0) CALL MAX_DELT
!
C********************************************************************************
c
         IF (JSTEADY .EQ. 1) THEN
          TIMEDAY=TIMETOT/86400.0d0
          WRITE (39,513) TIMEDAY,NBLOCK,NMINERAL
          WRITE (32,513) TIMEDAY,NBLOCK,NMINERAL
513         FORMAT(/' A chemical steady-state is reached until',
     +      ' time',E11.4,'days.',/4x,' Only update mineral phases',
     +      '(Block=',I5,',    Mineral=',I5,')')
          WRITE (39,*) '                     -----         '
c
         END IF
c
         IF (IFLOWSS.EQ.1) GOTO 499
         NKIT=100
         I_QSS=1      ! indicator of chemical quasi-steady state
         IF (TOLDC.EQ.0.0D0.OR.TOLDR.EQ.0.0D0) I_QSS=0
         IF (IEOS.EQ.2.OR.IEOS.EQ.12.OR.IEOS.EQ.13) I_QSS=0
         IF (IEOS.EQ.14.OR.IEOS.EQ.5)                  I_QSS=0
         IF ((MOPR(1).EQ.0.or.mopr(1).gt.2).AND.I_QSS.EQ.1 
     &      .AND. KIT.EQ.NKIT) THEN   ! steady-state flow
         IFLOWSS=1
            WRITE (34,57)
   57       FORMAT(' FOR 100 CONSECUTIVE TIME STEPS HAVE CONVERGENCE ON'
     X      ' ITER = 1'/' STEADY STATE IS REACHED, NO FURTHER FLOW'
     X      ' CALCULATION')
           GOTO 499
         END IF
c
      END IF
c
c-------------------------------------------------------------------------*
c      Move data saving here from above, and include chemistry
c      saving. At this point flow has converged, thus KON test
c      is not needed (KON=1 is no flow convergence; KON=2 for convergence)
c      
       if(mod(kcyc,ntsave).eq.0) then
        kcyc_last=kcyc
        call wrifi
c
C-------------------Write data to restart.dat for reactive transport
        IF (MOPR(1).ne.2.and.mopr(1).ne.1) THEN
          OPEN (UNIT=21,FILE='savechem',STATUS='UNKNOWN')
          CALL WRITE_RESTART
          CLOSE (UNIT=21)
        END IF
C-------------------------------------------------------------------
       end if
c--------------------------------------------------------------------------*
c
      IF(KON.EQ.2)  GOTO 30  ! go to the next time step (KCYC=KCYC+1)
C+++++++++++++++++++++++++End   addition for solute transport
C
      IF(MOP(15).EQ.0) RETURN
      REWIND 8
      write(34,13)
   13 FORMAT(' WRITE COEFFICIENTS FOR SEMI-ANALYTICAL HEAT LOSS CALCULA'
     X'TION ONTO FILE *TABLE*')
      WRITE(8,14) (AI(N),N=1,NELA)
c      REWIND 111
c      WRITE(111,14) (AI2(N),N=1,NELA)
c      REWIND 112
c      WRITE(112,14) (AI3(N),N=1,NELA)
   14 FORMAT(4E20.13)
C
699   CONTINUE
C
C*************************************************************************
C
C-----WRITE RESULTS ON FILE SAVE.
   10 CONTINUE
      KCYC=KCYC-1
      IF (KCYC.EQ.KCYC_LAST) GOTO 39
      CALL WRIFI
C
C*************************Write data to restart.dat for reactive transport
      IF (MOPR(1).ne.2) THEN
C         CALL WRITE_PLOT
         OPEN (UNIT=21,FILE='savechem',STATUS='UNKNOWN')
         CALL WRITE_RESTART
         CLOSE (UNIT=21)
      END IF
C
C*************************************************************************
C
39    CONTINUE
!
!----------
!.....Write water chemistry in format same as chemical.inp to runlog.out file
!.....for selected grid blocks
!----------
!
      IF (MOPR(1).ne.2) then
       if(mopr(1).ne.1.and.mopr(1).ne.3.and.mopr(1).ne.5)call waterchem
         write(39,"(/2x,'...done!   simulation finished',/)")
         write(32,"(/2x,'...done!   simulation finished',/)")
         close (39)  !unit 39 is the iteration file opened in readsolu
         close (32)  !unit 33 is the run log file opened in init
      end if
      write(*,"(/2x,'...done! TOUGHREACT v2.0 simulation finished',/)")
      RETURN
      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE Evaluate_DOLD
c
c***********************************************************************
c***********************************************************************
c*                                                                     *
c*     Evaluate mass and energy accumulation terms for the previous    *
c*        time step of flow simulation. This subroutine is called      *
c         when coupled with reactive transport (porosity changes)      *
c*                                                                     *
c*                  Version 1.0 - July 18, 2006                        *
c*                                                                     *
c***********************************************************************
c***********************************************************************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
C
      INCLUDE 'flowpar_v2.inc'
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)
C
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
      COMMON/P5/DOLD(MNEQ*MNEL)
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
C
C$$$$$$$$$ COMMON BLOCKS FOR ROCK PROPERTIES $$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE6/SLOLD(MNEL)          ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)          ! old gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)        ! porosity at previous time step
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/BC/NELA
      CHARACTER*5 ELEM,ELEM1,ELEM2,MAT,ELEG
C
C######### LOCAL ARRAYS ###############################################
C
      DIMENSION D(11,12)
c
      SAVE ICALL
      DATA ICALL/0/
c
      ICALL=ICALL+1
c
      DO 100 N=1,NELA
C     IDENTIFY MATRIX BLOCK.
      NLOC=(N-1)*NEQ
      NLOCP=(N-1)*NK1
C     IDENTIFY THE START OF SECONDARY VARIABLES FOR ELEMENT N.
      NLOC2=(N-1)*NSEC*NEQ1
C
C     ASSIGN QUANTITIES WHICH DEPEND ONLY UPON ELEMENT INDEX N.
c
      PHIN=PHIOLD(N)
c
      NMAT=MATX(N)
C-------------------------------------------------------------------
C      CD=SH(NMAT)*DM(NMAT)*(1.-POR(NMAT))
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
        D(K,M)=0.D0
 1011 continue
C
C-----COMPUTE CHANGE IN POROSITY.
      DPRES=0.d0
      IF(M.EQ.2) DPRES=DELX(NLOCP+1)
      PRES=X(NLOCP+1)+DX(NLOCP+1)+DPRES
c      DPHI=PHIN*(COM(NMAT)*(PRES-P(N))+EXPAN(NMAT)*(PAR(NLM2+NSEC-1)
c     A-T(N)))
      PHINN=PHIN            !  +DPHI
C
      DO 102 NP=1,NPH
C     IDENTIFY BEGINNING OF SECONDARY VARIABLES FOR PHASE NP.
      NL2NP=NLM2+(NP-1)*NBK
C     SATURATION.
c
      SNP=PAR(NL2NP+1)
c
      IF(SNP.EQ.0.d0) GOTO 102
C     DENSITY.
      RHONP=PAR(NL2NP+4)
      PHISRO=PHINN*SNP*RHONP
C
C-----SUM OVER COMPONENTS IN EACH PHASE.
      DO 103 K=1,NK
C     MASS FRACTION OF COMPONENT K IN PHASE NP.
      XNPKM=PAR(NL2NP+NB+K)
      D(K,M)=D(K,M)+XNPKM*PHISRO
c
  103 continue
C
C
C     INTERNAL ENERGY IN PHASE NP.
      ENNP=PHISRO*PAR(NL2NP+5)-PHINN*SNP*PRES
      D(NK1,M)=D(NK1,M)+ENNP
C
  102 CONTINUE
C
C     ADD ROCK ENERGY.
      ENR=CD*PAR(NLM2+NSEC-1)
      D(NK1,M)=D(NK1,M)+ENR
C
C     ASSIGN ACCUMULATION TERM AT BEGINNING OF TIME STEP.
C
c      IF(ITER.EQ.1.AND.M.EQ.1) THEN
          DO 105 K=1,NEQ
             DOLD(NLOC+K)=D(K,1)
 105      CONTINUE
C
  101 CONTINUE

100   continue
c
      return
      end
c
      SUBROUTINE COUPLE(next_tstep)
c
C
C********* COUPLE MULTIPRIONENT SOLUTE TRANSPORT WITH CHEMICAL REACTION ********
C
C
C***** N O T A T I O N ********************
C
C     NEQ IS THE NUMBER OF EQUATIONS, AND THE NUMBER OF PRIMARY
C         DEPENDENT VARIABLES (PER ELEMENT).
C
C     NPH IS THE NUMBER OF PHASES WHICH CAN BE PRESENT.
C
C     NK  IS THE NUMBER OF COMPONENTS.
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
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      include 'perm_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      COMMON/XYZ33/ZZZ(mnel)
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E5/P(MNEL)
c      COMMON/E55/P_ini(MNEL)                               ! For Statoil
      COMMON/E6/T(MNEL)
C----------------------------------------------------- For using EOS9
      common/TEM_EOS9/Tc_EOS9(MNEL)  ! initial temperature (oC)
C--------------------------------------------------------------------
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
c
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
c     Add tortuosity exponent (ptort) and critical porosity (phicrit)
      common/torpar/ptort(maxmat),phicrit(maxmat)
c     Multiphase tortuosity at each grid block
      common/tortmp/tortliq(mnel),tortgas(mnel)
c
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C5/AREA(MNCON)
C
C$$$$$$$$$ COMMON BLOCKS FOR RESIDUALS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL*NEQ
C
      COMMON/P4/R(MNEQ*MNEL+1)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/PARNP/NPL,NPG          ! specify in EOS module
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR LINEAR EQUATIONS $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/L4/WKAREA(MNEQ*MNEL+10)
      COMMON/L7/JVECT(niwork)
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
c
      COMMON/LDIM/LICN,LIRN
c
      COMMON/NEQUAS/NEQC                         ! number of equations
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT IN LIQUID PHASE $$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)  ! darcy velocity
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL)
      COMMON/SOLUTE6/SLOLD(MNEL)         ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)         ! old gas saturation
      COMMON/SOLUTE8/SL1(MNEL)           ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)           ! new gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)       ! porosity at previous time step
      COMMON/PRINTC/NOW                  ! print control
      COMMON/AMMISC/IABC,ISOLVC
c
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
C
      PARAMETER (mnz1=(mnel+2*mncon)*2)
      COMMON/LL1/IRNO(mnz1)
      COMMON/LL2/ICNO(mnz1)
      COMMON/LL3/COO(mnz1)    ! solving solute transport
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G4/ELEG(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G12/LCOM(MNOGN)
C
C$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT IN BOTH LIQUID AND GAS PHASES $$$$$$$$$$
C
      COMMON/TRANGAS1/PFUGOLD(NMNOD,NMGAS)  ! old partial pressure
      COMMON/TRANGAS2/PFUGOLD2(NMNOD)
      COMMON/TRANGAS4/RHANDG(NMNOD)         ! right-han side for a gaseous species
      COMMON/TRANGAS5/RSOURCEG(NMNOD,NMGAS) ! gas source terms due to reactions
      COMMON/TRANGAS6/EKGAS(NMGAS)          ! Gas equilibrium constants
      COMMON/TRANGAS7/EKGAS2,PFUGB2(100)         
      COMMON/TRANGAS8/dcfgas(mgas,mnel),DIFUNG  ! gaseous species diffusion coefficient
      COMMON/TRANGAS9/NGAS1                     ! number of gaseous species
      COMMON/GASCONS1/GP0(NMNOD,NMGAS)          ! initial gas conc.
C------------------- Index of gas species for use in matrixg
      common/gasindx/kgas
C
C-----For co2 generation by mineral phase using EOS2 module
      common/co2_gene/nco2
      COMMON/ICO2/ICO2H2O  ! CO2 and H2O reaction sources considered in the flow
c
c----------------------------------------Indicators from EOS module
C
      COMMON/EOS_INDICATOR/ IEOS          ! Indicate EOS module used
      COMMON/CO2M_TMVOC/ ico2m, iTMVOC    ! Indicate EOS module used
!
C----------------------------------------------------For using EOS2 Flow module
C----------------------------------And calculating co2 gas fugacity coefficient
c
      common/co2_gene1/nco2g
      COMMON/REACTh2o/Rh2o(NMNOD)      ! H2O REACTION SOURCES
      COMMON/REACTco2/Rco2(NMNOD)      ! CO2 REACTION SOURCES
      COMMON/SOLIDco2/SMco2(NMNOD)     ! CO2 TRAPPED in solid phase
!
!.....Extract fugacity coefficients from TMgas
!
      common/fugacity_coe/fug_coe(mnel,18)
      common/gas_index/ichem(18)            ! No-cond gas index in chemical input
!
!.....Extracting CO2 fugacity coefficient from ECO2N
!
      common/fuga_coe /FugCoeCO2(mnel)
!
!----------
!.....For H2 generation by mineral phase using EOS5 module
!----------
!
      common/h2_gene/ nh2
      common/h2_gene1/ nh2g
      common/h2_gene2/ ih2gt0     !=1: initial Ph2>0
!
!----------
!.......For gas reaction sources for TMgas module
!----------
!
      common/gas_gene1/ nch4, nh2s, nso2
      common/gas_gene2/ Rgases(MNEL,6)
      common/co2_gene2/ ico2gt0        !=1: initial Pco2>0

C----------------CO2 partial pressure    For use Steve's CO2 module
      COMMON/PCO2_ALL/PCO2A(MNEL)      ! calculated from ECO2
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
!
!.....Use aqueous complex to model surface complexation
      common/ns_index/ ns
!
c-----------------------------Addition for chemical qusi-steady state
      common/minkin4/rkin(mnod,mmin),amin3(mnod,mmin)
c
      COMMON/EQUDISO/DELTAP(MNOD,MMIN),RKIN0(MNOD,MMIN)
      COMMON/DISRATE/DRATE(MNOD,MMIN)
      COMMON/STEADY/IFLOWSS,JSTEADY
      COMMON/TOL_STEADY/TOLDC,TOLDR         ! Concentration and dis/pre changes
      COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
c
C-------------------------------------------for coupling with TRANQUI
C
      common/isarea/imflg2(mmin),imflag(mnel,mmin)
c... Dissolution kinetics
      common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
c... Precipitation kinetics
c        modify this common block
      common/iprkin/ideprec(mmin)
c----------- Added common block for rate law designations
      common/irtlaw/nplaw(mmin)
c
c            porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
      common/satgas2/sg2
      integer*8 ielem
c----------- Next_tstep variable defined in routine cycit
      integer*8 next_tstep
      integer*8 ICALLCH(MNOD)
      double precision  UTEM(MNOD,MPRI)
      COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
      common/oldtemp/tcold(mnod),tcmix(mnod)  ! tcmix(mnod): temperature of dryingout node after mixing of brine
c
C---------------------------------------Interface area reduction factor
      common/afactor/a_fm2(mncon)  ! advection area reduction (flow from F to M)
      common/afactord/a_fmd(mncon) ! diffusion area reduction (Both sides)
c
c.... Save modified active fracture area for reaction --
c     limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      COMMON/KC/KC
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/BIND/DIFF0,TEXP,BE
      COMMON/DG/WUP,WNR
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/DFM/TIMAX,REDLT
      CHARACTER*5 ELEM,ELEM1,ELEM2,MAT,ELEG
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
C
C-----------------------Common blocks for dryout grid blocks
C
      COMMON/DRYOUT/IDRY(MNOD),ADRY(MNOD,MPRI)
      COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)    ! residual in precipitates
      common/dry_salt/nsalt,isalt(mmin)
c
c---------------------------COMMON blocks for Kd adsorption and decay
      common/kddca3/kddp(mpri)    ! pointer to the primary species
c--------------------------------------------------------------------
      common/dtlim/max_chem_it,delt_conne,id_chem
      character*16 delt_conne
      character*5 id_chem
c
      double precision rerror(mpri)
      double precision prev_c(mnod,mpri)   ! concentration at previous iteration
      double precision rsourceq(mnod,mpri) ! eq. minerals/gases source terms
      integer*8 ielr
      double precision densw,deltat0
c
      double precision prev_cmg(mnod,mgas) ! gas concentration (mol/kgw) at previous iteration
      integer*8 inode
c
c-----When porosity close to a minimum value, only dissolution is allowed
      common/dispre/idispre(mmin)
      common/pseduo/ps_ads,con
      common/aqkin16/NoTrans(mpri)    ! >0: not subject to transport
      COMMON/tot_solid_aq/icon_nod(mnod,mpri),ttt(mpri),
     &                ttt_nod(mnod,mpri)          ! total concentraion including both aqueous and solid
!
!.....Extract rock density for geochemical calculations such as exchange and sorption
      common/rock_density1/denss(mnel)
      common/rock_density2/denss2, sl2, phisl2, a_fmr2
!
!....................................................................................
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) THEN
        WRITE(34,899)
  899   FORMAT(/' *COUPLE q2.1, 1999.4.14: COUPLE TRANSPORT WITH'
     X  ' REACTION (els and ns changes)**********'/)
      END IF
      IF(KCYC.EQ.1)  THEN
         NCOUNT=0
      END IF

!----------------------------------------------------
      IF (ICALL.EQ.1) THEN  ! initialization for t=0
         ico2m = 0                                         
         if (ieos.eq.2 .or. ieos.eq.13 .or. ieos.eq.14)   ico2m = 1
         iTMVOC = 0
         if (ieos.eq.15 .or. ieos.eq.16) iTMVOC = 1
         max_chem_it = 0    ! Max chemical iterations for entire printout interval
C
C************************************ Addition for quasi-steady state
C
         DO I=1,NNOD
           DO M=1,NMEQU
             DELTAP(I,M)=0.0D0
           END DO
           DO M=1,NMIN
             RKIN0(I,M)=0.0D0
             RKIN(I,M)=0.0D0
           END DO
         END DO
C
C****************************************************************
C
         if(next_tstep.ne.2) then
            IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14
     &     .OR. IEOS .EQ. 15 .OR. IEOS .EQ. 16)   THEN    ! for ECO2 module)
               CALL WRITE_PLOT_ECO2
            ELSE
               CALL WRITE_PLOT
            END IF
c
            IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14
     &     .OR. IEOS .EQ. 15 .OR. IEOS .EQ. 16)   THEN    ! for ECO2 module)
               CALL  WRITE_TIME_ECO2
            ELSE
               CALL  WRITE_TIME
            END IF
         end if
         IF (NGAS1 .GT. 0)  THEN
c
c GP is gas amount in moles per liter medium.  Here we define the initial GP
c based on PFUG, the gas partial pressure in bars (starting PFUG is input
c in routine INIT and initially assumed=fugacity (at 1 bar total pressure)
c Note: at high P, conversion from pfug to gp below should include Z (compr factor) but we
c  neglect it in the back-conversion from gp to pfug (in newtoneq) and P does not
c  change during chem calculations, so we can ignore it as far as calculating pfug
c  (but gp is sytematically off)
            DO I=1,NNOD
              rt=gc*tkelv(i)   ! gas RT
              phisgrt = PHIOLD(I)*(1.d0-SLOLD(I))/rt
              phioslo= phiold(i)*slold(i)
              DO M=1,NGAS
                  GP(I,M)=PFUG(I,M)*phisgrt      ! gas conc in mol_gas/l medium
c
                  prev_cmg(i,m)=gp(i,m)/phioslo  ! gas conc in mol/kgw
              END DO
            END DO
         END IF
C
         DO 721 I=1,NNOD
c
c..............initialize old temperature to current
            tcold(i)=tc(i)
C
c          Aqueous species source terms due to mineral precip/dissol
c          These are in moles per liter liquid
c          rsource: current source terms due to kin. and eq. mineral and gas ppt/dissolution
c          rscold: previous time step source term for aq. transport (all mins + gas)
c          utem: ?? not used in this version
           DO 722 K=1,NPRI
              RSOURCE(I,K)=0.D0  !current aq
              rsourceq(i,k)=0.d0 !ns98/10 added
              UTEM(I,K)=UTOLD(I,K)
722        CONTINUE
C
721      CONTINUE
c
c  RSOURCEG: gas source terms due to gas exolution/dissolution
c  These are in moles per liter liquid
         IF (NGAS1.EQ.0) GOTO 299
C
         DO 817 I=1,NNOD
            DO 822 K=1,NGAS
c-------- rsoldg: previous time step source term for gas transport
              RSOURCEG(I,K)=0.D0
822         CONTINUE
817      CONTINUE
c
299      CONTINUE
c
      END IF
c
C-------end first-time-through-only initialization-------------------
C
c------- Set timetot here equal to SUMTIM to avoid potential for mixup.
c See comments in CYCIT about TIMETOT.  It stores total current time for
c printout during chemistry calcs below.  SUMTIM is a tough2 variable and
c stores the same data (updated after convergence of flow by call to CONVER
c in CYCIT.
      IF (IFLOWSS .eq. 0)        then
            timetot=sumtim
      else
            TIMETOT=TIMETOT+DELTEX    ! for steady-state flow
            SUMTIM=TIMETOT
      end if
c
c-----------------------------------------------------------------------
c--Now we proceed to initializations at the beginning of each new time step
c-----------------------------------------------------------------------
c
c-------- skip these initializations if we are not at the next
c time step (next_tstep=0) (i.e. if we iterate between flow/transp/chem)
c
      if(next_tstep.eq.1.or.next_tstep.eq.2) then
! ----------
!.....Extract hydrological and thermal parameters from TOUGH including water density and phase 
!     saturations, temperature and CO2 partial pressure (for using EOS2, ECO2, and ECO2N modules)
! ----------
      CALL Get_TOUGH_HT_Variables
!     SUBROUTINE Get_TOUGH_HT_Variables, placed close to the bottom of treact_v2.f
! ----------
C-------------------------------Set zero for measuring number of iterations
      ITERTR   = 0        ! Counter of transport iterations
      MAXITCH  = 0        ! Maximum iterations of solving whole chemistry
      AVERITCH=0.d0
      COUNTCH=0.d0
!
         DO 737 I=1,NNOD
c
c..........reset icallch array for sequential partly-iterative approach(SPIA)
            icallch(I)= 1
c
c..........Add phi*sl1 and phi*sg1 to reduce calcs later
            phisl1(i) = phi(i)*sl1(i)
            phisg1(i) = phi(i)*sg1(i)
c
c------- Initialisation of mineral/gas source terms and precipitation
c terms at the beginning of each time step (at t+delta t).
c PRE, RSOURCE,GP were calculated at the last time step (at t) and stored in common
c blocks and unchanged to this point.  Here we save these values as PRE0 and
c RSCOLD, and GP0 because original values will be updated for the current time step.
c
c       source terms
c          aqueous species due to kin. and eq. minerals + gas ppt/dissol
            do n=1,npri
c initialize source terms to zero to avoid neg. concentrations at beginning
c of seq. iterations
               rsource(i,n)=0.d0
               rsourceq(i,n)=0.d0
            end do
c
c-------- gas species (for gas transport)
            do n=1,ngas
               rsourceg(i,n)=0.d0
            end do
c
c       precipitation terms.  PRE contains total moles of mineral per
c       liter medium, including starting rock composition
c       (note: initial pre's are set to pinit in routine init)
c
            do m=1,nmin
               pre0(i,m)=pre(i,m)
            end do
c
c       gas exolution terms GP (equivalent of pre for gases)
c       these are in moles/liter medium
c       (note: initial gp's are set in routine init from input
c        gas pressures)
c
            if(ngas.gt.0) then
               do m=1,ngas
                  gp0(i,m)=gp(i,m)
                  pfugold(i,m)=pfug(i,m)  !gas partial pressure in bars (not fugacity!!!)
               end do
            end if
!
!.........Calculate H2(aq)+H2(g) at beginning of the time step (for EOS5 flow module)
!
            if (ieos .eq. 5)   then
               rt = gc*(tc(i)+273.15d0)              ! Added (gas RT)
c             Gh2g= (PFUG(I,nh2g)/rt)*PHI(I)*SG1(I)     ! gas moles/l medium  ns98/3 modify gc
c             rco2(i) = Gh2g + c(i,nh2)*phi(i)*sl1(i)   ! In mol/l medium
               Gh2g= (PFUG(I,nh2g)/rt)*phisg1(i)     ! gas moles/l medium  ns98/3 modify gc
               rco2(i) = Gh2g + c(i,nh2)*phisl1(i)   ! In mol/l medium
!            Here use rco2 store H2
            end if
!
!
c       Total aqueous concentrations UT (at t+dt) and UTOLD (at t)
c       These are total aqueous concentrations in mole/l liquid.
c       UT will get overriden by transport.
c       Note: UTOLD was already set equal to UT at the end of this routine
c       for consistency with printout and other routines which use UTOLD
c       instead of UT.  For now, leave the loop below c'd out
c
c       Note: UT is also used temporarily further below to save the
c       complete transport right hand side in moles/l medium, but it
c       ends up storing the current total aq. concentration in
c       moles/liter liquid at the end of the transport calcs (do not get
c       confused!).
c
737      CONTINUE
c
c       Sedimentation - advects downwards sediments (pre's) and porewater (ut's, pfug's)
c       Note, need to do this here rather than after chemistry so that the correct info is stored
c       for restarts later on
         if(mopr(11).eq.1.and.vsedDF.gt.0.d0) call sedim_1D
c
      end if
c
c.... Calculate tortuosity for aqueous phase transport
      call tortcalc(nel,sl1,tortliq)
!
!-----------Construct coefficient matrix for solute transport (aq.species)-------------
      if (wtime .eq. 1.0d0)   then
!........When implicit time-weighting (wtime = 1.0d0), the matrix is simpler
         CALL MATRIXC_implicit
      else
!........For general time-weighting scheme
         CALL MATRIXC
      end if
         NZ0=NZ
      DO 97 INZ=1,NZ
         COO(INZ)=CO(INZ)
         IRNO(INZ)=IRN(INZ)
         ICNO(INZ)=ICN(INZ)
97    CONTINUE
C
C----------------Add reaction source/sink terms to the right-hand sides
C-----RSOURCE is the dissolution rate array
c
c              ********************************
c--------------Start of transport/geochem cycle-------------------
c              ********************************
c
c--------- to skip all transport
      if(isolvc.eq.0) go to 509
c
1100  ITERTR=ITERTR+1
c
      istop=0      ! counts number of chemical convergence failures
c
      DO 101 I=1,NEL
c       replaced by array        phisl1 = phi(i)*sl1(i)
         DO 100 K=1,NPRI
            if (k .eq. nw) goto 100      ! Do not move water
            if (k .eq. nd) goto 100      ! AD is not subject to transport
!
            IF(NoTrans(K) .ne. 0) GOTO 100   ! Not subject to transport
!
c   UT here in mole/l medium used for right hand side, before transport
c       !! Solution to transport is stored later in UT in moles/l liquid
c      Note: rhand was previously set in MATRIXC as phiold*slold*utold
           UT(I,K)=RHAND(I,K)
c
c       ut(i,k)=ut(i,k)+phisl1
        ut(i,k)=ut(i,k)+phisl1(i)*(rsource(i,k)+rsourceq(i,k))
c
100     CONTINUE
101   continue
C
      if(MOPR(3).ge.1) then
         write (32,*) '   ---rsource---'
         do nnn=1,nel
            write (32,'(7E12.4)') (rsource(nnn,ip),ip=1,npri)
         end do
c
         write (32,*) '   ---rhand in COUPLE---'
         do nnn=1,nel
            write (32,'(7E12.4)') (rhand(nnn,ip),ip=1,npri)
         end do
c
         write (32,*) '   ---ut in COUPLE before transport ---'
         do nnn=1,nel
            write (32,'(7E12.4)') (ut(nnn,ip),ip=1,npri)
         end do
      end if
c        ------------------------------------------------------
c--------Solve linear equations for aqueous transport (loop 120)
c        -------------------------------------------------------
!
      DO 120 K=1,NPRI
!
         if (k.EQ.ns)    then            ! Use aqueous complex to model surface complexation
            do i=1,nel
               ut(i,k) = utold(i,k)
            end do
            goto 120
         end if
!
c.........Total concentrations of surf complexes now updated within speciation loop
c         using amounts of mineral present
         if(k.eq.nw) goto 120           ! Not subject to transport
         if(NoTrans(K).ne.0) goto 120           ! Not subject to transport
c
         DO 140 I=1,NEL
c
c    UT here in mole/l medium used for right hand side.  Do not confuse
c    with solution to transport, which is stored later in UT
c    in moles/l liquid
c
            R(I)=UT(I,K)
            IF (EVOL(I).EQ.0.0D0) R(I)=0.0D0
c
c-------------------------------------------------
c
  140    CONTINUE
C
         NZ=NZ0
         IABC=0
         DO 98 INZ=1,NZ
            CO(INZ)=COO(INZ)
            IRN(INZ)=IRNO(INZ)
            ICN(INZ)=ICNO(INZ)
98       CONTINUE
c
C------------------------------Begin modification of species with Kd and decay
c
         IF (KDDP(K) .GT. 0) THEN
            CALL MATRIXC_KDD(K)
         END IF
c
c-----------------------------------------------------------------------------
c
!
         IF (NCON.EQ.0 .OR. NEL.EQ.2)    THEN
!
!...........Get direct solution without calling solver for one or 
!...........two grid blocks or without connections
!
            CALL Direct_Solution (nel, ncon, Mnz, mneq, mnel, nz, 
     &                            irn, icn, co, r)   ! 1D arrays   
!
!...........Calling iteration solver
         ELSE
            NEQC=NEL
            CALL LINEQC(neqc)               ! Iteration solver
         END IF
!
         DO 160 I=1,NEL
!
!........UT now contains new total concentrations in moles/liter liquid!
!........Reset ut to old concentration if dry
!
            if (sl1(i) .gt. sl1min)     then
               UT(I,K)=R(I)
!
               if (iCO2M .eq. 1)   then
                  sl_scaling = (0.25d0*slold(i)+0.75d0*sl1(i))/sl1(i)
                  ut(i,k) = ut(i,k)*sl_scaling 
               end if
!
            else
               ut(i,k) = utold(i,k)
            endif
!
  160    CONTINUE
  120 CONTINUE
c                  ****************************
c------------------End of aq. species transport--------------------
c                  ****************************
      if(MOPR(3).ge.1) then
         write(32,*) ' npri,naqx,nmin,ngas,nads,nexc,ne,nw,nh,NEQC'
         write(32,'(11I5)') npri,naqx,nmin,ngas,nads,nexc,ne,nw,nh,NEQC
c
         write(32,*) '------utold(i,n)----'
         do ii=1,nnod
            write(32,'(5e12.4)') (utold(ii,nnn),nnn=1,npri)
         end do
c
         write(32,*) '------ut(i,n)----'
         do ii=1,nnod
            write(32,'(5e12.4)') (ut(ii,nnn),nnn=1,npri)
         end do
c
         write(32,*) '------c(i,n)----'
         do ii=1,nnod
            write(32,'(5e12.4)') (c(ii,nnn),nnn=1,npri+naqx)
         end do
         if(nmin.ge.1) then
            write(32,*) '------pre(i,n)----'
            do ii=1,nnod
               write(32,'(5e12.4)') (pre(ii,nnn),nnn=1,nmin)
            end do
         end if
      end if
c
c         ****************************************************
c         Full evporation calculations for nodes that went dry
c         ****************************************************
c     Need this here, before gas transport such that new gas concentrations
c     after dryout (from evaporation speciation calcs) are transported
c
c                      **********************
c----------------------Start of gas transport ---------------------------
c                      **********************
c.....When fixing Pco2, give a large amount of gas to maintain
c
      IF(NGAS1.EQ.0) THEN
         DO M=1,NGAS
            DO I=1,NNOD
               GP(I,M)=1.0D+01
            END DO
         END DO
         GOTO 499  !skip  gas transport
c---------- skip gas transport (not gas conc) if all nodes are liquid saturated
      else
         slflag = 1.d0 - sl1min
         do i = 1, nnod
            slflag = min(sl1(i),slflag)
         enddo
         if(slflag.ge.(1.d0-sl1min)) goto 509
      END IF
c--------------------------------------------------------
!
      IF (IEOS.EQ.2 .OR. IEOS.EQ.13) GOTO 509    ! for EOS2 and ECO2  modules
      IF (IEOS.EQ.5 .OR. IEOS.EQ.14) GOTO 509    ! for EOS5 and ECO2N modules
      IF (IEOS.EQ.15.OR. IEOS.EQ.16) GOTO 509    ! for TMgas
!
c                      --------------------------------
c----------------------Loop 95 solves for gas transport
c                      --------------------------------
c... Calculate tortuosity for gas species transport
      call tortcalc(nel,sg1,tortgas)
c

      DO 95 K=1,NGAS             ! changes to ngas
c
         kgas = k
c
c        Save some gas specific variables
c---------- gas molecular wt and molecular diameter
         wtmolg = dmwgas(k)
         dmolsq = diamol(k)**2
c
         DO 115 I=1,NNOD
            pfugold2(i)=pfugold(i,k)
c
c not needed            phisl1= phi(i)*sl1(i)
c            prev_cmg(i,k)=gp(i,k)/phisl1      ! gas conc in mol/kgw
            prev_cmg(i,k)=gp(i,k)/phisl1(i)    ! gas conc in mol/kgw
c
c      Moved here and saved gas diffusion coefficient at each grid block
c---------- Calculate diffusion coefficient of gaseous species (tracer)
            if(difung.lt.0.d0)then
               tkk = tkelv(i)
               ptg = p(i)
               call gasdiffus(tkk,ptg,wtmolg,dmolsq,dcfg)
               dcfgas(k,i) = dcfg
            elseif(difung.gt.0.d0)then
c......Add P-T dependent formulation consistent with TOUGH2
c......T and P should be also passed here, Otherwise got zero values and was crushed 
c        dcfgas(k,i) = difung*(1.01325d6/ptg)*(tkk/273.15d0)**texp
!
               tkk = tkelv(i)
               ptg = p(i)
               dcfgas(k,i) = difung*(1.01325d5/ptg)*(tkk/273.15d0)**texp
!                           p0=1.0d5 pa was used in TOUGH2
!
            else
               dcfgas(k,i) = 0.d0
            endif
c
115      CONTINUE
!
         call get_normalized(nnod,pfugold2,rnorm,pfugmin)
!
!........Construction of coefficient matrix for transport of gas species
!
         CALL MATRIXG_implicit
C
C-----------------------------------------------Printout for debugging
C
         IF(MOPR(2).GE.1)   THEN
            IF(NOW.EQ.1)  THEN
               write (32,*)  '   iRN  ICN    CO-MATRIX for gaseous'
               do nnn=1,NZ
                  write (32,'(i6, i6,E12.4)') irn(nnn),icn(nnn),co(nnn)
               end do
c
               write (32,*)  '   right-hand side for gaseous species'
               do nnn=1,NEL
                  write (32,'(5E12.4)') rhandg(nnn)
               end do
            END IF
         END IF
c
C-------------------------------
C
         DO 125 I=1, NNOD
c--------- rhandg (moles/l_medium)was calculated in matrixg as function
c       of pfugold (bars).  rsourceg is in moles/l_liq.
c       To correctly use source gas terms, we need to first subtract
c       rscoldg (below) because it was already accounted for
c       in pfugold during previous chemistry calcs.
c
            r(i)=rhandg(i)+ phisl1(i)*rsourceg(i,k)
            if (ieos.eq.9) r(i)=rhandg(i)
            IF (EVOL(I) .EQ. 0.0D0)   R(I)=0.0D0
            IF (SG1(I) .le. SL1MIN)   R(I)=0.0D0
125      CONTINUE
C-------------------------------------------------------------------
         NEQC=NEL                          ! number of equations
         CALL LINEQC(neqc)                 ! Iteration solver

         call get_realized(nnod,r,rnorm,pfugmin)

         DO 135 I=1,NNOD
            PFUG(I,K)=R(I)
c---------- in case we iterate and pfug goes negative
            if(pfug(i,k).le.0.d0) pfug(i,k)=1.d-40
c
            IF(MOPR(2).GE.1)   THEN
               IF(NOW.EQ.1)  THEN
                  write(32,*) '    I,   sg1(i),phi(i),PFUG(I,1)'
                  write(32,'(I5,6E11.3)')  I,sg1(i),phi(i),PFUG(I,1)
               END IF
            END IF
c
135      CONTINUE
95    CONTINUE
c
509   CONTINUE     ! for EOS2 and ECO2
c                 ----------------------------------
c---------------- End of loop to solve gas transport-----------------
c                 ----------------------------------
c
C------------------Calculate the amount of gas present in medium
C
c Note: at high P, conversion from pfug to gp below should include Z (compr factor) but we
c  neglect it in the back-conversion from gp to pfug (in newtoneq) and P does not
c  change during chem calculations, so we can ignore it as far as calculating pfug
c  (but gp is sytematically off)
c
      DO I=1,NNOD
         rt = gc*tkelv(i)              ! added (gas RT)
         psrt = phisg1(i)/rt
         DO M=1,NGAS
            GP(I,M)= PFUG(I,M)*psrt  ! gas moles/l medium  
         END DO
      END DO
c                      ********************
c----------------------End of gas transport--------------------------
c                      ********************
499   CONTINUE
c
      ERRU  = 0.0D+00
      ERRU1 = 0.0D+00
      ERRK  = 0.0D+00
c
c         ***************************************************
c---------Start of geochemical calculation loop for each node (loop 1000) ------
c         ***************************************************   
      call DRY_MAP
c......Skips chemical computations
      if (mopr(1).eq.3.or.mopr(1).eq.5) then
         do i=1,nnod
            do j=1,npri
               ctot(i,j)=ut(i,j)/dwat(i)  ! ut in mol/l, ctot in mol/kg
            end do
         end do
         goto 1001
      end if
      DO 1000 I=1,NNOD ! loop through number of elements
         if(idry(i).ge.1) goto 1000    ! calculated in drychem
c
c---------- moved sg1 redefinition here, below 499 point
c    Redefine sg1 in case we did not go through matrixg, where sg1 is defined
c    (needed for case when no gas transport occurs (ngas1=0))
c  already defined       sg1(i) = 1.d0-sl1(i)
c
c  already defined         phisl1 = phi(i)*sl1(i)
c
c       Moved density, temperature, gas sat assigments here
c        and take density from common block
         dliq=dwat(i)               ! density in kg/l
         densw = dwat(i)*1000.d0    ! densw in kg/m3    dwat in kg/dm3
c        dliq=densw/1.d3           ! density in kg/l
         sg2 = sg1(i)
         tc2 = tc(i)                ! added for use in NEWTONEQ
         tk2 = tkelv(i)
         Pt = p(i)/1.0d+5           ! total pressure (in bar)
c.......Get the ref pressure = 1 bar below 100C
c       and water saturation p above 100 C
cc        if (tc2.le.100.d0 .or. mopr(17).eq.1) then
         if (tc2.le.100.d0) then
            p0bar=1.d0
         else
            call SAT(tc2,psat)
            p0bar=psat/1.d+5
         end if
c
c.......Save this since it is used in many places
c       need to go through code carefully for this        eadum = ((1.0d0/tk2) - (1.0d0/298.15d0))/gc
c
         if (kcpl.EQ.2) then        ! only monitoring porosity change
            phii=phim(i)           ! phim(i) store the current porosity
         else
            phii=phi(i)
         end if
!
! ------------
!.......Get grid rock density for geochemical calculations such as exchange
! ------------
!
         a_fmr2 = a_fmr(i)        ! Modified active fracture area for reaction (not used in this version)
         phi2   = phii            ! Porosity
         sl2    = sl1(i)          ! Liquid saturation
         phisl2 = phisl1(i)       ! phi2*sl2
         denss2 = denss(i)        ! Density of solid rock (kg/dm^3)
         slosl1 = slold(i)/sl1(i) ! Saturation scaling
!
! ------------
!
c---- Change to skip chemistry if liquid saturation goes below a given value
c     also skip if porosity is nearly zero
c
         if(sl1(i).le.sl1min.OR.EVOL(I).LE.1.0D-15.or.
     +      EVOL(I).ge.1.d20.or.phii.lt.1.d-10) go to 999
         IF(ISPIA.EQ.1 .AND. icallch(I).EQ.0) GO TO 999  ! the SPIA approach
c
C-----------The total amount of solute in the mass balance, including solute
C-------transported in solution, adsorbed as ionic exchange and precipitated
c      assign chemical system according to the base of primary
c      species previously known for the node
c
         if(nads.gt.0) then
c
c-------Calculates new total site concentrations for surface complexes
c         get the surface area (surfads) in m2/kg_h2o for current node
c         then computes site concentrations and trial values
            call ads_area (i,mnel,nmin,nsurf,m_index,pre,phi,
     &             sl1,a_fmr,densw,vmin,surfads,supadn)
            call surface_conc(i,densw)
         end if
         call assign
c     assign the chemical parameters of the current node I
c
c     ut is total solute in solution (modified by transport) in moles/l liquid
c     c is concentration of primary species in moles/l liquid
c     ut is total concentration (primary+secondary species) in moles/l liquid
c     tt is total moles in solution (for given volume of liquid)
c
         stion=0.d0        ! stoichiometric ionic strength
         do n=1,npri
            tt(n)=ut(i,n)  ! tt in moles (per liter liq)
            stion=stion+zsqi(n)*ut(i,n)
         end do
c
c      Add statements below to always reset water according to liquid saturation
c      (we ignore ut(i,n) even if changed by chemical reaction)
cc      tt(nw)=ut(i,nw)  !use this to consider water from chemical reac, but then need to feed back to flow
c
         tt(nw)=rmh2o*dwat(i)
         stion=0.5d0*stion
c
c       stimax defined in paramete.inc.  Maximum allowable stoichiom. ionic strength
c       Skip chemical speciation/reaction if stimax is exeeded
cpitz        if(mopr(17).eq.0.and.stion.gt.stimax) goto 999
         if(stion.gt.stimax) goto 999                                                        !!!!!!!!!!!!!!!!!!!!!!
c
c---------  If sequential iterations, we need to subtract source terms
c for equilibrium minerals so that total mass is conserved
c
         if (itertr.gt.1) then
            do n=1,npri
               tt(n) = tt(n) - rsourceq(i,n) - rsource(i,n)
            end do
         end if
c
c   gas contribution to total concentrations
c
         do m=1,ngas
            ncp=ncpg(m)
c             gpphisl1 = gP(i,m)/phisl1
            gpphisl1 = gP(i,m)/phisl1(i)
            do k=1,ncp
               n=icpg(m,k)
               tt(n)=tt(n)+stqg(m,k)*gpphisl1
c---------    If sequential iterations, we need to subtract source terms
c             for gases at equilibrium so that total mass is conserved
               if (itertr.gt.1) then
                  tt(n) = tt(n)-stqg(m,k)*rsourceg(i,m)
               end if
            end do
         end do
!
!.......Cation exchange contribution to total concentrations
!
         do isite=1,NXsites     ! Loop over multi-sites
c
c          slosl1 = slold(i)/sl1(i)            ! Move above
            do k=1,nexc
               xcadasl = xcads(i, isite, k)*slosl1
               do n=1,npri
                  tt(n) = tt(n) + stqx(k,n)*xcadasl
c                  tt(n) = tt(n) + stqx(k,n)*xcads(i, isite, k)*
c     &                        slold(i)/sl1(i)   ! Saturation scaling     
               end do
            end do
         end do
c
c----if icon(i) equal to 5, the total (both solid and aqueous phase) is know, so
c---- tt is replaced by ttt,
c
         do n=1,npaq            !npri          excludes surface species
            if (icon_nod(i,n).eq.5) tt(n) = ttt_nod(i,n)
         end do
c
c--------- Assigns guess concentration values for chemical
c NR iterations.  Added if-else case when node rewets.  In such case
c it is better to not to take the previous concentrations as
c guess since the node was dryed out.  Use UT values (total conc.)
         if (slold(i).le.sl1min) then
            do n=1,npaq       !npri    excludes surface species
               if (tt(n).le.1.d-30) then
                  cp(n)=1.d-30
                  prev_c(i,n)=1.d-30
               else
                  ccx=dmax1(ut(i,n),1.d-7)
                  cp(n)=ccx/dwat(i)
                  prev_c(i,n)=ccx/dwat(i)
               end if
            end do
            do n=1,naqx
               cs(n)=1.d-10
            end do
            cp(nw)=dwat(i)
         else
c
c---- Assigns old concentrations at current node as initial guess values
c     cp,cs,cm etc are in moles/kg water liq. c's are in moles/Kg h2o
            do n=1,npaq  !npri     excludes surface species (assigned in routine surface_conc)
               cp(n)=c(i,n)
c-------- save previous concentration (for fully iterative procedure)
               prev_c(i,n)=c(i,n)
            end do
            do n=1,naqx
               cs(n)=c(i,npri+n)
            end do
c   for water, cp(nw) is kilogram water
            cp(nw)=dwat(i)          !use current water density
         end if
c   adsorption contribution to total concentrations
         do m=1,nads
            ncp=ncpad(m)
            do k=1,ncp
               n=icpad(m,k)
c           skips total for primary surface species
c           the total for these species comes from the surface area & site density data computed earlier
               if (n.le.npaq)
     &            tt(n)=tt(n)+stqd(m,k)*d(i,m)*cp(nw)
            end do
         end do
c
c+++++++ Calculate mineral reactive surface areas
c... Call routine to calculate surface areas in m^2/kg H2O
c
         ielr = i
c
         call rsfarea(ielr,densw)
c
         do m=1,nmin
c         cm is moles mineral (per liter liquid); pre in moles per liter medium
c         now cm is only incremental change
            cm(m)=0.d0
         end do
c
c       Save area (amin3)
         do m = 1, nmkin
            amin3(i,m) = amin2(m)
         end do
c
c   Assigns old gas moles and gas fugacity values at current node
c   as guess values
c
c       moved block to define sg1 further up
         rt = gc*tk2              ! added (gas RT)
c
         fact_g = sl1(i)/sg1(i)*rt
c
         do m=1,ngas
c        cmg is moles gas (per liter liq); gP in moles per liter medium
c        cg in bars is gas partial pressure
!
            if (iCO2m.eq.1 .or. iTMVOC.eq.1) then
               cmg(m) = gP(i,m)/phisl2
               cg(m)  = pfug(i,m)          ! now we can use this again
            else
               cmg(m) = prev_cmg(i,m)      ! we solve for this in newtoneq
               cg(m)  = cmg(m)*fact_g      ! this is used in cmq_cp for sat index
            end if
         end do
!---------------------------------------------- Gas fugacity coefficients
! moved up   Pt = p(i)/1.0d+5       ! Total pressure (in bar)
         if (ieos .eq. 9)   then
            Pt = 0.0d0
            do ig=1,ngas
               Pt = Pt + cg(ig)     ! Total pressure of gaseous species
            end do
         end if
!........Move all fugacity coefficient stuff into separate routine
         inode=i
         call fugacomp(inode)
!------------------------------------------------------------------------------
!.......Multi-site exchanges
         do isite=1,NXsites     ! Loop over multi-sites
            do k=1,nexc
               cxM(isite,k) = xcads(i, isite, k)
            end do
            cecM(isite) = cec(i, isite)
         end do
         phi2=phi(i)
         if(MOPR(3).ge.1) then
            write(32,1267) I
1267        format('grid block',I5)
            write(32,1268) (cp(nnn),nnn=1,npri)
1268        format(' --- cp ---- ',5e12.4)
            write(32,1269) (U2(nnn),nnn=1,npri)
1269        format(' --- U2 ---- ',5e12.4)
            write(32,1270) (cs(nnn),nnn=1,naqx)
1270        format(' --- cs ---- ',e12.4)
         end if
C----------------------------------Call geochemical subroutine for each node
         ielem = i
         izero = 0
         CALL NEWTONEQ(ielem,densw)

c--------- abort chemistry if ionic str became too high in newtoneq
c       flag is ielem was set to zero
         if(ielem.eq.0) goto 999
c
C----------------------------------------- ****************
C
         IF (ITERCH.GT.MAXITCH) MAXITCH=ITERCH
c.......Saves id of grid block with max chem iterations
         if (iterch .gt. max_chem_it) then
            id_chem=elem(ielem)
            max_chem_it=iterch    ! max for entire printout interval
         end if
c
         AVERITCH=AVERITCH+DBLE(ITERCH)
         COUNTCH = COUNTCH + 1.d0
         if (ichdump.eq.1) call chdump(timetot,ielem,iterch)  ! for testing only
         if (ichdump.eq.2) then  ! Printing speciation for specified grid blocks
            if (mod(kcyc,nwti).eq.0)   then
               iprint=0
            do ino=1,nwnod
               nng=iwnod(ino)
               if (nng.eq.ielem)   then
                  iprint=1
                  goto 3099
               end if
            end do
3099        continue
            if (iprint.eq.1) call chdump(timetot,ielem,iterch)
            end if
         end if
C
         if (iterch.ge.maxitpch) call chdump(timetot,ielem,iterch)
         IF(IRETURN .gt. 0)  then  !case when chemical convergence not reached
            write(32,"(' Time (seconds) = ',e10.4)") timetot ! and istop scheme below
            istop = istop + 1
            if (ireturn.eq.1)then
               ireturn = 0
               goto 999
            end if
            ireturn = 0
         end if
c
C
c     Now we need to convert concentrations from chem module:
c     cp's are in moles/kg h2o (molal) (per original 1 liter liquid fed into chem. module)
c     cm, cmg, rkin etc are in moles (per original 1 liter liquid fed into chem. module)
c     cp(nw) contains kg of water liquid (per original 1 liter liquid fed into chem. module)
c     vliq is new volume of liquid from the original 1 liter we used to run the
c     chemical module.  For now leave to 1 until we implement a way to
c     add/remove water in tough if vliq changes from 1
c
         sumsalts=0.d0   ! sum of salts weights in kg per kg water (assume zero for now)
         dliq=densw/1000.d0  ! liquid density in g/cc (kg/l)
         vliq=1.d0
         factw=cp(nw)/vliq   ! conversion factor = kg h2o liq/liter liquid
c
c-------------------------------------------------------------
c Source terms due to prec-dis. and ion exch. (current node i)
c-------------------------------------------------------------
c note: rsource and rsourceg are in moles per liter liquid
c  calculate source terms only for iterative approach
c  so we keep consistent with previous version.  All source terms
c  are assumed to remain zero if non-iterative approach is used
c
         if (ispia.ne.2) then
c
c      Initialize precip/diss sources to zero before updating them
c      add the two loops below
            do n=1,npri
               rsource(i,n) = 0.d0   !aq.species sources due to kin. and eq. minerals + gases
               rsourceq(i,n)= 0.d0   !aq.species sources due to eq. minerals+gases
c        (include all minerals and gases in rsource)
            end do
            do n=1,ngas
               rsourceg(i,n) = 0.d0  !gas sources (due to gas) (at equil.)
            end do
c
c       Equilibrium mineral contribution
            do m=1,nmequ
c          now cm is the incremental change
               dum=cm(m)/vliq
               ncp=ncpm(m)
               do k=1,ncp
                  n=icpm(m,k)
                  rsourceq(i,n)=rsourceq(i,n)-stqm(m,k)*dum
               end do
            end do
c
c       Kinetic mineral contribution
c
            dxvl = deltex/vliq
            do m=1,nmkin
               nkkn = m + nmequ
               ncp=ncpm(nkkn)
               rkdxvl = rkin2(m)*dxvl
               do k=1,ncp
                  n=icpm(nkkn,k)
                  rsource(i,n)=rsource(i,n)+stqm(nkkn,k)*rkdxvl
               end do
c
c..........Addition for quasi-steady state
               rkin(i,m)=rkin2(m)
            end do
c
c       Gas contribution
c
            do m=1,ngas
               if (phisl1(i).gt.0.d0)then
                  dum=cmg(m)/vliq-gP(i,m)/phisl1(i)
               else
                  dum=cmg(m)/vliq
               end if
               rsourceg(i,m) = rsourceg(i,m) + dum   ! gas source
               ncp=ncpg(m)
               do k=1,ncp
                  n=icpg(m,k)
                  rsourceq(i,n)=rsourceq(i,n)-stqg(m,k)*dum
               end do
            end do
c
c       Adsorption contribution
c
            do m=1,nads
               dum=cd(m)/vliq-d(i,m)     
               ncp=ncpad(m)
               do k=1,ncp
                  n=icpad(m,k)
                  rsourceq(i,n)=rsourceq(i,n)-stqd(m,k)*dum
               end do
            end do
c
         endif
!
!
! ----------
!........CALL subroutine to obtain CO2 reaction source terms for feeding 
!........back to flow equations (only for EOS2, ECO2, and ECO2N flow modules)
! ----------
!
!
!        Fkg = phi(i)*sl1(i)*dliq      
        Fkg = phisl1(i)*dliq      !!!!!! Also used for other modules
!
!
          if (ico2m .eq. 1)   then
!
             CALL CO2H2O_ReactionSource (deltex, Fkg, nmequ, nmkin, 
     &               Mpri, Mmin, nco2, nw, cm, ncpm, icpm, rkin2, 
     &               stqm,  Rco2_i, Rh2o_i)
!
            Rco2 (i) = Rco2_i 
            Rh2o (i) = Rh2o_i
!
!...........Amount of CO2 trapped in minerals phases
!
            SMco2(i) = SMco2(i) - Rco2_i*deltex      
!
          end if
!
!
! ----------
!........CALL subroutine to obtain H2O reaction source terms for feeding 
!........back to flow equations (only for EOS1, EOS3, EOS4, and EOS7)   
! ----------
!
         ih2om = 0
         if (icall .eq. 1)   then     ! This block can be moved at the beginning
            ih2om = 0
            if (ieos.eq.1 .or. ieos.eq.3  .or. ieos.eq.4 .or.
     &         ieos.eq.7)   ih2om = 1
         end if
!
         if (ih2om .eq. 1 .and. Ico2h2o .gt. 0) then   ! This module not yet incorporated, QLOSS needs modified
!
            CALL H2O_ReactionSource (deltex, Fkg, nmequ, nmkin, Mpri, 
     &               Mmin, nw, cm, ncpm, icpm, rkin2, stqm, Rh2o_i)
!
            Rh2o (i) = Rh2o_i
!
         end if
!
! ----------
!........CALL subroutine to obtain H2 reaction source terms for feeding 
!........back to flow equations (only for EOS5 flow module)
! ----------
!
         if (ieos .eq. 5 .and. Ico2h2o .gt. 0) then       ! This module not yet incorporated
!
            CALL H2_ReactionSource (deltex,  Fkg, nmequ, nmkin, 
     &              Mpri, Mmin, namin, rkin2, Rco2_i)
!
            Rco2 (i) = Rco2_i         ! R_CO2 now store for H2
!
         end if
!
!
!---------
!.......Consider gas reaction sources (mass transfer with mineral phases
!...... for flow module TMVOCs
!---------
!
         if (iTMVOC .eq. 1)       then
!
!           Fkg=phi(i)*sl1(i)*dliq     ! Assigned before CO2 module
!
            do ig=1,ngas
!
               Rgases (i,ig) = 0.d0    ! Reaction source/sink terms
!
               Gmolw = dmwgas(ig)      ! Gas molecular weight, g/mol
!
               Ngas_basis = 0
               if (nagas(ig) .eq. 'CO2')   then
                  Ngas_basis = nco2    ! The related basis species number
               end if
               if (nagas(ig) .eq. 'CH4')  Ngas_basis = nch4
               if (nagas(ig) .eq. 'H2S')  Ngas_basis = nh2s
               if (nagas(ig) .eq. 'SO2')  Ngas_basis = nso2
               if (nagas(ig) .eq. 'O2')   Ngas_basis = no2aq  
               if (nagas(ig) .eq. 'H2')   Ngas_basis = nh2   
!
!..............Calculate gas reaction source/sink (mass transfer with mineral phases)
!
               if (Ngas_basis .gt. 0)   then
!
                  CALL Gas_ReactionSource
     &              (deltex,      ! Time step, s
     &               Ngas_basis,  ! The gas related basis species number
     &               Gmolw,       ! Gas molecular weight, g/mol
     &               Fkg,         ! Factor for unit converion (mol/kg h2o to mol/dm**3 medium)
     &               Rgas_i)      ! Reaction source/sink for grid block i, gas ig
!
                  Rgases(i,ig) = Rgas_i  ! Store the source terms in a 2-D array, in kg/m**3/s
!
               end if
!
               if (nagas(ig) .eq. 'CO2')   then
                  SMco2(i) = SMco2(i) - Rgas_i*deltex   ! CO2 trapped in solid phase
               end if
!
            end do  ! ig
!
!...........Calculate H2O reaction source/sink (mass transfer with mineral phases)
!
            Rh2o(i) = 0.0d0    
!
            if (Ico2h2o .eq. 2 .and. nw .gt. 0)  then
               Ngas_basis = nw
               Gmolw      = 18.0D0
               CALL Gas_ReactionSource (deltex,Ngas_basis,
     &                                  Gmolw, Fkg, Rgas_i)
               Rh2o(i) = Rgas_i
            end if
!
!
         end if
c------------------------------------------------------------
c   Assign the new concentration values to the current node
c------------------------------------------------------------
c       Change of units - save c's as molalities
c       cp's, cs's are in moles/kg h2o liq
         do n=1,npri
            c(i,n)=cp(n)
         end do
c       Note, here we save the actual moles of H2O (per Kg H2O)
         c(i,nw)=rmh2o
c
         do n=1,naqx
            c(i,npri+n)=cs(n)
         end do
         ph(i)=ph2
         phislv = phisl1(i)/vliq
         do m=1,nmequ
            pre(i,m)=pre0(i,m) + cm(m)*phislv
            if (pre(i,m).le.1.0d-30) pre(i,m)=0.0d0
         end do
c
         phisvd = phisl1(i)*factw*deltex
         do m=1,nmkin      ! for kinetic dissolution
            nkkn = nmequ+m
            pretmp=pre0(i,nkkn)-rkin2(m)*phisvd
            pre(i,nkkn) = pretmp
            if (pretmp.le.1.0d-30) pre(i,nkkn)=0.0d0
c
c----------- for minerals specified at equilibrium precip and kinetic
c  dissol move the amount calculated under equil to the amount calculated
c  under kinetics and reset pre of equil. mineral to zero
c
            if (kineq(m).ne.0) then
               pre(i,nkkn)=pre(i,nkkn) + pre(i,kineq(m))
               pre(i,kineq(m)) = 0.d0
            end if
c
c........Add for printout
            rkin(i,m)=rkin2(m)
c
         end do
c
         do m=1,ngas
            gP(i,m)=cmg(m)*phislv   !new amount of gas moles/L_medium
            if (gp(i,m).lt.0.0d0) gp(i,m)=0.0d0
c
            prev_cmg(i,m)=cmg(m)
C
            IF (ico2m.eq.1 .and. nagas(m).eq.'co2(g)') go to 799
            pfug(i,m)=cg(m)                     ! new partial pressure
799         continue
         end do
c
         if (npads.gt.0) then
            do k=1,nads
               d(i,k)=cd(k)     ! both d and cd are in mol/kgw
            end do
            do n=1,nsurf
               phip(i,n)= phip2(n)
            end do
         end if
!
!.......Multi-site exchanges
!
         do isite=1,NXsites     ! Loop over multi-sites
            do k=1,nexc
               xcads(i, isite, k) = cxM(isite,k)/vliq
            end do
         end do
c
c-----------------------------------------------------------------------------
c
c Here, we get residual error on concentrations (erru) for the fully-iterative
c transport-speciation procedure, for the current node i.  This is the maximum
c relative difference between current concentrations and those at the previous
c transport-speciation step (not the previous time step!).
c Use component concentrations rather than total concentrations to save
c the extra loops needed to compute total concentrations.
c
         if (maxitptr.gt.1) then     !skip if not fully iterative case
            do n=1,npri
               if (n.ne.nw) then
                  diff=dabs(prev_c(i,n)-c(i,n))
                  dum=dabs(prev_c(i,n)+c(i,n))*0.5d0
                  diff=diff/dum
                  rerror(n)=diff  !save the error for later
                  if (diff.gt.erru) then
                     erru=diff
                     nerror=n      !index of species with greatest error
                  end if
               end if
            end do
            if (erru .le. toltr) icallch(i)= 0    ! for the SPIA approach - flag to call chemistry
         end if
c
c note: if ISPIA=1, then chemical speciation will not be performed at those nodes were
c icallch is set to zero by the above statement.
C
c needed up here so that ut is updated only if chemistry is performed
c    Compute new ut's (total aqueous concentrations in mol/L soln)
         do n=1,npaq           ! npri  skip surf complexes
            ctot(i,n)=c(i,n)   ! ctot in mol/kgw
         end do
         do j=1,naqx
            ncp=ncps(j)
            do k=1,ncp
               n=icps(j,k)
               utemp=ctot(i,n)+stqs(j,k)*c(i,npri+j)
               ctot(i,n)=utemp        !total concentrations in mol/kg water
            end do
         end do
c........Need loop below to compute ut's in mol/L
         do n=1,npri
            ut(i,n)=ctot(i,n)*factw  ! total concentrations in mol/L sln (for transport equations)
         end do
c
c------------- need to reset source terms for gases to zero
c if we do sequential iteration but the node is dry
c
         goto 1000
 999     continue
         do n=1,ngas
           rsourceg(i,n) = 0.d0  !gas sources (due to gas) (at equil.)
         end do
1000  CONTINUE
1001  continue
c
c             *************************************
c-------------end of geochemical loop for each node---------------------
c             *************************************
c
c---------- move setting of phiold here (and only for the first iteration)
      if (itertr.eq.1)then
         do i = 1, nnod
            phiold(i) = phi(i)
         end do
      end if
c
c...... Add changes for porosity modification
c
      if(kcpl.GE.1) call phichg
      if(kcpl.GE.1) call permchg
!
!.......For monitoring porosity and perm. changes but not affecting flow
!
      if (kcpl.EQ.2) then
         do i = 1, nnod
            phim(i) = phi(i)                 ! save for printout
            phi(i) = phi0(i)
!
            do jiso=1,3
               permm(jiso,i) = perm(jiso,i)  ! save for printout
               perm(jiso,i) = perm0(jiso,i)
            end do
         end do
      end if
c----------------------------------------------------------------------------------------
c... Do Leverett scaling whether or not permeability is changed
      if(mopr(6).eq.1)call levscale
c
      IF(MAXITPTR .EQ. 1) GO TO 889         ! for SNIA approach (non-iterative)
c
c----------------------------------------------------
c  convergenc test for transport/chemistry iteration
c----------------------------------------------------
c
      IF(NMIN.EQ.0 .AND. NEXC.EQ.0 .AND. NADS.EQ.0) GO TO 889
      IF(ITERTR.gt.1.and.itertr.ge.MAXITPTR) THEN
         goto 889
      ENDIF
c
      if(erru.gt.toltr) then          ! if tolerance is not met
         if(ispia.lt.2) go to 1100   ! back to transport only
         if(ispia.eq.3) then
c             flag - do not increment time step and iterate with flow
            next_tstep=0
            return                    ! back to cycit to iterate with flow
         end if
      end if
c
 889  continue
c
c ---below this point we have converged iterations between
c    flow, transport, and chemistry (if selected) and we
c    proceed towards the next time step.
c
c             **********************************
c-------------End of transport/geochemical cycle---------------------
c             **********************************
c
C-------------------------------------------- Map dry-out grid block
c
c------------Calculate the amount of minerals formed from dry-out
c
      IF (NSALT .GT. 0)THEN
         CALL DRY_MIN
         do i = 1, nnod
            do m = 1, nmin
               pre(i,m) = pre(i,m) + drypre(i,m)
            end do
         end do
      END IF
c---------------------------------------------------------------------
      if(countch.ne.0.d0) then
         AVERITCH=AVERITCH/COUNTCH
      else
         averitch=0.d0
      end if
c
c------Flag to indicate we can move to next time step
c      and stop iterating with flow
      next_tstep=1
c
C**************************************************************************
c
      deltat0 = deltex
c
1962  CONTINUE
c
c---------- calculate porosity change and resulting effect on
c permeability and capillary pressure, if option is selected.
c
      if (kcpl.ge.1) call phichg
      if (kcpl.ge.1) call permchg
C
C------For monitoring porosity and perm. changes but not affecting flow
      if (kcpl.EQ.2) then
         do i = 1, nnod
            phim(i) = phi(i)        ! save for printout
            phi(i) = phi0(i)
c
            do jiso=1,3
               permm(jiso,i) = perm(jiso,i)  ! save for printout
               perm(jiso,i) = perm0(jiso,i)
            end do
         end do
      end if
c----------------------------------------------------------------------------------------
c
      if((kcpl.eq.1.or.kcpl.eq.3).and.mopr(6).eq.1) call levscale
c
      CALL COMPUTE_MASS(icall,deltat0)
      NWMAS=10*NWTI
      IF(MOPR(8).ge.4.and.MOD(KCYC,NWMAS).EQ.0)CALL WRITE_MASS
c
      DO 220 I=1,nnod
c
c        add just in case
         sg1(i) = 1.d0 - sl1(i)          
         SGOLD(I)=SG1(I)              ! gas saturation
         SLOLD(I)=SL1(I)              ! liquid saturation
         if(phi(i).gt.0.d0.and.sl1(i).lt.0.99999d0.and.kcpl.eq.3)then
            sl1(i) = dmin1((phiold(i)*sl1(i)/phi(i)),0.99999d0)
            sg1(i) = 1.d0-sl1(i)
            NLOC2=(i-1)*NSEC*NEQ1      ! start of sec. variables for N
            NP=NPL
            NL2NP=NLOC2+(NP-1)*NBK
            PAR(NL2NP+1) = sl1(i)
            NP=NPG
            NL2NP=NLOC2+(NP-1)*NBK
            PAR(NL2NP+1) = sg1(i)
         endif
c
         tcold(i)=tc(i)
c
         do n=1,npri
            utold(i,n)=ut(i,n)
            adryr0(i,n) = adryr(i,n)
         enddo
c
220   CONTINUE
c
      if(MOPR(3).ge.1) then
        write(32,*) 'new porosity'
        write(32,*) (phi(i),i=1,nnod)
        write(32,*) 'new permeability'
        write(32,*) (perm(1,i),i=1,nnod)
        write(32,*) 'Capillary Pressure factor'
        write(32,*) (pcfact(i),i=1,nnod)
      endif
c
C------------------------------------------------ Addition for qusi-steady state
C-----Relative concentration and dissolution rate change between two steps
C---------------------For monitoring attainment of quasi-stationary state
c
      IF(IFLOWSS.EQ.1 .AND. TOLDC.GT.0.0D0 .AND. TOLDR.GT.0.0D0) THEN
         DO I=1,NNOD
c
            if(sl1(i).le.sl1min.OR.EVOL(I).LE.1.0D-15.or.
     +         EVOL(I).gt.1.d20.or.phi(i).lt.1.d-10) goto 949
C----------------------------------------------------------------
c
            DO N=1,NPRI
               IF (N.EQ.NW)    GOTO 379
               IF (N.EQ.NH)    GOTO 379
               IF (N.EQ.No2aq) GOTO 379
               DIFF=DABS(UTOLD(I,N))-DABS(UTEM(I,N))
               DUM=DABS(UTEM(I,N))
               IF(DUM .LE. 1.0D-35) GOTO 379
               DIFF=DABS(DIFF/DUM)
               IF(DIFF .GT. ERRU1) ERRU1=DIFF
379            CONTINUE
               UTEM(I,N)=UTOLD(I,N)
            END DO
            DO N=1,NMKIN
C------------------------------------------------
               nmqp1 = NMEQU+N
               DELTAP(I,nmqp1)=PRE(I,nmqp1)-PRE0(I,nmqp1) ! precipitation in current time step
               RKIN(I,N)= -DELTAP(I,nmqp1)/DELTEX
C-------------------------------------------------
               DIFF=DABS(RKIN(I,N)-RKIN0(I,N))
               DUM=DABS(RKIN0(I,N))
               IF(DUM .LE. 1.0D-30) GOTO 389
               DIFF=DIFF/DUM
               IF(DIFF .GT. ERRK) ERRK=DIFF
389            CONTINUE
               RKIN0(I,N)=RKIN(I,N)
            END DO
            DO N=1,NMEQU
               N1=N+NMKIN
               DELTAP(I,N)=PRE(I,N)-PRE0(I,N) ! precipitation in current time step
               RKIN(I,N1)= -DELTAP(I,N)/DELTEX
               DIFF=DABS(RKIN(I,N1)-RKIN0(I,N1))
               DUM=DABS(RKIN0(I,N1))
               IF(DUM .LE. 1.0D-30) GOTO 489
                 DIFF=DIFF/DUM
               IF(DIFF .GT. ERRK) ERRK=DIFF
489            CONTINUE
               RKIN0(I,N1)=RKIN(I,N1)
            END DO
949         CONTINUE
         END DO
C
C----------------------Update mineral abunance when quasi steady-state is reached
C
         NUMS=3
         JSTEADY=0           ! chemical staedy-state flag
         IF(ERRU1 .GT. TOLDC  .OR. ERRK.GT.TOLDR)  THEN
            KITC=0
            TIME1=0.0D0        ! accumulative time
            ntmin = nmequ + nmkin
            DO I=1,NNOD
               DO M=1,ntmin
                  DRATE(I,M)=0.0D0    ! average dissolution rate
               END DO
            END DO
         END IF
c
         IF(ERRU1 .LE. TOLDC .AND. ERRK.LE.TOLDR)  THEN
            KITC=KITC+1
            TIME1=TIME1+DELTEX
            ntmin = nmequ + nmkin
            DO I=1,NNOD
               DO M=1,ntmin
                  DRATE(I,M)=DRATE(I,M)-DELTAP(I,M)
               END DO
            END DO
         END IF
c
         IF (KITC .GE. NUMS)     THEN
            JSTEADY=1
            NOW=1
            ntmin = nmequ + nmkin
            DO I=1,NNOD
               DO M=1,ntmin
                  DRATE(I,M)=DRATE(I,M)/TIME1  ! average dissolution rate
               END DO
            END DO
            CALL STEADYC(NNOD,NMEQU,NMKIN,TIMETOT)
         END IF
      END IF
C
C--------------------------- End monitoring attainment of quasi-stationary state
C
C
C-------------------------------------------------Write concentrations
c
      IF(NOW.EQ.1)  THEN
         IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14 .OR.
     &       IEOS .EQ. 15 .OR. IEOS .EQ. 16)     THEN    ! for ECO2 module)
            CALL WRITE_PLOT_ECO2
         ELSE
            CALL WRITE_PLOT
         END IF
      END IF
      IF(MOD(KCYC,NWTI).EQ.0.OR.NOW.EQ.1)  THEN
         IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14 .OR.
     &       IEOS .EQ. 15 .OR. IEOS .EQ. 16)     THEN    ! for ECO2 module)
            CALL  WRITE_TIME_ECO2
         ELSE
            CALL  WRITE_TIME
         END IF
c------------------------------
      END IF
C
      i8ten = 10
      IF(MOD(KCYC,i8ten).EQ.0 .OR. KCYC.LE.20
     +     .or.icall.le.20.or.NOW.eq.1)   THEN
         CALL WRITE_ITER (KCYC,TIMETOT)
      END IF
      RETURN
      END
!
      SUBROUTINE Modify_PoroPerm
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      include 'perm_v2.inc'
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*   Modify permeability and porosity according to "SEED" in flow.inp  *
!*                     Called at CYCIT                                 *
!*                                                                     *
!*                   Version 1.0 - August 15, 2006                     *
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/E4/PHI(MNEL)
      common/E7/pm(MNEL)
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
!--------
!.....Modify permeability
!--------
      if (mopr(15) .eq. 1)   then                 ! Option 1
         do n=1,nel
            do i=1,3
               perm(i,n) = perm(i,n)*pm(n)
            end do
            pm(n) = 1.0d0
         end do
      end if
!--------
!.....Calculate permeability from porosity
!--------
      if (mopr(15) .eq. 2)   then                ! Option 2
         do n=1,nel
            phi(n) = phi(n)*pm(n)
            dlog10k = -2.5d0 + 8.0d0*phi(n)   ! log10k = -2.5+8PHI
            xk     = 10.0d0**dlog10k          ! In darcy
            xk     = 1.0D-12*xk               ! m*2
            do i=1,2                          ! Kx = ky
               perm(i,n)     = xk
            end do
            perm(3,n) = 0.01d0*xk             ! Kz
            pm(n) = 1.0d0
         end do
      end if
      return
      end

