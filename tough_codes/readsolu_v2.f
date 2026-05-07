c
c
       SUBROUTINE READSOLU
c
C***************Read input data for solute transport part*************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      INCLUDE 'common_v2.inc'
      character*20 thermo_in,pitzdata  ! name of thermo database
      COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
      common/wricon2/ prt_com(mpri),prt_min(mmin),prt_aq(maqx),
     &                 prt_ads(mads),prt_exc(mexc)
      character*20 prt_com,prt_min,prt_aq,prt_ads,prt_exc
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
      COMMON/AMMISC/IABC,ISOLVC                    ! for solver
      COMMON/TRANGAS8/dcfgas(mgas,mnel),DIFUNG     ! gaseous species diffusion coefficient
      COMMON/TRANGAS9/NGAS1   ! Number of gaseous species for transport
C
      COMMON/WRICON1/ELEMW(200)
      COMMON/E1/ELEM(MNEL)
c     common block with routine readtherm to save
c     the name of the thermodynamic data base
      common/thermodat/thermo_in,pitzdata
c
      CHARACTER*5 EL,ELEM,ELEMW
      double precision TOLAD
c
      COMMON/TOL_STEADY/TOLDC,TOLDR    ! Concentration and dis/pre changes
c
c-----------------------------------COMMON blocks for Kd adsorption and decay
      common/Kddca5/izonekd(mnod) ! Kd zone code
c
c----------------------------------------------------For EOS2 or ECO2 modules
      COMMON/ICO2/ICO2H2O  ! CO2 and H2O reaction sources considered in the flow
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
C
C---------------------------------------------------------------------
      common/drmin/numdr    ! >0 numerical derivatives of mineral kinetics
c
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
      COMMON/STR_THRES/str_threshold  ! Ionic strength threshold for switch between pitzer and DH
cpitz      COMMON/brine_salt/tol_remove,rate_remove,rate_plus_1,nremove
!
      common/TDS_REACT1/ iTDS_REACT
!
      common/Print_Unit_Name/ Name_Conc, Name_Mine
      character*10 Name_Conc
      character*45 Name_Mine
      character*200 inprec
!
c----------------------------------------------------------------------------
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' READSOLU 1.0, 2008.2.6: Read input data for solute'
     X' transport')
c
C---------------------------------------------------------------------
c
      OPEN(UNIT=31,FILE='solute.inp',STATUS='OLD') 
      OPEN(UNIT=32,FILE='solute.out',STATUS='UNKNOWN')
c     open run log file
      open(unit=33,file='runlog.out',status='unknown')

      WRITE(*,*) '   --> reading solute transport input data'
      WRITE(32,*) ' ---> start reading solute.inp'
      WRITE(33,*) ' ---> start reading solute.inp'
c
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
C-------------------------------------------------Read the title card
      READ(inprec,"(a82)") chemtitle  ! 82 characters - defined in common_v2.inc
      WRITE(32,"(1x,a82)") chemtitle
      WRITE(32,15)
15    FORMAT(/' ****** Writing solute transport input ******'/)
C----------Read option variables for controlling chemical calculations
      READ(31,*)
      WRITE(32,352)
352   FORMAT(/' Options for controlling chemical calculations:'
     1        /' *******************************************')
c
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      READ(inprec,*,end=9999,err=9999) ISPIA,INIBOUND,ISOLVC,ngamm,
     &        NGAS1,ichdump,kcpl, Ico2h2o, numdr
c
!      numdr = 0
      iTDS_REACT = 0
c
c  ispia     Flag for chemical solver iterative approach:
c             =0 full iteration between transport and chemistry
c             =1 partial iteration between transport and chemistry (not all points)
c             =2 no iteration
c             =3 full iteration between flow, tranport, and chemistry
c  inibound  Flag for recharge water chemistry
c             =0 no
c             =1 yes
c  isolvc    Flag for solver type
c             =1 Default Biconjugent Gradient Squared Solver (5)
c             otherwise: iterative solver
c  ngamm     Number of iterations between activity coefficient calculations and chemical solver
c             =0 Default (1)
c             >0 
c  ngas1     Number of gases for transport calcs
c
c  ichdump  Flag: 0 = disabled
c                 1 = chemistry dumps at each time step, each node !!!!for test only!
c  kcpl     Flag: 0 = disabled
c                 1 = enables porosity and permeability changes due to mineral pre/dis
c                 2 = Monitor porosity and permeability changes but don't affect flow
c                 3 = enables porosity and permeability changes and recalculates liquid saturation
c
c  iTDS_REACT     1 = Pass TDS_REACT to TOUGH primary variable X
c  iTDS_REACT     2 = Pass TDS_REACT to density subroutine in ECO2N
c  numdr    Flag: 0 = calculate derivatives of minetl kinetic rates by analytical method
c                 1 =                                               by mnmerical  method
c
      IF (ISPIA.LT.0 .OR. ISPIA.GT.3) ISPIA=2  ! defaults to no iteration
!
!.....Always initialize boundary/injection waters
      INIBOUND = 1
!
      IF (ISOLVC.EQ.1 .OR. ISOLVC.GT.6) ISOLVC=3       ! Default
      IF (NGAS1.LT.0) NGAS1=0
C
      IF(ISPIA .EQ. 0) THEN
         WRITE(32,353)
353      FORMAT(' Using full sequential iterative approach between'
     X   ' transport and chemistry (not flow)')
      ELSE IF(ISPIA .EQ. 1) THEN
         WRITE(32,354)
354      FORMAT(' Using sequential partly-iterative approach between'
     X   ' transport and chemistry (not flow)')
      ELSE IF(ISPIA .EQ. 2) THEN
         WRITE(32,*) ' Using sequential non-iterative approach'
      else if(ispia .eq. 3) then
         write(32,*) ' Using full sequential iterative approach'
      END IF
C
      IF(INIBOUND .EQ. 0) THEN
         WRITE(32,356)
356   FORMAT(' not initializing boundary & recharge input solution')
      ELSE IF(inibound .EQ. 1) THEN
         WRITE(32,357)
357   FORMAT(' Initializing boundary & recharge input solution')
      END IF
C
      IF(ISOLVC .EQ. 1) THEN
         WRITE(32,358)
358      FORMAT(' Default BiConjugent Gradient Squared Solver (5)')
      ELSE
         WRITE(32,359)
359      FORMAT(' Iterative solver - See Tough2 V2 Manual')
      END IF
c     Added NGAMM here instead of RCOUR
      IF(ngamm .EQ. 0) THEN
         WRITE(32,460)
460      FORMAT(' Default 2 Iterations: DH and Chemical Solver')
      ELSE
         WRITE(32,461)
461      FORMAT(' Iterations: DH and Chemical Solver')
      END IF
      if(ngamm.le.0)ngamm = 1
         WRITE(32,367) NGAS1
367      FORMAT(' Number of gaseous species for transport=',I3)
c------------------------------------------------------------------------
c---------------------- For using EOS2 and ECO2 Flow module
c         IF (IEOS.EQ.2 .OR. IEOS.EQ.13)   THEN
c            NGAS1=0     ! Pco2 is from the flow, Pco2 is fixed in one time step
c         END IF
c------------------------------------------------------------------------
c
      if(ichdump.eq.0) then
        write(32,"(' Not print chemical speciation results'
     +   ' at each grid block and each time step')")
      else if(ichdump.eq.1) then
        write(32,"(' Print chemical speciation results'
     +   ' at each grid block and each time step')")
      else
        write(32,"(' Print chemical speciation results'
     +   ' at specified grid blocks and specified times')")
      end if
c
c------------------------------------------------------------------------
c
      if(kcpl.lt.0 .or. kcpl.gt.3)   kcpl=0
      if(kcpl.eq.0) then
         write(32,"(' Not calculating porosity and '
     +   'permeability changes due to mineral dis/pre')")
      else if(kcpl.eq.1) then
         write(32,"(' Calculate porosity and permeability changes '
     +   'and affect fluid flow')")
      else if(kcpl.eq.2) then
         write(32,"(' Monitor porosity and permeability changes '
     +   'but do not affect fluid flow')")
      else if(kcpl.eq.3) then
         write(32,"(' Calculate porosity and permeability changes '
     +   'and update liquid saturation due to porosity changes')")
       end if
C
C---------------------- Addition for EOS2 or ECO2 modules (TX,11/29/2001)
      if (Ico2h2o.ge.3 .or.Ico2h2o.lt.0)  Ico2h2o=0
      if(Ico2h2o.eq.0) then
        write(32,"(2x,'CO2 and H2O reaction sources are not considered'
     +   ' in the flow simulation (only for EOS2 or ECO2 module)')")
      else if(Ico2h2o.eq.1) then
        write(32,"(2x,'CO2 reaction sources are considered'
     +   ' in the flow simulation (only for EOS2 or ECO2 module)')")
      else if(Ico2h2o.eq.2) then
        write(32,"(2x,'CO2 and H2O reaction sources are considered'
     +   ' in the flow simulation (only for EOS2 or ECO2 module)')")
      end if
c
c------------------------------------------------------------------------
      if(numdr.eq.0) then
        write(32,"(' Calculate derivatives of mineral kinetic'
     +   ' rates by analytical method')")
      else
        write(32,"(' Calculate derivatives of mineral kinetic'
     +   ' rates by numerical method')")
      end if
c
c---------------------------------------------------------------------
c
        ISWITCH=0
        NJACOb=1
c--------------------------------------------
c---------  read-in constraints for chemical solver
c   sl1min - minimum liquid saturation
c   replaced by rcour   d1min  - minimum inter-block distance
c 
c  rcour     Flag and Courant number multiplier
c                if = 0. no Courant number limitation on time steps
c                if > 0.  = Time step limitation by RCOUR X Courant number for liq. and gas advection/diffusion
c                if < 0.  = Time step limitation by RCOUR X Courant number for liq. advection/diffusion
c
c   stimax - maximum stoichiometric ionic strength
c-------- add cnfact
c   cnfact - factor for reaction source terms
c            =1 fully implicit, =0 fully explicit
c            (for kin. minerals only)
c----------- cnfact not used anymore - fully implicit
c
      write(32,"(/' Constraints for chemistry solver:',
     &      /' *********************************')")
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      read(inprec,*,end=9999,err=9999) sl1min, rcour, stimax, cnfact
c
c...... Added for use in cs_cp
        dlstmx = dlog10(stimax)
c------- forces fully explicit source terms if non-iterative approach is
c selected - note, this is not really needed: source terms=0 if ispia=2
      if(ispia.eq.2) cnfact = 0.d0
      if(cnfact.ne.0.d0) cnfact=1.d0  ! no weighting allowed for now
      write(32,"(' Minimum liquid saturation = ',e10.4)") sl1min
      IF(rcour .EQ. 0.d0) THEN
         WRITE(32,361)
361   FORMAT(' No Courant number limitation on transport time steps')
      else
          write(32,363) dabs(rcour)
363   format(' Transport time steps for aq. species and gas limited by'
     +,f5.2,' x Courant number')
      END IF
C
      write(32,"(' Maximum stoichiometric ionic strength = ',e10.4)")
     & stimax
c
C------------------------------------------Read input and output file names
c
c... Thermodynamic database
c
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      thermo_in=inprec(1:20)                             
C
c... Iteration information file
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      OUTiter=inprec(1:20)
      WRITE(32,"(' Iteration info:',A20,A60)") OUTiter !,dftitle
      OPEN(UNIT=39,FILE=OUTiter,STATUS='UNKNOWN')

c... Aqueous conc. tecplot file
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      OUTplot=inprec(1:20)
C
c... Solid tecplot file
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      OUTsolid=inprec(1:20)
C
c... Gaseous species tecplot
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      OUTgas=inprec(1:20)
C
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      OUTtime=inprec(1:20)
C
c    changed so that mopr(8)=1 prints SI, =2 rctn rate, =3 sf area
c... Mineral saturation index (optional)
      IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14
     &     .OR. IEOS .EQ. 15 .OR. IEOS .EQ. 16)   THEN    ! for ECO2 module)
         IF (MOPR(8).ge.1)  THEN
            open(unit=77,file='min_SI.out',status='unknown')
            open(unit=80,file='CO2trap_tim.out',status='unknown')  ! For printing CO2 trapping over time 
         END IF
      END IF
c
c... Mineral reaction rate (optional)
      IF (MOPR(8).ge.2)  THEN
         open(unit=78,file='rctn_rate.out',status='unknown')
      END IF
c
c... Mineral reactive surface areas (optional)
      IF (MOPR(8).ge.3)  THEN
         open(unit=79,file='rct_sfarea.out',status='unknown')
      END IF
c
C----------------------------- Read weighting parameters
      WRITE(32,345)
345   FORMAT(/' Time/space weighting factors and diffusion coefficients'
     &':'/' *******************************************************')
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      READ(inprec,*,err=9999,end=9999) WTIME,WUPC,DIFUN,DIFUNG
      IF(WTIME.LT.0.d0 .OR. WTIME.GT.1.d0)  WTIME=1.d0      ! default
      IF(WUPC.LT.0.d0 .OR. WUPC.GT.1.d0)    WUPC=1.d0
      WRITE(32,346) WTIME,WUPC,DIFUN,DIFUNG
346   FORMAT('    WTIME     WUPC    DIFUN     DIFUNG'/,2F10.2,2E10.3)
c
C----------------------------Read data related to convergence criteria
c      READ(31,*)
      WRITE(32,350)
350   FORMAT(/' Convergence tolerance data:'/' ***********************')
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      READ(inprec,*,err=9999,end=9999) MAXITPTR,TOLTR,MAXITPCH,TOLCH,
     & MAXITPAD,TOLAD,TOLDC,TOLDR
C-TOLDC,TOLDR: Concentration and dis/pre changes for chemical steady-state
      if(ispia .eq. 2) maxitptr=1
      WRITE(32,371) MAXITPTR,TOLTR,MAXITPCH,TOLCH,TOLDC
371   FORMAT(' Maximum iterations for solving transport =',I5/
     &' Relative tolerance for solving transport =',E10.3/' Maximum'
     &' iterations for whole chemical system =',I5/' Relative tolerance'
     &' for whole chemical system =',E10.3/' (absolute mass action'
     &' tolerance set at min. 1.E-6)'/' Rel. con. change tolerance for'
     &' steady-state =',E10.3)
      write(32,372) MAXITPAD,TOLAD
372   format(' Maximum iterations for solving adsorption =',I5/
     &' Relative tolerance for solving adsorption =',E10.3/)
     
C--------------------------------Read data related to writing control
      WRITE(32,380)
380   FORMAT(' Printout/output control flags:'/' *********************')
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo

c     note we do not use (read in) IWCOMT anymore
      IWCOMT=0
      READ(inprec,*,err=9999,end=9999)
     &       NWTI,NWNOD,NWCOM,NWMIN,nwaq,nwads,nwexc,iconflag,minflag

      WRITE(32,390) NWTI
390   FORMAT(' Writing frequency in time =',I5)
      if(nwnod.ne.0) then
        write(32,391) nwnod
391     format(' Number of nodes to output in time file =',I5)
      end if
      if(nwcom.ne.0) then
        write(32,392) nwcom
392     format(' Number of primary species (total amounts) to output'
     &' time/plot files =',I5)
      end if
      if(nwmin.ne.0) then
        write(32,393) nwmin
393     format(' Number of minerals to output time/plot files =',I5)
        WRITE(32,"(' Plot minerals:',A20,A60)") OUTsolid  !,dftitle
        OPEN(UNIT=63,FILE=OUTsolid,STATUS='UNKNOWN')
      end if
      if(nwaq.ne.0) then
        write(32,394) nwaq
394     format(' Number of aqueous species to output time/plot files ='
     &,I5)
      end if
      if(nwads.ne.0) then
        write(32,395) nwads
395     format(' Number of adsorbed species to output time/plot files ='
     &,I5)
      end if
      if(nwexc.ne.0) then
        write(32,396) nwexc
396     format(' Number of exchange species to output time/plot files ='
     &,I5)
      end if
      write(32,397) iconflag,minflag
397   format(' ICONFLAG(0=mol/kg H2O,1=mol/L,2=g/L,3=mg/L(ppm))=',i5/
     &' MINFLAG(0=change in mol/m3, 1=change in volume fraction, 2='
     &'current volume fraction, 3=change in volume%)=',i5)
      WRITE(32,7)
7     FORMAT(/' Input and output file names:',
     1       /' ***************************')
      WRITE(32,"(' Thermo. database:',a20,a60)") thermo_in  !, dftitle
      WRITE(32,"(' Plot aqueous:',A20,A60)") OUTplot  !,dftitle
      WRITE(32,"(' Plot time series:',A20,A60)") OUTtime  !,dftitle
      OPEN(UNIT=61,FILE=OUTplot,STATUS='UNKNOWN')
      OPEN(UNIT=62,FILE=OUTtime,STATUS='UNKNOWN')
      if(ngas.ne.0) then
        OPEN(UNIT=66,FILE=OUTgas,STATUS='UNKNOWN')
        WRITE(32,"(' Plot gas:',A20,A60)") OUTgas  !,dftitle
      end if
c
c---------------------------------------------------------------
c
c.....Print unit name for concentrations
c
      if (iconflag .eq. 0)  Name_Conc =  'mol/kg h2o'
      if (iconflag .eq. 1)  Name_Conc =  'mol/L'
      if (iconflag .eq. 2)  Name_Conc =  'g/L'
      if (iconflag .eq. 3)  Name_Conc =  'mg/L(ppm)'
c
c.....Print unit name for minerals
c
      if (minflag .eq. 0)  then
         Name_Mine = 'Changes of abundance in  mol/m**3 medium'
      end if
c
      if (minflag .eq. 1)  then
       Name_Mine = 'Changes of abundance in volume fraction'
      end if
c
      if (minflag .eq. 2)  Name_Mine = 'Abundance in volume fraction'
c
      if (minflag .eq. 3)  then
       Name_Mine = 'Changes of abundance in volume fraction (%)'
      end if
c
c---------------------------------------------------------------
c
      
      if(nwnod.ne.0) then
       write(32,"(/' Elements for which to output time series:')")
       do
        read(31,"(a200)",err=9999) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo

       elemw(1)='     '
       if(nwnod.gt.0) then
        READ(inprec,'(40A5)',err=9999,end=9999) (ELEMW(I),I=1,NWNOD)
        WRITE(32,"(15(1x,A5))") (ELEMW(I),I=1,NWNOD)
       else
        i=0
        nwnod=0
        do
          if(i.ne.0) read(31,"(a200)",err=9999) inprec
          if(inprec.eq.'') exit
          i=i+1
          if(i.gt.200) go to 900
          read(inprec,"(a5)",err=9999) elemw(i)
          nwnod=i
        enddo
        write (32,"(15(1x,a5))") (elemw(i), i=1,nwnod) 
       endif
      endif
c 
      DO 402 I=1,NWNOD
         IWNOD(I)=0
         DO 405 J=1,NNOD
            IF (ELEM(J).EQ.ELEMW(I)) IWNOD(I)=J
405      CONTINUE
         IF (IWNOD(I).EQ.0)  THEN
            WRITE(32,406) ELEMW(I)
406         FORMAT(/1X,'Error: Grid block ',A5,' is not found in MESH')
         STOP
      END IF
402   CONTINUE

c     Aqueous species (total)
C---------------------------------------------------------------
      if(nwcom.ne.0) then
       write(32,"(' Aqueous components for which to output data: ')")
       do 
        read(31,"(a200)",err=9999) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       prt_com(1)='                    '
       if(nwcom.gt.0) then
        READ(inprec,*,err=9999,end=9999) (IWCOM(I),I=1,NWCOM)
        WRITE(32,"(15I5)") (IWCOM(I),I=1,NWCOM)
       else if(nwcom.lt.0) then
        i=0
        nwcom=0 
        do
          if(i.ne.0) read(31,"(a200)",err=9999) inprec
          if(inprec.eq.'') exit
          i=i+1
          if(i.gt.mpri) go to 905 
          prt_com(i)= inprec(1:20)
          nwcom=i
        enddo
        write(32,"(1x,a20)") (prt_com(i),i=1,nwcom)
       endif
      endif

c     Minerals
C---------------------------------------------------------------
      if(nwmin.ne.0) then
       write(32,"(/' Minerals for which to output data: ')")
       do 
        read(31,"(a200)",err=9999) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       prt_min(1)='                    '
       if(nwmin.gt.0) then
         IF(NWMIN.NE.0)  THEN
          READ(inprec,*,err=9999,end=9999) (IWMIN(I),I=1,NWMIN)
          WRITE(32,"(20I5)") (IWMIN(I),I=1,NWMIN)
       END IF

      else
        i=0
        nwmin=0 
        do
          if(i.ne.0) read(31,"(a200)",err=9999) inprec
          if(inprec.eq.'') exit
          i=i+1
          if(i.gt.mmin) go to 910 
          prt_min(i)= inprec(1:20)
          nwmin=i
        enddo
        write(32,"(1x,a20)") (prt_min(i),i=1,nwmin)
       endif
      endif

c      Individual aq. species
c----------------------------
      if(nwaq.ne.0) then
       write(32,"(' Individual aqueous species to output molalities:')")
       do 
         read(31,"(a200)",err=9999) inprec
         if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       prt_aq(1)='                    '
       if(nwaq.gt.0) then
         read(inprec,*,err=9999,end=9999) (iwaq(i),i=1,nwaq)
         write(32,"(20I5)") (iwaq(i),i=1,nwaq)
       else if(nwaq.lt.0) then
        i=0
        nwaq=0
        do
          if(i.ne.0) read(31,"(a200)",err=9999) inprec
          if(inprec.eq.'') exit
          i=i+1
          if(i.gt.maqx) go to 915 
          prt_aq(i)= inprec(1:20)
          nwaq=i
        enddo
         write(32,"(1x,a20)") (prt_aq(i),i=1,nwaq)
       endif
      endif

c      Individual sorbed species
c-------------------------------
      if(nwads.ne.0) then
       write(32,"(/' Surface complexes to output molalities:')")
       do 
         read(31,"(a200)",err=9999) inprec
         if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       prt_ads(1)='                    '
       if(nwads.gt.0) then
         read(inprec,*,err=9999,end=9999) (iwads(i),i=1,nwads)
         write(32,"(20I5)") (iwads(i),i=1,nwads)
       elseif(nwads.lt.0) then
        i=0
        nwads=0
        do
          if(i.ne.0) read(31,"(a200)",err=9999) inprec
          if(inprec.eq.'') exit
          i=i+1
          if(i.gt.mads) go to 920 
          prt_ads(i)= inprec(1:20)
          nwads=i
        enddo
         write(32,"(1x,a20)") (prt_ads(i),i=1,nwads)
       endif
      endif

c      Individual exchanged species
c----------------------------------
      if(nwexc.ne.0) then
       write(32,"(/' Exchanged species for which to output data:')")
       do 
         read(31,"(a200)",err=9999) inprec
         if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       prt_exc(1)='                    '
       if(nwexc.gt.0) then
         read(inprec,*,err=9999,end=9999) (iwexc(i),i=1,nwexc)
         write(32,"(20I5)") (iwexc(i),i=1,nwexc)
       elseif(nwexc.lt.0) then
        i=0
        nwexc=0
        do
          if(i.ne.0) read(31,"(a200)",err=9999) inprec
          if(inprec.eq.'') exit
          i=i+1
          if(i.gt.mexc) go to 925 
          prt_exc(i)= inprec(1:20)
          nwexc=i
        enddo
        write(32,"(1x,a20)") (prt_exc(i),i=1,nwexc)
       endif
      endif

C     Default types of chemical zones
c--------------------------------------
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      WRITE(32,295)
295   FORMAT(/' Default values for chemical zones:'/)
c
      if(mopr(11).eq.1) then
      WRITE(32,296)
296   FORMAT
     & (2X,'           IZIW IZBW IZMI IZGS IZAD IZEX IZPP IZKD VSED',
     & /2X,'           ---- ---- ---- ---- ---- ---- ---- ---- ----')
      else     
      WRITE(32,4296)
4296   FORMAT
     & (2X,'           IZIW IZBW IZMI IZGS IZAD IZEX IZPP IZKD',
     & /2X,'           ---- ---- ---- ---- ---- ---- ---- ----')
      endif
c
      vsedDF=0.d0   ! default sedimentation velocity
      vsedn=0.d0
      if(mopr(11).eq.1) then
c        read(31,"(a200)",err=9999) inprec
        READ(inprec,*,err=9999,end=9999) IZIWDF,IZBWDF,IZMIDF,IZGSDF,
     &   IZADDF,IZEXDF,izppdf,IZKDDF, VSEDDF
      else     
c        read(31,"(a200)",err=9999) inprec
        READ(inprec,*,err=9999,end=9999) IZIWDF,IZBWDF,IZMIDF,IZGSDF,
     &    IZADDF,IZEXDF,izppdf,IZKDDF
      endif
c
      if(mopr(11).eq.1) then
      WRITE(32,298) IZIWDF,IZBWDF,IZMIDF,IZGSDF,IZADDF,IZEXDF,
     1              izppdf,IZKDDF, vsedDF
298   FORMAT(12X,8I5,E10.3)
      else
      WRITE(32,4298) IZIWDF,IZBWDF,IZMIDF,IZGSDF,IZADDF,IZEXDF,
     1              izppdf,IZKDDF
4298   FORMAT(12X,8I5)
      endif
c
c     Types of chemical zones for individual grid blocks
C-------------------------------------------------------

      WRITE(32,300)
300   FORMAT(/' Data for chemical zones:'/)
      if(mopr(11).eq.1)then
        WRITE(32,310)
310     FORMAT
     &  (2X,' ELEM INDX IZIW IZBW IZMI IZGS IZAD IZEX IZPP IZKD  VSED',
     &  /2x,' ---- ---- ---- ---- ---- ---- ---- ---- ---- ----  ----')
      else
        write(32,4310)
4310      FORMAT
     &  (2X,' ELEM INDX IZIW IZBW IZMI IZGS IZAD IZEX IZPP IZKD',
     &  /2x,' ---- ---- ---- ---- ---- ---- ---- ---- ---- ----')
      endif
c       
      DO 490 J=1,NNOD       ! removed the C'out 
         IZONEIW(J)=0
         IZONEBW(J)=0
         IZONEM(J)=0
         IZONEG(J)=0
         IZONED(J)=0
         IZONEX(J)=0
         izonpp(j)=0       ! permeability-porosity relations
         IZONEKD(J)=0      ! line adsorption Kd
         if(mopr(11).eq.1)vsed(j)=vsedDF  !added 7/09
c
490   CONTINUE
c
c
499   CONTINUE
c
      do 
       read(31,"(a200)",err=9999) inprec
       if(inprec(1:3).eq.'end'.or.
     &       inprec(1:3).eq.'END') go to 599
       if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      el=inprec(1:5)
      IF(EL .EQ. '     ') GO TO 599
c     to read sedimentation velocity vsedn (in m2/s)
      if(mopr(11).eq.1) then
       READ(inprec(6:200),*,err=9999,end=9999) NSEQ,NADD,IZIW,
     &       IZBW,IZMI,IZGS,IZAD,IZEX,izpp,IZKD,vsedn
      else
       READ(inprec(6:200),*,err=9999,end=9999) NSEQ,NADD,IZIW,
     &       IZBW,IZMI,IZGS,IZAD,IZEX,izpp,IZKD
      endif
c
c     check if node exists and print warning if not
      iflag=0
      do j=1,nnod
         if (elem(j) .eq. el) iflag=1
      enddo 
      if(iflag.eq.0) then
       write(32,"(/' Warning: grid block ',A5,' is not found'/)") el
       write(33,"(/' Warning: grid block ',A5,' is not found'/)") el 
      endif
c
      DO 500 J=1,NNOD
         IF(ELEM(J).NE.EL) GOTO 500
         IZONEIW(J)=IZIW
         IZONEBW(J)=IZBW
         IZONEM(J)=IZMI
         IZONEG(J)=IZGS
         IZONED(J)=IZAD
         IZONEX(J)=IZEX
         izonpp(j)=izpp
         IZONEKD(J)=IZKD
         if(mopr(11).eq.1)vsed(j)=vsedn    ! sedimentation velocity m/s
         GOTO 509
500   CONTINUE
c
509   CONTINUE
      DO 1905 I=1,NSEQ
         N1=J+I*NADD
         IZONEIW(N1)=IZONEIW(J)
         IZONEBW(N1)=IZONEBW(J)
         IZONEM(N1) =IZONEM(J)
         IZONEG(N1) =IZONEG(J)
         IZONED(N1) =IZONED(J)
         IZONEX(N1) =IZONEX(J)
         izonpp(n1) =izonpp(j)
         IZONEKD(N1)=IZONEKD(J)
         if(mopr(11).eq.1)vsed(n1)=vsed(j) ! sedimentation velocity m/s
c
 1905 CONTINUE
c
      GOTO 499
C
C---------------------------------------------Get default values
599   CONTINUE
c
      DO 520 J=1,NNOD
         IF(IZONEIW(J) .EQ. 0) IZONEIW(J)=IZIWDF
         IF(IZONEBW(J) .EQ. 0) IZONEBW(J)=IZBWDF
         IF(IZONEM(J)  .EQ. 0) IZONEM(J)=IZMIDF
         IF(IZONEG(J)  .EQ. 0) IZONEG(J)=IZGSDF
         IF(IZONED(J)  .EQ. 0) IZONED(J)=IZADDF
         IF(IZONEX(J)  .EQ. 0) IZONEX(J)=IZEXDF
         IF(izonpp(j)  .EQ. 0) izonpp(j)=izppdf
         IF(IZONEKD(J) .EQ. 0) IZONEKD(J)=IZKDDF
c         vsed(j)=vsedDF !ns5/09  sedimentation velocity m/s
c
520   CONTINUE
C
      DO 540 J=1,NNOD
c       initialize all to 1 to avoid array bouncing when 0     
        if(IZONEIW(J).eq.0) IZONEIW(J)=1
        if(IZONEBW(J).eq.0) IZONEBW(J)=1
        if(IZONEM(J).eq.0) IZONEM(J)=1
        if(IZONEG(J).eq.0) IZONEG(J)=1
        if(IZONED(J).eq.0) IZONED(J)=1
        if(IZONEX(J).eq.0) IZONEX(J)=1
        if(izonpp(j).eq.0) izonpp(j)=1
        if(IZONEKD(J).eq.0) IZONEKD(J)=1
      if(mopr(11).eq.1) then
         WRITE(32,560) ELEM(J),J,IZONEIW(J),IZONEBW(J),IZONEM(J),
     &          IZONEG(J),IZONED(J),IZONEX(J),izonpp(j),IZONEKD(J),
     &          vsed(j)
560      FORMAT(2X,A5,9I5,E10.3)
      else
         WRITE(32,4560) ELEM(J),J,IZONEIW(J),IZONEBW(J),IZONEM(J),
     &          IZONEG(J),IZONED(J),IZONEX(J),izonpp(j),IZONEKD(J)
4560     FORMAT(2X,A5,9I5)
      endif
c
540   CONTINUE
c
      WRITE(32,*) 'end'
C
      write(33,*) ' ---> finished reading solute.inp '
      write(32,*) ' ---> finished reading solute.inp '
c
      close (31)
      close (32)
      close (33)
c---------  unit 32 will be reopen in init and couple run_log file
c
      RETURN

900   write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*) 'grid blocks to print exceed maximum number of 200'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

905   write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*) 'number of primary species to print exceeds maximum'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

910   write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*) 'number of minerals to print exceeds maximum'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

915   write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*)'number of individual species to print exceeds maximum'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

920   write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*)'number of sorption species to print exceeds maximum'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

925   write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*)'number of exchange species to print exceeds maximum'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

9999  write(33,*) 'error reading solute.inp - see solute.out'
      write(32,*) 'error on input'
      write(32,"('  last input record starts with',a40)") inprec(1:40)
      stop

      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE get_prt_dat
c
C***** To get indices of species/minerals to print if these were entered
c      as names in solute.inp
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      INCLUDE 'common_v2.inc'
      COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
      common/wricon2/ prt_com(mpri),prt_min(mmin),prt_aq(maqx),
     &                 prt_ads(mads),prt_exc(mexc)
      character*20 prt_com,prt_min,prt_aq,prt_ads,prt_exc,blank
      data blank/'                    '/
c
      if(prt_com(1).ne.blank)
     &  call getind(npri,nwcom,napri,prt_com,iwcom) 
      if(prt_min(1).ne.blank)
     &  call getind(nmin,nwmin,namin,prt_min,iwmin)
      if(prt_aq(1).ne.blank)
     &  call getind(npri+naqx,nwaq,naaqt,prt_aq,iwaq)
      if(prt_ads(1).ne.blank)
     &  call getind(nads,nwads,naads,prt_ads,iwads)
      if(prt_exc(1).ne.blank)
     &  call getind(nexc,nwexc,naexc,prt_exc,iwexc)

      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine getind (nx,nprt,xname,prtnam,indprt)
c
c     returns indices of species and minerals to print if names instead of
c      indices were entered in solute.inp   NS10/08
c
c     nx:     total number of species/minerals in simulation 
c     nprt:   number of species/minerals to print
c     xname:  name of species/minerals in simulation
c     prtnam: name of species/minerals to print
c     indprt: index of species/minerals to print (returned) 
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      character*20 xname(nx),prtnam(nprt),label
      integer*8 nx, nprt, indprt(nprt)
c
      n=1
      do i = 1,nprt
       imatch=0
       label = prtnam(i)
       call name_conv(label)   !change to lower case etc., as in init
       prtnam(i)=label
       do j = 1,nx
        if(xname(j).eq.prtnam(i)) then
         indprt(n)=j
         n=n+1
         imatch=1
         exit
        endif
       enddo 
        if(imatch.eq.0) then
         write(32,"(/5x,'species to print, ',a20,
     &   /7x,'is not part of the chemical system - stop')") prtnam(i)
         stop
        endif
      enddo
c
      return
      end
