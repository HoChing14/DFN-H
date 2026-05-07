c
      SUBROUTINE MATRIXC_implicit
C
C---------------- SET UP THE COEFFICIENT MATRIX FOR TRANSPORT OF SOLUTES
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
      include 'perm_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
C
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
c
c.....Add tortuosity exponent (ptort) and critical porosity (phicrit)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
c.....Multiphase tortuosity at each grid block
      common/tortmp/tortliq(mnel),tortgas(mnel)
c
      COMMON/SOCH/MAT(MAXMAT)
c
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)
      COMMON/C4/DEL2(MNCON)
      COMMON/C5/AREA(MNCON)
      COMMON/C9/ELEM1(MNCON)
      COMMON/C10/ELEM2(MNCON)
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
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
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
c
      COMMON/PORVEL/VEL(MNPH*MNCON)
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)    ! darcy velocity
      COMMON/SOLUTE6/SLOLD(MNEL)           ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)           ! old gas saturation
      COMMON/SOLUTE8/SL1(MNEL)             ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)             ! new gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)         ! old (initial) porosity
      COMMON/SOLUTE11/NPRI,npaq,npads      ! number of chemical component
      COMMON/PARNP/NPL,NPG                 ! specify in EOS module
C
      COMMON/AMMISC/IABC,ISOLVC
      COMMON/PRINTC/NOW                    ! print control
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT 
C
C.....Interface area reduction factor
      common/afactor/a_fm2(mncon)  ! advection area reduction (flow from F to M)
      common/afactord/a_fmd(mncon) ! diffusion area reduction (Both sides)
c.... Save modified active fracture area for reaction --
c.....Limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
c
c.....Indicators from EOS module
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
      common/vmineral/pre(mnod,mmin),pre0(mnod,mmin),
     +  pinit(mnod,mmin+1)
!
      common/chemgrid/c(mnod,maqt),utold(mnod,mpri),ut(mnod,mpri),
     & rhand(mnod,mpri),rsource(mnod,mpri),ph(mnod),gP(mnod,mgas),
     & aream(mnod,mmin),sads(mnod),psi(mnod),
     & d(mnod,mads),supadn(mnod,msurf),phip(mnod,msurf),
     & surfads(msurf),ub(mbound,mpri),ctot(mnod,mpri),cnfact
!
      common/names/napri(mpri),naaqx(maqx),naaqt(maqt),namin(mmin),
     +   nagas(mgas),naexc(mexc),naads(mads)
      character*20 napri,naaqx,naaqt,namin,nagas,naexc,naads
      common/transport/izoneiw(mnod),izonebw(mnod),
     +   izonem(mnod),izoneg(mnod),izoned(mnod),izonex(mnod),
     +   izonpp(mnod)
      common/constraints/sl1min,stimax,dlstmx
      common/gasprop/dmwgas(mgas),diamol(mgas)
C
c...................porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G4/ELEG(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G12/LCOM(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
C
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
      COMMON/GASLAW/RGAS,AMS,AMA,CVGAS
      COMMON/DFM/TIMAX,REDLT
      COMMON/BC/NELA
      CHARACTER*5 ELEM,ELEM1,ELEM2,MAT,ELEG
C
C-------------------------------------- For coupling with reactive transport
C
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
c
C--------------------- For time stepping use (subroutine MAX_DELT)
      COMMON/DIFUNT_L1/DIFUNT_L(MNCON)

      COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3) 
C------------------------------------------------------------------------------
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *****MATRIXC q2.1, 1999.4.12: ASSEMBLE MATRIX FOR SOLUTE'
     X' TRANSPORT**********')
C
C$$$$$FOR IABC=0 COLUMN INDICES WILL BE STORED IN ICN, OTHERWISE IN JVECT
C-----INITIALIZE COUNTER FOR MATRIX ELEMENTS.
      NZ=0
C
C******************* LOOP OVER ELEMENTS ******************************
C
C*****COMPUTE ALL QUANTITIES WHICH DEPEND ONLY UPON VARIABLES PERTAINING
C     TO ONE VOLUME ELEMENT.
C
      DO 100 N=1,NEL
!
C------------------------------------COMPUTE RIGHT-HAND SIDE TERMS
C-------------------CONTRIBUTED FROM INITIAL CONDITIONS
C
      phislo = PHIOLD(N)*SLOLD(N)
      DO 110 IPRI=1,NPRI
         RHAND(N,IPRI) = phislo*UTOLD(N,IPRI)
  110 CONTINUE
C
C-----------------------------COMPUTE COEFFICIENT MATRIX ELEMENT
      IRN(NZ+1)=N                     ! row index
      IF(IABC.EQ.0) ICN(NZ+1)=N       ! column index
      IF(IABC.NE.0) JVECT(NZ+1)=N
      CO(NZ+1)=phisl1(n)
      NZ=NZ+1
C
C+++++++++END OF ASSIGNMENT OF ONE-ELEMENT TERMS++++++++++++++++++++
C
  100 CONTINUE
C
C
C************** TAKE INTO ACCOUNT EXTERNAL SOURCE *********************
C
      DO 140 IOGN=1,NOGN
         J=NEXG(IOGN)
c
         densw = dwat(j)*1000.d0  !??? densw not needed here?
c
         EVOLJ=EVOL(J)*1000.D0
         IF (EVOL(J).EQ.0.0D0)  EVOLJ=1000.0D0
         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).GE.0.D0
     &      .and.izonebw(j).ne.0) THEN
               IZONEJ=IZONEBW(J)
               genchm = G(IOGN)*DELTEX/(EVOLJ)
            DO 160 IPRI=1,NPRI
               GCHEM=UB(IZONEJ,IPRI)*genchm
               RHAND(J,IPRI)=RHAND(J,IPRI)+GCHEM
  160       CONTINUE
         END IF
  140 CONTINUE
C
      DO 145 IOGN=1,NOGN
         J=NEXG(IOGN)
c
         vliqw = dwat(j)             ! dwat is in kg/l
c
         EVOLJ=EVOL(J)*1000.D0
         IF (EVOL(J).EQ.0.0D0)  EVOLJ=1000.0D0
         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).LT.0.D0) THEN
            CO(J)=CO(J)-G(IOGN)*DELTEX/EVOLJ !Modify matrix for pump
         END IF
!
!-----------
!.....For well on deliverability
!    (production occurs against specified bottomhole pressure)
!-----------
!
         gdelv=g(iogn)*FF((IOGN-1)*NPH+NPL)
         IF(LCOM(IOGN).EQ.(NK1+1) .AND. gdelv .LT.0.D0)

     +                                                  THEN
            FFI=FF((IOGN-1)*NPH+NPL)
            CO(J)=CO(J)-gdelv*DELTEX/EVOLJ
         END IF
!
!-----------------------------------------------------------------------
!
  145 CONTINUE
C
      IF(MOPR(2).GE.1)   THEN
      IF(NOW.EQ.1)  THEN
C
         write (32,*)
         write (32,*)
         write (32,33)
  33     format ('      PHIOLD         PHI       SLOLD          SL',
     +   '         VOL     izonebw')
         do n=1,nel
            write (32,'(5E12.4,I12)') phiold(n),phi(n),slold(n),SL1(n),
     +                            EVOL(N),izonebw(N)
         end do
C
         write (32,*)
         write (32,*) '       GRID        LCOM      G(source) '
         do IOGN=1,NOGN
            J=NEXG(IOGN)
            write (32,'(2I10, E12.4)') J,LCOM(IOGN),G(IOGN)
         end do
c
         write (32,*)
         write(32,*) '------utold(i,n)----'
         do ii=1,nel
            write(32,'(5e12.4)') (utold(ii,nnn),nnn=1,npri)
         end do
C
         write (32,*)
         write (32,*) '   ---rhand---'
         do n=1,nel
            write (32,'(7E12.4)') (rhand(n,ip),ip=1,npri)
         end do
C
         write (32,*)
         write (32,*) '   ----ub----'
         write (32,'(5E12.4)') (ub(1,ip),ip=1,npri)
C
         write (32,*)
         write (32,*) ' nel  NPRi      DELTEX      TIMETOT     vliqw  '
         write (32,'(2I5,3E12.4)') nel,npri,DELTEX,TIMETOT,vliqw
C
         write (32,*)
         write (32,*) 'dar  vel_gas      vel_liquid'
         do n=1,ncon
            write (32,'(2E12.4)') (veldar((N-1)*NPH+NP),np=1,2)
         end do
C
         write (32,*)
         write (32,*) 'Por  vel_gas      vel_liquid'
         do n=1,ncon
            write (32,'(2E12.4)') (vel((N-1)*NPH+NP),np=1,2)
         end do
C
      end if
      END IF
C
C********* LOOP OVER CONNECTIONS ***********************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND UPON VARIABLES FOR TWO VOLUME
C     ELEMENTS ("INTERFACE QUANTITIES").
C
      DO 200 N=1,NCON
C
      N1=NEX1(N)         ! the first element of a connection
      N2=NEX2(N)         ! the second element of a connection
C
C-----OBTAIN SOME QUANTITIES PERTAINING TO CONNECTION.
C
c***************************************For active fracture model
      AX1 = AREA(N)
      a_fm=a_fm2(n)
      AX=a_fm*AX1      ! advection interface area
c------------------Diffusion reduction factor is for flow in both directions
c--------------------The factor is calculated according to fracture (F) side
      a_fmdd=a_fmd(n)
      AXD=a_fmdd*AX1      ! diffusion interface area
c**************************************************************************
!
      D1=DEL1(N)
      D2=DEL2(N)
      DISTAN=D1+D2
      if (sl1(n1).le.sl1min.or.sl1(n2).le.sl1min)   then
         sla = 0.d0
                                                    else
         pslt1 = phisl1(n1)*tortliq(n1)
         pslt2 = phisl1(n2)*tortliq(n2)
         sla = (distan*pslt1*pslt2)/(d1*pslt2+d2*pslt1)
      end if
!
!.....Harmonic mean of product of porosity and liquid saturation
!.....Distance weighting for tortuosity
!.....Mass balance terms (phi * SL) not included in Millington-Quirk above
!
      DIFUNT_L(N)=DIFUN*sla
      DIFLIQ=DIFUNT_L(N)
!
c--------------------------------------for active fracture model
      FLUXD=AXD*DIFLIQ/DISTAN     ! a term derived from diffusion
c---------------------------------------------------------------
c
C+++++ASSIGN ALL INTERFACE TERMS+++++++++++++++++++++++++++++++++
c
      NP=NPL
      NI=(N-1)*NPH+NP
      VELN=VELDAR(NI)                    ! liquid Darcy velocity
      IF (VELN.GE.0.D0)    THEN
         FAC12=1.D0-WUPC                 ! up-stream weighting fatcor
         FAC21=WUPC
                          ELSE
         FAC12=WUPC
         FAC21=1.D0-WUPC
      END IF
C
      EVOLN1=EVOL(N1)
      IF (EVOL(N1).EQ.0.0D0)  EVOLN1=1000.0D0
      DUM1=DELTEX/EVOLN1
      DELTV1  = DUM1
C
         EVOLN2=EVOL(N2)
      IF (EVOL(N2).EQ.0.0D0)  EVOLN2=1000.0D0
      DUM2=DELTEX/EVOLN2
      DELTV2  = DUM2
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM N2
c
      IRN(NZ+1)=N1
      IF(IABC.EQ.0) ICN(NZ+1)=N2
      IF(IABC.NE.0) JVECT(NZ+1)=N2
      DUM12=(FAC12-1.D0)*AX*VELN
      CO(NZ+1)=DELTV1*DUM12                 !  advection contribution
      CO(NZ+1)=CO(NZ+1)-FLUXD*DELTV1        !  diffusion contribution
C
      IF (SL1(N1).le.sl1min) CO(NZ+1)=0.D0
      IF (EVOL(N1).EQ. 0.0D0) CO(NZ+1)=0.D0
C
      NZ=NZ+1
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FRO N1
c
      IRN(NZ+1)=N2
      IF(IABC.EQ.0) ICN(NZ+1)=N1
      IF(IABC.NE.0) JVECT(NZ+1)=N1
      DUM21=(FAC21-1.D0)*AX*(-1.D0)*VELN
      CO(NZ+1)=DELTV2*DUM21                 !  advection contribution
      CO(NZ+1)=CO(NZ+1)-FLUXD*DELTV2        !  diffusion contribution
C
      IF(SL1(N2).le.sl1min)   CO(NZ+1)=0.D0
      IF (EVOL(N2).EQ. 0.0D0) CO(NZ+1)=0.D0
C
      NZ=NZ+1
C
C-------------------------DIAGONAL TERM IN EQUATION FOR ELEMENT N1
c
      DUM11=FAC12*AX*VELN
      CO(N1)=CO(N1)-DELTV1*DUM11          !  advection contribution
      CO(N1)=CO(N1)+FLUXD*DELTV1          !  diffusion contribution
C
C-------------------------DIAGONAL TERM IN EQUATION FOR ELEMENT N2
c
      DUM22=FAC21*AX*VELN
      CO(N2)=CO(N2)+DELTV2*DUM22          !  advection contribution
      CO(N2)=CO(N2)+FLUXD*DELTV2          !  diffusion contribution
C
C+++++END OF ASSIGNMENT OF INTERFACE TERMS+++++++++++++++++++++++
C
  200 CONTINUE
C-------------------------------------Modify matrix for inactive element
c
       DO 225 N=1,NEL
         IF (EVOL(N) .EQ. 0.0D0) CO(N)=1.0D0
         IF (SL1(N)  .le.sl1min) CO(N)=1.D0
225    CONTINUE
c
C-----------------------------------------------------------------------
C
      IF(MOPR(2).GE.1)   THEN
      IF(NOW.EQ.1)  THEN
         write (32,*)  '   iRN  ICN      CO------MATRIX-----'
         do nnn=1,NZ
            write (32,'(i6, i6,E12.4)') irn(nnn),icn(nnn),co(nnn)
         end do
      END IF
      END IF
c
C-----------------------------------------------------------------------
C
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE MATRIXC
C
C---------------- SET UP THE COEFFICIENT MATRIX FOR TRANSPORT OF SOLUTES
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
      include 'perm_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
C
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
c
c.....Add tortuosity exponent (ptort) and critical porosity (phicrit)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
c.....Multiphase tortuosity at each grid block
      common/tortmp/tortliq(mnel),tortgas(mnel)
c
      COMMON/SOCH/MAT(MAXMAT)
c
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)
      COMMON/C4/DEL2(MNCON)
      COMMON/C5/AREA(MNCON)
      COMMON/C9/ELEM1(MNCON)
      COMMON/C10/ELEM2(MNCON)
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
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
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
c
      COMMON/PORVEL/VEL(MNPH*MNCON)
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)    ! darcy velocity
      COMMON/SOLUTE6/SLOLD(MNEL)           ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)           ! old gas saturation
      COMMON/SOLUTE8/SL1(MNEL)             ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)             ! new gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)         ! old (initial) porosity
      COMMON/SOLUTE11/NPRI,npaq,npads      ! number of chemical component
      COMMON/PARNP/NPL,NPG                 ! specify in EOS module
C
      COMMON/AMMISC/IABC,ISOLVC
      COMMON/PRINTC/NOW                    ! print control
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT 
C
C.....Interface area reduction factor
      common/afactor/a_fm2(mncon)  ! advection area reduction (flow from F to M)
      common/afactord/a_fmd(mncon) ! diffusion area reduction (Both sides)
c.... Save modified active fracture area for reaction --
c.....Limit is sl1min, not residual saturation
      common/afactorr/a_fmr(mnel)
c
c.....Indicators from EOS module
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
      common/vmineral/pre(mnod,mmin),pre0(mnod,mmin),
     +  pinit(mnod,mmin+1)
!
      common/chemgrid/c(mnod,maqt),utold(mnod,mpri),ut(mnod,mpri),
     & rhand(mnod,mpri),rsource(mnod,mpri),ph(mnod),gP(mnod,mgas),
     & aream(mnod,mmin),sads(mnod),psi(mnod),
     & d(mnod,mads),supadn(mnod,msurf),phip(mnod,msurf),
     & surfads(msurf),ub(mbound,mpri),ctot(mnod,mpri),cnfact
!
      common/names/napri(mpri),naaqx(maqx),naaqt(maqt),namin(mmin),
     +   nagas(mgas),naexc(mexc),naads(mads)
      character*20 napri,naaqx,naaqt,namin,nagas,naexc,naads
      common/transport/izoneiw(mnod),izonebw(mnod),
     +   izonem(mnod),izoneg(mnod),izoned(mnod),izonex(mnod),
     +   izonpp(mnod)
      common/constraints/sl1min,stimax,dlstmx
      common/gasprop/dmwgas(mgas),diamol(mgas)
C
c...................porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G4/ELEG(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G12/LCOM(MNOGN)
      COMMON/G26/FF(MNPH*MNOGN)
C
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
      COMMON/GASLAW/RGAS,AMS,AMA,CVGAS
      COMMON/DFM/TIMAX,REDLT
      COMMON/BC/NELA
      CHARACTER*5 ELEM,ELEM1,ELEM2,MAT,ELEG
C
C-------------------------------------- For coupling with reactive transport
C
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
c
C--------------------- For time stepping use (subroutine MAX_DELT)
      COMMON/DIFUNT_L1/DIFUNT_L(MNCON)

      COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3) 
C------------------------------------------------------------------------------
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *****MATRIXC q2.1, 1999.4.12: ASSEMBLE MATRIX FOR SOLUTE'
     X' TRANSPORT*************')
C
C$$$$$FOR IABC=0 COLUMN INDICES WILL BE STORED IN ICN, OTHERWISE IN JVECT
C-----INITIALIZE COUNTER FOR MATRIX ELEMENTS.
      NZ=0
C
C******************* LOOP OVER ELEMENTS ******************************
C
C*****COMPUTE ALL QUANTITIES WHICH DEPEND ONLY UPON VARIABLES PERTAINING
C     TO ONE VOLUME ELEMENT.
C
      DO 100 N=1,NEL
!
C------------------------------------COMPUTE RIGHT-HAND SIDE TERMS
C-------------------CONTRIBUTED FROM INITIAL CONDITIONS
C
      phislo = PHIOLD(N)*SLOLD(N)
      DO 110 IPRI=1,NPRI
         RHAND(N,IPRI) = phislo*UTOLD(N,IPRI)
  110 CONTINUE
C
C-----------------------------COMPUTE COEFFICIENT MATRIX ELEMENT
      IRN(NZ+1)=N                     ! row index
      IF(IABC.EQ.0) ICN(NZ+1)=N       ! column index
      IF(IABC.NE.0) JVECT(NZ+1)=N
      CO(NZ+1)=phisl1(n)
      NZ=NZ+1
C
C+++++++++END OF ASSIGNMENT OF ONE-ELEMENT TERMS++++++++++++++++++++
C
  100 CONTINUE
C
C
C************** TAKE INTO ACCOUNT EXTERNAL SOURCE *********************
C
      DO 140 IOGN=1,NOGN
         J=NEXG(IOGN)
c
         densw = dwat(j)*1000.d0  !??? densw not needed here?
c
         EVOLJ=EVOL(J)*1000.D0
         IF (EVOL(J).EQ.0.0D0)  EVOLJ=1000.0D0
         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).GE.0.D0
     &      .and.izonebw(j).ne.0) THEN
               IZONEJ=IZONEBW(J)
               genchm = G(IOGN)*DELTEX/(EVOLJ)
            DO 160 IPRI=1,NPRI
               GCHEM=UB(IZONEJ,IPRI)*genchm
               RHAND(J,IPRI)=RHAND(J,IPRI)+GCHEM
  160       CONTINUE
         END IF
  140 CONTINUE
C
      DO 145 IOGN=1,NOGN
         J=NEXG(IOGN)
c
         vliqw = dwat(j)             ! dwat is in kg/l
c
         EVOLJ=EVOL(J)*1000.D0
         IF (EVOL(J).EQ.0.0D0)  EVOLJ=1000.0D0
         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).LT.0.D0) THEN
            CO(J)=CO(J)-WTIME*G(IOGN)*DELTEX/EVOLJ !Modify matrix for pump
               genchmu = (1.D0-WTIME)*G(IOGN)*DELTEX/(EVOLJ*vliqw)
            DO IPRI=1,NPRI
               GCHEM=UTOLD(J,IPRI)*genchmu
               RHAND(J,IPRI)=RHAND(J,IPRI)+GCHEM
            END DO
         END IF
!
!-----------
!.....For well on deliverability
!    (production occurs against specified bottomhole pressure)
!-----------
!
         gdelv=g(iogn)*FF((IOGN-1)*NPH+NPL)
         IF(LCOM(IOGN).EQ.(NK1+1) .AND. gdelv .LT.0.D0)

     +                                                  THEN
            FFI=FF((IOGN-1)*NPH+NPL)
            CO(J)=CO(J)-WTIME*gdelv*DELTEX/EVOLJ

               genchmu = (1.D0-WTIME)*gdelv*DELTEX/(EVOLJ*vliqw)
            DO IPRI=1,NPRI
               GCHEM=UTOLD(J,IPRI)*genchmu
               RHAND(J,IPRI)=RHAND(J,IPRI)+GCHEM
            END DO
         END IF
!
!-----------------------------------------------------------------------
!
  145 CONTINUE
C
      IF(MOPR(2).GE.1)   THEN
      IF(NOW.EQ.1)  THEN
C
         write (32,*)
         write (32,*)
         write (32,33)
  33     format ('      PHIOLD         PHI       SLOLD          SL',
     +   '         VOL     izonebw')
         do n=1,nel
            write (32,'(5E12.4,I12)') phiold(n),phi(n),slold(n),SL1(n),
     +                            EVOL(N),izonebw(N)
         end do
C
         write (32,*)
         write (32,*) '       GRID        LCOM      G(source) '
         do IOGN=1,NOGN
            J=NEXG(IOGN)
            write (32,'(2I10, E12.4)') J,LCOM(IOGN),G(IOGN)
         end do
c
         write (32,*)
         write(32,*) '------utold(i,n)----'
         do ii=1,nel
            write(32,'(5e12.4)') (utold(ii,nnn),nnn=1,npri)
         end do
C
         write (32,*)
         write (32,*) '   ---rhand---'
         do n=1,nel
            write (32,'(7E12.4)') (rhand(n,ip),ip=1,npri)
         end do
C
         write (32,*)
         write (32,*) '   ----ub----'
         write (32,'(5E12.4)') (ub(1,ip),ip=1,npri)
C
         write (32,*)
         write (32,*) ' nel  NPRi      DELTEX      TIMETOT     vliqw  '
         write (32,'(2I5,3E12.4)') nel,npri,DELTEX,TIMETOT,vliqw
C
         write (32,*)
         write (32,*) 'dar  vel_gas      vel_liquid'
         do n=1,ncon
            write (32,'(2E12.4)') (veldar((N-1)*NPH+NP),np=1,2)
         end do
C
         write (32,*)
         write (32,*) 'Por  vel_gas      vel_liquid'
         do n=1,ncon
            write (32,'(2E12.4)') (vel((N-1)*NPH+NP),np=1,2)
         end do
C
      end if
      END IF
C
C********* LOOP OVER CONNECTIONS ***********************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND UPON VARIABLES FOR TWO VOLUME
C     ELEMENTS ("INTERFACE QUANTITIES").
C
      DO 200 N=1,NCON
C
      N1=NEX1(N)         ! the first element of a connection
      N2=NEX2(N)         ! the second element of a connection
C
C-----OBTAIN SOME QUANTITIES PERTAINING TO CONNECTION.
C
c***************************************For active fracture model
      AX1 = AREA(N)
      a_fm=a_fm2(n)
      AX=a_fm*AX1      ! advection interface area
c------------------Diffusion reduction factor is for flow in both directions
c--------------------The factor is calculated according to fracture (F) side
      a_fmdd=a_fmd(n)
      AXD=a_fmdd*AX1      ! diffusion interface area
c**************************************************************************
!
      D1=DEL1(N)
      D2=DEL2(N)
      DISTAN=D1+D2
      if (sl1(n1).le.sl1min.or.sl1(n2).le.sl1min)   then
         sla = 0.d0
                                                    else
         pslt1 = phisl1(n1)*tortliq(n1)
         pslt2 = phisl1(n2)*tortliq(n2)
         sla = (distan*pslt1*pslt2)/(d1*pslt2+d2*pslt1)
      end if
!
!.....Harmonic mean of product of porosity and liquid saturation
!.....Distance weighting for tortuosity
!.....Mass balance terms (phi * SL) not included in Millington-Quirk above
!
       DIFUNT_L(N)=DIFUN*sla
       DIFLIQ=DIFUNT_L(N)
!
c--------------------------------------for active fracture model
       FLUXD=AXD*DIFLIQ/DISTAN     ! a term derived from diffusion
c---------------------------------------------------------------
c
C+++++ASSIGN ALL INTERFACE TERMS+++++++++++++++++++++++++++++++++
c
      NP=NPL
      NI=(N-1)*NPH+NP
      VELN=VELDAR(NI)                    ! liquid Darcy velocity
      IF (VELN.GE.0.D0)    THEN
         FAC12=1.D0-WUPC                 ! up-stream weighting fatcor
         FAC21=WUPC
                          ELSE
         FAC12=WUPC
         FAC21=1.D0-WUPC
      END IF
C
      EVOLN1=EVOL(N1)
      IF (EVOL(N1).EQ.0.0D0)  EVOLN1=1000.0D0
      DUM1=DELTEX/EVOLN1
      DELTV1=WTIME*DUM1
      DELTV1R=(1.D0-WTIME)*DUM1
C
         EVOLN2=EVOL(N2)
      IF (EVOL(N2).EQ.0.0D0)  EVOLN2=1000.0D0
      DUM2=DELTEX/EVOLN2
      DELTV2=WTIME*DUM2
      DELTV2R=(1.D0-WTIME)*DUM2
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM N2
c
      IRN(NZ+1)=N1
      IF(IABC.EQ.0) ICN(NZ+1)=N2
      IF(IABC.NE.0) JVECT(NZ+1)=N2
      DUM12=(FAC12-1.D0)*AX*VELN
      CO(NZ+1)=DELTV1*DUM12                 !  advection contribution
      CO(NZ+1)=CO(NZ+1)-FLUXD*DELTV1        !  diffusion contribution
        dv1d12 = DELTV1R*DUM12
        dv1fld = DELTV1R*FLUXD
      DO IPRI=1,NPRI
         RHAND(N1,IPRI)= RHAND(N1,IPRI)-dv1d12*UTOLD(N2,IPRI)
         RHAND(N1,IPRI)= RHAND(N1,IPRI)+dv1fld*UTOLD(N2,IPRI)
      END DO
C
      IF (SL1(N1).le.sl1min) CO(NZ+1)=0.D0
      IF (EVOL(N1).EQ. 0.0D0) CO(NZ+1)=0.D0
C
      NZ=NZ+1
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FRO N1
c
      IRN(NZ+1)=N2
      IF(IABC.EQ.0) ICN(NZ+1)=N1
      IF(IABC.NE.0) JVECT(NZ+1)=N1
      DUM21=(FAC21-1.D0)*AX*(-1.D0)*VELN
      CO(NZ+1)=DELTV2*DUM21                 !  advection contribution
      CO(NZ+1)=CO(NZ+1)-FLUXD*DELTV2        !  diffusion contribution
        dv2d21 = DELTV2R*DUM21
        dv2fld = DELTV2R*FLUXD
      DO IPRI=1,NPRI
         RHAND(N2,IPRI)= RHAND(N2,IPRI)-dv2d21*UTOLD(N1,IPRI)
         RHAND(N2,IPRI)= RHAND(N2,IPRI)+dv2fld*UTOLD(N1,IPRI)
      END DO
C
      IF(SL1(N2).le.sl1min)   CO(NZ+1)=0.D0
      IF (EVOL(N2).EQ. 0.0D0) CO(NZ+1)=0.D0
C
      NZ=NZ+1
C
C-------------------------DIAGONAL TERM IN EQUATION FOR ELEMENT N1
c
      DUM11=FAC12*AX*VELN
      CO(N1)=CO(N1)-DELTV1*DUM11          !  advection contribution
      CO(N1)=CO(N1)+FLUXD*DELTV1          !  diffusion contribution
        dv1d11 = DELTV1R*DUM11
      DO IPRI=1,NPRI
         RHAND(N1,IPRI)= RHAND(N1,IPRI)+dv1d11*UTOLD(N1,IPRI)
         RHAND(N1,IPRI)= RHAND(N1,IPRI)-dv1fld*UTOLD(N1,IPRI)
      END DO
C
C-------------------------DIAGONAL TERM IN EQUATION FOR ELEMENT N2
c
      DUM22=FAC21*AX*VELN
      CO(N2)=CO(N2)+DELTV2*DUM22          !  advection contribution
      CO(N2)=CO(N2)+FLUXD*DELTV2          !  diffusion contribution
        dv2d22 = DELTV2R*DUM22
      DO IPRI=1,NPRI
         RHAND(N2,IPRI)= RHAND(N2,IPRI)+dv2d22*UTOLD(N2,IPRI)
         RHAND(N2,IPRI)= RHAND(N2,IPRI)-dv2fld*UTOLD(N2,IPRI)
      END DO
C
C+++++END OF ASSIGNMENT OF INTERFACE TERMS+++++++++++++++++++++++
C
  200 CONTINUE
C-------------------------------------Modify matrix for inactive element
c
       DO 225 N=1,NEL
         IF (EVOL(N) .EQ. 0.0D0) CO(N)=1.0D0
         IF (SL1(N).le.sl1min) CO(N)=1.D0
225    CONTINUE
c
C-----------------------------------------------------------------------
C
      IF(MOPR(2).GE.1)   THEN
      IF(NOW.EQ.1)  THEN
         write (32,*)  '   iRN  ICN      CO------MATRIX-----'
         do nnn=1,NZ
            write (32,'(i6, i6,E12.4)') irn(nnn),icn(nnn),co(nnn)
         end do
      END IF
      END IF
c
C-----------------------------------------------------------------------
C
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE MATRIXC_Kdd(IPRI)
C
C
C*** Modified the aqueous transport matrix for the species with Kd and decay
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
      include 'perm_v2.inc'
C
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/E4/PHI(MNEL)
      common/e6/T(MNEL)
      COMMON/L1/IRN(mnz+1)
      COMMON/L2/ICN(mnz+1)
      COMMON/L3/CO(mnz+1)
      COMMON/L4/WKAREA(MNEQ*MNEL+10)
      COMMON/L7/JVECT(niwork)
      COMMON/AMMIS/MA,IPIV,U,IAB,NZ
      COMMON/AHTRAN/AHT(MNEL),STIME(MNEL),MLAGNR(MNEL),AMTT(MNEL)
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE6/SLOLD(MNEL)        ! old liquid saturation
      COMMON/SOLUTE8/SL1(MNEL)          ! new liquid saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)      ! old (initial) porosity
C
c porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
!
      common/chemgrid/c(mnod,maqt),utold(mnod,mpri),ut(mnod,mpri),
     & rhand(mnod,mpri),rsource(mnod,mpri),ph(mnod),gP(mnod,mgas),
     & aream(mnod,mmin),sads(mnod),psi(mnod),
     & d(mnod,mads),supadn(mnod,msurf),phip(mnod,msurf),
     & surfads(msurf),ub(mbound,mpri),ctot(mnod,mpri),cnfact
!
c---------------------------------COMMON blocks for Kd adsorption and decay
      common/kddca2/decayc(mpri)      ! decay constants
      common/kddca21/a_TDecay(mpri),  ! Thermal decay parameter, a
     &               b_TDecay(mpri)   ! Thermal decay parameter, b
      common/kddca3/kddp(mpri)        ! pointer to the primary species
      common/Kddca4/vkd(30,mpri)      ! values of Kd in initial zones
      common/Kddca5/izonekd(mnod)     ! Kd zone code
      common/Kddca6/sden(30,mpri)     ! solid density
      common/kddca8/kdflag(30,mpri)
c
c---------------------------------------------------------------------------
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
C
c---------------------------------------------------------------------------
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *MATRIXC_Kdd q2.2, 1999.6.21: MODIFY TRANPORT MATRIX FOR'
     X' SPECIES WITH Kd AND DECAY')
C
      DO 100 N=1,NEL
c
C------------------------------------COMPUTE RIGHT-HAND SIDE TERMS
c--solid density (kg/dm**3), and Kd(l/kg=mass/kg solid / mass/l water)
C
         KDDS=KDDP(IPRI)       ! Number in the species list for Kd and decay
         KDDZONE=IZONEKD(N)    ! Kd zone code
         IF (KDDZONE .LE. 0)  THEN
            SDEN2=0.0D0
            VKD2=1.0D0
            GO TO 300
         END IF
         SDEN2=SDEN(KDDZONE,KDDS)  ! solid density
         VKD2=VKD(KDDZONE,KDDS)    ! Kd value; or r factor if solid density=0
         IFLKD2=KDFLAG(KDDZONE,KDDS)
300      CONTINUE
         Dlamda=DECAYC(KDDS)       ! decay constant
!
!
!---------------
!.....If Dlamda is negative value, Dlamda is a temperature (or other
!.....values) dependent parameter. In this case call subroutine
!.....variable_decay to get decay constant
!---------------
!
      if (dlamda .lt. 0.0d0)     then
!
         a4  = a_TDecay(kdds)    ! Thermal decay parameter, a
         b4  = b_TDecay(kdds)    ! Thermal decay parameter, b
         tk4 = t(n) + 273.15d0   ! Absolute temperature
!
         CALL Variable_DecayConstant(a4, b4, tk4, dlamda)
!
      end if
!
!--------------------------------------------------------------------------------
!
!
         IF (KDDZONE .LT. 0) Dlamda=0.0D0
c
         IF (SDEN2.EQ.0.0D0)  THEN
c----------------If density is zero vkd2 is retardation factor
         IF (IFLKD2.EQ.0 .AND. VKD2.GE.1.0D0) THEN
            RETARD1=VKD2-1.0D0
         ELSE IF (IFLKD2.GT.0) THEN
            RETARD1=VKD2*AHT(N)-1.0D0
         END IF
C---------------------- Right-hand side terms contributed from Kd and decay
            R(N)=R(N)+PHIOLD(N)*SLOLD(N)*UTOLD(N,IPRI)*RETARD1  !for Retardion
C-------------------------Coefficient matrix contributed from Kd and decay
c            CO(N)=CO(N)+PHI(N)*SL1(N)*
            CO(N)=CO(N)+phisl1(n)*
     +                  (Dlamda*DELTEX*(1.0d0+RETARD1)+RETARD1)
         END IF
         
         IF (SDEN2.GT.0.0D0 .AND. VKD2.GE.0.0D0)  THEN
            R(N)=R(N) + (1.0D0-PHIOLD(N))*UTOLD(N,IPRI)*SDEN2*VKD2 ! for Kd
c            CO(N)=CO(N)+PHI(N)*SL1(N)*Dlamda*DELTEX+
            CO(N)=CO(N)+phisl1(n)*Dlamda*DELTEX+
     +            (1.0D0-PHI(N))*SDEN2*VKD2*(1.0d0+Dlamda*DELTEX)
         END IF
C
200      CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE Variable_DecayConstant(a4, b4, tk4, dlamda)
c
c
c**** calculate decay constant as a function of temperature or (other variables)****
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
!
      dlamda_ln = a4 - b4/tk4           ! ln(K) = a -b/T, ln-linear equation
      dlamda    = dexp(dlamda_ln)       ! 1/day
      dlamda    = dlamda/86400.0d0      ! 1/s
!
!
      return
!
      end
!
!
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE MATRIXG_implicit
C
C------SET UP THE COEFFICIENT MATRIX FOR GAS SPECIES TRANSPORT
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
      IMPLICIT REAL*8 (A-H,O-Z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'perm_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/E4/PHI(MNEL)
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
c
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),CWET(MAXMAT),
     +                SH(MAXMAT)
c.....multiphase tortuosity at each grid block
      common/tortmp/tortliq(mnel),tortgas(mnel)
c
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)
      COMMON/C4/DEL2(MNCON)
      COMMON/C5/AREA(MNCON)
      COMMON/C9/ELEM1(MNCON)
      COMMON/C10/ELEM2(MNCON)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
C
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
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)   ! darcy velocity
      COMMON/SOLUTE6/SLOLD(MNEL)          ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)          ! old gas saturation
      COMMON/SOLUTE8/SL1(MNEL)            ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)            ! new gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)        ! old (initial) porosity
      COMMON/SOLUTE11/NPRI,npaq,npads     ! number of chemical component
      COMMON/PARNP/NPL,NPG                ! specify in EOS module
C
      COMMON/AMMISC/IABC,ISOLVC
      COMMON/PRINTC/NOW              ! print control
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
C
C$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT IN BOTH LIQUID AND GAS PHASES $$$$$$$$$$
C
      COMMON/TRANGAS2/PFUGOLD2(NMNOD)     !old partial pressure ??
      COMMON/TRANGAS4/RHANDG(NMNOD)       ! right-han side for a gaseous species
      COMMON/TRANGAS7/EKGAS2,PFUGB2(100)
c     Need to save per gas species per grid block    
      COMMON/TRANGAS8/dcfgas(mgas,mnel),DIFUNG   ! gaseous species diffusion coefficient
C----------------------------------------------------------------------------------
C
      common/trangas3/pfug(mnod,mgas)  !New partial pressure =pfug(  )
      common/names/napri(mpri),naaqx(maqx),naaqt(maqt),namin(mmin),
     +   nagas(mgas),naexc(mexc),naads(mads)
      character*20 napri,naaqx,naaqt,namin,nagas,naexc,naads
      common/vmineral/pre(mnod,mmin),pre0(mnod,mmin),
     +  pinit(mnod,mmin+1)
!
      common/chemgrid/c(mnod,maqt),utold(mnod,mpri),ut(mnod,mpri),
     & rhand(mnod,mpri),rsource(mnod,mpri),ph(mnod),gP(mnod,mgas),
     & aream(mnod,mmin),sads(mnod),psi(mnod),
     & d(mnod,mads),supadn(mnod,msurf),phip(mnod,msurf),
     & surfads(msurf),ub(mbound,mpri),ctot(mnod,mpri),cnfact
!
      common/transport/izoneiw(mnod),izonebw(mnod),
     +   izonem(mnod),izoneg(mnod),izoned(mnod),izonex(mnod),
     +   izonpp(mnod)
c.....Added log10 of stimax
      common/constraints/sl1min,stimax,dlstmx
      common/oldtemp/tcold(mnod),tcmix(mnod)
c.....Add temperature block
      common/heattran/tc(mnod),tkelv(mnod)
c----------------- gaseous species properties (mol wt and mol diam)
      common/gasprop/dmwgas(mgas),diamol(mgas)
C
c............porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
C$$$$$ COMMON BLOCKS FOR SINKS/SOURCES $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/G4/ELEG(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G12/LCOM(MNOGN)
C
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
      COMMON/GASLAW/RGAS,AMS,AMA,CVGAS
      COMMON/DFM/TIMAX,REDLT
      COMMON/BC/NELA
      common/gasindx/kgas
      CHARACTER*5 ELEM,ELEM1,ELEM2,ELEG
C
C---------------------------------- For coupling with reactive transport
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
C
c----------------------------------------Indicators from EOS module
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
C
C-------------------------------------------------------------------------
c
c-----Local T and P variables and gas diffusion coefficients
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **MATRIXG q2.1, 1999.4.12: ASSEMBLE MATRIX FOR TRANSPORT'
     X' IN LIQUID AND GAS PHASES**********')
C
c     EKGAS2 is used as gas log K to compute gas aqueous concentration
c     during gas transport.  Not working for now. Need to set to zero.
        EKGAS2=0.0D0
C-----------------------------------------------------------------------
C$$$$$FOR IABC=0 COLUMN INDICES WILL BE STORED IN ICN, OTHERWISE IN JVECT
C-----INITIALIZE COUNTER FOR MATRIX ELEMENTS.
      NZ=0
C
C******************* LOOP OVER ELEMENTS ******************************
C
C*****COMPUTE ALL QUANTITIES WHICH DEPEND ONLY UPON VARIABLES PERTAINING
C     TO ONE VOLUME ELEMENT.
C
      DO 100 N=1,NEL
         rt = gc*tkelv(n)
         rtold = gc*(tcold(n)+273.15d0)
C------------------------------------COMPUTE RIGHT-HAND SIDE TERMS
C-------------------CONTRIBUTED FROM INITIAL CONDITIONS
C
         DUM1=PHIOLD(N)*(SGOLD(N)/rtold)
         RHANDG(N)=DUM1*PFUGOLD2(N)
C
C-----------------------------COMPUTE COEFFICIENT MATRIX ELEMENT
      IRN(NZ+1)=N                     ! row index
      IF(IABC.EQ.0) ICN(NZ+1)=N       ! column index
      IF(IABC.NE.0) JVECT(NZ+1)=N
      CO(NZ+1)=phisg1(n)/rt
      NZ=NZ+1
C
C+++++++++END OF ASSIGNMENT OF ONE-ELEMENT TERMS++++++++++++++++++++
C
  100 CONTINUE
C
C
C********* LOOP OVER CONNECTIONS ***********************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND UPON VARIABLES FOR TWO VOLUME
C     ELEMENTS ("INTERFACE QUANTITIES").
C
      DO 200 N=1,NCON
C
      N1=NEX1(N)         ! the first element of a connection
      N2=NEX2(N)         ! the second element of a connection
      rt1=gc*tkelv(n1)   ! gas RT at n1
      rt2=gc*tkelv(n2)   ! gas RT at n2
         rt12=(rt1+rt2)/2.d0         ! gas RT connection
C
C-----OBTAIN SOME QUANTITIES PERTAINING TO CONNECTION.
C
      AX=AREA(N)         ! interface area
      D1=DEL1(N)
      D2=DEL2(N)
      DISTAN=D1+D2
!
!.....Added new formulation for diffusion weighting
!.....Harmonic weighting of irregular grids, including tortuosity -----
!
       if (sg1(n1).le.sl1min.or.sg1(n2).le.sl1min)   then
         difs = 0.d0
       else
          psgtd1 = phisg1(n1)*tortgas(n1)*dcfgas(kgas,n1)
          psgtd2 = phisg1(n2)*tortgas(n2)*dcfgas(kgas,n2)
          difs = (distan*psgtd1*psgtd2)/(d1*psgtd2+d2*psgtd1)
       end if
!
!.....Harmonic mean of product
!
!.....Distance weighting for tortuosity
!
      FLUXDG=AX*DIFS/DISTAN   ! a term derived from diffusion for gas
C
C+++++ASSIGN ALL INTERFACE TERMS+++++++++++++++++++++++++++++++++
!.....Liquid phase fluxes not used in gas transport
!.....Assign upstream weighting parameters for gas phase
!TX05/18/2010, for EOS9 gas velocity is zero
      if (ieos .eq. 9)   then
         VELNG = 0.0d0            ! gas Darcy velocity
                         else
         NP    = NPG
         NI    = (N-1)*NPH+NP
         VELNG = VELDAR(NI)       ! gas Darcy velocity
      end if
C---------------------------------------
      IF (IEOS.EQ.9)    VELNG=0.0D0       ! no advection for EOS9 flow module
C---------------------------------------
      IF (VELNG.GE.0.D0)    THEN
         FAC12G=1.D0-WUPC                 ! up-stream weighting fatcor
         FAC21G=WUPC
                          ELSE
         FAC12G=WUPC
         FAC21G=1.D0-WUPC
      END IF
599   CONTINUE
C
      WTIMEG=1.d0                         ! implicit
      EVOLN1=EVOL(N1)
      IF (EVOL(N1).EQ.0.0D0)  EVOLN1=1000.0D0
      DUM1=DELTEX/EVOLN1
      DELTV1=WTIMEG*DUM1
C
      EVOLN2=EVOL(N2)
      IF (EVOL(N2).EQ.0.0D0)  EVOLN2=1000.0D0
      DUM2=DELTEX/EVOLN2
      DELTV2=WTIMEG*DUM2
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N1, ARISING FROM N2
      IRN(NZ+1)=N1
      IF(IABC.EQ.0) ICN(NZ+1)=N2
      IF(IABC.NE.0) JVECT(NZ+1)=N2
c
c---------- skip if node is saturated or zero volume
      if(sg1(n1).ge.sl1min.and.evol(n1).gt.0.d0)then
        DUM12=AX*((FAC12G-1.D0)/rt12*VELNG)
        CO(NZ+1)=DELTV1*DUM12                  ! advection contribution
        CO(NZ+1)=CO(NZ+1)                      ! diffusion contribution
     +                -FLUXDG/rt12*DELTV1
      else
        CO(NZ+1)=0.D0
      endif
C
      NZ=NZ+1
C
C-----OFF-DIAGONAL TERM IN EQUATION FOR ELEMENT N2, ARISING FRO N1
      IRN(NZ+1)=N2
      IF(IABC.EQ.0) ICN(NZ+1)=N1
      IF(IABC.NE.0) JVECT(NZ+1)=N1
c
c---------- skip if node is saturated or zero volume
      if(sg1(n2).ge.sl1min.and.evol(n2).gt.0.d0)then
        DUM21=AX*((FAC21G-1.D0)/rt12*(-1.d0)*VELNG)
        CO(NZ+1)=DELTV2*DUM21       !  advection contribution
        CO(NZ+1)=CO(NZ+1)           !  diffusion contribution
     +                    -FLUXDG/rt12*DELTV2
      else
        CO(NZ+1)=0.D0
      endif
C
      NZ=NZ+1
C
C-------------------------DIAGONAL TERM IN EQUATION FOR ELEMENT N1
      DUM11=AX*(FAC12G/rt1*VELNG)   !  modified gc
      CO(N1)=CO(N1)-DELTV1*DUM11    !  advection contribution
      CO(N1)=CO(N1)                 !  diffusion contribution
     +               + FLUXDG/rt1*DELTV1
C
C-------------------------DIAGONAL TERM IN EQUATION FOR ELEMENT N2
      DUM22=AX*(FAC21G/rt2*VELNG)
      CO(N2)=CO(N2)+DELTV2*DUM22    !  advection contribution
      CO(N2)=CO(N2)                 !  diffusion contribution
     +               + FLUXDG/rt2*DELTV2      !Modified gc
C
C+++++END OF ASSIGNMENT OF INTERFACE TERMS+++++++++++++++++++++++
C
  200 CONTINUE
C
C
C----------------Modify matrix for inactive element and gas saturation=0
C
       DO 225 N=1,NEL
         IF (SG1(N) .le. SL1MIN)   CO(N)=1.D0
         IF (EVOL(N) .EQ. 0.0D0)   CO(N)=1.0D0
225    CONTINUE
c
C----------------------------------------------------------------------
C
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
c
      SUBROUTINE MAX_DELT
C
C********Calculate maximum time step imposed by Courant number *****
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
      implicit double precision (A-H,O-Z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/PARNP/NPL,NPG          ! specify in EOS module
      COMMON/PORVEL/VEL(MNPH*MNCON)
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)
      COMMON/C4/DEL2(MNCON)
      COMMON/TRANGAS8/dcfgas(mgas,mnel),DIFUNG     ! gaseous species diffusion coefficient
      COMMON/SOLUTE9/SG1(MNEL)           ! new gas saturation
      COMMON/E3/EVOL(MNEL)
      COMMON/SOLUTE8/SL1(MNEL)           ! new liquid saturation 
      COMMON/TRANGAS9/NGAS1              ! number of gaseous species
      COMMON/DIFUNT_L1/DIFUNT_L(MNCON)
      COMMON/E2/MATX(MNEL)
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
c----------------------------------------Indicators from EOS module
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
c.....Time-step limit data
      COMMON/E1/ELEM(MNEL)
      common/dtlim/max_chem_it,delt_conne,id_chem
      character*5 elem
      character*16 delt_conne
      character*5 id_chem
c
C
C######################################################################
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ****MAX_DELT 1.0, 2003.7.30: CALCULATE MAXIMUM TIME STEP'
     x' BASED ON COURANT NUMBER************')
C
C********* LOOP OVER CONNECTIONS ***********************************
C
         dt_gasliq = delt
c
      DO 200 N=1,NCON
c
         delt_aq = 1.d50     ! max dt for liquid
         delt_gas = 1.d50    ! max dt for gas
c
         N1=NEX1(N)          ! the first element of a connection
         N2=NEX2(N)          ! the second element of a connection
c
c chemistry because internodal distance < d1min
c         if(iskip1(n1).ne.1.and.iskip1(n2).ne.1)goto 200
c
c---------- moved block below from after the if statements
         D1=DEL1(N)
         D2=DEL2(N)
         IF (EVOL(N1).EQ.0.0D0) D1=D2
         IF (EVOL(N2).EQ.0.0D0) D2=D1
         DISTAN=D1+D2       ! distance between blocks N1 and N2
c
c---------- skip courant limitation for nodes where we skip
c
         if (sl1(n1).gt.sl1min.or.sl1(n2).gt.sl1min)then
c
            NP=NPL
            NI=(N-1)*NPH+NP
            VELN=VEL(NI)                ! liquid pore velocity
            AVELN=ABS(VELN)
c
c...  Aqueous phase only Courant limited deltat, D=0
c
            if (aveln.gt.1.0D-15)then
               delt_aq=(distan*dabs(rcour))/aveln
               if (delt_aq.lt.dt_gasliq)  then
                  delt_conne(1:5)=elem(n1)
                  delt_conne(6:10)=elem(n2)
                  delt_conne(11:16)=' V_liq'
                  dt_gasliq=delt_aq
               end if
            end if
         end if
c
c.....Gas phase test
c
      if (rcour.gt.0.0d0.and.ngas1.gt.0)then
c
c..... For gas-dominant system
c
        if(sg1(n1).gt.sl1min.and.sg1(n2).gt.sl1min)then
c
!TX05/18/2010; for EOS9, gas velocity is zero 
           if (ieos .eq. 9)   then
              AVELN = 0.0d0          
                              else
              NP    = NPG
              NI    = (N-1)*NPH+NP
              VELN  = VEL(NI)                    ! gas Darcy velocity
              AVELN = ABS(VELN)
           end if
c
c... Gas phase Courant-limited delta t only (gas velocity > 0, D=0)
           if (aveln.gt.1.0D-15)then
              delt_gas=(distan*dabs(rcour))/aveln
              if (delt_gas.lt.dt_gasliq) then
                 delt_conne(1:5)=elem(n1)
                 delt_conne(6:10)=elem(n2)
                 delt_conne(11:16)=' V_gas'
                 dt_gasliq=delt_gas
              end if
           end if
        end if
c
       end if
c
c
  200 CONTINUE
c
c... Find minimum dt
c
       if (deltmx.gt.0.d0)then
          delt = min(dt_gasliq,delt,deltmx)
       else
          delt = min(dt_gasliq,delt)
       end if
          delt = max(delt,0.01d0)
c
c......save type of time step limitation
c
       if (delt.eq.deltmx.and.deltmx.gt.0.d0)then
          delt_conne = 'Max_Delta_t     '
       else if(delt.eq.0.01d0)then
          delt_conne = 'Min_Delta_t     '
       end if
c
      RETURN
      END
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE DRY_MAP
C
C*********************** May dry grid block and calculate solid amount **********
C
C
C***** N O T A T I O N ********************
C
      implicit double precision (A-H,O-Z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
C
C$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     THESE BLOCKS HAVE A LENGTH OF NEL (= NUMBER OF ELEMENTS)
C
      COMMON/E3/EVOL(MNEL)
c
C$$$$$$$$$ COMMON BLOCKS FOR CONNECTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$C
C
      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C5/AREA(MNCON)
C
C$$$$$$$$$ COMMON BLOCK FOR SECONDARY VARIABLES $$$$$$$$$$$$$$$$$$$$
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C$$$$$$$$$ COMMON BLOCKS FOR SOLUTE TRANSPORT $$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON)   ! darcy velocity
      COMMON/SOLUTE6/SLOLD(MNEL)          ! old liquid saturation
      COMMON/SOLUTE8/SL1(MNEL)            ! new liquid saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)        ! porosity from previous time step
      COMMON/PARNP/NPL,NPG                ! specify in EOS module
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
C
C######################################################################
c
C.....Interface area reduction factor
      common/afactor/a_fm2(mncon)  ! advection area reduction (flow from F to M)
C
C-----------------------Common blocks for dryout grid blocks
C
       COMMON/DRYOUT/IDRY(MNOD),ADRY(MNOD,MPRI)
C
C######################################################################
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ****DRY_MAP q2.3, 2001.3.26: Map dry-out grid blocks***')
C
      DO N=1,NEL
         IDRY(N)=0
      END DO
c
c.....Note that ADRY is moles of species i per dm3 of medium
      DO N=1,NEL
      DO I=1,NPRI
         ADRY(N,I)=0.0D0
      END DO
      END DO

c********** LOOP OVER GRID BLOCKS ************************************
c Add separate loop over each grid block for storage terms, and remove these
c terms from the loop over each connection
      do n=1,nel
c
        if (slold(n).gt.sl1min.and.sl1(n).le.sl1min)  then
          idry(n)=1
          do i=1,npri
            adry(n,i)=adry(n,i)+(slold(n)-sl1(n))*
     +         phiold(n)*utold(n,i)       ! Use UTOLD - same but cleaner
          end do
        end if
      end do
c
C********* LOOP OVER CONNECTIONS ***********************************
C
C-----COMPUTE ALL QUANTITIES WHICH DEPEND UPON VARIABLES FOR TWO VOLUME
C     ELEMENTS ("INTERFACE QUANTITIES").
C
      DO 200 N=1,NCON
C
      N1=NEX1(N)         ! the first element of a connection
      N2=NEX2(N)         ! the second element of a connection
      IF (SL1(N1).gt.sl1min.AND.SL1(N2).gt.sl1min)  GO TO 200
      IF (SL1(N1).le.sl1min.AND.SL1(N2).le.sl1min)  GO TO 200
C
      IF (SL1(N1).le.sl1min)  IDRY(N1)=1
      IF (SL1(N2).le.sl1min)  IDRY(N2)=1
C
C-----OBTAIN SOME QUANTITIES PERTAINING TO CONNECTION.
C
c*************************************************For active fracture model
      AX1 = AREA(N)
      a_fm=a_fm2(n)
      AX=a_fm*AX1      ! advection interface area
c**************************************************************************
c
C+++++ASSIGN ALL INTERFACE TERMS+++++++++++++++++++++++++++++++++
      NP=NPL
      NI=(N-1)*NPH+NP
      VELN=VELDAR(NI)           ! liquid Darcy velocity
C
      FL=DELTEX*AX*VELN         ! water flux (equivalent to m3 of water)
      FL1=FL/EVOL(N1)           ! water flux/m*3 medium (equivalent to phi*sat)
      FL2=FL/EVOL(N2)           ! water flux/m*3 medium
C
        IF (VELN.GT.0.D0.and.sl1(n1).le.sl1min)  THEN
           DO I=1,NPRI
              ADRY(N1,I)=ADRY(N1,I)+FL1*UT(N2,I)
           END DO
        ELSE IF(VELN.LT.0.d0.and.sl1(n2).le.sl1min) THEN
           DO I=1,NPRI
              ADRY(N2,I)=ADRY(N2,I)-FL2*UT(N1,I)
           END DO
        END IF
C
  200 CONTINUE
c
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE DRY_MIN
C
C*********************** Calculate the amount of precipitation at dry grid-blocks
C----The order of mineral precipitation is assigned in CHEMICAL.INP
C
C
C***** N O T A T I O N ********************
C
      implicit double precision (A-H,O-Z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
C
C######################################################################
C
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
c
C-----------------------Common blocks for dryout grid blocks
C
      COMMON/DRYOUT/IDRY(MNOD),ADRY(MNOD,MPRI)
      COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)    ! residual in precipitates
      common/dry_salt/nsalt,isalt(mmin)
      double precision Acom(MPRI)
C
C######################################################################
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' DRY_MIN q2.4, 2001.4.17: Calculate mineral precipitation'
     X' at dry grid blocks**********')
c
C-----------------------------------------Initialization
c
c... Initialize adryr in beginning only
c
      IF(ICALL.EQ.1) THEN
         do n=1,nel
         do i=1,npri
c...........set adryr to previous value for restart
            adryr(n,i) = adryr0(n,i)
         end do
         end do
      endif
c
         Ih  = 0
         Ih2o = 0
         do i=1,npri
            if(napri(i).eq.'h+') Ih=i
            if(napri(i).eq.'h2o') Ih2o=i
         end do
c
C-----------------------------------------------------------------------
c
c... Update residual moles of species and initialize dry mineral moles
      do n=1,nel
         do i=1,npri
           adryr(n,i) = adryr0(n,i) + adry(n,i)
         end do
c------- Initialization of drypre
         do m = 1, nmin
           drypre(n,m) = 0.d0
         end do
      end do
c-----------------------------------------------------------------------
c
      DO 100 N=1, NEL
c         IF(IDRY(N).EQ.0) GO TO 100
         IF(IDRY(N).gt.0) then
c
         do j=1,npri
            Acom (j) = ADRY(N,j)
         end do
c
         do 200 m=1,nsalt
            m1=isalt(m)
c                     m1 should always be nonzero       if(m1.gt.0)then
             ncp=ncpm(m1)
c
             damin = 1.0d+10
             do 150 k=1,ncp
               j=icpm(m1,k)
               if (j .eq. ih)   go to 150
               if (j .eq. ih2o) go to 150
               if (Acom(j) .le. 0.0d0) go to 200
               Acmin = Acom(j)/stqm(m1,k)
               if (Acmin .lt. damin) damin = Acmin
150          continue
c
            if(damin.lt.9.9d9)drypre(N,m1) = drypre(N,m1) + damin
c
            if(damin.gt.0.d0.and.damin.lt.9.d9)then
              do 160 k=1,ncp
                 j=icpm(m1,k)
                 if (j .eq. ih)   go to 160
                 if (j .eq. ih2o) go to 160
                 Acom (j) = Acom (j) - damin*stqm(m1,k)
                 ADRYR(N,j) = ADRYR(N,j) - damin*stqm(m1,k)
160           continue
            endif
c
 200     continue
c
        endif
c
100   CONTINUE
c
C-----------------------------------------------------------------------
C
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE CONVER1
C
C-----THIS SUBROUTINE IS CALLED AFTER A STEADY STATE FLOW IS REACHED
C
C------for define next time step size
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      COMMON/KC/KC
      COMMON/DFM/TIMAX,REDLT
      COMMON/SVZ/NOITE,MOP(24)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/DLT/NDLT,DLT(100)
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
C      COMMON/TIMES/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
        COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX    ! TOUGH2 V2
      COMMON/ITERATION/ITERFL,ITERMOD,ITERTR,ITERCH,
     1         MAXITCH,ITERAD,MAXITAD,Iremove
      COMMON/RITERATION/AVERITCH,AVERITAD,COUNTAD,countch
      COMMON/ICONVERGENCE/ MAXITPFL,MAXITPTR,MAXITPCH,MAXITPAD
c
      COMMON/STEADY/IFLOWSS,JSTEADY
      common/mexhaust/nexh       ! number of minerals exhausted
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
C
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***CONVER1 1.0, 2003.7.4: UPDATE PRIMARY VARIABLES AFTER'
     X' CONVERGENCE IS ACHIEVED (Used after steady-state flow is'
     X' reached)**********')
C
      SUMTIM=SUMTIM+DELTEX
      IF(TIMAX.NE.0.D0.AND.TIMAX.EQ.SUMTIM) NOWTIM=1
C-----AFTER CONVERGENCE UPDATE TOTAL TIME AND ASSIGN NEW TIME STEP.
      IF(NDLT.EQ.0) GOTO20
       IF (KC+1.GT.8*NDLT) GO TO 20
       IF (KC.LE.99) THEN
          IF (DLT(KC+1).NE.0.0D0) THEN
             DELT=DLT(KC+1)
             GO TO 22
                                  ELSE
             GO TO 20
          END IF
       END IF
C-----IF NO FURTHER TIME STEP INSTRUCTIONS ARE PROVIDED, KEEP
C     GOING WITH LAST TIME STEP.
C-----COME HERE FOR NEW TIME STEP ASSIGNMENT.
      DELT=DLT(KC+1)
      GOTO 22
   20 DELT=DELTEX
c
c......Made time stepping for steady state chemistry consistent with CONVER3
       if(mopr(16).lt.1) then  
        if(maxitch.ge.maxitpch)delt=max(0.5d0*deltex,1.d-2)
        if(maxitch.ge.100.and.maxitch.lt.maxitpch)delt=
     +      max(0.60d0*deltex,1.d-2)
        if(maxitch.lt.100.and.maxitch.ge.75)delt=
     +      max(0.8d0*deltex,1.d-2)
        if(maxitch.lt.75.and.maxitch.ge.50)delt=
     +      max(deltex,1.d-2)
        if(maxitch.lt.50.and.maxitch.ge.30)delt=
     +      max(1.5d0*deltex,1.d-2)
        if(maxitch.lt.30)delt=max(2.0d0*deltex,1.d-2)
       end if 
c
      IF (JSTEADY.EQ.1)   DELT=0.01D0*DELTEX
c
c......Set minimum to 0.01 second
       delt = max(delt,0.01d0)
c
199   CONTINUE
c
   22 IF (TIMAX.NE.0.D0)   DELT=MIN(DELT,TIMAX-SUMTIM)
      IF (DELTMX.NE.0.D0)  DELT=MIN(DELT,DELTMX)
c
c......Set minimum to 0.01 second
       delt = max(delt,0.01d0)
C
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE CONVER2
C
C-----Time step for both flow and reactive transport
c
C
C-----THIS SUBROUTINE IS CALLED AFTER SUCCESFULL COMPLETION OF
C     A TIME STEP.
C     IT UPDATES PRIMARY VARIABLES, AND DEFINES THE NEXT TIME STEP.
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
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
     +            GK(MAXMAT)
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
      COMMON/ITERATION/ITERFL,ITERMOD,ITERTR,ITERCH,
     1         MAXITCH,ITERAD,MAXITAD,Iremove
      COMMON/RITERATION/AVERITCH,AVERITAD,COUNTAD,countch
      COMMON/ICONVERGENCE/ MAXITPFL,MAXITPTR,MAXITPCH,MAXITPAD
c
      common/mexhaust/nexh       ! number of minerals exhausted
C
      COMMON/Rdphi/RPHI(mnel)    ! porosity change rate (unit time)
C
      DIMENSION DXM(10)
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***CONVER2 1.0, 2003.7.4: UPDATE PRIMARY VARIABLES AFTER'
     X' CONVERGENCE IS ACHIEVED (Used when coupled with reactive'
     x' geochemistry)************')
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
c
C*******************************Porosity change due to mineral dis/pre
c
cc         PHI(N)=PHI(N)+RPHI(N)*DELTEX
c
C*********************************************************************
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
C
      SUMTIM=SUMTIM+DELTEX
      IF(TIMAX.NE.0.d0.AND.TIMAX.EQ.SUMTIM) NOWTIM=1
C-----AFTER CONVERGENCE UPDATE TOTAL TIME AND ASSIGN NEW TIME STEP.
      IF(NDLT.EQ.0) GOTO20
       IF (KC+1.GT.8*NDLT) GO TO 20
       IF (KC.LE.99) THEN
          IF (DLT(KC+1).NE.0.0D0) THEN
             DELT=DLT(KC+1)
             GO TO 22
                                  ELSE
             GO TO 20
          END IF
       END IF
C-----IF NO FURTHER TIME STEP INSTRUCTIONS ARE PROVIDED, KEEP
C     GOING WITH LAST TIME STEP.
C-----COME HERE FOR NEW TIME STEP ASSIGNMENT.
      DELT=DLT(KC+1)
      GOTO 22
   20 DELT=DELTEX
c--------------------------------------------------------
      IF(ITER.LE.MOP(16)) DELT=2.d0*DELTEX
   22 IF(TIMAX.NE.0.d0) DELT=MIN(DELT,TIMAX-SUMTIM)
      IF(DELTMX.NE.0.d0) DELT=MIN(DELT,DELTMX)
      delt = max(delt,0.01d0)
C
      RETURN
      END
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE CONVER3
C
C-----Time step for reactive transport
c
C
C-----THIS SUBROUTINE IS CALLED AFTER SUCCESFULL COMPLETION OF
C     A TIME STEP.
C     IT UPDATES PRIMARY VARIABLES, AND DEFINES THE NEXT TIME STEP.
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
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
      COMMON/ITERATION/ITERFL,ITERMOD,ITERTR,ITERCH,
     1         MAXITCH,ITERAD,MAXITAD,Iremove
      COMMON/RITERATION/AVERITCH,AVERITAD,COUNTAD,countch
      COMMON/ICONVERGENCE/ MAXITPFL,MAXITPTR,MAXITPCH,MAXITPAD
      COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
c
      integer*8 max_chem_it
      common/dtlim/max_chem_it,delt_conne,id_chem
      character*16 delt_conne
      character*5 id_chem
      double precision delt_ch
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ****CONVER3 1.0, 2006.9.15: Adjust timestep for reactive'
     x' transport***********')
C
      delt_ch = delt
!
cc      if (delt_conne .eq. 'Min_Delta_t     ') then
cc         mopr(16) = 1     ! No timne step reduction due to chemistry iteration
cc      end if
!
cc      if(iter.le.mop(16))then
      if(iter.le.mop(16).and.mopr(16).lt.1)then
        if(maxitch.ge.maxitpch)delt_ch=max(0.5d0*deltex,1.d-2)
        if(maxitch.ge.100.and.maxitch.lt.maxitpch)delt_ch=
     +      max(0.60d0*deltex,1.d-2)
        if(maxitch.lt.100.and.maxitch.ge.75)delt_ch=
     +      max(0.8d0*deltex,1.d-2)
        if(maxitch.lt.75.and.maxitch.ge.50)delt_ch=
     +      max(deltex,1.d-2)
        if(maxitch.lt.50.and.maxitch.ge.30)delt_ch=
     +      max(1.5d0*deltex,1.d-2)
        if(maxitch.lt.30)delt_ch=max(2.0d0*deltex,1.d-2)
      else if(iter.gt.mop(16))then
        delt_ch=max(deltex,1.d-2)
      end if
c--------------------------------------------------------
      if (iter.gt.mop(16).and.delt_ch.eq.delt)then
        delt_conne='Flo_Delta_t     '
      elseif(iter.le.mop(16).and.delt_ch.lt.delt)then
        delt_conne='Chm_Delta_t     '
      elseif(iter.le.mop(16).and.delt_ch.eq.delt)then
        delt_conne='Flo_Delta_t     '
      elseif(timax-sumtim.lt.delt.and.timax-sumtim.lt.delt_ch)then
        delt_conne='Tot_time_Dt     '
      end if
c
      IF (TIMAX.NE.0.d0)  DELT=MIN(DELT,delt_ch,TIMAX-SUMTIM)
      IF (DELTMX.NE.0.d0) DELT=MIN(DELT,DELTMX)
      delt = max(delt,0.01d0)
C
      if (delt.eq.deltmx) delt_conne = 'Max_Delta_t     '
      if (delt.eq.0.01d0) delt_conne = 'Min_Delta_t     '
c
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE TSTEPT1
C
C-----THIS ROUTINE MODIFIES TIME STEPS TO COINCIDE WITH USER-
C     INPUT VALUES AT WHICH PRINTOUT IS DESIRED.
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
c      COMMON/TIMES/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX
      COMMON/tim/ITI,DELAF,ITPR,TIS(100),ITCO,NOWTIM,DELTMX    ! TOUGH2 V2
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/MODELT/IDELT   ! IDELT=1 if DELT modified by TSTEPT1
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' TSTEPT 1.0, 1991.3.4: ADJUST TIME STEPS TO COINCIDE WITH'
     X' USER-DEFINED TARGET TIMES********')
C
      IF(TIS(ITI).LE.SUMTIM) GOTO 100
C
C
C-----FIND TIME FOLLOWING SUMTIM.
      DO1 I=1,ITI
      IF(TIS(I).EQ.SUMTIM) GOTO 4
      IF(TIS(I).LT.SUMTIM) GOTO 1
      GOTO 2
    1 CONTINUE
C
    2 IF(SUMTIM+DELT.LT.TIS(I)) GOTO 10
C
C-----COME HERE TO ADJUST DELT.
      DELT=TIS(I)-SUMTIM
cels4/8/09
      delt = max(delt,0.01d0)
c
      IDELT=1   ! DELT modified by TSTEPT1
      NOWTIM=1
      GOTO 10
C
C-----COME HERE AFTER NEXT TIME HAS BEEN REACHED.
    4 ITCO=I
C
      IF(DELAF.GT.0.D0)
     ADELT=MIN(DELT,DELAF)
      IDELT=1   ! DELT modified by TSTEPT1
C
      delt = max(delt,0.01d0)
c
      IF(TIS(I+1).GT.TIS(I)+DELT) RETURN
      DELT=MIN(DELT,TIS(I+1)-TIS(I))
      delt = max(delt,0.01d0)
c
      IDELT=1   ! DELT modified by TSTEPT1
      NOWTIM=1
C
      RETURN
   10 ITCO=I-1
      delt = max(delt,0.01d0)
      RETURN
C
  100 ITCO=ITI
      delt = max(delt,0.01d0)
      RETURN
C
      END
c
c
c
c-------------------------------------------------------------------------------
c
c
c
      SUBROUTINE LINEQC(neqc)
C
C-----------------------for TOUGH2 V2
C
      IMPLICIT double precision (A-H,O-Z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
C
C-----THIS SUBROUTINE CALLS THE LINEAR EQUATION SOLVER T2CG2.
C     IT HAS LOGIC TO HANDLE FAILURES IN LINEAR EQUATION SOLUTION.
C
C     AFTER SOLUTION, LATEST UPDATED ITERATES ARE OBTAINED FOR
C     ALL PRIMARY DEPENDENT VARIABLES.
C
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
      COMMON/AMMISC/IABC,ISOLVC
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
      SAVE ICALL,N,iteruc,iprpro,ichino
      DATA ICALL,iteruc,iprpro,ichino/0,0,0,0/
C
C=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>   Main body of LINEQ
C
      iiuunn = iunit
      nnvvcc = nvectr
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***LINEQC 2.00, 2001.5.26: Interface for linear equation'
     &' solvers T2CG2, can call a direct solver or a package of'
     &' conjugate gradient solvers')
C
      N = NEQC
c
C*********************************************************************
C*                                                                   *
C*             ACCOUNTING FOR ZEROs ON THE MAIN DIAGONAL             *
C*                                                                   *
C*********************************************************************
C
C
      ichino = 0
      ichino = 0
      izerod = 0
      IF(ISOLVC.NE.6) THEN
          NNF = NEQC
c
         DO 20 I=1,NNF
c    IF(IRN(I).EQ.ICN(I).AND.ABS(CO(I)).EQ.0.0e0) THEN
            IF(IRN(I).EQ.ICN(I).AND.ABS(CO(I)).EQ.0.0d0) THEN
               izerod = izerod+1
               ichind = 1
            ELSE
               ichind = 0
            END IF
   20    CONTINUE
      END IF
C
      IF(ichind.NE.ichino) WRITE(15,6003) kcyc,iter,izerod,zprocs
      ichino = ichind
C
         zertio=0.d0
      IF(izerod.GT.0) THEN
         IF(zprocs.EQ.'Z1') iprpro = 1
         IF(zprocs.EQ.'Z2') iprpro = 2
         IF(zprocs.EQ.'Z3') iprpro = 3
         IF(zprocs.EQ.'Z4') iprpro = 4
         IF(oprocs.NE.'O0') THEN
            IF(zprocs.EQ.'Z0'.OR.zprocs.EQ.'Z1') THEN
               zerod0 = 1.0d2*izerod
c               dilong = 1.0d0*NEQ*NELA
               dilong = 1.0d0*NEQC
               zertio = zerod0/dilong
               IF(zertio.GT.2.0d-1) THEN
                  zerpro=zertio*1.d2
                  WRITE (34,6004) zerpro,oprocs
                  oprocs='O0'
               END IF
            END IF
         END IF
      END IF
C
 6003 FORMAT(/,T2,'At KCYC=',i5,' and ITER=',i5,', IZEROD=',I5,
     &              ' and ZPROCS = ',a2)
c
c-------Can't pass microsoft develop studio and unix compiler
c-------------------------------------------------
 6004 FORMAT(//,' ',15('WARNING-'),//T14,F5.2,'% of the elements on ',
     &         'the main diagonal of the Jacobian matrix are zeros'/
     &,'The matrix preprocessor OPROCS = ',A2,' cannot be used'/
C     &,'Action taken: reset OPROCS to \'0\' (no O-preprocessing);',
     & ' continue execution.')
C
c      INUM  = 0
      IGOOD = 0
C
C
C*********************************************************************
C*                                                                   *
C*      DETERMINATION OF THE BANDWIDTHS OF THE U AND L MATRICES      *
C*      OR PLACEMENT OF THE ELEMENTS INTO THE CG SOLUTION ARRAY      *
C*                                                                   *
C*********************************************************************
C
      co(mnzp1) = 0.0d0
      r(mnetp1) = 0.0d0
C
  395 IF(ISOLVC.EQ.6.AND.icall.EQ.1) THEN
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
         navdia = nnnbig/matord
C
         IF(ntotd.GT.navdia) THEN
            WRITE (34,6020) ntotd,navdia
            STOP
         END IF
      END IF
C
C
C*********************************************************************
C*                                                                   *
C*                         MATRIX SOLUTION                           *
C*                                                                   *
C*********************************************************************
C
      IF(ISOLVC.EQ.6) THEN
         CALL LUBAND(matord,nz,nsubdg,nsupdg,ntotd,ab,
     &               matord,r,JVECT,info)
C
      ELSE
C
         IF(izerod.NE.0.AND.NEQ.GT.1.AND.zprocs.NE.'Z0') THEN
            CALL MTRXIN(iprpro)
            GO TO 415
         END IF
c
         NEQC=1       !!!!!!!!!!!!!!!
         IF(oprocs.NE.'O0'.AND.NEQC.EQ.1) GO TO 415
         IF(oprocs.eq.'O0') THEN
            GO TO 415
         ELSE IF(oprocs.eq.'O1') THEN
            ione = 1
            CALL MTRXPR(ione)
         ELSE IF(oprocs.eq.'O2') THEN
            itwo = 2
            CALL MTRXPR(itwo)
         ELSE IF(oprocs.eq.'O3') THEN
            ithree = 3
            CALL MTRXPR(ithree)
         ELSE IF(oprocs.eq.'O4') THEN
            ifour = 4
            CALL MTRXIN(ifour)
         END IF
C
  415    NNZ = 0
         DO 420 I=1,NZ
            NNZ      = NNZ+1
            co(NNZ)  = CO(I)
            irn(NNZ) = IRN(I)
            icn(NNZ) = ICN(I)
  420    CONTINUE
C
C
         DO 440 i=1,n
            wkarea(i) = 0.0d0
  440    CONTINUE
C
         IF(MOP(6).NE.0) CALL THYME(0,TS,TT)
C
         IF(ISOLVC.EQ.2) THEN
            CALL DSLUBC(N,r,wkarea,NNZ,irn,icn,co,
     &                  CLOSUR,NMAXIT,ITERU,ERR,IERR,IUNIT,
     &                  AB,LENW,jvect,LENIW)
         ELSE IF(ISOLVC.EQ.3.or.isolvc.eq.1) THEN
            CALL DSLUCS(N,r,wkarea,NNZ,irn,icn,co,
     &                  CLOSUR,NMAXIT,ITERU,ERR,IERR,IUNIT,
     &                  AB,LENW,jvect,LENIW)
         ELSE IF(ISOLVC.EQ.4) THEN
            CALL DSLUGM(N,r,wkarea,NNZ,irn,icn,co,NVECTR,
     &                  CLOSUR,NMAXIT,ITERU,ERR,
     &                  IERR,IUNIT,AB,LENW,jvect,LENIW)
        ELSE IF(ISOLVC.EQ.5) THEN
            CALL DLUSTB(N,r,wkarea,NNZ,irn,icn,co,
     &                  CLOSUR,NMAXIT,ITERU,ERR,IERR,IUNIT,
     &                  AB,LENW,jvect,LENIW)
         END IF

         DO 450 I=1,N
            R(I) = wkarea(I)
  450    CONTINUE
C
         IF(MOP(6).NE.0) THEN
            CALL THYME(1,TSS,TT)
            WRITE (34,6040) TSS
         END IF
C
         iteruc = iteruc+iteru
C
      END IF
C
C
 6020 FORMAT(//,20('ERROR-'),//,T33,
     &             '       S I M U L A T I O N   A B O R T E D',
     &       /,T40,'DECLARED BANDWIDTH SMALLER THAN NEEDED',
     &       /,T33,'        PLEASE CORRECT AND TRY AGAIN',
     &       //,T20,'THE NUMBER OF NEEDED DIAGONALS IS', I4,
     &             ' WHILE THE AVAILABLE NUMBER IS ', I4,
     &       //,20('ERROR-'))
 6025 FORMAT(/,' The number of non-zero matrix elements is NNZ = ',
     &         I6/)
C
 6040 FORMAT('      SOLUTION TIME = ',1PE12.4,'  SECONDS')
 6045 FORMAT(T5,' At [',I4,',',I3,']',' DELT=',E12.6,
     &       ' IERR=',I1,'& ERR=',1PE12.6,' IT=',I5,' ITC=',I10)
C
C
C*********************************************************************
C*                                                                   *
C*                          UPDATE CHANGES                           *
C*                                                                   *
C*********************************************************************
C
  455 IF(MOP(6).GE.5) WRITE (34, 6050)
 6050 FORMAT(/,' ===== INCREMENTS ====   MASS BALANCES FIRST,',
     &         ' ENERGY BALANCE LAST',/)
C
         NEQC=1       !!!!!!!!!!!!!!!
         NK1C=1       !!!!!!!!!!!!!!!
      DO 500 NN=1,N
         NLOC  = (NN-1)*NEQC
         NLOCP = (NN-1)*NK1C
C
         IF(MOP(6).GE.5) WRITE (34,6055) ELEM(NN),(R(NLOC+K),K=1,NEQC)
C
         DO 480 K=1,NEQC
            DX(NLOCP+K) = DX(NLOCP+K)+WNR*R(NLOC+K)
  480    CONTINUE
  500 CONTINUE
C
 6055 FORMAT('       AT ELEMENT *',A5,'*   ',8(1X,E12.6))
C
      RETURN
      END
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!    
      SUBROUTINE Direct_Solution(
     &     nel,       ! Number of grid blocks
     &     ncon,      ! Number of connections
     &     Mnz,       ! The maximum length for zon-zero terms in the matrix
     &     Mneq,      ! Maximum Number of equations
     &     Mnel,      ! Maximum Number of grid blocks
     &     nz,        ! Number of zon-zero terms in the matrix
     &     irn,       ! Row indice of zon-zero terms in the matrix
     &     icn,       ! Column indice of zon-zero terms in the matrix
     &     co,        ! 1D array for zon-zero terms of the matrix
     &     r)         ! 1D array for right-hand side terms
! 
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*       Direct solution without calling solver for one or two         * 
!*              grid blocks or without connections                     *
!*                                                                     *
!*                   Version 1.0 - July 08, 2008                       *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
      IMPLICIT NONE                                      
!
! -------
! ... Integer variables
! -------
! 
      INTEGER*8  nel, Mnz, Mneq, Mnel, ncon, iel, nz, inz
!
! -------
! ... Integer arrays
! -------
! 
      INTEGER*8  irn(Mnz+1), icn(Mnz+1) 
! 
! -------
! ... Double precision variables
! -------
! 
      REAL*8     c11, c12, c21, c22, dumb, x11, x22
! 
! -------
! ... Double precision arrays
! -------
! 
      REAL*8     co(Mnz+1)
      REAL*8     r (Mneq*Mnel + 1)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>
!
!
      if (ncon == 0 .and. nel /= 2)   then
!
         do iel=1,nel
!
            if (co(iel) == 0.0d0)   then
               write (32,146) 
               stop
            end if
!
           r(iel) = r(iel)/co(iel)
! 
         end do
!        
      end if
!
!
      if (nel == 2)   then 
!
         c12 = 0.0d0
         c21 = 0.0d0
!
         do inz=1,nz
           if (irn(inz)==1 .and. icn(inz)==1)  c11 = co(inz) 
           if (irn(inz)==2 .and. icn(inz)==2)  c22 = co(inz)
           if (irn(inz)==1 .and. icn(inz)==2)  c12 = co(inz)
           if (irn(inz)==2 .and. icn(inz)==1)  c21 = co(inz)
         end do
!
         dumb = c11*c22 - c12*c21
!
         if (c11==0.0d0 .or. c22==0.0d0 .or. dumb==0.0d0)   then
            write (32,146) 
            stop
         end if
!           
         x11  = (c22*r(1)-c12*r(2))/dumb
         x22  = (r(2)-c21*x11)/c22
         r(1) = x11
         r(2) = x22
!
      end if
!
146   format (2x,'Error: divided by zero in solving', 
     &           ' solute transport equations')
!
      return
!
      end
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Get_TOUGH_HT_Variables
!
!
      IMPLICIT REAL*8(A-H,O-Z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
      include 'common_v2.inc'
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         Extract hydrological and thermal conditions from            *
!*         multiphase fluid flow (TOUGH) simulation results            *
!*                                                                     *
!*                    Version 1.0 - June 16, 2009                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
!
      COMMON/EOS_INDICATOR/ IEOS  ! indicate EOS module used
      COMMON/PARNP/NPL,NPG        ! phase index
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/CONTST/RE1,RE2,RERM,NER,KER,DFAC
      COMMON/KONIT/KON,DELT,IGOOD
      COMMON/SOLUTE8/SL1(MNEL)    ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)    ! new gas saturation
      COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)   ! water density (kg/dm**3)
      COMMON/E5/P(MNEL)
      COMMON/E6/T(MNEL)
      common/TEM_EOS9/Tc_EOS9(MNEL) 
      COMMON/P1/X((MNK+1)*MNEL)
      COMMON/P2/DX((MNK+1)*MNEL)
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
!
!.....For ECO2 and ECO2N flow modules
!
      common/co2_gene/nco2
      common/co2_gene1/ nco2g
      common/co2_gene2/ ico2gt0        !=1: initial Pco2>0
      COMMON/PCO2_ALL/PCO2A(MNEL)      ! calculated from ECO2
!
!.....For H2 generation by mineral phase using EOS5 module
!
      common/h2_gene1/ nh2g
      common/h2_gene2/ ih2gt0          !=1: initial Ph2>0
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>
!
!
! ----------
!.....Temperature
! ----------
!
      do n=1,nel
         tc(n) = t(n)
         tkelv(n) = tc(n) + 273.15d0
      end do
!
      IF (IEOS .EQ. 9)   THEN
         DO N=1,NEL
            T(N)=Tc_EOS9(N)
         END DO
      END IF
!
!-----------------------------------------------------
!
      if (ieos .le. 14) then
!
      do n=1,nel
!
! ----------
!........Water density
! ----------
!
         dwat_old(N)=dwat(N)   
         if (dwat_old(N).le.0.0d0.or.dwat_old(N).gt.2.0d0)   then
            dwat_old(N)=1.0d0
         end if
!
cc         Dwat(N)=1.0D0          ! kg/dm**3,   For using EOS9
!
         NLOC2=(N-1)*NSEC*NEQ1
         IF (IEOS .NE. 9)         THEN
            NLOC2L=NLOC2+NBK
            Dwat(N) = PAR(NLOC2L+4)
                                  ELSE  ! for eos9 module
            Dwat(N) = PAR(NLOC2+4)
         END IF
         Dwat(N)=Dwat(N)/1000.0D0       ! kg/dm**3
! 
! ----------
!........Phase saturations
! ----------
!
         np     = npl
         nl2np  = nloc2+(np-1)*nbk
         sl1(n) = par(nl2np+1)          ! Current liquid saturation
         sg1(n) = 1.0d0 - sl1(n)        ! Current gas    saturation
!
      end do
!
! ----------
!.....For using EOS2, ECO2, ECO2N modules
!.....Take CO2 partial pressure from these modules
! ----------
!
      IF (IEOS.EQ.2 .OR. IEOS.EQ.13 .OR.
     &    IEOS.EQ.5 .OR. IEOS.EQ.14)      THEN
         DO N=1,NEL
            IF (IEOS.EQ.2) THEN
               NLOC=(N-1)*NEQ
               PFUG2=X(NLOC+3)
               IF(KON.EQ.2) GOTO 3021
               PFUG2=PFUG2+DX(NLOC+3)
 3021          CONTINUE
               PFUG2=PFUG2/1.0d+5      ! co2 partial pressure (in bar)
            END IF
!
            NLOC2=(N-1)*NSEC*NEQ1      ! start of sec. variables for N
            NP=NPL
            NL2NP=NLOC2+(NP-1)*NBK
!
            IF (IEOS.EQ.13) THEN
               PFUG2=PCO2A(N)/1.0d+5   ! PFUG in bar, PCO2A from ECO2 module
            END IF
!
            IF (IEOS.EQ.14)      THEN  ! For ECO2N
               IF (SG1(N) .GT. sl1min)    then
!
! ----------
!   extracts dissolved CO2 concentration from ECO2N
!   If co2(aq) is not primary species, we need to add co2(aq) according to
!   to stoichiometry in database (e.g., co2 + h2o = hco3- + h+). We look
!   at the difference between the "react" and the ECO2N co2 concentration, and only
!   change by the difference.  Otherwise we can't distinguish co2(aq)contribution 
!   to total ut's for h+ and h2o.
! ----------
!
                  NLOC2L=NLOC2+NBK
                  X_CO2=PAR(NLOC2L+NB+3)             ! dissolved CO2 mass fraction
                  X_Salt=PAR(NLOC2L+NB+2)            ! dissolved NaCl mass fraction
                  X_H2O=1. - X_CO2 - X_Salt          ! water mass fraction                                  
                  co2_M = 1000.d0*X_CO2/44.d0/X_H2O * dwat(N)  ! aqueous total co2 molarity (mol/L) from ECO2N
                  co2_Mold=utold(n,nco2)             ! total aqueous co2 molality in "react"
!                 difference between "react" and ECO2N concentrations, expressed as molarity (mol/L)
                  diff = ( co2_M - co2_Mold )

                  ncp=ncpg(nco2g)                    ! number of components in CO2(g) stoichiometry
                  do icp=1,ncp                       ! loops over components in gas
                    stoic=stqg(nco2g,icp)            ! stoichio coeff
                    ind=icpg(nco2g,icp)              ! index of species in stoichio
                    UTold(N,ind)= UTold(N,ind) + stoic*diff 
                    ut(N,ind)=utold(n,ind)           ! not needed ?
                  end do
!
!.................Extract gas mole fraction for partial pressure calculations
!
                  NL2NP=NLOC2
!
                  fsH2O = PAR(NL2NP + nb + 1)         ! H2O Mass fraction in the gas phase
                  fsCO2 = PAR(NL2NP + nb + 3)         ! CO2 Mass fraction in the gas phase
                  Tfw   = fsH2O/18.016d0 + fsCO2/44.0d0 ! Sum (mass fraction/molar weight
                  fmol  = fsCO2/44.0d0/Tfw            ! Molar fraction
!
                  PFUG2 = fmol*P(N)/1.0d+5            ! partial pressure of CO2
                                          else
                  pfug2 = PFUG(N,NCO2G)
               END IF
            END IF
!
            IF (IEOS.EQ.5)   THEN      ! For EOS5, H2 generation
               IF (SG1(N) .GT. sl1min)    then
                  PFUG2 = P(N)/1.0d+5
                                          else
                  pfug2 = 1.0d-35
               END IF
            END IF
!
            IF (IEOS.EQ.2 .OR. IEOS.EQ.13 .OR. IEOS.EQ.14)      THEN
             IF (PFUG2 .GT. PFUG(N,NCO2G) .OR. ICO2GT0.NE.1) 
     &         PFUG(N,NCO2G)=PFUG2
            END IF
!
            IF (IEOS.EQ.5)      THEN
             IF (PFUG2 .GT. PFUG(N,NH2G) .OR. IH2GT0.NE.1) 
     &         PFUG(N,NH2G)=PFUG2
            END IF
!
         END DO
      END IF
!
!
      end if   ! ieos <= 14
!
!
      return
!
      end
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      subroutine sedim_1D
c
c     Routine to compute sedimentation by simple advection  of
c     solids+porewater vertically down with constant velocity.
c     There is no compaction, and the composition of the top block
c     is assumed constant (top boundary condition).  We loose at
c     the bottom what is advected down.
c
c     Note: this routine will work only for a 1D column, with 
c     gridblocks and connections ordered with depth (1st block is top of
c     sediment column).  Connections must be in order from top to bottom
c     with constant area throughout.
c
c     The sedimentation velocity, vsed, is positive downwards, in m/s
c
c     Simple forward explicit method
c
c     !!!! Make sure maximum dt is < minimum dx / vsed   !!!
c
c******************************************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'

      COMMON/C1/NEX1(MNCON)
      COMMON/C2/NEX2(MNCON)
      COMMON/C3/DEL1(MNCON)
      COMMON/C4/DEL2(MNCON)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT      
      COMMON/TRANGAS1/PFUGOLD(NMNOD,NMGAS)  ! old partial pressure
      COMMON/GASCONS1/GP0(NMNOD,NMGAS)      ! initial gas conc.
c     pre's are in mol_min/dm3 = mol_min/L solution
c     ut's are in mol/L solution
c     vsed is in m/s so we need dx in m 



      do n=1,ncon      ! loop through connections

        n1=nex1(n)     ! the first element of a connection
        n2=nex2(n)     ! the second element of a connection
        d1=del1(n)
        d2=del2(n)
        delta_x=d1+d2
        factor=vsed(n1)*deltex/delta_x    ! take upstream velocity
c
c------ Advects mineral amounts in discrete finite steps (forward explicit)
        do m=1,nmin
          pre(n2,m) = pre0(n2,m) + factor*(pre0(n1,m)-pre0(n2,m))
        enddo 
c
c------ Advects total concentrations
        do i=1,npri
          ut(n2,i) = utold(n2,i) + factor*(utold(n1,i)-utold(n2,i))
        enddo 
c
c------ Advects fugacities
        if(ngas.gt.0) then
          do m=1,ngas
            pfug(n2,m) = pfugold(n2,m) + factor*(pfugold(n1,m)-
     &                pfugold(n2,m))
            gp(n2,m) = gp0(n2,m) + factor*(gp0(n1,m)-gp0(n2,m))
          enddo 
        endif
c
      enddo

c     Need to resave new variables
      do n=1,nnod
           do m=1,nmin
              pre0(n,m)=pre(n,m)
           end do
           do i=1,npri
             utold(n,i) = ut(n,i)
           enddo 
           if(ngas.gt.0) then
             do m=1,ngas
                   gp0(n,m)=gp(n,m)
                   pfugold(n,m)=pfug(n,m)  !gas partial pressure in bars
             end do
           end if
c
      enddo
c
c     Need to recompute porosity in case of different velocities

      if(kcpl.ge.1) call phichg
c    
      return
      end
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE get_normalized(n,a,rnorm,pfugmin)
c
c *****************************************************************************
c     Subroutine Norm: This subroutine Nomalize an vector using the largest member
c ******************************************************************************
c
      implicit double precision (A-H,O-Z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
!
      real*8 a(mnod)
!
      rnorm   = a(1)
      pfugmin = a(1)
!
      do i=1,n
         if (a(i) .gt. rnorm)    rnorm   = a(i)
         if (a(i) .lt. pfugmin)  pfugmin = a(i)
      end do
!
      do i=1,n
         a(i) = a(i)/rnorm
      enddo

      return
      end
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE get_realized(n,a,rnorm,pfugmin)
c
c *****************************************************************************
c     Subroutine Norm: This subroutine realize a normalized vector
c *****************************************************************************
c
      implicit double precision (A-H,O-Z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
c
      real*8 a(mnod)
c
      do i=1,n
         a(i)=a(i)*rnorm
         if (a(i).lt.0.0d0) a(i)=pfugmin
      end do
      return
      end
c

