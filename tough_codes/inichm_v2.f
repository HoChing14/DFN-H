c
c
        subroutine init
c
c********************** Read and initialize chemical input data **************************
c
c  Chemical reaction variables used here.  These are defined
c  in common.inc and, for the most part, set by calling readtherm
c
c   npri  number of component species
c   naqx  number of derived aq. species
c   nmin  number of minerals
c   ngas  number of gases
c   nads  number of adsorption "species"
c   napri(i) names of primary species i=1,npri
c   naaqx(i) names of derived aqu. species i=1,naqx
c   namin(i) names of minerals i=1,nmin
c   nagas(i) names of gases i=1,ngas
c   naads(i) names of adsorption "species" i=1,nads
c   ncps(i),ncpm(m),ncpg(m),ncpad(m) number of comp. species in the
c      stoichiometry of derived aq. species, minerals, gases, and ads. species
c   icps(i,j),icpm(i,j),icpg(i,j),icpad(i,j) index of component species j in
c      reaction stoichiometry i of derived aq. species, minerals, gases, and ads. species
c   a0(i) D-H radii of primary and secondary aq. species i=1,npri+naqx
c   z(i) ionic charge of primary and secondary  aq. species   i=1,npri+naqx
c   akcoes(j,5) log10(K) f(T) regression coefficient for reaction j=1,naqx+nmin+ngas+nads
c     (5 coefficients, base=1)
c   akcoem(j,5), akcoeg(j,5), akcoead(j,5) same as above for derived aq. species,
c       minerals, gases, and ads. species
c   stqs(j,i) stoichiometric coefficient of component i in derived aq. species j
c   stqm(j,i) stoichiometric coefficient of component i in mineral j
c   stqg(j,i) stoichiometric coefficient of component i in gas  j
c   stqd(j,i) stoichiometric coefficient of component i in ads. species j
c   vmin(i)  molar volume of mineral i=1,nmin (in l/mole) (input in cm3/mole)
c   zd(i)  surface charge (?) for adsorption "species" i=1,nads
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
        include 'perm_v2.inc'
        character*200 inprec
c
        common/initial/cguess(mpri),icon(mpri)
        common/isarea/imflg2(mmin),imflag(mnel,mmin)
c
c...... Dissolution kinetics
        common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
        common/iprkin/ideprec(mmin)
        common/dispre/ idispre(mmin)
c... ...Add block for rate ph dependence parameters
        common/phdep/aH1(mmin),aH2(mmin),aH1p(mmin),aH2p(mmin),
     +  aHexp(mmin),aHexpp(mmin),aOHexp(mmin),aOHexpp(mmin)
c
c------For species dependent rate law: H+,  k(h+), expo(h+). term.
c------------------------------------ oh-, k(oh-), expo(oh-). term.
c------------------------------------k=k(h+)*[H+]**expo +....
      common/rksd0/ndep  ! number of minerals with kinetic rates from mechanisms other than neutral
      common/rksd1/  ids(mmin,mechm,mechsp) ! dependent species pointer in CT(naqt)
c                       (mmin, No. mechanisms, No.species involved)
      common/rksd2/ rkds(mmin,mechm)   ! constant for the dependent species
      common/rksd21/ eads(mmin,mechm)  ! activation energy for the dependent species
      common/rksd3/expds(mmin,mechm,mechsp) ! exponential term for the dependent species
      common/rksd4/ndis(mmin)      ! number of additional mechanisms
      common/rksd41/ nspds(mmin,mechsp) ! number of speciess involved in one mechanism
      common/rksd1p/  idsp(mmin,mechm,mechsp)
      common/rksd2p/ rkdsp(mmin,mechm)
      common/rksd21p/ eadsp(mmin,mechm)
      common/rksd3p/expdsp(mmin,mechm,mechsp)
      common/rksd4p/npre(mmin)
      common/rksd41p/ nsppr(mmin,mechsp)
      common/deriv_mech/der_ds(mmin,mechm,mpri),der_pr(mmin,mechm,mpri)
      character*20 nadis(mmin,5,5),napre(mmin,5,5) ! name of species involved in rate constant
c
c-------Added common block for rate law designations
        common/irtlaw/nplaw(mmin)
        COMMON/SOLUTE7/SGOLD(MNEL)        ! old gas saturation
c-------Added 2 common blocks below to use por(matx(i))
        COMMON/SOLID/NM,DROK(maxmat),POR(maxmat),PER(3,maxmat),
     +                CWET(MAXMAT),SH(MAXMAT)
        COMMON/E2/MATX(MNEL)
c-------Solid solutions
        common/solsol/iss(mmin),ncpss(msol),icpss(msol,mcpss),nss
C----------For mineral precipitation at dry grid blocks
c       common/dry_salt/nsalt,isalt(0:mmin),rkamp(mmin)
        common/dry_salt/nsalt,isalt(mmin)
        integer*8 ksalt(mmin)
c
c------------------------------------------------------------------------
        common/aqxs/iaqxs ! =1: define secondary species in CHEMICAL.INP
c-------------------------------COMMON blocks for Kd adsorption and decay
c
        common/kddca1/nakdd(mpri)   ! name of species with Kd and decay
        common/kddca2/decayc(mpri)  ! decay constants
        common/kddca21/a_TDecay(mpri),  ! Thermal decay parameter, a
     &                 b_TDecay(mpri)   ! Thermal decay parameter, b
        common/kddca3/kddp(mpri)    ! pointer to the primary species
        common/Kddca4/vkd(30,mpri)  ! values of Kd in initial zones
        common/Kddca5/izonekd(mnod) ! Kd zone code
        common/Kddca6/sden(30,mpri) ! solid density
        common/Kddca7/nkdd          ! number of species with Kd adsorption
        common/kddca8/kdflag(30,mpri)
        CHARACTER*20  nakdd
c
c-------------------For co2 generation by mineral phase using EOS2 module
        common/co2_gene/ nco2
        common/co2_gene1/ nco2g
        common/co2_gene2/ ico2gt0     !=1: initial Pco2>0
c
c-------------------For H2 generation by mineral phase using EOS5 module
        common/h2_gene/ nh2
        common/h2_gene1/ nh2g
        common/h2_gene2/ ih2gt0     !=1: initial Ph2>0
!
!.......For gas reaction sources for TMVOC modules
        common/gas_gene1/ nch4, nh2s, nso2
        common/OtherGases_ix/ nch4g, nh2sg, nso2g
!
!.......Use aqueous complex to model surface complexation
        common/ns_index/ ns
!
!.......Surface complexation
        common/ion_str1/str              ! ionic strength
!
!.................................................................
        common/ion_str2/str_node(mnel)   ! ionic strength for all nodes
C
        COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
        COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
C
C-------------------------------------------------------------------------
C
        common/water_activity/aw(mnel)    !water activity del gxzh
        COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
        COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
        common/nwater/nibw,iwtype,niwtype  ! for chdump only
        COMMON/MOP_REACT/MOPR(20)  ! control parameters for reactive transport
        common/eqinit/sateq(mpri),mineq(mpri),igaseq(mpri)   !NS3/06
        character*20 nameq         !NS3/06
        character*100 label

        character*20 nadum2,nadum1
        integer*8 idum1(mexc),idum2(mexc)
        double precision vol2(mmin)
        integer*8 idxss(msol)
        integer*8 iterch
        integer*8 inode
c
c.......For aqueous kinetics
        common/aqkin1/nrx       ! number of redox pair
        common/aqkin2/ntrx      ! total number of redox pair
        common/aqkin16/NoTrans(mpri)    ! >0: not subject to transport
c
c.......For ion exchange under unsaturated conditions
c
        common/satgas2/sg2
        COMMON/SOLUTE9/SG1(MNEL)            ! new gas saturation
        common/rock_density/density_rock    ! rock density of current node
        COMMON/tot_solid_aq/icon_nod(mnod,mpri),ttt(mpri),
     &                ttt_nod(mnod,mpri)    ! total concentraion including both aqueous and solid
c            Store ikin4 for no reaction option 
        integer*8 ikin4m(mmin)
!
!.....Extract rock density for geochemical calculations such as exchange and sorption
      common/rock_density1/denss(mnel)
      common/rock_density2/denss2, sl2, phisl2, a_fmr2
!
!....................................................................................
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' init 1.0, 2008.2.6: Read and initialize chemical input'
     x' data')
c
c.......initialization of indexes
c
        nw    = 0       ! water species index
        nh    = 0       ! hydrogen species index
        noh   = 0       ! oh- species index
        ne    = 0       ! electron species index
        no2aq = 0       ! o2 aq. index
        nco2  = 0       ! aqueous carbon primary species index
        nco2g = 0       ! co2 gas index
        nch4g = 0       ! ch4 gas index
        nh2sg = 0       ! h2s gas index
        nso2g = 0       ! so2 gas index
        nh2   = 0       ! aqueous H2 species index among all species
        nh2g  = 0       ! H2 gas index
        nd    = 0       ! Primary surface species index
        ns    = 0       ! Use aqueous complex to model surface complexation
c
c.......Addition for TMVOCs 
c
        nch4  = 0       ! CH4 related primary species index
        nh2s  = 0       ! H2S related primary species index
        nso2  = 0       ! SO2 related primary species index
c
c.... Chemical input file (move from readsolu)
c
      OPEN (UNIT=41,FILE='chemical.inp',STATUS='OLD')
      OPEN (UNIT=42,FILE='chemical.out',STATUS='UNKNOWN')
c
c     run log file - we created this file in readsolu
c     reopen it under another unit no. and append to it
c
      open(unit=32,file='runlog.out',status='unknown',position='append')
c
c.... Starts reading the CHEMICAL.INP file (unit 41), and echo the data in the
c     CHEMICAL.OUT file (unit 42)
c
      write(*,*) '   --> reading chemical input data'
      write(32,"(/
     &  ' --> start reading chemical input data from chemical.inp'/)")
      write(42,"(/
     &  ' --> start reading chemical input data from chemical.inp'/)")
c  
c--Read title for the problem (3 lines)
c
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      end do
      read(inprec,'(a)',err=9001) label
      write (42,"(a,/,80('-'))") label
c
c  Read primary (component) species in the system (aqueous and surface species)
c  ----------------------------------------------------
c  List of species must be ordered with aqueous species first, and surface
c  primary species next. The surface primary species must be neutral!
c
c  read header label (needed)
c
      write(42,"(' DEFINITION OF THE GEOCHEMICAL SYSTEM')")
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      i=0        !increment index of all primary species, 1 to npri
      npads=0    !index up to total number of surface primary species
      nsurf=0    !index up to total number of minerals (surfaces) tied to surface complexes
      npot=0     !index of surfaces with non-zero potential terms
c
c.... Start of loop through primary species
130   continue
c
       if(i.ne.0) read(41,"(a200)",err=9002) inprec
       if(index(inprec(1:20),"*").ne.0) goto 120     ! end of primary species
c
       i=i+1
       if(i.eq.1) write(42,"(' PRIMARY AQUEOUS SPECIES')")
       IF(i.gt.mpri) THEN
          WRITE(32,35) mpri
35        FORMAT(/' Error: maximum number of component species (MPRI)='
     1           ,i3,'was exceeded.')
          STOP
       END IF
c
       read(inprec,*,end=9002, err=9002) nadum1, iflag      ! iflag=1: not subject to transport
c                                                           !      =2: surface complexes
       if(iflag.eq.2) then ! primary species is a surface complex, get the name of the associated surface mineral, site density and ads model type
c
          read(inprec,*,end=9002, err=9002) napri(i), NoTrans(i),
     &                                  nadum1, sdens, imod
c
c       sdens-->site_dens    site density in mol_sites/m2)
c         if nadum1 starts with 'surface', the surface area will be fixed and input in m2/kgw later
c         if nadum1 is 'no_mineral' a default mineral density will be used for surface area calcs
c       imod-->iadmod()   Model type: 0 no eletrostatic terms
c                                     1 constant capacitance
c                                     2 double diffuse layer, linear
c                                     3 double diffuse layer, Gouy-Chapman
c         Species with same model types MUST have same surface names and capacitance, but
c          not necessarily same site density
c
          label(1:20) = nadum1
          CALL name_conv(label)
          nadum1 = label(1:20)

          label(1:20) = napri(i)
          call name_conv(label)
          napri(i) = label(1:20)
c
          capac=0.d0
c         If model type 1, we read the same string again to get the
c           capacitance in F/m2 (or C/V/m2)
          if(imod.eq.1) then
             read(inprec,*,end=9002, err=9002) napri(i),
     &        NoTrans(i), nadum1, sdens, imod, capac
            if(capac.eq.0.d0) then
              write(42,"(5x,/'Adsorption species: ',A12,' Selected ',
     &        'sorption model 1 requires non-zero capacitance' )")
     &         napri(i)
              write(32,"(5x,/'Adsorption species: ',A12,' Selected ',
     &        'sorption model 1 requires non-zero capacitance' )")
     &         napri(i)
              stop
            endif
          endif
c
          npads=npads+1
          if(npads.gt.mpads) then
            write(32,"(5x,/'Number of primary surface species exceeds'
     &       ' the maximum allowed (',i2,')/')") mpads
            write(42,"(5x,/'Number of primary surface species exceeds'
     &       ' the maximum allowed (',i2,')/')") mpads
            stop
          endif
c
          site_dens(npads)=sdens    !site density
c
          if(npads.eq.1) then
            nsurf=1
          else
            ifound=0
            do n=1,npads
              if(nadum1.eq.naads_min(n)) ifound=1
            enddo
            if(ifound.eq.0) nsurf=nsurf+1
          endif

          if(nsurf.gt.msurf) then
            write(32,"(5x,/'Number of different surfaces exceeds the'
     &       ' maximum allowed (',i2,')/')") msurf
            write(42,"(5x,/'Number of different surfaces exceeds the'
     &       ' maximum allowed (',i2,')'/)") msurf
            stop
          endif
          if(nadum1.eq.'                    ') nadum1='no_mineral'
          naads_min(nsurf) = nadum1
          phip2(nsurf)=0.d0            !initialize potential term for that surface
          isurfp(npads) = nsurf        !points to index of surface, for each primary surface complex
          iadmod(nsurf)=imod           !adsorption model for each surface (can change during run-time)
          capacitance(nsurf)=capac     !fixed capacitance if any
c
       else
          read(inprec,*,end=9002, err=9002) napri(i), NoTrans(i)   ! NoTrans=1: not subject to transport
          label(1:20) = napri(i)
          call name_conv(label)
          napri(i) = label(1:20)

       endif
c
       write (42,*) napri(i), NoTrans(i)  ! write primary species to file chemical.out
       if(notrans(i).eq.2) then
         label(1:30) = 'Site density (mol_sites/m2)='
         if(iadmod(nsurf).ne.1) then
            write(42,"(5x,'Surface species linked to: ',A20,/5x,A30,
     &      e10.4,' Adsorption model=',i2)")
     &      naads_min(nsurf), label(1:30),site_dens(nsurf),iadmod(nsurf)
         else
            write(42,
     &      "(5x,'Surface species linked to: ',A20,/5x,A30,e10.4,
     &      ' Adsorption model=',i2,' Constant capacitance (F/m2)=',
     &      e10.4)") naads_min(nsurf), label(1:30),site_dens(nsurf),
     &      iadmod(nsurf), capacitance(nsurf)
         endif
       endif
c
      goto 130
c.... end of loop through primary species
120   continue
      npri=i  ! number of primary species
      npaq=npri-npads ! no absorb primary species
      if(npri.gt.0) write(42,918) npri,npaq
918   format(/' finished reading all primary species'/' total'
     x' primary species = ',i3,', aqueous = ',i3)
      if(npads.ne.0) write(42,917) npads
917   format(' surface primary species = ',i3,/)
      ntrx=0  ! total number of redox pair
c
c Read (optional) intra-aqueous kinetics and biodegradation
c ----------------------------------------------------------
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      if(index(inprec(1:20),"*").ne.0) goto 131     ! end of secondary species
c
         read(inprec,*,err=9022) ntrx
         if(ntrx.gt.mrx) go to 9903
         if(ntrx.ne.0) write(42,919) ntrx
919      format('AQUEOUS KINETICS'/,i3)
         do nrx=1,ntrx
            CALL READ_AQKIN
         end do
         read(41,'(a)',err=9022) label   !reads the star if kinetic data were provided
c
131   continue
      if(ntrx.ne.0) write(42,"(' finished reading aqueous kinetics,'
     x' ntrx = ',i3,/)") ntrx
c
c Reads (optional) secondary aqueous species
c ----------------------------------------------------------
c     
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      j=0
c--- Start of loop through specified secondary species
132   continue
cns10/09       read(41,"(a200)",err=9202) inprec
      if(j.ne.0) read(41,"(a200)",err=9202) inprec
      if(index(inprec(1:20),"*").ne.0) goto 122     !end of secondary species
c
      j=j+1
      if(j.eq.1) write(42,"(' AQUEOUS COMPLEXES')")
      IF(MAQX .LT. (J-1)) THEN
         write(42,45) MAQX
         WRITE(32,45) MAQX
45       FORMAT(/' Error: maximum number of secondary species (MAQX)'
     x   ' was exceeded. Current max = ',i3)
         STOP
      END IF
c
      read(inprec,*,end=9202,err=9202) naaqx(j)
c
      label(1:20) = naaqx(j)
      CALL name_conv(label)
      naaqx(j) = label(1:20)
c
      write (42,*) naaqx(j)
c
      goto 132
c--- end of loop through specified secondary species
122   continue
      naqx=j
      if(naqx.ne.0) write(42,920) naqx
920   format(' finished reading secondary species, naqx = ',i3)
c
c     when iaqxs=0, automatically picks up secondary species with given components
c     when iaqxs=1, define secondary species below
      iaqxs=1
      if(naqx.eq.0) then
         write(42,"(' AQUEOUS COMPLEXES will be automatically'
     &     ' selected from thermodynamic database',/)")
      iaqxs=0
      endif
c
c     reads minerals
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
c  read minerals in the system
c  ----------------------------
c  Note: minerals need a header label, it was read earlier
133   continue
c
cels10/12/09 initialize some arrays
      do m = 1, mmin
        isalt(m) = 0
        ksalt(m) = 0
      enddo
c
      m=0
      nsalt=0    ! number of possible mineral phases at dry block (TX, 05/15/01)
      nmkin=0             !Addition for kinetics
      ndep=0     ! number of minerals with species dependent dis/pre rate
c
c--- start of loop through minerals
140   continue
cns10/09       read(41,"(a200)",err=9003) inprec
      if(m.ne.0) read(41,"(a200)",err=9003) inprec
      if(index(inprec(1:20),"*").ne.0) goto 150     !end of minerals
c
      m=m+1
      if(m.eq.1) write(42,"('MINERALS')")
      IF(m.gt.mmin) THEN
         WRITE(32,55) mmin
55       FORMAT(/' Error: maximum number of minerals (mmin) was'
     1           ' exceeded.  Current max = ',i3)
         STOP
      END IF
c
      read(inprec,*,end=9003,err=9003) namin(m),ikin(m),
     &            idispre(m), iss(m), ksalt(m)
      label(1:20)=namin(m)                     
      call name_conv(label)
      namin(m)=label(1:20)                    
c
c----------Flag kinetic mineral for kicking out when dryingout
c
c Removing Pitzer lines
c       if(mopr(17).eq.0) then
c         if (m1.gt.nsalt) nsalt=m1
c         isalt(m1)=m                             ! for mineral precipitation at dry block
c       else
c        nsalt=nmin
c        if(m1.lt.0) isalt(m)=-1    !kicking out
c        if(m1.eq.0) isalt(m)=0     !nothing to do (no change)
c        if(m1.ge.1) isalt(m)=m1    !amplifying rate constant for 99, and the exponent for >=100
c       endif
c
c       Rewrite so arrays do not start at zero 
        if (ksalt(m).gt.0)then
           nsalt = nsalt + 1
        endif
c
c   ikin is flag for kinetics (=0 for no kinetics, =1 for kinetics)
c   idispre is flag = 1 for dissol only, 2 for precipit, 3 for both
c   iss is solid solution index (endmembers must have the same index)
c
       write (42,*) namin(m),ikin(m),idispre(m)
       if (ikin(m).ge.1)  nmkin=nmkin+1
c
c supersaturation constant for equil. minerals
c ssq0(m) is full supersaturation "window" in +log(K) units
c sst1(m) is temperature (C) above which that window starts to
c         decrease exponentially with temperature (typically 25 C)
c sst2(m) is temperature (C) at which ssq0 becomes one hundreth of
c         starting value (essentially no window anymore)
       if (ikin(m).eq.0) then
        read(41,*,err=9003) ssq0(m), sst1(m), sst2(m)
        if (ssq0(m).lt.0.d0) ssq0(m)=0.d0
        if (sst1(m).lt.0.d0) sst1(m)=0.d0
        if (sst2(m).le.0.d0) then
           sst2(m)=0.d0
        elseif (sst2(m).le.sst1(m)) then
           write(42,*) namin(m),
     +     ' ** error - sst2 is non zero and less or equal to sst1 '
           stop
        endif
        ssq(m)=ssq0(m)  !ssq is reset with temp. in subroutine assign
        write (42,"(3f10.3)") ssq0(m), sst1(m), sst2(m)
       else if (ikin(m).ge.1) then
c---------- initialize ck's and ckprec's to 1 to
c to avoid bombs later in cr_cp if undefined.
       ck1(nmkin)=1.d0
       ck2(nmkin)=1.d0
       ck1prec(nmkin)=1.d0
       ck2prec(nmkin)=1.d0
c--------- added supersaturation constant
       ssqk(nmkin) = 0.d0
c
c... Read dissolution rate constants, activ energy, exponents, coeff
       if (idispre(m).eq.2) go to 1089
       read(41,*,err=9003)rkf(nmkin),idep(nmkin),ck1(nmkin),ck2(nmkin),
     +    ea(nmkin),acfdiss(nmkin),bcfdiss(nmkin),ccfdiss(nmkin)
c----------- for dissolution set nucleation radius=0
       rnucl(nmkin) = 0.d0
c----------- reads ph dependence parameters if idep<>0
c aH1 is activity of H+ near which pH influence starts
c aHexp is aH+ exponent (slope of log(rate) with pH)
c aH2 is activity of H+ near which pOH influence starts
c aOHexp is aOH exponent (slope of log(rate) with pOH)
       if(idep(nmkin).eq.1) then
         read(41,*,err=9003) ph1tr,aHexp(nmkin),ph2tr,aOHexp(nmkin)
         aH1(nmkin)=10.d0**(-ph1tr)
         aH2(nmkin)=10.d0**(-ph2tr)
       end if
c
c------For species dependent rate law: H+,  k(h+), expo(h+). term.
c------------------------------------ oh-, k(oh-), expo(oh-). term.
c-------------------------------------nadum2, rkds,  expds
c------------------------------------k=k(h+)*[H+]**expo +....
c      common/rksd1/  ids(mmin,mechm,mechsp) ! dependent species pointer in CT(naqt)
c      common/rksd2/ rkds(mmin,mechm) ! constant for the dependent species
c      common/rksd3/expds(mmin,mechm,mechsp) ! exponatial term for the dependent species
c                   nadis(mmin) ! name of involved species
       if(idep(nmkin) .ge. 2)  then
        ndep=ndep+1 ! number of minerals with species dependent dis/pre rate
          read(41,*) ndis(nmkin) ! number of mechanisms
          if(ndis(nmkin).gt.mechm) go to 9901
          do ii=1,ndis(nmkin)
c            modified to ensure array dimensions are not exceeded
             read(41,"(a200)",err=9003) inprec
             read (inprec,*,err=9003) dum,dum,ndum   
c            reset dum so  it doesn't interfere with other uses of dum
             dum = 0.d0
             if(ndum.gt.mechsp) go to 9902
             read (inprec,*,err=9003)  rkds(nmkin,ii),
     &                             eads(nmkin,ii),nspds(nmkin,ii),
     &                  (nadis(nmkin,ii,isp),expds(nmkin,ii,isp),
     &                    isp=1,nspds(nmkin,ii))
!
             do isp=1,nspds(nmkin,ii)
                label(1:20) = nadis(nmkin,ii,isp)                     
                call name_conv(label)
                nadis(nmkin,ii,isp) = label(1:20)                    
             end do
!
          end do
       end if
c
c... If precipitation kinetics are allowed, read in rate law data
c
1089   continue
c
       if(idispre(m).eq.3.or.idispre(m).eq.2)then
        read(41,*,err=9003)rkprec(nmkin),ideprec(nmkin),
     +    ck1prec(nmkin),ck2prec(nmkin),eaprec(nmkin),acfprec(nmkin),
     +    bcfprec(nmkin),ccfprec(nmkin),rnucl(nmkin),nplaw(nmkin)
c----------- reads ph dependence parameters if ideprec<>0
c note: see dissolution parameters explanations above
        if(ideprec(nmkin).eq.1) then
         read(41,*,err=9003) ph1trp,aHexpp(nmkin),ph2trp,aOHexpp(nmkin)
         aH1p(nmkin)=10.d0**(-ph1trp)
         aH2p(nmkin)=10.d0**(-ph2trp)
        endif
c
        if (ideprec(nmkin) .ge. 2)  then
         ndep=ndep+1 ! number of minerals with species dependent dis/pre rate
          read(41,*) npre(nmkin) ! number of species involved in the rate constant
          do ii=1,npre(nmkin)
             read (41,*,err=9003) rkdsp(nmkin,ii),
     &                            eadsp(nmkin,ii),nsppr(nmkin,ii),
     &         (napre(nmkin,ii,isp),expdsp(nmkin,ii,isp),
     &                          isp=1,nsppr(nmkin,ii))
!
             do isp=1,nsppr(nmkin,ii)
                label(1:20) = napre(nmkin,ii,isp)                     
                call name_conv(label)
                napre(nmkin,ii,isp) = label(1:20)                    
             end do
!
          end do
        end if
c
c
c-- move this after endif if(rnucl(nmkin).eq.0.d0) rnucl(nmkin) = 1.e-5  ! volume fraction
c----------- add separate input line for supersaturation window
c similar parameters as defined for equilibrium minerals above
        read(41,*,err=9003) ssqk0(nmkin), sstk1(nmkin), sstk2(nmkin)
        if (ssqk0(nmkin).lt.0.d0) ssqk0(nmkin)=0.d0
        if (sstk1(nmkin).lt.0.d0) sstk1(nmkin)=0.d0
        if (sstk2(nmkin).le.0.d0) then
           sstk2(nmkin)=0.d0
        else if (sstk2(nmkin).le.sstk1(nmkin)) then
           write(42,*) namin(m),
     +     ' ** error - sstk2 is non zero and less or equal to sstk1 '
           stop
        endif
        ssqk(nmkin)=ssqk0(nmkin)  !ssq is reset with temp. in subroutine assign
       endif
c
c   rkf: coefficient a in expression k = a*exp(ea/RT)
c      rate=k*surf_area*(1-(q/k)**m)**n * act_H+
c   idep: flag=1 for pH dependence on rate, =0 for no pH dependence
c   ck1: exponent (n) in rate expression above
c   ck2: Q/K exponent (m) in rate expression above
c   activation energy
c
       if (idispre(m).eq.2) go to 1090
         write(42,'(e10.3,i5,6f10.3)') rkf(nmkin),idep(nmkin),
     +                       ck1(nmkin),ck2(nmkin),ea(nmkin),
     +           acfdiss(nmkin),bcfdiss(nmkin),ccfdiss(nmkin)
         if(idep(nmkin).eq.1) write(42,'(4f10.4)') ph1tr,aHexp(nmkin),
     &            ph2tr,aOHexp(nmkin)
         if(idep(nmkin).ge.2)  then
            write(42,'(i10)') ndis(nmkin)
            do ii=1,ndis(nmkin)
               write (42,'(9x,e14.4,f10.4,i5,4x,5(2x,a15,f10.4))')
     +           rkds(nmkin,ii),eads(nmkin,ii),nspds(nmkin,ii),
     +           (nadis(nmkin,ii,isp),expds(nmkin,ii,isp),
     +                            isp=1,nspds(nmkin,ii))
            end do
         end if
c
1090   continue
       if(idispre(m).eq.3.or.idispre(m).eq.2)then
         write(42,'(e10.3,i5,6f10.3,e10.4,i3)')
     +           rkprec(nmkin),ideprec(nmkin),
     +           ck1prec(nmkin),ck2prec(nmkin),eaprec(nmkin),
     +           acfprec(nmkin),bcfprec(nmkin),ccfprec(nmkin),
     +           rnucl(nmkin),nplaw(nmkin)
        if(ideprec(nmkin).eq.1) write(42,"(4f10.4)") ph1trp,
     &       aHexpp(nmkin),ph2trp,aOHexpp(nmkin)
        if(ideprec(nmkin).ge.2)  then
             write (42,'(i10)') npre(nmkin)
             do ii=1,npre(nmkin)
                write (42,'(9x,e14.4,f10.4,i5,4x,5(2x,a15,f10.4))')
     +             rkdsp(nmkin,ii),eadsp(nmkin,ii),nsppr(nmkin,ii),
     +             (napre(nmkin,ii,isp),expdsp(nmkin,ii,isp),
     +                              isp=1,nsppr(nmkin,ii))
             end do
        end if
        write (42,'(3f10.3)') ssqk0(nmkin), sstk1(nmkin), sstk2(nmkin)
        endif
       end if
c
      goto 140
c--- end of loop through minerals
150   continue
      nmin=m
      nmequ=nmin-nmkin
c
c     added number of salt minerals nsalt below
      if(nmin.ne.0) write(42,925) nmin
925   format(' finished reading minerals'/' nmin = ',i3)
      if(nmequ.ne.0) write(42,"(' nmequ = ',i3)") nmequ
      if(nmkin.ne.0) write(42,"(' nmkin = ',i3)") nmkin
      if(nsalt.ne.0) write(42,"(' nsalt = ',i3,/)") nsalt
c
c--- pointer for minerals specified under equilibrium
c    for precipitation and kinetic for dissolution
      do m=1,nmkin
         kineq(m)=0
c---- nm is already a tough2 variable     do nm=1,nmequ
         do nim=1,nmequ
           if(namin(nim).eq.namin(m+nmequ)) kineq(m)=nim
         enddo
      enddo
c
c       reassignment and reordering of salt minerals
        mslt = 1
        iflgsalt=0
      do n = 1, nsalt
        do m = 1, nmin
          if (ksalt(m).gt.nsalt)iflgsalt=1
          if (ksalt(m).eq.mslt)then
             isalt(mslt) = m
             mslt = mslt + 1
          endif
        enddo
      enddo
      If(iflgsalt.eq.1)write(*,*) 'Fix salt mineral order'
c
c   Precipitation rate amplifying factor due to rate constant, This factor amplifies the rate to 10e-07mol/m2/s
c
c
cpitz_dry       do m=1,nmkin                                         
cpitz_dry          rkamp(m)=0.0d0
cpitz_dry          if (rkprec(m).gt.0.0d0) then
cpitz_dry               rkampm=0.0d0-dlog10(rkprec(m))
cpitz_dry               if(rkampm.lt.0.0d0) rkampm=0.0d0
cpitz_dry               rkamp(m)=10.0d0**(rkampm-7.0d0)
cpitz_dry          endif
cpitz_dry       enddo
c
c----------- assigns various flags/pointers for solid solution
c nss:  total number of solid solutions
c iss(m): solid solution index for mineral m
c ncpss(n): number of endmembers in solid solution index n (n=iss(m))
c icpss(n,k): mineral index of endmember k, in solid solution index n
c idxss(id): index n of of solid solution ID (to allow input flexibility)
c
c     initialize first
      nss = 0
      do m=1,msol
        ncpss(m)=0
        idxss(m)=0
        do n=1,mcpss
         icpss(m,n)=0
        enddo
      enddo
      if(nmin.ne.0) write(42,926)
926   format(' The following minerals form ideal solid solutions:')
c
      do m=1,nmin
       isol = iss(m)
       if(isol.gt.msol) then
         write(*,*)
     +    ' Solid solution ID exceeds maximum number of solid solutions'
           write(42,*)
     +    ' Solid solution ID exceeds maximum number of solid solutions'
           write(32,*)
     +    ' Solid solution ID exceeds maximum number of solid solutions'
           stop
       endif

       if(isol.ne.0) then
         if(idxss(isol).eq.0)  then    ! increment only if sol sol not already found once
             nss=nss+1                 ! solid solution order index (1,2,3,4 etc)
             idxss(isol) = nss
         endif
         if(nss.gt.msol) then
           write(*,*) ' Maximum number of solid solutions exceeded'
           write(42,*) ' Maximum number of solid solutions exceeded'
           write(32,*) ' Maximum number of solid solutions exceeded'
           stop
         endif
c
         ncp=1
         do m2=1,nmin
          if(iss(m2).eq.isol) then
            if(ncp.gt.mcpss) then
             write(42,"(/5x,'Maximum number of ss endmembers exceeded'
     &       /5x,'Mineral: ',a8,' solid solution index: ',i3)")
     &       namin(m2), iss(m2)
             write(*,*) ' Maximum number of ss endmembers exceeded'
             write(32,*) ' Maximum number of ss endmembers exceeded'
             stop
            endif
            if(icpss(nss,ncp).eq.0) then
             icpss(nss,ncp)=m2
             ncpss(nss)=ncp
             ncp=ncp+1
             write(42,"(5x,'Mineral: ',a12,5x,'Solid solution ID:',
     &          i3,'   Index:', i3)") namin(m2), iss(m2), nss
            endif
          endif
         enddo
         write(42,*)
       endif
      enddo
      if(nss.ne.0) then
         write(42,*) ' Total number of solid solutions:', nss
         write(42,*)
      endif
c
c  read gases in the system
c  ------------------------

c     read header label (needed for gases)
c      read(41,'(a)',err=9004) label             !gases
c      write (42,'(a)') label   
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      k=0
c--- start of loop through gases
170   continue
c
      if(k.ne.0) read(41,"(a200)",err=9004) inprec
      if(index(inprec(1:20),"*").ne.0) goto 190     !end of minerals
c
      k=k+1
      if(k.eq.1) write(42,"('GASES')")
      IF(k.gt.mgas) THEN
          WRITE(32,65) mgas
65        FORMAT(/' Error: maximum number of gases (mgas) was exceeded.'
     1           ' Current max = ',i3)
          STOP
      END IF
c
      read(inprec,*,end=9004,err=9004) nagas(k)
c
      label(1:20)=nagas(k)                     ! name regulation pitz zh 12/16/04
      call name_conv(label)
      nagas(k)=label(1:20)                     ! name regulation pitz zh 12/16/04
c
      write(42,*) nagas(k)
c
      goto 170
c--- end of loop through gases
190   continue
      ngas=k
      if(ngas.ne.0) write(42,930) ngas
930   format(' finished reading gases, ngas = ',i3/)

c  read (optional) surface complexes in the system
c  -------------------------------------
c      read(41,'(a)',err=9005) label               !surface complexes, assciated mineral, neutral complex
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      j=0
c--- start of loop through surface complexes
230   continue
      if(j.ne.0) read(41,"(a200)",err=9005) inprec
      if(index(inprec(1:20),"*").ne.0) goto 220     !end of surface complexes
c
      j=j+1
      if(j.eq.1) write(42,"(' SURFACE COMPLEXES')")
      IF(j.gt.mads) THEN
         WRITE(32,75) mads
75       FORMAT(/' Error: maximum number of ads. species (mads) was'
     1   ' exceeded. Current max = ',i3)
         STOP
      END IF
c
      read(inprec,*,end=9005,err=9005) naads(j)
c
      label(1:20) = naads(j)
      CALL name_conv(label)
      naads(j) = label(1:20)
c
      write (42,*) naads(j)
      goto 230
c--- end of loop through surface complexes
220   continue
      nads = j
      if(nads.ne.0) write(42,932) nads
932   format(' finished reading surface complexes, nads = ',i3)
      if(nads.eq.0.and.npads.gt.0) then
         write(42,"(' Surface complexes will be automatically'
     &     ' selected from the thermodynamic database',/)")
      end if
c
c
c  read species (primary) with linear equilibrium Kd and decay
c  -------------------------------------------------------------
c       nakdd ( ) : species name
c       decayc ( ): decay constant
c       kddp ( )  : pointer for species with Kd and decay
c
c     read header label (needed)
c    read(41,'(a)',err=9026) label
c      write (42,'(a)') label
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      k=0
c
c--- start of loop through kd/decay species
172   continue
c
       if(k.ne.0) read(41,"(a200)",err=9026) inprec
       if(index(inprec(1:20),"*").ne.0) goto 192     !end of kd/decay species

       k=k+1
       if(k.eq.1) write(42,1077)
1077   format('SPECIES WITH KD AND DECAY'/'Species',10x,'Decay const.'
     &' (1/s)    Th. decay a Th. decay b')
       if(k.gt.mpri) then
          write (42,"(/1x,'Error: maximum number of species with',
     &      ' kd/decay was exceeded.',/1x,'Current max= ',i3)") mpri
          stop
       end if
c
       read  (inprec,*,end=9026,err=9026) nakdd(k),   ! Name of primary species with kd or/and decay
     &                    decayc(k),    ! Decay constant (1/s)
     &                    a_TDecay(k),  ! Thermal decay parameter, a
     &                    b_TDecay(k)   ! Thermal decay parameter, b
c
       label(1:20) = nakdd(k)
       CALL name_conv(label)
       nakdd(k) = label(1:20)
c
       write (42,"(a20,3(e12.5,5x))")  nakdd(k), decayc(k), 
     &     a_TDecay(k), b_TDecay(k)
c
       goto 172
c--- end of loop through kd/decay species
192   continue
      nkdd=k
      if(nkdd.ne.0) write(42,1078) nkdd
1078  format(' finished reading kd/decay species, nkdd  = ',i3/)
c
      do k=1,nkdd
         ndum=0
         do i=1,npri
            if (nakdd(k) .eq. napri(i)) ndum=1
         end do
         if (ndum .eq. 0) then
             write(42,213)
213          format(' Error: not found in the primary species list')
             stop
         end if
      end do
c
      do i=1,npri
         kddp(i)=0
         do k=1,nkdd
            if (napri(i) .eq. nakdd(k)) then
                kddp(i)=k
            end if
         end do
      end do
c
c     keep values of nmin, ngas to the end of waters initialization
      nminold=nmin
      ngasold=ngas
c
c  establish initial set of stoich. coef. and logK values for all
c  derived species, minerals, gases, etc.
c  ------------------------------------------------------
c
      write (42,"(/' Temporarily stop reading chemical.inp to read',
     &     ' thermodynamic database....')")
      write (32,"(/' Temporarily stop reading chemical.inp to read',
     &     ' thermodynamic database....')")
      call readtherm_HKF                                       
c
      write (42,"(' ....back to reading chemical.inp')")
      write (32,"(' ....back to reading chemical.inp')")
c
c  read exchange cations of the system
c-------------------------------------------
c     read header label (needed)
c      read  (41,'(a)', err=9006)    label   ! Exchange cations of the system
c      write (42,'(a)') label
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo

      nexc=0
c     reads first record
c 
      if(nexc.ne.0) read (41,"(a200)") inprec
      if(index(inprec(1:20),"*").ne.0) goto 323     !end of kd/decay species
c
c     reads number of exchangeable sites
c
      read(inprec,*,end=9006,err=9006) NXsites, Mod_Xsl      
!
      if (Mod_Xsl .le. 0 .or. Mod_Xsl .ge. 3) Mod_Xsl = 1     ! Default value
!
      write (42,"(5x,'Number of exchangeable sites: ',i3)") NXsites
!
      write (42,"(5x,'Model for exchange dependence on water',
     &               ' saturation: Mod_Xsl =',i3)") Mod_Xsl
      write (42,"(10x,'=1: Simply divide by water saturation ')")
      write (42,"(10x,'=2: Only the wetted proportion ')")
!
      IF (MXsites .LT. NXsites ) THEN
         write (42,84) MXsites
         WRITE (*,84)  MXsites
84       FORMAT(' Error: Maximum number of exchanged sites (MXsites)'
     &          ' was exceeded. Current MXsites = ',i3)
         STOP
      END IF
c
c     Need another label above the list of exchange species (needed)
c     read(41,'(a)',err=9006) label
c      write (42,'(a)') label
      write(42,"(15x,'master     convention    selectivity')")
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      i=0
      nx = 0   ! master exchanged species (TX, 11/20/2001)
c---Start of loop through exchange species
330   continue
c
      if(i.ne.0) read(41,"(a200)",err=9006) inprec
      if(index(inprec(1:20),"*").ne.0) goto 320     !end of exchange species
c
      i=i+1
      if(i.eq.1) write(42,"('EXCHANGEABLE CATIONS')")
      IF(i.gt.mexc) THEN
          write (42,85) MEXC
          WRITE (*,85)  MEXC
85        FORMAT(/' Error: maximum number of exchanged species (MEXC)'
     1    ' was exceeded.  Current MEXC= ',i3)
          STOP
      END IF
c
      read(inprec,*,end=9006,err=9006) naexc(i), idum1(i), idum2(i),
     &         (ekxM(isite, i), isite =1, NXsites)  ! Selectivity for multi-sites
c
      label(1:20) = naexc(i)
      CALL name_conv(label)
      naexc(i) = label(1:20)
c
      write (42,8009) naexc(i),idum1(i),idum2(i),(ekxM(isite, i),
     & isite=1,NXsites)  ! Selectivity for multi-sites
8009  format (1x,a17,i3,11x,i3,10x,5f6.2)
      do j=1,npri
         if (napri(j).eq.naexc(i)) then
            nbx(i)=j
            do k=1,npri
               stqx(i,k)=0.0d0
               if (j.eq.k) stqx(i,k)=1.0d0
            end do
            go to 335
         end if
      end do
c
335   continue
c
c      in case exch. species is not in list of primary species
      if (nbx(i).eq.0.0d0) then
        write (32,"('Exchange species not in list of components: ',
     &    a20)") naexc(i)  
        write (42,"('Exchange species not in list of components: ',
     &    a20)") naexc(i)  
        stop
      end if

      if (z(nbx(i)).le.0.0d0) then
         write (32,*) 'error selecting a non cation for exchange'
         write (42,*) 'error selecting a non cation for exchange'
         stop
      end if
c
      if (idum1(i).eq.1.and.nx.eq.0) then
         nx=i
         iex=idum2(i)
      end if
      goto 330
c--- end of loop through exchange species
320   continue
      nexc=i
c
      if (nx .eq. 0) then
         write (42,*) 'error not selecting a master exchanged species'
         stop
      end if
c
      do i=1,nexc
       if (idum2(i).ne.iex) write(42,*) 'convention index assumed= ',iex
      end do
323   continue
      if(nexc.ne.0) write(42,927) nexc
927   format(' finished reading exchanged species, nexc = ',i3/)
c
c-----------finished reading first part of chemical inp----------
c
      write(42,"(80('-')/' SUMMARY OF THERMODYNAMIC DATA FOR THE'
     &' CHEMICAL SYSTEM'/)")
      call echotherm_HFK       !echo data in CHEMICAL.OUT file
c
c  create a permanent order for the names of aq. species
c
      do i=1,npri
       naaqt(i)=napri(i)
      end do
      do j=1,naqx
       naaqt(npri+j)=naaqx(j)
      end do
c
c  indexing h2o,h,oh,xoh and e- species:
c
      naqt=npri+naqx
      do i=1,naqt
        if(naaqt(i).eq.'h2o' .or. naaqt(i).eq.'H2O') NoTrans(i)=1
        if(naaqt(i).eq.'h2o' .or. naaqt(i).eq.'H2O') nw  = i
        if(naaqt(i).eq.'h+'  .or. naaqt(i).eq.'H+')  nh  = i
        if(naaqt(i).eq.'oh-' .or. naaqt(i).eq.'OH-') noh = i
        if(naaqt(i).eq.'e-') ne=i
        if(naaqt(i).eq.'o2(aq)' .or. naaqt(i).eq.'O2(aq)') no2aq = i
        if(naaqt(i).eq.'xoh' .or. naaqt(i).eq.'XOH')       nd    = i
        if(naaqt(i).eq.'h2(aq)' .or. naaqt(i).eq.'H2(aq)') nh2   = i
c
c.......Use aqueous complex to model surface complexation
c
        if(naaqt(i).eq.'s_xoh' .or. naaqt(i).eq.'S_XOH')   ns = i
c
      end do
c
c-----For co2 consumption by mineral phase using EOS2 module
c------------------------For H2 generation using EOS5 module
c
      do i=1,npri
        if(napri(i).eq.'hco3-'   .or. napri(i).eq.'HCO3-')   nco2=i
        if(napri(i).eq.'co3-2'   .or. napri(i).eq.'CO3-2')   nco2=i
        if(napri(i).eq.'co2(aq)' .or. napri(i).eq.'CO2(aq)') nco2=i
      end do
c
      do i=1,ngas
        if(nagas(i).eq.'co2(g)' .or. nagas(i).eq.'CO2(g)') nco2g=i
        if(nagas(i).eq.'h2(g)' .or. nagas(i).eq.'H2(g)')   nh2g=i
c
        if(nagas(i).eq.'h2s(g)' .or. nagas(i).eq.'H2S(g)') nh2sg=i
        if(nagas(i).eq.'ch4(g)' .or. nagas(i).eq.'CH4(g)') nch4g=i
        if(nagas(i).eq.'so2(g)' .or. nagas(i).eq.'SO2(g)') nso2g=i
      end do
c
c.....For gas generation using TMVOCs
c
      do i=1,npri
        if(napri(i).eq.'ch4(aq)' .or. napri(i).eq.'CH4(aq)') nch4 = i
        if(napri(i).eq.'so4-2'   .or. napri(i).eq.'SO4-2')   nh2s = i    
        if(napri(i).eq.'hs-'     .or. napri(i).eq.'HS-')     nh2s = i    
        nso2 = nh2s                                                      
      end do
c
c.....Get stoichometric coefficients for aqueous kinetics
c
      CALL AQKIN_STOICHOMETRY
c
c
c--------------
c.....For mineral kinetic rate contributed from dependent species
c--------------
c
      if (ndep .gt. 0)   then
         do ik=1,nmkin
c
c--------------------------------For dissolution rate constant
c
            if (idep(ik) .ge. 2)  then
               do ii=1,ndis(ik)         ! loop over mechanism
                  do isp=1,nspds(ik,ii) ! loop over species involved in one mechanism
                     icount=0
                     do iii=1,naqt
                        if (naaqt(iii).eq.nadis(ik,ii,isp)) then
                           icount=1
                           ids(ik,ii,isp)=iii
                           go to 119
                        end if
                     end do
c
                     if (icount.eq.0) then
                        write(42,128) namin(nmequ+ik)
                        stop
                     end if
119                  continue
                  end do
               end do
            end if    ! end of dissolution
c
c--------------------------------For precipitation rate constant
c
            if (ideprec(ik) .ge. 2)  then
               do ii=1,npre(ik)
                  do isp=1,nsppr(ik,ii)
                     icount=0
                     do iii=1,naqt
                        if (naaqt(iii).eq.napre(ik,ii,isp)) then
                           icount=1
                           idsp(ik,ii,isp)=iii
                           go to 129
                        end if
                     end do
c
                     if (icount.eq.0) then
                      write(42,128) namin(nmequ+ik)
128                   format (//'Error: The species involved in the  ',
     &                'kinetic rate constant for mineral:  ',a15,
     &                 /'       is not found in aqueous species list.')
                      stop
                     end if
129                  continue
                  end do
               end do
            end if  ! end of precipitation
c
c----------------------------------------------------------------
c
         end do

c        initialize derivatives
         do ik=1,nmkin  ! loop again through minerals
           do iik=1,mechm
             do iiik=1,npri
                der_ds(ik,iik,iiik) = 0.d0
                der_pr(ik,iik,iiik) = 0.d0
             enddo
           enddo 
         enddo
c
      end if
c
c
c------------ For adsorbed species, we create various pointers and make various checks
      if(nads.gt.0) then
c
c      checks if surface mineral is included in the system
       do n=1,nsurf
         if(naads_min(n).eq.'no_mineral') then
           m_index(n)=0    ! pointer to mineral index 1 to nmin
c
         elseif(naads_min(n)(1:7).eq.'surface') then
           m_index(n)=-1
         else
           ifound=0
           do m=1,nmin
             if(namin(m).eq.naads_min(n)) then
               ifound=1
               m_index(n)=m   !pointer to mineral index 1 to nmin
             endif
           enddo
           if(ifound.eq.0) then
           write(42,
     &       "(/5x,'Cannot find the mineral surface: ',A20,
     &       ' in the list of minerals.  Check this name')")
     &          naads_min(n)
           write(32,
     &       "(/5x,'Cannot find the mineral surface',A20,
     &       ' in the list of minerals.  Check this name')")
     &          naads_min(n)
             stop
           endif
         endif
       enddo
c
c      pointer to mineral surface index for each surface species 1-nads
       do j=1,nads
         iad_surf(j)=0  !points to index of surfaces
         iad_neu(j)=0   !points to index of primary surface species
       enddo
c
       do j=1,nads            !loops over derived surface species
         ncp=ncpad(j)
         ifound=0
         do n=1,ncp
           do i=npaq+1,npri
             kk=i-npaq                     !index of primary surface species
             if(icpad(j,n).eq.i) then
               iad_surf(j)=isurfp(kk)      !points to the surface index
               iad_neu(j)=i                !points to the primary surface species
               ifound=1
             endif
           enddo
         enddo
         if(ifound.eq.0) then
c        note: we should never get here because we already checked
c          when reading the database that all components are there.
           write(42,
     &       "(5x,'Cannot find a primary surface species in the'
     &       ' stoichiometry of ', A20)") naads(j)
           write(32,
     &       "(5x,'Cannot find a primary surface species in the'
     &       ' stoichiometry of ', A20)") naads(j)
           stop
         endif
       enddo
c
      endif


c-------------------------------------------------- -----------------------------
c                                                   
c   NSEC = total number of secondary components (derived species, mineral,
c    gases, adsorption, and exchange species) in the system
c     Note: assigment of logKs (as F(t)), stoichio, etc for speciation
c     calculations are now made through a call to ASSIGN further below
c
c-----------------------------------------------------------------------------
c  Finished reading the components of the chemical system, their stoichiometries,
c  and logK values.  Now we start reading the water compositions (initial waters
c  boundary waters, and recharge waters).
c  For each water: data is read first, then speciation calcs; then data for
c  next water is read, then next speciation; and so on.
c
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      niwtype=0
      nbwtype=0
c
      if(index(inprec(1:20),"'*'").ne.0) goto 301         ! end of water types
      read(inprec,*,end=9007,err=9007) niwtype,nbwtype ! number of initial and boundary waters
      nnn=niwtype+nbwtype  ! total number of water types
      if(nnn.ne.0) write(42,1202)
1202  format(/80('-')/,' INITIAL AND BOUDARY WATERS')
      IF(IRESTART .EQ. 1) GOTO 989
      IF(MBOUND.LT.NBWTYPE) THEN
         write(42,95) NBWTYPE
         WRITE(32,95) NBWTYPE
95       FORMAT(/' Error: MBOUND should be great than',I4)
         STOP
      END IF
      MAXITCH  = 0        ! Maximum iterations for solving whole chemistry
      AVERITCH = 0.d0     ! average itreations for solving chemistry
      COUNTCH  = 0.d0     ! counts number of times we go through iteration process (?)
c
c ----------------------------------------------------------------
c  Loop 300 starts here: loops through each water type
c  --------------------
c  For each water:
c   Reads water type and composition
c   Calculate initial speciation of the water (w/o minerals)
c   Assign concentrations of aq. species at nodes of the flow model
c
      do 300 nibw=1,nnn
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
      Pt=1.d0     ! default initialization pressure (for speciation)
      tc2=25.d0   ! default initialization temperature (for speciation)
      iwtype=1    ! default water type
      read(inprec,*,err=285) iwtype, tc2, Pt
c
      if(Pt.le.0.d0) Pt=1.d0 
c
c     iwtype is water type, but it is really printed (assigned?) as nibw below
c     tc2 is temperature for initial speciation removed ITC2 (not used)
c
285   write (42,*)

c     pressure and ref pressure for initial speciaton
c     these are updated later for pressure correction on logK's
      call SAT(tc2,psat)
      p0bar=psat/1.d+5
      if(Pt.lt.p0bar) Pt=p0bar 
c
      if (nibw.le.niwtype) then
         write(42,8012) iwtype
8012     format(' Initial Water Type',i3)
      else
         write(42,8120) iwtype
8120     format(' Boundary (+ injection) Water Type',i3)
      end if
      write(42,8011) tc2, Pt
8011  format(' T =',f7.2,' degrees Celsius, P =',F7.2,' bar')
      write(42,"(8x,'icon',8x,'guess',8x,'ctot',8x,'constrain')") 
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      ncount=0
c
c--- start of loop through each analyte for the current water analysis
290   continue
c
      if(ncount.ne.0) read(41,"(a200)", err=9007) inprec
      if(index(inprec(1:20),"*").ne.0) goto 275     ! end of analysis
c
      read(inprec,*,end=9007,err=9007) nadum2,icon2,cguess2,ctot2,
     &    nameq,qksat
      label(1:20)=nadum2                     
      call name_conv(label)
      nadum2=label(1:20)                    
      if(nameq(1:1).ne.' ') Then
        label(1:20)=nameq                    
        call name_conv(label)
        nameq=label(1:20)                   
      endif
      write (42,8005) nadum2,icon2,cguess2,ctot2
8005  format (1x,a8,i1,8x,e10.4,8x,e10.4,8x)
c
c  nadum2   Name of primary species j
c  icon2 -> icon(j)  This is a constraint flag for concentrations:
c                =1 total concentration is read
c                =2 the species concentration will be fixed by gas/mineral equilibrium
c                =3 tt(j) is read as a fixed activity for that species
c                =4 the species concentration will be set by charge balance
c  ctot2 -> tt(j) -> u2(j) -> utold(node,j) Total concentration of j-containing species in the water or
c  concentration as defined by icon flag above, in moles/kg solution
c  !!for water, enter kg of solution (i.e. 1 kg). It is converted to moles later!
c  **** ns5/30/02 made changes to input kg water instead of liter water
c  cguess2 -> cguess(j) -> cp(j) -> cn(node,i) Estimated concentration of actual species j  (enter in moles/liter solution)
c  CTOT2 -> ictot(j) -> IUB(n.bound,j) and IUR(n.recharge,j)   Not used at present??
       do j=1,npri
        if (napri(j).eq.nadum2) then
         icon(j)=icon2
         cguess(j)=cguess2
         tt(j)=ctot2
         ncount=ncount+1
c
c........For fixing concentration with gas or mineral at T
         mineq(j)=0
         igaseq(j)=0
         iflag=0
         if(icon(j).eq.2) then    ! ignore gas/mineral data if icon <> 2
          do m=1,nmin
           if(nameq.eq.namin(m)) then
               mineq(j)= m
               sateq(j)= qksat
               ncp=ncpm(m)
               do n=1,ncp
                 i=icpm(m,n)
                 if(i.eq.j) iflag=1
               enddo
             endif
           enddo
c
          do m=1,mgas
            if(nameq.eq.nagas(m)) then
               igaseq(j)= m
               sateq(j)= qksat
               ncp=ncpg(m)
               do n=1,ncp
                 i=icpg(m,n)
                 if(i.eq.j) iflag=1
               enddo
             endif
           enddo
           if(igaseq(j).eq.0.and.mineq(j).eq.0) then
             write(32,*) 'water composition: ',napri(j), 'icon=2',
     &      ' cannot find mineral or gas specified to fix this species'
             write(42,*) 'water composition: ',napri(j), 'icon=2',
     &      ' cannot find mineral or gas specified to fix this species'
            stop
           else if(iflag.eq.0) then
             write(32,*) 'water composition: ',napri(j), 'fixed by ',
     & nameq,' the mineral or gas specified does not contain ',napri(j)
             write(42,*) 'water composition: ',napri(j), 'fixed by ',
     & nameq,' the mineral or gas specified does not contain ',napri(j)
            stop
           end if
         end if

         goto 290
        end if
      end do
c
      go to 290
c--- end of loop through components of water analysis
c
275   if (ncount.ne.npri) then
       write(32,*) 'error in number of primary species in zone=',iwtype
       write(42,*) 'error in number of primary species in zone=',iwtype
       stop
      end if
c
c --  Normalize tt to one kg of water
      do j=1,npri
        if(j.ne.nw.and.icon(j).eq.1)tt(j)=tt(j)/tt(nw)
      enddo
c--For water mass balance, reset tt(nw) as total moles water
c   rmh2o is a constant defined in parameter.inc = 55.5 moles water/kg water
c   dliq = liquid density in g/cc or kg/l (now assume 1).
c   sumsalt=total kg of salt (per 1 l solution), now assume 0.
c
      tt(nw)=rmh2o    !total moles water (for 1 kg water)
c
c--------------- always initialize u2 just in case (moved from below)
c
      do j=1,npri
         u2(j)=tt(j)
      end do
c
c     initialize gas fugacity coeffients
      do j=1,ngas
         gamg(j)=1.d0
      end do
c
C-------If not initializing boundary and recharge solutions
c inibound is a flag read in solute.inp for recharge chemistry (0=none, 1=yes)
c
c   Removed bunch of lines here for icon=4 and icon=5 cases (options to fix activities
c   of some aq. species with minerals or gases - never tested - problems identified)
c
c   Assign trial concentrations of primary species to arrays cp, cs, and ct
c   These are in moles per kg water (molal units) for chemistry calculations.
c   They will be reassigned later to the matrix c in moles/l liquid (molar units) at
c    each node of flow model for transport calculations.
c
        kk=0
        do i=1,npri
          cp(i)=cguess(i)
        end do
        cp(nw)=1.d0  ! water concentration
        do j=1,naqx
          cs(j)=0.0d0
        end do
c
        tc(1)=tc2      ! sets temperature
        i1 = 1
        tk2 = tc2+273.15d0
c
        inode = 0
        call fugacomp(inode)  
c
c        npri=npaq   ! to skip surface complexation calcs initially
        call assign
        izero = 0
        call nrinit(izero) ! zero flag is to skip surface speciation calcs
c
c.....Reset the tt to total concentration for correct balances in chdump
c     if concentration was calculated from input activity, min/gas equilibria,
c     or charge balance (only affects printout of balances in chdump)
c
      do i=1,npri
        if(icon(i).ne.1.and.icon(i).ne.5) tt(i)=u2(i)
      enddo
      do i=1,npri
        if(icon(i).eq.5) ttt(i)=tt(i)
      enddo
c
      izero = 0
      tsim0 = 0.d0
      call chdump(tsim0,izero,iterch)
c
c     moved this here from ninit
c     case where we reached the maximum iterations (no convergence)
      if(iterch.ge.maxitpch) then
        WRITE(32,"(/' ERROR: convergence problem in initialization of'
     &' water composition, Please adjust convergence criteria'   
     &' regarding chemical iteration and initial guess of concentration'
     &' of primary species')")
        stop
      endif
c
      if(str.gt.stimax) then
c       in case we exceed ionic strength (stimax) with initial water
        write(32,*) ' Maximum ionic strength (stimax) exceeded'
        stop
      endif
c
c
      IF(ITERCH .GT. MAXITCH) MAXITCH=ITERCH
      AVERITCH=AVERITCH+DBLE(ITERCH)
      COUNTCH=COUNTCH+1.d0
c
c     go to boundary waters
      if (nibw.gt.niwtype .and. nibw.le.niwtype+nbwtype) go to 111
      if (nibw.gt.niwtype+nbwtype) go to 112
c
c   Assign concentrations, pH, and starting temperature
c   to each node of flow model for each initial water type.
c   For boundary waters, see afer 111 label
c   ---------------------------------------------------------
c
c   c(n,j) = conc. of aq. species j (mol/kgw) at each node n
c   ctot(n,j) = total concentration of aq. species j (mol/kgw) at each node n
c   u(n,j) = total concentration of aq. species j (mol/liter) at each node n
c     note moles per liter = mol/dm3 = molar units
c   izoneiw flag for each node's initial water type (read in solute.f from solute.inp)
c
c   We need a conversion factor to convert concentrations
c   from chem module in moles/kg h2o liq to moles/vol liq
c   cp(nw) contains kg of water
c
      sumsalts=0.d0   !sum of salts weights in kg (assume zero for now)
c
      do n=1,nnod
c---------- get water density for each node in kg/m3
c Watch!! Dry nodes have zero densities ==> concentrations will become zero
c
        NLOC2=(N-1)*NSEC*NEQ1
        NLOC2L=NLOC2+NBK
        densw=PAR(NLOC2L+4)
        IF (IEOS.EQ.9) densw=PAR(NLOC2+4)
        if(densw .gt. 2.d3) densw=1.d3     !???
        if(densw.eq.0.d0)   densw=1.d3
        dliq=densw/1.d3   ! density in kg/l
        h2ofact=cp(nw)/(cp(nw)+sumsalts)*dliq   ! factor in kg h2o liq/vol liq
c
       if (izoneiw(n).eq.iwtype) then   ! if type matches the one read in solute.inp
         do i=1,npri                    ! for each element (node)
          if(i.le.npaq) then ! aqueous species
           c(n,i)=cp(i)
c          cp(nw) contains kg water after initial speciation
           if(i.eq.nw) c(n,i)=rmh2o
c          u2 was initially set for 1 kg water (later it is total moles per initial liter)
c           but water after speciation may change
c          add division by cp(nw) below to get true molalities and molarities
           ctot(n,i)=u2(i)/cp(nw)
           utold(n,i)=u2(i)*h2ofact
c--------- add mult by h2ofact below
           ut(n,i)=u2(i)*h2ofact   !ut will get overriden by transport - put this here in
c          !  case there are components we do not move (e.g. water)
           if(icon(i).eq.5) ttt_nod(n,i)=ttt(i)
           icon_nod(n,i)=icon(i)
          else  !surface complexes (excluded so far from calculations)
           c(n,i)=cguess(i)
c          total concentrations will be computed later from surf.area & site density
          endif
         end do
c
c--------- add loop to assign default gas fugacity as equilibrium value
         do i=1,ngas
           pfug(n,i) = cg(i) ! gas fugacity computed in cmq_cp (no gas saturated)
         enddo
c
         do j=1,naqx
          c(n,npri+j)=cs(j)
         end do
          ph(n)=ph2
          str_node(n) = str        ! Save ionic strength for all nodes
          aw(n)       = gamp(nw)   ! Save water activity for vapor pressure lowering calculation Pitz zh
c
       end if
      end do
      go to 300
c
c     boundary waters
c     only need to store total concentrations in ub
c  
c
111   nnb=nibw-niwtype
        do i=1,npri
           ub(nnb,i)=u2(i)      !u2 is moles (per initial kg water), ub is molar
        end do
      go to 300
c
  112  continue
       h2ofact= 1.d0   !factor in kg h2o liq/vol liq
c
300   continue
c   END OF LOOP 300 THROUGH EACH INPUT WATER AND INITIAL SPECIATION
c------------------------------------------------------------------
301   continue
      if(nnn.ne.0) write(42,1057) niwtype,nbwtype
1057  format(/' Number of initial and boundary waters:',i3,i3/)
c
      if(countch.ne.0.d0) then   !test for zero division
       AVERITCH=AVERITCH/COUNTCH
      else
       averitch=0.d0
      endif
c
c     restitute nmin, ngas of the system
      nmin=nminold
      ngas=ngasold
c
c----------------------------------------------------------------------
c
c Now reads initial mineral, gas, adsorption, and exchange zones in the
c  CHEMICAL.INP file, and assigns them to nodes of the flow model
c
c  initial mineral zones
c  ---------------------
c
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      nmtype=0
c     add initialization of noreact (flag for not reacting specific mineral)
      do m=1,nmin
         do n=1,nnod
           noreact(n,m)=0
         enddo
      enddo
c
      if(index(inprec(1:20),"*").ne.0) goto 498     ! end of mineral zones
c
      read(inprec,*,end=9008,err=9008) nmtype        ! number of mineral zones
      if(nmtype.ne.0) write(42,"('INITIAL MINERAL ZONES')")
c----loop 400 through each mineral zone
      do 400 i=1,nmtype
c
       do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       read(inprec,*,err=9008) imtype

       write (42,8006) imtype
8006   format (/1x, 'MINERAL ZONE= ',i3)
c       read(41,'(a)',err=9008) label
c       write (42,'(a)') label
       write(42,"('mineral         vol.frac.')")
c
c------add initialization of vol2, amin2, and rad2
       do m=1,nmequ
         vol2(m)=0.d0
       enddo
       do m=1,nmkin
         vol2(nmequ+m)=0.d0
         amin2(m)=0.d0
         rad2(m)=0.d0
       enddo
       sumvol=0.d0
       iread=0
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
c
c---start to loop through each mineral for current zone
620   continue
c
        if(iread.ne.0) read(41,"(a200)",err=9008) inprec
        if(index(inprec(1:20),"*").ne.0) goto 630     !end of current mineral zone
        iread=iread+1
c
        read(inprec,*,end=9008,err=9008) nadum2,voldum,ikin4   ! for kinetics
c
c      nadum2 -> namin(m) name of mineral m
c      voldum -> vol2(m)  volume fraction of mineral m in zone (-> pre(node,m) = vol2/vmin)
c---------- Note that the volume fraction definition now excludes porosity
c      ikin4  kinetic flag for mineral m (=0 no kinetics;  =1 kinetics)  Not used anymore
c      ikin4 = 2 flag for not reacting this mineral  
c
        label(1:20)=nadum2                      
        call name_conv(label)
        nadum2=label(1:20)                     
c
c     assigns the volume fraction to the mineral
        ifthere=0
        do m=1,nmequ
         if (namin(m).eq.nadum2) then
           vol2(m)=voldum
           ikin4m(m)=ikin4      
           sumvol=sumvol+voldum
           ifthere=1
         end if
        end do
c
c---------- moved assignment of vol2,rad2 and amin2 for kin minerals
c here to avoid potential major problem if mineral order is
c changed in input file. Add separate loop for kin minerals
c Note: variable ikin4 is now obsolete in mineral zones (not used)
c
        do m=1,nmkin
         if (namin(nmequ+m).eq.nadum2) then
           vol2(nmequ+m)=voldum
           ikin4m(nmequ+m)=ikin4    !ns09/09
           sumvol=sumvol+voldum
           read(41,*,err=9008) rad4,amin4,imflg
c
c          rad4 -> rad2(i) -> rad(node,i)  mineral grain radius
c          amin4 -> amin2(i) -> amin(node,i)  mineral surface area
c          imflg = 1, then area is in m2/m3 medium(use for fracture mineral area)
c----------Note that the surface area is now in terms of cm^2/g mineral
c          i is order index of minerals that have kinetics specified
c
           rad2(m)=rad4
           amin2(m)=amin4
           imflg2(m) = imflg
           ifthere=1
         end if
        end do
c
        if(ifthere.eq.0) then                            ! added if-block
         write(42,"(/5x,'Mineral ',A12,' is not included'
     +      ,' in the list of minerals for the system')") nadum2
         stop
        end if
       go to 620      ! loop back to read next mineral in current zone
c---end of loop through minerals for each zone
c
630    continue       ! input done for current mineral zone
c
c-------------- print mineral zone array variables (i.e. echo the
c  arrays themselves, not the dummy input variables) to make sure
c  the correct data was read and assigned to correct arrays.
c
       do m=1,nmequ
          write(42,"(1x,a15,e10.4,' equilibrium')") namin(m),vol2(m)
       enddo
       do m=1,nmkin
           nik = nmequ+m
        write(42,"(1x,a15,e10.4,' kinetic ')") namin(nik),vol2(nik)
        write(42,"(2E12.4,i5)") rad2(m),amin2(m),imflg2(m)
       enddo
       write(42,"(1x,'Sum of mineral volume fractions in zone: '
     +   ,e10.4)") sumvol       
c
c     assign mineral zone parameters to the nodes
c     si=saturation index; p= mol min/dm3 sol (mol/liter);
c
       do n=1,nnod
        if (izonem(n).eq.imtype) then
         do m=1,nmin
          pre(n,m) = (vol2(m)/vmin(m))*(1.d0-phi0(n))
          pinit(n,m)=pre(n,m)
c         add array for making the mineral unractive 
          if(ikin4m(m).eq.2) noreact(n,m)=1    !
c         minor fix
          pre0(n,m)=pre(n,m)
         end do
C
c.... Surface area conversions (3 types of input)
c     Changed numbering and added option 2
c Note, all input surface areas are converted to units of m2_mineral/m3_mineral
c
c imflg2(i) = 0 -> cm^2/g mineral
c imflg2(i) = 1 -> m^2/m^3 mineral
c imflg2(i) = 2 -> m^2 rock area/m^3 medium
c imflg2(i) = 3 -> m^2 rock area/m^3 medium_solids
cns6/10 imflg2(i) = 4 -> constant rate is input in mol/sec (no surface area)
c
         do iki=1,nmkin
c
             m=nmequ+iki
c
c----------- for minerals at equil.ppt and kin. dissol.
           if(kineq(iki).ne.0) then
             pre(n,kineq(iki))=0.d0
             pinit(n,kineq(iki))=0.d0
           endif
c
c... Set per grid block
           rad(n,iki)= rad2(iki)
c          add flag per grid block
           imflag(n,iki) = imflg2(iki)
c
c          mineral volume fraction in m3_miner/m3_medium(total)
c          vfrac=pre(n,m)*vmin(m)
           vfrtmp = 1.d0-phi(n) 
           vfrac=dmax1(vfrtmp*rnucl(iki),1.d-10,pre(n,m)*vmin(m))
c
c          Initialize number of grains (used in routine rsfarea)
           grains(n,iki)=1.d0
           if(rad2(iki).ne.0.d0) then
             vfrac_ini=dmax1(vfrtmp*rnucl(iki),1.d-10)
             grains(n,iki)=max(1.d0,vfrac_ini*0.125d0/        
     &             (rad2(iki)*rad2(iki)*rad2(iki)))        !number of mineral grains from rnucl and initial grain radius 
             if(vfrac.gt.0.d0) 
     &             grains(n,iki)=grains(n,iki)*vfrac_ini/vfrac        !number of mineral grains from rnucl and initial grain radius 
           endif
c
c---------- Surface areas in cm^2/g to m^2/m^3 mineral
           if(imflg2(iki).eq.0)then
              amin(n,iki)= amin2(iki)*dmolwm(nmequ+iki)*0.1d0/
     +           vmin(nmequ+iki)

c---------- Surface areas in m^2 mineral/m^3 mineral
           elseif(imflg2(iki).eq.1)then
              amin(n,iki)= amin2(iki)

cns6/09 add this option and changed order
c---------- Surface areas in m2_mineral/m3_medium(total)
           elseif(imflg2(iki).eq.2) then
cels6/10/09              if(vfrac.eq.0.d0) goto 9905  !abort - this option requires vfrac non zero
              amin(n,iki) = amin2(iki)/vfrac
c
c---------- (m^2/m^3 medium )
           elseif(imflg2(iki).eq.3)then
cels6/10/09              amin(n,iki)= amin2(iki)*1.5d0
              amin(n,iki)= amin2(iki)
c
cns6/10
           elseif(imflg2(iki).eq.4)then
              amin(n,iki)= 1.d0
c

           endif
         enddo
c
c-------------- store the proportion on unreactive minerals at each
c     node (i.e. when mineral volume fractions do not add up to 1).
c     This quantity is stored in pinit(node,nmin+1). WATCH UNITS!!!
c     This quantity is stored as:
c       (volume of unreactive minerals) / (volume of medium)
c     Do not confuse with other pinit values in moles/vol medium!!!
c
         pinit(n,nmin+1)=(1.d0-sumvol)*(1.d0-phi0(n))
c
        end if
       end do
400   continue
c--- end of loop through each mineral zones
c
498   continue
      if(nmtype.ne.0) write(42,1059) nmtype
1059  format(' Total number of mineral zones:',i3/)
c
c----------------------------------------------------------------------
c  initial gas zones
c  -----------------
c
      ico2gt0=0    ! =1, initial Pco2>0
      ih2gt0 =0    ! =1, initial Ph2 >0
c
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      ngtype=0
c
      if(index(inprec(1:20),"*").ne.0) goto 599           ! end of gas zones
      read(inprec,*,end=9018,err=9018) ngtype             ! number of gas zones
      if(ngtype.ne.0) write(42,"('INITIAL GAS ZONES')")
c--- start of loop through each gas zone
      do 403 i=1,ngtype
c       read(41,*,err=9018) igtype
       do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       read(inprec,*,err=9018) igtype

       write (42,905) igtype
905    format (/1x, 'GAS ZONE= ',i3)
c
       do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       write (42,"('gas      partial pressure')")
c
       iread=0
c
c--- start of loop through gases of current gas zone
621   continue
c
      if(iread.ne.0) read(41,"(a200)",err=9018) inprec
      if(index(inprec(1:20),"*").ne.0) goto 631     !end of gas zones
      iread=iread+1
c
      read(inprec,*,end=9018,err=9018) nadum2,voldum
c
      label(1:20)=nadum2                      
      call name_conv(label)
      nadum2=label(1:20)                     
c
c   nadum2 -> nagas(m) name of gas
c   voldum -> cg(m) -> pfug(node,m)   gas partial pressure in bars
c   Note: the partial pressure read in is an initial pressure.  It
c   is recomputed later.  You can input a very low P to compute a pressure
c   that reflects only the water chemistry.
c
      ifthere=0
      do m=1,ngas
         if (nagas(m).eq.nadum2) then
           cg(m)=voldum
           ifthere=1
c
           if (m.eq.nco2g .and. cg(m).gt. 0.0d0) then
             ico2gt0=ico2gt0+1    ! =1, initial Pco2>0
           end if
c
           if (m.eq.nh2g  .and. cg(m).gt. 0.0d0) then
             ih2gt0=ih2gt0+1      ! =1, initial Ph2>0
           end if
c
         end if
      end do
      if(ifthere.eq.0) then                    ! added if-block
         write(42,"(/5x,'Gas ',A12,' is not included'
     +      ,' in the list of gases for the system')") nadum2
         stop
      end if
      goto 621
c---end of loop through gases of current gas zone
631   continue
c
c
c------add printout of array with do loop below
      do m=1,ngas
         write(42,"(1x,a15,e10.4)") nagas(m),cg(m)
      enddo
c
c     Assign gas zones parameters to the nodes
c     For now, assume original gas temperature is 25 C
c     The gas constant gc is defined in parameter.inc
c
      rt=gc*298.15d0                !  added (gas RT)
c
c     gp= moles of gas per liter medium (per total volume rock+gas+liquid)
c     cg= initial gas partial pressure in bars.  It will change, especially
c         if the water chemistry does not reflect this partial pressure
c
      do n=1,nnod
        if (izoneg(n).eq.igtype) then
         do m=1,ngas
c
c--------------------------------------------  check against v q2.01
c  Only changes relative to v q2.0 are indicated
c  pfug was set equal to equilibrium pressure with initial water (after call
c  to nrinit).
c  Now we reset it to the input pressure if the later is specified > 0
c
c stop using this variable name - use pfug instead
           if(cg(m).gt.0.d0) pfug(n,m)=cg(m)
c-------- gp here is at 25 C - it is recalculated in couple
           GP(N,M)= CG(M)/rt*PHI(N)*SGOLD(N)        ! gas con./l medium
         end do
        end if
      end do
403   continue
c---end of loop through each gas zones
c
599   continue
      if(ngtype.ne.0) write(42,1062) ngtype
1062  format(' Total number of gas zones:',i3/)
c
c---------------------------------------------------------------------
c  permeability-porosity law zones
c
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      nppzon=0
c
      if(index(inprec(1:20),"*").ne.0) goto 499     ! end of permeability zones
c
      read(inprec,*,end=9012,err=9012) nppzon !number of permeability zones
      if(nppzon.ne.0) write(42,"('INITIAL PER-PORO ZONES')")
      do k=1,nppzon
        do 
         read(41,"(a200)",err=9001) inprec
         if(inprec(1:1).ne."#".and.inprec.ne.'') exit
        enddo
        read(inprec,*,err=9012) ippzon
c
        write(42,8019) ippzon
        write(42,"('perm law',8x,'a-par',17x,'b-par')")
c
        do 
         read(41,"(a200)",err=9001) inprec
         if(inprec(1:1).ne."#".and.inprec.ne.'') exit
        enddo
        read(inprec,*,err=9012) ipptyp,apppar,bpppar
c
        write (42,8020) ipptyp,apppar,bpppar
 8019  format (/' PERMEABILITY ZONE= ',i3)
 8020  format (1x,i3,11x,1e7.2,15x,1e7.2)
c
c...... Assign values to nodes
        do n = 1, nnod
          if(izonpp(n).eq.ippzon)then
            ikplaw(n) = ipptyp
            aparpp(n) = apppar
            bparpp(n) = bpppar
          endif
        enddo
c
      enddo
c
c     Add record marking the end (e.g., '*' or blank) for consistency
c     with input of other zones
      read(41,'(a)',err=9012) label
  499 continue
      if(nppzon.ne.0) write(42,1061) nppzon
 1061 format(' Total number of permeability zones=',i3/)
c
c----------------------------------------------------------------------
c  initial surface adsorption zones
c--------------------------------
c
      do i=1,nads
         cd(i)=0.d0
      enddo
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      ndzone=0
c
      if(index(inprec(1:20),"*").ne.0) goto 705 ! end of surface adsorption zones
c
      read(inprec,*,end=9009,err=9009) ndzone ! number of sorption zones
      if(ndzone.ne.0) write(42,"(' INITIAL SURFACE ADSORPTION ZONES')")
c--- start of loop through each adsorption zone
      do 700 k=1,ndzone
c        read(41,*,err=9009) idzone,  iequil   !zone number and flag for initial equilibration
        do 
          read(41,"(a200)",err=9001) inprec
          if(inprec(1:1).ne."#".and.inprec.ne.'') exit
        enddo
        read(inprec,*,err=9009) idzone,iequil   ! zone number and flag for initial equilibration
c
        write(42,"(/' ADSORPTION ZONE=',i3,' IEQUIL=',i3)")
     &        idzone, iequil
c
       do 
         read(41,"(a200)",err=9001) inprec
         if(inprec(1:1).ne."#".and.inprec.ne.'') exit
       enddo
       write(42,"(' Surface name',10x,'Units Flag  Surface Area')")  
       iread=0
c
c------start of loop through each ads species
        do i=1,npads
c         input: name primary surface species, surface area, site density (mol/m2)
c         iflag=0 surface area input in cm2/g_mineral
c         iflag=1 surface area input in m2/m3_mineral
c         iflag is ignored if surface species was defined tied
c           to naads_min(1:7)='surface', in which case input surface area is automatically
c           assumed to be constant and in units of m2/kgw.
c
c          read(41,*,err=9009) nadum2,iflag,sadsdum
c
        if(iread.ne.0) read(41,"(a200)",err=9008) inprec
        read(inprec,*,err=9009) nadum2,iflag,sadsdum
        iread=iread+1
c
          label(1:20) = nadum2
          CALL name_conv(label)
          nadum2 = label(1:20)
c
          write (42,"(3x,A20,5x,i3,5x,e10.4)") nadum2,iflag,sadsdum
c
          ifound=0
          do n=npaq+1,npri
            if(napri(n).eq.nadum2) then
              ifound=1
              i_phi=isurfp(n-npaq)  !index of mineral surface for this primary surface species
              m_min=m_index(i_phi)  !index in mineral list (1 to nmin) for this mineral surface
              exit
            endif
          enddo
          if(ifound.eq.0) then
           write(42,
     &       "(' Cannot find the adsorption zone species ',A20,
     &       ' in the list of primary species.  Check this name')")
     &          nadum2
             stop
          endif
c
          if(m_min.lt.0) then
c         case when we read in a constant surface area in m2/kgw
             s_area=sadsdum
c
          elseif(iflag.eq.0)then
c         convert from cm2/g_mineral to m2/m3_mineral
             if(m_min.gt.0) then
c            if surface is tied to a mineral
                  s_area = sadsdum*dmolwm(m_min)*0.1d0/    !dmolwm molecular weight in g/mol
     +             vmin(m_min)                             !vmin molar volume in dm3/mol
             elseif(m_min.eq.0) then
c            case when surface species is not tied to a specific mineral
c            assume average mineral density of 2.65 cm3/g
                  s_area = sadsdum/2.65d0*100.d0
             endif
c          elseif(iflag.eq.1)then
c         convert from m2/m3_fracture_medium_solids to m2/m3_mineral
c            factor of 1.5 accounts for areal fraction vs. volume fraction for spheres
c             s_area = sadsdum*1.5d0
c          elseif(iflag.eq.2)then
          elseif(iflag.eq.1)then
c            no conversion needed, already in m2/m3_mineral
             s_area = sadsdum
          endif
c
c         assign surface area (m2/m3_mineral) to matching node of zone
          do n=1,nnod
            if (izoned(n).eq.idzone) then
             supadn(n,i_phi)=s_area
             isurfeq(n)=iequil
            endif
          enddo
        enddo
        read(41,"(a200)",err=9009) inprec
  700 continue
c---  end of loop through each adsorption zone
c
c     initializes all potential terms
      do i=1,nnod
        do j=1,nsurf
          phip(i,j)=0.d0
        enddo
        do j=1,nads
          d(i,j)=0.d0   !mol/kgw secondary surface species
        enddo
      enddo
c
  705 continue
      if(ndzone.ne.0) write(42,1063) ndzone
 1063 format(' Total number of adsorption zones=',i3/)
c
c--- Note, computation of site concentrations are done later with
c    call to surfequil
c
c  initial Kd zones
c  -------------------------------------------------------------
c  nakdd ( ) : species name with Kd and decay
c
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      nkdtype=0
c 
      if(index(inprec(1:20),"*").ne.0) goto 879          ! end of Kd zones
c
      read(inprec,*,end=9011,err=9011) nkdtype           ! number of gas zones
      if(nkdtype.ne.0) write(42,"('INITIAL LINEAR.EQ. Kd ZONES')")
      if (nkdtype.gt.30) then
         write (42,*) 'Kd zone dimension (30) is not enough'
      end if
c
      do i=1,nkdtype
         do m=1,nkdd
            sden(i,m)=0.0d0    ! solid density
            vkd(i,m)=0.0d0     ! Kd values
            kdflag(i,m)=0
         end do
      end do
c
c-Read: Name,solid density (kg/dm**3), and Kd(l/kg=mass/kg solid / mass/l water)
c------------------------If density is zero vkd2 is retardation factor
c start looping through kd zones
c
      do i=1,nkdtype
        do 
          read(41,"(a200)",err=9001) inprec
          if(inprec(1:1).ne."#".and.inprec.ne.'') exit
        enddo

        read(inprec,*,err=9011) ikdtype

         write (42,907) ikdtype
907      format (/1x, 'Kd ZONE= ',i3)
c         read(41,'(a)',err=9011) label
c         write (42,'(a)') label

         do 
           read(41,"(a200)",err=9001) inprec
           if(inprec(1:1).ne."#".and.inprec.ne.'') exit
         enddo
         write(42,"(1x,'species   solid-density(Sden,kg/dm**3)',
     &   '  Kd(l/kg=mass/kg solid / mass/l')")  
         iread=0
c
c     start looping through kd species
623      continue
c
         if(iread.ne.0) read(41,"(a200)",err=9001) inprec
         if(index(inprec(1:20),"*").ne.0) goto 633                 !end of Kd zones
c
         read(inprec,*,err=9011) nadum2,sden2,vkd2,iflkd2
         iread=iread+1
c
         label(1:20) = nadum2
         CALL name_conv(label)
         nadum2 = label(1:20)
c
         write(42,'(1x,a15,e10.4,15x,e10.4,15x,I2)')
     &   nadum2,sden2,vkd2,iflkd2
c
         if (sden2.eq.0.0d0 .and. iflkd2.eq.0.and.vkd2.lt.1.0d0) then
            write (42,*) 'Error: Retardation factor mest be >=1'
            stop
         end if
         ndum=0
         do m=1,nkdd
            if (nakdd(m).eq.nadum2) then
               vkd(ikdtype,m)=vkd2
               sden(ikdtype,m)=sden2
               kdflag(ikdtype,m)=iflkd2
               ndum=1
            end if
         end do
         if (ndum .eq. 0) then
             write (42,*) 'Error: The Kd species name is not correct'
             stop
         end if
         go to 623
633      continue
      end do
c
      do i=1,nkdtype
        do m=1,nkdd
          if (sden(i,m).eq.0.0d0 .AND. vkd(i,m).eq.0.0d0) then
            vkd(i,m)=1.0D0
          end if
        end do
      end do
c
879   continue
      if(nkdtype.ne.0) write(42,1064) nkdtype
 1064 format(' Total number of Kd zones=',i3/)

c--initial zones of cation exchange
c-----------------------------------
      
      do 
        read(41,"(a200)",err=9001) inprec
        if(inprec(1:1).ne."#".and.inprec.ne.'') exit
      enddo
c
      nxtype=0
c 
      if(index(inprec(1:20),"*").ne.0) goto 805  ! end of cation exchange zones
c
      read(inprec,*,end=9010,err=9010) nxtype !number of exchange zones
      if(nxtype.ne.0) write(42,1065)
 1065 format('INITIAL ZONES OF CATION EXCHANGE'/'zone    ex. capacity'
     x'(meq/100 g solid or cmol/kg)')
c
      do 800 k=1,nxtype
        do 
          read(41,"(a200)",err=9001) inprec
          if(inprec(1:1).ne."#".and.inprec.ne.'') exit
        enddo
        read(inprec, *, err=9010) ixtype,(cecM(isite),isite=1,NXsites)

        write (42,8010)     ixtype, (cecM(isite),isite=1,NXsites)
8010    format (1x,i4,11x,5e10.4)
c
c       assign exchange parameters to each node
c
        do 795 n=1,nnod
           NMATIA=MATX(n)
           density_rock=drok(NMATIA)/1000.0d0
         if (izonex(n).eq.ixtype) then
c
           phi2 = phi(n)
           sg2  = sg1(n)
           sl2  = 1.0d0 - sg2
!
! ------------
!..........Get grid rock density for geochemical calculations such as exchange and sorption
! ------------
!
!...........Modified active fracture area for reaction
            a_fmr2 = (sl2-sl1min)/(1.0d0-sl1min) 
!
            phisl2 = phi2*sl2       ! phi2*sl2
            denss2 = denss(n)       ! Density of solid rock (kg/dm^3)
!
! ------------
!
           do i=1,naqt
              ct(i)=c(n,i)
           end do
!          Recalculate activity coefficient
           do i=1,npri
             cp(i)=c(n,i)
           enddo
           do i=npri+1,naqt
             cs(i-npri)=c(n,i)
           enddo
             izero = 0
             call dh_hkf81(no_ch)
c
           do i=1,npri
              gamt(i)=gamp(i)
           end do
           do i=1,naqx
              gamt(i+npri)=gams(i)
           end do
c
c........Loop over multi-sites
c
           do isite=1, NXsites
              cec2          = cecM(isite)
              cec(n, isite) = cecM(isite)
              if (cec2 .eq. 0.0d0 .or. sl2 .le. sl1min)   then
                do j=1,nexc
                  xcads(n, isite, j) = 0.0d0
               end do
               go to 794
              end if
c
              do j=1,nexc
                ekx(j) = ekxM(isite, j)    ! Selectivity
              end do
              call cx_ct
              do j=1,nexc
                xcads(n, isite, j) = cx(j)
              end do
794           continue
           end do    ! Multi-sites
c
         end if
795     continue
c
800   continue
c--- end loop through exchange zones
c----------------------------------------------------------------------
805   continue
      if(nxtype.ne.0) write(42,1066) nxtype
 1066 format(' Total number of exchange zones=',i4/)
      write (*,*) '   --> read input data complete'
      write (42,"(' --> Finished reading chemical.inp file')")
      write (32,"(' --> Finished reading chemical.inp file')")
c
c-----------------------------------------------------------
c initial speciation calcs including surface complexation (node by node)
c-----------------------------------------------------------
cns3/2010      if(nads.gt.0) then
      if(nads.gt.0.and.kcyc.lt.1) then
        write(42,"(/5x,'...start initial surface complexation calcs'/)")
        write(32,"(/5x,'...start initial surface complexation calcs'/)")
         call surfequil
        write(42,"(/5x,'...done initial surface complexation calcs'/)")
        write(32,"(/5x,'...done initial surface complexation calcs'/)")
      endif
c
989   CONTINUE
c
      close (41)
      close (42)
      return
c
9001  write (32,*) 'error reading title in chemical.inp'
      write (42,*) 'error reading title'
      write (42,"('last input record starts with: ', a40/)") label
      stop
9002  write (32,*) 'error reading primary species in chemical.inp'
      write (42,*) 'error reading primary species'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9102  write (32,"('error reading header of block following',
     &  ' primary species in chemical.inp')")
      write (42,"('error reading header of block following',
     &   ' primary species')")
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9202  write (32,*) 'error reading secondary species in chemical.inp'
      write (42,*) 'error reading secondary species'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9003  write (32,*) 'error reading minerals in chemical.inp'
      write (42,*) 'error reading minerals'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9004  write (32,*) 'error reading gases in chemical.inp'
      write (42,*) 'error reading gases'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9005  write (32,*) 'error reading surface in chemical.inp'
      write (42,*) 'error reading surface complexes'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
c
9026  write (32,"('error reading species with Kd and decay',
     &                ' in chemical.inp')")
      write (42,*) 'error reading species with Kd and decay'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9011  write (32,*) 'error reading Kd zone in chemical.inp'
      write (42,*) 'error reading Kd zone'
      stop
9022  write (32,*) 'error reading aqueous kinetics in chemical.inp'
      write (42,*) 'error reading aqueous kinetics'
      stop
C---------------------------
9705  write (32,"('error reading surface complexation model',
     &           ' parameters in chemical.inp')")
      write (42,"('error reading surface complexation model',
     &           ' parameters')")
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9006  write (32,*) 'error reading exchange cations in chemical.inp'
      write (42,*) 'error reading exchange cations'
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9007  write (32,*) 'error reading initial water zone in chemical.inp'
      write (42,*) 'error reading initial water zone=',iwtype
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9008  write (32,*) 'error reading initial mineral zone in chemical.inp'
      write (42,*) 'error reading initial mineral zone=',imtype
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9018  write (32,*) 'error reading initial gas zone in chemical.inp'
      write (42,*) 'error reading initial gas zone=',imtype
      write (42,"('last input record starts with: ',a40/)") inprec(1:40)
      stop
9009  write (32,"('error reading initial adsorption zone'
     &    ' in chemical.inp')")
      write (42,*) 'error reading initial adsorption zone=',
     &    idzone
      stop
9010  write (32,*)'error reading initial cation exchange zone=',ixtype
      write (42,*) 'error reading initial cation exchange zone=',
     &    ixtype
9012  write (32,*) 'error reading permeability-porosity in chemical.inp'
      write (42,*) 'error reading permeability-porosity zone=',
     &    ipptyp
      stop
9901  write (32,*) 'max number of rate mechanisms is exceeded, max=',
     &  mechm
      write (42,*) 'max number of rate mechanisms is exceeded, max=',
     &  mechm
      stop
9902  write (32,*) 'number of species in rate product exceeds max of,'
     &  , mechsp
      write (42,*) 'number of species in rate product exceeds max of,'
     &  , mechsp
      stop
9903  write (32,*) 'number of aqueous kinetic reactions exceeds max of,'
     &  , mrx
      write (32,*) 'number of aqueous kinetic reactions exceeds max of,'
     &  , mrx
      stop
c
      write(32,"(//2x,'Mineral: ',a8,' is initially absent',
     &  /5x,'Its surface area cannot be entered in m2/m3_medium'
     &  /5x,'Change the IMFLG flag in chemical.inp to another option'
     &       )") namin(m)        
      stop
      end
c
c-----------------------------------------------------------------------------
c
      subroutine nrinit(isurfcalc)
c
c*************** Solve initial chemical system by Newton-Raphson iteration *****************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/newton/amat(mtot,mtot),bmat(mtot),nmatc
      common/satgas1/isatg(mgas),moutg(mgas),nsatg !keep track gas saturation
      integer*8 isurfcalc,iterch,no_ch
      integer*8 indx(mtot)
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***** nrinit 1.0, 2003.7.30: Solve initial chemical'
     x' system by Newton-Raphson iteration')
c
c     isurfcalc is a flag = 0 to skip surface complexation calcs
c                         = 1 surf compl calcs, without surface equilibration
c                         = 2 surf compl calcs with surface equilibration (unchanging water composition)

c     mtot is the maximum matrix size (primary sp+derived sp+min+gas+ads+exch)
c     set in paramete.inc
      facmax=0.5d0
c
c  --- chemical iteration loop 800 (for first equilibration only)
c
      nsatg=0         !do not consider gases or minerals at initialization
      nsat=0
c
      do 800 k=1,maxitpch
        iterch=k
c
c       xh2o is kilograms of liquid water (we will solve for it)
        xh2o = cp(nw)   ! initial guess
c        xh2o = 1.d0    ! initial guess
        cp(nw)=1.d0     ! set cp(nw)=1. only during iterations
c
        no_ch=0
        call cheminit(no_ch,isurfcalc)
c
        call jacobinit(isurfcalc)
c
C-------------------------------------------------------Call LU solver
c
        i0 = 0
        call ludcmp (amat,nmatc,mtot,indx,dd)
        call lubksb (amat,nmatc,mtot,indx,bmat)
c
c       nmat=matrix dimension (squared matrix)
c       amat=jacobian matrix (partial derivatives)
c       mnr=maximum matrix dimension (defined in paramete.inc)
c       dd=defined in ludcmp (not used??)
c       indx=array used in slud/sbs or ludcmp/lubksb
c       bmat=the right hand side (-F's) on input, and the solution on output
c
        cp(nw) = xh2o   !stores kg water back in cp(nw)
c
c       actualizing unknowns and convergence criteria
        errmax=0.d0
        do 350 i=1,npri
          errx=bmat(i)                    !this is relative error
c         recover relative concentrations
          BMAT(I)=BMAT(I)*CP(I)
c
          if (dabs(errx).gt.errmax) errmax=dabs(errx)
350     continue
        if (errmax.gt.facmax) then
        do 360 i=1,npri
360     bmat(i)= bmat(i)*facmax/errmax
        endif
c
        do 370 i=1,npri
         cp(i)=cp(i)+bmat(i)
370     continue
c
        errmax2=0.d0
c
        if(isurfcalc.ne.0) then
c    actualize potential term phi for surface complex
c    -----------------
          if(npot.ne.0) then    ! npot=number of surfaces= number of potential terms
            errmax2=0.d0
            do m=1,npot
              kk=npri+m
              jj=ipoten(m)
              bmatphi=bmat(kk)
              if (errx.gt.errmax) errmax=errx
               IF (ABS(bmatphi) > 0.1d0) THEN
                   bmatphi = SIGN(0.1d0,bmatphi)
               END IF
               errx=dabs(fphi(m))
               if(errx.gt.errmax2) errmax2=errx
               phip2(jj)=phip2(jj)+bmatphi
            enddo
          endif

        endif
c
        if ((errmax.lt.tolch.and.errmax2.lt.1.d-6)
     &      .or. maxitpch.eq.1) then
          return
        end if
c
800   continue
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine jacobinit(isurfcalc)
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/newton/amat(mtot,mtot),bmat(mtot),nmatc
      common/initial/cguess(mpri),icon(mpri)
      common/eqinit/sateq(mpri),mineq(mpri),igaseq(mpri)   !NS3/06
      integer*8 isurfcalc
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***** jacobinit 1.0, 2003.7.30: Construct Jacobian'
     x' matrix for initial chemical system')
c
c       initialize amat and bmat terms equal zero
        do i=1,npri+npot
          do j=1,npri+npot
           amat(i,j)=0.d0
          end do
          bmat(i)=0.d0
        end do
c
c     isurf=0, skip surface complex speciation, initial speciation (icon is functional)
c     isurf=1, include surface complex speciation, without euqilibrating with surface
c        note: requires first a speciation with isurf=0 (icon is not functional)
c     isurf=2, include surface complex speciation but equilibrate with surface (leaving
c        concentration of aqueous species unchanged
c        note: requires first a speciation with isurf=0 (icon is not functional)
c
c     isurfcalc: flag to skip (=0) or include (>0) surface speciation calcs
c     for initial speciation, isurfcalc is always = 0
c
      if(isurfcalc.eq.0) then
c--------------------------------------------------------------------------
c--- Initial speciation, without surface calcs
c--------------------------------------------------------------------------
c
        nmatc=0
c
c       mass balance equations (npri)
        do 500 i=1,npri
           nmatc=nmatc+1
c          fixed total solute concentration
           if (icon(i).eq.1.or.icon(i).eq.5) then
             do j=1,npri
             amat(i,j)=du2(i,j)
             end do
             bmat(i)=tt(i)-u2(i)

c          Equilibration with gases or minerals
           else if (icon(i).eq.2) then
             if(mineq(i).ne.0) then

c             deriv of min. mass action w/respect to cp aq.
               m=mineq(i)
               ncp=ncpm(m)
               do n=1,ncp
                   j=icpm(m,n)
                   amat(i,j)=-stqm(m,n)             ! relative increment scheme
                   if(j.eq.nw) amat(i,j) = 0.d0     ! deriv w. resp to xh2o
               enddo
               satind=si2(mineq(i))
               if(mineq(i).gt.nmequ) satind=dlog10(si2k(mineq(i)-nmequ))
               bmat(i) = satind-sateq(i)  ! log form: bmat = log10(q/k)-sateq

               elseif (igaseq(i).ne.0) then
c
c             deriv of gas mass action w/respect to cp aq.
               m=igaseq(i)
               ncp=ncpg(m)
               do n=1,ncp
                   j=icpg(m,n)
                   amat(i,j)=-stqg(m,n)            ! relative increment
                   if(j.eq.nw) amat(i,j) = 0.d0    ! deriv w/ respect to xh2o
               enddo
               bmat(i) = dlog10(cg(igaseq(i)))-sateq(i)   !log form: bmat = log10(q/k/p)-sateq
             end if
c
c..........Fixed activity
           else if (icon(i).eq.3) then
             do j=1,npri
              amat(i,j)=0.d0
              if (i.eq.j) amat(i,j)=gamp(i)*cp(i)  !relative increment
             enddo
             bmat(i)=tt(i)-cp(i)*gamp(i)
C..........Charge balance if icon=4
           else if (icon(i).eq.4) then
             do j=1,npri
              amat(i,j)=cp(j)*z(j)    !simplify - only include derivative wrt primary species
              if (j.eq.nw) amat(i,j)=0.0d0
             enddo
c            add derivatives wrt secondary species
             do n=1,naqx            
               ncp=ncps(n)
               do m=1,ncp
                  k=icps(n,m)  !index of primary species k in sec. species n
                  if(k.ne.nw) then
                    dcs=stqs(n,m)*cs(n)*(z(npri+n))  !derivatives wrt secondary species
                    amat(i,k) = amat(i,k) + dcs
                  endif
               enddo
             enddo
c            calculates charge balance
             chbal=0.0d0
             do j=1,npri
               chbal=chbal+cp(j)*(z(j))
             enddo
             do j=1,naqx
               chbal=chbal+cs(j)*(z(npri+j))
             enddo
             bmat(i)=-chbal
c
           end if
c
500     continue
c
c       to ignore/skip speciation calcs of surface complexes
        if(npads.gt.0) then
          do i=1,npads
            n=npads+npaq
            amat(n,n)=1.d0
            bmat(n)=0.d0
          enddo
        end if
c
      else
c-------------------------------------------------------------------
c------initial speciation calcs including surface calcs
c-------------------------------------------------------------------
c
        do i = 1,npri
          do j=1,npri
             amat(i,j)=du2(i,j)
          end do
          bmat(i)=u2(i)-tt(i)
        enddo

        do m=1,nads
          ncp=ncpad(m)
          do n=1,ncp
            i=icpad(m,n)
            bmat(i)=bmat(i)+stqd(m,n)*cd(m)*xh2o
          enddo
        enddo
c
        do k=1,nads
          ncp=ncpad(k)
          do n=1,ncp
             i=icpad(k,n)
             do m=1,ncp
               j=icpad(k,m)
               amat(i,j)=amat(i,j)+stqd(k,n)*dcd(k,j)
             enddo
          enddo
        end do
c
        do i=1,npri
          bmat(i)=-bmat(i)   !negative sign because we add bmats later
        enddo
c
c-------Put surface complexation in the jabobian matrix
c       the alpha terms (potential term)
         if(npot.gt.0) then
           nmatc=npri+npot
           do i=npri+1,nmatc  ! deriv of potential equilibrium/respect to cp aq.
            kk=i-npri
            do j=1,npri
               amat(i,j)=dphi_dcp(kk,j)   ! relative increment scheme (mult by cp included in the derivative)
               if(j.eq.nw) amat(i,j) = 0.d0    ! deriv w/ respect to xh2o
            enddo
           enddo
c
           do j=npri+1,nmatc       ! deriv of mass bal w/respect to potential term
              kk=j-npri
              do i=1,npri
                amat(i,j)=dcp_dphi(i,kk)
              enddo
           enddo
c
           do j=npri+1,nmatc       ! deriv of potential equilibrium w/respect to potential term
              kj=j-npri
              do i=npri+1,nmatc
                 ki=i-npri
                 amat(i,j)=dphi_dphi(ki,kj)
              enddo
           enddo
c
c          the independent term (npot)
           do i=npri+1,nmatc
             kk=i-npri
             bmat(i) =-Fphi(kk)   ! negative sign, so we add bmat's later
           enddo
c
         endif
c
c         to equilibrate the surface keeping water composition unchanged
          if(isurfcalc.eq.2) then
           do i=1,npaq   ! not nmat
               do j=1,nmatc
                 if(j.ne.nw.and.i.ne.nw) amat(i,j)=0.d0
               enddo
               if(i.ne.nw) then   ! needs this for good water mass balance
                 bmat(i)=0.d0
c                any non-zero value for amat(i,i) should work since we set bmat(i)=0
                 amat(i,i) = 1.d0
               endif
           enddo
          endif
      endif
      return
      end
c
c-----------------------------------------------------------------------------
c

      subroutine cheminit(no_ch,isurfcalc)
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      integer*8 no_ch,isurfcalc,i1,i0
      common/ion_str1/str              !ionic strength
c
c     Note: repeats calls below for activity coefficient calcs
      i0 = 0
      do i=1,3
         call dh_hkf81(no_ch)
         call cs_cp
      end do
      i1 = 1
      call cmq_cp(i1)
      i0 = 0
      call cr_cp(i0)  !only to compute q/k for output
c
      if(isurfcalc.eq.0) return
c -------below only if surface complexation calcs, after first initialisation
c        with isurf=0
      if (nads.gt.0) call admodel
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine echotherm_HFK
c
c  Routine to echo all data read in the thermodynamic database
c  See routine readtherm for variable definitions
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      double precision coef(5)
      double precision stoic(mpri)
      data iunit2/42/
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ********echotherm 1.0, 2003.7.30: Echo all data read in'
     x' chemical thermodynamic database********')
c
      tc2=25.d0   ! for logK calculation and printout only
      tk2=tc2+273.15d0
      tkk = tk2
      tki = 1.d0/tkk
      tksqi = tki*tki
      tklog = dlog(tkk)
c
      write(iunit2,*) ' Components     a0    charge'
      do i=1,npri
        write(iunit2,"(1x,a10,2x,f8.3,2x,f5.1)") napri(i),a0(i),z(i)
      enddo
c-----Write out temperature interpolation coefficients
      write(iunit2,
     & "(/' Log K interpolation coefficients as a function', 
     &            ' of temperature (K) (a,b,c,d,e)')")
      write(iunit2,
     & "(' valid temperature range (deg.C.): ',f5.1,' to ',f5.1)")
     &  tmpmin,tmpmax
      write(iunit2,
     &  "(' ',15x,'        a*ln(TK)       b          c*TK'
     &  ,'     d*(TK)**-1  e*(TK)**-2')" )
      do i=1,naqx
        write(iunit2,1000) naaqx(i),(akcoes(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoes(i,n)
        enddo
        aks(i)=fak(coef,tkk,tklog,tki,tksqi)
      enddo
      do i=1,nmin
        write(iunit2,1000) namin(i),(akcoem(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoem(i,n)
        enddo
        akm(i)=fak(coef,tkk,tklog,tki,tksqi)
      enddo
      do i=1, ngas
        write(iunit2,1000) nagas(i),(akcoeg(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoeg(i,n)
        enddo
        akg(i)=fak(coef,tkk,tklog,tki,tksqi)
      enddo
      do i=1,nads
        write(iunit2,1000) naads(i),(akcoead(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoead(i,n)
        enddo
        akd(i)=fak(coef,tkk,tklog,tki,tksqi)
      enddo
c     Write out new dV data read from database if present
c     Write out pressure interpolation coefficients
      do i = 1,naqx
        if(akcops(i,1).ne.0.d0) then
         if(i.eq.1) write(iunit2,1003) tmpmin,tmpmax
         write(iunit2,1000) naaqx(i),(akcops(i,j),j=1,5)
        end if
      enddo
      do i = 1,nmin
        if(akcopm(i,1).ne.0.d0) 
     &  write(iunit2,1000) namin(i),(akcoem(i,j),j=1,5)
      enddo
      do i = 1,ngas
        if(akcopg(i,1).ne.0.d0) 
     &   write(iunit2,1000) nagas(i),(akcoeg(i,j),j=1,5)
      enddo
 1000 format(' ',a20,5(1pe12.4))
 1003 format(/'dV interpolation coefficients as a function of'
     &' temperature (K) (a,b,c,d,e) for log(K) pressure dependency'/
     &' valid temperature range (deg.C.): ',f5.1,' to ',f5.1,/24x,'a'
     &' '13x,'b*TK'7x,'c*TK**2     d*(TK)**-1  e*(TK)**-2')
c     Write stoichiometries and logK's at T=temp (defined above)
      if(naqx.gt.0) then
          write(iunit2,"(/' Derived Species Reactions')")
          write(iunit2,590) tc2, (napri(j),j=1,npri)
          do i=1,naqx
            ncp=ncps(i)
            do j=1,npri
              stoic(j)=0.d0
              do n=1,ncp
                if(icps(i,n).eq.j) stoic(j)=stqs(i,n)
              enddo
            enddo
            write(iunit2,600) naaqx(i),a0(npri+i),z(npri+i),aks(i),
     &           (stoic(j),j=1,npri)
          enddo
      endif
c
      if (nmin.gt.0) then
          write(iunit2,"(/' Mineral Reactions')")
          write(iunit2,595) tc2, (napri(j),j=1,npri)
          iflgstop = 0
          do i = 1, nmin
            ncp=ncpm(i)
            do j=1,npri
              stoic(j)=0.d0
              do n=1,ncp
                if(icpm(i,n).eq.j) stoic(j)=stqm(i,n)
              enddo
            enddo
            write(iunit2,605) namin(i), vmin(i), akm(i),
     &         (stoic(j),j=1,npri)
          enddo
      endif
c
      if (ngas .gt. 0) then
          write(iunit2,"(/' Gas Reactions')")
          write(iunit2,596) tc2,(napri(j),j=1,npri)
          do i = 1, ngas
            ncp=ncpg(i)
            do j=1,npri
              stoic(j)=0.d0
              do n=1,ncp
                if(icpg(i,n).eq.j) stoic(j)=stqg(i,n)
              enddo
            enddo
            write(iunit2,606) nagas(i),akg(i),
     &         (stoic(j),j=1,npri)
          enddo
      endif
c
      if (nads .gt. 0) then
          write(iunit2,"(/' Surface Adsorption Reactions')")
          write(iunit2,597) tc2,(napri(j),j=1,npri)
          do i = 1, nads
            ncp=ncpad(i)
            do j=1,npri
              stoic(j)=0.d0
              do n=1,ncp
                if(icpad(i,n).eq.j) stoic(j)=stqd(i,n)
              enddo
            enddo
            write(iunit2,610) naads(i),zd(i), akd(i),
     &        (stoic(j),j=1,npri)
          enddo
      endif
c
  590 format(2x,'species',8x,'a0',3x,'charge',2x,'logK(',f3.0,' C)',
     &   2x,50a7)
  595 format(2x,'minerals',4x,'m.vol(L/mol)',2x,'logK(',f3.0,' C)',
     &   2x,50a7)
  596 format(2x,'gases   ',18x,'logK(',f4.0,'C)',2x,50a7)
  597 format(2x,'surf.cmplex',15x,'logK(',f4.0,'C)',2x,50a7)
  600 format(' ',a12,f8.3,1x,f5.1,3x,f10.3,50(f7.2))
  605 format(' ',a12,2x,f10.4,2x,f10.3,1x,50(f7.2))
  606 format(' ',a12,14x,f10.3,1x,50(f7.2))
  610 format(' ',a12,7x,f5.1,3x,f10.3,50(f7.2))
      return
      end
c
c-------------------------------------------------------------------------------
c
c
c*******************************************************************
      SUBROUTINE surfequil

c********* initial speciation of chemical reactions for each node,
c          including equilibration with surface
c*******************************************************************
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      include 'perm_v2.inc'
!
!$$$$$$$$$ COMMON BLOCKS FOR ELEMENTS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/E6/T(MNEL)
      common/TEM_EOS9/Tc_EOS9(MNEL)  ! initial temperature (oC)
      COMMON/SOLUTE8/SL1(MNEL)           ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)           ! new gas saturation
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
      common/irtlaw/nplaw(mmin)
cels5/6/08      common/satgas2/sg2  ! current node gas saturation
cels10/5/09      common/satgas2/sg2,tmp2  ! Current node gas saturation and temp ns98/3 added temp
      common/satgas2/sg2
      common/afactorr/a_fmr(mnel)
      COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
cels10/5/09      double precision densw,deltat0
      double precision densw
      common/aqkin16/NoTrans(mpri)    ! >0: not subject to transport
      COMMON/tot_solid_aq/icon_nod(mnod,mpri),ttt(mpri),
     &                ttt_nod(mnod,mpri)          ! total concentraion including both aqueous and solid
cels6/24/08 need this common block to pass timetot correctly
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
cels10/5/09
      common/dtlim/max_chem_it,delt_conne,id_chem
      character*16 delt_conne
      character*5 id_chem
c
      integer*8 ielem

C-----------------------------------------------------------------------------
      IF (IEOS .EQ. 9)   THEN
         DO N=1,NEL
            T(N)=Tc_EOS9(N)
         END DO
      END IF
c
      DO N=1,NNOD
         TC(N)=T(N)
         tkelv(n) = tc(n) + 273.15d0
      END DO

      max_chem_it = 0    !ns 9/06 max chemical iterations for entire printout interval
c
      ITERCH=0
      MAXITCH=0        ! Maximum iterations of solving whole chemistry
      AVERITCH=0.d0
      COUNTCH=0.d0
c
c---------Start of geochemical calculation loop for each node (loop 1000) ------
c         ***************************************************

       DO 1000 I=1,NNOD

c      get water density for each node in kg/m3
c      Watch!! Dry nodes have zero densities ==> concentrations will become zero
        NLOC2=(I-1)*NSEC*NEQ1
        if (IEOS .NE. 9)   then     ! for other modules
           NLOC2L=NLOC2+NBK
           densw = PAR(NLOC2L+4)
        else                      ! for eos9 module
           densw=PAR(NLOC2+4)
        endif
        if(densw .gt. 2.d3) densw=1.d3     !???
        if(densw.eq.0.d0) densw=1.d3

        dliq=densw/1.d3  !density in kg/l


c--- Calculates total site concentrations for surface complexes
c       get the surface area (surfads) in m2/kg_h2o for current node
        call ads_area (i,mnel,nmin,nsurf,m_index,pre,phi,
     &            sl1,a_fmr,densw,vmin,surfads,supadn)

        call surface_conc(i,densw)

        sg2 = sg1(i)
cels4_28_08_should be defined already        sl1 = 1.d0-sg1(i)
        tc2 = tc(i)                ! added for use in NEWTONEQ
        tk2 = tkelv(i)

        call assign
c
c     assign the chemical parameters of the current node I
        naqt=npri+naqx
c     c is concentration of primary species in moles/l liquid
c     ut is total concentration (primary+secondary species) in moles/l liquid
c     tt is total moles in solution (for given volume of liquid)
        stion=0.d0   !stoichiometric ionic strength
        do n=1,npri
           tt(n)=ut(i,n)  !tt in moles (per liter liq)
           stion=stion+zsqi(n)*ut(i,n)
        enddo
c
c   adsorption contribution to total concentrations
c   (not really needed here at initialization because d(i,j)=0)
        do m=1,nads
          ncp=ncpad(m)
          do k=1,ncp
            n=icpad(m,k)
c           skips total for primary surface species - we already have these
c           the total for these species comes from the surface area & site density data
            if(n.le.npaq)
     &         tt(n)=tt(n)+stqd(m,k)*d(i,m)*cp(nw)
          enddo
        enddo

c----if icon(i) equal to 5, the total (both solid and aqueous phase) is known, so
c---- tt is replaced by ttt,
         do n=1,npri
          if(icon_nod(i,n).eq.5)tt(n) = ttt_nod(i,n)
         enddo

c
c---- Assigns old concentrations at current node as initial guess values
c     cp,cs,cm etc are in moles/kg water liq. c's are in moles/Kg h2o
        do n=1,npaq    !npri   exclude surf complexes
          cp(n)=c(i,n)
        end do
        do n=1,naqx
          cs(n)=c(i,npri+n)
        end do
c   for water, cp(nw) is kilogram water
        tt(nw)=rmh2o*dliq
        cp(nw)=dliq
c
c   adsorption contribution to total concentrations
c   (not really needed here at initialization because d(i,j)=0)
        do m=1,nads
          ncp=ncpad(m)
          do k=1,ncp
            n=icpad(m,k)
c           skips total for primary surface species - we already have these
c           the total for these species comes from the surface area & site density data
            if(n.le.npaq)
     &         tt(n)=tt(n)+stqd(m,k)*d(i,m)*cp(nw)
          enddo
        enddo

c----if icon(i) equal to 5, the total (both solid and aqueous phase) is known, so
c---- tt is replaced by ttt,
         do n=1,npri
          if(icon_nod(i,n).eq.5)tt(n) = ttt_nod(i,n)
         enddo
c
c        do m=1,nmin
c cm is moles mineral (per liter liquid); pre in moles per liter medium
c now cm is  incremental change
c          cm(m)=0.d0
c        end do
c
cels10/6/09 not used in this routine        rt = gc*tk2              !added (gas RT)
c
!------------------------------------------------------------------------------
!
        do k=1,nads
cels6/10/08 Nic revised dimensioning of cp 
cels6/10/08           cd(k)=d(i,k)    !d is in mol/kgw
cels6/10/08           if(cd(k).eq.0.d0) cd(k) = cp(npaq+k)/1.d-3
           if(d(i,k).eq.0.d0)then
              ik = iad_neu(k)
              cd(k) = cp(ik)/1.d-3
           else
              cd(k)=d(i,k)
           endif
        enddo

        phi2=phi(i)

C------------------------------------------------
C---    Call Newton-Raphson subroutine for each node
c-------------------------------------------------
        ielem = i
C       Flag isurf=1: inclusion of surface complexes, no surface equilibration
c            isurf=2: inclusion of surface complexes, with initial surface equilibration
        isurfcalc=1
        if(isurfeq(i).ne.0) isurfcalc=2

        CALL NRINIT(isurfcalc)

        if(ichdump.eq.1) call chdump(timetot,ielem,iterch)  !for testing only
C
        if (ichdump.eq.2) then  ! Printing speciation for specified grid blocks
            do ino=1,nwnod
              nng=iwnod(ino)
              if (nng.eq.ielem) then
                 call chdump(timetot,ielem,iterch)
                 go to   3099
              end if
            end do
3099    continue
        end if

c     case where we reached the maximum iterations (no convergence)
      if(iterch.ge.maxitpch) then
        WRITE(20,"
     &   (/2X,'ERROR: convergence problem in initial surface ',
     &   'complexation calcs', /2X,'   Please check/revise input',
     &   ' parameters',/)")
        WRITE(32,"
     &   (/2X,'ERROR: convergence problem in initial surface ',
     &   'complexation calcs', /2X,'   Please check/revise input',
     &   ' parameters',/)")
        stop
      endif
C
c     Now we need to convert concentrations from chem module:
c     cp's are in moles/kg h2o (molal) (per original 1 liter liquid fed into chem. module)
c     cm, cmg, rkin etc are in moles (per original 1 liter liquid fed into chem. module)
c     cp(nw) contains kg of water liquid (per original 1 liter liquid fed into chem. module)
c     vliq is new volume of liquid from the original 1 liter we used to run the
c     chemical module.  For now leave to 1 until we implement a way to
c     add/remove water in tough if vliq changes from 1
c
      sumsalts=0.d0   !sum of salts weights in kg per kg water (assume zero for now)
      vliq=1.d0       !1 liter initial
      factw=cp(nw)/vliq  !conversion factor = kg h2o liq/liter liquid
c------------------------------------------------------------
c   Assign the new concentration values to the current node
c------------------------------------------------------------
c NS4/06 change of units - save c's as molalities
c   cp's, cs's are in moles/kg h2o liq
        do n=1,npri
          c(i,n)=cp(n)
        enddo
c       Note, here we save the actual moles of H2O (per Kg H2O)
        c(i,nw)=rmh2o
c
        do n=1,naqx
          c(i,npri+n)=cs(n)
        enddo
        ph(i)=ph2
c
c
        do k=1,nads
          d(i,k)=cd(k)
        enddo
        do n=1,nsurf
        phip(i,n)= phip2(n)
       enddo
!
!.......Multi-site exchanges
!
c ns-lz 3/09  do isite=1,NXsites     ! Loop over multi-sites
c ns-lz 3/09  do k=1,nexc
c ns-lz 3/09      xcads(i, isite, k) = cxM(isite,k)/vliq
c ns-lz 3/09   end do
c ns-lz 3/09  end do
!
!    Compute new ut's (total aqueous concentrations in mol/L soln)
!
         do n=1,npaq       !npri skip surf complexes
            ctot(i,n)=c(i,n)   !ctot in mol/kgw
         enddo
         do j=1,naqx
            ncp=ncps(j)
            do k=1,ncp
              n=icps(j,k)
              utemp=ctot(i,n)+stqs(j,k)*c(i,npri+j)
              ctot(i,n)=utemp        !total concentrations in mol/kg water
            enddo
         enddo
!........Need loop below to compute ut's in mol/L
         do n=1,npri
c             If we equilibrated the surface, the total amount of water increased to
c              account for sorbed water.  For this reason, we normalize all concentrations
c              (and water amount) such that the total amount of water is consistent
c              with the initial water density (for consistency with future calcs in couple).
c use ut below              if(isurfeq(i).eq.1) ctot(i,n)=ctot(i,n)*dliq/factw
              ut(i,n)=ctot(i,n)*factw  !total concentrations in mol/L sln (for transport equations)
                if(isurfeq(i).eq.1) ut(i,n)=ut(i,n)*dliq/factw
              utold(i,n)=ut(i,n)
         enddo

         isurfeq(i) = 0    !reset to 0 for any further calcs




1000  CONTINUE

c
c-------------end of geochemical loop for each node---------------------
      RETURN
      END

      subroutine ads_area (ig,mnel,nmin,nsurf,m_index,pre,phi,
     &  sl1,a_fmr,densw,vmin,surfads,supadn)
c
c... Calculates surface area for sorption in m^2/kg water - NS 2/08
c    Written after routine rsfarea from ELS 6/17/01
c
c  input:
c     ig       index of current grid block
c     mnel     local dimension, number of gridblocks
c     nmin     local dimension, number of minerals
c     nsurf    local dimension, number of surfaces
c     m_index  pointer to mineral index in list of minerals
c     pre      mineral amount in mol_min/dm3_medium
c     phi      porosity
c     sl1      liquid saturation
c     a_fmr    active-fracture parameter
c     densw    water density in kg/m3
c     vmin     initial input mineral molar volume converted to dm3_mineral/mol_mineral
c     supadn   initial input mineral surface area (per node, converted to m2_mineral/m3_mineral)
c  output:
c     surfads  surface area (per mineral, for current node ig, in m2_mineral/kg_h2o)

      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      integer*8 m,ig,mnel,nmin,nsurf
      integer*8 m_index(nsurf)
      double precision pre(mnel,nmin)
      double precision supadn(mnel,nsurf)
      double precision surfads(nsurf),vmin(nmin)
      double precision phi(mnel),sl1(mnel),a_fmr(mnel)
      double precision densw
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ads_area 1.0, 2008.2.28: Calculate reactive surface area'
     x' for sorption')
c
c...........Modified to account for reduced surface area: Factor of
c...          S from Liu et al. (1998) causes S in denominator
c...          to drop out (saturated system, also)
         if(a_fmr(ig).lt.sl1(ig).and.sl1(ig).gt.0.d0.and.
     +      a_fmr(ig).gt.0.d0)then
c..... Factor based on active fracture model at low saturations,
c....... and for saturations above zero
           actfrc = a_fmr(ig)/(phi(ig)*densw*sl1(ig))
         else
c..... Factor for saturated system, or unsaturated to consider only
c........ the wetted proportion
           actfrc = 1.d0/(phi(ig)*densw)
         endif
c
c... Volume fraction of rock
         vfrtmp = 1.d0-phi(ig)
c
        do m = 1, nsurf
             m_min=m_index(m)
             if(m_min.eq.0) then
c              surface is not tied to a mineral
             surfads(m) = supadn(ig,m)*vfrtmp*actfrc
           elseif(m_min.lt.0) then
c              we entered a constant surface area already in m2/kgw
               surfads(m) = supadn(ig,m)
             elseif(m_min.gt.0) then
c              Assign mol/L medium based on mineral amount
               rkmol = pre(ig,m_min)
c              Calculate area in m^2/kgw
               surfads(m) = supadn(ig,m)*rkmol*vmin(m_min)*actfrc
             endif
        enddo
      return
      end



c*************************************************************************
      SUBROUTINE surface_conc(i,densw)
c*************************************************************************
c---  Calculates total site concentrations for surface complexes
c     and set various pointers for potentiual term calculations
c
c       i        current node index
c       densw    water density in kg/m3
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      include 'perm_v2.inc'
      double precision densw

      dliq=densw/1.d3    !water density in kg/L
      do n=1,npads
c         total site concentration in mol_sites/kg_h2o is stored in ctot
c         total site concentration in mol_sites/L_solution stored in ut and utold
          kk=npaq+n
          jj=isurfp(n)
          ctot(i,kk)  = surfads(jj)*site_dens(n)

c         if no surface mineral is present (yet), concentrations will be
c         zero so we reset to a small number
          if(ctot(i,kk).eq.0.d0) ctot(i,kk) = 1.d-20
          utold(i,kk) = ctot(i,kk)*dliq
          ut(i,kk)=utold(i,kk)
      enddo

c     counts the number of surfaces for which we need to solve potential
c     and creates pointers (so we can simulatneously run surfaces with and without
c     potential terms).
c     When a mineral is exhausted (pre < 1.d-10), we remove potential calculations to
c     avoid problems but keep the surface in, with small total concentration =1.d-20
      do j=1,nads
         iad_phi(j)=0   !points only to surfaces requiring potential
      enddo
      npot=0         !counts number of surfaces for which to run potential calcs
      do j=1,nsurf
        ipoten(j)=0
        phip2(j)=0.d0
        imodel=iadmod(j)   !adsorption model type 0,1,2, or 3
        kk=m_index(j)
        if(kk.gt.0.and.pre(i,kk).lt.1.d-10) imodel=0  !reverts to model 0 if mineral exhausted
        if(imodel.gt.0) then
          npot=npot+1
          ipoten(npot)=j   !points to index of surface array for each potential
          phip2(j)=phip(i,j)
          do k=1,nads
            if(iad_surf(k).eq.j)
     &            iad_phi(k)=npot  !points to index of array of non-zero potentials
          enddo
          if(phip2(j).eq.0.d0) then
            phip2(j)=-0.1d0
          else
            phip2(j)=phip(i,j)
          endif
        endif
      enddo


c     trial values for primary surface complexes
      do k=1,npads
         kk=npaq+k
         if(ctot(i,kk).gt.1.d-20) then
           cp(kk)=c(i,kk)
         else
           cp(kk)=1.d-20
         endif
      enddo

c Not needed, recomputed in cd_cp from the primary surface species
c      do k=1,nads
c         cd(k)=d(i,k)     !both cd and d are in mol/kgw
c      enddo

      return
      end
c
c
!-----------------------------------------------------------------------------
!
      SUBROUTINE fugacomp(inode)
c     Moved from couple, with minor (non-effective) changes 
c******************** Prepare calculations of gas fugacity coefficients *******************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
      INCLUDE 'common_v2.inc'
      COMMON/EOS_INDICATOR/ IEOS           ! indicate EOS module used
      common/fuga_coe /FugCoeCO2(mnel)     ! fugacity coefficients from eco2n
      common/fugacity_coe/fug_coe(mnel,18) ! fugacity coefficients from tmgas
      common/co2_gene1/     nco2g
      common/h2_gene1/      nh2g
      common/OtherGases_ix/ nch4g, nh2sg, nso2g
      common/gas_index/ichem(18)           ! no-cond gas index in chemical input
      double precision gamg_ig
      integer*8 inode,ig
c
       do ig=1,ngas
c
           gamg(ig) = 1.0d0   ! Gas fugacity coefficients

c..........For all modules upon initilisation
          if(inode.eq.0) then
              CALL Gas_Fuga_Coe (Pt, tk2, ig, gamg_ig)
              gamg(ig)      = gamg_ig
c
c           split loops
c
            elseif(inode.gt.0)then
c
c..........For ECO2N  module
c
              if (ieos .eq. 14) then
c
                if (ig .eq. nco2g)    then
                  gamg(ig) = FugCoeCO2(inode)   ! FugCoeCO2 is from ECO2N
                end if
c
c..........For EOS_TMgas module
c
              else if (ieos.eq.15.or.ieos.eq.16)  then
c
                iFugaC_React = 0
                icg = ichem(ig)      ! no-cond gas Index in chemical input file
                gamg(ig) = fug_coe(inode,icg)
!                        fug_coe(i,icg) is from EOS_TMgas in different gas order
                if (gamg(ig) .le. 1.0d-10)  gamg(ig) = 1.0d0
!
c
c..........For other fluid flow modules
c
              else
c
c...............Call subroutine to obtain gas fugacity coefficients
c
                CALL Gas_Fuga_Coe (Pt, tk2, ig, gamg_ig)
c
                gamg(ig)      = gamg_ig
                fug_coe(inode,ig) = gamg(ig)
                if (ig .eq. nco2g)   then          !!!! Not needed
                   fugCoeCO2(inode)   = gamg(ig)   ! Same order in chemical input file
                end if
c
              end if
c
          endif
c
       end do
c
      return
      end
!
!
      SUBROUTINE Gas_Fuga_Coe (Pt, Ta, ig, gamg_ig)
!
!******************** Calculate gas fugacity coefficients *******************
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!----------
!.....Gas index for gas fagacity calculations in REACT
!----------
      common/co2_gene1/     nco2g
      common/h2_gene1/      nh2g
      common/OtherGases_ix/ nch4g, nh2sg, nso2g
      double precision gamg_ig,pt,ta,aaa,bbb,ccc,ddd,eee,fff
      integer*8 ig
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' fuga_coe_CH4 1.0, 2003.7.30: Calculate gas fugacity'
     &' coefficients')
!
!------------------------------------------Assume real gases and ideal mixture
!---------------------------CH4 fugacity coefficient is calculated from eq.(14)
!----------------------------------of the paper of Nicolas Spycher et al., 1988
!-------------------For 16-350 0C and 0-500 bar
!
!----------
!.....Default fugacity coefficient
!----------
!
      gamg_ig = 1.0d0
!
      iFugaC = 0        ! Indicator for gas fugacity corrections in REACT
!
      if (ieos.eq.2  .or. ieos.eq.13)  iFugaC = 1
      if (ieos.eq.1)                   iFugaC = 1  
      if (ieos.eq.3)                   iFugaC = 1  
      if (ieos.eq.4)                   iFugaC = 1  
      if (ieos.eq.5)                   iFugaC = 1  
      if (ieos.eq.7)                   iFugaC = 1  
      if (ieos.eq.9)                   iFugaC = 1  
      if (ieos.eq.14)                  iFugaC = 1  ! if eco2n, comes here only on chemistry initialisation
!
!--------------------------------------------------------------------------------
      if (iFugaC .eq. 1)       then
!--------------------------------------------------------------------------------
!
!----------
!........CO2 gas fugacity
!----------
!
         if (ig .eq. nco2g)    then
!
            aaa = -1430.87d+0                ! Parameters for CO2 gas
            bbb =  3.598d+0
            ccc = -227.376d-5
            ddd =  347.644d-2
            eee = -1042.47d-5
            fff =  846.271d-8
!
            CALL GasFugFunction(Pt, Ta, gamg_ig,
     &                        aaa, bbb, ccc, ddd, eee, fff)
!
         end if  ! IF_CO2g
!
!----------
!........CH4 gas fugacity
!----------
!
         if (ig .eq. nch4g)    then
!
            aaa = -537.779d+0          ! Parameters for CH4 gas
            bbb =  1.54946d+0
            ccc = -92.7827d-5
            ddd =  120.861d-2
            eee = -370.814d-5
            fff =  333.804d-8
!
            CALL GasFugFunction(Pt, Ta, gamg_ig,
     &                        aaa, bbb, ccc, ddd, eee, fff)
!
         end if  ! IF_CH4g
!
!----------
!........H2  gas fugacity
!----------
!
         if (ig .eq. nh2g)    then
!
            aaa=-12.5908d+0             ! Parameters for H2 gas
            bbb= 0.259789d+0
            ccc=-7.24730d-5
            ddd= 0.471947d-2
            eee=-2.69962d-5
            fff= 2.15622d-8
!
            CALL GasFugFunction(Pt, Ta, gamg_ig,
     &                        aaa, bbb, ccc, ddd, eee, fff)
!
         end if  ! IF_H2g
!
!----------
!........H2S gas fugacity
!----------
!
         if (ig .eq. nh2sg)    then
!
            aaa = -1430.87d+0            ! Parameters for CO2 gas  !!!!
            bbb =  3.598d+0
            ccc = -227.376d-5
            ddd =  347.644d-2
            eee = -1042.47d-5
            fff =  846.271d-8
!
            CALL GasFugFunction(Pt, Ta, gamg_ig,
     &                        aaa, bbb, ccc, ddd, eee, fff)
!
         end if  ! IF_H2Sg
!
!----------
!........SO2 gas fugacity
!----------
!
         if (ig .eq. nso2g)    then
!
            CALL GasFugFunction(Pt, Ta, gamg_ig,
     &                        aaa, bbb, ccc, ddd, eee, fff)
!
         end if  ! IF_SO2g
!
!--------------------------------------------------------------------------------
      end if  !  IF_iFugaC
!--------------------------------------------------------------------------------
!
!
      return
!
      end
!
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE GasFugFunction(Pt, Ta, gamg_ig,
     &                        aaa, bbb, ccc, ddd, eee, fff)
!
!******************** Calculate CO2 gas fugacity coefficients *******************
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      double precision gamg_ig,pt,ta,aaa,bbb,ccc,ddd,eee,fff
!
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' fuga_coe_CO2 1.0, 2003.7.30: Calculate gas fugacity'
     &' coefficients')
!
!------------------------------------------Assume real gases and ideal mixture
!---------------------------CO2 fugacity coefficient is calculated from eq.(14)
!----------------------------------of the paper of Nicolas Spycher et al., 1988
!-------------------For 50-350 0C and 0-500 bar
!
      gamg_ig = (aaa/(Ta*Ta) + bbb/Ta + ccc)*Pt
      gamg_ig = gamg_ig + (ddd/(Ta*Ta) + eee/Ta + fff)*Pt*Pt/2.0d0
      gamg_ig = dexp(gamg_ig)
!
      if (gamg_ig .le. 0.0D0 .or. gamg_ig .gt. 1.0D0)   then
         write (32,*) 'There is a problem with fugacity coefficient'
         stop
      end if
      return
      end
!
!
!-------------------------------------------------------------------------------
!
!
!
      SUBROUTINE STEADYC(NNOD,NMEQU,NMKIN,TIMETOT)
C
C****** Update mineral phase when chemical quasi-steady state is reached *******
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
!
      common/chemgrid/c(mnod,maqt),utold(mnod,mpri),ut(mnod,mpri),
     & rhand(mnod,mpri),rsource(mnod,mpri),ph(mnod),gP(mnod,mgas),
     & aream(mnod,mmin),sads(mnod),psi(mnod),
     & d(mnod,mads),supadn(mnod,msurf),phip(mnod,msurf),
     & surfads(msurf),ub(mbound,mpri),ctot(mnod,mpri),cnfact
!
      common/vmineral/pre(mnod,mmin),pre0(mnod,mmin),
     +  pinit(mnod,mmin+1)
c
      COMMON /DISRATE/ DRATE(MNOD,MMIN)
C---------------------------Common blocks for change in porosity
      COMMON/CPOROSITY/DELPHI(NMNOD)   ! change porosity from time 0
      COMMON/CPOROSOLD/DPHIOLD(NMNOD)  ! last time step porosity change
      COMMON/SOLUTE8/SL1(MNEL)         ! new liquid saturation
      COMMON/E3/EVOL(MNEL)
      COMMON/BMNO/NBLOCK,NMINERAL ! Block and mineral numbers where mineral exhausted
      COMMON/DFM/TIMAX,REDLT
      COMMON/TIMESTEA/TIMESTEA
      common/constraints/sl1min,stimax,dlstmx
C
      integer*8 nnod,nmequ,nmkin
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' STEADYC 1.0, 2003.7.30: Update minerals when chemical'
     x' quasi-steady state is reached')
c
C----------------------Minimam time needed let one mineral or more exhausted
c
      TIMESTEA=1.0D+20
      DO I=1,NNOD
       IF(SL1(I).LE.sl1min.OR. EVOL(I).LE.1.0D-15) GO TO 949
       DO M=1,NMEQU
          RATE1=DRATE(I,M)
              IF (RATE1.GT.1.0D-20 .AND. PRE(I,M).GT.0.0D0) THEN
                TIMES=PRE(I,M)/RATE1
                 IF(TIMES.LT.TIMESTEA) THEN
                    TIMESTEA=TIMES
                    NBLOCK=I         ! Block number
                    NMINERAL=M       ! Mineral number
                 END IF
              END IF
       END DO
       DO M=1,NMKIN
          RATE2=DRATE(I,NMEQU+M)
              IF (RATE2.GT.1.0D-20 .AND. PRE(I,NMEQU+M).GT.0.0D0) THEN
                TIMES=PRE(I,NMEQU+M)/RATE2
                 IF(TIMES.LT.TIMESTEA) THEN
                    TIMESTEA=TIMES
                    NBLOCK=I             ! Block number
                    NMINERAL=M+NMEQU     ! Mineral number
                 END IF
              END IF
       END DO
c
949   CONTINUE
      END DO
c
c---------------Extrapolation can not over maximum simulation time (TIMAX)
      DTIMAX=TIMAX-TIMETOT
      IF (TIMESTEA .GT. DTIMAX)  TIMESTEA=DTIMAX
c-------------------------------------------------------------------------------
      TIMETOT=TIMETOT+TIMESTEA
c
C-------------------------------------------------Update mineral abundance
c
        ntmin = NMEQU+NMKIN
      DO I=1,NNOD
       IF(SL1(I).LE.sl1min.OR. EVOL(I).LE.1.0D-15) GO TO 959
c       DO M=1,NMEQU+NMKIN
       DO M=1,ntmin
            RATE2=DRATE(I,M)
            PRE(I,M)=PRE(I,M)-RATE2*TIMESTEA
       END DO
959   CONTINUE
      END DO
c
C--------------------------------------------------------------------------
C
      RETURN
      END
c

