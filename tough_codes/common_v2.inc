C------------------------ common_v2.inc ----------------------
c ------------- TOUGHREACT V2.0 COMMON BLOCKS -------------------
C
cels7/1/08      common/icontrol/iswitch,njacob,ichdump
      common/icontrol/iswitch,njacob,ichdump,ngamm
      common/rcontrol/rcour
      common/control2/ispia,inibound,kcpl
      common/names/napri(mpri),naaqx(maqx),naaqt(maqt),namin(mmin),
     +   nagas(mgas),naexc(mexc),naads(mads)
      common/solute11/npri,npaq,npads
      common/ichemcons/naqx,nmequ,nmin,ngas,ne,nw,nh,
     +   noh,no2aq,nd,nads,nexc,nx,iex
cels6/26/08 added azero
      common/chemcons/a0(maqt),z(maqt),zd(mads),zz2(maqt),azero(maqt)
cels5/5/06 added square and absolute value of charge
cels6/26/08 added azero cubed, omega, and zabterm (term in dh_hkf81)
      common/chemcons2/zsqi(maqt),zabsi(maqt),azero3(maqt),
     +  omegahkf(maqt),zabterm(maqt)
cels6/26/08 added salting out index
      common/iso_ind/indx_so(8),nis
      common/minprop/vmin(mmin),dmolwm(mmin)
      common/gasprop/dmwgas(mgas),diamol(mgas)
cels10/12/09 no ipcal used anywhere      common/ichemnode/ncps(maqx),
cels10/12/09     +   icps(maqx,mpri),ncpm(mmin),icpm(mmin,mpri),
cels10/12/09     +   ncpg(mgas),icpg(mgas,mpri),kineq(mmin),ipcal       !ns2/7/02 ns5/19/09
      common/ichemnode/ncps(maqx),icps(maqx,mpri),ncpm(mmin),
     +   icpm(mmin,mpri),ncpg(mgas),icpg(mgas,mpri),kineq(mmin)
      common/chemnode/cp(mpri),tt(mpri),
     +   cs(maqx),u2(mpri),zp(mpri),zs(maqx),gamp(mpri),gams(maqx),
     +   ph2,cm(mmin),cmg(mgas),cg(mgas),si2(mmin),
     +   sig2(mgas),stqs(maqx,mpri),stqm(mmin,mpri),stqg(mgas,mpri),
     +   aks(maqx),akm(mmin),akg(mgas),xh2o,
     +   akcoes(maqx,5),akcoem(mmin,5),akcoeg(mgas,5),
     +   si2_old(mmin),akcops(maqx,5),akcopm(mmin,5),akcopg(mgas,5)
      common/dchemnode/du2(mpri,mpri),dsi2(mmin,mpri),
     +   dsig2(mgas,mpri),dgamt(maqt,mpri),dgamp(mpri,mpri),
     +   dgams(maqx,mpri)
      common/iadsnode/ncpad(mads),icpad(mads,mpri)

      COMMON/FUGACOE/GAMG(MGAS)    ! gas fugacity coefficient

      common/adsnode/cd(mads),dcd(mads,mpri),
     + stqd(mads,mpri),akd(mads),phip2(mads),akcoead(mads,5),
     + site_dens(mpads),capacitance(msurf)

      common/surface1/naads_min(mads),naads_pri(mads),namesurface(mads)
      common/surface3/dcp_dphi(mpri,mads),dphi_dcp(mads,mpri)
     +                  ,dphi_dphi(mads,mads),Fphi(mads) 

      common/iadpoint/isurfp(mpads),iad_phi(mads),m_index(msurf),
     + nsurf,ipoten(msurf),npot,
     + iadmod(msurf),isurfeq(mnod),iad_surf(mads),iad_neu(mads)

      common/exchnode/cx(mexc),dcx(mexc,mpri),cec2,
     +     stqx(mexc,mpri),ekx(mexc)
      common/iexchnode/nbx(mexc)
!
      common/exchange1/phi2,ct(mtot),gamt(mtot)  
      common/tempvars/tc2,tk2,tmpmin,tmpmax,Pt,p0bar    !ns5/09 add Pt,p0bar pressure bar 
      common/vmineral/pre(mnod,mmin),pre0(mnod,mmin),
     +  pinit(mnod,mmin+1)
      common/imineral/noreact(mnod,mmin)   !ns09/09
      common/kincons/ikin(mmin),idep(mmin),ck1(mmin),ck2(mmin),
     +     rkf(mmin),ea(mmin),thresh(mmin),deph(mmin),depoh(mmin),
     +     ck1h(mmin),ck2h(mmin),rkh(mmin),eah(mmin),
     +     ck1oh(mmin),ck2oh(mmin),rkoh(mmin),eaoh(mmin),
     +     cr(mpri),dcm(mmin,mpri)


      common/kincons1/nmkin
      common/kincons2/paimk(mmin),si2k(mmin),akin(mmin),
     &      rkin2(mmin), amin2(mmin), rad2(mmin), 
     &      amin(mnod,mmin),rad(mnod,mmin),grains(mnod,mmin),
     &      sumqk(mmin)   !ns5/09 add this
      common/preckin/rkprec(mmin),acfprec(mmin),bcfprec(mmin),
     +     ccfprec(mmin),eaprec(mmin),ck1prec(mmin),ck2prec(mmin),
     +     rnucl(mmin),ssqk(mmin),ssq(mmin),
     +     ssqk0(mmin), sstk1(mmin), sstk2(mmin),
     +     ssq0(mmin), sst1(mmin), sst2(mmin)


!
!.....Multiple exchangeable sites
!
      common/exchange_sl/ Mod_Xsl                       ! Model for exchange dependance on water saturation 
      common/MultiXSite1/ NXsites                       ! Number of exchangeable sites
      common/MultiXSite2/ cecM(MXsites)                 ! Batch CEC for multi-sites
      common/MultiXSite3/ ekxM(MXsites, mexc)           ! Selectivity
      common/MultiXSite4/ cxM (MXsites, mexc)           ! Exchanged species concentrations
      common/MultiXSite5/ dcxM(MXsites, mexc, mpri)     ! Devivatives
      common/MultiXSite6/ cec (mnod, MXsites)           ! CEC for grid blocks
      common/MultiXSite7/ xcads(mnod, MXsites, mexc)    ! Exchanged species concentrations for grid blocks
!
      common/chemgrid/c(mnod,maqt),utold(mnod,mpri),ut(mnod,mpri),
     & rhand(mnod,mpri),rsource(mnod,mpri),ph(mnod),gP(mnod,mgas),
     & aream(mnod,mmin),sads(mnod),psi(mnod),
     & d(mnod,mads),supadn(mnod,msurf),phip(mnod,msurf),
     & surfads(msurf),ub(mbound,mpri),ctot(mnod,mpri),cnfact
cels6/18/08     & surfads(msurf),
cels6/18/08     & ub(mbound,mpri),ur(mrech,mpri),ctot(mnod,mpri),cnfact
!
      common/trangas3/pfug(mnod,mgas)  !New partial pressure =pfug(  )
      common/E4/phi(mnod)
      common/transpor1/nnod,iskip1(mnod)   !ns98/7/26 added iskip1
cels9/24/09      common/transport/izoneiw(mnod),izonebw(mnod),izonerw(mnod),
      common/transport/izoneiw(mnod),izonebw(mnod),
     +    izonem(mnod),izoneg(mnod),izoned(mnod),izonex(mnod),
     +    izonpp(mnod)
cels5/4/06  get rid of unused variables and add kelvin temp  common/heattran/tc(mnod),t0(mnod),tbound(mbound),trp(mrech)
      common/heattran/tc(mnod),tkelv(mnod)
      COMMON/ICONVERGENCE/ MAXITPFL,MAXITPTR,MAXITPCH,MAXITPAD
      COMMON/CONVERGENCE/TOLFL,TOLTR,TOLCH,TOLAD
cels9/24/09      COMMON/RESTART/  IRETURN,IRESTART, NBWTYPE,NRWTYPE        
      COMMON/RESTART/IRETURN,IRESTART,NBWTYPE        
      COMMON/ITERATION/ITERFL,ITERMOD,ITERTR,ITERCH,
     1         MAXITCH,ITERAD,MAXITAD,Iremove
cTX   COMMON/RITERATION/AVERITCH,AVERITAD,COUNTAD
      COMMON/RITERATION/AVERITCH,AVERITAD,COUNTAD,COUNTCH
cels6/26/08 added log10 of stimax
cels10/28/08 d1min unused      common/constraints/sl1min,d1min,stimax,dlstmx
      common/constraints/sl1min,stimax,dlstmx
      common/prtflags/minflag
      common/prtflags1/iconflag
      common/chfiles/OUTiter,OUTplot,OUTtime,OUTsolid,outgas,
     +  chemtitle
      common/ipplaw/ikplaw(mnod)
      common/dpplaw/aparpp(mnod),bparpp(mnod)
      common/sedim/vsed(mnod),vsedDF   !ns 4/09 
cels10/7/09 Add heats of reaction
cels10/7/09       common/heatsrctn/dhrctn(mnod),aks(maqx),akm(mmin),akg(mgas)
c
      CHARACTER*82 chemtitle
      CHARACTER*20 OUTiter,OUTplot,OUTtime,OUTsolid,outgas
      character*20 napri,naaqx,naaqt,namin,nagas,naexc,naads
      character*20 naads_min


