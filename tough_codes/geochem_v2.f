c
c
      subroutine hkfpar(T,tk)
c
c  This routine calculates parameters as a function of temperature
c  to calculate activity coefficents with equations in
c  Helgeson, Kirkham and Flowers, 1981, A.J.S. p.1249-1516
c  (see subroutine dh_hkf81 for details).
c
c   bhat NaCl (b NaCl = bhat/(2.303RT)) from Table 29 (bi here)
c   b Na+Cl- from Table 30 (bil here)
c   A and B Debye-Huckel (adh and bdh) from Table 1 (cols. 3 and 4)
c
c  Polynomial regression coefficients a, b, c, d, e, and f are
c  stored in data statements below to calculate parameters as:
c     parameter(T) = a + b*T + c*T + d*T**2 + e*T**3 + f*T**4
c  where T is temperature in C.
c
c  These are good for the range 0 to 300 C, BUT(!) the data
c  available to fit bi and bil did not go down to T = 0 C
c  (first point at 25 C) so the extrapolation may not be too good
c  below 25 C for bi and bil (A and B data cover the entire range),
c  but I checked that the extrapolated bi and bil values below
c  25 C vary smoothly down to 0 C.
c
c  Note the following units:
c    aft    yields adh in kg**0.5 mol**(-0.5)
c    bft    yields bdh in kg**0.5 mol**(-0.5) Anstrom**(-1)
c    bift   yields bihat(NaCl) in kg/mol * 1e+3
c    bilft  yields bil (Na+Cl-) in kg/mol * 1e+2
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      real*8 aft(5), bft(5), bift(5), bilft(5)
      real*8 t,tk
      common /dhparam/adh,bdh,bi,bil
c
      data aft/0.49276542d+00,0.31857945d-03,0.11628933d-04,
     &   -0.52038832d-07,0.12633045d-09/
      data bft/0.32476341d+00,0.12018502d-03,0.79530646d-06,
     &   -0.30410531d-08,0.56931304d-11/
      data bift/0.26538636d+01,-0.52889569d-02,-0.11009615d-03,
     &    0.46820513d-06,-0.11104895d-08/
      data bilft/-0.14769091d+02,0.22563951d+00,-0.10225385d-02,
     & 0.30349650d-05,-0.35468531d-08/
c
      func(a,b,c,d,e,T) = a+(b+(c+(d+e*T)*T)*T)*T
c
      adh=func(aft(1),aft(2),aft(3),aft(4),aft(5),T)
      bdh=func(bft(1),bft(2),bft(3),bft(4),bft(5),T)
      bihat=func(bift(1),bift(2),bift(3),bift(4),bift(5),T)*1.d-3
      bi=bihat/(4.57061d0*tk)  ![kg cal/mol]
      bil=func(bilft(1),bilft(2),bilft(3),bilft(4),bilft(5),T)*1.d-2
c
      return
      end
c
c
c-------------------------------------------------------------------------------
c
c
        subroutine dh_hkf81(no_ch)
C
C****************** Calculate activity coefficient of aqueous species ************
C   idryout is for dryout without update activity coefficients 0 and 1 normal, 2, skip
C
c  This routine computes activity coefficients of aqueous species
c  and the activity of water using pitzer model or an extended DH model according

c  to a user-specified ionic strength threshold and a user option.  The extended

c  DH model uses equations and parameters given in Helgeson, Kirkham and

c  Flowers, 1981, A.J.S. p.1249-1516. Also computes neutral species activities
c  from Setchenov equation (Langmuir 1997, Aqueous Environmental
c  Geochemistry, Prentice Hall, p. 144).
c  Charged species: individual ion activity coefficients calculated from:
c     equation 298 (which includes 121, 122, 129, 130, 169, 170 and 297)
c      - uses true ionic strenght.
c     use bhat NaCl (b NaCl = bhat/(2.303RT)) from Table 29 (bi here)
c     use b Na+Cl- from Table 30 (bil here)
c     use Rej from Table 3 (input in thermodynamic database as
c       a0 variable, but note that values are NOT a0 values)
c     use A and B Debye-Huckel (adh and bdh) from Table 1 (cols. 3 and 4)
c     caclulate a0 from input Rej values and eq. 125,
c       assuming other dominant anion (for cations) is Cl- (rej=1.81 A)
c       and cation (for anions) is Na+ (rej=1.91 A)
c
c  Activity of water calculated from:
c     osmotic coefficient using equation 190 and same parameters
c       as above, assuming similar simplifications as done for
c       the calculation of activity coefficients, but using
c       stoichiometric ionic strength.
c     equation 106 relating activity of water to the osm. coef.
c
c  Regression coefficients a, b, c, d, e to obtain A, B, bi abd bil
c  parameter as a function of temperature were calculated using
c  4th order polynomials as follow:
c        f(T) = a + b*T + c*T**2 + d*T**3 + e*T**4
c
c  Neutral species: assume gamma = 1 or,
c  for weak acids and dissolved gases:
c    log(activity coef)=sltout*ionic strength
c    where sltout is input in a0 variable as 100+sltout
c    For now, no temperature dependence on sltout is assumed as
c    the effect is small compared to variation of solubility (K) with temp.
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      integer*8 no_ch

      double precision lambda, lambd2, mstar, mchr
      double precision cpion(maqt)
c
      common /dhparam/adh,bdh,bi,bil  ! from routine hkfpar
      common/mop_react/mopr(20)       ! controling parameters for reactive transport
      common/str_thres/str_threshold  ! Ionic strength threshold for switch between pitzer and DH
      common/ion_str1/str             ! ionic strength
      common/savegam/ dgamt0(maqt,mpri),gamt0(maqt)
      common/aqkin16/NoTrans(mpri)    ! >0: not subject to transport
      parameter(CC=-1.0312d0,FF=0.0012806d0,GG=255.9d0,EE=0.4445d0,
     +   HH=-0.001606d0)
c
      SAVE ICALL
cpitz      SAVE gam
      DATA ICALL/0/
      data amwh2o/18.0152d0/

      ICALL=ICALL+1
      IF(ICALL.EQ.1) then
       WRITE(34,899)
  899 FORMAT(' dh_hkf81 1.0, 2003.7.30: Calculate activity coefficient'
     x' of aqueous species in geochem_v2.f')
c
c......Initializes gammas
       do i=1,npri
        gamp(i) = 1.d0
       end do
c
       do i=1,naqx
        gams(i) = 1.d0
       end do
c
      end if
C
c
c     hkf parameters at given temperature
c     call moved to routine assign to avoid recalc at all iterations
c       adh in kg**0.5 mol**(-0.5)  (this is A Debye-Huckel)
c       bdh in kg**0.5 mol**(-0.5) Angstrom**(-1) (this is B Debye-Huckel)
c       bi  in kg/mol * 1e+3 (this is bhat (NaCl) )
c       bil in kg/mol * 1e+2 (this is bil (Na+ Cl-) )
c
      no_ch=0
c
c ---calculate ionic strength and other global concentrations
c    all concentrations must be in molal scale (mol/kgH2O)
c
      sum=0.d0               !to compute true ionic strength
      sum2=0.d0              !to compute stoichiometric ionic str
      mstar=0.d0             !total solute in solution
      mchr=0.d0              !total solute excluding neutral species
c
c     primary species contribution
c
      do i=1,npri
        ct(i)=cp(i)
        if(ct(i).lt.1.0d-16) ct(i)=0.0d0
        if(i.ne.nw.and.i.ne.ne.and.NoTrans(i).ne.2) then
          zi2 = zsqi(i)
          sum=sum+zi2*cp(i)
          sum2=sum2+zi2*u2(i)/xh2o
          mstar=mstar+cp(i)
          if(z(i).ne.0.d0) mchr=mchr+cp(i)
          cpion(i)=cp(i)
         endif
      end do
c
c     secondary species contribution
c
      do i=1,naqx
         nprsec = npri+i
         ct(nprsec)=cs(i)
         if(ct(nprsec).lt.1.0d-16) ct(nprsec)=0.0d0
         sum=sum+zsqi(nprsec)*cs(i)
         mstar=mstar+cs(i)
         if(z(nprsec).ne.0.d0) then
          mchr=mchr+cs(i)
c         total concentrations excluding neutral species
          ncp=ncps(i)
           do n=1,ncp
             j=icps(i,n)
             cpion(j)=cpion(j)+stqs(i,n)*cs(i)
           enddo
         endif
      enddo
      if(sum.lt.0.d0) sum=0.d0
      str=0.5d0*sum
c
      if(str.gt.stimax) then
        no_ch=1
        return
      endif
c
      str2=str
      if(sum2.gt.0.d0)then
        str2=0.5d0*sum2
        stroot=dsqrt(str)
        stroo2=dsqrt(str2)
      else
        stroot=dsqrt(str)
        stroo2=stroot
      endif
c
c     Conversion factor for molality scale (from eq. 122, 169)
c     Note, this factor turns out equal to log(H2O mole fraction)
c     D-H yields gamma for mole fraction scale convention, we
c     use molality scale convention, therefore needs conversion
c
      capgam=-dlog10(1.d0+0.01801528d0*mstar)
c
c ---activity and osmotic coefficients
c
       summt = 0.d0
       bdhstrt = bdh*stroot
       bdhstr2 = bdh*stroo2
       bdh3str2 = bdh*bdh*bdh*str2
       cpgmdms = capgam/(0.0180153d0*mstar)
       bistr = bi*str
       bistr2 = bi*str2
       adhstrt = adh*stroot
c
c..... Primary species
c
cns3/4/2010       do i=1,npri
       do i=1,npaq
         gamp(i)=1.0d0
         zi2 = zsqi(i)
         lambda=1.d0+bdhstrt*azero(i)
         lambd2=1.d0+bdhstr2*azero(i)
c
        if(z(i).ne.0.d0) then
          gamlog=-adhstrt*zi2/lambda  +  capgam  +
     &    (omegahkf(i)*bistr + (bil - zabterm(i))*str)
c
            gamp(i)=10.d0**gamlog        ! activity coef.
            summt=summt + cpion(i) * (
     &       adh*zi2 / (azero3(i)*bdh3str2)
     &       * (lambd2- 1.d0/lambd2 - 2.d0*dlog(lambd2))
     &       + cpgmdms - 0.5d0*(omegahkf(i)*bistr2 +
     &             (bil-zabterm(i))*mchr*0.5d0)
     &       )
        else
          if(a0(i).gt.100.d0) then
            sltout=a0(i)-100.d0     !read in as sltout + 100 (flag)
            gamlog=sltout*str
            if(i.ne.nw) then
              gamp(i)=10.d0**gamlog        ! activity coef.
            endif
          endif
        endif
      enddo
c
c..... Secondary species
c
       do i=npri+1,npri+naqx
         gams(i-npri)=1.d0
        zi2 = zsqi(i)
        lambda=1.d0+bdhstrt*azero(i)
        lambd2=1.d0+bdhstr2*azero(i)
c
        if(z(i).ne.0.d0) then
          gamlog=-adhstrt*zi2/lambda  +  capgam  +
     &      (omegahkf(i)*bistr + (bil - zabterm(i))*str)
c
c           secondary charged species
            gams(i-npri)=10.d0**gamlog    ! activity coef.
c
c.......neutral species
        else
          if(a0(i).gt.100.d0) then
            sltout=a0(i)-100.d0     ! read in as sltout + 100 (flag)
            gamlog=sltout*str
            if(i.ne.nw) then
              gams(i-npri)=10.d0**gamlog    ! activity coef.
            endif
          endif
        endif
      enddo
c
c     activity of water (Note: value of mstar here will
c     strongly affect the value of the osmotic coefficient but
c     does not affect the activity of water because mstar
c     cancels out in the final expression)
c
      osmo = -2.303d0*summt/mstar
c
c     **note: for water, gamp is activity, NOT the coefficient!
      gamp(nw)=dexp(-osmo*mstar/55.50837d0)
c
c     activity coeficient of the electron
      if(ne.gt.0) gamp(ne)=1.d0
c
c     activity coefficient of the primary adsorbed species
      if(nd.gt.0) gamp(nd)=1.d0
c
c---------------------------
c     Revised salting out, so calculations are in readtherm_hkf
c---------------------------
c
      if(nis.gt.0)then
          gamln=(CC+FF*tk2+GG/tk2)*str - (EE+HH*tk2)*(str/(str+1.0d0))
          gsalt=dexp(gamln)
        do niis = 1, nis
          iis = indx_so(niis)
          if(iis.le.npri)gamp(iis)=gsalt
          if(iis.gt.npri)gams(iis-npri)=gsalt
        enddo
      endif
c
        do i=1,npri
          do j=1,npri
            dgamp(j,i)=0.0d0
            dgamt(j,i)=0.0d0
          enddo
        enddo
        do j=1,naqx
          npripj = npri+j
          do i=1,npri
           dgams(j,i)=0.0d0
           dgamt(npripj,i)=0.0d0
          enddo
        enddo
c
        goto 140
c
ct      endif
c
c------------------------------------------------------------------------------
c
       do i=1,npri+naqx
       if (gamt(i).gt.1.0d+06.or.gamt(i).lt.1.0d-06) goto 120
       enddo

       if(gamt(nw).gt.1.0d0.or.gamt(nw).lt.1.0d-02) then
         gamt(nw)=gamt0(nw)*xh2o
         gamt(nw)=dmax1(0.0d0,gamt(nw))
         gamt(nw)=dmin1(1.0d0,gamt(nw))
         gamp(nw)=gamt(nw)
       endif
c
       goto 140
c
 120   continue
c
       do i=1,npri
c
         gamp(i)=gamt0(i)
         gamt(i)=gamt0(i)
c
         do j=1,npri
           dgamp(j,i)=dgamt0(j,i)
           dgamt(j,i)=dgamt0(j,i)
         enddo
c
       enddo
c
       gamt(nw)=gamt(nw)*xh2o
       gamt(nw)=dmax1(0.0d0,gamt(nw))
       gamt(nw)=dmin1(1.0d0,gamt(nw))
       gamp(nw)=gamt(nw)
c
       do j=1,naqx
         jnpri=j+npri
         gams(j)=gamt0(jnpri)
         gamt(jnpri)=gamt0(jnpri)
c
         do i=1,npri
           dgams(j,i)=dgamt0(jnpri,i)
           dgamt(jnpri,i)=dgamt0(jnpri,i)
         enddo
c
       enddo
c
 140  continue

      return
      end
c
c-----------------------------------------------------------------------
c
        subroutine cs_cp
c
c  Calculates total solute in solution and partial derivatives
c  with respect to component species (mass balance in moles).
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(' cs_cp 1.0: 2003.7.30: Calculate total concentrations'
     X  ' and partial derivatives in geochem_v2.f')
c
c---------------------------
c      calculated now in readsolu        dlstmx = dlog10(stimax)
c---------------------------
c
c   calculation of the concentrations of secondary (derived) species
c   cp is concentration of primary species, cs of secondary species
c   (in moles per kg water!)
c
        do j=1,naqx
           cs(j)=-dlog10(gams(j))-aks(j)
           ncp=ncps(j)
           do n=1,ncp
              i=icps(j,n)
              cs(j)=cs(j)+stqs(j,n)*dlog10(cp(i)*gamp(i))
           end do
c
c......... add if/else cut-off to avoid convergence problems in rare cases
c          note: cs(j) above is first calculated in log form then converted below
c
           if(cs(j).gt.dlstmx) then
            cs(j)=cs(j)/10.d0
c..........to avoid underflow
           else if(cs(j).lt.-200.d0) then
                   cs(j)=1.d-200
           else
             cs(j)=10.d0**cs(j)
           end if
        end do
c
c   total solute in solution is u2 = total moles (primary + secondary species)
c   --first we get u2 contributions from component species and their derivatives
c   temp variables to limit calcs inside loop
c
        xh2rmh2 = xh2o*rmh2o
        do i=1,npri
          u2(i)=cp(i)*xh2o      !xh2o is unknown kg of water liquid
c.........du2 is derivative w/respect to cp (relative increment scheme)
          do k=1,npri
            du2(i,k)=0.0d0      !initialization
          end do
          du2(i,i)=u2(i)        !diagonal derivatives, using relative increment
        end do
        u2(nw)=xh2rmh2          !for water
        du2(nw,nw)=xh2rmh2      !also use relative increment
c
c......now we get u2 contributions from derived species, and their derivatives
       do j=1,naqx
          ncp=ncps(j)
            csjxh = cs(j)*xh2o
         do n=1,ncp
            i=icps(j,n)
            u2(i)=u2(i)+stqs(j,n)*csjxh
!
          do m=1,ncp
             k=icps(j,m)
            if(k.ne.nw) then
             du2(i,k)=du2(i,k)+stqs(j,m)*(csjxh*stqs(j,n))
            else
             du2(i,k)=du2(i,k)+csjxh*(stqs(j,n))
            endif
          enddo !m
         enddo !n
        enddo !j
c
        if (nh.ne.0) ph2=-dlog10(gamp(nh)*cp(nh))
c
        return
        end
c
c----------------------------------------------------------------------
c
        subroutine cmq_cp(iinit)
c
c  This routine computes saturation indexes of minerals at equilibrium
c  and fugacities of gases at equilibrium.
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
        common/satgas1/isatg(mgas),moutg(mgas),nsatg   !to keep track ofgas saturation
        common/satgas2/sg2
        COMMON/TRANGAS9/NGAS1        ! Number of gaseous species
        COMMON/min_SI/SIM(MNEL,MMIN) ! Mineral saturation index (log(Q/K)) for all nodes
        integer*8 iinit
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(' cmq_cp 1.0, 2003.7.30: Calculate saturation indexes of'
     X  ' minerals in geochem_v2.f')
c
c       Saturation index of minerals
c       si2 = log10 Q/K
c
        do m=1,nmequ
           paim=0.0d0
           ncp=ncpm(m)
           do n=1,ncp
             i=icpm(m,n)
             paim=paim+stqm(m,n)*(dlog10(cp(i)*gamp(i)))
           end do
c          save the old Q/K
           si2_old(m)=si2(m)
           si2(m)=paim-akm(m)
        end do
c
c       Saturation index of gases (cg= partial pressure)
c       sig2 = log10 Q/K/Pgas
c
        do j=1,ngas
          cg_gamg = cg(j)*gamg(j)
          if (cg_gamg .lt. 1.0d-100) cg_gamg = 1.0d-100
          paig=0.0d0
          ncp=ncpg(j)
          do n=1,ncp
            i=icpg(j,n)
            paig=paig+stqg(j,n)*(dlog10(cp(i)*gamp(i)))
          end do
c
         if(sg2.le.sl1min) then
c
c          no gas phase present, computes gas partial pressures directly from
c          chemical equilibrium (no need to compute sig2)
           if (ngas1.gt.0.or.iinit.eq.1) then
            cg(j)=paig-akg(j) !log f = log(q/k) always, even if no gas exsolved
            cg(j)=( 10.d0**cg(j) )/gamg(j) !cg is partial pressure (in bar), not fugacity
!
!...........No gas phase present, but still want a gas buffer
!
                                         else                                      
            sig2(j)=paig-dlog10(cg_gamg)-akg(j)
!
!............................................................
!
           endif
c
         else                  ! if sg2 .ge. sl1min
c
c          gas phase present, cg is computed iteratively in newtoneq, except for initial speciation
           if (iinit.eq.1) then
            cg(j)=paig-akg(j) !log f = log(q/k) always, even if no gas exsolved
            cg(j)=10.d0**cg(j)/gamg(j)   !cg is partial pressure (in bar), not fugacity
           endif
            sig2(j)=paig-dlog10(cg_gamg)-akg(j)
         endif

       end do
c
        return
        end
c
c
c
c-------------------------------------------------------------------------------
c
c
c
       SUBROUTINE CR_CP(ielem)
c
c  This routine calculates mineral reaction rates and their derivative
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        integer*8 ielem
c
c...... Dissolution kinetics
        common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
        common/iprkin/ideprec(mmin)
c...... Add block for rate ph dependence parameters
        common/phdep/aH1(mmin),aH2(mmin),aH1p(mmin),aH2p(mmin),
     +  aHexp(mmin),aHexpp(mmin),aOHexp(mmin),aOHexpp(mmin)
c...... Added common block for rate law designations
        common/irtlaw/nplaw(mmin)
        common/rksd0/ndep  ! number of minerals with species dependent dis/pre rate
        common/rksd5 /rkf_ds(mmin)
        common/rksd5p/rkprec_ds(mmin)
c
        double precision rkfdum(mmin),fdeltag(mmin),phterm(mmin),
     +     deriv(mmin),qkterm1(mmin),qkterm2(mmin),ssqk10(mmin)
        common/dispre/idispre(mmin)     !=1 only dissolution,=2 prec. =3 both
        COMMON/SOLUTE8/SL1(mnel)        ! new liquid saturation
        COMMON/DM/DELTEN,DELTEX,FOR,FORD
        common/minkin2/dr(mpri,mpri)
        COMMON/min_SI/SIM(MNEL,MMIN) ! Mineral saturation index (log(Q/K)) for all nodes
c.....  solid solutions
        common/solsol/iss(mmin),ncpss(msol),icpss(msol,mcpss),nss
c
        double precision skold(mmin), sumpre(mmin) !moved sumqk in common.inc
        integer*8 nflip(mmin)
        integer*8 isup,irtsp
c
        save skold, nflip
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(' cr_cp 1.0: 2003.7.30: Calculate mineral reaction rates'
     X  ' and their derivative in geochem_v2.f')
c
c-------------------------------------------------------------------
c...... for flip-flops convergence problems
        if(iterch.le.1) then
          do i=1,mmin
            nflip(i)=0
            skold(i)=1.d0
          end do
        end if
c
c...... Temperature in Kelvin
        eadum = ((1.0d0/tk2) - (1.0d0/298.15d0))/(gc1*1.d-3)
c
c---------------------------
c       Moved call to rate_sp after calculation of saturation indexes
c---------------------------
c
        sumsalts=0.d0   !sum of salts weights in kg per kg water (assume zero for now)
        vliq=1.d0
        factw=xh2o/vliq  !conversion factor = kg h2o liq/liter liquid
c
c       Calculates the saturation index for kinetic minerals
        do i=1,nmkin
           m=nmequ+i        !stoichiometries ordered 1 to nmin
           paimk(i)=0.d0
           ncp=ncpm(m)
           do n=1,ncp
             j=icpm(m,n)
             paimk(i)=paimk(i)+stqm(m,n)*dlog10(cp(j)*gamp(j))
           end do
           si2k(i)=paimk(i)-akin(i)
c          save mineral saturation index for all nodes
           si2k(i)=10.d0**si2k(i)
c
cc.........Limits should be treated with an exponent < 1 on the q/k term in the rate law
cc         to avoid overshoot of precipitation rate, give cutoff of Q/K
cc
cc         if (si2k(i) .gt. 1.0d4 .and. ikin(m) .eq. 1)   then          
cc              si2k(i) = 1.0d4
cc           end if
c
c---------------
c..........Added block below for sol.sol minor bug - need to calculate sum q/k here first
c---------------
c
c          initialize sum as si2k
           sumqk(i)=si2k(i)
c          sum of endmember amounts for dissolution mole fraction
           if(ielem.ne.0) sumpre(i)=pre(ielem,m)
        enddo
c
c...... To save Q/K only, go to end after calculating Q/K
        if(ielem.eq.0) go to 210
c
c.......Must be after the if statement for ielem=0
        ratet = deltex*phi(ielem)*sl1(ielem)*factw
c
c---------- For solid solutions: sum of saturation indices and ss saturation index
c
        do n=1,nss
           nem=ncpss(n)    !no. of endmembers in solid solution n
           smqk = 0.d0    !sum of q/k
           smpre= 0.d0    !sum of pre (mineral amounts)
           qssl = 0.d0     !log Q of existing sol sol
           akssl = 0.d0    !log K of existing sol sol 
           xi = 1.d0
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
              smqk=smqk+si2k(i)
c             sum of endmember amounts for dissolution mole fraction
              smpre=smpre+pre(ielem,m) 
           end do

c         If a solid solution is present, we calculate the saturation index
c         for that solid solution given its composition
           if (smpre.ne.0.d0) then
             do k=1,nem
               m=icpss(n,k)  !mineral index of endmember k in solid sol n
               i=m-nmequ
c              mole fraction for existing endmember (explicit, composition a start of time step)
               xi=pre(ielem,m)/smpre
               if(xi.le.0.d0) xi=1.d-30
c              Q and K for saturation index of solid solution assuming ai=xi
               qssl = qssl + xi*paimk(i) 
               akssl = akssl + xi*akin(i) + xi*dlog10(xi)    
             enddo 
             qkss = 10.d0**(qssl-akssl)

           endif
c
c          now needs second loop to assign same sumqk, qkss, etc. to all sol.sol endmembers
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
c             Save total amounts of endmember (will use later for dissolution mole fraction)
              sumpre(i)=smpre
              if(smpre.eq.0.d0) then
c               If the mineral is not present, we use the sum of q/k as the 
c               saturation index of of the solid solution (the sum corresponds 
c               to the maximum sat index of the ss, at equilibrium, for ideal ss)
                sumqk(i)=smqk    
              else 
c               If the mineral is present in rock, we use the actual saturation
c               index of the solid solution (calculated from its known composition)   
                sumqk(i)=qkss  
c               Only for dissolving solid solution: we reset the saturation index of
c               individual endmembers to values for the solid solution so that  
c               we get stoichiometric dissolution
                if(sumqk(i).lt.1.d0) si2k(i)=qkss   
              endif
           end do
        end do
c
c----- Calculate dissolution and precipitation rates -------
c
        do 200 i=1,nmkin
c
           m=nmequ+i        !stoichiometries ordered 1 to nmin
c
c..........Save mineral saturation index for all nodes (note: too many unnecessary calcs here)
           sim(ielem,m)=dlog10(si2k(i))
c
c..........Add skold and nflip stuff to avoid flip-flops with
c          small amounts of minerals if mineral not allowed to ppt.
c          and mineral with saturation window
         if (idispre(m).eq.1.or.ssqk(i).ne.0.d0)  then
           if (si2k(i).ge.1.0d0.and.skold(i).lt.1.0d0)
     &        nflip(i)=nflip(i)+1
           if (nflip(i).gt.5) then
              si2k(i)=skold(i)+(1.d0-skold(i))/1.1d0
           end if
              skold(i)=si2k(i)
         end if
c
c rate pH dependence - defaults for no dependence
         phterm(i)=1.d0    !dependence R = R*phterm
         deriv(i)=0.d0     !deriv of phterm w/resp H+ mol (cp(nh))
c arrays to limit use of power function
         qkterm1(i)=0.d0   ! ((q/k)**ck2 - 1)**ck1
         qkterm2(i)=0.d0   ! (q/k)**ck2
c
c------------------------------------------------------------------------------------------------
c----->Rate for Dissolution Case: set to positive value as R=rkfdum amin (1-q/k) ----------------
c........ Dissolution rate - mol/s/kgw (undersaturated mineral)
c        (final positive rate)  
c
c        did not taken into acct solid solutions if (si2k(i).lt.1.0d0)  then
c          note: by default si2k=sumqk for non solid solutions
c
         if (sumqk(i).lt.1.d0)  then
            if (idispre(m).ne.2)   then
c............ First calculate rate constant
c             rkfdum(i)= rkf(i)*dexp(-ea(i)*eadum)*(10.d0**
c      +         (acfdiss(i) + bcfdiss(i)*tk2 + ccfdiss(i)/tk2))
!
              if (rkf(i) .eq. -1.0d0)   then                              ! EPDC rate law
                 rkfdum(i)= acfdiss(i) + bcfdiss(i)*tk2 + ccfdiss(i)/tk2
!
              else if(ea(i).eq.0.d0.and.acfdiss(i).ne.0.d0)then
                 rkfdum(i)= rkf(i)*(10.d0**
     +             (acfdiss(i) + bcfdiss(i)*tk2 + ccfdiss(i)/tk2))
              elseif(ea(i).eq.0.d0.and.acfdiss(i).eq.0.d0.and.
     +             bcfdiss(i).eq.0.d0)then
                 rkfdum(i)= rkf(i)
              elseif(ea(i).gt.0.d0)then
                 rkfdum(i)= rkf(i)*dexp(-ea(i)*eadum)
              endif
              rkfdum_dis=rkfdum(i)
c
             if(idep(i).ge.2) then
c
c               Move call from above, for one mineral only, to save speed
                isup=0    !flag =1 for supersaturated
                irtsp = i
                call rate_species(eadum,irtsp,isup)
                rkfdum(i)=rkfdum(i)+rkf_ds(i) ! contributed from dependent species
              end if
            else if(idispre(m).eq.2)  then
              rkfdum(i)= 0.d0
            end if
!
!
!..... Functional dependence of rate on delta G-reaction
!..... use arrays to reduce use of power function later
!
!........Lasaga TST rate law
!
            if (ikin(m) .eq. 1)   then
!
c            Rewrote to reduce use of exponents
c            qkterm2(i) = si2k(i)**ck2(i)
             if(ck2(i).eq.1.d0)then
                qkterm2(i) = si2k(i)
             else
                qkterm2(i) = si2k(i)**ck2(i)
             endif
             fdeltag(i)= 1.0d0-qkterm2(i)
c
             if(ck1(i).eq.1.d0)then
                qkterm1(i) = fdeltag(i)
             else
                qkterm1(i) = fdeltag(i)**ck1(i)
             endif
!
!........Hellmann-Tisserand rate law
!........(Hellmann and Tisserand, Geochimica et Cosmochimica Acta, 70, 364-383, 2006)
!
            else if (ikin(m) .eq. 2)   then
!
             qkterm2(i) = dlog(si2k(i))
             qkterm2(i) = dabs(qkterm2(i))                         ! g, take absolute value
             fdeltag(i) = dexp(-1.0d0*ck1(i)*qkterm2(i)**ck2(i))   !     exp(-ng**m)
             qkterm1(i) = 1.0d0 - fdeltag(i)                       ! 1 - exp(-ng**m)
!
            end if
!
            rkin2(i)= rkfdum(i)*amin2(i)*qkterm1(i)
!
!
!...... rate dependence on aH+/aOH-
!
            if (idep(i).eq.1) then
              aH=gamp(nh)*cp(nh)
              if (aH.gt.aH1(i)) then
               phterm(i)=(aH/aH1(i))**aHexp(i)
               deriv(i)=aHexp(i)*phterm(i)   !relative increm
              else if(aH.lt.aH2(i)) then
               phterm(i)=(aH/aH2(i))**(-aOHexp(i))
               deriv(i)=-aOHexp(i)*phterm(i)   !relative increm
              end if
              rkin2(i)=rkin2(i)*phterm(i)
            end if
         end if
c
c-----------------------------------------------------------------------------------------------------
c---->Rate for Precipitation Case: set to negative value as R = -rkdum amin (q/k-1)----------------
c....... Precipitation rate (mol/s/kgw) (supersaturated mineral)
c    (negative final rate)
c
         ssqk10(i)=10.d0**ssqk(i)
c
c        did not acct for solid solutions         if (si2k(i).gt.ssqk10(i)) then
c        note: by default si2k=sumqk for non solid solutions
c
         if (sumqk(i).gt.ssqk10(i)) then
c
c.......... Will precipitate under kin. up to the supersaturation
c           window specified for equilibrium, then equilibrium takes over
            if (kineq(i).gt.0) then
              fdeltag(i) = 1.d0
              rkfdum(i)= 0.d0
              rkin2(i) = 0.d0
c
            else if (idispre(m).ne.1) then
c
c............ First calculate rate constant
              if (ideprec(i).eq.5) then
                rkfdum(i)=rkfdum_dis/10.0d0**akin(i) !for reversible (kpre=kdis/Keq)
                go to 199
              end if
c              rkfdum(i)= rkprec(i)*dexp(-eaprec(i)*eadum)*(10.d0**
c     +         (acfprec(i) + bcfprec(i)*tk2 + ccfprec(i)/tk2))
!
              if (rkf(i) .eq. -1.0d0)   then                           ! For EPDC rate law
               rkfdum(i)= acfprec(i) + bcfprec(i)*tk2 + ccfprec(i)/tk2
!
              else if(eaprec(i).eq.0.d0.and.acfprec(i).ne.0.d0)then
                 rkfdum(i)= rkprec(i)*(10.d0**
     +             (acfprec(i) + bcfprec(i)*tk2 + ccfprec(i)/tk2))
              elseif(eaprec(i).eq.0.d0.and.acfprec(i).eq.0.d0.and.
     +             bcfprec(i).eq.0.d0)then
                 rkfdum(i)= rkprec(i)
              elseif(eaprec(i).gt.0.d0)then
                 rkfdum(i)= rkprec(i)*dexp(-eaprec(i)*eadum)
              endif
c
              if(ideprec(i).ge.2) then
c
c               Move call from above, for one mineral only, to save speed
                isup=1    !flag =1 for supersaturated
                irtsp = i
                call rate_species(eadum,irtsp,isup)
                rkfdum(i)=rkfdum(i)+rkprec_ds(i) ! contributed from dependent species
              end if
199           continue
c
c Calculate difference between Q and K for bounding precipitation rate
c Allow for different precipitation rate laws
c Make sure precipitation gives a negative rate (negative sign)
c ns  fdeltag is negative.  We make it positive to avoid bombs at
c odd powers and for consistency with derivative further below
c
              if (nplaw(i).eq.0) then
c
c  Subtract supersaturation constant
c  use arrays to reduce use of power function later
!
!
!........Lasaga TST rate law
!
                 if (ikin(m) .eq. 1)   then
!
!                 Added if's to skip unit exponents, here and further down
                  if(ck2prec(i).ne.1.d0) then 
                    qkterm2(i) = si2k(i)**ck2prec(i)
                  else
                    qkterm2(i) = si2k(i)
                  endif
                  fdeltag(i) = qkterm2(i)-ssqk10(i)   !keep positive as (q/k-1) but use minus sign with rkin2 later
                  if(ck1prec(i).ne.1.d0) then 
                    qkterm1(i) = fdeltag(i)**ck1prec(i)
                  else
                    qkterm1(i) = fdeltag(i)
                  endif                    
!
!........Hellmann-Tisserand rate law
!
                 else if (ikin(m) .eq. 2)then
!
            qkterm2(i) = dlog(si2k(i))                                     ! g, Positive value for precipitation
                  if(ck2prec(i).ne.1.d0)then
            fdeltag(i) = dexp(-1.0d0*ck1prec(i)*qkterm2(i)**ck2prec(i))    !     exp(-ng**m)
                  else
            fdeltag(i) = dexp(-1.0d0*ck1prec(i)*qkterm2(i))                !     exp(-ng**m)
                  endif
            qkterm1(i) = 1.0d0 - fdeltag(i)                                ! 1 - exp(-ng**m)
!
                 end if
!
                 rkin2(i)= -rkfdum(i)*amin2(i)*qkterm1(i)
!
!
              else if(nplaw(i).eq.1) then
                  if(ck2prec(i).ne.1.d0) then
                fdeltag(i) = si2k(i)**ck2prec(i)
                  else
                fdeltag(i) = si2k(i)
                  endif
                qkterm1(i) = fdeltag(i)
c               Try something that drops off quickly as q/k -> 1
                rkin2(i)= -(rkfdum(i)*amin2(i)*(fdeltag(i)-
     +           (1.d0/fdeltag(i)**2)))
c
              end if
c
            else if (idispre(m).eq.1)then
               fdeltag(i) = 1.d0
               rkfdum(i)= 0.d0
               rkin2(i) = 0.d0
            end if
c
c  rate dependence on aH+/aOH-
c
            if (ideprec(i).eq.1) then
              aH=gamp(nh)*cp(nh)
              if (aH.gt.aH1p(i)) then
                phterm(i)=(aH/aH1p(i))**aHexpp(i)
                deriv(i)=aHexpp(i)*phterm(i)   !relative increm - mult by cp(nh)
              else if(aH.lt.aH2p(i)) then
                phterm(i)=(aH/aH2p(i))**(-aOHexpp(i))
                deriv(i)=-aOHexpp(i)*phterm(i)   !relative increm
              end if
                rkin2(i)=rkin2(i)*phterm(i)
            end if
c
        end if

c       Copy from cr_cp_num   Equilibrium case
        if (si2k(i).ge.(1.d0).and.si2k(i).le.ssqk10(i)) then
          rkfdum(i) = 0.d0
          rkin2(i) = 0.d0
        end if
c
200     continue
c
210     if(ielem.eq.0) return             !only calculates log Ks
c
c
c.....Solid solutions: we add a correction to the rates (computed above) of each endmember
c     Ideal solid solution is assumed for now (ai=xi)
c     The current scheme works only for kinetic minerals, and without exponents on q/k terms
c     Additional derivatives (for Jacobian) are ignored for now
c
        do n=1,nss
           nem=ncpss(n)    !no. of endmembers in solid solution n
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
c             Change for precipitation only
              if(sumqk(i).ge.1.d0)then  !precipitation  rkin2 = -kA(q/k - 1) 
                xmi=si2k(i)/sumqk(i)    !mole fraction of endmember (reflects fluid composition)
c               Correct the reaction rate to reflect solid solution, 
c               (need dabs for cases when q/k < 1, then qkterm1 < 0 while rkin2 > 0)
                rkin2(i)=rkin2(i)+dabs(rkin2(i)/qkterm1(i))*(xmi-1.d0)
c
              else if(sumpre(i).gt.0.d0) then   !dissolution  rkin2 = kA(1-q/k) 
c            For dissolution, the rate was computed above according to the saturation index of
c            of the solid solution (stoichiometric dissolution). We multiply by mole fraction.
                 xmi=pre(ielem,m)/sumpre(i)
                 if(xmi.le.0.d0) xmi=1.d-30
                 rkin2(i)=xmi*rkin2(i)
              endif 
           end do
        end do
c
c.......Moved from above with new loop
        do i=1,nmkin
          if (sumqk(i).ge.1.d0.and.sumqk(i).le.ssqk10(i)) then
            rkfdum(i) = 0.d0
            rkin2(i) = 0.d0
          end if
c         Set maximum rate so not so great to overshoot for dissolution
          m=nmequ+i
          ratetmp = pre(ielem,m)/ratet
          if(rkin2(i).gt.ratetmp) then 
             rkin2(i) = ratetmp
          endif 
c
c.........Need this if mineral is absent since c/out earlier
c         Note: ideally, we should reset rkfdum to be consistent with rkin2
c         for better convergence when a mineral runs out, even if non-zero
          if(rkin2(i).eq.0.d0) rkfdum(i)=0.d0
        enddo
c
c-----------------------------------------------------------------------------------------
c     Computes cr(j) (moles of j tied in kin minerals) and dr(j,i) (derivative of cr wrt i)
c
c     Initialize cr and dr
      do i=1,npri
           cr(i)=0.d0        !cr is moles tied up in kinetic minerals
           do j=1,npri
              dr(i,j)=0.d0    !dr is derivative of cr i with resprect to j
           end do
      end do
c
      do i=1,nmkin         !main loop through kinetic minerals
           m=nmequ+i
           ncp=ncpm(m)
c
c          Added this to avoid flip-flops near equilibrium at high rates
           if(dabs(qkterm1(i)).lt.1.d-12) rkin2(i) = 0.d0
c
           rk2xh = rkin2(i)*xh2o
c
c..........Precipitation
c
           if (si2k(i).gt.ssqk10(i).and.idispre(m).ne.1) then
!
              if (nplaw(i).ne.1) then
!
!................Lasaga TST rate law
!
                 if (ikin(m) .eq. 1)   then
!
                   fdrv1=-rkfdum(i)*amin2(i)*
     +                 si2k(i)*xh2o
     +                 * ( ck1prec(i)*qkterm1(i)/fdeltag(i) )
     +                 * ( ck2prec(i)*qkterm2(i)/si2k(i) )
!
!................Hellmann-Tisserand rate law
!
                 else if (ikin(m) .eq. 2)   then
                    fdrv1 = - rkfdum(i)*amin2(i)*xh2o               ! - K*A
     &                       *fdeltag(i)*ck1prec(i)*ck2prec(i)      ! *exp(-ng**m)*n*m
     &                       *qkterm2(i)**(ck2prec(i) - 1.0d0)      ! *g**(m-1)
!
                 end if
!
              else
                 fdrv1=-rkfdum(i)*amin2(i)*
     +              si2k(i)*xh2o *
     +              ( ck2prec(i)*fdeltag(i)/si2k(i) +
     +              2.d0*ck2prec(i)/(fdeltag(i)**2*si2k(i)) )
              end if
!
!
!..........Dissolution:
!
           else if(si2k(i).lt.(1.d0).and.idispre(m).ne.2) then

!.............Lasaga TST rate law
!
              if (ikin(m) .eq. 1)   then
!
                 fdrv1= - rkfdum(i)*amin2(i)*
     +            si2k(i)*xh2o
     +            * ( ck1(i)*qkterm1(i)/fdeltag(i) )
     +            * ( ck2(i)*qkterm2(i)/si2k(i) )
!
!.............Hellmann-Tisserand rate law
!
              else if (ikin(m) .eq. 2)   then
!
                 fdrv1 = - rkfdum(i)*amin2(i)*xh2o
     &                   * fdeltag(i)*ck1(i)*ck2(i)
     &                   * qkterm2(i)**(ck2(i) - 1.0d0)
!
              end if
!
!
           end if
!
!
           do n=1,ncp      !loop through components of kinetic mineral
             j=icpm(m,n)
c
c            Changes below to account for rates in moles/kgh2o/sec,
c            not just moles/sec.  Will not change much unless xh2o varies
c            much from 1.
c
             cr(j)=cr(j)+stqm(m,n)*rk2xh
             iph=0
             do k=1,ncp       !sub loop through components of kinetic mineral (for derivative wrt l)
                l=icpm(m,k)
                deriv1=0.d0       !added to store parts of dr(j,l)
                if (l.ne.nw) then
c                  Some modifications here - note that earlier, fdeltag was made always positive
c                  Derivatives for precipitation or dissolution, skipping suppressed phases
                   if (si2k(i).gt.ssqk10(i).and.idispre(m).ne.1) then
                      deriv1 = stqm(m,k)*stqm(m,n)*fdrv1
c                     pH dependent rates - addition for deriv. with respect to H+
                      if (l.eq.nh) then
                         deriv1=deriv1*phterm(i) + stqm(m,n)*rk2xh*
     &                       deriv(i)/phterm(i)
                         iph=1
                      end if
c
                   else if(si2k(i).lt.(1.d0).and.idispre(m).ne.2) then
                      deriv1 = stqm(m,k)*stqm(m,n)*fdrv1
c                     pH dependent rates - addition for deriv. with respect to H+
                      if (l.eq.nh) then
                         deriv1=deriv1*phterm(i) + stqm(m,n)*rk2xh*
     &                     deriv(i)/phterm(i)
                        iph=1
                      end if
                   end if
c                  Now we add the derivative parts to the total derivative so far
                   dr(j,l)=dr(j,l)+deriv1
c
c               Added next 3 statements for derivative w/ respect to water
                else
                   dr(j,l)=dr(j,l)+stqm(m,n)*rk2xh    !relative increment scheme
                end if
             end do  !k=1,ncp   with  l=icpm(m,k)
c
c            pH dependent rates - in case H+ is not in mineral stoichiometry
c            pH dependence introduces pH in equation
             if (iph.eq.0) then
              dr(j,nh)=dr(j,nh) + stqm(m,n)*rk2xh*
     &          deriv(i)/phterm(i)
             end if
           end do    !n=1,ncp   with  j=icpm(m,n)
c
      end do  !i=1,nmkin with m=nmequ+i
c
      return
      end
c
c
c
c-----------------------------------------------------------------------------
c
c
c
       SUBROUTINE CR_CP_Num(ielem)
c
c**************This routine calculates mineral reaction rates ********************
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        integer*8 ielem
c
c... Dissolution kinetics
        common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
        common/iprkin/ideprec(mmin)
c....Add block for rate ph dependence parameters
        common/phdep/aH1(mmin),aH2(mmin),aH1p(mmin),aH2p(mmin),
     +  aHexp(mmin),aHexpp(mmin),aOHexp(mmin),aOHexpp(mmin)
c... Added common block for rate law designations
        common/irtlaw/nplaw(mmin)
        common/rksd0/ndep  ! number of minerals with species dependent dis/pre rate
        common/rksd5 /rkf_ds(mmin)
        common/rksd5p/rkprec_ds(mmin)
c
        double precision rkfdum(mmin),fdeltag(mmin),phterm(mmin),
     +     qkterm1(mmin),qkterm2(mmin),ssqk10(mmin)
        common/dispre/idispre(mmin)     !=1 only dissolution,=2 prec. =3 both
        COMMON/SOLUTE8/SL1(mnel)        ! new liquid saturation
c
c                      porosity*saturation
        common/phisat/phisl1(mnel),phisg1(mnel)
c
        COMMON/DM/DELTEN,DELTEX,FOR,FORD
        common/minkin2/dr(mpri,mpri)
        COMMON/min_SI/SIM(MNEL,MMIN) ! Mineral saturation index (log(Q/K)) for all nodes
c... solid solutions
        common/solsol/iss(mmin),ncpss(msol),icpss(msol,mcpss),nss
        double precision skold(mmin), sumpre(mmin) !moved sumqk in common.inc
        integer*8 nflip(mmin)
        integer*8 isup,irtsp
c
        save skold, nflip
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(' cr_cp 1.0, 2003.7.30: Calculate mineral reaction rates'
     X  ' and their derivative in geochem_v2.f')
c
c
c-------------------------------------------------------------------
c... for flip-flops convergence problems
        if(iterch.le.1) then
          do i=1,mmin
            nflip(i)=0
            skold(i)=1.d0
          end do
        end if
c
c... Temperature in Kelvin
        eadum = ((1.0d0/tk2) - (1.0d0/298.15d0))/(gc1*1.d-3)
c
c-------------------------------------------------------------------
c
        sumsalts=0.d0   !sum of salts weights in kg per kg water (assume zero for now)
        vliq=1.d0
        factw=xh2o/vliq  !conversion factor = kg h2o liq/liter liquid
c
c       Calculates the saturation index for kinetic minerals
        do i=1,nmkin
           m=nmequ+i        !stoichiometries ordered 1 to nmin
           paimk(i)=0.d0
           ncp=ncpm(m)
           do n=1,ncp
             j=icpm(m,n)
             paimk(i)=paimk(i)+stqm(m,n)*dlog10(cp(j)*gamp(j))
           end do
           si2k(i)=paimk(i)-akin(i)
c-------save mineral saturation index for all nodes
           si2k(i)=10.d0**si2k(i)
c
c---------------
c..........Added block below for sol.sol minor bug - need to calculate sum q/k here first
c---------------
c          initialize sum as si2k
           sumqk(i)=si2k(i)
c          sum of endmember amounts for dissolution mole fraction
           if(ielem.ne.0) sumpre(i)=pre(ielem,m)
        enddo
c
c... To save Q/K only, go to end after calculating Q/K
        if(ielem.eq.0) go to 210
c
c.......Must be after the if statement for ielem=0
        ratet = deltex*phisl1(ielem)*factw
c
c---------- For solid solutions: sum of saturation indices and ss saturation index
c
        do n=1,nss
           nem=ncpss(n)    !no. of endmembers in solid solution n
           smqk = 0.d0    !sum of q/k
           smpre= 0.d0    !sum of pre (mineral amounts)
           qssl = 0.d0     !log Q of existing sol sol
           akssl = 0.d0    !log K of existing sol sol 
           xi = 1.d0
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
              smqk=smqk+si2k(i)
c             sum of endmember amounts for dissolution mole fraction
              smpre=smpre+pre(ielem,m) 
           end do
c
c         if a solid solution is present, we calculate the saturation index
c         for that solid solution given its composition
           if (smpre.ne.0.d0) then
             do k=1,nem
               m=icpss(n,k)  !mineral index of endmember k in solid sol n
               i=m-nmequ
c              mole fraction for existing endmember (explicit, composition a start of time step)
               xi=pre(ielem,m)/smpre
               if(xi.le.0.d0) xi=1.d-30
c              Q and K for saturation index of solid solution assuming ai=xi
               qssl = qssl + xi*paimk(i) 
               akssl = akssl + xi*akin(i) + xi*dlog10(xi)    
             enddo 
             qkss = 10.d0**(qssl-akssl)

           endif
c
c  now needs second loop to assign same sumqk, qkss, etc. to all sol.sol endmembers
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
c             Save total amounts of endmember (will use later for dissolution mole fraction)
              sumpre(i)=smpre
              if(smpre.eq.0.d0) then
c               If the mineral is not present, we use the sum of q/k as the 
c               saturation index of of the solid solution (the sum corresponds 
c               to the maximum sat index of the ss, at equilibrium, for ideal ss)
                sumqk(i)=smqk    
              else 
c               If the mineral is present in rock, we use the actual saturation
c               index of the solid solution (calculated from its known composition)   
                sumqk(i)=qkss  
c               Only for dissolving solid solution: we reset the saturation index of
c               individual endmembers to values for the solid solution so that  
c               we get stoichiometric dissolution
                if(sumqk(i).lt.1.d0) si2k(i)=qkss   
              endif
           end do
        end do
c
c----- Calculate dissolution and precipitation rates -------
c
        do 200 i=1,nmkin
c
           m=nmequ+i        !stoichiometries ordered 1 to nmin
c
c..........Save mineral saturation index for all nodes (note: too many unnecessary calcs here)
           sim(ielem,m)=dlog10(si2k(i))

c         if(dabs(sim(ielem,m)).le.0.001d0.and.iterch.gt.150) then
c           dummy=0.d0
c         endif
c
c ---------- add skold and nflip stuff to avoid flip-flops with
c  small amounts of minerals if mineral not allowed to ppt.
c and mineral with saturation window
         if (idispre(m).eq.1.or.ssqk(i).ne.0.d0)  then
           if (si2k(i).ge.1.0d0.and.skold(i).lt.1.0d0)
     &        nflip(i)=nflip(i)+1
           if (nflip(i).gt.5) then
              si2k(i)=skold(i)+(1.d0-skold(i))/1.1d0
           end if
              skold(i)=si2k(i)
         end if
c
c rate pH dependence - defaults for no dependence
         phterm(i)=1.d0    !dependence R = R*phterm
c        arrays to limit use of power function
         qkterm1(i)=0.d0   ! ((q/k)**ck2 - 1)**ck1
         qkterm2(i)=0.d0   ! (q/k)**ck2

c
c------------------------------------------------------------------------------------------------
c----->Rate for Dissolution Case: set to positive value as R=rkfdum amin (1-q/k) ----------------
c........ Dissolution rate - mol/s/kgw (undersaturated mineral)
c    (final positive rate)  
c
cels9/4/08 did not taken into acct solid solutions if (si2k(i).lt.1.0d0)  then
c          note: by default si2k=sumqk for non solid solutions
         if (sumqk(i).lt.1.d0)  then
            if (idispre(m).ne.2)   then
c... First calculate rate constant
cels9/3/08              rkfdum(i)= rkf(i)*dexp(-ea(i)*eadum)*(10.d0**
cels9/3/08      +         (acfdiss(i) + bcfdiss(i)*tk2 + ccfdiss(i)/tk2))
              if(ea(i).eq.0.d0.and.acfdiss(i).ne.0.d0)then
                 rkfdum(i)= rkf(i)*(10.d0**
     +             (acfdiss(i) + bcfdiss(i)*tk2 + ccfdiss(i)/tk2))
              elseif(ea(i).eq.0.d0.and.acfdiss(i).eq.0.d0.and.
     +             bcfdiss(i).eq.0.d0)then
                 rkfdum(i)= rkf(i)
              elseif(ea(i).gt.0.d0)then
                 rkfdum(i)= rkf(i)*dexp(-ea(i)*eadum)
              endif
              rkfdum_dis=rkfdum(i)
c
             if(idep(i).ge.2) then
c
c               Move call from above, for one mineral only, to save speed
                isup=0    !flag =1 for supersaturated
                irtsp = i
!
                if (idep(i).eq.2 .or. idep(i).eq.3)   then 
                   call rate_species(eadum,irtsp,isup)
                end if
!...............Labradorite rate law (Carroll and Knauss, 2005; Chemical Geology)
                if (idep(i).eq.4)   then 
                   RTdum = 2.303d0*(gc1*1.d-3)*tk2
                   call rate_species_Labradorite(RTdum,irtsp,isup)
                end if
!
                rkfdum(i)=rkfdum(i)+rkf_ds(i) ! contributed from dependent species
              end if
            else if(idispre(m).eq.2)  then
              rkfdum(i)= 0.d0
            end if
!
!
!... Functional dependence of rate on delta G-reaction
!...... use arrays to reduce use of power function later
!
!........Lasaga TST rate law
!
            if (ikin(m) .eq. 1)   then
!
             qkterm2(i) = si2k(i)**ck2(i)
             fdeltag(i)= 1.0d0-qkterm2(i)
             qkterm1(i) = fdeltag(i)**ck1(i)
!
!........Hellmann-Tisserand rate law
!........(Hellmann and Tisserand, Geochimica et Cosmochimica Acta, 70, 364-383, 2006)
!
            else if (ikin(m) .eq. 2)   then
!
             qkterm2(i) = dlog(si2k(i))
             qkterm2(i) = dabs(qkterm2(i))                           ! g, take absolute value
             fdeltag(i) = dexp(-1.0d0*ck1(i)*qkterm2(i)**ck2(i))     !     exp(-ng**m)
             qkterm1(i) = 1.0d0 - fdeltag(i)                         ! 1 - exp(-ng**m)
!
            end if
!
            rkin2(i)= rkfdum(i)*amin2(i)*qkterm1(i)
!
!
c...... rate dependence on aH+/aOH-
            if (idep(i).eq.1) then
              aH=gamp(nh)*cp(nh)
              if (aH.gt.aH1(i)) then
               phterm(i)=(aH/aH1(i))**aHexp(i)
              else if(aH.lt.aH2(i)) then
               phterm(i)=(aH/aH2(i))**(-aOHexp(i))
              end if
              rkin2(i)=rkin2(i)*phterm(i)
            end if
         end if
c
c-----------------------------------------------------------------------------------------------------
c---->Rate for Precipitation Case: set to negative value as R = -rkdum amin (q/k-1)----------------
c....... Precipitation rate (mol/s/kgw) (supersaturated mineral)
c    (negative final rate)
         ssqk10(i)=10.d0**ssqk(i)
c          note: by default si2k=sumqk for non solid solutions
         if (sumqk(i).gt.ssqk10(i)) then
c
c.. will precipitate under kin. up to the supersaturation
c     window specified for equilibrium, then equilibrium takes over
            if (kineq(i).gt.0) then
              fdeltag(i) = 1.d0
              rkfdum(i)= 0.d0
              rkin2(i) = 0.d0
c       if (idispre(nmequ+i).ne.1)then
c        elseif (idispre(nmequ+i).ne.1)then
            else if (idispre(m).ne.1) then
c
c... First calculate rate constant
              if (ideprec(i).eq.5) then
                rkfdum(i)=rkfdum_dis/10.0d0**akin(i) !for reversible (kpre=kdis/Keq)
                go to 199
              end if
c              rkfdum(i)= rkprec(i)*dexp(-eaprec(i)*eadum)*(10.d0**
c     +         (acfprec(i) + bcfprec(i)*tk2 + ccfprec(i)/tk2))
              if(eaprec(i).eq.0.d0.and.acfprec(i).ne.0.d0)then
                 rkfdum(i)= rkprec(i)*(10.d0**
     +             (acfprec(i) + bcfprec(i)*tk2 + ccfprec(i)/tk2))
              elseif(eaprec(i).eq.0.d0.and.acfprec(i).eq.0.d0.and.
     +             bcfprec(i).eq.0.d0)then
                 rkfdum(i)= rkprec(i)
              elseif(eaprec(i).gt.0.d0)then
                 rkfdum(i)= rkprec(i)*dexp(-eaprec(i)*eadum)
              endif
c
              if(ideprec(i).ge.2) then
c
c               Move call from above, for one mineral only, to save speed
                isup=1    !flag =1 for supersaturated
                irtsp = i
!
                if (ideprec(i).eq.2 .or. ideprec(i).eq.3) then
                   call rate_species(eadum,irtsp,isup)
                end if
!...............Labradorite rate law (Carroll and Knauss, 2005; Chemical Geology)
                if (ideprec(i).eq.4) then
                   RTdum = 2.303d0*(gc1*1.d-3)*tk2
                   call rate_species_Labradorite(RTdum,irtsp,isup)
                end if
!
                rkfdum(i)=rkfdum(i)+rkprec_ds(i) ! contributed from dependent species
              end if
199           continue
c
c Calculate difference between Q and K for bounding precipitation rate
c Allow for different precipitation rate laws
c Make sure precipitation gives a negative rate (negative sign)
c ns  fdeltag is negative.  We make it positive to avoid bombs at
c odd powers and for consistency with derivative further below
              if (nplaw(i).eq.0) then
c Subtract supersaturation constant
c  use arrays to reduce use of power function later
!
!
!........Lasaga TST rate law
!
                 if (ikin(m) .eq. 1)   then
!
                  qkterm2(i) = si2k(i)**ck2prec(i)
                  fdeltag(i) = qkterm2(i)-ssqk10(i)   !keep positive as (q/k-1) but use minus sign with rkin2 later
                  qkterm1(i) = fdeltag(i)**ck1prec(i)
!
!........Hellmann-Tisserand rate law
!
                 else if (ikin(m) .eq. 2)   then
!
            qkterm2(i) = dlog(si2k(i))       ! g, Positive value for precipitation
            fdeltag(i) = dexp(-1.0d0*ck1prec(i)*qkterm2(i)**ck2prec(i))    !     exp(-ng**m)
            qkterm1(i) = 1.0d0 - fdeltag(i)                                ! 1 - exp(-ng**m)
!
                 end if
!
                 rkin2(i)= -rkfdum(i)*amin2(i)*qkterm1(i)
!
!
              else if(nplaw(i).eq.1) then
                fdeltag(i) = si2k(i)**ck2prec(i)
c               Try something that drops off quickly as q/k -> 1
                rkin2(i)= -(rkfdum(i)*amin2(i)*(fdeltag(i)-
     +           (1.d0/fdeltag(i)**2)))
c
              end if

            else if (idispre(m).eq.1)then
               fdeltag(i) = 1.d0
               rkfdum(i)= 0.d0
               rkin2(i) = 0.d0
            end if
c
c           rate dependence on aH+/aOH-
            if (ideprec(i).eq.1) then
              aH=gamp(nh)*cp(nh)
              if (aH.gt.aH1p(i)) then
                phterm(i)=(aH/aH1p(i))**aHexpp(i)
              else if(aH.lt.aH2p(i)) then
                phterm(i)=(aH/aH2p(i))**(-aOHexpp(i))
              end if
                rkin2(i)=rkin2(i)*phterm(i)
            end if
c
        end if
c
c       Copy from cr_cp_num   Equilibrium case
        if (si2k(i).ge.(1.d0).and.si2k(i).le.ssqk10(i)) then
          rkfdum(i) = 0.d0
          rkin2(i) = 0.d0
        end if
c
200     continue
c
210     if(ielem.eq.0) return             !only calculates log Ks
c
c------ Solid solutions: we add a correction to the rates (computed above) of each endmember
c Ideal solid solution is assumed for now (ai=xi)
c !!! The current scheme works only for kinetic minerals, and without exponents on q/k terms
c !!! Additional derivatives (for Jacobian) are ignored for now
        do n=1,nss
           nem=ncpss(n)    !no. of endmembers in solid solution n
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
c             change for precipitation only
              if(sumqk(i).ge.1.d0)then    !precipitation  rkin2 = -kA(q/k - 1) 
                xmi=si2k(i)/sumqk(i)  !mole fraction of endmember (reflects fluid composition)
c               correct the reaction rate to reflect solid solution, 
c               (need dabs for cases when q/k < 1, then qkterm1 < 0 while rkin2 > 0)
                rkin2(i)=rkin2(i)+dabs(rkin2(i)/qkterm1(i))*(xmi-1.d0)

              else if(sumpre(i).gt.0.d0) then   !dissolution  rkin2 = kA(1-q/k) 
c            for dissolution, the rate was computed above according to the saturation index of
c            of the solid solution (stoichiometric dissolution). We multiply by mole fraction.
                 xmi=pre(ielem,m)/sumpre(i)
                 if(xmi.le.0.d0) xmi=1.d-30
                 rkin2(i)=xmi*rkin2(i)
              endif 
           end do
        end do
c
c.......Moved from above with new loop
        do i=1,nmkin
          if (sumqk(i).ge.1.d0.and.sumqk(i).le.ssqk10(i)) then
            rkfdum(i) = 0.d0
            rkin2(i) = 0.d0
          end if
c         Set maximum rate so not so great to overshoot for dissolution
          m=nmequ+i
          ratetmp = pre(ielem,m)/ratet
          if(rkin2(i).gt.ratetmp) then 
             rkin2(i) = ratetmp
          endif 
c
c........Need this if mineral is absent since c/out earlier
c  Note: ideally, we should reset rkfdum to be consistent with rkin2
c        for better convergence when a mineral runs out, even if non-zero
          if(rkin2(i).eq.0.d0) rkfdum(i)=0.d0
        enddo
c
c       initialize cr
        do i=1,npri
           cr(i)=0.d0        !cr is moles tied up in kinetic minerals
        end do
c       computes cr
        do i=1,nmkin
           m=nmequ+i
           ncp=ncpm(m)
c
c          add this to avoid flip-flops near equilibrium at high rates
           if(dabs(qkterm1(i)).lt.1.d-12) rkin2(i) = 0.d0
c
           rk2xh = rkin2(i)*xh2o
c
          do n=1,ncp
             j=icpm(m,n)
c Changes below to account for rates in moles/kgh2o/sec,
c not just moles/sec.  Will not change much unless xh2o varies
c much from 1.
             cr(j)=cr(j)+stqm(m,n)*rk2xh
          end do
       end do !i=1,nmkin with m=nmequ+i
c
      return
      end
c
c
c
c-------------------------------------------------------------------------------
c
c
c
        subroutine dcr_dcp_num(ielem)
c
C************** Calculates derive of reaction rates by numerical method ***********
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
        common/minkin2/dr(mpri,mpri)
        double precision crold(mmin)
        integer*8 ielem
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(' dcx_dcp 1.0: 2003.9.11: Calculates derivative of'
     x' reaction rates by numerical method in geochem_v2.f')
c
        do j=1,npri
           crold(j)=cr(j)
        end do
        do i=1,npri
           dd=cp(i)*1.0d-07
           cp(i)=cp(i)+dd
           call cr_cp_num(ielem)
           do j=1,npri
              dr(j,i)=(cr(j)-crold(j))/dd
              cr(j)=crold(j)
           end do
           cp(i)=cp(i)-dd
        end do
c
        do j=1,npri
           do i=1,npri
              dr(j,i)=dr(j,i)*cp(i)   ! convert to relative increment
           end do
        end do
c
        return
        end
c
c
c
c-------------------------------------------------------------------------------
c
c
c
        subroutine rate_species(eadum,irtsp,isup)
c
c******** Calculate part of rate constant contributed from dependent species *******
c
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/iprkin/ideprec(mmin)
c
      common/rksd1/  ids(mmin,mechm,mechsp) ! dependent species pointer in CT(naqt)
c                       (mmin, No. mechanisms, No.species involved)
      common/rksd2/ rkds(mmin,mechm)   ! constant for the dependent species
      common/rksd21/ eads(mmin,mechm)  ! activation energy for the dependent species
      common/rksd3/expds(mmin,mechm,mechsp) ! exponential term for the dependent species
      common/rksd4/ndis(mmin)           ! number of additional mechanisms
      common/rksd41/ nspds(mmin,mechsp) ! number of speciess involved in one mechanism
      common/rksd1p/  idsp(mmin,mechm,mechsp)
      common/rksd2p/ rkdsp(mmin,mechm)
      common/rksd21p/ eadsp(mmin,mechm)
      common/rksd3p/expdsp(mmin,mechm,mechsp)
      common/rksd4p/npre(mmin)
      common/rksd41p/ nsppr(mmin,mechsp)
      common/rksd5 /rkf_ds(mmin)
      common/rksd5p/rkprec_ds(mmin)
      integer*8 irtsp
C
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(11,899)
  899   FORMAT(6X,'rate_species 1.0  11 April     2004',6X,
     X  'Calculate part of rate constant contributed from',
     X                                ' dependent species')
c
c------For species dependent rate law: H+,  k(h+), expo(h+). term.
c------------------------------------ oh-, k(oh-), expo(oh-). term.
c-------------------------------------nadum2, rkds,  expds
c------------------------------------k=k(h+)*[H+]**expo +....
c
c        Passed i 
         i = irtsp
      if(isup.eq.0) then
c
c........For dissolution:
c
         rkf_ds(i)=0.0d0
         do j=1,ndis(i)        ! loop over mechanism
               prod_term=1.0d0
               nsp=nspds(i,j)  ! number of species involved in one mechanism
               do isp=1,nsp
                  js=ids(i,j,isp)
                   if(js.le.npri) then
                     cj=cp(js)
                     gamj=gamp(js)
                     aj=cj*gamj
                   else
                     cj=cs(js-npri)
                     gamj=gams(js-npri)
                     aj=cj*gamj
                   endif
!
!..................Power term or half concentration 
                   Power_Half = expds(i,j,isp)
!
                   if (idep(i) .eq. 2)  then
                      prod_term=prod_term*aj**Power_Half
!
!..................Inhibition rate law
!..................expds(i,j,isp) store half concentration
!
                   else if (idep(i) .eq. 3)  then    
                      prod_term=prod_term*(Power_Half/(Power_Half+cj))
                   end if                    
!
               end do
!
               cst_kin =  rkds(i,j)*dexp(-eads(i,j)*eadum)*prod_term
               rkf_ds(i)=rkf_ds(i) + cst_kin
!
         end do   !j=1,ndis(i)
!
!.....For precipitation:
!
      else if (isup.eq.1) then 

         rkprec_ds(i)=0.0d0
         do j=1,npre(i)
               prod_term=1.0d0
               nsp=nsppr(i,j)  ! number of species involved in one mechanism
               do isp=1,nsp
                  js=idsp(i,j,isp)
                  if(js.le.npri) then
                     cj=cp(js)
                     gamj=gamp(js)
                     aj=cj*gamj
                  else
                     cj=cs(js-npri)
                     gamj=gams(js-npri)
                     aj=cj*gamj
                  endif
!
!..................Power term or half concentration 
                   Power_Half = expdsp(i,j,isp)
!
                   if (ideprec(i) .eq. 2)  then
                      prod_term=prod_term*aj**Power_Half
!
!..................Inhibition rate law
!..................expds(i,j,isp) store half concentration
!
                   else if (ideprec(i) .eq. 3)  then    
                      prod_term=prod_term*(Power_Half/(Power_Half+cj))
                   end if                    
!
               end do
!
               cst_kin = rkdsp(i,j)*dexp(-eadsp(i,j)*eadum)*prod_term
               rkprec_ds(i)=rkprec_ds(i) + cst_kin
!
!
         end do  !j=1,npre(i)
!
      end if
!
      return
      end
c
c
c
c-------------------------------------------------------------------------------
c
c
c
        subroutine rate_species_Labradorite(RTdum,irtsp,isup)
c
!.......Labradorite rate law (Carroll and Knauss, 2005; Chemical Geology)
!.......Available only for numerical derivatives
!.......Under the option of IDEP=4 in chemical.inp
c******** Calculate part of rate constant contributed from dependent species *******
c
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
      common/iprkin/ideprec(mmin)
c
      common/rksd1/  ids(mmin,mechm,mechsp) ! dependent species pointer in CT(naqt)
c                       (mmin, No. mechanisms, No.species involved)
      common/rksd2/ rkds(mmin,mechm)   ! constant for the dependent species
      common/rksd21/ eads(mmin,mechm)  ! activation energy for the dependent species
      common/rksd3/expds(mmin,mechm,mechsp) ! exponential term for the dependent species
      common/rksd4/ndis(mmin)           ! number of additional mechanisms
      common/rksd41/ nspds(mmin,mechsp) ! number of speciess involved in one mechanism
      common/rksd1p/  idsp(mmin,mechm,mechsp)
      common/rksd2p/ rkdsp(mmin,mechm)
      common/rksd21p/ eadsp(mmin,mechm)
      common/rksd3p/expdsp(mmin,mechm,mechsp)
      common/rksd4p/npre(mmin)
      common/rksd41p/ nsppr(mmin,mechsp)
      common/rksd5 /rkf_ds(mmin)
      common/rksd5p/rkprec_ds(mmin)
      integer*8 irtsp
C
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(11,899)
  899   FORMAT(6X,'rate_species 1.0  11 April     2004',6X,
     X  'Calculate part of rate constant contributed from',
     X                                ' dependent species')
c
c------For species dependent rate law: H+,  k(h+), expo(h+). term.
c------------------------------------ oh-, k(oh-), expo(oh-). term.
c-------------------------------------nadum2, rkds,  expds
c------------------------------------k=k(h+)*[H+]**expo +....
c
      i = irtsp
      if(isup.eq.0) then
c
c........For dissolution:
c
         rkf_ds(i)=0.0d0
         do j=1,ndis(i)        ! loop over mechanism
               prod_term=1.0d0
               nsp=nspds(i,j)  ! number of species involved in one mechanism
               do isp=1,nsp
                  js=ids(i,j,isp)
                   if(js.le.npri) then
                     cj=cp(js)
                     gamj=gamp(js)
                     aj=cj*gamj
                   else
                     cj=cs(js-npri)
                     gamj=gams(js-npri)
                     aj=cj*gamj
                   endif
!
!..................Power term or half concentration 
                   Power_Half = expds(i,j,isp)
!
                   prod_term=prod_term*aj**Power_Half
!
               end do
!
!..............Labradorite rate law (Carroll and Knauss, 2005; Chemical Geology)
!             
               Kt = ccfdiss(i)    !  ccfdiss(i) stored Kt
               HAl_Term1 = prod_term*Kt
               HAl_Term  = HAl_Term1/(1.0d0+HAl_Term1)   
               cst_kin =  rkds(i,j)*10.0d0**(-eads(i,j)/RTdum)*HAl_Term
               rkf_ds(i)=rkf_ds(i) + cst_kin
!
         end do   !j=1,ndis(i)
!
!
!.....For precipitation:
!
      else if (isup.eq.1) then 

         rkprec_ds(i)=0.0d0
         do j=1,npre(i)
               prod_term=1.0d0
               nsp=nsppr(i,j)  ! number of species involved in one mechanism
               do isp=1,nsp
                  js=idsp(i,j,isp)
                  if(js.le.npri) then
                     cj=cp(js)
                     gamj=gamp(js)
                     aj=cj*gamj
                  else
                     cj=cs(js-npri)
                     gamj=gams(js-npri)
                     aj=cj*gamj
                  endif
!
!..................Power term or half concentration 
                   Power_Half = expdsp(i,j,isp)
!
                   prod_term=prod_term*aj**Power_Half
!
               end do
!
!..............Labradorite rate law (Carroll and Knauss, 2005; Chemical Geology)
!             
               Kt = ccfprec(i)    !  ccfprec(i) stored Kt
               HAl_Term1 = prod_term*Kt
               HAl_Term  = HAl_Term1/(1.0d0+HAl_Term1)   
               cst_kin = rkdsp(i,j)*10.0d0**(-eadsp(i,j)/RTdum)*HAl_term
               rkprec_ds(i)=rkprec_ds(i) + cst_kin
!
         end do  !j=1,npre(i)
!
      end if
!
      return
      end
c
c-------------------------------------------------------------------------------
c
        subroutine cx_ct
c
c*************** Calculate concentrations of exchanged cations ****************
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        double precision dum(mexc),bx(mexc)
!
        common/satgas2/sg2
!
!.....Extract rock density for geochemical calculations such as exchange and sorption
      common/rock_density2/denss2, sl2, phisl2, a_fmr2
!
!
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(1x,15('*'),' cx_ct 1.0, 2003.7.30: Calculate'
     X  ' concentrations of exchanged cations in geochem_v2.f',15('*'))
c
c       check for positive charge of exchangeble ions
c
        do 100 j=1,nexc
        if (z(nbx(j)).le.0.d0) then
           write (*,990)
990        format (/3x,'error: check charge of exchageable cations')
           stop
        endif
100     continue
c
c       the exchange coeficient for primary adsorbate should be 1.0
c
        if (nexc.gt.0.and.ekx(nx).ne.1.d0) then
           write (*,995)
995        format (/3x,'error:check exchange coef.for primary ads.')
           stop
        endif
c
c       the terms of addition of total eq. fraction as a funtion of bx(nx)
c
        a1=0.d0
        a2=0.d0
        a3=0.d0
c
c       gaines-thomas and vanselow conventions
c
        if (iex.eq.1.or.iex.eq.2) then
        do 200 j=1,nexc
        dum(j)= ekx(j)**(-z(nbx(j)))*ct(nbx(j))*gamt(nbx(j))*
     +          (ct(nbx(nx))*gamt(nbx(nx)))**(-z(nbx(j))/z(nbx(nx)))
        if ((z(nbx(j))/z(nbx(nx))).eq.1.d0) a1=a1+dum(j)
        if ((z(nbx(j))/z(nbx(nx))).eq.2.d0) a2=a2+dum(j)
        if ((z(nbx(j))/z(nbx(nx))).eq.3.d0) a3=a3+dum(j)
200     continue
c       resolution of: a1*bx(nx) + a2*bx(nx)**2 + a3*bx(nx)**3 - 1.0 = 0
        if (a3.eq.0.d0) then
           if (a2.eq.0.d0) then
              if(a1.ne.0.d0) bx(nx)=1.d0/a1
           else
              bx(nx)= (-a1+dsqrt(a1*a1+4.d0*a2))*0.5d0/a2
           endif
        else
           if (a1.eq.0.d0.and.a2.eq.0.d0) then
              bx(nx)= (1.d0/a3)**(1.d0/3.d0)
           else
              p1=a1/a3
              p2=a2/a3
              p0=-1.d0/a3
              call cubic (p2,p1,p0,z1,z2,z3)
              if(z1.gt.0.d0.and.z1.lt.1.d0.and.z2.gt.0.d0.and.
     +          z2.lt.1.d0.and.dabs(z1-z2)/z1.gt.1.d-3) go to 2000
              if(z1.gt.0.d0.and.z1.lt.1.d0.and.z3.gt.0.d0.and.z3.
     +          lt.1.d0.and.dabs(z1-z3)/z1.gt.1.d-3) go to 2000
              if(z3.gt.0.d0.and.z3.lt.1.d0.and.z2.gt.0.d0.and.z2.
     +          lt.1.d0.and.dabs(z3-z2)/z3.gt.1.d-3) go to 2000
              bx(nx)=z1
              if(z2.gt.0.d0.and.z2.lt.1.d0) bx(nx)=z2
              if(z3.gt.0.d0.and.z3.lt.1.d0) bx(nx)=z3
           endif
        endif
        endif
c
c       gapon convention
c
        if (iex.eq.3) then
        do 500 j=1,nexc
        dum(j)= (1.d0/ekx(j))
     +          *(ct(nbx(j))*gamt(nbx(j)))**(1.d0/z(nbx(j)))
     +          *(ct(nbx(nx))*gamt(nbx(nx)))**(-1.d0/z(nbx(nx)))
        a1=a1+dum(j)
500     continue
        bx(nx)=1.d0/a1
        endif
c
        if (bx(nx).gt.1.d0.or.bx(nx).lt.0.d0) go to 2000
c       the rest of eq. fractions as function of bx(nx)
        do 300 j=1,nexc
           if (iex.eq.1.or.iex.eq.2) then
              bx(j)= dum(j)*bx(nx)**(z(nbx(j))/z(nbx(nx)))
                                     else
              bx(j)= dum(j)*bx(nx)
           end if
300     continue
!
! -------------
!.....Conversion of bx (eq. fraction) into cx (mol solute ads/dm3 sol)
! -------------
!
!     Mod_Xsl = 1  ! Model for exchange dependence on water saturation
!
!.....Simply divide by water saturation
!
      if (Mod_Xsl .eq. 1)       then
!
         cecmol = cec2*denss2*(1.d0-phi2)*1.d-2/phisl2
!
!
!
!.....Modified to account for reduced surface sites: Factor of
!.....S from Liu et al. (1998) causes S in denominator
!.....to drop out (saturated system, also)
!
      else if (Mod_Xsl .eq. 2)   then
!
         if (a_fmr2 .lt. sl2 .and. sl2 .gt. 0.d0 .and.
     &       a_fmr2 .gt. 0.d0)   then
!
!...........Factor based on active fracture model at low saturations,
!...........and for saturations above zero
!
            cecmol = a_fmr2*cec2*denss2*(1.d0-phi2)*1.d-2/phisl2
!                      ------
!
                              else
!
!...........Factor for saturated system, or unsaturated to consider only
!...........the wetted proportion
!
            cecmol = cec2*denss2*(1.d0-phi2)*1.d-2/phi2
!
         end if
!
      end if  ! IF_Mod_Sl
!
! -------------
!
        do 400 j=1,nexc
c       gaines&thomas and gapon conventions
        if (iex.eq.1.or.iex.eq.3) cx(j)=bx(j)*cecmol/z(nbx(j))
c       vanselow convention
        if (iex.eq.2) cx(j)=bx(j)*cecmol
400     continue
c
        return
c
2000    write (*,2100) z1,z2,z3
2100    format(///,1x,3e14.3,'error or ambiguity in ion exchange calc.')
        stop
        end
c
c
c------------------------------------------------------------------------------
c
c
        subroutine cubic (a,b,c,z1,z2,z3)
c
c       real roots of a cubic equation (newton-raphson)
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        z=1.d0  !trial value
40      f=z*z*z+a*z*z+b*z+c
        h=f/(3.d0*z*z+2.d0*a*z+b)
        if(dabs(h/z).le.1.0d-4) go to 90
        z=z-h
        go to 40
90      z1=z
        z=0.01d0  !second root trial
110     f=z*z*z+a*z*z+b*z+c
        h=f/(3.d0*z*z+2.d0*a*z+b)
        if(dabs(h/z).le.1.0d-4) go to 160
        z=z-h
        go to 110
160     z2=z
170     z3=-a-z1-z2
        return
        end
c
c
c
c-----------------------------------------------------------------------------
c
c
c
        subroutine dcx_dcp
c
c****************** Evaluate derivatives for cation exchange ******************
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        double precision cxold(mexc)
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(1x,15('*'),' dcx_dcp  1.0, 2003.7.30: Evaluate'
     X  ' derivatives for cation exchange in geochem_v2.f',15('*'))
c
c       numerical derivative of cx respecte a cp
c
        do 100 j=1,nexc
100        cxold(j)=cx(j)
        do 200 i=1,npri
        dd=cp(i)*1.0d-07
        cp(i)=cp(i)+dd
        ct(i)=cp(i)
        call cx_ct
        do 150 j=1,nexc
        dcx(j,i)=(cx(j)-cxold(j))/dd
        cx(j)=cxold(j)
150     continue
        cp(i)=cp(i)-dd
        ct(i)=cp(i)
c
200     continue
c
        return
        end
c
c
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine ludcmp(a,n,np,indx,d)
c
      implicit double precision(a-h,o-z)
      implicit integer*8 (i-n)
      parameter (nmax=100,tiny=1.0d-20)
c     Added definitions for passed integers and fp vars
      integer*8 n,np
      integer*8 indx(np)
      double precision a(np,np),vv(nmax),d
c
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          aijabs = dabs(a(i,j))
          if (aijabs.gt.aamax) aamax=aijabs
11      continue
c       Added to stop program and write to output file
        if (aamax.eq.0.d0)then
c            write(32,*)'Singular Matrix in Chemical Solver, STOP'
c            write(34,*)'Singular Matrix in Chemical Solver, STOP'
c            call chdump(timetot,ielem,iterch)
c            stop
            d=1.0d+33                        !go return to Newton  gxzh  11/22/05
            return                           !go return to Newton  gxzh  11/22/05
        endif
c......
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.d0)a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.d0)a(n,n)=tiny
      return
      end
c
c
c-------------------------------------------------------------------------------
c
c
        subroutine lubksb(a,n,np,indx,b)
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      integer*8 n,np
      integer*8 indx(np)
      double precision a(np,np),b(np)

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end
c
c
c
c-------------------------------------------------------------------------------
c
c
c
c
       SUBROUTINE WRITE_PLOT
C
C*******************WRITE variables versus space in TECPLOT format************
C
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        include 'perm_v2.inc'
        COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
        COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
        COMMON/E6/T(MNEL)
        COMMON/E5/P(MNEL)
        COMMON/SOLUTE6/SLOLD(MNEL)         ! old liquid saturation
        COMMON/SOLUTE8/SL1(MNEL)           ! new liquid saturation
        COMMON/SOLUTE9/SG1(MNEL)           ! new gas saturation
        COMMON/SOLUTE10/PHIOLD(MNEL)       ! porosity at previous time step
        COMMON/C3/DEL1(MNCON)
        COMMON/C4/DEL2(MNCON)
        COMMON/P1/X((MNK+1)*MNEL)          ! print Rh
        COMMON/XYZ11/XXX(mnel)             ! for TECPLOT
        COMMON/XYZ22/YYY(mnel)             ! for TECPLOT
        COMMON/XYZ33/ZZZ(mnel)             ! for TECPLOT
        double precision PRECIP(mmin),DGP2(mgas)
        COMMON/TRANGAS9/NGAS1                ! Number of gaseous species transported
c
        COMMON/PARNP/NPL,NPG          ! specify in EOS module
        COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
        COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
        COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
        COMMON/DRYOUT/IDRY(MNOD),ADRY(MNOD,MPRI)
        COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)    !residual in precipitates
        double precision CDUM(MAQT)   ! working array for aqueous concentrations
        double precision XDUM(MXsites,MEXC)  ! working array for exchanged concentrations
        double precision AQDUM(MAQT)   ! working array for aqueous species concentrations
        double precision ADDUM(MADS)   ! working array for adsorbed species concentrations
        character*60 form1,form2,form3
c
c              porosity*saturation
        common/phisat/phisl1(mnel),phisg1(mnel)
c
        COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
        COMMON/MOP_REACT/MOPR(20)     ! controling parameters for reactive transport
        COMMON/min_SI/SIM(MNEL,MMIN)  ! Mineral saturation index (log(Q/K)) for all nodes
C
        common/clay_swell1/ iswell
        common/clay_swell2/ vmin_old(mnel,mmin)  ! previous mole volume for all node
        common/clay_swell3/ vmin0(mmin)          ! initial mole volume
        common/water_activity/aw(mnod)           ! water activity
!
      common/Print_Unit_Name/ Name_Conc, Name_Mine
      character*10 Name_Conc, name_unit
      character*45 Name_Mine,astring(msurf)
!
!
!.....Molecular weight of all species, g/mol
      common/molweight/wm_aqt(maqt)
!
      common/minkin4/rkin(mnel,mmin),amin3(mnel,mmin)
      common/afactorr/a_fmr(mnel)
c
c--------------------------------------------------------------------
c
        SAVE ICALL
        DATA ICALL/0/
c
        SAVE AMA,AMS                       !print Rh
        DATA AMS/18.016d0/,AMA/28.96D0/    !print Rh
c
        SAVE HC                            !print Rh
        DATA HC/1.d-10/                    !print Rh
c
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(' WRITE_PLOT 1.1, 2008.2.8: Write variables vs. block'
     X  ' for species, minerals and gases in geochem_v2.f')
c
C-----------------------------------------------------------Write title
        IF (TIMETOT .EQ. 0.0D0.or.icall.eq.1) THEN
C
C-----------------------for aqueous component or species, adsorbed and exch species
C
           WRITE(61,501)
501        FORMAT(/10X,'---- Aqueous concentrations vs. grid blocks',
     +                 ' at specified times ----'/)
           WRITE(61,*)'  Unit: '
c
           WRITE(61,"(4x,'- Aqueous species: Total concencentrations',
     &       ' in ', A10)")   Name_Conc
           name_unit=Name_conc
           if(iconflag.gt.1) name_unit ='mol/L'
           IF (NWAQ .GT. 0)
     &      WRITE(61,"(4x,'- Aqueous species: Individual',
     &        ' concentrations in ',A10)") name_unit
           IF (NWADS .GT. 0)
     &      WRITE(61,"(4x,'- Surface species: Concentrations in ',A10)")
     &         name_unit
           IF (nwexc.GT. 0)
     &      WRITE(61,"(4x,'- Exchanged species: Concentrations in ',
     &         A10)") name_unit
           WRITE(61,"(4x,'- Component in dry nodes in mol/L medium')")
c
           WRITE(61,*)
c--------------------------------------
           WRITE(61,720) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &     (naaqt(iwaq(j)),j=1,nwaq),
     &     (naads(iwads(j)),j=1,nwads),
     &     (('X_'//NAEXC(J), J=1,NWEXC), jj=1,NXsites)            ! Multi-site

           form1 = "(2F11.3,F10.3,3e12.4,F9.5,2F8.4,200E12.4)"    ! default

720   FORMAT('VARIABLES =X,        Y,       Z,  P(bar)  ,      Sg,   ',
     +      '    Sl  ,     T(C)  ,  aH2O,    pH     ',200(',',A11))

           IF (NWMIN .GT. 0)  THEN
              WRITE(63,503)
503           FORMAT(/10X,'---- Changes of mineral abundance (or/and ',
     +           'exchanged species concentrations) vs. grid blocks',
     +                 ' at specified times ----'/)
              WRITE(63,*)'  Unit: '
c
              WRITE(63,27)
27            FORMAT ('               ')
              WRITE(63,15)  Name_mine
15            FORMAT ('    - Mineral: ', A45)
!
              WRITE(63,*)
c--------------------------------------
c
              do j=1,nsurf
                astring(j)='  exp(phi )'
                write(astring(j)(10:10),"(i1)") j
              enddo
              WRITE(63,725) (astring(j), j=1,nsurf),
     &         (NAMIN(IWMIN(J)),J=1,NWMIN)      
725           FORMAT('VARIABLES =X,       Y,        Z,      T,',
     +            ' Porosity,   Permeabi., ',65(A11,','))
           END IF
C
C-----------------------for gases
C
           IF (NGAS .GT. 0)   THEN
              WRITE(66,506)
506           FORMAT(/10X,'---- Gas vs. grid blocks',
     +                 ' at specified times ----'/)
              WRITE(66,*)'  Unit: '
              IF (IEOS.EQ.9)   THEN
                 WRITE(66,1004)
1004             FORMAT ('    - Gas: partial pressure (bar)')
                 GOTO 799
              END IF
              IF (MINFLAG .GE. 1)  THEN
                 WRITE(66,1013)
1013             FORMAT ('    - Gas: volume fraction')
              ELSE
                 WRITE(66,1014)
1014             FORMAT ('    - Gas: partial pressure (bar)')
              END IF
799           CONTINUE
              WRITE(66,*)
              IF (IEOS.NE.9)  THEN
                 WRITE(66,727) (NAGAS(J),      J=1,NGAS)
727              FORMAT('VARIABLES =X,        Y,       Z,   T(C),',
     +      '    Sg ,    RH,    P(bar)   ',12(',',A11))
              ELSE
                 WRITE(66,728) (NAGAS(J),      J=1,NGAS)
728              FORMAT('VARIABLES =X,        Y,       Z,   T(C),',
     +      '    Sg ,    RH    ',20(',',A11))
              END IF
           END IF
C
C-----------------------for mineral saturation index (optional)
C          changed so that mopr ge 1 prints SI and all others dropped by one
c          then mass balance only printed for mopr(8) ge 4
           IF (MOPR(8).ge.1)   THEN
              WRITE(77,601)
601           FORMAT(/10X,'---- Mineral saturation index (log(Q/K))',
     +                   ' vs. grid blocks at specified times ----'/)
              WRITE(77,*)
              WRITE(77,603) (NAMIN(J),J=1,NMIN)
603           FORMAT('VARIABLES =X,       Y,      Z ',7X,60(',',A11))
           END IF
C-----------------------for mineral reaction rate (optional)
c
           IF (MOPR(8).ge.2)   THEN
              WRITE(78,605)
605          FORMAT(/10X,'---- Reaction Rate (mol/kg H2O/s)',
     +                   ' vs. grid blocks at specified times ----'/)
              WRITE(78,*)
              WRITE(78,607) (NAMIN(J),J=nmequ+1,nmin)
607           FORMAT('VARIABLES =X,       Y,      Z ',7X,60(',',A11))
           END IF
C
C---- active fracture area factor and reactive surface areas (optional)
           IF (MOPR(8).ge.3)   THEN
              WRITE(79,609)
609          FORMAT(/10X,'---- Reactive Surface Areas (m^2/ kg H2O)',
     +                   ' vs. grid blocks at specified times ----'/)
              WRITE(79,*)
              WRITE(79,611) (NAMIN(J),J=nmequ+1,nmin)
611    FORMAT('VARIABLES =X,       Y,      Z,     Sl,      A_fmr,',
     +     '     A_Factor  ',60(',',A11))
           END IF
C
        END IF
C
        timeday = timetot/8.6400d4
        WRITE (61,730) TIMEDAY
        IF (NWMIN .GT. 0) WRITE (63,730) TIMEDAY
        IF (NGAS .GT. 0) WRITE (66,730) TIMEDAY
730     FORMAT('ZONE T= "',E12.6,' d"','  F=POINT')
C
C       Changed so that mbalance moves to 4, all others dropped by one
c       also should print at time zero
        IF (MOPR(8).ge.1 .AND. TIMETOT.GT.0.0D0)   THEN
           WRITE (77,739) TIMEDAY
739        FORMAT('ZONE T= "',E12.6,' d"','  F=POINT')
        END IF

        IF (MOPR(8).ge.2 .AND. TIMETOT.GT.0.0D0)   THEN
           WRITE (78,730) TIMEDAY
        END IF
C
        IF (MOPR(8).ge.3 .AND. TIMETOT.GT.0.0D0)   THEN
           WRITE (79,730) TIMEDAY
        END IF
C
C------Write results (Sl,pH,logfo2,dissolved con.,precip.,exchange and adsorption)
C

       if(mopr(7).eq.0.or.mopr(7).eq.4)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E12.4)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E12.4)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E12.4)"
       elseif(mopr(7).eq.1)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E9.1)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E9.1)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E9.1)"
       elseif(mopr(7).eq.2)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E10.2)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E10.2)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E10.2)"
       elseif(mopr(7).eq.3)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E11.3)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E11.3)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E11.3)"
       elseif(mopr(7).eq.5)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E12.5)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E13.5)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E12.5)"
       elseif(mopr(7).eq.6)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E13.6)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E14.6)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E13.6)"
       elseif(mopr(7).eq.7)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E14.7)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E15.7)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E14.7)"
       elseif(mopr(7).ge.8)then
           form1 = "(2F11.3,F10.3,3e12.4,F8.3,2F8.4,200E15.8)"
           form2 = "(2F11.3,F10.3,F8.3,F10.5,e13.5,60E16.8)"
           form3 = "(2F11.3,F10.3,F8.3,E12.4,F8.3,20E15.8)"
       endif
c
          DO 740 I=1,NNOD
c
c-----calculate current relative humidity (rlwr)-------------------
c
          nloc=(i-1)*nk1                  ! print Rh
          nloc2=(i-1)*neq1*nsec
          nloc2l=nloc2+nbk
          pres=x(nloc+1)
          tx=par(nloc2+nsec-1)

          sgwr=par(nloc2+1)
          slwr=1.0d0-sgwr
c
c..........initialize rlwr
          rlwr = 1.0d0
c
        if(nk.eq.2) then
           xairl=par(nloc2l+nb+2)
           xairg=par(nloc2+nb+2)
           call sat(tx,ps)
          if(sgwr.gt.0.d0.and.slwr.gt.0.d0)then
            xmol=(xairl/ama)/(xairl/ama+(1.d0-xairl)/ams)
            pa=xmol/hc
            rlwr=(pres-pa)/ps
          else                                    !if no liquid water present
            xairv=1.d0-xairg                      !h2o gas mass fraction
            xmolv=xairv*ama/(ama*xairv+ams*xairg)  !h2o gas mole fraction
            rlwr=pres*xmolv/ps
          endif
         endif                                 !!print Rh
c
         densw=dwat(i)                ! water density (kg/l)
c
         IF (ICONFLAG .GE. 1 )               THEN
c     ... output in mol, g, or mg/L (for cdum only)
c
            IF(IDRY(I).ge.2) THEN
             DO J=1,NPRI
                CDUM(J)=ADRYR(I,J)
             END DO
            ELSE
             DO J=1,NPRI
                CDUM(J)=ctot(i,j)*dwat(i)
             END DO
            ENDIF
c
c           ---------------------
c           Use additional separate arrays for printing of secondary species,
c           sorbed species, and exchange species
            do j=1,nwaq
               aqdum(j)=c(i,iwaq(j))*densw
            enddo
            do j=1,nwads
               addum(j)=d(i,iwads(j))*densw
            enddo
            do j=1,nwexc
                do isite=1,NXsites
                  XDUM(isite,J)=XCADS(I,isite,iwexc(j))*densw
              enddo
            enddo
c           ---------------------
            if (iconflag .ge. 2)   then           ! converts cdum only
               do j=1,npri
                 if (napri(j) .ne. 'tracer') then
                   CDUM(J)=CDUM(J)*wm_aqt(j)   ! g/L
                   if (iconflag .eq. 3)   then
                     CDUM(J)=CDUM(J)*1.0d3    ! mg/L (ppm)
                   end if
                 endif
               end do
            end if

         ELSE
c     ..... output in mol/kg
            IF(IDRY(I).ge.2) THEN
               DO J=1,NPRI
                  CDUM(J)=ADRYR(I,J)
               END DO
            ELSE
               DO J=1,NPRI
                  CDUM(J)=ctot(i,j)       !UTOLD(I,J)*dwinv   ctot is in mol/kgw
               END DO
            ENDIF
c
c-------------------
c           NS3/08 use separate arrays for printing of secondary species,
c           sorbed species, and exchange species
            do j=1,nwaq
               aqdum(j)=c(i,iwaq(j))
            enddo
            do j=1,nwads
                addum(j)=d(i,iwads(j))
            enddo
            do j=1,nwexc
                do isite=1,NXsites
                  XDUM(isite,J)=XCADS(I,isite,iwexc(j))
                enddo
            enddo

         END IF
c--------------------
c
         if (iswell.eq.1) then      ! iswell=1   for clay mineral swelling
            do m = 1, nmin
               vmin(m)=vmin_old(i,m)
            end do
         end if
c
         if(minflag.eq.1 .or. minflag.eq.3) then            
           DO J=1,NWMIN
             IF (ISWELL .NE. 1)    THEN
                PRECIP(J)=PRE(I,IWMIN(J))-PINIT(I,IWMIN(J))  ! pre in moles per liter medium
                PRECIP(J)=PRECIP(J)*vmin(iwmin(j))           ! change in volume fraction
                                 ELSE
                PRECIP(J)=PRE(I,IWMIN(J))*vmin(iwmin(j))-
     +                  PINIT(I,IWMIN(J))*vmin0(iwmin(j))    ! for clay swelling
             END IF
!
             if (minflag .eq.3 ) PRECIP(J)=PRECIP(J)*100.0d0 ! Change in Vf %
!
           END DO
         else if(minflag.eq.2)then
           DO J=1,NWMIN
             PRECIP(J)=PRE(I,IWMIN(J))*vmin(iwmin(j))        ! total volume fraction
           END DO
         else
           DO J=1,NWMIN
             PRECIP(J)=PRE(I,IWMIN(J))-PINIT(I,IWMIN(J))  !pre in moles per liter medium
             PRECIP(J)=PRECIP(J)*1000.0D0    ! change in mol/m**3 medium
           END DO
         endif
C
         DO 330 J=1,NGAS
c....... for minflag >= 1, volume fraction, partial pressure (bar) otherwise
         if(minflag.ge.1)then
           dgp2(j)=pfug(i,j)*1.d5/p(i)
         else
           dgp2(j)=pfug(i,j)
         endif
         IF (ieos.eq.9)   dgp2(j)=pfug(i,j)
330          CONTINUE
c
              WRITE (61,form1)XXX(I),YYY(I),ZZZ(I),P(I)/1.0D+5,sg1(I),
     &              sl1(I),T(I),AW(I),PH(I),(CDUM(IWCOM(J)),J=1,NWCOM),
     &                    (aqdum(J), j=1,nwaq),
     &                    (addum(J), j=1,nwads),
     &                    ((XDUM(isite,J), J=1,NWEXC),isite=1,NXsites)
           if (nwmin.ne.0) then
             if (kcpl.ne.2)   then
              WRITE (63,form2)XXX(I),YYY(I),ZZZ(I),T(I),
     &          phi(i),perm(2,i),(dexp(phip(i,j)),j=1,nsurf),
     &        (PRECIP(J),J=1,NWMIN)       !,(D(I,J), J=1,NADS)
c     &                    ((XDUM(isite,J), J=1,NWEXC),isite=1,NXsites)
                            else         ! only monitor porosity and Per. changes
              WRITE (63,form2)XXX(I),YYY(I),ZZZ(I),T(I),
     &        phim(i),permm(2,i),(dexp(phip(i,j)),j=1,nsurf),
     &        (PRECIP(J),J=1,NWMIN)       !,(D(I,J),J=1,NADS)
c     &                    ((XDUM(isite,J), J=1,NWEXC),isite=1,NXsites)
             end if
           end if
c
        IF (IEOS.NE.9)  THEN
           IF (NGAS1 .NE. 0)  THEN
c              WRITE (66,form3)XXX(I),YYY(I),ZZZ(I),T(I),sg1(I),rlwr,
c     +                     P(I)/1.0D+5,(PFUG(I,J),          J=1,NGAS)
c                           ELSE
              WRITE (66,form3)XXX(I),YYY(I),ZZZ(I),T(I),sg1(I),rlwr,
     +                     P(I)/1.0D+5,(DGP2(J),          J=1,NGAS)
           END IF
                        ELSE
              WRITE (66,form3)XXX(I),YYY(I),ZZZ(I),T(I),sg1(I),
     +                     rlwr,(PFUG(I,J),          J=1,NGAS)
c
        END IF
C
C-----------------------for mineral saturation index (optional)
       IF (MOPR(8).ge.1 .AND. TIMETOT.GT.0.0D0)   THEN
          WRITE (77,'(3F11.3,60F12.6)') XXX(I),YYY(I),ZZZ(I),
     +                              (SIM(I,J),J=1,NMIN)
       END IF
C--------------------------------------------------------------
C-----------------------for Reaction rate (optional)
       IF (MOPR(8).ge.2 .AND. TIMETOT.GT.0.0D0)   THEN
          WRITE (78,'(3F11.3,60e12.4)') XXX(I),YYY(I),ZZZ(I),
     +                   (rkin(I,J),J=1,nmkin)
       END IF
C-----------------------for Reactive surface areas (optional)
       IF (MOPR(8).ge.3 .AND. TIMETOT.GT.0.0D0)   THEN
c
c         Modified to account for reduced surface area: Factor of
c...          S from Liu et al. (1998) causes S in denominator
c...          to drop out (saturated system, also)
         if(a_fmr(i).lt.sl1(i).and.sl1(i).gt.0.d0.and.
     +      a_fmr(i).gt.0.d0)then
c..... Factor based on active fracture model at low saturations,
c....... and for saturations above zero
           actfrc = a_fmr(i)/(phisl1(i)*densw*1.d3)
c        Need to skip for dry nodes
         elseif(sl1(i).eq.1.d0.or.a_fmr(i).gt.sl1(i).and.
     +     sl1(i).ne.0.d0)then
c..... Factor for saturated system, or unsaturated to consider only
c........ the wetted proportion
           actfrc = 1.d0/(phi(i)*densw*1.d3)
         elseif(sl1(i).eq.0.d0)then
           actfrc=0.d0
         endif
c
          WRITE (79,'(3F11.3,3e12.5,60e12.4)') XXX(I),YYY(I),ZZZ(I),
     +         sl1(i),a_fmr(i),actfrc,(amin3(I,J),J=1,nmkin)
       END IF
C
740      CONTINUE
c
       RETURN
       END
c
c-------------------------------------------------------------------------------
c
       SUBROUTINE WRITE_PLOT_ECO2
C
C                                  For ECO2N flow module
C
C*******************WRITE varibles versus space in TECPLOT format************
C
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        include 'perm_v2.inc'
        COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
        COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
        COMMON/E3/EVOL(MNEL)
        COMMON/E6/T(MNEL)
        COMMON/E5/P(MNEL)
        COMMON/E1/ELEM(MNEL)
        CHARACTER*5 ELEM
        COMMON/SOLUTE6/SLOLD(mnel)      ! old liquid saturation
        COMMON/SOLUTE8/SL1(mnel)        ! new liquid saturation
        COMMON/SOLUTE9/SG1(mnel)        ! new gas saturation
        COMMON/SOLUTE10/PHIOLD(mnel)    ! porosity at previous time step
        COMMON/C3/DEL1(MNCON)
        COMMON/C4/DEL2(MNCON)
        COMMON/XYZ11/XXX(mnel)          ! for TECPLOT
        COMMON/XYZ22/YYY(mnel)          ! for TECPLOT
        COMMON/XYZ33/ZZZ(mnel)          ! for TECPLOT
        double precision PRECIP(mmin),DGP2(mgas)
        COMMON/TRANGAS9/NGAS1           ! Number of gaseous species transported
        COMMON/PARNP/NPL,NPG            ! specify in EOS module
        COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
        COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
!
        double precision CDUM(MAQT)   ! working array for aqueous concentrations
        double precision AQDUM(MAQT), ADDUM(MADS)  ! working array for aq and ads concentrations
        double precision XDUM(MXsites,MEXC)        ! working array for exchanged concentrations
!
        common/clay_swell1/ iswell
        common/clay_swell2/ vmin_old(mnel,mmin)  ! previous mole volume for all node
        common/clay_swell3/ vmin0(mmin)          ! initial mole volume
!
!              porosity*saturation
        common/phisat/phisl1(mnel),phisg1(mnel)
c
c----------------------------------For co2 sequstration
        common/co2_gene/ nco2
        common/co2_gene1/ nco2g
        COMMON/SOLIDco2/SMco2(NMNOD)     ! CO2 TRAPPED in solid phase
c
c.....................Added arrays for printing surface areas and factors
c
        common/minkin4/rkin(mnel,mmin),amin3(mnel,mmin)
        common/afactorr/a_fmr(mnel)
c
C----------------------------------------------------------------
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
!.....Extract fugacity coefficients dissolved gas concentrations from TMgas
      common/fugacity_coe/fug_coe(mnel,18)
      common/gas_index/ichem(18)            ! No-cond gas index in chemical input
!
!.....Extracting CO2 fugacity coefficient from ECO2N
      common/fuga_coe /FugCoeCO2(mnel)
!
!.....in the BMW array molar weights are ordered as in PAR !
      COMMON /MOLWtbio/ BMW(20)
!
      double precision CmolKg(20)       ! Dissolved gas concentration from TMgas
      double precision Xgg(20)          ! Gas mass frcation in the gas phase for TMVOCs
c
c----------------------------------------------------------------
c
        COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
        COMMON/min_SI/SIM(MNEL,MMIN) ! Mineral saturation index (log(Q/K)) for all nodes
        COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
!.....Molecular weight of all species, g/mol
      common/molweight/wm_aqt(maqt)!
!
      common/dissolved_solid/ TDS(mnel)             !!!! for use in EOS
      common/Print_Unit_Name/ Name_Conc, Name_Mine
      character*10 Name_Conc, name_unit
      character*45 Name_Mine
!
C----------------------------------------------------------------
c
        character*42 form1,form2,form3,astring(msurf)
C
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(11,899)
  899   FORMAT(' WRITE_PLOT_ECO2 1.0 30 July 2003',6X,
     X  'Write variables vs. block for species, minerals and gases')
C-----------------------------------------------------------Write title
        IF (TIMETOT .EQ. 0.0D0.or.icall.eq.1) THEN
c
c-----------------------for aqueous component or species, adsorbed and exch species
c         
           WRITE(61,501)
501        FORMAT(/10X,'---- Aqueous concentrations vs. grid blocks',
     +                 ' at specified times ----'/)
           WRITE(61,*)
           WRITE(61,*)'  Unit: '

           WRITE(61,"(4x,'- Aqueous species: Total concencentrations',
     &       ' in ', A10)")   Name_Conc
           name_unit=Name_conc
           if(iconflag.gt.1) name_unit ='mol/L'
           IF (NWAQ .GT. 0)
     &      WRITE(61,"(4x,'- Aqueous species: Individual',
     &        ' concentrations in ',A10)") name_unit
           IF (NWADS .GT. 0)
     &      WRITE(61,"(4x,'- Surface species: Concentrations in ',A10)")
     &         name_unit
           IF (nwexc.GT. 0)
     &      WRITE(61,"(4x,'- Exchanged species: Concentrations in ',
     &         A10)") name_unit
           WRITE(61,"(4x,'- Component in dry nodes in mol/L medium')")
           WRITE(61,*)
c--------------------------------------
           WRITE(61,720) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &        (naaqt(iwaq(j)),j=1,nwaq),
     &        (naads(iwads(j)),j=1,nwads),
     &        (('X_'//NAEXC(J), J=1,NWEXC), jj=1,NXsites)                ! Multi-site
720           FORMAT('VARIABLES =X,       Y,      Z,      Sg,',
     +             '     Sl,         T,        pH,   Density_L',
     +          '       ',200(A11,','))
C
C-----------------------for minerals
C

c
           IF (NWMIN .GT. 0)  THEN
              WRITE(63,503)
503           FORMAT(/10X,'---- Changes of mineral abundance (or/and ',
     +           'exchanged species concentrations) vs. grid blocks',
     +                 ' at specified times ----'/)
              WRITE(63,*)
              WRITE(63,*)'  Unit: '
C
              WRITE(63,26)
26            FORMAT ('    - SMco2: Total CO2 sequstered in',
     +                ' mineral phases in  kg/m**3 medium')
              WRITE(63,27)
27            FORMAT ('    - Permeability in  m**2')
C
              WRITE(63,"('    - Minerals: ', A45)")  Name_mine
              WRITE(63,*)
c--------------------------------------
c
              do j=1,nsurf
                astring(j)='  exp(phi )'
                write(astring(j)(10:10),"(i1)") j
              enddo
              WRITE(63,725) (astring(j), j=1,nsurf),
     &            (NAMIN(IWMIN(J)),J=1,NWMIN)
725           FORMAT('VARIABLES =X,    Y,        Z,      T,   SMco2 ,',
     +        ' Porosity,   Permeabi., ',65(A11,','))
           END IF
C-----------------------for gases
           IF (NGAS .GT. 0)   THEN
              WRITE(66,506)
506           FORMAT(/10X,'---- Gas pressures (bar), fugacity',
     &        ' coefficients, and  dissolved concentrations',
     &        ' (mol/kg) ',
     &        /'              vs. grid blocksat specified times ----'/)
           end if
           if (ieos .eq. 15 .or. ieos .eq. 16)   then      ! For TMVOCs
              WRITE(66,728) (NAGAS(J),J=1,NGAS)
728           FORMAT('VARIABLES =X,        Y,       Z,      T,',
     &         '  ',20(A11,'Fcoe         Conc       Xgg         ,'))
                                                 else      ! For ECO2
727           FORMAT('VARIABLES =X,        Y,       Z,      T,',
     &         '  ',20(A11,'Fcoe         Conc       ,'))
           end if
c
c-----------------------for mineral saturation index (optional)
c
           IF (MOPR(8).ge.1)   THEN
              WRITE(77,601)
601           FORMAT(/10X,'---- Mineral saturation index (log(Q/K))',
     +                   ' vs. grid blocks at specified times ----'/)
              WRITE(77,*)
              WRITE(77,603) (NAMIN(J),J=1,NMIN)
603           FORMAT('VARIABLES =X,       Y,      Z,',7X,60(A11,','))
           END IF
C
C-----------------------for mineral reaction rate (optional)
c
           IF (MOPR(8).ge.2)   THEN
              WRITE(78,605)
605          FORMAT(/10X,'---- Reaction Rate (mol/kg H2O/s)',
     +                   ' vs. grid blocks at specified times ----'/)
              WRITE(78,*)
              WRITE(78,607) (NAMIN(J),J=nmequ+1,nmin)
607           FORMAT('VARIABLES =X,       Y,      Z,',7X,60(A11,','))
           END IF
C
C---- active fracture area factor and reactive surface areas (optional)
           IF (MOPR(8).ge.3)   THEN
              WRITE(79,609)
609          FORMAT(/10X,'---- Reactive Surface Areas (m^2/ kg H2O)',
     +                   ' vs. grid blocks at specified times ----'/)
              WRITE(79,*)
              WRITE(79,611) (NAMIN(J),J=nmequ+1,nmin)
611    FORMAT('VARIABLES =X,       Y,      Z,     Sl,      A_fmr,',
     +     '     A_Factor,  ',60(A11,','))
           END IF
C
        END IF
!
!       Initialize all printing arrays to 0
!
        do j=1,npri
             cdum(j)=0.d0
        enddo
        do j=1,nwaq
             aqdum(j)=0.d0
        enddo
        do j=1,nwads
             addum(j)=0.d0
        enddo
        do j=1,nwexc
          do isite=1,NXsites
                XDUM(isite,J)=0.d0
          enddo
        enddo
C
C----------------For co2 disposal
          TAco2=0.0d0
          TGco2=0.0d0
          TSco2=0.0d0
          aco2 = 0.d0
          tco2 = 0.d0
c
C-----------------------------------------------------
c
        timeday = timetot/8.6400d4
C
        WRITE (61,730) TIMEDAY
        IF (NWMIN .GT. 0) WRITE (63,730) TIMEDAY
        IF (NGAS .GT. 0) WRITE (66,730) TIMEDAY
730     FORMAT('ZONE T= "',E12.6,' d"','  F=POINT')
C
        IF (MOPR(8).ge.1 .AND. TIMETOT.GT.0.0D0)   THEN
           WRITE (77,739) TIMEDAY
739        FORMAT('ZONE T= "',E12.6,' d"','  F=POINT')
        END IF
C
        IF (MOPR(8).ge.2 .AND. TIMETOT.GT.0.0D0)   THEN
           WRITE (78,730) TIMEDAY
        END IF
C
        IF (MOPR(8).ge.3 .AND. TIMETOT.GT.0.0D0)   THEN
           WRITE (79,730) TIMEDAY
        END IF
C
C------Write results (Sl,pH,logfo2,dissolved con.,precip.,exchange and adsorption)
C
       if(mopr(7).eq.0.or.mopr(7).eq.4)then
           form1 = "(2F11.3,F10.3,2F8.3,200E12.4)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E12.4)"
           form3 = "(2F11.3,F10.3,F8.3,20E12.4)"
       elseif(mopr(7).eq.1)then
           form1 = "(2F11.3,F10.3,2F8.3,200E9.1)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E9.1)"
           form3 = "(2F11.3,F10.3,F8.3,20E9.1)"
       elseif(mopr(7).eq.2)then
           form1 = "(2F11.3,F10.3,2F8.3,200E10.2)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E10.2)"
           form3 = "(2F11.3,F10.3,F8.3,20E10.2)"
       elseif(mopr(7).eq.3)then
           form1 = "(2F11.3,F10.3,2F8.3,200E11.3)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E11.3)"
           form3 = "(2F11.3,F10.3,F8.3,20E11.3)"
       elseif(mopr(7).eq.5)then
           form1 = "(2F11.3,F10.3,2F8.3,200E12.5)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E13.5)"
           form3 = "(2F11.3,F10.3,F8.3,20E12.5)"
       elseif(mopr(7).eq.6)then
           form1 = "(2F11.3,F10.3,2F8.3,200E13.6)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E14.6)"
           form3 = "(2F11.3,F10.3,F8.3,20E13.6)"
       elseif(mopr(7).eq.7)then
           form1 = "(2F11.3,F10.3,2F8.3,200E14.7)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E15.7)"
           form3 = "(2F11.3,F10.3,F8.3,20E14.7)"
       elseif(mopr(7).ge.8)then
           form1 = "(2F11.3,F10.3,2F8.3,200E15.8)"
           form2 = "(2F11.3,F10.3,F8.3,2F10.5,E13.5,60E16.8)"
           form3 = "(2F11.3,F10.3,F8.3,20E15.8)"
       endif
C
          DO 740 I=1,NNOD
c---------------EPA-CO2, Jaffe
cc          if (elem(i)(1:3).ne.'A11')   go to 740                                            
cc            if (elem(i)(4:5).eq.'50')   go to 740                       
c
C---------------------Aqueous CO2 sequstrated
           Aco2=Aco2*1000.0D0        !--->  mol/m**3
           Aco2=Aco2*44.0D0          !--->  g/m**3
           Aco2=Aco2/1000.0D0        !--->  kg/m**3
C---------------------Total (aqueous+solid) CO2 sequstrated
           Tco2=Tco2*1000.0D0        !--->  mol/m**3
           Tco2=Tco2*44.0D0          !--->  g/m**3
           Tco2=Tco2/1000.0D0        !--->  kg/m**3
c---------------------Gaseous CO2 sequestrated
c
             NLOC2=(I-1)*NSEC*NEQ1
             NLOC2L=NLOC2+NBK
c
             Dco2=PAR(NLOC2+4)           ! CO2 density
             Sco2=1.0D0-SLOLD(I)         ! gas co2 saturation
             Gco2=Dco2*PHI(I)*Sco2
c
c---------------------aqueous CO2 sequestrated
c
             Dliqc2=PAR(NLOC2L+4)           ! CO2 density
             qM=Dliqc2*PHI(I)*SLOLD(I)      ! liquid mass in Kg/m3
             Xco2l=par(nloc2+nbk+nb+3)      ! mass fraction of CO2 in liquid phase
             Aco2=Xco2l*qM                  ! aqueous co2 sequestrated
c----------------------------------------------------------
             TGco2=TGco2+Gco2*evol(i)
             TAco2=TAco2+Aco2*evol(i)
             TSco2=TSco2+SMco2(I)*evol(i)
C
c........Modified for printout of aqueous concentrations
c
         NLOC2=(I-1)*NSEC*NEQ1
         NL2NP=NLOC2+(NPL-1)*NBK
         dlkgl = dwat(i)
c
         if(dlkgl.eq.0.d0)then
            ph22  = 0.0d0
         else
            ph22  = ph(i)
         endif
c
         IF (ICONFLAG.ge.1)  THEN
c    ...... output in mol, g, or mg/L

            DO J=1,NPRI
               if(sl1(i).ne.0.d0) CDUM(J)=ctot(I,J)*dlkgl
            END DO
c-------------------
c           Use additional separate arrays for printing of secondary species,
c           sorbed species, and exchange species (all in molalities)
            do j=1,nwaq
                if(sl1(i).ne.0.d0) aqdum(j)=c(i,iwaq(j))*dlkgl
            enddo
            do j=1,nwads
                if(sl1(i).ne.0.d0) addum(j)=d(i,iwads(j))*dlkgl
            enddo
            do j=1,nwexc
                 do isite=1,NXsites
                   if(sl1(i).ne.0.d0)
     &                  XDUM(isite,J)=XCADS(I,isite,iwexc(j))*dlkgl
                 enddo
            enddo
c           --------
            if (iconflag .ge. 2)   then        ! converts cdum only
               do j=1,npri
                 if (napri(j) .ne. 'tracer') then
                   CDUM(J)=CDUM(J)*wm_aqt(j)   ! g/L
                   if (iconflag .eq. 3)   then
                     CDUM(J)=CDUM(J)*1.0d3     ! mg/L (ppm)
                   end if
                 endif
               end do
            end if

         ELSE
c    ... output in mol/kg
c
            DO J=1,NPRI
               if(sl1(i).ne.0.d0) CDUM(J)=ctot(I,J)       ! ctot is already in mol/kgw
            END DO
c
c           -------------------
c           NS3/08 use additional separate arrays for printing of secondary species,
c           sorbed species, and exchange species (all in molalities)
            do j=1,nwaq
               if(sl1(i).ne.0.d0) aqdum(j)=c(i,iwaq(j))
            enddo
            do j=1,nwads
               if(sl1(i).ne.0.d0) addum(j)=d(i,iwads(j))
            enddo
            do j=1,nwexc
                 do isite=1,NXsites
                  if(sl1(i).ne.0.d0)
     &                  XDUM(isite,J)=XCADS(I,isite,iwexc(j))
                 enddo
            enddo

         END IF
!
!-----------
!........Extract fugacity coefficients
!........Extract dissolved CO2 from flow modules
!-----------
!
         do ig = 1, ngas
!
!..........For TMVOC1 (regular TMVOC) module
!
           if (ieos .eq. 15)   then    ! Now for TMVOC
!
              icg = ichem(ig)        ! No-cond gas Index in chemical input file
              gamg(ig) = fug_coe(i,icg)
!                        fug_coe(i,icg) is from EOS_TMgas in different gas order
!
           end if
!
!
!..........For TMVOC2 module
!
!
           if (ieos .eq. 16)   then
!
              icg = ichem(ig)        ! No-cond gas Index in chemical input file
              gamg(ig) = fug_coe(i,icg)
!                        fug_coe(i,icg) is from EOS_TMgas in different gas order
!
              Xh2oL = PAR(NL2NP + nb + 1)    ! H2O mass fraction in aqueous phase
              icg        = ichem(ig)         ! No-cond gas Index in chemical input file
              XicgL      = PAR(NL2NP + nb+1 +icg)       ! Mass fraction (First is H2O)
              CmolKg(ig) = 1.0d03*XicgL/bmw(icg)/Xh2oL  ! Dissolved gas concentration (mol/kg)
              if (dlkgl .eq. 0.d0)  CmolKg(ig) = 0.0d0
!
           end if
!
!..........For EOS2, ECO2, and ECO2N  modules
!
           if (ieos .eq. 2  .or. ieos .eq. 13
     &        .or. ieos .eq. 14      )              then
!
              gamg(ig) = fugCoeCO2(i)
!
              Xh2oL      = PAR(NL2NP + nb + 1) ! H2O mass fraction in the aqueous phase
              Xco2L      = PAR(NL2NP + nb + 3) ! CO2 Mass fraction
              CmolKg(ig) = 1.0d03*Xco2L/44.0d0/Xh2oL  ! Dissolved gas concentration (mol/kg)
              if (dlkgl .eq. 0.d0)  CmolKg(ig) = 0.0d0
!
           end if
!
         end do
!
!-----------
!........Extract gas mass fraction in the gas phase
!-----------
!
         if (ieos .eq. 15 .or. ieos .eq. 16)   then    ! for TMVOCs
            NLOC2g = (I-1) * nsec * NEQ1 + nb                
            do ig = 1, ngas
               icg = ichem(ig)      
               Xgg(ig) = PAR(NLOC2g+1+icg)
            end do
         end if
!
!-----------
!........Calculate for printout mineral options
!-----------
!
         if(minflag.eq.1 .or. minflag.eq.3) then           
           DO J=1,NWMIN
             IF (ISWELL .NE. 1)    THEN
                PRECIP(J)=PRE(I,IWMIN(J))-PINIT(I,IWMIN(J))  ! pre in moles per liter medium
                PRECIP(J)=PRECIP(J)*vmin(iwmin(j))           ! change in volume fraction
                                 ELSE
                PRECIP(J)=PRE(I,IWMIN(J))*vmin(iwmin(j))-
     +                  PINIT(I,IWMIN(J))*vmin0(iwmin(j))    ! for clay swelling
             END IF
!
             if (minflag .eq.3 ) PRECIP(J)=PRECIP(J)*100.0d0 ! change in Vf %
!
           END DO
         else if(minflag.eq.2)then
           DO J=1,NWMIN
             PRECIP(J)=PRE(I,IWMIN(J))*vmin(iwmin(j))        ! total volume fraction
           END DO
         else
           DO J=1,NWMIN
             PRECIP(J)=PRE(I,IWMIN(J))-PINIT(I,IWMIN(J))  ! pre in moles per liter medium
             PRECIP(J)=PRECIP(J)*1000.0D0                 ! change in mol/m**3 medium
           END DO
         end if
!
!-----------
!.............................................................................
!-----------
!
             DO 330 J=1,NGAS
             dgp2(j)=pfug(i,j)
330          CONTINUE
C
c              if(yyy(i).eq.5.0d0) then
              WRITE (61,form1)XXX(I),YYY(I),ZZZ(I),sg1(I),sl1(I),
c     +                 T(I),PH(I),Pco2,(CDUM(IWCOM(J)),J=1,NWCOM)
     +                 T(I),pH22,dwat(i), (CDUM(IWCOM(J)),J=1,NWCOM),
     &                    (aqdum(J), j=1,nwaq),
     &                    (addum(J), j=1,nwads),
     &                    ((XDUM(isite,J), J=1,NWEXC),isite=1,NXsites)
c
           if (kcpl.ne.2)   then
              WRITE (63,form2)XXX(I),YYY(I),ZZZ(I),T(I),SMco2(I),
     +          phi(i),perm(2,i),(dexp(phip(i,j)),j=1,nsurf),
     +           (PRECIP(J),J=1,NWMIN)        !,(D(I,J),J=1,NADS),
c     &                    ((XDUM(isite,J), J=1,NEXC),isite=1,NXsites)
                             else     ! only monitor porosity and perm. changes
              WRITE (63,form2)XXX(I),YYY(I),ZZZ(I),T(I),SMco2(I),
     +         phim(i),permm(2,i),(dexp(phip(i,j)),j=1,nsurf),
     +         (PRECIP(J),J=1,NWMIN)         !,(D(I,J),J=1,NADS),
c     &                    ((XDUM(isite,J), J=1,NEXC),isite=1,NXsites)
            end if
c            endif
!
              if (ieos .eq. 15 .or. ieos .eq. 16)   then          ! For TMVOCs
                 WRITE (66,form3)XXX(I),YYY(I),ZZZ(I),T(I),       
     &                   (DGP2(J), gamg(j), CmolKg(j),Xgg(j), J=1,NGAS)
                                                    else
                 WRITE (66,form3)XXX(I),YYY(I),ZZZ(I),T(I),       
     &                   (DGP2(J), gamg(j), CmolKg(j), J=1,NGAS)
              end if
!
!-----------------------for mineral saturation index (optional)
!
       IF (MOPR(8).ge.1 .AND. TIMETOT.GT.0.0D0)   THEN
          WRITE (77,'(3F11.3,60F12.6)') XXX(I),YYY(I),ZZZ(I),
     +                   (SIM(I,J),J=1,NMIN)
       END IF
C-----------------------for Reaction rate (optional)
       IF (MOPR(8).ge.2 .AND. TIMETOT.GT.0.0D0)   THEN
          WRITE (78,'(3F11.3,60e12.4)') XXX(I),YYY(I),ZZZ(I),
     +                   (rkin(I,J),J=1,nmkin)
       END IF
C-----------------------for Reactive surface areas (optional)
       IF (MOPR(8).ge.3 .AND. TIMETOT.GT.0.0D0)   THEN
c
c         Modified to account for reduced surface area: Factor of
c...          S from Liu et al. (1998) causes S in denominator
c...          to drop out (saturated system, also)
         if(a_fmr(i).lt.sl1(i).and.sl1(i).gt.0.d0.and.
     +      a_fmr(i).gt.0.d0)then
c..... Factor based on active fracture model at low saturations,
c....... and for saturations above zero
           actfrc = a_fmr(i)/(phisl1(i)*dlkgl*1.d3)
c        Need to skip for dry nodes
         elseif(sl1(i).eq.1.d0.or.a_fmr(i).gt.sl1(i).and.
     +     sl1(i).ne.0.d0)then
c..... Factor for saturated system, or unsaturated to consider only
c........ the wetted proportion
           actfrc = 1.d0/(phi(i)*dlkgl*1.d3)
         elseif(sl1(i).eq.0.d0)then
           actfrc=0.d0
         endif
c
          WRITE (79,'(3F11.3,3e12.5,60e12.4)') XXX(I),YYY(I),ZZZ(I),
     +         sl1(i),a_fmr(i),actfrc,(amin3(I,J),J=1,nmkin)
       END IF
c
C--------------------------------------------------------------
C
740      CONTINUE
c
c
       RETURN
       END
c
c
c
c-------------------------------------------------------------------------------
c
c
c
       SUBROUTINE WRITE_TIME
C
C*****************WRITE variables versus time at specified emements**********
C
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        include 'perm_v2.inc'
        COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
        COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
        COMMON/E6/T(MNEL)
        COMMON/SOLUTE6/SLOLD(MNEL)        ! old liquid saturation
        COMMON/SOLUTE8/SL1(MNEL)          ! new liquid saturation
        COMMON/SOLUTE9/SG1(MNEL)          ! new gas saturation
        COMMON/SOLUTE10/PHIOLD(MNEL)      ! porosity at previous time step
        COMMON/WRICON1/ELEMW(200)
        COMMON/E1/ELEM(MNEL)
        COMMON/E5/P(MNEL)
        COMMON/DRYOUT/IDRY(MNOD),ADRY(MNOD,MPRI)
        COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)                ! residual in precipitates
        common/water_activity/aw(mnod)    ! water activity
        common/ion_str2/str_node(mnel)    ! ionic strength for all nodes
        double precision PRECIP(mmin),DGP2(mgas)
        COMMON/TRANGAS9/NGAS1             ! Number of gaseous species
        CHARACTER*5 ELEM,ELEMW
        COMMON/PARNP/NPL,NPG           ! specify in EOS module
        COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
        COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
        double precision CDUM(MAQT)    ! working array for aqueous concentrations
        double precision AQDUM(MAQT)   ! working array for aqueous species concentrations
        double precision ADDUM(MADS)   ! working array for adsorbed species concentrations
        double precision XDUM(MXsites,MEXC)  ! working array for exchanged concentrations
        COMMON/MOP_REACT/MOPR(20)      ! controlling parameters for reactive transport
        COMMON/min_SI/SIM(MNEL,MMIN)   ! Mineral saturation index (log(Q/K)) for all nodes
        COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
        character*100 form1
        common/molweight/wm_aqt(maqt)  ! Molecular weight of all species, g/mol
        common/Print_Unit_Name/ Name_Conc, Name_Mine
        character*10 Name_Conc, name_unit
        character*45 Name_Mine
!
        common/clay_swell1/ iswell
        common/clay_swell2/ vmin_old(mnel,mmin)  ! previous mole volume for all node
        common/clay_swell3/ vmin0(mmin)          ! initial mole volume
!
        COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
c
c--------------------------------------------------------------------
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(/1X,15('*'),'WRITE_TIME 1.0, 2003.7.30: Write variables'
     x' vs. time for species, minerals and gases',15('*'))
c
C--------------------------------------------------------------Write title
c
        IF (TIMETOT .EQ. 0.0D0.or.icall.eq.1) THEN
           WRITE(62,600)
600        FORMAT(/10X,'----Time evolution at specified elements----'/)
           WRITE(62,*) '  Unit: '
           WRITE(62,"(4x,'- Aqueous species: Total concencentrations',
     &       ' in ', A10)")   Name_Conc
           name_unit=Name_conc
           if(iconflag.gt.1) name_unit ='mol/L'
           IF (NWAQ .GT. 0)
     &      WRITE(62,"(4x,'- Aqueous species: Individual',
     &        ' concentrations in ',A10)") name_unit
           IF (NWADS .GT. 0)
     &      WRITE(62,"(4x,'- Surface species: Concentrations in ',A10)")
     &         name_unit
           IF (nwexc.GT. 0)
     &      WRITE(62,"(4x,'- Exchanged species: Concentrations in ',
     &         A10)") name_unit
           WRITE(62,"(4x,'- Component in dry nodes in mol/L medium')")
c
           IF (NWMIN .GT. 0)   THEN
              WRITE(62,"('    - Minerals: ', A45)")  Name_mine

           END IF
c
           IF (NGAS .GT. 0)   THEN
              IF (IEOS.EQ.9)   THEN
                 WRITE(62,1004)
1004             FORMAT ('    - Gas: partial pressure (bar)')
                 GO TO 799
              END IF
              IF (MINFLAG .GE. 1)  THEN
                 WRITE(62,1013)
1013             FORMAT ('    - Gas: volume fraction')
                                   ELSE
                 WRITE(62,1014)
1014             FORMAT ('    - Gas: partial pressure (bar)')
              END IF
           END IF
799        CONTINUE
c
           WRITE(62,*)
c--------------------------------------
C
       if(mopr(7).eq.0.or.mopr(7).eq.4)then
              WRITE(62,721) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).eq.1)then
              WRITE(62,722) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).eq.2)then
              WRITE(62,723) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).eq.3)then
              WRITE(62,724) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).eq.5)then
              WRITE(62,725) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).eq.6)then
              WRITE(62,726) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).eq.7)then
              WRITE(62,727) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       elseif(mopr(7).ge.8)then
              WRITE(62,728) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
       endif
c
721      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A11))
722      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A8))
723      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A9))
724      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A10))
725      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A12))
726      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A13))
727      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A14))
728      FORMAT('ELEM    Time(day)      P(bar)        Sg         Sl',
     +      '      pH    I.STR         aH2O',
     +      '   Porosity  Perm(m^2)  ',200(1x,A15))
        END IF
c
C--------------------------------------------------------------Write results
c
        DO 300 INO=1,NWNOD
         NNG=IWNOD(INO)
         dlkgl=dwat(nng)              ! water density (kg/l)
c
         IF (ICONFLAG.ge.1) THEN   ! in mol/l
              cb = 0.d0
              cbtot = 0.d0
               DO J=1,NPRI
                CDUM(J)=ctot(nng,j)*dlkgl
c---------------- add for charge balance
                  cb=cdum(j)*z(j)+cb
                  cbtot = dabs(cdum(j)*z(j)) + cbtot
               END DO
c              need to include charged surface in chg balance
               do j=1,nads
                  cb=d(nng,j)*zd(j)+cb
                  cbtot = dabs(d(nng,j)*zd(j)) + cbtot
               enddo

c              -------------------
c              Use additional separate arrays for printing of secondary species,
c              sorbed species, and exchange species (all in molalities)
               do j=1,nwaq
                aqdum(j)=c(nng,iwaq(j))*dlkgl
               enddo
               do j=1,nwads
                addum(j)=d(nng,iwads(j))*dlkgl
               enddo
               do j=1,nwexc
                 do isite=1,NXsites
                   XDUM(isite,J)=XCADS(nng,isite,iwexc(j))*dlkgl
                 enddo
               enddo
C              -------------------
               if (iconflag .ge. 2)   then      !converts cdum only
                do j=1,npri
c                 if (napri(j) .ne. 'tracer') then
                   CDUM(J)=CDUM(J)*wm_aqt(j)    ! g/L
                   if (iconflag .eq. 3)   then
                     CDUM(J)=CDUM(J)*1.0d3      ! mg/L (ppm)
c                   end if
                 endif
                end do
               end if

         ELSE      ! in mol/kg
            cb = 0.d0
            cbtot = 0.d0
               DO J=1,NPRI
                  CDUM(J)=ctot(nng,j)
c                 add for charge balance
                  cb=cdum(j)*z(j)+cb
                  cbtot = dabs(cdum(j)*z(j)) + cbtot
               END DO
c              need to include charged surface in chg balance
               do j=1,nads
                  cb=d(nng,j)*zd(j)+cb
                  cbtot = dabs(d(nng,j)*zd(j)) + cbtot
               enddo
c
c              -------------------
c              NS3/08 use additional separate arrays for printing of secondary species,
c              sorbed species, and exchange species (all in molalities)
               do j=1,nwaq
                aqdum(j)=c(nng,iwaq(j))
               enddo
               do j=1,nwads
                addum(j)=d(nng,iwads(j))
               enddo
               do j=1,nwexc
                 do isite=1,NXsites
                   XDUM(isite,J)=XCADS(nng,isite,iwexc(j))
                 enddo
               enddo
c              -------------------
c
         END IF
!
!........Calculate for printout mineral options (minflag=3, change of volume fraction in %)
!
         if(minflag.eq.1 .or. minflag.eq.3) then           
           DO J=1,NWMIN
             IF (ISWELL .NE. 1)    THEN
                PRECIP(J)=PRE(NNG,IWMIN(J))-PINIT(NNG,IWMIN(J))  ! pre in moles per liter medium
                PRECIP(J)=PRECIP(J)*vmin(iwmin(j))               ! change in volume fraction
                                 ELSE
                PRECIP(J)=PRE(NNG,IWMIN(J))*vmin(iwmin(j))-
     +                  PINIT(NNG,IWMIN(J))*vmin0(iwmin(j))      ! for clay swelling
             END IF
!
             if (minflag .eq.3 ) PRECIP(J)=PRECIP(J)*100.0d0     ! Change in Vf %
!
           END DO
         else if(minflag.eq.2)then
           DO J=1,NWMIN
             PRECIP(J)=PRE(NNG,IWMIN(J))*vmin(iwmin(j))       ! total volume fraction
           END DO
         else
           DO J=1,NWMIN
             PRECIP(J)=PRE(NNG,IWMIN(J))-PINIT(NNG,IWMIN(J))  ! pre in moles per liter medium
             PRECIP(J)=PRECIP(J)*1000.0D0                     ! change in mol/m**3 medium
           END DO
         endif
!
!..................
!
         DO 330 J=1,NGAS
c.... for minflag >= 1, volume fraction, partial pressure (bar) otherwise
            if(minflag.ge.1)then
              dgp2(j)=pfug(nng,j)*1.d5/p(nng)
            else
              dgp2(j)=pfug(nng,j)
            endif
            if (ieos.eq.9) dgp2(j)=pfug(nng,j)
330      CONTINUE
C
        if (kcpl.ne.2)   then
           PHInng=phi(nng)
           PERMnng=perm(2,nng)

                          else     ! only monitor porosity and perm. changes
           PHInng=phim(nng)
           PERMnng=permm(2,nng)
         end if
C
        TIMEDAY=TIMETOT/86400.d0   ! Time in day
c
       if(mopr(7).eq.0.or.mopr(7).eq.4)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E12.4)"
       elseif(mopr(7).eq.1)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E9.1)"
       elseif(mopr(7).eq.2)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E10.2)"
       elseif(mopr(7).eq.3)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E11.3)"
       elseif(mopr(7).eq.5)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E13.5)"
       elseif(mopr(7).eq.6)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E14.6)"
       elseif(mopr(7).eq.7)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E15.7)"
       elseif(mopr(7).ge.8)then
        form1= "(A5,1x,E12.6,e12.4,2E11.4,F8.3,2E12.4,F9.5,E12.4,
     &200E16.8)"
       endif
c
       WRITE(62,form1) ELEM(NNG),TIMEDAY,P(NNG)/1.0D+5,
     +               sg1(NNG),sl1(NNG),PH(NNG),
     +                              str_node(nng),aw(nng),
     +                               PHInng,PERMnng,
     +                 (CDUM(IWCOM(J)),J=1,NWCOM),
     &                 (aqdum(J), j=1,nwaq),
     &                 (addum(J), j=1,nwads),
     &                 ((XDUM(isite,J), J=1,NWEXC),isite=1,NXsites),
     +                 (PRECIP(J),          J=1,NWMIN),
     +                 (dgp2(J),            J=1,NGAS), cb
c     +                 100.d0*(cb/cbtot)
300    CONTINUE
C---------------------------------
C
500      CONTINUE
C
       RETURN
       END
c
c-------------------------------------------------------------------------------
c
       SUBROUTINE WRITE_TIME_ECO2
C
C                                  For ECO2N flow module
C
C*****************WRITE varibles versus time at specified emements**********
C
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        INCLUDE 'common_v2.inc'
        include 'perm_v2.inc'
        COMMON/WRICON/ NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                 IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                 nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                 iwexc(mexc)
        COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
        COMMON/E6/T(MNEL)
        COMMON/SOLUTE6/SLOLD(MNEL)        ! old liquid saturation
        COMMON/SOLUTE8/SL1(MNEL)          ! new liquid saturation
        COMMON/SOLUTE9/SG1(MNEL)          ! new gas saturation
        COMMON/SOLUTE10/PHIOLD(MNEL)      ! porosity at previous time step
        COMMON/WRICON1/ELEMW(200)
        COMMON/E1/ELEM(MNEL)
        COMMON/E5/P(MNEL)
        double precision PRECIP(mmin),DGP2(mgas)
        COMMON/TRANGAS9/NGAS1          ! Number of gaseous species
        CHARACTER*5 ELEM,ELEMW
        COMMON/PARNP/NPL,NPG           ! specify in EOS module
        COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
        COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
!
        double precision CDUM(MAQT)    ! working array for aqueous concentrations
        double precision AQDUM(MAQT)   ! working array for aqueous species concentrations
        double precision ADDUM(MADS)   ! working array for adsorbed species concentrations
        double precision XDUM(MXsites,MEXC)  ! working array for exchanged concentrations
!
        common/clay_swell1/ iswell
        common/clay_swell2/ vmin_old(mnel,mmin)  ! previous mole volume for all node
        common/clay_swell3/ vmin0(mmin)          ! initial mole volume
!
        COMMON/SOLIDco2/SMco2(NMNOD)   ! CO2 TRAPPED in solid phase
        COMMON/MOP_REACT/MOPR(20)      ! controlling parameters for reactive transport
        COMMON/min_SI/SIM(MNEL,MMIN)   ! Mineral saturation index (log(Q/K)) for all nodes
        COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
        character*100 form1
        common/molweight/wm_aqt(maqt)  ! molecular weight of all species, g/mol
        common/Print_Unit_Name/ Name_Conc, Name_Mine
        character*10 Name_Conc,name_unit
        character*45 Name_Mine
!
        common/co2_gene/nco2
        COMMON/EOS_INDICATOR/ IEOS       ! Indicate EOS module used
        COMMON/E3/EVOL(MNEL)
!
c--------------------------------------------------------------------
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(11,899)
  899   FORMAT(6X,'WRITE_TIME_ECO2 1.0 30 July    2003',6X,
     X  'Write variables vs. time for species, minerals and gases')
c
C--------------------------------------------------------------Write title
c
        IF (TIMETOT .EQ. 0.0D0.or.icall.eq.1) THEN
           WRITE(62,600)
600        FORMAT(/10X,'----Time evolution at specified elements----'/)
           WRITE(62,*)
           WRITE(62,*)'  Unit: '
           WRITE(62,"(4x,'- Aqueous species: Total concencentrations',
     &       ' in ', A10)")   Name_Conc
           name_unit=Name_conc
           if(iconflag.gt.1) name_unit ='mol/L'
           IF (NWAQ .GT. 0)
     &      WRITE(62,"(4x,'- Aqueous species: Individual',
     &        ' concentrations in ',A10)") name_unit
           IF (NWADS .GT. 0)
     &      WRITE(62,"(4x,'- Surface species: Concentrations in ',A10)")
     &         name_unit
           IF (nwexc.GT. 0)
     &      WRITE(61,"(4x,'- Exchanged species: Concentrations in ',
     &         A10)") name_unit
           WRITE(62,"(4x,'- Component in dry nodes in mol/L medium')")
c
           IF (NWMIN .GT. 0)   THEN
              WRITE(62,26)
26            FORMAT ('    - SMco2: Total CO2 sequstered in',
     +                ' mineral phases in  kg/m**3 medium')
              WRITE(62,27)
27            FORMAT ('    - Permeability in  m**2')
C
              WRITE(62,"('    - Minerals: ', A45)")  Name_mine
           END IF
c
           IF (NGAS .GT. 0)  THEN
              WRITE(62,1014)
1014            FORMAT ('    - Gas: Partial pressure (bar)')
           END IF
c
           WRITE(62,*)
c--------------------------------------
C
           WRITE(62,720) ('t_'//NAPRI(IWCOM(J)),J=1,NWCOM),
     &                     (naaqt(iwaq(j)),j=1,nwaq),
     &                     (naads(iwads(j)),j=1,nwads),
     &                     (('X_'//NAEXC(J), J=1,NWEXC),jj=1,NXsites),                ! Multi-site
     &                     (NAMIN(IWMIN(J)),J=1,NWMIN),
     &                     (NAGAS(J),       J=1,NGAS),
     &                     'chg.bal.(equiv)'
720        FORMAT('ELEM    Time(day)        Sg         Sl',
     &         '      T(C)     pH       SMco2',
     &         '     Porosity    Perm(m^2) ',200(1x,A11))
!
!..........Write title for CO2 trapping file
!
           WRITE(80,"('    Time(day)   SMco2_tot(kg)   AQco2_tot(kg)')")
!
        END IF
c
c       Initialize all printing arrays to 0
        do j=1,npri
             cdum(j)=0.d0
        enddo
        do j=1,nwaq
             aqdum(j)=0.d0
        enddo
        do j=1,nwads
             addum(j)=0.d0
        enddo
        do j=1,nwexc
          do isite=1,NXsites
                XDUM(isite,J)=0.d0
          enddo
        enddo
c
C--------------------------------------------------------------Write results
!
!
c        IF(TIMETOT .EQ. 0.d0)  GO TO 500
        DO 300 INO=1,NWNOD
         NNG=IWNOD(INO)
         dlkgl=dwat(nng)              ! water density (kg/l)
c
         cb=0.d0
         cbtot=0.d0
         IF (ICONFLAG.ge.1) THEN   ! in mol/l
               DO J=1,NPRI
                if(sl1(nng).ne.0.d0) CDUM(J)=ctot(NNG,J)*dlkgl
                  cb=cdum(j)*z(j)+cb
                  cbtot = dabs(cdum(j)*z(j)) + cbtot
               END DO
c              need to include charged surface in chg balance
               do j=1,nads
                  cb=d(nng,j)*zd(j)+cb
                  cbtot = dabs(d(nng,j)*zd(j)) + cbtot
               enddo

c              -------------------
c              Use additional separate arrays for printing of secondary species,
c              sorbed species, and exchange species (all in molalities)
               do j=1,nwaq
                if(sl1(nng).ne.0.d0) aqdum(j)=c(nng,iwaq(j))*dlkgl
               enddo
               do j=1,nwads
                if(sl1(nng).ne.0.d0) addum(j)=d(nng,iwads(j))*dlkgl
               enddo
               do j=1,nwexc
                 do isite=1,NXsites
                   if(sl1(nng).ne.0.d0) 
     &                  XDUM(isite,J)=XCADS(nng,isite,iwexc(j))*dlkgl
                 enddo
               enddo
c              -------------------
               if (iconflag .ge. 2)   then           ! converts cdum only
                do j=1,npri
                 if (napri(j) .ne. 'tracer') then
                   CDUM(J)=CDUM(J)*wm_aqt(j)   ! g/L
                   if (iconflag .eq. 3)   then
                     CDUM(J)=CDUM(J)*1.0d3     ! mg/L (ppm)
                   end if
                 endif
                end do
               end if

         ELSE          ! in mol/kg
               DO J=1,NPRI
                if(sl1(nng).ne.0.d0)  CDUM(J)=ctot(nng,j)
                  cb=cdum(j)*z(j)+cb
                  cbtot = dabs(cdum(j)*z(j)) + cbtot
               END DO
c              need to include charged surface in chg balance
               do j=1,nads
                  cb=d(nng,j)*zd(j)+cb
                  cbtot = dabs(d(nng,j)*zd(j)) + cbtot
               enddo
c              -------------------
c              Use additional separate arrays for printing of secondary species,
c              sorbed species, and exchange species (all in molalities)
               do j=1,nwaq
                if(sl1(nng).ne.0.d0) aqdum(j)=c(nng,iwaq(j))
               enddo
               do j=1,nwads
                addum(j)=d(nng,iwads(j))
               enddo
               do j=1,nwexc
                 do isite=1,NXsites
                   if(sl1(nng).ne.0.d0)
     &                   XDUM(isite,J)=XCADS(nng,isite,iwexc(j))
                 enddo
               enddo
c              -------------------
         END IF
!
!........Calculate printout mineral options (minflag=3, change of volume fraction in %)
!
         if(minflag.eq.1 .or. minflag.eq.3) then           
           DO J=1,NWMIN
             IF (ISWELL .NE. 1)    THEN
                PRECIP(J)=PRE(NNG,IWMIN(J))-PINIT(NNG,IWMIN(J))  ! pre in moles per liter medium
                PRECIP(J)=PRECIP(J)*vmin(iwmin(j))               ! change in volume fraction
                                 ELSE
                PRECIP(J)=PRE(NNG,IWMIN(J))*vmin(iwmin(j))-
     +                  PINIT(NNG,IWMIN(J))*vmin0(iwmin(j))    ! for clay swelling
             END IF
!
             if (minflag .eq.3 ) PRECIP(J)=PRECIP(J)*100.0d0   ! Change in Vf %
!
           END DO
         else if(minflag.eq.2)then
           DO J=1,NWMIN
             PRECIP(J)=PRE(NNG,IWMIN(J))*vmin(iwmin(j))        ! total volume fraction
           END DO
         else
           DO J=1,NWMIN
             PRECIP(J)=PRE(NNG,IWMIN(J))-PINIT(NNG,IWMIN(J))   ! pre in moles per liter medium
             PRECIP(J)=PRECIP(J)*1000.0D0                      ! change in mol/m**3 medium
           END DO
         end if
!
!.............................
!
         DO 330 J=1,NGAS
c.... for minflag = 1, print vol fraction,fugacity otherwise
c             if(minflag.ge.1)then
c                dgp2(j)=pfug(nng,j)*1.d5/p(nng)
c             else
                dgp2(j)=pfug(nng,j)
c             endif
330     CONTINUE
C
        if (kcpl.ne.2)   then
           PHInng=phi(nng)
           PERMnng=perm(2,nng)

                          else     ! only monitor porosity and perm. changes
           PHInng=phim(nng)
           PERMnng=permm(2,nng)
         end if
C
         TIMEDAY=TIMETOT/8.6400d4    ! Time in year
c
       if(mopr(7).eq.0.or.mopr(7).eq.4)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E12.4)"
       elseif(mopr(7).eq.1)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E9.1)"
       elseif(mopr(7).eq.2)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E10.2)"
       elseif(mopr(7).eq.3)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E11.3)"
       elseif(mopr(7).eq.5)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E13.5)"
       elseif(mopr(7).eq.6)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E14.6)"
       elseif(mopr(7).eq.7)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E15.7)"
       elseif(mopr(7).ge.8)then
           form1 = "(A5,1x,E12.6,2e11.4,2F8.3,200E16.8)"
       endif
!
       WRITE(62,form1) ELEM(NNG),TIMEDAY,sg1(NNG),sl1(NNG),
     &                           T(NNG),PH(NNG),
     &                  SMco2(NNG),PHInng,PERMnng,
     &                  (CDUM(IWCOM(J)),J=1,NWCOM),
     &                  (aqdum(J), j=1,nwaq),
     &                  (addum(J), j=1,nwads),
     &                  ((XDUM(isite,J), J=1,NWEXC),isite=1,NXsites),
     &                  (PRECIP(J),          J=1,NWMIN),
     &                  (DGP2(J),            J=1,NGAS),
     &                  cb
300    CONTINUE
C---------------------------------
C
cc500    CONTINUE
!
!--------
!.....Calculate trapped in different phases and then write ..........
!--------
!
      if (mopr(8) .ge.1 )   then 
!
         SMco2_tot = 0.0d0
         AQco2_tot = 0.0d0
!
         do i=1,nnod
!
!...........In solid phase
!
           if (evol(i) .lt. 1.0d20)   then
            if (SMco2(i) .gt. 0.01d0)   then
               SMco2_tot = SMco2_tot + evol(i)*SMco2(i)   ! Kg
            end if
!
!...........In aqueous phase
!
            if (UT(i,nco2) .gt. 1.0d-6)   then
               AQco2_tot = AQco2_tot + evol(i)*phi(i)*SL1(i)
     &                                 *UT(i,nco2)*44.0d0         ! Kg
            end if
           end if
!
         end do
!
!.......Write
!
      end if
!
      WRITE(80,'(E12.4,2E16.4)') TIMEDAY, SMco2_tot, AQco2_tot
c
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
      SUBROUTINE COMPUTE_MASS(icall,deltat0)
C
C**** Compute mass input to and output from the system for water and solutes ****
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      INCLUDE 'common_v2.inc'
      COMMON/BALAN1/SOLUTINP(MPRI),SOLUTOUT(MPRI),
     1              SOLUTINI(MPRI),SOLIDINI(MPRI),
     2              SOLUTNOW(MPRI),SOLIDNOW(MPRI)
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E3/EVOL(MNEL)
C
      COMMON/SOLUTE1/VELDAR(MNPH*MNCON) ! darcy velocity
      COMMON/SOLUTE6/SLOLD(MNEL)        ! old liquid saturation
      COMMON/SOLUTE7/SGOLD(MNEL)        ! old gas saturation
      COMMON/SOLUTE8/SL1(MNEL)          ! new liquid saturation
      COMMON/SOLUTE9/SG1(MNEL)          ! new gas saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)      ! porosity at previous time step
      COMMON/PARNP/NPL,NPG              ! specify in EOS module
C
      COMMON/TRANGAS9/NGAS1          ! Number of gaseous species FOR TRANSPORT
C
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
C
      COMMON/G4/ELEG(MNOGN)
      COMMON/G7/G(MNOGN)
      COMMON/G8/EG(MNOGN)
      COMMON/G9/NEXG(MNOGN)
      COMMON/G12/LCOM(MNOGN)
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
c
c--------------------------------COMMON blocks for Kd adsorption and decay
c
      common/kddca3/kddp(mpri)    ! pointer to the primary species
      common/Kddca4/vkd(30,mpri)  ! values of Kd in initial zones
      common/Kddca5/izonekd(mnod) ! Kd zone code
      common/Kddca6/sden(30,mpri) ! solid density
      common/Kddca7/nkdd          ! number of species with Kd adsorption
C--------------------------------------------------------------------------
C
      COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
      character*5 elem,eleg
      integer*8 icall
      double precision deltat0
c
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(1x,15('*'),' COMPUTE_MASS 1.0, 2003.7.30: Compute mass'
     X' input to and output from the system',15('*'))
c
C---------------------------------------Set variables to zero at the beginning
      if (icall.EQ.1) then
         DO 20 J=1,NPRI
            SOLUTINP(J)= 0.0D0
            SOLUTOUT(J)= 0.0D0
C
            SOLUTINI(J)= 0.0D0
            SOLIDINI(J)= 0.0D0
20       CONTINUE
C-----------------------Initial amount of component presented in aqueous phase
        DO 100 I=1,NNOD
           IF (EVOL(I).GE.1.0D+20)  GO TO 100
           EVOLI=EVOL(I)*1000.D0
           phislov = EVOLI*PHIOLD(I)*SLOLD(I)
           DO 200 J=1,NPRI
              SOLUTINI(J)=SOLUTINI(J) + phislov
     +                  *UTOLD(I,J)
200        CONTINUE
C-------------------------Initial amount of component presented in solid phase
           do m=1,nmin
            ncp=ncpm(m)
            do k=1,ncp
              j=icpm(m,k)
                 solidini(j)=solidini(j)+stqm(m,k)*pre0(i,m)
     +                       *evoli
            enddo
           enddo
c
           do m=1,ngas
            ncp=ncpg(m)
            do k=1,ncp
              j=icpg(m,k)
                 solidini(j)=solidini(j)+stqg(m,k)*gp(i,m)
     +                       *evoli
            enddo
           enddo
C
C*************************************************** Addition for Kd adsorption
C
      IF (NKDD .GT. 0)    THEN
         DO J=1,NPRI
c
c--solid density (kg/dm**3), and Kd(l/kg=mass/kg solid / mass/l water)
C
            KDDS=KDDP(J)       ! Number in the species list for Kd and decay
            IF (KDDS .EQ. 0)  GO TO 399
            KDDZONE=IZONEKD(I)    ! Kd zone code
            IF (KDDZONE .LE. 0)  THEN
               SDEN2=0.0D0
               VKD2=1.0D0
               GO TO 300
            END IF
            SDEN2=SDEN(KDDZONE,KDDS)  ! solid density
            VKD2=VKD(KDDZONE,KDDS)    ! Kd value; or r factor if solid density=0
300         CONTINUE
C
            IF (SDEN2.EQ.0.0D0 .AND. VKD2.GE.1.0D0)  THEN
c----------------If density is zero vkd2 is retardation factor
               RETARD1=VKD2-1.0D0       ! R1=R-1
            END IF
            IF (SDEN2.GT.0.0D0 .AND. VKD2.GE.0.0D0)  THEN
               RETARD1=(1.0D0-PHIOLD(I))*SDEN2*VKD2/(PHIOLD(I)*SLOLD(I))
            END IF
C
            SOLIDINI(J)=SOLIDINI(J) + phislov*UTOLD(I,J)*RETARD1
C
399         CONTINUE
         END DO
      END IF
C
C**********************************************************************
c
           do m=1,nads
            ncp=ncpad(m)
            do k=1,ncp
              j=icpad(m,k)
                 solidini(j)=solidini(j)+stqd(m,k)*d(i,m)*
     +                       phislov
            enddo
           enddo
C
C************************************************ Addition for exchange
C
      do isite=1,NXsites
      do k=1,nexc
         do j=1,npri
            solidini(j)=solidini(j)+stqx(k,j)*xcads(i,isite,k)*phislov
         end do
      end do
      end do
C
C**********************************************************************
c
100     continue
         go to 299
      END IF
C
C-------------------------------Calculate mass entering and leaving the system
      DO 140 IOGN=1,NOGN
         J=NEXG(IOGN)
c
         vliqw = dwat(j)
c........Need to keep vliqw greater than zero
         if(dwat(j).eq.0.d0)vliqw=1.d0
         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).GE.0.D0) THEN
            IZONEJ=IZONEBW(J)
            DO 160 IPRI=1,NPRI
               SOLUTINP(IPRI)= SOLUTINP(IPRI) + G(IOGN)
     +                      *UB(IZONEJ,IPRI)*deltat0/vliqw
160         CONTINUE
         END IF
cns13/6/08         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).GE.0.D0) THEN
         IF(LCOM(IOGN).EQ.1 .AND. G(IOGN).GE.0.D0
     &         .and.IZONEBW(J).ne.0) THEN
            DO 180 IPRI=1,NPRI
               SOLUTOUT(IPRI)= SOLUTOUT(IPRI) - G(IOGN)*UT(J,IPRI)
     +                      *deltat0
180         CONTINUE
         END IF
140   CONTINUE
C
  299   CONTINUE
C
        RETURN
        END
c
c
c-------------------------------------------------------------------------------
c
c
       SUBROUTINE WRITE_MASS
C
C**** Write mass input to and output from the system for water and solutes ****
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      INCLUDE 'common_v2.inc'
      COMMON/WRICON/NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1              IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2              nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3              iwexc(mexc)
      COMMON/BALAN1/SOLUTINP(MPRI),SOLUTOUT(MPRI),
     1              SOLUTINI(MPRI),SOLIDINI(MPRI),
     2              SOLUTNOW(MPRI),SOLIDNOW(MPRI)
      double precision BALANCE1(MPRI),BALANCE2(MPRI),
     1         BDIFF(MPRI),RELD(MPRI),relinp(mpri)
C
      COMMON/E1/ELEM(MNEL)
      COMMON/E3/EVOL(MNEL)
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
      COMMON/SOLUTE8/SL1(MNEL)          ! new liquid saturation
      COMMON/SOLUTE10/PHIOLD(MNEL)      ! porosity at previous time step
      COMMON/TRANGAS9/NGAS1             ! Number of gaseous species FOR TRANSPORT
      COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)    ! residual in precipitates
      common/dry_salt/nsalt,isalt(mmin)
      common/aqkin16/NoTrans(mpri)      ! >0: not subject to transport
      double precision tdryc(mpri),adryrmax(mpri),tdryfrc(mpri)
      double precision voldry
      CHARACTER*5 ELEM
c
c--------------------------------COMMON blocks for Kd adsorption and decay
      common/kddca3/kddp(mpri)    ! pointer to the primary species
      common/Kddca4/vkd(30,mpri)  ! values of Kd in initial zones
      common/Kddca5/izonekd(mnod) ! Kd zone code
      common/Kddca6/sden(30,mpri) ! solid density
      common/Kddca7/nkdd          ! number of species with Kd adsorption
!
!............................................................................
!
      common/co2_gene/nco2
      COMMON/EOS_INDICATOR/ IEOS       ! Indicate EOS module used
      COMMON/SOLIDco2/SMco2(NMNOD)     ! CO2 TRAPPED in solid phase
!
!............................................................................
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(1X,15('*'),'WRITE_MASS 1.0, 2003.7.30: Write mass input to'
     X' and output from the system in geochem_v2.f',15('*'))
c
      IF(ICALL.EQ.1) OPEN(UNIT=37,FILE='mbalance.out',STATUS='UNKNOWN')
c
C-----------------------------------Write mass balance information
C
c... Save balances from previous call to routine and reset to zero
      DO J=1,NPRI
         SOLUTNOW(J)= 0.0D0
         SOLIDNOW(J)= 0.0D0
         tdryc(j) = 0.0d0
      end do
C---------------Current amount of component presented in aqueous phase
        DO 400 I=1,NNOD
          IF (EVOL(I).lt.1.0D+20)then
            EVOLI=EVOL(I)*1000.D0
            phislv = evoli*phiold(i)*sl1(i)
            DO 500 J=1,NPRI
              SOLUTNOW(J)=SOLUTNOW(J) +  phislv*UT(I,J)
500        CONTINUE
c
C-----------------Current amount of component presented in solid phase
          do m=1,nmin
            ncp=ncpm(m)
            do k=1,ncp
              j=icpm(m,k)
              solidnow(j)=solidnow(j)+stqm(m,k)
     +                            *pre(i,m)*evoli
            end do
          end do
c
          do m=1,ngas
            ncp=ncpg(m)
            do k=1,ncp
              j=icpg(m,k)
                 solidnow(j)=solidnow(j)+stqg(m,k)*gp(i,m)
     +                       *evoli
            enddo
          enddo
C
C******************* Addition for Kd adsorption (Tianfu Xu, 11/28/2001)
C
      IF (NKDD .GT. 0)    THEN
         DO J=1,NPRI
c
c--solid density (kg/dm**3), and Kd(l/kg=mass/kg solid / mass/l water)
C
            KDDS=KDDP(J)       ! Number in the species list for Kd and decay
            IF (KDDS .EQ. 0)  GO TO 399
            KDDZONE=IZONEKD(I)    ! Kd zone code
            IF (KDDZONE .LE. 0)  THEN
               SDEN2=0.0D0
               VKD2=1.0D0
               GO TO 300
            END IF
            SDEN2=SDEN(KDDZONE,KDDS)  ! solid density
            VKD2=VKD(KDDZONE,KDDS)    ! Kd value; or r factor if solid density=0
300         CONTINUE
C
            IF (SDEN2.EQ.0.0D0 .AND. VKD2.GE.1.0D0)  THEN
c----------------If density is zero vkd2 is retardation factor
               RETARD1=VKD2-1.0D0       ! R1=R-1
            END IF
            IF (SDEN2.GT.0.0D0 .AND. VKD2.GE.0.0D0)  THEN
               RETARD1=(1.0D0-PHI(I))*SDEN2*VKD2/(PHI(I)*SL1(I))
            END IF
C
            SOLIDNOW(J)=SOLIDNOW(J) + phislv*UT(I,J)*RETARD1
C
399         CONTINUE
         END DO
      END IF
C
C**********************************************************************
c
          do m=1,nads
            ncp=ncpad(m)
            do k=1,ncp
              j=icpad(m,k)
                 solidnow(j)=solidnow(j)+stqd(m,k)*d(i,m)*
     +                       phislv
            enddo
          enddo
C
C************************ Addition for exchange (Tianfu Xu, 11/27/2001)
C
      do isite=1,NXsites
      do k=1,nexc
         do j=1,npri
           solidnow(j)=solidnow(j)+stqx(k,j)*xcads(i,isite,k)*phislv
         end do
      end do
      end do
C
C**********************************************************************
c
c
        endif
c
400     continue
C
C----------Total of each component in in residual solids at dry blocks
C
        do j = 1,npri
           adryrmax(j) = 0.d0
           tdryc(j) = 0.0d0
           voldry = 0.d0
           tdryfrc(j) = 0.d0
           do i = 1,nnod
c............Aadd maximum adryr (mol/m^3)
             if(evol(i).lt.1.d20)then
               if(adryr(i,j).gt.0.d0)voldry= evol(i) + voldry
                 adryrmax(j) = max(1000.D0*adryr(i,j),adryrmax(j))
                 tdryc(j) = tdryc(j) + EVOL(I)*1000.D0*adryr(i,j)
             endif
           end do
            if(voldry.gt.0.d0)tdryfrc(j) = tdryc(j)/voldry
        end do
C
        TIMEDAY=TIMETOT/86400.0D0
        WRITE(37,680) TIMEDAY
680     FORMAT(//1X,'Mass (mol) balance for the whole system at time',
     +          E12.5,' days:'/)
        write(37,*)'                       Total Moles    '
        WRITE(37,685)
685     FORMAT(1X,'--------------------------------------------',
     1    '-------------------------------------------------')
        WRITE(37,700)
700     FORMAT(1X,' Component     Input    Output',
     2    '        Initial           Current             Residual'
     3/35X,'Aqueous Solid+gas Aqueous  Solid+gas   Tot(mol) Max(mol/m3)'
     4        /1X,'--------------------------------------------',
     5    '-------------------------------------------------')
        DO 740 J=1,NPRI
        IF (J.EQ.NW.OR.J.EQ.NH.OR.J.EQ.NE.OR.NoTrans(j).NE.0) GO TO 740
           WRITE(37,760)  NAPRI(J),SOLUTINP(J),SOLUTOUT(J),
     1                             SOLUTINI(J),SOLIDINI(J),
     2                             SOLUTNOW(J),SOLIDNOW(J),
     3                             tdryc(j),adryrmax(j)
760        FORMAT(2X,A10,8E10.3)
740     CONTINUE
c
C-------------------------------------------Calculate mass balance
C
        DO 785 J=1,NPRI
           BALANCE1(J)=SOLUTINP(J) +
     1                 SOLUTINI(J) + SOLIDINI(J)
c
           BALANCE2(J)= SOLUTOUT(J) +
     1                 SOLUTNOW(J) + SOLIDNOW(J)
           BDIFF(J) = BALANCE2(J) - BALANCE1(J)
           RELD(J)=(BDIFF(J)*100.d0)/BALANCE1(J)
           relinp(j) = (SOLIDNOW(J)-SOLIDINI(J))/
     +        (SOLUTINP(J) - (SOLUTNOW(J)- SOLUTINI(J)))
           relinp(j) = (relinp(j)-1.d0)*1.d2
785     CONTINUE
C
        write(37,*)''
        write(37,*)'            Total Moles and Mass Balances     '
c
        WRITE(37,800)
800     FORMAT(14X,' Input+Initial  Output+Current',
     1      ' Difference   RelDiff%   DSol/liq%  Tdry/Tvol(mol/m3)')
        WRITE(37,685)
        DO 840 J=1,NPRI
        IF (J.EQ.NW.OR.J.EQ.NH.OR.J.EQ.NE.OR.NoTrans(J).NE.0) GO TO 840
           WRITE(37,860)  NAPRI(J),BALANCE1(J),BALANCE2(J),
     1                BDIFF(J),RELD(J),relinp(j),tdryfrc(j)
860        FORMAT(1X,A10,2X,E12.4,4X,E12.4,4E13.4)
840     CONTINUE
        WRITE(37,685)
!
!
        RETURN
        END
c
c
c-------------------------------------------------------------------------------
c
c
       SUBROUTINE WRITE_ITER (istep,TIME)
C
C*********************** WRITE iteration messages *****************************
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      INCLUDE 'common_v2.inc'
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/STEADY/IFLOWSS,JSTEADY
      COMMON/TIMESTEA/TIMESTEA
      COMMON/E1/ELEM(MNEL)
      common/dtlim/max_chem_it,delt_conne,id_chem
      character*5 elem
      character*16 delt_conne
      character*5 id_chem
      double precision time,timeyr,timedum
      integer*8 istep
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(1x,15('*'),' *WRITE_ITER 1.0, 2003.7.30: Write iteration'
     x' information in geochem_v2.f',15('*'))
c
c---------------------------------------------------------------WRITE titles
c
      IF (istep.EQ.0.or.icall.eq.1) THEN
        max_chem_it=maxitch   !ns 9/06 first time only
        WRITE (39,100)
100     FORMAT(/10X,'---- Modelling Progress and Iteration',
     1                 ' Messages ----',
     1   //'   Description of column-headers and variables:',
     2     /'   ITERFL  = Iterations to solve flow',
     4     /'   ITERSQ  = Sequential transport-chemistry iterations',
     5     /'   MAXITCH = Maximum iterations to solve chemistry',
     7     /'   AVGITCH = Average iterations to solve chemistry',
     8     /'   DELTA-T = Current time step (days)',
     9     /'   Dt limit = Criterion limiting timestep',
     +     /'   Maxitch ID = Grid block: maximum chemical iterations')
        WRITE (39,120)
120     FORMAT(/5x,'STEP   Time(day)   ITERFL ITERSQ MAXITCH',
     1             ' AVGITCH Delta-t(day)  Dt limit',
     1  '       Maxitch ID',
     2  /1X,'--------------------------------------------------',
     3  '--------------------------------------------------------')

      write(32,"(/' see output file: ',a20,
     &    ' for time and time step information'/)") outiter

      END IF
C---------------------------------------------------WRITE iteration messages
      TIMEDUM=TIME
      IF (JSTEADY.EQ.1) TIMEDUM=TIME-TIMESTEA       ! for qusi steady-state
      timeday = timedum/8.6400d4
      dday=deltex/86400.d0
      WRITE (39,140)istep,timeday,ITERFL,ITERTR,max_chem_it,AVERITCH,
     1           dday, delt_conne, id_chem
140   FORMAT(1x,I8,1x,E13.6,I6,I6,I7,F9.3,2x,e11.5,2x,a16,
     1    2x,a5)
c
      max_chem_it = 0
c
      RETURN
      END
c
c
c
c-------------------------------------------------------------------------------
c
c
c
      subroutine chdump(time,node,itchem)
c
C******** This subroutine outputs all chemical data for the current node *********
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      COMMON/E1/ELEM(MNEL)
      COMMON/E5/P(MNEL)
      common/solute8/sl1(MNEL)           ! new liquid saturation
      common/solute9/sg1(MNEL)           ! new gas saturation
      common/dm/delten,deltex,for,ford
      common/nwater/nibw,iwtype,niwtype  ! for chdump only
      common/aqkin2/ntrx     ! total number of redox pair
      common/aqkin11/crx(mpri)
c
      character*5 elem
      double precision solid(mpri),gas(mpri),aq2(mpri),time
      double precision sorbed(mpri), aqkin_mol(mpri)
      integer*8 node,itchem
c
      data iout/20/
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(1x,15('*'),' CHDUMP 1.01: 2006.3.23: Write chemical'
     x' speciation data in geochem_v2.f',15('*'))
c
      if(icall.eq.1) then
        open(unit=iout,file='chdump.out',status='unknown')
      endif
c
c return if we exceeded 1000 dumps (to limit file size)
!      if(icall.gt.1000) then
!       write(32,"(' Calls to chdump reached limit of 1000',
!     &   ' at node: ',a5,' - printout is skipped')")
!     &   elem(node)
!       return
!      endif
c
      do i=1,npri
         solid(i)=0.d0
         aq2(i)=0.d0
         gas(i)=0.d0
         sorbed(i)=0.d0
      enddo
c
      if(node.eq.0) then
        if(nibw.le.niwtype) then
          write(iout,"(2x,'INITIAL WATER NUMBER:',I5)") iwtype
        else
          write(iout,"(2x,'BOUNDARY WATER NUMBER:',I5)") iwtype
        endif
        write(iout,"(2x,'TIME t = 0')")
        write(iout,"(2x,'SPECIATION TEMPERATURE (C): ',f10.3,
     &     '    PRESSURE (bar): ',f10.3)") tc2, Pt
      else
        write(iout,"(2x,'ELEMENT:         ',a5)") elem(node)
        timed=time/86400.d0
        write(iout,"(2x,'TIME (s/d/yr)    ',2(e12.5,'/'),e12.5)") time,
     &    timed, timed/365.24d0
        write(iout,"(2x,'TIME STEP (s) :  ',e10.5)") deltex
        write(iout,"(2x,'TEMPERATURE (C): ',f10.3)") tc2
        write(iout,"(2x,'PRESSURE (Pa):   ',e10.5)") p(node)
        write(iout,"(2x,'LIQ SATURATION:  ',f10.6)") sl1(node)
        write(iout,"(2x,'GAS SATURATION:  ',f10.6)") sg1(node)
        sumsat=sg1(node)+sl1(node)
        write(iout,"(2x,'SUM SATURATIONS: ',f10.6)") sumsat
        write(iout,"(2x,'POROSITY:        ',e10.5)") phi(node)
      endif
      write(iout,"(2x,'CHEM.ITERATIONS: ',i5)") itchem
c
      sum=0.d0
      sstr=0.d0        ! add stoichiometric i  ns 1/31/02
cns3/4/2010      do i=1,npri
      do i=1,npaq
        sum=sum+zsqi(i)*cp(i)
        sstr=sstr+zsqi(i)*u2(i)/u2(nw)*rmh2o
      enddo
      sstr=0.5d0*sstr
c
c--first we recalculate the amount of components in derived species,
c  minerals and gases for the balance check below:
c
c     derived aq.species contribution to total concentrations
      do j=1,naqx
        ncp=ncps(j)
        do n=1,ncp
          i=icps(j,n)
          aq2(i)=aq2(i)+stqs(j,n)*cs(j)
        enddo
        sum=sum+zsqi(npri+j)*cs(j)
      enddo
c
      str=0.5d0*sum
      write(iout,"('  IONIC STRENGTH:  ',e10.4,10x,
     & 'STOIC. IONIC STRENGTH:  ',e10.4)") str, sstr
c
c     sorbed contributions to total moles
      do j=1,nads
        ncp=ncpad(j)
        do n=1,ncp
          i=icpad(j,n)
          sorbed(i)=sorbed(i)+stqd(j,n)*cd(j)*cp(nw)
       enddo
      enddo

      if(node.ne.0) then  !skip initialization
c      If we equilibrated the surface species without changing the
c      input water composition (except water), need to add the sorbed
c      amounts to the total
       if(nads.gt.0.and.isurfeq(node).ne.0.and.time.eq.0.d0) then
        do i=1,npaq
          if(i.ne.nw) tt(i) = tt(i) + sorbed(i)
        enddo
       endif
      endif
c
c     equilibrated mineral contribution to total moles
      do j=1,nmequ
        ncp=ncpm(j)
        do n=1,ncp
          i=icpm(j,n)
          solid(i)=solid(i)+stqm(j,n)*cm(j)
       enddo
      enddo
c     kinetic minerals contribution to total moles
        dxh2o = deltex*xh2o
      do j=1,nmkin
        nmqj = nmequ+j
        ncp=ncpm(nmqj)
         rkdxh2 = rkin2(j)*dxh2o
        do n=1,ncp
           i=icpm(nmqj,n)
          solid(i)=solid(i) - stqm(nmqj,n)*rkdxh2
       enddo
      enddo
c     gas contribution to total moles
      do j=1,ngas
        ncp=ncpg(j)
        do n=1,ncp
          i=icpg(j,n)
          gas(i)=gas(i)+stqg(j,n)*cmg(j)
       enddo
      enddo
c
cns6/10  aqueous kinetic reactions contributions
      if (ntrx.gt.0) then
        do i=1,npri
           aqkin_mol(i)=crx(i)*deltex
        end do
      end if
c
c--primary species concentrations
      write(iout,"(/145x,'input')")
      write(iout,"( 5x,' species  ',11x,' molality ',2x,
     & '      gamma    ',
     & 'log activ.  tot.aqueous  tot.aqueous  tot.solid    tot.gas  ',
     & '    tot.sorbed   ','aq.kinetics ',
     & ' tot.system  tot.balance  aq.balance',
     & /25x,'(water in kg)',1x,'(H2O activity)',14x,
     & 'molality',6x,'moles',7x,'moles',7x,'moles',10x,'moles',8x,
     & 'moles',8x,'moles'/5x,174('-'))")
c
      do i=1,npri
        if(i.eq.nw) then
          bal2=(rmh2o+aq2(i))*cp(nw)-u2(i)  !water aqueous balance (cp(nw) is kg of water)
        else
          bal2=(cp(i)+aq2(i))*cp(nw)-u2(i)  !aqueous balance
        endif
        if(i.gt.npaq) then
            sorbed(i)=sorbed(i)+u2(i)
            u2(i)=0.d0
        endif
cns6/10        bal1=u2(i)+solid(i)+gas(i)+sorbed(i)-tt(i)    !total mass balance
        bal1=u2(i)+solid(i)+gas(i)+sorbed(i)-tt(i) - aqkin_mol(i)    !total mass balance
        actlog=dlog10(gamp(i)*cp(i))
c       water activity is stored in gamp(nw)  ! NS3/06
        if(i.eq.nw) actlog=dlog10(gamp(i))    ! NS3/06
        if(i.gt.npaq.and.node.eq.0) then
          write(iout,"(1x,i3,1x,a20,
     &   ' (Not part of initial speciation calculations)')") i,napri(i)
        else
          write(iout,"(1x,i3,1x,a20,e12.5,e15.5,f10.3,2x,9(e12.5,1x))")
     &    i,napri(i),cp(i),gamp(i),actlog,u2(i)/cp(nw),u2(i),
     &    solid(i),gas(i),sorbed(i),aqkin_mol(i),tt(i), bal1, bal2
        endif
      enddo
c--derived aq. species
      do i=1,naqx
        actlog=dlog10(gams(i)*cs(i))
        write(iout,"(1x,i3,1x,a20,e12.5,e15.5,f10.3,2x)")
     &       i+npri,naaqx(i),cs(i),gams(i),actlog
      enddo
c
c  surface complexes info
c  output boltzman terms in place of activity coefficients
      if (nads.gt.0)  then

        do k=1,nads
         if(node.ne.0) then
          gammads=dexp(phip2(iad_surf(k))*(-zd(k)))  !phi2p = exp(-psi F /RT)
          actlog=dlog10(cd(k)*gammads)
          write(iout,"(1x,i3,1x,a20,e12.5,e15.5,f10.3,2x)")
     &        k,naads(k),cd(k),gammads,actlog
          iprint_note=0
         else
          write(iout,"(1x,i3,1x,a20,' See note below - '
     &   'linked to surface of: ',10a20)")
     &        k,naads(k), naads_min(iad_surf(k))
          iprint_note=1
         endif
        enddo
        if(iprint_note.eq.0) then
          write(iout,"(/1x,'For surface species, activity coefficients',
     &     ' are effective values equal to 1/(exp{-F*psi/RT})**z ',
     &     ' for use in checking mass action')")
        elseif(iprint_note.eq.1) then
          write(iout,"(/1x,'Surface complexes are included but not yet',
     &    ' equilibrated - use the ichdump option to see equilibration',
     &    ' results on a gridblock-by-gridblock basis')")
        endif
c
        if(npot.ne.0.and.node.ne.0) then
         do n=1,nsurf
            write(iout,"(1x,'Potential term (-F*psi/RT) for surface',
     &         i2,' (',A12,'): ',e10.4)") n, naads_min(n),phip2(n)
         enddo
        endif
      endif
c
c.........Add GX changes for printing charge balance
        zzt_abs=0.0d0        ! for charge balance 08/09/05
        zzt_bls=0.0d0        
        do i=1,npri
          zzt_abs=zzt_abs+cp(i)*zabsi(i)
          zzt_bls=zzt_bls+cp(i)*z(i)
         end do
        do j=1,naqx
          zzt_abs=zzt_abs+cs(j)*zabsi(npri+j)
          zzt_bls=zzt_bls+cs(j)*z(npri+j)
         end do
        do j=1,nads
          zzt_abs=zzt_abs+cd(j)*zd(j)
          zzt_bls=zzt_bls+cd(j)*zd(j)
         end do

         pct_bal=zzt_bls*100.0d0/zzt_abs    
      write(iout,"(/5x,'Solution charge imbalance = ',
     +  e10.4,' equiv',2x,'(',e10.4,'%)')") zzt_bls, pct_bal
c
c--minerals
      if(nmequ.gt.0)
     & write(iout,"(/5x,'eq.minerals             moles  ',6x,
     & ' log(Q/K)   v.frac solid',/5x,60('-'))")
      do i=1,nmequ
        if(node.ne.0) then
         if(phi(node).ne.1.d0) then
           vfra=pre(node,i)/(1.d0-phi(node))*vmin(i)
           write(iout,"(1x,i3,1x,a20,e12.5,e15.5,2x,e10.4)")
     &       i,namin(i),cm(i),si2(i),vfra
         endif
        endif
        if(node.eq.0) then
        write(iout,"(1x,i3,1x,a20,e12.5,f15.5)")
     &       i,namin(i),cm(i),si2(i)
        endif
      enddo
      if(nmkin.gt.0)
     & write(iout,"(/5x,'kin.minerals            moles  ',6x,
     & ' log(Q/K)   v.frac solid rate(mol/sec)',/5x,75('-'))")
      do i=1,nmkin
        cmk=-rkin2(i)*deltex*xh2o
        sik=dlog10(si2k(i))
c-- print volume fractions for debug
       if(node.ne.0) then
         if(phi(node).ne.1.d0) then
         vfra=pre(node,nmequ+i)/(1.d0-phi(node))*vmin(nmequ+i)
         write(iout,"(1x,i3,1x,a20,e12.5,f15.5,2x,e10.4,2xe10.4)")
     &       i,namin(nmequ+i),cmk,sik,vfra,rkin2(i)
          endif
       endif
       if(node.eq.0) then
         write(iout,"(1x,i3,1x,a20,e12.5,f15.5)")
     &       i,namin(nmequ+i),cmk,sik
        endif
      enddo
c--gases
c     add printing of fugacity coeffients
      if(ngas.gt.0)
     & write(iout,"(/5x,' gases    ',2x,'            moles  ',3x,
     & ' log(Q/K/fuga)    log(fuga)   phi    '/5x,72('-'))")
      do i=1,ngas
        fug=dlog10(cg(i))
        write(iout,"(1x,i3,1x,a20,e12.5,2f15.5,f10.5)")
     &       i,nagas(i),cmg(i),sig2(i),fug, gamg(i)
      enddo
c
      write(iout,"(//)")
c
      return
      end
c
c-------------------------------------------------------------------------------
c
        subroutine waterchem
c
C*** Write water chemistry in format same as chemical.inp to chdump.out file ******
c              for selected grid blocks
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
        COMMON/E1/ELEM(MNEL)
        COMMON/E5/P(MNEL)
        COMMON/E6/T(MNEL)
        COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
        COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
        COMMON/WRICON/NWXY,NWDIM,NWTI,NWNOD,NWCOM,NWMIN,IWNOD(200),
     1                IWCOM(maqx),IWMIN(mmin),NWTS,NWTT,IWCOMT,
     2                nwaq,nwads,nwexc,iwaq(maqx),iwads(mads),
     3                iwexc(mexc)
        COMMON/SECPAR/PAR((MNEQ+1)*(MNPH*(MNB+MNK)+2)*MNEL)
        COMMON/PARNP/NPL,NPG          ! specify in EOS module
        COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
        character*5 elem
        data iout1/20/
c
        SAVE ICALL
        DATA ICALL/0/
        ICALL=ICALL+1
        IF(ICALL.EQ.1) WRITE(34,899)
  899   FORMAT(1x,15('*'),'Waterchem 1.0, 2004.3.30: Write chemical'
     x  ' speciation like the format in chemical.inp in geochem_v2.f')
c
        do iw=1,4
           write (iout1,*)
        end do
        write (iout1,909)
909     format (6X,'--- Write chemical speciation like the format',
     +  ' in chemical.inp -- ')
c
        timeyr = timetot/3.15569d7
c
        WRITE (iout1,910) TIMEYR
910     FORMAT(/18X,'  T= "',E12.6,' yr'//)
c
       do ino=1,nwnod
         i=iwnod(ino)
         write (iout1,912) elem(i)
912      format (/'  Grid block:  ',a5)
          write (iout1,'(i4,f10.3,2x,f10.4)') ino,t(i),p(i)/1.d5
         write (iout1,914)
914   format ('# component  flag    guess      ctotal    fix log(Q/K)')
         do j=1,npri
          cpp=c(i,j)
          utt=ctot(i,j)
          if(j.eq.nw) then
             cpp=1.0d0
             utt=1.0d0
          end if
          write (iout1,"(
     &     '''', a11,'''',1x,'1',1x,2E12.4,'  ''*''   0.0 ')")
     &      napri(j),cpp,utt
         end do
          write(iout1,916)
916       format("'*'")
!
      end do
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE READ_RESTART
C
C*************** Read restart data for reactive chemical transport **********************
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      COMMON/TRANGAS9/NGAS1                ! Number of gaseous species transport
      COMMON/TRANGAS1/PFUGOLD(NMNOD,NMGAS) ! Old partial pressure
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT
      COMMON/POV6/TSTART
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/SOLUTE6/SLOLD(MNEL)
      COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)    ! residual in precipitates
      common/water_activity/aw(mnod)    !water activity
c
c----------------------------------------Indicators from EOS module
C
      COMMON/SOLIDco2/SMco2(NMNOD)     ! CO2 TRAPPED in solid phase
      COMMON/EOS_INDICATOR/ IEOS       ! Indicate EOS module used
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' READ_RESTART 1.0, 2003.7.30: Read restart data for'
     X' reactive chemical transport')
C
      write(*,*) '   --> reading restart data for reactive transport'
C
      READ(21,*)     ! Geochemical state variables at time (s)'
      READ(21,'(E14.8)') TIMETOT
      if(timetot.ne.timin) then
        write(*,"(' Restart time in INCHEM file is not the',
     + ' same as'/' starting time for flow calculations')")
        write(*,*) 'timetot= ',timetot,' tstart=',timin
      endif
C
      READ(21,*)           ! Concentration of primary species
      DO 20 I=1,NNOD
         READ(21,25)  (C(I,N),N=1,NPRI)
25       FORMAT(6e15.8)
20    CONTINUE
C
c... pH
      read(21,*)
      DO I=1,NNOD
         read(21,444) ph(i)
      enddo
444      FORMAT(f10.6)
c
      NAQT=NPRI + NAQX
      IF (NAQX .EQ. 0)  GO TO 59
      READ(21,*)           ! Concentration of secondary species
      DO 27 I=1,NNOD
         READ(21,28)  (C(I,N),N=NPRI+1,NAQT)
28       FORMAT(6e15.8)
27    CONTINUE
59    CONTINUE
c
c... Read in water activities
      read(21,*)
      DO I=1,NNOD
         read(21,555) aw(i)
      enddo
555      FORMAT(f12.8)
c
c
C-----------------------Calculate Total concentration of primary species
c          ctot are total molalities (mol/kwg), we read UT further below
c
      do i=1,nnod
       do n=1,npri
         ctot(i,n)=c(i,n)
       enddo
        do j=1,naqx
            ncp=ncps(j)
            do k=1,ncp
              n=icps(j,k)
              ctot(i,n)=ctot(i,n)+stqs(j,k)*c(i,npri+j)
            enddo
        enddo
      enddo
c
C---------------------------------------------------------------Mineral amount
c
      IF(NMIN .EQ. 0)  GOTO 99
      READ(21,*)                      ! Initial mineral amount
      DO 40 I=1,NNOD
         READ(21,45)  (PINIT(I,N),N=1,NMIN+1)
45       FORMAT(4e21.14)
40    CONTINUE
c
      READ(21,*)                      ! Current mineral amount
      DO 50 I=1,NNOD
         READ(21,55)  (PRE(I,N),N=1,NMIN)
50    CONTINUE
55       FORMAT(4e21.14)
c
C-------------------------------------------------------Gas partial pressure
c
99    CONTINUE
       IF(NGAS .EQ. 0)  GOTO 109
C
      READ(21,*)                      ! gas partial pressure
      DO 150 I=1,NNOD
         READ(21,155)  (PFUG(I,M),M=1,NGAS)
150    CONTINUE
155      FORMAT(6e15.8)
!
!
!.................................................Multi-sites cation exchange
!
109   CONTINUE
!
      IF(NEXC .EQ. 0)  GOTO 199
!
      READ(21,*)                      ! CEC of ion exchange
      READ(21,65)  ((CEC(I,isite), isite=1,NXsites), I=1,NNOD)
!
      READ(21,*)                      ! Exchange amount
      DO 60 I=1,NNOD
         READ(21,65)  ((XCADS(I, isite, N), N=1,NEXC), isite=1,NXsites)
60    CONTINUE
65    FORMAT(6e15.8)
!
!
C-----------------------------------------------------------------Adsorption
199   CONTINUE
c
c--------- residual moles of species after dryout
            read(21,*)
         do i=1,nnod
            read(21,75)(adryr0(i,n),n=1,npri)
         enddo
  75     FORMAT(6e15.8)
         do i=1,nnod
           do n=1,npri
             adryr(i,n) = adryr0(i,n)
           enddo
         enddo
c
c     Add reading UT's separately
      read(21,*)           ! UT  concentration of primary species mol/L
      do i=1,nnod
         read(21,85)  (ut(i,n),n=1,npri)
         do n=1,npri
           utold(i,n)=ut(i,n)
         enddo
      enddo
 85     FORMAT(6e15.8)
!
!-------------------------------------------- CO2 trapped in the solid phase
!
      IF (   IEOS.EQ.2  .OR. IEOS.EQ.13
     &  .OR. IEOS.EQ.14 .OR. IEOS.EQ.15 .OR. IEOS.EQ.16)   THEN
         READ(21,*)
         DO I=1,NNOD
            READ(21,401)  SMco2(I)
       END DO
401     FORMAT(e21.14)
      END IF
c
c     Reads number of grains (used for surface area calculation)
c     only read in if nmkin > 0
      if(nmkin.gt.0)then
        read(21,*)                      ! Number of grains
        do  i=1,nnod
           read(21,"(4e21.14)")  (grains(i,n),n=1,nmkin)
        enddo
      endif
!
!----------------------------------------------------------------------------
!
      RETURN
      END
c
c
c-------------------------------------------------------------------------------
c
c
       SUBROUTINE WRITE_RESTART
C
C*************** Write restart data for reactive chemical transport **********************
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      COMMON/TRANGAS1/PFUGOLD(NMNOD,NMGAS)  ! Old partial pressure
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT      
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/SOLUTE6/SLOLD(MNEL)        ! old liquid saturation
      COMMON/DRYOUT1/adryr(MNOD,MPRI),adryr0(mnod,mpri),
     +   drypre(mnod,mmin)              ! residual in precipitates
      common/water_activity/aw(mnod)    !water activity
         COMMON/Wdensity/Dwat(NMNOD),Dwat_old(NMNOD)    ! Water density (kg/dm**3)
!
      COMMON/SOLIDco2/SMco2(NMNOD)      ! CO2 TRAPPED in solid phase
      COMMON/EOS_INDICATOR/ IEOS        ! Indicate EOS module used
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(15('*'),' WRITE_RESTART 1.0, 2003.7.30: Write restart data'
     X' for reactive chemical transport',15('*'))
C
C-------------------------------------------------------------------------------
c
      WRITE(21,"('    Geochemical state variables at restart time'
     + ' (s):')")
      WRITE(21,'(E14.8)') SUMTIM
C
      WRITE(21,*)  '    Concentration of primary species mol/kgw :'
      DO 20 I=1,NNOD
         WRITE(21,25)  (C(I,N),N=1,NPRI)
25       FORMAT(6e15.8)
20    CONTINUE
C
c... pH
c
      WRITE(21,*)'    pH :'
      DO I=1,NNOD
         WRITE(21,444) ph(i)
      enddo
444      FORMAT(f10.6)
c
      NAQT=NPRI + NAQX
      IF (NAQX .EQ. 0)  GO TO 59
      WRITE(21,*)  '    Concentration of secondary species mol/kgw :'
      DO 27 I=1,NNOD
         WRITE(21,28)  (C(I,N),N=NPRI+1,NAQT)
28       FORMAT(6e15.8)
27    CONTINUE
59    CONTINUE
c
c... Water activities
      write(21,*)  '    Water Activity :'
      DO I=1,NNOD
         write(21,555) aw(i)
      enddo
555      FORMAT(f12.8)
c
C---------------------------------------------------------------Mineral amount
c
      IF(NMIN .EQ. 0)  GOTO 99
      WRITE(21,*)  '    Initial mineral amount mol/dm^3 :'
      DO 40 I=1,NNOD
         WRITE(21,45)  (PINIT(I,N),N=1,NMIN+1)
45       FORMAT(4e21.14)
40    CONTINUE
C
      WRITE(21,*)  '    Current mineral amount mol/dm^3 :'
      DO 50 I=1,NNOD
         WRITE(21,55)  (PRE(I,N),N=1,NMIN)
55       FORMAT(4e21.14)
50    CONTINUE
C-------------------------------------------------------Gas partial pressure
99    CONTINUE
      IF(NGAS .EQ. 0)  GOTO 109
      WRITE(21,*)  '     Gas partial pressure :'
      DO 150 I=1,NNOD
         WRITE(21,155)  (PFUG(I,M),M=1,NGAS)   ! TX on 20-Aug-1999
155       FORMAT(6e15.8)
150    CONTINUE
!
!.................................................Multi-sites cation exchange
!
109   CONTINUE
!
      IF(NEXC .EQ. 0)  GOTO 199
!
      WRITE(21,*)  '   CEC of ion exchange :'
      WRITE(21,65)  ((CEC(I,isite), isite=1,NXsites), I=1,NNOD)
      WRITE(21,*)  '    Exchange amount :'
      DO 60 I=1,NNOD
         WRITE(21,65)  ((XCADS(I, isite, N), N=1,NEXC), isite=1,NXsites)
65       FORMAT(6E15.8)
60    CONTINUE
c
C-----------------------------------------------------------------Adsorption
c
199   CONTINUE
c
c------------ residual concentrations
            write(21,*)  '  Residual Solute after Dryout moles/dm^3 :'
         do i=1,nnod
           write(21,75)(adryr(i,n),n=1,npri)
         enddo
 75    FORMAT(6e15.8)

c       Add saving of UT here, since we changed units of C (mol/kgw)
c and we want to avoid using densities for conversions on restart (in case of dry nodes)
      write(21,*)  ' UT   Concentration of primary species mol/L :'
      do i=1,nnod
         write(21,85) (ut(i,n),n=1,npri)
      enddo
 85     FORMAT(6e15.8)
!
!-------------------------------------------- CO2 trapped in the solid phase
!
      IF (   IEOS.EQ.2  .OR. IEOS.EQ.13
     &  .OR. IEOS.EQ.14 .OR. IEOS.EQ.15 .OR. IEOS.EQ.16)   THEN
         WRITE(21,*)  '  CO2 trapped in the solid phase :'
         DO I=1,NNOD
            WRITE(21,401)  SMco2(I)
       END DO
401      FORMAT(e21.14)
      END IF
c
c     Writes number of grains (used for surface area calculation)
c     only write out if nmkin > 0
      if(nmkin.gt.0)then
        write(21,*)  ' Number of Grains :'
         do  i=1,nnod
         write(21,"(4e21.14)")  (grains(i,n),n=1,nmkin)
        enddo
      endif
!
!
      RETURN
      END
c
c
c
c-------------------------------------------------------------------------------
c
c
c
       subroutine gasdiffus(tkk,ptg,wtmolg,dmolsq,dcfgaspt)
c
c... Gaseous species diffusion coefficient (after Lasaga, 1998,
c....    Kinetic Theory in the Earth Sciences, pg. 322)
c
c      tkk = temperature (Kelvin)
c      p = pressure (kg/m/s^2)
c      wtmolg = molecular weight (g/mol)
c      dmolsq = molecular diameter (m) **2
c      dcfgas = diffusion coefficient of gaseous species (m^2/s)
c      rgas = molar gas constant (m^2 kg/s^2/mol/K)
c      avog = Avogadro's number (1/mol)
c      ppi = pi
c ag1 = rgas/(3*sqrt(2)*pi*avog)
c ag2 = 8*rgas*1000/pi
c
       implicit double precision (a-h,o-z)
       implicit integer*8 (i-n)
       parameter(ag1=1.03585728603761d-24)
       parameter(ag2=2.11727258541919d4)
       double precision tkk,ptg,wtmolg,dmolsq,dcfgaspt
c
c
       SAVE ICALL
       DATA ICALL/0/
       ICALL=ICALL+1
       IF(ICALL.EQ.1) WRITE(11,899)
  899  FORMAT(6X,'gasdiffus 1.0     30 July      2003',6X,
     X 'Calculate gaseous species diffusion coefficient')
c
       dcfgaspt = (ag1*tkk/(ptg*dmolsq))*dsqrt(ag2*tkk/wtmolg)
c
       return
       end
c
c
c-------------------------------------------------------------------------------
c
c
c
c
