c
c
      subroutine phichg
c
C************* This routine calculates porosity changes from mineral ***************
c              precipitation/dissolution
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
c
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +              CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      CHARACTER*5 MAT
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
C.....For clay swelling
      common/ion_str2/str_node(mnel)    !ionic strength for all nodes
      common/clay_swell1/ iswell
      common/clay_swell2/ vmin_old(mnel,mmin)  ! previous mole volume for all node
      common/clay_swell3/ vmin0(mmin)          ! initial mole volume
      double precision rtnode(mnel,mmin)
      double precision str_old(mnel)
c-----------------------------------------
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *phichg 1.0, 2003.7.30: Calculate porosity changes from'
     X' mineral precipitation/dissolution**********')
c
C-----------------------------------------Get clay swelling parameters
c
      if (icall.eq.1)  then
         iswell=0
         do n=1,nm
            if (mat(n).eq.'SWELL' .or. mat(n).eq.'swell') then
               iswell=1
               strm=dm(n)     ! Minimum ionic strength to maintain the density
               Denr=por(n)    ! Density reduction in percentage when I = 0
               Dtim=per(1,n)  ! Time needed for the density reduction (in second)
c              never used               rtim=1.0d0/Dtim
             end if
         end do
       if (iswell.eq.1) then
            do m = 1, nmin
               vmin0(m)=vmin(m)
            end do
c
            do i = 1, nnod
               do m = 1, nmin
                  vmin_old(i,m)=vmin0(m)
                str_old(i)=str_node(i)
                rtnode(i,m)=0.0d0
               end do
            end do
         end if
      end if
c
c.... Gives new porosity after mineral dissolution/precipitation
      do i = 1, nnod
c----------------------------------- for clay swelling
      if (iswell.eq.1) then
       str2=str_node(i)      ! ionic strength
         if (str2.ge.strm)  go to 199
       str_old2=str_old(i)
         if (str_old2 .gt. strm)   str_old2=strm
         dtr2=Denr*(str_old2-str2)/strm   ! density reduction factor during one DT
         do m = 1, nmin
            if (namin(m)(1:5) .eq. 'swell'
     +                  .or. namin(m)(1:5) .eq. 'SWELL')   then
                rtnode(i,m)=rtnode(i,m)+dtr2
                if (rtnode(i,m) .gt. Denr) rtnode(i,m)=Denr
                fden=1.0d0-rtnode(i,m)
                vmin(m)=vmin0(m)/fden
                vmin_old(i,m)=vmin(m)
                str_old(i)=str2
          end if
         end do
199      continue
      end if
c-----------------------------------------------------
        vfmtx = 0.d0
        do m = 1, nmin
          vfmtx = pre(i,m)*vmin(m) + vfmtx
        end do
c
c       In case the sum of initial volume do not add up to 1,
c       the proportion of initial unreactive minerals were
c       stored in pinit(node,nmin+1) as vol.unreactive/vol.medium
c       (this was done in routine INIT)
         phi(i) = 1.d0 - vfmtx - pinit(i,nmin+1)
c----------- Keep porosity from going below zero
         phi(i) = dmax1(phi(i),0.d0)
c
      end do
c
      return
      end
c
c
c-------------------------------------------------------------------------------
c
c
      subroutine permchg
c
C************** Calculate permeability changes due to mineral dis./prec. **********
c
c     Modified back to k goes to zero
c     Modified so that permeability only can go down to 1e-20
c     Removed comments on old changes
c     Added two new laws for fractures and matrix
c     Removed dependence on tortuosity and replaced by flag (ikplaw)
c     Added matrix permeability change for porosity
c      less than 0.8 (otherwise cubic law for fractures)
c     matrix permeability change based on simplified Carman-Kozeny
c      Lasaga (1998, pg. 236) Kinetic Theory in the Earth Sciences
c     modifications to catch a zero permeability and porosity
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'perm_v2.inc'
      include 'common_v2.inc'

c      parameter(ppi=3.1415926536d0)
      parameter(pi3=9.42477796076938d0,pi4=12.5663706143592d0)
      parameter(pi6=18.8495559215388d0,pid128=2.45436926061703d-2)
      parameter(tfdpi=0.238732414637843d0)

      COMMON/E1/ELEM(MNEL)
      COMMON/E2/MATX(MNEL)
      COMMON/SOLID/NM,DM(MAXMAT),POR(MAXMAT),PER(3,MAXMAT),
     +              CWET(MAXMAT),SH(MAXMAT)
      COMMON/SOCH/MAT(MAXMAT)
      CHARACTER*5 MAT
      CHARACTER*5 ELEM
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' *permchg 1.0, 2003.7.30: Calculate permeability changes'
     X' due to mineral dis./prec.**********')
c
c... Calculate permeability changes at each grid element
c
      do i = 1, nnod
c
        if(perm(1,i).eq.0.d0.and.perm(2,i).eq.0.d0.and.perm(3,i).
     +    eq.0.d0)goto 100
c
c...... Save the ratio of new to initial porosity as phifact
        if(ikplaw(i).eq.1.or.ikplaw(i).eq.3) then
          if(phi0(i).gt.0.d0)then
             phifact=phi(i)/phi0(i)
             phif3 = phifact**3
          else
             phifact=1.d0
             phif3 = 1.d0
          end if
        end if
c
c........Permeability changes for planar fractures and porous matrix
c........Adjust 3 permeability components equally (k)
c
c........Fracture: cubic law for constant length, width,
c...................................and frequency) - Type 3
         if(ikplaw(i).eq.3) then
           do k = 1, 3
c............Modify only nonzero permeabilities
             if(perm0(k,i).gt.0.d0) then
               perm(k,i) = dmax1(perm0(k,i)*phif3,1.d-30)
             end if
           end do
         end if
c
c........Matrix: simplified Carman-Kozeny relation - Type 1
         if(ikplaw(i).eq.1.and.phi(i).ne.1.d0) then
             phirat = ((1.d0-phi0(i))**2)/((1.d0-phi(i))**2)
           do k = 1, 3
c............Modify only nonzero permeabilities
             if(perm0(k,i).gt.0.d0)then
               perm1 = perm0(k,i)*(phif3)
               perm(k,i) = dmax1(perm1*phirat,1.d-30)
             end if
           end do
         end if
c
c........Matrix:     Verma and Pruess (1988) ------ Type 5
c
         if (ikplaw(i).eq.5.and.phi(i).ne.1.d0)then
!
            phic=aparpp(i)
            power=bparpp(i)
!
            if (phi(i) .gt. phic)    then
!
               phirat = (phi(i)-phic)/(phi0(i)-phic)
               phirat = phirat**power
               do k = 1, 3
                  perm1 = perm0(k,i)
                  perm(k,i) = dmax1(perm1*phirat,1.d-30)
               end do
                                    else
               do k = 1, 3
                  perm(k,i) = 1.0d-30
               end do
!
c              IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14)   THEN       ! for ECO2 module)
c                 CALL WRITE_PLOT_ECO2
c                                  ELSE
c                 CALL WRITE_PLOT
c              END IF
!
c              write (32, 999)   elem(i)
c999           format(/2x, 'Warning: Permeability in grid block: ', a6,
c     &                   ', reaches zero, Execution was stopped')
c              stop
!
            end if
!
         end if
!
!........For UNOCAL scaling cleanup by HF to recover porosity and permeability
!........This is not a general option
!
         if(ikplaw(i).eq.15.and.phi(i).ne.1.d0)then
           phi00=0.5d0            ! initial porosity
           perm00=4.3d-12         ! initial K
           phic=aparpp(i)
           power=bparpp(i)
           phirat = (phi(i)-phic)/(phi00-phic)
           phirat = phirat**power
           do k = 1, 3
             perm1 = perm00
             perm(k,i) = dmax1(perm1*phirat,1.d-30)
           end do
         end if
c--------------------------------------------------------------------------
c
c........Matrix: Fractal pore-space geometry (Pape et al., 1999) - Type 6
c
         if (ikplaw(i).eq.6.and.phi(i).ne.1.d0)then
            df1=aparpp(i)         ! for porosity<=0.1
            df2=bparpp(i)         ! for porosity> 0.1
            phirat = phi(i)/phi0(i)
            if (phi(i).le.0.1d0) phirat = phirat**df1
            if (phi(i).gt.0.1d0) phirat = phirat**df2
            do k = 1, 3
              perm1 = perm0(k,i)
              perm(k,i) = dmax1(perm1*phirat,1.d-30)
            end do
         end if
c
c--------------------------------------------------------------------------
C
c... Matrix: Pore throat model - Type 2 
c
c Cubic packed spheres used to calculate pore throat diameters
c Permeability then calculated using Hagen-Poiseulle equation (modified from Ehrlich et al. 1991)
c     All areas and volumes are in m^2 or m^3
c aparpp = number of effective throats per pore (~2 to 3)
c bparpp = number of pores per unit area
c diam0  = initial pore throat diameter
c areapt = initial area of all pore throats
c phipt  = initial porosity occupied by pore throats
c phipor = initial porosity occupied by pores
c dnpore = total number of pores
c radpor = initial pore radius
c areapr = initial area of all pores
c areatp = initial area of overlap of throats and pores
c atotal = total area of throats and pores
c ddiam  = change in pore throat diameter
c
         if (ikplaw(i).eq.2.and.phi0(i).gt.0.d0)  then
             diamrat = pid128*aparpp(i)*bparpp(i)
             pi3bp = pi3*bparpp(i)
             dnpore = bparpp(i)**1.5d0
             pi4dnp = pi4*dnpore
             pi6dnp = pi6*dnpore
             pfact = aparpp(i)*bparpp(i)*pid128
             phidif = phi(i) - phi0(i)
           do k = 1, 3
c............Modify only nonzero permeabilities
             if (perm0(k,i).gt.0.d0) then
               diam0 = (perm0(k,i)/diamrat)**0.25d0
               areapt = pi3bp*diam0
               phipt  = pi3bp*(diam0*0.5d0)**2
               phipor = phi0(i) - phipt
               radpor = ((phipor/dnpore)*tfdpi)**(0.333333333333333d0)
               areapr = pi4dnp*(radpor**2)
               areatp = pi6dnp*(diam0*0.5d0)**2
               atotal = areapr + areapt - areatp
               ddiam  = phidif/atotal
               deldiam = diam0 + ddiam
               if (deldiam.gt.0.d0) then  
                 perm(k,i) = dmax1(pfact*(deldiam**4),1.d-30)
               else
                 perm(k,i) = 1.d-30
               endif
             endif
           enddo
         endif
c
c.....Fracture: Hydraulic Aperture Model - Type 4 
c
c     Hydraulic aperture model based on cubic law permeability relation
c     All areas and volumes are in m^2 or m^3
c      aparpp = true fracture porosity / fracture-matrix area
c      bparpp = true fracture spacing
c      aper0 = initial calculated hydraulic aperture
c      daper = change in hydraulic aperture
c
         if (ikplaw(i).eq.4.and.phi0(i).gt.0.d0) then
             bp12 = bparpp(i)*12.d0
             phidif = (phi(i) - phi0(i))/phi0(i)
             daper = phidif*aparpp(i)
             nmatx = matx(i)
           do k = 1, 3
             if(per(k,nmatx).gt.0.d0.and.perm0(k,i).gt.0.d0) then
               aper0 = (bp12*perm0(k,i))**(0.333333333333333d0)
c..............One more check so apertures can't go negative
               delaper = aper0 + daper
               if(delaper.gt.0.d0)then  
                 perm(k,i) = dmax1((delaper**3)/bp12,1.d-30)
               else
                 perm(k,i) = 1.d-30
               endif
c
c_only_for v311 ....Minimum changed to six orders of magnitude below mean
c...............fracture permeability (from the rock properties file)
c_only_for v311                perm(k,i) = ((aper0 + daper)**3)/bp12
c_only_for v311                cxfct=perm(k,i)/per(k,matx(i))
c              write(*,*) elem(i),perm0(k,i),per(k,matx(i)),cxfct
c_only_for v311              if(cxfct.lt.1.d-6) then
c_only_for v311                 perm(k,i)=1.d-6*per(k,nmatx)
c_only_for v311              endif
c
             end if
           end do
         end if
c
100     continue
c
      end do
c
      return
      end
c
c
c-------------------------------------------------------------------------------
c
c
      subroutine rsfarea(ig,densw)
c
c... .Calculates reactive surface area in m^2/kg water 
c.....Now uses grain radius to calculate surface area
c.....cut out of treact22.f and modified to use a mineral volume
c.....fraction no smaller than rnucl
c.....added water density dependence 
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/afactorr/a_fmr(mnel)
      common/isarea/imflg2(mmin),imflag(mnel,mmin)
      COMMON/SOLUTE8/SL1(MNEL)           ! new liquid saturation
c            porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c
      double precision densw
      integer*8 ig
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' rsfarea 1.1, 2009.6.10: Calculate reactive surface area')
c
c...........Modified to account for reduced surface area: Factor of
c...          S from Liu et al. (1998) causes S in denominator
c...          to drop out (saturated system, also)
         if(a_fmr(ig).lt.sl1(ig).and.sl1(ig).gt.0.d0.and.
     +      a_fmr(ig).gt.0.d0)  then
c..... Factor based on active fracture model at low saturations,
c....... and for saturations above zero
           actfrc = a_fmr(ig)/(phisl1(ig)*densw)
         else
c..... Factor for saturated system, or unsaturated to consider only
c........ the wetted proportion
           actfrc = 1.d0/(phi(ig)*densw)
         end if
c
c... Volume fraction of rock (V_solid/V_medium)
          vfrtmp = 1.d0-phi(ig)
c
        do m = 1, nmkin
c
c..... Index for kinetic minerals
          imk = nmequ + m
c
c        add unreactive option
         if(noreact(ig,imk).ne.0) then
            amin2(m)=0.d0
c
cns6/10         else
c   add constant rate option
         elseif(imflg2(m).eq.4) then
            amin2(m)=1.d0 
         else
c
c... Set grain radius
          rad2(m) = rad(ig,m)
c
c... Calculate area (m^2/kg H2O) based on grain radius if above zero,
c..... and pre nearly zero
c      if(rad2(m).gt.0.d0.and.pre(ig,imk).lt.1.d-30)then
c        anucl = ppi*0.5d0/rad2(m)
c        amin2(m) = anucl*actfrc  
c        Above assumes anucl in m2_min/m3_med, or in m2_min/m3_min with 
c            rnucl*(1-phi)=1 (thus entire gridblock filled with mineral).
c      else
c... Assign mol/L medium based on mineral amount or by preset minimum volume
c..... fraction (rnucl)
c        rkmol = dmax1(pre(ig,imk),(rnucl(m)*vfrtmp/
c     +    vmin(imk)))
c... Calculate area in m^2/kg H2O
c        amin2(m) = amin(ig,m)*rkmol*vmin(imk)*actfrc
c      end if
c
cxxxxxxxxxxx
c       Changed.  If the initial grain radius (rad) is specified, it is used
c       to compute initial surface area, then increases as the mineral precipitates,
c       resulting in surface area decrease, up to the point when 
c       the surface area equals the input surface area.  
c
c       amin(node,mineral) is input surface area (converted in routine init to m2_mineral/m3_mineral)
c       amin2(node) is current surface area for rate equation (m2_mineral/kg_water)
c       rad(mineral) input initial grain radius (m)
c       rnucl(mineral) input assumed initial volume fraction (vol_mineral/vol_solid) if mineral absent
c       pre(node,mineral) current mineral amount (mol_mineral/L_medium)
c       vmin(mineral)  input molar volume of mineral (converted in routine init to L_mineral/mol_mineral)
c       vfrac  calculated volume fraction of mineral (vol_mineral/vol_medium)  
c       afrc is conversion factor equal to V_medium/kgw
c
c       General relationship is: amin2(m2_miner/kgw) = vfrac(m3_miner/m3_med) * amin(m2_miner/m3_miner) * afrc(m3_med/kgw)
c
          anucl=0.d0
          vfrac0=vfrtmp*rnucl(m)
          vfrac=dmax1(pre(ig,imk)*vmin(imk),vfrac0,1.d-10)

c         Option to calculate a "nucleation" surface area - !!! this is NOT a nucleation model !!!
c         This option is turned on by rad2 non zero - we switch it off when number of grains (.le. 1)
c         as anucl will be quite small
          if(rad2(m).gt.0.d0.and.grains(ig,m).gt.1.d0)then             !flag to turn this option on - use radius to calculate s.area
c
c            Note: grains is number of grains calculated from vfrac and radius, assuming spheres.
c            grains is initially set in routine init.  We reset it here in case a mineral diappeared and reappeared 
             if(pre(ig,imk).eq.0.d0) then
               grains(ig,m)=vfrac*0.125d0/(rad2(m)*rad2(m)*rad2(m)) !reset number of initial mineral grains 
             endif
c
c            Calculate new radius from volume fraction (grains is never zero - default set in init)
             radius= 
     &          ( vfrac*0.125d0/grains(ig,m))**(1.d0/3.d0)   !assumes spheres 

c            S.area in m2_mineral/m3_mineral for cubic packing of spheres; 8 spheres per cube of volume (4*radius)**3 
!............Changed back according to Nic's email of 05/26/10
             anucl = ppi*0.5d0/radius      
c             anucl = 0.5d0/radius
c
c            Decrease number of grains proportionally to cumulative volume fraction change (only on growth)
c            for next go round (this is an arbitrary, smooth decrease yielding desired trend). 
             if(sumqk(m).gt.1.d0) then
                 vfrac_ini=dmax1(vfrac0,1.d-10)
                 grains(ig,m)=grains(ig,m)*vfrac_ini/vfrac
             endif

          endif
c
c         Sum both input and "nucleation" surfaces areas - anucl becomes small so summing both is
c         numerically much smoother than stwiching from anucl to amin when amin becomes greater
c         Note, on dissolution, grains remains constant and anucl decreases with decreasing r
c
          if(imflag(ig,m).eq.3) then
            amin2(m) = ( anucl*vfrac + amin(ig,m))*actfrc
          else
            amin2(m) = ( anucl+amin(ig,m) )*vfrac*actfrc
          endif
c
         endif    !end of else for noreact option
c
cxxxxxxx
c
        end do
c
c.....Add endmember surface areas for solid solutions
      call add_surf(ig)  
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine levscale
c
c     Leverett scaling separated out from propchg
c     pcfact is used for capillary pressure modification
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'perm_v2.inc'
      include 'common_v2.inc'
      COMMON/SOLID/NM,DROK(maxmat),POR(maxmat),PER(3,maxmat),
     +   CWET(maxmat),SH(maxmat)
      COMMON/E2/MATX(MNEL)
c 	Passed tortuosity so scaling can be done for fractures
c 	differently than for matrix. Fractures have a very small tortuosity
c 	because the tortuosity is multiplied by the true porosity. The
c 	scaling for fractures has a 1/3 power, while capillaries have a
c	 1/2 power     
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' levscale 1.0, 2003.7.30: Leverett scaling due to'
     X' permeability changes')
c
      do i = 1, nnod
c
c... Material index at each node
        nmat=matx(i)
c
c... Save the ratio of new to old porosity as phifact
c
        if(por(nmat).gt.0.d0)  then
           phifact=phi(i)/por(nmat)
        else
           phifact=0.d0
        end if
c
c... Save the ratio of new to old permeabilities as perfact
c    (note: 3 values for 3 permeability components - all values
c    above should be the same if non-zero
         perfact=0.d0
       do k = 1, 3
         if(dabs(per(k,nmat)).gt.0.d0)then
            ratio = perm(k,i)/per(k,nmat)
         else
            ratio = 1.d0
         endif
            ratio = dmax1(ratio,perfact)
            perfact = ratio
       enddo
c
c... Capillary pressure correction based
c    on function by Slider (1976) pg 280.
       if (phifact.gt.0.d0.and.perfact.gt.0.d0)  then
         pcfact(i) = dsqrt(phifact/perfact)
       else
         pcfact(i) = 1.d0
       end if
c
      end do
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine add_surf(ig)
c     Sum the surface areas for endmembers of solid solutions
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/solsol/iss(mmin),ncpss(msol),icpss(msol,mcpss),nss
!
      do n=1,nss
           nem=ncpss(n)    !no. of endmembers in solid solution n
           asum = 0.d0
c          changes here and below to sum only areas of endmembers present
           asum1 = 0.d0 
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
              asum=asum+amin2(i)                            !sum regardles if present or not
              if(pre(ig,m).ge.1.d-30) asum1=asum1+amin2(i)  !sum only mineral present, use 1d-30 for consistency with rsfarea
           enddo
           do k=1,nem
              m=icpss(n,k)  !mineral index of endmember k in solid sol n
              i=m-nmequ
              if(asum1.ne.0.d0) then
                 amin2(i)=asum1
              else
                 amin2(i)=asum
              endif
           enddo
      enddo
      return
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
      SUBROUTINE TDS_From_REACT 
! 
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      This routine computes total dissolved solid from -REACT        *
!*                                                                     *
!*                  Version 1.0 - August 30, 2007                      *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
!
      common/co2_gene/ nco2
      common/TDS_REACT1/ iTDS_REACT
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/P1/X((MNK+1)*MNEL)
!
!.....Total dissolved solid mass fraction and CO2 concentration
      common/dissolved_solid/ TDS(mnel)              
      common/dissolved_CO2/ ctot_co2(mnel)                 
!
!.....Molecular weight of primary sopecies, g/mol                     
      common/molweight/wm_aqt(maqt)     ! Read from the database          
!
! -------
! ... Integer variables
! -------
! 
      integer*8  i, n, nel, npri, nskip, nloc
!
! -------
! ... Real variables
! -------
! 
      REAL*8 TDS2
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>>=>  Main body
!
!
!.....Loop over grid blocks
!
      do i=1,nel
!
         nloc=(i-1)*nk1
         TDS2 = 0.0d0
!
!........Loop over components
!
         do n =1, npri+naqx
!
            nskip = 0
            if (n .eq. nw   )  nskip = 1
            if (n .eq. nh   )  nskip = 1
            if (n .eq. noh  )  nskip = 1
            if (n .eq. ne   )  nskip = 1
            if (n .eq. no2aq)  nskip = 1
!
            if (nskip .eq. 0)   then
               TDS2 = TDS2 + c(i, n)*Wm_Aqt(n)    ! In g/kg H2O
            end if
!
         end do
!
!........Mass fraction of total dissolved solid
!
         TDS(i) = TDS2/(TDS2 + 1000.0d0)   
!
!.....Pass the total dissolved solid to TOUGH2 primary variables
!
         if (iTDS_REACT .eq. 1)   then
            x(nloc+2) = TDS(i)
         end if
!
!........Total dissolved CO2 concentrations
!
         ctot_CO2(i) = ctot(i, nco2)
!
!
      end do
!
!
      return

      end SUBROUTINE TDS_From_REACT
!      
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE Henry_Constant_From_REACT(
     &              nel,       ! Actual  number of grid blocks
     &              nagas,     ! Name of gaseous species 
     &              utot,      ! Dissolved component, mol/kg
     &              pfug)      ! Gas partial pressure, bar
! 
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*         This routine computes Henry constants from -REACT           *
!*                                                                     *
!*                  Version 1.0 - July 03, 2007                        *     
!*                                                                     *
!***********************************************************************
!***********************************************************************
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
!
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
!.....Gases related basis species number
!
      common/co2_gene/ nco2
      common/h2_gene/ nh2
      common/gas_gene1/ nch4, nh2s, nso2
      common/ichemcons/naqx,nmequ,nmin,ngas,ne,nw,nh,
     +   noh,no2aq,nd,nads,nexc,nx,iex           ! only : no2aq, ngas
!
      common/Henry_React/Kh_REACT(mnel,mgas)  ! Henry constants calculated from -REACT
      common/gas_index/ichem(18)       ! No-cond gas index in chemical input  
! 
! -------
! ... Integer variables
! -------
! 
      integer*8  i, ig, nel       
!
! -------
! ... Real variable
! -------
! 
c     not used      REAL*8 TDS2
!
! -------
! ... Real array
! -------
! 
      REAL*8 pfug(mnel,mgas)  ! New partial pressure =pfug(  )
!
! -------
! ... Character array
! -------
! 
      character*20  nagas(mgas)
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>>=>  Main body
!
!
!----------------
!.....For ECO2 or ECO2N flow module
!----------------
!
      if (ieos .eq. 13 .or. ieos .eq. 14)        then
!
!........Loop over grid blocks
!
         do i=1,nel
!
!...........Calculate Henry constants
!
            Kh_REACT(i,1) = utot(i,nco2)/pfug(i,1)
!
         end do  ! Loop over grid blocks
!
!
!----------------
!.....For TMVOC1 or TMVOC2 flow module
!----------------
!
      else if (ieos .eq. 15 .or. ieos .eq. 16)   then
!
!........Loop over grid blocks
!
         do i=1,nel
!
!...........Loop over gases
!
            do ig=1,ngas
!
               icg = ichem(ig)             ! Non-cond gas Index
!
               Ngas_basis = 0
!
               if (nagas(ig) .eq. 'CO2'    .or.
     &             nagas(ig) .eq. 'CO2(g)' .or.
     &             nagas(ig) .eq. 'co2(g)'      )      then 
                      Ngas_basis = nco2     ! The related basis species number
               end if
!
               if (nagas(ig) .eq. 'O2')   Ngas_basis = no2aq    
               if (nagas(ig) .eq. 'CH4')  Ngas_basis = nch4 
               if (nagas(ig) .eq. 'H2S')  Ngas_basis = nh2s 
               if (nagas(ig) .eq. 'SO2')  Ngas_basis = nso2 
!
!..............Calculate Henry constants
!
               Kh_REACT(i,icg) = utot(i,Ngas_basis)/pfug(i,ig)                  
!
            end do  ! Loop over gases
!
         end do  ! Loop over grid blocks
!
!
      end if     ! ieos
!
!
      return

      end SUBROUTINE Henry_Constant_From_REACT
!
!
!
!***********************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!
!
!
      SUBROUTINE CO2H2O_ReactionSource 
     &    (deltex, Fkg, nmequ, nmkin, Mpri, Mmin, nco2, nw, cm, 
     &     ncpm, icpm, rkin2, stqm,  Rco2_i, Rh2o_i)
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Obtain CO2 reaction source terms for feeding back to flow      *
!*           equations (only flow EOS2 and ECO2 flow modules)          *
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
      INTEGER*8  nco2, nw, nmequ, nmkin, Mpri, Mmin
      INTEGER*8  m, n, k, ncp, nmqp1
!
! -------
! ... Integer arrays
! -------
! 
      INTEGER*8  ncpm(mmin)
      INTEGER*8  icpm(mmin,mpri)
!
! -------
! ... Double precision variables
! -------
! 
      REAL*8     deltex, Fkg, dum, Rco2_i, Rh2o_i
! 
! -------
! ... Double precision arrays
! -------
! 
      REAL*8     cm(mmin), rkin2(mmin) 
      REAL*8     stqm(mmin,mpri)
!
! 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>
!
!
      Rco2_i = 0.0D0
      Rh2o_i = 0.0D0
!
! ----------
!.....Equilibrium mineral contribution        
! ----------
!
      do m=1,nmequ
!
         dum = cm(m)         ! cm in mol/kg h2o
         ncp = ncpm(m)
!
         do k=1,ncp
            n = icpm(m,k)
            if (n .eq. nco2)  Rco2_i = Rco2_i - stqm(m,k)*dum
            if (n .eq. nw)    Rh2o_i = Rh2o_i - stqm(m,k)*dum
         end do

      end do
!
! ----------
!.....Kinetic mineral contribution        
! ----------
!
      do m=1,nmkin
!
         dum   = rkin2(m)*deltex
         nmqp1 = m + nmequ
         ncp   = ncpm(nmqp1)
!
         do k=1,ncp
            n = icpm(nmqp1,k)
            if (n .eq. nco2)  Rco2_i = Rco2_i + stqm(nmqp1,k)*dum
            if (n .eq. nw)    Rh2o_i = Rh2o_i + stqm(nmqp1,k)*dum
         end do
!
      end do
!
! ----------
!.....For CO2 reaction source terms
! ----------
!
      Rco2_i = Rco2_i*Fkg       ! mol/kg h2o --> mol/dm**3 medium
      Rco2_i = Rco2_i*44.4D0    ! mol/dm**3 ----> kg/m**3 medium
      Rco2_i = Rco2_i/deltex    ! ---> kg/m**3/s
!
! ----------
!.....For H2O reaction source terms
! ----------
!
      Rh2o_i = Rh2o_i*Fkg       ! mol/kg h2o --> mol/dm**3 medium  
      Rh2o_i = Rh2o_i*18.0D0    ! mol/dm**3 ----> kg/m**3 medium
      Rh2o_i = Rh2o_i/deltex    ! ---> kg/m**3/s
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
      SUBROUTINE H2O_ReactionSource 
     &    (deltex, Fkg, nmequ, nmkin, Mpri, Mmin, nw, cm, ncpm,
     &     icpm, rkin2, stqm, Rh2o_i)
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Obtain H2O reaction source terms for feeding back to flow      *
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
      INTEGER*8  nw, nmequ, nmkin, Mpri, Mmin
      INTEGER*8  m, n, k, ncp, nmqp1
!
! -------
! ... Integer arrays
! -------
! 
      INTEGER*8  ncpm(mmin)
      INTEGER*8  icpm(mmin,mpri)
!
! -------
! ... Double precision variables
! -------
! 
      REAL*8     deltex, Fkg, dum, Rh2o_i
! 
! -------
! ... Double precision arrays
! -------
! 
      REAL*8     cm(mmin), rkin2(mmin) 
      REAL*8     stqm(mmin,mpri)
!
! 
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>
!
!
      Rh2o_i = 0.0D0
!
! ----------
!.....Equilibrium mineral contribution        
! ----------
!
      do m=1,nmequ
!
         dum = cm(m)         ! cm in mol/kg h2o
         ncp = ncpm(m)
!
         do k=1,ncp
            n = icpm(m,k)
            if (n .eq. nw)    Rh2o_i = Rh2o_i - stqm(m,k)*dum
         end do

      end do
!
! ----------
!.....Kinetic mineral contribution        
! ----------
!
      do m=1,nmkin
!
         dum   = rkin2(m)*deltex
         nmqp1 = m + nmequ
         ncp   = ncpm(nmqp1)
!
         do k=1,ncp
            n = icpm(nmqp1,k)
            if (n .eq. nw)    Rh2o_i = Rh2o_i + stqm(nmqp1,k)*dum
         end do
!
      end do
!
! ----------
!.....For H2O reaction source terms
! ----------
!
      Rh2o_i = Rh2o_i*Fkg       ! mol/kg h2o --> mol/dm**3 medium  
      Rh2o_i = Rh2o_i*18.0D0    ! mol/dm**3 ----> kg/m**3 medium
      Rh2o_i = Rh2o_i/deltex    ! ---> kg/m**3/s
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
      SUBROUTINE H2_ReactionSource (deltex,  Fkg, nmequ, nmkin, 
     &                              Mpri, Mmin, namin, rkin2, Rco2_i)
!
!
!***********************************************************************
!***********************************************************************
!*                                                                     *
!*      Obtain H2 reaction source terms for feeding back to flow       *
!*              equations (only for EOS5 flow modules)                 *
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
      INTEGER*8  nmequ, nmkin, Mpri, Mmin
      INTEGER*8  i, m, nkkn, m_iron, n_iron
      INTEGER*8  icall
      SAVE       icall
      DATA       icall/0/
!
! -------
! ... Double precision variables
! -------
!
      REAL*8     deltex, Fkg, Rco2_i
! 
! -------
! ... Double precision arrays
! -------
! 
      REAL*8     rkin2(mmin) 
!
! -------
! ... Character array
! -------
! 
      CHARACTER*20   namin(mmin)
!
!
!  =>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>>=>=>=>=>=>=>
!
!
      icall = icall + 1
!
      if (icall .eq. 1)     then

         do m=1,nmkin      
!      
            nkkn   = nmequ + m
!
            if (namin(nkkn)(1:4) .eq. 'iron')   then
               m_iron = m
               n_iron = nmequ + m_iron
            end if
!
         end do
!
      end if
!
! ----------
!.....Current H2 amount in mol/l medium
! ----------
!
      Rco2_i = rkin2(m_iron)         ! Fe, mol/kg h2o, dissolution                   
      Rco2_i = Rco2_i*Fkg            ! Fe, mol/l medium                   
      Rco2_i = Rco2_i*4.0d0/3.0d0    ! H2, mol/l medium                  
      Rco2_i = Rco2_i*2.016D0        ! H2, mol/dm**3 ----> kg/m**3 medium                  
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
      SUBROUTINE Gas_ReactionSource
     &              (deltex,      ! Time step, s
     &               Ngas_basis,  ! The gas related basis species number
     &               Gmolw,       ! Gas molecular weight, g/mol
     &               Fkg,         ! Factor for unit converion (mol/kg h2o to mol/dm**3 medium)
     &               Rgas_i)      ! Reaction source/sink for grid block i, gas ig
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
      include 'common_v2.inc'
!
!
!--------
!.....Equlibrium minerals
!--------
!
      Rgas_i = 0.0D0
!
      do m=1,nmequ
         dum = cm(m)          ! cm in mol/kg h2o
         ncp = ncpm(m)
         do k=1,ncp
            n = icpm(m,k)
            if (n.eq.Ngas_basis) Rgas_i = Rgas_i - stqm(m,k)*dum
         end do
      end do
!
!--------
!.....Kinetic minerals
!--------
!
      do m=1,nmkin
         dum   = rkin2(m)*deltex
         nmqp1 = m + nmequ
         ncp = ncpm(nmqp1)
         do k=1,ncp
            n = icpm(nmqp1,k)
            if (n.eq.Ngas_basis) Rgas_i = Rgas_i + stqm(nmqp1,k)*dum
         end do
      end do
!
!--------
!.....Unit conversion
!--------
!
      Rgas_i = Rgas_i*Fkg       ! mol/kg h2o --> mol/dm**3 medium
      Rgas_i = Rgas_i*Gmolw     ! mol/dm**3 ----> kg/m**3 medium
      Rgas_i = Rgas_i/deltex    ! ---> kg/m**3/s
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
      SUBROUTINE QLOSS_Rgases
!
!
!....Modify residual terms due to gas generation by reactions
!....This subroutine is called when TMgas flow module is used
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
      include 'common_v2.inc'
!
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
!
      COMMON/ICO2/ICO2H2O  ! CO2 and H2O reaction sources considered in the flow
      COMMON/REACTh2o/Rh2o(NMNOD)      ! H2O REACTION SOURCES
!
!.....For gas reaction sources for TMgas module
!
      common/gas_index/ichem(18)       ! No-cond gas index in chemical input
      common/gas_gene2/ Rgases(MNEL,6)
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(6X,'QLOSS_Rgases 1.0    22 April   1999',6X,
     &'Modify residual terms for CO2 reaction sources')
!
      do N=1,NEL
!
         NLOC=(N-1)*NEQ
!
         do ig =1,ngas
!
            icg = ichem(ig)          ! Index in TMgas module
!
            NLco2 = NLOC + icg + 1   ! The first component is water vapor
!
            R(NLco2) = R(NLco2) - Rgases(N,ig)*DELTEX  ! soucre, kg/m**3 medium
!
            if (Ico2h2o .eq. 2)  then
               NLh2o = NLOC+1                          ! H2O is component 1
               R(NLh2o) = R(NLh2o) - Rh2o(N)*DELTEX    ! soucre, kg/m**3 medium
            end if
!
         end do
!
      end do
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
         SUBROUTINE QLOSS_Rco2
C
C
C-------Modify residual terms due to CO2 and H2 generation by reactions
C       This subroutine is called when EOS2, ECO2, ECO2N, and EOS5 modules is used
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
C
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
C
c----------------------------------------Indicators from EOS module
      COMMON/EOS_INDICATOR/ IEOS  
      COMMON/CO2M_TMVOC/ ico2m, iTMVOC   
C
C------------------------------------------------------------------
      COMMON/ICO2/ICO2H2O  ! CO2 and H2O reaction sources considered in the flow
      COMMON/REACTco2/Rco2(NMNOD)    ! CO2 REACTION SOURCES
      COMMON/REACTh2o/Rh2o(NMNOD)    ! H2O REACTION SOURCES
C
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(6X,'QLOSS_Rco2  1.0     22 April   1999',6X,
     X'Modify residual terms for CO2 reaction sources')
!
!
      DO 5 N=1,NEL
         NLOC=(N-1)*NEQ
!
         IF (IEOS .EQ. 2 .OR.         ! for EOS2
     &       IEOS .EQ. 5)     THEN    ! for EOS5
            NLco2=NLOC+2
         ELSE IF (IEOS .EQ. 13 .OR. IEOS .EQ. 14)   THEN
            NLco2=NLOC+3     ! In ECO2   CO2 is component 3
         END IF
!
         R(NLco2)=R(NLco2)-Rco2(N)*DELTEX       ! soucre, kg/m**3 medium
!        R(NLco2) store H2 generation when EOS5 is used
!
!........Iron corrsion generates one mol of H2 and
!........consume one mole of H2O (2 over -18)
         if (ieos .eq. 5)    then
             Rh2o(N) = -9.0d0*Rco2(N)           ! Rco2 stores H2 here
         end if
!
         IF (ico2m .eq. 1 .and. Ico2h2o .EQ. 2 .or.
     &       ieos  .eq. 5)                             THEN
            NLh2o=NLOC+1                       ! H2O is component 1
            R(NLh2o)=R(NLh2o)-Rh2o(N)*DELTEX   ! soucre, kg/m**3 medium
         END IF
    5 CONTINUE
!
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
         SUBROUTINE QLOSS_Rh2o
!
!........Modify residual terms due to H2O generation by reactions
!........This subroutine can be called when EOS1, or EOS3, EOS4, EOS7 
!........modules is used
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      INCLUDE 'flowpar_v2.inc'
      INCLUDE 'chempar_v2.inc'
!
      COMMON/P4/R(MNEQ*MNEL+1)
      COMMON/NN/NEL,NCON,NOGN,NK,NEQ,NPH,NB,NK1,NEQ1,NBK,NSEC,NFLUX
      COMMON/BC/NELA
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/REACTh2o/Rh2o(NMNOD)      ! H2O REACTION SOURCES
!
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(11,899)
  899 FORMAT(6X,'QLOSS_Rh2o  1.0     22 April   2008',6X,
     &'Modify residual terms for H2O reaction sources')
!
      DO N=1,NEL
         NLOC=(N-1)*NEQ
         NLh2o=NLOC+1                       ! H2O is component 1
         R(NLh2o)=R(NLh2o)-Rh2o(N)*DELTEX   ! soucre, kg/m**3 medium
      END DO
!
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
      subroutine tortcalc(nel,satphas,tortphas)
c
c..... Calculation of tortuosity for gas or liquid (satphas = SL or SG)
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'perm_v2.inc'
c
      COMMON/SOLI/COM(maxmat),EXPAN(maxmat),CDRY(maxmat),
     +    TORT(maxmat),GK(maxmat)
      COMMON/E2/MATX(MNEL)
      common/E4/phi(mnel)
      common/torpar/ptort(maxmat),phicrit(maxmat)
      double precision satphas(mnel),tortphas(mnel)
      double precision tortn
      integer*8 nel,nmat
c
c... Original Millington - Quirk (TORT = 0.0) and new generalized 
c..... power law (after Burnol & Claret, 2009; Lagneau, 2002) with 
c      saturation power added (els 11/18/09)
c
      do n = 1, nel
        nmat = matx(n)
        tortn = tort(nmat)
c
c.... Unsaturated
        if(satphas(n).gt.0.d0.and.satphas(n).lt.1.d0)then
c
c.... Millington - Quirk: tort = 0.0
          if(tortn.eq.0.d0)then
            tortphas(n) = (phi(n)**0.333333D0)*(satphas(n)**2.333333D0)
c
c.... Generalized power law: tort < 0.0
          elseif(tortn.lt.0.d0)then
            tortphas(n) = -tortn*(satphas(n)**2.333333D0)*
     +        (((phi(n)-phicrit(nmat))/(phi0(n)-phicrit(nmat)))
     +        **ptort(nmat))
          else
            tortphas(n) = tortn
          endif
c
c.... Fully saturated (gas or liquid)
        elseif(satphas(n).eq.1.d0)then
c
          if(tortn.eq.0.d0)then
            tortphas(n) = phi(n)**0.333333D0
          elseif(tortn.lt.0.d0)then
            tortphas(n) = -tortn*
     +        (((phi(n)-phicrit(nmat))/(phi0(n)-phicrit(nmat)))
     +        **ptort(nmat))
          else
            tortphas(n) = tortn
          endif
c
c.... Zero Saturation (gas or liquid)
        else
            tortphas(n) = 0.d0
c
        endif
c
      enddo
c
      return
      end
c
c
