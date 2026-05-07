c
c-------------------------------------------------------------------------------
c
      subroutine assign
c
C************* This subroutine assigns logK values to each reaction **************
c***************** also assigns temperature at current node *********
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      integer*8 n
      double precision coef(5)
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **********assign 1.0, 2003.7.30: Assigns logK values to'
     x' each reaction**********')
c
c------------ Calculates activity coefficient parameters at node temperature
        call hkfpar(tc2,tk2)
c
        degc=tc2
        if(degc.gt.tmpmax) then
            degc=tmpmax                 ! add fit constraints
        elseif(degc.lt.tmpmin) then
            degc=tmpmin                 ! add fit constraints
        endif
        tkk = degc + 273.15d0
        tki = 1.d0/tkk
        tksqi = tki*tki
        tksq = 1.d0/tksqi   !ns5/09
        tklog = dlog(tkk)
c-----
c.....calculate temperature factors for T + dT (not yet)
c        dTrct = 1.d-3
c        tkrct = tkk + dTrct
c        tkirct = 1.d0/tkrct
c        tksqirct = tkirct*tkirct
c        tklogrct = dlog(tkrct)
c-----
        do j=1,naqx
          do k=1,5
            coef(k)=akcoes(j,k)
          end do
          aks(j)=fak(coef,tkk,tklog,tki,tksqi)
c
c.....calculate numerical derivative of K for enthalpy of reaction
c          aksrct=fak(coef,tkrct,tklogrct,tkirct,tksqirct)
c          htrctnaks(j) = (aksrct-aks(j))/dTrct
c
c.........Pressure correction
          if(akcops(j,1).ne.0.d0) then
           do k=1,5
             coef(k)=akcops(j,k)
           end do
           aks(j)= aks(j)
     &      + pkcor(coef,tkk,tksq,tki,tksqi,Pt,p0bar)
          endif
        end do
c
        do j=1,nmequ
          do k=1,5
            coef(k)=akcoem(j,k)
          end do
          akm(j)=fak(coef,tkk,tklog,tki,tksqi)
c.........pressure correction
          if(akcopm(j,1).ne.0.d0) then
           do k=1,5
             coef(k)=akcopm(j,k)
           end do
           akm(j)= akm(j)
     &      + pkcor(coef,tkk,tksq,tki,tksqi,Pt,p0bar)
          endif
c
c.........Add temperature dependent supersaturation window ssq0 (in log(K) units)
c         exponential decrease from temp = sst1 to temp = sst2 (1/100 of initial value)
c
          ssq(j)=ssq0(j)
          if(tc2.gt.sst1(j).and.
     +      sst1(j).ne.0.d0.and.sst2(j).ne.0.d0) ssq(j)=ssq0(j)*
     +      dexp( -4.61d0/(sst2(j)-sst1(j)) * (tc2-sst1(j)) )
        end do
c
        nmkin=nmin-nmequ
        do j=1,nmkin
          jkin = nmequ+j
          do k=1,5
            coef(k)=akcoem(jkin,k)
          end do
          akin(j)=fak(coef,tkk,tklog,tki,tksqi)
c.........pressure correction
          if(akcopm(jkin,1).ne.0.d0) then
           do k=1,5
             coef(k)=akcopm(jkin,k)
           end do
           akin(j)= akin(j)
     &      + pkcor(coef,tkk,tksq,tki,tksqi,Pt,p0bar)
          endif
c  
c.........Add temperature dependent supersaturation window ssqk0 (in log(K) units)
c         exponential decrease from temp = sstk1 to temp = sstk2 (1/100 of initial value)
c
          ssqk(j)=ssqk0(j)
          if(tc2.gt.sstk1(j).and.
     +      sstk1(j).ne.0.d0.and.sstk2(j).ne.0.d0) ssqk(j)=ssqk0(j)*
     +      dexp( -4.61d0/(sstk2(j)-sstk1(j)) * (tc2-sstk1(j)) )
        end do
c
        do j=1,ngas
          do k=1,5
            coef(k)=akcoeg(j,k)
          end do
          akg(j)=fak(coef,tkk,tklog,tki,tksqi)
c........ pressure correction
          if(akcopg(j,1).ne.0.d0) then
           do k=1,5
             coef(k)=akcopg(j,k)
           end do
           akg(j)= akg(j)
     &      + pkcor(coef,tkk,tksq,tki,tksqi,Pt,p0bar)
          endif
        end do
c
        do j=1,nads
          do k=1,5
            coef(k)=akcoead(j,k)
          end do
          akd(j)=fak(coef,tkk,tklog,tki,tksqi)
        end do
c
        return
        end
c
c
c-------------------------------------------------------------------------------
c
c
        double precision function fak(coef,tkk,tklog,tki,tksqi)
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        double precision coef(5)
        double precision tkk,tklog,tki,tksqi
c
        fak=coef(1)*tklog+coef(2)+coef(3)*tkk+coef(4)*tki+
     +     coef(5)*tksqi
c
        return
        end

        double precision function pkcor(coef,tkk,tksq,tki,tksqi,
     &           pbar,p0bar)
c*********************************    
c 
c       Returns pressure correction to log(K), which is the second term of 
c       the equation below (i.e., Poynting factor in log10 form):
c
c       log(K)T,P = log(K)T,P0 + (1/2.302585)*(-dV*(P-P0)/(gc*1000*Tk))  
c
c       dV is the volume change (in cm3) for the reaction (Vproducts-Vreactants)
c       therefore negative sign in front of dV !!
c       For gases this corresponds to partial molal volume of pure condensed phase
c       (i.e. for reaction X(gas) <==> X(aq), with dV(gas) = 0 by convention)
c
c       dV = coef(1) + coef(2))Tk + coef(3)*Tk**2 + coef(4)/Tk + coef(5)/Tk**2          
c       dV in cm3!!!
c       gc = gas constant in liter*bar/K/mole (from chempar.inc)
c       gc*1000 is gas constant in cm3*bar/K/mole
c
c       coef regression coefficients read in thermo database
c       tkk = temperature in Kelvin 
c       tki = 1.d0/tkk
c       tksq = tkk*tkk
c       tksqi = tki*tki
c       pbar = current pressure in bar
c       p0bar = reference pressure in bar (sat P of pure water for current logK data) 
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        double precision coef(5)
        double precision tkk,tksq,tki,tksqi
c
        pkcor= -(
     &   coef(1)+coef(2)*tkk+coef(3)*tksq+coef(4)*tki+coef(5)*tksqi
     &   )  *  (pbar-p0bar)/gc/1.d3/tkk/2.302585d0
c
        return
        end
c
c-------------------------------------------------------------------------------
c
      subroutine newtoneq(ielem,densw)
c
C*********** Solve equations of chemical system by Newton-raphson iteration *********
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      integer*8 ielem
      common/matrix/amat(mnr,mnr),bmat(mnr),isat(mmin),
     +    mout(mmin),nmtc,nsat
      common/satgas1/isatg(mgas),moutg(mgas),nsatg !keep track gas saturation
      common/satgas2/sg2
c              porosity*saturation
      common/phisat/phisl1(mnel),phisg1(mnel)
c... ...Dissolution kinetics
      common/disskin/acfdiss(mmin),bcfdiss(mmin),ccfdiss(mmin)
      common/iprkin/ideprec(mmin)
c...... Added common block for rate law designations
      common/irtlaw/nplaw(mmin)
c
      common/minkin3/ cr0(mpri)
      common/drmin/numdr    ! >0 numerical derivatives of mineral kinetics
      COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
      COMMON/E1/ELEM(MNEL)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT  ! for passing time to slud
      COMMON/SOLUTE8/SL1(MNEL)                    ! new liquid saturation
      COMMON/SVZ/NOITE,MOP(24)
C-----------------------------------------------------------------------------
      COMMON/MOP_REACT/MOPR(20) ! controling parameters for reactive transport
c---------------------------------------------------Indicators from EOS module
      COMMON/EOS_INDICATOR/IEOS    ! Indicate EOS module used
C------------------------------------------------------------------
C.......for aqueous kinetics
      common/aqkin2/ntrx      ! total number of redox pair
C
      character*5 elem        ! for node error printout
      character*20 errname    ! save name of non-converging species
      character*8 errtype     ! absolute or relative !NS3/06 chge from 10 to 8
      integer*8 indx(mnr)
      integer*8 no_ch,ntot
      integer*8 nflip(mmin)
      double precision convm(mmin), convg(mgas),convphi ! absolute convergence criterion
      double precision densw
      double precision balm(mpri)    ! chemical mass balance for each component
c
      COMMON/TRANGAS9/NGAS1          ! Number of gaseous species
c----------------------------------------------------------------------------------
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' **********newtoneq 1.0, 2003.7.30: Solve equations of'
     x' chemical system by Newton-Raphson method**********')
c
      no_ch=0                  ! flag=1 to skip chemistry
      errtype = 'relative'
      rt=gc*tk2
      facmax=0.5d0
      IF(NAQX.EQ.0 .AND. NMIN.EQ.0 .AND.
     1  NADS.EQ.0 .AND. NGAS.EQ.0 .AND. NTRX.EQ.0)  FACMAX=1000.d0
      kjacob=0
c.......define absolute convergence criteria for mineral/gas mass action
      do m=1,nmequ
        aklog=dabs(akm(m))
        convm(m)=dmax1(1.d-6*aklog,1.d-6)    ! akm is in log10 form
      end do
      do m=1,ngas
        aklog=dabs(akg(m))
        convg(m)=dmax1(1.d-6*aklog,1.d-6)    ! akg is in log10 form
      end do
      convphi=1.d-6
c
c Determine if gas is saturated here for rare
c cases when the gas become unsaturated during transport
c (i.e. if gas sturation becomes less than sl1min, we do
c not form a separate gas phase).  If no test here and
c amount of gas is zero, it will bomb in cmq_cp
c   (changed such that fixed fugacity works with fully saturated nodes)
c
      if(sg2.le.sl1min.and.ngas1.ne.0) nsatg=0
c
c Now we keep minerals at equilibrium  in the matrix at
c all times, and solve for cm(i) as the amount of mineral dissolving
c (negative) or precipitating (positive) for the current time step.
c If a mineral is not present originally (PRE=0), we flag it to kick it
c out of the jacobian as long as its log(Q/K) is less than 0.
c
      do m=1,nmequ
c         Counter of "flip-flops" for equilibrium minerals with supersat
c         window repetitively going from undersat into the window
         nflip(m)=0
         si2_old(m)=0.d0   ! to save previous iteration logQ/K
         mout(m)=0
         if(pre(ielem,m).le.0.d0) mout(m)=1
c          noreact option
         if(noreact(ielem,m).ne.0) mout(m)=2
      enddo
c  ------------------------------------
c  Start of chemical iteration loop 800
c  ------------------------------------
      do 800 k=1,maxitpch
         iterch=k
         kjacob=kjacob+1
c
c    xh2o variable = kilograms of liquid water
c    see INIT for definition of xh2o.  During iterations, we solve
c    for xh2o and set cp(nw)=1.  After iterations, xh2o is stored
c    in cp(nw).
c
        xh2o = cp(nw)
        cp(nw) = 1.d0
c
c       actualize jacobian
        it=k
        call chemeq(it,no_ch,ielem)
        if(no_ch.ne.0) goto 900  ! abort chemistry
c
c.......Limit flip-flopping problems in/out of supersat window
c       Note: could put loop into existing one at top of jacobeq but cleaner here (?)
c
        do m=1,nmequ
c         logQ/K sign reversal when the division (flip) is negative
          flip=0.d0
c
          if(dabs(si2(m)).gt.convm(m)) flip=si2_old(m)/si2(m)
c
          if(flip.lt.0.d0.and.ssq(m).gt.0.d0) nflip(m)=nflip(m) + 1
          if(nflip(m).gt.5) mout(m)=1
        end do
c
c.......kinetics of mineral dissolution and precipitation
        if (nmkin.gt.0)   then
           if (numdr.eq.0)  then
            call cr_cp(ielem)       ! Analytical derivatives
           else
            call cr_cp_num(ielem)
            call dcr_dcp_num(ielem) ! Numerical derivatives
           end if
        end if
c
        call jacobeq
        kjacob=0
c.......Save mass balances before call to solver for use in convergence criteria
        do i=1,npri
          balm(i)=bmat(i)
        enddo
C-----------------------------------------------------------Call LU solver
        ntot = mnr
        call ludcmp (amat,nmtc,ntot,indx,dd)
        call lubksb (amat,nmtc,ntot,indx,bmat)
        cp(nw) = xh2o     !stores xh2o back into cp(nw)
c--Relative convergence criteria
c  -----------------------------
c  Also need bmat multiplication by cp for increment in newton raphson
        errmax=0.d0
        do i=1,npri                ! for primary species
          errx=bmat(i)
          BMAT(I)=BMAT(I)*CP(I)
c skip convergence test if both tt and balm(=tt-u2) <1e-50
          if(dabs(tt(i)).lt.1.d-50.and.dabs(balm(i)).lt.1.d-50) then
            bmat(i)=0.d0
            errx=0.d0
          endif
          if(dabs(errx).gt.errmax) then
            errmax=dabs(errx)
            errname=napri(i)       !ns98/3 added to save species name
            error=errmax
            errconv=tolch
          end if
        end do
C
c  This is an important if-block to constrain divergence!
c  It is most needed at the first equilibration when all minerals
c  are introduced at once into the system
c
        if (errmax.gt.facmax) then
          fcermx = facmax/errmax
          do i=1,nmtc
            bmat(i)= bmat(i)*fcermx
          end do
        end if
c
c...Moved the entire convergence criteria bloc here, before the
c   update of variables (and change goto 800 to new 750)
c
c--Chemical convergence criteria
c  -----------------------------
c
        if(errmax.lt.tolch) then
c
c  At this point, we passed the relative convergence criteria
c  for aq. species.  Now we make sure mineral/gas mass actions
c  pass the absolute convergence criteria set earlier.
c
          errtype='absolute'
        do m=1,nsat
          lm=isat(m)
c---------Add if-else clause (for supersat "gap")
          if(si2(lm).le.0.d0) then
            if(dabs(si2(lm)).gt.convm(lm)) then
              error=si2(lm)
              errname=namin(lm)
              errconv=convm(lm)
              goto 750    !iterate more
            end if
c---------Added else block below
          else
            if(dabs(si2(lm)-ssq(lm)).gt.convm(lm)) then
              error=si2(lm)-ssq(lm)
              errname=namin(lm)
              errconv=convm(lm)
              goto 750    !iterate more
            end if
          end if
        end do
c
        do m=1,nsatg
          lg=isatg(m)
          if(dabs(sig2(lg)).gt.convg(lg)) then
            error=sig2(lg)
            errname=nagas(lg)
            errconv=convg(lg)
            goto 750   !iterate more
          end if
        end do
c
        if(npot.ne.0) then
          do m=1,npot
            if(dabs(Fphi(m)).gt.convphi) then
              error=bmat(nmtc-nsurf+m)
              errname='potential'
              errconv=convphi
              goto 750   !iterate more
            end if
          end do
        endif

c.......Add last check to make sure chemical mass balance is within acceptable limits
c       Note, in 99.99% of cases, previous constraints will already yield accurate mass balance.
c
        errtype='mass bal'
        error=0.d0
        do i=1,npri
          errconv=max(dabs(tt(i))*tolch,dabs(cr(i))*deltex*tolch,1.d-12)
          diff=dabs(balm(i))
c at very high rates, diff can't catch up to errconv, so we exit if the relative error is near machine precision       
         if(diff.gt.errconv.and.errmax.gt.1.d-12) then
            error=diff
            errname=napri(i)
            goto 750   ! iterate more
          endif
        enddo

        goto 900  !finally converged!
      endif

  750 continue
c
c    actualize mol primary species
c    -----------------------------
        do i=1,npri
           cp(i)=cp(i)+bmat(i)
        enddo
c    actualize mol mineral
c    ---------------------
c---------- Need vliq for mineral calculation
        sumsalts=0.d0   !sum of salts weights in kg (assume zero for now)
        dliq=densw/1000.d0       !liquid density in g/cc (kg/l) (assume 1. for now)
        vliq=1.d0
        gpressure=rt*sl1(ielem)/(vliq*sg2)   !converter from cmg(mol/kg) to pressure (bar)
c
        vphisl1 = vliq/phisl1(ielem)
c
        do m=1,nsat
          bmatm=bmat(npri+m)
          n=isat(m)
          cm(n)=cm(n)+bmatm
c
          premn=pre(ielem,n)    ! NS3/06
          dum=premn*vphisl1
c        if( (cm(n)+dum).le.0.d0 ) then
        if( (cm(n)+dum).le.0.d0.and.si2(n).le.0.0d0 ) then     ! to remove only when both are smaller than zero  gxzh
c---------- Dissolves what is left and remove the mineral
          cm(n) = -dum
          mout(n) = 2  ! removes mineral from jacobian
          errmax=1.d0  ! force one more iteration
        end if
c
        if(cm(n).le.0.0d0.and.premn.le.1.0d-30) cm(n)=0.0d0    ! NS3/06 don't dissolve if there is not solid phase gxzh

        enddo
c
c    actualize mol gas
c    -----------------
c
c     conversion factor to calculate gas pressure
c
        do m=1,nsatg
         bmatg=bmat(npri+nsat+m)*cmg(isatg(m))   ! modified for solving for delta_X/X
         cmg(isatg(m))=cmg(isatg(m))+bmatg
         if(cmg(isatg(m)).le.0.d0)   cmg(isatg(m))=1.d-35  !cmg cannot be zero if nsatg > 0
c
c  cmg = mol_gas/kgw
c  cg = gas partial pressure (bar) for current node (later pfug is set to cg for all nodes) 
c  We convert below cmg in mol_gas/kgw to cg (pfug) in bar.  At high P, this conversion
c  should include Z (compr factor) but we neglect Z in routine couple when back-converting
c  from gp mol_gas/L_medium to pfug in bar, and since P does not change during chem calculations,
c  we are ok (cg and pfug are ok, which is what we save/transport, but gp is sytematically off)

c        skip if ngas1 =0, for fixed pressure at all nodes for all gases
         if (ngas1.gt.0) then
          cg(isatg(m))=cmg(isatg(m))*gpressure   ! pressure calculation, moved /sg2 into gpressure
         end if
c
        enddo
c
c    unsaturated gases - cg was set to gas partial pressure =log(q/k)/gamg in routine cmq_cp
c    note: gases are removed from matrix in jacobeq if gas saturation less than sl1min
c
c    actualize potential term phi
c    -----------------
      if(npot.gt.0) then    !npot=number of surfaces= number of potential terms
        do m=1,nsurf
          kk=npri+nsat+nsatg+m
          jj=ipoten(m)
          errx=dabs(bmat(kk)/phip2(jj))     !no relative increment
          bmatphi=bmat(kk)
          IF (ABS(bmatphi) > 0.1d0) THEN
            bmatphi = SIGN(0.1d0,bmatphi)
          END IF
          phip2(jj)=phip2(jj)+bmatphi
        enddo
      endif
c
c
      if (mopr(3).ge.1) then
        if(ielem.eq.1 .or. ielem.eq.7) then
          write(32,*) 'grid block in Newtoneq,iterch'
          write(32,'(2I5)') ielem, iterch
          write(32,*) 'ngas1 npri nmtc nsatg         sg2      CG(1) '
          write(32,'(4I5,e12.4,2E11.3)') ngas1,npri,nmtc,nsatg,sg2,CG(1)
          write(32,*) '      cmg(1)      errmax        dliq        vliq'
          write(32,'(4E11.3)') cmg(1),errmax,dliq,vliq
          write(32,*) '        --- cp ----  '
          write(32,'(5e12.4)') (cp(nnn),nnn=1,npri)
          write(32,*) '        --- U2 ----  '
          write(32,'(5e12.4)') (U2(nnn),nnn=1,npri)
          write(32,*) '        ---tt ----  '
          write(32,'(5e12.4)') (tt(nnn),nnn=1,npri)
          write(32,*) '        --- cs ----  '
          write(32,'(5e12.4)') (cs(nnn),nnn=1,naqx)
        end if
      end if
c
c.......Moved the entire convergence criteria block from here to further up - cleaner
c
800   continue   !keep iterating
      IF (IEOS.EQ.9 .AND. NO2AQ.GT.0)  GOTO 896
C---------------------------------------------------
c   can only get here if do 800 loop reaches maxitpch
       write(32,"(
     &  /2x,'Warning: chemistry did not converge at node ',a5,
     &   ' (routine NEWTONEQ), Non-convergence type: ',a10,
     &   '  Node temperature (C): ',F10.2,'  Liq.sat.: ',e10.5)")
     &    elem(ielem),errtype,tc2,sl1(ielem)
c
       if(errtype.eq.'relative') then
         write(32,
     &   "( 2x,'Species: ',A10,' Relative error =',e10.4,
     &       '  Tolerance= ',e10.4,
     &       '  Program execution was not aborted.  Check results!')
     &   ") errname,error,errconv
        ireturn=1
c
       else if(errtype.eq.'absolute') then
         write(32,
     &   "( 2x,'Mineral: ',A10,' Mass action error =',e10.4,
     &       '  Tolerance= ',e10.4,
     &       '  Program execution was not aborted.  Check results!')
     &   ") errname,error,errconv
        ireturn=2
c
       else if(errtype.eq.'mass bal') then
         write(32,
     &   "( 2x,'Species: ',A10,' Mass balance error =',e10.4,
     &       '  Tolerance= ',e10.4,
     &       '  Program execution was not aborted.  Check results!')
     &   ") errname,error,errconv
        ireturn=1
       endif
       return
c
896    CONTINUE
       IRETURN=1
       if(errtype.eq.'absolute')IRETURN=2
       RETURN
c
c--We have converged below this point
c
900    continue
       if(no_ch.ne.0) then ! ionic str (calculated in dh_hkf81) too large,through call to chemeq
          ielem=0
          cp(nw) = xh2o     !stores xh2o back into cp(nw)
          return
       endif
c
       do m=1,nmequ
          mout(m)=0
       enddo
       return
       end
c-------------------------------------------------------------------------------
      subroutine jacobeq
c
C********************** Construct Jacobian matrix for chemical system *****************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
      common/matrix/amat(mnr,mnr),bmat(mnr),isat(mmin),
     +  mout(mmin),nmtc,nsat
      common/satgas1/isatg(mgas),moutg(mgas),nsatg !keep track gas saturation
      common/minkin2/ dr(mpri,mpri)
      common/minkin3/ cr0(mpri)
      COMMON/DM/DELTEN,DELTEX,FOR,FORD
      common/satgas2/sg2
!
      COMMON/TRANGAS9/NGAS1      ! Number of gaseous species
      common/co2_gene2/ ico2gt0  ! =1: initial Pco2>0
!
!........For H2 generation by mineral phase using EOS5 module
      common/h2_gene2/ ih2gt0       ! =1: initial Ph2>0
!
      COMMON/EOS_INDICATOR/ IEOS    ! Indicate EOS module used
!
!.......For aqueous kinetics
      common/aqkin2/ntrx     ! total number of redox pair
      common/aqkin11/crx(mpri)
      common/aqkin12/drx(mpri,mpri)
      common/initial/cguess(mpri),icon(mpri)
c
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ***** jacobeq 1.0, 2003.7.30: Construct Jacobian matrix'
     x' for chemical system')
c
c.......extension of the matrix with saturated minerals
      nmtc=npri
      nsat=0
c
c Now we keep minerals at equilibrium in the matrix at
c all times, and solve for cm(i) as the amount of mineral dissolving
c (negative) or precipitating (positive) for the current time step.
c Minerals get kicked off only if within a certain "precipitation gap"
c (between log(Q/K)=0 and some positive log(Q/K) input as ssq), or if
c exhausted (if PRE=0 and log(q/k) < 0), or if very small
c
      do 850 m=1,nmequ
         if(mout(m).eq.2) goto 850
         if((mout(m).eq.1.and.si2(m).le.0.d0 ).or.
     +      (si2(m).gt.0.d0.and.si2(m).le.ssq(m)
     +      .and.cm(m).le.0.d0 ))  then
            cm(m)=0.d0
         else
           nmtc=nmtc+1
           nsat=nsat+1
           isat(nsat)=m
         endif
850   continue
C
      nsatg=0
c---------------------------------------test for saturated gases
c (changed such that fixed fugacity works with fully saturated nodes)
c
      itestgas=0
      if (ieos.eq.2.or.ieos.eq.13.or.
     &    ieos.eq.14.or.ieos.eq.15.or.ieos.eq.16) then
          if(sg2.gt.sl1min) itestgas=1
          if(ico2gt0 .eq. 1)itestgas=1    ! take initial background Pco2
      else
          if(sg2.gt.sl1min.or.ngas1.eq.0) itestgas=1
      end if
!
      if (ieos.eq.5)    then                ! H2 generation for EOS5
         if(sg2.gt.sl1min)  itestgas=1
         if(ih2gt0 .eq. 1)  itestgas=1   ! take initial background Ph2
      else
         if(sg2.gt.sl1min.or.ngas1.eq.0) itestgas=1
      end if
!
       if (itestgas.eq.1)  then
          DO 860 M=1,NGAS
            IF(CG(M).lt.1.0d-100)  GO TO 860
            nmtc=nmtc+1
            nsatg=nsatg+1
            isatg(nsatg)=m
860       continue
        end if
c----------------------------------------------------------------------
c
        ntot=npri+nmequ+ngas+nsurf                    
c
c       initialize amat and bmat terms equal zero
        do 200 i=1,ntot
        do 100 j=1,ntot
           amat(i,j)=0.d0
100     continue
           bmat(i)=0.d0
200     continue
c
c       mass balance equations (npri)
c       the alpha terms (partial derivatives)
c
        do  i=1,npri
          do j=1,npri
            amat(i,j)=du2(i,j)
          end do
        end do
c
        if (nexc .eq. 0) go to 109
        do  i=1,npri
          do j=1,npri
            do k=1,nexc
            do isite=1,NXsites     ! Loop over multi-sites
               amat(i,j)=amat(i,j)+stqx(k,i)*dcxM(isite,k,j)*cp(j)
            end do
          end do
          end do
        end do
c
109     continue
c
        if (nads .eq. 0) go to 119
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
119     continue
c
c.......aqueous kinetic rate
c
        if (ntrx.gt.0) then
           do i=1,npri
c            add if to limit depletion of an aq. kinetic reactant
c            if(dabs(tt(i)).gt.1.d-50) then
             do j=1,npri
               amat(i,j)=amat(i,j)-drx(i,j)*deltex*cp(j)
             end do
c            endif
           end do
        end if
c...........................
c
c.......the alpha terms contributed from kinetic dissolution
c
        if (nmkin.gt.0) then
           do i=1,npri
             do j=1,npri
               amat(i,j)=amat(i,j)-dr(i,j)*deltex
             end do
           end do
        end if
c
c       the independent term (npri)
c
        do i=1,npri
          bmat(i)=u2(i)-tt(i)
        enddo
c
        do m=1,nmequ
          ncp=ncpm(m)
          do n=1,ncp
            i=icpm(m,n)
            bmat(i)=bmat(i)+stqm(m,n)*cm(m)
          enddo
        enddo
c
        do m=1,nsatg
          ncp=ncpg(m)
          do n=1,ncp
            i=icpg(m,n)
            bmat(i)=bmat(i)+stqg(m,n)*cmg(m)
          enddo
        enddo
c
        do m=1,nads
          ncp=ncpad(m)
          do n=1,ncp
            i=icpad(m,n)
            bmat(i)=bmat(i)+stqd(m,n)*cd(m)*xh2o
          enddo
        enddo
!
        do k=1,nexc
           do i=1,npri
           do isite=1,NXsites        ! Loop over multi-sites
              bmat(i) = bmat(i)+stqx(k,i)*cxM(isite,k)
           end do
           end do
        end do
!
        do i=1,npri
          bmat(i)=-bmat(i)
        end do
c
c.......aqueous kinetic rate
c
        if (ntrx.gt.0) then
          do i=1,npri
c            add if to limit depletion of an aq. kinetic reactant 
c            if(dabs(tt(i)).gt.1.d-50) then
             dum=crx(i)*deltex
             bmat(i)=bmat(i)+dum
c            endif
          end do
        end if
c
c...........................
c
c.......the independent term contributed from kinetic dissolution
c
        if (nmkin.gt.0) then
          do i=1,npri
c             bmat(i)=bmat(i)+(cr(i)-cr0(i))*deltex
             dum=cr(i)*deltex
             bmat(i)=bmat(i)+dum
             cr0(i)=cr(i)
          end do
        end if
c
c       the alpha terms (nsat)
        do i=npri+1,npri+nsat       !deriv of min. mass action w/respect to cp aq.
          m=isat(i-npri)
          ncp=ncpm(m)
          do n=1,ncp
            j=icpm(m,n)
c
cpitz            dgamck=0.0d0
cpitz           do k=1,ncp
cpitz             ik=icpm(m,k)
cpitz             dgamck=dgamck+stqm(m,k)*dgamp(ik,j)/gamp(ik)  ! add gamma derivative term pitz
cpitz           enddo
cpitz            amat(i,j)=-stqm(m,n)-dgamck*cp(j)              ! relative increment scheme
c
            amat(i,j)=-stqm(m,n)
c            if(j.eq.nw) amat(i,j) = 0.d0     ! deriv w. resp to xh2o
          enddo
            amat(i,nw)=0.d0
        enddo
c
        do j=npri+1,npri+nsat        ! deriv of mass bal w/respect to n minerals
          m=isat(j-npri)
          ncp=ncpm(m)
          do n=1,ncp
            i=icpm(m,n)
            amat(i,j)=stqm(m,n)
          enddo
        enddo
c
c       the alpha terms (ngas)        ! deriv of gas mass action w/respect to cp aq.
        do i=npri+nsat+1,npri+nsat+nsatg
          m=isatg(i-(npri+nsat))
          ncp=ncpg(m)
          do n=1,ncp
            j=icpg(m,n)
c
cpitz            dgamck=0.0d0
cpitz           do k=1,ncp
cpitz             ik=icpg(m,k)
cpitz             dgamck=dgamck+stqg(m,k)*dgamp(ik,j)/gamp(ik)  !add gamma derivative term pitz Zh 11/04
cpitz           enddo
cpitz            amat(i,j)=-stqg(m,n)-dgamck*cp(j)              !relative increment add gamma derivative term
c
            amat(i,j)=-stqg(m,n)  ! relative increment add gamma derivative term
c            if(j.eq.nw) amat(i,j)=0.d0     ! deriv w/ respect to xh2o
          enddo
            amat(i,nw)=0.d0    ! deriv w/ respect to xh2o
        enddo
c
        do j=npri+nsat+1,npri+nsat+nsatg ! deriv of mass bal w/respect to n gas
          m=isatg(j-(npri+nsat))
          ncp=ncpg(m)
          do n=1,ncp
            i=icpg(m,n)
            amat(i,j)=stqg(m,n)*cmg(m)  !modified for solving for delta_X/X
          enddo
        enddo
c
        if (ngas1.gt.0) then
           do i=npri+nsat+1,npri+nsat+nsatg     ! deriv gas mass action w/ respect to ngas
              m=isatg(i-(npri+nsat))
              amat(i,i)=1.d0                    ! modified for solving for delta_X/X
           end do 
        endif
c
c       the independent term (nsat)
        do i=npri+1,npri+nsat
          m=isat(i-npri)
c
c ssq is the amount of supersaturation allowed (in log(q/k) units)
          if(si2(m).le.0.d0) then
            bmat(i) = si2(m)            ! log form: bmat = log10(q/k)
          else
            bmat(i) = si2(m)-ssq(m)
          endif
        enddo
c
c       the independent term (ngas)
        do i=npri+nsat+1,npri+nsat+nsatg
          m=isatg(i-(npri+nsat))
          bmat(i) = sig2(m)           ! log form: bmat = log10(q/k/p)
        enddo
c
c.....Put surface complexation in the jabobian matrix
c     the alpha terms (potential term)
c
      if(npot.gt.0) then
         nmtc=nmtc+npot
         kkk=npri+nsat+nsatg

         do i=kkk+1,nmtc      ! deriv of potential equilibrium/respect to cp aq.
          kk=i-kkk
          do j=1,npri
            amat(i,j)=dphi_dcp(kk,j)    ! relative increment scheme (deriv includes mult by cp)
c            if(j.eq.nw) amat(i,j) = 0.d0    !deriv w/ respect to xh2o
         enddo
            amat(i,nw) = 0.d0           ! deriv w/ respect to xh2o
        enddo
c
        do j=kkk+1,nmtc  ! deriv of mass bal w/respect to potential term
          kk=j-kkk
          do i=1,npri
            amat(i,j)=dcp_dphi(i,kk)
          enddo
        enddo
c
        do j=kkk+1,nmtc  ! deriv of potential equilibrium w/respect to potential term
         kj=j-kkk
         do i=kkk+1,nmtc
           ki=i-kkk
            amat(i,j)=dphi_dphi(ki,kj)
          enddo
        enddo
c
c       the independent term (nsurf)
        do i=kkk+1,nmtc
          kk=i-kkk
          bmat(i) =-Fphi(kk)      !negative sign, so we add bmats later
        enddo
c
      endif
c.....end addition of potential term in jacobian
c
        return
        end
c
c
c-------------------------------------------------------------------------------
c
c
        subroutine chemeq(it,no_ch,ielem)
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
c-------------------------
        common/satgas2/sg2
        COMMON/MOP_REACT/MOPR(20)  ! controling parameters for reactive transport
        COMMON/min_SI/SIM(MNEL,MMIN) ! Mineral saturation index (log(Q/K)) for all nodes
        integer*8 it,ielem,no_ch
c---------------------
        common/ion_str1/str              !ionic strength
        common/ion_str2/str_node(mnel)    !ionic strength for all nodes
        common/water_activity/aw(mnel)    !water activity
C.......For aqueous kinetics
        common/aqkin2/ntrx     ! total number of redox pair
c----------------------
c Now read in NGAMM       NGAMM=1
c Doesn't seem necessary  IF(IT .LE. 4)  NGAMM=3
c Some difficult problems need this, but fewer is much faster
c        do i=1,NGAMM
        NGMM=1
        IF(IT .LE. 3.or.it.gt.80)NGMM=NGAMM
        do i=1,NGMM
           call dh_hkf81(no_ch)
!
cpitz           if (mopr(9) .eq. 0)   then
              call cs_cp        ! Without derivatives of activity coefficient for DH   
cpitz                                 else
cpitz              call cs_cp_gam    ! Consider derivtives of activity coefficient for Pitzer
cpitz           end if
!
        enddo
!
        str_node(ielem)=str     ! Save ionic strength for all nodes
        aw(ielem)=gamp(nw)      ! Save water activity for vapor pressure lowering calculation Pitz zh
!
        if(no_ch.ne.0) return
!
!.......For intra-aqueous kinetics and biodegradation
!
        if (ntrx .gt. 0)   then
           call cr_cp_rx
           call dcr_dcp_rx
        end if
!
        if (nmequ.gt.0 .or. ngas.gt.0) then
           izero = 0
           call cmq_cp(izero)
c..........save mineral saturation index for all nodes
c          if (mopr(8).eq.1)  then
c          dropped all by one, and moved mbalance to 4
           if (mopr(8).ge.1)  then
              do m=1,nmequ
c              need to print SI before equilibration   sim(ielem,m)=si2(m)
               sim(ielem,m)=si2_old(m)
              end do
           end if
c--------------------------
        end if
!
!.......Cation exchange
!
        if (nexc.gt.0) then
!
            do i=1,npri
              ct(i)=cp(i)
              gamt(i)=gamp(i)
            end do
            do i=1,naqx
              ct(i+npri)=cs(i)
              gamt(i+npri)=gams(i)
            end do
!
!
!...........Loop over multi-sites
!
            do isite=1, NXsites
!
               cec2 = cecM(isite)
               if (cec2 .eq. 0.0d0) then
                  do k=1,nexc
                     cxM(isite,k) = 0.0d0
                  end do
                  go to 794
               end if
!
               do j=1,nexc
                  ekx(j) = ekxM(isite, j)    ! Selectivity
               end do
!
               call cx_ct
               call dcx_dcp
!
               do k=1,nexc
                  cxM(isite,k) = cx(k)
                  do jp=1,npri
                     dcxM(isite,k,jp) = dcx(k,jp)
                  end do
             end do
!
794            continue
            end do  ! multi-sites
!
        end if
!
795     continue
!
!.......................................
!
        if (nads.gt.0)  call admodel
c
        return
        end
c
c
c-------------------------------------------------------------------------------
c
c
        subroutine write_mat(ielem)
c
c *****************************************************************************
c     Subroutine Norm: This subroutine write the jacobin matrix in case of singularity
c                  :
c ******************************************************************************
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
c
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
c
        common/matrix/amat(mnr,mnr),bmat(mnr),isat(mmin),
     +    mout(mmin),nmtc,nsat
        common/satgas1/isatg(mgas),moutg(mgas),nsatg !keep track gas saturation
        common/satgas2/sg2
        COMMON/e1/elem(mnel)
        COMMON/KK1/NPRINT,WTIME,WUPC,DIFUN,TIMETOT  !for passing time to slud
        COMMON/SOLUTE8/SL1(mnel)           ! new liquid saturation
        COMMON/CYC/KCYC,ITER,ITERC,TIMIN,SUMTIM,GF,TIMOUT
        character*5 elem
c
        open(unit=69,file='jaconbin',status='unknown')
c
        write(69,*)
        write(69,*)
        write(69,*) elem(ielem),'        Time Step=',KCYC
        write(69,*) elem(ielem),'        Time(sec)=',timetot
        write(69,*) elem(ielem),'        Liq. Satr=',sl1(ielem)
        write(69,*)
        write(69,*)
c
        write(69,'(15x,100(1x,A15))')
     +   (napri(j),j=1,npri),
     +   (namin(isat(j)),j=1,nsat),
     +   (nagas(isatg(j)),j=1,nsatg),
     +   ' Righthand:    ',' cp_cm_cmg:    '
        write(69,*)
c
        do i=1,npri
        write(69,'(A15,100(1x,E15.9))') napri(i),
     +           (amat(i,j),j=1,nmtc),bmat(i),cp(i)
        enddo
        write(69,*)
c
        do i=npri+1,npri+nsat
        write(69,'(A15,100(1x,E15.9))') namin(isat(i-npri)),
     +           (amat(i,j),j=1,nmtc),bmat(i),
     +                         cm(isat(i-npri))
        enddo
        write(69,*)
c
        do i=npri+nsat+1,nmtc
        write(69,'(A15,100(1x,E15.9))') nagas(isatg(i-npri-nsat)),
     +           (amat(i,j),j=1,nmtc),bmat(i),
     +           cg(isatg(i-npri-nsat))
        enddo
c
        write(69,*)
        write(69,*)
c
        write(69,*) 'Activity coefficents of aqueous species: '
        write(69,'(100(1x,A15))')
     +   (napri(j),j=1,npri),
     +   (naaqx(j),j=1,naqx)
        write(69,'(100(1x,E15.9))') (gamt(j),j=1,npri+naqx)
        write(69,*)
        close(unit=69)
        return
        end
c
c***************************************************************************
       SUBROUTINE  admodel
c***************************************************************************
c      Setup surface complexation model constraints and partial derivatives
c      NS 2/08 Combine routines admodel0, admodel2 and admodel3 from LZ and make
c      use of pointers to run different surfaces and model types simultaneously (with
c      and/or without potential terms)
c
c      Sorption model type for each surface n is:
c        iadmod(n) = 0  no electrical double layer (no Boltzman terms)
c        iadmod(n) = 1  constant capacitance
c        iadmod(n) = 2  diffuse double layer with linear potential function
c        iadmod(n) = 3  diffuse double layer with Gouy-Champan potential
c
c      Note on arrays:
c      There are two array types, one for surfaces, one for non-zero potential calcs:
c        Index 1 to nsurf (total surfaces): s.areas (surfads), model type (iadmod),
c           and all potential terms (phip2) (zero and non-zero)
c        Index 1 to npot (total non-zero potentials): used for all arrays
c           related to the jacobian setup
c      Pointers between arrays are as follows:
c        iad_surf(i=1,nads)  points to index of surface 1 to nsurf
c        iad_phi(i=1,nads)   points to index of non-zero potential terms 1 to npot
c        ipoten(i=1,npot)    points to index of surface 1 to nsurf

        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
        common/ion_str1/str              !ionic strength
!
C-------------------------------------------------Double layer thickness (dm)
c     gc1 is gas constant in J/K/mol, Faraday is Fadaray constant  (chempar_v2.inc)
c        farsq = faraday**2
c     epsi is the multipication of the water dielec constant at 25C=78.5 and
c        the absolute vacuum permittivity=8.854e-12 (defined in chempar_v2.inc)
c        epsigc1 = epsi*gc1
c
C--------------------------------------------------------Calculate alpha term
C----Call CD_CP to obtain the amount of surface complexes
      CALL  CD_CP

      do j=1,npot     !main loop through each surface with non-zero potential term

        jj=ipoten(j)  !points to index in surface arrays

        if(IADMOD(jj).eq.1) then
cels4/30/08          ADFACTOR=FARADAY*FARADAY/(gc1*tk2)/capacitance(jj)   !!!ALFA without specific surface
          ADFACTOR=farsq/(gc1*tk2)/capacitance(jj)   !!!ALFA without specific surface
          dphi_dphi(j,j)=-1.0d0
          chg_term1=-phip2(jj)   !charge density (C/m2) divided by adfactor

        elseif(iadmod(jj).eq.2) then
c         F/capacitance term for model 2 (linear potential function)
cels4/30/08          CAPPAINV=DSQRT(epsi*gc1*tk2/(FARADAY*FARADAY*STR*2.0d3))
cels4/30/08          ADFACTOR=CAPPAINV*FARADAY*FARADAY/(EPSI*gc1*tk2)    !!!ALFA without specific surface
cels5/12/08 added temp variable epsgctk2
            epsgctk2 = epsigc1*tk2
          CAPPAINV=DSQRT(epsgctk2/(farsq*STR*2.0d3))
          ADFACTOR=CAPPAINV*farsq/epsgctk2       !!!ALFA without specific surface
          chg_term1=-phip2(jj)                   !charge density (C/m2) divided by adfactor
          dphi_dphi(j,j)=-1.0d0

        elseif(iadmod(jj).eq.3) then
c         F/capacitance term for model 3 (Gouy-Chapman)
cels4/30/08          adfactor=faraday/DSQRT(8.0d3*epsi*gc1*tk2*STR)
          adfactor=faraday/DSQRT(8.0d3*epsigc1*tk2*STR)
          chg_term1=dsinh(-phip2(jj)*0.5d0)  !charge density (C/m2) divided by adfactor
          dphi_dphi(j,j)=-0.5d0*dcosh(-phip2(jj)*0.5d0)

        endif

c       charge-density term from sum of surface complexes
        chg_term2=0.0d0
        do k=1,nads
          if(iad_phi(k).eq.j) then
             chg_term2=chg_term2+ADFACTOR/surfads(jj)*ZD(K)*CD(K)
          endif
        enddo

cns3/4/2010 add contribution from ads primary species for special cases
c       when the ads primary species is charged
        do kk=1,npads
         kkk=kk+npaq  !primary species index
         if(z(kkk).ne.0.d0) then
          if(isurfp(kk).eq.ipoten(j))
     &       chg_term2=chg_term2+ADFACTOR/surfads(jj)*cp(kkk)*z(kkk)
         endif 
        enddo


c       function to minimize (charge density balance)
        Fphi(j)=chg_term1-chg_term2


c----initialize derivatives
         do i=1,npri
           dcp_dphi(i,j)=0.0d0
           dphi_dcp(j,i)=0.0d0
         enddo
         do i=1,npot
           if(i.ne.j) dphi_dphi(j,i)=0.0d0    !(i,i) terms were defined above for each model
         enddo

C----------------------------calculate dcp_dphi----------------------------
         do k=1,nads
            if(iad_phi(k).eq.j) then
            ncp=ncpad(k)
            zdcdx = zd(k)*cd(k)*xh2o
            do n=1,ncp
              i=icpad(k,n)
cels4/30/08              dcp_dphi(i,j)=dcp_dphi(i,j)+stqd(k,n)*zd(k)*cd(k)*xh2o
              dcp_dphi(i,j)=dcp_dphi(i,j)+stqd(k,n)*zdcdx
            enddo
          endif
         enddo

C----------------------------calculate dphi_dcp----------------------------
         do k=1,nads
         if(iad_phi(k).eq.j) then
            ncp=ncpad(k)
            aszc = ADFACTOR*zd(k)*cd(k)/surfads(jj)
            do n=1,ncp
              i=icpad(k,n)
cels4/30/08              dphi_dcp(j,i)=dphi_dcp(j,i)-ADFACTOR/surfads(jj)
cels4/30/08     1                      *stqd(k,n)*zd(k)*cd(k)   !mult by cp for relative increment!
              dphi_dcp(j,i)=dphi_dcp(j,i)-aszc
     1                      *stqd(k,n)   !mult by cp for relative increment!
            enddo
         endif
         enddo

cns3/4/2010 add derivatives for case of charged ads primary species
         adsur = ADFACTOR/surfads(jj)
         do kk=1,npads
          kkk=kk+npaq  !primary species index
          if(z(kkk).ne.0.d0) then
           if(isurfp(kk).eq.ipoten(j))
     &       dphi_dcp(j,kkk)=dphi_dcp(j,kkk)-
     &          adsur*cp(kkk)*z(kkk)  !mult by cp for rel. incr.
          endif 
         enddo

C----------------------------calculate dphi_dphi----------------------------
         do i=1,npot
           ii=ipoten(i)  !points to surface index
           do k=1,nads
            if(i.eq.j.and.iad_phi(k).eq.j) then
            dphi_dphi(j,i)=dphi_dphi(j,i)-ADFACTOR/surfads(ii)
     &                    *zd(k)*zd(k)*cd(k)
           endif
           enddo
         enddo

      enddo  !loop j through each surface

      CALL DCD_DCP
C
        RETURN
        END
c
c
c********************************************************************************
        subroutine cd_cp
C********************************************************************************
c       Computes concentrations of derived surface complexes (cd, in moles/kgw)
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
c
cc         write (32,*) 'nads'
cc         write (32,'(i5)') nads
!
C------------------------------calculate the amount of surface complexs
      do k=1,nads
c
         cd(k) = -akd(k)
         ncp   = ncpad(k)
c
         do n=1,ncp
            i = icpad(k,n)
            cd(k)=cd(k)+stqd(k,n)*(dlog10(cp(i)*gamp(i)))      !)+dlog10(gamp(i))) saves cpu
         end do

c        include electrostatic effects for models if npot>0
c         if(iadmod(iad_surf(k)).ne.0) then
         do j=1,npot
           jj=ipoten(j)
           if(iad_phi(k).eq.j) then
             cd(k)=cd(k)+zd(k)*phip2(jj)/2.3026d0              !dlog(10.0d0)  saves cpu
           endif
         enddo
c         endif

         cd(k)=10.d0**(cd(k))
c
      end do
C
        return
        end
c
c
c**********************************************************
        subroutine dcd_dcp
c**********************************************************
c       the derivative of cd with regards to pri. sp.

c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
c
      do k=1,nads
         do i=1,npri
           dcd(k,i)=0.d0   !initialization
         enddo
         ncp   = ncpad(k)
cels5/12/08 add temp variable
         cdkxh2o = cd(k)*xh2o
         do n=1,ncp
            i = icpad(k,n)
            if(i.ne.nw) then
cels5/12/08              dcd(k,i)=cd(k)*stqd(k,n)*xh2o        !with relative increment
              dcd(k,i)=cdkxh2o*stqd(k,n)        !with relative increment
            else
cels5/12/08              dcd(k,i)=cd(k)*xh2o                  !with relative increment
              dcd(k,i)=cdkxh2o                  !with relative increment
            endif
         end do
      end do

      return
      end
c
