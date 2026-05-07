c
c
      subroutine readtherm_hkf
c
c Subroutine to replace database.f
c Reads the thermodynamic data base.  For each derived species, mineral,
c and gas, the first record contains soichiometry and logK's, and the second
c contains regression coefficients of logK's a a function ot temperature.
c
c Main variables (these are defined in common.inc and paramete.inc)
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
c   icps(j,i),icpm(j,i),icpg(j,i),icpad(j,i) index of component species (stored order i)  
c      in j derived aq. species, minerals, gases, and ads. species
c   stqs(j,i) stoichiometric coefficient of component (in order i) in derived aq. species j
c   stqm(j,i) stoichiometric coefficient of component (in order i) in mineral j
c   stqg(j,i) stoichiometric coefficient of component (in order i) in gas  j
c   stqd(j,i) stoichiometric coefficient of component (in order i) in ads. species j
c   a0(i) D-H radii of primary and secondary aq. species i=1,npri+naqx
c   z(i) ionic charge of primary and secondary  aq. species   i=1,npri+naqx
c   akcoes(j,5) log10(K) f(T) regression coefficient for reaction j=1,naqx+nmin+ngas+nads
c     (5 coefficients, base=1)
c   akcoem(j,5), akcoeg(j,5), akcoead(j,5) same as above for derived aq. species, 
c       minerals, gases, and ads. species
c   vmin(i)  molar volume of mineral i=1,nmin
c   zd(i)  surface charge (?) for adsorption "species" i=1,nads
c
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc' 

      common/aqxs/iaqxs ! =1: define secondary species in CHEMICAL.INP
      common/thermodat/thermo_in,pitzdata
!
!.....Molecular weight of all species, g/mol                     
      common/molweight/Wm_Aqt(maqt)           
!
!.....LogK pressure dependence
      dimension akcoes2(5),iaqxf(maqx)
      dimension akcops2(5), coefp(5)  
      dimension tempa(ntmp,ntmp),w(ntmp,ntmp),vb(ntmp)
!
      dimension tempc(ntmp),nam(mpri),s(mpri),coef(5),icheck(maqx),   
     &  stoic(mpri),alogk(ntmp),indx(ntmp)
c---------------------------------------------------------------------------
      character*20 name, nam, dummy,namexp(1000)
      character*20 naprit(mpri)    ! list of primary species in input file
      character*20 thermo_in,pitzdata
      character*20 label
      character*1000 inprec        ! input record
c************************************      
c
      data iunit2,iunit5/42,45/
      data nfit/5/
      SAVE ICALL
      DATA ICALL/0/
      ICALL=ICALL+1
      IF(ICALL.EQ.1) WRITE(34,899)
  899 FORMAT(' ****readtherm_hkf, 2008.2.2: Read chemical thermodynamic'
     X' database*************')
      iflgstop=0  ! set to 1 to stop program 
c
c.....Initialize arrays
c
      do i=1,mmin
        icheck(i)=0    
      enddo
      do i=1,maqx
        do j=1,npri
          stqs(i,j)=0.d0
        enddo
      enddo  
      do i=1,nmin
        do j=1,npri
          stqm(i,j)=0.d0
        enddo
      enddo  
      do i=1,ngas
        do j=1,npri
          stqg(i,j)=0.d0
        enddo
      enddo  
      do i=1,nads
        do j=1,npri
          stqd(i,j)=0.d0
        enddo
      enddo  
c
       do j=1, npri
          label(1:20)=napri(j)
          call name_conv(label)
          napri(j)=label(1:20)      
       enddo
       do j=1, naqx
          label(1:20)=naaqx(j)
          call name_conv(label)
          naaqx(j)=label(1:20)      
       enddo
c
       do j=1, nmin
          label(1:20)=namin(j)
          call name_conv(label)
          namin(j)=label(1:20)      
       enddo
       do j=1, ngas
          label(1:20)=nagas(j)
          call name_conv(label)
          nagas(j)=label(1:20)      
       enddo
c
c  The name of the thermo database is read in file SOLUTE.INP
c--Opens the thermodynamic database
c      write(32,*) '   --> Thermodynamic database: thermok.dat'
c      write(iunit2,*) '   --> Thermodynamic database: thermok.dat'
c
       open(unit=iunit5,file=thermo_in, status='old',err=500)
       write(32,*) '   --> Start reading thermodynamic database: ', 
     &    thermo_in
       write(iunit2,*) '   --> Start reading thermodynamic database: ',
     &    thermo_in
c
c--Skips records until we find end-of-header record
2      read (iunit5,"(a20)",end=1) dummy
       if(dummy(1:14).ne.'!end-of-header') goto 2
       goto 3
1      continue
       write(*,*) 
     & ' Cannot find end-of-header record in thermo database'
       write(iunit2,*) 
     & ' Cannot find end-of-header record in thermo database'
       write(32,*) 
     & ' Cannot find end-of-header record in thermo database'
       stop

3      continue
c
c--Reads temperature data.  These are the temperatures for
c  -----------------------
c     which the logK data are given, and over which the log K
c     data were regressed.  Calculation of logKs from input regression
c     coefficients will be restrained to the min and max temperature
c     listed here.
      read(iunit5,*,err = 6001) name,ntemp,(tempc(l),l=1,ntemp)  !name is dummy here
      if (ntemp .gt. ntmp) then
        write(32,*) '  Too many temperature points in database!'
        write(32,*) '  Maximum allowed is: ', ntmp 
        write(iunit2,*) '  Too many temperature points in database!'
        write(iunit2,*) '  Maximum allowed is: ', ntmp 
        stop
      endif
c     stores min and max temperature
      tmpmin = tempc(1)
      tmpmax = tempc(ntemp)
c
      call pre_reg(tempc,ntemp,nfit,ntmp,indx,tempa,w)
c   
c--Reads component species
c  -----------------------
      if(npri.gt.mpri) then
       write(32,*) '  Too many components are specified in this run'
       write(32,*) '  Maximum allowed is:', mpri
       write(iunit2,*) 'Too many components are specified in this run'
       write(iunit2,*) 'Maximum allowed is:', mpri
       stop
      endif
      iprit=0     ! total primary species counter
c
c--- start of loop through primary species  
   10 continue
c
       read(iunit5,"(a1000)",err=6002) inprec
c
c       if (inprec(2:5).eq.'null'.or. 
c     &      inprec(2:4).eq.'end'
       if (index(inprec(1:20),"'null'").ne.0.or. 
     &     index(inprec(1:20),"'end'").ne.0
     &    ) goto 20                  ! end of components
c
       read(inprec,*,err=6002,end=6002) name,aa0,zz,dmolwt
c       name -> napri component species name
c       a0       D-H parameter
c       zz -> z  charge
c
c             modified logic/edits down to line "naqxp=0" to avoid 
c             unecessary use of local arrays and remove these arrays
       label(1:20)=name(1:20)
       call name_conv(label)
       name(1:20)=label(1:20)      
c
       nis = 0
       do i=1,npri
        if(name.eq.napri(i)) then 
           iprit=iprit+1
           z(i)=zz     
           zsqi(i)=zz*zz
           zabsi(i)=dabs(zz)
           zabterm(i) = 0.19d0*(zabsi(i)-1.d0)
           a0(i)=aa0          
           omegahkf(i) = 1.66027d+5*zsqi(i)/a0(i)
           rex=1.81d0
           if(z(i).lt.0.d0) rex=1.91d0
           azero(i) = 2.d0*(a0(i) + zabsi(i)*rex)/(zabsi(i)+1.d0)
           azero3(i) = azero(i)**3
!
           Wm_Aqt(i) = dmolwt
           if (dmolwt .eq. 999.999d0)   then
               Wm_Aqt(i) = 0.0d0
           end if
!
           naprit(iprit)=name     ! array to check later if we found all names
              if   (napri(i).eq.'co2(aq)' .or. napri(i).eq.'CO2(aq)'
     &           .or. napri(i).eq.'ch4(aq)' .or. napri(i).eq.'CH4(aq)'
     &           .or. napri(i).eq. 'h2(aq)' .or. napri(i).eq. 'H2(aq)'
     &           .or. napri(i).eq.'h2s(aq)' .or. napri(i).eq.'H2S(aq)'
     &           .or. napri(i).eq. 'o2(aq)' .or. napri(i).eq. 'O2(aq)'
     &           .or. napri(i).eq.'so2(aq)' .or. napri(i).eq.'SO2(aq)'
     &               )   then
                   nis = nis+1
                   indx_so(nis)=i
              end if
           go to 10
        endif
       enddo            
c
       go to 10     
c--- end of loop through primary species
c
20    continue
c
       nerror=0
       do j = 1, npri
          iyes=0
          do jj=1,npri
             if (napri(j).eq.naprit(jj)) iyes=1
          end do  
c		 
          if (iyes.eq.0)  then
             write(32,*) 'Error: The primary species:    ',
     +          napri(j),'    is not found in the thermo database'                    
             nerror=nerror+1
          end if   
       end do
       if (nerror.ge.1)  then
          write(32,*) 'Execution aborted - stop   '
          stop          
       end if   
c
c--Reads derived species (automatically selects species that contain all components)
c  ---------------------
      naqxp=0      ! all possible secondary species
      if (iaqxs.eq.0) then      ! for automatic slection of secondary species
         naqx=1
                      else
         do ix=1,naqx
            iaqxf(ix)=0  ! the selected species in CHEMICAL.INP is not found in the database
         end do   
      end if
c
c*******************************************
c--- loop through derived species starts here
   30 continue

c      initialize stoichiometry coefficientsreads  
       do j=1,npri
         stoic(j)=0.d0       
       enddo 
c      reads  a0, z, and stoichiometries 
       read(iunit5,"(a1000)",err=6004) inprec

c        if (inprec(2:5).eq.'null'.or. 
c     &      inprec(2:4).eq.'end'
       if (index(inprec(1:20),"'null'").ne.0.or. 
     &     index(inprec(1:20),"'end'").ne.0
     &     ) goto 50        ! end of derived species list
c
        if (inprec(1:2).eq.'* '
     &       .or.inprec(1:1).eq.'#'.or.inprec.eq.'') go to 30    ! skips record and continues cycling

       read(inprec,*,err=6003,end=6003) name,Xmwt,aa0,zz,ncp,
     &                 (s(n),nam(n),n=1,ncp)
c       write(32,*) ' Read: ', name,'stoichiometry'
C-----------xmwt: molecular weight
c
c------------ check if it is a possible secondary species
       label(1:20)=name(1:20)
       call name_conv(label)
       name(1:20)=label(1:20)      

       do i=1,ncp
        label(1:20)=nam(i)
        call name_conv(label)
        nam(i)=label(1:20)      
       enddo

       do n=1,ncp
         iyes=0
         do j=1,npri
            if (nam(n) .eq. napri(j)) iyes=1
         end do
         if (iyes .eq. 0) then 
            read(iunit5,*)
            read(iunit5,*)
c            write(32,*) '  - Skip this species, because', 
c     +                       ' one basis species is not found'
            go to 30   ! skip and reads next species if one component is not found
         end if
       end do
c
       naqxp=naqxp+1        ! all possible secondary species
       namexp(naqxp)=name
c
c     --------------------------------------------------------------------------------
       if (iaqxs.eq.0) then      ! for automatic slection of secondary species
         naaqx(naqx)=name
         a0(npri+naqx)=aa0
         z(npri+naqx)=zz
!
         Wm_aqt(npri+naqx) = Xmwt                    
         if (Xmwt .eq. 999.999d0)   then                     
            Wm_aqt(npri+naqx) = 0.0d0
         end if
!
         zsqi(npri+naqx)=zz*zz
         zabsi(npri+naqx)=dabs(zz)
         ncps(naqx)=ncp
           zabterm(npri+naqx) = 0.19d0*(zabsi(npri+naqx)-1.d0)
           omegahkf(npri+naqx) = 1.66027d+5*zsqi(npri+naqx)/
     +       a0(npri+naqx)
           rex=1.81d0
           if(z(npri+naqx).lt.0.d0) rex=1.91d0
           azero(npri+naqx) = 2.d0*(a0(npri+naqx) + 
     +        zabsi(npri+naqx)*rex)/(zabsi(npri+naqx)+1.d0)
           azero3(npri+naqx) = azero(npri+naqx)**3
           if   (naaqx(naqx).eq.'co2(aq)' .or. naaqx(naqx).eq.'CO2(aq)'
     &      .or. naaqx(naqx).eq.'ch4(aq)' .or. naaqx(naqx).eq.'CH4(aq)'
     &      .or. naaqx(naqx).eq. 'h2(aq)' .or. naaqx(naqx).eq. 'H2(aq)'
     &      .or. naaqx(naqx).eq.'h2s(aq)' .or. naaqx(naqx).eq.'H2S(aq)'
     &      .or. naaqx(naqx).eq. 'o2(aq)' .or. naaqx(naqx).eq. 'O2(aq)'
     &      .or. naaqx(naqx).eq.'so2(aq)' .or. naaqx(naqx).eq.'SO2(aq)'
     &          )   then
              nis = nis+1
              indx_so(nis)=npri+naqx
           end if
       end if   
c
       read(iunit5,*,err=6003) dummy,(alogk(nn),nn=1,ntemp)   
c       write(32,*) ' Reading: ', name,'logK record'
c
c ns5/05 add coefficients akcop to calculate logK as function of pressure
c        akcop is used to input dV in cm3!! 
c        dv = akcop(1) + akcop(2)*Tk + akcop(3)*Tk**2 + akcop(4)/Tk +akcop(5)*/Tk**2  
       do n=1,5
        akcoes2(n)=0.d0
        akcops2(n)=0.d0
       enddo
       read(iunit5,"(a1000)",err=6004) inprec
       read(inprec,*,end=41,err=6003) dummy,(akcoes2(nn),nn=1,5),
     &       (akcops2(nn),nn=1,5)
c
c keep input regression coefficients if given, but regress
c log K values if all coefficients are zero   
c 
   41  ireg=0    ! flag=0 if all input regression coeffs = 0
       do n=1,5
        if(akcoes2(n).ne.0.d0) ireg=1
       enddo
c
       if(ireg.eq.0) then    ! computes logK regression coefficients
c
        call evaluate_logk(alogk,ntmp,ntemp,tempc,name)
        call regression_logK(alogk,tempa,w,indx,nfit,ntemp,ntmp,vb)
        call evaluate_regression(vb,ntmp,ntemp,tempc,name)

       endif
c
       if (iaqxs.eq.0) then        ! case when all secondary species are automatically picked up     
          if(ireg.eq.0) then       ! computed regression coeffs
             do n=1,5
               akcoes(naqx,n)=vb(n)
             end do
          else                     ! input regression coeffs
             do n=1,5
               akcoes(naqx,n)=akcoes2(n)
             enddo
          endif
c
c.........for p correction
          do n=1,5
            akcops(naqx,n)=akcops2(n)
          enddo
c
        do 40 n=1,ncp  ! loops over components n in the stoichiometry of species naqx
          do j=1,npri
            if(nam(n).eq.napri(j)) then 
               stqs(naqx,n)=s(n)            ! saves component soichiometry coef 
               icps(naqx,n)=j               ! species index in stoichiometry
               goto 40           
            endif                            
          enddo
   40   continue      ! loop 
c
c       we get here only if all components were found
        naqx=naqx+1                         ! increments naqx   
        if(naqx.gt.maqx) then
         write(32,*) '  Too many derived aq. species in THERMOK.DAT'
         write(32,*) '  Maximum allowed is:', maqx
         write(iunit2,*) 'Too many derived aq. species in THERMOK.DAT'
         write(iunit2,*) 'Maximum allowed is:', maqx
         stop
        endif
       end if
c 
c     -------------------------------------------------------------------------------
       if (iaqxs.eq.1) then  ! case when secondary species defined in CHEMICAL.INP
         do jx=1,naqx
            if (name .eq. naaqx(jx)) then 
               iaqxf(jx)=1       ! the selected species is found in the database
               a0(npri+jx)=aa0
               z(npri+jx)=zz
               zsqi(npri+jx)=zz*zz
               zabsi(npri+jx)=dabs(zz)
               ncps(jx)=ncp
               zabterm(npri+jx) = 0.19d0*(zabsi(npri+jx)-1.d0)
               omegahkf(npri+jx) = 1.66027d+5*zsqi(npri+jx)/
     +            a0(npri+jx)
               rex=1.81d0
              if(z(npri+jx).lt.0.d0) rex=1.91d0
              azero(npri+jx) = 2.d0*(a0(npri+jx) + 
     +          zabsi(npri+jx)*rex)/(zabsi(npri+jx)+1.d0)
              azero3(npri+jx) = azero(npri+jx)**3
           if   (naaqx(jx).eq.'co2(aq)' .or. naaqx(jx).eq.'CO2(aq)'
     &      .or. naaqx(jx).eq.'ch4(aq)' .or. naaqx(jx).eq.'CH4(aq)'
     &      .or. naaqx(jx).eq. 'h2(aq)' .or. naaqx(jx).eq. 'H2(aq)'
     &      .or. naaqx(jx).eq.'h2s(aq)' .or. naaqx(jx).eq.'H2S(aq)'
     &      .or. naaqx(jx).eq. 'o2(aq)' .or. naaqx(jx).eq. 'O2(aq)'
     &      .or. naaqx(jx).eq.'so2(aq)' .or. naaqx(jx).eq.'SO2(aq)'
     &          )   then
              nis = nis+1
              indx_so(nis)=npri + jx 
           end if
c
               if(ireg.eq.0.d0) then  ! use computed reg. coeffs if not given
                  do n=1,5
                    akcoes(jx,n)=vb(n)
                  end do
               else                   ! input regression coeffs
                   do n=1,5
                    akcoes(jx,n)=akcoes2(n)
                   enddo
               endif
c..............for p correction
               do n=1,5
                 akcops(jx,n)=akcops2(n)
               enddo
c
               do 48 n=1,ncp                 ! loops over components n in the stoichiometry of species naqx
                  do j=1,npri
                   if(nam(n).eq.napri(j)) then 
                    stqs(jx,n)=s(n)          ! saves component soichiometry coef 
                    icps(jx,n)=j             ! species index in stoichiometry
                    goto 48           
                   end if                            
                  end do
   48          continue      ! loop 
               go to 30
            end if   
         end do
       end if   
c
      goto 30   
c--- end of loop through secondary species
c
   50 continue
c
      if (iaqxs.eq.0) then     
         naqx=naqx-1    ! number of derived aq. species
                      else
         do jx=1,naqx
            if (iaqxf(jx) .eq. 0)  then
               write(32,*)  
               write(32,*) 'Error: The defined secondary species:    ',
     +          naaqx(jx),'               is not found in the database'
               stop
            end if
         end do   
      end if
c
c     create a permanent order for the names of aq. species 
      do i=1,npri
       naaqt(i)=napri(i)
      end do
      do j=1,naqx
       naaqt(npri+j)=naaqx(j)
      end do
c
      write(32,*)
      write(32,*) ' --- List of all possible aqueous complexes --- '
      do ix=1,naqxp
         write(32,'(10x,a20)') namexp(ix)
      end do
c
      if (iaqxs.eq.1) then      
       write(32,*)
       write(32,*) ' --- Aqueous complexes selected -------------- '
       do ix=1,naqx
          write(32,'(10x,a20)') naaqx(ix)
       end do
       write(32,*) ' ---------------------------------------------- '
       write(32,*)
      endif
c
c***********************************************************************
c
c--Reads minerals
c  --------------
c
      if(nmin.gt.mmin) then
        write(32,*) '  Too many minerals are specified in this run'
        write(32,*) '  Maximum allowed is:', mmin
        write(iunit2,*) 'Too many minerals are specified in this run'
        write(iunit2,*) 'Maximum allowed is:', mmin
        stop
      endif
c
c---  loop for minerals starts here   
   70 continue
       read(iunit5,"(a1000)",err=6004) inprec
c
c       if (inprec(2:5).eq.'null'.or. 
c     &      inprec(2:4).eq.'end'
       if (index(inprec(1:20),"'null'").ne.0.or. 
     &     index(inprec(1:20),"'end'").ne.0
     &     ) goto 80                     ! end of mineral list

       if (inprec(1:2).eq.'* '
     &       .or.inprec(1:1).eq.'#'.or.inprec.eq.'') go to 70    ! skips record and continues cycling
c
       read(inprec,*,err=6004,end=6004) name,dmwm,vol,ncp,
     +     (s(n),nam(n),n=1,ncp)
c       name -> namin   mineral name
c       vmin   molar volume in cm3/mole (converted later to l/mole)
c       ncp->ncpm    number of stoichiometric components
c       s -> stqm      stoichiometric coefficient
c       nam    name of component species in stoichiometry  
c
       label(1:20)=name(1:20)
       call name_conv(label)
       name(1:20)=label(1:20)      
c
       do j=1,ncp
        label(1:20)=nam(j)
        call name_conv(label)
        nam(j)=label(1:20)      
       enddo
c
       read(iunit5,*,err=6004) dummy,(alogk(nn),nn=1,ntemp)   ! read logK record
c      reads regression coefficients for logK as f(T)     
       do n=1,5
         coef(n)=0.d0
         coefp(n)=0.d0  
       enddo
       read(iunit5,"(a1000)",err=6004) inprec
       read(inprec,*,end=112,err=6004) dummy,(coef(n),n=1,5),  
     &       (coefp(nn),nn=1,5)  
c       
c keep input regression coefficients if given, but regress
c log K values if all coefficients are zero    
  112   ireg=0    !flag=0 if all input regression coeffs = 0
        do n=1,5
         if(coef(n).ne.0.d0) ireg=1
        enddo
c
c       checks if mineral is in the list of minerals specified for this run; if yes store data
        do m = 1,nmin              ! loops over specified minerals in the system
          if(name.eq.namin(m)) then       ! found the mineral
            ncpm(m)=ncp 
            do 60 n=1,ncp                 ! loops over components in mineral
              do j=1,npri                 ! checks to see if all components are in the system
                if(nam(n).eq.napri(j)) then
                  stqm(m,n) = s(n)            ! number of stoichio components
                  icpm(m,n)=j                 ! species index
                  go to 60                    ! resumes looping over components
                endif
              enddo
              write(32,*) '  Mineral:', name  ! gets to this point only if no match is found - skip and read next mineral
              write(32,*) '  Components missing:', nam(n), 'STOP'  
              write(iunit2,*) '  Mineral:', name  ! gets to this point only if no match is found - read next mineral
              write(iunit2,*) '  Component missing:', nam(n), 'STOP'  
              stop
   60       continue   ! loop 
c
            if(ireg.eq.0) then   ! if no regression coefficients were given
             call evaluate_logk(alogk,ntmp,ntemp,tempc,name)
             call regression_logK(alogk,tempa,w,indx,nfit,ntemp,ntmp,vb)
             call evaluate_regression(vb,ntmp,ntemp,tempc,name)
             do n=1,5           
                akcoem(m,n)=vb(n)
             enddo
            else                 ! if regression coefficients  are input
             do n=1,5  
                akcoem(m,n)=coef(n)
             enddo
            endif
c
c...........coefs for logK pressure dependency
            do n=1,5  
                akcopm(m,n)=coefp(n)
            enddo
c
            vmin(m)=vol*1.d-3       ! convert from cm3/mole to l/mole
            dmolwm(m)=dmwm
            if (vmin(m).eq.0.d0) then
              write(32,*) 'zero molar volume for mineral ',namin(m)
              write(32,*) 'check the database'
              write(iunit2,*) 'zero molar volume for mineral ',namin(m)
              iflgstop=1     !to stop later
            endif
            write(32,*) namin(m), vmin(m)
            icheck(m)=1
          endif
        enddo              
      goto 70 
c--- end of loop for minerals
c 
   80 continue
c
      do m=1,nmin                 ! cheks which minerals were not found     
        if(icheck(m).ne.1) then
          write(32,*) ' Mineral:', namin(m)
          write(32,*) ' Not found in database - stop'
          write(iunit2,*) ' Mineral:', namin(m)
          write(iunit2,*) ' Not found in database - stop'
          stop
        endif
      enddo
c
c--Reads gas species
c-------------------
      if(ngas.gt.mgas) then
        write(32,*) '  Too many gases specified in this run'
        write(32,*) '  Maximum allowed is:', mgas
        write(iunit2,*) '  Too many gases specified in this run'
        write(iunit2,*) '  Maximum allowed is:', mgas
        stop
      endif  
c
c--- start of loop through gases    
  200 continue
c
       read(iunit5,"(a1000)",err=6004) inprec
c
c       if (inprec(2:5).eq.'null'.or. 
c     &      inprec(2:4).eq.'end'
       if (index(inprec(1:20),"'null'").ne.0.or. 
     &     index(inprec(1:20),"'end'").ne.0
     &     ) goto 100                     !end of mineral list

       if (inprec(1:2).eq.'* '
     &       .or.inprec(1:1).eq.'#'.or.inprec.eq.'') go to 200     
c
       read(inprec,*,err=6005,end=6005) name,dmolwt,dmdiam,ncp,
     +     (s(n),nam(n),n=1,ncp)
c       name -> nagas   gas name
c       vbargas   molar volume (not used - calculated independently later)
c       ncp    number of stoichiometric components
c       s -> stqg      stoichiometric coefficient
c       nam    name of component species in stoichiometry  
c
       label(1:20)=name(1:20)
       call name_conv(label)
       name(1:20)=label(1:20)      
c
       do j=1,ncp
         label(1:20)=nam(j)
         call name_conv(label)
         nam(j)=label(1:20)      
       enddo
c
       read(iunit5,*,err=6004) dummy,(alogk(nn),nn=1,ntemp)   !read logK record
c       reads regression coefficients for logK as f(T)     
       do n=1,5
         coef(n)=1.d0
       enddo
       read(iunit5,"(a1000)",err=6005) inprec
       read(inprec,*,end=201,err=6005) dummy,(coef(n),n=1,5),  
     &     (coefp(n),n=1,5)      !ns5/09
c
c keep input regression coefficients if given, but regress
c log K values if all coefficients are zero    
  201   ireg=0    !flag=0 if all input regression coeffs = 0
        do n=1,5
         if(coef(n).ne.0.d0) ireg=1
        enddo
c       
c       checks if gas is in the list of gases specified for this run; if yes store data
        do m = 1,ngas                      ! loops over specified gases in the system
          if(name.eq.nagas(m)) then        ! found the gas
            ncpg(m)=ncp                    ! number of components in stoichio
            dmwgas(m) = dmolwt
            diamol(m) = dmdiam
            do 90 n=1,ncp                  ! loops over components in gas
              do j=1,npri                      
                if(nam(n).eq.napri(j)) then
                  stqg(m,n) = s(n)         ! stoichio coeff
                  icpg(m,n) = j            ! index of species in stoichio
                  go to 90                 ! resumes looping over components
                endif
              enddo
              write(32,*) '  Gas:', name   ! gets to this point only if no match is found - read next mineral
              write(32,*) '  One or more components missing - stop'  
              write(iunit2,*) '  Gas:', name  ! gets to this point only if no match is found - read next mineral
              write(iunit2,*) '  One or more components missing - stop'  
              stop
   90       continue   ! loop 
c
            if(ireg.eq.0) then  ! if no regression coeffs are input, calculate them
             call evaluate_logk(alogk,ntmp,ntemp,tempc,name)
             call regression_logK(alogk,tempa,w,indx,nfit,ntemp,ntmp,vb)
             call evaluate_regression(vb,ntmp,ntemp,tempc,name)
             do n=1,5  
                 akcoeg(m,n)=vb(n)
             enddo
            else                ! use input coefficients
             do n=1,5  
                 akcoeg(m,n)=coef(n)
             enddo
            endif
c
c...........coefs for logK pressure dependency
            do n=1,5  
                akcopg(m,n)=coefp(n)
            enddo
c
            icheck(m)=1
          endif
        enddo              
      goto 200 
c--- end of loop through gases
c
  100 continue 
      do m=1,ngas      
        if(icheck(m).ne.1) then
          write(32,*) ' Gas:', nagas(m)
          write(32,*) ' Not found in database - stop'
          write(iunit2,*) ' Gas:', nagas(m)
          write(iunit2,*) ' Not found in database - stop'
          stop
        endif
      enddo
c
c--Reads adsorbed species
c------------------------
c     skips reading surface complexes if not needed
      if(npads.eq.0.and.nads.eq.0) go to 600

      iauto_read=0                ! default, assumes we already read species in chemical.inp
      if(nads.eq.0) iauto_read=1  ! we will automatically read species in the database 
      ns_tot=0 
c--- start of loop through adsorbed species (surface complexes)
  300 continue
c
       read(iunit5,"(a1000)",err=6004) inprec
c
c       if (inprec(2:5).eq.'null'
c     &   .or.inprec(2:4).eq.'end') goto 120                    ! end of mineral list
       if (index(inprec(1:20),"'null'").ne.0.or. 
     &     index(inprec(1:20),"'end'").ne.0) goto 120

       if (inprec(1:2).eq.'* '
     &       .or.inprec(1:1).eq.'#'.or.inprec.eq.'') go to 300  ! skips record and continues cycling

       read(inprec,*,err=6006,end=6006) name,zz,ncp,
     &               (s(n),nam(n),n=1,ncp)
c       name -> naads   ads. species name
c       ncp    number of stoichiometric components
c       s -> stqd      stoichiometric coefficient
c       nam    name of component species in stoichiometry  
c       zd     charge
c
        label(1:20) = name
        CALL name_conv(label)
        name = label(1:20)      
c
         do n=1,ncp
            label(1:20) = nam(n)
            CALL name_conv(label)
            nam(n) = label(1:20)      
         end do
c
       label(1:20)=name(1:20)
       call name_conv(label)
       name(1:20)=label(1:20)      
c
       do j=1,ncp
       label(1:20)=nam(j)
       call name_conv(label)
       nam(j)=label(1:20)      
       enddo
c
       read(iunit5,*,err=6004) dummy,(alogk(nn),nn=1,ntemp)   ! read logK record
c
c      reads regression coefficients for logK as f(T)     
       do n=1,5
        coef(n)=0.d0
       enddo
       read(iunit5,"(a1000)",err=6006) inprec
       read(inprec,*,end=301,err=6006) dummy,(coef(n),n=1,5)  ! ordered in terms of nsec!!!
c
c keep input regression coefficients if given, but regress
c log K values if all coefficients are zero    
  301   ireg=0    !flag=0 if all input regression coeffs = 0
        do n=1,5
         if(coef(n).ne.0.d0) ireg=1
        enddo
c
c       checks if all components in species are part of the chemical system
        do n=1,ncp
           iyes=0
           do j=1,npri
            if (nam(n) .eq. napri(j)) iyes=1
           end do
           if (iyes .eq. 0) go to 300  ! if not all components found, skips to next species 
        end do
c
        ns_tot=ns_tot+1
        namexp(ns_tot)=name
c
        if(iauto_read.eq.1) then  ! if we did not list surface complexes and pick all of them automatically
c       ------------------------------------------------------------------------------
c
         nads=ns_tot
         if(nads.gt.mads) then
             write(32,*) '  Too many adsorption species specified in',
     &            ' this run'
             write(32,*) '  Maximum allowed is:', mads
             write(iunit2,*) ' Too many adsorption species specified',
     &           ' in this run'
             write(iunit2,*) ' Maximum allowed is:', mads
             stop
         endif
c
         naads(nads)=name
         ncpad(nads)=ncp
         do n=1,ncp                   ! loops over components in the ads.species
           do j=1,npri                      
            if(nam(n).eq.napri(j)) then
              stqd(nads,n) = s(n)
              icpad(nads,n) = j
            endif 
           enddo
         enddo
         zd(nads)=zz
         if(ireg.eq.0) then ! if no regression coefficients are given, calculate them
             call evaluate_logk(alogk,ntmp,ntemp,tempc,name)
             call regression_logK(alogk,tempa,w,indx,nfit,ntemp,ntmp,vb)
             call evaluate_regression(vb,ntmp,ntemp,tempc,name)
             do n=1,5  
                akcoead(nads,n)=vb(n)
             enddo
         else                 ! if regression coefficients  are input
             do n=1,5  
                akcoead(nads,n)=coef(n)
             enddo
         endif
         icheck(nads)=1
c
        else    ! if we listed all surface complexes in chemical.inp file 
c       ----------------------------------------------------------------------------
c       
c       checks if ads.sp. is in the list of ads.sp. specified for this run; if yes store data
         do m = 1,nads               ! loops over specified ads.sp. in the system
          if(name.eq.naads(m)) then  ! found the ads.species
            ncpad(m)=ncp
            do 110 n=1,ncp           ! loops over components in the ads.species
              do j=1,npri                      
                if(nam(n).eq.napri(j)) then
                  stqd(m,n) = s(n)
                  icpad(m,n) = j
                  go to 110          ! resumes looping over components
                endif
              enddo
              write(32,*) '  Ads.species:', name      ! gets to this point only if no match is found - read next mineral
              write(32,*) '  One or more components missing - stop'  
              write(iunit2,*) '  Ads.species:', name  ! gets to this point only if no match is found - read next mineral
              write(iunit2,*) '  One or more components missing - stop'  
              stop
  110       continue   ! loop 
            zd(m)=zz
            if(ireg.eq.0) then ! if no regression coefficients are given, calculate them
             call evaluate_logk(alogk,ntmp,ntemp,tempc,name)
             call regression_logK(alogk,tempa,w,indx,nfit,ntemp,ntmp,vb)
             call evaluate_regression(vb,ntmp,ntemp,tempc,name)
             do n=1,5  
                akcoead(m,n)=vb(n)
             enddo
            else                 ! if regression coefficients  are input
             do n=1,5  
                akcoead(m,n)=coef(n)
             enddo
            endif
            icheck(m)=1
          endif
         enddo              
        endif
c
      goto 300
c--- end of loop through adsorbed species
c
  120 continue 
      do m=1,nads      
        if(icheck(m).ne.1) then
          write(32,*) ' Ads.species:', naads(m)
          write(32,*) ' Not found in database - stop'
          write(iunit2,*) ' Ads.species:', naads(m)
          write(iunit2,*) ' Not found in database - stop'
          stop
        endif
      enddo
c
      write(32,*)
      write(32,*) ' --- List of all possible surface complexes --- '
      do ix=1,ns_tot
         write(32,'(10x,a20)') namexp(ix)
      end do
c
      if (iauto_read.eq.0) then
       write(32,*)
       write(32,*) ' --- Surface complexes selected -------------- '
       do ix=1,nads
          write(32,'(10x,a20)') naads(ix)
       end do
       write(32,*) ' ---------------------------------------------- '
       write(32,*)
      endif
c
c
c--------done reading the database----------------------------
c     close data file
  600 close(iunit5)
      write(32,*) '   --> Finished reading thermodynamic database'
      write(iunit2,*) '   --> Finished reading thermodynamic database'
c
      if(iflgstop.eq.1) stop
      return 
c
6001  write(32,*) 'Error reading temperature data in database: stop'
      stop
6002  write(32,*) 'Error reading primary species in database: stop'
      write(32,*) ' Input record starting with:      ', inprec(1:40)
      write(32,*) ' Current or previous record name: ', name
      stop
6003  write(32,*) 'Error reading secondary species: stop'
      write(32,*) ' Input record starting with:      ', inprec(1:40)
      write(32,*) ' Current or previous record name: ', name
      stop
6004  write(32,*) 'Error reading minerals: stop'
      write(32,*) ' Input record starting with:      ', inprec(1:40)
      write(32,*) ' Current or previous record name: ', name
      stop
6005  write(32,*) 'Error reading gases: stop'
      write(32,*) ' Input record starting with:      ', inprec(1:40)
      write(32,*) ' Current or previous record name: ', name
      stop
6006  write(32,*) 'Error reading adsorbed species: stop'
      write(32,*) ' Input record starting with:      ', inprec(1:40)
      write(32,*) ' Current or previous record name: ', name
c
      stop
500   write(32,*) 'Error in opening database file: ',thermo_in, 'stop'
      stop
c
      end
c
c
c
       subroutine name_conv(label)
c
c *****************************************************************************
c     Subroutine ptz: This subroutine regulates the chemical names
c
c *****************************************************************************
c
      Character*20 label
c
      i=1
 10   if (label(i:i).eq.' ') then
         do j=1,19
          label(j:j)=label(j+1:j+1)
        enddo
        goto 10
      endif
c
      do i=1,20
c
      if (label(i:i).eq.'A') label(i:i)='a'
      if (label(i:i).eq.'B') label(i:i)='b'
      if (label(i:i).eq.'C') label(i:i)='c'
      if (label(i:i).eq.'D') label(i:i)='d'
      if (label(i:i).eq.'E') label(i:i)='e'
      if (label(i:i).eq.'F') label(i:i)='f'
      if (label(i:i).eq.'G') label(i:i)='g'
      if (label(i:i).eq.'H') label(i:i)='h'
      if (label(i:i).eq.'I') label(i:i)='i'
      if (label(i:i).eq.'J') label(i:i)='j'
      if (label(i:i).eq.'K') label(i:i)='k'
      if (label(i:i).eq.'L') label(i:i)='l'
      if (label(i:i).eq.'M') label(i:i)='m'
      if (label(i:i).eq.'N') label(i:i)='n'
      if (label(i:i).eq.'O') label(i:i)='o'
      if (label(i:i).eq.'P') label(i:i)='p'
      if (label(i:i).eq.'Q') label(i:i)='q'
      if (label(i:i).eq.'R') label(i:i)='r'
      if (label(i:i).eq.'S') label(i:i)='s'
      if (label(i:i).eq.'T') label(i:i)='t'
      if (label(i:i).eq.'U') label(i:i)='u'
      if (label(i:i).eq.'V') label(i:i)='v'
      if (label(i:i).eq.'W') label(i:i)='w'
      if (label(i:i).eq.'X') label(i:i)='x'
      if (label(i:i).eq.'Y') label(i:i)='y'
      if (label(i:i).eq.'Z') label(i:i)='z'
      enddo
c
      do i=2,15
      if (label(i:i+5).eq.'++++++') then
            label(i:i+1)='+6'
            label(i+2:i+5)='    '
      endif
      if (label(i:i+4).eq.'+++++') then
            label(i:i+1)='+5'
            label(i+2:i+4)='   '
      endif
      if (label(i:i+3).eq.'++++') then
            label(i:i+1)='+4'
            label(i+2:i+3)='  '
      endif
      if (label(i:i+2).eq.'+++') then
            label(i:i+1)='+3'
            label(i+2:i+2)=' '
      endif
      if (label(i:i+1).eq.'++') label(i:i+1)='+2'
c
      if (label(i:i+5).eq.'------') then
          label(i:i+1)='-6'
            label(i+2:i+5)='    '
      endif
      if (label(i:i+4).eq.'-----') then
          label(i:i+1)='-5'
            label(i+2:i+4)='   '
      endif
      if (label(i:i+3).eq.'----') then
          label(i:i+1)='-4'
            label(i+2:i+3)='  '
      endif
      if (label(i:i+2).eq.'---') then
          label(i:i+1)='-3'
            label(i+2:i+2)=' '
      endif
      if (label(i:i+1).eq.'--') label(i:i+1)='-2'
      enddo
c
      return
      end
c
c
c
c *****************************************************************************
c
       subroutine evaluate_logk(alogk,ntmp,ntemp,tempc,name)
c
c *****************************************************************************
c     Subroutine evaluate_logk: This subroutine checks logk values, if there is any
c                  :blank (500.00) put values based on interpolation/extrapolation
c
c *****************************************************************************
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      common/evaluate_logk_regression/islop0
      double precision alogk(ntmp),aalogk(ntmp),tk(ntmp),
     +     tempa(ntmp,ntmp),tempc(ntmp),slop(ntmp)
      integer*8 indx(ntmp),ipoint(ntmp)
      character*20 name
      character*1 number
      isign=0
      do i=1,ntemp
      ipoint(i)=0
      if (dabs(alogk(i)-500.0d0).lt.1.0d-20) then
      ipoint(i)=1
      isign=isign+1
      endif
      enddo
c
      if (isign.eq.0) then
       islop=0
       istrenge=1
       slop0=0.0d0
       do i=2,ntemp-1    ! calculate slop
        slop(i)=(alogk(i)-alogk(i-1))/(tempc(i)-tempc(i-1))
        slop(i+1)=(alogk(i+1)-alogk(i))/(tempc(i+1)-tempc(i))
        if(slop(i)*slop(i+1).lt.-1.0d-20) islop=islop+1
       enddo
c
       do i=2,ntemp
        do j=2,ntemp
          if (j.eq.i)slop0=slop0+dabs(slop(j))
        enddo
        slop0=slop0/(ntemp-2)
        if (slop0.lt.1.0d-20) return
        if (istrenge.lt.slop(i)/slop0/2)
     +     istrenge=slop(i)/slop0/2
       enddo
c
       islop0=islop
       if (istrenge.gt.1.or.islop.gt.1) then
c        write(*,*) '****** Warning! LogK values scatter: ',name
        write(32,'(a36)') '****** Warning! LogK values acatter!'
c        write(*,'(a8,8f9.3)') 'Temp:   ',(tempc(i),i=1,ntemp)
c        write(*,'(a8,8f9.3)') 'LogK:   ',(alogk(i),i=1,ntemp)
c        write(*,'(a18)') 'Ignore this (y/n)?'
c        read(*,*) number
c         if (number.eq.'Y'.or.number.eq.'y') return
c         if (number.eq.'N'.or.number.eq.'n') then
c          write(*,*) ' Stop, check database: data0.ypf!'
c          stop
c         endif
       endif
      return
      endif
c
      if(isign.eq.ntemp) then
        write(32,*) ' Stop! Logk values are not available'
        return
      endif
c
      if(isign.eq.(ntemp-1)) then
       write(32,*)
       write(32,'(a38)') '****** Caution! Constant logk is used:'
       do i=1,ntemp
        if (dabs(alogk(i)-500.0d0).gt.1.0d-20) then
         do j=1,ntemp
          alogk(j)=alogk(i)
         enddo
         write(32,'(5x,8f10.4)')
     +    (alogk(j),j=1,ntemp)
         write(32,*)
         return
        endif
       enddo
      endif
c
      write(32,'(a60)')
     +'****** Caution! logk values short,interpolated/extrapolated:'
c
       ii=0
       do i=1,ntemp
        if (ipoint(i).ne.1) then
           ii=ii+1
           tk(ii)=tempc(i)+273.15d0
           aalogk(ii)=alogk(i)
        endif
       enddo
c
      do i=1,ii
        tkk=tk(i)
        tempa(i,1)=dlog(tkk)
        tempa(i,2)=1.0d0
        tempa(i,3)=1.0d0/tkk
        tempa(i,4)=tkk
        tempa(i,5)=1.0d0/tkk**2
        tempa(i,6)=tkk**2
        tempa(i,7)=1.0d0/tkk**3
        tempa(i,8)=tkk**3
      enddo
c
      call ludcmp(tempa,ii,ntmp,indx,dd)
c
      call lubksb(tempa,ii,ntmp,indx,aalogk)
c
      do i=ii+1,ntemp
       aalogk(i)=0.0
      enddo
c
      do i=1,ntemp
       tkk=tempc(i)+273.15d0
       aaalogk=
     + aalogk(1)*dlog(tkk)+
     + aalogk(2)+
     + aalogk(3)*1.0d0/tkk+
     + aalogk(4)*tkk+
     + aalogk(5)*1.0d0/tkk**2+
     + aalogk(6)*tkk**2+
     + aalogk(7)*1.0d0/tkk**3+
     + aalogk(8)*tkk**3
       alogk(i)=aaalogk
      enddo
      write(32,'(5x,8f10.4)')
     +    (alogk(i),i=1,ntemp)
c
       islop=0
       istrenge=1
       slop0=0.0d0
       do i=2,ntemp-1    ! calculate slop
        slop(i)=(alogk(i)-alogk(i-1))/(tempc(i)-tempc(i-1))
        slop(i+1)=(alogk(i+1)-alogk(i))/(tempc(i+1)-tempc(i))
        if(slop(i)*slop(i+1).lt.-1.0d-20) islop=islop+1
       enddo

       do i=2,ntemp
        do j=2,ntemp
          if (j.eq.i)slop0=slop0+dabs(slop(j))
        enddo
        slop0=slop0/(ntemp-2)
        if (slop0.lt.1.0d-20) return
        if (istrenge.lt.slop(i)/slop0/2)
     +     istrenge=slop(i)/slop0/2
       enddo
c
       islop0=islop
       if (istrenge.gt.1.or.islop.gt.1) then
c        write(*,*) '****** Warning! LogK values scatter: ',name
        write(32,'(a36)') '****** Warning! LogK values acatter!'
c        write(*,'(a8,8f9.3)') 'Temp:   ',(tempc(i),i=1,ntemp)
c        write(*,'(a8,8f9.3)') 'LogK:   ',(alogk(i),i=1,ntemp)
c        write(*,'(a18)') 'Ignore this (y/n)?'
c        read(*,*) number
c         if (number.eq.'Y'.or.number.eq.'y') return
c         if (number.eq.'N'.or.number.eq.'n') then
c          write(*,*) ' Stop, check database: data0.ypf!'
c          stop
c         endif
       endif
c
      return
c
      end
c
c
c
       subroutine evaluate_regression(vb,ntmp,ntemp,tempc,name)
c
c *****************************************************************************
c     Subroutine evaluate_regression: This subroutine checks regression equation,
c                  :if there is any oscillations between two adjacent temperature
c                  :points send a warning
c
c *****************************************************************************
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      character*20 name
      character*1 number
      common/evaluate_logk_regression/islop0
      double precision vb(ntmp),aalogk(100),
     +          tempc(ntmp),slop(ntmp),aaa(5000)
c
c=====checking second order derivative term, too rigorous!===========
c
      goto 40
      tkk=tempc(1)+273.15d0
      ii=0
 30   tkk=tkk+0.1d0
      ii=ii+1
      if (tkk.gt.(tempc(ntemp)+273.15d0)) then
      write(32,*) name
      write(32,'(10(e8.2,1x))') (aaa(i),i=1,ii-1)
      return
      endif
       aaa(ii)=
     +      -vb(1)/tkk**2+
     + 2.0d0*vb(4)/tkk**3+         ! alogk(2) term =0.0d0
     + 6.0d0*vb(5)/tkk**4          ! alogk(3) term =0.0d0
      if (ii.lt.2) goto 30
      if (aaa(ii)*aaa(ii-1).lt.-1.0d-50) then
       write(32,'(a39,f10.3,a3)')
     +  'Logk regression equation oscillates at:',tkk-273.15d0, '(C)'
c       write(*,'(a20,a39,f10.3,a3)')
c     +  name,'Logk regression equation oscillates at:',
c     +  tkk-273.15d0, '(C)'
c       write(*,'(a18)') 'Ignore this (y/n)?'
c       read(*,*) number
c         if (number.eq.'Y'.or.number.eq.'y') goto 30
c         if (number.eq.'N'.or.number.eq.'n') then
c          write(*,*) ' Stop, check database: data0.ypf!'
c          stop
c         endif
      endif
c
      goto 30
c
c====check slop (first order derivative term) between original data points======
c
 40   do i=1,ntemp
        tkk=tempc(i)+273.15d0
        slop(i)=vb(1)/tkk
     +         +vb(3)                   ! alogk(2) term =0.0d0
     +         -vb(4)/tkk**2
     +         -2.0d0*vb(5)/tkk**3
      enddo
c
      islop=0
      do i=2,ntemp
       if (slop(i)*slop(i-1).lt.-1.0d-20) islop=islop+1
      enddo
c
c      write(2,'(a8,8f10.3)') 'Slop = ',(slop(i),i=1,ntemp)
c
      if (islop.gt.1.and.islop0.eq.1) then
       write(32,*) 'islop0= ',islop0
       write(32,'(a47)')
     +'****** Oscillations in LogK regression equation'
c       write(*,'(a49,a20)')
c     +'****** Oscillations in LogK regression equation: ',name
c
      do i=1,ntemp
       aalogk(i)=fak2(vb,tempc(i))
      enddo
      write(32,'(5x,8f10.4)')
     +    (aalogk(i),i=1,ntemp)
      write(32,*)
c
c       write(*,'(a18)') 'Ignore this (y/n)?'
c       read(*,*) number
c         if (number.eq.'Y'.or.number.eq.'y') return
c         if (number.eq.'N'.or.number.eq.'n') then
c          write(*,*) ' Stop, check database: data0.ypf!'
c          stop
c         endif
      endif
c
c====check slop (first order derivative term) between intervals with 1/100 of original data points======
c
      do i=2,ntemp
        islop=0
        do j=1,100
        tkk1=tempc(i-1)+dble(j-1)*(tempc(i)-tempc(i-1))*0.01d0
        tkk1=tkk1+273.15d0

        tkk2=tempc(i-1)+dble(j)*(tempc(i)-tempc(i-1))*0.01d0
        tkk2=tkk2+273.15d0
c
        slop1=vb(1)/tkk1
     +       +vb(3)
     +       -vb(4)/tkk1**2        ! alogk(2) term =0.0d0
     +       -2.0d0*vb(5)/tkk1**3
c
        slop2=vb(1)/tkk2
     +       +vb(3)
     +       -vb(4)/tkk2**2        ! alogk(2) term =0.0d0
     +       -2.0d0*vb(5)/tkk2**3
c
        if(slop1*slop2.lt.-1.0d-50) islop=islop+1
        enddo
c
      if (islop.gt.1) then !.and.islop0.eq.1
       write(32,'(a47)')
     +'****** Oscillations in LogK regression equation'
c       write(*,'(a49,a20)')
c     +'****** Oscillations in LogK regression equation: ',name
c
      do j=1,100
        tkk=tempc(i-1)+dble(j-1)*(tempc(i)-tempc(i-1))*0.01d0
       aalogk(i)=fak2(vb,tkk)
c       write(2,'(3f20.4)')
c     +      tempc(i),aalogk(i),alogk(i)-aalogk(i)
      enddo
c
       write(32,'(a14,2f10.3)') 'Between temp:', tempc(i-1),tempc(i)
       write(32, '(8f10.3)') (aalogk(j),j=1,100)
c
c       write(*,'(a18)') 'Ignore this (y/n)?'
c       read(*,*) number
c         if (number.eq.'Y'.or.number.eq.'y') return
c         if (number.eq.'N'.or.number.eq.'n') then
c          write(*,*) ' Stop, check database: data0.ypf!'
c          stop
c         endif
         endif
c
      enddo
c
      return
      end
c
c
c
       subroutine pre_reg(tempc,ntemp,nfit,ntmp,indx,tempa,w)
c
c *****************************************************************************
c     Subroutine pre_reg: This subroutine prepares the matrix for logK regression
c                  :
c     Arguments    :tempc(ntemp): temperature points; ntemp # of temperature points;
c                  :w(nfit,nfit): regression matrix; nfit: # of fit parameters
c *****************************************************************************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      double precision tempa(ntmp,ntmp),w(ntmp,ntmp),tempc(ntmp)
      integer*8 indx(ntmp)
      do i=1,ntemp
        tkk=tempc(i)+273.15d0
        tempa(1,i)=dlog(tkk)
        tempa(2,i)=1.0d0
        tempa(3,i)=tkk
        tempa(4,i)=1.0d0/tkk
        tempa(5,i)=1.0d0/tkk**2
      enddo
c
      do i=1,nfit
        do j=i,nfit
          w(i,j)=0.0d0
            do k=1,ntemp
              w(i,j)=w(i,j)+tempa(i,k)*tempa(j,k)
            enddo
          if(i.ne.j) w(j,i)=w(i,j)
        enddo
      enddo
c
      call ludcmp(w,nfit,ntmp,indx,dd)
c
      return
      end
c
c
c
       subroutine regression_logK(alogk,tempa,w,indx,nfit,ntemp,ntmp,vb)
c
c *****************************************************************************
c     Subroutine regression_logK: This subroutine calculates the regression
c                  :coefficients for logK
c
c     Arguments    :tempc(ntemp): temperature points; ntemp # of temperature points;
c                  :w(nfit,nfit): regression matrix; nfit: # of fit parameters
c                  :vb(nfit): the coefficients returned.
c *****************************************************************************
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      double precision tempa(ntmp,ntmp),w(ntmp,ntmp),vb(ntmp),
     +   alogk(ntmp)
      integer*8 indx(ntmp)
      do j=1,nfit
        vb(j)=0.0d0
        do i=1,ntemp
          vb(j)=vb(j)+alogk(i)*tempa(j,i)
        enddo
      enddo
c
      call lubksb(w,nfit,ntmp,indx,vb)
c
      return
      end
c
c
c
      subroutine echotherm
c
c  Routine to echo all data read in the thermodynamic database
c  See routine readtherm for variable definitions
c  By: N.S.  LBL 5/98
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
c---------------------------------------------------------------------------

      double precision coef(5)
      double precision stoic(mpri)
      data iunit2/42/
c
      temp=25.d0   ! for logK calculation and printout only
c
      write(iunit2,*) '  Components     a0    charge'
      do i=1,npri
        write(iunit2,"(2x,a10,2x,f8.3,2x,f5.1)") napri(i),a0(i),z(i)
      enddo
c.....Write out temperature interpolation coefficients
      write(iunit2,
     & "(//' Log K interpolation coefficients (a,b,c,d,e)')")
      write(iunit2,
     & "(' valid temperature range (deg.C.): ',f5.1,' to ',f5.1)")
     &  tmpmin,tmpmax
      write(iunit2,
     &  "(' ',15x,'        a*ln(TK)       b          c*TK'
     &  ,'     d*(TK)**-1  e*(TK)**-2')" )
      do i = 1, naqx
        write(iunit2,1000) naaqx(i),(akcoes(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoes(i,n)
        enddo
        aks(i)=fak2(coef,temp)
      enddo
      do i = 1, nmin
        write(iunit2,1000) namin(i),(akcoem(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoem(i,n)
        enddo
        akm(i)=fak2(coef,temp)
      enddo
      do i = 1, ngas
        write(iunit2,1000) nagas(i),(akcoeg(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoeg(i,n)
        enddo
        akg(i)=fak2(coef,temp)
      enddo
      do i = 1, nads
        write(iunit2,1000) naads(i),(akcoead(i,j),j=1,5)
        do n=1,5
         coef(n)=akcoead(i,n)
        enddo
        akd(i)=fak2(coef,temp)
      enddo
 1000 format(' ',a20,5(1pe12.4))
c
c.....Write stoichiometries and logK's at T=temp (defined above)
      if(naqx.gt.0) then
          write(iunit2,"(/' Derived Species Reactions')")
          write(iunit2,590) temp, (napri(j),j=1,npri)
          do i = 1, naqx
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
          write(iunit2,595) temp, (napri(j),j=1,npri)
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
          write(iunit2,596) temp,(napri(j),j=1,npri)
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
          write(iunit2,597) temp,(napri(j),j=1,npri)
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
     &   2x,30a7)
  595 format(2x,'minerals',4x,'m.vol(L/mol)',2x,
     &   'logK(',f3.0,' C)',2x,30a7)
  596 format(2x,'gases   ',18x,'logK(',f4.0,'C)',2x,30a7)
  597 format(2x,'surf.cmplex',15x,'logK(',f4.0,'C)',2x,30a7)
  600 format(' ',a12,f8.3,1x,f5.1,3x,f10.3,30(f7.2))
  605 format(' ',a12,2x,f10.4,2x,f10.3,1x,30(f7.2))
  606 format(' ',a12,14x,f10.3,1x,30(f7.2))
  610 format(' ',a12,7x,f5.1,3x,f10.3,30(f7.2))
c
      return
      end
c
c
c
        double precision function fak2(coef,temp)
c
        implicit double precision (a-h,o-z)
        implicit integer*8 (i-n)
        include 'flowpar_v2.inc'
        include 'chempar_v2.inc'
        include 'common_v2.inc'
        double precision coef(5)
c
        degc=temp
        if(degc.gt.tmpmax) then
c           write(*,*) 'Temperature is outside of log K fit range'
c           write(*,*) 'Temperature:',temp,' reset down to:',tmpmax
            degc=tmpmax                 ! add fit constraints
        elseif(degc.lt.tmpmin) then
            degc=tmpmin                 ! add fit constraints
        endif
        tk=degc+273.15d0                ! use tc instead of tc2
        fak2=coef(1)*dlog(tk)+coef(2)+coef(3)*tk+coef(4)/tk+
     +     coef(5)/(tk**2)
        return
        end
c
c
       subroutine char_number(cc,rnumb)
c
c *****************************************************************************
c     Subroutine char_number: This subroutine convert the character into number
c                  
c *****************************************************************************
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      character*1 cc
      rnumb=0.d0
      if (cc.eq.'1') rnumb=1.d0
      if (cc.eq.'2') rnumb=2.d0
      if (cc.eq.'3') rnumb=3.d0
      if (cc.eq.'4') rnumb=4.d0
      if (cc.eq.'5') rnumb=5.d0
      if (cc.eq.'6') rnumb=6.d0
      if (cc.eq.'7') rnumb=7.d0
      if (cc.eq.'8') rnumb=8.d0
      if (cc.eq.'9') rnumb=9.d0
c
      return
      end
c
