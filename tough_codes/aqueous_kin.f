C
C
      SUBROUTINE READ_AQKIN
C
c******************** Read data for aqueous kinetics
c
C
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
C
      common/aqkin1/nrx
      common/aqkin21/i_mod(mrx)
      common/aqkin22/n_mech(mrx)     ! Number of mechanisms
      common/aqkin3/rkaq(mrx,mechaq)
      common/aqkin4/rkaq_abc(mrx,mechaq,3)   ! Parameters for calculating rkaq
      common/aqkin5/ncp_rx(mrx),s_rx(mrx,maqsp),
     +              coef_rx(mrx,mechaq),nam_rx(mrx,maqsp)
      common/aqkin13/ncp_rx1(mrx,mechaq),s_rx1(mrx,mechaq,maqsp),
     +               ia1(mrx,mechaq,maqsp),nam_rx1(mrx,mechaq,maqsp)   ! Product terms
      common/aqkin14/ncp_rx2(mrx,mechaq),s_rx2(mrx,mechaq,maqsp),
     +               ia2(mrx,mechaq,maqsp),nam_rx2(mrx,mechaq,maqsp)   ! Monod terms
      common/aqkin15/ncp_rx3(mrx,mechaq),s_rx3(mrx,mechaq,maqsp),
     +               ia3(mrx,mechaq,maqsp),nam_rx3(mrx,mechaq,maqsp)   ! Inhibition terms
      character*20 nam_rx,nam_rx1,nam_rx2,nam_rx3
      character*100 dummy
      character*80 label
      character*200 inprec 
C
      read(41,*,err=9022) irx
      write (42,'(i3)')   irx
      read(41,*,err=9022) ncp_rx(nrx),(s_rx(nrx,n),
     +             nam_rx(nrx,n),n=1,ncp_rx(nrx))
!
      do n=1,ncp_rx(nrx)
         label(1:20)=nam_rx(nrx,n)
         call name_conv(label)
         nam_rx(nrx,n)=label(1:20)
      end do
!
      write(42,'(i5,10(f8.4,2x,a13))') ncp_rx(nrx),
     +     (s_rx(nrx,n),nam_rx(nrx,n),n=1,ncp_rx(nrx))
      read  (41,*,err=9022)   i_mod(nrx),n_mech(nrx)  ! reaction model,No. mechanisms
      write (42,'(2i5)')   i_mod(nrx),n_mech(nrx)
!
!.....f number of mechanisms exceeds max
      if(n_mech(nrx).gt.mechaq) then
         write (32,"('maximum number of aqueous kinetic mechanisms is',
     &     ' exceeded'/' max allowed =',i5)") mechaq
         write (42,"('maximum number of aqueous kinetic mechanisms is',
     &     ' exceeded'/' max allowed =',i5)") mechaq
         stop
      endif
!
      do im=1,n_mech(nrx)
         read  (41,*,err=9022)   rkaq(nrx,im)         ! forward rate constant
         write (42,'(4x,e12.4)')    rkaq(nrx,im)
!
!........If rkaq is negative, the rate constant will be calculated from variables such as T
!
         if (rkaq(nrx,im) .eq. -1.0d0)   then
            read  (41,*)              rk_a, rk_b      !, rk_c
            write (42,'(2x,e12.4)')   rk_a, rk_b      !, rk_c
!
            rkaq_abc(nrx,im,1) = rk_a                 ! Store k25  
            rkaq_abc(nrx,im,2) = rk_b                 ! Store Ea 
         end if           
!
!...........................................................................................
!
         read(41,"(a200)",err=9022) inprec
         read(inprec,*,err=9022) ndum  

!........If number of species in product term exceeds max
         if(ndum.gt.maqsp) then
          write (32,"('maximum number of aqueous kinetic mechanisms is',
     &     ' exceeded'/' max allowed =',i5)") maqsp
          write (42,"('maximum number of aqueous kinetic mechanisms is',
     &     ' exceeded'/' max allowed =',i5)") maqsp
          stop
         endif
!
         read(inprec,*,err=9022) ncp_rx1(nrx,im),(nam_rx1(nrx,im,n),  ! product term
     +      ia1(nrx,im,n),s_rx1(nrx,im,n),n=1,ncp_rx1(nrx,im))
!
         do n=1,ncp_rx1(nrx,im)
            label(1:20)=nam_rx1(nrx,im,n)
            call name_conv(label)
            nam_rx1(nrx,im,n)=label(1:20)
         end do
!
         write(42,'(2x,i5,10(15x,a13,i5,f6.2))') ncp_rx1(nrx,im),
     +                      (nam_rx1(nrx,im,n),ia1(nrx,im,n),
     +                   s_rx1(nrx,im,n),n=1,ncp_rx1(nrx,im))
!
         read(41,*,err=9022) ncp_rx2(nrx,im),(nam_rx2(nrx,im,n),  ! Monod term
     +      ia2(nrx,im,n),s_rx2(nrx,im,n),n=1,ncp_rx2(nrx,im))
!
         do n=1,ncp_rx2(nrx,im)
            label(1:20)=nam_rx2(nrx,im,n)
            call name_conv(label)
            nam_rx2(nrx,im,n)=label(1:20)
         end do
!
         write(42,'(2x,i5,10(15x,a13,i5,f6.2))') ncp_rx2(nrx,im),
     +                      (nam_rx2(nrx,im,n),ia2(nrx,im,n),
     +                   s_rx2(nrx,im,n),n=1,ncp_rx2(nrx,im))
!
         read(41,*,err=9022) ncp_rx3(nrx,im),(nam_rx3(nrx,im,n),  ! Inhibition term
     +      ia3(nrx,im,n),s_rx3(nrx,im,n),n=1,ncp_rx3(nrx,im))
!
         do n=1,ncp_rx3(nrx,im)
            label(1:20)=nam_rx3(nrx,im,n)
            call name_conv(label)
            nam_rx3(nrx,im,n)=label(1:20)
         end do
!
         write(42,'(2x,i5,10(15x,a13,i5,f6.2))') ncp_rx3(nrx,im),
     +                      (nam_rx3(nrx,im,n),ia3(nrx,im,n),
     +                   s_rx3(nrx,im,n),n=1,ncp_rx3(nrx,im))
!
         if (i_mod(nrx).eq.2)     then      ! i_mod=2 for reversible reaction
            read(41,'(a100)') dummy         ! skips logK record
            write (42,'(a100)') dummy
            read(41,*,err=9022)  (coef_rx(nrx,n),n=1,5)
            write(42,'(5e17.8)')  (coef_rx(nrx,n),n=1,5)
         end if
!
      end do
!
      return
!
9022  write (42,*) 'error reading aqueous kinetics of the system'
      write (32,*) 'error reading aqueous kinetics of the system'
      stop
!
      end
!
!
!
      SUBROUTINE AQKIN_STOICHOMETRY
!
!************ Obtain stoichometry coefficients for aqueous redox reactions ******
!
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
!
      common/aqkin2/ntrx              ! total number of redox pair
      common/aqkin22/n_mech(mrx)     ! Number of mechanisms
      common/aqkin5/ncp_rx(mrx),s_rx(mrx,maqsp),
     +              coef_rx(mrx,mechaq),nam_rx(mrx,maqsp)
      common/aqkin6/stqrx(mrx,maqsp),icprx(mrx,maqsp)
      common/aqkin13/ncp_rx1(mrx,mechaq),s_rx1(mrx,mechaq,maqsp),
     +               ia1(mrx,mechaq,maqsp),nam_rx1(mrx,mechaq,maqsp)   ! Product terms
      common/aqkin14/ncp_rx2(mrx,mechaq),s_rx2(mrx,mechaq,maqsp),
     +               ia2(mrx,mechaq,maqsp),nam_rx2(mrx,mechaq,maqsp)   ! Monod terms
      common/aqkin15/ncp_rx3(mrx,mechaq),s_rx3(mrx,mechaq,maqsp),
     +               ia3(mrx,mechaq,maqsp),nam_rx3(mrx,mechaq,maqsp)   ! Inhibition terms
      common/aqkin33/icprx1(mrx,mechaq,maqsp),icprx2(mrx,mechaq,maqsp),
     +               icprx3(mrx,mechaq,maqsp)                      ! Species index
      character*20 nam_rx,nam_rx1,nam_rx2,nam_rx3
!
      naqt=npri+naqx
!
!-------------------------- Obtain stoichometric coefficients for mass balance term
!
      do m = 1,ntrx
         ncp=ncp_rx(m)
         do n=1,ncp
           NotFound  = 1
           do j=1,npri
             if(nam_rx(m,n).eq.napri(j)) then
               stqrx(m,n) = s_rx(m,n)       ! number of stoichio components
               icprx(m,n)=j                 ! species index
               NotFound = 0
             end if
           end do
           if (NotFound .eq. 1)   then
            write (32, 894)   nam_rx(m,n)
              stop
         end if
         end do
      end do
!
!-------------Obtain stoichometric coefficients for power terms of the rate expression
!
      do m = 1,ntrx          ! loop over reactions
       do im = 1,n_mech(m)   !      over mechanisms
!-------------------------------------------------------Product term
         ncp=ncp_rx1(m,im)
         do n=1,ncp
           NotFound  = 1
           do j=1,naqt
             if(nam_rx1(m,im,n).eq.naaqt(j)) then
               icprx1(m,im,n)=j                 ! Species index
               NotFound = 0
             end if
           end do
           if (NotFound .eq. 1)   then
            write (32, 895)   nam_rx1(m,im,n)
              stop
         end if
!
           if (ia1(m,im,n) .eq. 3)  then  ! Have to be in the list of primary species
              IPrimary = 0
                do ip=1,npri
                 if (nam_rx1(m,im,n) .eq. napri(ip))   then
                    IPrimary = 1
                 end if
              end do
              if (IPrimary .eq. 0)   then
                 write (32, 896)   nam_rx1(m,im,n)
                 stop
            end if
           end if
!
         end do
!-------------------------------------------------------Monod term
         ncp=ncp_rx2(m,im)
         do n=1,ncp
           NotFound  = 1
           do j=1,naqt
             if(nam_rx2(m,im,n).eq.naaqt(j)) then
               icprx2(m,im,n)=j                 ! Species index
               NotFound = 0
             end if
           end do
           if (NotFound .eq. 1)   then
            write (32, 895)   nam_rx2(m,im,n)
              stop
         end if
!
           if (ia2(m,im,n) .eq. 3)  then  ! Have to be in the list of primary species
              IPrimary = 0
                do ip=1,npri
                 if (nam_rx2(m,im,n) .eq. napri(ip))   then
                    IPrimary = 1
                 end if
              end do
              if (IPrimary .eq. 0)   then
                 write (32, 896)   nam_rx2(m,im,n)
                 stop
            end if
           end if
!
         end do
!-------------------------------------------------------Inhibition term
         ncp=ncp_rx3(m,im)
         do n=1,ncp
           NotFound  = 1
           do j=1,naqt
             if(nam_rx3(m,im,n).eq.naaqt(j)) then
               icprx3(m,im,n)=j                 ! Species index
               NotFound = 0
             end if
           end do
           if (NotFound .eq. 1)   then
            write (32, 895)   nam_rx3(m,im,n)
              stop
         end if
!
           if (ia3(m,im,n) .eq. 3)  then  ! Have to be in the list of primary species
              IPrimary = 0

                do ip=1,npri
                 if (nam_rx3(m,im,n) .eq. napri(ip))   then
                    IPrimary = 1
                 end if
              end do
              if (IPrimary .eq. 0)   then
                 write (32, 896)   nam_rx3(m,im,n)
                 stop
            end if
           end if
!
         end do
!---------------------------------------------------------------------
       end do
      end do
!
894   format (//2X,'Error: Species for mass balance terms',
     &            ' in AQUEOUS KINETICS block: ',a15,
     &    /6x, ' was not found in the list of primary species')
!
895   format (//2X,'Error: Species for Product/Monod/Inhibition terms',
     &             ' at AQUEOUS KINETICS block: ',a15,
     &    /6x, ' was not found in the list of all aqueous species')
!
896   format (//2X,'Error: When using the total concentration',
     &            ' in AQUEOUS KINETICS block: the species, ',a15,
     &    /6x, ', must be in the list of primary species')
C
      return
      end
C
C
C
c----------------------------------------------------------------------
c
      SUBROUTINE CR_CP_RX
c
c******** This routine calculates aqueous kinetic reaction rates ***********
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
C
      common/aqkin2/ntrx              ! total number of redox pair
      common/aqkin21/i_mod(mrx)
      common/aqkin22/n_mech(mrx)     ! Number of mechanisms
      common/aqkin3/rkaq(mrx,mechaq)
      common/aqkin4/rkaq_abc(mrx,mechaq,3)   ! Parameters for calculating rkaq
      common/aqkin5/ncp_rx(mrx),s_rx(mrx,maqsp),
     +              coef_rx(mrx,mechaq),nam_rx(mrx,maqsp)
      common/aqkin13/ncp_rx1(mrx,mechaq),s_rx1(mrx,mechaq,maqsp),
     +               ia1(mrx,mechaq,maqsp),nam_rx1(mrx,mechaq,maqsp)   ! Product terms
      common/aqkin14/ncp_rx2(mrx,mechaq),s_rx2(mrx,mechaq,maqsp),
     +               ia2(mrx,mechaq,maqsp),nam_rx2(mrx,mechaq,maqsp)   ! Monod terms
      common/aqkin15/ncp_rx3(mrx,mechaq),s_rx3(mrx,mechaq,maqsp),
     +               ia3(mrx,mechaq,maqsp),nam_rx3(mrx,mechaq,maqsp)   ! Inhibition terms
      common/aqkin33/icprx1(mrx,mechaq,maqsp),icprx2(mrx,mechaq,maqsp),
     +               icprx3(mrx,mechaq,maqsp)                      ! Species index
      common/aqkin6/stqrx(mrx,maqsp),icprx(mrx,maqsp)
      double precision rkin2rx(mrx)
      common/aqkin11/crx(mpri)
!
      common/satgas2/sg2
c
      character*20 nam_rx,nam_rx1,nam_rx2,nam_rx3
      double precision at(maqt)
c
c--------------------------------Take concentration
c
      do j=1,npri
         ct(j)=cp(j)
      end do
c
      do j=1,naqx
         ct(npri+j)=cs(j)
      end do
c
c--------------------------------Take activity of primary species
c
      do j=1,npri
         at(j)=cp(j)*gamp(j)
      end do
c
      do j=1,naqx
         at(npri+j)=cs(j)*gams(j)
      end do
c
C-------------------------------------------------Initialize crx
      do j=1,npri
        crx(j)=0.0d0     ! crx is moles tied up in aqueous kinetics
      end do
C---------------------------------------------------------------
      do 200 i=1,ntrx          ! Over kinetic reaction
c
         rkin2rx(i)=0.0d0
         do im=1,n_mech(i)     ! Over mechanisms
!
           rdum= rkaq(i,im)
!
!..........If rkaq is negative, the rate constant will be calculated from variables such as T
!
           if (rkaq(i,im) .eq. -1.0d0)   then               
               eadum = ((1.0d0/tk2) - (1.0d0/298.15d0))/(gc1*1.d-3)
               p_k25 = rkaq_abc(i,im,1)
               p_Ea  = rkaq_abc(i,im,2)
               rdum  = p_k25*dexp(-1.0d0*p_Ea*eadum)
           end if 
!
c---------------------------------------Account product term
           ncp1=ncp_rx1(i,im)
           do n=1,ncp1
              j=icprx1(i,im,n)
              if (ia1(i,im,n).eq.1) cc=at(j)     ! Activity
              if (ia1(i,im,n).eq.2) cc=ct(j)     ! Concentration
              if (ia1(i,im,n).eq.3) cc=u2(j)     ! Total concentration
              d1=cc**s_rx1(i,im,n)
              rdum=rdum*d1
           end do
c---------------------------------------Account monod term
           ncp2=ncp_rx2(i,im)
           do n=1,ncp2
              j=icprx2(i,im,n)
              if (ia2(i,im,n).eq.1) cc=at(j)     ! Activity
              if (ia2(i,im,n).eq.2) cc=ct(j)     ! Concentration
              if (ia2(i,im,n).eq.3) cc=u2(j)     ! Total concentration
              d2=cc/(cc+s_rx2(i,im,n))
              rdum=rdum*d2
           end do
c---------------------------------------Account inhibition term
           ncp3=ncp_rx3(i,im)
           do n=1,ncp3
              j=icprx3(i,im,n)
              if (ia3(i,im,n).eq.1) cc=at(j)     ! Activity
              if (ia3(i,im,n).eq.2) cc=ct(j)     ! Concentration
              if (ia3(i,im,n).eq.3) cc=u2(j)     ! Total concentration
              d3=s_rx3(i,im,n)/(cc+s_rx3(i,im,n))
              rdum=rdum*d3
           end do
c--------------------------------------------------------------
           rkin2rx(i)=rkin2rx(i)+rdum
         end do
c
C-------------------------------------------------------------------------
C
         sl2 = 1.0d0 - sg2        ! Liquid saturation
!
         ncp=ncp_rx(i)
         do n=1,ncp
            j=icprx(i,n)
            f_sl2 = 1.d0           ! Stoichmetric correction factor
!                                for adsorbed species in unsaturated zone
            if (nam_rx(i,n) .eq. 'fe(ads)') f_sl2 = sl2
            crx(j)=crx(j)+stqrx(i,n)*rkin2rx(i)*xh2o*f_sl2
         end do
C
200     continue
C
      return
      end
C
C
C
c----------------------------------------------------------------------
c
      SUBROUTINE DCR_DCP_RX
c
c******** This routine calculates aqueous kinetic rate derivatives ***********
c
      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)
      include 'flowpar_v2.inc'
      include 'chempar_v2.inc'
      include 'common_v2.inc'
C
      common/aqkin2/ntrx     ! Total number of redox pair
      common/aqkin6/stqrx(mrx,maqsp),icprx(mrx,maqsp)
      common/aqkin11/crx(mpri)
      common/aqkin12/drx(mpri,mpri)
C
      double precision crxold(mpri)
c
      do j=1,npri
         crxold(j)=crx(j)
      end do
C
      do 200 i=1,npri
         dd=cp(i)*1.0d-07
         cp(i)=cp(i)+dd
c
         call cs_cp
         call cr_cp_rx
c
         do j=1,npri
            drx(j,i)=(crx(j)-crxold(j))/dd
            crx(j)=crxold(j)
         end do
         cp(i)=cp(i)-dd
C
200   continue
c
c
      return
      end
C
C
c----------------------------------------------------------------------
c
c