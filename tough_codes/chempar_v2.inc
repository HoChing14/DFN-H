C------------------------ chempar_v2.inc ----------------------
c... Include file for Chemical Parameters in TOUGHREACT V2.0
C----------------------------------------------------------------------
      parameter (mpri=35,maqx=200,mmin=45,mgas=8,mexc=12,mads=60,
     +   MXsites = 4, msurf=5, mpads=10,mbas=1,ntmp=8,
     +   mechm=5,   !number of mineral kinetic rate mechanisms
     +   mechsp=8,  !number of species in product term for mineral kinetics
     +   mrx=30,    !number of aqueous kinetic reactions     
     +   mechaq=8,  !number of aqueous kinetic mechanisms
     +   maqsp=8,   !number of species in aqueous kinetic rate eq.  
     +   nc=mpri,ng=mgas,nmgas=mgas,mspec=maqx,
     +   msecd=maqx+mmin+mgas+mexc+mads,
     +   mtot=mpri+msecd,ndim2=mtot,
     +   mnr=mpri+mmin+mgas+msurf+2,
     +   maqt=mpri+maqx,
     +   maxnpr=mpri,maxnse=maqx,maxnmi=mmin,maxnex=mexc,maxnad=mads)
c
cels6/26/08      parameter (mnod=mnel,nmnod=mnel,mbound=200,mrech=1,
cels6/26/08     +   maxb=1,maxqm=1,mpair=1000)
cels6/4/09      parameter (mnod=mnel,nmnod=mnel,mbound=200,maxb=1,
cels6/4/09         + maxqm=1,mpair=1000)
      parameter (mnod=mnel,nmnod=mnel,mbound=20)
C
c... Solid solution array dimensions
c            msol: maximum number of solid solutions
c            mcpss: maximum number of components in each
      parameter (msol=5, mcpss=6)
c
C---------------Explanation of the above array dimension variables-----------
C
C        MNEL=MNOD        Maximum number of grid blocks
C        MPRI=NC          Maximum number of primary species
C        MAQX=MAXNSE      Maximum number of secondary species
C        MMIN=NM          Maximum number of minerals
C        MGAS=NG          Maximum number of gases
!
C        MXsites          Maximum number of exchangeable sites         
!
C        MEXC             Maximum number of exchangeable species
C        MADS             Maximum number of secondary adsorbed species
C        MPADS            Maximum number of primary adsorbed species 
c        NSURF            Maximum number of sorption surfaces (minerals)
C        MBOUND           Maximum number of boundary (+ injecting) water types
C
C----------------------------------------------------------------------
c------ Constants
c... Gas constant in liter*bar/K/mole and in J/K/mol
        parameter (gc=0.083144d0, gc1=8.3144d0)
c... Moles of water per kg water
        parameter (rmh2o=55.50868d0)
c...Faraday constant, gas constant in J/K/mol
c     epsi is water dielec constant at 25C=78.5 times absolute vacuum permittivity =8.854e-12, epsi=6.95d-10
        parameter (faraday=96485.d0, dielec=78.5d0, permitt=8.854d-12,
     +     epsi=dielec*permitt,farsq=faraday*faraday,epsigc1=epsi*gc1)
c
      parameter (ppi=3.1415926536d0)
