C------------------------ PERM_v2.INC ----------------------
c... Include file for Permeability-Porosity-Capillary Pressure
c           Parameters in TOUGHREACT V2.0
C----------------------------------------------------------------------
      common/permrct/perm(3,mnel),perm0(3,mnel)
      common/phirct/phi0(mnel),pcfact(mnel)
C-For monitoring porosity and perm. changes but not affecting flow, Tianfu Xu, 04/23/2002 
C    Save for printout
      common/permsav/phim(mnel),permm(3,mnel)
c
