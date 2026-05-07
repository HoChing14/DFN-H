C*********************************************************************
C*******  file 'T2' for module t2cg2.f of TOUGH2, Version 2.0  *******
C*********************************************************************
C*                                                                   *
C*    ASSIGN PARAMETERS FOR FLEXIBLE DIMENSIONING OF LARGE ARRAYS    *
C*                                                                   *
C*********************************************************************
C
C
C---  MNEL: Number of elements
cc      PARAMETER (MNEL = 60000)
cc      PARAMETER (MNEL = 50000)
      PARAMETER (MNEL = 50000)
C
C---  MNCON: Number of connections
cc       PARAMETER (MNCON = 180000)
cc       PARAMETER (MNCON = 150000)
       PARAMETER (MNCON = 180000)
C
C---  MNEQ: Number of equations per element (# of Primary variables)
cc      PARAMETER (MNEQ = 9)
      PARAMETER (MNEQ = 4)
C
C---  MNK: Number of components
cc      PARAMETER (MNK = 8)
      PARAMETER (MNK = 3)
C
C---  MNPH: Number of phases which can be present
cc      PARAMETER (MNPH = 3)
      PARAMETER (MNPH = 2)
C
C---  MNB: Number of phase-dependent secondary variables
      PARAMETER (MNB = 12)
C
C---  MNBK = Total number of phase-dependent secondary variables
      PARAMETER (MNBK = MNB+MNK)
C
C---  MNOGN: Number of sources/sinks
      PARAMETER (MNOGN = 200)
C
C---  MGSTAB: Number of data point pairs in tables of time-variable
C             injection and/or production rates
      PARAMETER (MGTAB = 10000)
C
C---  nredm: The total number of equations (order of the full matrix)
      PARAMETER (NREDM = MNEL*MNEQ)
C
C---  MNZ: The total number of non-zero entries in the full matrix
      PARAMETER (MNZ = (MNEL+2*MNCON)*MNEQ*MNEQ)
C
C---  NRWORK: The size of the REAL work array for the CG solvers
      PARAMETER (NRWORK = 1000+MNZ+38*nredm)
C
C---  NRWORK: The size of the integer*8 work array for the CG solvers
      PARAMETER (NIWORK = 32+mnz+5*nredm)
C
C---  NDIM: The number of active dimensions (1,2,or 3)
C
cels4/28/06 not used      PARAMETER (NDIM = 3)
C
C---  NCONUP: The total number of connections of elements with
C---          other elements of larger element numbers
C
      PARAMETER (NCONUP = 11*MNCON/10)
C
C---  NCONDN: The total number of connections of elements with
C---          other elements of smaller element numbers
C
      PARAMETER (NCONDN = 11*MNCON/10)
C
      COMMON/D4ARAY3/NX0(MNEQ,MNEL+1)
      COMMON/D4ARAY4/NB0(MNEQ,MNEQ,MNEL+1)
      COMMON/D4ARAY5/NA0(NCONDN,MNEQ,MNEQ)
      COMMON/D4ARAY6/NC0(NCONUP,MNEQ,MNEQ)
C
      COMMON/D4COM1/IIJJKK(MNEL+1)
      COMMON/D4COM2/NO(MNEL+2*MNCON+10)
      COMMON/D4COM3/IJKMM(MNEL+1)
      COMMON/D4COM4/IJKPP(MNEL+1)
      COMMON/D4COM5/ISUMPP(0:MNEL+1)
      COMMON/D4COM7/ISUMMM(0:MNEL+1)
C
      integer*8 BNDWTH
      COMMON/D4COM6/BNDWTH,NCRTM,NEQRM,NEQUH,IDIM,IRELSZ
C
C
C*********************************************************************
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C*********************************************************************
C---------------------------------------------------------------------------
        parameter(maxmat = 100)      ! in TOUGH2 V2, maxmat = 27
C
cels10/8/98 number of time steps to skip to write save file
        parameter (ntsave =1)
c
