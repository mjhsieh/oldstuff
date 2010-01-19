*  Main include file for tranal v.5.0 programs
      implicit real*8 (A-H,O-Z)
* max number of mol.types, sites, molecules  and atoms
      PARAMETER (NTPS=10, NS=1000, NPART  = 10000, NTOT=50000)
* max number of bonds
      PARAMETER (NBM=20000)
* max number pof trajectory files
      PARAMETER (MAXFIL=1000)
*
*  Math and Physical constants
*
      PARAMETER (AVSNO=6.02252D23,BOLTZ=1.38054D-23,FKCALJ=4.1868)
      PARAMETER (PI=3.1415926536D0,ELECHG=1.602D-19,EPS0=8.854D-12)
      PARAMETER (TOTPI=1.128379169D0,PI2=2.D0*PI)
      PARAMETER (TODGR=180.D0/PI)
*                                 
      real*4 DTN
      real*8 MASS,MASSD,MASSDI
      logical LXMOL,LVEL,LFILTR,LBOND
      character*128  TAKESTR,STR,NAMOL,FILNAM,FILXM,FXMOL
      character NAME*6,NFORM*4,NM*4
      COMMON /COORD/ SX(NTOT),SY(NTOT),SZ(NTOT),
     +               OX(NTOT),OY(NTOT),OZ(NTOT),
     +               X(NPART),Y(NPART),Z(NPART),
     +               WX(NTOT),WY(NTOT),WZ(NTOT),
     +               PX(NTOT),PY(NTOT),PZ(NTOT),
     +               QX(NTOT),QY(NTOT),QZ(NTOT),
     +               VX(NTOT),VY(NTOT),VZ(NTOT)
      COMMON /TPSTPS/ NSPEC(NTPS),NSITS(NTPS),ITS(NS),ITYPE(NTOT)
     X               ,IADDR(NTPS+1),ISADR(NTPS+1),ISADDR(NTPS+1)
     +               ,NNUM(NTOT),NSITE(NTOT),LIST(NTPS)
     +               ,NUMR(NTOT),NUMO(NTOT),ITM(NPART)
      COMMON/INT/ ISTEP,IBREAK,NF
      COMMON /MOLPAR/ SIGMA(NS),EPSIL(NS),CHARGE(NS),R(NS,3)
     X               ,MASS(NS),MASSD(NS),SUMMAS(NTPS)
     +               ,MASSDI(NTOT),Q(NTOT),TIMIN(MAXFIL)
      COMMON/MISC/BREAKM,FULTIM1,VFAC,TEMP
      COMMON /CONTRL/ NOP,NSITES,NTYPES,NSTOT,IAN,IPRINT,NUMTR,DTN
      common /NAMES/ NAMOL(NTPS),NAME(NTPS),NM(NS),FILNAM(MAXFIL)
     +  ,LLNX,NFORM,FILXM,FXMOL
      COMMON /SIM/ BOXL,BOYL,BOZL,HBOXL,HBOYL,HBOZL,TOTMAS,FULTIM,VOL
      common /BONDS/ID(NBM),IB(NBM),JB(NBM),NRB(NTPS),IADB(NTPS),NRBB
      common/LOGLOG/LXMOL,LF,LVEL,LFILTR,LBOND

