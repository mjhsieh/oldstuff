C  This is the include file containing static dimensioning
C  of the working fields. Practically all modules depend somehow on it, 
C  so you must recompile the sources each time you change this file.
C
      IMPLICIT REAL*8 (A-H,O-Z) 
*
*   1. Define maximum sizes of working arrays
*   -----------------------------------------
C   This part is now in a separate file dimpar.h
      include "dimpar.h"
*
*
*   2. Constants
*   ------------
C   Specfied here are constants available in all program modules which
C   have this file included
C
C   Physical constants:
C   Avagadro number, Boltzmann constant, J/cal conversion factor:
      PARAMETER (AVSNO=6.02252D23,BOLTZ=1.38054D-23,FKCALJ=4.1868)
C   Electron charge, Vaccuum "dielectric" constant:  
      PARAMETER (ELECHG=1.602D-19,EPS0=8.854D-12)
C   Parameters to calculate erfc function      
      PARAMETER (A1=0.254829592D0,A2=-0.284496736D0,A3=1.421413741D0)
      PARAMETER (A4=-1.453152027D0,A5=1.061405429D0,B1=0.3275911D0)
C   Pi and angles      
      PARAMETER (PI=3.1415926536D0,TOTPI=1.128379169D0,PI2=2.D0*PI)
      PARAMETER (TODGR=180.D0/PI,ODIN=1.d0)
      parameter (cos30=0.8660254037844386,BOXYC=0.5/cos30,SQ3=2.*cos30)
*      
*   3. Variables
*   ------------
C   Specfied here are common variables available in all program modules which
C   have this file included
*
*   3.1. Names (character)
C
C   NM     - name of sites
C   NAMOL  - name of molecules;  NAME - short molecule name
C   NMGR   - name of RDF
C   PATHDB - path to molecular database
C   Lenght of NAMOL should be sinchronized with p.5.2, input.f
      CHARACTER NM*4,NAMOL*32,PATHDB*112
      CHARACTER*6 NAME
      CHARACTER*8 NMGR
C   label -  program version
C   File names:
C        fdump  - restart file
C        frdf   - RDF file
C        ftcf   - Time correlation functions (TCF) file
C        ftrk   - Trajectory file
C        fmol   - Coordinate file for molecular viewers
C        filref - Some atoms may be bind to some locations, specified in 
C                 this file 
      character*24 label,fdump,frdf,ftcf,ftrk,fmol,filref,fname*20
      common /NAMES/ label,fdump,frdf,ftrk,fname,ftcf,fmol,PATHDB,filref
C       masses: MASS - at. units; MASSD - dimensional; MASSDI - inverse
*
*    3.2 Variables
*
C    More information on used variables are given in appropriate parts of program
C
      REAL*8 MASS,MASSD,MASSDI
C   Integers  
C   Some important numbers:
C         NOP    - total molecules
C         NSITES - total sites
C         NTYPES - total molecule types
C         other - see description in appropriate part of the source
      COMMON /INTEG/ NOP,NSITES,NTYPES,NSTOT,NFREQ,MOLINT,IPRINN
     X               ,NSTEP,NSTEPS,MSTEP,NSTTOT,NNAV,NAVT,NHIST,IHIST
     X               ,NAVER,IPRINT,NTRAJ,NSTEPA,NDUMP,IAVER,ICHECK
     +               ,ITREK,NTREK,NTRKF,LTREK,NTRKZ,ICFN,ICHNB,KMAX  
     +	     ,NRBT,NRAT,NRRT,NRIT,NRAA,NRBB,NRTT,NRII,NHDM,IO,
     +        NNAD,NRBND,NTYP1,ICELL,NPULL,IEXT,ISTAV,ICMOM
C   energy/temperature related variables     
      COMMON /ENERGY/ PE,TEMPF,ENERF,COULF,UNITE,UNITM,UNITT
     X               ,UNITP,TOTMAS,DT,TSTEP,TSTEB,HTSTEP,HTSTEB,RHO
     X               ,FNST,TKE,TRTEMP,TDELT,TRPRES,FEXP,TOLER,BDTOL
     +               ,TRYCK,TRID,TRYCKM,PINT,CONVET,RTP
     X               ,TEMP,FTIMES,FTSCL,FNOP,RCUTT,EXAMPL,EXFREQ
     +               ,EABSS,EABS,EFIELD,POTEN
     +               ,TTR,TROT,TINT,FNSTI,FNSTR,RCUT,TEMPV,FACVT,RLR
     +               ,PE1,PELS1,PE2,PELS2,SPE,TTREK,EFACT,PERMOL
     +		     ,TRYCKX,TRYCKY,TRYCKZ,RFF,RFF2
C   Virials
      common /VIRRR/ VIRB,VIRA(3),VIRAN(3),WIRS,WIRSS,VIR1,VIR2,VIRD,
     +               VIRX,VIRY,VIRZ,VIRFX,VIRFY,VIRFZ
C   Control NVT/NPT ensembles
      COMMON /NVTNPT/ QT,QPR,DQT,SC,SCV,DQP,SCL,SCLX,SCLY,SCLZ,
     +                SCM(NTPS)
C   Size related quantities
      COMMON /SIZES/  UNITL,BOXL,BOYL,BOZL,HBOXL,HBOYL,HBOZL,ALPHA
     X               ,ALPHAD,RSSQ,SHORT,SHORT2,SHORTA,BOXY3,RDNA,VOL
C   These numbers keep CPU time for different part of the program     
      common /TIMESZ/ TIM,TIMA,time0,timel,timeg,timef,timee,timen,
     +timeb,timea,timet,timest,timev,times  
C   Logical constants      
      LOGICAL        LSCLT,LSCFT,LNVT,LRST,LDMP,LCHL,LFORMI,LFORMO    
     X               ,LNPT,LMOVE,LRR,LSEP,LOCT,LHEX,LBONDL,LCMOM
     X               ,LINPUT,L0AVS,L0VS,L0CPU,LSHEJK,LCHP,LCHT,LEQ
     X               ,L14NB,L15NB,LXMOL,LMNBL,LGR,LGDMP,LGRST,LGEOM
     +               ,LMOL1,LMOL2,LCOMD,LKONG,LGATHER,LSCTD,LVISUAL
      LOGICAL         LCF,LCFDMP,LCFRST              ! in common TCFTCF
      COMMON /LOGLOG/ LSCLT,LSCFT,LNVT,LRST,LDMP,LCHL,LFORMI,LFORMO
     X               ,LNPT,LXMOL,LRR,LSEP,LOCT,LHEX,LBONDL,LCMOM
     X               ,LINPUT,L0AVS,L0VS,L0CPU,LSHEJK,LCHP,LCHT,LEQ
     X               ,L14NB(NTPS),L15NB(NTPS),LMOVE(NTPS),LMNBL,LKONG        
     +               ,LMOL1,LMOL2,LCOMD,LGEOM,LGATHER,LSCTD,LVISUAL
*
*   5. Arrays
*   ---------
C   Specfied here are common arrays available in all program modules which
C   have this file included. Dimensions of these arrays are defined above
C   in Sec.1
C
*     Molecules:
C   X,Y,Z              - molecular center of mass coordinates
C   PX,PY,PZ           - momentum of the molecule
C   QX,QY,QZ           - angular molecular momemtum
C   DPX,...            - dipole momentum
C   UPX,...            - unit vector fixed by two selected atoms
C   XPL,...            - average coordinates of "linked" atoms
        COMMON /COMCOM/  X(NPART), Y(NPART), Z(NPART)
     X               ,PX(NPART),PY(NPART),PZ(NPART)
     X               ,QX(NPART),QY(NPART),QZ(NPART)
     +               ,DPX(NPART),DPY(NPART),DPZ(NPART)
     +               ,UPX(NPART),UPY(NPART),UPZ(NPART)
     +               ,XPL(NPLM),YPL(NPLM),ZPL(NPLM)   
C     Atoms:
C   SX,SY,SZ           - atom coordinates
C   VX,VY,VZ           - atom momenta 
C   FX,FY,FZ           - forces (sometimes used as temporary arrays)
C   GX,GY,GZ           - slow forces
C   HX,HY,HZ           - fast forces
C   WX,WY,WZ           - local atom coordinates (relatively molecular center-of-mass)
C   OX,OY,OZ           - temporary array for atom coordinates or velocities
C   MASSDI             - (mass)**(-1)  (internal units)
C   Q                  - charge
      COMMON /SITSIT/ SX(NTOT) ,SY(NTOT) ,SZ(NTOT)
     X               ,VX(NTOT) ,VY(NTOT) ,VZ(NTOT)
     X               ,FX(NTOT) ,FY(NTOT) ,FZ(NTOT)
     X               ,GX(NTOT) ,GY(NTOT) ,GZ(NTOT)
     X               ,HX(NTOT) ,HY(NTOT) ,HZ(NTOT)
     X               ,WX(NTOT) ,WY(NTOT) ,WZ(NTOT)
     X               ,OX(NTOT) ,OY(NTOT) ,OZ(NTOT)
     +               ,MASSDI(NTOT), Q(NTOT)
C   This arrays handle lists of neighbours    
	COMMON /LISTNB/MBSH,MBLN,NBS1(NBSMAX),NBS2(NBSMAX),
     +NBL1(NBLMAX),NBL2(NBLMAX),NNBB(NTOT),INBB(NTOT,NBDMAX)
C   Addresses and references
C   See file units.f for details
	common /ADDREFS/ IADDR(NTPS+1),ISADR(NTPS+1),ISADDR(NTPS+1),
     +                 ITM(NPART),ITYPE(NTOT),NNUM(NTOT),NSITE(NTOT),
     +                 MDX(NTPS,NTPS),ITS(NS),INPL(NPLM)
C   Molecule type characteristics   
C         NSPEC  - total molecules of this type
C         NSITS  - total atoms in a molecule of this type
C         ISHEJK - which molecules (1) are constrained
      COMMON /MOLTYP/ NSPEC(NTPS),NSITS(NTPS),IINIT(NTPS)              
     X               ,LIST(NTPS),ISHEJK(NTPS)
     X               ,NRB(NTPS),NRA(NTPS),NRT(NTPS),NRI(NTPS)
     X               ,IPOT(NTPS),NRTSC(NTPS),NAME(NTPS),NAMOL(NTPS)
C    some arrays for intermediate data     
      common /SOMMOL/ POTES(NTPP),POTLJ(NTPP),POTE1(NTPP),POTL1(NTPP)
     X               ,PELRC(NTPP),VRLRC(NTPP),SPEE(NTPS)
     X               ,EKIN (NTPS+1),TEMPR(NTPS+1),SELFPE(NTPS)
     X               ,SUMMAS(NTPS),C14LJ(NTPS),C14EL(NTPS)
     X               ,TFACT(NTPS),QSUM(NTPS)
     +           ,PES14(NTPS),PSR14(NTPS),PES141(NTPS),PSR141(NTPS)
C   Site characteristics
C        R  - initial local site coordinates from .mol file
C        NM - names of atom sites
      COMMON /SITTES/ SIGMA(NS),EPSIL(NS),CHARGE(NS),R(NS,3)
     +               ,EPAD(NNADM),SIGAD(NNADM),ILJ(NNADM)
     X               ,MASS(NS),MASSD(NS),NM(NS)         
C   Everything on bonds      
      COMMON /BONBON/ BB(NB),EB(NB),RB(NB),FB(NB),EBOND(NTPS),BR(NB)
     +	       ,BE(NB),DB(NB),RO(NB),IADB(NTPS),IB(NB),JB(NB),ID(NB)
     +	       ,IBI(NBO),JBJ(NBO),IBUK(NBO)
C   Everything on angles
      COMMON /ANGANG/ AA(NA),EA(NA),RA(NA),FA(NA),EANGL(NTPS),AR(NA)
     X               ,AE(NA),IADA(NTPS),IA(NA),JA(NA),KA(NA)
     +	       ,IAI(NAO),JAJ(NAO),KAK(NAO),IAUK(NAO)
C   Everything on torsions and impropers     
      COMMON /TORTOR/ TT (NT),ET (NT),ETORS(NTPS),FT1(NT),FT2(NT)
     X               ,FT3(NT),FT4(NT),FT5 (NT),TR(NT),TE(NT)
     X               ,RT (NT),FT (NT),NMUL(NT)
     X               ,IADT(NTPS),IT(NT),JT(NT),KT(NT),LT(NT)
     X               ,ITORS(NT),ITI(NTO),JTJ(NTO),KTK(NTO),LTL(NTO),
     +                ITUK(NTO)
C   Rotation matrix for each molecule     
      COMMON /ROTROT/ RMX(NPART,3,3)
C   Accumulators for averages and history of averages     
      COMMON /AVRAWR/ AV(NRQS),AW(NRQS),HIST(NRQS,LHIST)
C   RDF-calculations related arrays/variables      
      COMMON /RDFRDF/RDFINC,DRDFI,RDFCUT,RDFCUT2,LGR,LGDMP,LGRST
     X               ,MAXRDF,MAX,IGR(MAXGR),JGR(MAXGR),IIN(MAXGR)
     X               ,NRDF,IRDF(MAXS,MAXGR),NMGR(0:MAXGR)
     +		         ,NGRI(NS),IGRI(NS,MXGRS),MGRI(NS,MXGRS)
C   TCF-related arrays/variables
      CHARACTER*6 NAMECF
      PARAMETER (MCF=MAXCF*NTPS) 
	real*4 CC1,CC2,CC3,CC4,CC5,CC6
      COMMON /TCFTCF/ CC1(MAXTCF),CFV(MCF)
     X               ,CC2(MAXTCF),CFA(MCF)
     X               ,CC3(MAXTCF),CP1(MCF),CP2(MCF)
     X               ,CC4(MAXTCF),CU1(MCF),CU2(MCF)
     +               ,CC5(MAXTCF),CC6(MAXTCF)
     X               ,CVX(MCF),CVY(MCF),CVZ(MCF)
     X               ,CAX(MCF),CAY(MCF),CAZ(MCF)
     X               ,CCF(MCF)
     X               ,LCF,LCFDMP,LCFRST
     X               ,INX(MAXCF),NCALLS,NSTEG,NSTCF,NOM,JUMP
     X               ,N1(NTPS),N2(NTPS),ITCF(NTCF+6)
     X               ,NAMECF(NTCF)
*
*  6. Parameters for parallel execution (MPI)
*  ------------------------------------------
      PARAMETER (NBUFF  = 3*NTOT+2*NA+2*NB+2*NT)        
      integer NUMTASK,TASKID 
	integer NAB(0:NODES),NAE(0:NODES),NABS(0:NODES),NAP(0:NODES),
     +NAP3(0:NODES),NABS3(0:NODES)
      common /PARA/ NUMTASK,TASKID,IDG,LENMP,LENTT,LENIT,NBUF,NAP3,
     +NAP,NAB,NAE,NABS,NORD(NTPS),LENMB,LENMA,LENMT,LENMI,MAST,NAPC,
     +NABS3
*    for safety 
	save




