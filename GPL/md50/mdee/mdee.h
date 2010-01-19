      implicit real*8 (A-H,O-Z)
      PARAMETER (NBASE  =   7  , NPART  = 4*NBASE**3       )
      PARAMETER (NTPS   =   5                              )
	PARAMETER (NS     = 100  , NTOT   = 4000              )
	parameter (NBSMAX = 64*NTOT, NBLMAX=800*NTOT         )
      PARAMETER (NB     = 200  , NA     = 200              ) 
      PARAMETER (NT     = 200  , NI     = 20               )
      PARAMETER (MAXS   = 200  , MAXGR  = NS*(NS-1)/2      )
      PARAMETER (NRQS   = 5000 , NKVMAX = 10000            )
      PARAMETER (MMX    = 2000 ,NNX    = 153, NBDMAX=40   )
      PARAMETER (LHIST  = 1000, NEMAX=100,   NRH=4+NEMAX, MXGRS=40)
*
      PARAMETER (AVSNO=6.02252D23,BOLTZ=1.38054D-23,FKCALJ=4.1868)
      PARAMETER (PI=3.1415926535897932D0,TODGR=180.D0/PI)
	parameter (ELECHG=1.602D-19,EPS0=8.854D-12)
      PARAMETER (A1=0.254829592D0,A2=-0.284496736D0,A3=1.421413741D0)
      PARAMETER (A4=-1.453152027D0,A5=1.061405429D0,B1=0.3275911D0)
      PARAMETER (TOTPI=1.128379169D0,PI2=2.D0*PI)
*
      COMMON /COMCOM/  X(NPART), Y(NPART), Z(NPART)
     X               ,PX(NPART),PY(NPART),PZ(NPART)
     X               ,TX(NPART),TY(NPART),TZ(NPART)
     X               ,QX(NPART),QY(NPART),QZ(NPART)
      COMMON /SITSIT/ SX(NTOT) ,SY(NTOT) ,SZ(NTOT)
     X               ,VX(NTOT) ,VY(NTOT) ,VZ(NTOT)
     X               ,FX(NTOT) ,FY(NTOT) ,FZ(NTOT)
     X               ,GX(NTOT) ,GY(NTOT) ,GZ(NTOT)
     X               ,HX(NTOT) ,HY(NTOT) ,HZ(NTOT)
     X               ,WX(NTOT) ,WY(NTOT) ,WZ(NTOT)
     X               ,OX(NTOT) ,OY(NTOT) ,OZ(NTOT)
	COMMON /EE/NE,ME,NDEL,NNH,EE(NEMAX),EC(NEMAX),
     +           IWALK(NEMAX,NEMAX)
	COMMON /LISTNB/ MBSH,MBLN,NBS1(NBSMAX),NBS2(NBSMAX),
     +NBL1(NBLMAX),NBL2(NBLMAX)
      DOUBLE PRECISION MASS,MASSD,MASSDI
      CHARACTER*4 NM,NAM,NAMOL*12
      CHARACTER*6 NAME
      CHARACTER*8 NMGR
      character*24 label,fdump,flst,frdf,ftcf,ftrk,fbio,fname*20
      integer TASKID
      common /NAMES/ label,fdump,flst,frdf,ftrk,fname,ftcf,fbio
      common /TIMESZ/ time0,timel,timeg,timef,timee,timen,timeb,timea,
     +               timet,timest,timev,times,time2
      COMMON /NBRNBR/ MNB(NTPS)    ,NNB(NTPS,MMX)
     X               ,INB(NTPS,MMX),JNB(NTPS,MMX,NNX)
      COMMON /TPSTPS/ NSPEC(NTPS),NSITS(NTPS),IEE(NTPS)
     X               ,IADDR(NTPS+1),ISADR(NTPS+1),ISADDR(NTPS+1)
     X               ,MDX(NTPS,NTPS),NRNOD(NTPS)
     X               ,NRB(NTPS),NRA(NTPS),NRT(NTPS),NRI(NTPS)
     X               ,IPOT(NTPS),NRTSC(NTPS),NAME(NTPS),NAMOL(NTPS)
      COMMON /MOLMOL/ SIGMA(NS),EPSIL(NS),CHARGE(NS),R(NS,3)
     X               ,SIX(NS,NS),TWL(NS,NS),ONE(NS,NS)
     X               ,EPS(NS,NS),SIG(NS,NS),ITS(NS)
     X               ,MASS(NS),MASSD(NS),NM(NS)
      COMMON /SITTES/ MASSDI(NTOT),Q(NTOT),ITYPE(NTOT)
     +                ,NAM(NTOT),NNUM(NTOT),NSITE(NTOT)
      COMMON /CONTRL/ TIM,TIMA,ICHE,ICHU,ME0
     +               ,NOP,NSITES,NTYPES,NSTOT,NFREQ,MOLINT
     X               ,NSTEP,NSTEPS,MSTEP,NSTTOT,NNAV,NAVT,NHIST,IHIST
     X               ,NAVER,IPRINT,NTRAJ,NSTEPA,NDUMP,IAVER
     +               ,ITREK,NTREK,NTRKF,LTREK,NTRKZ,ICFN,ICHNB,NCEE
     +               ,NUMTASK,TASKID
      COMMON /PHYSIC/ PE,TEMPF,ENERF,COULF,UNITL,UNITE,UNITM,UNITT
     X               ,UNITP,TOTMAS,DT,TSTEP,TSTEB,RHO,RCUTT,TOLER,BDTOL
     X               ,RF,FNST,TKE,TRTEMP,TDELT,AVB,AVA,AVT,TRPRES,TPRES
     +               ,TRYCK,TRID,TRYCKM,VOL,PLJ,PINT,PEL,PEST,PELS,PELR
     X               ,PES,PQE,POTEN,EKINE,ETOT,TEMP,FTIMES,FTSCL,FNOP
     +               ,TTR,TROT,TINT,FNSTI,FNSTR,RCUT,ELFAC,BETA,SPE
     X               ,POTES(NTPS*(NTPS+1)/2),POTLJ(NTPS*(NTPS+1)/2)
     X               ,PELRC(NTPS*(NTPS+1)/2),VRLRC(NTPS*(NTPS+1)/2)
     X               ,EKIN (NTPS),TEMPR(NTPS),SELFPE(NTPS),QTOT
     X               ,SUMMAS(NTPS),VOLF,BOXF,BOXYC,BOXY3 
     X               ,TFACT(NTPS),QSUM(NTPS),ENNO,PLOW
      common /VIRRR/ VIR,VIRL,VIRB,VIRA,VIRT,VIRN,VIRO,WIRF,WIRS,VIRS,
     +               VIRD,VIRLS,WIRSS
      COMMON /NVTNPT/ DQT,SC,SCV,DQP,SCL,SCLV
      COMMON /SISTEM/ BOXL,BOYL,BOZL,HBOXL,HBOYL,HBOZL
     X               ,BOXLI,BOYLI,BOZLI,RSSQ,SHORT,SHORT2,SHORTA,CLN
      COMMON /NBNBNB/ PES14(NTPS),PSR14(NTPS),C14LJ(NTPS),C14EL(NTPS),
     +                NNBB(NTOT),INBB(NTOT,NBDMAX)
      COMMON /BONBON/ BB(NB),EB(NB),RB(NB),FB(NB),EBOND(NTPS)
     X               ,DB(NB),RO(NB),IADB(NTPS),IB(NB),JB(NB),ID(NB)
      COMMON /ANGANG/ AA(NA),EA(NA),RA(NA),FA(NA),EANGL(NTPS)
     X               ,IADA(NTPS),IA(NA),JA(NA),KA(NA)
      COMMON /TORTOR/ TT (NT),ET (NT),ETORS(NTPS)
     X               ,FT1(NT),FT2(NT),FT3 (NT)
     X               ,RT (NT),FT (NT),NMUL(NT)
     X               ,IADT(NTPS),IT(NT),JT(NT),KT(NT),LT(NT)
     X               ,ITORS(NTPS)
      COMMON /IMPIMP/ TI(NI),EI(NI),RI(NI),FI(NI),EIMP(NTPS)
     X               ,IADI(NTPS),IM(NI),JM(NI),KM(NI),LM(NI)
      COMMON /DIPDIP/ DPX(NPART),DPY(NPART),DPZ(NPART)
      COMMON /VECVEC/ UPX(NPART),UPY(NPART),UPZ(NPART)
     X               ,APX(NPART),APY(NPART),APZ(NPART)
     X               ,BPX(NPART),BPY(NPART),BPZ(NPART)
     X               ,CPX(NPART),CPY(NPART),CPZ(NPART)
      COMMON /ROTROT/ RMX(NPART,3,3)
      COMMON /AVRAWR/ AV(NRQS,NEMAX),AW(NRQS,NEMAX),ITE(NEMAX,LHIST),
     +                HIST(NRH,LHIST)
      COMMON /AUXAUX/ XR(NTOT),  YR(NTOT),  ZR(NTOT)
     X               , A(NTOT),   B(NTOT),   C(NTOT), D(NTOT)
     X               ,RSQI(NTOT),FSRN(NTOT),FESN(NTOT)
      LOGICAL        LSCLT,LSCFT,LNVT,LRST,LDMP,LCHL
     X               ,LNPT,LXE,LCSRDF
     X               ,LCHECK,L0AVS,L0VS,LSHEJK,LCHP,LCHT,LEQ
     X               ,L14NB,L15NB,LXMOL,LECH,LMNBL,LHEX,LOCT
      LOGICAL         LGR,LGDMP,LGRST,ALLSTS
      COMMON /LOGLOG/ LSCLT,LSCFT,LNVT,LRST,LDMP,LPLT,LCHL
     X               ,LXE,LNPT,LXMOL,LCSRDF,LMNBL
     X               ,LCHECK,L0AVS,L0VS,LSHEJK,LCHP,LCHT,LEQ
     X               ,L14NB(NTPS),L15NB(NTPS),LECH,LNEWE,LHEX,LOCT
      COMMON /RDFRDF/ RDFINC,DRDFI,RDFCUT,LGR,LGDMP,LGRST,ALLSTS
     X               ,MAXRDF,MAX,IGR(MAXGR),JGR(MAXGR),IIN(MAXGR)
     X               ,NRDF,IDX(NS,NS),IRDF(MAXS,MAXGR),NMGR(0:MAXGR)
     +		         ,NGRI(NS),IGRI(NS,MXGRS),MGRI(NS,MXGRS)
	common /FURFUR/ AK(NKVMAX),RX(NKVMAX),RY(NKVMAX),RZ(NKVMAX)
     +               ,FEXP,ALPHA,ALPHAD,NKV
       save
