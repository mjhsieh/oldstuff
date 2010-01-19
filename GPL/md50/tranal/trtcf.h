*   Tranal v 5.0
*   header file for TCF computations (tcf.f) 
*
      parameter (MAXCF=1000,NTCF=12)
      parameter (MAXTCF=(MAXCF+2)*3*NPART,MCF=MAXCF*NTPS)
      real*4 CC1,CC2,CC3,CC4,CC5,CC6
      character*6 NAMECF, NMM*2
      common /MOLEC/ DPX(NPART),DPY(NPART),DPZ(NPART)
     +              ,UPX(NPART),UPY(NPART),UPZ(NPART)
     +              ,RMX(NPART,3,3)
      COMMON /TCFTCF/ CC1(MAXTCF),CFV(MCF)
     X               ,CC2(MAXTCF),CFA(MCF)
     X               ,CC3(MAXTCF),CP1(MCF),CP2(MCF)
     X               ,CC4(MAXTCF),CU1(MCF),CU2(MCF)
     +               ,CC5(MAXTCF),CC6(MAXTCF)
     X               ,CVX(MCF),CVY(MCF),CVZ(MCF)
     X               ,CAX(MCF),CAY(MCF),CAZ(MCF)
     X               ,CCF(MCF)
     X               ,LCF,LCFDMP,LCFRST
     X               ,INX(MAXCF),NCALLS,NSTEG,NSTCF,NOM
     X               ,N1(NTPS),N2(NTPS),ITCF(18)
     X               ,NAMECF(12)
 
