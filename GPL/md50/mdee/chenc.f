*======================= CHENS ================================
*
	subroutine CHENS
	include "mdee.h" 
	real*4 rand 
	external rand
	if(rand().lt.0.5)then
	  ME1=ME-1
	else
	  ME1=ME+1
	end if
C	if(EC(ME).lt.1.d-7.and.LNPT)then
C	  do ITYP=1,NTYPES
C	    if(IEE(ITYP).eq.0)then
c	      do I=ISADDR(ITYP)+1,ISADDR(ITYP+1)
c		if(abs(SX(I)).gt.HBOXL)return
c		if(abs(SY(I)).gt.HBOYL)return
c		if(abs(SZ(I)).gt.HBOZL)return
c	      end do
c	    end if
c	  end do
c	end if  
c	if(EC(ME1).lt.1.d-7.and.LNPT)then
c	  do ITYP=1,NTYPES
c	    if(IEE(ITYP).eq.0)then
c	      do I=ISADDR(ITYP)+1,ISADDR(ITYP+1)
c		if(abs(SX(I)).gt.BOXF)return
c		if(abs(SY(I)).gt.BOYF)return
c		if(abs(SZ(I)).gt.BOZF)return
c	      end do
c	    end if
c	  end do
c	end if  
	ICHE=ICHE+1
	if(ME1.lt.1.or.ME1.gt.NE)return
	time20=cputime(0.d0)     
	call DIFFUR(ME1,EFN,EFO)
	call DIFEN(ME1,EEN,EEO)
	DEF=EFN-EFO
	DEE=EEN-EEO  
	DE=(DEE+DEF)*BETA+EE(ME1)-EE(ME)
	FACTR=  0.001*ENERF/FNOP
	ENNN=(EEN+EFN)*FACTR
	ENNO=(EEO+EFO)*FACTR
	if(IPRINT.ge.6)write(*,'(a5,I5,a4,I5,a4,f12.4,3x,3f11.3)')
     + 'from ',ME,' to ',ME1,' DE=',DE,ENNN,ENNO,EE(ME1)-EE(ME)
	if(DE.gt.0)then
	  if(DE.gt.20)go to 100
	  if(exp(-DE).lt.rand())go to 100
	end if
	ME=ME1
	ICHU=ICHU+1
	TRPRES=TPRES+(1.d0-EC(ME))*NDEL*BOLTZ*TEMP*1.d-5/(VOL*UNITL**3)
	call RECINT
      do I=1,ME-1
        if(IWALK(I,ME).le.IWALK(ME,I))
     +  IWALK(I,ME)=IWALK(I,ME)+1
      end do
      do I=ME+1,NE
        if(IWALK(I,ME).lt.IWALK(ME,I))
     +  IWALK(I,ME)=IWALK(I,ME)+1
      end do
 100	time2=time2+cputime(time20)
	return
	end                     
*	
*======================= RECINT ===============================
*
      subroutine RECINT
      include "mdee.h"
      NDEL=0
      do ITYP=1,NTYPES
	if(IEE(ITYP).eq.0)then
	  ESCI=EC(ME)
	  NDEL=NDEL+NSPEC(ITYP) 
	else
	  ESCI=1.d0
	end if   
        DO I       =  ISADR(ITYP)+1,ISADR(ITYP+1)
	  do JTYP=1,NTYPES
	    if(IEE(JTYP).eq.0)then
	      ESCJ=EC(ME) 
	    else
	      ESCJ=1.d0
	    end if   
	    ESCL=ESCI*ESCJ
            DO J       =  ISADR(JTYP)+1,ISADR(JTYP+1)
              EPSI       =   ESCL*DSQRT(EPSIL(I))*DSQRT(EPSIL(J))
     X        *1.D3/AVSNO/UNITE
              SIGM       =  (SIGMA(I)+SIGMA(J))/2.d0
              TSIX       = -4.D0*EPSI*SIGM**6
              TTWL       =  4.D0*EPSI*SIGM**12
              ONE(I,J)   =  ESCL*CHARGE(I)*CHARGE(J)*COULF
              SIX(I,J)   =  TSIX
              TWL(I,J)   =  TTWL
              EPS(I,J)   =  DSQRT(EPSIL(I))*DSQRT(EPSIL(J))
     X        *1.D3/AVSNO/UNITE
              SIG(I,J)   =  SIGM
            END DO! OF J
          END DO
	end do
      end do
      return
      end
*
*===================== DIFFUR ===========================================
*
	subroutine DIFFUR(ME1,EFN,EFO)
	include "mdee.h"
      DATA TWOPI/6.2831853072D0/
*
	IF(NTYPES.EQ.1.AND.NSPEC(1).EQ.1) RETURN
      timef0	= cputime(0.d0)
*
      QSS         = 0.D0
      CL	= BOXL
      CFUR=4.d0*PI*COULF/VOL     
      DO ITYP     = 1,NTYPES
        QSS         = QSS+QSUM(ITYP)
      END DO! OF ITYP
      IF(QSS.EQ.0.D0) RETURN
*
      EFN     = 0.D0
      EFO	= 0.d0 
      ESFF=0.d0
*   reciprocal space Evald
      do IKV=1,NKV
        SCSO=0.
        SSNO=0.
	  SCSN=0.
	  SSNN=0.
 	  RKX=RX(IKV)
	  RKY=RY(IKV)
	  RKZ=RZ(IKV)
C$DIR NO_RECURRENCE
        do I=1,NSTOT
	    IC=NSITE(I)
	    ITYP=ITYPE(I)
	    SCP=RKX*SX(I)+RKY*SY(I)+RKZ*SZ(I)
            SSC=sin(SCP)
	    CSC=cos(SCP)
	    if(IEE(ITYP).eq.0)then
	      QO=CHARGE(IC)*EC(ME)
	      QN=CHARGE(IC)*EC(ME1) 
	    else
	      QO=CHARGE(IC)
	      QN=QO
	    end if
          SSNO=SSNO+QO*SSC
          SCSO=SCSO+QO*CSC
          SSNN=SSNN+QN*SSC
          SCSN=SCSN+QN*CSC
        end do
	  AKC = CFUR*AK(IKV)
        EFO=EFO+(SSNO**2+SCSO**2)*AKC
        EFN=EFN+(SSNN**2+SCSN**2)*AKC
      end do
*	FACTR=  0.001*ENERF/FNOP
*	if(IPRINT.ge.7)write(*,'(a,4f12.5)')' FUR:',FACTR*EFN,FACTR*EFO
	EFF1=EFO
	EFF2=0.
	EFF3=0.      
      DIF	= EFN-EFO 
*  From eterms 
	ABC = ALPHAD/DSQRT(PI) 
	do ITYP=1,NTYPES
*
        ISP         = ISADDR (ITYP)
        ISBEG       = ISADR  (ITYP)+1
        ISEND       = ISADR  (ITYP+1)
* molecules
        DO 300 I    = 1,NSPEC(ITYP)
* first atom
          DO 200 IS   = ISBEG,ISEND
            ISP         = ISP+1
            QIIO         = ONE(IS,IS)
	      if(IEE(ITYP).eq.0)then
		QIIN	= (CHARGE(IS)*EC(ME1))**2*COULF
	      else
		QIIN	= QIIO
	      end if 
	      EFN=EFN-ABC*QIIN
	      EFF2=EFF2-ABC*QIIN
	      EFO=EFO-ABC*QIIO
            JSP         = ISP
            ISB         = IS+1
*
            IF(ISB.LE.ISEND) THEN
* second atom in the same molecule
              DO 100 JS   = ISB,ISEND
      	  JSP         = JSP+1
      	  QIJO        = ONE(IS,JS)
	          QIJN	= CHARGE(IS)*CHARGE(JS)*COULF
	 	  if(IEE(ITYP).eq.0)QIJN=QIJN*EC(ME1)**2
      	  DX          = SX(ISP)-SX(JSP)
      	  DY          = SY(ISP)-SY(JSP)
      	  DZ          = SZ(ISP)-SZ(JSP)
       	  if(DX.gt. HBOXL)DX=DX-BOXL
        	  if(DX.lt.-HBOXL)DX=DX+BOXL
        	  if(DY.gt. HBOYL)DY=DY-BOYL
        	  if(DY.lt.-HBOYL)DY=DY+BOYL
        	  if(DZ.gt. HBOZL)DZ=DZ-BOZL
        	  if(DZ.lt.-HBOZL)DZ=DZ+BOZL
*	    	  call PBC(RX,RY,RZ)
      	  R2          = DX**2+DY**2+DZ**2
      	  R1          = DSQRT(R2)
      	  R1I         = 1.d0/R1
      	  ALPHAR      = ALPHAD*R1
      	  TTT         = 1.D0/(1.D0+B1*ALPHAR)
      	  EXP2A       = DEXP(-ALPHAR**2)
      	  ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A
      	  EESO  = QIJO*(1.0D0-ERFC)*R1I
      	  EESN  = QIJN*(1.0D0-ERFC)*R1I
      	  EFN=EFN-EESN
	          EFF3=EFF3-EESO
	          EFO=EFO-EESO
 100  	  CONTINUE  ! end do of JS
            END IF
  200     CONTINUE
  300   CONTINUE
	END DO ! of ITYP 
*	FACTR=  0.001*ENERF/FNOP
*	if(IPRINT.ge.7)write(*,'(a,4f12.5)')' FURtot:',FACTR*EFN,FACTR*EFO
      return
      end  
*
*====================== DIFEN ========================================
*                                                                     
	subroutine DIFEN(ME1,EEN,EEO)
	include "mdee.h"
	EEN=0.d0
	EEO=0.d0         
	ELJO=0.
	ESTO=0.
	ELJN=0.
	ESTN=0.
	DO INP      = 1,MBLN
	  ISP=iabs(NBL1(INP))
	  JSP=NBL2(INP)
	  ITYP=ITYPE(ISP)
	  JTYP=ITYPE(JSP)
	  ISB =NSITE(ISP)
  	  JSB =NSITE(JSP)
	  IMOL=NNUM(ISP)
	  JMOL=NNUM(JSP)
	  if(IMOL.ne.JMOL)then
  	  A6O =SIX(ISB,JSB)
	  B12O=TWL(ISB,JSB)
	  QIJO=ONE(ISB,JSB)
	  if(IEE(ITYP).eq.0)then
	    ESCI=EC(ME1) 
	  else
	    ESCI=1.d0
	  end if
	  if(IEE(JTYP).eq.0)then
	    ESCJ=EC(ME1) 
	  else
	    ESCJ=1.d0
	  end if
	  ESCL=ESCI*ESCJ
	  SS6	= (SIG(ISB,JSB))**6
	  EPSI	= ESCL*EPS(ISB,JSB)
        A6N        = -4.D0*EPSI*SS6
        B12N       =  4.D0*EPSI*SS6**2
        QIJN       =  ESCL*CHARGE(ISB)*CHARGE(JSB)*COULF
        DX           = SX(ISP)-SX(JSP)
        DY           = SY(ISP)-SY(JSP)
        DZ           = SZ(ISP)-SZ(JSP)
*                  call PBC(DX,DY,DZ)
        if(DX.gt. HBOXL)DX=DX-BOXL
        if(DX.lt.-HBOXL)DX=DX+BOXL
        if(DY.gt. HBOYL)DY=DY-BOYL
        if(DY.lt.-HBOYL)DY=DY+BOYL
        if(DZ.gt. HBOZL)DZ=DZ-BOZL
        if(DZ.lt.-HBOZL)DZ=DZ+BOZL
*  
	RR           = DX**2+DY**2+DZ**2
	R2I	=  1.d0/RR
        R1           = DSQRT(RR)
*  Electrostatic interaction
        ALPHAR       = ALPHAD*R1
	  TTT          = 1.D0/(1.D0+B1*ALPHAR)
        EXP2A        = DEXP(-ALPHAR**2)
        ERFCR = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A/R1
        EESO          = QIJO*ERFCR
        EESN          = QIJN*ERFCR
*  LJ interaction
        R6I          = R2I*R2I*R2I
        ESRO          = (     A6O+      B12O*R6I)*R6I 
        ESRN          = (     A6N+      B12N*R6I)*R6I
        EEN=EEN+EESN+ESRN
	EEO=EEO+EESO+ESRO
	ELJN=ELJN+ESRN
	ESTN=ESTN+EESN
	ELJO=ELJO+ESRO
	ESTO=ESTO+EESO
	end if !  IMOL.ne.JMOL
*	if(IPRINT.ge.7.and.ESCL.ne.1.d0)write(*,'(2i5,2f12.5,2e14.6)')
*     +ISP,JSP,ENERF*0.001*(ESRO+EESO),ENERF*0.001*(ESRN+EESN),A6O,QIJO
*
      END DO! OF INP					
*    Local forces 
	DO INP      = 1,MBSH
	  ISP=iabs(NBS1(INP))
	  JSP=NBS2(INP)
	  IMOL=NNUM(ISP)
	  JMOL=NNUM(JSP)
	  if(IMOL.ne.JMOL)then
	  ITYP=ITYPE(ISP)
	  JTYP=ITYPE(JSP)
	  ISB =NSITE(ISP)
  	  JSB =NSITE(JSP)
  	  A6O =SIX(ISB,JSB)
	  B12O=TWL(ISB,JSB)
	  QIJO=ONE(ISB,JSB)
	  if(IEE(ITYP).eq.0)then
	    ESCI=EC(ME1) 
	  else
	    ESCI=1.d0
	  end if
	  if(IEE(JTYP).eq.0)then
	    ESCJ=EC(ME1) 
	  else
	    ESCJ=1.d0
	  end if
	  ESCL=ESCI*ESCJ
	  SS6	= SIG(ISB,JSB)**6
	  EPSI	= ESCL*EPS(ISB,JSB)
        A6N        = -4.D0*EPSI*SS6
        B12N       =  4.D0*EPSI*SS6**2
        QIJN       =  ESCL*CHARGE(ISB)*CHARGE(JSB)*COULF
        DX           = SX(ISP)-SX(JSP)
        DY           = SY(ISP)-SY(JSP)
        DZ           = SZ(ISP)-SZ(JSP)
*                  call PBC(DX,DY,DZ)
        if(DX.gt. HBOXL)DX=DX-BOXL
        if(DX.lt.-HBOXL)DX=DX+BOXL
        if(DY.gt. HBOYL)DY=DY-BOYL
        if(DY.lt.-HBOYL)DY=DY+BOYL
        if(DZ.gt. HBOZL)DZ=DZ-BOZL
        if(DZ.lt.-HBOZL)DZ=DZ+BOZL
*  
	  RR           = DX**2+DY**2+DZ**2
	  R2I	= 1.d0/RR
        R1           = DSQRT(RR)
*  Electrostatic interaction
        ALPHAR       = ALPHAD*R1
	  TTT          = 1.D0/(1.D0+B1*ALPHAR)
        EXP2A        = DEXP(-ALPHAR**2)
        ERFCR = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A/R1
        EESO          = QIJO*ERFCR
        EESN          = QIJN*ERFCR
*  LJ interaction
        R6I          = R2I*R2I*R2I
        ESRO          = (     A6O+      B12O*R6I)*R6I 
        ESRN          = (     A6N+      B12N*R6I)*R6I
        EEN=EEN+EESN+ESRN
	  EEO=EEO+EESO+ESRO 
	ELJN=ELJN+ESRN
	ESTN=ESTN+EESN
	ELJO=ELJO+ESRO
	ESTO=ESTO+EESO
	end if   ! of IMOL.ne.JMOL
*
      END DO! OF INP	     
	FACTR=  0.001*ENERF/FNOP
	if(IPRINT.ge.7)write(*,'(a,2f12.5,a,2f12.5)')
     +' LJ:',FACTR*ELJN,FACTR*ELJO,' El:',FACTR*ESTN,
     +FACTR*ESTO
	return
	end 
*
