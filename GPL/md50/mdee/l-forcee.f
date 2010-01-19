*=============== SLOW FORCES =============================================
*
      SUBROUTINE FORCES
*
*  LJ and Re-part of electrostatic interactions
*
      include "mdee.h"
*
      timel0	= cputime(0.d0)
      DO N         = 1,NSTOT
        FX(N)        = 0.D0
        FY(N)        = 0.D0
        FZ(N)        = 0.D0
        OX(N)        = 0.D0
        OY(N)        = 0.D0
        OZ(N)        = 0.D0
      END DO! OF N 
      VIR1=0.
*
      IF(NTYPES.EQ.1.AND.NSPEC(1).EQ.1) RETURN
*
*
* Sum over atom pairs
      DO INP      = 1,MBLN
 10	ISP0=NBL1(INP)
	JSP=NBL2(INP)
        ISP=iabs(ISP0)
	ITYP=ITYPE(ISP)
	JTYP=ITYPE(JSP)
	ISB =NSITE(ISP)
  	JSB =NSITE(JSP)
        IMOL=NNUM(ISP)
        JMOL=NNUM(JSP)
        if(IMOL.eq.JMOL)then
          EPSI       =   DSQRT(EPSIL(ISB))*DSQRT(EPSIL(JSB))
     X        *1.D3/AVSNO/UNITE
          SIGM       =  (SIGMA(ISB)+SIGMA(JSB))/2.d0
          A6         = -4.D0*EPSI*SIGM**6
          B12        =  4.D0*EPSI*SIGM**12
          QIJ        =  CHARGE(ISB)*CHARGE(JSB)*COULF
        else
  	  A6=SIX(ISB,JSB)
	  B12=TWL(ISB,JSB)
	  QIJ=ONE(ISB,JSB)
        end if
	MTR=MDX(ITYP,JTYP)
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
        R1           = DSQRT(RR)
C        if(R1.le.3..and.(ISP.eq.12.or.JSP.eq.12))
C     +write(*,*)'F- long ',R1,' atoms ',ISP,JSP
        R2I          = 1.0D0/RR
	R1I	   = R1*R2I
*  Electrostatic interaction
        ALPHAR       = ALPHAD*R1
	TTT          = 1.D0/(1.D0+B1*ALPHAR)
        EXP2A        = DEXP(-ALPHAR**2)
        ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A 
	EE1	= QIJ*R1I
*  LJ interaction
        R6I          = R2I*R2I*R2I
        ESR          = (     A6+      B12*R6I)*R6I
        FSR          = (6.D0*A6+12.D0*B12*R6I)*R6I*R2I
	if(IMOL.eq.JMOL)then 
C   This can not happen for different molecules... 
C   Scaling constant from the input 
          if(ISP0.le.0)then
	    ESR=ESR*C14LJ(ITYP)
            FSR=FSR*C14LJ(ITYP)
	    EE1=EE1*C14EL(ITYP)
	    EES=EE1
            FES=EES*R2I
          else
            FES          = QIJ*(TOTPI*ALPHAR*EXP2A+ERFC)*R1I*R2I
            EES          = EE1*ERFC
          end if
        else   ! normal case
          FES          = QIJ*(TOTPI*ALPHAR*EXP2A+ERFC)*R1I*R2I
          EES          = EE1*ERFC
	end if
        FORCE        = FES+FSR
*  Forces
          FX(ISP)      = FX(ISP)+FORCE*DX
          FY(ISP)      = FY(ISP)+FORCE*DY
          FZ(ISP)      = FZ(ISP)+FORCE*DZ
          FX(JSP)      = FX(JSP)-FORCE*DX
          FY(JSP)      = FY(JSP)-FORCE*DY
          FZ(JSP)      = FZ(JSP)-FORCE*DZ
          VIR1	 = VIR1+FSR*RR
	  POTES(MTR)   = POTES(MTR)+EE1
	  POTLJ(MTR)   = POTLJ(MTR)+ESR
	  PEST=PEST+EES
	  PE=PE+EES+ESR
*
      END DO! OF INP					
*
*     print *,'potes(1)',potes(1)*enerf*1.d-3/256.
*     print *,'potlj(1)',potlj(1)*enerf*1.d-3/256.
*     print *,'potes(2)',potes(2)*enerf*1.d-3/256.
*     print *,'potlj(2)',potlj(2)*enerf*1.d-3/256.
*     print *,'potes(3)',potes(3)*enerf*1.d-3/256.
*     print *,'potlj(3)',potlj(3)*enerf*1.d-3/256.
*
*Accumulate total forces
*
      DO N         = 1,NSTOT
	M	= NNUM(N)
        GX(N)        = GX(N)+FX(N)
        GY(N)        = GY(N)+FY(N)
        GZ(N)        = GZ(N)+FZ(N)
	  WIRS	= WIRS-FX(N)*(SX(N)-X(M))-FY(N)*(SY(N)-Y(M))-
     +               FZ(N)*(SZ(N)-Z(M))
      END DO! OF N
      VIR=VIR+VIR1
      VIRL=VIRL+VIR1
      timel	= timel+cputime(timel0)
*
      RETURN
      END
*
*======================= LOCALF ==========================================
*
      SUBROUTINE LOCALF
*
*  LJ and electrostatic interactions of closest neigbours
*
      include "mdee.h"
*
      times0	= cputime(0.d0)
      DO N         = 1,NSTOT
        FX(N)        = 0.D0
        FY(N)        = 0.D0
        FZ(N)        = 0.D0
      END DO! OF N
      VIR1=0.  
*
      IF(NTYPES.EQ.1.AND.NSPEC(1).EQ.1) RETURN
*
* Sum over atom pairs
      EEE=0.
      ELJ=0.
      DO INP      = 1,MBSH
	ISP0=NBS1(INP)
	JSP=NBS2(INP)
        ISP=iabs(ISP0)
	ITYP=ITYPE(ISP)
	JTYP=ITYPE(JSP)
	ISB =NSITE(ISP)
  	JSB =NSITE(JSP)
        IMOL=NNUM(ISP)
        JMOL=NNUM(JSP)
        if(IMOL.eq.JMOL)then
          EPSI       =   DSQRT(EPSIL(ISB))*DSQRT(EPSIL(JSB))
     X        *1.D3/AVSNO/UNITE
          SIGM       =  (SIGMA(ISB)+SIGMA(JSB))/2.d0
          A6         = -4.D0*EPSI*SIGM**6
          B12        =  4.D0*EPSI*SIGM**12
          QIJ        =  CHARGE(ISB)*CHARGE(JSB)*COULF
        else
  	  A6=SIX(ISB,JSB)
	  B12=TWL(ISB,JSB)
	  QIJ=ONE(ISB,JSB)
        end if
	MTR=MDX(ITYP,JTYP)
*
*   9.2 Calculate interaction parameters for the given atom pair 
*   ------------------------------------------------------------
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
        R1           = DSQRT(RR)
C        if(R1.le.3..and.(ISP.eq.12.or.JSP.eq.12))
C     +write(*,*)'F- short ',R1,' atoms ',ISP,JSP
          R2I          = 1.0D0/RR
	  R1I	   = R1*R2I
*  Electrostatic interaction
          ALPHAR       = ALPHAD*R1
	  TTT          = 1.D0/(1.D0+B1*ALPHAR)
          EXP2A        = DEXP(-ALPHAR**2)
          ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A 
	  EE1	= QIJ*R1I
	  if(ISP0.le.0)then 
C   This can not happen for different molecules... 
C   Scaling constant from the input 
	    EE1=EE1*C14EL(ITYP)
            EES=EE1
            FORCE=EES*R2I
          else
            FORCE        = QIJ*(TOTPI*ALPHAR*EXP2A+ERFC)*R1I*R2I
            EES          = EE1*ERFC
	  end if
	  POTES(MTR)    = POTES(MTR)+EE1
	  PELS=PELS+EES
          PE=PE+EES
*  LJ interaction
	  if(A6.ne.0.d0.or.B12.ne.0.d0)then
            R6I          = R2I*R2I*R2I
            ESR          = (     A6+      B12*R6I)*R6I
            FSR          = (6.D0*A6+12.D0*B12*R6I)*R6I*R2I
	    if(ISP0.le.0)then 
C   This can not happen for different molecules... 
	      if(IMOL.ne.JMOL)write(*,*)
     +'!!! different mol. types for 1-4 interaction:',
     +ITYP,JTYP,ISP,JSP,' SF'
C   Scaling constant from the input 
	        ESR=ESR*C14LJ(ITYP)
                FSR=FSR*C14LJ(ITYP)
	      end if
	      FORCE        = FORCE+FSR
              VIR1		= VIR1+FSR*RR
	      POTLJ(MTR)    = POTLJ(MTR)+ESR
	      PE		= PE+ESR
	     end if     
*        write(*,'(4i5,5e11.4)')ISP0,JSP,ISB,JSB,R1,ESR,EES,PE
*  Forces
        FX(ISP)      = FX(ISP)+FORCE*DX
        FY(ISP)      = FY(ISP)+FORCE*DY
        FZ(ISP)      = FZ(ISP)+FORCE*DZ
        FX(JSP)      = FX(JSP)-FORCE*DX
        FY(JSP)      = FY(JSP)-FORCE*DY
        FZ(JSP)      = FZ(JSP)-FORCE*DZ
*
      END DO! OF INP
      PERMOL=0.001*ENERF / NOP
*
*Accumulate total forces
*
      DO N         = 1,NSTOT
	M	= NNUM(N)
        HX(N)        = HX(N)+FX(N)
        HY(N)        = HY(N)+FY(N)
        HZ(N)        = HZ(N)+FZ(N)
        WIRSS	= WIRSS-FX(N)*(SX(N)-X(M))-FY(N)*(SY(N)-Y(M))-
     +               FZ(N)*(SZ(N)-Z(M))
      END DO! OF N
      VIRS=VIRS+VIR1
      VIRLS=VIRLS+VIR1
      times	= times+cputime(times0)
*
      RETURN
      END 
*
*=============== FURIR =======================================
*
      SUBROUTINE FURIR
*
      include "mdee.h"
      real*8 TSIN(NTOT),TCOS(NTOT)
*                                       
	timef0	= cputime(0.d0)
*
      call INIFUR
*
      IF(NTYPES.EQ.1.AND.NSPEC(1).EQ.1) RETURN
	if(NKV.eq.0)return
*
      QSS         = 0.D0
      DO ITYP     = 1,NTYPES
        QSS         = QSS+QSUM(ITYP)
      END DO! OF ITYP
      IF(QSS.EQ.0.D0) RETURN
*
      QPE         = 0.D0
      WIR1        = 0.D0
      CFUR=4.d0*PI*COULF/VOL     
*
      DO I        = 1,NSTOT
	  FX(I)=0.d0
	  FY(I)=0.d0
	  FZ(I)=0.d0
      END DO! OF I
*   reciprocal space Evald
      do IKV=1,NKV
        RXV=RX(IKV)
        RYV=RY(IKV)
        RZV=RZ(IKV)
        SCS=0.
        SSN=0.
C$DIR NO_RECURRENCE
        do I=1,NSTOT
          SCP=RXV*SX(I)+RYV*SY(I)+RZV*SZ(I)
          ITYP=ITYPE(I)
          IC=NSITE(I)
	  if(IEE(ITYP).eq.0)then
	      QO=CHARGE(IC)*EC(ME)
	  else
	      QO=CHARGE(IC)
	  end if
*  sin and cos from tables
          SINSC=QO*sin(SCP)
          COSSC=QO*cos(SCP)
          TSIN(I)=SINSC
          TCOS(I)=COSSC
          SSN=SSN+SINSC
          SCS=SCS+COSSC
        end do
        AKC=CFUR*AK(IKV)
        QPE=QPE+AKC*(SSN**2+SCS**2)
        AK2=2.d0*AKC
        do I=1,NSTOT
          QFOR=AK2*(SCS*TSIN(I)-SSN*TCOS(I))
          FX(I)=FX(I)+RXV*QFOR
          FY(I)=FY(I)+RYV*QFOR
          FZ(I)=FZ(I)+RZV*QFOR
        end do
      end do      
      DO I     = 1,NSTOT
        GX(I)    = GX(I)+FX(I)
        GY(I)    = GY(I)+FY(I)
        GZ(I)    = GZ(I)+FZ(I)
        IMOL	= NNUM(I)
        WIR1	= WIR1-FX(I)*(SX(I)-X(IMOL))-FY(I)*(SY(I)-Y(IMOL))-
     +               FZ(I)*(SZ(I)-Z(IMOL))
      END DO! OF I
*
      PE        = PE +QPE
      PEST	= PEST+QPE
      PQE       = QPE         
      WIRF	= WIRF+WIR1
      FACTR=  0.001*ENERF/FNOP
*	if(IPRINT.ge.7)write(*,'(a,4f12.5)')
*     +' FUR_TRUE:',FACTR*QPE
      timef	= timef+cputime(timef0)
      return
      end
*
*================= INIFUR =========================================
*
      subroutine INIFUR
	include "mdee.h"
      data INIT/0/
      save INIT
      if(INIT.ne.0)return 
        INIT=1
	if(ALPHAD.eq.0.d0)then
	  NKV=0
	  return
      end if    
      ALP2=ALPHAD**2
      CUTK2=4.d0*ALP2*FEXP
	CUTK=sqrt(CUTK2)
      FKX=2.d0*PI/BOXL
      FKY=2.d0*PI/BOYL
      FKZ=2.d0*PI/BOZL
      IKV=0
	KMAX=CUTK/FKX+1
	KMAY=CUTK/FKY+1
	KMAZ=CUTK/FKZ+1
      do IKZ=-KMAZ,KMAZ
        do IKX=0,KMAX
          if(IKX.eq.0)then
            if(IKZ.gt.0)then
              KY1=0
            else
              KY1=1
            end if
            KY2=KMAY
          else
            KY1=-KMAY
            KY2=KMAY
          end if
          do IKY=KY1,KY2
            RK2=(IKX*FKX)**2+(IKY*FKY)**2+(IKZ*FKZ)**2
            if(RK2.le.CUTK2)then
      	      IKV=IKV+1
              if(IKV.gt.NKVMAX)stop 'increase NKVMAX'
              RX(IKV)=IKX*FKX
              RY(IKV)=IKY*FKY
              RZ(IKV)=IKZ*FKZ
              AK(IKV)=exp(-0.25d0*RK2/ALP2)/RK2
            end if
          end do
        end do
      end do      
      NKV=IKV
      write(*,*)'Ewald sum:  Kmax=',KMAX,'   Num. of k-vectors ',NKV
      return	
      end 
*   12.  Calculation of self-interaction term in the Ewald sum
*   ----------------------------------------------------------
C
C   Really this procedure calculates two terms:
C   1) Standard self-interaction Ewald term, giving only constant 
C   contribution to the energy
C   2) "Molecular correction" which is due to the fact, that bound
C   atoms in a molecule do not interact electrostatically, but 
C   reciprocal Ewald term includes these interactions also.
*
*============= ETERMS =================================================
*
      SUBROUTINE ETERMS
*  12.1 Front matter
      include "mdee.h"         
      real*8 SPEEI(NTPS)
      data INIT/0/  
      save SPEEI,SPEI,ABC,INIT
      timee0	= cputime(0.d0)
      if(LMNBL)return
      WIR1=0.
*
*  12.2  Electrostatic self-interaction
*  ------------------------------------
	INIT=1
        ABC = ALPHAD/DSQRT(PI)
	SPEI=0.
        do ITYP=1,NTYPES 
	  SPET=0.
	  do I=1,NSITS(ITYP)
	    IS=ISADR(ITYP)+I
  	    QII=ONE(IS,IS)
            SPET         = SPET-NSPEC(ITYP)*ABC*QII
	  end do         
	  SPET=SPET/NUMTASK
	  SPEEI(ITYP)=SPET
	  SPEI=SPEI+SPET
	end do
      do I=1,NTYPES
	SELFPE(I)=SPEEI(I)
      end do            
      SPE=SPEI
*
*  12.3  Calculate intramolecular correction
*  -----------------------------------------
*  12.3.1 Start cycle over bound atoms 
C  list of bound atoms is prepeared in BNDLST (setup.f)   
      IBEG=NUMTASK-TASKID 
      do ISP=IBEG,NSTOT,NUMTASK
	NSP = NNBB(ISP)
	do J=1,NSP
	  JSP=    INBB(ISP,J)
	  if(JSP.le.0)JSP=-JSP
	  if(ISP.gt.JSP)then  ! avoid double count!
	    IS	= NSITE(ISP)
	    JS	= NSITE(JSP)
	    ITYP  = ITYPE(ISP)
  	    QIJ = ONE(IS,JS)
            DX          = SX(ISP)-SX(JSP)
            DY          = SY(ISP)-SY(JSP)
            DZ          = SZ(ISP)-SZ(JSP)
*  12.3.2  Apply periodic boundary conditions
            if(DX.gt. HBOXL)DX=DX-BOXL
            if(DX.lt.-HBOXL)DX=DX+BOXL
	    if(LHEX)then
C  hexagonal periodic cell)
	      XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
              DY=DY-BOYL
              DX=DX-HBOXL
	      XY=BOXYC*DX
            end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
              DY=DY+BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.gt.BOXY3+XY)then
              DY=DY-BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3+XY)then
              DY=DY+BOYL
              DX=DX-HBOXL
            end if
	    else                   
C  rectangular cell
              if(DY.gt. HBOYL)DY=DY-BOYL
              if(DY.lt.-HBOYL)DY=DY+BOYL
	    end if
            if(DZ.gt. HBOZL)DZ=DZ-BOZL
            if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	    if(LOCT)then
	      CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	      DX=DX-sign(CORSS,DX)
	      DY=DY-sign(CORSS,DY)
	      DZ=DZ-sign(CORSS,DZ)
	    end if
*	  call PBC(DX,DY,DZ)
*   12.3.3 Calculate contributions to energy and forces
            R2          = DX**2+DY**2+DZ**2
            R2I         = 1.D0/R2
            R1          = DSQRT(R2)
            R1I         = R1*R2I
            ALPHAR      = ALPHAD*R1
*  	          call RERFC(ALPHAR,ERFC,EXP2A)
            TTT         = 1.D0/(1.D0+B1*ALPHAR)
            EXP2A       = DEXP(-ALPHAR**2)
            ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A
            EES  = QIJ*(ERFC-1.d0)*R1I
            SPE         = SPE+EES
	    SELFPE(ITYP)  = SELFPE(ITYP)+EES
              FES  = QIJ*((TOTPI*ALPHAR*EXP2A+ERFC)-1.d0)*R1I*R2I
              GX(ISP)     = GX(ISP)+FES*DX
              GY(ISP)     = GY(ISP)+FES*DY
              GZ(ISP)     = GZ(ISP)+FES*DZ
              GX(JSP)     = GX(JSP)-FES*DX
              GY(JSP)     = GY(JSP)-FES*DY
              GZ(JSP)     = GZ(JSP)-FES*DZ
              IMOL	= NNUM(ISP)
        WIR1	= WIR1-FES*(DX*(SX(ISP)-X(IMOL))+DY*(SY(ISP)-Y(IMOL))+
     +               DZ*(SZ(ISP)-Z(IMOL)))
              JMOL	= NNUM(JSP)
        WIR1	= WIR1+FES*(DX*(SX(JSP)-X(JMOL))+DY*(SY(JSP)-Y(JMOL))+
     +               DZ*(SZ(JSP)-Z(JMOL)))
	  end if
	end do 
      END DO! OF I
      PE = PE+SPE
      PEST=PEST+SPE
      WIRF=WIRF+WIR1
 100  timee	= timee+cputime(timee0)
      RETURN
      END
*                 
*============= ELRCLJ =================================================
*
      SUBROUTINE ELRCLJ
*
*   long-range corrections to LJ energy and pressure
*
      include"mdee.h"
      data INIT/0/
*
      DO I        = 1,MOLINT
      PELRC(I)    = 0.D0
      VRLRC(I)    = 0.D0
      END DO! OF I
*
      RC          = DSQRT(RSSQ)
      DO 400 ITYP = 1    ,NTYPES
      ISBEG       = ISADR(ITYP)+1
      ISEND       = ISADR(ITYP +1)
      NSPI        = NSPEC(ITYP)
      DO 400   IS = ISBEG,ISEND
      DO 300 JTYP = ITYP ,NTYPES
      MT          = MDX  (ITYP,JTYP)
      JSBEG       = ISADR(JTYP)+1
      JSEND       = ISADR(JTYP +1)
      NSPJ        = NSPEC(JTYP)
      FNOPIJ      = DFLOAT(NSPI*NSPJ)
      DO 300   JS = JSBEG,JSEND
	  EP          = dsqrt(EPSIL(IS)*EPSIL(JS))*1.d3/AVSNO/UNITE
        SI          = 0.5*(SIGMA(IS)+SIGMA(JS))   
        PELRC(MT)   = PELRC(MT)+EP*SI**3*(( 1.D0/9.D0)*(SI/RC)**9
     X                                - ( 1.D0/3.D0)*(SI/RC)**3)*FNOPIJ
        VRLRC(MT)   = VRLRC(MT)+EP*SI**3*((-4.D0/3.D0)*(SI/RC)**9
     X                                + (      2.D0)*(SI/RC)**3)*FNOPIJ
  300   CONTINUE
  400 CONTINUE
*
      FAC       = 1.d-3*ENERF/DFLOAT(NOP)
      FCV	= UNITP/(3.*VOL)
      PELR	= 0.d0
      VIRD	= 0.d0     
      DO I      = 1,MOLINT
        PELRC(I)  = PELRC(I)/VOL
        VRLRC(I)  = -VRLRC(I)/VOL
        PELR=PELR+PELRC(I)
        VIRD=VIRD+VRLRC(I)    
        if(INIT.eq.0.and.IPRINT.ge.6)then
        PRINT "('*** PE  LRC ',I2,': ',F12.8,' kJ/mol')",I,PELRC(I)*FAC
        PRINT "('*** PR  LRC ',I2,': ',F12.4,' atm   ')",I,VRLRC(I)*FCV
        end if
      END DO
      INIT=1
*
      RETURN
      END
*
C===========================================================================
*
*  15. Recalculation of list of neighbours
*  ---------------------------------------
C  This procedure calculate lists of closest and far neighbours
C  for subroutines LOCALF and FORCES
C  These lists must be recalculated after some time
C  Each node calculate goes through its own list of atom pairs 
C  and creates own lists, which are used in subroutines for force calculation
C
C================ CHCNB =======================================================
C     
      subroutine CHCNB(LRDF)
*  15.1 Front matter     
      include "mdee.h"
      logical LGRC,LRDF 
      parameter (RMX2=5.**2,MMOL=5)
      data INIT/0/ 
      timev0	= cputime(0.d0)
      MBLN=0
      MBSH=0
      IOVER=0
      IESH=0
      LGRC=LGR.and.LRDF
      RDFCUT2=RDFCUT**2
*
*  15.2  Organize cycle over ALL atom pairs
*  ----------------------------------------
      IBEG=NUMTASK-TASKID 
C  cycle over I,J is organized is follows:
C  it takes I>J if I and J are of the same parity and I<J if I and J have
C  different parity. This allows nearly uniformly distribute atom pairs 
C  between the nodes
      do I=1,NSTOT
	IMOL=NNUM(I) 
	ITYP=ITYPE(I)
*   15.2.1 Even strips
*        do J=I-2,1,-2
        do J=1,I-1
	  JMOL=NNUM(J)
	  DX=SX(I)-SX(J)
	  DY=SY(I)-SY(J)
	  DZ=SZ(I)-SZ(J)
*   15.2.1.1 Periodic boundary conditions
C	    call PBC (DX,DY,DZ)
	  if(DX.gt. HBOXL)DX=DX-BOXL
          if(DX.lt.-HBOXL)DX=DX+BOXL
          if(DZ.gt. HBOZL)DZ=DZ-BOZL
          if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	    if(LHEX)then
C  hexagonal periodic cell)
	      XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
              DY=DY-BOYL
              DX=DX-HBOXL
	      XY=BOXYC*DX
            end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
              DY=DY+BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.gt.BOXY3+XY)then
              DY=DY-BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3+XY)then
              DY=DY+BOYL
              DX=DX-HBOXL
            end if
	    else                   
C  rectangular cell
            if(DY.gt. HBOYL)DY=DY-BOYL
            if(DY.lt.-HBOYL)DY=DY+BOYL
	    end if
	    if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	      CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	      DX=DX-sign(CORSS,DX)
	      DY=DY-sign(CORSS,DY)
	      DZ=DZ-sign(CORSS,DZ)
	    end if                  
	    RR2=DX**2+DY**2+DZ**2  
*  15.2.1.2 Case of "molecular" cutoff - calculate COM distance
C  If we do not use Ewald summation (LMNBL=.true.), 
C  it is important that at least
C  small molecules (less than MMOL atoms) where inside or outside
C  cut-off radius as a whole 
	    if(LMNBL)then
	      NMI=NSITS(ITYP) 
	      NMJ=NSITS(ITYPE(J))
	      if(NMI.gt.MMOL.and.NMJ.gt.MMOL)then
		RRM=RR2
		go to 110
	      else if(NMI.gt.MMOL.and.NMJ.le.MMOL)then
	        DX=SX(I)-X(JMOL)
	        DY=SY(I)-Y(JMOL)
	        DZ=SZ(I)-Z(JMOL)
	      else if(NMJ.gt.MMOL)then
	        DX=X(IMOL)-SX(J)
	        DY=Y(IMOL)-SY(J)
	        DZ=Z(IMOL)-SZ(J)
	      else 
	        DX=X(IMOL)-X(JMOL)
	        DY=Y(IMOL)-Y(JMOL)
	        DZ=Z(IMOL)-Z(JMOL)
	      end if
*	      call PBC (DX,DY,DZ)
	      if(DX.gt. HBOXL)DX=DX-BOXL
              if(DX.lt.-HBOXL)DX=DX+BOXL
              if(DZ.gt. HBOZL)DZ=DZ-BOZL
              if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	      if(LHEX)then
C  hexagonal periodic cell)
	        XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
                DY=DY-BOYL
                DX=DX-HBOXL
	        XY=BOXYC*DX
              end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
                DY=DY+BOYL
                DX=DX+HBOXL
	        XY=BOXYC*DX
              end if
              if(DY.gt.BOXY3+XY)then
                DY=DY-BOYL
                DX=DX+HBOXL
	        XY=BOXYC*DX
              end if
              if(DY.lt.-BOXY3+XY)then
                DY=DY+BOYL
                DX=DX-HBOXL
              end if
	      else                   
C  rectangular cell
                if(DY.gt. HBOYL)DY=DY-BOYL
                if(DY.lt.-HBOYL)DY=DY+BOYL
	      end if
	      if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	        CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	        DX=DX-sign(CORSS,DX)
	        DY=DY-sign(CORSS,DY)
	        DZ=DZ-sign(CORSS,DZ)
	      end if                  
C  square distance betwen COM of small molecules
	      RRM=DX**2+DY**2+DZ**2 
	    else   ! .not.LMNBL
	      RRM=RR2              
	    end if
*  15.2.1.3 Check if the atoms are bound
 110	    if(IMOL.eq.JMOL)then    ! bound atoms should not be in the list
	      if(RR2.le.RMX2)then
C   Check that the atoms are non-bound. Check using binding list of the
C   atom which has less bound atoms
	        if(NNBB(I).le.NNBB(J))then
		  do K=1,NNBB(I) 
C   1-3 bound:  not in the list
		    if(INBB(I,K).eq.J)go to 125
                    if(INBB(I,K).eq.-J.and.L14NB(ITYP))then
C   1-4 bound:  put into list with negative first index
C   (these atom pairs may have a special way of force calculations)
  	              MBSH=MBSH+1
	              if(MBSH.gt.NBSMAX)IESH=1
                      if(IESH.eq.0)then
	                NBS1(MBSH)=-I
		        NBS2(MBSH)=J
                      end if
*            write(*,'(a,2i5,f12.5,2i8)')'1-4s',I,J,sqrt(RR2),MBSH,MBLN
	              go to 125
		    end if  
	  	  end do
	        else
		  do K=1,NNBB(J) 
C   1-3 bound
		    if(INBB(J,K).eq.I)go to 125
                    if(INBB(J,K).eq.-I.and.L14NB(ITYP))then
C   1-4 bound                                     
	              MBSH=MBSH+1
	              if(MBSH.gt.NBSMAX)IESH=1
                      if(IESH.eq.0)then
	                NBS1(MBSH)=-I
		        NBS2(MBSH)=J
                      end if
*            write(*,'(a,2i5,f12.5,2i8)')'1-4s',I,J,sqrt(RR2),MBSH,MBLN
                      go to 125
		    end if  
	  	  end do
	        end if
	      end if ! (RR2.le.RMAX2
	    end if   ! (IMOL.eq.JMOL
*   15.2.1.4 Pick up non-bound atoms to the list of neigbours
	    if(RRM.le.SHORT2)then
	      MBSH=MBSH+1
	      if(MBSH.gt.NBSMAX)IESH=1
              if(IESH.eq.0)then
	        NBS1(MBSH)=I
	        NBS2(MBSH)=J
              end if
*            write(*,'(a,2i5,f12.5,2i8)')'no_s',I,J,sqrt(RR2),MBSH,MBLN
	    else if(RRM.le.RSSQ)then
	      MBLN=MBLN+1
	      if(MBLN.gt.NBLMAX)then
                if(IOVER.eq.0)write(*,*)
     +     '!! list of far neigbours overfilled'
                IOVER=1
              end if
              if(IOVER.eq.0)then
	        NBL1(MBLN)=I
	        NBL2(MBLN)=J
*            write(*,'(a,2i5,f12.5,2i8)')'no_l',I,J,sqrt(RR2),MBSH,MBLN
              end if
	    end if
*   15.2.1.5 RDF calculation
C  omitted            
C   Only at this point we go through all the atom pairs
C   RDF calculation
 125        continue
          IF(LGRC.and.RR2.le.RDFCUT2.and.ME.eq.1.and.NHIST.ge.IHIST)THEN	!=====> LGR
	      ISB	= NSITE(J)
	      JSB	= NSITE(I)
*              write(*,*)ISB,JSB,RR2
	      if(NGRI(ISB).le.NGRI(JSB))then
		do JS=1,NGRI(ISB) 
		  if(IGRI(ISB,JS).eq.JSB)then 
		    IP = MGRI(ISB,JS)
	            go to 137  
		  end if
		end do
	      else
		do JS=1,NGRI(JSB)
		  if(IGRI(JSB,JS).eq.ISB)then
		    IP = MGRI(JSB,JS)
	            go to 137
		  end if  
		end do         
	      end if  
	      go to 138
 137	      R1	= sqrt(RR2)
              MM        = R1*DRDFI+1
            if(MM.le.MAX)IRDF(MM,IP) = IRDF(MM,IP)+1
*            write(*,'(2i4,f10.4,4i5)')I,J,R1,MM,IP,NRDF,IRDF(MM,IP)
	    end if ! of LGRC
 138        continue
        end do
        end do
        if(ME.eq.1.and.NHIST.ge.IHIST)NRDF=NRDF+1
        if(IOVER.ne.0)then
          write(*,*)' Insrease NBLMAX in prcm.h for this job to ',MBLN
          I=MBLN/NTOT+1
          write(*,*)' (change  to:  NBLMAX =',I,'* NTOT)'
          stop
        end if
        timev	= timev+cputime(timev0)
	if(IPRINT.lt.7)return
	write(*,*)' Task ',TASKID,
     +'total num.of interactions: ',MBSH,' and ',MBLN
        if(IESH.ne.0)then
          write(*,*)'!!! list of closest neigbours overfilled ',MBSH
          stop
        end if
	return
	end
