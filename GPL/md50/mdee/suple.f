*================ CUBIC =================================================
*
      SUBROUTINE CUBIC(X0,Y0,Z0,X,Y,Z,BOXL,NOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION X0(*),Y0(*),Z0(*),X(*),Y(*),Z(*)
*
      FNOP     = DFLOAT(NOP)
      NR       = AINT(FNOP**(1.D0/3.D0) +.99999D0)
      NRCUBE   = NR**3
      NEMPTY   = NRCUBE-NOP
      SHIFT    = BOXL/DFLOAT(NR)
      SKIFT    = 0.5*BOXL*(1.D0 + 1.D0/DFLOAT(NR))
      L        = 0
      DO 100 I = 1,NR
      DO 100 J = 1,NR
      DO 100 K = 1,NR
      L        = L+1
      X(L)     = DFLOAT(I)*SHIFT
      Y(L)     = DFLOAT(J)*SHIFT
      Z(L)     = DFLOAT(K)*SHIFT
  100 CONTINUE
      L        = 0
      NSKIP    = 0
      IF(NEMPTY.NE.0) NSKIP=NOP/NEMPTY+1
      PRINT "(/1X,'*** LATTICE: ',I2,' SITES/EDGE  -> ',I4,
     X' SITES TOTALLY ( ',I3,' VACANT  - WITH INTERVALL OF ',I3,' )'/)"
     X,NR,NRCUBE,NEMPTY,NSKIP
*
      IF(NEMPTY.NE.0) THEN
*
      DO 200 I = 1,NRCUBE
      IF(MOD(I,NSKIP).EQ.0) GO TO 200
      L        = L+1
      X0(L)    = X(I)-SKIFT
      Y0(L)    = Y(I)-SKIFT
      Z0(L)    = Z(I)-SKIFT
      IF(L.EQ.NOP) GO TO 300
  200 CONTINUE
*
  300 CONTINUE
*
      ELSE
*
      DO 400 I = 1,NRCUBE
      X0(I)    = X(I)-SKIFT
      Y0(I)    = Y(I)-SKIFT
      Z0(I)    = Z(I)-SKIFT
  400 CONTINUE
*
      END IF
*
      DO 500 I = 1,NOP
        call PBC(X0(I),Y0(I),Z0(I))
  500 CONTINUE
*
      RETURN
      END
*
*============ FMELT ================================================
*
      SUBROUTINE FMELT(X0,Y0,Z0,NSP)
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DIMENSION X0(*),Y0(*),Z0(*)
      DATA PI/3.1415926536D0/
      CN=PI*(2*NSP)**(1.D0/3.D0)
      AB=0.D0
      BC=0.D0
      CD=0.D0
      DE=0.D0
      EF=0.D0
      FG=0.D0
      DO 10 I=1,NSP
      AB=AB+DCOS(CN*X0(I))
      BC=BC+DSIN(CN*X0(I))
      CD=CD+DCOS(CN*Y0(I))
      DE=DE+DSIN(CN*Y0(I))
      EF=EF+DCOS(CN*Z0(I))
      FG=FG+DSIN(CN*Z0(I))
   10 CONTINUE
      FMLT=(AB*AB+BC*BC+CD*CD+DE*DE+EF*EF+FG*FG)/DFLOAT(3*NSP)
      PRINT "(/1X,'  *****   MELTING FACTOR: ',F12.4/)",FMLT
      RETURN
      END
*
*=============== FCC ===================================================
*
      SUBROUTINE FCC(X,Y,Z,BOXL,NSP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NSP),Y(NSP),Z(NSP)
      DIMENSION XB(4),YB(4),ZB(4)
      DATA XB/0.5,-.5,0.5,-.5/,YB/0.5,-.5,-.5,0.5/,ZB/0.5,0.5,-.5,-.5/
C     SET UP STARTING FCC LATTICE 
      HBOXL	= 0.5*BOXL
      NC     = (NSP/4)**(1.D0/3.D0)+0.9999
      DRINK	= 0.01*HBOXL/NC
      FACT   = HBOXL/DFLOAT(NC)
      SHIFT  = BOXL/DFLOAT(NC)
      BASE   = HBOXL*(-1.D0+1./DFLOAT(NC))
      IL     = 0
      DO 20 IB = 1,4
        ZS       = BASE+ZB(IB)*FACT
        DO 19 IZ = 1,NC
          YS       = BASE+YB(IB)*FACT
          DO 18 IY = 1,NC
            XS       = BASE+XB(IB)*FACT
            DO 17 IX = 1,NC
              IL       = IL+1
              X(IL)    = XS
              Y(IL)    = YS
              Z(IL)    = ZS
              if(IL.ge.NSP)go to 24
              XS       = XS+SHIFT
   17       CONTINUE
            YS       = YS+SHIFT
   18     CONTINUE
          ZS       = ZS+SHIFT
   19   CONTINUE
   20 CONTINUE
*
   24 DO I     = 1,NSP
        FACT     = (2.D0*(RANF()-0.5D0))*DRINK
        X(I)     = X(I)+FACT*X(I)
        FACT     = (2.D0*(RANF()-0.5D0))*DRINK
        Y(I)     = Y(I)+FACT*Y(I)
        FACT     = (2.D0*(RANF()-0.5D0))*DRINK
        Z(I)     = Z(I)+FACT*Z(I)
      END DO
*
      do I=1,NSP
        call PBC(X(I),Y(I),Z(I))
      end do
*
      PRINT *,'    FCC LATTICE CREATED FOR ',NSP,' POINTS'
*
      RETURN
      END
*
*=============== RANQ =========================================
*
      SUBROUTINE RANQ(QP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QP(4)
*
      PI=3.14159D0
      A=PI*RANF()/2.D0
      B=PI*RANF()
      C=PI*RANF()
      QP(1)=DSIN(A)*DSIN(B-C)
      QP(2)=DSIN(A)*DCOS(B-C)
      QP(3)=DCOS(A)*DSIN(B+C)
      QP(4)=DCOS(A)*DCOS(B+C)
      RETURN
      END
*
*=============== ROTATE ==============================================
*
      SUBROUTINE ROTATE(Q,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Q(4)
      SIGN = -1.D0
      AA   = Q(1)*Q(1)
      BB   = Q(2)*Q(2)
      CC   = Q(3)*Q(3)
      DD   = Q(4)*Q(4)
      AB   = Q(1)*Q(2)
      AC   = Q(1)*Q(3)
      AD   = Q(1)*Q(4)*SIGN
      BC   = Q(2)*Q(3)
      BD   = Q(2)*Q(4)*SIGN
      CD   = Q(3)*Q(4)*SIGN
      RX   = X*(-AA+BB-CC+DD)+Y*(+CD-AB)*2.D0 +Z*(+BC+AD)*2.D0
      RY   = X*(-AB-CD)*2.D0 +Y*(+AA-BB-CC+DD)+Z*(+BD-AC)*2.D0
      RZ   = X*(+BC-AD)*2.D0 +Y*(-AC-BD)*2.D0 +Z*(-AA-BB+CC+DD)
      X    = RX
      Y    = RY
      Z    = RZ
      RETURN
      END
*
*=============== COMVEL ========================================
*
      SUBROUTINE COMVEL(VX,VY,VZ,TEMP,NTOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VX(*),VY(*),VZ(*)
*
      RTEMP   = DSQRT(TEMP)
      DO N    = 1,NTOT
      VX(N)   = RTEMP*GAUSS(DUMMY)
      VY(N)   = RTEMP*GAUSS(DUMMY)
      VZ(N)   = RTEMP*GAUSS(DUMMY)
      END DO! OF N
*
      SUMX    = 0.D0
      SUMY    = 0.D0
      SUMZ    = 0.D0
      DO N    = 1,NTOT
      SUMX    = SUMX+VX(N)
      SUMY    = SUMY+VY(N)
      SUMZ    = SUMZ+VZ(N)
      END DO! OF N
      SUMX    = SUMX/DFLOAT(NTOT)
      SUMY    = SUMY/DFLOAT(NTOT)
      SUMZ    = SUMZ/DFLOAT(NTOT)
*
      DO N    = 1,NTOT
      VX(N)   = VX(N)-SUMX
      VY(N)   = VY(N)-SUMY
      VZ(N)   = VZ(N)-SUMZ
      END DO! OF N
*
      RETURN
      END
*
*============ GAUSS ===============================================
*
      DOUBLE PRECISION FUNCTION GAUSS(DUMMY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (A1=3.949846138,A3=0.252408784)
      PARAMETER (A5=0.076542912,A7=0.008355968)
      PARAMETER (A9=0.029899776)
      SUM   = 0.0
      DO I  = 1,12
      SUM   = SUM + RANF()
      END DO
      R     = (SUM-6.0)/4.0
      R2    = R*R
      GAUSS = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)
      RETURN
      END
*
*============================== RANF =========================
*
      DOUBLE PRECISION FUNCTION RANF()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER L,C,M,SEED
      PARAMETER (L=1029,C=221591,M=1048576)
      SAVE SEED
      DATA SEED /0/
      SEED    = MOD(SEED*L+C,M)
      RANF    = DFLOAT(SEED)/M
      RETURN	
      END
*
*=============== SCALEV ==========================================
*
      SUBROUTINE SCALEV(VX,VY,VZ,FTMP,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION VX(*),VY(*),VZ(*)
      DO I  = 1,N
      VX(I) = VX(I)*FTMP
      VY(I) = VY(I)*FTMP
      VZ(I) = VZ(I)*FTMP
      END DO
      END
*
*=============== SCLFRC ===================================
*
*    4. Cut large forces
*    -------------------
C
C    This subroutine scale down large forces arising e.g. if
C    initial conformation was too bad. Input parameter is
C    NNN, which says which forces should be scaled:
C    1 - HX (fast forces in double time step algorithm) 
C    2 - GX (slow forces in double time step algorithm) 
C    3 - FX (forces in shake algorithm)
C    Output parameter DONE=.true. if scaling of forces have been done
C    otherwise it returns .false.
C
      SUBROUTINE SCLFRC(DONE,NNN)
      include "mdee.h"
*
      LOGICAL DONE
      DONE      = .FALSE.

      if(.not.LSCFT)return
*
      GO TO (1000,2000,3000),NNN
*
      STOP '!!! WRONG OPTION IN SCLFRC !!!'
*
*  4.1 Cut large fast forces in double time step algorithm
*
 1000 CONTINUE
      DO I      = 1,NSTOT
        FSTERM    = HX(I)**2+HY(I)**2+HZ(I)**2
        FSTERM    = DSQRT(FSTERM)*FTSCL
        IF(FSTERM.GT.1.D0) THEN
          if(IPRINT.ge.8)write(*,*)' Force on atom ',I,' scaled by ',
     +    1./FSTERM
          HX(I)     = HX(I)/FSTERM
          HY(I)     = HY(I)/FSTERM
          HZ(I)     = HZ(I)/FSTERM
          DONE      = .TRUE.
        END IF
      END DO! OF I
      if(DONE)write(*,*)' fast forces have been cut'
      RETURN
*
*  4.2 Cut large slow forces in double time step algorithm
*
 2000 CONTINUE
      DO I      = 1,NSTOT
        FSTERM    = GX(I)**2+GY(I)**2+GZ(I)**2
        FSTERM    = DSQRT(FSTERM)*FTSCL
        IF(FSTERM.GT.1.D0) THEN
          GX(I)     = GX(I)/FSTERM
          GY(I)     = GY(I)/FSTERM
          GZ(I)     = GZ(I)/FSTERM
          if(IPRINT.ge.8)write(*,*)' Force on atom ',I,' scaled by ',
     +    1./FSTERM
          DONE      = .TRUE.
        END IF
      END DO! OF I
      if(DONE)write(*,*)' slow forces have been cut'
      RETURN
*
*  4.3 Cut large forces in double time step algorithm
*
 3000 CONTINUE
      DO I      = 1,NSTOT
        FSTERM    = FX(I)**2+FY(I)**2+FZ(I)**2
        FSTERM    = DSQRT(FSTERM)*FTSCL
        IF(FSTERM.GT.1.D0) THEN
          FX(I)     = FX(I)/FSTERM
          FY(I)     = FY(I)/FSTERM
          FZ(I)     = FZ(I)/FSTERM
          if(IPRINT.ge.8)write(*,*)' Force on atom ',I,' scaled by ',
     +    1./FSTERM
          DONE      = .TRUE.
        END IF
      END DO! OF I
      if(DONE)write(*,*)' forces have been cut'
*
      RETURN
      END
*
*=============== DIPOLE ===========================================
*
      SUBROUTINE DIPOLE
      include "mdee.h"
*
      FAK	  = 1.0/0.52917
      FAC         = 2.5418
      N           = 0
      DO ITYP     = 1,NTYPES
      NSP         = NSPEC(ITYP)
      FNSPI       = 1.D0/DFLOAT(NSP)
      DIPMOM      = 0.D0
      ISB         = ISADR(ITYP)+1
      ISE         = ISADR(ITYP +1)
      IBE         = IADDR(ITYP)+1
      IEN         = IADDR(ITYP +1)
      DO I        = IBE,IEN
      DMX         = 0.D0
      DMY         = 0.D0
      DMZ         = 0.D0
      DO IS       = ISB,ISE
      N           = N+1
Calculate permanent dipole moment (magnitude)
      DMX         = DMX+CHARGE(IS)*WX(N)
      DMY         = DMY+CHARGE(IS)*WY(N)
      DMZ         = DMZ+CHARGE(IS)*WZ(N)
      END DO! OF IS
      DNORM       = DSQRT(DMX**2+DMY**2+DMZ**2)
      DIPMOM      = DIPMOM+DNORM*FAC*FAK*RF
      IF(DNORM.GT.0.D0) DNORM = 1.D0/DNORM
      DPX(I)      = DMX*DNORM
      DPY(I)      = DMY*DNORM
      DPZ(I)      = DMZ*DNORM
      END DO! OF I
      DIPMOM      = DIPMOM*FNSPI
      IF(MOD(NSTEP,10).EQ.0.and.IPRINT.ge.8) THEN
      IF(DABS(DIPMOM).GT.0.00001D0) THEN
      PRINT *,'--> --> --> --> --> --> --> --> --> --> --> --> --> -->'
      PRINT "(' DIPOLE MOMENT(D):',F8.3,' FOR TYPE  ',A6)"
     X,DIPMOM,NAME(ITYP)
      END IF
      END IF
      END DO! OF ITYP
      RETURN
      END
*
*=============== ZEROFS =============================================
*
      SUBROUTINE ZEROFS(NNN)
      include "mdee.h"

*Zero forces:
*
      GO TO (1000,2000),NNN
*
 1000 CONTINUE
      PELS=0.
	VIRS=0.
	VIRB=0.
	VIRO=0.
	VIRT=0.
	VIRA=0.
	WIRSS=0.
	VIRLS=0.
      PINT	= 0. 
      PE         = 0.D0
      DO I       = 1,MOLINT
      POTES(I)   = 0.D0
      POTLJ(I)   = 0.D0
      END DO
	do I=1,NTYPES
	SELFPE(I)      = 0.D0
      PES14 (I)      = 0.D0
      PSR14 (I)      = 0.D0
      end do
      DO I       = 1,NSTOT
      HX(I)      = 0.D0
      HY(I)      = 0.D0
      HZ(I)      = 0.D0
      END DO
      RETURN
*
 2000 CONTINUE
*
	PEST=0.
      VIR=0.
      VIRL=0.
      VIRN=0.
      WIRF=0.
	WIRS=0.
	DO I       = 1,NSTOT
      GX(I)      = 0.D0
      GY(I)      = 0.D0
      GZ(I)      = 0.D0
      END DO
      RETURN
      END
*
*=============== GETCOM ============================================
*
      SUBROUTINE GETCOM
      include "mdee.h"
*
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES
      NSBEG       = ISADR  (ITYP)+1
      NSEND       = ISADR  (ITYP +1)
      SUMM        = SUMMAS (ITYP)
      DO J        = 1,NSPEC(ITYP)
      I           = I+1
      X (I)       = 0.D0
      Y (I)       = 0.D0
      Z (I)       = 0.D0 
	PX (I)       = 0.D0
      PY (I)       = 0.D0
      PZ (I)       = 0.D0
Calculate C.O.M. vectors:
      DO IS       = NSBEG,NSEND
      N           = N+1
      X (I)       = X(I)+MASS(IS)*SX(N)
      Y (I)       = Y(I)+MASS(IS)*SY(N)
      Z (I)       = Z(I)+MASS(IS)*SZ(N)
      PX (I)       = PX(I)+VX(N)
      PY (I)       = PY(I)+VY(N)
      PZ (I)       = PZ(I)+VZ(N)
      END DO! OF IS
      X (I)       = X(I)/SUMM
      Y (I)       = Y(I)/SUMM
      Z (I)       = Z(I)/SUMM
      END DO! OF I
      END DO! OF ITYP
*
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES
      NSBEG       = ISADR  (ITYP)+1
      NSEND       = ISADR  (ITYP +1)
      DO J        = 1,NSPEC(ITYP)
      I           = I+1
      DO IS       = NSBEG,NSEND
      N           = N+1
Calculate short site vectors
      WX(N)       = SX(N)-X(I)
      WY(N)       = SY(N)-Y(I)
      WZ(N)       = SZ(N)-Z(I)
      END DO! OF IS
      END DO! OF I
      END DO! OF ITYP
*Apply periodic boundaries on COM points
	do I=1,NOP
	  call PBC(X(I),Y(I),Z(I))
	end do
c	do ITYP	= 1,NTYPES
c	  if(LNPT.and.IEE(ITYP).eq.0.and.EC(ME).lt.1.d-7)then
c         DO I        = IADDR(ITYP)+1,IADDR(ITYP+1)
c            if(X(I).gt. BOXF)X(I)=X(I)-2.*BOXF
c            if(X(I).lt.-BOXF)X(I)=X(I)+2.*BOXF
c            if(Y(I).gt. BOXF)Y(I)=Y(I)-2.*BOXF
c            if(Y(I).lt.-BOXF)Y(I)=Y(I)+2.*BOXF
c            if(Z(I).gt. BOXF)Z(I)=Z(I)-2.*BOXF
c            if(Z(I).lt.-BOXF)Z(I)=Z(I)+2.*BOXF
c	    end do
c	  else
c          DO I        = IADDR(ITYP)+1,IADDR(ITYP+1)
c            call PBC(X(I),Y(I),Z(I))
c	    end do
c	  end if
c      END DO! OF I
*
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES
      NSBEG       = ISADR  (ITYP)+1
      NSEND       = ISADR  (ITYP +1)
      DO J        = 1,NSPEC(ITYP)
      I           = I+1
      DO IS       = NSBEG,NSEND
      N           = N+1
*Reconstruct long site vectors
      SX(N)       = WX(N)+X(I)
      SY(N)       = WY(N)+Y(I)
      SZ(N)       = WZ(N)+Z(I)
      END DO! OF IS
      END DO! OF I
      END DO! OF ITYP
*
      RETURN
      END
*
*=============== MPOLES =============================================
*
      SUBROUTINE MPOLES(SITE,Q,DM,QM,OM,HM,NSB,NSE,NS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION Q(NS),SITE(NS,3)
     X,DM(3),QM(3,3),OM(3,3,3),HM(3,3,3,3)
*
      FAK		= 1.0/0.52917
      DO 10 IS		= NSB,NSE
      DO 10 J		= 1,3
      DM(J) 		= 0.D0
      DO 10 K		= 1,3
      QM(J,K)		= 0.D0
      DO 10 L		= 1,3
      OM(J,K,L)		= 0.D0
      DO 10 M		= 1,3
      HM(J,K,L,M)	= 0.D0
   10 CONTINUE
*
      DO 20 IS		= NSB,NSE
      DO 20 J		= 1,3
      DM(J) 		= DM(J)+Q(IS)*SITE(IS,J)*FAK
      DO 20 K		= 1,3
      QM(J,K)		= QM(J,K)+Q(IS)*SITE(IS,J)*SITE(IS,K)*FAK**2
      DO 20 L		= 1,3
      OM(J,K,L)		= OM(J,K,L)+Q(IS)*SITE(IS,J)*
     X                    SITE(IS,K)*SITE(IS,L)*FAK**3
      DO 20 M		= 1,3
      HM(J,K,L,M)	= HM(J,K,L,M)+Q(IS)*SITE(IS,J)*
     X                    SITE(IS,K)*SITE(IS,L)*SITE(IS,M)*FAK**4
   20 CONTINUE
      RETURN
      END

*....
*======== XINTGR =======================================================
*....
      FUNCTION XINTGR(H,Y,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(*)
      HT=H/3.0
      IF(NDIM.LT.6)RETURN
      SUM1=4.00*Y(2)
      SUM1=HT*(Y(1)+SUM1+Y(3))
      AUX1=4.00*Y(4)
      AUX1=SUM1+HT*(Y(3)+AUX1+Y(5))
      AUX2=HT*(Y(1)+3.8750*(Y(2)+Y(5))+2.6250*(Y(3)+Y(4))+Y(6))
      SUM2=4.00*Y(5)
      SUM2=AUX2-HT*(Y(4)+SUM2+Y(6))
      AUX=4.00*Y(3)
      IF(NDIM-6)5,5,2
    2 DO 4 I=7,NDIM,2
      SUM1=AUX1
      SUM2=AUX2
      AUX1=4.00*Y(I-1)
      AUX1=SUM1+HT*(Y(I-2)+AUX1+Y(I))
      IF(I-NDIM)3,6,6
    3 AUX2=4.00*Y(I)
      AUX2=SUM2+HT*(Y(I-1)+AUX2+Y(I+1))
    4 CONTINUE
    5 CONTINUE
      XINTGR=AUX2
      RETURN
    6 CONTINUE
      XINTGR=AUX1
      RETURN
      END
*
*=================== CHTERM ===================================
*                                                              
	subroutine CHTERM
	include "mdee.h"
	PTX=0.
	PTY=0.
	PTZ=0.
*  Check exess momenta
	do I=1,NSTOT
	  PTX=PTX+VX(I)
	  PTY=PTY+VY(I)
	  PTZ=PTZ+VZ(I)
	end do	
	call GETEMP
	TEMPO=TEMP
	SUMASD=TOTMAS/UNITM
	DVX=PTX/SUMASD
	DVY=PTY/SUMASD
	DVZ=PTZ/SUMASD
	do I=1,NSTOT
	  VX(I)=VX(I)-DVX/MASSDI(I)
	  VY(I)=VY(I)-DVY/MASSDI(I)
	  VZ(I)=VZ(I)-DVZ/MASSDI(I)
	end do
      call GETEMP
	  if(dabs(TEMPO-TEMP)/TEMP.gt.0.01)then
	  write(*,*)' anomalous excess momenta removed'
	  write(*,*)'old temperature ',TEMPO
        write(*,*)'new temperature ',TEMP
      end if
*  Set new temperature
	if(LCHT)then
	  SCT=sqrt(TRTEMP/TEMP)
	  do I=1,NSTOT
	    VX(I)=VX(I)*SCT
          VY(I)=VY(I)*SCT
	    VZ(I)=VZ(I)*SCT
        end do
	  call GETEMP
	  write(*,*)' New temperature set to',TEMP
	end if
*  Set new volume/density
	if(LCHP)then
	  IF(RHO.gt.1.d-20) THEN
          IF(CLN.le.1.d-20) THEN
	      VOL = TOTMAS*1.d27/RHO
            CL        = VOL**(1.0/3.0)
	      write(*,'(2(a,f10.4))')' New density: ',RHO,' BOXL=',CL
          else
            CL=CLN
	      VOL=CL**3
	      DENS	= TOTMAS*1.d-3/(VOL*UNITL**3)
	      write(*,'(2(a,f10.4))')' New BOX size ',CL,' Density',DENS
          ENDIF
	    SCX=CL/BOXL
        ELSE
*Vacuum simulation - arbitrary box size:
          CL        = 1000.D0
	    SCX=1.
          PRINT "(/,'*** VACUUM SIMULATION ***',/)"
        END IF
*  Scaling coordinates
	  N=0
	  DO ITYP  = 1,NTYPES                    ! over types
          ISB      = ISADR (ITYP)+1
          ISE      = ISADR (ITYP +1)
          JBE      = IADDR (ITYP)+1
          JEN      = IADDR (ITYP +1)
          DO J     = JBE,JEN                  ! over molecules
	      DX=(SCX-1.d0)*X(J)
	      DY=(SCX-1.d0)*Y(J)
	      DZ=(SCX-1.d0)*Z(J)
	      X(J)=X(J)+DX
            Y(J)=Y(J)+DY
            Z(J)=Z(J)+DZ
            DO IS    = ISB,ISE                ! over atoms
              N        = N+1
              SX(N)=SX(N)+DX
              SY(N)=SY(N)+DY
              SZ(N)=SZ(N)+DZ
            END DO! OF IS
          end do ! of J
	  end do ! of ITYP 
	  write(*,*)' Scaling distances by ',SCX
	  call RECLEN(CL)
	end if
	return
	end
*	
*=============== GETKIN ================================================
*
      SUBROUTINE GETKIN
      include "mdee.h"
*
      DIMENSION XX(NPART),XY(NPART),XZ(NPART)
      DIMENSION UX(NTOT),UY(NTOT),UZ(NTOT)
      DIMENSION RMOI(3,3,NPART),RMOII(3,3,NPART)
     X         ,FMOI(9)      ,FMOII(9)
	real*8 TOR(9),TORI(9) 
      DOUBLE PRECISION MASSA
*
      N        = 0
	EKTR	= 0.
	EKROT	= 0.
	TKET	= 0.
      DO ITYP  = 1,NTYPES                    ! over types
	  SUMMD	= SUMMAS(ITYP)/(AVSNO*1000.0*UNITM)
        NSS      = NSITS (ITYP)
        NSP      = NSPEC (ITYP)
        ISB      = ISADR (ITYP)+1
        ISE      = ISADR (ITYP +1)
        JBE      = IADDR (ITYP)+1
        JEN      = IADDR (ITYP +1)  
        DO J     = JBE,JEN                  ! over molecules
          PX(J)    = 0.D0
          PY(J)    = 0.D0
          PZ(J)    = 0.D0
          QX(J)    = 0.D0
          QY(J)    = 0.D0
          QZ(J)    = 0.D0
	    do II=1,9
	      TOR(II)=0.d0
	    end do
          DO IS    = ISB,ISE                ! over atoms
            N      = N+1
	      if(LSHEJK)then
	        UX(N)=0.5d0*(VX(N)+HX(N))
	        UY(N)=0.5d0*(VY(N)+HY(N))
	        UZ(N)=0.5d0*(VZ(N)+HZ(N))
	      else
	        UX(N)=VX(N)
	        UY(N)=VY(N)
	        UZ(N)=VZ(N)
	      end if   
*
*	Molecular linear momenta:
*
            PX(J)  = PX(J)+ UX(N)
            PY(J)  = PY(J)+ UY(N)
            PZ(J)  = PZ(J)+ UZ(N)
	      TKET = TKET+0.5*(UX(N)**2+UY(N)**2+UZ(N)**2)*MASSDI(N)
*
*	Molecular angular momenta:
*
            QX(J)    = QX(J)+(UY(N)*WZ(N)-WY(N)*UZ(N))
            QY(J)    = QY(J)+(UZ(N)*WX(N)-WZ(N)*UX(N))
            QZ(J)    = QZ(J)+(UX(N)*WY(N)-WX(N)*UY(N))
	      TOR(1)=TOR(1)+MASSD(IS)*(WY(N)**2+WZ(N)**2)
	      TOR(5)=TOR(5)+MASSD(IS)*(WX(N)**2+WZ(N)**2)
	      TOR(9)=TOR(9)+MASSD(IS)*(WX(N)**2+WY(N)**2)  
	      TOR(2)=TOR(2)-MASSD(IS)*WX(N)*WY(N)
	      TOR(3)=TOR(3)-MASSD(IS)*WX(N)*WZ(N)
	      TOR(6)=TOR(6)-MASSD(IS)*WY(N)*WZ(N)
          END DO! OF IS 
	    TOR(4)=TOR(2)
	    TOR(7)=TOR(3)
	    TOR(8)=TOR(6) 
	    call INV3B3(TOR,TORI,DET,IER)
	    EKTR = EKTR+0.5*(PX(J)**2+PY(J)**2+PZ(J)**2)/SUMMD
          RXJ    = 0.D0
          RYJ    = 0.D0
          RZJ    = 0.D0 
	    N	= N-NSS
          DO IS    = ISB,ISE
            N      = N+1
            FAC=MASS(IS)/SUMMAS(ITYP)
* local atomic momenta
            VTX = UX(N)-PX(J)*FAC
            VTY = UY(N)-PY(J)*FAC
            VTZ = UZ(N)-PZ(J)*FAC
* local angular molecular momenta
            RXJ    = RXJ+(VTY*WZ(N)-WY(N)*VTZ)
            RYJ    = RYJ+(VTZ*WX(N)-WZ(N)*VTX)
            RZJ    = RZJ+(VTX*WY(N)-WX(N)*VTY)
          END DO! OF IS    
* rotational kinetic energy
	    if(IER.eq.0)then
	      EROTM=0.5*(TORI(1)*RXJ**2+(TORI(2)+TORI(4))*RXJ*RYJ+
     +                 TORI(5)*RYJ**2+(TORI(3)+TORI(7))*RXJ*RZJ+
     +                 TORI(9)*RZJ**2+(TORI(6)+TORI(8))*RYJ*RZJ)
	    else
	      ROTM=0.
            PP=RXJ**2+RYJ**2+RZJ**2
            if(PP.lt.1.d-50)go to 100
	      N	= N-NSS
            DO IS    = ISB,ISE
              N        = N+1
              SS=WX(N)**2+WY(N)**2+WZ(N)**2
              SP=WX(N)*RXJ+WY(N)*RYJ+WZ(N)*RZJ
              ROTM=ROTM+MASSD(IS)*(SS*PP-SP**2)
            end do
            if(ROTM.ge.1.d-50)then
	        EROTM=0.5d0*PP**2/ROTM
	      else
	        EROTM=0.d0
	      end if
	    end if
	    EKROT=EKROT+EROTM
          TROTM        = EROTM*ENERF*TEMPF
	  if(IPRINT.ge.10)
     +  write(*,'(2i5,e14.5,f10.3)')J,IER,EROTM,TROTM
 100	  continue
	end do ! of J
      end do ! of ITYP
      EKINT=TKET-EKTR-EKROT
      TTR         = EKTR*ENERF*TEMPF/(FNOP-1.d0)
      if(FNSTR.gt.0.01)then
          TROT        = 3.*EKROT*ENERF*TEMPF/FNSTR
	else
	  TROT=TEMP
	end if 
      if(FNSTI.gt.0.01)then
	  TINT	= 3.d0*EKINT*ENERF*TEMPF/FNSTI
	else
	  TINT=TEMP
	end if
*	if(ME.le.3)then
*	TTTT= 3.d0*TKET*ENERF*TEMPF/FNST
*	WW= 0.001*ENERF/FNOP 
*	write(*,'(a5,4f12.4)')'kien:',TKET*WW,EKTR*WW,EKROT*WW,EKINT*WW
*	write(*,'(I7,i4,a6,4f12.4)')MSTEP,ME,' temp:',TTTT,TTR,TROT,TINT 
*	end if
	RETURN
      END
*
*=============== GETMOI ========================================
*

      SUBROUTINE GETMOI(X,Y,Z,M,MOMENT,NSS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION M,MOMENT
      DIMENSION X(NSS),Y(NSS),Z(NSS),M(NSS)
      DIMENSION MOMENT(3,3),R(3)
*
      DO INERT           = 1,3
      DO IA              = 1,3
      TT                 = 0.D0
      MOMENT(INERT,IA)   = 0.D0
      IF(INERT.EQ.IA) TT = 1.D0
      DO I               = 1,NSS
      R(1)               =  X(I)
      R(2)               =  Y(I)
      R(3)               =  Z(I)
      SQS                =  X(I)**2+Y(I)**2+Z(I)**2
      MOMENT(INERT,IA)   = MOMENT(INERT,IA)+M(I)*(TT*SQS-R(INERT)*R(IA))
      END DO! OF I
      END DO! OF IA
      END DO! OF INERT
*
      RETURN
      END
*
*===================== INV3B3 ==========================================
*
      SUBROUTINE INV3B3(A,AI,DET,IER)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION A(9),AI(9)
*
      AI(1) = A(5)*A(9)-A(6)*A(8)
      AI(2) = A(3)*A(8)-A(2)*A(9)
      AI(3) = A(2)*A(6)-A(3)*A(5)
      AI(4) = A(6)*A(7)-A(4)*A(9)
      AI(5) = A(1)*A(9)-A(3)*A(7)
      AI(6) = A(3)*A(4)-A(1)*A(6)
      AI(7) = A(4)*A(8)-A(5)*A(7)
      AI(8) = A(2)*A(7)-A(1)*A(8)
      AI(9) = A(1)*A(5)-A(2)*A(4)
*
      DET   = A(1)*AI(1)+A(4)*AI(2)+A(7)*AI(3)
	SPUR=A(1)+A(5)+A(9)   
	IER=0
	if(dabs(DET).lt.1.d-6*dabs(SPUR**3))IER=1
      IF(DABS(DET).gt.0.D0) then
	  R=1.D0/DET
	else
	  R=0.d0
	  IER=2
	end if
*
      AI(1) = R*AI(1)
      AI(2) = R*AI(2)
      AI(3) = R*AI(3)
      AI(4) = R*AI(4)
      AI(5) = R*AI(5)
      AI(6) = R*AI(6)
      AI(7) = R*AI(7)
      AI(8) = R*AI(8)
      AI(9) = R*AI(9)
*
      RETURN
      END
*
*================= COPY =================================================
*
      SUBROUTINE COPY(A,B,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M),B(M)
      DO I   = 1,M
       A(I)  = B(I)
      END DO! OF I
      RETURN
      END
*
*========== GETMSS =======================================================
*
      DOUBLE PRECISION FUNCTION GETMSS(N)
      PARAMETER (N1=9,N2=16)
      CHARACTER*4 N
      CHARACTER*2 B
      CHARACTER*1 C
      DIMENSION B(N2),BM(N2)
      DIMENSION C(N1),CM(N1)
      DATA B /'LI','NA','MG','CL','CA','CR'
     X       ,'MN','FE','CO','NI','BR','AG','PC','XE','ME','LP'/
      DATA BM/6.939D0,22.9898D0,24.312D0,35.453D0,40.08D0,51.996D0
     X       ,54.938D0,55.847D0,58.9332D0,58.71D0,79.904D0,107.868D0
     X       ,0.D0,131.3D0,15.033D0,10.0D0/
      DATA C /'H','C','N','O','F','P','S','K','I'/
      DATA CM/1.008D0,12.011D0,14.007D0,15.999D0,18.998D0,30.974D0
     X       ,32.064D0,39.102D0,126.944D0/
*
      A       = 0.D0
      DO 100 I=1,N1
      IF(N(1:1).EQ.C(I)) A=CM(I)
  100 CONTINUE
*
      DO 200 I=1,N2
      IF(N(1:2).EQ.B(I)) A=BM(I)
  200 CONTINUE
*
      GETMSS=A
*
      RETURN
      END
*
*================= RDFINP =============================================
*
      SUBROUTINE RDFINP
      include"mdee.h" 
	character*128 str , TAKESTR
        integer IAUX(NS)
*     
        IO=5
        PRINT *
        PRINT *,'----------------------------------------'
        PRINT *,'*** PREPARATION FOR RDF CALCULATIONS ***'
        PRINT *,'----------------------------------------'
	NMGR(0)    = 'NONE    '
	do IS=1,NSITES
	  NGRI(IS)=0
	  do J=1,MXGRS
	    IGRI(IS,J)=0
	    MGRI(IS,J)=0
	  end do
	end do
      IF(ALLSTS) THEN
        N          = 0
        DO I       = 1,NSITES
          DO J       = I,NSITES
            N          = N+1
	      NGRI (I) = NGRI(I)+1
	      if(J.ne.I)NGRI(J)=NGRI(J)+1
	      if(NGRI(I).gt.MXGRS.or.NGRI(J).gt.MXGRS.or.N.gt.MAXGR)then
 	        write(*,*)' List of RDF limited to ',N,' pairs'
	        write(*,*)
     +' To calculate other RDF, specify list of RDF-s explicitly'           
	        go to 21
	      end if 
	      IGRI(I,NGRI(I))=J
	      MGRI(I,NGRI(I))=N
	      IGRI(J,NGRI(J))=I
	      MGRI(J,NGRI(J))=N
          END DO  
        END DO 
*
      ELSE
*     
        N          = 0    
	IE=-1 
        STR=TAKESTR(IO,IE)
        read(str,*)MAXRDF
        write(*,*)' Reading ', MAXRDF, ' pairs of RDFs' 
        DO K       = 1,MAXRDF
          N          = N+1  
	    IE=-1
	    STR=TAKESTR(IO,IE)
	    if(IE.lt.0)go to 19
	    if(str(1:1).eq.'&')then
            READ(str(2:24),*)IREP
	      do IR=1,IREP
	        STR=TAKESTR(IO,IE)
	        read(STR,*,err=99)I,J  
	        NGRI(I)=NGRI(I)+1 
	        if(J.ne.I)NGRI(J)=NGRI(J)+1
	        if(NGRI(I).gt.MXGRS)then
 	          write(*,*)' List of RDF exceeded for site ',I,' ',NM(I)
	          write(*,*)' increase MXGRS in prcm.h'
		  write(*,*)' RDF-s input interrupted'          
	          go to 19
	        end if 
	        if(NGRI(J).gt.MXGRS)then
 	          write(*,*)' List of RDF exceeded for site ',J,' ',NM(J)
	          write(*,*)' increase MXGRS in prcm.h'
		  write(*,*)' RDF-s input interrupted'          
	          go to 19
	        end if 
*  NGRI - number of RDF-pairs for this site
*  IGRI - 2-nd site for this RDF
*  MGRI - RDF number
	        IGRI(I,NGRI(I))=J
	        MGRI(I,NGRI(I))=N
	        IGRI(J,NGRI(J))=I
	        MGRI(J,NGRI(J))=N
	      end do
	    else 
	      IPER=1
	      read(str,*,err=99)I,J
	      NGRI(I)=NGRI(I)+1 
	      if(J.ne.I)NGRI(J)=NGRI(J)+1
	      if(NGRI(I).gt.MXGRS)then
 	        write(*,*)' List of RDF exceeded for site ',I,' ',NM(I)
	        write(*,*)' increase MXGRS in prcm.h'
		write(*,*)' RDF-s input interrupted'          
	        go to 19
	      end if 
	      if(NGRI(J).gt.MXGRS)then
 	        write(*,*)' List of RDF exceeded for site ',J,' ',NM(J)
	        write(*,*)' increase MXGRS in prcm.h'
	        write(*,*)' RDF-s input interrupted'          
	        go to 19
	      end if 
	      IGRI(I,NGRI(I))=J
	      MGRI(I,NGRI(I))=N
	      IGRI(J,NGRI(J))=I
	      MGRI(J,NGRI(J))=N
	    end if
        END DO
	end if      
 19	continue
*
  21  MAXRDF     = N
      if(TASKID.eq.MAST)
     +PRINT "('NR OF PAIRS IN RDF CALCULATIONS: ',I5)",N
      IF(MAXRDF.GT.MAXGR) STOP '!!! MAXRDF.GT.MAXGR IN RDFINP'
*
      if(TASKID.eq.MAST)
     +PRINT "('*** NR OF RADIAL DISTRIBUTION FUNCTIONS: ',I4)",N
*
      NRDF       = 0
      DO I       = 1,MAXGR
	  IIN(I)=0
	  IGR(I)=0
	  JGR(I)=0
        DO J       = 1,MAX
          IRDF(J,I)  = 0
        END DO
      END DO
*     
*  Check RDF list and calculate normalising factors
      DO ITYP    = 1,NTYPES
        NSPI       = NSPEC(ITYP)
        ISB        = ISADR(ITYP)+1
        ISE        = ISADR(ITYP +1)
        DO IS      = ISB,ISE 
	    NNR=NGRI(IS)
	    do K=1,NNR
	      JS = IGRI(IS,K)
	      if(IS.le.JS)then
		N    = MGRI(IS,K)
		JTYP = ITS(JS)
              NSPJ = NSPEC(JTYP)
              NMGR(N)    = NM(IS)//NM(JS)
              IGR (N)    = NSPI
              JGR (N)    = NSPJ
	        if(IS.eq.JS.and.NSPI.ne.1)then
	          IIN(N)=IIN(N)+(NSPI-1)*NSPI/2
	        else
	          IIN(N)=IIN(N)+NSPI*NSPJ
	        end if     
	      end if
	    end do
*
        END DO! OF IS
      END DO! OF ITYP
*
      DRDF     = RDFCUT/DFLOAT(MAX)
      DRDFI    =  1.D0/DRDF
      if(TASKID.eq.MAST)then
        PRINT "('*** ONE RDF SLICE IN ANGSTROM:  ',F8.5)",DRDF
        PRINT *
	end if 
*  Check RDF table
	if(IPRINT.ge.7)then   
	  write(*,*)' Calculate RDF:'  
	  do IS=1,NSITES
	    NNR=NGRI(IS) 
          do KS=1,NNR
	      JS=IGRI(IS,KS)            
	      K=MGRI(IS,KS)
	    write(*,'(a,2i5,3x,a4,a1,a4,i6,3i9)')
     +' sites ',IS,JS,NM(IS),'-',NM(JS),K,IGR(K),JGR(K),IIN(K)/2
	    end do                                     
	  end do
	end if
*
 22   RETURN      
 99   IE=99
	STR=TAKESTR(IO,IE)
      END                       
*
*================= RDFOUT =============================================
*
      SUBROUTINE RDFOUT(IOUT)
      include"mdee.h"
*
      CHARACTER*12 LABRDF
      DIMENSION GR(MAXS,MAXGR),GRO(MAXS,MAXGR),DF(MAXS)
      integer IRDF1(MAXS,MAXGR)
      DATA JINIT/0/, LABRDF /'rdf:        '/
      save GRO,FNRNR,JINIT,NRNR
*     
*      write(*,*)' RDFOUT:',NRDF,LGDMP,IHIST,NHIST
      IF(NRDF.EQ.0.and.(.not.LGDMP).and.(IHIST.lt.NHIST)) RETURN
CM
CM   collecting RDF
CM    
      if(IOUT.ne.0)then
      PRINT *
      PRINT *,'+++++++++++++++++++++++++++++++++++++++++++++++'
      PRINT *,'***      OUTPUT FROM RDF CALCULATIONS       ***'
      PRINT *,'+++++++++++++++++++++++++++++++++++++++++++++++'
      PRINT *
      PRINT "('*** RDF  COLLECTED  DURING:  ',I6,'   STEPS')",NRDF
      PRINT "('*** NR  OF  RDF TYPES  CALCULATED       ',I6)",MAXRDF
      PRINT "('*** CUT OFF RADIUS DIVIDED TO ',I6,' SLICES')",MAX
      PRINT *
*
      PRINT "(1X,130('='))"
      end if
      FNRDF        = DFLOAT(NRDF)
      if(JINIT.eq.0)then
      NRNR         = 0
      FNRNR        = 0.D0
      DO I         = 1,MAXRDF
      DO J         = 1,MAX
      GRO(J,I)      = 0.D0
      END DO! OF J
      END DO! OF I
      end if
*
      IF(LGRST.and.JINIT.eq.0) THEN
        OPEN (UNIT=20,file=frdf,status='old',FORM='UNFORMATTED',
     +        err=111)
        READ(20,err=111) NRNR
        DO I         = 1,MAXRDF
          READ(20,err=112,end=112) (GRO(J,I),J=1,MAX) 
        END DO
        PRINT "('*** RDF RESTART FOR ',I8,' OLD CONFIGURARTIONS')",NRNR
        FNRNR        = DFLOAT(NRNR)
        close (20)
        JINIT=1
      END IF 
	go to 113
 111  write(*,*)' RDF restart failed. File ',frdf
	write(*,*)' RDF are collected only from this run'
	go to 113
 112	write(*,*)' Invalid RDF restart file ',frdf
	frdf(1:1)='?'
	if(LGDMP)write(*,*)' RDF dump will be in file',frdf 
 113  continue
*
      FNRINV       =  1.D0/(FNRDF+FNRNR)
      if(NRDF.eq.0)then
        FNRDF=1.
      end if
      RDFINC       =  RDFCUT/DFLOAT(MAX)
*
      DO 200 I     = 1,MAXRDF
*
      FNI          = DFLOAT(IGR(I))
      FNJ          = DFLOAT(JGR(I))
      FIJ	   = dfloat(IIN(I)) 
*
      if(IOUT.ne.0)then
        write(*,*)
        write(*,'(a,a8,a,i8)')'#*** THIS PAIR: ',NMGR(I),
     +'   --->  TOTAL NR OF CONFIGURATIONS',NRNR+NRDF
        write(*,*)
        write(*,'(3X,5HR(A) ,5X,3HRDF,3X,7HINT RDF)')
*
      end if
      AVRNRI       =  FNRDF
      DNSNRI       = VOL/FNI
      DNSNRJ       = VOL/FNJ 
      RDFINT       = 0.D0 
	LABRDF(5:12)=NMGR(I)
*
      M            = 0
      DO 100 J     = 1,MAX
        M            = M+1    
        DIST         = RDFINC*(DFLOAT(J)-0.5)
        SHELLV       = 4.D0*PI*RDFINC*(DIST**2+RDFINC**2/12.D0)
        FIRDF        = DFLOAT(IRDF(J,I))
        RDFINT       = RDFINT+FIRDF/ (AVRNRI*FNI)
        RDFIJ        = VOL*FIRDF/(FIJ*AVRNRI*SHELLV)
*
        GR(J,I)      = (FNRDF*RDFIJ+FNRNR*GRO(J,I))*FNRINV
        RDFIJ        = GR(J,I)
	RDFINT	     = RDFINT+RDFIJ*SHELLV/DNSNRJ
        DF(M)        =  RDFIJ*DIST**2
        FINT         =  4.D0*PI*XINTGR(RDFINC,DF,M)
        FINTI        =  FINT/DNSNRI
        FINTJ        =  FINT/DNSNRJ
*
      IF(RDFINT.NE.0.D0.AND.DIST.LE.RDFCUT.and.IOUT.ne.0)PRINT 
     X"(4(1X,F8.4),10x,a12,3e12.5)",
     +DIST,RDFIJ,FINTI,FINTJ,LABRDF
  100 CONTINUE
  200 CONTINUE
*
      if(IOUT.ne.0)PRINT "(1X,130('='))"       
	NRDFT=NRDF+NRNR
      IF(LGDMP.and.(.not.lcheck).and.NRDFT.gt.0) THEN
      OPEN (UNIT   = 21,file=frdf,status='unknown',FORM='UNFORMATTED')
      WRITE(21) NRDFT
      DO I         = 1,MAXRDF
      WRITE(21) (GR(J,I),J=1,MAX)
      END DO! OF I
      PRINT "('*** RDF DUMP FOR ',I8,' CONFIGURARTIONS')",NRNR+NRDF
      close (21)
      END IF
*
      RETURN
      END
*
*=============== GETROT ==========================================
*
      SUBROUTINE GETROT
*
      include"mdee.h"
*
      DIMENSION XX(NS),XY(NS),XZ(NS),FMASS(NS)
      DIMENSION XS(NS),YS(NS),ZS(NS)
      DIMENSION FMOI(3,3)      
      DIMENSION  ROT(3,3),G(3),H(3)
*
      RF       = 1.
      N        = 0
      DO ITYP  = 1,NTYPES
      NSS      = NSITS (ITYP)
      NSP      = NSPEC (ITYP)
      ISB      = ISADR (ITYP)+1
      ISE      = ISADR (ITYP +1)
      JBE      = IADDR (ITYP)+1
      JEN      = IADDR (ITYP +1)
      DO J     = JBE,JEN
      I        = 0
      DO IS    = ISB,ISE
      I        = I+1
      N        = N+1
      XX   (I) = WX(N)
      XY   (I) = WY(N)
      XZ   (I) = WZ(N)
      FMASS(I) = MASSD(IS)
      END DO! OF IS
*
      CALL GETMOI(XX,XY,XZ,FMASS,FMOI,NSS)
*
      CALL HH3BY3(FMOI,G,H)
*
      DO K       = 1,3
      DO L       = 1,3
      RMX(J,K,L) = FMOI(K,L)
      END DO! OF K
      END DO! OF L
*
      END DO! OF J
*
      END DO! OF ITYP
*
      RETURN
      END
*
*========== HH3BY3 ============================================
*
      SUBROUTINE HH3BY3(A,D,E)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      data ICOUNT/0/
*
      DIMENSION A(3,3),D(3),E(3)
*
      H       = 0.
      SCALE   = ABS(A(3,1))+ABS(A(3,2))
      IF(SCALE.EQ.0) THEN
      E(3)    = A(3,2)
      ELSE
      A(3,1)  = A(3,1)/SCALE
      A(3,2)  = A(3,2)/SCALE
      H       = A(3,1)**2+A(3,2)**2+H
      F       = A(3,2)
      G       =-SIGN(SQRT(H),F)
      E(3)    = SCALE*G
      H       = H-F*G
      A(3,2)  =   F-G
*
      A(1,3)  = A(3,1)/H
      G       = A(1,1)*A(3,1)+A(2,1)*A(3,2)
      E(1)    = G/H
      F       = E(1)*A(3,1)
*
      A(2,3)  = A(3,2)/H
      G       = A(2,1)*A(3,1)+A(2,2)*A(3,2)
      E(2)    = G/H
      F       = F+E(2)*A(3,2)
*
      HH      = F/(H+H)
      F       = A(3,1)
      G       = E(1)-HH*F
      E(1)    = G
      A(1,1)  = A(1,1)-F*E(1)-G*A(3,1)
      F       = A(3,2)
      G       = E(2)-HH*F
      E(2)    = G
      A(2,1)  = A(2,1)-F*E(1)-G*A(3,1)
      A(2,2)  = A(2,2)-F*E(2)-G*A(3,2)
      END IF
      D(3)    = H
*
      SCALE   = 0.
      E(2)    = A(2,1)
      D(2)    = 0.
*
      E(1)    = 0.
      D(1)    = A(1,1)
      A(1,1)  = 1.
      IF(D(2).NE.0) THEN
      G       = A(2,1)*A(1,1)
      A(1,1)  = A(1,1)-A(1,2)*G
      END IF
      D(2)    = A(2,2)
      A(2,2)  = 1.
      A(2,1)  = 0.
      A(1,2)  = 0.
      IF(D(3).NE.0) THEN
      G       = A(3,1)*A(1,1)+A(3,2)*A(2,1)
      A(1,1)  = A(1,1)-A(1,3)*G
      A(2,1)  = A(2,1)-A(2,3)*G
      G       = A(3,1)*A(1,2)+A(3,2)*A(2,2)
      A(1,2)  = A(1,2)-A(1,3)*G
      A(2,2)  = A(2,2)-A(2,3)*G
      END IF
      D(3)    = A(3,3)
      A(3,3)  = 1.
      A(3,1)  = 0.
      A(1,3)  = 0.
      A(3,2)  = 0.
      A(2,3)  = 0.
*
*Diagonalize it:
*
      E(1)    = E(2)
      E(2)    = E(3)
      E(3)    = 0.
      DO 15 L = 1,3
      ITER    = 0
    1 CONTINUE
      DO 12 M = L,2
      DD      = ABS(D(M))+ABS(D(M+1))
      IF(ABS(E(M))+DD.EQ.DD) GO TO 2
   12 CONTINUE
      M       = 3
    2 CONTINUE
      IF(M.NE.L) THEN
      IF(ITER.EQ.30) then
        write(*,*) '!!! TOO MANY ITERATIONS!    in HH3BY3'
        ICOUNT=ICOUNT+1
        if(ICOUNT.gt.1000)stop
        return
      end if
      ITER    = ITER+1
      G       = (D(L+1)-D(L))/(2.*E(L))
      R       = SQRT(G**2+1.)
      G       = D(M)-D(L)+E(L)/(G+SIGN(R,G))
      S       = 1.
      C       = 1.
      P       = 0.
      DO 14 I = M-1,L,-1
      F       = S*E(I)
      B       = C*E(I)
      IF(ABS(F).GE.ABS(G)) THEN
      C       = G/F
      R       = SQRT(C**2+1.)
      E(I+1)  = F*R
      S       = 1./R
      C       = C*S
      ELSE
      S       = F/G
      R       = SQRT(S**2+1.)
      E(I+1)  = G*R
      C       = 1./R
      S       = S*C
      END IF
      G       = D(I+1)-P
      R       =(D(I)-G)*S+2.*C*B
      P       = S*R
      D(I+1)  = G+P
      G       = C*R-B
      F       = A(1,I+1)
      A(1,I+1)= S*A(1,I)+C*F
      A(1,I)  = C*A(1,I)-S*F
      F       = A(2,I+1)
      A(2,I+1)= S*A(2,I)+C*F
      A(2,I)  = C*A(2,I)-S*F
      F       = A(3,I+1)
      A(3,I+1)= S*A(3,I)+C*F
      A(3,I)  = C*A(3,I)-S*F
   14 CONTINUE
      D(L)    = D(L)-P
      E(L)    = G
      E(M)    = 0.
      GO TO 1
      END IF
   15 CONTINUE
      RETURN
      END
*
*========= ROTMAT =================================================
*
      SUBROUTINE ROTMAT(XX,YY,ZZ,IBG,N,NBTOP)
*
      include "mdee.h"
*
      DIMENSION XX(*),YY(*),ZZ(*)
      LOGICAL NBTOP
*
      IBEG  = IBG
      IEND  = IBG+N-1
*
      IF(NBTOP) THEN	! from box to principal:
      DO I  = IBEG,IEND
      AQ    = XX(I)
      BQ    = YY(I)
      CQ    = ZZ(I)
      XX(I) = RMX(I,1,1)*AQ+RMX(I,2,1)*BQ+RMX(I,3,1)*CQ
      YY(I) = RMX(I,1,2)*AQ+RMX(I,2,2)*BQ+RMX(I,3,2)*CQ
      ZZ(I) = RMX(I,1,3)*AQ+RMX(I,2,3)*BQ+RMX(I,3,3)*CQ
      END DO! OF I
      ELSE 		! from principal to box:
      DO I  = IBEG,IEND
      AQ    = XX(I)
      BQ    = YY(I)
      CQ    = ZZ(I)
      XX(I) = RMX(I,1,1)*AQ+RMX(I,1,2)*BQ+RMX(I,1,3)*CQ
      YY(I) = RMX(I,2,1)*AQ+RMX(I,2,2)*BQ+RMX(I,2,3)*CQ
      ZZ(I) = RMX(I,3,1)*AQ+RMX(I,3,2)*BQ+RMX(I,3,3)*CQ
      END DO! OF I
      END IF
*
      RETURN
      END
*
*======================= PBC ========================================
*
      subroutine PBC(XX,YY,ZZ)
      include "mdee.h"
      if(XX.gt. HBOXL)XX=XX-BOXL
      if(XX.lt.-HBOXL)XX=XX+BOXL
      if(YY.gt. HBOYL)YY=YY-BOYL
      if(YY.lt.-HBOYL)YY=YY+BOYL
      if(ZZ.gt. HBOZL)ZZ=ZZ-BOZL
      if(ZZ.lt.-HBOZL)ZZ=ZZ+BOZL
      return
      end 
*
**================ DISP ============================================
*
      SUBROUTINE DISP(FT,FM,FD,NI0,INN)  
      implicit real*8 (A-H,O-Z)
      DIMENSION FT(INN)                  
*
*  Dispersion calculation
*
      if(INN.le.NI0)then
        FM=FT(INN)
        FD=2.*FT(INN)
        return
      end if
      FM        = 0.0 
      DO INU    = NI0,INN
        FM        = FM+FT(INU) 
      END DO
      FM        = FM/(INN-NI0+1)
      FD        = 0.0
      DO INU    = NI0,INN
        FD        = FD+(FT(INU)-(FM))**2  
      END DO
      FD        = DSQRT(FD)/(INN-NI0)
      RETURN                     
      END              
*
*======================= TESTCONF ===================================
*   
	subroutine TESTCONF
*  Check configuration
	include "mdee.h"  
	logical OK/.true./
	do I=1,NSTOT-1
	  IMOL=NNUM(I)
	  ITYP=ITYPE(I)
	  ISB =NSITE(I) 
	  do J=I+1,NSTOT
	    JMOL=NNUM(J)
	    JTYP=ITYPE(J)
	    JSB=NSITE(J)
	    if(IMOL.ne.JMOL)then
	      DX=SX(I)-SX(J)
	      DY=SY(I)-SY(J)
	      DZ=SZ(I)-SZ(J)
	      call PBC (DX,DY,DZ)
	      RR=sqrt(DX**2+DY**2+DZ**2) 
	      if(RR.le.0.3*(SIGMA(ISB)+SIGMA(JSB)).or.RR.le.1.)then
		if(OK)then
		  OK=.false.
		  write(*,*)'!!! Bad configuration '
		  write(*,*)'The following atom pairs are too close:'
	 	  write(*,*)
     +' N1 site1  mol1  N2  site2  mol2      R      Eel      Elj'
		end if
  	        QIJ=ONE(ISB,JSB)
              A6=SIX(ISB,JSB)
	        B12=TWL(ISB,JSB)
		EES=0.001*ENERF*QIJ/RR
		ELJ=0.001*ENERF*(A6/RR**6+B12/RR**12)
	     	write(*,'(2(I5,1x,a4,1x,a6),f10.4,2e12.4)')
     +I,NM(ISB),NAME(ITYP),J,NM(JSB),NAME(JTYP),RR,EES,ELJ 
	      end if
	    end if
	  end do
	end do
	if(OK)write(*,*)' Configuration OK'
	return
	end 
*
*===================== XMOLDUMP ======================================
*
*   4. Dump configuration in "XMOL" format
*   --------------------------------------
	subroutine XMOLDUMP
	include "mdee.h"     
	open(unit=21,file=fbio,status='unknown') 
      NSDP=0
	do ITYP=1,NTYPES
	    NSDP=NSDP+NSPEC(ITYP)*NSITS(ITYP)
	end do
	write(21,*)NSDP 
	fultim=TIM+NSTEP*DT
	write(21,*)' after ',fultim*1.d12,' ps'
	do I=1,NSTOT 
        ITYP=ITYPE(I)
	  IS=NSITE(I)
	  write(21,'(a2,3(2x,f12.5))')NM(IS)(1:2),SX(I),SY(I),SZ(I)
	end do  
	close(21)
	return
	end
*       
*==================== RAND =============================================
*
*	subroutine RANDU(IX,IY,W)
	real*4 function RAND()
	implicit real*8 (A-H,O-Z)
	data ITWO/2/,M2/0/,IX/14371/,IA,IC,MIC/3*0/,S/0.d0/
	save IX,IY
	if(M2.ne.0)go to 20
	  M=1
 10	  M2=M
	  M=ITWO*M2
	  if(M.gt.M2)go to 10
	  HALFM=M2
	  IA0=HALFM*datan(1.d0)/8.d0
	  IC0=HALFM*(0.5d0-dsqrt(3.d0)/6.d0)
	  IA=8*IA0+5
	  IC=2*IC0+1
	  MIC=(M2-IC)+M2
	  S=sngl(0.5d0/HALFM)
 20	IY=IX*IA
	if(IY.gt.MIC)IY=(IY-M2)-M2
	IY=IY+IC
	if(IY/2.gt.M2)IY=(IY-M2)-M2
	if(IY.lt.0)IY=(IY+M2)+M2
	IX=IY
	RAND=IY*S
	return
	end 
* 









