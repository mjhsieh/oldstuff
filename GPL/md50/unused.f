*
*=============== MPOLES =============================================
*
      SUBROUTINE MPOLES(SITE,Q,DM,QM,OM,HM,NSB,NSE,NS)
      IMPLICIT real*8 (A-H,O-Z)
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

*
*=============== GETMOI ========================================
*

      SUBROUTINE GETMOI(X,Y,Z,M,MOMENT,NSS)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 M,MOMENT
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
*=============== GETROT ==========================================
*
      SUBROUTINE GETROT
*
      include"prcm.h"
*
      DIMENSION RX(NS),RY(NS),RZ(NS),FMASS(NS)
      DIMENSION FMOI(3,3)      
      DIMENSION G(3),H(3)
*
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
      RX   (I) = WX(N)
      RY   (I) = WY(N)
      RZ   (I) = WZ(N)
      FMASS(I) = MASSD(IS)
      END DO! OF IS
*
      CALL GETMOI(RX,RY,RZ,FMASS,FMOI,NSS)
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
      IMPLICIT real*8 (A-H,O-Z)
*
      DIMENSION A(3,3),D(3),E(3)
      data ICOUNT/0/
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
      include "prcm.h"
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
*============ FMELT ================================================
*
      SUBROUTINE FMELT(X0,Y0,Z0,NSP)
      IMPLICIT real*8  (A-H,O-Z)
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
*   Do read(STR,*)(NAME(I),I=1,N) where STR and NAME are characters
*   -------------------------------------------------------------------
C   (only because some stupid compilers do not understand this)
C========= PARSES ============================
C
	subroutine PARSES(STR,CH,N,L,ie)
	character STR(L),CHR
	character*32 CH(N)
	logical LI
C  LCH=32 is length of NAMOL in prcm.h
	LCH=32
	do I=1,N
	   do J=1,LCH
	      CH(I)(J:J)=' '
	   end do
	end do
	LI=.false.
	M=0
	ie=0
	do I=1,L
	  CHR=STR(I)
	  ICHR=ichar(CHR)
	  if(LI)then
C  Delimiters - space, null and tab
	    if(CHR.eq.' '.or.ICHR.eq.0.or.ICHR.eq.9)then		
	      LI=.false.
	    else
	      IC=IC+1
	      if(IC.le.LCH)CH(M)(IC:IC)=CHR
	    end if
	  else
	     if(CHR.ne.' '.and.ICHR.ne.0.and.ICHR.ne.9)then
		M=M+1
		if(M.gt.N)return
		CH(M)(1:1)=CHR
		IC=1
		LI=.true.
	     end if
	  end if
	end do
	if(M.lt.N)ie=1
	return
	end
*
*=============== NEXTITEM =================================
*
	function NEXTITEM(STR)
	character*128 STR
	do I=1,128
          if(STR(I:I).ne.' ')go to 10
	end do
	NEXTITEM=0
	return
 10	do J=I,128
          if(STR(J:J).eq.' ')go to 20
	end do
	NEXTITEM=128
	return
 20	do I=J,128
          if(STR(I:I).ne.' ')go to 30
	end do   
 30	NEXTITEM=I
	return
	end

