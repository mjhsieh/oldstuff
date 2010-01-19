C========================================
C
C     MDynaMix v.4.3
C
*     PART 11
*
*     File util.f
*     -------------
*
*============ GAUSS ===============================================
*
      real*8 FUNCTION GAUSS(DUMMY)
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (A1=3.949846138,A3=0.252408784)
      PARAMETER (A5=0.076542912,A7=0.008355968)
      PARAMETER (A9=0.029899776)
      SUM   = 0.0
      DO I  = 1,12
      SUM   = SUM + RANF()
      END DO                   
      R     = (SUM-6.0)/4.0
      R2    = R*R
      GAUSS = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R
      RETURN
      END
*
*============================== RANF =========================
*
      real*8 FUNCTION RANF()
      IMPLICIT real*8 (A-H,O-Z)
      INTEGER L,C,M,SEED
      PARAMETER (L=1029,C=221591,M=1048576)
      SAVE SEED
      DATA SEED /0/
      SEED    = MOD(SEED*L+C,M)
      RANF    = DFLOAT(SEED)/M
      RETURN	
      END                                            
*
*=============== RANQ =========================================
*
      SUBROUTINE RANQ(QP)
      IMPLICIT real*8 (A-H,O-Z)
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
*=============== ROTATE ==============================================
*
      SUBROUTINE ROTATE(Q,X,Y,Z)
      IMPLICIT real*8 (A-H,O-Z)
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
*================== JACOBI =============================================
*
*   rotation matrix
*
      SUBROUTINE JACOBI(A,V)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION A(3,3),V(3,3)
      N=3
      RHO=1.0D-8
      TES=0.D0
      SCL=0.D0
      DO 10 I=1,N
      SCL=SCL+A(I,I)*A(I,I)
   10 CONTINUE
      SCL=DSQRT(DABS(SCL))/FLOAT(N)
      DO 20 I=1,N
      DO 20 J=1,N
   20 A(I,J)=A(I,J)/SCL
      DO 30 I=1,N
      DO 30 J=1,N
      V(I,J)=0.D0
      IF(I.EQ.J)V(I,J)=1.D0
   30 CONTINUE
      DO 100 I=2,N
      IU1=I-1
      DO 100 J=1,IU1
  100 TES=TES+2.D0*A(I,J)*A(I,J)
      TES=DSQRT(DABS(TES))
      M=0
  105 TES=TES/FLOAT(N)
      IF(TES.LT.RHO)TES=RHO
  110 DO 165 I=2,N
      IUU1=I-1
      DO 165 J=1,IUU1
      IF(DABS(A(I,J))-TES)165,115,115
  115 M=1
      V1=A(J,J)
      V2=A(I,J)
      V3=A(I,I)
      U=0.5D0*(V1-V3)
      IF(DABS(U)-RHO)120,125,125
  120 OMG=-1.D0
      GO TO 130
  125 SQVVUU=DSQRT(DABS(V2*V2+U*U))
      OMG=-V2/SQVVUU
      IF(U.LT.0.D0)OMG=-OMG
  130 SQOOO=DSQRT(DABS(2.0*(1.00+DSQRT(DABS(1.00-OMG*OMG)))))
      S=OMG/SQOOO
      C=DSQRT(DABS(1.00-S*S))
      DO 160 K=1,N
      IF(K-I)140,135,135
  135 TEM=A(K,J)*C-A(K,I)*S
      A(K,I)=A(K,J)*S+A(K,I)*C
      A(K,J)=TEM
      GO TO 155
  140 IF(K-J)145,150,150
  145 TEM=A(J,K)*C-A(I,K)*S
      A(I,K)=A(J,K)*S+A(I,K)*C
      A(J,K)=TEM
      GO TO 155
  150 TEM=A(K,J)*C-A(I,K)*S
      A(I,K)=A(K,J)*S+A(I,K)*C
      A(K,J)=TEM
  155 TEM=V(K,J)*C-V(K,I)*S
      V(K,I)=V(K,J)*S+V(K,I)*C
      V(K,J)=TEM
  160 CONTINUE
      A(J,J)=V1*C*C+V3*S*S-2.D0*V2*S*C
      A(I,I)=V1*S*S+V3*C*C+2.D0*V2*S*C
      A(I,J)=(V1-V3)*S*C+V2*(C*C-S*S)
  165 CONTINUE
      IF(M-1)175,170,170
  170 M=0
      GO TO 110
  175 IF(TES-RHO)180,180,105
  180 DO 190 I=1,N
      DO 190 J=1,N
  190 A(I,J)=SCL*A(I,J)
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
*========== ORDER3x3=================================================
*
C   Set values of vector G in decreasing order and matrix A correspondingly
C
      subroutine ORDER3X3(A,G)
      real*8 A(3,3),G(3),B(3,3),H
      integer IO(3)
      IO(1)=1
      IO(2)=2
      IO(3)=3
      if(G(3).gt.G(2))then
         IO(2)=3                  !   G
         IO(3)=2
         H=G(3)
         G(3)=G(2)
         G(2)=H
         if(G(2).gt.G(1))then
           G(2)=G(1)
           G(1)=H
           IO(1)=3
           IO(2)=1
           if(G(3).gt.G(2))then
             H=G(3)
             G(3)=G(2)
             G(2)=H
             IO(2)=2
             IO(3)=1
           end if
         end if
       else
         if(G(2).gt.G(1))then
           H=G(1)
           G(1)=G(2)
           G(2)=H
           IO(1)=2
           IO(2)=1
           if(G(3).gt.G(2))then
             G(2)=G(3)
             G(3)=H
             IO(2)=3
             IO(3)=1
           end if
         end if
       end if
       if(G(1).lt.G(2).or.G(2).lt.G(3))write(*,*)' wrong G order:',G
       do I=1,3
         do J=1,3
           JJ=IO(J)
           B(I,J)=A(I,JJ)
         end do
       end do
       do I=1,3
          do J=1,3
             A(I,J)=B(I,J)
           end do
       end do
       return
       end
*....
*======== XINTGR =======================================================
*....
      FUNCTION XINTGR(H,Y,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(*)
      HT=H/3.0 
	XINTGR=0.d0
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
*================ LENG ===============================================
*
      function LENG(STR,LS)
      character*1 STR(LS)
      IL=0 
      do I=1,LS
	ISTR=ichar(STR(I))
	if(ISTR.le.15)then
	  LENG=I-1
	  return
	end if
        if(STR(I).ne.' '.and.ISTR.gt.16)IL=I
      end do 
      LENG=IL
      return
      end    

