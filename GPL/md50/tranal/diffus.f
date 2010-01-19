C    Tranal v.5.0
C    Computation of MSD and diffusion
C    
      program DIFFUS
      include "tranal.h"
*  NTCD - num of time stamps to compute MSD
      parameter (NTCD=1002)
      real*8 RCUM(NTCD),RCUMX(NTCD),RCUMY(NTCD),RCUMZ(NTCD),TTIM(NTCD),
     +TRCX(NTCD,NPART),TRCY(NTCD,NPART),TRCZ(NTCD,NPART)
      integer NCUM(NTCD),ISHX(NPART),ISHY(NPART),ISHZ(NPART)
      character*64 FILDIF
      logical LCOM
      data IREC/1/,IINIT/0/,JINIT/0/
      namelist /DIFF/IDF,NTT,DTT,IAT,LCOM,BREAKM,FILDIF,FBEG
*
* Parameters:
*    IDF - molecule time for which diffusion is calculated
*    DTT - time interval for MSD calculations in s
*    NTT - number of steps for calculation of MSD and diffusion
*    IAT - atom number in molecule motion of which is followed to compute MSD
*          (if 0, COM is followed)
*    LCOM - whether to correct for COM motion
*    BREAKM - time interval considered as "break" in trajectory
*    FILDIF - name of file for output
*    FBEG - start of linear fitting of MSD, as a fraction of the whole
*           interval  DTT*NTT
      MTR=NTCD-1
      ISTEP=1
      LCOM=.false.
      FBEG=0.2
      call SETTRAJ
      read(*,DIFF)
      if(NTT+1.gt.NTCD)stop ' increase NTCD'
      do IT=1,NTT
	RCUM(IT)=0.d0
	NCUM(IT)=0
      end do   
      do I=1,NPART
	ISHX(I)=0
	ISHY(I)=0
	ISHZ(I)=0
      end do
      IEND=0
      do while(IEND.eq.0)
        do III=1,ISTEP
          call READCONF(IEND)
          if(IEND.ne.0)go to 34
        end do
        ITYP		= IDF
        N           = ISADDR(ITYP)
        I           = IADDR(ITYP)
        IML=NSPEC(ITYP)
	XCM=0.
	YCM=0.
	ZCM=0.
        call GETCOM
        DO J        = 1,IML   ! sum over molecules
          I           = I+1
          if(IAT.gt.0)then
	    NSADR=ISADDR(ITYP)+NSITS(ITYP)*(J-1)+IAT
            XX       = SX(NSADR)
            YY       = SY(NSADR)
            ZZ       = SZ(NSADR)
          else
            XX=X(I)
            YY=Y(I)
            ZZ=Z(I)
          end if
          XX       = XX+ISHX(J)*BOXL
          YY       = YY+ISHY(J)*BOYL
          ZZ       = ZZ+ISHZ(J)*BOZL
*  PBC**(-1)
	  if(IINIT.ne.0.or.IREC.ne.1)then
	    IREC1=IREC-1
	    if(IREC1.le.0)IREC1=IREC1+MTR 
	    if(XX-TRCX(IREC1,J).gt. HBOXL)then
	      ISHX(J)=ISHX(J)-1
	      XX=XX-BOXL
	    end if
	    if(XX-TRCX(IREC1,J).lt.-HBOXL)then
	      ISHX(J)=ISHX(J)+1
	      XX=XX+BOXL
	    end if
	    if(YY-TRCY(IREC1,J).gt. HBOYL)then
	      ISHY(J)=ISHY(J)-1
	      YY=YY-BOYL
	    end if
	    if(YY-TRCY(IREC1,J).lt.-HBOYL)then
	      ISHY(J)=ISHY(J)+1
	      YY=YY+BOYL
	    end if
	    if(ZZ-TRCZ(IREC1,J).gt. HBOZL)then
	      ISHZ(J)=ISHZ(J)-1
	      ZZ=ZZ-BOZL
	    end if
	    if(ZZ-TRCZ(IREC1,J).lt.-HBOZL)then
	      ISHZ(J)=ISHZ(J)+1
	      ZZ=ZZ+BOZL
	    end if
	  end if	  
*  remember 
	  TRCX(IREC,J)=XX	  
	  TRCY(IREC,J)=YY	  
	  TRCZ(IREC,J)=ZZ  
	  JM=J
	  XCM=XCM+XX
	  YCM=YCM+YY
	  ZCM=ZCM+ZZ
        END DO! OF J
        XCM=XCM/IML
        YCM=YCM/IML
        ZCM=ZCM/IML
        if(JINIT.eq.0)then
	  XCM0=XCM
	  YCM0=YCM
	  ZCM0=ZCM
	  JINIT=1
        end if
        if(LCOM)then
          TRCX(IREC,IML+1)=XCM
          TRCY(IREC,IML+1)=YCM
          TRCZ(IREC,IML+1)=ZCM
        else
          TRCX(IREC,IML+1)=0.
          TRCY(IREC,IML+1)=0.
          TRCZ(IREC,IML+1)=0.
        end if
        TTIM(IREC)=FULTIM*1.d-12            ! in  s
*      if(IPRINT.ge.6)write(*,'(4f10.3,I7,2x,3I3)')
*     +FULTIM,TRCX(IREC,ICNT),TRCY(IREC,ICNT),TRCZ(IREC,ICNT),
*     +IREC,ISHX(ICNT),ISHY(ICNT),ISHZ(ICNT)
	if(IPRINT.ge.6)write(*,'(f12.1,3f12.6)')
     +    FULTIM,XCM-XCM0,YCM-YCM0,ZCM-ZCM0
*  calculate <R(t)-R(0)>
	IRB=IREC  
	IRB1=IRB+1
	if(IRB1.gt.MTR)IRB1=1
	IROUND=0
	do ITIM=1,NTT
	  TBACK=FULTIM*1.d-12-ITIM*DTT
*	  write(*,*)ITIM,FULTIM,TBACK*1.d12,IRB,TTIM(IRB)*1.d12,IRB1,
*     + TTIM(IRB1)*1.d12
 10	  if(TTIM(IRB).le.TBACK)then
	    if(TTIM(IRB1)-TTIM(IRB).lt.5.*DTT)then
	    AL=(TBACK-TTIM(IRB))/(TTIM(IRB1)-TTIM(IRB))      
	    if(AL.lt.0.d0.or.AL.gt.1.d0)write(*,'(a,f9.5,I7,3e14.5)')
     +  ' strange ALPHA: ',AL,IRB,TTIM(IRB),TBACK,TTIM(IRB1)  
          DO J        = 1,IML   ! sum over molecules
	      X0=(TRCX(IRB,J)-TRCX(IRB,IML+1))*(1.-AL)+
     +           (TRCX(IRB1,J)-TRCX(IRB1,IML+1))*AL
	      Y0=(TRCY(IRB,J)-TRCY(IRB,IML+1))*(1.-AL)+
     +           (TRCY(IRB1,J)-TRCY(IRB1,IML+1))*AL
	      Z0=(TRCZ(IRB,J)-TRCZ(IRB,IML+1))*(1.-AL)+
     +           (TRCZ(IRB1,J)-TRCZ(IRB1,IML+1))*AL
	      RAX=(X0-TRCX(IREC,J)+TRCX(IREC,IML+1))**2
	      RAY=(Y0-TRCY(IREC,J)+TRCY(IREC,IML+1))**2
	      RAZ=(Z0-TRCZ(IREC,J)+TRCZ(IREC,IML+1))**2
	      RR=RAX+RAY+RAZ
	      RCUMX(ITIM)=RCUMX(ITIM)+RAX
	      RCUMY(ITIM)=RCUMY(ITIM)+RAY
	      RCUMZ(ITIM)=RCUMZ(ITIM)+RAZ
	      RCUM(ITIM)=RCUM(ITIM)+RR
	      NCUM(ITIM)=NCUM(ITIM)+1  
	      if(IPRINT.ge.8)write(*,'(2I6,2f12.6)')
     +      ITIM,J,ITIM*DTT*1.d12,sqrt(RR)
	    end do   
	    end if
	  else
	    IRB1=IRB
	    IRB=IRB-1
	    if(IRB.le.0)then
	      if(IINIT.eq.0)go to 30
	      IRB=MTR
	      IROUND=1
	    end if
	    if(IRB.le.IREC.and.IROUND.eq.1)go to 30
	    go to 10
	  end if
	end do ! of ITIM 
 30	IREC=IREC+1
	if(IREC.gt.MTR)then
	  IREC=1
	  IINIT=1
	end if
      end do
 34   write(*,*)IAN,' configurations analysed'
	open(unit=9,file=fildif,status='unknown')
	write(*,'(a,i5)')
     +'# translational time -correlation function of type',IDF
	write(9,'(a)')'#  D units is 10**-9 m**2/s '
	write(9,'(a,I3)')
     +'#  translational time -correlation function of type',IDF
	write(9,'(a)')
     +'#  t(ps)  Dav    Ddif    sqrt(MSD)(A)     Dx     Dy     Dz'
	write(*,*)NTT,DTT 
	TSHO=0.
	TSHO=0.
        XYS=0.
        YS=0.
        XS=0.
        X2S=0.
        S=0
        TLAST=0.
	do I=1,NTT-1
	  TIM=I*DTT*1.d12
	  if(NCUM(I).gt.0)then
	    TCOR=10.*RCUM(I)/(6.d0*TIM*NCUM(I))
	    TCX=10.*RCUMX(I)/(2.d0*TIM*NCUM(I))
	    TCY=10.*RCUMY(I)/(2.d0*TIM*NCUM(I))
	    TCZ=10.*RCUMZ(I)/(2.d0*TIM*NCUM(I))
	    TSH2=RCUM(I)/NCUM(I)
	    TSH=sqrt(TSH2)
	    TSD=1.d-11*TSH*(TSH-TSHO)/(3.*DTT) 
	    TSHO=TSH
            write(9,'(f10.3,2f10.5,f10.3,3f10.5,i9)')
     +TIM,TCOR,TSD,TSH,TCX,TCY,TCZ,NCUM(I)
      write(*,'(f10.3,2f10.5,f10.3,3f10.5,i9)')
     +TIM,TCOR,TSD,TSH,TCX,TCY,TCZ,NCUM(I)
            if(I*1./NTT.gt.FBEG)then
              S=S+1.
              XS=XS+TIM
              X2S=X2S+TIM**2
              YS=YS+TSH2
              XYS=XYS+TIM*TSH2
            end if
            TLAST=I*DTT*1.d12
	  end if
	end do
        DDD=10.*(S*XYS-YS*XS)/(6.*(S*X2S-XS**2))
      write(*,'(a,f12.6,a)')
     +' Diffusion coefficient ',DDD,' 10**(-5) cm/S**2'
      write(*,'(a,f10.3,a,f10.3,a)')'#  for MSD between ',
     +FBEG*NTT*DTT*1.d12,' and ',TLAST,' ps'
      write(9,'(a,f12.6,a)')
     +'#  Diffusion coefficient ',DDD,' 10**(-5) cm/S**2'
      write(9,'(a,f10.3,a,f10.3,a)')'#  for MSD between ',
     +FBEG*NTT*DTT*1.d12,' and ',TLAST,' ps'
      stop
      end
