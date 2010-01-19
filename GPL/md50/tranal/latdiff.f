C    Tranal v.5.0
C    Computation of lateral MSD and diffusion in lipid bilayers
C    
      program LATDIFF
      include "tranal.h"
*  NTCD - num of time stamps to compute MSD
      parameter (NTCD=2002)
      real*8 RCUM(NTCD),RCUMX(NTCD),RCUMY(NTCD),RCUMZ(NTCD),TTIM(NTCD),
     +TRCX(NTCD,NPART),TRCY(NTCD,NPART),DRZ(NPART)
      real*8 XT(NPART),YT(NPART),ZT(NPART)
      integer NCUM(NTCD),ISHX(NPART),ISHY(NPART)
      character*64 FILDIF
      logical LCOM
      data IREC/1/,IINIT/0/,JINIT/0/
      namelist /DIFF/IDF,NTT,DTT,IAT,ISTEP,LCOM,BREAKM,FILDIF,IRT,FBEG
*
* Parameters:
*    IDF - molecule type for which diffusion is calculated
*    DTT - time interval for MSD calculations (in s)
*    NTT - number of steps for calculation of MSD 
*    IAT - atom number in molecule motion of which is followed to compute MSD
*          (if 0, molecule COM is used)
*    ISTEP - how often take configurations for analysis
*    IRT  - type number of lipid molecules 
*    LCOM - whether to correct for monolayer COM motion 
*    BREAKM - time interval considered as "break" in trajectory
*    FILDIF - name of file for output
*    FBEG - start of linear fitting of MSD, as a fraction of the whole
*           interval  DTT*NTT
      MTR=NTCD-2
      LCOM=.true.
      FBEG=0.2
      call SETTRAJ
      read(*,DIFF)
      if(NTT+2.gt.NTCD)stop ' increase NTCD'
      if(IRT.le.0.or.IRT.gt.NTYPES)then
        write(*,*)' wrong lipid type number (IRT):',IRT
        stop
      end if
      if(IDF.le.0.or.IDF.gt.NTYPES)then
        write(*,*)' wrong type number (IDF):',IDF
        stop
      end if
      do IT=1,NTT
	RCUM(IT)=0.d0
	NCUM(IT)=0
      end do   
      do I=1,NPART
	ISHX(I)=0
	ISHY(I)=0
      end do
      IEND=0
      do while(IEND.eq.0)
        do JJ=1,ISTEP
          call READCONF(IEND)
          if(IEND.ne.0)go to 34
        end do
        ITYP		= IDF
        N           = ISADDR(IDF)
        I           = IADDR(IDF)
        IML=NSPEC(IDF)
        call GETCOM
*  Computation of bilayer mid-plane
        NRFM=0
	ZRS=0.
	do I=1,NSPEC(IRT)
	  IMOL=I+IADDR(IRT)
	  if(NRFM.eq.0)then
	    ZRS=Z(IMOL)
	    ZRF=ZRS
	    NRFM=1
	  else
	    ZZ=Z(IMOL)
	    DZ=Z(IMOL)-ZRF
	    if(DZ.gt.HBOZL)ZZ=ZZ-BOZL
	    if(DZ.lt.-HBOZL)ZZ=ZZ+BOZL
	    ZRS=ZRS+ZZ
	    NRFM=NRFM+1
	    ZRF=ZRS/NRFM
	  end if
	  if(IPRINT.ge.8)write(*,'(i5,3f12.4)')I,Z(IMOL),ZZ,ZRF
	end do
	if(ZRF.gt.HBOZL)ZRF=ZRF-BOZL
	if(ZRF.lt.-HBOZL)ZRF=ZRF+BOZL
	if(IPRINT.ge.6)write(*,*)'Time=',FULTIM,' Zref = ',ZRF
        NN1=0
        NN2=0
        DO J        = 1,IML   ! sum over molecules of type IDF
          I           = J+IADDR(IDF)   ! total mol number
          if(IAT.gt.0)then
	    NSADR=ISADDR(IDF)+NSITS(IDF)*(J-1)+IAT
            XX       = SX(NSADR)
            YY       = SY(NSADR)
            ZZ       = SZ(NSADR)
          else
            XX=X(I)
            YY=Y(I)
            ZZ=Z(I)
          end if
          XX       = XX+ISHX(I)*BOXL
          YY       = YY+ISHY(I)*BOYL
          call PBC(XX,YY,ZZ)
*  PBC**(-1)
	  if(IINIT.ne.0.or.IREC.ne.1)then
	    IREC1=IREC-1
	    if(IREC1.le.0)IREC1=IREC1+MTR 
	    if(XX-TRCX(IREC1,J).gt. HBOXL)then
	      ISHX(I)=ISHX(I)-1
	      XX=XX-BOXL
	    end if
	    if(XX-TRCX(IREC1,J).lt.-HBOXL)then
	      ISHX(I)=ISHX(I)+1
	      XX=XX+BOXL
	    end if
	    if(YY-TRCY(IREC1,J).gt. HBOYL)then
	      ISHY(I)=ISHY(I)-1
	      YY=YY-BOYL
	    end if
	    if(YY-TRCY(IREC1,J).lt.-HBOYL)then
	      ISHY(I)=ISHY(I)+1
	      YY=YY+BOYL
	    end if
	  end if	  
*  remember 
	  TRCX(IREC,J)=XX	  
	  TRCY(IREC,J)=YY	  
	  DZ=ZZ-ZRF
	  if(DZ.gt.HBOZL)DZ=DZ-BOZL
	  if(DZ.lt.-HBOZL)DZ=DZ+BOZL
          DRZ(J)=DZ
        END DO! OF J
*  Correction for COM motion of lipids
        if(LCOM)then
	  XCM1=0.
	  YCM1=0.
	  ZCM1=0.
	  XCM2=0.
	  YCM2=0.
	  ZCM2=0.
          DO J        = 1,NSPEC(IRT)   ! sum over molecules of type IRT
            I           = J+IADDR(IRT)
            if(IRT.eq.IDF)then
              XX=TRCX(IREC,J)
              YY=TRCY(IREC,J)
              ZZ=Z(I)
            else
              XX=X(I)
              YY=Y(I)
              ZZ=Z(I)
              call PBC(XX,YY,ZZ)
              XX       = XX+ISHX(I)*BOXL
              YY       = YY+ISHY(I)*BOYL
*  PBC**(-1)
	      if(IINIT.ne.0.or.IREC.ne.1)then
	        if(XX-XT(I).gt. HBOXL)then
	          ISHX(I)=ISHX(I)-1
	          XX=XX-BOXL
	        end if
	        if(XX-XT(I).lt.-HBOXL)then
	          ISHX(I)=ISHX(I)+1
	          XX=XX+BOXL
	        end if
	        if(YY-YT(I).gt. HBOYL)then
	          ISHY(I)=ISHY(I)-1
	          YY=YY-BOYL
	        end if
	        if(YY-YT(I).lt.-HBOYL)then
	          ISHY(I)=ISHY(I)+1
	          YY=YY+BOYL
	        end if
	      end if	
              XT(I)=XX
              YT(I)=YY
            end if
	    DZ=ZZ-ZRF
	    if(DZ.gt.HBOZL)DZ=DZ-BOZL
	    if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	    if(DZ.gt.0)then
	      XCM1=XCM1+XX
	      YCM1=YCM1+YY
	      ZCM1=ZCM1+ZZ
	      NN1=NN1+1
	    else
	      XCM2=XCM2+XX
	      YCM2=YCM2+YY
	      ZCM2=ZCM2+ZZ
	      NN2=NN2+1
	    end if
          end do
          if(NN1.ne.NN2)write(*,*)' different layers: ',NN1,NN2
          XCM1=XCM1/NN1
          YCM1=YCM1/NN1
          ZCM1=ZCM1/NN1
          XCM2=XCM2/NN2
          YCM2=YCM2/NN2
          ZCM2=ZCM2/NN2
	  TRCX(IREC,IML+1)=XCM1
	  TRCY(IREC,IML+1)=YCM1
	  TRCX(IREC,IML+2)=XCM2
	  TRCY(IREC,IML+2)=YCM2
        else 
	  TRCX(IREC,IML+1)=0.
	  TRCY(IREC,IML+1)=0.
	  TRCX(IREC,IML+2)=0.
	  TRCY(IREC,IML+2)=0.
        end if
        TTIM(IREC)=FULTIM*1.d-12            ! in  s
*      if(IPRINT.ge.6)write(*,'(5f10.3,I7,2x,3I3)')
*     +FULTIM,X(11),Y(11),XT(11),YT(11),
*     +IREC,ISHX(11),ISHY(11)
	if(IPRINT.ge.6.and.LCOM)write(*,'(a,4f12.5)')
     +    '  COM shift: ',XCM1,YCM1,XCM2,YCM2
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
	    if(DRZ(J).gt.0)then
	      X0=(TRCX(IRB,J)-TRCX(IRB,IML+1))*(1.-AL)+
     +            (TRCX(IRB1,J)-TRCX(IRB1,IML+1))*AL
	      Y0=(TRCY(IRB,J)-TRCY(IRB,IML+1))*(1.-AL)+
     +            (TRCY(IRB1,J)-TRCY(IRB1,IML+1))*AL
	      RAX=(X0-TRCX(IREC,J)+TRCX(IREC,IML+1))**2
	      RAY=(Y0-TRCY(IREC,J)+TRCY(IREC,IML+1))**2
	    else
	      X0=(TRCX(IRB,J)-TRCX(IRB,IML+2))*(1.-AL)+
     +            (TRCX(IRB1,J)-TRCX(IRB1,IML+2))*AL
	      Y0=(TRCY(IRB,J)-TRCY(IRB,IML+2))*(1.-AL)+
     +            (TRCY(IRB1,J)-TRCY(IRB1,IML+2))*AL
	      RAX=(X0-TRCX(IREC,J)+TRCX(IREC,IML+2))**2
	      RAY=(Y0-TRCY(IREC,J)+TRCY(IREC,IML+2))**2
	    end if
	      RR=RAX+RAY
	      RCUMX(ITIM)=RCUMX(ITIM)+RAX
	      RCUMY(ITIM)=RCUMY(ITIM)+RAY
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
     +'# lateral MSD of type',IDF
	write(9,'(a)')'#  D units is 10**-9 m**2/s '
	write(9,'(a,I3)')
     +'#  lateral MSD of type',IDF
	if(LCOM)write(9,'(a)')
     +'#  correction for COM motion of each layear is made'
	write(9,'(a)')
     +'#  t(ps)  Dav    Ddif    sqrt(MSD)(A)     Dx     Dy      Nav'
	TSHO=0.
        XYS=0.
        YS=0.
        XS=0.
        X2S=0.
        S=0
	write(*,*)NTT,DTT 
        TLAST=0.
	do I=1,NTT-1
	  TIM=I*DTT*1.d12
	  if(NCUM(I).gt.0)then
	    TCOR=10.*RCUM(I)/(4.d0*TIM*NCUM(I))
	    TCX=10.*RCUMX(I)/(2.d0*TIM*NCUM(I))
	    TCY=10.*RCUMY(I)/(2.d0*TIM*NCUM(I))
	    TSH2=RCUM(I)/NCUM(I)
	    TSH=sqrt(TSH2)
	    TSD=1.d-11*TSH*(TSH-TSHO)/(2.*DTT) 
	    TSHO=TSH
            write(9,'(f10.3,2f10.5,f10.3,2f10.5,i9)')
     +TIM,TCOR,TSD,TSH,TCX,TCY,NCUM(I)
      write(*,'(f10.3,2f10.5,f10.3,2f10.5,i9)')
     +TIM,TCOR,TSD,TSH,TCX,TCY,NCUM(I)
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
      DDD=2.5*(S*XYS-YS*XS)/(S*X2S-XS**2)
      write(*,'(a,f12.6,a)')
     +' Lateral diffusion ',DDD,' 10**(-5) cm/S**2'
      write(*,'(a,f10.3,a,f10.3,a)')'#  for MSD between ',
     +FBEG*NTT*DTT*1.d12,' and ',TLAST,' ps'
      write(9,'(a,f12.6,a)')
     +'#  Lateral diffusion ',DDD,' 10**(-5) cm/S**2'
      write(9,'(a,f10.3,a,f10.3,a)')'#  for MSD between ',
     +FBEG*NTT*DTT*1.d12,' and ',TLAST,' ps'
      stop
      end
