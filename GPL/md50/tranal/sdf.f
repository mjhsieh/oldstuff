C    Tranal v.5.0
C    Computation of Spatial Distribution Functions
C    Output file is suitable for visualization using "GOpenMol" program
C    
 	program SDF
	include "tranal.h"
C  max size of SDF grid
	parameter (NOMXM=200,NOMYM=200,NOMZM=200)
	parameter (NOIM=100,MSOR=100)
	dimension E0(3),D1(3),D2(3),EX(3),EY(3),EZ(3),RV(3)
	dimension XAV(NS),YAV(NS),ZAV(NS),XAV2(NS),YAV2(NS),ZAV2(NS)
	dimension IAA(NS),ISM(NS,NOIM),ISREP(NS)
	dimension ISOR(MSOR),IO1(NOIM),IO2(NOIM),IO3(NOIM)
	dimension ID3(0:NOMXM,0:NOMYM,0:NOMZM)
	dimension TMP(0:NOMXM)
	integer*8 NCOR
	character*2 NMM(NS)
	character*64 FILORI,FILCRD,AUX
	namelist /SDFIN/ NSOR,ISOR,NOI,IO1,IO2,IO3,FILORI,FILCRD,
     + NOMX,NOMY,NOMZ,RXMAX,RYMAX,RZMAX 
*    Recuired parameters:
*    NOI - number of equivalent local coord. systems 
*    IO1, IO2, IO3 (NOI) - define local coordinate system
*    FSET - file with atom numbers for calculation of average structure 
*
	NOMX=50
	NOMY=50
	NOMZ=50
	RXMAX=5.
	RYMAX=5.
	RZMAX=5.
	NOI=1
	NSOR=1
        call SETTRAJ
	read(*,SDFIN)
	if(NOI.gt.NOIM)stop 'increase NOIM !!!'
	if(NSOR.gt.MSOR)stop 'increase MSOR !!!'
	if(NOMX.gt.NOMXM)stop 'increase NOMXM !!!'
	if(NOMY.gt.NOMYM)stop 'increase NOMYM !!!'
	if(NOMZ.gt.NOMZM)stop 'increase NOMZM !!!'
	do IGR=1,NSOR   ! equivalent atoms "to"
	  if(ISOR(IGR).le.0.or.ISOR(IGR).gt.NSITES)then
	     write(*,*)' Wrong value ',IGR,' of array ISOR ',ISOR(IGR)
	     stop
	  end if
	end do
	do II=1,NOI   ! equivalent atom groups "from"
	  if(IO1(II).le.0.or.IO1(II).gt.NSITES)then
	     write(*,*)' Wrong value ',II,' of array IO1 ',IO1(II)
	     stop
	  end if
	  if(IO2(II).le.0.or.IO2(II).gt.NSITES)then
	     write(*,*)' Wrong value ',II,' of array IO2 ',IO2(II)
	     stop
	  end if
	  if(IO2(II).le.0.or.IO3(II).gt.NSITES)then
	     write(*,*)' Wrong value ',II,' of array IO3 ',IO3(II)
	     stop
	  end if
	end do
*  reading parameters for drawing the central molecule
	read(*,'(a)')AUX
	do I=1,61
	  if(AUX(I:I+2).eq.'ALL')go to 7
	end do
	read(AUX,*)NMSTR
	do I=1,NMSTR
	  read(*,*,err=902,end=902)ISREP(I),(ISM(I,II),II=1,NOI)
	end do
	write(*,*)' read ',NMSTR,' reference atom coordinates'
	go to 8
 902	  write(*,*)' input error in file ',FSET,' line ',I
	  stop
 7	NMSTR=NSITS(ITS(IO1(1)))
	do I=1,NMSTR
	  II=I+ISADR(ITS(IO1(1)))
	  ISREP(II)=1
	  ISM(II,1)=II
	end do
 8	continue
	DDX=2.d0*RXMAX/NOMX
	DDY=2.d0*RYMAX/NOMY
	DDZ=2.d0*RZMAX/NOMZ
	ROMAX2=(RXMAX**2+RYMAX**2+RZMAX**2)
	DDO=ROMAX/NOM     
*  Zero tables
	do IX=0,NOMXM
	   do IY=0,NOMYM
	      do IZ=0,NOMZM
		 ID3(IX,IY,IZ)=0
	      end do
	   end do
	end do
	do I=1,NS
	   XAV(I)=0.
	   YAV(I)=0.
	   ZAV(I)=0.
	   XAV2(I)=0.
	   YAV2(I)=0.
	   ZAV2(I)=0.
	   IAA(I)=0
	end do
      VOLS=0.
      NCOR=0
      IEND=0
      do while(IEND.eq.0)
        call READCONF(IEND)
	VOLS=VOLS+VOL
        if(IEND.ne.0)go to 140
	do IGR=1,NSOR   ! equivalent atoms "to"
	  do II=1,NOI   ! equivalent atom groups "from"
	    ITL=ITS(IO1(II))    ! type of "from" molecule
	    do IML=1,NSPEC(ITL)  ! cycle over "from" molecules
*   3 atom difining local coord.syst:  I, I1, I2 
	      ISHIFT=ISADDR(ITL)+(IML-1)*NSITS(ITL)-ISADR(ITL)
	      IA=ISHIFT+IO1(II)
	      IA1=ISHIFT+IO2(II)
	      IA2=ISHIFT+IO3(II)
              D1(1)=SX(IA1)-SX(IA)
              D1(2)=SY(IA1)-SY(IA)
              D1(3)=SZ(IA1)-SZ(IA)
              D2(1)=SX(IA2)-SX(IA)
              D2(2)=SY(IA2)-SY(IA)
              D2(3)=SZ(IA2)-SZ(IA)
	      E0(1)=SX(IA)               ! center of local coord. syst.
	      E0(2)=SY(IA)
	      E0(3)=SZ(IA)
	      call PBC(D1(1),D1(2),D1(3))
	      call PBC(D2(1),D2(2),D2(3))
	      if(IPRINT.ge.7)write(*,*)IA,IA1,IA2,' ',NM(IO1(II)),
     +  NM(IO2(II)),NM(IO3(II))
	      call VECT(D1,D2,RV)            ! normal v. to 123 plane  (Z axis)
	      call NORM(RV,EZ,RN,3)
	      RV(1)=D1(1)+D2(1)              ! mediana of 213 angle  ( X axis)
	      RV(2)=D1(2)+D2(2)
	      RV(3)=D1(3)+D2(3)
	      call NORM(RV,EX,RN2,3)            !  X axis
	      call VECT(EZ,EX,EY)               !  Y axis
*  Compute average structure
CC              call MEANSTR(II,IML,E0,EX,EY,EZ,IML,IML)
*  Sum over atoms which mean positions should be defined
**	do I=1,NMSTR
	      I=1
 205	      if(I.gt.NMSTR)go to 220
	      ISS=ISM(I,II)           ! site
	      ISHIFT=ISADDR(ITL)+(IML-1)*NSITS(ITL)-ISADR(ITL)
	      IS=ISS+ISHIFT
	      NMM(I)=NM(NSITE(IS))(1:2)
	      DX=SX(IS)-E0(1)
	      DY=SY(IS)-E0(2)
	      DZ=SZ(IS)-E0(3)
	      call PBC(DX,DY,DZ)
	      XL=DX*EX(1)+DY*EX(2)+DZ*EX(3)
	      YL=DX*EY(1)+DY*EY(2)+DZ*EY(3)
	      ZL=DX*EZ(1)+DY*EZ(2)+DZ*EZ(3)
	   if(ISREP(I).gt.1.and.IAA(I).gt.0)then
*   check the next atom and accumulated average
	     I1=I+1
	     ISS2=ISM(I1,II)+ISHIFT
	     DX2=SX(ISS2)-E0(1)
	     DY2=SY(ISS2)-E0(2)
	     DZ2=SZ(ISS2)-E0(3)
	     call PBC(DX2,DY2,DZ2)
	     XL2=DX2*EX(1)+DY2*EX(2)+DZ2*EX(3)
	     YL2=DX2*EY(1)+DY2*EY(2)+DZ2*EY(3)
	     ZL2=DX2*EZ(1)+DY2*EZ(2)+DZ2*EZ(3)
	     AX0=XAV(I)/IAA(I)
	     AY0=YAV(I)/IAA(I)
	     AZ0=ZAV(I)/IAA(I)
	     AX1=XAV(I1)/IAA(I1)
	     AY1=YAV(I1)/IAA(I1)
	     AZ1=ZAV(I1)/IAA(I1)
	     RR1=(XL-AX0)**2+(YL-AY0)**2+(ZL-AZ0)**2           !  1->1
	     RR2=(XL2-AX0)**2+(YL2-AY0)**2+(ZL2-AZ0)**2        !  2->1
	     RR11=(XL-AX1)**2+(YL-AY1)**2+(ZL-AZ1)**2          !  1->2
	     RR21=(XL2-AX1)**2+(YL2-AY1)**2+(ZL2-AZ1)**2       !  2->2
	     if(ISREP(I).eq.3)then
	       I2=I+2
	       ISS3=ISM(I2,II)+ISHIFT
	       DX3=SX(ISS3)-E0(1)
	       DY3=SY(ISS3)-E0(2)
	       DZ3=SZ(ISS3)-E0(3)
               if(IPROC.ne.4)call PBC(DX3,DY3,DZ3)
	       AX2=XAV(I2)/IAA(I2)
	       AY2=YAV(I2)/IAA(I2)
	       AZ2=ZAV(I2)/IAA(I2)
	       XL3=DX3*EX(1)+DY3*EX(2)+DZ3*EX(3)
	       YL3=DX3*EY(1)+DY3*EY(2)+DZ3*EY(3)
	       ZL3=DX3*EZ(1)+DY3*EZ(2)+DZ3*EZ(3)
	       RR3 =(XL3-AX0)**2+(YL3-AY0)**2+(ZL3-AZ0)**2       !  3->1
	       RR31=(XL3-AX1)**2+(YL3-AY1)**2+(ZL3-AZ1)**2       !  3->2
	       RR32=(XL3-AX2)**2+(YL3-AY2)**2+(ZL3-AZ2)**2       !  3->3
	       RR12=(XL -AX2)**2+(YL -AY2)**2+(ZL -AZ2)**2       !  1->3
	       RR22=(XL2-AX2)**2+(YL2-AY2)**2+(ZL2-AZ2)**2       !  2->3
*  evaluation of different transpositions
	       ICS=1
	       R123=RR1+RR21+RR32
	       RMIN=R123
	       R213=RR2+RR11+RR32
	       if(R213.lt.RMIN)then
		  ICS=2
		  RMIN=R213
	       end if
	       R132=RR1+RR31+RR22
	       if(R132.lt.RMIN)then
		  ICS=3
		  RMIN=R132
	       end if
	       R231=RR11+RR31+RR12
	       if(R231.lt.RMIN)then
		  ICS=4
		  RMIN=R231
	       end if
	       R312=RR3+RR11+RR22
	       if(R312.lt.RMIN)then
		  ICS=5
		  RMIN=R312
	       end if
	       R321=RR3+RR21+RR12
	       if(R321.lt.RMIN)ICS=6
	       if(ICS.eq.2.or.ICS.eq.4)then       !   2 -> 1
		  XAV(I)=XAV(I)+XL2
		  YAV(I)=YAV(I)+YL2
		  ZAV(I)=ZAV(I)+ZL2
	          XAV2(I)=XAV2(I)+XL2**2
	          YAV2(I)=YAV2(I)+YL2**2
	          ZAV2(I)=ZAV2(I)+ZL2**2
	if(IPRINT.ge.7)write(*,'(3I5,a4,3f9.3)')
     +  I,IS,ICS,NM(NSITE(IS)),XL2,YL2,ZL2
	          IAA(I)=IAA(I)+1
		  if(ICS.eq.2)then              !  1->2  3
		     I1=I+1
		     I2=I+2
		  else       !  ICS=4              3->2; 1->3
		     I1=I+2
		     I2=I+1
		  end if
		  XAV(I1)=XAV(I1)+XL
		  YAV(I1)=YAV(I1)+YL
		  ZAV(I1)=ZAV(I1)+ZL
	if(IPRINT.ge.7)write(*,'(3I5,a4,3f9.3)')
     +  I,ISS2,ICS,NM(NSITE(ISS2)),XL,YL,ZL
	          XAV2(I1)=XAV2(I1)+XL**2
	          YAV2(I1)=YAV2(I1)+YL**2
	          ZAV2(I1)=ZAV2(I1)+ZL**2
		  XAV(I2)=XAV(I2)+XL3
		  YAV(I2)=YAV(I2)+YL3
		  ZAV(I2)=ZAV(I2)+ZL3
	if(IPRINT.ge.7)
     +write(*,'(3I5,a4,3f9.3)')I,ISS3,ICS,NM(NSITE(ISS3)),XL3,YL3,ZL3
	          XAV2(I2)=XAV2(I2)+XL3**2
	          YAV2(I2)=YAV2(I2)+YL3**2
	          ZAV2(I2)=ZAV2(I2)+ZL3**2
		  IAA(I1)=IAA(I1)+1
		  IAA(I2)=IAA(I2)+1
		  I=I+2
		  go to 210
	       else if (ICS.ge.5)then   !        3 -> 1
		  XAV(I)=XAV(I)+XL3
		  YAV(I)=YAV(I)+YL3
		  ZAV(I)=ZAV(I)+ZL3
	          XAV2(I)=XAV2(I)+XL3**2
	          YAV2(I)=YAV2(I)+YL3**2
	          ZAV2(I)=ZAV2(I)+ZL3**2
	if(IPRINT.ge.7)write(*,'(3I5,a4,3f9.3)')
     +I,IS,ICS,NM(NSITE(IS)),XL3,YL3,ZL3
	          IAA(I)=IAA(I)+1
		  if(ICS.eq.5)then                !    312
		     I1=I+1
		     I2=I+2
		  else   !   ICS=6                     321
		     I1=I+2
		     I2=I+1
		  end if
		  XAV(I1)=XAV(I1)+XL
		  YAV(I1)=YAV(I1)+YL
		  ZAV(I1)=ZAV(I1)+ZL
	if(IPRINT.ge.7)write(*,'(3I5,a4,3f9.3)')
     +I,ISS2,ICS,NM(NSITE(ISS2)),XL,YL,ZL
	          XAV2(I1)=XAV2(I1)+XL**2
	          YAV2(I1)=YAV2(I1)+YL**2
	          ZAV2(I1)=ZAV2(I1)+ZL**2
		  XAV(I2)=XAV(I2)+XL2
		  YAV(I2)=YAV(I2)+YL2
		  ZAV(I2)=ZAV(I2)+ZL2
	if(IPRINT.ge.7)
     +  write(*,'(3I5,a4,3f9.3)')I,ISS3,ICS,NM(NSITE(ISS3)),XL2,YL2,ZL2
	          XAV2(I2)=XAV2(I2)+XL2**2
	          YAV2(I2)=YAV2(I2)+YL2**2
	          ZAV2(I2)=ZAV2(I2)+ZL2**2
		  IAA(I1)=IAA(I1)+1
		  IAA(I2)=IAA(I2)+1
		  I=I+2
		  go to 210
	       end if     ! (RR2.lt.RR1.and.RR2.lt.RR3
	     else if((RR1+RR21).gt.(RR2+RR11)) then  !   ISREP=2,  1 <-> 2
	       I1=I+1
	       XAV(I1)=XAV(I1)+XL
	       YAV(I1)=YAV(I1)+YL
	       ZAV(I1)=ZAV(I1)+ZL
	       if(IPRINT.ge.7)
     + write(*,'(2I5,a4,3f9.3)')I1,ISS2,NM(NSITE(ISS2)),XL,YL,ZL
	       XAV2(I1)=XAV2(I1)+XL**2
	       YAV2(I1)=YAV2(I1)+YL**2
	       ZAV2(I1)=ZAV2(I1)+ZL**2
	       XAV(I)=XAV(I)+XL2
	       YAV(I)=YAV(I)+YL2
	       ZAV(I)=ZAV(I)+ZL2
	if(IPRINT.ge.7)write(*,'(2I5,a4,3f9.3)')
     +I,IS,NM(NSITE(IS)),XL2,YL2,ZL2
	       XAV2(I)=XAV2(I)+XL2**2
	       YAV2(I)=YAV2(I)+YL2**2
	       ZAV2(I)=ZAV2(I)+ZL2**2
	       IAA(I)=IAA(I)+1
	       IAA(I1)=IAA(I1)+1
	       I=I+1
	       go to 210
	     end if   !  ISEP.eq.3
	   end if     !  ISEP.ne.1
	   XAV(I)=XAV(I)+XL
	   YAV(I)=YAV(I)+YL
	   ZAV(I)=ZAV(I)+ZL
	   if(IPRINT.ge.7)write(*,'(2I5,a4,3f9.3)')
     +I,IS,NM(NSITE(IS)),XL,YL,ZL
	   XAV2(I)=XAV2(I)+XL**2
	   YAV2(I)=YAV2(I)+YL**2
	   ZAV2(I)=ZAV2(I)+ZL**2
	   IAA(I)=IAA(I)+1
 210	   I=I+1
	   go to 205
 220	   continue
C   End of calculation of the mean structure   
	      JTYP=ITS(ISOR(IGR)) ! type of "to" molecule
	      NBG=ISADDR(JTYP)+ISOR(IGR)-ISADR(JTYP)
	      NEN=ISADDR(JTYP+1)
	      NST=NSITS(JTYP)
	      do J=NBG,NEN,NST        ! over molecules with given site
	        if(NNUM(J).eq.NNUM(IA))go to 120
	        DX=SX(J)-E0(1)
	        DY=SY(J)-E0(2)
	        DZ=SZ(J)-E0(3)
	        call PBC(DX,DY,DZ)
	        RR2=DX**2+DY**2+DZ**2
                NCOR=NCOR+1
	        if(RR2.lt.ROMAX2)then
	          RR=sqrt(RR2)
                  XM=EX(1)*DX+EX(2)*DY+EX(3)*DZ       ! X-coord in molec.c.s. 
                  YM=EY(1)*DX+EY(2)*DY+EY(3)*DZ       ! Y-coord in molec.c.s. 
                  ZM=EZ(1)*DX+EZ(2)*DY+EZ(3)*DZ       ! Z-coord in molec.c.s.
*  calculate 3D density 
	          NXM=(XM+RXMAX)/DDX
	          NYM=(YM+RYMAX)/DDY
	          NZM=(ZM+RZMAX)/DDZ
	          if(NXM.ge.0.and.NXM.le.NOMX.and.NYM.ge.0.
     + and.NYM.le.NOMY.and.NZM.ge.0.and.NZM.le.NOMZ)then
     	            ID3(NXM,NYM,NZM)=ID3(NXM,NYM,NZM)+1
		    if(IPRINT.ge.7)write(*,*)I,J,NXM,NYM,NZM,XM,YM,ZM
	          end if
	        end if     !  RR2.lt.ROMAX
                if(NCOR.le.0)then
	          write(*,*)' Integer overflow!!!'
	          stop
	        end if
 120	        continue
	      end do     !  J=NBG    
 	    end do      !  IML=
	  end do       ! II
	end do       ! IGR
	if(IPRINT.ge.6)write(*,*)' time=',FULTIM,' NCOR:',NCOR
      end do !  trajectory
 140  continue
      VOL=VOLS/IAN
C  SDF output	
      open(unit=31,file=filori,status='unknown')
      X0=-RXMAX+DDX			! SDF box coord
      X1= RXMAX
      Y0=-RYMAX+DDY
      Y1= RYMAX
      Z0=-RZMAX+DDZ
      Z1= RZMAX
C  This is "formatted" gOpenMol PLT-file 
      N3=3
      N200=200
      write(31,*)N3,N200
      write(31,*)NOMZ,NOMY,NOMX
      write(31,'(6e13.5)')Z0,Z1,Y0,Y1,X0,X1
      FACORI=VOL/(DDX*DDY*DDZ*NCOR)
      write(*,*)' writing 3D distribution ...'
      write(*,'(i14,3f8.2,3f8.4,e13.5)')
     +NCOR,RXMAX,RYMAX,RZMAX,DDX,DDY,DDZ,FACORI
      do IZ=1,NOMZ
	do IY=1,NOMY
          do IX=1,NOMX
            TMP(IX)=ID3(IX,IY,IZ)*FACORI
	  end do
	  write(31,'(6e13.5)')(TMP(IX),IX=1,NOMX)
	end do
      end do
      close(31)
C  Mean structure output (XMOL format)
      open(unit=32,file=filcrd,status='unknown')
      write(32,'(i5,a)')NMSTR,
     &'      average structure from MolDynamix'
      write(32,'(a)')
     +'At        x       y        z     ...        RMS'
      do I=1,NMSTR
	XX=XAV(I)/IAA(I)
	XX2=abs(XAV2(I)/IAA(I)-XX**2)
	YY=YAV(I)/IAA(I)
        YY2=abs(YAV2(I)/IAA(I)-YY**2)
	ZZ=ZAV(I)/IAA(I)
        ZZ2=abs(ZAV2(I)/IAA(I)-ZZ**2)
	RR=sqrt(XX2+YY2+ZZ2)
	write(32,'(a5,3(1x,f10.4),4x,f9.5)')NMM(I),XX,YY,ZZ,RR
      end do
      close(32)
      write(*,*)' Done!   '
      write(*,*)' Use pltfile utility to convert file ',FILORI,
     +' to binary .plt format'
      stop
      end
