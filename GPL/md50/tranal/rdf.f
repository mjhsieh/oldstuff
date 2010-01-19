C    Package tranal v.5.0
C    Calculations of radial distribution functions
      program RDF
      include "tranal.h"
      character*64 FOUTRDF
      parameter (MAXRDF=50,MAXGR=50,NAAM=1000)
      integer IREP(MAXRDF),IRD(MAXRDF,MAXGR),JRD(MAXRDF,MAXGR),
     +IRDF(MAXRDF,NAAM),IRDFI(MAXRDF,NAAM)
      dimension TMP(MAXRDF)
      namelist /RDFIN/ FOUTRDF,RDFCUT,NA,NAI,RMI,RMAX,NRDF
*  defaults
      RDFCUT=10.
      NA=200
      NAI=100
      RMI=0
      RMAX=5.
*  read and setup trajectory parameters
      call SETTRAJ
*  read RDF computation parameters
      read(*,RDFIN)
      if(NRDF.gt.MAXRDF)stop 'Increase MAXRDF!!!'
      if(NA.gt.NAAM.or.NAI.gt.NAAM)stop 'Increase NAAM!'
      DRI=RMAX-RMI
*  initialization
      BOXS=0.           ! to accumulate average box sizes
      BOYS=0.
      BOZS=0.
*  read RDFs
      N=0
      write(*,*)' Trying to read ',NRDF,' groups for RDF calculations'
      DO K       = 1,NRDF
 	IE=-1
	STR=TAKESTR(5,IE)
	if(str(1:1).eq.'&')then
          READ(str(2:24),*)IRP
	  if(IRP.gt.MAXGR) stop 'Increase MAXGR!!!'
	  IREP(K)=IRP
	  do IR=1,IRP
	    STR=TAKESTR(5,IE)
	    read(STR,*,err=99)I,J  
	    IRD(K,IR)=I
	    JRD(K,IR)=J
	  end do
	else 
	  IREP(K)=1
	  read(STR,*,err=99)I,J  
	  IRD(K,1)=I
	  JRD(K,1)=J
	end if
	write(*,*)' RDF  No ',K,' equiv ',IREP(K)
	do I=1,IREP(K)
	  write(*,'(2I6,2(3x,a4))')
     +    IRD(K,I),JRD(K,I),NM(IRD(K,I)),NM(JRD(K,I))
	end do
      END DO
      go to 101
 99   STR=TAKESTR(5,99)
      stop
 101  continue
      do I=1,NA
	do J=1,NRDF
	  IRDF(J,I)=0
	end do
      end do
      do I=1,NAI
	do J=1,NRDF
	  IRDFI(J,I)=0
	end do
      end do
*  Starting trajectory analysis
      IEND=0
      do while(IEND.eq.0)              ! over configurations
        call READCONF(IEND)
        if(IEND.ne.0)go to 200
        if(IPRINT.ge.6)write(*,'(a,f12.3,a,3f10.3)')
     +  ' Processing time',FULTIM
        do I=1,NRDF                    !  over RDF groups
	  do K=1,IREP(I)               !  over site pairs of the same group
	    IS=IRD(I,K)
	    JS=JRD(I,K)
	    ITYP=ITS(IS)
	    JTYP=ITS(JS)
	    NIS=NSPEC(ITYP)
	    NJS=NSPEC(JTYP)
	    do IM=1,NIS                ! over first molecules
	      JBEG=1
	      if(IS.eq.JS)JBEG=IM+1
	      ISHIFT=ISADDR(ITYP)+IS-ISADR(ITYP)-NSITS(ITYP)
	      do JM=JBEG,NJS           ! over second molecules
	        JSHIFT=ISADDR(JTYP)+JS-ISADR(JTYP)-NSITS(JTYP)
	        IST=ISHIFT+IM*NSITS(ITYP)
	        JST=JSHIFT+JM*NSITS(JTYP)
	        IMOL=NNUM(IST)
	        JMOL=NNUM(JST)
*		  write(*,*)I,K,IS,JS,IST,JST,NM(IS),NM(JS)
	        DX=SX(IST)-SX(JST)
	        DY=SY(IST)-SY(JST)
	        DZ=SZ(IST)-SZ(JST)
	        call PBC(DX,DY,DZ)
	        RR=sqrt(DX**2+DY**2+DZ**2)
*  This distinguishes sites on the same molecule
	        if(IMOL.ne.JMOL)then
		  NR=NA*RR/RDFCUT+1
		  if(NR.gt.0.and.NR.le.NA)IRDF(I,NR)=IRDF(I,NR)+1
	        else
		  NR=NAI*(RR-RMI)/DRI+1
*		write(*,'(7i6,f8.4)')I,K,IS,JS,IST,JST,NR,RR
		  if(NR.gt.0.and.NR.le.NAI)IRDFI(I,NR)=IRDFI(I,NR)+1
	        end if
	      end do
	    end do
	  end do
        end do
        BOXS=BOXS+BOXL
        BOYS=BOYS+BOYL
        BOZS=BOZS+BOZL
      end do
 200  write(*,*)IAN,' configurations analysed'
*  RDF output
      open(unit=16,file=FOUTRDF,status='unknown')
      write(16,'(a)')'   Radial distribution functions '
*   Averagre box sizes and volume
      BOXS=BOXS/IAN
      BOYS=BOYS/IAN
      BOZS=BOZS/IAN
      VOLS = BOXS*BOYS*BOZS
      write(16,'(a,3f12.4)')'  Average box sizes: ',BOXS,BOYS,BOZS
      IFST=NA
      do K=1,NRDF
	CONU=0.
	NPAIR=0
	NPRI=0                   ! intramolecular. pairs
	do KK=1,IREP(K)
	  IS=IRD(K,KK)
	  JS=JRD(K,KK)
	  ITYP=ITS(IS)
	  JTYP=ITS(JS)
	  NN1 = NSPEC(ITYP)
	  NN2 = NSPEC(JTYP)
	  if(ITYP.eq.JTYP)then
	    if(IS.eq.JS)then
	      NPAIR=NPAIR+NN1*(NN1-1)/2
	    else
	      NPAIR=NPAIR+NN1*(NN1-1)
	    end if
	  else
	    NPAIR=NPAIR+NN1*NN2
	  end if
	  if(ITYP.eq.JTYP.and.IS.ne.JS)NPRI=NPRI+1
	end do
        write(16,*)
	FACT=VOLS/(NPAIR*IAN)
	write(16,*)' This pair: ',NM(IS),' - ',NM(JS),
     +  '   Npairs=',NPAIR,'  Nint=',NPRI
	write(16,*)'------------------------------------'
        CONU=0.
	do I=1,NA
	  RDFV=IRDF(K,I)
	  RR=(I-0.5)*RDFCUT/NA
	  DRR=RDFCUT/NA
   	  SH12=4*PI*DRR*(RR**2+DRR**2/12.d0)
	  CONU=CONU+RDFV*FACT/VOLS
	  CONU1=CONU*NSPEC(ITS(IS))
	  CONU2=CONU*NSPEC(ITS(JS))
  	  if(CONU.gt.1.d-8)then
	    write(16,'(f8.3,3f12.5,5x,a4,2x,a4,a1,a4,I15)')
     +RR,RDFV*FACT/SH12,CONU1,CONU2,'rdf:',NM(IS),'-',NM(JS),IRDF(K,I)
	    if(I.lt.IFST)IFST=I
	   end if
	 end do
	 if(NPRI.ge.1)then
	   write(16,*)
	   write(16,*)' Intramolecular RDF: '
	   write(16,*)' ------------------- '
 	   FACT=1.d0/(NSPEC(ITS(IS))*NPRI*IAN)
	   do I=1,NAI
	     RR=RMI+(I-0.5)*DRI/NAI
	     DRR=DRI/NAI
	     RDFV=IRDFI(K,I)*FACT/DRR
	     if(RDFV.gt.1.d-8)then
	      write(16,'(f8.3,f12.5,25x,a8,2x,a4,a1,a4)')
     +     RR,RDFV,'rdf_int:',NM(IS),'-',NM(JS)
	     end if
	   end do
	 end if
      end do
      write(16,*)' SUMMARY:'
      write(16,'(a1,7x,50(1x,a4,a1,a4))')
     +  '#',(NM(IRD(K,1)),'-',NM(JRD(K,1)),K=1,NRDF)
      do I=IFST,NA
	RR=(I-0.5)*RDFCUT/NA
	DRR=RDFCUT/NA
  	SH12=4*PI*DRR*(RR**2+DRR**2/12.d0)
        do K=1,NRDF
	 NPAIR=0
	 do KK=1,IREP(K)
	   IS=IRD(K,KK)
	   JS=JRD(K,KK)
	   if(IS.eq.JS)then
	     NPAIR=NPAIR+NSPEC(ITS(IS))*(NSPEC(ITS(JS))-1)/2
	   else
	     NPAIR=NPAIR+NSPEC(ITS(IS))*NSPEC(ITS(JS))
	   end if
	 end do
	 FACT=VOLS/(NPAIR*1.d0*IAN)
	 TMP(K)=IRDF(K,I)*FACT/SH12
	end do
	write(16,'(f8.3,50f10.5)')RR,(TMP(K),K=1,NRDF)
      end do
      close(16)
      stop
      end
