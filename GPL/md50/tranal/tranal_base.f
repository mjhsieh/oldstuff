*  
*   Part of MDynaMix package v.5.0 
*
*   Block TRANAL_BASE - analysis of trajectorisproceedies
*
*   subroutines:
*        SETTRAJ - prepare list of trajectories and data structures
*        READCONF - reads next configuration
*        READMOL - read information about molecules from .mmol files
*        GETCOM - compute molecular centres of mass
*        PBC - apply periodic boundary conditions
*
*   This subroutines prepare data structures 
*   for subsequent analysis
*
 	subroutine SETTRAJ
	include "tranal.h"
	character*128 FFILTR,PATHDB,FNAME,cdum*8,chr4*4
	real*4 TX(NTOT),TY(NTOT),TZ(NTOT)
	save
	namelist /TRAJ/NFORM,FNAME,PATHDB,NTYPES,NAMOL,NFBEG,NFEND,
     + ISTEP,BREAKM,LVEL,LFILTR,FFILTR,LXMOL,FXMOL,IPRINT,LBOND,NSPEC,
     + BOXL,BOYL,BOZL,VFAC
*  requested parameters:
*
*    NFORM - trajectory format
*        MDYN     - binary MDynaMix
*        XMOL     - XMOL  (xyz)
*        PDBT     - PDB trajectory from GROMACS
*        DCDT     - DCD binary trajectory from NAMD
*    
*    PATHDB - path to molecular database (.mmol files)
*    NTYPES - number of molecular types
*    NAMOL - names of molecules to analyse
*    NFBEG -  number of the first trajectory file 
*    NFEND - number of the last trajectory file
*    NSPEC - number of molecules of each type (not necessary in MDYN format)
*
*  other things:
*       
*  BOXL,BOYL,BOZL - box sizes. They are used if trajectory file does not 
*                   contain information on them
*
*  IPRINT - output level
*  BREAKM - which break in trajectory file (in ps ) is allowed
*
*   LVEL - read velocities from the trajectory files
*
*   LXMOL - eventually rewrite trajectory to XMOL format
*   FXMOL - base name of XMOL file
*   ISTEP - how often take configurations
*
*   LFILTR - analysis for particular atoms
*   FFILTR - file name specifying which atoms take for the analysis
*
*   LBOND - information about bonds is read from .mmol file
*
*   VFAC - factor to multiply velocities to bring them to A/fs
*
*  This is default trajectory format
	NFORM='MDYN'
	LVEL=.false.
	LFILTR=.false.
	LXMOL=.false.
	LBOND=.false.
	ISTEP=1
	IPRINT=5
	BREAKM=11.
	VFAC=1.
*  Read input data
	read(*,TRAJ)
	if(NFORM.eq.'PDBT'.or.NFORM.eq.'DCDT')then
	  if(LVEL)write(*,*)' no velocities in this trajectory format'
          LVEL=.false.
	end if
*  defining the length of the file name
	IL=0 
	do I=1,120
	  ISTR=ichar(FNAME(I:I))
	  if(ISTR.le.15)then
	    LF=I-1
	    go to 11
	  end if
        if(FNAME(I:I).ne.' ')IL=I
	end do 
	LF=IL
*  Prepare for eventual rewriting trajectory files in XMOL format
 	if(LXMOL)then
 	  LLNX=0 
 	  do I=1,128
 	    ISTR=ichar(FXMOL(I:I))
	    if(ISTR.le.15)then
	      LLNX=I-1
	      go to 11
	    end if
            if(FXMOL(I:I).ne.' ')LLNX=I
	  end do 
	end if
*  Check and reconstruct file list
 11	continue
	write(*,*)FNAME(1:LF)
	NF=NFEND-NFBEG+1
	J=0
 	do I=1,NF
	  NFL=NFBEG+I-1
	  filnam(I)(1:LF)=FNAME(1:LF)
	  filnam(I)(LF+1:LF+4)='.000'
	  do M=LF+5,128
	    filnam(I)(M:M)=' '
	  end do
	  if(NFL.lt.10)then
	    write(filnam(I)(LF+4:LF+4),'(i1)')NFL
	  else if(NFL.lt.100)then
	    write(filnam(I)(LF+3:LF+4),'(i2)')NFL
	  else
	    write(filnam(I)(LF+2:LF+4),'(i3)')NFL
	  end if 
	  write(*,*)' looking for file ',filnam(I)
	  if(NFORM.eq.'PDBT')then
            open (unit=12,file=filnam(I),status='old',err=5)
  	    if(IPRINT.ge.6)write(*,'(a6,a64,a8)')
     +                 ' file ',filnam(I)(1:64),' is open'
 9	    read(12,'(a)',err=10)STR
*     PDB trajectory format
	    if(STR(1:6).ne.'HEADER'.and.STR(1:5).ne.'TITLE')go to 9
	    do K=7,128
	      if(STR(K:K+1).eq.'t=')then
	        read(str(K+2:K+20),*)FULTIM
  	        if(IPRINT.ge.5)write(*,*)' FULTIM ',FULTIM,' ps' 
		go to 13
	      end if
	    end do
	    write(*,*)' time is not defined '
 13	    continue
	    go to 15 
*  XMOL trajectory format
	  else if(NFORM.eq.'XMOL')then
            open (unit=12,file=filnam(I),status='old',err=5)
  	    if(IPRINT.ge.6)write(*,'(a6,a64,a8)')
     +                 ' file ',filnam(I)(1:64),' is open'
	    read(12,*,err=10)NSTOT
	    read(12,*,err=14)cdum,FULTIM
	    FULTIM=0.001*FULTIM            ! in ps
	    if(IPRINT.ge.5)write(*,*)' FULTIM ',FULTIM,' ps'  
 	    go to 15
 14	    write(*,*)'time is not defined'
	    FULTIM=I
	    go to 15
*  Binary MDynaMix format
	  else if(NFORM.eq.'MDYN')then
            open (unit=12,file=filnam(I),status='old',
     +    form='unformatted',err=5)
  	    if(IPRINT.ge.6)write(*,'(a6,a64,a8)')
     +                 ' file ',filnam(I)(1:64),' is open'
            read(12,err=10,end=10)IA,dt,boxl,boyl,bozl,unitl,ntypes
	    if(IA.eq.1)LVEL=.true.
	    if(IPRINT.ge.5)write(*,*)' num.of types ',NTYPES  
	    if(I.gt.1)then
	      if(NTP.ne.NTYPES)then
	        write(*,*)
     +      ' uncompatible NTYPES in ',filnam(I),' and ',filnam(I-1) 
	        stop
	      end if
	    else
	      rewind (12)
              read(12,err=10) IA, dt, boxl, boyl, bozl, unitl, ntypes,
     x                (nspec(ii),nsits(ii),logo,ii=1,ntypes)
	    end if 
	    NTP=NTYPES  
	    do ITYP=1,NTYPES
	      read(12,err=10)
	      if(IPRINT.ge.8)write(*,*)ITYP,NSPEC(ITYP),NSITS(ITYP)
	    end do
	    read(12)IA,FULTIM
	    FULTIM=FULTIM*1.d12         ! in ps
	    if(IPRINT.ge.5)write(*,*)' FULTIM ',FULTIM,' ps'  
	    go to 15
*  NAMD DCD format
	  else if(NFORM.eq.'DCDT')then
            open (unit=12,file=filnam(I),status='old',
     +    form='unformatted',err=5)
  	    if(IPRINT.ge.6)write(*,'(a6,a64,a8)')
     +                 ' file ',filnam(I)(1:64),' is open'
            read(12)CHR4,NCF,NST1,NFREQ,NSTF,I1,I2,I3,I4,I5,DTN
	    write(*,*) CHR4,NCF
	    FULTIM=NST1*DTN*0.04888            ! in ps
	    if(IPRINT.ge.5)write(*,*)' FULTIM ',FULTIM,' ps'  
	    go to 15 
	  else
	    write(*,*)' Unsupported trajectory format ',NFORM
	    stop
	  end if
  5	  write(*,*)' File not found: ',filnam(I)
	  go to 20
 10	  write(*,*)' File ',filnam(I),': input error'
	  go to 20
 15	  J=J+1
	  TIMIN(J)=FULTIM 
	  filnam(J)=filnam(I)
 20	  close(12)
	end do
        write(*,*)' input of the file list is completed:',J,' files'	  
 30	NF=J
	FULTIM0=TIMIN(1)
 	TTER=TTER*1.d-12    
	if(NTYPES.le.0)then
	   write(*,*)'!!! Wrong number of molecule types : ',NTYPES
	   stop
	end if
*  Adresses and referenses
	IAN=0                  ! number of analysed configurations
	ISADR(1)=0
	ISADDR(1)=0    
	IADDR(1)=0
	NOP=0 
	NSITES=0
	NSTOT=0 
	FNSTR=0.
	IATOM=0
	IMOL=0
	IST=0   
	NCOR=0    
	JBREAK=1                ! indicate break of trajectory
	NX=0
*  Read molecular information
	do ITYP=1,NTYPES
	  LIST(ITYP)=1
	  if(NFORM.eq.'MDYN')NSTS=NSITS(ITYP)
          CALL READMOL (SUMM,NX,ITYP,PATHDB,NAMOL(ITYP))
	  if(NFORM.eq.'MDYN'.and.NSTS.ne.NSITS(ITYP))then
	    write(*,*)' Molecule of type ',ITYP,' has ',NSTS,
     &' sites in the trajectory but ',NSITS(ITYP),' in the input file'
	    stop
	  end if
          NOP     = NOP    + NSPEC(ITYP)
          TOTMAS       = TOTMAS+SUMM*NSPEC(ITYP)
          SUMMAS(ITYP) = SUMM 
	  NAME(ITYP)=NAMOL(ITYP)(1:6)
        end do           
	if(NOP.le.0)then
	   write(*,*)'!!! Wrong number of molecules :      ',NOP
	   stop
	end if
        FACTM=AVSNO*1.d3
        do ityp=1, ntypes
 	  ISADR(ITYP+1)  = ISADR(ITYP)+ NSITS(ITYP)          ! site addr
	  ISADDR(ITYP+1) = ISADDR(ITYP)+NSPEC(ITYP)*NSITS(ITYP)  ! atom adr
          IADDR (ITYP+1) = IADDR (ITYP)+NSPEC(ITYP)         ! mol. addr
          NSITES  = NSITES + NSITS(ITYP)
          if(NSITS(ITYP).eq.2)FNSTR=FNSTR+2.d0*NSPEC(ITYP)
          if(NSITS(ITYP).ge.3)FNSTR=FNSTR+3.d0*NSPEC(ITYP)
          NSTOT   = NSTOT  + NSPEC(ITYP)*NSITS(ITYP)
          FNSS3   = 3.D0*DFLOAT(NSITES)
          write(*,'(a6,i4,3x,a5,i5,4x,a6,i4)')
     +' type ',ITYP,' Nmol',NSPEC(ITYP),
     +'Nsites',NSITS(ITYP)
	  do J=ISADR(ITYP)+1,ISADR(ITYP+1)
	    MASSD(J)=MASS(J)/TOTMAS/FACTM
	    ITS(J)=ITYP
	  end do                     !    Nsite -> type
	  do IS=1,NSPEC(ITYP)
	    IMOL=IMOL+1
            do I=1,NSITS(ITYP)
	      ISITE=ISADR(ITYP)+I
	      IATOM=IATOM+1
	      NNUM(IATOM)=IMOL        !   Natom -> Nmol
	      NSITE(IATOM)=ISITE      !   Natom -> Nsite
	      ITYPE(IATOM)=ITYP       !   Natom -> type
	      if(IPRINT.ge.9)write(*,'(4I6,2x,2a4,a6)')
     +        IATOM,IMOL,ISITE,ITYP,NM(ISITE)
	    end do
	    ITM(IMOL)=ITYP            !  mol -> type
	  end do
        end do                   
* Renumeration array
	if(LFILTR)then
	   do I=1,NUMTR
	     NUMR(I)=0
	   end do
	   write(*,*)' renumerating atoms from file ',FFILTR
	   open(unit=14,file=FFILTR,status='old')
	   do I=1,NSTOT
	      read(14,*,err=921,end=921)II,NUMO(II)
	      NUMR(NUMO(II))=II
	      if(IPRINT.ge.7)write(*,*)I,II,NUMO(II)
	   end do
	   close(14)
	else
	   NUMTR=NSTOT
	   do I=1,NSTOT
	      NUMO(I)=I
	      NUMR(I)=I
	   end do	   
	end if
	FULTIM1=TIMIN(1)
	write(*,*)' FULTIM1=',FULTIM1
	return
 921	write(*,*)' error reading renumeration file ',FFILTR
	stop
	end
*
*============= READCONF ===========================
*
*  read configuration
*    IEND signals the end of trajectories
	subroutine READCONF(IEND)
	include "tranal.h"
	real*4 TX(NTOT),TY(NTOT),TZ(NTOT)
	character*32 cdum,chr4*4
	data IINIT/0/ 
	save
	if(IINIT.eq.0)then
	  IINIT=1
	  IFL=1    ! number of the first file to read
 	  IIN=0     ! Signals that we read headers
	end if
	IEND=0    ! signals end of the analysis
       	ITRI=0    ! signals about double point
*  Begin reading of the next trajectory file
*  Reading HEADERs
*  Case of PDB trajactory
 35	if(NFORM.eq.'PDBT')then
	  if(IIN.eq.0)then
            open (unit=12,file=filnam(IFL),status='old',
     +      form='formatted') 
	    if(IPRINT.ge.7)write(*,*)' file ',IFL,' open'
	  end if
 31	  read(12,'(a)',end=45) STR
	  if(STR(1:6).ne.'HEADER'.and.STR(1:5).ne.'TITLE')go to 31 
	  do K=7,128
	    if(STR(K:K+1).eq.'t=')then
	      read(str(K+2:K+20),*)FULTIM
	      go to 33
	    end if
	  end do
	  write(*,*)' time is not defined '
 33	  IFLR=0
* Line 2 
 32	  read( 12,'(a)',end=45 ) STR
	  if(STR(1:4).eq.'ATOM')then
	    if(IPRINT.ge.6)write(*,*)' box sizes undetermined'
	    IFLR=1
	    go to 34
	  end if
	  if(STR(1:4).ne.'CRYS')go to 32
	  read(str,*)cdum,BOXL,BOYL,BOZL
	  if(IPRINT.ge.7)write(*,'(a,3f12.3)')' Box: ',BOXL,BOYL,BOZL
	else if(NFORM.eq.'XMOL')then
*  XMOL trajectory
	  if(IIN.eq.0)then
            open (unit=12,file=filnam(IFL),status='old',
     +    form='formatted') 
	    if(IPRINT.ge.7)write(*,*)' file ',IFL,' open'
	  end if
 	  read(12,*,err=40,end=45) NUMTR
          read( 12,'(a)',err=40,end=40 )STR
          read(STR,*,end=40,err=36)cdum,FULTIM
	  do J=1,100
            if(STR(J:J+2).eq.'BOX')then
              read(STR(J+4:128),*,end=40,err=40)BOXL,BOYL,BOZL
	      if(IPRINT.ge.7)write(*,'(a,3f12.3)')
     &        ' Box: ',BOXL,BOYL,BOZL
 	      go to 37
	    end if
          end do
	  if(IPRINT.ge.7)write(*,*)' box sizes undetermined'
	  go to 37
 36	  if(IPRINT.ge.6)write(*,*)' time is undetermined'
	  FULTIM=IAN
 37	  FULTIM=0.001*FULTIM	! in ps
	else if(NFORM.eq.'DCDT')then
*  NAMD DCD trajectory
	  if(IIN.eq.0)then
            open (unit=12,file=filnam(IFL),status='old',
     +    form='unformatted')
            read(12)CHR4,NCF,NST1,NFREQ,NSTF,I1,I2,I3,I4,I5,DTN
	    FULTIM=NST1*DTN*0.04888            ! in ps
 	    read(12)
	    read(12)NUMTR
	  else
	    FULTIM=FULTIM+NFREQ*DTN*0.04888
	    end if
	    read(12,err=40,end=45)BOXL,Y1,BOYL,Y2,Y3,BOZL
	    if(IPRINT.ge.7)write(*,'(a,3f12.3)')' Box: ',BOXL,BOYL,BOZL
	    if(LXMOL)then
	      filxm(1:LLNX)=FXMOL(1:LLNX)
	      filxm(LLNX+1:LLNX+4)=filnam(IFL)(LF+1:LF+4)
	      open(unit=19,file=filxm,status='unknown')
	    end if
	  else
*  Binary MDynamix trajectory
	    if(IIN.eq.0)then
              open (unit=12,file=filnam(IFL),status='old',
     +        form='unformatted',err=45) 
              read(12) IVEL, dt, boxl, boyl, bozl, unitl, ntypes,
     x                (nspec(i),nsits(i),logo,i=1,ntypes)
	      if(LXMOL)then
	        filxm(1:LLNX)=FXMOL(1:LLNX)
	        filxm(LLNX+1:LLNX+4)=filnam(IFL)(LF+1:LF+4)
	        open(unit=19,file=filxm,status='unknown')
	      end if
	      if(IPRINT.ge.7)write(*,*)
     &          'read first line.',NTYPES,' types'
              do ityp=1, ntypes
                NSB = ISADR (ITYP)+1
                NSE = ISADR (ITYP+1)
* Lines from 2 to NTYPES+1
                read( 12 ) (mass(j), j=NSB,NSE),LIST(ITYP)
		if(IPRINT.ge.8)write(*,*)ityp,NSB,NSE,LIST(ITYP)
              end do                    
	      if(IPRINT.ge.7)write(*,*)' read masses   '
	    end if
* First line of a configureation
            read(12,err=40,end=45)IA,fultim,cumtim,
     &      TEMP,PRES,EP,BOXL,BOYL,BOZL,(LIST(LL),LL=1,NTYPES)
	    FULTIM=FULTIM*1.d12                 ! in ps
	    if(IPRINT.ge.7)write(*,'(2(a,3f10.2))')
     &      'Box: ',BOXL,BOYL,BOZL,' term:',TEMP,PRES,EP
	  end if    ! if(NFORM...
34	  if(IIN.eq.0)then
 	    write(*,'(2a/a,f14.3,a8,6I3)')
     +  ' begin analys of file  ',filnam(IFL), 
     +  ' init. T = ',fultim
	    IIN=1
	  end if  
	  IBREAK=0
	  DIFT=FULTIM-FULTIM1
  	  HBOXL=0.5*BOXL
	  HBOYL=0.5*BOYL
	  HBOZL=0.5*BOZL
	  VOL=BOXL*BOYL*BOZL 
	  if(DIFT.gt.BREAKM)then
	    write(*,*)
     +    'break in trajectory file. from ',FULTIM1,' to ',
     +    FULTIM,'ps'
	    IBREAK=1
	  end if    
	  if(IFL.lt.NF.and.fultim.ge.TIMIN(IFL+1))go to 45
	  FULTIM1=FULTIM 
	  if(IPRINT.ge.7)write(*,*)' reading conf at time ',FULTIM
	  if(NFORM.eq.'MDYN')then
c - for each type
            do ityp=1,ntypes
              NSB          = ISADDR(ITYP)+1
              NSE         = ISADDR(ITYP+1) 
* Other NTYPES lines of each point
              if(LIST(ITYP).ne.0)then
	        if(IVEL.le.0)then
	          read(12,err=40,end=45) (TX(j),j=NSB,NSE),
     *                (TY(j),j=NSB,NSE),
     *                (TZ(j),j=NSB,NSE)           
 	          do J=NSB,NSE
	            SX(J)=TX(J)
	            SY(J)=TY(J)
	            SZ(J)=TZ(J)
	          end do
	        else    !  IVEL=1
	          read(12,err=40,end=45) (TX(j),j=NSB,NSE),
     *                (TY(j),j=NSB,NSE),
     *                (TZ(j),j=NSB,NSE)
 	          do J=NSB,NSE
	            SX(J)=TX(J)
	            SY(J)=TY(J)
	            SZ(J)=TZ(J)
	           end do
	           read(12,err=40,end=45) (TX(j),j=NSB,NSE),
     *                (TY(j),j=NSB,NSE),
     *                (TZ(j),j=NSB,NSE)
 	           do J=NSB,NSE
	             VX(J)=TX(J)
	             VY(J)=TY(J)
	             VZ(J)=TZ(J)
	           end do
	         end if    !  IVEL
	       end if     ! LIST
             end do    !  ITYP
	  else if(NFORM.eq.'DCDT')then
	    read(12,err=40,end=45)(TX(I),I=1,NUMTR)
	    read(12,err=40,end=45)(TY(I),I=1,NUMTR)
	    read(12,err=40,end=45)(TZ(I),I=1,NUMTR)
            do i=1,NUMTR
	      OX(I)=SX(I)
	      OY(I)=SY(I)
	      OZ(I)=SZ(I)
	      SX(I)=TX(I)
	      SY(I)=TY(I)
	      SZ(I)=TZ(I)
	    end do
	  else     ! ASCII trajectories
            do i=1,NUMTR
	      OX(I)=SX(I)
	      OY(I)=SY(I)
	      OZ(I)=SZ(I)
	      IN=NUMR(I)
	      if(LVEL.and.NFORM.eq.'XMOL')then
	        read(12,*,err=40,end=40)CDUM,XX,YY,ZZ,VVX,VVY,VVZ
		if(IN.ne.0)then
		   SX(IN)=XX
		   SY(IN)=YY
		   SZ(IN)=ZZ
		   VX(IN)=VVX
		   VY(IN)=VVY
		   VZ(IN)=VVZ
		   if(IPRINT.ge.8)write(*,'(2i5,3f10.4)')IN,I,XX,YY,ZZ
		end if
              else
*   GROMACS PDB trajectory
		if(NFORM.eq.'PDBT')then
 53		  read(12,'(a)',err=40)STR
		  if(STR(1:6).eq.'ENDMDL'.or.STR(1:6).eq.'HEADER')then
		    write(*,*)' Number of atoms in the trajectory ',IN-1,
     &		    ' < than in the input file ',NSTOT
		    stop
		  end if
		  if(STR(1:4).ne.'ATOM')go to 53
		  read(str(31:128),*)XX,YY,ZZ
		else        ! XMOL case
	          read(12,*,err=40,end=40)CDUM,XX,YY,ZZ
		end if		
		if(IN.ne.0)then
		   SX(IN)=XX
		   SY(IN)=YY
		   SZ(IN)=ZZ
		   if(IPRINT.ge.8)write(*,'(2i5,3f10.4)')IN,I,XX,YY,ZZ
		end if
	      end if
	  end do
	end if
	JBREAK=0
	if(ITRI.ne.0.and.DIFT.le.0.00001)then
	   write(*,*)' double point in the trajectory file at ',FULTIM
	   go to 35
	end if
	ITRI=1
	DT=DIFT
	if(DIFT.le.0.0001)DIFT=0.0001
	EKIN=0.
	do I=1,NSTOT
	  DX=SX(I)-OX(I)
	  DY=SY(I)-OY(I)
	  DZ=SZ(I)-OZ(I)
 	  if(DX.gt.HBOXL)DX=DX-BOXL
 	  if(DX.lt.-HBOXL)DX=DX+BOXL
 	  if(DY.gt.HBOYL)DY=DY-BOZL
 	  if(DY.lt.-HBOYL)DY=DY+BOZL
 	  if(DZ.gt.HBOZL)DZ=DZ-BOZL
 	  if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	  IS=NSITE(I)
	  if(.not.LVEL)then
	    VX(I)=DX*VFAC*1.d-3/DIFT            !  in A/fs
	    VY(I)=DY*VFAC*1.d-3/DIFT
	    VZ(I)=DZ*VFAC*1.d-3/DIFT
	  end if
	  EKIN=EKIN+MASS(IS)*(VX(I)**2+VY(I)**2+VZ(I)**2)   ! 2 Ekin
	  OX(I)=SX(I)
	  OY(I)=SY(I)
	  OZ(I)=SZ(I)
	end do
	TEMP=EKIN/(3.d-7*AVSNO*NSTOT*BOLTZ)
	IEND=0
	if(mod(IAN,ISTEP).eq.0.and.LXMOL)then
	   NSAT=0
	   do ITYP=1,NTYPES
	     if(LIST(ITYP).ne.0)NSAT=NSAT+NSITS(ITYP)*NSPEC(ITYP)
	   end do
	   write(19,*)NSAT
	   write(19,'(a,f15.2,a,3f11.5)')
     +' after ',FULTIM*1.d3,' fs,   BOX:',BOXL,BOYL,BOZL
	   do I=1,NSTOT
	     ITYP=ITYPE(I)
	     IS=NSITE(I)
	     if(LIST(ITYP).ne.0)write(19,'(a2,3(1x,f9.4))')
     +       NM(IS)(1:2),SX(I),SY(I),SZ(I)
	   end do
	end if
	IAN=IAN+1
	if(BOXL.le.0..or.BOYL.le.0..or.BOZL.le.0.)write(*,*)
     +'!!! problem with BOX size:',BOXL,BOYL,BOZL
	return
*  signal reading error
 40	write(*,*)' input error in file ',filnam(IFL),' T=',FULTIM
 45	write(*,'(2a/a,f10.3,a,i10)')
     +' finish analys of file ',filnam(IFL), 
     +  ' final T = ',fultim,'      N conf.  ',IAN
	close(12)
	IFL=IFL+1
	IIN=0
	if(IFL.gt.NF)then
	  IEND=1
	  return
	end if
	go to 35
	end
*
*================ READMOL =============================================
*
*    Define molecule from DATABASE
*
	subroutine READMOL(SUMM,NX,NTYP,PATHDB,NAMN)
	include "tranal.h"
*
        DIMENSION CM(3)
	character PATHDB*128,namn*128,FIL*256
        FIL=' '
	NX0=NX             
	NRBB = 0
	IE=0
	IL=0 
	do I=1,128
	  ISTR=ichar(PATHDB(I:I))
	  if(ISTR.le.15)then
	    LD=I-1
	    go to 11
	  end if
        if(PATHDB(I:I).ne.' ')IL=I
	end do 
	LD=IL
 11	continue
	FIL(1:LD)=PATHDB(1:LD)
        FIL(LD+1:LD+1)='/'
	IL=0 
	do I=1,32
	  ISTR=ichar(NAMN(I:I))
	  if(ISTR.le.15)then
	    LN=I-1
	    go to 11
	  end if
        if(NAMN(I:I).ne.' ')IL=I
	end do 
	LN=IL
 12	continue
	FIL(LD+2:LD+LN+1)=NAMN(1:LN)
	FIL(LD+LN+2:LD+LN+6)='.mmol'
	open(unit=12,file=fil,status='old',err=999)
	if(IPRINT.ge.5)then
	PRINT "(80('-'))"
        PRINT *
        write(*,'(a,i4,a)')'*** MOLECULE OF TYPE No.',NTYP,NAMN(1:LN)
	end if
	STR=TAKESTR(12,IE)
	read(STR,*)NSTS
	if(IPRINT.ge.6)PRINT * ,'*** NR OF SITES:                ',NSTS
	if(NX+NSTS.gt.NS)
     & stop '!!! Increase NS (number of sites) in tranal.h'
	NSITS(NTYP)=NSTS
	NSBEG=NX+1
	do I=1,NSTS
	  STR=TAKESTR(12,IE)
	  IX=NX+I
	  NM(IX)=STR(1:4)
	  read(STR(5:128),*,err=903)(R(IX,J),J=1,3),
     +  MASS(IX),CHARGE(IX)
          if(IPRINT.ge.6)
     +      write(*,'(a5,a4,a4,3f7.3,a3,f8.4,a3,f6.3,a5,f8.4,a4,f9.4)')
     +      'atom ',NM(IX),'  R=',(R(IX,J),J=1,3),
     +      ' M=',MASS(IX),' Q=',CHARGE(IX)
	end do
	NX=NX+NSTS
	NSEND=NX
*   Reading eventual bond list
        if(LBOND)then
          STR=TAKESTR(12,IE)
          read(STR,*,err=904,end=904)NSR
          do I=1,NSR
	    STR=TAKESTR(12,IE)
	    if(IPRINT.ge.6)write(*,'(a128)')STR
          end do
          STR=TAKESTR(12,IE)
C  Number of bonds
          read(STR,*,err=902,end=902)NBD
          NRB(NTYP)=NBD
          do K=NRBB+1,NRBB+NBD
	    STR=TAKESTR(12,IE)
	    read(STR,*,err=902,end=902)IDUM,IB(K),JB(K)
          end do
	  NRBB = NRBB+NBD
	  IADB(NTYP)=NRBB-NRB(NTYP)
        end if
* COM coordinates 
      SUMM       = 0.D0
      DO IS      = NSBEG,NSEND
        SUMM       = SUMM+MASS(IS)
      END DO! OF IS
      SUMMAS(NTYP)=SUMM
*
      CM(1)      = 0.D0
      CM(2)      = 0.D0
      CM(3)      = 0.D0
      DO IS      = NSBEG,NSEND
        CM(1)      = CM(1)+R(IS,1)*MASS(IS)
        CM(2)      = CM(2)+R(IS,2)*MASS(IS)
        CM(3)      = CM(3)+R(IS,3)*MASS(IS)
      END DO! OF IS
*
      CM(1)      = CM(1)/SUMM
      CM(2)      = CM(2)/SUMM
      CM(3)      = CM(3)/SUMM
	if(IPRINT.ge.10)then
        PRINT *,'*** Molecular GEOMETRY (C.O.M. IN ORIGIN): '
        DO IS      = NSBEG,NSEND
          PRINT "('*** ',I4,3x,A4,3X,3(1X,F7.3))",
     +    IS,NM(IS),R(IS,1),R(IS,2),R(IS,3)
        END DO! OF IS
      end if
      return
 903	write(*,*)' Error in the list of atoms, line ',I
	stop
 902  write(*,*)' Error in the list of bonds, line ',K
	stop
 904	write(*,*)' Error in the reference of .mol file, line '
	stop
 999  write(*,*)'!!! Molecule ',NAMN,' not found in the data base'
      write(*,*)'    File not found: ',FIL
	stop
	end 
*
*============== TAKESTR =========================================
*
      character*128 function TAKESTR(KAN,IE)
	character*128 AUX, fn*32 
	integer LCOUNT(256)   
	data LCOUNT /256*0/
	save AUX,LCOUNT           
	if(IE.eq.99)go to 10 
 1	LCOUNT(KAN)=LCOUNT(KAN)+1
 	read(KAN,'(a128)',err=10,end=20)AUX
	if(AUX(1:1).eq.'#')go to 1
	TAKESTR=AUX 
	IE=0
	return
 10	write(*,*)'!!! error in input file '   
	if(KAN.eq.5)then
	  write(*,*)' Standard input ' 
	else
	  inquire(unit=KAN,name=fn) 
	  write(*,*)' File ',fn
	end if
	write(*,*)' in line ',LCOUNT(KAN)
	write(*,*)AUX
	stop
 20	if(IE.ge.0)then
 	  write(*,*)'!!! end of file reached' 
	  if(KAN.eq.5)then
	    write(*,*)' Standard input ' 
	  else
	    inquire(unit=KAN,name=fn) 
	    write(*,*)' File ',fn
	  end if
	  stop
	end if  
	TAKESTR='  ' 
	IE=-1
	return
	end 
*
*===================== UTILITIES =================================
*
*
*=============== GETCOM ============================================
*                                      
*    Compute molecular centres of mass
*  --------------------------------------------------------------------
C     input: SX,SY,SZ coordinates
C     output:  X,Y,Z - molecular COMc
C              WX,WY,WZ - atom coordinates relative to molecular COMs
C              PX,PY,PZ - molecular momenta
C              QX,QY,QZ - molecular angular momenta
C
      SUBROUTINE GETCOM
      include "tranal.h"
*
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES                ! over types
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        SUMM        = SUMMAS (ITYP)
        DO J        = 1,NSPEC(ITYP)         ! over molecules
          I           = I+1
          X (I)       = 0.D0
          Y (I)       = 0.D0
          Z (I)       = 0.D0 
	  PX (I)       = 0.D0
          PY (I)       = 0.D0
          PZ (I)       = 0.D0
C    Calculate C.O.M. vectors relative to the first atom
          N0 = N+1
          DO IS       = NSBEG,NSEND          ! over sites
            N         = N+1
            DX =SX(N)-SX(N0) 
            DY =SY(N)-SY(N0) 
            DZ =SZ(N)-SZ(N0)
            call PBC(DX,DY,DZ)
            X (I)     = X(I)+MASS(IS)*DX
            Y (I)     = Y(I)+MASS(IS)*DY
            Z (I)     = Z(I)+MASS(IS)*DZ
            PX (I)    = PX(I)+VX(N)*MASS(IS)
            PY (I)    = PY(I)+VY(N)*MASS(IS)
            PZ (I)    = PZ(I)+VZ(N)*MASS(IS)
CD            if(NNUM(N).ne.I)
CD     +write(*,*)' wrong atom/mol in GETCOM -1: ',N,NNUM(N),I,IS,TASKID
          END DO! OF IS
          X (I)       = X(I)/SUMM+SX(N0)
          Y (I)       = Y(I)/SUMM+SY(N0)
          Z (I)       = Z(I)/SUMM+SZ(N0)
	  call PBC(X(I),Y(I),Z(I))
        END DO! OF I
      END DO! OF ITYP
*
*   7.2 Calculate local atom coordinates relative to molecular COM:
*   --------------------------------------------------------------
      N           = 0
      I           = 0
      DO ITYP     = 1,NTYPES
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        DO J        = 1,NSPEC(ITYP)
          I           = I+1
          QX(I)    = 0.D0
          QY(I)    = 0.D0
          QZ(I)    = 0.D0
          DO IS       = NSBEG,NSEND
            N           = N+1
Calculate short site vectors
            WX(N)       = SX(N)-X(I)
            WY(N)       = SY(N)-Y(I)
            WZ(N)       = SZ(N)-Z(I)
            call PBC(WX(N),WY(N),WZ(N))
            QX(I)    = QX(I)+(VY(N)*WZ(N)-WY(N)*VZ(N))*MASS(IS)
            QY(I)    = QY(I)+(VZ(N)*WX(N)-WZ(N)*VX(N))*MASS(IS)
            QZ(I)    = QZ(I)+(VX(N)*WY(N)-WX(N)*VY(N))*MASS(IS)
CD            if(NNUM(N).ne.I)
CD     +write(*,*)' wrong atom/mol in GETCOM -2: ',N,NNUM(N),I,IS,TASKID
          END DO! OF IS
        END DO! OF I
      END DO! OF ITYP
*  reset molecular to COM
*
      do I=1,NSTOT
	IMOL=NNUM(I)
	SX(I)=X(IMOL)+WX(I)
	SY(I)=Y(IMOL)+WY(I)
	SZ(I)=Z(IMOL)+WZ(I)
      end do
      RETURN
      END
*
*======================= PBC ========================================
*
      subroutine PBC(XX,YY,ZZ)
      include "tranal.h"
      if(XX.gt. HBOXL)XX=XX-BOXL
      if(XX.lt.-HBOXL)XX=XX+BOXL
      if(YY.gt. HBOYL)YY=YY-BOYL
      if(YY.lt.-HBOYL)YY=YY+BOYL
      if(ZZ.gt. HBOZL)ZZ=ZZ-BOZL
      if(ZZ.lt.-HBOZL)ZZ=ZZ+BOZL
      return
      end 
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
        if(STR(I).ne.' ')IL=I
	end do 
	LENG=IL
	return
      end   
*
*
*============= VECT =========================================
*	                                                       
	subroutine VECT(R1,R2,R3)
	real*8 R1(3),R2(3),R3(3)
	R3(1)= R1(2)*R2(3)-R1(3)*R2(2)
	R3(2)=-R1(1)*R2(3)+R1(3)*R2(1)
	R3(3)= R1(1)*R2(2)-R1(2)*R2(1) 
	return
	end
*
*================ NORM =======================================
*                                                             
	subroutine NORM(A,B,R,N)
	real*8 R,A(*),B(*) 
	R=0.d0
	do I=1,N
        R=R+A(I)**2
	end do       
	R=sqrt(R) 
	if(R.eq.0.d0)then
	  do I=1,N
	    B(I)=0.d0
	  end do
	else
	  do I=1,N
	    B(I)=A(I)/R
	  end do 
	end if
	return
	end
*....
*======== XINTGR =======================================================
*....
      FUNCTION XINTGR(H,Y,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(0:NDIM)
      HT=H/3.0
      IF(NDIM.LT.5)RETURN
      SUM1=4.00*Y(1)
      SUM1=HT*(Y(0)+SUM1+Y(2))
      AUX1=4.00*Y(3)
      AUX1=SUM1+HT*(Y(2)+AUX1+Y(4))
      AUX2=HT*(Y(0)+3.8750*(Y(1)+Y(4))+2.6250*(Y(2)+Y(3))+Y(5))
      SUM2=4.00*Y(4)
      SUM2=AUX2-HT*(Y(3)+SUM2+Y(5))
      AUX=4.00*Y(2)
      IF(NDIM-6)5,5,2
    2 DO 4 I=6,NDIM,2
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
*===================== INV3B3 ==========================================
*
*  Invert matrix 3x3
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
