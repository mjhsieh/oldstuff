*=============== RESTRT ==========================================
*
      SUBROUTINE RESTRT(NNN)
*
      include "mdee.h"
      character*12 LAB
	save NSTTOT0,NSTEPA0
	data NSTTOT0/0/,NSTEPA0/0/
*
      GO TO (1000,2000),NNN
*
      STOP ' NO VALID OPTION IN RESTRT !!!'
 1000 CONTINUE
      open (unit=7,file=fdump,status="old",form='unformatted')
      read(7)LAB
      if(LAB.ne.LABEL)then
        write(*,*)"!!! Wrong dump file: label=",LAB
*        stop ' '
      end if 
      MOP=NOP
      MS=NSITES
      read(7) NOP,NSITES,NSTOT,NNAV,DTT
      READ (7) MSTEP,NSTTOT,NSTEPA,NAVT,NHIST,TIM,TIMA,NNH,
     +ME,NEOLD
      IF(MOP.NE.NOP)then
        write(*,*)'NOP =',NOP,'...    CERTAINLY WRONG RESTART FILE !!!'
        stop
      end if
      IF(MS .NE.NSITES )then
      write(*,*)'Nsite =',NSITES,'...  CERTAINLY WRONG RESTART FILE !!!'
      stop
      end if
      READ (7) SC,SCL,BOXL,VOL 
      if(LECH)then
	  read(7)
          if(NE.ne.NEOLD)then
             write(*,*)' Number of subensembles has changed to ',NE
             ME=ME0
             write(*,*)' Now in ',ME0,' subensemble'
          end if
	  write(*,*)' New EE parameters. Iteration ',NHIST 
	  if(IPRINT.ge.6)write(*,'(10(10f8.3))')(EC(IC),IC=1,NE)
	  if(IPRINT.ge.6)write(*,'(10(10f8.3))')(EE(IC),IC=1,NE)
      else 
         NE=NEOLD
	 read(7)(EE(I),EC(I),I=1,NE)
      end if
      read(7) timest,timea,timeb,timee,timef,timeg,timel,
     +timen,timet,timev,times,time2
  12  read(7) (SX(I),SY(I),SZ(I),I=1,NSTOT)
      read(7) (VX(I),VY(I),VZ(I),I=1,NSTOT)
	do J=1,NEOLD
      READ (7) (AV(I,J),AW(I,J),I=1,NNAV)
	end do
      do J=1,NHIST
        read(7)(HIST(I,J),I=1,NNH),(ITE(I,J),I=1,NEOLD)
      end do
	if(LECH)HIST(3,NHIST)=1.d0
        read(7,end=11,err=11)((IWALK(I,K),I=1,NEOLD),K=1,NEOLD)
	read(7,end=11,err=11)ICHU,ICHE
 11   close (7)
	NSTTOT0=NSTTOT
	NSTEPA0=NSTEPA
      PRINT *
      PRINT "('*** restart file read in from ',a)",fdump
      PRINT *
      PRINT "('*** NR OF PARTICLES              ',I6)",MOP
      PRINT "('*** NR OF SITES                  ',I6)",MS
      PRINT "('*** NR OF PREVIOUS STEPS       ',I8)",NSTTOT
      print "('*** NR of steps with averaging ',I8)",NSTEPA
      PRINT "('*** TIME STEP USED(in fs)              ',f10.4)",
     +DTT*1.d15
	IF(DTT.NE.DT)then
        PRINT "('... NOT THE SAME AS IN PRESENT RUN --> ',f10.4)",
     +  DT*1.d15
        PRINT *, "Velocities are scaled"
	  do I=1,NSTOT
	    VX(I)=VX(I)*DT/DTT
	    VY(I)=VY(I)*DT/DTT
	    VZ(I)=VZ(I)*DT/DTT
	  end do
	end if
	PRINT "('*** short time step (in fs)            ',f10.4)",
     +DT*1.d15/NFREQ
*
	write(*,*)' Num of subensembles ',NE,' now in ',ME
	if(L0AVS)then
	  TIMEST=0.
	  TIMEA=0.
	  TIMEB=0.
	  TIMEE=0.
	  TIMEF=0.
	  TIMEG=0.
	  TIMEL=0.
	  TIMEN=0.
	  TIMET=0.
	  TIMEV=0.
	  TIMES=0.
	  TIME2=0.
	end if
      IF(L0VS) THEN    
      DO N  = 1,NSTOT
      VX(N) = 0.D0
      VY(N) = 0.D0
      VZ(N) = 0.D0
      END DO! OF N 
*
      PRINT "(/,'*** OLD VELOCITIES AND FORCES ZEROED ! ***',/)"
*
      END IF
*
      DO N  = 1,NSTOT
      OX(N) = SX(N)
      OY(N) = SY(N)
      OZ(N) = SZ(N)
      END DO! OF N
*
      RETURN
 2000 CONTINUE
      NSTTOT1=NSTTOT0+NSTEP
      TIM0	= TIM+NSTEP*DT
      TIMA0	= TIMA+NSTEP*DT
      time	= timest+cputime(time0)
      open (unit=8,file=fdump,status="unknown",form='unformatted')
      write(8)LABEL
      write(8) NOP,NSITES,NSTOT,NNAV,DT
      write(8) NSTEP,NSTTOT1,NSTEPA,NAVT,NHIST,TIM0,TIMA0
     +,NNH,ME,NE
      WRITE(8) SC,SCL,BOXL,VOL
	write(8)(EE(I),EC(I),I=1,NE)
      write(8) time,timea,timeb,timee,timef,timeg,timel,timen,timet,
     +timev,times,time2      
	WRITE(8) (SX(I),SY(I),SZ(I),I=1,NSTOT)
      WRITE(8) (VX(I),VY(I),VZ(I),I=1,NSTOT)
	do J=1,NE
        write(8) (AV(I,J),AW(I,J),I=1,NNAV)
	end do
      do J=1,NHIST
        write(8)(HIST(I,J),I=1,NNH),(ITE(I,J),I=1,NE)
      end do
      write(8)((IWALK(I,K),I=1,NE),K=1,NE)
      write(8)ICHU,ICHE,A1,A2,A3,A4,A5
      close(8)
      PRINT *
      write(6,'(2a)')'*** restart file dumped on ',fdump
      PRINT *
      if(LGR.and.LGDMP)call RDFOUT(0)
*
      RETURN
      END 
*
*======================= TRACE ========================================
*
      subroutine TRACE
c Dump to file all interesting values
      include "mdee.h"
c
	logical lfirst /.true./
c
      integer chan  /55/
c
      double precision fdir(3)/3*0.0d0/
	if(ITREK.le.0)return
	NTRKZ=NTRKZ+1
	if(NTRKZ.gt.LTREK)then
	  if(.not.lfirst)close(chan)
	  lfirst=.true.
	end if
c
      if( lfirst ) then ! initialization of trk file
        lfirst = .false.
c
c - delete old file (useful on Convex)
        NTRKZ=0
	  NTRKF=NTRKF+1
        write(ftrk(ICFN:ICFN+2),'(i3)')NTRKF
	  if(ftrk(ICFN:ICFN).eq.' ')FTRK(ICFN:ICFN)='0'
	  if(ftrk(ICFN+1:ICFN+1).eq.' ')FTRK(ICFN+1:ICFN+1)='0'
	  open( chan, file=ftrk, iostat=i ) ! open file as text-file
        if( i.ne.0 ) then
          write( *,* )'---Can not open coordinate file'
          return
        endif
        rewind chan      !\ - goto beginning
        write( chan, * ) !| - write 1 byte <LF>
        close( chan )    !/ - close == leave file contained 1 byte
c
c - start to write binary file:
        open( chan, form='unformatted', file=ftrk )
c
c - 0 (null) = nomer shaga (not in use)       I4
c - deltat - shag integrirovaniya [s]         R8
c - box* - [Angstr]                         3*R8
c - unitl - [m]                               R8
c - ntypes - kol-vo tipov                     I4
c - nspecs, nsites - ponyatno               2*I4
c - .true. - ostatok ot LMOVE(I) (not in use) L4
c
        write( chan ) 0, dt, boxl, boyl, bozl, unitl, ntypes,
     x                (nspec(i),nsits(i),.true.,i=1,ntypes)
c
c - dmass2(INDEX(ityp,iatom)) - atomic mass in [AMU]  R8
c - dlya kazhdogo tipa - otdel'naya FORTRAN'naya zapis
        do ityp=1, ntypes
          NSB            =  ISADR (ITYP)+1
          NSE            =  ISADR (ITYP+1)
          write( chan ) (mass(j), j=NSB,NSE)
        end do
	  if(ITREK.ne.0)write(*,*)' start of dumping trace in file',ftrk
      endif ! first call
c
c ---- novyj shag:
c
c - 0 (null) - nomer shaga
c - fultim - polnoe vremya s nachala MD-modelirovaniya [s]  R8
c - cumtim - vremya s nachala analiza MD-modelirovaniya [s] R8
c - fdir(3) - massiv orientacii direktora zh.k.           3*R8
	fultim=TIM+NSTEP*DT
	cumtim=HIST(1,IHIST)*1.d-12
      write( chan ) 0, fultim, cumtim, fdir,BOXL,BOYL,BOZL,ME
c
c - dlya kazhdogo tipa!
      do ityp=1,ntypes
           NSB          = ISADDR(ITYP)+1
           NSE         = ISADDR(ITYP+1)
        write( chan ) (real(SX(j)),j=NSB,NSE),
     *                (real(SY(j)),j=NSB,NSE),
     *                (real(SZ(j)),j=NSB,NSE)
      end do
      end 
*
*==================== BIODUMP =======================================
* 
	subroutine BIODUMP
* System configuration for BIOSYM program
	include "mdee.h" 
	ANG=90.
	open(unit=19,file=fbio,status='unknown')
	write(19,'(a)')'!BIOSYM archive 2'
      write(19,'(a)')'PBC=ON'
      write(19,'(a)')' system configuration'
      write(19,'(a)')'!DATE'
      write(19,'(a3,6f10.4,a5)')
     +'PBC',BOXL,BOYL,BOZL,ANG,ANG,ANG,' (P1)'
	N=0
	I=0
	do ITYP=1,NTYPES
        NSBEG       = ISADR  (ITYP)+1
        NSEND       = ISADR  (ITYP +1)
        DO J        = 1,NSPEC(ITYP)
	    I=I+1
          DO IS       = NSBEG,NSEND
            N           = N+1
            write(19,'(a4,1x,3f15.7,1x,a4,i4,2x,a2,1x,a2,1x,f8.4,i4)')
     + NM(IS),SX(N),SY(N),SZ(N),NAME(ITYP)(1:4),ITYP,'h ',NM(IS)(1:2),
     + CHARGE(IS),N
	    end do
	    write(19,'(a3)')'end'
	  end do
	end do
      write(19,'(a3)')'end'
	write(*,*)' Configuration file dumped to ',fbio
	return
	end  
