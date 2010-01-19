*     PART 7
*
*     File restart.f
*     -------------
*
C     This file contains subroutines responsible for writing and
C     reading restart file, dumping trajectories and other operations
C     with input/output files
C
C     Subroutines:
C
C   1. RESTRT     Read or write restart file 
C   2. TRACE      Dumping trajectory
C   3. BIODUMP    Dump configuration for "BIOSIM" program
C   4. XMOLDUMP   Dump configuration in "XMOL" format
C   5. CONFDUMP   Dump coordinates and velocities for all atoms
C
*
*     1. Read or write restart file   
*     -----------------------------
C=============== RESTRT ==========================================
C
      SUBROUTINE RESTRT(NNN)
      include "prcm.h"
      character*12 LAB
C  reserved place in the .dmp file   
      real*8 SHUTA(10)
      save NSTTOT0,NSTEPA0
      data NSTTOT0/0/,NSTEPA0/0/,SHUTA/10*0.d0/
*
      GO TO (1000,2000),NNN
      STOP ' RESTRT: NO VALID OPTION !!!'
*
*    1.1 Read restart file
*    ---------------------
 1000 CONTINUE
      I43=1
*    1.1.1 Read formatted restart file
C   (this option appeared in v. 4.3)
      if(LFORMI)then
        open(unit=7,file=fdump,status="old",form='formatted',err=1001)
        read(7,'(a)')LAB
        if(LAB.ne.LABEL)then
          if(I43.ne.0.and.TASKID.eq.MAST)
     +    write(*,*)"!!! Wrong label of dump file: label=",LAB
        end if 
        MOP=NOP
        MS=NTYPES
        read(7,'(4i9,1x,d16.9)') NOP,NTYPES,NSTOT,NNAV,DTT
        READ(7,'(5i9,2d16.9,2i9)',err=19) 
     +  MSTEP,NSTTOT,NSTEPA,NAVT,NHIST,TIM,TIMA,NTRKF,NTRKZ
        if(IPRINT.ge.7)write(*,*)' read integer data '
C  checking some parameters 
        IF(MOP.NE.NOP)then
           write(*,*)
     +     '!!! NOP ! Number of particles wrong in RESTART FILE !!!'
        end if
        IF(MS.NE.NTYPES)then
           write(*,*)
     + '!!! NTYPES ! Number of molecule types wrong in RESTART FILE !!!'
           write(*,*)' use now as in the input file'
           NTYPES=MS
        end if
        READ (7,'(5d16.9)',err=19) 
     +  SC,SCL,BOXL,VOL,CUTRDF,BOYL,BOZL,SCLX,SCLY,SCLZ,
     +  (SCM(I),I=1,NTYPES) 
        if(IPRINT.ge.7)write(*,*)' read BOX sizes '
        read(7,'(5d16.9)',err=19) 
     +timest,timea,timeb,timee,timef,timeg,timel,timen,timet,timev,times
        if(IPRINT.ge.7)write(*,*)' read times'
        read(7,'(5d16.9)',err=19) (SX(I),SY(I),SZ(I),I=1,NSTOT)
        if(IPRINT.ge.7)write(*,*)' read coordinates'
        read(7,'(5d16.9)',err=19) (VX(I),VY(I),VZ(I),I=1,NSTOT)
        if(IPRINT.ge.7)write(*,*)' read velocities'
        READ (7,'(5d16.9)',err=28) (AV(I),AW(I),I=1,NNAV)
        if(IPRINT.ge.7)write(*,*)' read averages'
        do J=1,NHIST
          read(7,'(5d16.9)',err=29)(HIST(I,J),I=1,NNAV)
          if(IPRINT.ge.7)write(*,*)' read history ',J
        end do
        go to 30
 28     write(*,*)' error reading averages'
        go to 30
 29     write(*,*)' error reading history'
 30     close(7)
*    1.1.2 Read unformatted restart file 
C         (default)
      else
        open(unit=7,file=fdump,status="old",form='unformatted',err=1001)
C  First line - program label
    5   read(7)LAB
        if(LAB.ne.LABEL)then
          if(I43.ne.0.and.TASKID.eq.MAST)
     +    write(*,*)"!!! Wrong dump file: label=",LAB
        end if 
        MOP=NOP
        MS=NTYPES
C   Second line - number of particles, sites, atoms, averages and time step
        read(7) NOP,NTYPES,NSTOT,NNAV,DTT
C   Third line:   MSTEP - MD steps made in previous run (not really used)
C                 NSTTOT - steps made from the beginning of the simulation
C                 NSTEPA - num. of steps with time accumulation
C                 NAVT   - total averaged configurations
C                 NHIST  - number of averaging "window"
C                 TIM    - real simulation time (in s)
C                 NTRKF  - number of the last trajectory file
C                 NTRKZ  - number of config. in the last trajectory file  
        READ (7) MSTEP,NSTTOT,NSTEPA,NAVT,NHIST,TIM,TIMA,NTRKF,NTRKZ
C  checking some parameters 
        IF(MOP.NE.NOP)then
           write(*,*)
     +     '!!! NOP ! Number of particles wrong in RESTART FILE !!!'
        end if
        IF(MS.NE.NTYPES)then
           write(*,*)
     + '!!! NTYPES ! Number of molecule types wrong in RESTART FILE !!!'
           write(*,*)' use now as in the input file'
           NTYPES=MS
        end if
C  for compatibility with v<4.3
        if(I43.ne.0)then 
          READ (7,err=8) SC,SCL,BOXL,VOL,CUTRDF,BOYL,BOZL,SCLX,SCLY,
     +    SCLZ,(SCM(I),I=1,NTYPES) 
        else
          READ (7,err=14) SC,SCL,BOXL,VOL,CUTRDF,BOYL,BOZL,SCLX,SCLY,
     +    SCLZ 
C  for compatibility with v.<4.2
C         READ (7,err=14) SC,SCL,BOXL,VOL,CUTRDF,BOYL,BOZL 
          do I=1,NTYPES
            SCM(I)=0.
          end do
	  go to 13
 14 	  LSEP=.false.
        end if   ! if(I43
        go to 13
  8     I43=0 
        rewind (7)
	go to 5
C    forth line - CPU time wasted
  13    read(7,err=12) timest,timea,timeb,timee,timef,timeg,timel,
     +  timen,timet,timev,times
C    5-th line - coordinates
  12    read(7) (SX(I),SY(I),SZ(I),I=1,NSTOT)
C    6-th line - momenta (velocities)
        read(7) (VX(I),VY(I),VZ(I),I=1,NSTOT)
C    7-th line - accumulated averages in the last window
        READ (7) (AV(I),AW(I),I=1,NNAV)
        do J=1,NHIST
C    8-th - 7+NHIST lines : averages from all previous windows
          read(7,err=11)(HIST(I,J),I=1,NNAV)
        end do
 11     close (7)    
      end if
      if(.not.LSEP)then                             
   	SCLX=SCL/3.
   	SCLY=SCL/3.
	SCLZ=SCL/3.    
      end if
      if(LGRST.and.CUTRDF.gt.0.d0)RDFCUT=CUTRDF
      NSTTOT0=NSTTOT
      NSTEPA0=NSTEPA
*  1.1.3   Report status of previous run 
      if(TASKID.eq.MAST)then
        PRINT *
        PRINT "('*** restart file read in from ',a)",fdump
        PRINT *
        PRINT "('*** NR OF PARTICLES            ',I6)",MOP
        PRINT "('*** NR OF SITES                ',I6)",MS
        PRINT "('*** NR OF PREVIOUS STEPS       ',I8)",NSTTOT
        print "('*** NR of steps with time ac.  ',I8)",NSTEPA
        PRINT "('*** TIME STEP USED(in fs)              ',f10.4)",
     +  DTT*1.d15          
      end if
*  1.1.4 Scale reduce velocities if time step changed
      IF(DTT.NE.DT)then
        if(TASKID.eq.MAST)PRINT 
     +  "('... NOT THE SAME AS IN PRESENT RUN --> ',f10.4)",DT*1.d15
        if(TASKID.eq.MAST)PRINT *, "Velocities are scaled"
	  do I=1,NSTOT
	    VX(I)=VX(I)*DT/DTT
	    VY(I)=VY(I)*DT/DTT
	    VZ(I)=VZ(I)*DT/DTT
	  end do  
	  SCL=SCL*DT/DTT     
	  SCY=SCY*DT/DTT     
	  SCZ=SCZ*DT/DTT     
	  SC=SC*DT/DTT     
	  call GETEMP(.true.)
	  write(*,*)'Temperature is ',TEMP
	end if
	if(TASKID.eq.MAST)
     +PRINT "('*** short time step (in fs)            ',f10.4)",
     +DT*1.d15/NFREQ
*   1.1.5 Zero cpu time counter (if specified)
	if(L0CPU)then
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
	  NSTEPA0=0
	end if
*   1.1.6 Zero velocities (if specified)
      IF(L0VS) THEN    
        DO N  = 1,NSTOT
          VX(N) = 0.D0
          VY(N) = 0.D0
          VZ(N) = 0.D0
        END DO! OF N
      end if
      RETURN
*
*   1.2 Dump restart file
*   ---------------------
 2000 CONTINUE
*   1.2.1 Prepare
      if(TASKID.eq.MAST)then
        NSTTOT1=NSTTOT0+NSTEP
        NSTEPA1=NSTEPA0+NSTEP
        TIM0	= TIM+NSTEP*DT
        TIMA0	= TIMA+NSTEP*DT
        time	= timest+cputime(time0)
	if(.not.LSEP)then                             
   	  SCLX=SCL/3.
   	  SCLY=SCL/3.
	  SCLZ=SCL/3.    
	end if
        write(*,*)
*   1.2.2 Dump in formatted form
        if(LFORMO)then
          open (unit=8,file=fdump,status="unknown",form='formatted')
          write(8,'(a)')LABEL
          write(8,'(4i9,1x,d16.9)') NOP,NTYPES,NSTOT,NNAV,DT
          write(8,'(5i9,2d16.9,2i9)') 
     +NSTEP,NSTTOT1,NSTEPA1,NAVT,NHIST,TIM0,TIMA0,NTRKF,NTRKZ
          WRITE(8,'(5d16.9)') SC,SCL,BOXL,VOL,RDFCUT,BOYL,BOZL,
     +SCLX,SCLY,SCLZ,(SCM(I),I=1,NTYPES)
          write(8,'(5d16.9)') 
     +time,timea,timeb,timee,timef,timeg,timel,timen,timet,timev,times      
	  WRITE(8,'(5d16.9)') (SX(I),SY(I),SZ(I),I=1,NSTOT)
          WRITE(8,'(5d16.9)') (VX(I),VY(I),VZ(I),I=1,NSTOT)
          write(8,'(5d16.9)') (AV(I),AW(I),I=1,NNAV)
          do J=1,NHIST
            write(8,'(5d16.9)')(HIST(I,J),I=1,NNAV)
          end do
          write(6,'(2a)')'*** formatted restart file dumped on ',fdump
        else
*   1.2.3  Dump in unformatted form (default)  
C          Comments on format of the restart file see in p.1.1.2
          open (unit=8,file=fdump,status="unknown",form='unformatted')
          write(8)LABEL
          write(8) NOP,NTYPES,NSTOT,NNAV,DT
          write(8) NSTEP,NSTTOT1,NSTEPA1,NAVT,NHIST,TIM0,TIMA0
     +    ,NTRKF,NTRKZ
          WRITE(8) SC,SCL,BOXL,VOL,RDFCUT,BOYL,BOZL,SCLX,SCLY,SCLZ,
     +    (SCM(I),I=1,NTYPES),SHUTA
          write(8) time,timea,timeb,timee,timef,timeg,timel,timen,timet,
     +    timev,times      
	  WRITE(8) (SX(I),SY(I),SZ(I),I=1,NSTOT)
          WRITE(8) (VX(I),VY(I),VZ(I),I=1,NSTOT)
          write(8) (AV(I),AW(I),I=1,NNAV)
          do J=1,NHIST
            write(8)(HIST(I,J),I=1,NNAV)
          end do
          write(6,'(2a)')'*** restart file dumped on ',fdump
        end if
        close(8)
      end if
*   1.2.4  Dump RDF
      if(LGR.and.LGDMP.and.NHIST.ge.IHIST)call RDFOUT(0)      ! aver.f
      if(TASKID.ne.MAST)return 
*   1.2.5  Dump TCF
      if(LCF.and.LCFDMP.and.NHIST.ge.IHIST)call TCFOUT(0) ! tcf.f
      RETURN       
*   1.2.6  I/O error diagnostics
 1001 write(*,*)' File ',fdump,' not found on node ',TASKID
      call FINAL
 19   write(*,*)
     +' Input error trying to read formatted restart file ',fdump
      call FINAL
      END 
*
*======================= TRACE ========================================
*
*   2. Dump trajectory
*   ------------------
*

      subroutine TRACE
C  Front matter
      include "prcm.h"
      logical lfirst 
      real*8 fdir(3)
      data lfirst /.true./,fdir/3*0.d0/
*
*   2.1 Check number of configurations
*   ----------------------------------
C   Dump trajectory only from master node, and if parameter ITREK > 0
      if(ITREK.le.0.or.TASKID.ne.MAST)return
C   NTRKZ is number of current configuration, written in the trajectory file 
      NTRKZ=NTRKZ+1
C   If number of configurations exceed given level, close trajectory file
C   and prepare to start new (set "lfirst" to .true.)
      if(NTRKZ.gt.LTREK)then
	if(.not.lfirst)then
	  close(55)
	  write(*,*)' close file ',FTRK
	end if
	lfirst=.true.
      end if
C   IVEL = 0/1  dump or no dump velocities
      if(ITREK.le.2)then
         IVEL=0
      else
         IVEL=1
      end if
*
*    2.2 Start a new trajectory file 
*    -------------------------------
*    2.2.1. Give name to a new file
      if( lfirst ) then ! initialization of trk file
        lfirst = .false.
        NTRKZ=1
        if(NTRKF.eq.0)NTRKZ=0
C    NTRKF is number of current trajectory file
	NTRKF=NTRKF+1
        if(NTRKF.ge.1000)NTRKF=NTRKF-1000
C    Create name of new trajectory file  (ftrk)
        write(ftrk(ICFN:ICFN+2),'(i3)')NTRKF
	if(ftrk(ICFN:ICFN).eq.' ')FTRK(ICFN:ICFN)='0'
	if(ftrk(ICFN+1:ICFN+1).eq.' ')FTRK(ICFN+1:ICFN+1)='0'
*    2.2.2  Open new trajectory file
	if(mod(ITREK,2).ne.0)then
C - start to write binary file:
          open(unit=55, form='unformatted',file=ftrk)
	else
C - start to write ASCII file (XMOL format)
	  open(unit=55, file=ftrk)
	end if
*    2.2.3  Write heading for the new trajectory file
C   This is only for binary file   
C   First line:
c
c - 0/1 - no velocites / dump velocities in Å/fs     I4
c - deltat - time step            [s]                R8
c - box size - [Angstr]                            3*R8
c - unitl - [m]                                      R8
c - ntypes - number of mol. types                    I4
c - nspecs, nsites, .true.      ntypes*(2*I4+L4)
c
	if(mod(ITREK,2).ne.0)then
          write( 55 ) IVEL, dt, boxl, boyl, bozl, unitl, ntypes,
     x                (nspec(i),nsits(i),.true.,i=1,ntypes)
c
c - for each type: masses of atoms and LIST = 0 or 1
C   LIST = 0 - no coordinates for this type
C          1 - coordinates are present  
c - 
          do ityp=1, ntypes
            NSB            =  ISADR (ITYP)+1
            NSE            =  ISADR (ITYP+1)
            write( 55 ) (mass(j), j=NSB,NSE),LIST(ITYP)
          end do
	end if
	if(ITREK.ne.0)write(*,*)' start of dumping trace in file ',ftrk
      endif ! first call
*
*    2.3  Dumping a configuration:
*
	fultim=TIM+NSTEP*DT
*    2.3.1  ASCII format (XMOL)
	if(mod(ITREK,2).eq.0)then
	  NSDP=0 
	  do ITYP=1,NTYPES
	    if(LIST(ITYP).ne.0)NSDP=NSDP+NSPEC(ITYP)*NSITS(ITYP)
	  end do
	  write(55,*)NSDP
	  if(IPRINT.ge.5)
     +    write(*,*)'add ',NTRKZ,' config. after ',fultim*1.d15,' fs'
C  Second commentary line of XMOL format contains time in fs
	  write(55,'(a,f15.2,a,3f12.4)')
     +' after ',fultim*1.d15,' fs, BOX: ',BOXL,BOYL,BOZL
	  do I=1,NSTOT
	    ISB=NSITE(I)
	    ITYP=ITYPE(I)
	    if(LIST(ITYP).ne.0)then
               if(IVEL.le.0)then
                  write(55,'(a2,3(1x,f8.3))')
     +            NM(ISB)(1:2),SX(I),SY(I),SZ(I)
               else
                  VVX=VX(I)*MASSDI(I)*1.d-15/UNITT
                  VVY=VY(I)*MASSDI(I)*1.d-15/UNITT
                  VVZ=VZ(I)*MASSDI(I)*1.d-15/UNITT
                  write(55,'(a2,3(1x,f8.3),3(1x,f11.5))')
     +            NM(ISB)(1:2),SX(I),SY(I),SZ(I),VVX,VVY,VVZ
               end if
            end if
            end do
*    2.3.2  Binary format
	else                    
  	  cumtim=HIST(1,IHIST)*1.d-12
	  FDIR(1)=TEMP
          FDIR(2)=TRYCK
          FDIR(3)=POTEN
C  First line of a record:
C   null (not used)
c - fultim - full time from the beginning                [s]  R8
c - cumtim - full time from begining collecting averages  [s] R8
c - fdir - temperature,pressure,pot.energy                  3*R8
C   BOXL,BOYL,BOZL                                          3*R8
C   LIST                                                      I4
          write( 55 ) IVEL, fultim, cumtim, fdir,BOXL,BOYL,BOZL,LIST
	  if(IPRINT.ge.5)
     +    write(*,*)'add ',NTRKZ,' config. after ',fultim*1.d15,' fs'
c
c - one record for each type
        do ityp=1,ntypes
           NSB          = ISADDR(ITYP)+1
           NSE         = ISADDR(ITYP+1)
           if(LIST(ITYP).ne.0)then
              if(IVEL.le.0)then
                 write( 55 ) (real(SX(j)),j=NSB,NSE),
     *                (real(SY(j)),j=NSB,NSE),
     *                (real(SZ(j)),j=NSB,NSE)
              else
                 write( 55 ) (real(SX(j)),j=NSB,NSE),
     *                (real(SY(j)),j=NSB,NSE),
     *                (real(SZ(j)),j=NSB,NSE)
         write( 55 ) (real(VX(j)*MASSDI(j)*1.d-15/UNITT),j=NSB,NSE),
     *               (real(VY(j)*MASSDI(j)*1.d-15/UNITT),j=NSB,NSE),
     *               (real(VZ(j)*MASSDI(j)*1.d-15/UNITT),j=NSB,NSE)
              end if
           end if
        end do
	end if
      end 
*
*==================== BIODUMP =======================================
*
*   3. Dump configuration for the BIOSYM or CERIUS program
*   ------------------------------------------------------
      subroutine BIODUMP
* System configuration for BIOSYM program
      include "prcm.h" 
      ANG=90.
      if(TASKID.ne.MAST)return
      open(unit=19,file=fmol,status='unknown')
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
      write(*,*)' Configuration file dumped to ',fmol
      return
      end               
*
*===================== XMOLDUMP ======================================
*
*   4. Dump configuration in "XMOL" format
*   --------------------------------------
      subroutine XMOLDUMP
      include "prcm.h"     
      if(TASKID.ne.MAST)return
      open(unit=21,file=fmol,status='unknown') 
*
*  4.1 Dump only COM coordinates of molecules
* 
      if(LCOMD)then
        do I=1,NOP
          ITYP=ITM(I)
          write(21,'(3f10.4,3x,a6,2x,i6)')X(I),Y(I),Z(I),NAME(ITYP),I
        end do
*
*  4.2 Dump coordinates in XMOL format
*
      else
        NSDP=0
	do ITYP=1,NTYPES
	    if(LIST(ITYP).ne.0)NSDP=NSDP+NSPEC(ITYP)*NSITS(ITYP)
	end do
	write(21,*)NSDP 
	fultim=TIM+NSTEP*DT
	write(21,'(a,f15.2,a,3f12.4)')
     +' after ',fultim*1.d15,' fs, BOX: ',BOXL,BOYL,BOZL
	do I=1,NSTOT 
          ITYP=ITYPE(I)
          if(LIST(ITYP).ne.0)then
	    IS=NSITE(I)
            if(IVEL.le.0)then
              write(21,'(a2,3(1x,f8.3))')
     +        NM(IS)(1:2),SX(I),SY(I),SZ(I)
            else
              VVX=VX(I)*MASSDI(I)*1.d-15/UNITT
              VVY=VY(I)*MASSDI(I)*1.d-15/UNITT
              VVZ=VZ(I)*MASSDI(I)*1.d-15/UNITT
              write(21,'(a2,3(1x,f8.3),3(1x,f11.5))')
     +        NM(IS)(1:2),SX(I),SY(I),SZ(I),VVX,VVY,VVZ
            end if
          end if
	end do
      end if
      close(21)
      return
      end
*
*====================== CONFDUMP ==================================
*
*   5. Dump coordinates and velocities for all atoms
*   ------------------------------------------------
*
C   This usually called if something goes wrong
C
      subroutine CONFDUMP(IOU)
      include "prcm.h"
C  write down coordinates and velocities of all atoms
C  (used only for debug purposes)
      write(IOU,*)' Coordinates and velocities of all particles:'
      write(IOU,*)' node ',TASKID
      do I=1,NSTOT
        VFC=MASSDI(I)*1.d-12/UNITT
        AT=1.5*(VX(I)**2+VY(I)**2+VZ(I)**2)*MASSDI(I)*ENERF*TEMPF
	write(IOU,'(I5,2x,a4,3(1x,f9.3),2x,3(1x,f8.4),1x,F8.0)')
     +I,NM(NSITE(I)),SX(I),SY(I),SZ(I),VX(I)*VFC,VY(I)*VFC,VZ(I)*VFC,AT
      end do
      return
      end
