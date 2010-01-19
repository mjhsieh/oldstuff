*     PART 6  
*
*     File mpi.f
*     -------------
*
*     Subroutines for parallel execution (MPI)
*     Use single precision (R4) buffer to transfer the data
*     which saves time for communications but may loose
*     some precision
*
*    1. Parallel initialization
*    --------------------------          
C    This subroutine initialize MPI interface and get two numbers:
C      NUMTASK  - total number of nodes available
C      TASKID   - current node  (TASKID=0,1,...,NUMTASK-1)  
C
C======= PARAINI ===================================================
C
      subroutine PARAINI(numtask,taskid)   
      include "mpif.h"
      integer numtask,taskid
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numtask, ierr ) 
      print *, "Process ", taskid, " of ", numtask, " is alive"
      return
      end 
* 
*    2. Broadcast atom positions to all nodes
*    ----------------------------------------
C    This subroutine broadcast positions of atoms to all the nodes
C    Each node operate with atoms from NAB(TASKID) to NAE(TASKID),
C    but for force evaluation each node should know exact positions 
C    of all the atoms.
C    The are two regimes of this subroutine defined by parameter IND:
C      IND=1  - atoms distributed non-equally between the nodes
C               (case of constrained dunamics)
C      IND=2  - each node (except the last one) has equal number of atoms
C               This is used in DOUBLE subroutine 
C
C============ SINHR ===================================
C 
	subroutine SINHR(IND)
	include "prcm.h"
C  Making buffers real*4 will not affect precision noticaebly, but
C  save communication time greatly
	real*4 BUFF(NBUFF),BUF2(NBUFF)
	include "mpif.h"      
        timen0	= cputime(0.d0)
C   This normally not necessary
*      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
*
*   2.1 Case of non-equal number of atoms on processors
*   ---------------------------------------------------
	if(IND.eq.1)then   !    (used in constrain dynamics)
*   2.1.1 Put coordinates into buffer
C   This is number atoms on the processor    
          NPC=NAP(TASKID)          
	  NPC3=NPC*3
	  do K=1,NPC    
	    I=NABS(TASKID)+K
	    BUFF(K)=SX(I)
	    BUFF(K+NPC)=SY(I)
	    BUFF(K+2*NPC)=SZ(I)
	  end do
*    2.1.2 Communicate
C    This subroutine gather coordinates from all the nodes 
C    and then broadcast them
          call MPI_ALLGATHERV(BUFF,NPC3,MPI_REAL,BUF2,NAP3,
     +NABS3,MPI_REAL,MPI_COMM_WORLD,ierr)
C   This normally not necessary
          call MPI_BARRIER(MPI_COMM_WORLD,IERR)
*    2.1.3 Reconstruct coordinates from the buffer 
	  do K=0,NUMTASK-1
	      NPC=NAP(K)
              do J=NAB(K),NAE(K)  
	        I=J+2*NABS(K)
	        SX(J)=BUF2(I)
	        SY(J)=BUF2(I+NPC)
	        SZ(J)=BUF2(I+2*NPC)
	      end do
	  end do
	end if
*
*    2.2  Case of equal number of atoms on processors
*    ------------------------------------------------
	if(IND.eq.2)then
*  in DOUBLE - each short time step 
	  NAPC3=NAPC*3
	  ISH=NABS(TASKID)                
	  if(NBUFF.le.NAPC3)then
	    write(*,*)' Increase NBUFF in prcm.h to ',NAPC3
	    call FINAL
	  end if
*    2.2.1 Write in buffer
	  do I=1,NAPC    
	    BUFF(I)=SX(I+ISH)
	    BUFF(I+NAPC)=SY(I+ISH)
	    BUFF(I+2*NAPC)=SZ(I+ISH)
	  END DO 
*    2.2.2 Communicate
        call MPI_ALLGATHER(BUFF,NAPC3,MPI_REAL,BUF2,NAPC3,
     +  MPI_REAL,MPI_COMM_WORLD,ierr)
C   This normally not necessary
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
*    2.2.3 Restore coordinates from the buffer
	  do K=0,NUMTASK-1
            do J=NAB(K),NAE(K)  
	      I=J+2*NABS(K)
	      SX(J)=BUF2(I)
	      SY(J)=BUF2(I+NAPC)
	      SZ(J)=BUF2(I+2*NAPC)
	    end do
	end do
      end if    ! IF(IAL...
C
      timen	= timen+cputime(timen0) 
      return
      end            
* 
*    3. Communicate atom velocities to the master node
*    ----------------------------------------
C    This subroutine communicate velocities of atoms to the master node
C    Each node operate with atoms from NAB(TASKID) to NAE(TASKID),
C    but for evaluation of some averages the master node should know 
C    velocities of all the atoms.
C    The are two regimes of this subroutine defined by parameter IND:
C      IND=1  - atoms distributed non-equally between the nodes
C               (case of constrained dunamics)
C      IND=2  - each node (except the last one) has equal number of atoms
C               This is used in DOUBLE subroutine 
C
C============ SINHV ===================================
C 
	subroutine SINHV(IND)
	include "prcm.h"
	include "mpif.h"      
	real*4 BUFF(8),BUF2(8)
        timen0	= cputime(0.d0)
	ISH=NAB(TASKID) 
*
*   3.1 Case of non-equal number of atoms on processors
*   ---------------------------------------------------
	if(IND.eq.1)then 
	  NAPP=NAP(TASKID)             
	  call MPI_GATHERV(VX(ISH),NAPP,MPI_DOUBLE_PRECISION,VX,
     +NAP,NABS,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie)
	  call MPI_GATHERV(VY(ISH),NAPP,MPI_DOUBLE_PRECISION,VY,
     +NAP,NABS,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie)
	  call MPI_GATHERV(VZ(ISH),NAPP,MPI_DOUBLE_PRECISION,VZ,
     +NAP,NABS,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie)
*
*   3.2 Case of equal number of atoms on processors
*   -----------------------------------------------
	else   ! IND=2
	  call MPI_GATHER(VX(ISH),NAPC,MPI_DOUBLE_PRECISION,VX,
     +NAPC,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie)
	  call MPI_GATHER(VY(ISH),NAPC,MPI_DOUBLE_PRECISION,VY,
     +NAPC,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie)
	  call MPI_GATHER(VZ(ISH),NAPC,MPI_DOUBLE_PRECISION,VZ,
     +NAPC,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie) 
	end if
* 
*   3.3 Synchronize all the distancies and thermostat parameters  
*   ------------------------------------------------------------
C   In principle, each node solve the same set of Nose-Hoover equations 
C   with the same parameters, and should synchroniously change them.
C   This is only to ensure that they are the same on all the nodes
	if(TASKID.eq.MAST)then
	  BUFF(1)=SC
	  BUFF(2)=BOXL
	  BUFF(3)=BOYL
	  BUFF(4)=BOZL
	  BUFF(5)=SCL
	  BUFF(6)=SCLX
	  BUFF(7)=SCLY
	  BUFF(8)=SCLZ
	end if
	call MPI_BCAST(BUFF,8,MPI_REAL,MAST,
     +MPI_COMM_WORLD,ie)	
	SC=BUFF(1)
	BOXL=BUFF(2)
	BOYL=BUFF(3)
	BOZL=BUFF(4)
	SCL=BUFF(5)
	SCLX=BUFF(6)
	SCLY=BUFF(7)
	SCLZ=BUFF(8)
	call RECLEN(ODIN,ODIN,ODIN)
        timen	= timen+cputime(timen0) 
      return
      end            
*
*    4.  Sum up fast forces
*    ----------------------
C   This procedure sum up fast forces (HX) acting on each particle
C   The result is send on "home" node for each particle
C   Contributions to energies and virials from all the nodes are 
C   also summed up
C
C=================== CM_ADDF ======================================
C
	subroutine CM_ADDF(LAVER)
	include "prcm.h" 
	include "mpif.h"
	real*4 BUFF(NBUFF),BUF2(NBUFF) 
	logical LAVER
        timen0	= cputime(0.d0)
*
*   4.1 Sum up energies and virials
* 
C  This ie done only at the last short time step, because all these
C  variables are not used during intermidiate time steps          
	if(LAVER)then 
	  BUFF(1)=PE2             !  LJ energy
	  BUFF(2)=PELS2           !  electrostatic energy
	  BUFF(3)=PINT            !  bond/angle intramolec. energy
	  BUFF(4)=VIR2            !  LJ virial
	  BUFF(5)=WIRSS           !  correction to molec. virial
	  BUFF(6)=VIRB            !  bond virial
	  BUFF(7)=VIRFX           !  intramolec. virial projections
	  BUFF(8)=VIRFY
	  BUFF(9)=VIRFZ
	  ICNT=9+MOLINT           
	  do INT=1,MOLINT
	    BUFF(9+INT)=POTE1(INT)
	    BUFF(ICNT+INT)=POTL1(INT)
	  end do
	  ICNT=ICNT+MOLINT
	  ICNT2=ICNT+NTYPES
	  do I=1,NTYPES
	    BUFF(ICNT+I)=PES141(I)
	    BUFF(ICNT2+I)=PSR141(I)
	  end do
	  ICNT=ICNT2+NTYPES
          if(ISTAV.gt.0)then
	  if(NRBB.gt.0)then 
	    ICNT2=ICNT+NRBB
	    do I=1,NRBB
	      BUFF(ICNT+I)=BR(I)
	      BUFF(ICNT2+I)=BE(I)
	    end do
	    ICNT=ICNT2+NRBB
	  end if
	  if(NRAA.gt.0)then 
	    ICNT2=ICNT+NRAA
	    do I=1,NRAA
	      BUFF(ICNT+I)=AR(I)
	      BUFF(ICNT2+I)=AE(I)
	    end do
	    ICNT=ICNT2+NRAA
	  end if
	  if(NRTT.gt.0)then 
	    ICNT2=ICNT+NRTT
	    do I=1,NRTT
	      BUFF(ICNT+I)=TR(I)
	      BUFF(ICNT2+I)=TE(I)
	    end do
	    ICNT=ICNT2+NRTT
	  end if
          end if
	  ICNTS=ICNT  
	  if(ICNT.gt.NBUFF)then
	    write(*,*) 'increase NBUFF to ',ICNT
	    call FINAL
	  end if       
          call MPI_REDUCE(BUFF,BUF2,ICNT,MPI_REAL,MPI_SUM,
     $    MAST,MPI_COMM_WORLD,ierr)
 	  if(TASKID.eq.MAST)then 
	    PE2=BUF2(1)
	    PELS2 = BUF2(2)
	    PINT = BUF2(3)
            VIR2= BUF2(4)
	    WIRSS=BUF2(5) 
	    VIRB=BUF2(6)
	    VIRFX=BUF2(7)        
	    VIRFY=BUF2(8)        
	    VIRFZ=BUF2(9)        
	    ICNT2=9+MOLINT
	    do INT=1,MOLINT
	      POTES(INT)=POTES(INT)+BUF2(9+INT)
	      POTLJ(INT)=POTLJ(INT)+BUF2(ICNT2+INT)
	    end do                  
	    ICNT=ICNT2+MOLINT
	    ICNT2=ICNT+NTYPES
	    do I=1,NTYPES
	      PES14(I)=PES14(I)+BUF2(ICNT+I)
	      PSR14(I)=PSR14(I)+BUF2(ICNT2+I)
	    end do   
	    ICNT=ICNT2+NTYPES
            if(ISTAV.gt.0)then
	    if(NRBB.gt.0)then
	      ICNT2=ICNT+NRBB
	      do I=1,NRBB
	        BB(I)=BB(I)+BUF2(ICNT+I)
	        EB(I)=EB(I)+BUF2(ICNT2+I)
	      end do                   
	      ICNT=ICNT2+NRBB
	    end if
	    if(NRAA.gt.0)then
	      ICNT2=ICNT+NRAA
	      do I=1,NRAA
	        AA(I)=AA(I)+BUF2(ICNT+I)
	        EA(I)=EA(I)+BUF2(ICNT2+I)
	      end do                   
	      ICNT=ICNT2+NRAA
	    end if
	    if(NRTT.gt.0)then
	      ICNT2=ICNT+NRTT
	      do I=1,NRTT
	        TT(I)=TT(I)+BUF2(ICNT+I)
	        ET(I)=ET(I)+BUF2(ICNT2+I)
	      end do                   
	      ICNT=ICNT2+NRTT
	    end if
            end if
	    if(ICNT.ne.ICNTS)write(*,*)' sent ',ICNTS,' received ',ICNT
	  end if
        end if !  LAVER
*
*   4.2  Sum up forces
*   ------------------     
 	do J=0,NUMTASK-1
	  I0=3*NABS(J)
	  I1=I0+NAP(J)
	  I2=I1+NAP(J)
	  do K=1,NAP(J)    
	    I=NABS(J)+K
	    BUFF(I0+K)=HX(I)
	    BUFF(I1+K)=HY(I)
	    BUFF(I2+K)=HZ(I)
	  END DO
	end do
      call MPI_REDUCE_SCATTER(BUFF,BUF2,NAP3,MPI_REAL,
     +MPI_SUM,MPI_COMM_WORLD,ierr)
      do J=NAB(TASKID),NAE(TASKID)  
	  I=J-NABS(TASKID)
	  HX(J)=BUF2(I)
	  HY(J)=BUF2(I+NAP(TASKID))
	  HZ(J)=BUF2(I+2*NAP(TASKID))
      end do
      timen	= timen+cputime(timen0)      
      return
      end
*
*    5.  Sum up slow forces
*    ----------------------
C   This procedure sum up slow forces (GX) acting on each particle
C   The result is send on "home" node for each particle
C   Contributions to energies and virials from all the nodes are 
C   also summed up
CM
CM=================== CM_ADDLF ======================================
CM
      subroutine CM_ADDLF
      include "prcm.h" 
      real*4 BUFF(NBUFF),BUF2(NBUFF)
      include "mpif.h"
      timen0	= cputime(0.d0) 
*
*   5.1 Sum up forces
*   -----------------
      do J=0,NUMTASK-1
	I0=3*NABS(J)
	I1=I0+NAP(J)
	I2=I1+NAP(J)
C
C   NAB(I) is the first atom for which node I is responsible
C   NAE(I) is the last atom for which node I is responsible
C   NABS(I)=NAB(I)-1
C   NAP(I)=NAE(I)-NABS(I) is the number of atoms for node I
C   (these are set in module UNITS, file setup.f ) 
C
	do K=1,NAP(J)    
	  I=NABS(J)+K
	  BUFF(I0+K)=GX(I)
	  BUFF(I1+K)=GY(I)
	  BUFF(I2+K)=GZ(I)
	END DO
      end do
      call MPI_REDUCE_SCATTER(BUFF,BUF2,NAP3,MPI_REAL,
     +MPI_SUM,MPI_COMM_WORLD,ierr)
      do J=NAB(TASKID),NAE(TASKID)  
	I=J-NABS(TASKID)
	GX(J)=BUF2(I)
	GY(J)=BUF2(I+NAP(TASKID))
	GZ(J)=BUF2(I+2*NAP(TASKID))
	M	= NNUM(J)
	WIRS	= WIRS-GX(J)*(SX(J)-X(M))-GY(J)*(SY(J)-Y(M))-
     +               GZ(J)*(SZ(J)-Z(M))
      end do
*
*   5.2  Sum up other data
*   ----------------------
      BUFF(1) = VIRX
      BUFF(2) = VIRY
      BUFF(3) = VIRZ
      BUFF(4) = VIR1
      BUFF(5) = PE1
      BUFF(6) = PELS1
      BUFF(7) = SPE
      BUFF(8) = WIRS
      ICNT=8
      do INT=1,MOLINT
	BUFF(INT+ICNT)=POTE1(INT)
	BUFF(INT+MOLINT+ICNT)=POTL1(INT)
      end do
      ICNT=ICNT+2*MOLINT
      ICNT2=ICNT+NTYPES
      ICNT3=ICNT2+NTYPES
      do I=1,NTYPES
	BUFF(I+ICNT)=PES141(I)
	BUFF(I+ICNT2)=PSR141(I)
	BUFF(I+ICNT3)=SPEE(I)
      end do 
      ICNT=ICNT3+NTYPES
      call MPI_REDUCE(BUFF,BUF2,ICNT,MPI_REAL,MPI_SUM,
     $     MAST,MPI_COMM_WORLD,ierr)
      if(TASKID.eq.MAST)then 
	VIRX = BUF2(1)
	VIRY = BUF2(2)
	VIRZ = BUF2(3)
	VIR1 = BUF2(4)
	PE1  = BUF2(5)
	PELS1= BUF2(6)
	SPE  = BUF2(7)
        WIRS = BUF2(8)
	ICNT=8
	do INT=1,MOLINT
	  POTES(INT)=POTES(INT)+BUF2(INT+ICNT)
	  POTLJ(INT)=POTLJ(INT)+BUF2(INT+MOLINT+ICNT)
	end do                                 
	ICNT=ICNT+2*MOLINT
	ICNT2=ICNT+NTYPES
	ICNT3=ICNT2+NTYPES
	do I=1,NTYPES
	  PES14(I)=PES14(I)+BUF2(I+ICNT)
	  PSR14(I)=PSR14(I)+BUF2(I+ICNT2)
	  SELFPE(I)=SELFPE(I)+BUF2(I+ICNT3)
	end do
      end if       
      timen	= timen+cputime(timen0) 
      return
      end         
*
*    6.  Sum up all forces
*    ---------------------
C   This procedure sum up all forces (FX) acting on each particle
C   The result is send on "home" node for each particle
C   Contributions to energies and virials from all the nodes are 
C   also summed up
CM
CM=================== CM_ADDLA ======================================
CM
      subroutine CM_ADDLA
      include "prcm.h" 
      real*4 BUFF(NBUFF),BUF2(NBUFF)
      include "mpif.h"
      timen0	= cputime(0.d0) 
      do J=0,NUMTASK-1
	I0=3*NABS(J)
	I1=I0+NAP(J)
	I2=I1+NAP(J)
	do K=1,NAP(J)    
	  I=NABS(J)+K
	  BUFF(I0+K)=GX(I)+HX(I)
	  BUFF(I1+K)=GY(I)+HY(I)
	  BUFF(I2+K)=GZ(I)+HZ(I)
	  M	= NNUM(I)
	  WIRS	= WIRS-GX(I)*(SX(I)-X(M))-GY(I)*(SY(I)-Y(M))-
     +               GZ(I)*(SZ(I)-Z(M))
	END DO
      end do
      call MPI_REDUCE_SCATTER(BUFF,BUF2,NAP3,MPI_REAL,
     +MPI_SUM,MPI_COMM_WORLD,ierr)
      do J=NAB(TASKID),NAE(TASKID)  
	I=J-NABS(TASKID)
	FX(J)=BUF2(I)
	FY(J)=BUF2(I+NAP(TASKID))
	FZ(J)=BUF2(I+2*NAP(TASKID))
      end do
      BUFF(1) = VIRX
      BUFF(2) = VIRY
      BUFF(3) = VIRZ
      BUFF(4) = VIR1
      BUFF(5) = VIR2
      BUFF(6)= VIRB
      BUFF(7)= VIRFX
      BUFF(8)= VIRFY
      BUFF(9)= VIRFZ
      BUFF(10)=WIRSS
      BUFF(11) = VIRA(1)
      BUFF(12) = VIRA(2)
      BUFF(13) = VIRA(3)
      BUFF(14) = PE1
      BUFF(15) = PELS1
      BUFF(16) = PE2
      BUFF(17) = PELS2
      BUFF(18) = PINT
      BUFF(19) = SPE
      BUFF(20) = WIRS
      ICNT=20
      do INT=1,MOLINT
	BUFF(INT+ICNT)=POTE1(INT)
	BUFF(INT+MOLINT+ICNT)=POTL1(INT)
      end do
      ICNT=ICNT+2*MOLINT
      ICNT2=ICNT+NTYPES
      ICNT3=ICNT2+NTYPES
      do I=1,NTYPES
	BUFF(I+ICNT)=PES141(I)
	BUFF(I+ICNT2)=PSR141(I)
	BUFF(I+ICNT3)=SPEE(I)
      end do 
      ICNT=ICNT3+NTYPES
      if(ISTAV.gt.0)then
      if(NRBB.gt.0)then 
	ICNT2=ICNT+NRBB
	do I=1,NRBB
	  BUFF(ICNT+I)=BR(I)
	  BUFF(ICNT2+I)=BE(I)
	end do
	ICNT=ICNT2+NRBB
      end if
      if(NRAA.gt.0)then 
	    ICNT2=ICNT+NRAA
	    do I=1,NRAA
	      BUFF(ICNT+I)=AR(I)
	      BUFF(ICNT2+I)=AE(I)
	    end do
	    ICNT=ICNT2+NRAA
      end if
      if(NRTT.gt.0)then 
	    ICNT2=ICNT+NRTT
	    do I=1,NRTT
	      BUFF(ICNT+I)=TR(I)
	      BUFF(ICNT2+I)=TE(I)
	    end do
	    ICNT=ICNT2+NRTT
      end if
      end if
      ICNTS=ICNT  
      if(ICNT.gt.NBUFF)then
	    write(*,*) 'increase NBUFF to ',ICNT
	    call FINAL
      end if       
      call MPI_REDUCE(BUFF,BUF2,ICNT,MPI_REAL,MPI_SUM,
     $     MAST,MPI_COMM_WORLD,ierr)
      if(TASKID.eq.MAST)then 
	VIRX = BUF2(1)
	VIRY = BUF2(2)
	VIRZ = BUF2(3)
	VIR1 = BUF2(4)
	VIR2 = BUF2(5)
	VIRB = BUF2(6)
	VIRFX= BUF2(7)        
	VIRFY= BUF2(8)        
	VIRFZ= BUF2(9)
        WIRSS= BUF2(10)
        VIRA(1)=BUF2(11)
        VIRA(2)=BUF2(12)
        VIRA(3)=BUF2(13)
	PE1  = BUF2(14)
	PELS1= BUF2(15)
	PE2  = BUF2(16)
	PELS2= BUF2(17)
	PINT = BUF2(18)
	SPE  = BUF2(19)
        WIRS = BUF2(20)
	ICNT=20
	do INT=1,MOLINT
	  POTES(INT)=POTES(INT)+BUF2(INT+ICNT)
	  POTLJ(INT)=POTLJ(INT)+BUF2(INT+MOLINT+ICNT)
	end do                                 
	ICNT=ICNT+2*MOLINT
	ICNT2=ICNT+NTYPES
	ICNT3=ICNT2+NTYPES
	do I=1,NTYPES
	  PES14(I)=PES14(I)+BUF2(I+ICNT)
	  PSR14(I)=PSR14(I)+BUF2(I+ICNT2)
	  SELFPE(I)=SELFPE(I)+BUF2(I+ICNT3)
	end do
        ICNT=ICNT3+NTYPES
        if(ISTAV.gt.0)then
	if(NRBB.gt.0)then
	    ICNT2=ICNT+NRBB
	    do I=1,NRBB
	      BB(I)=BB(I)+BUF2(ICNT+I)
	      EB(I)=EB(I)+BUF2(ICNT2+I)
	    end do                   
	    ICNT=ICNT2+NRBB
	end if
	if(NRAA.gt.0)then
	    ICNT2=ICNT+NRAA
	    do I=1,NRAA
	      AA(I)=AA(I)+BUF2(ICNT+I)
	      EA(I)=EA(I)+BUF2(ICNT2+I)
	    end do                   
	    ICNT=ICNT2+NRAA
	end if
	if(NRTT.gt.0)then
	    ICNT2=ICNT+NRTT
	    do I=1,NRTT
	      TT(I)=TT(I)+BUF2(ICNT+I)
	      ET(I)=ET(I)+BUF2(ICNT2+I)
	    end do                   
	    ICNT=ICNT2+NRTT
	end if
        end if
        if(ICNT.ne.ICNTS)write(*,*)' sent ',ICNTS,' received ',ICNT
      end if       
      timen	= timen+cputime(timen0) 
      return
      end         
*
*==================== ADD_RDF ======================================
*                                                                   
	subroutine ADD_RDF(IRDF,IRDF1,MAXS,MAXGR,MAST)
	include "mpif.h"
	integer IRDF(MAXS,MAXGR),IRDF1(MAXS,MAXGR) 
	LEN=MAXS*MAXGR
      call MPI_REDUCE(IRDF,IRDF1,LEN,MPI_INTEGER,MPI_SUM,MAST,
     $     MPI_COMM_WORLD,ierr)                 
	return
	end
*
*===================== GATHER =======================================
*                                                                    
	subroutine GATHER(NRNOD)
	include "prcm.h"
	include "mpif.h" 
	integer AUX(NBDMAX*NTOT),NRNOD(NTPS)
	do ITYP=1,NTYPES	
	  IABEG=ISADDR(ITYP)+1
	  IAEND=ISADDR(ITYP+1) 
	  NRN  =NRNOD(ITYP)       
*  prepare to send       
	  if(NRN.eq.TASKID)then
	    J=0
	    do I=IABEG,IAEND
	      J=J+1
	      AUX(J)=NNBB(I)
	      do INB=1,NNBB(I)
	        J=J+1
	        AUX(J)=INBB(I,INB)
	      end do
	    end do  
	    JBUF=J
	  end if
*  transfer buffer size         
        call MPI_BCAST(JBUF,1,MPI_INTEGER,NRN,MPI_COMM_WORLD,ie) 
*  transfer data
	  call MPI_BCAST(AUX,JBUF,MPI_INTEGER,NRN,MPI_COMM_WORLD,ie)
*  unpack data
	  if(NRN.ne.TASKID)then
	    J=0
	    do I=IABEG,IAEND
	      J=J+1
	      NNBB(I)=AUX(J)
	      do INB=1,NNBB(I)
	        J=J+1
	        INBB(I,INB)=AUX(J)
	      end do
	    end do  
	  end if
	  if(J.ne.JBUF)write(*,*)' in GATHER: send ',JBUF,' received ',J,
     +' for type ',ITYP,' node ',TASKID
      end do
	return
	end 
*
*======================== TRANSRDF ========================
*
 	subroutine TRANSRDF(NS,MXGRS,MAST,TASKID,NGRI,IGRI,MGRI,IAUX)
	include "mpif.h"
	integer NGRI(NS),IGRI(NS,MXGRS),MGRI(NS,MXGRS),IAUX(NS),TASKID
	call MPI_BCAST(NGRI,NS,MPI_INTEGER,MAST,MPI_COMM_WORLD,ie)
	do I=1,MXGRS
	  if(TASKID.eq.MAST)then 
	    do J=1,NS
	      IAUX(J)=IGRI(J,I)
	    end do
	  end if  
	  call MPI_BCAST(IAUX,NS,MPI_INTEGER,MAST,MPI_COMM_WORLD,ie) 
	  if(TASKID.eq.MAST)then 
	    do J=1,NS
	      IAUX(J)=MGRI(J,I)
	    end do
	  else
	    do J=1,NS
	      IGRI(J,I)=IAUX(J)
	    end do
	  end if  
	  call MPI_BCAST(IAUX,NS,MPI_INTEGER,MAST,MPI_COMM_WORLD,ie)  
	  if(TASKID.ne.MAST)then 
	    do J=1,NS
	      MGRI(J,I)=IAUX(J)
	    end do
	  end if
	end do  
	return  
	end                 
*
*============================= ALLSUM ==================================
*     
	subroutine ALLSUM(A,B,N)
	real*8 A(N),B(N)
	include "mpif.h"
	  call MPI_ALLREDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     MPI_COMM_WORLD,ierr)     
	return
	end
*
*============================ SUMMA ===============================
*                                                                  
      subroutine SUMMA(A,B,N,MAST)
	real*8 A(N),B(N)
	include "mpif.h"
      call MPI_REDUCE(A,B,N,MPI_DOUBLE_PRECISION,MPI_SUM,
     $     MAST,MPI_COMM_WORLD,ierr)
	return
	end   
*
*============================ BCAST ===============================
*                                                                  
      subroutine BCAST(A,N,MAST)
      real*8 A(N)
      include "mpif.h"
      call MPI_BCAST(A,N,MPI_DOUBLE_PRECISION,MAST,MPI_COMM_WORLD,ie)	
      return
      end   
*
*========================== BARRIER ================================
*
      subroutine BARRIER
      include "mpif.h"
      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      return
      end
*
*============================= FINAL ========================
*                                                            
	subroutine FINAL
	include "mpif.h"
        call MPI_BARRIER(MPI_COMM_WORLD,IERR)
	call MPI_FINALIZE(ie)
C        write(*,*)' Program finalized'
        close(6)
	stop 
	end
*
*============================= ABORT ========================
*                                                            
	subroutine ABORT
	include "mpif.h"
        call MPI_ABORT(MPI_COMM_WORLD,IERR)
	call MPI_FINALIZE(ie)
        write(*,*)' Program aborted'
        close(6)
	stop 'MPI finalized' 
	end
