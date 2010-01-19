*     PART 4
*
*     File mdstep.f
*     -------------
*
C     This file contains subroutines responsible for integration of equation 
C     of motion:
C
C     SULTAN     - reversible Verlet algorithm with optiohnal 
C                  constrained dynamics for each molecule type
C
C     DOUBLE     - double time step reversible algorithm of Martyna et al 
C                  (Mol.Phys., v.87(5), pp.1117-1157)
C
C     SCALING    - realises specified statistical ensemble (NVT, NPT or NVE 
C                  with periodical velocity scaling)
C
C     SHEJK      - SHAKE algorithm for constrained dynamics
C
*
*  
*   1. Leap frog algorithm with constraint dynamics (SULTAN)
*   --------------------------------------------------------
*
C   This subroutine realises leap-frog algorithm for MD-simulation
C   input: x(t), v(t-1/2)
C   The unit is called from main() after completion of initialisation work
C
      SUBROUTINE SULTAN
      include "prcm.h"  
*
*   1.1 Local definitions
*
      LOGICAL DONE,LRDF
      real*8 VIRR(3)
      DATA LRDF/.false./
*
*   1.2 Preparation
*   ---------------
*   
*   1.2.1  Compose list of neigbours       
      call CHCNB(LRDF)        ! l-forces.f  
*   1.2.2  Calculate temperature
      call GETEMP(.true.)
*   1.2.3  Calculate force for the first time
C   (no differense in firces in this case)
      call ALLFORCE                ! forces.f
*   1.2.4  Organize cycle for integration equation of motion	
      LRDF=.true.
C   NSTEP  - step number in this run
C   NSTTOT - previously made steps
C   MSTEP  - step number fron the begining of simulation
      NSTEP=0
 1    NSTEP=NSTEP+1
      MSTEP=NSTTOT+NSTEP
*
*   1.3  Integration equations of motion
*   ------------------------------------
C   The reversible integration scheme applied:
C   First - integration thermostat-related variables to time t+1/2
C       (including correction of positions to t due to thermostat)
C   Then - Verlet integration of positions and velocoties to t+1:
C       a) velocities at t+1/2 from forces at t
C       b) positions at t from velocities at t+1/2 (including constrain)
C       c) forces at t
C       d) velocities at t from forces at t 
C   Third - integration of thermostat-related variables to t+1
      timeg0	= cputime(0.d0)
*   1.3.1  Remember old positions for constrained dynamics
*   ------------------------------------------------------
C   In the parallel algorithm each node is responsible for motion
C   of atoms from NAB(TASKID) to NAE(TASKID). Arrays NAB and NAE are 
C   defined in units.f 
      DO I         = NAB(TASKID),NAE(TASKID)
	ITYP	= ITYPE(I)
C   LMOVE define which molecules can move (normally, LMOVE=.true.)
	if(LMOVE(ITYP))then 
          OX(I)        = SX(I)
          OY(I)        = SY(I)
          OZ(I)        = SZ(I)
        end if
      end do
*
*   1.3.2 Verlet step
*   -----------------
      DO I         = NAB(TASKID),NAE(TASKID)
	ITYP	= ITYPE(I)
	if(LMOVE(ITYP))then 
*   1.3.2.1 Unconstrained velocities at t+1/2 
          VX(I)        = VX(I)+TSTEP*FX(I)
          VY(I)        = VY(I)+TSTEP*FY(I)
          VZ(I)        = VZ(I)+TSTEP*FZ(I)
*   1.3.2.2  Unconstrained postions  at t+1 
          SX(I)        = SX(I)+TSTEP*VX(I)*MASSDI(I)
          SY(I)        = SY(I)+TSTEP*VY(I)*MASSDI(I)
          SZ(I)        = SZ(I)+TSTEP*VZ(I)*MASSDI(I)
        end if
      end do
*
*   1.3.2.3 Apply constraints for positions and velocities
C  IMB,...,IME - these numbers defined which molecules this node 
C                is responsible for
C  Arrays of forces FX,FY,FZ are used in SHEJK as temporarly storage
C  
      IMB=NNUM(NAB(TASKID))
      IME=NNUM(NAE(TASKID))
      VIRAN(1)=0.
      VIRAN(2)=0.
      VIRAN(3)=0.
      do IMOL=IMB,IME  
	CALL SHEJK(IMOL)
      end do
*  first time - correct velocities
      if(MSTEP.eq.1)then
        call GETEMP(.true.)
        FTMP        = SQRT(TRTEMP/TEMP)
        DO     N    = 1,NSTOT
          VX(N)       = VX(N)*FTMP
          VY(N)       = VY(N)*FTMP
          VZ(N)       = VZ(N)*FTMP
        end do
      end if
      timeg	= timeg+cputime(timeg0)
C  transfer virials from constraints to the master node
      call SUMMA(VIRAN,VIRA,3,MAST)
      timeg0	= cputime(0.d0)
      call SCALING(3)
      timeg	= timeg+cputime(timeg0)
*
*   1.3.2.4 Broadcast new coordinates and velocities
C   MPI - only
      call SINHR(1)             !  mpi.f or scalar.f
C   Velocities are collected only on the master node for analysis
C   (of course, each node keeps velocities of own atoms)	
      call SINHV(1)        ! mpi.f or scalar.f
*   1.3.2.5 Apply PBC and calculate new centre of masses
*   --------------------------------------------------      
      timeg0	= cputime(0.d0)
      CALL GETCOM     ! supl.f
      timeg	= timeg+cputime(timeg0)
*   1.3.2.6 After a given number of steps recalculate list of neigbours      
*   -----------------------------------------------------------------
      if(mod(MSTEP,ICHNB).eq.0)then
        call CHCNB(LRDF)
      end if 
*   1.3.3 Recalculate forces for the next step
      call ALLFORCE                ! forces.f
*      
*   1.4 Average collection
*   ----------------------
      timeg0	= cputime(0.d0)
* remove excess COM momenta if specified
      if(ICMOM.gt.0.and.mod(MSTEP,ICMOM).eq.0)call CHKMOM
      CALL GETKIN     ! supl.f
*   1.4.1 Calculate dipole moments
      if(LCF)CALL DIPOLE
*   1.4.2 Calculate tcf      
      if(NHIST+1.ge.IHIST)CALL GETTCF
*   1.4.3 Collect other averages     
C   electrostatic energy: sum up "fast" and "slow" contributions      
        PEST = PELS1+PELS2+SPE
        PE = PE1+PE2+PEST
	if(mod(MSTEP,IAVER).eq.0)CALL GETAVR(1)
*   1.4.4 Dump trajectory	
	if(mod(MSTEP,NTREK).eq.0)CALL TRACE
*   1.4.5 Calculate some quantities for current output  
C   Total intermolecular energy
        PENL=PE*PERMOL
C   Electrostatic potential energy      
        PEL=PEST*PERMOL
C   Intramolecular potential energy      
        PENS=PINT*PERMOL
        PEKI=TKE*PERMOL                 !   kinetic energy
        ETOT=PENL+PENS+PEKI                  !   total energy
        POTEN	= PENL+PENS                   !  total potential energy
        WIRSUM=WIRS+WIRSS
C   These are diffrent contribution to pressure:	
        TRYCKE	= PEST*UNITP/(3.*VOL)    ! electrostatic contr.
        TRYCKL	= (VIR1+VIR2)*UNITP/(3.*VOL)     ! LJ
        TRYCKB	= VIRB*UNITP/(3.*VOL)    ! bond (flexible)
        TRYCKA	= (VIRA(1)+VIRA(2)+VIRA(3))*UNITP/(3.*VOL) ! shejk-correction 
        TRYCKS	= WIRSUM*UNITP/(3.*VOL) ! orientational correction
        TRYCKD	= VIRD*UNITP/(3.*VOL) ! long-range correction
        DENS	= TOTMAS*1.d-3/(VOL*UNITL**3) !  current density
*   1.4.7  Current output
C   Made only by master node 
      if(TASKID.eq.MAST)then
C OUTPUT for visualisation
        if(LVISUAL)then
          if(mod(MSTEP,IAVER).eq.0)then
            write(*,'(a)')'@mm coord.start'
            do I=1,NSTOT
              write(*,'(a,3f12.4)')'@mm ',SX(I),SY(I),SZ(I)
            end do
            write(*,'(a)')'@mm coord.end'
          end if
          WRITE(*,'(a,I7,4(1X,F9.3),2(1X,f9.2),f9.4,3f10.4)') 
     +    '@mm E ',MSTEP,POTEN,PENS,PEL,ETOT,TEMP,TRYCKM,DENS,
     +     BOXL,BOYL,BOZL
        end if  ! LVISUAL
        if(IPRINT.ge.5.and.mod(MSTEP,IAVER).eq.0)
     +   WRITE(6,'(I8,4(1X,F9.3),2(1X,f9.2),f9.4)') 
     +    MSTEP,POTEN,PENS,PEL,ETOT,TEMP,TRYCKM,DENS
        if(IPRINT.ge.6)write(*,'(2f9.6,6f10.2)')
     +   SC,SCL,TRYCK,TTR,TROT,TINT,TEMPV
	if(IPRINT.ge.6.and.LSEP)write(*,'(4f11.1,3f11.3)')
     + TRYCKX,TRYCKY,TRYCKZ,(TRYCKX+TRYCKY+TRYCKZ)/3.,BOXL,BOYL,BOZL
          if(IPRINT.ge.6.and.IEXT.eq.1)write(*,'(2(a,e14.6))')
     +    ' E abs: ',(EABSS+EABS)*PERMOL
 	if(IPRINT.ge.7.and.mod(MSTEP,IAVER).eq.0)
     +  WRITE(6,'(8f10.2)')TRYCKL,TRYCKE,TRYCKB,TRYCKA,TRYCKS,TRYCKD
C   temperature and termostat parameter for each species
          if(IPRINT.ge.6.and.LSCTD)
     +    write(*,'(8f10.5)')(TEMPR(I),SCM(I),I=1,NTYPES)
*** Put you own subroutine for evaluation of whatever you want here:
***     call USER
***
      end if
*   1.4.8  Intermediate output after "NAVER" steps   
      if(mod(MSTEP,NAVER).eq.0)call GETAVR(2)
*   1.4.9  Each "NDUMP" steps dump restart file
C   Important: this should be called by ALL nodes!
      if(LDMP.and.mod(NSTTOT+NSTEP,NDUMP).eq.0)call RESTRT(2)
      timeg	= timeg+cputime(timeg0)
*
*   1.4.10 Return to p. 1.3.1
*   Check STOP file
      if(LVISUAL)then
        open(unit=77,file='MD_STOP',status='old',err=10)
        close(77)
        return
      end if
 10   if(NSTEP.lt.NSTEPS)go to 1
*<---------------------
*
      RETURN
      END
*
*====================================================================
*
*   2. DOUBLE (su) - Double time step algorithm for flexible molecules
*   ------------------------------------------------------------------
C   This subroutine realises reversible double time step algorithm 
C   (M.Tuckerman, B.J.Berne, G.J.Martyna, J.Chem.Phys. v.97, p.1990 (1992))
C   All the forces are divided into two groups: "fast" and "slow"
C   Fast forces are: 
C      - due to covalent bond, angles and dihedral angles;
C      - Lennard-Jones and real-space electrostatic within cutoff defined
C   by SHORT parameter (normally, 5A)                 
C
C   Slow forces are:
C      - Lennard-Jones and real-space electrostatic on distances
C   between SHORT and CUTOFF 
C      - Reciprocal-space electrostatics              
C   Fast forces are recalculated each short time step (TSETP/NFREQ)
C   and slow forces are recalculated each long time step
C  
      SUBROUTINE DOUBLE
      include "prcm.h"
*
*   2.1 Local definitions
*
C   LRDF defines wheter to calculate RDFs
      LOGICAL LRDF,DONE  
C   this is temporary array to keep long-range contributions to virial
      real*8 TRLL(3)
C   RDF should not be calculated during initialisation stage (2.2)
      LRDF=.false.
      NFREQ2=NFREQ/2+1
      VIRA(1)=0.
      VIRA(2)=0.
      VIRA(3)=0.
*
*   2.2 Preparation
*   ---------------
*   
*   2.2.1  Compose list of neigbours       
      call CHCNB(LRDF)      ! mdstep.f 
*   2.2.2  Calculate centre-of-masses
      CALL GETCOM                ! supl.f 
*   2.2.3  Calculate all forces (should be done before entering the cycle)
      LRDF=.true.  
C   slow forces  (important: slow forces calculated before fast ones!)
      call SLFORCES
C   Fast forces
      call FFORCES(.false.)
*   2.2.4  Calculate different contributions to pressure and energy
C     (this is done mostly for control purposes)
      timeg0	= cputime(0.d0)
      if(IPRINT.ge.7)then
	PELL=PELS1+PELS2+SPE
        PENL=(PE1+PE2+PELL)*PERMOL
        PEL=PELL*PERMOL
        PENS=PINT*PERMOL
        PEKI=TKE*PERMOL
        ETOT=PENL+PENS+PEKI
        WRITE(6,'(4(a6,f12.5))')
     +' Eex= ',PENL,' Ein= ',PENS,' Eel= ',PEL,' Etot=',ETOT 
      end if
*
*   2.3  Double time step algorithm
*   -------------------------------
*   2.3.1 Organize the cycle (till 2.4.9)
      NSTEP=0
 1    NSTEP=NSTEP+1
      MSTEP=NSTTOT+NSTEP
*   2.3.2 First thermostat integration: to half of the long time step
      call SCALING(1)
*   2.3.3 First integration of slow forces
C     Velocities are corrected as due to slow forces acting during 
C     half of the long time step
      DO I           =  NAB(TASKID),NAE(TASKID)
	ITYP	= ITYPE(I)
C   LMOVE define which molecules can move (normally, LMOVE=.true.)
	if(LMOVE(ITYP))then 
          VX(I)          =   VX(I)+HTSTEP*GX(I)
          VY(I)          =   VY(I)+HTSTEP*GY(I)
          VZ(I)          =   VZ(I)+HTSTEP*GZ(I)
	end if
      END DO! OF I
*
*    2.3.4  Integration of the fast forces
*    ------------------------------------- 
*    2.3.4.1  Cycle over short time steps
C  This set zeros for current averages bond/angles and their energies
C  These averages are taken over the molecules and over short time
C  steps within this large step 
      call ZEROAV
      DO NSTEB       =  1,NFREQ
        DO I           =  NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then 
*    2.3.4.3  Velocities at t(short)+1/2
            VX(I) =   VX(I)+HTSTEB*HX(I)         ! V(t+1/2)
            VY(I) =   VY(I)+HTSTEB*HY(I)
            VZ(I) =   VZ(I)+HTSTEB*HZ(I)
*    2.3.4.4  Coordinates at t(short)+1
            SX(I) =   VX(I)*TSTEB*MASSDI(I)+SX(I)     ! X(t+1)
            SY(I) =   VY(I)*TSTEB*MASSDI(I)+SY(I)
            SZ(I) =   VZ(I)*TSTEB*MASSDI(I)+SZ(I)
	  end if
        END DO! OF I
        timeg	= timeg+cputime(timeg0)   
*  2.3.4.5  Broadcast coordinates over nodes  
C   Each node broadcast coordinates of own atoms to all nodes  
	call SINHR(2)  
*  2.3.4.6  Recalculate centre-of-mass coordinates
	CALL GETCOM         ! services.f
*
*  2.3.5 Recalculate forces
*  ------------------------
*  2.3.5.1 After "ICHNB" long time steps - recalculate list of neighbours
C  This is done in the middle of long time step  
	if(NSTEB.eq.NFREQ2.and.mod(MSTEP,ICHNB).eq.0)call CHCNB(LRDF)     
	if(NSTEB.eq.NFREQ)then
*  2.3.5.2 Each long time step - calculate slow forces
          CALL SLFORCES
*  2.3.5.3 Claculate fast forces
C         The logical parameter signals whether to calculate averages 
          call FFORCES(.true.)
        else
          call FFORCES(.false.)
        end if
*   
*   2.3.6 Conclude integration of fast forces
*   -----------------------------------------
C   (continued in the beginning of the cycle, see p.2.3.2    
        timeg0	= cputime(0.d0)
        DO I           =  NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then 
            VX(I)          =  VX(I)+HTSTEB*HX(I)     !  V(t+1)
            VY(I)          =  VY(I)+HTSTEB*HY(I)
            VZ(I)          =  VZ(I)+HTSTEB*HZ(I)
          end if
	END DO! OF I      
        call GETEMP(.true.)
	if(IPRINT.ge.8)write(*,'(a,I3,f12.3)')' sh.st.',NSTEB,TEMP
*
      END DO! OF NSTEB 
*
*   2.3.7 Conclude integration of slow forces
*   -----------------------------------------  
*   2.3.7.1 Velocities at t(long)+1
        DO I           = NAB(TASKID),NAE(TASKID)
	  ITYP	= ITYPE(I)
	  if(LMOVE(ITYP))then 
            VX(I)          =  VX(I)+HTSTEP*GX(I)
            VY(I)          =  VY(I)+HTSTEP*GY(I)
            VZ(I)          =  VZ(I)+HTSTEP*GZ(I)
	  end if
        END DO! OF I
      timeg	= timeg+cputime(timeg0)
*   2.3.7.2 Conclude integration of NHC
      call SCALING(2)
*   2.3.6.3 Report velocities to the master node
      call SINHV(2)
      timeg0	= cputime(0.d0)
*
*   2.4 calculate averages for intermediate output
*   ----------------------------------------------
* remove excess COM momenta if specified
      if(ICMOM.gt.0.and.mod(MSTEP,ICMOM).eq.0)call CHKMOM
      if(TASKID.eq.MAST)then
*   2.4.1 Different contribution to the kinetic energy 
        CALL GETKIN       ! services.f
*   2.4.2 Dipole moment
        if(LCF)then
          CALL DIPOLE       ! serivices.f
*   2.4.3 Time correlation functons
          if(NHIST+1.ge.IHIST)CALL GETTCF          ! tcf.f
        end if
*   2.4.4 Acccumulate averages 
        PEST = PELS1+PELS2+SPE
        PE = PE1+PE2+PEST
        if(mod(MSTEP,IAVER).eq.0)CALL GETAVR(1)  ! aver.f
*   2.4.5 Dump trajectory
	if(mod(MSTEP,NTREK).eq.0)CALL TRACE      ! restart.f
*   2.4.6 Calculate different contributions to pressure and energy
        PENL = PE*PERMOL                ! intermolecular energy
        PEL= PEST*PERMOL                ! electrostatic energy          
        PENS=PINT*PERMOL                ! intramolecular energy
        PEKI=TKE*PERMOL                 ! kinetic energy
        POTEN=PENL+PENS                 ! total potential energy
        ETOT =POTEN+PEKI                ! total energy   
        WIRSUM=WIRS+WIRSS
        TRYCKE	= PEL*UNITP/(3.*VOL*PERMOL)     ! electrostatic
        TRYCKL	= (VIR1+VIR2)*UNITP/(3.*VOL)    ! LJ
        TRYCKB	= VIRB*UNITP/(3.*VOL)           ! bonds
        TRYCKS	= WIRSUM*UNITP/(3.*VOL)         ! molecular
*        TRYCKA	= (VIRA(1)+VIRA(2)+VIRA(3))*UNITP/(3.*VOL) ! molecular
        TRYCKD	= VIRD*UNITP/(3.*VOL)           ! long-range corr.
        DENS	= TOTMAS*1.d-3/(VOL*UNITL**3)
*   2.4.7 Current output
C OUTPUT for visualisation
        if(LVISUAL)then
          if(mod(MSTEP,IAVER).eq.0)then
            write(*,'(a)')'@mm coord.start'
            do I=1,NSTOT
              write(*,'(a,3f12.4)')'@mm ',SX(I),SY(I),SZ(I)
            end do
            write(*,'(a)')'@mm coord.end'
          end if
          WRITE(*,'(a,I7,4(1X,F9.3),2(1X,f9.2),f9.4,3f10.4)') 
     +    '@mm E ',MSTEP,POTEN,PENS,PEL,ETOT,TEMP,TRYCKM,DENS,
     +     BOXL,BOYL,BOZL
        end if   ! if(LVISUAL
        if(mod(MSTEP,IAVER).eq.0)then
          if(IPRINT.ge.5)
     +      WRITE(6,'(I8,4(1X,F9.3),2(1X,f9.2),f9.4)')
     +      MSTEP,POTEN,PENS,PEL,ETOT,TEMP,TRYCK,DENS
          if(IPRINT.ge.6)write(*,'(2f7.4,2f9.2,4f7.0,I7,I9)')
     + SC,SCL,TRYCKM,TTR,TROT,TINT,TEMPV,MBSH,MBLN   
   	  if(IPRINT.ge.6.and.LSEP)write(*,'(4f11.1,3f11.3)')
     + TRYCKX,TRYCKY,TRYCKZ,(TRYCKX+TRYCKY+TRYCKZ)/3.,BOXL,BOYL,BOZL
          if(IPRINT.ge.6.and.IEXT.eq.1)write(*,'(2(a,e14.6))')
     +    ' E abs: ',(EABS+EABSS)*PERMOL,' Ext.f. ',EFIELD
 	  if(IPRINT.ge.7.and.mod(MSTEP,IAVER).eq.0)
     +  WRITE(*,'(8f10.1)')TRYCKL,TRYCKE,TRYCKB,TRYCKS,TRYCKD,TRID
          if(IPRINT.ge.6.and.LSCTD)
     +    write(*,'(8f10.5)')(TEMPR(I),SCM(I),I=1,NTYPES)
*** Put you own subroutine for evaluation of whatever you want here:
***     call USER
***
        end if  ! if(mod(MSTEP
      end if    ! if(TASKID ...
*   2.4.8  Intermediate output after "NAVER" steps   
      if(mod(MSTEP,NAVER).eq.0)CALL GETAVR(2)
*   2.4.9  Dump restart file   
C   Important - all nodes
      if(LDMP.and.mod(MSTEP,NDUMP).eq.0)call RESTRT(2)
      timeg	= timeg+cputime(timeg0)
*
*   2.4.10 Conclude MD step
*   Check STOP file (visual regime only)
      if(LVISUAL)then
        open(unit=77,file='MD_STOP',status='old',err=10)
        close(77)
        return
      end if
 10   if(NSTEP.lt.NSTEPS)go to 1
*<---------------------
*
      RETURN
      END
*
*=============== SCALING ==============================================
*
      SUBROUTINE SCALING(IAL)
*
*  3  Realisation of the specified ensemble
*  ----------------------------------------
C  Case NVT/NPT - Nose-Hoover thermostat, "modular-invariant" version 
C  of G.J.Martyna, D.J.Tobais and M.L.Klein,
C  J.Chem.Ohys., v.101(5), p.4177 (1994)
C  
C  NVE - if specified, do simple velocity scaling if temperature 
C        deviate from reference value more that TDELT 
C
C  IAL = 1 - first half-step in double time step algorithm
C  IAL = 2 - second half-step in double time step algorithm
C  IAL = 3 - Leap-frog algorithm
C
*
*  3.1 Definitions
*  ---------------
      include "prcm.h"   
      parameter (FP3=0.5/3.)
      real*8 DPE(4),DTE(NTPS)
      data DKEO/0./,ISKK/0/
*
*  3.2 Case of Leap-frog - correct velocities due to barostat
*
      if(IAL.eq.3.and.LNVT)then
*   3.2.1    Separate thermostat on each molecule type
*            (velocity scaled before thermostat parameter correction)
        if(LSCTD)then
          do I=NAB(TASKID),NAE(TASKID)
	    ITYP = ITYPE(I)
            SCV= 1.-SCM(ITYP)*TSTEP
	    if(LMOVE(ITYP))then
              VX(I) = VX(I)*SCV
              VY(I) = VY(I)*SCV
              VZ(I) = VZ(I)*SCV
            end if
          end do
*   3.2.2  Common thermostat
        else
          SCV = 1.-SC*TSTEP
          do I=NAB(TASKID),NAE(TASKID)
	    ITYP = ITYPE(I)
	    if(LMOVE(ITYP))then 
              VX(I) = VX(I)*SCV
              VY(I) = VY(I)*SCV
              VZ(I) = VZ(I)*SCV
            end if
          end do
        end if
      end if
*
*  3.3 Recalculate temperatures and pressures
*  ---------------
*  3.3.1 Calculate temperature    
      call GETEMP(.true.)
*  3.3.2 Calculate pressure  
C        True only for MAST node
C        This collect virial for "atomic" pressure
      VIRSUM	= VIR1+VIR2+VIRB+PELS1+PELS2+SPE+VIRD+
     +            VIRA(1)+VIRA(2)+VIRA(3)
C        This collect virial for "molecular" pressure
      WIRSUM	= VIR1+VIR2+PELS1+PELS2+SPE+VIRD+WIRS+WIRSS
CD      write(*,'(5e13.5)')VIR1,VIR2,PELS1,PELS2,SPE,VIRD,WIRS,WIRSS
C        Kinetic (ideal gas) contributions to pressure
      TRID	= 2.*TKE*UNITP/(3.*VOL)
      TRIDM	= NOP*BOLTZ*TEMP*1.d-5/(VOL*UNITL**3)
C        Pressure in atm and its projections
      TRYCK	= VIRSUM*UNITP/(3.*VOL)+TRID
      TRYCKM	= WIRSUM*UNITP/(3.*VOL)+TRIDM 
C        Projections are calculated for "atomic" pressure
      TRYCKX	= (VIRX+VIRFX+VIRD/3.+VIRA(1))*UNITP/VOL+TRID
      TRYCKY	= (VIRY+VIRFY+VIRD/3.+VIRA(2))*UNITP/VOL+TRID
      TRYCKZ	= (VIRZ+VIRFZ+VIRD/3.+VIRA(3))*UNITP/VOL+TRID
CD	if(TASKID.eq.MAST)write(*,'(5f13.2)')
CD     +TRYCK,TRYCKX,TRYCKY,TRYCKZ,(TRYCKX+TRYCKY+TRYCKZ)/3.
C
C     Now temperature TEMP and pressure TRYCK are defined for the 
C     current configuration, as well as their components
*
*   3.3.3 Calculate deviation from thermostat T and P
      DKE         = TEMP/TRTEMP-1.d0               !  total temperature 
      do I=1,NTYPES
CD         write(*,*)' temp ',I,TEMPR(I),TASKID
         DTE(I)   = TEMPR(I)/TRTEMP-1.d0           !  for each species
      end do
C    in the case of constrained dynamics, pressure is defined by 
C    "molecule" algorithm - it corresponds scaling of molecular COM
C    for flexible molecules, scaling of atom positions is employed,
C    which corresponds to the "atomic" pressure 
C    Exception: case of separate pressure control in each direction 
C    (LSEP=.true.). Then pressure in each direction is determined from
C    "atomic" pressure. 
      if(LSHEJK)then
	 DPE(1) = TRYCKM-TRPRES
      else
         DPE(1) = TRYCK -TRPRES
      end if
      DPE(2) = TRYCKX-TRPRES
      DPE(3) = TRYCKY-TRPRES
      if(LHEX)then
        DPE(2)=0.5*(DPE(2)+DPE(3))
        DPE(3)=DPE(2)
      end if
      DPE(4) = TRYCKZ-TRPRES
C    This broadcast pressure to all nodes
      call BCAST(DPE,4,MAST)         
*
*   3.4  First integration of Nose-Hoover equations  (t->t+1/2)
*   -----------------------------------------------------------
      if(IAL.eq.1)then
        if(LNVT)then
*   3.4.1 Correction ksi due to thermostat
C         coefficients DQT,DQP and RTP were defined in main.f, p.6.2
C
C    Nose eqn is:    dP/dt   -> F - P*ksi/Q
C                    dksi/dt -> 2*Ekin - Nf*kT
C    Here ksi/Q = SC
C             Q = Nf*kT*tau**2  ; tau -> QT (i.u.)
C                                 Q/(Nf*kT) -> 1./DQT 
C    So     d(SC)/DT = (2*Ekin/Nf*kT-1)*DQT = DKE*DQT
C
C    This is separate thermostat for each species 
          if(LSCTD)then
            SCA = 0.
            do I=1,NTYPES
              SCM(I) = SCM(I) + DTE(I)*DQT*HTSTEP     !  DQT -thermostat mass
              SCA = SCA + SCM(I)*EKIN(I) 
            end do
            SC = SCA/TKE
C   common thermostat
          else
 	    SC = SC + DKE*DQT*HTSTEP
          end if
          if(LNPT)then
*   3.4.2 Barostat corrections (NPT)
*   -------------------------------
*   3.4,2,1 Correction ksi due to barostat
C
C   Additional correction to ksi:
C       d(ksi)=ita**2/W-kT
C       ita/W -> SCL
C     so d(ksi)  ->   SCL**2*(W/Q) - kT/Q 
C                       
            if(LSCTD)then      
              SCA = 0.
              do I=1,NTYPES
                SCM(I) = SCM(I) + (SCL**2*RTP-DQT/FNST)*HTSTEP
                SCA = SCA + SCM(I)*EKIN(I) 
              end do
              SC = SCA/TKE
            else
              SC = SC + (SCL**2*RTP-DQT/FNST)*HTSTEP    !  RTP=DQT/DQP
            end if
*   3.4.2.2 Correction ita 
C         TKE is kinetic energy
            if(LSEP)then
              SCLX = SCLX+DQP*(VOL*DPE(2)+2.*TKE/FNST-SCLX*SC)*HTSTEP
              SCLY = SCLY+DQP*(VOL*DPE(3)+2.*TKE/FNST-SCLY*SC)*HTSTEP
              if(LHEX)then
                 SCLX=0.5*(SCLX+SCLY)
                 SCLY=SCLX
              end if
              SCLZ = SCLZ+DQP*(VOL*DPE(4)+2.*TKE/FNST-SCLZ*SC)*HTSTEP
              SCL = SCLX+SCLY+SCLZ   ! trace(P)
            else
              SCL = SCL + DQP*(3.*VOL*DPE(1)+6.*TKE/FNST-SCL*SC)*HTSTEP
              SCLX=SCL
              SCLY=SCL
              SCLZ=SCL
C    Control fluctuations
	      if(abs(SCL*HTSTEP).ge.0.1)then
	        write(*,*)' too strong fluctuations in NPT algorithm'
                if(NUMTASK.gt.1)write(*,*)' at node ',TASKID
	        write(*,*)' scaling factor ',SCL*HTSTEP
   	        if(IPRINT.ge.7)then
	          write(*,*)'------------------->'
	          write(*,*)SCL,TRYCK,TRYCKM
	          write(*,*)VIR1*UNITP/(3.*VOL),VIR2*UNITP/(3.*VOL),
     + PELS1*UNITP/(3.*VOL),PELS2*UNITP/(3.*VOL),VIRB*UNITP/(3.*VOL)
	          write(*,*)'<-------------------'
	        end if
                ISKK=ISKK+1
                if(ISKK.gt.100)then 
	write(*,*)' !!! repeated failure of NPT-algorithm' 
	write(*,*)' !!! restart program with constant volume'
	write(*,*)' !!! or increase thermostat parameter for pressure'
                  call FINAL
                end if
              end if    ! if(LSEP   
            end if
*   3.4.2.3 Calculate scaling coefficients  
            SCVX = 1.-((1.-3./FNST)*SCL+SC)*HTSTEP
            SCVY = SCVX
            SCVZ = SCVX
            DRCX = TSTEP*SCLX
            DRCY = TSTEP*SCLY
            DRCZ = TSTEP*SCLZ
            DRC  = (DRCX+DRCY+DRCZ)/3.
            DRCV = DRC*HTSTEP           !  These are second order corrections 
            DRCX = 1.+DRCX + 0.5*DRCX**2 !  and (may be?) ommitted
            DRCY = 1.+DRCY + 0.5*DRCY**2 !  and (may be?) ommitted
            DRCZ = 1.+DRCZ + 0.5*DRCZ**2 !  and (may be?) ommitted
*
*   3.4.2.4 Correct coordinates and velocities - flexible molecules
            do I=NAB(TASKID),NAE(TASKID)
	        ITYP = ITYPE(I)
	        if(LMOVE(ITYP))then 
                  VX(I) = VX(I)*SCVX
                  VY(I) = VY(I)*SCVY
                  VZ(I) = VZ(I)*SCVZ
                  SX(I) = SX(I)*DRCX + VX(I)*DRCV*MASSDI(I)
                  SY(I) = SY(I)*DRCY + VY(I)*DRCV*MASSDI(I)
                  SZ(I) = SZ(I)*DRCZ + VZ(I)*DRCV*MASSDI(I)
                end if
            end do
*   3.4.2.5  Scale simulation box
            call RECLEN(DRCX,DRCY,DRCZ)
         else      !  .not.LNPT
*  3.4.3 Velocities corrections in the NVT ensemble 
           call GETEMP(.true.)
           if(LSCTD)then
             do I=NAB(TASKID),NAE(TASKID)
	       ITYP = ITYPE(I)
               SCV = 1.-SCM(ITYP)*HTSTEP
	       if(LMOVE(ITYP))then
                 VX(I) = VX(I)*SCV
                 VY(I) = VY(I)*SCV
                 VZ(I) = VZ(I)*SCV
               end if
             end do
           else
             SCV = 1.-SC*HTSTEP
             do I=NAB(TASKID),NAE(TASKID)
	       ITYP = ITYPE(I)
	       if(LMOVE(ITYP))then
                 VX(I) = VX(I)*SCV
                 VY(I) = VY(I)*SCV
                 VZ(I) = VZ(I)*SCV
               end if
             end do
           end if   ! if(LSCTD)
         end if     ! if(LNPT)
*   3.4.4  Absorbed energy
         EABS = EABS - (SC+SCL*(1.-3./FNST))*TKE*TSTEP
       end if       ! if(LNVT)
*
*   3.5  Second integration of Nose-Hoover equations (t+1/2->t+1)
*   -------------------------------------------------------------
      else if(IAL.eq.2)then
        if(LNVT)then
*   3.5.1 Correct velocities due to thermostat
          if(LSCTD)then
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
              SCV = 1.-SCM(ITYP)*HTSTEP
	      if(LMOVE(ITYP))then
                VX(I) = VX(I)*SCV
                VY(I) = VY(I)*SCV
                VZ(I) = VZ(I)*SCV
              end if
            end do
          else
            SCV = 1.-SC*HTSTEP
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
	      if(LMOVE(ITYP))then 
                VX(I) = VX(I)*SCV
                VY(I) = VY(I)*SCV
                VZ(I) = VZ(I)*SCV
              end if
            end do
          end if
*   3.5.2 Correct velocities due to barostat
          if(LNPT)then
            SCVX = 1.-(1.-3./FNOP)*SCL*HTSTEP
            SCVY = SCVX
            SCVZ = SCVX
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
	      if(LMOVE(ITYP))then 
                VX(I) = VX(I)*SCVX
                VY(I) = VY(I)*SCVY
                VZ(I) = VZ(I)*SCVZ
              end if
            end do
*   3.5.3  Recalculate temperature
            call GETEMP(.true.) 
*   3.5.4  Correct ita
            if(LSEP)then
              SCLX = SCLX+DQP*(VOL*DPE(2)+2.*TKE/FNST-SCLX*SC)*HTSTEP
              SCLY = SCLY+DQP*(VOL*DPE(3)+2.*TKE/FNST-SCLY*SC)*HTSTEP
              if(LHEX)then
                 SCLX=0.5*(SCLX+SCLY)
                 SCLY=SCLX
              end if
              SCLZ = SCLZ+DQP*(VOL*DPE(4)+2.*TKE/FNST-SCLZ*SC)*HTSTEP
              SCL = SCLX+SCLY+SCLZ   ! trace(P)
            else
              SCL = SCL + DQP*(3.*VOL*DPE(1)+6.*TKE/FNST-SCL*SC)*HTSTEP
              SCLX=SCL
              SCLY=SCL
              SCLZ=SCL
C    Control fluctuations
	      if(abs(SCL*HTSTEP).ge.0.1)then
	        write(*,*)' too strong fluctuations in NPT algorithm'
                if(NUMTASK.gt.1)write(*,*)' at node ',TASKID
	        write(*,*)' scaling factor ',SCL*HTSTEP
   	        if(IPRINT.ge.7)then
	          write(*,*)'------------------->'
	          write(*,*)SCL,TRYCK,TRYCKM
	          write(*,*)VIR1*UNITP/(3.*VOL),VIR2*UNITP/(3.*VOL),
     + PELS1*UNITP/(3.*VOL),PELS2*UNITP/(3.*VOL),VIRB*UNITP/(3.*VOL)
	          write(*,*)'<-------------------'
	        end if
                ISKK=ISKK+1
                if(ISKK.gt.100)then 
	write(*,*)' !!! repeated failure of NPT-algorithm' 
	write(*,*)' !!! restart program with constant volume'
	write(*,*)' !!! or increase thermostat parameter for pressure'
                  call FINAL
                end if
              end if    ! if(LSEP   
            end if
*   3.5.5  Correct ksi due to barostat
            SC = SC + (SCL**2*RTP-DQT/FNST)*HTSTEP
            if(LSCTD)then
              do I=1,NTYPES
                SCM(I) = SCM(I) + (SCL**2*RTP-DQT/FNST)*HTSTEP
              end do
            end if
          else      !  .not.LNPT
*   recalculate temperature
            call GETEMP(.true.)
          end if    ! if(LNPT
*   3.5.6 Correction ksi due to thermostat
          if(LSCTD)then
            SCA = 0.
            do I=1,NTYPES
              DTE(I)   = TEMPR(I)/TRTEMP-1.d0          
              SCM(I) = SCM(I) + DTE(I)*DQT*HTSTEP
              SCA = SCA + SCM(I)*EKIN(I) 
            end do
            SC = SCA/TKE
          else 
            DKE  = TEMP/TRTEMP-1.d0
	    SC = SC + DKE*DQT*HTSTEP
          end if
*   3.5.7  Absorbed energy
         EABS = EABS - (SC+SCL*(1.-3./FNST))*TKE*TSTEP
        end if ! if (LNVT
*
*   3.6  Case of single-step SHAKE algorithm
*   ----------------------------------------
      else if(IAL.eq.3)then
*   3.6.1  Correct ksi(t+1) from velocities and temperature at t+1/2
        if(LNVT)then
	  SC = SC + DKE*DQT*TSTEP      !  Ksi(t+1)
C    This is separate thermostat for each species (may be not work)
          if(LSCTD)then
            SCA = 0.
            do I=1,NTYPES
              SCM(I) = SCM(I) + DTE(I)*DQT*TSTEP     !  DQT -thermostat mass
              SCA = SCA + SCM(I)*EKIN(I) 
            end do
            SC = SCA/TKE
          end if
          if(LNPT)then
*   3.6.2 Barostat corrections (NPT)
*   -------------------------------
*   3.6.2.1 Correction ita  t-1/2  -> t+1/2 
C         TKE is kinetic energy
            if(LSEP)then
              SCLX = SCLX+DQP*(VOL*DPE(2)+2.*TKE/FNST-SCLX*SC)*TSTEP
              SCLY = SCLY+DQP*(VOL*DPE(3)+2.*TKE/FNST-SCLY*SC)*TSTEP
              if(LHEX)then
                 SCLX=0.5*(SCLX+SCLY)
                 SCLY=SCLX
              end if
              SCLZ = SCLZ+DQP*(VOL*DPE(4)+2.*TKE/FNST-SCLZ*SC)*TSTEP
              SCL = SCLX+SCLY+SCLZ   ! trace(P)
            else
              SCL = SCL + DQP*(3.*VOL*DPE(1)+6.*TKE/FNST-SCL*SC)*TSTEP
              SCLX=SCL
              SCLY=SCL
              SCLZ=SCL
*   3.6.2.2    Control fluctuations
	      if(abs(SCL*TSTEP).ge.0.1)then
	        write(*,*)' too strong fluctuations in NPT algorithm'
                if(NUMTASK.gt.1)write(*,*)' at node ',TASKID
	        write(*,*)' scaling factor ',SCL*TSTEP
   	        if(IPRINT.ge.7)then
	          write(*,*)'------------------->'
	          write(*,*)SCL,TRYCK,TRYCKM
	          write(*,*)VIR1*UNITP/(3.*VOL),VIR2*UNITP/(3.*VOL),
     + PELS1*UNITP/(3.*VOL),PELS2*UNITP/(3.*VOL),VIRB*UNITP/(3.*VOL)
	          write(*,*)'<-------------------'
	        end if
                ISKK=ISKK+1
                if(ISKK.gt.100)then 
	write(*,*)' !!! repeated failure of NPT-algorithm' 
	write(*,*)' !!! restart program with constant volume'
	write(*,*)' !!! or increase thermostat parameter for pressure'
                  call FINAL
                end if
              end if    ! if(LSEP   
            end if
*  this is for correction of velocities
            SCV = ((1.-3./FNST)*SCL)*TSTEP
*   3.6.2.3 Calculate scaling coefficients  
            DRCX = TSTEP*SCLX       ! eta*dt
            DRCY = TSTEP*SCLY
            DRCZ = TSTEP*SCLZ
*   3.6.2.4 Correct coordinates and velosities 
*   3.6.2.4.1   Calculate centre of mass and molecular momenta
	    IMB = NNUM(NAB(TASKID))
	    IME = NNUM(NAE(TASKID))
C    Cycle over molecules
	    do IMOL = IMB,IME
              ITYP = ITM(IMOL)
              SUMM        = SUMMAS (ITYP)
	      if(LMOVE(ITYP))then 
		NSBEG = ISADR(ITYP)+1
		NSEND = ISADR(ITYP+1)
		XC       = 0.D0
      	        YC       = 0.D0
                ZC       = 0.D0
	        PXC	= 0.
                PYC	= 0.
                PZC	= 0.
C   Cycle over atoms within the molecule                  
		IBEG=ISADDR(ITYP)+1+NSITS(ITYP)*(IMOL-IADDR(ITYP)-1)
		IEND=ISADDR(ITYP)+NSITS(ITYP)*(IMOL-IADDR(ITYP))
                DO I       = IBEG,IEND
                  IS	     = NSITE(I)
                  XC      = XC+MASS(IS)*SX(I)
                  YC      = YC+MASS(IS)*SY(I)
                  ZC      = ZC+MASS(IS)*SZ(I)
C   Molecular momenta
                  PXC     = PXC+VX(I)
                  PYC     = PYC+VY(I)
                  PZC     = PZC+VZ(I)
                END DO! OF N
                XC       = XC/SUMM    
      	        YC       = YC/SUMM
      	        ZC       = ZC/SUMM
*   3.6.2.4.2  Correction to COM momenta
                DVX	= PXC*SCV
                DVY	= PYC*SCV
                DVZ	= PZC*SCV 
*   3.6.2.4.3  COM displacements
                DXC = XC*DRCX
                DYC = YC*DRCY
                DZC = ZC*DRCZ
*   3.6.2.4.4  Correct atom velocities
C    DVX - molecular momentum correction
C    DVX/SUMM - COM velocity correction
C    MASS(IS)*DVX/SUMM - atom momentum correction
      	        DO I       = IBEG,IEND
      	          IS	= NSITE(I)
                  FAC	= MASS(IS)/SUMM
                  VX(I)	= VX(I)-DVX*FAC
                  VY(I)	= VY(I)-DVY*FAC
                  VZ(I)	= VZ(I)-DVZ*FAC
*   3.6.2.4.5  Correct atom positions
                  SX(I) = SX(I) + DXC 
                  SY(I) = SY(I) + DYC 
                  SZ(I) = SZ(I) + DZC 
                end do
	      end if
            end do ! of ITYP
*   3.6.2.5  Scale simulation box
            DRCX=DRCX+1.
            DRCY=DRCY+1.
            DRCZ=DRCZ+1.
            call RECLEN(DRCX,DRCY,DRCZ)
*   3,6,2,6 Correction ksi due to barostat
C
C   Additional correction to ksi:
C       d(ksi)=ita**2/W-kT
C       ita/W -> SCL
C     so d(ksi)  ->   SCL**2*(W/Q) - kT/Q 
C                       
            SC = SC + (SCL**2*RTP-DQT/FNST)*TSTEP    !  RTP=DQT/DQP
            if(LSCTD)then    
              SCA = 0.
              do I=1,NTYPES
                SCM(I) = SCM(I) + (SCL**2*RTP-DQT/FNST)*TSTEP
                SCA = SCA + SCM(I)*EKIN(I) 
              end do
              SC = SCA/TKE
            end if
          end if     ! if(LNPT)
*   3.6.4  Absorbed energy
          EABS = EABS - (SC+SCL*(1.-3./FNST))*TKE*TSTEP
        end if       ! if(LNVT)
      end if   ! if(IAL
*
*  3.7 Forcible adjusting of velocities to given temperature
*  ---------------------------------------------------------
      IF(LSCLT.and.(IAL.eq.2.or.IAL.eq.3)) THEN 
        if(LSCTD)then
          DO ITYP     = 1,NTYPES
            IF(DABS(TRTEMP-TEMPR(ITYP)).GT.TDELT.and.LMOVE(ITYP))THEN
              SCV = DSQRT(TRTEMP/TEMPR(ITYP))
	      IBEG=ISADDR(ITYP)+1
	      IEND=ISADDR(ITYP+1) 
	      if(IBEG.le.NAB(TASKID))IBEG=NAB(TASKID)
	      if(IEND.gt.NAE(TASKID))IEND=NAE(TASKID)
	      do I=IBEG,IEND
                VX(I)       = SCV*VX(I)
                VY(I)       = SCV*VY(I)
                VZ(I)       = SCV*VZ(I) 
	      end do  
              if(IPRINT.ge.5)write(*,*)
     +      'velocities scaled by ',SCV,' for type',ITYP
              NRTSC(ITYP) =  NRTSC(ITYP)+1
	    end if
	  end do
        else
          IF(DABS(TRTEMP-TEMP).GT.TDELT)THEN
            SCV = dsqrt(TRTEMP/TEMP)
            do I=NAB(TASKID),NAE(TASKID)
	      ITYP = ITYPE(I)
	      if(LMOVE(ITYP))then 
                VX(I) = VX(I)*SCV
                VY(I) = VY(I)*SCV
                VZ(I) = VZ(I)*SCV
              end if
            end do
            if(IPRINT.ge.5)write(*,*)
     +      'velocities scaled by ',SCV
            NRTSC(ITYP) =  NRTSC(ITYP)+1
          end if
        end if
        call GETEMP(.true.)
      END IF!(LSC...
*
      RETURN
      END    
*
*=============== SHEJK ===========================================
*
*   4. SHAKE algorithm for rigid molecules
*   --------------------------------------
C   This algorithm provides integration of equations of motions
C   with constraints: All bond lengths in the given molecule (IMOL) are fixed
C   The algorithm is based on Lagrange formalism of solution of 
C   equations of motion with constrains
C
      SUBROUTINE SHEJK(IMOL)
*
*   4.1  Definitions 
*   ----------------
      include "prcm.h"
      LOGICAL MOVING(NTOT),MOVED (NTOT),DONE
      PARAMETER (MAXIT=2000)
*
*   4.2  Prepare
*   ------------
*   4.2.1 Calculate some references
      ITYP      = ITM(IMOL)              ! type of the molecule
      NRBS      = NRB(ITYP)              ! number ob bonds (constraints)
      if(ISHEJK(ITYP).le.0.or.NRBS.le.0)return
      TOL2      = TOLER**2               ! tolerance parameter
      NBBEG     = IADB(ITYP)             ! first bond in the molecule
      NBEND     = NBBEG+NRBS-1           ! last bond
      ISHF      = ISADDR(ITYP)+NSITS(ITYP)*(IMOL-IADDR(ITYP)-1)
      ISBEG     = ISHF+1                 ! first atom in the molecule
      ISEND     = ISHF+NSITS(ITYP)       ! last atom
*   4.2.2 Remember unconstrained positions
C   Arrays for forces FX,... are used as a temporary storage
C   Forces are not needed at this step
      do I=ISBEG,ISEND
        FX(I)=SX(I)
        FY(I)=SY(I)
        FZ(I)=SZ(I)
      end do
*   4.2.3 Set up logical arrays controlling constraints
      DO M      = NBBEG,NBEND
	LI	= IB(M)+ ISHF
	LJ	= JB(M)+ ISHF
        MOVING(LI) = .FALSE.
        MOVED (LI) = .TRUE.
        MOVING(LJ) = .FALSE.
        MOVED (LJ) = .TRUE.
      END DO! OF M
*
*  4.2.4  Begin iterations:
*  -----------------   
      ITER      = 0
      DONE      = .FALSE.
C  Point to returm:
    1 CONTINUE
C  Maximum iterations exceed
      IF(ITER.GT.MAXIT) GO TO 1000
*
*  4.3  Adjustment of constraints 
*  ------------------------------
*  4.3.1 Begin cycle over the bonds
      IF(.NOT.DONE) THEN
C   "DONE" will be .true. at the end of iteration if all constrains satisfied
        DONE      = .TRUE.
        DO M      = NBBEG,NBEND
          I       = IB(M)+ISHF
          J       = JB(M)+ISHF
C  MOVED=.true. means that the atom should be checked:
C  either if the bond was moved after the constrain was adjusted 
C  or the constrain was not adjusted at all      
          if(MOVED(I).OR.MOVED(J)) then
*  4.3.2 Calculate deviation
            AMAS      = 1.D0/MASSDI(I)
            BMAS      = 1.D0/MASSDI(J)
            BNDSQ     = RB(M)**2
            BXN       = SX(I)-SX(J)
            BYN       = SY(I)-SY(J)
            BZN       = SZ(I)-SZ(J)
	    call PBC(BXN,BYN,BZN)
            BSQ       = BXN**2+BYN**2+BZN**2
            DIFFSQ    = BNDSQ-BSQ
            if(IPRINT.ge.10)write(*,'(2i5,f9.6,4(2x,l1),i5)')
     +I,J,sqrt(BSQ),moved(I),moved(j),moving(I),moving(j),ITER
*   4.3.3 Deviation is still large - adjust the atoms
            IF(DABS(DIFFSQ).GT.(BNDSQ*TOL2))THEN
C   OX.. Constrained positions from step n
              DXN       = OX(I)-OX(J)
              DYN       = OY(I)-OY(J)
              DZN       = OZ(I)-OZ(J)
	      call PBC(DXN,DYN,DZN)
              VDOT      = BXN*DXN+BYN*DYN+BZN*DZN
              IF(VDOT.LT.BNDSQ*BDTOL) THEN
C  This happens if the deviation increased....
      write(*,*)' --> constrained failure for atoms ',I,J
      WRITE(*,66) 'old     :',OX(I),OY(I),OZ(I),OX(J),OY(J),OZ(J)
      WRITE(*,66) 'unconstr:',FX(I),FY(I),FZ(I),FX(J),FY(J),FZ(J)
      WRITE(*,66) 'current :',SX(I),SY(I),SZ(I),SX(J),SY(J),SZ(J)
      write(*,*)'OLD DISTANCE: ',DSQRT((OX(I)-OX(J))**2
     X                             +(OY(I)-OY(J))**2
     X                             +(OZ(I)-OZ(J))**2)
      write(*,*)'unconstrained:',DSQRT((FX(I)-FX(J))**2
     X                             +(FY(I)-FY(J))**2
     X                             +(FZ(I)-FZ(J))**2)
      write(*,*)'NEW DISTANCE: ',DSQRT((SX(I)-SX(J))**2
     X                             +(SY(I)-SY(J))**2
     X                             +(SZ(I)-SZ(J))**2)
*
   66 FORMAT(1X,a,3(1X,F10.6),3X,3(1X,F10.6))
CD    STOP 'CONSTRAINED FAILURE'
              END IF
              RMA       = MASSDI(I)
              RMB       = MASSDI(J)
              GAB       = DIFFSQ/(2.D0*(RMA+RMB)*VDOT)
              DXN       = DXN*GAB
              DYN       = DYN*GAB
              DZN       = DZN*GAB
C    Moving the atoms
              SX(I)     = SX(I)+RMA*DXN
              SY(I)     = SY(I)+RMA*DYN
              SZ(I)     = SZ(I)+RMA*DZN
              SX(J)     = SX(J)-RMB*DXN
              SY(J)     = SY(J)-RMB*DYN
              SZ(J)     = SZ(J)-RMB*DZN
C    Mark atom as MOVING during the iteration
              MOVING(I) = .TRUE.
              MOVING(J) = .TRUE.
C    Next iteration required because some atoms were moved
              DONE      = .FALSE.
            END IF
          END IF
        END DO! OF M - all bonds
*
*  4.4  Check that all constraints are satisfied
*  ---------------------------------------------
*  Obs!    These are separate loops ( Reccurence possible )
C  Set OK to all the bonds  
        DO M      = NBBEG,NBEND
          I       = IB(M)+ISHF
	  J	  = JB(M)+ISHF
	  MOVED(I)  = .false.
          MOVED(J)  = .false.
        END DO! OF M
C  If the atom was moved during the last iteration, mark it as requiring check
C  at the next iteration
        DO M      = NBBEG,NBEND
          I       = IB(M)+ISHF
	  J	  = JB(M)+ISHF
          if(MOVING(I))MOVED(I)=.true.
          if(MOVING(J))MOVED(J)=.true.
C  ... and set MOVING as false for the next iteratio
          MOVING(I) = .FALSE.
          MOVING(J) = .FALSE.
        END DO! OF M
        ITER      = ITER + 1
C  Continue iterations
        GO TO 1
      END IF  !(.not.DONE
*
*   4.5 Calculate corrections to the velocities and contribution to virial
*   ----------------------------------------------------------------------
 1000 CONTINUE
*
      IF(ITER.GE.MAXIT) THEN
        write(*,*)'??? NO CONVERGENCE AFTER MAXIT INTERATIONS'
        write(*,'(a,i6,a,i6,a,f12.6,a,f12.6/)')
     +'??? BOND BETWEEN ATOMS ',I,' AND',J,'  Reqv = ',
     +dsqrt(BNDSQ),' Rcalc = ',dsqrt(BSQ)
      END IF
C   velocity correction
      DO I      = ISBEG,ISEND
        DIX       = (SX(I)-FX(I))/(MASSDI(I)*TSTEP)   
        DIY       = (SY(I)-FY(I))/(MASSDI(I)*TSTEP)
        DIZ       = (SZ(I)-FZ(I))/(MASSDI(I)*TSTEP)
        VIRAN(1)=VIRAN(1)+DIX*OX(I)/TSTEP
        VIRAN(2)=VIRAN(2)+DIY*OY(I)/TSTEP
        VIRAN(3)=VIRAN(3)+DIZ*OZ(I)/TSTEP
        VX(I)     = VX(I)+DIX
        VY(I)     = VY(I)+DIY
        VZ(I)     = VZ(I)+DIZ
      END DO! OF N
      RETURN
      END
