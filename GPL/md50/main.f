      PROGRAM MD
C  MDynamix v.5.0
C
C  This is the main program unit taking care of all what follows. 
C  
C   This program unit is organized as follows:
C 
C   1. Definitions
C   2. Introduction
C      2.1 Set version information and data
C      2.2 Parallel initialisation
C      2.3 Some control parameters 
C   3. Program input 
C      3.1 Read in input data 
C      3.2 Define file names
C      3.3 Setting initial values for some variables / arrays
C      3.4 Input of molecule structures and force field from data base
C   4. Calculate some total values
C      4.1 Set up initial values to zero and organise cycle over molecule types
C      4.2 General total values
C      4.3 Calculate the number of degrees of freedom:
C      4.4 Report total values
C      4.5 Check array boundaries 
C      4.6 Zero average accumulators
C      4.7 miscellaneous
C   5. Setup internal units. Convert all data to internal units
C   6. Set up parameters for NVT and NPT ensembles:
C      6.1  NVT
C      6.2  NPT
C   7. Set up coordinates and velocities of atoms.    Choice:
C      7.1 Read from restart file (continue old run)   
C      7.2 Prepare a new run:
C         7.2.1  Initial zero-s
C         7.2.2  Set up initial center-of-mass coordinates of molecules
C                Choice:
C            7.2.2.1  FCC lattice 
C            7.2.2.2  Molecule(s) in a cylindrical hole with solvent around
C            7.2.2.3  Molecule(s) in a spherical hole with solvent around
C            7.2.2.4  Cubic lattice
C            7.2.2.5  read atom coord. from input file
C            7.2.2.6  read molecular center-of-mass coord. from input file
C            7.2.2.7  read atom coord. from input file in XMOL format
C            7.2.2.8  Non-implemented way of setting up coordinates 
C         7.2.3  Set up orientations of the molecules and atom coordinates
C         7.2.4  Initialize velocities
C            7.2.4.1  Maxwell distribution of velocities
C            7.2.4.2  Calculate temperature and scale velosities
C            7.2.4.3  Recalculate temperature after velocity scaling
C   8. Calculate parameters for electrostatics.     
C      8.1 Define minimum box size
C      8.2 Reaction field method
C      8.3 Ewald method 
C      8.4 Calculate long-range corrections from Lennard-Jones interactions 
C   9. Report basic simulation parameters
C   10. Finalize preparation
C      10.1  Prepare RDF calculations
C      10.2  Prepare TCF calculations
C      10.3  Prepare list of bound atoms
C      10.4  Change temperature and/or density if nesessary
C      10.5  Bind some atoms to specified locations (if LRR=.true.)
C      10.6  Check configuration on non-overlapping
c      10.7  Verbose: Dump start configuration
C   11. Case if we do not really wanted to run MD, but only (Choice): 
C      11.1  get results of previously made run
C      11.2  check that input OK
C   12. Run MD.           Choice:  
C      12.1  SHAKE algorithm
C      12.2  Double time step algorithm
C   13. Report final results
C      13.1  Report calculated averages
C      13.2  Report RDF-s 
C      13.3  Report TCF-s 
C      13.4  Dump final configuration in "Xmol" format   
C   14. Finalise (close up MPI connections)
* 
*-----------------------------------------------
*
*   1. Definitions
*   -------------- 
C  prcm.h is the include file containing static dimensioning
C  of the working fields. Practically all modules depend somehow on it, 
C  so you must recompile the sources each time you change this file.
      include "prcm.h"
C   local arrays:
C   QP is used for quaternions (to set molecule's orientations)  
      DIMENSION  QP(4)
      character*128 AUX,TAKESTR
      character*20 FNM, DATE*8,FINP*24, CDUB*4, CD1*1
      logical LMA
*   2. Introduction
*   ---------------
*  2.1. Set version information and data
      label = "MDYN_v5.0"
      DATE  = '02.01.07'
*  2.2 Parallel initialisation
*  -----------------------
C  Paralel version requires initialisation of MPI system. |NUMTASK| 
C  reports the number of processes the MD program is running. The 
C  identification number of this instance of process is in variable
C  |taskid|. One of the processes is selected as the master. 
C  By definition here, his |taskid| is |NUMTASK-1|. This process
C  takes care about user oriented I/O.
	call PARAINI(numtask,taskid) 
	MAST=NUMTASK-1
	if(TASKID.eq.MAST)then
          write(*,*)'*** MOLECULAR DYNAMICS',
     X' PROGRAM FOR MIXTURES OF RIGID or FLEXIBLE MOLECULES'
	  write(*,*)
	  if(NUMTASK.gt.1)then
	    write(*,*)' Parallel version for connected machines '
            write(*,*)' Running on ',NUMTASK,' processors'
            write(*,*)   
	  end if
	  write(*,*)'Version ',label,'           from ',DATE
          write(*,*)
        end if
*
	time0	= cputime(0.d0)
*
*  2.3 Some control parameters 
*  ---------------------------
      BDTOL      = 0.0001   ! for constrained dynamics    
      MSTEP      = 0        ! Num. of MD steps in this run 
      DO I       = 1,NTPS
        NRA(I)   = 0
        NRB(I)   = 0
        NRT(I)   = 0
        NRI(I)   = 0
      END DO  
      LFORMI = .false.     ! control format of the restart file
      LFORMO = .false.     
*   
*   3. Program input 
*   ----------------
*   3.1. Read in input data 
*   -----------------------
      call INPUT                  ! see file readmol.f
      IPRINN=IPRINT               ! auxilary input control parameter
      if(TASKID.ne.MAST)IPRINT=0  ! set output level to 0 on non-master nodes 
C  This will set averages accumulators to 0 for a new run
      if(.not.LRST) L0AVS=.true.  
C  Check arrays boundaries
      IF(NTYPES.GT.NTPS) STOP '!!! increase NTPS in dimpar.h!'
      if(MAX.gt.MAXS)stop '!!! increase MAXS in dimpar.h!'
*
*   3.2 Define file names
*   ---------------------  
C  This is how fortran gets the full names of 
C  restart (dump) file (*.dmp)
C  radial distribution function file (*.rdf), 
C  time correlation function file (*.tcf),
C  cartesian output coordinate file (*.car),
C  zero-th trajectory file (*.000),
C  initial coordinates *.start file by extension of the base filename |FNAME|.
      IC=LENG(FNAME,20)
      IC1=IC+1
      IC2=IC+4
      ICFN=IC+2
      FNM=  '                    '
      FDUMP='                        '
      FRDF=FDUMP
      FTCF=FDUMP
      FTRK=FDUMP
      FMOL=FDUMP
      FINP=FDUMP
      FNM(1:IC)=FNAME(1:IC)
      FDUMP(1:IC)=FNM
      FDUMP(IC1:IC2)='.dmp'         !
      FRDF(1:IC)=FNM
      FRDF(IC1:IC2)='.rdf'
      FTCF(1:IC)=FNM
      FTCF(IC1:IC2)='.tcf'
      FMOL(1:IC)=FNM
      FMOL(IC1:IC2+1)='.xmol'
      FTRK(1:IC)=FNM
      FTRK(IC1:IC2)='.000'  
      FINP(1:IC)=FNM
      FINP(IC1:IC2+2)='.start'  
      if(IPRINT.ge.8)then
	write(*,*)' file names:'
        write(*,*)FDUMP
	write(*,*)FRDF
	write(*,*)FMOL
	write(*,*)FTRK
      end if
*
*  3.3 Setting initial values for some variables / arrays
*  -------------------------------------------------
      NTREK	= TTREK/DT+0.5   ! num. of configurations in a trajectory file
      if(NTREK.le.0)NTREK=1
C  remember size/density from input 
C  ( they will be owerwritten by the restart file)
      BOXLN=BOXL
      BOYLN=BOYL
      if(LHEX)BOYLN=BOXLN*cos30
      BOZLN=BOZL
      RHON=RHO  
      RRLL=RLR
      NRBB         = 0
      NRAA         = 0
      NRTT         = 0
      IADB(1)      = 1
      IADA(1)      = 1
      IADT(1)      = 1
      TOTMAS       = 0.
      NNAD=0
C  Input error flag
      IE	= 0
*
*  3.4 Input of molecule structures and force field from data base
*  ---------------------------------------------------------------   
	if(TASKID.eq.MAST)write(*,*)
     +'--------------------------------------------------------------'
	do ITYP=1,NTYPES
          IF(NSPEC(ITYP).LE.0)then
             write(*,*)'!!! Zero molecules of type ',ITYP
             write(*,*)'    Check the input file'
             stop
          end if
          CALL READMOL (ITYP)
          TOTMAS       = TOTMAS+SUMMAS(ITYP)*NSPEC(ITYP)
	  NAME(ITYP)=NAMOL(ITYP)(1:6)
          if(IPOT(ITYP).eq.2.and.TASKID.eq.MAST)write(*,*)
     +' Internal motion: flexible SPC water'          
        end do
      if(IPRINT.ge.5)write(*,*)
     +'--------------------------------------------------------------'
      IF(LSHEJK) THEN
C  Case of constrained dynamics. Simple verlet algorithm in this case.
        NFREQ = 1 
	if(TASKID.eq.MAST)then
           write(*,*)
          write(*,*) '*** CONSTRAINED DYNAMICS WILL BE USED! ***'
          do I=1,NTYPES
            if(ISHEJK(I).le.0)write(*,*)'  exception: ',NAME(I)
          end do
        end if
      END IF 
*
*  4. Calculate some total values
*  ------------------------------
*  
C  The composition of the system is described in terms of
C atoms, molecules, sites and types. 
C "Atoms" run over all the atoms in the system
C "molecules" run over all molecules
C "Types" are different molecule types
C "Sites" are different atoms on a molecule. They run over each molecule type
C and, within this type, over all atoms on a molecule of this type. 
C
C  The molecular internal structure is described in term of "bonds",  
C  "angles" and "dihedral angles" (torsions and impropers) in according
C  with AMBER or CHARMM force field formalism
C
C Total values:
C	NSTOT		total atoms
C	NOP		total molecules
C	NSITES	 	total sites
C       NRBND           total bonds
C       NRANG           total angles
C       NRTOR           total dehedral angles
C
*  4.1 Set up initial values to zero and organise cycle over molecule types
*  ------------------------------------------------------------------------
      NOP     = 0
      NSITES  = 0
      NSTOT   = 0
      NRBND   = 0
      NRANG   = 0
      NRTOR   = 0
      FNST    = 0.            
      FTTR    = 0.
      FNSTR   = 0.
      LMA=.true.
      if(IPRINT.ge.10)then
	  write(*,*)
	  write(*,*)' Check references:'
	  write(*,*)' Iatom    Imol    Isit   Ityp     Name'
      end if  
      DO ITYP = 1,NTYPES
*  4.2 General total values
*  ------------------------
        NOP     = NOP    + NSPEC(ITYP)       ! Num. of molecules
        NSITES  = NSITES + NSITS(ITYP)       ! Num. of sites
        NSTOT   = NSTOT  + NSPEC(ITYP)*NSITS(ITYP)  ! Num. of atoms
	NRBND=NRBND+NRB(ITYP)*NSPEC(ITYP)    ! Num. of bonds
	NRANG=NRANG+NRA(ITYP)*NSPEC(ITYP)    ! Num. of angles
	NRTOR=NRTOR+NRT(ITYP)*NSPEC(ITYP)    ! Num. of dihedrals
*  4.3 Calculate the number of degrees of freedom:
*  -----------------------------------------------
C  FNST - total
C  FTTR - translational
C  FNSTR - rotational
	if(LMOVE(ITYP))then                               
	  FNST  = FNST+3*NSPEC(ITYP)*NSITS(ITYP)   
	  FTTR  = FTTR+3*NSPEC(ITYP)               
	  if(LSHEJK.and.ISHEJK(ITYP).gt.0)
     +      FNST=FNST-NRB(ITYP)*NSPEC(ITYP)
          if(NSITS(ITYP).eq.2)FNSTR=FNSTR+2.*NSPEC(ITYP)
          if(NSITS(ITYP).ge.3)FNSTR=FNSTR+3.*NSPEC(ITYP)
	else
C  if there are no fixed molecules ( LMOVE=.t. for all types), subtract 3
C  from the number degrees of freedom because of translational invariancy
C  ( will be done below if LMA=.true.)
	  LMA=.false.
	end if
      end do
*
*   4.4 Report total values
*   -----------------------
      if(TASKID.eq.MAST)then
        PRINT "(/'*** NR OF SITES/MOL (NSITES)  ',I6)" ,NSITES
        PRINT "( '*** NR OF SITES/TOT (NSTOT)   ',I6)" ,NSTOT
        PRINT "( '*** NR OF PARTICLES (NOP)     ',I6)",NOP
        PRINT "( '*** NR OF BONDS               ',I6)",NRBND
	PRINT "( '*** NR OF ANGLES              ',I6)",NRANG
	PRINT "( '*** NR OF TORSIONS            ',I6)",NRTOR
      end if
*  FNSTI - number of internal degrees of freedom 
      FNSTI	= FNST-FTTR-FNSTR
      if(LMA.and.NOP.gt.1)then
         FNST      = FNST-3.
         FTTR      = FTTR-3.
      end if
      if(TASKID.eq.MAST)then
        write(*,'(a,f8.0)')'*** NR of DEGREES of FREEDOM ',FNST
        write(*,'(a,f8.0)')'***        translationalnal: ',FTTR
        write(*,'(a,f8.0)')'***        rotational:       ',FNSTR
        write(*,'(a,f8.0)')'***        internal:         ',FNSTI
	write(*,*)
      end if
*        
*    4.5 Check array boundaries  
*    --------------------------
      IF(NRBB.GT.NB  ) STOP '!!! INCREASE NB in dimpar.h  !'
      IF(NRAA.GT.NA  ) STOP '!!! INCREASE NA in dimpar.h  !'
      IF(NRTT.GT.NT  ) STOP '!!! INCREASE NT in dimpar.h  !'
      IF(NOP   .GT.NPART) STOP '!!! INCREASE NPART in dimpar.h !'
      IF(NSITES.GT.NS   ) STOP '!!! INCREASE NS in dimpar.h   !'
      IF(NSTOT .GT.NTOT ) STOP '!!! INCREASE NTOT in dimpar.h !'
      IF(TOTMAS.EQ.0 ) STOP '!!! NO PARTICLES WITH MASSES DEFINED !'
*
*  4.6 Zero average accumulators
*  -----------------------------   
      if(L0AVS)then
        NHIST = 0              ! set series counter to 0
	NAVT	= 0            ! total averaged configurations
	if(TASKID.eq.MAST)write(*,*)' New averaging'
        NSTEPA	= 0            ! num. of MD steps with averaging
        TIMA	= 0.           ! time of averaging
C  these are average and av**2 accumulators 
        DO I        = 1,NRQS
          AV(I)       = 0.D0
          AW(I)       = 0.D0
        END DO! OF I
C   These are for temporary averages  
	call ZEROAV             ! see service.f
      end if
*
*   4.7 miscellaneous
*   -----------------
*      IF(NOP.EQ.1) RHO = 0.
      FTIMES      = 1./DFLOAT(NFREQ) 
      FNOP        = DFLOAT(NOP)
*
*  5. Setup internal units. Convert all data to internal units
*  -------------------------------------------------------- 
      CALL IUNITS               !  file units.f
*
*  6. Set up parameters for NVT and NPT ensembles:
*  -----------------------------------------------
*  6.1  Temperature controll
*  -------------------------
      IF(LNVT) THEN
        DQT	= (UNITT*1.d15/QT)**2
	if(TASKID.eq.MAST)then
	  write(*,*)
          write(*,'(a,f10.3)')'*** Constant temperature MD: T = ',TRTEMP
          write(*,*)'   NOSE-NVT RELAXATION TIME      ',QT,' fs'
	end if
      END IF!(LNVT...
      if(LSCLT.and.TASKID.eq.MAST)then
        write(*,*)' Temperature will be adjusted to ', TRTEMP,
     + ' if its deviation exceeds ',TDELT
      end if
      if(LSCTD.and.TASKID.eq.MAST)write(*,*)
     +' Separate temperature controll for each',
     +' molecular species is enforced '
*
*   6.2 NPT
*   -------
      IF(LNPT) THEN
        if(.not.LNVT.and.TASKID.eq.MAST)
     +  write(*,*)'!!! Attention: NPE-ensemble???'
        DQT	= (UNITT*1.d15/QT)**2
        DQP	= 3.*UNITE*UNITT**2/(UNITP*(FNST+3)*BOLTZ*TRTEMP*
     +            (QPR*1.d-15)**2)
        RTP     = (QPR/QT)**2*(FNST+3.)/FNST
        if(TASKID.eq.MAST)write(*,'(a,f10.3,a)')
     +  '*** Constant-pressure MD: P = ',TRPRES,' atm'
        if(TASKID.eq.MAST)
     +  write(*,*)'   NOSE-NPT RELAXATION TIME     ',QPR,' fs'
      END IF!(LNVT...                   
      FACVT = FNST*TRTEMP*(QPR*1.d-15/UNITT)**2/3.
      if(LSEP)then
        FACVT=FACVT/3.
        if(TASKID.eq.MAST)write(*,*)
     +' Separate pressure controll in different directions'
      end if
*  6.3  Other notes
      if(LRR.and.TASKID.eq.MAST)then
        write(*,*)
        write(*,*)
     +'!!! Forcible keeping of the structure defined in ',filref,
     +'  DR =',RRLL
      end if
      if(LKONG.and.TASKID.eq.MAST)write(*,*)
     +'   Kong combinatioin rules are used for cross-interactions'
      if(LGEOM.and.TASKID.eq.MAST)write(*,*)
     +'   Geometrical average for Sigma is used for cross-interactions'
*
*  7. Set up coordinates and velocities of atoms
*  ---------------------------------------------
      IF(LRST) THEN
*  7.1 Read from restart file (continue old run)
*  ---------------------------------------------
        CALL RESTRT(1)                  ! resfart.f
C  setup lengths                        
        call RECLEN(ODIN,ODIN,ODIN)     ! setup.f
C  calculate centre-of-masses
        CALL GETCOM                     ! service.f
        if(TASKID.eq.MAST)write(*,'(a,f12.4,a2)')
     +  ' The system was simulated ',TIM*1.d12,'ps'
        if(TASKID.eq.MAST)write(*,'(a,f12.4,a2)')
     +  ' Averages was collected   ',TIMA*1.d12,'ps'
*
      ELSE
*
* 7.2. Prepare a new run:
* ----------------------
* 7.2.1. Initial zero-s
* ---------------------
        NHIST	= 0
        TIM	= 0.
        timest	= 0.
	NTRKF=0           ! num. of trace file
	NTRKZ=0           ! num. of points in the trace file
* 7.2.2. Set up initial center-of-mass coordinates of molecules
* -------------------------------------------------------------
C   Different ways of distributing the molecules in the cell  
*  7.2.2.1. FCC lattice 
        IF(ICHECK.eq.1) THEN 
            CALL FCC(X,Y,Z,BOXL,BOYL,BOZL,NOP,ICELL)    ! setup.f
*  7.2.2.2. Molecule(s) in a cylindrical hole with solvent around
          ELSE if(ICHECK.eq.2)then
            call DNASOL                     ! setup.f
*  7.2.2.3. Molecule(s) in a spherical hole with solvent around
          ELSE if(ICHECK.eq.3)then
            call SPHOLE                     ! setup.f
*  7.2.2.4. Cubic lattice
	  else if(ICHECK.eq.4)then
            CALL CUBIC(X,Y,Z,BOXL,BOYL,BOZL,NOP,ICELL,IPRINT)   ! setup.f
	  else if(ICHECK.eq.-2)then
*  7.2.2.5  read atom coord. from input file
	    open(unit=13,file=FINP,status='old',err=298)
C    reading from input file
            do I=1,NSTOT                                
	      AUX=TAKESTR(13,IE)
	      read(AUX,*,err=297,end=297)SX(I),SY(I),SZ(I) 
	    end do          
	    if(TASKID.eq.MAST)write(*,*)
     +'read initial coordinetes from file ',FINP
	    close(13)  
*  7.2.2.6  read molecular center-of-mass coord. from input file
	  else if(ICHECK.eq.-1)then
	      open(unit=13,file=FINP,status='old',err=398)
	      go to 399         
C   something goes wrong:
C   can not read input file...         
 397	      if(TASKID.eq.MAST)then
              write(*,*)'too few molecules in the input file'
	        IE=99
	        AUX=TAKESTR(13,IE)
	        write(*,*)' 1.e1000 thanks to your system administrator'
	      end if
	      stop
C   can not open input file                    
 398	      if(TASKID.eq.MAST)
     +write(*,*)' input file ',FINP,' not found (ICHECK=-1)'
	      stop
C    reading...
 399          do I=1,NOP                                
		AUX=TAKESTR(13,IE)
	        read(AUX,*,err=397,end=397)X(I),Y(I),Z(I)
		call PBC(X(I),Y(I),Z(I))
	      end do          
	      if(TASKID.eq.MAST)write(*,*)
     +'read initial COM from file ',FINP
	      close(13)  
	  else if(ICHECK.eq.0)then
*    reading input coordinates from XMOL file
	      open(unit=13,file=FINP,status='old',err=298)
              read(13,*)NAINI
              if(NAINI.lt.NSTOT)then
                if(TASKID.eq.MAST)write(*,*)
     &' Input coordinate file ',FINP,' contains too few atoms ',
     & NAINI,' <',NSTOT
                stop
              else if(NAINI.gt.NSTOT)then
                if(TASKID.eq.MAST)write(*,*)
     &'Warning: Input coordinate file ',FINP,
     &' contains more atoms than required '
              end if
              read(13,'(a128)')AUX
              do I=1,80
                if(AUX(I:I+2).eq.'BOX')then
                  IFOUND=1
                  read(AUX(I+4:128),*,end=437,err=437)BOXL,BOYL,BOZL
                  if(IPRINT.ge.5.and.TASKID.eq.MAST)
     &              write(*,'(a,3f10.2)')
     & 'Box sizes are taken from initial structure:',BOXL,BOYL,BOZL
                  go to 438
                end if
              end do
 437          IFOUND=0
              if(IPRINT.ge.5.and.TASKID.eq.MAST)write(*,*)
     &' Box sizes not found in the initial structure'
 438          continue
*    reading coordinates
              do I=1,NSTOT
                read(13,*,err=296,end=297)CD1,SX(I),SY(I),SZ(I)
              end do
              if(TASKID.eq.MAST)write(*,*)
     & ' Coordinates are read from file ',FINP
	  else
*  7.2.2.8. Non-implemented way of setting up coordinates 
	      write(*,*)' available values of ICHECK from -2 to 4'
	      write(*,*)' now ICHECK =',ICHECK
	      stop
        END IF!(ICHE...
        go to 299
*  7.2.2.9
C   something goes wrong:
C   can not read input file...         
 296        if(TASKID.eq.MAST)then
              write(*,*)'read error in the input file'
	        IE=99
	        AUX=TAKESTR(13,IE)
	        write(*,*)' fatal processor error'
	      end if
	      stop       
C   too short input file
 297	      if(TASKID.eq.MAST)then
              write(*,*)'too few atoms in the input file'
	        IE=99
	        AUX=TAKESTR(13,IE)
	        write(*,*)' fatal processor error'
	      end if
	      stop       
C   can not open input file                    
 298	      if(TASKID.eq.MAST)
     +write(*,*)' input file ',FINP,' not found (ICHECK=0)'
	      stop
 299    continue
*
* 7.2.3. Set up orientations of the molecules and atom coordinates
* ----------------------------------------------------------------
C   (exclude the case if we already have atom coordinates)
	if(ICHECK.ne.0.and.ICHECK.ne.-2)then  
          N              = 0
          I              = 0 
C   cycle over types 
          DO ITYP        = 1,NTYPES
            IBES           = ISADR (ITYP)+1
            IENS           = ISADR (ITYP +1)
C   cycle over molecules of this type 
            DO J           = 1,NSPEC(ITYP)
              I              = I+1
C   generate random quaternion
              CALL RANQ(QP)            !  supl.f
C   rotate randomly the molecule and setup atom coordinates
              DO IS          = IBES,IENS
                N              = N+1
                RXI            = R(IS,1)
                RYI            = R(IS,2)
                RZI            = R(IS,3)
C   do not rotate, if IINIT(ITYP)=1 
                if(IINIT(ITYP).ne.1)CALL ROTATE  (QP,RXI,RYI,RZI)
                SX(N)          = RXI+X(I)
                SY(N)          = RYI+Y(I)
                SZ(N)          = RZI+Z(I)
              END DO! OD IS
            END DO! OF I
          END DO! OF ITYP
          if(TASKID.eq.MAST)
     +    PRINT "('*** MOLECULAR ORIENTATIONS INITIALIZED')"
	end if
C    recalculate centre-of-mass
	call GETCOM
*
*  7.2.4. Initialize velocities
*  ----------------------------
        IF(L0VS) THEN    
*  7.2.4.1  Set zero velocitries
          DO N  = 1,NSTOT
            VX(N) = 0.D0
            VY(N) = 0.D0
            VZ(N) = 0.D0
          END DO! OF N
        else

*  7.2.4.2.  Maxwell distribution of velocities
          CALL COMVEL                      !  service.f 
*
          VSQ         = 0.0
          DO N        = 1,NSTOT  
	    ITYP=ITYPE(N)
            if(LMOVE(ITYP))VSQ = VSQ+
     +     (VX(N)**2+VY(N)**2+VZ(N)**2)*MASSDI(N)
          END DO! OF N
          TEMP        = 1.50*VSQ*ENERF*TEMPF/FNST
C  Rescale velocities on factor:
          FTMP        = SQRT(TRTEMP/TEMP)
          DO     N    = 1,NSTOT
	    ITYP=ITYPE(N)
	    if(LMOVE(ITYP))then
              VX(N)       = VX(N)*FTMP
              VY(N)       = VY(N)*FTMP
              VZ(N)       = VZ(N)*FTMP
C  Case if the molecule is fixed (LMOVE=.f.) 
	    else
              VX(N)       = 0.
              VY(N)       = 0.
              VZ(N)       = 0. 
	    end if
          END DO! OF N
        end if
*  7.2.4.3.  Recalculate temperature after velocity scaling and report
        VSQ         = 0.0
        DO N        = 1,NSTOT
          VSQ         = VSQ+(VX(N)**2+VY(N)**2+VZ(N)**2)*MASSDI(N)
        END DO! OF N
        TEMP        = 1.5*VSQ*ENERF*TEMPF/FNST
C  Report       
	if(TASKID.eq.MAST)then
          write(*,'(a)')'*** VELOCITIES INITIATED'
          write(*,*)'   INITIAL TEMPERATURE CALCULATED',TEMP
          write(*,*)
	end if
*  7.2.5  Initial thermostat/barostat parameters
        SC	 = 0.      ! d(ksi)/dt in Nose-Hoover thermostat
        SCL	= 0.       ! ita and its projections
        SCLX    = 0.
        SCLY    = 0.
        SCLZ    = 0.
        do I=1,NTYPES
          SCM(I)=0.        ! d/ksi)/dt for each mol. species
        end do
C  Number of previously steps made
        NSTTOT	= 0      
      END IF!(LRST...
*
* 8. Calculate parameters for electrostatics      
* ------------------------------------------
* 8.1 Define some parameters 
C  minimal box length
      CL=BOXL
      if(CL.gt.BOYL)CL=BOYL
      if(CL.gt.BOZL)CL=BOZL
*
* 8.2 Reaction field method
* -------------------------
C  If ALPHA from input data is < 0, then reaction field method is applied
      if(ALPHA.le.0.d0)then 
C  In this case -ALPHA defines dielectric constant DIEL 
	 DIEL=-ALPHA
C  This parameter sets up special procedure for choose of neighbours
C  See details in subr. CHCNB, file l-forces.f
	 LMNBL=.true.  
C  Cut-off distance must be < L/2 (spherical cut-off essential)
	 if(RCUTT.gt.0.5*CL)then
	    RCUTT=0.5*CL
	    RSSQ= RCUTT**2
	    if(TASKID.eq.MAST)write(*,*)' Cut off distance set to ',RCUTT
	 end if
C   Reaction field method add some terms to energy (RFF) and force (RFF2)
C   If standard reaction field method, RFF and RFF2 will be calculated by:
C        RFF = (1/DIEL-1)/RCUTT
C	 RFF2=2.*(DIEL-1.d0)/(2.d0*DIEL+1)/RCUTT**3
C   Generalised reaction field method take into account ion concentration:
C   see Tironi et al, J.Chem.Phys. v.102,p.5451 (1995)) 
C   Expressions depends on ion concentration through Debue radius
C   and are more complicated
C   FEXP define Debue length in Å.
	 RDEB=FEXP            ! now it is Debue length in Å    
	 if(RDEB.gt.1.d-3)then
C   For convenience, the case of absence salt ( formally RDEB=inf.)
C   can be specified also by FEXP=FDEB=0.
C   RCK is ratio of cutoff to Debue radius
	    RCK=RCUTT/RDEB
	 else
	    RCK=0.
	 end if
C   In the case RCK=0, RFF is defined as in standard reaction field
C   method (see abive)
	 RFF2=-(2.*(DIEL-1.)*(RCK+1.)+DIEL*RCK**2)/
     &      (((2.*DIEL+1.)*(RCK+1.)+DIEL*RCK**2)*RCUTT**3)  
         RFF = (1./(DIEL*(1.+RCK))-1.)/RCUTT
C   This is case of absence of any accounting of long-range electrostatics
	 if(DIEL.le.1.0000001d0)then
            RFF=0.
            RFF2=0.
	    if(TASKID.eq.MAST)then
	      write(*,*)
	      write(*,*)
     +' no electrostatics out Rcutoff'
	      write(*,*) 
	   end if
         else
	  if(TASKID.eq.MAST)then
	    write(*,*)
	    write(*,*)
     +' reaction field method for electrostatic interactions'
	    write(*,*) 
	    write(*,*)' macroscopic dielectric constant ',DIEL 
	    if(RDEB.gt.1.d-3)write(*,*)' Debue length ',RDEB
	    write(*,*)
	  end if
         end if
C   Set Ewald's alpha to zero
	 ALPHAD=0.d0   
      else
* 8.3. Ewald method 
* ----------------- 
C Define "alpha" parameter - from cut-off distance or from the box size
C ALPHAD is "alpha" which is used in Ewald formula  
C "Kmax" parameter will be define latter in subr. INIFUR (l-forces.f) 
	if(RCUTT.le.0.5*CL)then
C  Define alpha from cut-off
          ALPHAD	= ALPHA/RCUTT
	else
C  Define alpha from the box size
	  ALPHAD	= 2.d0*ALPHA/CL 
	end if  
	LMNBL=.false.
      end if
*
* 8.4 Calculate long-range corrections from Lennard-Jones interactions 
* --------------------------------------------------------------------
      CALL ELRCLJ       !  l-forces.f
*
* 9. Report basic simulation parameters
* -------------------------------------    
      if(TASKID.eq.MAST)then
        write(*,'(a)')'   Concentrations of different molecules:'
	do ITYP=1,NTYPES
	  CONC=NSPEC(ITYP)/(AVSNO*1.d-27*VOL)
	  write (*,'(4x,a4,f9.4,a4)')NAME(ITYP)(1:4),CONC,' M/l'  
	end do
*     
	PRINT "(/'*** DENSITY           ',D12.6,' g/cm**3 ' )",RHO
        PRINT "( '*** BOX VOLUME        ',D12.6,' A**3    ' )",VOL
        PRINT "( '*** BOX LENGTHS       ',3(D12.6,' A ') )"
     +  ,BOXL,BOYL,BOZL
        PRINT "( '*** CUT-OFF           ',D12.6,' A       ' )",RCUTT
	PRINT "( '*** Closest R         ',D12.6,' A       ' )",SHORT
	PRINT "( '*** RDF cut-off       ',D12.6,' A       ' )",RDFCUT
	PRINT "( '*** ALPHA (Ew.sum)    ',D12.6)",ALPHA  
	PRINT "( '*** Fexp  (Ew.sum)    ',D12.6)",FEXP
        PRINT *
      if(.not.LRST)then
	PRINT "(/'*** LONG  TIME STEP   ',F12.6,' fs ' )",DT*1.d15
        if(.not.LSHEJK)PRINT 
     +"('*** SHORT TIME STEP   ',F12.6,' fs' )",DT*1.d15/(NFREQ+1.d-14)    
	end if
        if(LSCFT)write(*,*)' Large forces will be cut on level ',FTSCL
        write(*,*)
      end if  
      FTSCL=1./FTSCL
*
*   10.  Finalise preparation
*   -------------------------
*   10.1  Prepare RDF calculations
*   ------------------------------
      IF(LGR)   CALL RDFINP            !  supl.f
*   10.2  Prepare TCF calculations
*   ------------------------------
      IF(LCF)   CALL TCFINP            !   tcf.f
*   10.3  Prepare list of bound atoms
*   ---------------------------------
      CALL BNDLST                      !   units.f
*   10.4  Change temperature and/or density if nesessary
*   ----------------------------------------------------
C   This really works in the case of continuation of an old run,
C   if one of parameters LCHT or LCHP is .true. 
C   One can change temperature or density (by scaling vecocities or
C   distances) after restart and continue simulations with new TEMP or RHO
C   specified in the input file   
      if(LCHP)RHO=RHON
      call CHTERM(BOXLN,BOYLN,BOZLN)       !  service.f
*   10.5  Bind some atoms to specified locations (if LRR=.true.)
*   -----------------------------------------------------------
      if(LRR)call PULLINI                         ! setup.f
*   10.6  Check configuration on non-overlapping
*   --------------------------------------------
      call TESTCONF                        !  service.f
*   10.7  Verbose: Dump start configuration
*   ---------------------------------------
C   6 here is standard output channel     
      if(IPRINT.ge.8)call CONFDUMP(6)       !  service.f
*   10.8  Gather molecule in one cell
      if(LGATHER)call GROUP
*
*   11.   We do not really wanted to run MD, but only: (LINPUT=.t.)
*   --------------------------------------------------
      IF(LINPUT)then
*   11.1. ... get results of previously made run
*   --------------------------------------------  
	if(LRST)then
        if(IHIST.gt.NHIST)IHIST=NHIST
	  call GETAVR(3)
	  IF(LGR.and.LGRST)   CALL RDFOUT(1)
        IF(LCF.and.LCFRST)  CALL TCFOUT(1)
      else  
*   11.2. ... check that input OK
*   -----------------------------
	  if(TASKID.eq.MAST)WRITE(*,*)' STOP after cheking the input' 
	end if
	if(LXMOL)call XMOLDUMP
	go to 9999
      end if
*  
*   12.  Run MD
*   -----------              
      if(TASKID.eq.MAST)
     +PRINT "('*** MD STEPS ABOUT TO START FOR',I8,' STEPS:')",NSTEPS
      if(TASKID.eq.MAST)PRINT * 
      if(IPRINT.ge.5)WRITE(*,*)'  STEP     POTE      PE(int)   E(el)'
     +,'   Etot       TEMP     PRES      Dens '
C   start dumping trajectory file
      if(.not.LRST.and.TASKID.eq.MAST)call TRACE
C   wait for all nodes to be ready to start
      call BARRIER
*  12.1.  SHAKE algorithm
*  ---------------------- 
      IF(LSHEJK) THEN
        CALL SULTAN
      ELSE
*  12.2.  Double time step algorithm  
*  ---------------------------------
        CALL DOUBLE
      END IF
      write(*,*) ' main MD loop finished'
*
*  13. Report final results
*  ------------------------      
      if(IHIST.gt.NHIST)IHIST=NHIST
*  13.1. Report calculated averages 
      call GETAVR(3)                    !  aver.f
*  13.2. Report RDF-s   
      IF(LGR)   CALL RDFOUT(1)          ! aver.f
      if(TASKID.ne.MAST)go to 9999
*  13.3. Report TCF-s 
      IF(LCF)   CALL TCFOUT(1)          !  tcf.f
*  13.4. Dump final configuration in "Xmol" format   
      if(LXMOL) call XMOLDUMP           !  restart.f
*  14. Finalise (close up MPI connections) 
 9999 if(LVISUAL)then
        do I=1,100
          write(*,*)'@mm'
        end do
      end if
      call FINAL
*
      END
*
