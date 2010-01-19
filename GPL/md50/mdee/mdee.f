      PROGRAM MDEE
*
      include "mdee.h"
*-----------------------------------------------
      DIMENSION  QP(4)
*
      NAMELIST /INPUT/ NTYPES,NSPEC,NSITS,FNAME,TPRES
     X         ,TRTEMP,RHO,DT,RCUT,TDELT,FTSCL,QT,QPR
     X         ,BOXL,BOYL,BOZL,BDTOL,TOLER,C14LJ,C14EL
     X         ,NSTEPS,NFREQ,IPOT,NAVER,IPRINT,NAMOL
     X         ,ITORS,NDUMP,IAVER,ICHECK,IHIST,ALPHA,FEXP
     X         ,LRST,LDMP,LSCLT,LSCFT,LNVT,LNPT,LCHECK,L14NB
     X         ,L0AVS,L0VS,LSHEJK,LCHP,LCHT,LXMOL,LECH,L15NB
     X         ,SHORT,ICHNB,LGR,LGDMP,LGRST,ALLSTS,RDFCUT,PATHDB
     X         ,LCFDMP,LCFRST,NSTEG,JUMP,N1,N2,MAX
     +         ,ME0,NE,IEE,NCEE,PLOW
*  
      character*20 FNM, FINP*24,PATHDB*64,DATE*8
      character*128 TAKESTR,AUX
      label = "MDEE_v3.9   "
	DATE  = '10.02.00'
      PRINT *,'*** MOLECULAR DYNAMICS',
     X' PROGRAM FOR MIXTURES OF RIGID or FLEXIBLE MOLECULES'
      PRINT *,'*** ',
     X '         with removing/inserting of molecules '
	write(*,*)'     (expanded ensemble method)'
      PRINT *
	write(*,*)'Version ',label,'           from ',DATE
	write(*,*)
	time0	= cputime(0.d0)
*
*     Default parameters
*
*Physical parameters (default) for a many-particle system:
*
      RHO        = 1.00              ! density   (if BOXL=0)
      TRTEMP     = 300.0             ! temperature
      TPRES      = 1.                ! pressure (atm)
      NOP        = NPART             ! number of particles
      NFREQ      = 10                ! num. of short time steps in DT
      NAVER      = 100               ! num. of steps for averaging 
      DT         = 1.D-15            ! time step
      RCUT       = 10.D0             ! cut-off
      TDELT      = 50.D0             ! ??
      RDFCUT     = 10.D0             ! cut off for RDF  
      ALPHA	= 2.5                  ! Ewald: erfc(alpha)=0
      FEXP	= 8.                   !        exp(-fexp)=0
      NUMTASK=1
      TASKID=0
      LHEX=.false.
      LOCT=.false.
      FTSCL      = 0.001                ! cutoff large forces 
      QT	 = 10.D0
      QPR	 = 300.d0
      JUMP       = 1
      BOXL       = 0.D0
      BOYL       = 0.D0
      BOZL       = 0.D0
	SHORT	= 4.                   ! R for closest neigbours   
      TOLER      = 1.D-4
      BDTOL      = 1.D-4 
	MAX	= MAXS
      ICHECK	= 1
      LMNBL=.false.
	IHIST	= 1
	ICHNB	= 10    ! cheking neigbours   list
	NCEE	= 1     ! num. of steps followed by change of subensemble
*
*Some control parameters
*
      FNAME	 = "md"
      LCHECK     = .FALSE.              ! true - ony input checking
      LNVT       = .FALSE.
      LNPT       = .FALSE.
      LSCLT      = .TRUE.
      LSCFT      = .FALSE.
      L0AVS      = .FALSE.
      L0VS       = .FALSE.
      LRST       = .FALSE.
      LDMP       = .FALSE.
      LCHL       = .FALSE.
      LCHP       = .FALSE.
      LCHT       = .FALSE.
	LECH	= .false.
      MSTEP      =  0
*     
      DO I       = 1,NTPS
        ITORS(I)   = 0
        NRB(I)     = 0
        C14LJ(I)   = 0.D0
        C14EL(I)   = 0.D0
        L14NB(I)=.true.
        L15NB(I)=.true.
      END DO
*
*Read in new values
*
      READ(*,INPUT)
      if(.not.LRST)then
	  L0AVS=.true.
	  ME=ME0
      end if
      do IE=1,NE
	read(*,*)EC(IE),EE(IE) 
	if(IPRINT.ge.7)write(*,*)EC(IE),EE(IE)
*  some historical reasons...
        EE(IE)=-EE(IE)
      end do    
      if(LCHECK)L0AVS=.false.
      if(MAX.gt.MAXS)stop "increase MAXS"
      FTSCL=1.d0/FTSCL
*
*   file names
*
      IC=1
      do I=1,20
        if(FNAME(I:I).ne.' ')IC=I
      end do
      IC1=IC+1
      IC2=IC+4
	ICFN=IC+2
      FNM='                    '
	write(FDUMP,'(24(1h  ))')
	write(FLST,'(24(1h  ))')
	write(FRDF,'(24(1h  ))')
	write(FTRK,'(24(1h  ))')
      write(FBIO,'(24(1h  ))')
      write(FINP,'(24(1h  ))')
      FNM(1:IC)=FNAME(1:IC)
      FDUMP(1:IC)=FNM
      FDUMP(IC1:IC2)='.dmp'
      FLST(1:IC)=FNM
      FLST(IC1:IC2)='.lst'
      FRDF(1:IC)=FNM
      FRDF(IC1:IC2)='.rdf'
      FBIO(1:IC)=FNM
      FBIO(IC1:IC2)='.car'
      FINP(1:IC)=FNM
      FINP(IC1:IC2)='.inp'
	CLN=BOXL
*
*---------------- DEFINE MOLECULE TYPES -------------------------------
*
      SUMM         = 0.D0
      NRBB         = 0
      NRAA         = 0
      NRTT         = 0
      NRII         = 0
      NX           = 0
      IADB(1)      = 1
      IADA(1)      = 1
      IADT(1)      = 1
      TOTMAS       = 0.D0
	NEFAC=0
*
*  Input of molecules from data base
        NTYP=0
*     
	do ITYP=1,NTYPES
        CALL READMOL (SUMM,NRBB,NRAA,NRTT,NRII,NX,NTYP,PATHDB,
     +  NAMOL(ITYP))
        TOTMAS       = TOTMAS+SUMM*NSPEC(ITYP)
        SUMMAS(NTYP) = SUMM
	  NAME(ITYP)=NAMOL(ITYP)(1:6)
        if(IEE(ITYP).eq.0)NEFAC=NEFAC+NSPEC(ITYP)
      end do
*     
      IF(NX  .GT.NS  ) STOP '!!! INCREASE NS   !'
      IF(NTYPES.GT.NTPS) STOP '!!! INCREASE NTPS !'
      IF(NRBB.GT.NB  ) STOP '!!! INCREASE NB   !'
      IF(NRAA.GT.NA  ) STOP '!!! INCREASE NA   !'
      IF(NRTT.GT.NT  ) STOP '!!! INCREASE NT   !'
      IF(NRII.GT.NI  ) STOP '!!! INCREASE NI   !'
      IF(TOTMAS.EQ.0 ) STOP '!!! NO PARTICLES WITH MASSES DEFINED !'
*
      IF(LSHEJK) THEN
      NFREQ = 1
      PRINT "(/'*** CONSTRAINED DYNAMICS WILL BE USED! ***'/)"
      END IF
	NRBT=0.
	do ITYP=1,NTYPES
	  NRBT=NRBT+NRB(ITYP)*NSPEC(ITYP)
	end do
*
      NOP     = 0
      NSITES  = 0
      NSTOT   = 0
      IATOM=0
      IMOL=0
      FNSTR	= 0.
      IST=0
      TOTCH=0.
      if(IPRINT.ge.9)then
	  write(*,*)
	  write(*,*)' Check references:'
	  write(*,*)' Iatom    Imol    Isit   Ityp     Name      Nb  '
	end if
	DO ITYP = 1,NTYPES
        NOP     = NOP    + NSPEC(ITYP)
        NSITES  = NSITES + NSITS(ITYP)
        if(NSITS(ITYP).eq.2)FNSTR=FNSTR+2.d0*NSPEC(ITYP)
        if(NSITS(ITYP).ge.3)FNSTR=FNSTR+3.d0*NSPEC(ITYP)
        NSTOT   = NSTOT  + NSPEC(ITYP)*NSITS(ITYP)
        FNSS3   = 3.D0*DFLOAT(NSITES)
        IF(NSPEC(ITYP).LE.0) STOP '!!! CHECK NSPEC !!!'
                 TFACT(ITYP) = 1.D0
        IF(LSHEJK) TFACT(ITYP) = FNSS3/(FNSS3-NRB(ITYP))
        NRTSC(ITYP) = 0
        write(*,'(a6,i4,a3,a6,3x,a5,i5,4x,a6,i4)')
     +' type ',ITYP,' - ',NAME(ITYP),' Nmol',NSPEC(ITYP),
     +'Nsites',NSITS(ITYP)
	  do IS=1,NSPEC(ITYP)
	    IMOL=IMOL+1
            do I=1,NSITS(ITYP)
	      ISITE=IST+I
	      IATOM=IATOM+1
	      NNUM(IATOM)=IMOL
              ITS(ISITE)=ITYP
	      NSITE(IATOM)=ISITE
	      ITYPE(IATOM)=ITYP
	      TOTCH	= TOTCH+CHARGE(ISITE)
	      if(IPRINT.ge.9)write(*,'(4I6,2x,2a4,a6,2x,I5)')
     +      IATOM,IMOL,ISITE,ITYP,NM(ISITE),' on ',NAME(ITYP),NRB(ITYP)
	    end do
	  end do
	  IST=IST+NSITS(ITYP)
      END DO! OF ITYP
      write(*,*)
      write(*,'(a,f8.3)')'*** Total charge            Q =  ',TOTCH
      IF(NOP.EQ.1) RHO = 0.D0
*
      PRINT "(/'*** NR OF SITES/MOL (NSITES)  ',I6)" ,NSITES
      PRINT "( '*** NR OF SITES/TOT (NSTOT)   ',I6)" ,NSTOT
      PRINT "( '*** NR OF PARTICLES (NOP)     ',I6)",NOP
	PRINT "( '*** NR OF BONDS (NRBT)        ',I6)",NRBT
*
      IF(NOP   .GT.NPART) STOP '!!! INCREASE NPART !'
      IF(NSITES.GT.NS   ) STOP '!!! INCREASE NS    !'
      IF(NSTOT .GT.NTOT ) STOP '!!! INCREASE NTOT  !'
*
*  FNST - number of degrees of freedom 
*  Obs!  total momenta is fixed! 
      FNST        = 3.d0*DFLOAT(NSTOT)-3.d0
      if(LSHEJK)FNST=FNST-NRBT
	FTIMES      = 1.D0/DFLOAT(NFREQ) 
      FNOP        = DFLOAT(NOP)
      FNSTI	= FNST-3.d0*FNOP-FNSTR+3.d0
      write(*,'(a,f6.0)')'*** NR of DEGREES of FREEDOM ',FNST
      write(*,'(a,f6.0)')'***              rotational: ',FNSTR
      write(*,'(a,f6.0)')'***              internal:   ',FNSTI
*     
	write(*,*)
      CALL IUNITS
      CALL ELRCLJ
	write(*,*)
*
*  concentrations in M/l
*     
	write(*,'(a)')'   Concentrations of different molecules:'
	do ITYP=1,NTYPES
	  CONC=NSPEC(ITYP)/(AVSNO*1.d-27*VOL)
	  write (*,'(4x,a4,f9.4,a4)')
     +  NAME(ITYP)(1:4),CONC,' M/l'  
	end do
*  Zeros
	do I=1,NRH
	  do J=1,LHIST
	    HIST(I,J)=0.
	  end do
	end do
*
*	Parameters for NVT-ensemble:
*
      IF(LNVT) THEN
        DQT	= (UNITT*1.d15/QT)**2
	  write(*,*)
        write(*,'(a,f10.3)')'*** Constant temperature MD: T = ',TRTEMP
        PRINT *,'   NOSE-NVT RELAXATION TIME      ',QT,' fs'
        LSCLT      = .FALSE.
      END IF!(LNVT...
*
*	Parameters for NPT-ensemble:
*
      IF(LNPT) THEN
        if(.not.LNVT)write(*,*)'!!! Attention: NPE-ensemble???'
        DQT	= (UNITT*1.d15/QT)**2
        DQP	= UNITE*UNITT**2/(UNITP*FNOP*BOLTZ*TRTEMP*(QPR*1.d-15)**2)
        write(*,'(a,f10.3,a)')
     +  '*** Constant-pressure MD: P = ',TPRES,' atm'
        PRINT *,'   NOSE-NPT RELAXATION TIME     ',QPR,' fs'
      END IF!(LNVT...
        SC	 = 0.D0
        SCL	= 0.D0
*
Continue an old run:
*
      IF(LRST) THEN
*
      CL=BOXL
      CALL RESTRT(1)
      if(LCHP)BOXL=CL
      call CHTERM
      call RECLEN(BOXL)
      write(*,'(a,3f12.4)')' box: ',BOXL,BOYL,BOZL
      call RECINT
*     
	if(LECH)then
	  do I=1,NE
	    EE(I)=EE(I)*NDEL
	  end do
	end if 
      CALL GETCOM
      write(*,'(a,f12.4,a2)')
     +' The system was simulated ',TIM*1.d12,'ps'
      write(*,'(a,f12.4,a2)')
     +' Averages was collected   ',TIMA*1.d12,'ps'
*
      ELSE
*
*Prepare a new run:
*
*Generate start coordinates for C.O.M. points:
*
      NHIST	= 1
      TIM	= 0.
      timest	= 0.
	ICHE=0
	ICHU=0
	NTRKF=0           ! num. of trace file
	NTRKZ=0           ! num. of points in the trace file
      IF(NOP.EQ.1) THEN
      X(1)   = 0.D0
      Y(1)   = 0.D0
      Z(1)   = 0.D0
*
      ELSE
*
	  if(ICHECK.eq.0)then
*  7.2.2.5  read atom coord. from input file
	    open(unit=13,file=FINP,status='old',err=298)
	    go to 299
C   something goes wrong:
C   can not read input file...         
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
C    reading...
 299        do I=1,NSTOT                                
		AUX=TAKESTR(13,IE)
	        read(AUX,*,err=297,end=297)SX(I),SY(I),SZ(I) 
	      end do          
	      if(TASKID.eq.MAST)write(*,*)
     +'read initial coordinetes from file ',FINP
	      close(13)  
              call GETCOM
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
            else IF(ICHECK.eq.1) THEN
              CALL FCC(X,Y,Z,BOXL,NOP)
            ELSE
              CALL CUBIC(X,Y,Z,FX,FY,FZ,BOXL,NOP)
            END IF!(ICHE...
*
      END IF!(NOP...
*
*Initialize velocities
*
      TMPRED     = EPSIL(1)*504.D0/TRTEMP
      CALL COMVEL(VX,VY,VZ,TMPRED,NSTOT)
*
      PRINT *,'*** RANDOM VELOCITIES FOR   ',NOP,'  PARTICLES'
*
      VSQ         = 0.0
      DO N        = 1,NSTOT
      VSQ         = VSQ+(VX(N)**2+VY(N)**2+VZ(N)**2)*MASSDI(N)
      END DO! OF N
*
      TEMP        = 1.50*VSQ*ENERF*TEMPF/FNST
      FTMP        = SQRT(TRTEMP/TEMP)
      DO     N    = 1,NSTOT
      VX(N)       = VX(N)*FTMP
      VY(N)       = VY(N)*FTMP
      VZ(N)       = VZ(N)*FTMP
      END DO! OF N
*
      VSQ         = 0.0
      DO N        = 1,NSTOT
      VSQ         = VSQ+(VX(N)**2+VY(N)**2+VZ(N)**2)*MASSDI(N)
      END DO! OF N
*
      TEMP        = 1.5*VSQ*ENERF*TEMPF/FNST
*
      PRINT *,'*** INITIAL TEMPERATURE CALCULATED',TEMP
      PRINT *,'*** VELOCITIES INITIATED'
      PRINT *
*
*Orientations of the molecules:
*
      if(ICHECK.ne.0)then
      N              = 0
      I              = 0
      DO ITYP        = 1,NTYPES
      IBES           = ISADR (ITYP)+1
      IENS           = ISADR (ITYP +1)
      DO J           = 1,NSPEC(ITYP)
      I              = I+1
      CALL RANQ(QP)
      DO IS          = IBES,IENS
      N              = N+1
      RXI            = R(IS,1)
      RYI            = R(IS,2)
      RZI            = R(IS,3)
      CALL ROTATE  (QP,RXI,RYI,RZI)
      SX(N)          = RXI+X(I)
      SY(N)          = RYI+Y(I)
      SZ(N)          = RZI+Z(I)
      WX(N)          = RXI
      WY(N)          = RYI
      WZ(N)          = RZI
      GX(N)          = 0.D0
      GY(N)          = 0.D0
      GZ(N)          = 0.D0
      HX(N)          = 0.D0
      HY(N)          = 0.D0
      HZ(N)          = 0.D0
      OX(N)          = SX(N)
      OY(N)          = SY(N)
      OZ(N)          = SZ(N)
      END DO! OD IS
      END DO! OF I
      END DO! OF ITYP
*
*
      PRINT "('*** MOLECULAR ORIENTATIONS INITIALIZED')" 
      end if
      NSTTOT	= 0      
      CL=BOXL
      call RECLEN(CL)	    
	call RECINT
	  do I=1,NE
	    EE(I)=EE(I)*NDEL
	  end do
*
      END IF!(LRST...
* Parameter for Ewald sum
	if(RCUTT.le.0.5*BOXL)then
        ALPHAD	= ALPHA/RCUTT
	else
	  ALPHAD	= 2.d0*ALPHA/BOXL 
	end if 
*
*Zero average accumulators
*     
	if(L0AVS)then
	  NAVT	= 0
	  write(*,*)' New averaging'
          NSTEPA	= 0
          TIMA	= 0.
	  call ZEROAV 
          NSTEPA=0
          DO I        = 1,NRQS
	    do J=1,NE
              AV(I,J)       = 0.D0
              AW(I,J)       = 0.D0
	    end do
          END DO! OF I
          do I=1,NE
            do J=1,NE
              IWALK(I,J)=0
            end do
          end do
	end if
	TRPRES=TPRES
	PRINT "(/'*** DENSITY           ',D12.6,' g/cm**3 ' )",RHO
      PRINT "( '*** BOX VOLUME        ',D12.6,' A**3    ' )",VOL
      PRINT "( '*** BOX LENGTH        ',D12.6,' A       ' )",BOXL
      PRINT "( '*** CUT-OFF           ',D12.6,' A       ' )",RCUTT
	PRINT "( '*** Closest R         ',D12.6,' A       ' )",SHORT
	PRINT "( '*** ALPHA (Ew.sum)    ',D12.6)",ALPHA  
	PRINT "( '*** Fexp  (Ew.sum)    ',D12.6)",FEXP
      PRINT *       
*
      IF(LGR)   CALL RDFINP
*
      DO ITYP = 1,NTYPES
      CALL BNDLST(ITYP)
      END DO! OF ITYP
* Check for electroneutrality

*
      IF(LCHECK)then
	  if(LRST)then
	    call GETAVR(3)
	    IF(LGR.and.LGRST)   CALL RDFOUT(1)
        else  
	    WRITE(*,*)' STOP after cheking the input' 
	  end if
	  if(LXMOL)call XMOLDUMP
	  go to 9999
      end if
	PRINT "('*** MD STEPS ABOUT TO START FOR',I8,' STEPS:')",NSTEPS
      PRINT * 
*
      if(IPRINT.ge.5)WRITE(6,900)
* A MD step 
      IF(LSHEJK) THEN
        CALL SULTAN
      ELSE
        CALL DOUBLE
      END IF
*
      if(LXMOL)call XMOLDUMP
      if(IHIST.gt.NHIST)IHIST=NHIST
      call GETAVR(3)
      IF(LGR.and.LGRST)   CALL RDFOUT(1)
*
  900 FORMAT(/1X,'  STEP     POTE      PE(int)   E(el)'
     X,'   Etot       TEMP     PRES      Dens '/)
*
      IF(.NOT.LNVT) THEN
      PRINT *,'... TEMP SCALED:'
      DO ITYP = 1,NTYPES
      PRINT *,'...',NRTSC(ITYP),'  TIMES IN'
     X,NSTEPS,'  STEPS - FOR TYPE: ',ITYP
      END DO! OF ITYP
      END IF
 9999 stop
*
      END
*
