*=====================================================================
*     MDynaMix v.5.0
*     Part 3
*
*     Simulation setup
*     ----------------
C
C   This file contains subroutines setting up simulation:
C
C     1. IUNITS  - set up units, addresses, references
C     2. BNDLST  - create list of bound, nonbound and
C                  1-4 neighbouring atoms
C     3. CUBIC   - set up molecular COM on cubic lattice
C     4. FCC     - set up molecular COM on FCC lattice
C     5. DNASOL  - put some molecule(s) in a cylindrical hole and
C                  distribute other around
C     6. SPHOLE  - put some molecule(s) in a spherical hole and
C                  distribute other around
C     7. MIXTURE - mix molecules randomly
C     8. PULLINI - initialization of restraints
C
*=============== IUNITM ===============================================
* 
*
      SUBROUTINE IUNITS
      include "prcm.h"
*
*   1.1 Calculate box dimensions:
*   ----------------------------
*  1.1.1  Some constants 
      FACTM       = AVSNO*1000.0
      TOTMAS      = TOTMAS/FACTM             ! in kg
      FNOP        = DFLOAT(NOP)
      AVOFNR      = AVSNO/FNOP
      IF(RHO.NE.0.D0) THEN
        VOLM        = AVOFNR*(TOTMAS/(RHO*1000.0))
*  1.1.2 Define box size from density
        IF(BOXL*BOYL*BOZL.lt.1.D-20) THEN
	   if(TASKID.eq.MAST)then
              write(*,*)' Define box size from the density'
              write(*,*)
           end if
           VOL = VOLM*FNOP*1.d30/AVSNO    
	   if(ICELL.eq.1)VOL=2.*VOL
	   if(ICELL.eq.2)VOL=VOL/cos30
           CL = VOL**(1.0/3.0)       
	   BOXL=CL
	   BOYL=CL
	   BOZL=CL
C  Truncated octahedron geometry
	   if(ICELL.eq.1)VOL=0.5*VOL
C  Hexagonal cell geometry
	   if(ICELL.eq.2)then
	     BOYL=BOXL*cos30
	     VOL=BOXL*BOYL*BOZL
	   end if
*  1.1.3  Box size difened from the input file    
        else
          CL=BOXL
        ENDIF
      ELSE
*  1.1.4 Case of vacuum simulation
        BOXL        = 1000.D0
        BOYL        = 1000.D0
        BOZL        = 1000.D0
C  No treatment of long-range electrostatics  
	ALPHA	      = 0.d0
        if(TASKID.eq.MAST)PRINT "(/,'*** VACUUM SIMULATION ***',/)"
      END IF          
      if(ICELL.eq.0.and.TASKID.eq.MAST)
     +write(*,*)' Rectangular geometry is used'
      if(ICELL.eq.1.and.TASKID.eq.MAST)
     +write(*,*)' Truncated octahedron geometry is used'
      if(ICELL.eq.2.and.TASKID.eq.MAST)write(*,*)
     +' Hexagonal periodic cell is used'
*  1.1.5 Calculate size-dependent variables
      call RECLEN(ODIN,ODIN,ODIN) 
      VOLM        = 1.d-30*VOL*AVOFNR
      DT = DT*1.d-15          ! in sec
      TTREK=TTREK*1.d-15
*
*  1.2 Define different constants
*  ----------------------------    
*  1.2.1 Internal units:
      UNITL       = 1.d-10                !  m
      UNITM       = TOTMAS             ! mass of the system in kg
      UNITT       = DT                 ! time step
      UNITE       = UNITM*(UNITL/UNITT)**2
      UNITP       = UNITE*1.e-5/UNITL**3 
* 1.2.2 Dimensionless quantities and auxilarly constants
      TSTEP       = UNITT/SQRT((UNITM/UNITE)*UNITL**2)
      HTSTEP      = 0.5*TSTEP
      TSTEB       = TSTEP*FTIMES
      HTSTEB      = 0.5*TSTEB
      BETA	  = UNITE/(BOLTZ*TRTEMP)
      ENERF       = AVSNO*UNITE
      TEMPF       = 2.0/(3.0*AVSNO*BOLTZ)
      COULF       = (ELECHG/EPS0)*(ELECHG/UNITE)/(4.D0*PI*UNITL)
      PERMOL      = ENERF*1.D-3/FNOP 
      SHORT2	    = SHORT**2
      SHORTA	    = SHORT-1.d-11/UNITL !  -0.1A
      RLR	    = 1.d0/(BETA*RLR**2) 
      if(.not.LGRST.and.RDFCUT.le.0.d0)RDFCUT=HBOXL 
      RDFCUT2=RDFCUT**2
* 1.2.3 Some output
	if(IPRINT.ge.7)then
        PRINT "(/'*** UNIT MASS         ',D12.6,' kg      ' )",UNITM
        PRINT "( '*** UNIT TIME         ',D12.6,' S       ' )",UNITT
        PRINT "( '*** UNIT ENERGY       ',D12.6,' J       ' )",UNITE
        PRINT "( '*** UNIT LENGTH       ',D12.6,' m       ' )",UNITL
        PRINT"('*** UNIT PRESSURE     ',D12.6,' 10**5 N/m**2  '/)",UNITP
        PRINT "( '*** ENERGY  FACTOR    ',D12.6,' J/mol   ' )",ENERF
        PRINT "( '*** TEMP    FACTOR    ',D12.6,' K/J     ' )",TEMPF
        PRINT "( '*** COULOMB FACTOR    ',D12.6,'         '/)",COULF
        PRINT "( '*** 1/kT (beta)       ',D12.6,' u.e.    ' )",BETA
        PRINT "( '*** LONG  TIME STEP   ',F12.6,' .I.U.   ' )",TSTEP
        if(.not.LSHEJK)PRINT 
     +"( '*** SHORT TIME STEP   ',F12.6,' .I.U.   ' )",TSTEB
	end if
*
*  1.2.4  Calculate energy/temperature convertion factors
*  ---------------------------------------------------
C  These factors are calculated for each molecule type taking into
C  account possible constraints. 
      CONVET = 3.*ENERF*TEMPF/FNST
      do ITYP=1,NTYPES
        FNSS3   = 3.*DFLOAT(NSITS(ITYP))
        FNSS1   = NSTOT/(NSTOT-1.)
        TFACT(ITYP) = ENERF*TEMPF*FNSS1/(NSITS(ITYP)*NSPEC(ITYP))
        if(LSHEJK.and.ISHEJK(ITYP).ne.0)
     +    TFACT(ITYP)=TFACT(ITYP)*FNSS3/(FNSS3-NRB(ITYP))
        NRTSC(ITYP) = 0       ! will report num. of temperature scalings
      end do
*
*  1.2.5 dimensionless mass parameters:
      DO IS      =  1,NSITES
        MASSD (IS)  = MASS (IS)/FACTM/UNITM  ! fraction of the total mass 
      END DO! OF I
*
*  1.3 Addresses and references
*  ----------------------------
*  1.3.1 Addresses arrays
C  See comment to p.4 in main.f
C  The addresses provide atom number corresponding 
C  to each new molecule, site ,type
      IATOM   = 0
      IMOL    = 0 
      IST     = 0
      IADDR (1)      = 0
      ISADR (1)      = 0
      ISADDR(1)      = 0
      TOTCH	= 0. 
      NTYP1 = NTYPES + 1
      DO ITYP        = 1,NTYPES
        IADDR (ITYP+1) = IADDR (ITYP)+NSPEC(ITYP)              ! molecules
        ISADR (ITYP+1) = ISADR (ITYP)+NSITS(ITYP)              ! sites
        ISADDR(ITYP+1) = ISADDR(ITYP)+NSITS(ITYP)*NSPEC(ITYP)  ! atoms
*
*  1.3.2 Setup reference arrays
*  --------------------------    
C Now, we can create list of references
C
C Reference arrays
C	NNUM          atom -> molec
C       NSITE         atom -> site
C       ITYPE         atom -> type
C       ITS           site -> type
C       ITM           molec -> type
	if(TASKID.eq.MAST)
     +  write(*,'(a6,i4,a3,a16,3x,a5,i5,4x,a6,i4)')
     +  '  type ',ITYP,' - ',NAME(ITYP),' Nmol',NSPEC(ITYP),
     +  'Nsites',NSITS(ITYP)
	do IS=1,NSPEC(ITYP)
	  IMOL=IMOL+1
          do I=1,NSITS(ITYP)
	      ISITE=IST+I
	      IATOM=IATOM+1
	      NNUM(IATOM)=IMOL                  ! atom -> molec
	      NSITE(IATOM)=ISITE                ! atom -> site
	      ITYPE(IATOM)=ITYP                 ! atom -> type
	      ITS(ISITE)=ITYP                   ! site -> type
	      if(IPRINT.ge.10)write(*,'(4I7,3x,2a4,a6)')
     +        IATOM,IMOL,ISITE,ITYP,NM(ISITE),' on ',NAME(ITYP)  
	      TOTCH	= TOTCH+CHARGE(ISITE)
	  end do  
	  ITM(IMOL)=ITYP                      ! molec -> type
	end do
	IST=IST+NSITS(ITYP)
      END DO! OF ITYP
*
*  1.4  Report some quantities
*  ---------------------------
      if(TASKID.eq.MAST)then
        write(*,*)
        write(*,'(a,f8.3)')'*** Total charge            Q =  ',TOTCH
      end if
*     
      if(IPRINT.ge.10)then            
        write(*,*)
        PRINT "(25X,'*** ADDRESSES: ***')"
        DO ITYP        = 1,NTYP1
          PRINT "(10X,4(1X,I8))",
     +    ITYP,ISADR(ITYP),IADDR(ITYP),ISADDR(ITYP)
        END DO! OF ITYP
      end if
*
*   1.5 Setup values of different working arrays
*   -------------------------------------------
      N              = 0
      MT             = 0
*   1.5.1 Inverse masses, charges
      DO ITYP        = 1,NTYPES
        QS             = 0.D0
        NSB            =  ISADR (ITYP)+1
        NSE            =  ISADR (ITYP+1)
        DO I           = 1,NSPEC(ITYP)
          DO IS          = NSB,NSE
            N              = N+1
            MASSDI(N)      = 1.0/MASSD(IS)
            Q     (N)      =    CHARGE(IS)
            QS             = QS+DABS(Q(N))
          END DO! OF I
        END DO! OF IS
        QSUM(ITYP)     = QS
*    1.5.2 Type pair indeces; energy components
        DO JTYP        = ITYP,NTYPES
          MT             = MT+1
          MDX(ITYP,JTYP) = MT
          MDX(JTYP,ITYP) = MT
          POTES(MT)      = 0.D0
          POTLJ(MT)      = 0.D0
        END DO! OF JTYP
      END DO! OF ITYP      
      MOLINT         = MT+1
      POTES(MOLINT)=0.
      POTLJ(MOLINT)=0.
C This counts energy absorbed from the thermostat 
      EABS = 0.
      EABSS = 0.
*    1.5.3 External field parameters
C  Frequence and amplitude of the external field in int.units
      if(IEXT.eq.1)then
        write(*,'(a,e12.5,a,e12.5,a)')' External filed: ',EXAMPL,
     +' V/cm;  Frequence ',EXFREQ,' hz'
        EXAMPL=EXAMPL*UNITL*1.d2*ELECHG/UNITE           ! V/cm -> i.u.
        write(*,*)' In internal units: ',EXAMPL
      end if
*
      if(IPRINT.ge.9)
     +PRINT "(/' *** ',I2,' MOLECULAR INTERACTIONS ***'/)",MT
*
      DO ITYP         = 1,NTYPES
        DO JTYP         = ITYP,NTYPES
          MT              = MDX(ITYP,JTYP)
          if(IPRINT.ge.9)
     +PRINT "(' NR: ',I2,5X,A6,'  <-->  ',A6)",MT,NAME(ITYP),NAME(JTYP)
        END DO! OF JTYP
      END DO! OF ITYP
      LMOL1=.false.
      LMOL2=.false.
*
*   1.6 Creating bonds, angles, torsion lists - for each node 
*   ---------------------------------------------------------
*   1.6.1 Loop over molecule types 
      NRBT=0    
      NRAT=0    
      NRRT=0
      NRIT=0    
      ICNB=0
      ICNA=0
      ICNT=0
      if(LBONDL)then
         open(unit=26,file="bond.list",status='unknown')
         write(26,'(a,a)')'# Bond list for job ',FNAME
      end if
      do ITYP=1,NTYPES
	NSS=NSPEC(ITYP) 
*   1.6.2 Create bond list
	if(IPOT(ITYP).eq.0)then    ! otherwise special subroutine used (STRBND)
	  NBS=NRB(ITYP) 
	  JSFB=IADB(ITYP)-1
	  NAS=NRA(ITYP) 
	  JSFA=IADA(ITYP)-1
	  NTS=NRT(ITYP) 
	  JSFT=IADT(ITYP)-1
 	  do I=1,NSS                 ! over molecules
            IML=ISADDR(ITYP)+(I-1)*NSITS(ITYP)
	    do J=1,NBS               ! over bonds
	      IS=JSFB+J
	      IBD=ICNB+(I-1)*NBS+J
	      if(mod(IBD,NUMTASK).eq.TASKID)then
		NRBT=NRBT+1
	        if(NRBT.gt.NBO)stop "increase NBO"
		IBI(NRBT)=IB(IS)+IML
		JBJ(NRBT)=JB(IS)+IML
	        IBUK(NRBT)=IS
       	      end if
C  write bond list to file "bonds.list" if specified
              if(LBONDL)then
                IS1=IB(IS)+ISADR(ITYP)
                IS2=JB(IS)+ISADR(ITYP)
                write(26,'(2I6,6x,A4,A3,2A4,A6,2x,3I5,a)')
     +  IB(IS)+IML,JB(IS)+IML,NM(IS1),' - ',NM(IS2),' at ',NAME(ITYP),I
              end if
	    end do 
*    1.6.3  Create covalent angle list
	    do J=1,NAS   
	      IS=JSFA+J
	      IBD=ICNA+(I-1)*NAS+J
	      if(mod(IBD,NUMTASK).eq.TASKID)then
		NRAT=NRAT+1                
	        if(NRAT.gt.NAO)stop "increase NAO"
		IAI(NRAT)=IA(IS)+IML
		JAJ(NRAT)=JA(IS)+IML
		KAK(NRAT)=KA(IS)+IML
	        IAUK(NRAT)=IS
	      end if
	    end do 
*    1.6.4  Create dihedral angle list
	    do J=1,NTS   
	      IS=JSFT+J
	      IBD=ICNT+(I-1)*NTS+J
	      if(mod(IBD,NUMTASK).eq.TASKID)then
		NRRT=NRRT+1
	        if(NRRT.gt.NTO)stop "increase NTO"
		ITI(NRRT)=IT(IS)+IML
		JTJ(NRRT)=JT(IS)+IML
		KTK(NRRT)=KT(IS)+IML
		LTL(NRRT)=LT(IS)+IML 
	 	ITUK(NRRT)=IS
	      end if
	    end do 
	  end do
	  ICNB=ICNB+NSS*NBS
	  ICNA=ICNA+NSS*NAS
	  ICNT=ICNT+NSS*NTS
        else
*   1.6.5  Generate list of bonds for special potential model
C          (only for output; not used in MD)   
          if(LBONDL)then
	    NBS=NRB(ITYP) 
	    JSFB=IADB(ITYP)-1
 	    do I=1,NSS                 ! over molecules
              IML=ISADDR(ITYP)+(I-1)*NSITS(ITYP)
	      do J=1,NBS               ! over bonds
	        IS=JSFB+J
                IS1=IB(IS)+ISADR(ITYP)
                IS2=JB(IS)+ISADR(ITYP)
                write(26,'(2I6,6x,A4,A3,2A4,A6,2x,3I5,a)')
     + IB(IS)+IML,JB(IS)+IML,NM(IS1),' - ',NM(IS2),' at ',NAME(ITYP),I
	      end do 
            end do
          end if
	end if     
      end do
C      write(*,*)' node ',TASKID,' angles ',NRAT  
C	do I=1,NRAT
C	  write(*,*)TASKID,I,IAI(I),JAJ(I),KAK(I),IAUK(I)
C	end do
C      write(*,*)' node ',TASKID,' torsions ',NRRT
C	do I=1,NRRT
C	  write(*,*)TASKID,I,ITI(I),JTJ(I),KTK(I),LTL(I),ITUK(I)
C	end do
*  1.6.5 Calculate lengths of transfer arrays
	LENMB=8*NRBB
	LENMA=8*NRAA
	LENMT=8*NRTT
        LENMP=8*NSTOT
        LENIT=8*MOLINT
        LENTT=8*NTYPES
C  In unused.f :
C  Additional list of bonds
C       call ADDBONDS
*
*  1.7	CONVERSION TO INTERNAL UNITS:
*  ----------------------------------
*  1.7.1 Energy factors
      FACTOR     = 1.D3/AVSNO/UNITE
*      EFACT      = FKCALJ*FACTOR      !  kcal
      EFACT      = FACTOR              ! kJ
*  1.7.2 Bonds force cnstants
C   Harmonic bonds:
      DO I       = 1,NB
        FB(I)      = FB(I)*EFACT
      END DO! OF I
C Morse potentials:
      DO I       = 1,NB
        DB(I)      = DB(I)*EFACT
      END DO! OF I
*                  
* 1.7.3  Angle bending:
*
      DO I         = 1,NA
        FA(I)        = FA(I)*EFACT
        RA(I)        = RA(I)/TODGR
      END DO! OF I
*
*  1.7.4 Torsional angles
*
      DO  I        = 1,NT
        RT (I)       = RT (I)/TODGR
        FT (I)       = FT (I)*EFACT
        FT1(I)       = FT1(I)*EFACT
        FT2(I)       = FT2(I)*EFACT
        FT3(I)       = FT3(I)*EFACT
        FT4(I)       = FT4(I)*EFACT
        FT5(I)       = FT5(I)*EFACT
      END DO! OF I
*  1.7.5 Lennard-Jones parameters:
	do IS=1,NS
	  EPSIL(IS)=2.d0*dsqrt(EPSIL(IS)*EFACT)
          if(LGEOM)then
            SIGMA(IS)=sqrt(SIGMA(IS))
          else
	    SIGMA(IS)=0.5d0*SIGMA(IS)
          end if
	end do
	do IS=1,NNAD
	  EPAD(IS)=2.d0*dsqrt(EPAD(IS)*EFACT)
          if(LGEOM)then
            SIGAD(IS)=sqrt(SIGMA(IS))
          else
	    SIGAD(IS)=0.5d0*SIGAD(IS)
          end if
	end do
*                 
*  1.8 Distribution of atoms over procesors
*  ----------------------------------------   
C  (Relevant only for parallel execution) 
	NAPC=NSTOT/NUMTASK+1
	if(LSHEJK)then
*  1.8.1 distribution of molecules - trying to give equal work for SHEJK
	  I=0 
	  NABS(0)=0
	  NAB(0)=1
	  IWRK=0  
	  AWRK=NRBND*1./NUMTASK
	  do ITYP=1,NTYPES
	    IBEG=IADDR(ITYP)+1
	    IEND=IADDR(ITYP+1) 
	    NSS =NSITS(ITYP)
	    NWRK=NRB(ITYP)
	    if(ITYP.lt.NTYPES)NWRK1=NRB(ITYP+1)
	    do IMOL=IBEG,IEND
            IWRK=IWRK+NWRK   
	      DIFF=IWRK-AWRK*(I+1)
	      if(DIFF.ge.0)then
	        NAE(I)=ISADDR(ITYP)+(IMOL-IBEG)*NSS  
	      if(DIFF.le.0.5*NWRK.or.NAE(I).eq.NABS(I))NAE(I)=NAE(I)+NSS
	        NAP(I)=NAE(I)-NABS(I)
	        NAP3(I)=3*NAP(I)
	        I=I+1
	        if(I.ge.NUMTASK)go to 65
	        NABS(I)=NAE(I-1)
	        NAB(I)=NABS(I)+1
	      end if
	    end do
	  end do
 65	  NAE(NUMTASK-1)=NSTOT
	  NAP(NUMTASK-1)=NSTOT-NABS(NUMTASK-1) 
	  NAP3(NUMTASK-1)=3*NAP(NUMTASK-1)              
	else  ! not.SHEJK
*  1.8.2 Equilibrium distribution of atoms over processors
	  do I=0,NUMTASK-1
	    NABS(I)=NAPC*I
	    NAB(I)=NABS(I)+1
	    NAE(I)=NABS(I)+NAPC 
	    NAP(I)=NAPC 
	    NAP3(I)=3*NAPC
	    NABS3(I)=3*NABS(I)
	  end do   
	  if(NAE(NUMTASK-1).gt.NTOT)stop 'increase NTOT for overhead' 
	  NAE(NUMTASK-1)=NSTOT
	  NAP(NUMTASK-1)=NSTOT-NABS(NUMTASK-1) 
	  NAP3(NUMTASK-1)=3*NAP(NUMTASK-1)              
	end if
	do I=0,NUMTASK-1
	  NABS3(I)=3*NABS(I)
	end do   
	if(IPRINT.ge.7)then
	  write(*,*)' node   Nat   from    to   NABS '
	  do I=0,NUMTASK-1
	     write(*,'(5i7)')I,NAP(I),NAB(I),NAE(I),NABS(I)
	  end do
	end if
      RETURN
      END
*
*=============== BNDLST ===========================================
*
*   2. Prepare lists of bonds, angles, dihedrals 
C      and corresponding reference arrays
      SUBROUTINE BNDLST
*
      include "prcm.h" 
      integer NRNOD(NTPS)
*                                 
      do I=1,NSTOT
	NNBB(I)=0                ! number of atoms bound to this
	do IS=1,NBDMAX
          INBB(I,IS)=0            ! which atoms are bound
	end do
      end do       
      IBEG=TASKID+1 
      NRN=NUMTASK  
      if(NRN.gt.NTYPES)NRN=NTYPES
      do ITYP=1,NTYPES
	NRNOD(ITYP)=mod(ITYP-1,NUMTASK)
      end do  
      if(IBEG.gt.NTYPES)go to 101
      do ITYP=IBEG,NTYPES,NUMTASK
*           
        NSP       = NSPEC(ITYP)     ! num of mol of this type
        NSS       = NSITS(ITYP)     ! num of sites on a mol of this type
        NRBS      = NRB(ITYP)       ! num of bonds 
        NBBEG     = IADB(ITYP)      ! first bond
        NBEND     = NBBEG+NRBS-1    ! last bond             
        ISBEG     = ISADR(ITYP)+1   ! first site
        ISEND     = ISADR(ITYP+1)   ! last site
        ISHF      = ISADR(ITYP)  ! shift of global site num. relatevely local
*
*  First pass: 1-2 neigbours = bonds 
*
        do I=NBBEG,NBEND
          IBB=IB(I)+ISHF
          JBB=JB(I)+ISHF
          if(ID(I).ne.2)call BNDREG(IBB,JBB,1)
        end do
*
*   1 - 3 neigbours: atom + its bonds
*
        do IS=ISBEG,ISEND
          do J=NBBEG,NBEND
            IBB=IB(J)+ISHF
            JBB=JB(J)+ISHF
            if(IS.eq.IBB.and.ID(J).ne.2)then
              do K=J+1,NBEND
                IBK=IB(K)+ISHF
                JBK=JB(K)+ISHF
                if(IS.eq.IBK.and.ID(K).ne.2)then
                  call BNDREG(JBB,JBK,1)
                else if(IS.eq.JBK.and.ID(K).ne.2)then
                  call BNDREG(JBB,IBK,1)
                end if
              end do
            else if(IS.eq.JBB.and.ID(J).ne.2)then
              do K=J+1,NBEND
                IBK=IB(K)+ISHF
                JBK=JB(K)+ISHF
                if(IS.eq.IBK.and.ID(K).ne.2)then
                  call BNDREG(IBB,JBK,1)
                else if(IS.eq.JBK.and.ID(K).ne.2)then
                  call BNDREG(IBB,IBK,1)
                end if
              end do
            end if
          end do
        end do
*
*  1 - 4 neighbours: bonds and its bonds
*  1 - 4 neigbours are marked by minus in the list
* 
        if(.not.L14NB(ITYP))go to 100
        do I=NBBEG,NBEND
          IS=IB(I)
          JS=JB(I)
          if(ID(I).ne.2)then
            do J=NBBEG,NBEND
              if(IS.eq.JB(J).and.ID(J).ne.2)then
                do K=NBBEG,NBEND
                  if(JS.eq.IB(K).and.ID(K).ne.2)then
                    call BNDREG(IB(J)+ISHF,JB(K)+ISHF,-1)
                  else if(JS.eq.JB(K).and.ID(K).ne.2)then
                    call BNDREG(IB(J)+ISHF,IB(K)+ISHF,-1)
                  end if
                end do
              else if(IS.eq.IB(J).and.ID(J).ne.2)then
                do K=NBBEG,NBEND
                  if(JS.eq.IB(K).and.ID(K).ne.2)then
                    call BNDREG(JB(J)+ISHF,JB(K)+ISHF,-1)
                  else if(JS.eq.JB(K).and.ID(K).ne.2)then
                    call BNDREG(JB(J)+ISHF,IB(K)+ISHF,-1)
                  end if
                end do
              end if
            end do
          end if
        end do
*----------------------------
 100    write(*,*)
     +  ' node ',TASKID,'  list of bound atoms for ',NAME(ITYP)
      end do   ! of ITYP      
 101  call  GATHER(NRNOD)
CM
      if(IPRINT.ge.8.and.TASKID.eq.MAST)then
	write(*,*)' list of bound atoms: '
	do I=1,NSTOT
	  write(*,'(i5,i3,2x,40i5)')I,NNBB(I),(INBB(I,J),J=1,NNBB(I))
	end do 
      end if
*----------------------------------------
      RETURN
      END   
*
*============ BNDREG ==========================================================
*
*   2.2 registration of bound sites
      subroutine BNDREG(IBB,JBB,ISIG)
      include "prcm.h"
      if(IBB.eq.JBB)return
      ITYP=ITS(IBB)
      NSP=NSITS(ITYP)
      if(ISIG.ge.0)then
        ISG=1
      else
        ISG=-1
      end if
      if(ITYP.ne.ITS(JBB))then
        write(*,*)'   Internal error: different types in BNDREG: ',
     +IBB,JBB,':',ITYP,ITS(JBB)
        stop
      end if
      ISHA = ISADDR(ITYP)-ISADR(ITYP)
      do IM=1,NSPEC(ITYP)
        IBS=IBB+ISHA+(IM-1)*NSP
        JBS=JBB+ISHA+(IM-1)*NSP
*  check that this pair is not already registered
        do I=1,NNBB(IBS)
          if(iabs(INBB(IBS,I)).eq.JBS)go to 10 
        end do
	NNBB(IBS)=NNBB(IBS)+1
	INBB(IBS,NNBB(IBS))=ISG*JBS
	      if(NNBB(IBS).gt.NBDMAX)then
	        write(*,*)
     +' list of bound atoms exceeded for atom ',IBS,' ',NM(IBS)
        write(*,*)' Possibly, non-bonded interaction set off ???' 
	        stop
	      end if
 10     do I=1,NNBB(JBS)
          if(iabs(INBB(JBS,I)).eq.IBS)go to 20 
        end do

	NNBB(JBS)=NNBB(JBS)+1
	INBB(JBS,NNBB(JBS))=ISG*IBS
	      if(NNBB(JBS).gt.NBDMAX)then
	        write(*,*)
     +' list of bound atoms exceeded for atom ',JBS,' ',NM(JBS)
        write(*,*)' Possibly, non-bonded interaction set off ???' 
	        stop
	      end if
 20     continue
      end do
      return
      end
*
*================ CUBIC =================================================
*                                    
*   3. Set molecules on a cubic lattice
*
      SUBROUTINE CUBIC(X,Y,Z,BOXL,BOYL,BOZL,NOP,ICELL,IPRINT)
      IMPLICIT real*8 (A-H,O-Z)
*
      DIMENSION X(*),Y(*),Z(*)
*
      if(NOP.eq.1)then
         X(1)=0.
         Y(1)=0.
         Z(1)=0.
         write(*,*)' Put molecule COM to 0 '
         return
      end if
      FNOP     = DFLOAT(NOP) 
	if(ICELL.eq.1)then
	  FAC=2.
	else
	  FAC=1.
	end if          
	SIZE	 = (BOXL*BOYL*BOZL/FNOP/FAC)**(1.d0/3.d0)
      NRX      = BOXL/SIZE
      NRY      = BOYL/SIZE
      NRZ      = BOZL/SIZE
 10 	if(IPRINT.ge.7)write(*,*)NRX,NRY,NRZ
	if(NRX*NRY*NRZ.lt.NOP)NRX=NRX+1
	if(NRX*NRY*NRZ.lt.NOP)NRY=NRY+1
	if(NRX*NRY*NRZ.lt.NOP)then
        NRZ=NRZ+1
	  go to 10
	end if
      NRCUBE   = NRX*NRY*NRZ
      NEMPTY   = NRCUBE-NOP
      L        = 0
      DO 100 I = 1,NRX
      DO 100 J = 1,NRY
      DO 100 K = 1,NRZ
      L        = L+1
      X(L)     = (DFLOAT(I)-0.5)*BOXL/NRX-0.5*BOXL
      Y(L)     = (DFLOAT(J)-0.5)*BOYL/NRY-0.5*BOYL
      Z(L)     = (DFLOAT(K)-0.5)*BOZL/NRZ-0.5*BOZL 
	if(ICELL.eq.1.and.(abs(X(L))+abs(Y(L))+abs(Z(L))).gt.0.75*BOXL)
     +   L=L-1
	if(L.ge.NOP)go to 110
  100 CONTINUE 
	if(L.lt.NOP)then
	  NRX=NRX+1
	  go to 10
	end if
  110 NSKIP    = 0
      IF(NEMPTY.NE.0) NSKIP=NOP/NEMPTY+1
      PRINT "(/1X,'*** CUBIC LATTICE: ',3I4,' SITES/EDGE  -> ',I4/
     X' SITES TOTALLY ( ',I3,' VACANT  - WITH INTERVALL OF ',I3,' )'/)"
     X,NRX,NRY,NRZ,NRCUBE,NEMPTY,NSKIP
*
      DO I = 1,NOP
        call PBC(X(I),Y(I),Z(I))
	end do
*                           
	call MIXTURE(X,Y,Z,NOP)
*
      RETURN
      END
*
*=============== FCC ===================================================
*                                   
*  4. Set molecules on the FCC lattice
*  -----------------------------------
      SUBROUTINE FCC(X,Y,Z,BOXL,BOYL,BOZL,NSP,ICELL)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION X(NSP),Y(NSP),Z(NSP)
      DIMENSION XB(4),YB(4),ZB(4) 
      DATA XB/0.5,-.5,0.5,-.5/,YB/0.5,-.5,-.5,0.5/,ZB/0.5,0.5,-.5,-.5/
C     SET UP STARTING FCC LATTICE 
      if(NSP.eq.1)then
         X(1)=0.
         Y(1)=0.
         Z(1)=0.
         write(*,*)' Put molecule COM to 0 '
         return
      end if
      HBOXL	= 0.5*BOXL 
      VOL=BOXL*BOYL*BOZL
      if(ICELL.eq.1)then
	FAC=2.
      else
	FAC=1.
      end if  
      SCELL=(4.*VOL/NSP)**0.333333333333
      NCX=BOXL/SCELL+1        
      NCY=BOYL/SCELL+1        
      NCZ=BOZL/SCELL+1 
 5    continue
      if(4*NCX*NCY*NCZ.lt.NSP)then
	NCX=NCX+1
        if(4*NCX*NCY*NCZ.lt.NSP)then
	  NCY=NCY+1
          if(4*NCX*NCY*NCZ.lt.NSP)then
	    NCZ=NCZ+1
	    go to 5
          end if
        end if
      end if        
 10   FACT   = HBOXL/DFLOAT(NCX)
      SHIFX  = BOXL/DFLOAT(NCX)
      SHIFY  = BOYL/DFLOAT(NCY)
      SHIFZ  = BOZL/DFLOAT(NCZ)
      BASE   = HBOXL*(-1.D0+1./DFLOAT(NCX))
      IL     = 0
      DO 20 IB = 1,4
        ZS       = BASE+ZB(IB)*FACT
        DO 19 IZ = 1,NCZ
          YS       = BASE+YB(IB)*FACT
          DO 18 IY = 1,NCY
            XS       = BASE+XB(IB)*FACT
            DO 17 IX = 1,NCX
              IL       = IL+1
              X(IL)    = XS
              Y(IL)    = YS
              Z(IL)    = ZS
		if(ICELL.eq.1.and.(abs(X(IL))+abs(Y(IL))+abs(Z(IL))).gt.
     +          0.75*BOXL)IL=IL-1
              if(IL.ge.NSP)go to 24
              XS       = XS+SHIFX
   17       CONTINUE
            YS       = YS+SHIFY
   18     CONTINUE
          ZS       = ZS+SHIFZ
   19   CONTINUE
   20 CONTINUE  
	if(IL.lt.NSP)then
	  IL=0
	  NCY=NCY+1
	  go to 10
	end if
*
   24 do I=1,NSP
        call PBC(X(I),Y(I),Z(I))
      end do
*
      write(*,*)'   FCC LATTICE CREATED FOR ',NSP,' POINTS'
*                           
      call MIXTURE(X,Y,Z,NSP)
*
      RETURN
      END
*
*=================== DNASOL ====================================
* 
*  5. Put one large molecule in a cylindrical hole along Z axis and
*  distribute other molecules around it
*
	subroutine DNASOL
	include "prcm.h"
	dimension XRS(NPART),YRS(NPART),ZRS(NPART)
	NOPSM=0
	do ITYP=1,NTYPES
	  if(IINIT(ITYP).ne.1)NOPSM=NOPSM+NSPEC(ITYP)
	end do
	VOLAV=VOL-BOZL*PI*RDNA**2
	if(TASKID.eq.MAST)
     +write(*,*)' Distribute ',NOPSM,' molecules in ',VOLAV,' A**3'
	SHAG=(VOLAV/NOPSM)**0.333333333 
	IAT=0 
	NSX=BOXL/SHAG+0.7
	NSY=BOYL/SHAG+0.7
	NSZ=BOZL/SHAG+0.5
 10	NCOUNT=0 
	do IX=1,NSX
	  XX=-HBOXL+BOXL*(dfloat(IX)-0.6)/NSX
	  do IY=1,NSY
	    YY=-HBOYL+BOYL*(dfloat(IY)-0.45)/NSY 
	    RR=sqrt(XX**2+YY**2)
	    if(RR.gt.RDNA)then
	      do IZ=1,NSZ
	        NCOUNT=NCOUNT+1
	        XRS(NCOUNT)=XX
	 	YRS(NCOUNT)=YY
	        ZRS(NCOUNT)=-HBOZL+BOZL*(dfloat(IZ)-0.45)/NSZ
	        if(NCOUNT.ge.NOPSM)go to 20
            end do
	    end if
	  end do
	end do
	if(NCOUNT.le.NOPSM)then
	  IAT=IAT+1
	  if(mod(IAT,3).eq.0)then
	    NSX=NSX+1
	    if(TASKID.eq.MAST)write(*,*)' increase NSX to ',NSX
	  else if(mod(IAT,3).eq.1)then
	    NSY=NSY+1
	    if(TASKID.eq.MAST)write(*,*)' increase NSY to ',NSY
	  else
	    NSZ=NSZ+1
	    if(TASKID.eq.MAST)write(*,*)' increase NSZ to ',NSZ
	  end if
	  go to 10
	end if
*  peremeshivanie
 20	call MIXTURE(XRS,YRS,ZRS,NOPSM)
	J=0
	JMOL=0
	do ITYP=1,NTYPES 
	  IBEG=IADDR(ITYP)+1
	  IEND=IADDR(ITYP+1)
	  if(IINIT(ITYP).eq.1)then
	    if(NSPEC(ITYP).gt.1)stop 
     +' only 1 molecule of such type is allowed'
	    JMOL=JMOL+1
	    X(JMOL)=0.d0
	    Y(JMOL)=0.d0
	    Z(JMOL)=0.d0
	    if(IPRINT.ge.8)write(*,*)
     +' put molecule ',JMOL,' at zero ref.point' 
	  else
	    do I=IBEG,IEND 
	      J=J+1
	      JMOL=JMOL+1
	      X(JMOL)=XRS(J)	
	      Y(JMOL)=YRS(J)	
	      Z(JMOL)=ZRS(J)	
	    end do
	  end if
	end do	      	          
 	if(TASKID.eq.MAST)write(*,*)JMOL,' molecules are distributed' 
	if(IPRINT.ge.8)then
	  WRITE(*,*)' COM of molecules:'
	  DO ityp=1,ntypes	
	    IBEG=IADDR(ITYP)+1
	    IEND=IADDR(ITYP+1)
          do I=IBEG,IEND
	       write(*,'(2i5,3f14.5)')ITYP,I,X(I),Y(I),Z(I)
	    end do
	  end do
	end if
	return
	end
*
*=================== SPHOLE ====================================
*
*  6. Put one large molecule in a spherical hole and
*  distribute other molecules around it
*
	subroutine SPHOLE
	include "prcm.h"
	dimension XRS(NPART),YRS(NPART),ZRS(NPART)
	NOPSM=0
	do ITYP=1,NTYPES
	  if(IINIT(ITYP).ne.1)NOPSM=NOPSM+NSPEC(ITYP)
	end do
	VOLAV=VOL-4.d0*PI*RDNA**3/3.d0
	if(TASKID.eq.MAST)
     +write(*,*)' Distribute ',NOPSM,' molecules in ',VOLAV,' A**3'
	SHAG=(VOLAV/NOPSM)**0.333333333 
	IAT=0 
	NSX=BOXL/SHAG+0.7
	NSY=BOYL/SHAG+0.7
	NSZ=BOZL/SHAG+0.5
 10	NCOUNT=0 
	do IX=1,NSX
	  XX=-HBOXL+BOXL*(dfloat(IX)-0.6)/NSX
	  do IY=1,NSY
	    YY=-HBOYL+BOYL*(dfloat(IY)-0.45)/NSY
	    do IZ=1,NSZ 
	      ZZ=-HBOZL+BOZL*(dfloat(IZ)-0.55)/NSZ
	      RR=sqrt(XX**2+YY**2+ZZ**2)
	      if(RR.gt.RDNA.and.(.not.LOCT.or.
     +abs(XX)+abs(YY)+abs(ZZ).lt.0.75*BOXL))then
	        NCOUNT=NCOUNT+1
	        XRS(NCOUNT)=XX
	 	YRS(NCOUNT)=YY
	        ZRS(NCOUNT)=ZZ
	        if(NCOUNT.ge.NOPSM)go to 20
            end if
	    end do
	  end do
	end do
	if(NCOUNT.le.NOPSM)then
	  IAT=IAT+1
	  if(mod(IAT,3).eq.0)then
	    NSX=NSX+1
	    if(TASKID.eq.MAST)write(*,*)' increase NSX to ',NSX
	  else if(mod(IAT,3).eq.1)then
	    NSY=NSY+1
	    if(TASKID.eq.MAST)write(*,*)' increase NSY to ',NSY
	  else
	    NSZ=NSZ+1
	    if(TASKID.eq.MAST)write(*,*)' increase NSZ to ',NSZ
	  end if
	  go to 10
	end if
*  blandning
 20	call MIXTURE(XRS,YRS,ZRS,NOPSM)
*  put fixed molecules on their places
	J=0
	JMOL=0
	do ITYP=1,NTYPES 
	  IBEG=IADDR(ITYP)+1
	  IEND=IADDR(ITYP+1)
	  if(IINIT(ITYP).eq.1)then
	    if(NSPEC(ITYP).gt.1)stop 
     +' only 1 molecule of such type is allowed'
	    JMOL=JMOL+1
	    X(JMOL)=0.d0
	    Y(JMOL)=0.d0
	    Z(JMOL)=0.d0
	    if(IPRINT.ge.8)write(*,*)
     +' put molecule ',JMOL,' at zero ref.point' 
	  else
	    do I=IBEG,IEND 
	      J=J+1
	      JMOL=JMOL+1
	      X(JMOL)=XRS(J)	
	      Y(JMOL)=YRS(J)	
	      Z(JMOL)=ZRS(J)	
	    end do
	  end if
	end do	      	          
 	if(TASKID.eq.MAST)write(*,*)JMOL,' molecules are distributed' 
	if(IPRINT.ge.8)then
	  WRITE(*,*)' COM of molecules:'
	  DO ityp=1,ntypes	
	    IBEG=IADDR(ITYP)+1
	    IEND=IADDR(ITYP+1)
          do I=IBEG,IEND
	       write(*,'(2i5,3f14.5)')ITYP,I,X(I),Y(I),Z(I)
	    end do
	  end do
	end if
	return
	end     
*
*===================== MIXTURE ===================================
* 
*  7. Mix molecules to avoid artificial correlations
*
	subroutine MIXTURE(X,Y,Z,NOP)
	real*8 X(NOP),Y(NOP),Z(NOP)
*
*  making mixture
*
 	do ICT=2,7
	  if(mod(ICT,2).eq.0)then
	    J=NOP
	    do I=1,NOP,ICT
	      XX=X(I)
	      YY=Y(I)
	      ZZ=Z(I)
	      X(I)=X(J)
	      Y(I)=Y(J)
	      Z(I)=Z(J)
	      X(J)=XX
	      Y(J)=YY
	      Z(J)=ZZ
	      J=J-1
	    end do
	  else
	    J=1
	    do I=1,NOP,ICT
	      XX=X(I)
	      YY=Y(I)
	      ZZ=Z(I)
	      X(I)=X(J)
	      Y(I)=Y(J)
	      Z(I)=Z(J)
	      X(J)=XX
	      Y(J)=YY
	      Z(J)=ZZ
	      J=J+1
	    end do
	  end if
      end do
	return
	end
*
*================ PULLINI =======================================
*
      subroutine PULLINI
C
C   8. This subroutine read coordinates of points to which some
C   atoms will be bound by harmonical forces (see subr PULLBACK)
C   Normally, it is not used (LRR=.false.)
C
      include "prcm.h"
      character*128 STR,TAKESTR  
      open(unit=31,file=filref,status='old',err=98)
      STR=TAKESTR(31,IE)
      read(STR,*,end=99,err=99)NPULL
      if(NPULL.gt.NPLM)stop '!!! increase NPLM!!!'
      if(IPRINT.ge.6)write(*,*)' Linked atoms: '
      do I=1,NPULL
C  IS - atom number
        STR=TAKESTR(31,IE) 
        read(STR,*,end=99,err=99)IS,XX,YY,ZZ
        INPL(I)=IS
        XPL(I)=XX
        YPL(I)=YY
        ZPL(I)=ZZ
        if(IPRINT.ge.6)write(*,'(i5,2x,a,3f14.6)')
     +  IS,NM(NSITE(IS)),XX,YY,ZZ
      end do
      close(31)
      return
 98   if(TASKID.eq.MAST)write(*,*)
     +' LRR=.true. and file ',filref,' not found'
      call FINAL
 99   IE=99
      STR=TAKESTR(31,IE)
      end
