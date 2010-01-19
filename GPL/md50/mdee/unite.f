*=============== IUNITS ===============================================
*
      SUBROUTINE IUNITS
*
      include "mdee.h"
*
      DIMENSION  DM(3),QM(3,3),OM(3,3,3),HM(3,3,3,3)
     X          ,DI(3),QI(3,3),OI(3,3,3),HI(3,3,3,3)
*
*Calculate box dimensions etc:
*
      FACTM       = AVSNO*1000.0
      TOTMAS      = TOTMAS/FACTM
      AVOFNR      = AVSNO/FNOP
      IF(RHO.NE.0.D0) THEN
        VOLM        = AVOFNR*(TOTMAS/(RHO*1000.0))
        IF(BOXL.EQ.0.D0) THEN
          VOL      = VOLM*FNOP*1.d30/AVSNO
          CL        = VOL**(1.0/3.0)
        else
          CL=BOXL
        ENDIF
      ELSE
*Vacuum simulation - arbitrary box size:
        CL        = 1000.D0
        PRINT "(/,'*** VACUUM SIMULATION ***',/)"
      END IF
      call RECLEN(CL)
       VOLM        = 1.d-30*VOL*AVOFNR 
	VOLF	= VOL  
	BOXF	= 0.5*CL
*
      
*Generate internal units:
*      
      UNITL       = 1.d-10
      UNITM       = TOTMAS
      UNITT       = DT
      UNITE       = UNITM*(UNITL/UNITT)**2
      UNITP       = UNITE*1.e-5/UNITL**3
*
*Define dimensionless quantities and conversion factors
*
      TSTEP       = UNITT/SQRT((UNITM/UNITE)*UNITL**2)
      TSTEB       = TSTEP*FTIMES
      ENERF       = AVSNO*UNITE
      TEMPF       = 2.0/(3.0*AVSNO*BOLTZ)
      BETA	    = UNITE/(TRTEMP*BOLTZ)
      COULF       = (ELECHG/EPS0)*(ELECHG/UNITE)/(4.D0*PI*UNITL)
      ELFAC	    = COULF*BETA
      PERMOL      = ENERF*1.D-3/FNOP 
      SHORT2	    = SHORT**2
      SHORTA	    = SHORT-1.d-11/UNITL   !  -0.1A
*
	if(IPRINT.ge.7)then
        PRINT "(/'*** UNIT MASS         ',D12.6,' kg      ' )",UNITM
        PRINT "( '*** UNIT TIME         ',D12.6,' S       ' )",UNITT
        PRINT "( '*** UNIT ENERGY       ',D12.6,' J       ' )",UNITE
        PRINT "( '*** UNIT LENGTH       ',D12.6,' m       ' )",UNITL
        PRINT"('*** UNIT PRESSURE     ',D12.6,' 10**5 N/m**2  '/)",UNITP
        PRINT "('*** LONG  TIME STEP   ',F12.6,' .I.U.   ' )",TSTEP
        PRINT "( '*** SHORT TIME STEP   ',F12.6,' .I.U.   ' )",TSTEB
        PRINT "( '*** ENERGY  FACTOR    ',D12.6,' J/mol   ' )",ENERF
        PRINT "( '*** TEMP    FACTOR    ',D12.6,' K/J     ' )",TEMPF
        PRINT "( '*** COULOMB FACTOR    ',D12.6,'         '/)",COULF
	end if
*
*dimensionless mass parameters:
*
      DO IS      =  1,NSITES
      MASSD (IS)  = MASS (IS)/FACTM/UNITM
      END DO! OF I
*
*Specific data for molecule types:
*
      IADDR (1)      = 0
      ISADR (1)      = 0
      ISADDR(1)      = 0
      DO ITYP        = 1,NTYPES
      IADDR (ITYP+1) = IADDR (ITYP)+NSPEC(ITYP)
      ISADR (ITYP+1) = ISADR (ITYP)+NSITS(ITYP)
      ISADDR(ITYP+1) = ISADDR(ITYP)+NSITS(ITYP)*NSPEC(ITYP)
      END DO! OF ITYP
*     
	if(IPRINT.ge.9)then
      PRINT "(25X,'*** ADDRESSES: ***'/)"
      DO ITYP        = 1,NTYPES+1
      PRINT "(10X,4(1X,I8))",ITYP,ISADR(ITYP),IADDR(ITYP),ISADDR(ITYP)
      END DO! OF ITYP
	end if
*
*** 	call RECINT
      N              = 0
      MT             = 0
      DO ITYP        = 1,NTYPES
      QS             = 0.D0
      if(IPRINT.ge.10)
     +PRINT "(/'*** POTENTIAL MODEL FOR :',A6)",NAME(ITYP)
      NSB            =  ISADR (ITYP)+1
      NSE            =  ISADR (ITYP+1)
      IF(IPRINT.GT.10) PRINT "(/,' SOME SITE DATA: ')"
      DO I           = 1,NSPEC(ITYP)
      DO IS          = NSB,NSE
      N              = N+1
      MASSDI(N)      = 1.0/MASSD(IS)
      Q     (N)      =    CHARGE(IS)
      NAM   (N)      = NM(IS)
      QS             = QS+DABS(Q(N))
      IF(IPRINT.GT.10) 
     XPRINT "(I4,3X,I4,3X,A4,6X,'Q: ',F6.3,5X,'MI:  ',F12.6)",
     XN,IS,NAM(N),Q(N),MASSDI(N)
      END DO! OF I
      END DO! OF IS
*
      QSUM(ITYP)     = QS
*
      IF(IPRINT.GT.9) THEN
      PRINT *
      DO IS          = NSB,NSE
      PRINT "('# SITE:',I3,':  ',A4
     X,' EPS:',F7.3,'  SIG:',F7.3,'  CHARGE:',F7.3,'  MASS:',F7.3)"
     X,IS,NM(IS),EPSIL(IS),SIGMA(IS),CHARGE(IS),MASS(IS)
      END DO! OF IS
      END IF!(IPR...
*
      NNN            = NS
      CALL MPOLES(R,CHARGE,DM,QM,OM,HM,NSB,NSE,NNN)
*
      DO I		= 1,3
      DI(I)		= DM(I)*2.5418
      DO J		= 1,3
      QI(I,J)		= QM(I,J)*1.345
      DO K		= 1,3
      OI(I,J,K)		= OM(I,J,K)*0.7118
      DO L		= 1,3
      HI(I,J,K,L)	= HM(I,J,K,L)*0.3767
      END DO! OF L
      END DO! OF K
      END DO! OF J
      END DO! OF I
*	
      DII		= DSQRT(DI(1)**2+DI(2)**2+DI(3)**2)
      DMM		= DSQRT(DM(1)**2+DM(2)**2+DM(3)**2)
*
      IF(DABS(DMM).GT.0.0001D0) THEN
      if(IPRINT.ge.9)WRITE(6,160) DMM,DII
  160 FORMAT(/1X,'*** DIPOLE MOMENT    (MAGNITUDE)   ',F10.5,
     +' A.U.  -> ',F10.5,'  10**-18 esu '/)
      END IF
*
      DO JTYP        = ITYP,NTYPES
      MT             = MT+1
      MDX(ITYP,JTYP) = MT
      MDX(JTYP,ITYP) = MT
      POTES(MT)      = 0.D0
      POTLJ(MT)      = 0.D0
      END DO! OF JTYP
      END DO! OF ITYP
      MOLINT         = MT
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
*
*	CONVERSION TO INTERNAL UNITS:
*
*	Bond stretching 
*
*Harmonic potentials
*
      FACTOR     = 1.D3/AVSNO/UNITE
      EFACT      = FACTOR
      EFACB      = EFACT
      DO I       = 1,NB
        FB(I)      = FB(I)*EFACB
      END DO! OF I
*
* Morse potentials:
*
      DO I       = 1,NB
        DB(I)      = DB(I)*FACTOR
      END DO! OF I
*                  
*	Angle bending:
*
      DO I         = 1,NA
        FA(I)        = FA(I)*EFACT
        RA(I)        = RA(I)/TODGR
      END DO! OF I
*
*	Torsional angles
*
      DO  I        = 1,NT
        FT (I)       = FT (I)*EFACT
        RT (I)       = RT (I)/TODGR
        FT1(I)       = FT1(I)*EFACT
        FT2(I)       = FT2(I)*EFACT
        FT3(I)       = FT3(I)*EFACT
      END DO! OF I
*
*
*	Improper torsional angles
*
      DO I         = 1,NI
      FI(I)        = FI(I)*EFACT
      RI(I)        = RI(I)/TODGR
      END DO! OF I
*
      RETURN
      END
*
*=============== RECLEN ===============================================
*
      SUBROUTINE RECLEN(CL)
*
      include "mdee.h"
*
      if(CL.ne.0.d0)then
        BOXL=CL
        BOYL=BOXL
        BOZL=BOXL
      end if
      HBOXL=0.5*BOXL
      HBOYL=0.5*BOYL
      HBOZL=0.5*BOZL
      RF          = HBOXL
      RCUTT=RCUT
      IF(RCUT.LE.0.D0)RCUTT = HBOXL
      RSSQ        = RCUTT**2
      BOXLI       =   1.0/BOXL
      BOYLI       =   1.0/BOYL
      BOZLI       =   1.0/BOZL
      VOL	= BOXL*BOYL*BOZL
      ALPHAD	= ALPHA/RCUTT
      RHO	= TOTMAS*1.d27/VOL
      return
      end
