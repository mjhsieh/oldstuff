*     Program MDynaMix,  v.4.3
*
*     Part 5
*
*     File forces.f
*     -------------
*
C     This file contains subroutines responsible for calculations 
C     of different forces. Forces are divided into two groups:
C     slow and fast forces. This division is essential only if
C     double time step algorithm is used
C     
C     Subroutines:
C
C     1. FFORCES - collect fast forces  
C     2. SLFORCES - collect slow forces (double time step)
C     3. ALLFORCE - collect all forces (used in simple verlet)
C
C     4. GETBND - calculate covalent bonds forces
C     5. GETANG - calculate covalent angle forces
C     6. GETTRS - calculate forces due to dihedral angles (all types)
C     7. STRBND - a special procedure for SPC water model
C       
C     8. LOCALF - non-bonded forces of closest neighbours
C     9. FORCES - non-bonded forces of far neighbours
C     10. FURIR - reciprical space Ewald contribution
C     11. INIFUR - initialisation of FURIR subroutine
C     12. ETERMS - self-interaction term of Ewald sum and intramolecular corrections
C     13. ELRCLJ - long-range corrections from LJ forces
C     14. PULLBACK - fixed harmonic potential for selected atoms

C     15. CHCNB - update list of neighbours 
*
*  1. Collect fast forces
*  ----------------------
      subroutine FFORCES(LAVER)
      logical LAVER,DONE
*  1.1   Set zero-s 
      CALL ZEROFS(1)  ! service.f
*  1.2   Calculate different contributions to the fast forces
C        (see details in descriptions of corresponding subroutines below)   
      call LOCALF    
      CALL GETBND  
      CALL GETANG  
      CALL GETTRS  
*  1.3   Sum up fast forces from different nodes
      call CM_ADDF(LAVER)  ! mpi.f or scalar.f
*  1.4   Cut large forces of specified
      call SCLFRC(DONE,1)
      return
      end
* 
*  2. Collect slow forces
*  ----------------------
      subroutine SLFORCES
      logical DONE
*  2.1   Set zero-s 
      call ZEROAV     ! service.f
      CALL ZEROFS(2)  ! service.f
*  2.2   Calculate different contributions to the slow forces
C        (see details in descriptions of corresponding subroutines below)   
      CALL FURIR      ! 
      CALL ETERMS     ! 
      CALL FORCES     ! l-forces.f   
*   2.3  Sum up "slow" forces from different nodes
      CALL CM_ADDLF    ! mpi.f
C   This may take essential time but seems has no effect
*      if(LNPT)call ELRCLJ     ! l-forces.f 
C   If a molecule is softly fixed
      call PULLBACK
*   2.4   Cut large forces of specified
      call SCLFRC(DONE,2)
      return
      end
*
*  3. Collect all forces
*  ----------------------
      subroutine ALLFORCE
      logical DONE
*  3.1   Set zero-s 
      call ZEROAV     ! service.f
      CALL ZEROFS(1)  ! service.f
      CALL ZEROFS(2)  ! service.f
*  3.2   Calculate different contributions to the forces
C        (see details in descriptions of corresponding subroutines below)   
      call LOCALF    
      CALL GETBND  
      CALL GETANG  
      CALL GETTRS  
      CALL FURIR    
      CALL ETERMS   
      CALL FORCES   
*   3.3  Sum up all forces from different nodes
      CALL CM_ADDLA    ! mpi.f
C   This may take essential time but seems has no effect
*      if(LNPT)call ELRCLJ     ! l-forces.f 
C   If a molecule is softly fixed
      call PULLBACK
*   3.4   Cut large forces of specified
      call SCLFRC(DONE,3)
      return
      end
*
*=============== INTRAMOLECUULAR FORCES ==========================
*
*    4. Covalent bond contribution
*    -----------------------------
C
C    This subroutine calculate contribution to forces from 
C    covalent bonds. Each node has own list of bonds. Distribution
C    of bonds overthe  nodes takes place in subroutine UNITS (setup.f)
C 
C=============== GETBND ===========================================
C
      SUBROUTINE GETBND
*
*   4.1 Front matter
*   ----------------
      include "prcm.h"
      timeb0	= cputime(0.d0) 
C  This set zeroth in arrays calculation average bond
      DO M      = 1,NRBB
        BR(M)     = 0.D0
        BE(M)     = 0.D0
      END DO
C  This call special procedure to calculate bonds in the case if 
C  parameter IPOT for a molecule is 1 or 2
      do ITYP=1,NTYPES
	if(IPOT(ITYP).eq.1.or.IPOT(ITYP).eq.2)call STRBND(ITYP)
      end do
      if(NRBB.le.0)go to 200        
*
*   4.2 Organize cycle over bonds and define parameters of each bond
*   ----------------------------------------------------------------
C  Cycle over the bonds for this node
      do N=1,NRBT
        I     = IBI(N)              !  first atom
	ITYP  = ITYPE(I)            !  type
	M     = IBUK(N)             !  global bond number
        J     = JBJ(N)              !  second atom
	NSP   = NSPEC(ITYP)         !  num. of molecules of this type       
        FNSPI     = 1.D0/DFLOAT(NSP)
C   Equilibrium length of the bond
        REQ       = RB(M)
*
*    4.3 Calculate bond length
*    -------------------------
        BX        = SX(J)-SX(I)
        BY        = SY(J)-SY(I)
        BZ        = SZ(J)-SZ(I)
        call PBC(BX,BY,BZ)
        BSQ       = BX*BX+BY*BY+BZ*BZ
        BD        = DSQRT(BSQ)
        R1I       = 1.D0/BD
C    Deviation fron equilibrium
        BBB       =  BD-REQ
*
*    4.4 Calculate energy and force
*    ------------------------------
*    4.4.1 Case of harmonic bond
	if(ID(M).ne.1.or.BBB.gt.REQ)then  
C    Potential parameter
          FCN       = FB(M)
          DBND      = BBB*FCN
C    Energy of the bond
          EBDD	    = DBND*BBB
C    Force (abs. value)
          FORB      =  -2.0D0*DBND*R1I
*    4.4.2 Case of Morse potential 
	else  
	  DDIS = DB(M)
	  ROH  = RO(M)
	  EXM  = exp(-ROH*BBB)
          EXMP = 1-EXM
C    Energy of the bond
	  EBDD = DDIS*EXMP**2
C    Force (abs. value)
	  FORB = -2.d0*DDIS*ROH*EXM*EXMP*R1I       
	end if
*
*     4.5 Add results to accumullating arrays
*     ---------------------------------------
        BR(M)     = BR(M)+BD*FNSPI            ! average length
        BE(M)     = BE(M)+EBDD*FNSPI          ! average energy
        PINT	  = PINT+EBDD                 ! total intramolec. energy
        if(LMOVE(ITYP))then
C    Forces
          HX(I)     = HX(I)-BX*FORB
          HY(I)     = HY(I)-BY*FORB
          HZ(I)     = HZ(I)-BZ*FORB
          HX(J)     = HX(J)+BX*FORB
          HY(J)     = HY(J)+BY*FORB
          HZ(J)     = HZ(J)+BZ*FORB
C    Contribution to the virial
          VIRB	    = VIRB+FORB*BSQ
	  VIRFX	    = VIRFX+FORB*BX**2  
	  VIRFY	    = VIRFY+FORB*BY**2  
	  VIRFZ	    = VIRFZ+FORB*BZ**2
        end if
	if(IPRINT.ge.9)write(*,*)I,J,BD,EBDD/EFACT
*       
 100	continue
      END DO! OF M       
*
 200  timeb	= timeb+cputime(timeb0)
      RETURN
      END
*
*    5. Covalent angle contribution
*    -------------------------------
C
C    This subroutine calculate contribution to forces from 
C    covalent angles. Each node has own list of angles. Distribution
C    of angles over the nodes takes place in subroutine UNITS (setup.f)
C 
C
C=============== GETANG ===========================================
C
      SUBROUTINE GETANG
*
*   5.1 Front matter
*   ----------------
      include "prcm.h"
      IF(NRAA.LE.0) RETURN
      timea0	= cputime(0.d0)
*     
      DO M      = 1,NRAA
        AR(M)     = 0.D0
        AE(M)     = 0.D0
      END DO
*
*   5.2 Organize cycle over the angles
*   ----------------------------------
      do N=1,NRAT
	I=IAI(N)                            ! first atom
	ITYP=ITYPE(I)                       ! type
        J         = JAJ(N)                  ! second atom
        K	  = KAK(N)                  ! third atom
	M	  = IAUK(N)                 ! global angle number
        FCN       = FA(M)                   ! force constant
        AEQ       = RA(M)                   ! equilibrium angle
        NSP       = NSPEC (ITYP)       
        FNSPI     = 1.D0/DFLOAT(NSP)
*
*   5.3 Calculate the angle and energy 
*   ----------------------------------
*   5.3.1  I->J  vector
          AX        = SX(I)-SX(J)
          AY        = SY(I)-SY(J)
          AZ        = SZ(I)-SZ(J)
          call PBC(AX,AY,AZ)
          ASQ       = AX*AX+AY*AY+AZ*AZ
          ASQI      = 1.D0/ASQ
*   5.3.2  K->J  vector
          BX        = SX(K)-SX(J)
          BY        = SY(K)-SY(J)
          BZ        = SZ(K)-SZ(J)
          call PBC(BX,BY,BZ)
          BSQ       = BX*BX+BY*BY+BZ*BZ
          BSQI      = 1.D0/BSQ
*   5.3.3 Calculate Angle
          AB        = DSQRT(ASQ*BSQ)
          ABI       = 1.D0/AB
          COSA      = AX*BX+AY*BY+AZ*BZ
          COSA      = COSA*ABI
          COSA      = DMIN1(COSA, 1.D0)
          COSA      = DMAX1(COSA,-1.D0)
          ALF       = DACOS(COSA)
          SINA      = DSIN(ALF)
          SINAI     =  1.D0/SINA
*   5.3.4 Calculate energy
          DALF      =  ALF-AEQ          ! deviation from equilibrium
          ABC       = DALF*FCN          ! derivative of the energy
          EANG	    = ABC*DALF	        ! energy
          AR(M)     = AR(M)+ALF*FNSPI
          AE(M)     = AE(M)+EANG*FNSPI 
          PINT	    = PINT+EANG 
	  if(IPRINT.ge.9)write(*,*)I,J,K,ALF*TODGR,EANG/EFACT
*
*   5.4 Force calculation
*   ---------------------
          if(LMOVE(ITYP))then
            SAB       =-2.D0*ABC*SINAI
            CAB       = SAB*COSA
            FAB       = SAB*ABI
            FAA       = CAB*ASQI
            FBB       = CAB*BSQI
            FAX       = FAB*BX-FAA*AX
            FAY       = FAB*BY-FAA*AY
            FAZ       = FAB*BZ-FAA*AZ
            FBX       = FAB*AX-FBB*BX
            FBY       = FAB*AY-FBB*BY
            FBZ       = FAB*AZ-FBB*BZ
            HX(I)     = HX(I)-FAX
            HY(I)     = HY(I)-FAY
            HZ(I)     = HZ(I)-FAZ
            HX(J)     = HX(J)+FAX+FBX
            HY(J)     = HY(J)+FAY+FBY
            HZ(J)     = HZ(J)+FAZ+FBZ
            HX(K)     = HX(K)-FBX
            HY(K)     = HY(K)-FBY
            HZ(K)     = HZ(K)-FBZ     
C   Contributions to virials (only to projections!)
	    VIRFX      = VIRFX-AX*FAX-BX*FBX
	    VIRFY      = VIRFY-AY*FAY-BY*FBY
	    VIRFZ      = VIRFZ-AZ*FAZ-BZ*FBZ
          END IF
      END DO! OF N
CD	write(*,*)TASKID,PINT*0.001*ENERF/FNOP,PINT*0.001*ENERF/FNOP
      timea	= timea+cputime(timea0)
      RETURN
      END
*
*    6. Dihedral angle contribution
*    -------------------------------
C
C    This subroutine calculate contribution to forces from 
C    dihedral angles. Each node has own list of angles. Distribution
C    of angles over the nodes takes place in subroutine UNITS (setup.f)
C
C    This subroutine was initially written by Andrei Komolkin 
C
C
C=============== GETTRS ============================================== 
C
      SUBROUTINE GETTRS
*
*   6.1 front matter
*   ----------------
      include "prcm.h"
      IF(NRTT.LE.0) RETURN
      timet0	= cputime(0.d0)
      DO M      = 1,NRTT
        TR(M)     = 0.D0
        TE(M)     = 0.D0
      END DO
*
*   6.2 Organize cycle over dihedral angles
*   ---------------------------------------              
      do N=1,NRRT
	I	= ITI(N)           !  first atom (I)
	ITYP	= ITYPE(I)         !  type
        J         = JTJ(N)         !  second atom (J)
        K         = KTK(N)         !  third atom (K)
        L         = LTL(N)         !  forth atom (L)
	M	    = ITUK(N)      !  global dihedral angle number
        NSP       = NSPEC (ITYP)
        FNSPI     = 1.D0/DFLOAT(NSP)
*
*  6.3 Calculate dihedral angle
*  ----------------------------
*  6.3.1 Calculate vectors IJ,JK,KL
        ax        = sx(j)-sx(i)
        ay        = sy(j)-sy(i)
        az        = sz(j)-sz(i)
        bx        = sx(k)-sx(j)
        by        = sy(k)-sy(j)
        bz        = sz(k)-sz(j)
        cx        = sx(l)-sx(k)
        cy        = sy(l)-sy(k)
        cz        = sz(l)-sz(k)
        call PBC(AX,AY,AZ)
        call PBC(BX,BY,BZ)
        call PBC(CX,CY,CZ)
*  6.3.2 Scalar products
        ab        = ax*bx + ay*by + az*bz
        bc        = bx*cx + by*cy + bz*cz
        ac        = ax*cx + ay*cy + az*cz
        at        = ax*ax + ay*ay + az*az
        bt        = bx*bx + by*by + bz*bz
        ct        = cx*cx + cy*cy + cz*cz
*  6.3.3 Vector products
        axb       = (at*bt)-(ab*ab)
        bxc       = (bt*ct)-(bc*bc)
        fnum      = (ab*bc)-(ac*bt)
        den       = axb*bxc
*  6.3.4 Check that any three atoms is not on one line
C       (Otherwise contribution is zero)
        if(den.gt.1.0d-12) then 
          den       = dsqrt(den)         
C   Cosine of angle:
          co        = fnum/den
          CO        = DMIN1(CO, 1.D0)
          CO        = DMAX1(CO,-1.D0)
C    Sign of angle:
      signum    = ax*(by*cz-cy*bz)+ay*(bz*cx-cz*bx)+az*(bx*cy-cx*by)
*   6.3.5   Define angle:
          arg       = dsign(dacos(co),signum)
          SIN1F        = dsin(arg)
	  si	= SIN1F
          if( dabs(si).lt.1.0d-12 ) si = dsign(1.0d-12,si)
*
*   6.4 Choice from several dihedral angle types: 
*   ---------------------------------------------
*   6.4.1  Normal AMBER-type torsion angle 
          if(ITORS(M).eq.0)then
            FCN       = FT(M)
            TEQ       = RT(M)
            MUL       = NMUL(M)
            EARG      = MUL*ARG-TEQ
            POTT      = FCN*(1.D0+DCOS(EARG))  
            DERI      = -FCN*MUL*DSIN(EARG)
*   6.4.2 MM3 force field torsional angle
	  else if(ITORS(M).eq.1)then
            FN1   = FT1(M)
            FN2   = FT2(M)
            FN3   = FT3(M)  
	    COS1F = dcos(ARG)
            COS2F =  2.D0*COS1F**2-1.D0
            COS3F = COS1F*(2.D0*COS2F-1.D0)
	    SIN2F = 2.d0*SIN1F*COS1F
	    SIN3F = SIN1F*COS2F+SIN2F*COS1F
            POTT  = FN1 * (1.D0 + COS1F)
     X             +FN2 * (1.D0 - COS2F)
     X             +FN3 * (1.D0 + COS3F)  
	    DERI	= -FN1*SIN1F+2.d0*FN2*SIN2F-3.d0*FN3*SIN3F
*   6.4.3            torsional angle 
	  else if(ITORS(M).eq.5)then  
	    FN0   = FT(M)
            FN1   = FT1(M)
            FN2   = FT2(M)
            FN3   = FT3(M)  
            FN4   = FT4(M)
            FN5   = FT5(M)
            TEQ   = RT(M)
            EARG  = ARG-TEQ
            COSI      =  DCOS(EARG)
            COSI2     =  COSI*COSI
            COSI3     =  COSI*COSI2
            COSI4     =  COSI*COSI3
            COSI5     =  COSI*COSI4
	    SINI	=  dsin(EARG)
            POTT =    FN0+ FN1*COSI+FN2*COSI2+     FN3*COSI3+
     +                              FN4*COSI4+     FN5*COSI5 
	    DERI =  -SINI*(FN1+2.d0*FN2*COSI +3.d0*FN3*COSI2+
     +                         4.d0*FN4*COSI3+5.d0*FN5*COSI4)
*   6.4.4 Improper dihedral angle
	  else if(ITORS(M).eq.-1)then                        
            FCN       = FT(M)
            TEQ       = RT(M)
            EARG      = ARG-TEQ
            POTT      = FCN*EARG**2  
            DERI      = 2.d0*FCN*EARG
	  else
	    write(*,*)'  Torsion type ',ITORS(M),' not supported'
	    write(*,*)'  Torsion ',I,J,K,L,' on ',NAME(ITYP)
	    call FINAL
          end if
          TR(M)     = TR(M)+abs(ARG*FNSPI)
          TE(M)     = TE(M)+POTT*FNSPI
          PINT	= PINT+POTT
	  if(IPRINT.ge.9)write(*,*)I,J,K,L,ARG*TODGR,POTT/EFACT
*
*   6.5 Calculate Forces:
*   ---------------------
          if(LMOVE(ITYP))then
          de1       = DERI/den/si
          axb       = axb/den*co
          bxc       = bxc/den*co
*   6.5.1 X-components
          dnum      =  cx*bt - bx*bc
          dden      = (ab*bx - ax*bt)*bxc
          FFI1      = (dnum  - dden) * de1
          dnum      = ((bx-ax)*bc - ab*cx ) + (2.0*ac*bx - cx*bt)
          dden      = axb*(bc*cx-bx*ct) + (ax*bt-at*bx-ab*(bx-ax))*bxc
          FFJ1      = (dnum - dden) * de1
          dnum      = ab*bx - ax*bt
          dden      = axb*( bt*cx - bc*bx )
          FFL1      = (dnum - dden) * de1
          FFK1      = -(ffi1+ffj1+ffl1)
C   Forces:
          HX(I)     = HX(I)+ffi1
          HX(J)     = HX(J)+ffj1
          HX(K)     = HX(K)+ffk1
          HX(L)     = HX(L)+ffl1 
*   6.5.2 Y-components
          dnum      =  cy*bt - by*bc
          dden      = (ab*by - ay*bt)*bxc
          FFI2      = (dnum  - dden) * de1
          dnum      = ((by-ay)*bc - ab*cy ) + (2.0*ac*by - cy*bt)
          dden      = axb*(bc*cy-by*ct) + (ay*bt-at*by-ab*(by-ay))*bxc
          FFJ2      = (dnum - dden) * de1
          dnum      = ab*by - ay*bt
          dden      = axb*( bt*cy - bc*by )
          FFL2      = (dnum - dden) * de1
          FFK2      = -(ffi2+ffj2+ffl2)
          HY(I)     = HY(I)+ffi2
          HY(J)     = HY(J)+ffj2
          HY(K)     = HY(K)+ffk2
          HY(L)     = HY(L)+ffl2
*    6.5.3  Z-components
          dnum      =  cz*bt - bz*bc
          dden      = (ab*bz - az*bt)*bxc
          FFI3      = (dnum  - dden) * de1
          dnum      = ((bz-az)*bc - ab*cz ) + (2.0*ac*bz - cz*bt)
          dden      = axb*(bc*cz-bz*ct) + (az*bt-at*bz-ab*(bz-az))*bxc
          FFJ3      = (dnum - dden) * de1
          dnum      = ab*bz - az*bt
          dden      = axb*( bt*cz - bc*bz )
          FFL3      = (dnum - dden) * de1
          FFK3      = -(ffi3+ffj3+ffl3)
          HZ(I)     = HZ(I)+ffi3
          HZ(J)     = HZ(J)+ffj3
          HZ(K)     = HZ(K)+ffk3
          HZ(L)     = HZ(L)+ffl3
*   6.5.4 Contributions to projections of virial
	  VIRFX  = VIRFX-(AX+BX)*FFI1-BX*FFJ1+CX*FFL1
	  VIRFY  = VIRFY-(AY+BY)*FFI2-BY*FFJ2+CY*FFL2
	  VIRFZ  = VIRFZ-(AZ+BZ)*FFI3-BZ*FFJ3+CZ*FFL3
          end if  !  if(LMOVE
	end if    !  den > 0
*              
      END DO! OF N
      timet	= timet+cputime(timet0)
      RETURN
      END
*
*    7. Special procedure for SPC water model
*    ----------------------------------------
C
C    This subroutine calculate intramolecular forces in 
C    flexible SPC water model
C    Depending on parameter IPOT for the water, it employ 
C    harmonic intramolecular potential (IPOT=1) or
C    anharmonic Morse potential (IPOT=2)
C    In the both case potential bond strength parameter defined 
C    in .mol file have no meaning 
C    (this is mostly because the model includes some cross-interaction 
C    terms not defined in the AMBER-like potential)
C
C=============== STRBND ===========================================
C
      SUBROUTINE STRBND(ITYP)
*
*    7.1  Front matter
*
      include "prcm.h"
C   Parameters of intramolecular SPC potential:
      PARAMETER (AAP  = 2809.56D0 , BBP = 2809.56D0)
      PARAMETER (CCP  =  687.41D0 , DDP =   2.566D0)
      PARAMETER (EEP  = -884.63D0 , FFP = -884.63D0)
C   BBSA : if OH bond for water exceed BBSA, warning message is issued
      PARAMETER (GGP  =  467.31D0 ,BBMA=1.3,BBSA=1.4  )
C   Check that it is really water molecule
	if(NAME(ITYP)(1:3).ne.'H2O'.and.NAME(ITYP)(1:3).ne.'h2o')then
	  write(*,*)
     +'Potential of type 1 or 2 can be used only for water: ',
     +' mol. name H2O*'
	  call FINAL
	end if
* 
*   7.2 Prepare parameters for the calculations 
*   ------------------------------------------
C  bond number for waters
      I1       = IADB(ITYP)
      I2       = I1+1
      I3       = I2+1
C  Initial and last water molecules
      IBEG     = IADDR(ITYP)+1
      IEND     = IADDR(ITYP+1)
      NSP      = NSPEC(ITYP)
      NSS      = NSITS(ITYP)   ! should be allways 3
      FNSPI    = 1.D0/DFLOAT(NSP)
C  Equilibrium distances
      ROH1     = RB(I1)
      ROH2     = RB(I2)
      RH1H2    = RB(I3)
C   Convert potential parameters to internal units
      FACTOR   = 1.D3/AVSNO/UNITE
      AAA      = AAP*FACTOR
      BBB      = BBP*FACTOR
      DDD      = DDP
      IPT      = IPOT(ITYP)
      IF(IPT.EQ.2) THEN
        AAA      = AAA/DDD**2
        BBB      = BBB/DDD**2
      END IF
      CCC      = CCP*FACTOR
      EEE      = EEP*FACTOR
      FFF      = FFP*FACTOR
      GGG      = GGP*FACTOR
      BBM	= BBMA
      BBS	= BBSA
      AADD     = AAA*DDD
      BBDD     = BBB*DDD
C   These are atom numbers for the first water molecule
      ISHF     = ISADDR(ITYP)
      IBB      = ISHF+1
      JBB      = ISHF+2
      KBB      = ISHF+3
*
*    7.3  Case of anharmonic Morse potential
*    ---------------------------------------
      IF(IPT.EQ.2) THEN
*  7.3.1  Cycle over molecules
C   (each node takes care of its own set of molecules)
      DO N     = TASKID+1,NSP,NUMTASK
C   Atom numbers
        I        = (N-1)*NSS+IBB
        J        = (N-1)*NSS+JBB
        K        = (N-1)*NSS+KBB
*   7.3.2 Bond lengths calculation
        AX       = SX(J)-SX(I)
        AY       = SY(J)-SY(I)
        AZ       = SZ(J)-SZ(I)
        call PBC(AX,AY,AZ)
C
        BX       = SX(K)-SX(I)
        BY       = SY(K)-SY(I)
        BZ       = SZ(K)-SZ(I)
        call PBC(BX,BY,BZ)
C
        CX       = SX(K)-SX(J)
        CY       = SY(K)-SY(J)
        CZ       = SZ(K)-SZ(J)
        call PBC(CX,CY,CZ)
C
        RRA      = AX*AX+AY*AY+AZ*AZ
        RRB      = BX*BX+BY*BY+BZ*BZ
        RRC      = CX*CX+CY*CY+CZ*CZ
        ASS      = DSQRT(RRA)
        BSS      = DSQRT(RRB)
        CSS      = DSQRT(RRC)
        if(ASS.gt.BBS)then
          write(*,'(a,f14.6,2i7,i10,i5)')
     +    ' !!! OH-bond length ',ASS,I,J,MSTEP,TASKID
          write(*,*)SX(I),SY(I),SZ(I)
          write(*,*)SX(J),SY(J),SZ(J)
CD         call XMOLDUMP
CD         stop
        end if  
        if(BSS.gt.BBS)then
          write(*,'(a,f14.6,2i7,i10,i5)')
     +    ' !!! OH-bond length ',BSS,I,K,MSTEP,TASKID
          write(*,*)SX(I),SY(I),SZ(I)
          write(*,*)SX(K),SY(K),SZ(K)
CD         call XMOLDUMP
CD         stop
        end if
*  7.3.3 Energy calculations
        SA       = ASS-ROH1
        SB       = BSS-ROH2
        SH       = CSS-RH1H2
C
        EXPA     = DEXP(-DDD*SA)
        EXPAM    = 1.D0-EXPA
        EXPB     = DEXP(-DDD*SB)
        EXPBM    = 1.D0-EXPB
C
        WB1       = AAA*EXPAM**2+EEE*SA*SH
        WB2       = BBB*EXPBM**2+FFF*SB*SH
        WB3       = CCC*SH*SH   +GGG*SA*SB
        WB        = WB1+WB2+WB3
*   7.3.4 Fill arrays for average energy and bond length
        BR(I1)    = BR(I1)+ASS*FNSPI
        BR(I2)    = BR(I2)+BSS*FNSPI
        BR(I3)    = BR(I3)+CSS*FNSPI
C
        BE(I1)    = BE(I1)+WB1*FNSPI
        BE(I2)    = BE(I2)+WB2*FNSPI
        BE(I3)    = BE(I3)+WB3*FNSPI
C   total intramolec. energy
        PINT	= PINT+WB1+WB2+WB3
*   7.3.5 Force calculations
        if(LMOVE(ITYP))then
        if(ASS.gt.BBM)then
C   We do not have intention to dissociate water...
          FB3A	= -2.*SA*AADD*DDD
        else
          FB3A	= -2.*AADD*EXPA*EXPAM
        end if
        FB3A     = FB3A   - EEE*SH-GGG*SB
        if(BSS.gt.BBM)then
          FB3B	= -2.*SB*BBDD*DDD
        else
          FB3B	= -2.*BBDD*EXPB*EXPBM
        end if
        FB3B     = FB3B   - FFF*SH-GGG*SA
        FB3C     = -2.0*CCC*SH            - EEE*SA-FFF*SB
*
        FB3A     = FB3A/ASS
        FB3B     = FB3B/BSS
        FB3C     = FB3C/CSS
*
        AXX       = FB3A*AX
        AYY       = FB3A*AY
        AZZ       = FB3A*AZ
*
        BXX       = FB3B*BX
        BYY       = FB3B*BY
        BZZ       = FB3B*BZ
*
        CXX       = FB3C*CX
        CYY       = FB3C*CY
        CZZ       = FB3C*CZ
*
        HX(I)    = HX(I)-AXX-BXX
        HY(I)    = HY(I)-AYY-BYY
        HZ(I)    = HZ(I)-AZZ-BZZ
*
        HX(J)    = HX(J)+AXX   -CXX
        HY(J)    = HY(J)+AYY   -CYY
        HZ(J)    = HZ(J)+AZZ   -CZZ
*
        HX(K)    = HX(K)   +BXX+CXX
        HY(K)    = HY(K)   +BYY+CYY
        HZ(K)    = HZ(K)   +BZZ+CZZ
*  7.3.6  Virial calculations
	VIX	= AX*AXX + BX*BXX + CX*CXX
	VIY	= AY*AYY + BY*BYY + CY*CYY
	VIZ	= AZ*AZZ + BZ*BZZ + CZ*CZZ
	VIRB	= VIRB + VIX+VIY+VIZ
	VIRFX	= VIRFX+VIX
	VIRFY	= VIRFY+VIY
	VIRFZ	= VIRFZ+VIZ
        end if
      END DO! OF N
*
      ELSE
*
*  7.4   Harmonic flexible SPC model
*  ---------------------------------
*  7.4.1  Cycle over water molecules
      DO N     = TASKID+1,NSP,NUMTASK
        I        = (N-1)*NSS+IBB
        J        = (N-1)*NSS+JBB
        K        = (N-1)*NSS+KBB
*  7.4.2  Bond length calculations
        AX       = SX(J)-SX(I)
        AY       = SY(J)-SY(I)
        AZ       = SZ(J)-SZ(I)
        BX       = SX(K)-SX(I)
        BY       = SY(K)-SY(I)
        BZ       = SZ(K)-SZ(I)
        CX       = SX(K)-SX(J)
        CY       = SY(K)-SY(J)
        CZ       = SZ(K)-SZ(J)
        RRA      = AX*AX+AY*AY+AZ*AZ
        RRB      = BX*BX+BY*BY+BZ*BZ
        RRC      = CX*CX+CY*CY+CZ*CZ
        ASS      = DSQRT(RRA)
        BSS      = DSQRT(RRB)
        CSS      = DSQRT(RRC)
*  7.4.3  Deviations
        SA       = ASS-ROH1
        SB       = BSS-ROH2
        SH       = CSS-RH1H2
*  7.4.4 Energy calculations
        WB1      = AAA*SA*SA+EEE*SA*SH
        WB2      = BBB*SB*SB+FFF*SB*SH 
        WB3      = CCC*SH*SH+GGG*SA*SB
        WB       = WB1+WB2+WB3
        BR(I1)    = BR(I1)+ASS*FNSPI
        BR(I2)    = BR(I2)+BSS*FNSPI
        BR(I3)    = BR(I3)+CSS*FNSPI
        BE(I1)    = BE(I1)+WB1*FNSPI
        BE(I2)    = BE(I2)+WB2*FNSPI
        BE(I3)    = BE(I3)+WB3*FNSPI
        PINT	= PINT+WB1+WB2+WB3
*  7.4.5  Force calculations
        if(LMOVE(ITYP))then
        FB3A     = -2.0*AAA*SA - EEE*SH-GGG*SB
        FB3B     = -2.0*BBB*SB - FFF*SH-GGG*SA
        FB3C     = -2.0*CCC*SH - EEE*SA-FFF*SB
        FB3A     = FB3A/ASS
        FB3B     = FB3B/BSS
        FB3C     = FB3C/CSS
        AXX       = FB3A*AX
        AYY       = FB3A*AY
        AZZ       = FB3A*AZ
        BXX       = FB3B*BX
        BYY       = FB3B*BY
        BZZ       = FB3B*BZ
        CXX       = FB3C*CX
        CYY       = FB3C*CY
        CZZ       = FB3C*CZ
        HX(I)    = HX(I)-AXX-BXX                           
        HY(I)    = HY(I)-AYY-BYY
        HZ(I)    = HZ(I)-AZZ-BZZ
        HX(J)    = HX(J)+AXX   -CXX
        HY(J)    = HY(J)+AYY   -CYY
        HZ(J)    = HZ(J)+AZZ   -CZZ
        HX(K)    = HX(K)   +BXX+CXX
        HY(K)    = HY(K)   +BYY+CYY
        HZ(K)    = HZ(K)   +BZZ+CZZ
*  7.4.6 Virial calculations                      
	VIX	= AX*AXX + BX*BXX + CX*CXX
	VIY	= AY*AYY + BY*BYY + CY*CYY
	VIZ	= AZ*AZZ + BZ*BZZ + CZ*CZZ
	VIRB	= VIRB + VIX+VIY+VIZ
	VIRFX	= VIRFX+VIX
	VIRFY	= VIRFY+VIY
	VIRFZ	= VIRFZ+VIZ
        end if   !  if(LMOVE
      END DO! OF N
      END IF
      RETURN
      END
*
*
*    8. Calculation of LJ and electrostatic forces - small distances 
*    ----------------------------------------------------------------
C    This procedure calculates forces due to Lennard-Jones and real-space
C    electrostatics interactions for atoms on small distances:
C    closer than SHORT. These forces are calculated each short time step.
C    Interactions between atoms with distance more than SHORT are 
C    calculated each short time step (see subr. FORCES, p.8) 
C
C======================= LOCALF ==========================================
C
      SUBROUTINE LOCALF
      include "prcm.h"  
      times0	= cputime(0.d0)
      IF(MBSH.le.0) RETURN
*
*   8.1  Organize cycle over atom pairs 
*   -----------------------------------
      DO INP      = 1,MBSH
	ISP0=NBS1(INP) 
        ISP=iabs(ISP0)           !   atom
	JSP=NBS2(INP)
	ITYP=ITYPE(ISP)          !   type
	JTYP=ITYPE(JSP)
	ISB =NSITE(ISP)          !   site
  	JSB =NSITE(JSP)
	IMOL=NNUM(ISP)           !   molecule
	JMOL=NNUM(JSP)
	MTR=MDX(ITYP,JTYP)       !  type-pair index
*
*   8.2  define interaction parameters for the given pair
*   -----------------------------------------------------
	if(ISP0.le.0)then  
* signals that this pair is 1-4 neigbours
	  if(ITYP.ne.JTYP)write(*,*)
     +'!!! different mol. types for 1-4 interaction:',
     +ITYP,JTYP,ISP,JSP,' FF'
	  FCTLJ=C14LJ(ITYP)
*   8.2.1 special 1-4 interactions
c   Look through the list of special 1-4 interactions
C   and check whether this pair belong to it
C   This procedure may seem cumbersobe, but I didnot find better way
C   without using (NSITES,NSITES) array. Really this is not very
C   time-consuming; saving memory is preferable
	  do I=1,NNAD 
	    ILA=ILJ(I)
	    if(ILA.eq.ISB)then  
	      do J=1,NNAD 
	        JLA=ILJ(J)
	        if(JLA.eq.JSB)then
                  EPSI       =  EPAD(I)*EPAD(J)
                  if(LGEOM)then
                    SIGM       =  SIGAD(I)*SIGAD(J)
                  else
                    SIGM       =  SIGAD(I)+SIGAD(J)
                  end if
                  go to 195
	        end if
	      end do
                EPSI       =  EPAD(I)*EPSIL(JSB)
                if(LGEOM)then
                  SIGM       =  SIGAD(I)*SIGMA(JSB)
                else
                  SIGM       =  SIGAD(I)+SIGMA(JSB)
                end if
              go to 195
	    end if
	  end do
	  do J=1,NNAD
	    JLA=ILJ(J)
	    if(JLA.eq.JSB)then
              EPSI       =  EPSIL(ISB)*EPAD(J)
              if(LGEOM)then
                SIGM       =  SIGMA(ISB)*SIGAD(J)
              else
                SIGM       =  SIGMA(ISB)+SIGAD(J)
              end if
              go to 195
	    end if
	  end do
*   8.2.2  Normal 1-4 connected pairs
C
C   Note: EPSIL here, above and below is really sqrt(epsil) in LJ potential
C  (see also setup.f )
C
          EPSI       =  EPSIL(ISB)*EPSIL(JSB)
          if(LGEOM)then
            SIGM       =  SIGMA(ISB)*SIGMA(JSB)
          else
            SIGM       =  SIGMA(ISB)+SIGMA(JSB)
          end if
 195      A6         = -EPSI*SIGM**6*FCTLJ    
          B12       =  EPSI*SIGM**12*FCTLJ    
  	  QIJ=CHARGE(ISB)*CHARGE(JSB)*COULF*C14EL(ITYP)
*   8.2.3 Normal nonbonded pairs 
	else  ! if(ISP0   
*
*  Note:  SIGMA here is SIGMA/2 or (if LGEOM) sqrt(sigma) 
*     and EPSIL here is 2*sqrt(EPSIL)
*
*   Geometrical rule for sigma
          if(LGEOM)then
            EPSI       =  EPSIL(ISB)*EPSIL(JSB)
            SIGM       =  SIGMA(ISB)*SIGMA(JSB)
            A6         = -EPSI*SIGM**6
            B12       =  EPSI*SIGM**12
          else if(LKONG)then
*   Kong rules
            if(EPSIL(ISB).gt.0.d0)then
              A6=-EPSIL(ISB)*EPSIL(JSB)*(SIGMA(ISB)*SIGMA(JSB)*4.d0)**3
              TTWLI = (EPSIL(ISB)*SIGMA(ISB)**6)**2 
              TTWLJ = (EPSIL(JSB)*SIGMA(JSB)**6)**2
              B12 = 0.5*(1.d0+(TTWLJ/TTWLI)**(1./13.))**13*TTWLI
            else
              B12=0.
              A6=0.
            end if
*         use Lorentz-Berthelot rules
          else
            EPSI       =  EPSIL(ISB)*EPSIL(JSB)
            SIGM       =  SIGMA(ISB)+SIGMA(JSB)
            A6         = -EPSI*SIGM**6
            B12       =  EPSI*SIGM**12
          end if
  	  QIJ   = CHARGE(ISB)*CHARGE(JSB)*COULF
	end if
*
*   8.3  Calculate distances between the atoms
*   ------------------------------------------
*   8.3.1  Vector between the atom pair
        DX           = SX(ISP)-SX(JSP)
        DY           = SY(ISP)-SY(JSP)
        DZ           = SZ(ISP)-SZ(JSP)
*   8.3.2  Apply periodic boundary conditions
*                  call PBC(DX,DY,DZ)
        if(DX.gt. HBOXL)DX=DX-BOXL
        if(DX.lt.-HBOXL)DX=DX+BOXL
        if(DZ.gt. HBOZL)DZ=DZ-BOZL
        if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	if(LHEX)then
C  hexagonal periodic cell
	  XY=BOXYC*DX
          if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
            DY=DY-BOYL
            DX=DX-HBOXL
	    XY=BOXYC*DX
          end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
            DY=DY+BOYL
            DX=DX+HBOXL
	    XY=BOXYC*DX
          end if
          if(DY.gt.BOXY3+XY)then
            DY=DY-BOYL
            DX=DX+HBOXL
	    XY=BOXYC*DX
          end if
          if(DY.lt.-BOXY3+XY)then
            DY=DY+BOYL
            DX=DX-HBOXL
          end if
	else                   
C  rectangular cell
          if(DY.gt. HBOYL)DY=DY-BOYL
          if(DY.lt.-HBOYL)DY=DY+BOYL
	end if
	if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	  CORSS = HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	  DX=DX-sign(CORSS,DX)
	  DY=DY-sign(CORSS,DY)
	  DZ=DZ-sign(CORSS,DZ)
	end if
*   8.3.3  Calculate the distance and its exponents
	RR           = DX**2+DY**2+DZ**2
        R1           = DSQRT(RR)
        R2I          = 1.0D0/RR
	R1I	   = R1*R2I
*
*  8.4  Electrostatic interaction
*  ------------------------------
*  8.4.1   1-4 connected pairs - direct Coulomb
	if(ISP0.le.0)then 
	  EE1 = QIJ*R1I
	  if(LMOVE(ITYP))then
	    FORCE = EE1*R2I
	  else
	    FORCE=0.
	  end if
	  EES = EE1
	else
*  8.4.2 Reaction field method
	  if(LMNBL)then
	    EE1 = QIJ*(R1I-0.5*RFF2*RR+RFF)
	    EES = EE1
	    FORCE = QIJ*(R1I*R2I+RFF2)
C  correction to virial from the reaction field is added to LJ virial 	      
	    VIR2	= VIR2 + QIJ*(1.5*RFF2*RR-RFF)
*  8.4.3  Ewald method   
	  else
            ALPHAR       = ALPHAD*R1
	    TTT          = 1.D0/(1.D0+B1*ALPHAR)
            EXP2A        = DEXP(-ALPHAR**2)
            ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A 
	    EE1	= QIJ*R1I
            EES          = EE1*ERFC
            FORCE        = EE1*(TOTPI*ALPHAR*EXP2A+ERFC)*R2I
	  end if
	end if
*   8.4.4  Contributions to the electrostatic energy
	if(IMOL.eq.JMOL)then
	  PES141(ITYP)=PES141(ITYP)+EE1
	else
	  POTE1(MTR)    = POTE1(MTR)+EE1
	end if 
	PELS2=PELS2+EES
*
*   8.5  LJ interaction
*   -------------------
	if(B12.ne.0.d0)then
	  R6I	  = R2I*R2I*R2I
          ESR   = (     A6+      B12*R6I)*R6I
          FSR   = (6.D0*A6+12.D0*B12*R6I)*R6I*R2I
	  FORCE = FORCE+FSR
 	  if(IMOL.eq.JMOL)then
	    if(LMOVE(ITYP))then
              PSR141(ITYP) = PSR141(ITYP)+ESR
              VIR2        = VIR2+FSR*RR
	    else
	      ESR=0.
	      FORCE=0.
 	    end if
	  else
	    POTL1(MTR)    = POTL1(MTR)+ESR
            VIR2		= VIR2+FSR*RR
	  end if
	  PE2		= PE2+ESR
	end if
*
*  8.6  Calculate forces and virials
*  --------------------------------
        HX(ISP)      = HX(ISP)+FORCE*DX
        HY(ISP)      = HY(ISP)+FORCE*DY
        HZ(ISP)      = HZ(ISP)+FORCE*DZ
        HX(JSP)      = HX(JSP)-FORCE*DX
        HY(JSP)      = HY(JSP)-FORCE*DY
        HZ(JSP)      = HZ(JSP)-FORCE*DZ 
	VIRFX	= VIRFX+FORCE*DX**2  
	VIRFY	= VIRFY+FORCE*DY**2  
	VIRFZ	= VIRFZ+FORCE*DZ**2  
      END DO! OF INP					
*
*   8.7  Intramolecular correction for "molecular" virial
*   -----------------------------------------------------
      DO N         = 1,NSTOT
	M	= NNUM(N)
        WIRSS	= WIRSS-HX(N)*(SX(N)-X(M))-HY(N)*(SY(N)-Y(M))-
     +               HZ(N)*(SZ(N)-Z(M))
      END DO! OF N
      times	= times+cputime(times0)
      RETURN
      END
*
*    9. Calculation of LJ and electrostatic forces - medium distances 
*    ----------------------------------------------------------------
C    This procedure calculates forces due to Lennard-Jones and real-space
C    electrostatics interactions for atoms on "medium" distances:
C    between SHORT and RCUT. These forces are calculated each long time step.
C    Interactions between atoms with distance less than SHORT are 
C    calculated each short time step (see subr. LOCALF, p.8) 
C    Interactions out RCUT are assumed insignificant and do not calculated
C    at all (long-range electrostatic forces are calculated by Ewald
C    method, see subr. FURIER)
C    Program flow is mostly identical to the previous section
C
C=============== SLOW FORCES =============================================
C
      SUBROUTINE FORCES
      include "prcm.h"
      timel0	= cputime(0.d0)
C  if nothing to do ....
      IF(MBLN.le.0) go to 100
*
*   9.1 Organize cycle over atom pairs
*   ----------------------------------
C   List of atom pairs is determined in subroutine CHCNB 
c   This list is local for each node
      DO INP      = 1,MBLN           !  pair number
 	  ISP0=NBL1(INP)               !  atom number
	  JSP=NBL2(INP)
	  ISP=iabs(ISP0)
	  ITYP=ITYPE(ISP)              !  type number
	  JTYP=ITYPE(JSP)
	  ISB =NSITE(ISP)              !  site number
  	  JSB =NSITE(JSP)   
	  IMOL=NNUM(ISP)               !  molecule number
	  JMOL=NNUM(JSP)
	  MTR=MDX(ITYP,JTYP)
*
*   9.2 Calculate interaction parameters for the given atom pair 
*   ------------------------------------------------------------
	  if(ISP0.le.0)then 
C   This can not happen for different molecules... 
	    if(ITYP.ne.JTYP)write(*,*)
     +'!!! different mol. types for 1-4 interaction:',
     +ITYP,JTYP,ISP,JSP,' SF'
C   Scaling constant from the input 
	    FCTLJ=C14LJ(ITYP)
*   9.2.1 Case of special 1-4 interactions
C   Special LJ parameters for 1-4 interactions may be defined in the .mol file
	  do I=1,NNAD 
	    ILA=ILJ(I)
	    if(ILA.eq.ISB)then  
	      do J=1,NNAD 
	        JLA=ILJ(J)
	        if(JLA.eq.JSB)then
                  EPSI       =  EPAD(I)*EPAD(J)
                  if(LGEOM)then
                    SIGM       =  SIGAD(I)*SIGAD(J)
                  else
                    SIGM       =  SIGAD(I)+SIGAD(J)
                  end if
                  go to 95
	        end if
	      end do
                EPSI       =  EPAD(I)*EPSIL(JSB)
                if(LGEOM)then
                  SIGM       =  SIGAD(I)*SIGMA(JSB)
                else
                  SIGM       =  SIGAD(I)+SIGMA(JSB)
                end if
              go to 95
	    end if
	  end do
	  do J=1,NNAD
	    JLA=ILJ(J)
	    if(JLA.eq.JSB)then
              EPSI       =  EPSIL(ISB)*EPAD(J)
              if(LGEOM)then
                SIGM       =  SIGMA(ISB)*SIGAD(J)
              else
                SIGM       =  SIGMA(ISB)+SIGAD(J)
              end if
              go to 95
	    end if
	  end do
*   9.2.2  Normal 1-4 interactions
C       (defined by normal Lorentz-Berthelot combination rules,
C       with scaling factor) 
            EPSI       =  EPSIL(ISB)*EPSIL(JSB)
            if(LGEOM)then
              SIGM       =  SIGMA(ISB)*SIGMA(JSB)
            else
              SIGM       =  SIGMA(ISB)+SIGMA(JSB)
            end if
  95        A6         = -EPSI*SIGM**6*FCTLJ    
            B12        =  EPSI*SIGM**12*FCTLJ    
  	    QIJ=CHARGE(ISB)*CHARGE(JSB)*COULF*C14EL(ITYP)
*
	  else    ! normal non-bonded pairs
*
*   9.2.3 Normal LJ parameters 
*
*  Note:  SIGMA here is SIGMA/2 or (if LGEOM) sqrt(sigma) 
*     and EPSIL here is 2*sqrt(EPSIL)
*
*   Geometrical rule for sigma
            if(LGEOM)then
              EPSI       =  EPSIL(ISB)*EPSIL(JSB)
              SIGM       =  SIGMA(ISB)*SIGMA(JSB)
              A6         = -EPSI*SIGM**6
              B12       =  EPSI*SIGM**12
            else if(LKONG)then
*   Kong rules
              if(EPSIL(ISB).gt.0.d0)then
            A6=-EPSIL(ISB)*EPSIL(JSB)*(SIGMA(ISB)*SIGMA(JSB)*4.d0)**3
                TTWLI = (EPSIL(ISB)*SIGMA(ISB)**6)**2 
                TTWLJ = (EPSIL(JSB)*SIGMA(JSB)**6)**2
                B12 = 0.5*(1.d0+(TTWLJ/TTWLI)**(1./13.))**13*TTWLI
              else
                B12=0.
                A6=0.
              end if
*         use Lorentz-Berthelot rules
            else
              EPSI       =  EPSIL(ISB)*EPSIL(JSB)
              SIGM       =  SIGMA(ISB)+SIGMA(JSB)
              A6         = -EPSI*SIGM**6
              B12       =  EPSI*SIGM**12
            end if
  	    QIJ=CHARGE(ISB)*CHARGE(JSB)*COULF
	  end if
*
*   9.3 Calculate the distance between the atoms
*   --------------------------------------------
*   9.3.1 Define the vector 
          DX           = SX(ISP)-SX(JSP)
          DY           = SY(ISP)-SY(JSP)
          DZ           = SZ(ISP)-SZ(JSP)
*   9.3.2 Apply periodic boundary conditions
*                  call PBC(DX,DY,DZ)
C  Inline periodic boundary conditions
          if(DX.gt. HBOXL)DX=DX-BOXL
          if(DX.lt.-HBOXL)DX=DX+BOXL
          if(DZ.gt. HBOZL)DZ=DZ-BOZL
          if(DZ.lt.-HBOZL)DZ=DZ+BOZL   
	  if(LHEX)then
C  hexagonal periodic cell)
	    XY=BOXYC*DX
          if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
              DY=DY-BOYL
              DX=DX-HBOXL
	      XY=BOXYC*DX
            end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
              DY=DY+BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.gt.BOXY3+XY)then
              DY=DY-BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3+XY)then
              DY=DY+BOYL
              DX=DX-HBOXL
            end if
	  else    ! not.LHEX                   
C  rectangular cell
            if(DY.gt. HBOYL)DY=DY-BOYL
            if(DY.lt.-HBOYL)DY=DY+BOYL
	  end if
	  if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	    CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	    DX=DX-sign(CORSS,DX)
	    DY=DY-sign(CORSS,DY)
	    DZ=DZ-sign(CORSS,DZ)
	  end if
*  9.3.3 Calculate the distance
	  RR           = DX**2+DY**2+DZ**2
          R1           = DSQRT(RR)
          R2I          = 1.0D0/RR
	  R1I	   = R1*R2I
*
*  9.4 Electrostatic interaction
*  ---------------------------------------
C  EE1 is contribution to type-type energy (direct Coulomb term)
C  EES is exact contribution to the total el-st energy
C  FORCE is |force|/r
C
*  9.4.1  1-4 connected atoms
C  Direct Coulomb forces
	  if(ISP0.le.0)then                             
	    EE1 = QIJ*R1I
	    if(LMOVE(ITYP))then
	      FORCE = EE1*R2I
	    else
	      FORCE=0.
	    end if 
	    EES = EE1
	  else
	    if(LMNBL)then
*  9.4.2 Reaction field method
	      EE1 = QIJ*(R1I-0.5*RFF2*RR+RFF)
	      EES = EE1
	      FORCE = QIJ*(R1I*R2I+RFF2)
C  correction to virial from the reaction field is added to LJ virial 	      
	      VIR1	= VIR1 + QIJ*(1.5*RFF2*RR-RFF)
	    else
*  9.4.3  Ewald method (real space term)
C  Polinom expansion for erfc is used
              ALPHAR       = ALPHAD*R1
	      TTT          = 1.D0/(1.D0+B1*ALPHAR)
              EXP2A        = DEXP(-ALPHAR**2)
              ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A
	      EE1	= QIJ*R1I
              EES       = EE1*ERFC
              FORCE     = EE1*(TOTPI*ALPHAR*EXP2A+ERFC)*R2I 
	    end if
	  end if
	  if(IMOL.eq.JMOL)then
	    PES141(ITYP)=PES141(ITYP)+EE1
	  else
	    POTE1(MTR)    = POTE1(MTR)+EE1
	  end if 
	  PELS1=PELS1+EES
* 
*  9.5  LJ interactions
*  --------------------
	if(B12.ne.0.d0)then
	  R6I          = R2I*R2I*R2I
          ESR          = (     A6+      B12*R6I)*R6I
          FSR          = (6.D0*A6+12.D0*B12*R6I)*R6I*R2I
	  FORCE	= FORCE+FSR   
	  if(IMOL.eq.JMOL)then
	    if(LMOVE(ITYP))VIR1	= VIR1+FSR*RR
            PSR141(ITYP) = PSR141(ITYP)+ESR
	  else			
	    POTL1(MTR)   = POTL1(MTR)+ESR
	    VIR1	= VIR1+FSR*RR
          end if
          PE1 = PE1 + ESR
        end if
*  
*  9.6 Contributions to virial projections
*  -----------------------------------------
	if(IMOL.ne.JMOL.or.LMOVE(ITYP))then
	  VIRX	= VIRX+FORCE*DX**2  
	  VIRY	= VIRY+FORCE*DY**2  
	  VIRZ	= VIRZ+FORCE*DZ**2  
	end if
CD	   write(*,*)ISP,JSP,IMOL,JMOL,R1,ESR*PERMOL,EES*PERMOL
*
*  9.7  Calculate forces
*  ---------------------
        GX(ISP)      = GX(ISP)+FORCE*DX     
        GY(ISP)      = GY(ISP)+FORCE*DY
        GZ(ISP)      = GZ(ISP)+FORCE*DZ
        GX(JSP)      = GX(JSP)-FORCE*DX
        GY(JSP)      = GY(JSP)-FORCE*DY
        GZ(JSP)      = GZ(JSP)-FORCE*DZ
      END DO! OF INP					
 100  timel	= timel+cputime(timel0) 
      RETURN
      END
*
*=============================================================
*
*   10. Ewald method. Reciprocal space contribution
*   -----------------------------------------------
C   This procedure calculates contribution to the energy and forces 
C   from long-range part of electrostatic interactions, calculated 
C   in the reciprocal space by the Furie transform
C
C   Subroutines:
C   FURIR - the main procedure for Ewald sum calculation
C   INIFUR - initialisation of arrays for Ewald sum calculation
C=============== FURIR =======================================
C
      SUBROUTINE FURIR
*
*   10.1  Definitions     
*   -----------------
      include "prcm.h"
      parameter (NKVM=50000) ! see also 2 next subroutines
      parameter (FKTX=1.05,FKTS=0.95)
      parameter(NPUSH=10000)
      DIMENSION   TSIN(NTOT),TCOS(NTOT)
      dimension RX(NKVM),RY(NKVM),RZ(NKVM),AK(NKVM)
      data INIT/0/
*
*   10.2  Initialization
*   --------------------
      timef0	= cputime(0.d0)
      if(INIT.eq.0)then 
        INIT=1 
        NFPUSH=NPUSH/NSTOT+1
*  10.2.1  remember box size
	BOXLO=BOXL
	BOYLO=BOYL
	BOZLO=BOZL
*  10.2.2  Initialize arrays for the calculations
	if(ALPHAD.gt.0.d0)then
          call INIFUR(RX,RY,RZ,NKV,BOXL,BOYL,BOZL,FEXP,ALPHAD,ICELL)
	  ALP2=ALPHAD**2
          if(TASKID.eq.MAST.and.IPRINT.ge.4)write(*,*)
     +'Ewald sum:  Num. of k-vectors ',NKV
*  10.2.3  Case when Ewald method is not used
        else
	  NKV=0
          if(TASKID.eq.MAST.and..not.LMNBL)write(*,*)
     +'No Ewald'
	  end if
	end if     
        IF(NKV.eq.0.or.NTYPES.EQ.1.AND.NSPEC(1).EQ.1) RETURN
*  10.2.4  If box size changed too strong, reinitialize
	DBX=BOXL/BOXLO
	DBY=BOYL/BOYLO
	DBZ=BOZL/BOZLO
	if(DBX.gt.FKTX.or.DBX.lt.FKTS.or.DBY.gt.FKTX.or.
     *   DBY.lt.FKTS.or.DBZ.gt.FKTX.or.DBZ.lt.FKTS)then
	  BOXLO=BOXL
	  BOYLO=BOYL
	  BOZLO=BOZL
          call INIFUR(RX,RY,RZ,NKV,BOXL,BOYL,BOZL,FEXP,ALPHAD,ICELL)
	  ALP2=ALPHAD**2
          if(TASKID.eq.MAST.and.IPRINT.ge.4)write(*,*)
     +'Box size changed. New num. of k-vectors ',NKV
      end if
*
*  10.3  Reciprocal space Ewald sum
*  --------------------------------
      QPE         = 0.D0
      CFUR=4.d0*PI*COULF/VOL     
*  10.3.1  Start cycle over vectors in the reciprocal space
      IBEG=NUMTASK-TASKID
        if(LVISUAL)then
          if(mod(IKV,NFPUSH).eq.0)write(*,'(a)')'@mm '
        end if
      do IKV=IBEG,NKV,NUMTASK
C   These are coordinates of reciprocal vectors
        RXV=RX(IKV)/BOXL
        RYV=RY(IKV)/BOYL
        RZV=RZ(IKV)/BOZL 
        RK2=RXV**2+RYV**2+RZV**2
        SCS=0.
        SSN=0.
*  10.3.2 Calculate Sin(rk) and Cos(rk) and their sums for all the atoms
        do I=1,NSTOT       
C  Scalar product r.k 
          SCP=RXV*SX(I)+RYV*SY(I)+RZV*SZ(I)
C  sin and cos are taken from look-up tables
          QQ=Q(I)
          SINSC=QQ*sin(SCP)
          COSSC=QQ*cos(SCP)
C   they were calculated at initialization time
C	  if(SCP.le.0.d0)then
C	    XT=-SCP*DTAB
C	    IID=1
C          else
C	    XT=SCP*DTAB
C	    IID=0
C	  end if
C	  IOT=idint(XT)
C	  IU=mod(IOT,LTAB)                             
C	  IU1=IU+1
C	  DD=XT-dfloat(IOT) 
C	  DD1=1.d0-DD
C          SINSC=Q(I)*(DD1*SINT(IU)+DD*SINT(IU1))
C          COSSC=Q(I)*(DD1*COST(IU)+DD*COST(IU1))
C	  if(IID.ne.0)SINSC=-SINSC
*  10.3.3 remember sin and cos in temporarly arrays (for forces)
          TSIN(I)=SINSC
          TCOS(I)=COSSC
C  These are sums
          SSN=SSN+SINSC
          SCS=SCS+COSSC
        end do
*  10.3.4 Calculate energy 
        AKC=CFUR*exp(-0.25d0*RK2/ALP2)/RK2
        QPE=QPE+AKC*(SSN**2+SCS**2)
*  10.3.5 Calculate contribution to forces
        AK2=2.d0*AKC
        do I=1,NSTOT
          QFOR=AK2*(SCS*TSIN(I)-SSN*TCOS(I))
          GX(I)=GX(I)+RXV*QFOR
          GY(I)=GY(I)+RYV*QFOR
          GZ(I)=GZ(I)+RZV*QFOR
        end do
*  10.3.6  Calculate contribution to virial projections
	AKV=0.5*AKC*(4.*ALP2+RK2)/(ALP2*RK2) 
	VIRX=VIRX-RXV**2*AKV*(SSN**2+SCS**2) 
	VIRY=VIRY-RYV**2*AKV*(SSN**2+SCS**2)
	VIRZ=VIRZ-RZV**2*AKV*(SSN**2+SCS**2)
      end do  
      VIRX=VIRX+QPE
      VIRY=VIRY+QPE
      VIRZ=VIRZ+QPE
      PELS1=PELS1+QPE
      timef	= timef+cputime(timef0)
      return
      end
C=================================================================
*
*   11. Initialisation of Ewald sum calculations
*   --------------------------------------------
*
C  This subroutine defines a set of vectors in the reciprocal space
C  which will be used for Ewald sum calculations in dependence of geometry
C  Basis vectors in the reciprocal space are determined from vector
C  products:  Kx = Ly x Lz and so on (Lxyz are basis vectros in real space) 
C  For hexagonal cell:  Kx=(1,1/sqrt(3),0), Ky=(0,2/sqrt(3),0), Kz=(0,0,1)
C  For octahedron cell: Kx=(1,0,-1), Ky=(0,1,-1), Kz=(0,0,2)  
C  (multiplied by 2*Pi/Lbox in each direction; note that for hexagonal cell
C  BOYL defined as BOYL=sqrt(3)*Ly/2)
*
*================= INIFUR =========================================
*
      subroutine INIFUR
     +(RX,RY,RZ,NKV,BOXL,BOYL,BOZL,FEXP,ALPHAD,ICELL) 
      implicit real*8 (A-H,O-Z)
      parameter (NKVM=50000,PI=3.1415926536)
      dimension RX(NKVM),RY(NKVM),RZ(NKVM)
*
*  11.1  Calculate look-up tables for sin and cos
*  ----------------------------------------------
C      call SINTAB
*
*  11.2  Calculate parameters which set up cut-off in the reciprocal space
*  -----------------------------------------------------------------------
      ALP2=ALPHAD**2
      CUTK2=4.d0*ALP2*FEXP
      CUTK=sqrt(CUTK2)
      TWOPI=2.d0*PI
      IKV=0
      KMAX=CUTK*BOXL/TWOPI+1
      KMAY=CUTK*BOYL/TWOPI+1
      KMAZ=CUTK*BOZL/TWOPI+1
*
*   11.3  Calculate set of reciprocal vectors
*   -----------------------------------------
*   11.3.1 Cycle over reciprocal vectors in a cube
      do IKZ=-KMAZ,KMAZ
        do IKX=0,KMAX
          if(IKX.eq.0)then
            if(IKZ.gt.0)then
              KY1=0
            else
              KY1=1
            end if
            KY2=KMAY
          else
            KY1=-KMAY
            KY2=KMAY
          end if
          do IKY=KY1,KY2
*   11.3.2 Calculate cartesian coordinates of reciprocal vectors
            if(ICELL.eq.1)then
C  Truncated octahedron
C  These are cartesian coordinates of reciprocal vectors multiplied by 
C  the box size (which may change in NPT-dynamics)
              RXV = TWOPI*IKX
              RYV = TWOPI*IKY
              RZV = TWOPI*(-IKX+IKY+2*IKZ)
            else if(ICELL.eq.2)then
C  Hexagonal
              RXV = TWOPI*IKX
              RYV = TWOPI*IKY+PI*IKX
              RZV = TWOPI*IKZ
            else  ! ICELL=0 
C  Rectangular cell
              RXV = TWOPI*IKX
              RYV = TWOPI*IKY
              RZV = TWOPI*IKZ
            end if
*   11.3.3 Remember reciprocal vectors which are within cutoff
            RK2 = (RXV/BOXL)**2+(RYV/BOYL)**2+(RZV/BOZL)**2
            if(RK2.le.CUTK2)then
      	      IKV=IKV+1
              if(IKV.gt.NKVM)then
                write(*,*)'max number of recirpocal vectors exceeded'
                write(*,*)'check Ewald summation parameters or '
                write(*,*)'increase NKVM in forces.f/FURIER and INIFUR'
                stop 'in INIFUR'
              end if
              RX(IKV)=RXV
              RY(IKV)=RYV
              RZ(IKV)=RZV
            end if
          end do
        end do
      end do      
      NKV=IKV
      return	
      end
*
*  11.4  Create look-up tables for sin and cos
*  -------------------------------------------
*=========================== SINTAB =================================
*                                                                    
C	subroutine SINTAB
C	implicit real*8 (A-H,O-Z)
C	parameter (PI= 3.1415926535897932d0, PI2=2.d0*PI, LTAB=50000)
C	common /TSNCS/ DTAB,SINT(0:LTAB),COST(0:LTAB)
C	do I=0,LTAB 
C	  X=dfloat(I)*PI2/LTAB         
C	  sint(I)=sin(X)
C	  cost(I)=cos(X)
C	end do 
C	DTAB=dfloat(LTAB)/PI2
C	return
C	end
C========================================================================
*
*   12.  Calculation of self-interaction term in the Ewald sum
*   ----------------------------------------------------------
C
C   Really this procedure calculates two terms:
C   1) Standard self-interaction Ewald term, giving only constant 
C   contribution to the energy
C   2) "Molecular correction" which is due to the fact, that bound
C   atoms in a molecule do not interact electrostatically, but 
C   reciprocal Ewald term includes these interactions also.
*
*============= ETERMS =================================================
*
      SUBROUTINE ETERMS
*  12.1 Front matter
      include"prcm.h"         
      real*8 SPEEI(NTPS)
      data INIT/0/  
      save SPEEI,SPEI,ABC,INIT
      timee0	= cputime(0.d0)
      if(LMNBL)return
*
*  12.2  Electrostatic self-interaction
*  ------------------------------------
C  This constant term may be calculated only once
      if(INIT.eq.0)then
	INIT=1
        ABC = ALPHAD/DSQRT(PI)
	SPEI=0.
        do ITYP=1,NTYPES 
	  SPET=0.
	  do I=1,NSITS(ITYP)
	    IS=ISADR(ITYP)+I
  	    QII=COULF*CHARGE(IS)**2
            SPET         = SPET-NSPEC(ITYP)*ABC*QII
	  end do         
	  SPET=SPET/NUMTASK
	  SPEEI(ITYP)=SPET
	  SPEI=SPEI+SPET
	end do
      end if       
      do I=1,NTYPES
	SPEE(I)=SPEEI(I)
      end do            
      SPE=SPEI
*
*  12.3  Calculate intramolecular correction
*  -----------------------------------------
*  12.3.1 Start cycle over bound atoms 
C  list of bound atoms is prepeared in BNDLST (setup.f)   
      IBEG=NUMTASK-TASKID 
      do ISP=IBEG,NSTOT,NUMTASK
	NSP = NNBB(ISP)
	do J=1,NSP
	  JSP=    INBB(ISP,J)
	  if(JSP.le.0)JSP=-JSP
	  if(ISP.gt.JSP)then  ! avoid double count!
	    IS	= NSITE(ISP)
	    JS	= NSITE(JSP)
	    ITYP  = ITYPE(ISP)
  	    QIJ=COULF*CHARGE(IS)*CHARGE(JS)
            DX          = SX(ISP)-SX(JSP)
            DY          = SY(ISP)-SY(JSP)
            DZ          = SZ(ISP)-SZ(JSP)
*  12.3.2  Apply periodic boundary conditions
            if(DX.gt. HBOXL)DX=DX-BOXL
            if(DX.lt.-HBOXL)DX=DX+BOXL
	    if(LHEX)then
C  hexagonal periodic cell)
	      XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
              DY=DY-BOYL
              DX=DX-HBOXL
	      XY=BOXYC*DX
            end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
              DY=DY+BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.gt.BOXY3+XY)then
              DY=DY-BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3+XY)then
              DY=DY+BOYL
              DX=DX-HBOXL
            end if
	    else                   
C  rectangular cell
              if(DY.gt. HBOYL)DY=DY-BOYL
              if(DY.lt.-HBOYL)DY=DY+BOYL
	    end if
            if(DZ.gt. HBOZL)DZ=DZ-BOZL
            if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	    if(LOCT)then
	      CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	      DX=DX-sign(CORSS,DX)
	      DY=DY-sign(CORSS,DY)
	      DZ=DZ-sign(CORSS,DZ)
	    end if
*	  call PBC(DX,DY,DZ)
*   12.3.3 Calculate contributions to energy and forces
            R2          = DX**2+DY**2+DZ**2
            R2I         = 1.D0/R2
            R1          = DSQRT(R2)
            R1I         = R1*R2I
            ALPHAR      = ALPHAD*R1
*  	          call RERFC(ALPHAR,ERFC,EXP2A)
            TTT         = 1.D0/(1.D0+B1*ALPHAR)
            EXP2A       = DEXP(-ALPHAR**2)
            ERFC = ((((A5*TTT+A4)*TTT+A3)*TTT+A2)*TTT+A1)*TTT*EXP2A
            EES  = QIJ*(ERFC-1.d0)*R1I
            SPE         = SPE+EES
	    SPEE(ITYP)  = SPEE(ITYP)+EES
            if(LMOVE(ITYP))then 
              FES  = QIJ*((TOTPI*ALPHAR*EXP2A+ERFC)-1.d0)*R1I*R2I
              GX(ISP)     = GX(ISP)+FES*DX
              GY(ISP)     = GY(ISP)+FES*DY
              GZ(ISP)     = GZ(ISP)+FES*DZ
              GX(JSP)     = GX(JSP)-FES*DX
              GY(JSP)     = GY(JSP)-FES*DY
              GZ(JSP)     = GZ(JSP)-FES*DZ
	      VIRX	= VIRX+FES*DX**2  
	      VIRY	= VIRY+FES*DY**2  
	      VIRZ	= VIRZ+FES*DZ**2
            end if
	  end if
	end do 
      END DO! OF I
 100  timee	= timee+cputime(timee0)
      RETURN
      END
*=====================================================================
*                 
*   13. Calculate out-cut-off corrections to energy and pressure from
*       LJ interactions
*
*============= ELRCLJ =================================================
*
      SUBROUTINE ELRCLJ
      include"prcm.h" 
      dimension PELRC0(NTPP),VRLRC0(NTPP)
      data INIT/0/
      if(INIT.eq.0)then                                
        DO I        = 1,MOLINT
          PELRC0(I)    = 0.D0
          VRLRC0(I)    = 0.D0
        END DO! OF I
        RC          = DSQRT(RSSQ)
C  Cycle over all site pairs
        DO 400 ITYP = 1    ,NTYPES
          ISBEG       = ISADR(ITYP)+1
          ISEND       = ISADR(ITYP +1)
          NSPI        = NSPEC(ITYP)
          DO 400   IS = ISBEG,ISEND
            DO 300 JTYP = ITYP ,NTYPES
              MT          = MDX  (ITYP,JTYP)
              JSBEG       = ISADR(JTYP)+1
              JSEND       = ISADR(JTYP +1)
              NSPJ        = NSPEC(JTYP)
              FNOPIJ      = DFLOAT(NSPI*NSPJ)
              if(ITYP.eq.JTYP)FNOPIJ=0.5*FNOPIJ
              DO 300   JS = JSBEG,JSEND  
                EP       =   EPSIL(IS)*EPSIL(JS)
                if(EP.gt.1.d-20)then
	          SI3	= (SIGMA(IS)+SIGMA(JS))**3     ! SI**3      
      PELRC0(MT)   = PELRC0(MT)+EP*SI3*(( 1.D0/9.D0)*(SI3/RC**3)**3
     X                                - ( 1.D0/3.D0)*(SI3/RC**3))*FNOPIJ
      VRLRC0(MT)   = VRLRC0(MT)+EP*SI3*((-4.D0/3.D0)*(SI3/RC**3)**3
     X                                + (      2.D0)*(SI3/RC**3))*FNOPIJ
	        end if
  300       CONTINUE
  400   CONTINUE         
      end if
      FCV       = UNITP/(3.*VOL)
      VIRD	= 0.
      DO I      = 1,MOLINT
        PELRC(I)  = 4.*PI*PELRC0(I)/VOL
        VRLRC(I)  = -4.*PI*VRLRC0(I)/VOL
        VIRD=VIRD+VRLRC(I)
        if(INIT.eq.0.and.IPRINT.ge.5)then
          PRINT "('*** PE  LRC ',I2,': ',F12.8,' kJ/mol')",
     +    I,PELRC(I)*PERMOL
        PRINT "('*** PR  LRC ',I2,': ',F12.4,' atm   ')",I,VRLRC(I)*FCV
          INIT=1
        end if
      END DO
*
      RETURN
      END
C===================================================================
*
*   14. Harmonic force aroun fixed position
*   ---------------------------------------
C   
C   Each atom  may be put in a harmonic well around 
C   some position defined in "filref" file (if LRR=.true)
*====================== PULLBACK ====================================
* 
	subroutine PULLBACK
	include "prcm.h"     
        if((.not.LRR.or.NPULL.eq.0).and.IEXT.ne.1)return
        if(LRR)then
	  ELINK=0.          
          do IL=1,NPULL             ! num of linkage point
            I = INPL(IL)            ! num of atom
            if(I.le.0)write(*,*)' wrong index in PULLBACK!!!'
	    DX=SX(I)-XPL(IL)
	    DY=SY(I)-YPL(IL)
	    DZ=SZ(I)-ZPL(IL)
C   Normally, PBC are not used here
C	  call PBC(DX,DY,DZ)
	    RR2=DX**2+DY**2+DZ**2
	    ELINK=ELINK+RR2
	    GX(I)=GX(I)-RLR*DX
	    GY(I)=GY(I)-RLR*DY
	    GZ(I)=GZ(I)-RLR*DZ
	  end do
	  ELNK=ELINK*PERMOL*RLR 
	  if(IPRINT.ge.6)write(*,*)
     +'Deviation from correct structure ',sqrt(ELINK/NPULL),ELNK
        end if
*   External electrostatic field
        if(IEXT.eq.1)then
          FULTIM=TIM+NSTEP*DT
          EFIELD=EXAMPL*cos(EXFREQ*FULTIM*0.5/PI)
          do I=NAB(TASKID),NAE(TASKID)
             IS=NSITE(I)
             GZ(I)=GZ(I)+CHARGE(IS)*EFIELD
          end do
        end if
	return
	end      
C
C===========================================================================
*
*  15. Recalculation of list of neighbours
*  ---------------------------------------
C  This procedure calculate lists of closest and far neighbours
C  for subroutines LOCALF and FORCES
C  These lists must be recalculated after some time
C  Each node calculate goes through its own list of atom pairs 
C  and creates own lists, which are used in subroutines for force calculation
C
C================ CHCNB =======================================================
C     
      subroutine CHCNB(LRDF)
*  15.1 Front matter     
      include "prcm.h"
      logical LGRC,LRDF 
      parameter (RMX2=5.**2,MMOL=5)
      data INIT/0/ 
      timev0	= cputime(0.d0)
      MBLN=0
      MBSH=0
      IOVER=0
      LGRC=LGR.and.LRDF.and.(NHIST+1.ge.IHIST)
*
*  15.2  Organize cycle over ALL atom pairs
*  ----------------------------------------
      IBEG=NUMTASK-TASKID 
C  cycle over I,J is organized is follows:
C  it takes I>J if I and J are of the same parity and I<J if I and J have
C  different parity. This allows nearly uniformly distribute atom pairs 
C  between the nodes
      do I=IBEG,NSTOT,NUMTASK
	IMOL=NNUM(I) 
	ITYP=ITYPE(I)
*   15.2.1 Even strips
        do J=I-2,1,-2
	  JMOL=NNUM(J)
	  DX=SX(I)-SX(J)
	  DY=SY(I)-SY(J)
	  DZ=SZ(I)-SZ(J)
*   15.2.1.1 Periodic boundary conditions
C	    call PBC (DX,DY,DZ)
	  if(DX.gt. HBOXL)DX=DX-BOXL
          if(DX.lt.-HBOXL)DX=DX+BOXL
          if(DZ.gt. HBOZL)DZ=DZ-BOZL
          if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	  if(LHEX)then
C  hexagonal periodic cell)
	    XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
              DY=DY-BOYL
              DX=DX-HBOXL
	      XY=BOXYC*DX
            end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
              DY=DY+BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.gt.BOXY3+XY)then
              DY=DY-BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3+XY)then
              DY=DY+BOYL
              DX=DX-HBOXL
            end if
	  else                   
C  rectangular cell
            if(DY.gt. HBOYL)DY=DY-BOYL
            if(DY.lt.-HBOYL)DY=DY+BOYL
	  end if
	  if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	      CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	      DX=DX-sign(CORSS,DX)
	      DY=DY-sign(CORSS,DY)
	      DZ=DZ-sign(CORSS,DZ)
	    end if                  
	    RR2=DX**2+DY**2+DZ**2  
*  15.2.1.2 Case of "molecular" cutoff - calculate COM distance
C  If we do not use Ewald summation (LMNBL=.true.), 
C  it is important that at least
C  small molecules (less than MMOL atoms) where inside or outside
C  cut-off radius as a whole 
	    if(LMNBL)then
	      NMI=NSITS(ITYP) 
	      NMJ=NSITS(ITYPE(J))
	      if(NMI.gt.MMOL.and.NMJ.gt.MMOL)then
		RRM=RR2
		go to 110
	      else if(NMI.gt.MMOL.and.NMJ.le.MMOL)then
	        DX=SX(I)-X(JMOL)
	        DY=SY(I)-Y(JMOL)
	        DZ=SZ(I)-Z(JMOL)
	      else if(NMJ.gt.MMOL)then
	        DX=X(IMOL)-SX(J)
	        DY=Y(IMOL)-SY(J)
	        DZ=Z(IMOL)-SZ(J)
	      else 
	        DX=X(IMOL)-X(JMOL)
	        DY=Y(IMOL)-Y(JMOL)
	        DZ=Z(IMOL)-Z(JMOL)
	      end if
*	      call PBC (DX,DY,DZ)
	      if(DX.gt. HBOXL)DX=DX-BOXL
              if(DX.lt.-HBOXL)DX=DX+BOXL
              if(DZ.gt. HBOZL)DZ=DZ-BOZL
              if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	      if(LHEX)then
C  hexagonal periodic cell)
	        XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
                DY=DY-BOYL
                DX=DX-HBOXL
	        XY=BOXYC*DX
              end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
                DY=DY+BOYL
                DX=DX+HBOXL
	        XY=BOXYC*DX
              end if
              if(DY.gt.BOXY3+XY)then
                DY=DY-BOYL
                DX=DX+HBOXL
	        XY=BOXYC*DX
              end if
              if(DY.lt.-BOXY3+XY)then
                DY=DY+BOYL
                DX=DX-HBOXL
              end if
	      else                   
C  rectangular cell
                if(DY.gt. HBOYL)DY=DY-BOYL
                if(DY.lt.-HBOYL)DY=DY+BOYL
	      end if
	      if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	        CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	        DX=DX-sign(CORSS,DX)
	        DY=DY-sign(CORSS,DY)
	        DZ=DZ-sign(CORSS,DZ)
	      end if                  
C  square distance betwen COM of small molecules
	      RRM=DX**2+DY**2+DZ**2 
	    else   ! .not.LMNBL
	      RRM=RR2              
	    end if  ! if(LMNBL
*  15.2.1.3 Check if the atoms are bound
 110	    if(IMOL.eq.JMOL)then    ! bound atoms should not be in the list
              if(.not.L15NB(ITYP))go to 125
	      if(RR2.le.RMX2)then
C   Check that the atoms are non-bound. Check using binding list of the
C   atom which has less bound atoms
	        if(NNBB(I).le.NNBB(J))then
		  do K=1,NNBB(I) 
C   1-3 bound:  not in the list
		    if(INBB(I,K).eq.J)go to 125
                    if(INBB(I,K).eq.-J.and.L14NB(ITYP))then
C   1-4 bound:  put into list with negative first index
C   (these atom pairs may have a special way of force calculations)
  	              MBSH=MBSH+1
	              if(MBSH.gt.NBSMAX)write(*,*)
     +        '!!! list of closest neigbours overfilled. atoms ',I,J
	              NBS1(MBSH)=-I
		      NBS2(MBSH)=J
	              go to 125
		    end if  
	  	  end do
	        else
		  do K=1,NNBB(J) 
C   1-3 bound
		    if(INBB(J,K).eq.I)go to 125
                    if(INBB(J,K).eq.-I.and.L14NB(ITYP))then
C   1-4 bound                                     
	              MBSH=MBSH+1
	              if(MBSH.gt.NBSMAX)write(*,*)
     +        '!!! list of closest neigbours overfilled. atoms ',I,J
	              NBS1(MBSH)=-I
		      NBS2(MBSH)=J 
	              go to 125
		    end if  
	  	  end do
	        end if
	      end if ! (RR2.le.RMAX2
	    end if   ! (IMOL.eq.JMOL
*   15.2.1.4 Pick up non-bound atoms to the list of neigbours
	    if(RRM.le.SHORT2)then
	      MBSH=MBSH+1
	      if(MBSH.gt.NBSMAX)write(*,*)
     +      '!!! list of closest neigbours overfilled. atoms ',I,J
	      NBS1(MBSH)=I
	      NBS2(MBSH)=J
	    else if(RRM.le.RSSQ)then
	      MBLN=MBLN+1
	      if(MBLN.gt.NBLMAX)then
                if(IOVER.eq.0)write(*,*)
     +     '!! list of far neigbours overfilled'
                IOVER=1
              end if
              if(IOVER.eq.0)then
	        NBL1(MBLN)=I
	        NBL2(MBLN)=J
              end if
	    end if
*   15.2.1.5 RDF calculation
C   Only at this point we go through all the atom pairs
 125        IF(LGRC.and.RR2.le.RDFCUT2) THEN			!=====> LGR
	      ISB	= NSITE(J)
	      JSB	= NSITE(I) 
C  find RDF pair
	      if(NGRI(ISB).le.NGRI(JSB))then
		do JS=1,NGRI(ISB) 
		  if(IGRI(ISB,JS).eq.JSB)then 
		    IP = MGRI(ISB,JS)
	            go to 127  
		  end if
		end do
	      else
		do JS=1,NGRI(JSB)
		  if(IGRI(JSB,JS).eq.ISB)then
		    IP = MGRI(JSB,JS)
	            go to 127
		  end if  
		end do         
	      end if  
	      go to 128
 127	      R1	= sqrt(RR2)
              GG        = R1*DRDFI
	      MM        = IDINT(GG)+1
            if(MM.le.MAX)IRDF(MM,IP) = IRDF(MM,IP)+1
	    end if ! of LGRC
 128	    continue
	  end do
*  15.2.2   Odd strips
C  Everythig here is a complete analogue to p 15.2.1    
        do J=I+1,NSTOT,2
	  JMOL=NNUM(J)
	  DX=SX(I)-SX(J)
	  DY=SY(I)-SY(J)
          DZ=SZ(I)-SZ(J)
C	    call PBC (DX,DY,DZ)
	  if(DX.gt. HBOXL)DX=DX-BOXL
          if(DX.lt.-HBOXL)DX=DX+BOXL
          if(DZ.gt. HBOZL)DZ=DZ-BOZL
          if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	  if(LHEX)then
C  hexagonal periodic cell)
	    XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
              DY=DY-BOYL
              DX=DX-HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
              DY=DY+BOYL
              DX=DX+HBOXL
	      XY=BOXYC*DX
            end if
            if(DY.gt.BOXY3+XY)then
              DY=DY-BOYL
              DX=DX+HBOXL
              XY=BOXYC*DX
            end if
            if(DY.lt.-BOXY3+XY)then
              DY=DY+BOYL
              DX=DX-HBOXL
            end if
	  else                   
C  rectangular cell
            if(DY.gt. HBOYL)DY=DY-BOYL
            if(DY.lt.-HBOYL)DY=DY+BOYL
	  end if
	  if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	    CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	    DX=DX-sign(CORSS,DX)
	    DY=DY-sign(CORSS,DY)
	    DZ=DZ-sign(CORSS,DZ)
	  end if
          RR2=DX**2+DY**2+DZ**2
	  if(LMNBL)then
	    NMI=NSITS(ITYP) 
	    NMJ=NSITS(ITYPE(J))
	    if(NMI.gt.MMOL.and.NMJ.gt.MMOL)then
	      RRM=RR2
	      go to 130
	    else if(NMI.gt.MMOL.and.NMJ.le.MMOL)then
	        DX=SX(I)-X(JMOL)
	        DY=SY(I)-Y(JMOL)
	        DZ=SZ(I)-Z(JMOL)
	      else if(NMJ.gt.MMOL)then
	        DX=X(IMOL)-SX(J)
	        DY=Y(IMOL)-SY(J)
	        DZ=Z(IMOL)-SZ(J)
	      else 
	        DX=X(IMOL)-X(JMOL)
	        DY=Y(IMOL)-Y(JMOL)
	        DZ=Z(IMOL)-Z(JMOL)
	      end if
*	      call PBC (DX,DY,DZ)
	      if(DX.gt. HBOXL)DX=DX-BOXL
            if(DX.lt.-HBOXL)DX=DX+BOXL
            if(DZ.gt. HBOZL)DZ=DZ-BOZL
            if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	      if(LHEX)then
C  hexagonal periodic cell)
	        XY=BOXYC*DX
            if(DY.gt.BOXY3-XY.and.(DX.gt.0..or.DY.gt.2.*BOXY3-XY))then
                DY=DY-BOYL
                DX=DX-HBOXL
	        XY=BOXYC*DX
              end if
          if(DY.lt.-BOXY3-XY.and.(DX.lt.0..or.DY.lt.-2.*BOXY3-XY))then
                DY=DY+BOYL
                DX=DX+HBOXL
	        XY=BOXYC*DX
              end if
              if(DY.gt.BOXY3+XY)then
                DY=DY-BOYL
                DX=DX+HBOXL
	        XY=BOXYC*DX
              end if
              if(DY.lt.-BOXY3+XY)then
                DY=DY+BOYL
                DX=DX-HBOXL
              end if
	      else                   
C  rectangular cell
              if(DY.gt. HBOYL)DY=DY-BOYL
              if(DY.lt.-HBOYL)DY=DY+BOYL
	      end if
	      if(LOCT)then
C  truncated octahedron abs(x)+abs(y)+abs(z) < 0.75*BOXL
	        CORSS=HBOXL*int((4./3.)*(abs(DX)+abs(DY)+abs(DZ))/BOXL)
	        DX=DX-sign(CORSS,DX)
	        DY=DY-sign(CORSS,DY)
	        DZ=DZ-sign(CORSS,DZ)
	      end if                  
C  square distance betwen COM of small molecules
	      RRM=DX**2+DY**2+DZ**2 
	    else
	      RRM=RR2              
	    end if
 130	    if(IMOL.eq.JMOL)then 
              if(.not.L15NB(ITYP))go to 135
	      if(RR2.le.RMX2)then
C   Check that the atoms are non-bound
	      if(NNBB(I).le.NNBB(J))then
		do K=1,NNBB(I) 
C   1-3 bound
		  if(INBB(I,K).eq.J)go to 135
                if(INBB(I,K).eq.-J.and.L14NB(ITYP))then
C   1-4 bound                                     
	            MBSH=MBSH+1
	            if(MBSH.gt.NBSMAX)write(*,*)
     +        '!!! list of closest neigbours overfilled. atoms ',I,J
	            NBS1(MBSH)=-I
		    NBS2(MBSH)=J
	            go to 135
		  end if  
	  	end do
	      else      
		do K=1,NNBB(J) 
C   1-3 bound
		  if(INBB(J,K).eq.I)go to 135
                if(INBB(J,K).eq.-I.and.L14NB(ITYP))then
C   1-4 bound                                     
	            MBSH=MBSH+1
	            if(MBSH.gt.NBSMAX)write(*,*)
     +        '!!! list of closest neigbours overfilled. atoms ',I,J
	            NBS1(MBSH)=-I
		    NBS2(MBSH)=J
	            go to 135
		  end if  
	  	end do
	      end if   ! if(NNB
	      end if   ! RR2.le
	    end if     ! IMOL.eq.JMOL
	    if(RRM.le.SHORT2)then
	      MBSH=MBSH+1
	      if(MBSH.gt.NBSMAX)then
              write(*,*)
     +'!!! list of closest neigbours overfilled, atoms ',I,J,TASKID   
	        call FINAL
	      end if
	      NBS1(MBSH)=J
	      NBS2(MBSH)=I
	    else if(RRM.le.RSSQ)then
	      MBLN=MBLN+1
	      if(MBLN.gt.NBLMAX)then
                if(IOVER.eq.0)write(*,*)
     +     '!!! list of far neigbours overfilled',TASKID
                IOVER=1
              end if
              if(IOVER.eq.0)then
	        NBL1(MBLN)=J
	        NBL2(MBLN)=I
              end if
	    end if
C   RDF calculation
 135      IF(LGRC.and.RR2.le.RDFCUT2) THEN			!=====> LGR
	      ISB	= NSITE(J)
	      JSB	= NSITE(I) 
	      if(NGRI(ISB).le.NGRI(JSB))then
		do JS=1,NGRI(ISB) 
		  if(IGRI(ISB,JS).eq.JSB)then 
		    IP = MGRI(ISB,JS)
	            go to 137  
		  end if
		end do
	      else
		do JS=1,NGRI(JSB)
		  if(IGRI(JSB,JS).eq.ISB)then
		    IP = MGRI(JSB,JS)
	            go to 137
		  end if  
		end do         
	      end if  
	      go to 138
 137	      R1	= sqrt(RR2)
            GG        = R1*DRDFI
	      MM        = IDINT(GG)+1
            if(MM.le.MAX)IRDF(MM,IP) = IRDF(MM,IP)+1
	    end if ! of LGRC
 138	    continue 
	  end do
	end do 
	if(LGRC)NRDF=NRDF+1 
*
*  15.3 Check that list of neighbours is not overfilled 
*
        if(IOVER.ne.0)then
          write(*,*)' Insrease NBLMAX in prcm.h for this job to ',MBLN
          I=MBLN/NTOT+1
          write(*,*)' (change  to:  NBLMAX =',I,'* NTOT)'
          call FINAL
        end if
        timev	= timev+cputime(timev0)
	if(IPRINN.lt.7)return
	write(*,*)' Task ',TASKID,
     +'total num.of interactions: ',MBSH,' and ',MBLN
	return
	end
	
