C========================================
C
C     MDynaMix v.4.3
C
*     PART 8
*
*     File aver.f
*     -------------
*
C     This file contains subroutines responsible for collecting averages 
C
*============= GETAVR =================================================
*
      SUBROUTINE GETAVR(NNN)
*
      include "prcm.h"
      external SHIFT
      character*32000 STR
      dimension AUX(NRQS),AUX2(4,LHIST)     
      data INIT/0/
CC
*
* !!! Attention! Any corrections should be done accordingly in 
*                parts 1,2 and 3 (labels 1000,2000,3000)
*  Loacl constants
      TOKJ           = ENERF*1.D-3
      FNOPI          = 1.D0/DFLOAT(NOP)
*
      GO TO (1000,2000,3000),NNN
*
      STOP '!!! NNN OUT OF ORDER IN GETAVR !'
*
*    1.1  Accumulate current values of physical properties
*    -----------------------------------------------------
 1000 CONTINUE
      if(TASKID.ne.MAST)return
      IAV            = 2
      if(ISTAV.gt.0)then
*
*    1.1.1  Covalent bonds:
*    (miss if ISTAV = 0 )
        M              = 0
        DO ITYP        = 1,NTYPES
          NSP            = NSPEC(ITYP)
          FNSPI          = 1.D0/DFLOAT(NSP)
          NRBS           = NRB(ITYP)
          IF(NRBS.NE.0)  THEN
            FNBI         = 1.D0/DFLOAT(NRBS)
            DO I         = 1,NRBS
              M          = M+1
              IAV        = IAV+1
              BBM        = BB(M) 
              AV(IAV)    = AV(IAV)+BBM
              AW(IAV)    = AW(IAV)+BBM**2
              IAV        = IAV+1
              EBM        = EB(M)*TOKJ
              AV(IAV)    = AV(IAV)+EBM
              AW(IAV)    = AW(IAV)+EBM**2
              EBOND(ITYP)= EBOND(ITYP)+EBM
            END DO! OF I
          END IF
        END DO! OF ITYP
*Angles:
        M              = 0
        DO ITYP        = 1,NTYPES
          NSP            = NSPEC(ITYP)
          FNSPI          = 1.D0/DFLOAT(NSP)
          NRAS           = NRA(ITYP)
          IF(NRAS.NE.0)  THEN
            FNAI         = 1.D0/DFLOAT(NRAS)
            DO I         = 1,NRAS
              M              = M+1
              IAV            = IAV+1
              AAM            = AA(M)*TODGR
              AV(IAV)        = AV(IAV)+AAM
              AW(IAV)        = AW(IAV)+AAM**2
              IAV            = IAV+1
              EAM            = EA(M)*TOKJ
              AV(IAV)        = AV(IAV)+EAM
              AW(IAV)        = AW(IAV)+EAM**2
              EANGL(ITYP)    = EANGL(ITYP)+EAM
            END DO! OF I
          END IF
        END DO! OF ITYP
*Torsions:
      M              = 0
      DO ITYP        = 1,NTYPES
      NSP            = NSPEC(ITYP)
      FNSPI          = 1.D0/DFLOAT(NSP)
      NRTS           = NRT(ITYP)
      IF(NRTS.NE.0)  THEN
      FNTI           = 1.D0/DFLOAT(NRTS)
      DO I           = 1,NRTS
      M              = M+1
      IAV            = IAV+1
      TTM            = TT(M)*TODGR
      AV(IAV)        = AV(IAV)+TTM
      AW(IAV)        = AW(IAV)+TTM**2
      IAV            = IAV+1
      ETM            = ET(M)*TOKJ
      AV(IAV)        = AV(IAV)+ETM
      AW(IAV)        = AW(IAV)+ETM**2
      ETORS(ITYP)    = ETORS(ITYP)+ETM
      END DO! OF I
      END IF
      END DO! OF ITYP
      end if ! ISTAV.gt.0
*Other energies
      PEL=0.
      DO I           = 1,NTYPES
          NSP          = NSPEC(I)
          FNSPI        = 1.D0/DFLOAT(NSP)
          SELFPE(I)    = SELFPE(I)*TOKJ*FNSPI
          PES14 (I)    = PES14 (I)*TOKJ*FNSPI
          PSR14 (I)    = PSR14 (I)*TOKJ*FNSPI
          IAV          = IAV+1
          AV(IAV)      = AV(IAV)+SELFPE(I)
          AW(IAV)      = AW(IAV)+SELFPE(I)**2
          IAV          = IAV+1
          AV(IAV)      = AV(IAV)+PES14(I)
          AW(IAV)      = AW(IAV)+PES14(I)**2
          IAV          = IAV+1
          AV(IAV)      = AV(IAV)+PSR14(I)
          AW(IAV)      = AW(IAV)+PSR14(I)**2
	  PEL	= PEL+SELFPE(I)+PES14(I)
      END DO
*
*Intermolecular non-bonded
      DO ITYP        =     1,NTYPES
        NSPI           = NSPEC(ITYP)
        DO JTYP        =  ITYP,NTYPES
          NSPJ           = NSPEC(JTYP)
          I              =   MDX(ITYP,JTYP)
                       FNSPI = 1.D0/DFLOAT(NSPI)
          IF(NSPJ.lt.NSPI)FNSPI = 1.D0/DFLOAT(NSPJ)
          POTSEL         = PELRC(I)*TOKJ*FNSPI
          POTESE         = POTES(I)*TOKJ*FNSPI
          POTLJE         = POTLJ(I)*TOKJ*FNSPI
          PEL	= PEL+POTESE
          IAV            = IAV+1
          AV(IAV)        = AV(IAV)+POTSEL
          AW(IAV)        = AW(IAV)+POTSEL**2
          IAV            = IAV+1
          AV(IAV)        = AV(IAV)+POTESE
          AW(IAV)        = AW(IAV)+POTESE**2
          IAV            = IAV+1
          AV(IAV)        = AV(IAV)+POTLJE
          AW(IAV)        = AW(IAV)+POTLJE**2
          POTES(I)       = 0.D0
          POTLJ(I)       = 0.D0
        END DO! OF JTYP
      END DO! OF ITYP
* temperature
      IAV            = IAV+1
      AV(IAV)        = AV(IAV)+TEMP
      AW(IAV)        = AW(IAV)+TEMP**2
* pot.energy
      IAV            = IAV+1
      PEX	= PE*TOKJ*FNOPI
      AV(IAV)        = AV(IAV)+PEX
      AW(IAV)        = AW(IAV)+PEX**2
*  intramolecular energy
      IAV            = IAV+1
      PINTT	= PINT*TOKJ*FNOPI
      AV(IAV)        = AV(IAV)+PINTT
      AW(IAV)        = AW(IAV)+PINTT**2
*  pressure
      IAV	= IAV+1
      PRES= TRYCK
      AV(IAV)	= AV(IAV)+PRES
      AW(IAV)	= AW(IAV)+PRES**2     
      IAV	= IAV+1
      AV(IAV)	= AV(IAV)+TRYCKM
      AW(IAV)	= AW(IAV)+TRYCKM**2
*   Pressure X,Y,Z components
      IAV     = IAV+1
      AV(IAV) = AV(IAV)+TRYCKX
      AW(IAV) = AW(IAV)+TRYCKX**2
      IAV     = IAV+1
      AV(IAV)=AV(IAV)+TRYCKY
      AW(IAV) = AW(IAV)+TRYCKY**2
      IAV     = IAV+1
      AV(IAV)=AV(IAV)+TRYCKZ
      AW(IAV) = AW(IAV)+TRYCKZ**2
*  box sizes
      IAV	= IAV+1
      AV(IAV)	= AV(IAV)+BOXL
      AW(IAV)	= AW(IAV)+BOXL**2
      IAV	= IAV+1
      AV(IAV)	= AV(IAV)+BOYL
      AW(IAV)	= AW(IAV)+BOYL**2
      IAV	= IAV+1
      AV(IAV)	= AV(IAV)+BOZL
      AW(IAV)	= AW(IAV)+BOZL**2
*  translational and rotational kinetic energies
	IAV	= IAV+1
      AV(IAV)	= AV(IAV)+TTR
      AW(IAV)	= AW(IAV)+TTR**2
*
	IAV	= IAV+1
      AV(IAV)	= AV(IAV)+TROT
      AW(IAV)	= AW(IAV)+TROT**2  
*
	IAV	= IAV+1
      AV(IAV)	= AV(IAV)+TINT
      AW(IAV)	= AW(IAV)+TINT**2
*  Temperatures of different species
	do JTYP=1,NTYPES 
          IAV            = IAV+1
          AV(IAV)        = AV(IAV)+TEMPR(JTYP)
          AW(IAV)        = AW(IAV)+TEMPR(JTYP)**2
	end do
*  "Temperature" of volume fluctuations   
      if(LSEP)then
        if(LHEX)then
          TEMPV = FACVT*(SCLX**2+SCLZ**2)/2.
        else
          TEMPV = FACVT*(SCLX**2+SCLY**2+SCLZ**2)/3.
        end if
      else            ! isotropic cell fluctuations
	  TEMPV   = FACVT*SCL**2
      end if
      IAV	= IAV+1
      AV(IAV) = AV(IAV)+TEMPV
      AW(IAV) = AW(IAV)+TEMPV**2
*  Energy absorbed from the thermostat
      IAV	= IAV+1
      AV(IAV) = EABS                    ! energy absorbed on this iteration
      if(INIT.eq.0)then
        EABSS=AW(IAV)
        INIT=1
      end if
*  total averaged quantities
      NNAV	= IAV+1
*  num of configurations 
        NAVT	= NAVT+1
      IF(IAV.GT.NRQS) THEN
        PRINT *,'!!! INCREASE NRQS TO ',IAV
        call FINAL 
      END IF
*
Calculate total kinetic energy:
*
      EPOTB          = 0.D0
      EKINE          = 0.D0
      DO I           = 1,NTYPES
        EKINE          = EKINE+EKIN (I)
        EPOTB          = EPOTB+EBOND(I)+EANGL(I)+ETORS(I)
        EBOND(I)       = 0.D0
        EANGL(I)       = 0.D0
        ETORS(I)       = 0.D0                                    
      END DO
      EKINE          = EKINE*TOKJ*FNOPI
*
      RETURN
*
*    1.2  Calculate intermediate averages and remember
*    -------------------------------------------------
 2000 CONTINUE
*
      NHIST	= NHIST+1
      if(TASKID.ne.MAST)return
      MNSTEP         = MSTEP       
      FNAVSI         = 1.D0/DFLOAT(NAVT)
      TTIM	= (TIM+NSTEP*DT)*1.d12
*
      PRINT "(80('-'))"
      PRINT 
     +"('*** INTERMEDIATE AVERAGE QUANTITIES AFTER ',I8,' STEPS: ')",
     +NAVT*IAVER
      write(*,*)' (after ',TTIM,'ps of simulation, point ',NHIST+1,')'
      write(*,*)
*
      if(IPRINT.ge.6)then
      write(*,*)' CPU time for some procedures:'
      time=timest+cputime(time0)
      write(*,*)' medium-range intermolecular forces:',timel      
      write(*,*)' short-range intermolecular forces: ',times
	write(*,*)' choice of neuigbours:              ',timev   
      write(*,*)' FURIR:                             ',timef      
      write(*,*)' ETERM:                             ',timee      
      write(*,*)' covalent bonds:                    ',timeb      
      write(*,*)' covalent angles:                   ',timea      
      write(*,*)' torsions:                          ',timet      
      write(*,*)' moving:                            ',timeg      
      write(*,*)' data transfer:                     ',timen      
      write(*,*)'-------------------------------------------------'
      write(*,*)' Total time:                        ',time      
      end if
      if(NHIST.gt.LHIST)then
        call SHIFT(HIST,LHIST,NNAV,NRQS)
        NHIST=LHIST
      end if  
      HIST(1,NHIST)=TTIM
      HIST(2,NHIST)=dfloat(NAVT*IAVER)
      IAV            = 2
      if(ISTAV.gt.0)then
      DO ITYP        = 1,NTYPES
      EBONDS	= 0.
      NRBS           = NRB(ITYP)
      IF(NRBS.NE.0)  THEN
        ISHF           = ISADR(ITYP)
        NBBEG          = IADB(ITYP)
        NBEND          = NBBEG+NRBS-1
	  if(IPRINT.ge.5.and.LMOVE(ITYP))then
        PRINT "(/'*** BOND LENGTHS AND ENERGIES FOR ',A6)",NAME(ITYP) 
	write(*,*)'                       bond length            energy' 
	  end if
        DO I           = NBBEG,NBEND
          II             = IB(I)+ISHF
          JJ             = JB(I)+ISHF
          IAV            = IAV+1
          AV1            = AV(IAV)*FNAVSI
          AW1            = DSQRT(DABS(AW(IAV)*FNAVSI-AV1**2))
          IAV            = IAV+1
          AV2            = AV(IAV)*FNAVSI
          EBONDS=EBONDS+AV2
          AW2            = DSQRT(DABS(AW(IAV)*FNAVSI-AV2**2))
	    if(IPRINT.ge.5.and.LMOVE(ITYP))
     +    write(*,'(a4,a3,a4,2(f12.4,a3,f9.4,3x))')
     +    NM(II),' - ',NM(JJ),AV1,'+/-',AW1,AV2,'+/-',AW2
        END DO! OF I
        write(*,'(a,f12.4,a)')'  Total bond energy ',EBONDS,' kj/mol'    
      else
	  IAV=IAV+2*NRBS
      END IF 
      END DO! OF ITYP
*     
      DO ITYP        = 1,NTYPES
	EANGLES=0.
      NRAS           = NRA(ITYP)
      IF(NRAS.NE.0)  THEN
	if(IPRINT.ge.5.and.LMOVE(ITYP))then
      PRINT "(/'***      ANGLES  AND ENERGIES FOR ',A6)",NAME(ITYP)
	write(*,*)'                         angles          energies'
	end if
      NABEG          = IADA(ITYP)
      NAEND          = NABEG+NRAS-1
      ISHF           = ISADR(ITYP)
      DO I           = NABEG,NAEND
      II             = IA(I)+ISHF
      JJ             = JA(I)+ISHF
      KK             = KA(I)+ISHF
      IAV            = IAV+1
      AV1            = AV(IAV)*FNAVSI
      AW1            = DSQRT(DABS(AW(IAV)*FNAVSI-AV1**2))
      IAV            = IAV+1
      AV2            = AV(IAV)*FNAVSI
      AW2            = DSQRT(DABS(AW(IAV)*FNAVSI-AV2**2))    
	EANGLES	=EANGLES+AV2
	if(IPRINT.ge.5.and.LMOVE(ITYP))
     +write(*,'(a4,a1,a4,a1,a4,3x,2(f12.4,a3,f9.4,3x))')
     +NM(II),'-',NM(JJ),'-',NM(KK),AV1,'+/-',AW1,AV2,'+/-',AW2
      END DO! OF I
	write(*,'(a,a,f12.4,a)')'Totan angles energy of typ ',NAME(ITYP),
     +EANGLES,' kJ/M' 
	else
	  IAV=IAV+2*NRAS
      END IF              
      END DO! OF ITYP
*   torsions
      DO ITYP        = 1,NTYPES                       
	ETORT	= 0.
      NRTS           = NRT(ITYP)
      IF(NRTS.NE.0)  THEN
	if(IPRINT.ge.5.and.LMOVE(ITYP))
     +PRINT "(/'*** TORS. ANGLES AND ENERGIES FOR ',A6)",NAME(ITYP)
      NTBEG          = IADT(ITYP)
      NTEND          = NTBEG+NRTS-1
      ISHF           = ISADR(ITYP)
      DO I           = NTBEG,NTEND
      II             = IT(I)+ISHF
      JJ             = JT(I)+ISHF
      KK             = KT(I)+ISHF
      LL             = LT(I)+ISHF
      IAV            = IAV+1
      AV1            = AV(IAV)*FNAVSI
      AW1            = DSQRT(DABS(AW(IAV)*FNAVSI-AV1**2))
      IAV            = IAV+1
      AV2            = AV(IAV)*FNAVSI
      AW2            = DSQRT(DABS(AW(IAV)*FNAVSI-AV2**2))  
	ETORT	= ETORT+AV2
	if(IPRINT.ge.5.and.LMOVE(ITYP))
     +write(*,'(a4,3(a1,a4),3x,2(f12.4,a3,f9.4,3x))')
     +NM(II),'-',NM(JJ),'-',NM(KK),'-',NM(LL),AV1,'+/-',AW1,
     +AV2,'+/-',AW2
      END DO! OF I
	write(*,'(a,a,f12.4,a)')' total torsion energy of type ',
     +NAME(ITYP),ETORT,' kJ/M' 
	else
	  IAV=IAV+2*NRTS
      END IF             
      END DO! OF ITYP   
      end if   ! ISTAV.gt.0
	write(*,*)
*  Intramolecular
      DO I           = 1,NTYPES
      write(*,*)
     +  '*** NON-BONDED INTRAMOLECULAR CONTRIBUTIONS FOR ',NAME(I)
      IAV            = IAV+1
      AV1            = AV(IAV)*FNAVSI
      AW1            = DSQRT(DABS(AW(IAV)*FNAVSI-AV1**2))
      IAV            = IAV+1
      AV2            = AV(IAV)*FNAVSI
      AW2            = DSQRT(DABS(AW(IAV)*FNAVSI-AV2**2))
      IAV            = IAV+1
      AV3            = AV(IAV)*FNAVSI
      AW3            = DSQRT(DABS(AW(IAV)*FNAVSI-AV3**2))
      PRINT "(15X,' Electrostatic :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV2,AW2
      PRINT "(15X,' Lennard-Jones :',F12.4,' kj/mol  +/- ',F12.4/)",
     +AV3,AW3
*      PRINT "(15X,' Self-electrost:',F12.4,' kj/mol  +/- ',F12.4/)" ,
*     +AV1,AW1
      END DO
*
      DO ITYP        =    1,NTYPES
      DO JTYP        = ITYP,NTYPES
      I              =  MDX(ITYP,JTYP)
      PRINT "('*** INTERMOLECULAR ENERGIES: ',A6,' - ',A6,' ***')"
     X,NAME(ITYP),NAME(JTYP)
      IAV            = IAV+1
      AV1            = AV(IAV)*FNAVSI
      AW1            = DSQRT(DABS(AW(IAV)*FNAVSI-AV1**2))
      IAV            = IAV+1
      AV2            = AV(IAV)*FNAVSI
      AW2            = DSQRT(DABS(AW(IAV)*FNAVSI-AV2**2))
      IAV            = IAV+1
      AV3            = AV(IAV)*FNAVSI
      AW3            = DSQRT(DABS(AW(IAV)*FNAVSI-AV3**2))
      PRINT "(15X,'Electrostatic  :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV2,AW2
      PRINT "(15X,'Lennard-Jones  :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV3,AW3
      PRINT "(15X,'Out cutoff corr:',F12.4,' kj/mol  +/- ',F12.4/)" ,
     +AV1,AW1
      END DO! OF JTYP
      END DO! OF ITYP
      IAV            = IAV+1
      AV4            = AV(IAV)*FNAVSI
      AW4            = DSQRT(DABS(AW(IAV)*FNAVSI-AV4**2))
      PRINT "(22X,' TEMP  :',F12.4,' K       +/- ',F12.4/)",AV4,AW4
      IAV            = IAV+1
      AV4            = AV(IAV)*FNAVSI
      AW4            = DSQRT(DABS(AW(IAV)*FNAVSI-AV4**2))
      PRINT "(22X,' EPnon-bond:',F12.4,' kj/mol  +/- ',F12.4)",AV4,AW4
      IAV            = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV4**2))
      PRINT "(22X,' EPbonded :',F12.4,' kj/mol  +/- ',F12.4)",AV5,AW5
      AV3	= AV4+AV5
      AW3	= sqrt(AW5**2+AW4**2)
      PRINT "(22X,' EPtotal  :',F12.4,' kj/mol  +/- ',F12.4/)",AV3,AW3  
      IAV            = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' PRES(A) :',F12.4,' ATM     +/- ',F12.4)",AV5,AW5
      IAV            = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' PRES(M) :',F12.4,' ATM     +/- ',F12.4/)",AV5,AW5
      IAV      = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' PRES(X) :',F12.4,' ATM     +/- ',F12.4)",AV5,AW5
      IAV      = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' PRES(Y) :',F12.4,' ATM     +/- ',F12.4)",AV5,AW5
      IAV      = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' PRES(Z) :',F12.4,' ATM     +/- ',F12.4/)",AV5,AW5
      IAV            = IAV+1
      AVX            = AV(IAV)*FNAVSI
      AWX            = DSQRT(DABS(AW(IAV)*FNAVSI-AVX**2))
      IAV            = IAV+1
      AVY            = AV(IAV)*FNAVSI
      AWY            = DSQRT(DABS(AW(IAV)*FNAVSI-AVY**2))
      IAV            = IAV+1
      AVZ            = AV(IAV)*FNAVSI
      AWZ            = DSQRT(DABS(AW(IAV)*FNAVSI-AVZ**2))
      AV5            = AVX*AVY*AVZ
      AW5            = AWX*AVY*AVZ+AWY*AVX*AVZ+AWZ*AVY*AVZ   
	if(ICELL.eq.1)then
	  AV5=AV5/2
	  AW5=AW5/2
	end if
      DENS	= TOTMAS*1.d-3/(AV5*UNITL**3)
	DENSE	= TOTMAS*AW5*1.d-3/(AV5**2*UNITL**3)
	PRINT "(22X,' DENSITY :',F12.4,' g/cm**3 +/- ',F12.4)",DENS,DENSE
      PRINT "(22X,' BOX-X   :',F12.4,' A       +/- ',F12.4/)",AVX,AWX
      PRINT "(22X,' BOX-Y   :',F12.4,' A       +/- ',F12.4/)",AVY,AWY
      PRINT "(22X,' BOX-Z   :',F12.4,' A       +/- ',F12.4/)",AVZ,AWZ
      write(*,*)
	IAV            = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' Trans.temp:',F12.4,' K   +/- ',F12.4)",AV5,AW5
      IAV            = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' Rot.temp:  ',F12.4,' K   +/- ',F12.4)",AV5,AW5
      IAV            = IAV+1
      AV5            = AV(IAV)*FNAVSI
      AW5            = DSQRT(DABS(AW(IAV)*FNAVSI-AV5**2))
      PRINT "(22X,' Int.temp.: ',F12.4,' K   +/- ',F12.4)",AV5,AW5
      do ITYP=1,NTYPES
        IAV            = IAV+1
        AV1             = AV(IAV)*FNAVSI
        AW1            = DSQRT(DABS(AW(IAV)*FNAVSI-AV1**2)*FNAVSI)
	  write(*,'(22x,a7,a6,a3,f12.4,a,f12.4)')
     +  ' Temp (',NAME(ITYP),')  ',AV1,' K  +/-',AW1 
      end do                  
      IAV	= IAV + 1
      AV2	= AV(IAV)*FNAVSI
      AW2            = DSQRT(DABS(AW(IAV)*FNAVSI-AV2**2))
      PRINT "(22X,' Vol.temp.: ',F12.4,' K   +/- ',F12.4)",AV2,AW2
      write(*,*)
      IAV	= IAV + 1
      IAVA = IAV
*   Absorbed energy
      AV2	= AV(IAV)
      if(NHIST.gt.1)then
        FTIM=HIST(1,NHIST-1)
      else
        FTIM=0.
      end if
      EABSB = AV2*PERMOL/(TTIM-FTIM)
      if(INIT.eq.0)then
        EABSS=AW(IAV)
        INIT=1
      end if
      if(IEXT.ne.0)write(*,'(22x,a,g13.6,a)')
     +' Absorbed energy ',EABSB,' kJ/M/ps'
      PRINT "(80('-'))"  
	IAV=IAV+2
      do I=3,IAV
        HIST(I,NHIST)=AV(I)*FNAVSI
      end do
      HIST(IAVA,NHIST)=EABSB
      EABSS = EABSS + EABS
      AW(IAVA) = EABSS
      EABS = 0.
      DO I        = 1,NRQS
        AV(I)       = 0.D0
        AW(I)       = 0.D0
      END DO! OF I
      NAVT=0
*
      RETURN     
 3000 continue
      if(IHIST.gt.NHIST.or.NHIST.le.0.or.IPRINT.lt.2)return
      if(TASKID.ne.MAST)go to 3010
      IBEG=1
      if(IPRINT.lt.7)IBEG=IHIST
      SUMT=0.
      if(NHIST.eq.IHIST)then
        FACD=1.d0/sqrt(NHIST-IHIST+1.d0)
      else
        FACD=1.d0/sqrt(1.*NHIST-IHIST)
      end if
* Calculation of averages
      do I=IHIST,NHIST
        SUMT=SUMT+HIST(2,I)
      end do
      do I=1,3
	do J=1,NHIST
	  AUX2(I,J)=0.
	end do
      end do
      do J=3,NNAV
        AV(J)=0.
        AW(J)=0.
        do I=IHIST,NHIST
          AV(J)=AV(J)+HIST(J,I)*HIST(2,I)
        end do
        AV(J)=AV(J)/SUMT
        do I=IHIST,NHIST
          AW(J)=AW(J)+(AV(J)-HIST(J,I))**2*HIST(2,I)
        end do
        AW(J)=sqrt(AW(J)/SUMT)*FACD
      end do
*
      NST=SUMT
      if(IHIST.gt.1)then
        TM0=HIST(1,IHIST-1)
      else
        TM0=0.
      end if
      PRINT "(80('-'))"
      PRINT 
     +"('#** FINAL AVERAGE QUANTITIES AFTER ',I8,' STEPS: ')",NST
      write(*,*)'  from ',TM0,'ps   to ',HIST(1,NHIST),'ps'
      write(*,*)
*     
 3010 time=timest+cputime(time0)
      if(IPRINN.ge.7)then
        write(*,*)' CPU time for some procedures:'
        write(*,*)' medium-range intermolecular forces:',timel,TASKID      
        write(*,*)' short-range intermolecular forces: ',times,TASKID
	write(*,*)' choice of neuigbours:              ',timev,TASKID   
        write(*,*)' FURIR:                             ',timef,TASKID      
        write(*,*)' ETERM:                             ',timee,TASKID      
        write(*,*)' covalent bonds:                    ',timeb,TASKID      
        write(*,*)' covalent angles:                   ',timea,TASKID      
        write(*,*)' torsions:                          ',timet,TASKID      
        write(*,*)' moving:                            ',timeg,TASKID      
        write(*,*)' data transfer:                     ',timen,TASKID      
        write(*,*)'-------------------------------------------------'
        write(*,*)' Total time:                        ',time,TASKID
	write(*,*)
      else if(IPRINT.ge.5)then
        write(*,*)' CPU time for some procedures:'
        write(*,*)' medium-range intermolecular forces:',timel      
        write(*,*)' short-range intermolecular forces: ',times
	write(*,*)' choice of neuigbours:              ',timev   
        write(*,*)' FURIR:                             ',timef     
        write(*,*)' ETERM:                             ',timee      
        write(*,*)' covalent bonds:                    ',timeb      
        write(*,*)' covalent angles:                   ',timea      
        write(*,*)' torsions:                          ',timet      
        write(*,*)' moving:                            ',timeg      
        write(*,*)' data transfer:                     ',timen      
        write(*,*)'-------------------------------------------------'
        write(*,*)' Total time:                        ',time
        write(*,*)
      end if
      if(TASKID.ne.MAST)return
      IAV            = 2
      if(ISTAV.gt.0)then
	write(str,'(16(500x,500x,500x,500x))')   
      STR(1:24)='#  Time     Naver       '
      DO ITYP        = 1,NTYPES
        EBONDS	= 0.
        NRBS           = NRB(ITYP)
        IF(NRBS.NE.0)  THEN
          ISHF           = ISADR(ITYP)
          NBBEG          = IADB(ITYP)
          NBEND          = NBBEG+NRBS-1
          IC0=24
          IAV0	= IAV+1 
	    if(IPRINT.ge.5)then
	    write(*,*)'-----------------------------------------------'
          PRINT "('*** BOND LENGTHS AND ENERGIES FOR ',A6)",NAME(ITYP)
	    end if
            DO I           = NBBEG,NBEND
            II             = IB(I)+ISHF
            JJ             = JB(I)+ISHF 
            IAV            = IAV+1
	      AV1	= AV(IAV)
	      AW1	= AW(IAV)
            IAV            = IAV+1
            AV2            = AV(IAV)
            EBONDS=EBONDS+AV2
            AW2            = AW(IAV)
	      if(IPRINT.ge.5)write(*,'(a4,a3,a4,2(f12.4,a3,f9.4,3x))')
     +      NM(II),' - ',NM(JJ),AV1,'+/-',AW1,AV2,'+/-',AW2
            if(IC0.le.31990)STR(IC0:IC0+5)=NM(II)(1:2)//'-'//NM(JJ)(1:2)
            IC0=IC0+17
          END DO! OF I
          if(IPRINT.ge.4)write(*,'(a,a,2x,f12.4,a)')
     +    '###  Total bond energy for ',NAME(ITYP),EBONDS,' kj/mol'
          if(IPRINT.ge.6.and.LMOVE(ITYP))then
	      IF(IC0.GT.31990)IC0=31990
            write(*,*)  
            write(*,*)' History of bonds'
	      ISSB=1
	      ISSE=2488
 3112	      if(ISSE.gt.IC0)ISSE=IC0
	      write(*,*)STR(ISSB:ISSE)
	      if(ISSE.lt.IC0)then
	        ISSB=ISSE+1
	  	ISSE=ISSB+2464
	        go to 3112
	      end if     
	      do J=IBEG,NHIST 
          if(J.eq.IHIST)write(*,*)'-----------------------------'
	        write(*,3333)HIST(1,J),HIST(2,J),(HIST(I,J),I=IAV0,IAV)
	      end do             
	    end if
	    do J=IBEG,NHIST 
	      do I=IAV0+1,IAV,2
	        AUX2(1,J)=AUX2(1,J)+HIST(I,J)*NSPEC(ITYP)
	      end do
	    end do
        END IF                                        

      END DO! OF ITYP
 3333 format(f9.2,f7.0,145(1x,f8.4,f8.3))
*
	write(str,'(16(500x,500x,500x,500x))')   
      STR(1:24)='#  Time     Naver       '
	DO ITYP        = 1,NTYPES
	EBONDS=0.
      NRAS           = NRA(ITYP)
      IF(NRAS.NE.0)  THEN
      NABEG          = IADA(ITYP)
      NAEND          = NABEG+NRAS-1
      ISHF           = ISADR(ITYP)
      IC0=22
      IAV0	= IAV+1 
	if(IPRINT.ge.5)then
      PRINT "(/'***      ANGLES  AND ENERGIES FOR ',A6)",NAME(ITYP)
      write(*,*)'-----------------------------------------------'   
	end if
      DO I           = NABEG,NAEND
      II             = IA(I)+ISHF
      JJ             = JA(I)+ISHF
      KK             = KA(I)+ISHF
      IAV            = IAV+1
      AV1            = AV(IAV)
      AW1            = AW(IAV)
      IAV            = IAV+1
      AV2            = AV(IAV)
      AW2            = AW(IAV)
      EBONDS=EBONDS+AV2
	if(IPRINT.ge.5)
     +write(*,'(a4,a1,a4,a1,a4,3x,2(f12.4,a3,f9.4,3x))')
     +NM(II),'-',NM(JJ),'-',NM(KK),AV1,'+/-',AW1,AV2,'+/-',AW2
      if(IC0.le.31985)STR(IC0:IC0+14)=NM(II)//'-'//NM(JJ)//'-'//NM(KK)
      IC0=IC0+17
      END DO! OF I
      if(IPRINT.ge.4)write(*,'(a,a,2x,f12.4,a)')
     +    '###  Total angles energy for ',NAME(ITYP),EBONDS,' kj/mol'
      if(IPRINT.ge.6.and.LMOVE(ITYP))then
	      IF(IC0.GT.31990)IC0=31990
            write(*,*)
            write(*,*)' History of angless'
	      ISSB=1
	      ISSE=2488
 3113	      if(ISSE.ge.IC0)ISSE=IC0
	      write(*,*)STR(ISSB:ISSE)
	      if(ISSE.lt.IC0)then
	        ISSB=ISSE+1
	        ISSE=ISSB+2464
	        go to 3113
	      end if
	      do J=IBEG,NHIST
              if(J.eq.IHIST)write(*,*)'-----------------------------'
	        write(*,3334)HIST(1,J),HIST(2,J),(HIST(I,J),I=IAV0,IAV)
	      end do
	    end if
	    do J=IBEG,NHIST
	      do I=IAV0+1,IAV,2
		AUX2(2,J)=AUX2(2,J)+HIST(I,J)*NSPEC(ITYP)
	      end do
	    end do
 	  END IF
      END DO! OF ITYP
 3334 format(f9.2,f7.0,145(1x,f8.2,f8.3))
*
	write(str,'(16(500x,500x,500x,500x))')   
      STR(1:24)='#  Time     Naver       '
	DO ITYP        = 1,NTYPES
      NRTS           = NRT(ITYP)
	EBONDS=0.
      IF(NRTS.NE.0)  THEN 
	if(IPRINT.ge.5)then
        PRINT "(/'*** DIHEDRAL ANGLES AND ENERGIES FOR ',A6)",NAME(ITYP)
        write(*,*)'-----------------------------------------------'
      end if
	NTBEG          = IADT(ITYP)
      NTEND          = NTBEG+NRTS-1
      ISHF           = ISADR(ITYP)
      IC0=22
      IAV0	= IAV+1
      DO I           = NTBEG,NTEND
      II             = IT(I)+ISHF
      JJ             = JT(I)+ISHF
      KK             = KT(I)+ISHF
      LL             = LT(I)+ISHF
      IAV            = IAV+1
      AV1            = AV(IAV)
      AW1            = AW(IAV)
      IAV            = IAV+1
      AV2            = AV(IAV)
      AW2            = AW(IAV)
      EBONDS=EBONDS+AV2
      if(IC0.le.31985)STR(IC0:IC0+14)=
     +NM(II)(1:3)//'-'//NM(JJ)(1:3)//'-'//NM(KK)(1:3)//NM(LL)(1:3)
	if(IPRINT.ge.5)write(*,'(a4,3(a1,a4),3x,2(f12.4,a3,f9.4,3x))')
     +NM(II),'-',NM(JJ),'-',NM(KK),'-',NM(LL),AV1,'+/-',AW1,
     +AV2,'+/-',AW2
      IC0=IC0+17
      END DO! OF I
      if(IPRINT.ge.4)write(*,'(a,a,2x,f12.4,a)')
     +    '###  Total dihedral energy for ',NAME(ITYP),EBONDS,' kj/mol'
      if(IPRINT.ge.6.and.LMOVE(ITYP))then
	      IF(IC0.GT.31982)IC0=31982
            write(*,*)
            write(*,*)' History of dihedrals'
	      ISSB=1
 	      ISSE=2486
 3114	      if(ISSE.ge.IC0)ISSE=IC0
	      write(*,*)STR(ISSB:ISSE)
	      if(ISSE.lt.IC0)then
	        ISSB=ISSE+1
		ISSE=ISSB+2464
	        go to 3114
	      end if
	      do J=IBEG,NHIST
          if(J.eq.IHIST)write(*,*)'-----------------------------'
	        write(*,3334)HIST(1,J),HIST(2,J),(HIST(I,J),I=IAV0,IAV)
	      end do
	    end if
	    do J=IBEG,NHIST
	      do I=IAV0+1,IAV,2
		AUX2(3,J)=AUX2(3,J)+HIST(I,J)*NSPEC(ITYP)
	      end do
	    end do
 	  END IF
        END DO! OF ITYP 
*
	do I=1,3
	  do J=IBEG,NHIST
	    AUX2(I,J)=AUX2(I,J)/NOP
	  end do
	end do
      end if   ! if (ISTAV.gt.0
*
      write(str,'(16(500x,500x,500x,500x))')   
      STR(1:24)='#  Time     Naver       '
	IC0=18
      IAV0	= IAV+1
      write(*,*)'---------------------------------------------------'
      write(*,*)
      ENONB=0.
      DO I           = 1,NTYPES
      write(*,*)
     +'*** NON-BONDED INTRAMOLECULAR CONTRIBUTIONS FOR ',NAME(I)
      IAV            = IAV+1
      AV1            = AV(IAV)
      AW1            = AW(IAV)
      IAV            = IAV+1
      AV2            = AV(IAV)
      AW2            = AW(IAV)
      IAV            = IAV+1
      AV3            = AV(IAV)
      AW3            = AW(IAV)
      PRINT "(15X,'Electrostatic  :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV2,AW2
      PRINT "(15X,'Lennard-Jones  :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV3,AW3
      ENONB=ENONB+(AV2+AV3)*NSPEC(I)
      if(IPRINT.ge.7)PRINT 
     +"(15X,'Self-inter.Ew:',F12.4,' kj/mol  +/- ',F12.4/)" ,AV1,AW1
	write(STR(IC0:IC0+1),'(I2)')I
      if(IC0.le.2475)STR(IC0+2:IC0+18)=': Self  Electro   LJ  '
      IC0=IC0+25
      END DO
	if(ISTAV.gt.0.and.IPRINT.ge.7)then
	      IF(IC0.GT.2490)IC0=2490
            write(*,*)
            write(*,*)' History of intramolecular energy contributions'
	      write(*,*)STR(1:IC0)
	      do J=IBEG,NHIST
          if(J.eq.IHIST)write(*,*)'-----------------------------'
	        write(*,3335)HIST(1,J),HIST(2,J),(HIST(I,J),I=IAV0,IAV)
	      end do
	end if
 3335 format(f9.2,f7.0,200(1x,3f8.2))
*
	write(str(25:2496),'(100x,100x,100x,100x,71x),4(250x,250x)')   
	IC0=22
      IAV0	= IAV+1
	write(*,*)
      write(*,*)'-----------------------------------------------'
      DO ITYP        =    1,NTYPES
      DO JTYP        = ITYP,NTYPES
      I              =  MDX(ITYP,JTYP)
      PRINT "('*** INTERMOLECULAR ENERGIES: ',A6,' - ',A6,' ***')"
     X,NAME(ITYP),NAME(JTYP)
      IAV            = IAV+1
      AV1            = AV(IAV)
      AW1            = AW(IAV)
      IAV            = IAV+1
      AV2            = AV(IAV)
      AW2            = AW(IAV)
      IAV            = IAV+1
      AV3            = AV(IAV)
      AW3            = AW(IAV)
      PRINT "(15X,'Electrostatic  :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV2,AW2
      PRINT "(15X,'Lennard-Jones  :',F12.4,' kj/mol  +/- ',F12.4)",
     +AV3,AW3
      PRINT "(15X,'Out cutoff corr:',F12.4,' kj/mol  +/- ',F12.4/)" ,
     +AV1,AW1
      if(IC0.le.2485)STR(IC0:IC0+12)=NAME(ITYP)//'-'//NAME(JTYP)
      IC0=IC0+17
      END DO! OF JTYP
      END DO! OF ITYP
      if(IPRINT.ge.6)then
	      IF(IC0.GT.2490)IC0=2490
            write(*,*)
            write(*,*)' History of intermolecular energy contributions'
	      write(*,*)STR(1:IC0)
	      do J=IBEG,NHIST
          if(J.eq.IHIST)write(*,*)'-----------------------------'
                IIC=0
                do I=IAV0,IAV
                   if(mod(I-IAV0,3).ne.0)then
                      IIC=IIC+1
                      AUX(IIC)=HIST(I,J)
                   end if
                end do
	        write(*,3336)HIST(1,J),HIST(2,J),(AUX(I),I=1,IIC)
	      end do
	end if
 3336 format(f9.2,f7.0,100(1x,f8.1,f8.2))
*
	IAVF=IAV+4
        if(IPRINT.ge.5.and.NRAA.gt.0.and.ISTAV.ge.0)then  
            write(*,*)
            write(*,*)' History of some energy contributions'
	      STR(18:68)=
     +'   Ebond     Eangl     Etors     '
	      write(*,*)STR(1:70)
	      do J=IBEG,NHIST   
          if(J.eq.IHIST)write(*,*)'-----------------------------'
	    if(J.ge.IHIST.or.IPRINT.ge.6)
     +        write(*,'(f9.2,f7.0,4f10.4,f10.3)')
     +        HIST(1,J),HIST(2,J),(AUX2(I,J),I=1,3)
	      end do
	end if
*
	write(str(25:2496),'(100x,100x,100x,100x,71x),4(250x,250x)')   
	IAV0	= IAV+1
      write(*,*)'-----------------------------------------------'
      write(*,*)' FINAL thermodynamics properties:'
	write(*,*)
      IAV            = IAV+1
      AV1=AV(IAV)
	PRINT "(22X,' TEMP    :',F12.4,' K       +/- ',F12.4/)",
     +AV1,AW(IAV)
	STR(18:30)='     TEMP    '
      IAV            = IAV+1
      AV2=AV(IAV)
      PRINT "(22X,' EPnon-bd:',F12.4,' kj/mol  +/- ',F12.4)",
     +AV2,AW(IAV)
	STR(31:39)=' EPnon-bd'
      IAV            = IAV+1
      AV3=AV(IAV)
      EINTER = AV2 - ENONB/FNOP
      PRINT "(22X,' EPbonded:',F12.4,' kj/mol  +/- ',F12.4)",
     +AV3,AW(IAV)
      PRINT "(22X,' EP-inter:',F12.4,' kj/mol  +/- ',F12.4)",
     +EINTER
      AV4=AV2+AV3
      AW4=sqrt(AW(IAV)**2+AW(IAV-1)**2)
      PRINT "(22X,' EPtotal :',F12.4,' kj/mol  +/- ',F12.4/)",
     +AV4,AW4
	STR(40:48)=' EPintra'
	STR(49:56)=' EPtotal'
      IAV            = IAV+1
      AV5=AV(IAV)
      PRINT "(22X,' PRES(At):',F12.4,' ATM     +/- ',F12.4)",
     +AV5,AW(IAV)
	STR(57:65)=' PRES(A) '
      IAV            = IAV+1
      AV6=AV(IAV)
      PRINT "(22X,' PRES(Ml):',F12.4,' ATM     +/- ',F12.4/)",
     +AV6,AW(IAV)
	STR(66:75)=' PRES(M)  '
      IAV            = IAV+1
      AV10=AV(IAV)
      PRINT "(22X,' PRES(X) :',F12.4,' ATM     +/- ',F12.4)",
     +AV10,AW(IAV)
      IAV            = IAV+1
      AV11=AV(IAV)
      PRINT "(22X,' PRES(Y) :',F12.4,' ATM     +/- ',F12.4)",
     +AV11,AW(IAV)
      IAV            = IAV+1
      AV12=AV(IAV)
      PRINT "(22X,' PRES(Z) :',F12.4,' ATM     +/- ',F12.4/)",
     +AV12,AW(IAV)
      IAV            = IAV+1
      AVX            = AV(IAV)
      AWX            = AW(IAV)
      IAV            = IAV+1
      AVY            = AV(IAV)
      AWY            = AW(IAV)
      IAV            = IAV+1
      AVZ            = AV(IAV)
      AWZ            = AW(IAV)
      AV7            = AVX*AVY*AVZ
      AW7            = AWX*AVY*AVZ+AWY*AVX*AVZ+AWZ*AVX*AVY
	if(ICELL.eq.1)then
	  AV7=AV7/2.
	  AW7=AW7/2.
	end if
      DENS	= TOTMAS*1.d-3/(AV7*UNITL**3)
      DENSE	= TOTMAS*AW7*1.d-3/(AV7**2*UNITL**3)
      PRINT "(22X,' DENSITY :',F12.4,' g/cm**3 +/- ',F12.4/)",DENS,DENSE
      PRINT "(22X,' Box-X   :',F12.4,' A       +/- ',F12.4)",AVX,AWX
      PRINT "(22X,' Box-Y   :',F12.4,' A       +/- ',F12.4)",AVY,AWY
      PRINT "(22X,' Box-Z   :',F12.4,' A       +/- ',F12.4)",AVZ,AWZ
      write(*,*)
      PRINT "(22X,' Trans.temp.:',F12.4,' K  +/- ',F12.4)",
     +AV(IAV+1),AW(IAV+1)
      PRINT "(22X,' Rot.temp.:  ',F12.4,' K  +/- ',F12.4)",
     +AV(IAV+2),AW(IAV+2)
      PRINT "(22X,' Inter.temp.:',F12.4,' K  +/- ',F12.4/)",
     +AV(IAV+3),AW(IAV+3)
      IAV=IAV+3
	do ITYP=1,NTYPES
        IAV            = IAV+1
	  write(*,'(22x,a7,a6,a3,f12.4,a,f12.4)')
     +  ' Temp (',NAME(ITYP),')  ',AV(IAV),' K  +/-',AW(IAV) 
	end do
	PRINT "(80('-'))" 
	if(IPRINT.lt.2)go to 9999
      write(*,*)
      write(*,*)' History of global thermodynamic properties'
	write(*,*)STR(1:80)
	do J=1,NHIST            
          if(J.eq.IHIST)write(*,*)'-----------------------------'
	  if(J.ge.IHIST.or.IPRINT.ge.5)
     +  write(*,'(i3,f9.2,f8.0,f10.2,3f9.3,2f9.1)')J,
     +  HIST(1,J),HIST(2,J),HIST(IAV0,J),HIST(IAV0+1,J),HIST(IAV0+2,J),
     +  HIST(IAV0+1,J)+HIST(IAV0+2,J), HIST(IAV0+3,J),HIST(IAV0+4,J)  
	end do             
      write(*,*)'----------------------------------------------------'
      write(*,'(21x,3f9.3,2f9.2,2f9.1)')AV1,AV2,AV3,AV4,AV5,AV6
*     
	if(IPRINT.lt.4)go to 9999
      write(*,*)' Pressure - Density - Box Size'
	STR(23:79)=
     +' Dens    BoxX     BoxY     BoxZ    Px      Py     Pz'
	write(*,*)STR(1:79)
	do J=1,NHIST
          if(J.eq.IHIST)write(*,*)'-----------------------------'
          IAV=IAV0+5 
	  VOL=HIST(IAV+3,J)*HIST(IAV+4,J)*HIST(IAV+5,J)      
	  if(ICELL.eq.1)VOL=VOL/2.
	  if(J.ge.IHIST.or.IPRINT.ge.6)
     +  write(*,'(I3,f9.2,f8.0,f9.4,3f8.2,3f8.0)')
     +  J,HIST(1,J),HIST(2,J),TOTMAS*1.d-3/(VOL*UNITL**3),
     +  (HIST(I,J),I=IAV+3,IAV+5),(HIST(I,J),I=IAV,IAV+2)
	end do 
      write(*,*)'---------------------------------------------------'
      write(*,'(19x,f9.4,3f8.2,3f8.0)')DENS,AVX,AVY,AVZ,AV10,AV11,AV12
      VOL = AV7
      if(ICELL.eq.1)VOL=VOL/2.
*
      write(*,*)
      write(*,*)' Fractional temperatures'
	STR(23:44)=' Ttrn    Trot     Tint      '
	do ITYP=1,NTYPES
	  IC1=40+9*ITYP
	  STR(IC1:IC1+5)=NAME(ITYP)
	end do    
	IAV=IAV+6
	STR(IC1+6:IC1+12)='    vol'
	write(*,*)STR(1:IC1+12)
	do J=IBEG,NHIST
          if(J.eq.IHIST)write(*,*)'-----------------------------'
	  if(J.ge.IHIST.or.IPRINT.ge.6)
     +  write(*,'(I3,f9.2,f8.0,f9.4,f9.2,20f9.3)')
     +  J,HIST(1,J),HIST(2,J),(HIST(I,J),I=IAV,IAV+3+NTYPES)
	end do
        write(*,*)'----------------------------------------------------'
	write(*,'(20x,f9.4,f9.2,20f9.3)')(AV(J),J=IAV,IAV+3+NTYPES)
*    Absorbed energy
      IAV = IAV+4+NTYPES
      EABSB = AV(IAV)
      DEAB  = AW(IAV)
      if(IEXT.ne.0)write(*,'(a,g13.6,a,g13.6)')
     &' Absorbed energy ',EABSB,' kJ/M/ps    +/-',DEAB
*      
      DO I        = 1,NRQS
        AV(I)       = 0.D0
        AW(I)       = 0.D0
      END DO! OF I
      NAVT=0
*  
      IF(.NOT.LNVT) THEN
        PRINT *,'... TEMP SCALED:'
        DO ITYP = 1,NTYPES
          PRINT *,'...',NRTSC(ITYP),'  TIMES IN'
     X,NSTEPS,'  STEPS - FOR TYPE: ',ITYP
        END DO! OF ITYP
      END IF
*
 9999 RETURN
      END 
*
*================== SHIFT =====================================
*
	subroutine SHIFT(HIST,LHIST,NNAV,IN)
	real*8 HIST(IN,LHIST)
	do I=2,LHIST
	  do J=1,NNAV 
	    HIST(J,I-1)=HIST(J,I)
	  end do
	end do
	return
	end 
*
*
*================= RDFINP =============================================
*
      SUBROUTINE RDFINP
      include"prcm.h" 
      character*128 str , TAKESTR
      integer IAUX(NS)
*     
      if(TASKID.eq.MAST)then
        PRINT *
        PRINT *,'----------------------------------------'
        PRINT *,'*** PREPARATION FOR RDF CALCULATIONS ***'
        PRINT *,'----------------------------------------'
      end if
      if(TASKID.eq.MAST)
     +PRINT "('NR OF PAIRS IN RDF CALCULATIONS: ',I5)",MAXRDF
      IF(MAXRDF.GT.MAXGR) STOP '!!! MAXRDF.GT.MAXGR IN RDFINP'
*
      N=MAXRDF
      if(TASKID.eq.MAST)
     +PRINT "('*** NR OF RADIAL DISTRIBUTION FUNCTIONS: ',I4)",N
*
      NRDF       = 0
      DO I       = 1,MAXGR
	  IIN(I)=0
	  IGR(I)=0
	  JGR(I)=0
        DO J       = 1,MAX
          IRDF(J,I)  = 0
        END DO
      END DO
*     
*  Check RDF list and calculate normalising factors
      DO ITYP    = 1,NTYPES
        NSPI       = NSPEC(ITYP)
        ISB        = ISADR(ITYP)+1
        ISE        = ISADR(ITYP +1)
        DO IS      = ISB,ISE 
	    NNR=NGRI(IS)
	    do K=1,NNR
	      JS = IGRI(IS,K)
	      if(IS.le.JS)then
		N    = MGRI(IS,K)
		JTYP = ITS(JS)
              NSPJ = NSPEC(JTYP)
              NMGR(N)    = NM(IS)//NM(JS)
              IGR (N)    = NSPI
              JGR (N)    = NSPJ
	        if(IS.eq.JS.and.NSPI.ne.1)then
	          IIN(N)=IIN(N)+(NSPI-1)*NSPI/2
	        else
	          IIN(N)=IIN(N)+NSPI*NSPJ
	        end if     
	      end if
	    end do
*
        END DO! OF IS
      END DO! OF ITYP
*
      DRDF     = RDFCUT/DFLOAT(MAX)
      DRDFI    =  1.D0/DRDF
      if(TASKID.eq.MAST)then
        PRINT "('*** ONE RDF SLICE IN ANGSTROM:  ',F8.5)",DRDF
        PRINT *
	end if 
*  Check RDF table
	if(IPRINT.ge.7)then   
	  write(*,*)' Calculate RDF:'  
	  do IS=1,NSITES
	    NNR=NGRI(IS) 
          do KS=1,NNR
	      JS=IGRI(IS,KS)            
	      K=MGRI(IS,KS)
	    write(*,'(a,2i5,3x,a4,a1,a4,i6,3i9)')
     +' sites ',IS,JS,NM(IS),'-',NM(JS),K,IGR(K),JGR(K),IIN(K)/2
	    end do                                     
	  end do
	end if
*
 22   RETURN      
      END                       
*
*================= RDFOUT =============================================
*
      SUBROUTINE RDFOUT(IOUT)
      include"prcm.h"
*
      CHARACTER*12 LABRDF
      DIMENSION GR(MAXS,MAXGR),GRO(MAXS,MAXGR),DF(MAXS)
      integer IRDF1(MAXS,MAXGR)
      DATA JINIT/0/, LABRDF /'rdf:        '/
      save GRO,FNRNR,JINIT,NRNR
*     
      IF(NRDF.EQ.0.and.(.not.LGRST)) RETURN
CM
CM   collecting RDF
CM    
      call ADD_RDF(IRDF,IRDF1,MAXS,MAXGR,MAST)  
      if(TASKID.ne.MAST)RETURN
      if(IOUT.ne.0)then
      PRINT *
      PRINT *,'+++++++++++++++++++++++++++++++++++++++++++++++'
      PRINT *,'***      OUTPUT FROM RDF CALCULATIONS       ***'
      PRINT *,'+++++++++++++++++++++++++++++++++++++++++++++++'
      PRINT *
      PRINT "('*** RDF  COLLECTED  DURING:  ',I6,'   STEPS')",NRDF
      PRINT "('*** NR  OF  RDF TYPES  CALCULATED       ',I6)",MAXRDF
      PRINT "('*** CUT OFF RADIUS DIVIDED TO ',I6,' SLICES')",MAX
      PRINT *
*
      PRINT "(1X,130('='))"
      end if
      FNRDF        = DFLOAT(NRDF)
      if(JINIT.eq.0)then
      NRNR         = 0
      FNRNR        = 0.D0
      DO I         = 1,MAXRDF
      DO J         = 1,MAX
      GRO(J,I)      = 0.D0
      END DO! OF J
      END DO! OF I
      end if
*
      IF(LGRST.and.JINIT.eq.0) THEN
        OPEN (UNIT=20,file=frdf,status='old',FORM='UNFORMATTED',
     +        err=111)
        READ(20,err=111) NRNR
        DO I         = 1,MAXRDF
          READ(20,err=112,end=112) (GRO(J,I),J=1,MAX) 
        END DO
        PRINT "('*** RDF RESTART FOR ',I8,' OLD CONFIGURARTIONS')",NRNR
        FNRNR        = DFLOAT(NRNR)
        close (20)
        JINIT=1
      END IF 
	go to 113
 111  write(*,*)' RDF restart failed. File ',frdf
	write(*,*)' RDF are collected only from this run'
	go to 113
 112	write(*,*)' Invalid RDF restart file ',frdf
	frdf(1:1)='?'
	if(LGDMP)write(*,*)' RDF dump will be in file',frdf 
 113  continue
*
      FNRINV       =  1.D0/(FNRDF+FNRNR)
      if(NRDF.eq.0)then
        FNRDF=1.
      end if
      RDFINC       =  RDFCUT/DFLOAT(MAX)
*
      DO 200 I     = 1,MAXRDF
*
        FNI          = DFLOAT(IGR(I))
        FNJ          = DFLOAT(JGR(I))
	FIJ	     = dfloat(IIN(I)) 
*
      if(IOUT.ne.0)then
        write(*,*)
        write(*,'(a,a8,a,i8)')'#*** THIS PAIR: ',NMGR(I),
     +'   --->  TOTAL NR OF CONFIGURATIONS',NRNR+NRDF
        write(*,*)
        write(*,'(3X,5HR(A) ,5X,3HRDF,3X,7HINT RDF)')
*
      end if
      AVRNRI       =  FNRDF
      DNSNRI       = VOL/FNI
      DNSNRJ       = VOL/FNJ 
      RDFINT       = 0.D0 
      LABRDF(5:12)=NMGR(I)
*
      M            = 0
      DO 100 J     = 1,MAX
        M            = M+1    
        DIST         = RDFINC*(DFLOAT(J)-0.5)
        SHELLV       = 4.D0*PI*RDFINC*(DIST**2+RDFINC**2/12.D0)
        FIRDF        = DFLOAT(IRDF1(J,I))
        RDFINT       = RDFINT+FIRDF/ (AVRNRI*FNI)
        RDFIJ        = VOL*FIRDF/(FIJ*AVRNRI*SHELLV)
*
        GR(J,I)      = (FNRDF*RDFIJ+FNRNR*GRO(J,I))*FNRINV
        RDFIJ        = GR(J,I)
	RDFINT	     = RDFINT+RDFIJ*SHELLV/DNSNRJ
        DF(M)        =  RDFIJ*DIST**2
        FINT         =  4.D0*PI*XINTGR(RDFINC,DF,M)
        FINTI        =  FINT/DNSNRI
        FINTJ        =  FINT/DNSNRJ
*
        IF(RDFINT.NE.0.D0.AND.DIST.LE.RDFCUT.and.IOUT.ne.0)PRINT 
     X"(4(1X,F8.4),10x,a12,I9)",
     +DIST,RDFIJ,FINTI,FINTJ,LABRDF
  100   CONTINUE
  200 CONTINUE
*
      if(IOUT.ne.0)PRINT "(1X,130('='))"       
	NRDFT=NRDF+NRNR
      IF(LGDMP.and.(.not.linput).and.NRDFT.gt.0) THEN
      OPEN (UNIT   = 21,file=frdf,status='unknown',FORM='UNFORMATTED')
      WRITE(21) NRDFT
      DO I         = 1,MAXRDF
      WRITE(21) (GR(J,I),J=1,MAX)
      END DO! OF I
      PRINT "('*** RDF DUMP FOR ',I8,' CONFIGURARTIONS')",NRNR+NRDF
      close (21)
      END IF
*
      RETURN
      END
*
