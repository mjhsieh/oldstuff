*
*   Averaging procedures
*
*============ GETEMP =========================================
*
*  Temperature calculations
      SUBROUTINE GETEMP
      include "mdee.h"       
      TKE         = 0.D0
      TEMP        = 0.D0 
      do ITYP=1,NTYPES
	TEMPR(ITYP)=0.
        I           = ISADDR(ITYP)
        NSS         = NSITS (ITYP)
        NSP         = NSPEC (ITYP)
        NSPNS       = NSP*NSS
        VSQ         = 0.D0
	if(LMOVE(ITYP))then
          DO N      = 1,NSPNS
            I       = I+1
            VSQ     = VSQ+(VX(I)**2+VY(I)**2+VZ(I)**2)*MASSDI(I)
          END DO! OF N
          TKE         = TKE+ 0.5D0*VSQ  
          EKIN(ITYP)  =      0.5D0*VSQ
          TEMPR(ITYP) = EKIN(ITYP)*TFACT(ITYP) 
	end if
      END DO! OF ITYP   
C  This is the main result:
      TEMP        = TKE*CONVET
      RETURN
      END
*
*============= GETAVR =================================================
*
      SUBROUTINE GETAVR(NNN)
*
      include "mdee.h"
	character*500 STR
	dimension PEE(NEMAX),FEN(NEMAX)
     +,TMP1(LHIST),TMP2(LHIST),TINN(NTPS*(NTPS+1)/2)
*
* !!! Attention! Any corrections should be done accordingly in 
*                parts 1 and 3 (labels 1000,3000)
*
      TOKJ           = ENERF*1.D-3
      FNOPI          = 1.D0/DFLOAT(NOP)
*
      GO TO (1000,2000,2000),NNN
*
      STOP '!!! NNN OUT OF ORDER IN GETAVR !'
*
 1000 CONTINUE 
      AV(1,ME)=AV(1,ME)+1
      IAV            = 2
*
*Bonds:
      AVB            = 0.D0
      M              = 0
      DO ITYP        = 1,NTYPES
      NSP            = NSPEC(ITYP)
      FNSPI          = 1.D0/DFLOAT(NSP)
      NRBS           = NRB(ITYP)
      IF(NRBS.NE.0)  THEN
      FNBI           = 1.D0/DFLOAT(NRBS)
      DO I           = 1,NRBS
      M              = M+1
      IAV            = IAV+1
      BBM            = BB(M)*FTIMES 
      AV(IAV,ME)        = AV(IAV,ME)+BBM
      AW(IAV,ME)        = AW(IAV,ME)+BBM**2
      IAV            = IAV+1
      EBM            = EB(M)*FTIMES*TOKJ
      AV(IAV,ME)        = AV(IAV,ME)+EBM
      AW(IAV,ME)        = AW(IAV,ME)+EBM**2
      AVB            = AVB+EBM*FNBI
      EBOND(ITYP)    = EBOND(ITYP)+EBM
      END DO! OF I
      END IF
      END DO! OF ITYP
*Angles:
      AVA            = 0.D0
      M              = 0
      DO ITYP        = 1,NTYPES
      NSP            = NSPEC(ITYP)
      FNSPI          = 1.D0/DFLOAT(NSP)
      NRAS           = NRA(ITYP)
      IF(NRAS.NE.0)  THEN
      FNAI           = 1.D0/DFLOAT(NRAS)
      DO I           = 1,NRAS
      M              = M+1
      IAV            = IAV+1
      AAM            = AA(M)*FTIMES*TODGR
      AV(IAV,ME)        = AV(IAV,ME)+AAM
      AW(IAV,ME)        = AW(IAV,ME)+AAM**2
      IAV            = IAV+1
      EAM            = EA(M)*FTIMES*TOKJ
      AV(IAV,ME)        = AV(IAV,ME)+EAM
      AW(IAV,ME)        = AW(IAV,ME)+EAM**2
      EANGL(ITYP)    = EANGL(ITYP)+EAM
      AVA            = AVA+EAM*FNAI
      END DO! OF I
      END IF
      END DO! OF ITYP
*Torsions:
      AVT            = 0.D0
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
      TTM            = TT(M)*FTIMES*TODGR
      AV(IAV,ME)        = AV(IAV,ME)+TTM
      AW(IAV,ME)        = AW(IAV,ME)+TTM**2
      IAV            = IAV+1
      ETM            = ET(M)*FTIMES*TOKJ
      AV(IAV,ME)        = AV(IAV,ME)+ETM
      AW(IAV,ME)        = AW(IAV,ME)+ETM**2
      ETORS(ITYP)    = ETORS(ITYP)+ETM
      AVT            = AVT+ETM*FNTI
      END DO! OF I
      END IF
      END DO! OF ITYP
*Other energies
	PEL=0.
      DO I           = 1,NTYPES
      NSP            = NSPEC(I)
      FNSPI          = 1.D0/DFLOAT(NSP)
      SELFPE(I)      = SELFPE(I)*FTIMES*TOKJ*FNSPI
      PES14 (I)      = PES14 (I)*FTIMES*TOKJ*FNSPI
      PSR14 (I)      = PSR14 (I)*FTIMES*TOKJ*FNSPI
      IAV            = IAV+1
      AV(IAV,ME)        = AV(IAV,ME)+SELFPE(I)
      AW(IAV,ME)        = AW(IAV,ME)+SELFPE(I)**2
      IAV            = IAV+1
      AV(IAV,ME)        = AV(IAV,ME)+PES14(I)
      AW(IAV,ME)        = AW(IAV,ME)+PES14(I)**2
      IAV            = IAV+1
      AV(IAV,ME)        = AV(IAV,ME)+PSR14(I)
      AW(IAV,ME)        = AW(IAV,ME)+PSR14(I)**2
	PEL	= PEL+SELFPE(I)+PES14(I)
      END DO
*
Calculate total potential energy:
*
      POTEN          = 0.D0
*Intermolecular non-bonded
      DO ITYP        =     1,NTYPES
        NSPI           = NSPEC(ITYP)
        DO JTYP        =  ITYP,NTYPES
          NSPJ           = NSPEC(JTYP)
          I              =   MDX(ITYP,JTYP)
          FNSPI = 1.D0/      DFLOAT(NSPI)
          IF(ITYP.NE.JTYP)FNSPI = 1.D0/DMIN1(DFLOAT(NSPI),DFLOAT(NSPJ))
          POTEN          = POTEN+PELRC(I)+POTES(I)+POTLJ(I)
          POTLRC         = PELRC(I)*TOKJ*FNSPI
          POTESE         = POTES(I)*TOKJ*FNSPI
          POTLJE         = POTLJ(I)*TOKJ*FNSPI
          PEL	= PEL+POTESE
          IAV            = IAV+1
          AV(IAV,ME)        = AV(IAV,ME)+POTESE
          AW(IAV,ME)        = AW(IAV,ME)+POTESE**2
          IAV            = IAV+1
          AV(IAV,ME)        = AV(IAV,ME)+POTLJE
          AW(IAV,ME)        = AW(IAV,ME)+POTLJE**2
          POTES(I)       = 0.D0
          POTLJ(I)       = 0.D0
        END DO! OF JTYP
      END DO! OF ITYP
* temperature
      IAV            = IAV+1
      AV(IAV,ME)        = AV(IAV,ME)+TEMP
      AW(IAV,ME)        = AW(IAV,ME)+TEMP**2
* pot.energy
      IAV            = IAV+1
      PEX	= PE*TOKJ*FNOPI
      AV(IAV,ME)        = AV(IAV,ME)+PEX
      AW(IAV,ME)        = AW(IAV,ME)+PEX**2
*  intramolecular energy
      IAV            = IAV+1
      PINTT	= PINT*TOKJ*FNOPI
      AV(IAV,ME)        = AV(IAV,ME)+PINTT
      AW(IAV,ME)        = AW(IAV,ME)+PINTT**2
      POTEN          = POTEN*TOKJ*FNOPI+PQE*TOKJ*FNOPI
	PEL	= PEL+PQE*TOKJ*FNOPI
*  pressure
      IAV	= IAV+1
      PRES= TRYCK
      AV(IAV,ME)	= AV(IAV,ME)+PRES
      AW(IAV,ME)	= AW(IAV,ME)+PRES**2     
      IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+TRYCKM
      AW(IAV,ME)	= AW(IAV,ME)+TRYCKM**2
*  box sizes
      IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+BOXL
      AW(IAV,ME)	= AW(IAV,ME)+BOXL**2
      IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+BOYL
      AW(IAV,ME)	= AW(IAV,ME)+BOYL**2
      IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+BOZL
      AW(IAV,ME)	= AW(IAV,ME)+BOZL**2
*  translational and rotational kinetic energies
	IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+TTR
      AW(IAV,ME)	= AW(IAV,ME)+TTR**2
*
	IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+TROT
      AW(IAV,ME)	= AW(IAV,ME)+TROT**2  
*
	IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+TINT
      AW(IAV,ME)	= AW(IAV,ME)+TINT**2
*
	IAV	= IAV+1
      AV(IAV,ME)	= AV(IAV,ME)+ENNO
      AW(IAV,ME)	= AW(IAV,ME)+ENNO**2
*
      NNAV	= IAV
	NAVT	= NAVT+1
*
      IF(IAV.GT.NRQS) THEN
      PRINT *,'!!! INCREASE NRQS TO ',IAV
      STOP 
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
      POTEN	= PE+PINT
*
      ETOT           = POTEN+EKINE
	ITE(ME,NHIST)=ITE(ME,NHIST)+IAVER
	HIST(2,NHIST)=HIST(2,NHIST)+dfloat(IAVER)	
*
      RETURN
*
 2000 CONTINUE
*     
	if(NNN.eq.3.and.HIST(2,NHIST).lt.0.5d0)then
	  NHIST=NHIST-1
	  go to 3000
	end if
	do I=1,NE
	  PEE(I)=dfloat(ITE(I,NHIST))/HIST(2,NHIST)
	end do
      if(IPRINT.ge.4)then
        write(*,'(i4,10(10f7.4/))')NHIST,(PEE(I),I=1,NE)
        write(*,'(a,3f11.3,a,f10.5)')' Box: ',BOXL,BOYL,BOZL,
     +             '   Dens: ',TOTMAS*1.d-3/(VOL*UNITL**3)

      end if
	do I=1,NE
	  if(NDEL.gt.0)HIST(I+3,NHIST)=EE(I)/NDEL
	end do
*
	MNSTEP         = MSTEP       
      FNAVSI         = 1.D0/DFLOAT(NAVT)
      TTIM	= (TIM+NSTEP*DT)*1.d12
	HIST(1,NHIST)=TTIM
	NNH=NE+3
*     
	if(NNN.eq.3)go to 3000
        NHIST	= NHIST+1
	do I=1,NNH
	  HIST(I,NHIST)=0.d0
	end do 
	write(*,*)' begin ',NHIST,' iteration'
	if(NHIST.eq.IHIST)then
	  NAVT	= 0
	  write(*,*)' New averaging'
          TIMA	= 0.
	call ZEROAV
        DO I        = 1,NRQS 
	    do J=1,NE
            AV(I,J)       = 0.D0
            AW(I,J)       = 0.D0
	    end do
        END DO! OF I
        NSTEPA=0
        do I=1,NE
          do J=1,NE
            IWALK(I,J)=0
          end do
        end do
      end if
      if(NHIST.gt.LHIST)then
        call ISHIFT(ITE,LHIST,NE,NEMAX)
        call SHIFT(HIST,LHIST,NNH,NRH)
        NHIST=LHIST
      end if
	return  
 3000 continue                                                     
      if(IHIST.gt.NHIST)IHIST=NHIST 
	IBEG=1
	if(IPRINT.lt.6)IBEG=IHIST
      SUMT=0.
      FACD=1.d0/sqrt(NHIST-IHIST+1.d0)
* Calculation of averages
      do I=IHIST,NHIST
        SUMT=SUMT+HIST(2,I)
      end do
*
      NST=SUMT
      if(IHIST.gt.1)then
        TM0=HIST(1,IHIST-1)
      else
        TM0=0.
      end if
	do I=1,NE
	  if(AV(1,I).lt.1.d0)AV(1,I)=1.d0
	  do J=1,NNAV 
	    AW(J,I)=AV(J,I)/AV(1,I)
	  end do
	end do
      PRINT "(80('-'))"
      PRINT 
     +"('#** FINAL AVERAGE QUANTITIES AFTER ',I8,' STEPS: ')",NST
      write(*,*)'  from ',TM0,'ps   to ',HIST(1,NHIST),'ps'
      write(*,*)
*     
	if(IPRINT.ge.6)then
      write(*,*)' CPU time for some procedures:'
      time=timest+cputime(time0)
      write(*,*)' medium-range intermolecular forces:',timel      
      write(*,*)' short-range intermolecular forces: ',times
	write(*,*)' change of subensemble:             ',time2
	write(*,*)' choice of neuigbours:              ',timev   
      write(*,*)' FURIR:                             ',timef      
      write(*,*)' ETERM:                             ',timee      
      write(*,*)' atom-atom intramolecular forces:   ',timen      
      write(*,*)' covalent bonds:                    ',timeb      
      write(*,*)' covalent angles:                   ',timea      
      write(*,*)' torsions:                          ',timet      
      write(*,*)' moving:                            ',timeg      
      write(*,*)'-------------------------------------------------'
      write(*,*)' Total time:                        ',time
	write(*,*)
	end if 
	write(*,*)
	write(*,*)'Acceptance ratio ',float(ICHU)/(float(ICHE)+0.001) 
        write(*,*)'==================================================='
	write(*,*)' DISTRIBUTION OVER SUBENSEMBLES'  
	write(*,*)' #  alpha     ita       prob.      bF     F(kJ/M)',
     +'   P-       P+'
	if(NDEL.eq.0)NDEL=1 
        if(NST.le.0)NST=1
        if(PLOW.le.1.d0/NST)PLOW=1.d0/NST
	do I=1,NE
	  PEE(I)=0.
	  do J=IHIST,NHIST
	    PEE(I)=PEE(I)+dfloat(ITE(I,J))/SUMT
	  end do
	end do
	do I=1,NE
          if(PEE(I).le.PLOW)then
            if(PEE(NE).le.PLOW)then
	      FEN(I)=-(EE(I)-EE(NE))/NDEL                        ! Beta*F
            else
	      FEN(I)=-((EE(I)-EE(NE)+dlog(PLOW/PEE(NE)))/NDEL)    ! Beta*F
            end if
          else
            if(PEE(NE).le.PLOW)then
	      FEN(I)=-((EE(I)-EE(NE)+dlog(PEE(I)/PLOW))/NDEL)    ! Beta*F
            else
	      FEN(I)=-((EE(I)-EE(NE)+dlog(PEE(I)/PEE(NE)))/NDEL)    ! Beta*F
            end if
          end if
          FF=FEN(I)*0.001*ENERF/BETA                         ! F (kJ/M)
          IF(PEE(I).NE.0.) THEN
            if(I.gt.1) PTR0=IWALK(I,I-1)*NCEE/(PEE(I)*NSTEPA)
            if(I.lt.NE)PTR1=IWALK(I,I+1)*NCEE/(PEE(I)*NSTEPA)
          ELSE
            PTR0      = 0.
            PTR1      = 0.
          END IF
	  write(*,'(I3,f9.5,f9.4,f9.5,4f9.4)')
     +  I,EC(I),-EE(I)/NDEL,PEE(I),FEN(I),FF,PTR0,PTR1
	end do
        write(*,*)'==================================================='
	write(*,*)' Distribution history:'
	do I=IHIST,NHIST
          if(HIST(3,I).eq.1.d0)write(*,'(a10,10x,80f8.3)')
     +  ' new n_ : ',(HIST(J,I),J=4,NE+3)
      write(*,'(i4,f8.2,f8.0,80f8.5)')
     +  I,HIST(1,I),HIST(2,I),(ITE(J,I)/HIST(2,I),J=1,NE)
	  TMP1(I)=ITE(1,I)/HIST(2,I)
	  TMP2(I)=ITE(NE,I)/HIST(2,I)
	end do                       
	call DISP(TMP1,PM1,PD1,IHIST,NHIST)
	call DISP(TMP2,PM2,PD2,IHIST,NHIST)
      write(*,*)'-----------------------------------------------'
	write(*,'(a,f10.7,a3,f10.7)')'  P(T) = ',PM1,'+/-',PD1
	write(*,'(a,f10.7,a3,f10.7)')'  P(0) = ',PM2,'+/-',PD2 
	if(PM1.gt.0.d0.and.PM2.gt.0.d0)then
	  FRM=(EE(1)-EE(NE)+dlog(PM1/PM2))/NDEL    ! -Beta*F
        FRM=-FRM*0.001*ENERF/BETA                         ! F (kJ/M)
	  FRD=sqrt((PD1/PM1)**2+(PD2/PM2)**2)*0.001*ENERF/(BETA*NDEL)
	  write(*,'(a,f10.4,a3,f10.4)')'  F    = ',FRM,'+/-',FRD 
	end if
      IAV            = 2
      EBONDS	= 0.
      IC0=0
      STR(1:20)='#     alpha              '
      write(str(21:500),'(100x,100x,100x,100x,80x)')   
      DO ITYP        = 1,NTYPES
        NRBS           = NRB(ITYP)
        IF(NRBS.NE.0)  THEN
          ISHF           = ISADR(ITYP)
          NBBEG          = IADB(ITYP)
          NBEND          = NBBEG+NRBS-1
          IC0=17
          IAV0	= IAV+1
            DO I           = NBBEG,NBEND
            II             = IB(I)+ISHF
            JJ             = JB(I)+ISHF
            IAV            = IAV+2
            if(IC0.le.490)STR(IC0:IC0+5)=NM(II)(1:2)//'-'//NM(JJ)(1:2)
            IC0=IC0+17
          END DO! OF I
          if(IPRINT.ge.6)then
	      IF(IC0.GT.500)IC0=500
            write(*,*)
            write(*,*)' Distribution of bonds.     Type ',ITYP
	      write(*,*)STR(1:IC0)
	      do J=1,NE
	        write(*,3333)J,EC(J),(AW(I,J),I=IAV0,IAV)
	      end do
	    end if                
      END IF
      END DO! OF ITYP
 3333 format(I5,f7.4,100(1x,f8.4,f8.3))
*
      write(str(21:500),'(100x,100x,100x,100x,80x)')   
      DO ITYP        = 1,NTYPES
      NRAS           = NRA(ITYP)
      IF(NRAS.NE.0)  THEN
      NABEG          = IADA(ITYP)
      NAEND          = NABEG+NRAS-1
      ISHF           = ISADR(ITYP)
      IC0=17
      IAV0	= IAV+1
      write(*,*)'-----------------------------------------------'
      DO I           = NABEG,NAEND
      II             = IA(I)+ISHF
      JJ             = JA(I)+ISHF
      KK             = KA(I)+ISHF
      IAV            = IAV+2
      if(IC0.le.485)STR(IC0:IC0+14)=NM(II)//'-'//NM(JJ)//'-'//NM(KK)
      IC0=IC0+17
      END DO! OF I
          if(IPRINT.ge.6)then
	      IF(IC0.GT.500)IC0=500
            write(*,*)
            write(*,*)' Distribution of angles'
	      write(*,*)STR(1:IC0)
	      do J=1,NE
	        write(*,3334)J,EC(J),(AW(I,J),I=IAV0,IAV)
	      end do
	    end if                
 	END IF
      END DO! OF ITYP
 3334 format(I5,f7.4,100(1x,f8.2,f8.3))
*
	write(str(21:500),'(100x,100x,100x,100x,80x)')  
	DO ITYP        = 1,NTYPES
      NRTS           = NRT(ITYP)
      IF(NRTS.NE.0)  THEN
      NTBEG          = IADT(ITYP)
      NTEND          = NTBEG+NRTS-1
      ISHF           = ISADR(ITYP)
      IC0=17
      IAV0	= IAV+1
          write(*,*)'-----------------------------------------------'
      DO I           = NTBEG,NTEND
      II             = IT(I)+ISHF
      JJ             = JT(I)+ISHF
      KK             = KT(I)+ISHF
      LL             = LT(I)+ISHF
      IAV            = IAV+2
      if(IC0.le.485)STR(IC0:IC0+14)=NM(II)//'-'//NM(JJ)//'-'//NM(LL)
      IC0=IC0+17
      END DO! OF I
          if(IPRINT.ge.6)then
	      IF(IC0.GT.500)IC0=500
            write(*,*)
            write(*,*)' Distribution of torsions'
	      write(*,*)STR(1:IC0)
	      do J=1,NE
	        write(*,3333)J,EC(J),(AW(I,J),I=IAV0,IAV)
	      end do
	    end if                
 	END IF
      END DO! OF ITYP
*
      write(str(21:500),'(100x,100x,100x,100x,80x)')   
	IC0=17
      IAV0	= IAV+1
      if(IPRINT.ge.7)write(*,*)'------------------------------------'
      DO I           = 1,NTYPES
        IAV            = IAV+3
        if(IPRINT.ge.7)then
          PRINT "('*** SOME ENERGY CONTRIBUTIONS: *** - TYPE:',I2)",I
	    write(STR(IC0:IC0+1),'(I2)')I
          if(IC0.le.485)STR(IC0+2:IC0+23)=' :SelfPE  PES14  PLJ14'
          IC0=IC0+25
        end if
	END DO
	IF(IC0.GT.500)IC0=500
      if(IPRINT.ge.7)then
        write(*,*)
        write(*,*)' Distribution of some intramolecular energies'
	  write(*,*)STR(1:IC0)
	  do J=1,NE
	    write(*,3335)J,EC(J),(AW(I,J),I=IAV0,IAV)
	  end do
	end if
 3335 format(I5,f7.4,50(1x,3f8.3))
*
	write(str(21:500),'(100x,100x,100x,100x,80x)')  
	IC0=17
      IAV0	= IAV+1 
	if(IPRINT.ge.5)then
	  write(*,*)
        write(*,*)'-----------------------------------------------'
      end if
	DO ITYP        =    1,NTYPES
        DO JTYP        = ITYP,NTYPES
          I              =  MDX(ITYP,JTYP)
          IAV            = IAV+2
          if(IC0.le.485)STR(IC0:IC0+12)=NAME(ITYP)//'-'//NAME(JTYP)
          IC0=IC0+17
        END DO! OF JTYP
      END DO! OF ITYP
	IF(IC0.GT.500)IC0=500
      if(IPRINT.ge.5)then
        write(*,*)
        write(*,*)' Distribution of intermolecular energies'
	  write(*,*)STR(1:IC0)
	  NP=NTYPES*(NTYPES+1)/2
	  do IP=1,NP
	    TINN(IP)=0.
	  end do
	  do J=1,NE
	    write(*,3336)J,EC(J),(AW(I,J),I=IAV0,IAV)
	    if(J.gt.1.and.EC(J).ne.0.d0)then
	      do IP=1,NP
	        TINN(IP)=TINN(IP)+
     +        (sqrt(EC(J-1))/sqrt(EC(J))-1.d0)*AW(IAV0-1+2*IP,J)+
     +        (sqrt(EC(J-1))/sqrt(EC(J))-1.d0)*AW(IAV0-2+2*IP,J)-
     -        (sqrt(EC(J))/sqrt(EC(J-1))-1.d0)*AW(IAV0-1+2*IP,J-1)-
     -        (sqrt(EC(J))/sqrt(EC(J-1))-1.d0)*AW(IAV0-2+2*IP,J-1)
	      end do
	    end if
	    if(EC(J).eq.0.d0)then
	      do IP=1,NP
	        TINN(IP)=TINN(IP)-
     -        (sqrt(EC(J))/sqrt(EC(J-1))-1.d0)*AW(IAV0-1+2*IP,J-1)+
     -        (sqrt(EC(J))/sqrt(EC(J-1))-1.d0)*AW(IAV0-2+2*IP,J-1)
	      end do
	    end if
	  end do           
        write(*,*)'-----------------------------------------------'
	  write(*,'(a10,2x,50(9x,f8.2))')
     +  ' integral ',(TINN(IP),IP=1,NP)
	end if                
 3336 format(I5,f7.4,50(1x,2f8.2))
*
      write(str(19:500),'(100x,100x,100x,100x,82x)')   
	IAV0	= IAV+1
      write(*,*)'-----------------------------------------------'
	write(*,*)
      IAV            = IAV+6
	STR(13:19)=' TEMP  '
	STR(20:28)=' EPinter '
	STR(29:37)=' EPintra '
	STR(38:46)=' EPtotal '
	STR(47:55)=' PRES    '
	STR(56:64)=' Esolv   '
	STR(65:73)=' TS      '
      write(*,*)
      write(*,*)' Distribution of global thermodynamic properties'
	write(*,*)STR(1:80)
	ET0=AW(IAV0+1,NE)+AW(IAV0+2,NE)
        if(LSHEJK)then
           IAVP=IAV0+4
        else
           IAVP=IAV0+3
        end if
	do J=1,NE                                                   
	  ETOT=AW(IAV0+1,J)+AW(IAV0+2,J)
	  ESOLV=NOP*(ETOT-ET0)/NDEL
          FF=FEN(J)*0.001*ENERF/BETA                         ! F (kJ/M)
          SS=ESOLV-FF
	  write(*,'(I3,f8.4,f8.2,3f9.3,5f9.2)')
     +  J,EC(J),AW(IAV0,J),AW(IAV0+1,J),AW(IAV0+2,J),
     +  ETOT,AW(IAVP,J),ESOLV,SS  
	end do
	PP1=AW(IAVP,1)
	PPM=AW(IAVP,NE)
	if(LNPT)then
          VOL1            = AW(IAV,1)
          VOLM            = AW(IAV,NE)
	  PVTERM1	= -TRPRES*VOL1*AVSNO*1.d-28/NOP
	  PVTERM2	= -TRPRES*VOLM*AVSNO*1.d-28/(NOP-NDEL)
	else
	  PVTERM1	= -PP1*VOL*AVSNO*1.d-28/NOP
	  PVTERM2	= -PPM*VOL*AVSNO*1.d-28/(NOP-NDEL)
	end if 
	PVTERM=0.5*(PVTERM1+PVTERM2)
*	write(*,*)' Pressure= ',PRESA, 'atm'
*	if(.NOT.LNPT)then
          PKT=0.001*ENERF/BETA
          PKT2=0.001*ENERF*NOP/(BETA*(NOP-NDEL))
*	  write(*,*)' PV-term 1 = ',PVTERM1,' kJ/M' 
*	  write(*,*)' PV-term 2 = ',PVTERM2,' kJ/M' 
	  write(*,*)' kT-term   = ',PKT,' kJ/M'
*	  write(*,*)' kT-term 2 = ',PKT2,' kJ/M'
	  write(*,671)' Free energy:',FRM,'+/-',FRD,' kJ/M' 
*	  write(*,671)' Alt. Fr.en.:',FRM+PVTERM2+PKT2,'+/-',FRD,' kJ/M' 
*	else
*	  PKT=-0.001*ENERF*NDEL*dlog(TPRES*VOL*BETA/(UNITP*NDEL))/BETA
*	  write(*,*)' kT-term = ',PKT,' kJ/M'
*	  write(*,671)' Free energy:',FRM+PKT,'+/-',FRD,' kJ/M'
*	end if 
 671	format(a,f12.4,a,f9.4,a)
	write(*,*)'----------------------------------------------------'
      write(*,'(17x,3f9.3,4f9.2)')AV1,AV2,AV3,AV4,AV5,AV6
	if(IPRINT.lt.4)return
*
      write(*,*)
      write(*,*)' Fractional temperatures'
      STR(14:80)=
     +'   Dens     boxX    boxY   boxZ    Ttrans    Trot     Tint'    
	write(*,*)STR(1:75)
	do J=1,NE
          AVX = AW(IAV,J)
          AVY = AW(IAV+1,J)
          AVZ = AW(IAV+2,J)
        AV7            = AVX*AVY*AVZ
        DENS	= TOTMAS*1.d-3/(AV7*UNITL**3)
	  write(*,'(I5,f8.4,f9.4,3f9.2,3f9.3)')
     +  J,EC(J),DENS,
     +  AVX,AVY,AVZ,(AW(I,J),I=IAV+3,IAV+5)
	end do
      write(*,*)'----------------------------------------------------'
*     
      write(*,*)' Table of transitions'
      write(*,'(4x,50i4)')(I,I=1,NE)
      do J=1,NE
        write(*,'(51i4)')J,(IWALK(J,I),I=1,NE)
      end do          
      RETURN
      END 
*
*================== ISHIFT =====================================
*
	subroutine ISHIFT(IST,LHIST,NNAV,IN)
	integer*4 IST(IN,LHIST)
	do I=2,LHIST
	  do J=1,NNAV 
	    IST(J,I-1)=IST(J,I)
	  end do
	end do
	return
	end 
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
*=============== ZEROAV =============================================
*
      SUBROUTINE ZEROAV
      include "mdee.h"
*
*Zero averages
*
      DO I       = 1,NB
        BB(I)      = 0.D0
        EB(I)      = 0.D0
      END DO
*
      DO I       = 1,NA
        AA(I)      = 0.D0
        EA(I)      = 0.D0
      END DO
*
      DO I       = 1,NT
        TT(I)      = 0.D0
        ET(I)      = 0.D0
      END DO 
*
      RETURN
      END
