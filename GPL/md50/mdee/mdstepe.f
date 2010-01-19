*=============== SULTAN ================================================
*
      SUBROUTINE SULTAN
*
*  Constrained MD
*
      include "mdee.h"
*
      LOGICAL LRDF
      DATA LRDF/.false./
*
*>----------------------------
      CALL GETCOM     ! supl.f      
      call CHCNB(LRDF)       ! mdstep.f
      LRDF=.true.
      NSTEP=0
 1    NSTEP=NSTEP+1
      MSTEP=NSTTOT+NSTEP
*>----------------------------
*
        CALL ZEROFS(1)
	call ZEROAV
	LCSRDF=.true.
*       
      if(mod(MSTEP,ICHNB).eq.0)call CHCNB(LRDF)
*
      DO ITYP = 1,NTYPES
      CALL GETBND(ITYP)   ! i-forces.f
      CALL GETANG(ITYP)   ! i-forces.f
      CALL GETTRS(ITYP)   ! i-forces.f
      CALL GETIMP(ITYP)   ! i-forces.f
      END DO!  OF ITYP 
      call LOCALF         ! l-forces.f
*

      CALL ZEROFS(2)  ! supl.f
      CALL FURIR      ! l-forces.f
      CALL ETERMS     ! l-forces.f
      CALL FORCES     ! l-forces.f
      if(LNPT)call ELRCLJ     ! l-forces.f
*
      timeg0	= cputime(0.d0)
      GSUM           = 0.D0
      HSUM           = 0.D0
*  Forces summation and pressure evaluation
      DO I           = 1,NSTOT
        FX(I)          = HX(I)+GX(I)+OX(I)
        FY(I)          = HY(I)+GY(I)+OY(I)
        FZ(I)          = HZ(I)+GZ(I)+OZ(I)
        HSUM           = HSUM+HX(I)+HY(I)+HZ(I)
        GSUM           = GSUM+GX(I)+GY(I)+GZ(I)
      END DO! OF I
*
      CALL SCLFRC(DONE,3)
*
      DO I         = 1,NSTOT
Constrained velocities at n-1/2
      HX(I)        = VX(I)
      HY(I)        = VY(I)
      HZ(I)        = VZ(I)
Constrained positions at n
      OX(I)        = SX(I)
      OY(I)        = SY(I)
      OZ(I)        = SZ(I)
      END DO! OF I
*
      DO ITYP      = 1,NTYPES
      NSP          = NSPEC(ITYP)
      NSS          = NSITS(ITYP)
      NSNSP        = NSS*NSP
*
      I            = ISADDR(ITYP)
      DO L         = 1,NSNSP
      I            = I+1
*Unconstrained velocities at step n+1/2
      VX(I)        = VX(I)+TSTEP*FX(I)
      VY(I)        = VY(I)+TSTEP*FY(I)
      VZ(I)        = VZ(I)+TSTEP*FZ(I)
*Unconstrained  postions  at step n+1
      SX(I)        = SX(I)+TSTEP*VX(I)*MASSDI(I)
      SY(I)        = SY(I)+TSTEP*VY(I)*MASSDI(I)
      SZ(I)        = SZ(I)+TSTEP*VZ(I)*MASSDI(I)
*fx is used as a temporary array for unconstrained positions (now OX...)
      FX(I)       = SX(I)
      FY(I)       = SY(I)
      FZ(I)       = SZ(I)
      END DO! OF L
*
      END DO! OF ITYP...
*  Constrain dynamics                           
	do ITYP=1,NTYPES
	  CALL SHEJK(ITYP)
	end do
*  Realization of specified ensemble
      CALL SCALING(1)
*  Change of ensemble
      if(mod(MSTEP,NCEE).eq.0)call CHENS
*  Collecting averages
      CALL GETKIN     ! supl.f
      CALL GETCOM     ! supl.f
*      CALL GETROT     ! supl.f
*      CALL DIPOLE
      if(mod(MSTEP,IAVER).eq.0)CALL GETAVR(1)
*
      NSTEPA=NSTEPA+1
      VIRL=VIRL+VIRLS
      PEST=PEST+PELS
      WIRS=WIRS+WIRSS+WIRF
      TRYCKE	= PEST*UNITP/(3.*VOL)    ! electrostatic contr.
      TRYCKL	= VIRL*UNITP/(3.*VOL)    ! LJ
      TRYCKB	= VIRB*UNITP/(3.*VOL)    ! bond (flexible)
      TRYCKT	= VIRT*UNITP/(3.*VOL)    ! imp.
      TRYCKN	= VIRN*UNITP/(3.*VOL)    ! interaction inside a molecule
      TRYCKF	= WIRF*UNITP/(3.*VOL)    ! self-int. correction (rigid)
      TRYCKA	= VIRA*UNITP/(3.*VOL)    ! shejk-correction 
      TRYCKS	= WIRS*UNITP/(3.*VOL)    ! orientational correction
      TRYCKD	= VIRD*UNITP/(3.*VOL)    ! long-range correction
      PENL=PE*0.001*ENERF/FNOP
      PEL=PEST*0.001*ENERF/FNOP
      PENS=PINT*0.001*ENERF/FNOP
      PEKI=TKE*0.001*ENERF
      ETOT=PENL+PENS+PEKI
      POTEN	= PENL+PENS
      DENS	= TOTMAS*1.d-3/(VOL*UNITL**3)
*
      if(IPRINT.ge.5.and.mod(MSTEP,IAVER).eq.0)
     +WRITE(6,991) MSTEP,ME,POTEN,PENS,PEL,ETOT,TEMP,TRYCKM,DENS
*  990 FORMAT('%',A1,I5,2(1X,F10.3),1X,F9.2,1X,F9.3,1X,F10.3,1X,F8.4)
      if(IPRINT.ge.7)write(*,'(2f10.6,f10.2,3f10.2)')
     +   SC,SCL,TRYCK,TTR,TROT,TINT
 	if(IPRINT.ge.7.and.mod(MSTEP,IAVER).eq.0)
     +WRITE(6,'(8f10.3)')TRYCKL,TRYCKE,TRYCKB,TRYCKA,TRYCKN,
     +TRYCKF,TRYCKS,TRYCKD
  991 FORMAT(I7,I4,4F9.3,2(1X,f9.2),f9.4)
**
      IF(mod(MSTEP,NAVER).eq.0)call GETAVR(2)
      if(LDMP.and.mod(NSTTOT+NSTEP,NDUMP).eq.0)call RESTRT(2)
      timeg	= timeg+cputime(timeg0)
*
*<---------------------
      if(NSTEP.lt.NSTEPS)go to 1
*<---------------------
*
      RETURN
      END
*
*=============== DOUBLE ================================================
*
      SUBROUTINE DOUBLE
*
      include "mdee.h"
*
      LOGICAL LRDF
      DATA LRDF/.false./
*
*>----------------------------
*    Restore forces, COM and list of neigbours
*
      PERMOL = 0.001*ENERF/NOP
      CALL ZEROFS(1)
      CALL GETCOM     ! supl.f 
      call CHCNB(LRDF)      ! mdstep.f
      LRDF=.true.
      NFR2	= NFREQ/2+1
**       
      LCSRDF=.false.
      DO ITYP = 1,NTYPES
        CALL GETBND(ITYP)   ! i-forces.f
        CALL GETANG(ITYP)   ! i-forces.f
        CALL GETTRS(ITYP)   ! i-forces.f
      END DO!  OF ITYP
      call LOCALF         ! l-forces.f
      CALL SCLFRC(DONE,1)
*
*
      CALL ZEROFS(2)  ! supl.f 
      CALL FURIR      ! l-forces.f
      CALL ETERMS     ! l-forces.f
      CALL FORCES     ! l-forces.f
      call ELRCLJ     ! l-forces.f
      CALL SCLFRC(DONE,2)
*
      if(IPRINT.ge.7)then
         write(*,'(3(a,f12.4))')
     +' Eel= ',(PEST+PELS)*PERMOL,' Eint=',PINT*PERMOL,' Epot',PE*PERMOL
      end if
*
*   Begin MD
*
      NSTEP=0
 1    NSTEP=NSTEP+1
	MSTEP=NSTTOT+NSTEP
*>----------------------------
                         
*  velocity correction due to slow forces
      DO I           =  1,NSTOT
        VX(I)          =   VX(I)+TSTEP*GX(I)*0.5D0
        VY(I)          =   VY(I)+TSTEP*GY(I)*0.5D0
        VZ(I)          =   VZ(I)+TSTEP*GZ(I)*0.5D0
      END DO! OF I
      if(LNPT)then
	  TRL=(VIR+PEST+VIRD)*UNITP/(3.*VOL)
	  SCL=SCL+0.5d0*(TRL-TRPRES)*VOL*DQP*TSTEP
	end if 
*
*- - - - - - - FAST MOTION - - - - - - - -
*     
      LCSRDF=.false.
      call ZEROAV
      DO NSTEB       =  1,NFREQ 
	if(NSTEB.eq.NFREQ)LCSRDF=.true.
        DO I           =  1,NSTOT
	  FX(I) = VX(I)
          FY(I) = VY(I)
          FZ(I) = VZ(I)
          OX(I) = SX(I)
          OY(I) = SY(I)
          OZ(I) = SZ(I)
          VX(I) =   VX(I)+TSTEB*HX(I)*0.5D0         ! V(t+1/2)
          VY(I) =   VY(I)+TSTEB*HY(I)*0.5D0
          VZ(I) =   VZ(I)+TSTEB*HZ(I)*0.5D0
          SX(I) =   VX(I)*TSTEB*MASSDI(I)+SX(I)     ! X(t+1)
          SY(I) =   VY(I)*TSTEB*MASSDI(I)+SY(I)
          SZ(I) =   VZ(I)*TSTEB*MASSDI(I)+SZ(I)
        END DO! OF I
*  Realization of the specified ensemble
        CALL SCALING(2)       ! mdstep.f  
        CALL ZEROFS(1)     ! supl.f
	CALL GETCOM 
	if(NSTEB.eq.NFR2.and.mod(MSTEP,ICHNB).eq.0)then
           call CHCNB(LRDF)
        end if
*  Recalculation of slow forces
	if(mod(NSTEB,NFREQ).eq.0)then
          CALL ZEROFS(2)      ! supl.f
          CALL FURIR          ! l-forces.f
          CALL ETERMS         ! l-forces.f
          CALL FORCES         ! l-forces.f
*          if(LNPT)call ELRCLJ         ! l-forces.f
          CALL SCLFRC(DONE,2)   ! supl.f
        end if
	DO ITYP = 1,NTYPES
          CALL GETBND(ITYP)     ! i-forces.f           F(t+1)
          CALL GETANG(ITYP)     ! i-forces.f
          CALL GETTRS(ITYP)     ! i-forces.f
        END DO!  OF ITYP
	call LOCALF            ! l-forces.f
*
        CALL SCLFRC(DONE,1)    ! supl.f
*         
	  if(NSTEB.ne.NFREQ)then
          DO I           =  1,NSTOT
            VX(I)          =  VX(I)+TSTEB*HX(I)*0.5D0     !  V(t+1)
            VY(I)          =  VY(I)+TSTEB*HY(I)*0.5D0
            VZ(I)          =  VZ(I)+TSTEB*HZ(I)*0.5D0
          END DO! OF I        
	  else
          DO I           =  1,NSTOT
            VX(I)          =  VX(I)+TSTEB*(HX(I)-OX(I))*0.5D0     !  V(t+1)
            VY(I)          =  VY(I)+TSTEB*(HY(I)-OY(I))*0.5D0
            VZ(I)          =  VZ(I)+TSTEB*(HZ(I)-OZ(I))*0.5D0
          END DO! OF I        
	  end if
*
*  Realization of the specified ensemble
        CALL SCALING(3)       ! mdstep.f
*  Change of subensemble
	  if(NSTEB.eq.NFR2.and.mod(MSTEP,NCEE).eq.0)call CHENS
*
        if(IPRINT.ge.8)then
          CALL GETEMP     ! aver.f
  990 FORMAT('%',A1,I5,2(1X,F10.3),1X,F9.2,1X,F10.3,1X,F9.2,1X,F9.2)
          WRITE(6,990) WHICH,NSTEB,POTEN,ETOT,TEMP,VIR,TRL,TRYCK
        end if
      END DO! OF NSTEB
*
*- - - - - - - SLOW MOTION - - - - - - - -
*     
      timeg0	= cputime(0.d0)
*
      DO I           = 1,NSTOT
        VX(I)          =  VX(I)+TSTEP*(GX(I)+OX(I))*0.5D0
        VY(I)          =  VY(I)+TSTEP*(GY(I)+OY(I))*0.5D0
        VZ(I)          =  VZ(I)+TSTEP*(GZ(I)+OZ(I))*0.5D0
      END DO! OF I
      if(LNPT)then
	TRL=(VIR+PEST+VIRD)*UNITP/(3.*VOL)
	SCL=SCL+0.5d0*(TRL-TRPRES)*VOL*DQP*TSTEP
      end if 
*
*  Collecting averages
      NSTEPA=NSTEPA+1
      VIRLL=VIRL+VIRLS
      PELL=PEST+PELS
      WIRSUM=WIRS+WIRSS+WIRF
      TRYCKE	= PELL*UNITP/(3.*VOL)
      TRYCKL	= VIRLL*UNITP/(3.*VOL )
      TRYCKB	= VIRB*UNITP/(3.*VOL)
      TRYCKT	= VIRT*UNITP/(3.*VOL)
      TRYCKN	= VIRN*UNITP/(3.*VOL)
      TRYCKF	= WIRF*UNITP/(3.*VOL)
      TRYCKA	= VIRA*UNITP/(3.*VOL)
      TRYCKS	= WIRSUM*UNITP/(3.*VOL)
      TRYCKD	= VIRD*UNITP/(3.*VOL)
*      T1	= VIR*UNITP/(3.*VOL)
*	T2	= VIRS*UNITP/(3.*VOL)
*	write(*,'(a4,3f12.3)')' P= ',T1,T2,T1+T2
      CALL GETROT       ! supl.f
*      CALL GETKIN       ! supl.f
*      CALL DIPOLE       ! supl.f
      if(mod(MSTEP,IAVER).eq.0)CALL GETAVR(1)
      PENL=PE*0.001*ENERF/FNOP
      PEL=PELL*0.001*ENERF/FNOP
      PENS=PINT*0.001*ENERF/FNOP
      PEKI=TKE*0.001*ENERF
      ETOT=PENL+PENS+PEKI
      POTEN	= PENL+PENS
      DENS	= TOTMAS*1.d-3/(VOL*UNITL**3)
*
*      if(IPRINT.ge.6)WRITE(6,'(4(a6,f12.5))')
*     +' Eex= ',PENL,' Ein= ',PENS,' Eki= ',PEKI,' Etot=',ETOT 
      if(IPRINT.ge.5.and.mod(MSTEP,IAVER).eq.0)
     +WRITE(6,991) MSTEP,ME,POTEN,PENS,PEL,ETOT,TEMP,TRYCK,DENS
      if(IPRINT.ge.7)write(*,'(2f10.6,4f10.2)')
     +     SC,SCL,TRYCKM,TTR,TROT,TINT
 	if(IPRINT.ge.7.and.mod(MSTEP,IAVER).eq.0)
     +WRITE(6,'(8f10.3)')TRYCKL,TRYCKE,TRYCKB,TRYCKA,TRYCKN,
     +TRYCKF,TRYCKS,TRYCKD
  991 FORMAT(I7,I4,4F9.3,2(1X,f9.2),f9.4)
*
      if(mod(MSTEP,NAVER).eq.0)CALL GETAVR(2)
      if(LDMP.and.mod(MSTEP,NDUMP).eq.0)call RESTRT(2)
      timeg	= timeg+cputime(timeg0)
*
*<---------------------
      if(NSTEP.lt.NSTEPS)go to 1
*<---------------------
*
      RETURN
      END
*
*=============== SCALING ==============================================
*
      SUBROUTINE SCALING(IAL)
*
*  NVT, NPT or NVE with temp-scaling ensembles
*    IAL=1      leap frog algorithm 
*        2      double time scale, before moves 
*        3                         after
      include "mdee.h"
*
      DIMENSION TTT(NTPS)
	data DKEO/0./
*
      call GETEMP
* Pressure calculation
      VIRSUM	= VIR+VIRS+PEST+PELS+VIRD
      WIRSUM	= VIRL+VIRLS+PEST+PELS+VIRD+WIRS+WIRSS+WIRF
      TRID	= 2.*TKE*UNITP/(3.*VOL)
      TRIDM	= (NOP-1)*BOLTZ*TEMP*1.d-5/(VOL*UNITL**3)
      TRYCK	= VIRSUM*UNITP/(3.*VOL)+TRID
      TRYCKM	= WIRSUM*UNITP/(3.*VOL)+TRIDM 
      TRL=(VIR+PEST+VIRD)*UNITP/(3.*VOL)       ! slow virial 
      TRS=(VIRS+PELS)*UNITP/(3.*VOL)+TRID      ! fast virial
*
*    NPT molecular dynamics
*    
*   new thermostat parameters
      if(LNVT)then
	  DKE         = TEMP/TRTEMP-1.d0
*
*    1 step leap-frog algorithm
*
	if(IAL.eq.1)then     
	  if(LNPT)then
	    DPE = TRYCKM-TRPRES
	    SCLO= SCL
	    SCL = (SCL+DPE*VOL*DQP*TSTEP)/(1.d0-SCL*TSTEP) ! ita (t+1/2)
	    SCLM= 0.5d0*(SCLO+SCL)                         ! ita (t)
	    N   = 0
            N0  = 0
            I         = 0
	    SCC	= 0.5d0*SCL*TSTEP
	    DCC	= SCC/(1.d0-SCC)
	    SCCV	= 0.5d0*SCLM*TSTEP
	    DCV	= SCCV/(1.d0-SCCV)
            DO ITYP   = 1,NTYPES              ! types
              NSBEG       = ISADR  (ITYP)+1
              NSEND       = ISADR  (ITYP +1)
              SUMM        = SUMMAS (ITYP)
              DO J        = 1,NSPEC(ITYP)     ! molecules
                I           = I+1
		XC       = 0.D0
      	        YC       = 0.D0
      	        ZC       = 0.D0
	        PXC	= 0.
                PYC	= 0.
                PZC	= 0.
Calculate C.O.M. vectors:
      	        DO IS       = NSBEG,NSEND
      	          N           = N+1
      	          XC       = XC+MASS(IS)*SX(N)
      	          YC       = YC+MASS(IS)*SY(N)
      	          ZC       = ZC+MASS(IS)*SZ(N)
                  PXC	= PXC+VX(N)
                  PYC	= PYC+VY(N)
                  PZC	= PZC+VZ(N)
                END DO! OF IS
                XC       = XC/SUMM    ! COM after intergation velocities
      	        YC       = YC/SUMM
      	        ZC       = ZC/SUMM
*    NPT-scaling  - shift of COM and transl. momenta
		DX       = (XC+X(I))*DCC
                DY       = (YC+Y(I))*DCC
                DZ       = (ZC+Z(I))*DCC
                DVX	= (PXC+PX(I))*DCP
                DVY	= (PYC+PY(I))*DCP
                DVZ	= (PZC+PZ(I))*DCP
                DO IS       = NSBEG,NSEND
                  N0           = N0+1
                  SX (N0)       = SX(N0)+DX
                  SY (N0)       = SY(N0)+DY
                  SZ (N0)       = SZ(N0)+DZ
                  FAC	= MASS(IS)/SUMM
                  VX(N0)	= VX(N0)-DVX*FAC
                  VY(N0)	= VY(N0)-DVY*FAC
                  VZ(N0)	= VZ(N0)-DVZ*FAC
                end do
	        if(N.ne.N0)write(*,*)
     +        '!!! count error in SCALING',ITYP,J,I,N,N0
	      end do   ! of J (I)
            end do ! of ITYP 
            CL = BOXL*(1.d0+0.5d0*SCL*TSTEP)/(1.d0-0.5d0*SCL*TSTEP)
            call RECLEN(CL)
	    end if   ! if(LNPT...
*   NVT-constraints dynamics
	    SCC = 0.5d0*SC*TSTEP
	    do I=1,NSTOT
	      VX(I)=(VX(I)-HX(I)*SCC)/(1.d0+SCC)
	      VY(I)=(VY(I)-HY(I)*SCC)/(1.d0+SCC)
	      VZ(I)=(VZ(I)-HZ(I)*SCC)/(1.d0+SCC)
	    end do
          if(LNPT)SC = SC+(0.5d0*SCL**2-DQT/FNST)*TSTEP
	    SC=SC+DKE*DQT*TSTEP                   !  zeta(t+1)
*
*    multiple time step; t -> t+1/2
*
	else if(IAL.eq.2)then  
	  SCC = 0.5d0*SC*TSTEB
          if(LNPT)then
	    SCC = SCC+0.5d0*SCL*TSTEB          
	    SCL = (SCL+0.5d0*TRS*VOL*DQP*TSTEB)/(1.d0-0.5d0*SCL*TSTEB) !ita(t+1/2)
	    DCL = 0.5d0*SCL*TSTEB
	    do I=1,NSTOT
	      SX(I)=(SX(I)+OX(I)*DCL)/(1.d0-DCL)    !  t+1
	      SY(I)=(SY(I)+OY(I)*DCL)/(1.d0-DCL)
	      SZ(I)=(SZ(I)+OZ(I)*DCL)/(1.d0-DCL)
	    end do
	    CL = BOXL*(1.d0+0.5d0*SCL*TSTEB)/(1.d0-0.5d0*SCL*TSTEB)
            call RECLEN(CL)
          end if
	  do I=1,NSTOT
	    VX(I) = VX(I) - SCC*FX(I)           ! t+1/2
	    VY(I) = VY(I) - SCC*FY(I)
	    VZ(I) = VZ(I) - SCC*FZ(I)
	  end do
	  DKEO=DKE       ! (t+1/2)
*
*   multiple time step; t+1/2 -> t+1
*
        else if(IAL.eq.3)then 
	    if(LNPT)then
	      SC  = SC+(0.5d0*SCL**2-DQT/FNST)*TSTEB     
	      SCL = (SCL+0.5d0*TRS*VOL*DQP*TSTEB)/(1.d0-0.5d0*SCL*TSTEB) !ita(t+1)
	      SCC = 0.5d0*SCL*TSTEB
	    else
	      SCC = 0.d0
	    end if
	    SC  = SC+DKEO*DQT*TSTEB                   !  zeta(t+1)
	    SCC = SCC+0.5d0*SC*TSTEB
	    DCC = 1.d0/(1.d0+SCC)  
	    do I=1,NSTOT
	      VX(I) = VX(I)*DCC
	      VY(I) = VY(I)*DCC
	      VZ(I) = VZ(I)*DCC
	    end do
	  end if
      END IF!(LNVT...
*
* scale velocities if necessary
*
      IF(LSCLT) THEN
        DO ITYP     = 1,NTYPES
          I           = ISADDR(ITYP)
          IF(DABS(TRTEMP-TTT(ITYP)).GT.TDELT) THEN
            SCV         = DSQRT(TRTEMP/TTT(ITYP))
            NSPNS       = NSPEC(ITYP)*NSITS(ITYP)
            DO N        = 1,NSPNS
              I           = I+1
              VX(I)       = SCV*VX(I)
              VY(I)       = SCV*VY(I)
              VZ(I)       = SCV*VZ(I)
            END DO! OF I
            if(IPRINT.ge.5)write(*,*)
     +      'velocities scaled by ',SCV,' for type',ITYP
            NRTSC(ITYP) =  NRTSC(ITYP)+1
          END IF!(DABS..
        END DO! OF ITYP
        call GETEMP
      END IF!(LSC...
      TKE	= TKE/FNOP
*
      RETURN
      END
