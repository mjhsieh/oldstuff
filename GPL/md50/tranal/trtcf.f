C========================================
C
C     MDynaMix package
C     Trajectory analysis block
*
C     Calculation of time correlation functions 
C     
*
*============ TRTCF ===================================================
*
      program TRTCF
      include "tranal.h"
      include "trtcf.h"
*
      DIMENSION M1(NTPS),M2(NTPS),JTCF(NTCF)
      character*32 FILTCF
      namelist /TCF/ NSTEG,DTCF,ITCF,N1,N2,FILTCF
*
*  1.1 TCF types
*  
      NAMECF(1)  = 'VACF  '
      NAMECF(2)  = 'WACF  '
      NAMECF(3)  = 'P1(dp)'
      NAMECF(4)  = 'P2(dp)'
      NAMECF(5)  = 'unitv1' 
      NAMECF(6)  = 'unitv2'
      NAMECF(7)  = 'VACF-X'
      NAMECF(8)  = 'VACF-Y'
      NAMECF(9)  = 'VACF-Z'
      NAMECF(10) = 'WACF-X'
      NAMECF(11) = 'WACF-Y'
      NAMECF(12) = 'WACF-Z'
      NSTCF=0
      ICTCF=0
*  1.2 Setup trajectory
      call SETTRAJ
      IEND=0
* read input
      read(*,TCF)
      DTCF=DTCF*0.001    !   in ps
      do I=1,NTCF
	if(ITCF(I).ne.0)NSTCF=NSTCF+1
      end do
      if(ITCF(1).ne.0.or.ITCF(7).ne.0.or.ITCF(8).ne.0.or.ITCF(9).ne.0)
     +ITCF(NTCF+1)=1                                           ! velocity
      if(ITCF(2).ne.0.or.ITCF(10).ne.0.or.ITCF(11).ne.0.or.
     +ITCF(12).ne.0)ITCF(NTCF+2)=1                             ! ang. vel.
      if(ITCF(3).ne.0.or.ITCF(4).ne.0)ITCF(NTCF+3)=1           ! dip.moment
      if(ITCF(5).ne.0.or.ITCF(6).ne.0)then
        ITCF(NTCF+4)=1           ! unit vec.
*                                                        
        do ITYP=1,NTYPES
	  NSS = NSITS(ITYP)
          IF(N1(ITYP).LT.1  .OR.N2(ITYP).LT.1.or.
     +       N1(ITYP).GT.NSS.OR.N2(ITYP).GT.NSS)then
	  write(*,*)' Fault unit vector for tcf calculation. Type ',
     +ITYP,' Sites ',N1(ITYP),N2(ITYP),' of ',NSS
	    stop      
	  end if
	end do
      end if
*  Case of TCF in local coord.system
      if(ITCF(7).eq.2.or.ITCF(8).eq.2.or.ITCF(9).eq.2)then
        ITCF(NTCF+5)=1
        LMOL1=.true.
      else
        LMOL1=.false.
      end if
      if(ITCF(10).eq.2.or.ITCF(11).eq.2.or.ITCF(12).eq.2)then
        ITCF(NTCF+6)=1
        LMOL2=.true.
      else
        LMOL2=.false.
      end if
*
      NOM        = NOP*3
      NNTCF	= (MAXCF+2)*NOM
      NCALLS     = 0           
      MEM=NSTCF*MCF*8 
      do I=NTCF+1,NTCF+6
	if(ITCF(I).ne.0)MEM=MEM+4*NNTCF
      end do
      DMEM=0.001*MEM
*     
        PRINT "(/'*** TCF CALCULATION WILL BE PERFORMED ***')"
        PRINT "('   .... FOR ',I3,' TYPES')",NSTCF
        PRINT "(/'*** CC ARRAYS WILL TAKE ',f9.3
     X       ,'  KB')",DMEM 
*
        write(*,*)
        write(*,*)
     +  '*** FOLLOWING CORRELATION FUNCTIONS WILL BE CALCULATED:'
        write(*,*)
        DO I       = 1,NTCF
          if(ITCF(I).ne.0)
     +    PRINT "(10X,'*** TYPE:',I3,'     --> ',A6)",I,NAMECF(I)
        END DO! OF I
        PRINT *
*
      IF(NSTEG.GT.MAXCF)            STOP '!!! NSTEG.GT.MAXCF'
      IF((NOM*(NSTEG+2)).GT.MAXTCF) STOP '!!! INCREASE MAXTCF !'
*
*  1.3 Zero arrays
*
      DO I       = 1,NSTEG
        INX(I)     = 0
      END DO! OF I
*
      DO  I      = 1,MAXCF*NTPS
      CFV(I)     = 0.D0
      CFA(I)     = 0.D0
      CP1(I)     = 0.D0
      CP2(I)     = 0.D0
      CU1(I)     = 0.D0
      CU2(I)     = 0.D0
      CVX(I)     = 0.D0
      CVY(I)     = 0.D0
      CVZ(I)     = 0.D0
      CAX(I)     = 0.D0
      CAY(I)     = 0.D0
      CAZ(I)     = 0.D0
      END DO! OF I
*
      DO  I      = 1,NNTCF
        CC1(I)     = 0.D0
        CC2(I)     = 0.D0
        CC3(I)     = 0.D0
        CC4(I)     = 0.D0
        CC5(I)     = 0.D0
        CC6(I)     = 0.D0
      END DO! OF I
      FULTIM0 = TIMIN(1)
      TRTEMP=0.
*
*  1.4 Start trajectory analysis
*
      do while(IEND.eq.0)
        call READCONF(IEND)
        if(IEND.eq.0)then
	  if((FULTIM-FULTIM0)+1.d-6.ge.DTCF*ICTCF)then
 237	    ICTCF=ICTCF+1
	    call GETCOM
	    call DIPOLE
	    call GETROT
	    call GETTCF
	    TRTEMP=TRTEMP+TEMP
	    if(IPRINT.ge.6)write(*,'(a,I6,f11.5,3x,f10.2)')
     +  ' tcf: ',ICTCF,FULTIM,TEMP
	    if(FULTIM-FULTIM0.gt.DTCF*ICTCF)then
	      write(*,*)' no config. at t= ',DTCF*ICTCF+FULTIM0
	      go to 237
	    end if
	    if(IPRINT.ge.8)then
	      write(*,*)' local atom coordinates in principal s'
	      do I=1,NSTOT
		 IMOL=NNUM(I)
		 QX(IMOL)=WX(I)
		 QY(IMOL)=WY(I)
		 QZ(IMOL)=WZ(I)
		 call ROTMAT(QX,QY,QZ,IMOL,1,.TRUE.)
             write(*,'(2I5,3f10.4)')I,IMOL,QX(IMOL),QY(IMOL),QZ(IMOL)
	      end do
	    end if
	  end if  ! if/FULTIM
        end if
      end do
      write(*,*)IAN,' configurations analysed'
      TRTEMP=TRTEMP/IAN
*   2. TCF output
*   --------------
*    2.1. Linear momenta:
*    --------------------
      open(unit=25,file=filtcf,status='unknown')
      DEL      = DTCF
      if(ITCF(1).ne.0)then
        write(25,*)'----------------------------------------------'
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(1),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
*
        write(25,*)' Temperature is ',TRTEMP
        write(25,*)
        DO ITYP  = 1,NTYPES
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      NSP      = NSPEC(ITYP)
      FNSPI3   = 1.D0/DFLOAT(NSP)/3.D0
      ISB      = ISADR(ITYP)+1
      ISE      = ISADR(ITYP +1)
      FMASS    = 0.D0
      DO IS    = ISB,ISE
      FMASS    = FMASS+MASS(IS)
      END DO! OF IS
      JCF      = (ITYP-1)*NSTEG
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI3
      CNRM     = CFV(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  write(25,*)'!!! EMPTY TCF ARRAY !!!'
	  go to 201
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI3
      T1       = CFV(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
      DV       = BOLTZ*TRTEMP*AVSNO/FMASS*TAU*1.E-5
      write(25,'(a,e12.4,a)')' *** Diffusion coefficient:  ',DV,
     + ' cm**2/s'  
*                          
 201  continue
      END DO! OF ITYP                                   
	end if
*
*    2.2 Angular momenta:
*    --------------------
      if(ITCF(2).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(2),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	write(25,*)' Angular velocity autocorrelation function'
*
        DO ITYP = 1,NTYPES
*     
	  if(NSITS(ITYP).gt.1)then
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI3   = 1.D0/DFLOAT(NSP)/3.D0
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI3
      CNRM     = CFA(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  write(25,*)'!!! EMPTY TCF ARRAY !!!'  
	  go to 202
	  end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI3
      T1       = CFA(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.4,a)')' *** Correlation time:  ',TAU,'  PS'
*     
 202  continue
	end if
      END DO! OF ITYP
	end if
*
*    2.3 Reorientational motion of dipole axis I order Legendre:
*    -----------------------------------------------------------
 	if(ITCF(3).ne.0)then
          write(25,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(3),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
*
      DO ITYP = 1,NTYPES   
	if(NSITS(ITYP).gt.1)then
*
      write(25,*)'*** Reorientation TCF for dipole moment for ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CP1(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 203
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CP1(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*
 203  continue
	end if
      END DO! OF ITYP
	end if
*
*    2.4 Reorientational motion of dipole axis II order Legendre:
*    ------------------------------------------------------------
 	if(ITCF(4).ne.0)then
          write(25,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(4),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	  write(25,*)
     +'Reorientational motion of dipole axis II order Legendre'
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      write(25,*)'***     Time autocorrelation function for  ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI2   = 1.D0/DFLOAT(NSP)/2.D0	! (3xx-1) divided by 2
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI2
      CNRM     = CP2(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 204
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI2
      T1       = CP2(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*                         
 204      continue
	  end if
        END DO! OF ITYP
      end if
*
*    2.5 Reorientational motion of spec. unit vector I order Legendre:
*    -----------------------------------------------------------------
      if(ITCF(5).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(5),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	write(25,*)
     +'Reorientational motion of spec. unit vector I order Legandre'
*
        DO ITYP = 1,NTYPES
	  if(N1(ITYP).ne.N2(ITYP))then
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
          DO I     = 1,NSTEG
            CCF(I)   = 0.D0
          END DO! OF I
          JCF      = (ITYP-1)*NSTEG
          NSP      = NSPEC(ITYP)
          FNSPI    = 1.D0/DFLOAT(NSP)
          FNF      = 1.D0/FLOAT(INX(1))*FNSPI
          CNRM     = CU1(JCF+1)*FNF
          IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
          IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 205
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CU1(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*                    
 205  continue
	end if
      END DO! OF ITYP
	end if
*
*      2.6 Reorientational motion of spec. unit vector II order Legendre:
*      ------------------------------------------------------------------
 	if(ITCF(6).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(6),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	write(25,*)
     +'Reorientational motion of spec. unit vector II order Legandre'
*
      DO ITYP = 1,NTYPES
	if(N1(ITYP).ne.N2(ITYP))then
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI2   = 1.D0/DFLOAT(NSP)/2.D0 ! (3xx-1) divided by 2
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI2
      CNRM     = CU2(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 206
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI2
      T1       = CU2(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*                                 
 206  continue
	end if
      END DO! OF ITYP
      end if
*
*    2.7 Linear momenta X-component:
*    -------------------------------
      if(ITCF(7).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(7),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL1)then
	  write(25,*)' Linear momenta in principal system X-component'
        else
	  write(25,*)' Linear momenta in laboratory system X-component'
        end if
*
      DO ITYP = 1,NTYPES
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CVX(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 207
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CVX(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*
 207	continue
      END DO! OF ITYP
	end if
*
*   2.8 Linear momenta Y-component:
*   -------------------------------
      if(ITCF(8).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(8),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL1)then
	  write(25,*)' Linear momenta in principal system Y-component'
        else
	  write(25,*)' Linear momenta in laboratory system Y-component'
        end if
*
      DO ITYP = 1,NTYPES
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CVY(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 208
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CVY(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*
 208	continue
      END DO! OF ITYP
	end if
*
*    2.9 Linear momenta Z-component:
*    -------------------------------
      if(ITCF(9).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(9),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL1)then
	  write(25,*)' Linear momenta in principal system Z-component'
        else
	  write(25,*)' Linear momenta in laboratory system Z-component'
        end if
*
      DO ITYP = 1,NTYPES
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CVZ(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')" 
	  go to 209
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CVZ(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*
 209  continue 
      END DO! OF ITYP
	end if
*
*    2.10  Angular momenta X-component:
*    ----------------------------------
 	if(ITCF(10).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(10),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL2)then
	  write(25,*)' Angular momenta in principal system X-component'
        else
	  write(25,*)' Angular momenta in laboratory system X-component'
        end if
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CAX(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')" 
	  go to 210
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CAX(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*             
 210  continue
	end if
      END DO! OF ITYP
	end if
*
*    2.11 Angular momenta Y-component:
*    ---------------------------------
 	if(ITCF(11).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(11),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL2)then
	  write(25,*)' Angular momenta in principal system Y-component'
        else
	  write(25,*)' Angular momenta in laboratory system Y-component'
        end if
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CAY(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 211
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CAY(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*              
 211	continue
	end if
      END DO! OF ITYP
	end if
*
*    2.12  Angular momenta Z-component:
*    ----------------------------------
      if(ITCF(12).ne.0)then
        write(25,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(12),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL2)then
	  write(25,*)' Angular momenta in principal system Z-component'
        else
	  write(25,*)' Angular momenta in laboratory system Z-component'
        end if
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      write(25,*)'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',
     +NAME(ITYP)
      write(25,*)
      write(25,*)'     TIME       TCF          UNNORM.'
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      JCF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CAZ(JCF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')" 
	  go to 212
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CAZ(JCF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      write(25,'(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)')TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      write(25,'(a,f12.6,a)')' *** CORRELATION TIME:  ',TAU,'  PS'
*
 212	continue
	end if
      END DO! OF ITYP
      end if 
      stop
      END
*
*============ GETTCF ===================================================
*
*   3.  Put data in arrays for tcf calculations
*   -------------------------------------------
      SUBROUTINE GETTCF
      include "tranal.h"
      include "trtcf.h"
      character*16 filtcf
*
      if(NSTCF.le.0)return
*
      NOP3    = NOP*3
*
*   3.1 Linear momenta:
*     
	if(ITCF(NTCF+1).ne.0)then
	  IXCF=1
          DO ITYP = 1,NTYPES
            NSB     = IADDR(ITYP)
            NSP     = NSPEC(ITYP)
	    IXCC    = NSB+1
            if(ITCF(1).ne.3)then
              CALL COPY2(CC1(IXCF),PX(IXCC),NSP)
              IYCF    = IXCF+NSP
              CALL COPY2(CC1(IYCF),PY(IXCC),NSP)
              IZCF    = IYCF+NSP
              CALL COPY2(CC1(IZCF),PZ(IXCC),NSP)
            else
              CALL COPY2(CC1(IXCF),PX(IXCC),NSP)
              IYCF    = IXCF+NSP
              CALL COPY2(CC1(IYCF),PY(IXCC),NSP)
              IZCF    = IYCF+NSP
              CALL COPY2(CC1(IZCF),PZ(IXCC),NSP)
            end if
            if(IPRINT.ge.7)then
               do II=1,NSP
               III=IXCC+II
               write(*,'(2I5,3f12.8)')IXCC,II,PX(III),PY(III),PZ(III)
               end do
            end if
            if(LMOL1)then
              CALL ROTMAT(PX,PY,PZ,IXCC,NSP,.TRUE.) ! Rotate to princip coords
              CALL COPY2(CC5(IXCF),PX(IXCC),NSP)
              CALL COPY2(CC5(IYCF),PY(IXCC),NSP)
              CALL COPY2(CC5(IZCF),PZ(IXCC),NSP)
            if(IPRINT.ge.7)then
               do II=0,NSP-1
               III=IXCC+II
             write(*,'(a,I5,3f12.8)')' rot',III,PX(III),PY(III),PZ(III)
               end do
            end if
            end if
	    IXCF    = IXCF+NSP*3
          END DO! OF ITYP  
	end if
*
*	3.2 Angular momenta:
*       --------------------
      if(ITCF(NTCF+2).ne.0)then
        IXCF=1
        DO ITYP = 1,NTYPES
          NSB     = IADDR(ITYP)
          NSP     = NSPEC(ITYP)
          IXCC    = NSB +1
          CALL COPY2(CC2(IXCF),QX(IXCC),NSP)
          IYCF    = IXCF+NSP
          CALL COPY2(CC2(IYCF),QY(IXCC),NSP)
          IZCF    = IYCF+NSP
          CALL COPY2(CC2(IZCF),QZ(IXCC),NSP)
          if(LMOL2)then
            CALL ROTMAT(QX,QY,QZ,IXCC,NSP,.TRUE.) ! Rotate to princip coords
            CALL COPY2(CC6(IXCF),QX(IXCC),NSP)
            CALL COPY2(CC6(IYCF),QY(IXCC),NSP)
            CALL COPY2(CC6(IZCF),QZ(IXCC),NSP)
          end if
          IXCF=IXCF+NSP*3
        END DO! OF ITYP
      end if
*
*    3.3 Reorientational motion - dipole moment
*    ------------------------------------------
      if(ITCF(NTCF+3).ne.0)then
        IXCF=1
        DO ITYP = 1,NTYPES
          NSB     = IADDR(ITYP)
          NSP     = NSPEC(ITYP)
          IXCC    = NSB +1
          CALL COPY2(CC3(IXCF),DPX(IXCC),NSP)
          IYCF    = IXCF+NSP
          CALL COPY2(CC3(IYCF),DPY(IXCC),NSP)
          IZCF    = IYCF+NSP
          CALL COPY2(CC3(IZCF),DPZ(IXCC),NSP)
	  IXCF=IXCF+3*NSP
*
        END DO! OF ITYP
      end if
*
*	3.4 Reorientational motion - spec. unit vector
*       ----------------------------------------------
	if(ITCF(NTCF+4).ne.0)then
      CALL GETVEC
*
      IXCF=1
      DO ITYP = 1,NTYPES
      NSB     = IADDR(ITYP)
      NSP     = NSPEC(ITYP)
      IXCC    = NSB +1
      CALL COPY2(CC4(IXCF),UPX(IXCC),NSP)
      IYCF    = IXCF+NSP
      CALL COPY2(CC4(IYCF),UPY(IXCC),NSP)
      IZCF    = IYCF+NSP
      CALL COPY2(CC4(IZCF),UPZ(IXCC),NSP)
	IXCF=IXCF+3*NSP
*
      END DO! OF ITYP
	end if
*
      CALL ADDTCF
*
      RETURN
      END
*
*============ ADDTCF ====================================================
*
      SUBROUTINE ADDTCF
      include "tranal.h"
      include "trtcf.h"
*
      NCALLS     = NCALLS+1
      LIS        = MOD(NCALLS-1,NSTEG)+1
      LST        = NOM+(NSTEG-LIS)*NOM
*
      if(ITCF(NTCF+1).ne.0)CALL COPY1(CC1(LST+1),CC1(1),NOM)
      if(ITCF(NTCF+2).ne.0)CALL COPY1(CC2(LST+1),CC2(1),NOM)
      if(ITCF(NTCF+3).ne.0)CALL COPY1(CC3(LST+1),CC3(1),NOM)
      if(ITCF(NTCF+4).ne.0)CALL COPY1(CC4(LST+1),CC4(1),NOM)
      if(ITCF(NTCF+5).ne.0)CALL COPY1(CC5(LST+1),CC5(1),NOM)
      if(ITCF(NTCF+6).ne.0)CALL COPY1(CC6(LST+1),CC6(1),NOM)
*
      DO J       = 1,LIS
        INX(J)     = INX(J)+1
      END DO! OF J
*
* Vectors in box frame:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      JCF        = NSTEG*(ITYP-1)+1
      NSP        = NSPEC (ITYP  )
      NSB        = IXCF+1
      NSE        = IXCF+NSP*3
      DO I       = NSB,NSE
      L          = I+LST
      if(ITCF(1).ne.0)
     +CALL SUMCF(LIS,CC1(I),CC1(L),NOM,CFV(JCF),NSTEG)
      if(ITCF(2).ne.0)
     +CALL SUMCF(LIS,CC2(I),CC2(L),NOM,CFA(JCF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* Specific unit vectors:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      JCF        = NSTEG*(ITYP-1)+1
      NSP        = NSPEC (ITYP  )
      NSB        = IXCF+1
      NSE        = IXCF+NSP
      DO I       = NSB,NSE
      J          = I+NSP
      K          = J+NSP
      L          = I+LST
      M          = J+LST
      N          = K+LST
      if(ITCF(NTCF+3).ne.0)
     +CALL SUMPL(LIS,CC3(I),CC3(J),CC3(K),CC3(L),CC3(M),CC3(N),NOM
     X          ,CP1(JCF),CP2(JCF),NSTEG)
      if(ITCF(NTCF+4).ne.0)
     +CALL SUMPL(LIS,CC4(I),CC4(J),CC4(K),CC4(L),CC4(M),CC4(N),NOM
     X          ,CU1(JCF),CU2(JCF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* X Y Z projections:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      JCF        = NSTEG*(ITYP-1)+1
      NSP        = NSPEC (ITYP  )
      NSB        = IXCF+1
      NSE        = IXCF+NSP
      DO I       = NSB,NSE
      J          = I+NSP
      K          = J+NSP
      L          = I+LST
      M          = J+LST
      N          = K+LST
      if(LMOL1)then
       CALL SUMPS(LIS,CC5(I),CC5(J),CC5(K),CC5(L),CC5(M),CC5(N),NOM
     X          ,CVX(JCF),CVY(JCF),CVZ(JCF),NSTEG)
      else
        if(ITCF(NTCF+1).ne.0)
     +  CALL SUMPS(LIS,CC1(I),CC1(J),CC1(K),CC1(L),CC1(M),CC1(N),NOM
     X          ,CVX(JCF),CVY(JCF),CVZ(JCF),NSTEG)
      end if
      if(LMOL2)then
       CALL SUMPS(LIS,CC6(I),CC6(J),CC6(K),CC6(L),CC6(M),CC6(N),NOM
     X          ,CAX(JCF),CAY(JCF),CAZ(JCF),NSTEG)
      else
        if(ITCF(NTCF+2).ne.0)
     +  CALL SUMPS(LIS,CC2(I),CC2(J),CC2(K),CC2(L),CC2(M),CC2(N),NOM
     X          ,CAX(JCF),CAY(JCF),CAZ(JCF),NSTEG)
      end if
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
      IF (NCALLS.LE.NSTEG) RETURN
*
      LEFT       = NSTEG-LIS
*
      IF(LEFT.LE.0) RETURN
*
      DO J       = 1,LEFT
      INX(LIS+J) = INX(LIS+J)+1  
      END DO! OF J
*
* Vectors in box frame:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      JCF        = (ITYP-1)*NSTEG+LIS+1
      NSP        = NSPEC(ITYP)
      NSB        = IXCF+1
      NSE        = IXCF+NSP*3
      DO I       = NSB,NSE
      L          = I+NOM
      if(ITCF(1).ne.0)CALL SUMCF(LEFT,CC1(I),CC1(L),NOM,CFV(JCF),NSTEG)
      if(ITCF(2).ne.0)CALL SUMCF(LEFT,CC2(I),CC2(L),NOM,CFA(JCF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* Specific unit vectors:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      JCF        = (ITYP-1)*NSTEG+LIS+1
      NSP        = NSPEC(ITYP)
      NSB        = IXCF+1
      NSE        = IXCF+NSP
      DO I       = NSB,NSE
      J          = I+NSP
      K          = J+NSP
      L          = I+NOM
      M          = J+NOM
      N          = K+NOM 
	if(ITCF(NTCF+3).ne.0)
     +CALL SUMPL(LEFT,CC3(I),CC3(J),CC3(K),CC3(L),CC3(M),CC3(N),NOM
     X          ,CP1(JCF),CP2(JCF),NSTEG)
	if(ITCF(NTCF+4).ne.0)
     +CALL SUMPL(LEFT,CC4(I),CC4(J),CC4(K),CC4(L),CC4(M),CC4(N),NOM
     X          ,CU1(JCF),CU2(JCF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* X Y Z projections
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
        JCF        = (ITYP-1)*NSTEG+LIS+1
        NSP        = NSPEC(ITYP)
        NSB        = IXCF+1
        NSE        = IXCF+NSP
        DO I       = NSB,NSE
          J          = I+NSP
          K          = J+NSP
          L          = I+NOM
          M          = J+NOM
          N          = K+NOM
          if(LMOL1)then
       CALL SUMPS(LEFT,CC5(I),CC5(J),CC5(K),CC5(L),CC5(M),CC5(N),NOM
     X          ,CVX(JCF),CVY(JCF),CVZ(JCF),NSTEG)
          else
            if(ITCF(NTCF+1).ne.0)
     +  CALL SUMPS(LEFT,CC1(I),CC1(J),CC1(K),CC1(L),CC1(M),CC1(N),NOM
     X          ,CVX(JCF),CVY(JCF),CVZ(JCF),NSTEG)
          end if
          if(LMOL2)then
       CALL SUMPS(LEFT,CC6(I),CC6(J),CC6(K),CC6(L),CC6(M),CC6(N),NOM
     X          ,CAX(JCF),CAY(JCF),CAZ(JCF),NSTEG)
          else
            if(ITCF(NTCF+2).ne.0)
     +  CALL SUMPS(LEFT,CC2(I),CC2(J),CC2(K),CC2(L),CC2(M),CC2(N),NOM
     X          ,CAX(JCF),CAY(JCF),CAZ(JCF),NSTEG)
          end if
        end do
        IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
      RETURN
      END
*
*========= SUMCF ===============================================
*
      SUBROUTINE SUMCF(N,X,AX,INC,CF,NSTEG)
      IMPLICIT real*8 (A-H,O-Z)
	real*4 X,AX
      DIMENSION   AX(*),CF(*)
*
      II       = 1
      DO 100 I = 1,N
      CF(I)    = CF(I)+X*AX(II)
      II       = II + INC
  100 CONTINUE
      RETURN
      END
*
*========= SUMPL ===============================================
*
      SUBROUTINE SUMPL(N,X,Y,Z,AX,AY,AZ,INC,CF1,CF2,NSTEG)
      IMPLICIT real*8 (A-H,O-Z)
	real*4 X,Y,Z,AX,AY,AZ
      DIMENSION   AX(*),AY(*),AZ(*),CF1(*),CF2(*)
*
      II       = 1
      DO 100 I = 1,N
      T        = X*AX(II)+Y*AY(II)+Z*AZ(II)
      CF1(I)   = CF1(I)+T
      CF2(I)   = CF2(I)+3.D0*T*T-1.D0
      II       = II + INC
  100 CONTINUE
      RETURN
      END
*
*========= SUMPS ===============================================
*
      SUBROUTINE SUMPS(N,X,Y,Z,AX,AY,AZ,INC,CF1,CF2,CF3,NSTEG)
      IMPLICIT real*8 (A-H,O-Z)
      real*4 X,Y,Z,AX,AY,AZ
      DIMENSION   AX(*),AY(*),AZ(*),CF1(*),CF2(*),CF3(*)
*
      II       = 1
      DO 100 I = 1,N
      CF1(I)   = CF1(I)+X*AX(II)
      CF2(I)   = CF2(I)+Y*AY(II)
      CF3(I)   = CF3(I)+Z*AZ(II)
      II       = II + INC
  100 CONTINUE
      RETURN
      END

*
*============ GETVEC =============================================
*
      SUBROUTINE GETVEC
      include "tranal.h"
      include "trtcf.h"
*
      DO ITYP = 1,NTYPES
      I1      = N1(ITYP)
      I2      = N2(ITYP)
      NSB     = IADDR(ITYP)+1
      NSE     = IADDR(ITYP +1)
      NSS     = NSITS(ITYP)
      DO I    = NSB,NSE      
	IAT     = ISADDR(ITYP)+(I-NSB)*NSS
         IS      = IAT+I1
         JS      = IAT+I2
         UTX     = WX(IS)-WX(JS)
         UTY     = WY(IS)-WY(JS)
         UTZ     = WZ(IS)-WZ(JS)
         UNORM   = DSQRT(UTX**2+UTY**2+UTZ**2)
	if(UNORM.gt.1.d-10)then
         UPX(I)  = UTX/UNORM
         UPY(I)  = UTY/UNORM
         UPZ(I)  = UTZ/UNORM 
	else
	UPX(I)=0.d0 
	UPY(I)=0.d0 
	UPZ(I)=0.d0
	end if
      END DO! OF I
      END DO! OF ITYP
*
      RETURN
      END
*
*================= COPY1 =================================================
*
      SUBROUTINE COPY1(A,B,M)
      real*4 A(M),B(M)
      DO I   = 1,M
       A(I)  = B(I)
      END DO! OF I
      RETURN
      END
*
*================= COPY2 =================================================
*
      SUBROUTINE COPY2(A,B,M)
      real*4 A(M)
	real*8 B(M)
      DO I   = 1,M
       A(I)  = B(I)
      END DO! OF I
      RETURN
      END
*========= ROTMAT =================================================
*
*   12. Recalculate coordinates in principal coord. system
*   ------------------------------------------------------
      SUBROUTINE ROTMAT(XX,YY,ZZ,IBG,N,NBTOP)
*
      include "tranal.h"
      include "trtcf.h"
*
      DIMENSION XX(*),YY(*),ZZ(*)
      LOGICAL NBTOP
*
      IBEG  = IBG
      IEND  = IBG+N-1
*
      IF(NBTOP) THEN	! from box to principal:
        DO I  = IBEG,IEND
          AQ    = XX(I)
          BQ    = YY(I)
          CQ    = ZZ(I)
          XX(I) = RMX(I,1,1)*AQ+RMX(I,2,1)*BQ+RMX(I,3,1)*CQ
          YY(I) = RMX(I,1,2)*AQ+RMX(I,2,2)*BQ+RMX(I,3,2)*CQ
          ZZ(I) = RMX(I,1,3)*AQ+RMX(I,2,3)*BQ+RMX(I,3,3)*CQ
        END DO! OF I
      ELSE 		! from principal to box:
        DO I  = IBEG,IEND
          AQ    = XX(I)
          BQ    = YY(I)
          CQ    = ZZ(I)
          XX(I) = RMX(I,1,1)*AQ+RMX(I,1,2)*BQ+RMX(I,1,3)*CQ
          YY(I) = RMX(I,2,1)*AQ+RMX(I,2,2)*BQ+RMX(I,2,3)*CQ
          ZZ(I) = RMX(I,3,1)*AQ+RMX(I,3,2)*BQ+RMX(I,3,3)*CQ
        END DO! OF I
      END IF
*
      RETURN
      END
*
*=============== GETROT ==========================================
*
*  11. Calculate rotation matrix to principal coord. system
*  --------------------------------------------------------
      SUBROUTINE GETROT
*
      include "tranal.h"
      include "trtcf.h"
      DIMENSION RX(NS),RY(NS),RZ(NS),FMASS(NS)
      DIMENSION FMOI(3,3)      
      DIMENSION G(3),H(3)
*
      N        = 0
      DO ITYP  = 1,NTYPES
        NSS      = NSITS (ITYP)
        NSP      = NSPEC (ITYP)
        ISB      = ISADR (ITYP)+1
        ISE      = ISADR (ITYP +1)
        JBE      = IADDR (ITYP)+1
        JEN      = IADDR (ITYP +1)
        DO J     = JBE,JEN
          I        = 0
          DO IS    = ISB,ISE
            I        = I+1
            N        = N+1
            RX   (I) = WX(N)
            RY   (I) = WY(N)
            RZ   (I) = WZ(N)
            FMASS(I) = MASS(IS)
          END DO! OF IS
*
          CALL GETMOI(RX,RY,RZ,FMASS,FMOI,NSS)
*
          CALL HH3BY3(FMOI,G,H)
*
* ordering eigen values G1>G2>G3
          call ORDER3X3(FMOI,G)
*
          DO K       = 1,3
            DO L       = 1,3
              RMX(J,K,L) = FMOI(K,L)
            END DO! OF K
          END DO! OF L
*
        END DO! OF J
*
      END DO! OF ITYP
*
      RETURN
      END
*
*========== ORDER3x3
*
      subroutine ORDER3X3(A,G)
      real*8 A(3,3),G(3),B(3,3),H
      integer IO(3)
      IO(1)=1
      IO(2)=2
      IO(3)=3
      if(G(3).gt.G(2))then
         IO(2)=3                  !   G
         IO(3)=2
         H=G(3)
         G(3)=G(2)
         G(2)=H
         if(G(2).gt.G(1))then
           G(2)=G(1)
           G(1)=H
           IO(1)=3
           IO(2)=1
           if(G(3).gt.G(2))then
             H=G(3)
             G(3)=G(2)
             G(2)=H
             IO(2)=2
             IO(3)=1
           end if
         end if
       else
         if(G(2).gt.G(1))then
           H=G(1)
           G(1)=G(2)
           G(2)=H
           IO(1)=2
           IO(2)=1
           if(G(3).gt.G(2))then
             G(2)=G(3)
             G(3)=H
             IO(2)=3
             IO(3)=1
           end if
         end if
       end if
       if(G(1).lt.G(2).or.G(2).lt.G(3))write(*,*)' wrong G order:',G
       do I=1,3
         do J=1,3
           JJ=IO(J)
           B(I,J)=A(I,JJ)
         end do
       end do
       do I=1,3
          do J=1,3
             A(I,J)=B(I,J)
           end do
       end do
       return
       end
*
*=============== GETMOI ========================================
*
      SUBROUTINE GETMOI(X,Y,Z,M,MOMENT,NSS)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 M,MOMENT
      DIMENSION X(NSS),Y(NSS),Z(NSS),M(NSS)
      DIMENSION MOMENT(3,3),R(3)
*
      DO INERT           = 1,3
        DO IA              = 1,3
          TT                 = 0.D0
          MOMENT(INERT,IA)   = 0.D0
          IF(INERT.EQ.IA) TT = 1.D0
          DO I               = 1,NSS
            R(1)               =  X(I)
            R(2)               =  Y(I)
            R(3)               =  Z(I)
            SQS                =  X(I)**2+Y(I)**2+Z(I)**2
      MOMENT(INERT,IA)   = MOMENT(INERT,IA)+M(I)*(TT*SQS-R(INERT)*R(IA))
          END DO! OF I
        END DO! OF IA
      END DO! OF INERT
*
      RETURN
      END
*
*========== HH3BY3 ============================================
*
      SUBROUTINE HH3BY3(A,D,E)
*
      IMPLICIT real*8 (A-H,O-Z)
*
      DIMENSION A(3,3),D(3),E(3)
      data ICOUNT/0/
*
      H       = 0.
      SCALE   = ABS(A(3,1))+ABS(A(3,2))
      IF(SCALE.EQ.0) THEN
        E(3)    = A(3,2)
      ELSE
        A(3,1)  = A(3,1)/SCALE
        A(3,2)  = A(3,2)/SCALE
        H       = A(3,1)**2+A(3,2)**2+H
        F       = A(3,2)
        G       =-SIGN(SQRT(H),F)
        E(3)    = SCALE*G
        H       = H-F*G
        A(3,2)  =   F-G
*
        A(1,3)  = A(3,1)/H
        G       = A(1,1)*A(3,1)+A(2,1)*A(3,2)
        E(1)    = G/H
        F       = E(1)*A(3,1)
*
        A(2,3)  = A(3,2)/H
        G       = A(2,1)*A(3,1)+A(2,2)*A(3,2)
        E(2)    = G/H
        F       = F+E(2)*A(3,2)
*
        HH      = F/(H+H)
        F       = A(3,1)
        G       = E(1)-HH*F
        E(1)    = G
        A(1,1)  = A(1,1)-F*E(1)-G*A(3,1)
        F       = A(3,2)
        G       = E(2)-HH*F
        E(2)    = G
        A(2,1)  = A(2,1)-F*E(1)-G*A(3,1)
        A(2,2)  = A(2,2)-F*E(2)-G*A(3,2)
      END IF
      D(3)    = H
*
      SCALE   = 0.
      E(2)    = A(2,1)
      D(2)    = 0.
*
      E(1)    = 0.
      D(1)    = A(1,1)
      A(1,1)  = 1.
      IF(D(2).NE.0) THEN
        G       = A(2,1)*A(1,1)
        A(1,1)  = A(1,1)-A(1,2)*G
      END IF
      D(2)    = A(2,2)
      A(2,2)  = 1.
      A(2,1)  = 0.
      A(1,2)  = 0.
      IF(D(3).NE.0) THEN
        G       = A(3,1)*A(1,1)+A(3,2)*A(2,1)
        A(1,1)  = A(1,1)-A(1,3)*G
        A(2,1)  = A(2,1)-A(2,3)*G
        G       = A(3,1)*A(1,2)+A(3,2)*A(2,2)
        A(1,2)  = A(1,2)-A(1,3)*G
        A(2,2)  = A(2,2)-A(2,3)*G
      END IF
      D(3)    = A(3,3)
      A(3,3)  = 1.
      A(3,1)  = 0.
      A(1,3)  = 0.
      A(3,2)  = 0.
      A(2,3)  = 0.
*
*Diagonalize it:
*
      E(1)    = E(2)
      E(2)    = E(3)
      E(3)    = 0.
      DO 15 L = 1,3
      ITER    = 0
    1 CONTINUE
      DO 12 M = L,2
      DD      = ABS(D(M))+ABS(D(M+1))
      IF(ABS(E(M))+DD.EQ.DD) GO TO 2
   12 CONTINUE
      M       = 3
    2 CONTINUE
      IF(M.NE.L) THEN
      IF(ITER.EQ.30) then
        write(*,*) '!!! TOO MANY ITERATIONS!    in HH3BY3'
        ICOUNT=ICOUNT+1
        if(ICOUNT.gt.1000)stop
        return
      end if
      ITER    = ITER+1
      G       = (D(L+1)-D(L))/(2.*E(L))
      R       = SQRT(G**2+1.)
      G       = D(M)-D(L)+E(L)/(G+SIGN(R,G))
      S       = 1.
      C       = 1.
      P       = 0.
      DO 14 I = M-1,L,-1
      F       = S*E(I)
      B       = C*E(I)
      IF(ABS(F).GE.ABS(G)) THEN
      C       = G/F
      R       = SQRT(C**2+1.)
      E(I+1)  = F*R
      S       = 1./R
      C       = C*S
      ELSE
      S       = F/G
      R       = SQRT(S**2+1.)
      E(I+1)  = G*R
      C       = 1./R
      S       = S*C
      END IF
      G       = D(I+1)-P
      R       =(D(I)-G)*S+2.*C*B
      P       = S*R
      D(I+1)  = G+P
      G       = C*R-B
      F       = A(1,I+1)
      A(1,I+1)= S*A(1,I)+C*F
      A(1,I)  = C*A(1,I)-S*F
      F       = A(2,I+1)
      A(2,I+1)= S*A(2,I)+C*F
      A(2,I)  = C*A(2,I)-S*F
      F       = A(3,I+1)
      A(3,I+1)= S*A(3,I)+C*F
      A(3,I)  = C*A(3,I)-S*F
   14 CONTINUE
      D(L)    = D(L)-P
      E(L)    = G
      E(M)    = 0.
      GO TO 1
      END IF
   15 CONTINUE
      RETURN
      END
*
*=============== ROTATE ==============================================
*
      SUBROUTINE ROTATE(Q,X,Y,Z)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION Q(4)
      SIGN = -1.D0
      AA   = Q(1)*Q(1)
      BB   = Q(2)*Q(2)
      CC   = Q(3)*Q(3)
      DD   = Q(4)*Q(4)
      AB   = Q(1)*Q(2)
      AC   = Q(1)*Q(3)
      AD   = Q(1)*Q(4)*SIGN
      BC   = Q(2)*Q(3)
      BD   = Q(2)*Q(4)*SIGN
      CD   = Q(3)*Q(4)*SIGN
      RX   = X*(-AA+BB-CC+DD)+Y*(+CD-AB)*2.D0 +Z*(+BC+AD)*2.D0
      RY   = X*(-AB-CD)*2.D0 +Y*(+AA-BB-CC+DD)+Z*(+BD-AC)*2.D0
      RZ   = X*(+BC-AD)*2.D0 +Y*(-AC-BD)*2.D0 +Z*(-AA-BB+CC+DD)
      X    = RX
      Y    = RY
      Z    = RZ
      RETURN
      END
*
*================== JACOBI =============================================
*
*   rotation matrix
*
      SUBROUTINE JACOBI(A,V)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION A(3,3),V(3,3)
      N=3
      RHO=1.0D-8
      TES=0.D0
      SCL=0.D0
      DO 10 I=1,N
      SCL=SCL+A(I,I)*A(I,I)
   10 CONTINUE
      SCL=DSQRT(DABS(SCL))/FLOAT(N)
      DO 20 I=1,N
      DO 20 J=1,N
   20 A(I,J)=A(I,J)/SCL
      DO 30 I=1,N
      DO 30 J=1,N
      V(I,J)=0.D0
      IF(I.EQ.J)V(I,J)=1.D0
   30 CONTINUE
      DO 100 I=2,N
      IU1=I-1
      DO 100 J=1,IU1
  100 TES=TES+2.D0*A(I,J)*A(I,J)
      TES=DSQRT(DABS(TES))
      M=0
  105 TES=TES/FLOAT(N)
      IF(TES.LT.RHO)TES=RHO
  110 DO 165 I=2,N
      IUU1=I-1
      DO 165 J=1,IUU1
      IF(DABS(A(I,J))-TES)165,115,115
  115 M=1
      V1=A(J,J)
      V2=A(I,J)
      V3=A(I,I)
      U=0.5D0*(V1-V3)
      IF(DABS(U)-RHO)120,125,125
  120 OMG=-1.D0
      GO TO 130
  125 SQVVUU=DSQRT(DABS(V2*V2+U*U))
      OMG=-V2/SQVVUU
      IF(U.LT.0.D0)OMG=-OMG
  130 SQOOO=DSQRT(DABS(2.0*(1.00+DSQRT(DABS(1.00-OMG*OMG)))))
      S=OMG/SQOOO
      C=DSQRT(DABS(1.00-S*S))
      DO 160 K=1,N
      IF(K-I)140,135,135
  135 TEM=A(K,J)*C-A(K,I)*S
      A(K,I)=A(K,J)*S+A(K,I)*C
      A(K,J)=TEM
      GO TO 155
  140 IF(K-J)145,150,150
  145 TEM=A(J,K)*C-A(I,K)*S
      A(I,K)=A(J,K)*S+A(I,K)*C
      A(J,K)=TEM
      GO TO 155
  150 TEM=A(K,J)*C-A(I,K)*S
      A(I,K)=A(K,J)*S+A(I,K)*C
      A(K,J)=TEM
  155 TEM=V(K,J)*C-V(K,I)*S
      V(K,I)=V(K,J)*S+V(K,I)*C
      V(K,J)=TEM
  160 CONTINUE
      A(J,J)=V1*C*C+V3*S*S-2.D0*V2*S*C
      A(I,I)=V1*S*S+V3*C*C+2.D0*V2*S*C
      A(I,J)=(V1-V3)*S*C+V2*(C*C-S*S)
  165 CONTINUE
      IF(M-1)175,170,170
  170 M=0
      GO TO 110
  175 IF(TES-RHO)180,180,105
  180 DO 190 I=1,N
      DO 190 J=1,N
  190 A(I,J)=SCL*A(I,J)
      RETURN
      END
*
*=============== DIPOLE ===========================================
*                                                                
*   9. Calculate dipole moment
*   --------------------------
      SUBROUTINE DIPOLE
      include "tranal.h"
      include "trtcf.h"
*
      FAK	  = 1.0/0.52917
      FAC         = 2.5418
      N           = 0
      DO ITYP     = 1,NTYPES
      NSP         = NSPEC(ITYP)
	NSS	    = NSITS(ITYP) 
      FNSPI       = 1.D0/DFLOAT(NSP)
      DIPMOM      = 0.D0
      ISB         = ISADR(ITYP)+1
      ISE         = ISADR(ITYP +1)
      IBE         = IADDR(ITYP)+1
      IEN         = IADDR(ITYP +1)
      DO I        = IBE,IEN
	if(NSS.gt.1)then
      DMX         = 0.D0
      DMY         = 0.D0
      DMZ         = 0.D0
      DO IS       = ISB,ISE
      N           = N+1
Calculate permanent dipole moment (magnitude)
      DMX         = DMX+CHARGE(IS)*WX(N)
      DMY         = DMY+CHARGE(IS)*WY(N)
      DMZ         = DMZ+CHARGE(IS)*WZ(N)
      END DO! OF IS
      DNORM       = DSQRT(DMX**2+DMY**2+DMZ**2)
      DIPMOM      = DIPMOM+DNORM*FAC*FAK
      IF(DNORM.GT.1.d-10) then
      DPX(I)      = DMX/DNORM
      DPY(I)      = DMY/DNORM
      DPZ(I)      = DMZ/DNORM
	else
	DPX(I)=0.d0
	DPY(I)=0.d0
	DPZ(I)=0.d0
	end if   
	else
	DIPMOM=0.d0
	end if
      END DO! OF I
      DIPMOM      = DIPMOM*FNSPI
      IF(MOD(NSTEP,10).EQ.0.and.IPRINT.ge.8) THEN
      IF(DABS(DIPMOM).GT.0.00001D0) THEN
      PRINT "(' DIPOLE MOMENT(D):',F8.3,' FOR TYPE  ',A6)"
     X,DIPMOM,NAME(ITYP)
      END IF
      END IF          
      END DO! OF ITYP
      RETURN
      END   

