C========================================
C
C     MDynaMix v.4.3
C
*     PART 9
*
*     File tcf.f
*     -------------
*
C     This file contains subroutines responsible for calculation 
C     of time correlation functions
*   
*
*============ TCFINP ===================================================
*
*    1.1 Initialize working arrays for tcf calculations
*
      SUBROUTINE TCFINP
      include "prcm.h"
*
      DIMENSION M1(NTPS),M2(NTPS),JTCF(NTCF)
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
      if(TASKID.eq.MAST)then
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
      end if
*
      IF(NSTEG.GT.MAXCF)            STOP '!!! NSTEG.GT.MAXCF'
      IF((NOM*(NSTEG+2)).GT.MAXTCF) STOP '!!! INCREASE MAXTCF !'
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
*
*----------------------------------------------------------------------
*
      IF(.NOT.LCFRST) RETURN
*
      LUIN=50
      OPEN (UNIT=LUIN,file=ftcf,FORM='UNFORMATTED',STATUS='OLD',err=80)
*
      READ(LUIN) CT,MSTEG,MSTCF,IUMP,MTYPES,MOP,MOM,JTCF
*
                           ICHECK = 0
      IF(CT    .NE.DT    ) ICHECK = ICHECK+1
      IF(MSTEG .NE.NSTEG ) ICHECK = ICHECK+1
      IF(MSTCF .NE.NSTCF ) ICHECK = ICHECK+1
      IF(IUMP  .NE.JUMP  ) ICHECK = ICHECK+1
      IF(MOM   .NE.NOM   ) ICHECK = ICHECK+1
      IF(MTYPES.NE.NTYPES) ICHECK = ICHECK+1
	do I=1,NTCF
	  if(ITCF(I).ne.JTCF(I)) then
	    WRITE(*,*)'!!! Wrong tcf num.',I
	    ICHECK=ICHECK+1                 
	  end if
	end do
*             
	if(TASKID.eq.MAST)then
      PRINT *
      PRINT "('... RESTART DT: ',D10.3
     X      ,' ... CURRENT DT: ',D10.3)",CT,DT
      PRINT "('... RESTART NSTEG: ',I6
     X      ,' ... CURRENT NSTEG: ',I6)",MSTEG,NSTEG
      PRINT "('... RESTART NSTCF  ',I6
     X      ,' ... CURRENT NSTCF: ',I6)",MSTCF,NSTCF
      PRINT "('... RESTART JUMP:  ',I6
     X      ,' ... CURRENT JUMP:  ',I6)",IUMP,JUMP
      PRINT "('... RESTART NOM:   ',I6
     X      ,' ... CURRENT NOM:   ',I6)",MOM,NOM
      PRINT "('... RESTART NOP:   ',I6
     X      ,' ... CURRENT NOP:   ',I6)",MOP,NOP
      PRINT "('... RESTART NTYPES ',I6
     X      ,' ... CURRENT NTYPES ',I6)",MTYPES,NTYPES
      PRINT *            
	end if
*
      IF(ICHECK.NE.0) STOP '!!! WRONG TCF RESTART FILE ! '
*     
      READ (LUIN) M1,M2
      READ (LUIN) NCALLS,INX
      if(ITCF(1).ne.0)READ (LUIN) CFV
      if(ITCF(2).ne.0)READ (LUIN) CFA
      if(ITCF(3).ne.0)READ (LUIN) CP1
      if(ITCF(4).ne.0)READ (LUIN) CP2
      if(ITCF(5).ne.0)READ (LUIN) CU1
      if(ITCF(6).ne.0)READ (LUIN) CU2
      if(ITCF(7).ne.0)READ (LUIN) CVX
      if(ITCF(8).ne.0)READ (LUIN) CVY
      if(ITCF(9).ne.0)READ (LUIN) CVZ
      if(ITCF(10).ne.0)READ (LUIN) CAX
      if(ITCF(11).ne.0)READ (LUIN) CAY
      if(ITCF(12).ne.0)READ (LUIN) CAZ
      if(ITCF(NTCF+1).ne.0)READ (LUIN) (CC1(I),I=1,NNTCF)
      if(ITCF(NTCF+2).ne.0)READ (LUIN) (CC2(I),I=1,NNTCF)
      if(ITCF(NTCF+3).ne.0)READ (LUIN) (CC3(I),I=1,NNTCF)
      if(ITCF(NTCF+4).ne.0)READ (LUIN) (CC4(I),I=1,NNTCF)
      if(ITCF(NTCF+5).ne.0)READ (LUIN) (CC5(I),I=1,NNTCF)
      if(ITCF(NTCF+6).ne.0)READ (LUIN) (CC6(I),I=1,NNTCF)
      CLOSE(LUIN)
*     
	if(ITCF(NTCF+1).ne.0)then 
	  IXCF=0
        DO ITYP  = 1,NTYPES
          NSP      = NSPEC(ITYP)
          NSB      = IADDR(ITYP)+1
          NSE      = IADDR(ITYP +1)
          VSQ      = 0.D0
          DO I     = NSB,NSE
	    IXC = IXCF+I
            VSQ = VSQ+CC1(IXC)**2+CC1(IXC+NSP)**2+CC1(IXC+2*NSP)**2
          END DO! OF I
	    IXCF=IXCF+3*NSP
          VSQ = VSQ/(SUMMAS(ITYP)/UNITM/(AVSNO*1000.D0))/DFLOAT(NSP)
          TTT = 0.5*VSQ*ENERF*TEMPF
      PRINT "(/'*** TEMPERATURE FROM TCF RESTART FILE:',F12.4,' K')",
     +TTT
*
        END DO! OF ITYP
	end if
*
      PRINT "(/,'*** TCF RESTART FILE READ FROM UNIT ',I3)",LUIN
      PRINT "('*** CC ARRAYS HAVE DIMENSIONS: ',I10 )",(MSTEG+2)*MOM
      PRINT "('*** CF ARRAYS HAVE DIMENSIONS: ',I10/)", MSTEG*MTYPES
*
      RETURN 
 80	write(*,*)' restart TCF file ',FTCF,' not found'
      END
*
*========== TCFOUT ==================================================
*
*   2.  TCF output
*   --------------
      SUBROUTINE TCFOUT(IMODE)
      include"prcm.h"   
      if(TASKID.ne.MAST)return
      if(IMODE.eq.0)go to 1000
*
*    2.1. Linear momenta:
*    -------------------- 
      if(ITCF(1).ne.0)then
        print "(/1X,33('-=-'))"
        DEL      = FLOAT(JUMP)*DT*1.D12
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(1),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
*
        DO ITYP  = 1,NTYPES
*
      PRINT "(/1X,'*** VELOCITY AUTOCORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
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
      FMASS    = FMASS+MASSD(IS)
      END DO! OF IS
      ICF      = (ITYP-1)*NSTEG
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI3
      CNRM     = CFV(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 201
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI3
      T1       = CFV(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(1X,' *** CORRELATION TIME:  ',F12.6,'  PS')",TAU
      DV       = BOLTZ*TRTEMP/(FMASS*UNITM)*TAU*1.E-12
      PRINT "(1X,' *** DIFFUSION COEFFICIENT: ',E12.4,' m**2/s'/)",DV
*                          
 201  continue
      END DO! OF ITYP                                   
	end if
*
*    2.2 Angular momenta:
*    --------------------
      if(ITCF(2).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(2),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	write(*,*)' Angular velocity autocorrelation function'
*
        DO ITYP = 1,NTYPES
*     
	  if(NSITS(ITYP).gt.1)then
      PRINT "(/1X,'*** ANGULAR VELOCITY AUTOCORRELATION FUNCTION FOR ',
     +A6)",NAME(ITYP)
*
            PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI3   = 1.D0/DFLOAT(NSP)/3.D0
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI3
      CNRM     = CFA(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 202
	  end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI3
      T1       = CFA(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'   PS'/)",TAU
*     
 202  continue
	end if
      END DO! OF ITYP
	end if
*
*    2.3 Reorientational motion of dipole axis I order Legendre:
*    -----------------------------------------------------------
 	if(ITCF(3).ne.0)then
          write(*,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(3),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
*
      DO ITYP = 1,NTYPES   
	if(NSITS(ITYP).gt.1)then
*
      PRINT"(/1X,
     +'*** REORINTATIONAL DIPOLE CORRELATION FUNCTION 1 FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CP1(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 203
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CP1(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*
 203  continue
	end if
      END DO! OF ITYP
	end if
*
*    2.4 Reorientational motion of dipole axis II order Legendre:
*    ------------------------------------------------------------
 	if(ITCF(4).ne.0)then
          write(*,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(4),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	  write(*,*)
     +'Reorientational motion of dipole axis II order Legendre'
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      PRINT"(/1X,'*** REORIENTATIONAL CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI2   = 1.D0/DFLOAT(NSP)/2.D0	! (3xx-1) divided by 2
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI2
      CNRM     = CP2(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 204
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI2
      T1       = CP2(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*                         
 204      continue
	  end if
        END DO! OF ITYP
      end if
*
*    2.5 Reorientational motion of spec. unit vector I order Legendre:
*    -----------------------------------------------------------------
      if(ITCF(5).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(5),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	write(*,*)
     +'Reorientational motion of spec. unit vector I order Legandre'
*
        DO ITYP = 1,NTYPES
	  if(N1(ITYP).ne.N2(ITYP))then
*
      PRINT "(/1X,'*** REORIENTATIONAL CORRELATION FUNCTIONS FOR '
     +,A6)",NAME(ITYP)
*
          PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
          DO I     = 1,NSTEG
            CCF(I)   = 0.D0
          END DO! OF I
          ICF      = (ITYP-1)*NSTEG
          NSP      = NSPEC(ITYP)
          FNSPI    = 1.D0/DFLOAT(NSP)
          FNF      = 1.D0/FLOAT(INX(1))*FNSPI
          CNRM     = CU1(ICF+1)*FNF
          IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
          IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 205
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CU1(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*                    
 205  continue
	end if
      END DO! OF ITYP
	end if
*
*      2.6 Reorientational motion of spec. unit vector II order Legendre:
*      ------------------------------------------------------------------
 	if(ITCF(6).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)')'*** ',NAMECF(6),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
	write(*,*)
     +'Reorientational motion of spec. unit vector II order Legandre'
*
      DO ITYP = 1,NTYPES
	if(N1(ITYP).ne.N2(ITYP))then
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI2   = 1.D0/DFLOAT(NSP)/2.D0 ! (3xx-1) divided by 2
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI2
      CNRM     = CU2(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 206
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI2
      T1       = CU2(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*                                 
 206  continue
	end if
      END DO! OF ITYP
      end if
*
*    2.7 Linear momenta X-component:
*    -------------------------------
      if(ITCF(7).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(7),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL1)then
	  write(*,*)' Linear momenta in principal system X-component'
        else
	  write(*,*)' Linear momenta in laboratory system X-component'
        end if
*
      DO ITYP = 1,NTYPES
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CVX(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 207
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CVX(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*
 207	continue
      END DO! OF ITYP
	end if
*
*   2.8 Linear momenta Y-component:
*   -------------------------------
      if(ITCF(8).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(8),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL1)then
	  write(*,*)' Linear momenta in principal system Y-component'
        else
	  write(*,*)' Linear momenta in laboratory system Y-component'
        end if
*
      DO ITYP = 1,NTYPES
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CVY(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"  
	  go to 208
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CVY(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*
 208	continue
      END DO! OF ITYP
	end if
*
*    2.9 Linear momenta Z-component:
*    -------------------------------
      if(ITCF(9).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(9),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL1)then
	  write(*,*)' Linear momenta in principal system Z-component'
        else
	  write(*,*)' Linear momenta in laboratory system Z-component'
        end if
*
      DO ITYP = 1,NTYPES
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CVZ(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')" 
	  go to 209
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CVZ(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*
 209  continue 
      END DO! OF ITYP
	end if
*
*    2.10  Angular momenta X-component:
*    ----------------------------------
 	if(ITCF(10).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(10),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL2)then
	  write(*,*)' Angular momenta in principal system X-component'
        else
	  write(*,*)' Angular momenta in laboratory system X-component'
        end if
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CAX(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')" 
	  go to 210
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CAX(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*             
 210  continue
	end if
      END DO! OF ITYP
	end if
*
*    2.11 Angular momenta Y-component:
*    ---------------------------------
 	if(ITCF(11).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(11),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL2)then
	  write(*,*)' Angular momenta in principal system Y-component'
        else
	  write(*,*)' Angular momenta in laboratory system Y-component'
        end if
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(//1X,'     TIME       TCF          UNNORM.'//)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CAY(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')"
	  go to 211
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CAY(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(//1X,' *** CORRELATION TIME:  ',F12.6,'  PS'//)",TAU
*              
 211	continue
	end if
      END DO! OF ITYP
	end if
*
*    2.12  Angular momenta Z-component:
*    ----------------------------------
      if(ITCF(12).ne.0)then
        write(*,'(/1x,a,a6,a,i5,a,f12.6,a)') '*** ',NAMECF(12),
     +    '    TCF DURING ',INX(1),' STEPS,    -> DELTA:',DEL,' ps'
        if(LMOL2)then
	  write(*,*)' Angular momenta in principal system Z-component'
        else
	  write(*,*)' Angular momenta in laboratory system Z-component'
        end if
*
      DO ITYP = 1,NTYPES
	if(NSITS(ITYP).gt.1)then
*
      PRINT "(/1X,'*** TIME CORRELATION FUNCTIONS FOR ',A6)",
     +NAME(ITYP)
*
      PRINT "(/1X,'     TIME       TCF          UNNORM.'/)"
*
      DO I     = 1,NSTEG
      CCF(I)   = 0.D0
      END DO! OF I
      ICF      = (ITYP-1)*NSTEG
      NSP      = NSPEC(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FNF      = 1.D0/FLOAT(INX(1))*FNSPI
      CNRM     = CAZ(ICF+1)*FNF
      IF(CNRM.NE.0.D0) CNRM=1.D0/CNRM
      IF(CNRM.EQ.0.D0) then
	  PRINT "(/1X,'!!! EMPTY TCF ARRAY !!!')" 
	  go to 212
	end if
      DO J     = 1,NSTEG
      IF(INX(J).GT.0) THEN
      TIME     = DFLOAT(J-1)*DEL
      FAC      = 1.D0/DFLOAT(INX(J))*FNSPI
      T1       = CAZ(ICF+J)*FAC
      TCF1     = CNRM*T1
      CCF(J)   = TCF1
      PRINT "(1X,F12.6,1X,F9.6,5X,E12.4,5X,I8)",TIME,TCF1,T1,INX(J)
      END IF!(INX(J...
      END DO! OF J
*
      TAU      = XINTGR(DEL,CCF,NSTEG)
      PRINT "(/1X,' *** CORRELATION TIME:  ',F12.6,'  PS'/)",TAU
*
 212	continue
	end if
      END DO! OF ITYP
	end if 
 	return
*
*   2.13  Dump TCF restart file
*   --------------------------- 
 1000 IF(.NOT.LCFDMP) RETURN
	NNTCF=(MAXCF+2)*3*NOP
*
                LUT=51
      OPEN(UNIT=LUT,file=ftcf,FORM='UNFORMATTED',STATUS='unknown')
      WRITE (LUT) DT,NSTEG,NSTCF,JUMP,NTYPES,NOP,NOM,ITCF
      WRITE (LUT) N1,N2
      WRITE (LUT) NCALLS,INX
      if(ITCF(1).ne.0)WRITE (LUT) CFV
      if(ITCF(2).ne.0)WRITE (LUT) CFA
      if(ITCF(3).ne.0)WRITE (LUT) CP1
      if(ITCF(4).ne.0)WRITE (LUT) CP2
      if(ITCF(5).ne.0)WRITE (LUT) CU1
      if(ITCF(6).ne.0)WRITE (LUT) CU2
      if(ITCF(7).ne.0)WRITE (LUT) CVX
      if(ITCF(8).ne.0)WRITE (LUT) CVY
      if(ITCF(9).ne.0)WRITE (LUT) CVZ
      if(ITCF(10).ne.0)WRITE (LUT) CAX
      if(ITCF(11).ne.0)WRITE (LUT) CAY
      if(ITCF(12).ne.0)WRITE (LUT) CAZ
      if(ITCF(NTCF+1).ne.0)WRITE (LUT) (CC1(I),I=1,NNTCF)
      if(ITCF(NTCF+2).ne.0)WRITE (LUT) (CC2(I),I=1,NNTCF)
      if(ITCF(NTCF+3).ne.0)WRITE (LUT) (CC3(I),I=1,NNTCF)
      if(ITCF(NTCF+4).ne.0)WRITE (LUT) (CC4(I),I=1,NNTCF)
      if(ITCF(NTCF+5).ne.0)WRITE (LUT) (CC5(I),I=1,NNTCF)
      if(ITCF(NTCF+6).ne.0)WRITE (LUT) (CC6(I),I=1,NNTCF)
      CLOSE (LUT)
*
      PRINT "(/,'*** TCF RESTART FILE DUMPED in file ',a/)",ftcf
*
      RETURN
      END
*
*============ GETTCF ===================================================
*
*   3.  Put data in arrays for tcf calculations
*   -------------------------------------------
      SUBROUTINE GETTCF
      include "prcm.h"
*
      IF(.NOT.LCF) RETURN
      if(NSTCF.le.0) return
      IF(MOD(NSTEP,JUMP).NE.0) RETURN
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
            CALL COPY2(CC1(IXCF),PX(IXCC),NSP)
            IYCF    = IXCF+NSP
            CALL COPY2(CC1(IYCF),PY(IXCC),NSP)
            IZCF    = IYCF+NSP
            CALL COPY2(CC1(IZCF),PZ(IXCC),NSP)
            if(LMOL1)then
              CALL ROTMAT(PX,PY,PZ,IXCC,NSP,.TRUE.) ! Rotate to princip coords
              CALL COPY2(CC5(IXCF),PX(IXCC),NSP)
              CALL COPY2(CC5(IYCF),PY(IXCC),NSP)
              CALL COPY2(CC5(IZCF),PZ(IXCC),NSP)
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
      include "prcm.h"
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
      ICF        = NSTEG*(ITYP-1)+1
      NSP        = NSPEC (ITYP  )
      NSB        = IXCF+1
      NSE        = IXCF+NSP*3
      DO I       = NSB,NSE
      L          = I+LST
      if(ITCF(1).ne.0)
     +CALL SUMCF(LIS,CC1(I),CC1(L),NOM,CFV(ICF),NSTEG)
      if(ITCF(2).ne.0)
     +CALL SUMCF(LIS,CC2(I),CC2(L),NOM,CFA(ICF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* Specific unit vectors:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      ICF        = NSTEG*(ITYP-1)+1
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
     X          ,CP1(ICF),CP2(ICF),NSTEG)
      if(ITCF(NTCF+4).ne.0)
     +CALL SUMPL(LIS,CC4(I),CC4(J),CC4(K),CC4(L),CC4(M),CC4(N),NOM
     X          ,CU1(ICF),CU2(ICF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* X Y Z projections:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      ICF        = NSTEG*(ITYP-1)+1
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
     X          ,CVX(ICF),CVY(ICF),CVZ(ICF),NSTEG)
      else
        if(ITCF(NTCF+1).ne.0)
     +  CALL SUMPS(LIS,CC1(I),CC1(J),CC1(K),CC1(L),CC1(M),CC1(N),NOM
     X          ,CVX(ICF),CVY(ICF),CVZ(ICF),NSTEG)
      end if
      if(LMOL2)then
       CALL SUMPS(LIS,CC6(I),CC6(J),CC6(K),CC6(L),CC6(M),CC6(N),NOM
     X          ,CAX(ICF),CAY(ICF),CAZ(ICF),NSTEG)
      else
        if(ITCF(NTCF+2).ne.0)
     +  CALL SUMPS(LIS,CC2(I),CC2(J),CC2(K),CC2(L),CC2(M),CC2(N),NOM
     X          ,CAX(ICF),CAY(ICF),CAZ(ICF),NSTEG)
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
      ICF        = (ITYP-1)*NSTEG+LIS+1
      NSP        = NSPEC(ITYP)
      NSB        = IXCF+1
      NSE        = IXCF+NSP*3
      DO I       = NSB,NSE
      L          = I+NOM
      if(ITCF(1).ne.0)CALL SUMCF(LEFT,CC1(I),CC1(L),NOM,CFV(ICF),NSTEG)
      if(ITCF(2).ne.0)CALL SUMCF(LEFT,CC2(I),CC2(L),NOM,CFA(ICF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* Specific unit vectors:
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
      ICF        = (ITYP-1)*NSTEG+LIS+1
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
     X          ,CP1(ICF),CP2(ICF),NSTEG)
	if(ITCF(NTCF+4).ne.0)
     +CALL SUMPL(LEFT,CC4(I),CC4(J),CC4(K),CC4(L),CC4(M),CC4(N),NOM
     X          ,CU1(ICF),CU2(ICF),NSTEG)
      END DO! I
      IXCF       = IXCF+NSP*3
      END DO! OF ITYP
*
* X Y Z projections
*
      IXCF       = 0
      DO ITYP    = 1,NTYPES
        ICF        = (ITYP-1)*NSTEG+LIS+1
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
     X          ,CVX(ICF),CVY(ICF),CVZ(ICF),NSTEG)
          else
            if(ITCF(NTCF+1).ne.0)
     +  CALL SUMPS(LEFT,CC1(I),CC1(J),CC1(K),CC1(L),CC1(M),CC1(N),NOM
     X          ,CVX(ICF),CVY(ICF),CVZ(ICF),NSTEG)
          end if
          if(LMOL2)then
       CALL SUMPS(LEFT,CC6(I),CC6(J),CC6(K),CC6(L),CC6(M),CC6(N),NOM
     X          ,CAX(ICF),CAY(ICF),CAZ(ICF),NSTEG)
          else
            if(ITCF(NTCF+2).ne.0)
     +  CALL SUMPS(LEFT,CC2(I),CC2(J),CC2(K),CC2(L),CC2(M),CC2(N),NOM
     X          ,CAX(ICF),CAY(ICF),CAZ(ICF),NSTEG)
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
      include "prcm.h"
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
