*=============== SHEJK ===========================================
*
      SUBROUTINE SHEJK(ITYP)
*
      include "mdee.h"
*
      LOGICAL MOVING(NTOT),MOVED (NTOT),DONE
      PARAMETER (MAXIT=2000)
*
      VIR1	= 0.
      NRBS      = NRB(ITYP)
      IF(NRBS.LE.0)  RETURN
      TOL2      = TOLER**2
*
      NBBEG     = IADB(ITYP)
      NBEND     = NBBEG+NRBS-1
*
      NSP       = NSPEC (ITYP)
      NSS       = NSITS (ITYP)
      ISHF      = ISADDR(ITYP)
      ISBEG     = ISHF+1
      ISEND     = ISHF+NSP*NSS
*>    
      DO K      = 1,NSP
*>
        NTS       = (K-1)*NSS
        DO M      = NBBEG,NBEND
          IBB       = IB(M)+ISHF
          JBB	      = JB(M)+ISHF
	    LI        = NTS+IBB
	    LJ	      = NTS+JBB
          MOVING(LI) = .FALSE.
          MOVED (LI) = .TRUE.
          MOVING(LJ) = .FALSE.
          MOVED (LJ) = .TRUE.
        END DO! OF M
*
* Begin iterations:
*     
        ITER      = 0
        DONE      = .FALSE.
*<
    1   CONTINUE
*<
        IF(ITER.GT.MAXIT) GO TO 1000
        IF(.NOT.DONE) THEN
*     
        DONE      = .TRUE.
        DO M      = NBBEG,NBEND
          IBB       = IB(M)+ISHF
          JBB       = JB(M)+ISHF
          I         = NTS+IBB
          J         = NTS+JBB
          IF(MOVED(I).OR.MOVED(J)) THEN
            AMAS      = 1.D0/MASSDI(I)
            BMAS      = 1.D0/MASSDI(J)
            BNDSQ     = RB(M)**2
*Unconstrained positions at step n+1
            BXN       = SX(I)-SX(J)
            BYN       = SY(I)-SY(J)
            BZN       = SZ(I)-SZ(J)
            BSQ       = BXN**2+BYN**2+BZN**2
            DIFFSQ    = BNDSQ-BSQ
      if(IPRINT.ge.10)then
      write(*,'(3i5,f9.6,4(2x,l1),i5)')
     +K,I,J,sqrt(BSQ),moved(I),moved(j),moving(I),moving(j),ITER
      end if
*
            IF(DABS(DIFFSQ).GT.(BNDSQ*TOL2))THEN
*
Constrained positions from step n
              DXN       = OX(I)-OX(J)
              DYN       = OY(I)-OY(J)
              DZN       = OZ(I)-OZ(J)
              VDOT      = BXN*DXN+BYN*DYN+BZN*DZN
              IF(VDOT.LT.BNDSQ*BDTOL) THEN
      WRITE(6,66) I,SX(I),SY(I),SZ(I),OX(I),OY(I),OZ(I)
      WRITE(6,66) J,SX(J),SY(J),SZ(J),OX(J),OY(J),OZ(J)
*
      PRINT *,'NEW DISTANCE:',DSQRT((SX(I)-SX(J))**2
     X                             +(SY(I)-SY(J))**2
     X                             +(SZ(I)-SZ(J))**2)
      PRINT *,'OLD DISTANCE:',DSQRT((OX(I)-OX(J))**2
     X                             +(OY(I)-OY(J))**2
     X                             +(OZ(I)-OZ(J))**2)
*
   66 FORMAT(1X,I5,3(1X,F10.6),3X,3(1X,F10.6))
      call XMOLDUMP
      STOP 'CONSTRAINED FAILURE'
*      print *,' ... constraint failure!'
              END IF
              RMA       = MASSDI(I)
              RMB       = MASSDI(J)
              GAB       = DIFFSQ/(2.D0*(RMA+RMB)*VDOT)
              DXN       =   DXN*GAB
              DYN       =   DYN*GAB
              DZN       =   DZN*GAB
              SX(I)     = SX(I)+RMA*DXN
              SY(I)     = SY(I)+RMA*DYN
              SZ(I)     = SZ(I)+RMA*DZN
              SX(J)     = SX(J)-RMB*DXN
              SY(J)     = SY(J)-RMB*DYN
              SZ(J)     = SZ(J)-RMB*DZN
              MOVING(I) = .TRUE.
              MOVING(J) = .TRUE.
              DONE      = .FALSE.
            END IF
          END IF
        END DO! OF M - all bonds
* 
*  Obs!    These are separate loops ( Reccurence possible )
*
        DO M      = NBBEG,NBEND
          IBB       = IB(M)+ISHF
	  JBB	      = JB(M)+ISHF
          I         = NTS+IBB
	  J	      = NTS+JBB
	  MOVED(I)  = .false.
          MOVED(J)  = .false.
        END DO! OF M
        DO M      = NBBEG,NBEND
          IBB       = IB(M)+ISHF
	  JBB	      = JB(M)+ISHF
          I         = NTS+IBB
	  J	      = NTS+JBB
          if(MOVING(I))MOVED(I)=.true.
          MOVING(I) = .FALSE.
          if(MOVING(J))MOVED(J)=.true.
          MOVING(J) = .FALSE.
        END DO! OF M
*
        ITER      = ITER + 1
*
        GO TO 1
*
        END IF
*
 1000   CONTINUE
*
        IF(ITER.GE.MAXIT) THEN
      PRINT "(/1X,'??? NO CONVERGENCE AFTER MAXIT INTERATIONS')"
      PRINT "(1X,'??? BOND BETWEEN ATOMS',I4,'  AND',I4,' '
     X      ,/'  EQ. BOND = ',F12.6,' A --> CALC.= ',F12.6,' A'/)"
     X      ,I,J,DSQRT(BNDSQ),DSQRT(BSQ)
        END IF
*<
      END DO! OF K - all molecules
C   velocity correction
      N         = ISADDR(ITYP)
      NSPNS     = NSPEC(ITYP)*NSITS(ITYP)
      DO I      = 1,NSPNS
        N         = N+1
        DIX       = (SX(N)-FX(N))/(MASSDI(N)*TSTEP)
        DIY       = (SY(N)-FY(N))/(MASSDI(N)*TSTEP)
        DIZ       = (SZ(N)-FZ(N))/(MASSDI(N)*TSTEP)
        VX(N)     = VX(N)+DIX
        VY(N)     = VY(N)+DIY
        VZ(N)     = VZ(N)+DIZ
	  VIR1	= VIR1+(OX(N)*DIX+OY(N)*DIY+OZ(N)*DIZ)/TSTEP
      END DO! OF N
      VIR	= VIR+VIR1
      VIRA	= VIRA+VIR1
      RETURN
      END
*
*=============== BNDLST ===========================================
*
      SUBROUTINE BNDLST
*
      include "mdee.h" 
	parameter (RMX2=5.**2)
*                                 
	do I=1,NSTOT
	  NNBB(I)=0
	  do IS=1,NBDMAX
           INBB(I,IS)=0
	  end do
	end do       
*	IBEG=TASKID+1 
*	NRN=NUMTASK  
*	if(NRN.gt.NTYPES)NRN=NTYPES
	do ITYP=1,NTYPES
	  NRNOD(ITYP)=mod(ITYP-1,NUMTASK)
	end do  
*	if(IBEG.gt.NTYPES)go to 101
	do ITYP=1,NTYPES
*           
        NNMAX     = -99999
        NSP       = NSPEC(ITYP)     ! num of mol of this type
        NSS       = NSITS(ITYP)     ! num of sites on a mol of this type
        NRBS      = NRB(ITYP)       ! num of bonds 
        NRAS      = NRA (ITYP)      ! num of angles 
        NRTS      = NRT(ITYP)       ! num of torsions
*
        NBBEG     = IADB(ITYP)
        NBEND     = NBBEG+NRBS-1
        NABEG     = IADA(ITYP)
        NAEND     = NABEG+NRAS-1
        NTBEG     = IADT(ITYP)
        NTEND     = NTBEG+NRTS-1                       
*             
*
        ISBEG     = ISADR(ITYP)+1
        ISEND     = ISADR(ITYP+1)
        ISHF      = ISADR(ITYP)  ! shift of global site num. relatevely local
	ISHA	    = ISADDR(ITYP)-ISHF
*
*  Double sum over sites 
*  
	 DO IS=ISBEG,ISEND-1
	    DO JS=IS+1,ISEND   
*  IS, JS - site number
*  Bonds    
	      if(.not.L15NB(ITYP))go to 100
	      I=ISHA+IS
	      J=ISHA+JS                
	      DX=SX(I)-SX(J)
	      DY=SY(I)-SY(J)
	      DZ=SZ(I)-SZ(J)
	      call PBC(DX,DY,DZ)                
*   atoms are too far
              RRR=DX**2+DY**2+DZ**2
              if(RRR.gt.RMX2)go to 150
	      do IIS=NBBEG,NBEND
		II=IB(IIS)+ISHF
		JJ=JB(IIS)+ISHF
		if(II.eq.IS.and.JJ.eq.JS)go to 100
		if(II.eq.JS.and.JJ.eq.IS)go to 100
	      end do
*  Angles
	      do IIS=NABEG,NAEND
		II=IA(IIS)+ISHF
		KK=KA(IIS)+ISHF
		if(II.eq.IS.and.KK.eq.JS)go to 100
		if(II.eq.JS.and.KK.eq.IS)go to 100
	      end do          
*  Torsions      
	      do IIS=NTBEG,NTEND
		II=IT(IIS)+ISHF
		LL=LT(IIS)+ISHF 
		ITR=ITORS(ITYP)
		if((II.eq.IS.and.LL.eq.JS).or.(II.eq.JS.and.LL.eq.IS))
     +	then
		  if(ITR.eq.-1.or..not.L14NB(ITYP))then
		    go to 100
		  else
	  	    go to 90
	 	  end if
		end if
	      end do
*  Pass all checks - non-bound 
            go to 150
*  These atoms are 1-4 bound (Scaling LJ parameters may be used)
 90         do IIM=0,NSP-1
	        ISHIFT=ISADDR(ITYP)-ISADR(ITYP)+IIM*NSS 
		ISA=IS+ISHIFT 
		JSA=JS+ISHIFT 
	 	NNBB(ISA)=NNBB(ISA)+1      ! number of bound neigbours for this atom
		INBB(ISA,NNBB(ISA))=-JSA   ! bound neigbours ( - for 1-4 interactions)
		if(NNBB(ISA).gt.NBDMAX)then
		  write(*,*)
     +' list of bound atoms exceeded for atom ',ISA,' ',NM(ISA)
		  stop
		end if 
	 	NNBB(JSA)=NNBB(JSA)+1
		INBB(JSA,NNBB(JSA))=-ISA  
		if(NNBB(JSA).gt.NBDMAX)then
		  write(*,*)
     +' list of bound atoms exceeded for atom ',JSA,' ',NM(JSA)
		  stop
		end if 
	      end do 
	      go to 150
*   Atoms are 1-3 or 1-2 bounded. 
 100        do IIM=0,NSP-1
	        ISHIFT=ISADDR(ITYP)-ISADR(ITYP)+IIM*NSS
		ISA=IS+ISHIFT 
		JSA=JS+ISHIFT 
	 	NNBB(ISA)=NNBB(ISA)+1
		INBB(ISA,NNBB(ISA))=JSA
		if(NNBB(ISA).gt.NBDMAX)then
		  write(*,*)
     +' list of bound atoms exceeded for atom ',ISA,' ',NM(ISA)
		  stop
		end if 
	 	NNBB(JSA)=NNBB(JSA)+1
		INBB(JSA,NNBB(JSA))=ISA  
		if(NNBB(JSA).gt.NBDMAX)then
		  write(*,*)
     +' list of bound atoms exceeded for atom ',JSA,' ',NM(JSA)
		  stop
		end if 
	      end do
 150        continue
	    end do
	  end do
*----------------------------
 	  write(*,*)
     +  ' node ',TASKID,'  list of bound atoms for ',NAME(ITYP)
	end do   ! of ITYP      
 101	continue
CM call  GATHER(NRNOD)
CM
	if(IPRINT.ge.8)then
	  write(*,*)' list of bound atoms: '
	  do I=1,NSTOT
	    write(*,'(i5,i3,2x,40i5)')I,NNBB(I),(INBB(I,J),J=1,NNBB(I))
	  end do 
	end if
*----------------------------------------
      RETURN
      END   
