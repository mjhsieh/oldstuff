*
*================ READMOL =============================================
*
*    Define molecule from DATABASE
	subroutine READMOL(SUMM,NRBB,NRAA,NRTT,NRII,NX,NTYP,PATHDB,NAMN)
	include "mdee.h"
*
        DIMENSION CM(3)
	character PATHDB*62,namn*12,STR*128,FIL*128,TAKESTR*128
        FIL=' '
	NTYP        = NTYP+1
      LD=LENG(PATHDB,62)
      FIL(1:LD)=PATHDB(1:LD)
        FIL(LD+1:LD+1)='/'
	LN=LENG(NAMN,12)
	FIL(LD+2:LD+LN+1)=NAMN(1:LN)
	FIL(LD+LN+2:LD+LN+6)='.mmol'
	IE=0
	open(unit=12,file=fil,status='old',err=999)
	if(IPRINT.ge.6)then
	PRINT "(80('-'))"
      PRINT *
      PRINT * ,'*** MODEL ONE MOLECULE OF TYPE: ',NAMN,'  No.',NTYP
	end if
	STR=TAKESTR(12,IE)
	read(STR,*)NSTS
	if(IPRINT.ge.6)PRINT * ,'*** NR OF SITES:                ',NSTS
	NSITS(NTYP)=NSTS
	NSBEG=NX+1
	do I=1,NSTS
	  STR=TAKESTR(12,IE)
	  IX=NX+I
	  NM(IX)=STR(1:4)
	  read(STR(5:80),*,err=901)(R(IX,J),J=1,3),
     +  MASS(IX),CHARGE(IX),SIGMA(IX),EPSIL(IX)
          if(IPRINT.ge.6)
     +      write(*,'(a5,a4,a3,f9.4,a3,f6.3,a5,f9.4,a5,f9.4)')
     +      'atom ',NM(IX),' M=',MASS(IX),' Q=',CHARGE(IX),' sig=',
     +      SIGMA(IX),' eps=',EPSIL(IX)
	end do
	NX=NX+NSTS
	NSEND=NX
* COM coordinates 
	SUMM       = 0.D0
      DO IS      = NSBEG,NSEND
        MASS(IS)   = MASS(IS)
        SUMM       = SUMM+MASS(IS)
      END DO! OF IS
*
      CM(1)      = 0.D0
      CM(2)      = 0.D0
      CM(3)      = 0.D0
      DO IS      = NSBEG,NSEND
        CM(1)      = CM(1)+R(IS,1)*MASS(IS)
        CM(2)      = CM(2)+R(IS,2)*MASS(IS)
        CM(3)      = CM(3)+R(IS,3)*MASS(IS)
      END DO! OF IS
*
      CM(1)      = CM(1)/SUMM
      CM(2)      = CM(2)/SUMM
      CM(3)      = CM(3)/SUMM
*
      DO K       = 1,3
        DO IS      = NSBEG,NSEND
          R(IS,K)    = R(IS,K)-CM(K)
        END DO! OF IS
      END DO! OF I
*  Output of the reference
	STR=TAKESTR(12,IE)
	read(STR,*)NSR
	do I=1,NSR
	  read(12,'(a80)')STR
	  if(IPRINT.ge.8)write(*,'(a80)')STR
	end do
	if(IPRINT.ge.8)then
        PRINT *,'*** Molecular GEOMETRY (C.O.M. IN ORIGIN): '
        DO IS      = NSBEG,NSEND
          PRINT "('*** ',A4,3X,3(1X,F7.3))",
     +    NM(IS),R(IS,1),R(IS,2),R(IS,3)
        END DO! OF IS
      end if
*
	if(IPRINT.ge.8)then
        PRINT *
        PRINT *,'*** INTER ATOMIC DISTANCES '
        do IS=NSBEG,NSEND
	    do JS=NSBEG+1,NSEND
	      RR2=0.
	      do K=1,3
	        RR2=RR2+(R(IS,K)-R(JS,K))**2
	      end do
 	      RR=sqrt(RR2)
            PRINT "('*** ',A4,'-  ',A4,'  --> ',F7.3)",
     +      NM(IS  ),NM(JS),RR
	    end do
	  end do
	end if
*  Bounds
	STR=TAKESTR(12,IE)
	read(STR,*)NBD
	NRB(NTYP)=NBD
	do K=NRBB+1,NRBB+NBD
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=702)ID(K),IB(K),JB(K),RB(K),FB(K),DB(K),RO(K)
	  go to 703
 702	  read(STR,*,err=902)ID(K),IB(K),JB(K),RB(K),FB(K)
	  ID(K)=0
 703	  if(RB(K).le.0.d0)RB(K)=sqrt((R(IB(K),1)-R(JB(K),1))**2+
     +  (R(IB(K),2)-R(JB(K),2))**2+(R(IB(K),3)-R(JB(K),3))**2)
	end do
	if(IPRINT.ge.8)then
        PRINT *
        PRINT "('*** PARAMETERS FOR ',I3,'  BONDS')",NBD
        PRINT *,'*** HARMONIC POTENTIALS:'
        IDTOT      = 0
        DO M       = NRBB+1,NRBB+NBD
          if(ID(M).ne.0)IDTOT=IDTOT+1
          PRINT"('*** BOND: ',I3,'  - ',I3,'  R(eq): ',F7.3,'  FORCE ',
     X    'CONSTANT:',F9.3)",IB(M),JB(M),RB(M),FB(M)
        END DO! OF M
*
        IF(IDTOT.GT.0) THEN
          PRINT *
          PRINT *,'*** MORSE POTENTIALS:'
          DO M       = NRBB+1,NRBB+NRB(NTYP)
            if(ID(M).ne.0)
     +      PRINT "('*** BOND: ',I3,'  - ',I3,'  RO: ',F7.3,'  DISS  ',
     X      'ENERGY  :',F9.3)",IB(M),JB(M),RO(M),DB(M)
          END DO! OF M
        END IF
	end if
*  Angles
	STR=TAKESTR(12,IE)
	read(STR,*)NAN
	NRA(NTYP)=NAN
	do K=NRAA+1,NRAA+NAN
	  STR=TAKESTR(12,IE)
	  read(STR,*,err=903)IA(K),JA(K),KA(K),RA(K),FA(K)
	  if(RA(K).le.0)RA(K)=ANGLG(R,IA(K),JA(K),KA(K),NS)*TODGR
	end do
	if(IPRINT.ge.8)then
	  PRINT *
        PRINT *,'*** COVALENT ANGLES '
	  do K=NRAA+1,NRAA+NAN
          PRINT "('***',3(1X,A4),'  ->',F9.3,'  Force c.',f9.3)",
     +    NM(IA(K)),NM(JA(K)),NM(KA(K)),RA(K),FA(K)        
	  end do
	end if
*  Dihedrals
      STR=TAKESTR(12,IE)
	read(STR,*)NTT
	NRT(NTYP)=NTT
	if(NTT.gt.0)then
*	  STR=TAKESTR(12,IE)
*	  read(STR,*)NTT0
	  NTT0=ITORS(NTYP)
	  IF(NTT0.EQ.0) THEN
	  do K=NRTT+1,NRTT+NTT
	    STR=TAKESTR(12,IE)
	    read(STR,*,err=904)
     +	    IT(K),JT(K),KT(K),LT(K),RT(K),FT(K),NMUL(K)
	    FT1(K)=0.d0
	    FT2(K)=0.d0
	    FT3(K)=0.d0
	  end do
	  else
	  do K=NRTT+1,NRTT+NTT
	    STR=TAKESTR(12,IE)
	    read(STR,*,err=904)
     +	    IT(K),JT(K),KT(K),LT(K),FT1(K),FT2(K),FT3(K)
	    NMUL(K)=0
	    RT(K)=0.d0
	    FT(K)=0.d0
	  end do
	  end if
	  if(IPRINT.ge.8)then
	    PRINT *
            PRINT *,'*** DIHEDRALS '
	    if(NTT0.eq.0)then
	      write(*,*)'  1+cos(nf) type:'
	      do K=NRTT+1,NRTT+NTT
                PRINT"(4I6,4(1X,A4),' ->',I2,F8.2,'  F.c.',f9.3)"
     +          ,   IT(K),JT(K),KT(K),LT(K)
     +          ,NM(IT(K)),NM(JT(K)),NM(KT(K)),NM(LT(K)),
     +          NMUL(K),RT(K),FT(K)        
	      end do
	    end if
	    if(NTT0.gt.0)then
	      write(*,*)'  other type:'
	      do K=NRTT+1,NRTT+NTT
                PRINT"('***',4(1X,A4),'  ->','  Force c.',3f9.3)",
     +          NM(IT(K)),NM(JT(K)),NM(KT(K)),NM(LT(K)),
     +          FT1(K),FT2(K),FT3(K)        
	      end do
	    end if
	  end if
        end if	
*
	NRBB         = NRBB+NRB(NTYP)
      NRAA         = NRAA+NRA(NTYP)
      NRTT         = NRTT+NRT(NTYP)
      IADB(NTYP)   = NRBB-NRB(NTYP)+1
      IADA(NTYP)   = NRAA-NRA(NTYP)+1
      IADT(NTYP)   = NRTT-NRT(NTYP)+1
*     
	close(12)
      return
 901  write(*,*)' Error in the list of atoms, line ',I
	stop
 902  write(*,*)' Error in the list of bonds, line ',K
	stop
 903  write(*,*)' Error in the list of angless, line ',K
	stop
 904  write(*,*)' Error in the list of dihedralss, line ',K
	stop
 999  write(*,*)'!!! Molecule ',NAMN,' not found in the data base'
      write(*,*)'    File not found: ',FIL
	stop
	end 
*
*=============== ANGLG =============================================
*
	FUNCTION ANGLG(R,I,J,K,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(N,3)
	RIJ=0.
	RJK=0.
	SC=0.
	DO M=1,3
	  RIJ=RIJ+(R(I,M)-R(J,M))**2
        RJK=RJK+(R(J,M)-R(K,M))**2
        SC =SC +(R(I,M)-R(J,M))*(R(K,M)-R(J,M))
	end do
      ANGLG = DACOS(SC/sqrt(RIJ*RJK))
      RETURN
      END
*================================================================
*
*    3. Read and check a line from input file
*    ----------------------------------------
*
C    This subroutine reads a line from Fortran input file 
C    defined by "KAN" (i.e. opened with parameter unit=KAN)
C    It returns the line if it is not started with # (commentaries); 
C    otherwise it reads and check the next line.
C    Parameter IE:
C    If input value of IE is NOT 99, the subroutine performs normal action
C    (see above) and return IE=0 (normal exit) or IE=-1 (read error)
C    If input value of IE=99, the subroutine is used for diagnostics:
C    it reports the file name and the line number where input error occur.
C    
C
C============== TAKESTR =========================================
C
*   3.1 Definitions
      character*128 function TAKESTR(KAN,IE)
      character*128 AUX, fn*32 
      integer LCOUNT(256)   
      data LCOUNT /256*0/
      save AUX,LCOUNT
*   3.2 Normal action           
      if(IE.eq.99)go to 10
C  LCOUNT variable count line number (including commentaries) 
 1	LCOUNT(KAN)=LCOUNT(KAN)+1
 	read(KAN,'(a128)',err=10,end=20)AUX
	if(AUX(1:1).eq.'#')go to 1
	TAKESTR=AUX 
	IE=0
	return
*   3.3 Error diagnostic
 10	write(*,*)'!!! error in input file '   
	if(KAN.eq.5)then
	  write(*,*)' Standard input ' 
	else
	  inquire(unit=KAN,name=fn) 
	  write(*,*)' File ',fn
	end if
	write(*,*)' in line ',LCOUNT(KAN)
	write(*,*)AUX
	stop
*    3.4 End of file case (is not signaled if input IE < 0 )
 20	if(IE.ge.0)then
 	  write(*,*)'!!! end of file reached' 
	  if(KAN.eq.5)then
	    write(*,*)' Standard input ' 
	  else
	    inquire(unit=KAN,name=fn) 
	    write(*,*)' File ',fn
	  end if
	  stop
	end if  
	TAKESTR='  ' 
	IE=-1
	return
	end 
*
*================ LENG ===============================================
*
	function LENG(STR,LS)
	character*1 STR(LS)
	IL=0 
	do I=1,LS
	  ISTR=ichar(STR(I))
	  if(ISTR.le.15)then
	    LENG=I-1
	    return
	  end if
        if(STR(I).ne.' ')IL=I
	end do 
	LENG=IL
	return
      end
