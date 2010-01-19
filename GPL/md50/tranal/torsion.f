C    Tranal - package for rajectory analysis, MDynaMix v.5.0
C    Compute distribution over torsion angles
C    
      program TORSION
      include "tranal.h"
      parameter (NTOR=100,NAAM=1000,NT=1000)
      integer MANGL(NTOR),IADT(NTPS),IT(NT),JT(NT),KT(NT),LT(NT)
      integer IGS(NTOR),ITRAN(NTOR),NPOP(NTOR,NAAM)
      character FTOR*64
      call SETTRAJ
*  Input of data - after "TRAJ" section, 
*  positional input with optional commentaries
*  
*  1) File name of the output file
      IE=0
      STR=TAKESTR(5,IE)
      read(STR,*)FTOR
*  2) number of different torsion types and number of bins in [-180:180]
      STR=TAKESTR(5,IE)
      read(STR,*)NAG,NAA
      if(NAG.gt.NTOR)stop ' increase NTT'
      if(NAA.gt.NAAM)stop ' increase NAAM'
*   3) how many torsions of each type
      STR=TAKESTR(5,IE)
      read(STR,*)(MANGL(I),I=1,NAG)    ! torsions of each type
      write(*,*)(I,MANGL(I),I=1,NAG)
      IADT(1)=0
      do I=1,NAG
	IADT(I+1)=IADT(I)+MANGL(I)
      end do
      NADT=IADT(NAG+1)                   ! total  
      if(NADT.gt.NT)stop 'increase NT'
*   4) List of torsions - 4 site numbers for each torsion
*      1st group, then 2-nd ..
      do I=1,NADT                   
	STR=TAKESTR(5,IE)
  	READ(STR,*)IT(I),JT(I),KT(I),LT(I)    ! which
      end do
      do I=1,NAG
	MANGL(I)=MANGL(I)*NSPEC(ITS(IT(IADT(I)+1)))
      end do
      if(IPRINT.ge.6)then 
	write(*,*)' Compute distributions over torsion angles:'
	do I=1,NAG
        write(*,*)I,' type'  
	do J=IADT(I)+1,IADT(I+1)  
	  write(*,*)NM(NSITE(IT(J))),NM(NSITE(JT(J))),
     + NM(NSITE(KT(J))),NM(NSITE(LT(J)))
	  end do
	end do               
      end if
      do I=1,NAG
	IGS(I)=0
	ITRAN(I)=0
	do J=1,NAA
	  NPOP(I,J)=0 
	end do
      end do
*  start analysis
      IEND=0
      do while(IEND.eq.0)
        call READCONF(IEND)
	do ITYP=1,NAG     !  type of torsion
	  do N=IADT(ITYP)+1,IADT(ITYP+1)
            I         = IT(N)  
            J         = JT(N)
            K         = KT(N)
            L         = LT(N)
	    MTYP = ITS(I)                 !  type of molecule
	    do IMOL=1,NSPEC(MTYP)         !  over molecules of this type
	      IST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+I
	      JST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+J
	      KST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+K
	      LST = ISADDR(MTYP)+(IMOL-1)*NSITS(MTYP)-ISADR(MTYP)+L
        ax        = sx(jst)-sx(ist)
        ay        = sy(jst)-sy(ist)
        az        = sz(jst)-sz(ist)
        bx        = sx(kst)-sx(jst)
        by        = sy(kst)-sy(jst)
        bz        = sz(kst)-sz(jst)
        cx        = sx(lst)-sx(kst)
        cy        = sy(lst)-sy(kst)
        cz        = sz(lst)-sz(kst)
*
        call PBC(AX,AY,AZ)
        call PBC(BX,BY,BZ)
        call PBC(CX,CY,CZ)
*
        ab        = ax*bx + ay*by + az*bz
        bc        = bx*cx + by*cy + bz*cz
        ac        = ax*cx + ay*cy + az*cz
        at        = ax*ax + ay*ay + az*az
        bt        = bx*bx + by*by + bz*bz
        ct        = cx*cx + cy*cy + cz*cz

        axb       = (at*bt)-(ab*ab)
        bxc       = (bt*ct)-(bc*bc)
*
        fnum      = (ab*bc)-(ac*bt)
        den       = axb*bxc
        if(den.gt.0.0d0) then 
*
          den       = dsqrt(den)         
Cosine of angle:
          co        = fnum/den
          CO        = DMIN1(CO, 1.D0)
          CO        = DMAX1(CO,-1.D0)
* Sign of angle:
      signum    = ax*(by*cz-cy*bz)+ay*(bz*cx-cz*bx)+az*(bx*cy-cx*by)
* Value of angle:
          arg       = dsign(dacos(co),signum)
	  else
	    arg=0.d0
	  end if    
	  NARG=0.5*(arg+PI)*NAA/PI+1
	  if(abs(arg*180/PI).le.120)then
	    IGS(ITYP)=IGS(ITYP)+1
	  else
	    ITRAN(ITYP)=ITRAN(ITYP)+1
	  end if
	  NPOP(ITYP,NARG)=NPOP(ITYP,NARG)+1 
	  if(IPRINT.ge.8)write(*,'(6i5,f10.3)')
     +ITYP,IMOL,IST,JST,KST,LST,ARG*180/PI
	  end do
	  end do
	end do
      end do
      write(*,*)IAN,' configurations analysed'
      write(*,*)' Results are written in file ',FTOR
	  open(unit=17,file=ftor,status='unknown')
	  write(17,'(a)')'# distribution over torsion angles'
	  write(17,'(a1,10x,50(a4,5x))')'#',
     +(NM(NSITE(IT(IADT(I)+1))),I=1,NAG)
	  write(17,'(a1,10x,50(a4,5x))')'#',
     +(NM(NSITE(JT(IADT(I)+1))),I=1,NAG)
	  write(17,'(a1,10x,50(a4,5x))')'#',
     +(NM(NSITE(KT(IADT(I)+1))),I=1,NAG)
	  write(17,'(a1,10x,50(a4,5x))')'#',
     +(NM(NSITE(LT(IADT(I)+1))),I=1,NAG)
	  do M=1,NAA
	    ANG=-180+(M-0.5)*360/NAA 
	    FACT=NAA*1.d0/IAN
	write(17,'(f7.1,50f9.4)')ANG,
     + (NPOP(I,M)*FACT/MANGL(I),I=1,NAG)
	  end do
	  write(17,'(a,50f9.6)')'# Gaushe ',
     +(IGS(I)*1./(IGS(I)+ITRAN(I)),I=1,NAG)
	  write(17,'(a,50f9.6)')'# Trans  ',
     +(ITRAN(I)*1./(IGS(I)+ITRAN(I)),I=1,NAG)
	  close(17)
      stop
      end
