C    Tranal v.5.0
C    
C    Computation of Bilayer electron density
*
      program BILELDENS
      include "tranal.h"
      parameter (MAXDIV=1000)
      dimension DENE(MAXDIV),DENM(MAXDIV),DECH(MAXDIV),EPOT(MAXDIV)
      dimension EPOT2(MAXDIV)
      integer IATE(NS)
      character elm,fdens*64
      logical LCOM,LSYM
      namelist /ELDEN/ FDENS,IRT,NA,ZMAX,IATE,LSYM
*  FDENS - output file with results
*  IRT - type of lipid molecules
*  NA - resolution of density
*  ZMAX - maximum Z coordinate
*  IATE - for which sites compute density (1/0) (default = all)
*  LSYM  - forcible symmetrization of the distributions
      LSYM=.false.
      call SETTRAJ
      do I=1,NSITES
        IATE(I)=1
      end do
*  Enter analysis parameters
      read(*,ELDEN)
*  Prepare
      do I=1,NA
        DENE(I)=0.
        DENM(I)=0
        DECH(I)=0
      end do
      IEND=0
      BOXS=0.
      BOYS=0.
      BOZS=0.
      do while(IEND.eq.0)
        call READCONF(IEND)
        if(IEND.ne.0)go to 101
* compute Z-coordimate of membrane midplane
        call GETCOM
        NRFM=0
	ZRS=0.
	do IML=1,NSPEC(IRT)
	  IMOL=IML+IADDR(IRT)
	  if(NRFM.eq.0)then
	    ZRS=Z(IMOL)
	    ZRF=ZRS
	    NRFM=1
	  else
	    ZZ=Z(IMOL)
	    DZ=Z(IMOL)-ZRF
	    if(DZ.gt.HBOZL)ZZ=ZZ-BOZL
	    if(DZ.lt.-HBOZL)ZZ=ZZ+BOZL
	    ZRS=ZRS+ZZ
	    NRFM=NRFM+1
	    ZRF=ZRS/NRFM
	  end if
	  if(IPRINT.ge.8)write(*,'(i5,3f12.4)')IML,Z(IMOL),ZZ,ZRF
	end do
	if(ZRF.gt.HBOZL)ZRF=ZRF-BOZL
	if(ZRF.lt.-HBOZL)ZRF=ZRF+BOZL
        EU=0.
        ED=0.
*  Cycle over atoms 
        do I=1,NSTOT
	  IS=NSITE(I)
          if(IATE(IS).ne.0)then
	    DZ=SZ(I)-ZRF
	    if(DZ.gt.HBOZL)DZ=DZ-BOZL
	    if(DZ.lt.-HBOZL)DZ=DZ+BOZL
	    if(DZ.gt.HBOZL)DZ=DZ-BOZL
	    if(DZ.lt.-HBOZL)DZ=DZ+BOZL
**            write(*,'(i5,4f12.3)')I,ZRF,SZ(I),DZ,Z(NNUM(I))
	    NR=NA*((DZ+ZMAX)/2.)/ZMAX+1
            elm=NM(IS)(1:1)
*  Add here nuclear charges of atoms
	  if(elm.eq.'H')then
	    DED=1.-CHARGE(IS)
          else if(elm.eq.'C')then
	    DED=6.-CHARGE(IS)
	  else if(elm.eq.'N')then
	    DED=7.-CHARGE(IS)
	  else if(elm.eq.'O')then
	    DED=8.-CHARGE(IS)
	  else if(elm.eq.'P')then
	    DED=15.-CHARGE(IS)
	  else if(elm.eq.'S')then
	    DED=16.-CHARGE(IS)
	  else 
	    write(*,*)' Atom type ',elm,' not supported ',IS,NM(IS)
	    DED=0.
	  end if
*  nuclear density
	    DEM=MASS(IS)
	    if(NR.gt.0.and.NR.le.NA)then
              DENE(NR)=DENE(NR)+DED
              DENM(NR)=DENM(NR)+DEM
              DECH(NR)=DECH(NR)+CHARGE(IS)
            end if
            if(DZ.ge.0.)then
              EU=EU+CHARGE(IS)
            else
               ED=ED+CHARGE(IS)
            end if
          end if
        end do
        if(IPRINT.ge.6)write(*,'(a,f11.3,a,f9.4,a,2f12.5)')
     +'Time: ',FULTIM,' Zref = ',ZRF,' Q+/-: ',EU,ED
*  computing average BOX
        BOXS=BOXS+BOXL
        BOYS=BOYS+BOYL
        BOZS=BOZS+BOZL
      end do
 101  BOXS=BOXS/IAN
      BOYS=BOYS/IAN
      BOZS=BOZS/IAN
      write(*,*)IAN,' configurations analysed'
      write(*,*)' Average BOX:',BOXS,BOYS,BOZS
*  symmetrization
      EU=0.
      ED=0.
      if(LSYM)then
        do I=1,NA/2
          AEL=0.5*(DENE(I)+DENE(NA-I+1))
          DENE(I)=AEL
          DENE(NA-I+1)=AEL
          ANM=0.5*(DENM(I)+DENM(NA-I+1))
          DENM(I)=ANM
          DENM(NA-I+1)=ANM
          ACH=0.5*(DECH(I)+DECH(NA-I+1))
          DECH(I)=ACH
          DECH(NA-I+1)=ACH
        end do
      end if
*  potential calculation
      FAC=NA/(BOXS*BOYS*2.*ZMAX*IAN)
      FACM=1.d24/AVSNO                 ! to g/cm**3
      EFACT=ELECHG*1.d10*FAC/EPS0      ! 1d10:  Å -> m
      do I=1,NA
        POT=0.
        do J=1,NA
          RZ=abs((J-I)*2.)*ZMAX/NA
          POT=POT-RZ*DECH(J)
        end do
        if(I.le.NA/2)then
          EU=EU+DECH(I)
        else
          ED=ED+DECH(NA-I+1)
        end if
        EPOT(I)=EFACT*POT*ZMAX/NA
      end do
*  potential by double integration relative Z=0 (Tieleman&Berendsen 1996)
      NA2=NA/2
      do I=1,NA
        POT=0.
        if(I.gt.NA2)then
          do J=NA2+1,I
            QS=0.
            do K=NA2+1,J
              QS=QS+DECH(K)
            end do
            POT=POT+QS
          end do
        else
          do J=NA2,I,-1
            QS=0.
            do K=NA2,J,-1
              QS=QS+DECH(K)
            end do
            POT=POT+QS
          end do
        end if
        EPOT2(I)=-EFACT*POT*(2*ZMAX/NA)**2
      end do
*  output  (units - e/A^3)
      open(unit=8,file=FDENS,status='unknown')
      write(8,'(a)')
     +'#    Z        el.dens     mass dens  charge den   el.pot(V)',
     +' e.pot. DI'
      EU=EU/IAN
      ED=ED/IAN
      do I=1,NA
        RR=(I-0.5)*2.*ZMAX/NA-ZMAX
        DE=ECOR*(I-0.5*(NA-1))
        write(8,'(f12.3,6f12.6)')
     +  RR,DENE(I)*FAC,DENM(I)*FAC*FACM,DECH(I)*FAC,EPOT(I),EPOT2(I)
      end do
      if(LSYM)write(8,'(a)')'  Densities are symmetrized'
      write(*,*)' Q_up = ',EU,'  Q_down=',ED,' Q total ',EU+ED
      write(8,*)'# Q_up = ',EU,'  Q_down=',ED,' Q total ',EU+ED
      stop
      end
