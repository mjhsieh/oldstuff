C    Tranal - package for rajectory analysis, MDynaMix v.5.0
*
*  Compute order parameter / dipole couplings
*
*
      program ORDER
      include "tranal.h"
      parameter (NEQHM=100,IODM=500)
      integer NC(NEQHM),NH(NEQHM),IOR(NEQHM),ORF1(NEQHM),ORF2(NEQHM)
      real*8 CSHL(NEQHM),DCP(NEQHM),SORD(NEQHM),ORDR(IODM,NEQHM)
      logical LREFL
      character FORD*64
      call SETTRAJ
*  Input of data
*  1) Name of output file
      STR=TAKESTR(5,IE)
      read(STR,*)FORD
*  2) Type number of lipid molecules
      STR=TAKESTR(5,IE)
      read(STR,*)IRT
      if(IRT.le.0)then
        IRT=1
        LREFL=.false.
      else
        LREFL=.true.
      end if
      if(IRT.gt.NTYPES)then
        write(*,*)
     +' parameter IRT is greater than the number of molecule types:',
     + IRT
        stop
      end if
*  3) number of order parameters
      STR=TAKESTR(5,IE)
      read(STR,*)NEQH
      if(NEQH.le.0)then
        write(*,*)' Unacceptable number of CH bonds in order parameter',
     &' computations: ',NEQH
        stop
      end if
      if(NEQH.gt.NEQHM)then
        write(*,*)' Increase NEQHM'
        stop
      end if
      J=1
*   4) Sites of fisrt and second atoms as well as 0 or 1
 19   STR=TAKESTR(5,IE)
        read(STR,*,end=20)NC(J),NH(J),IOR(J)
        if(IOR(J).gt.0)then
          if(J.eq.NEQH)then
            write(*,*)' Warning: IOR changed to 0 for the last pair ' 
            write(*,*)' in order parameters calculations'
            IOR(J)=0
            go to 15
          end if
          read(STR,*,end=20)NC(J),NH(J),IOR(J),ORF1(J),ORF2(J)
          J=J+1
          STR=TAKESTR(5,IE)
*   4a) If third index in previous line was "1":
          read(STR,*,err=20)NC(J),NH(J)
          ORF1(J)=ORF1(J-1)
          ORF2(J)=ORF2(J-1)
          IOR(J)=0
        else
          ORF1(J)=0
          ORF2(J)=0
          IOR(J)=0
        end if
        J=J+1
        if(J.le.NEQH)go to 19
 15   write(*,*)' Order parameter will be calculated for:'
      do J=1,NEQH
        write(*,'(a4,a1,a4,a,i5,a,i5)')
     +  NM(NC(J)),'-',NM(NH(J)),', or ',NC(J),' - ',NH(J)
        if(IOR(J).eq.1)write(*,*)' in pair with '
        CSHL(J)=0.              ! length
        DCP(J)=0.               ! dipole coupling parameter
	SORD(J)=0               ! order parameter 
      end do
*   5) Last line - of not zero, compute distribution of Cos(theta) angles
      STR=TAKESTR(5,IE)
      read(STR,*)IOD         ! resolution of the distributions
      if(IOD.lt.0)IOD=0
      if(IOD.gt.IODM)IOD=IODM
      close(29)
      if(IOD.gt.0)then
        do J=1,NEQH
          do I=1,IOD
            ORDR(I,J)=0.
          end do
        end do
      end if
* Begin analysis
      IEND=0
      do while(IEND.eq.0)
        call READCONF(IEND)
        J=0
 5      J=J+1
        if(J.gt.NEQH)go to 50
*  compute reference Z-coord 
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
	end do
	if(IPRINT.ge.6)write(*,*)' Zref = ',ZRF
        ITYP=ITS(NC(J))
        NML=NSPEC(ITYP)
        if(IOR(J).eq.1)then
* If we want to distinguish enantiomers (flag IOR=1)
          J1=J+1
          do N=1,NML
	    IC = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+NC(J)-ISADR(ITYP)
	    IH = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+NH(J)-ISADR(ITYP)
            DX1 = SX(IH)-SX(IC)
            DY1 = SY(IH)-SY(IC)
            DZ1 = SZ(IH)-SZ(IC)
            call PBC(DX1,DY1,DZ1)
            RR1 = sqrt(DX1**2+DY1**2+DZ1**2)
	    IC2 = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+NC(J1)-ISADR(ITYP)
	    IH2 = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+NH(J1)-ISADR(ITYP)
            DX2 = SX(IH2)-SX(IC2)
            DY2 = SY(IH2)-SY(IC2)
            DZ2 = SZ(IH2)-SZ(IC2)
            call PBC(DX2,DY2,DZ2)
            RR2 = sqrt(DX2**2+DY2**2+DZ2**2)
            CORD1 = DZ1/RR1                  ! Cos(theta)
            CORD2 = DZ2/RR2
            ORD1 = 0.5*(3.*CORD1**2-1.)
            ORD2 = 0.5*(3.*CORD2**2-1.)
            DCP1 = ORD1/RR1**3
            DCP2 = ORD2/RR2**3
*  reference vector
	    I3 = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+ORF1(J)-ISADR(ITYP)
	    I4 = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+ORF2(J)-ISADR(ITYP)
            DX3 = SX(I3)-SX(I4)
            DY3 = SY(I3)-SY(I4)
            DZ3 = SZ(I3)-SZ(I4)
            call PBC(DX3,DY3,DZ3)
            RR3 = sqrt(DX3**2+DY3**2+DZ3**2)
*  Checking scalar product
            SC1 = (DX1*DX3+DY1*DY3+DZ1*DZ3)/(RR1*RR3)
            SC2 = (DX2*DX3+DY2*DY3+DZ2*DZ3)/(RR2*RR3)
            if(SC1.le.SC2)then     ! the same order   
              CSHL(J)=CSHL(J)+RR1
              SORD(J)=SORD(J)+ORD1
              DCP(J)=DCP(J)+DCP1
              CSHL(J1)=CSHL(J1)+RR2
              SORD(J1)=SORD(J1)+ORD2
              DCP(J1)=DCP(J1)+DCP2
              if(IOD.gt.0)then
                DZ=SZ(IC)-ZRF
                call PBC(DX,DY,DZ)
                if(DZ.le.0.and.LREFL)then
                  CORD1=-CORD1
                  CORD1=-CORD2
                end if
                IN1=(CORD1+1.)*0.5*IOD+1
                IN2=(CORD2+1.)*0.5*IOD+1
                ORDR(IN1,J)=ORDR(IN1,J)+1.
                ORDR(IN2,J1)=ORDR(IN2,J1)+1.
              end if
            else                   ! change order
              CSHL(J)=CSHL(J)+RR2
              SORD(J)=SORD(J)+ORD2
              DCP(J)=DCP(J)+DCP2
              CSHL(J1)=CSHL(J1)+RR1
              SORD(J1)=SORD(J1)+ORD1
              DCP(J1)=DCP(J1)+DCP1
              if(IOD.gt.0)then
                DZ=SZ(IC)-ZRF
                call PBC(DX,DY,DZ)
                if(DZ.le.0.and.LREFL)then
                  CORD1=-CORD1
                  CORD1=-CORD2
                end if
                IN1=(CORD1+1.)*0.5*IOD+1
                IN2=(CORD2+1.)*0.5*IOD+1
                ORDR(IN1,J1)=ORDR(IN1,J1)+1.
                ORDR(IN2,J)=ORDR(IN2,J)+1.
              end if
            end if
          end do
          J=J1
        else             !  IOR = 0 (just this vector)
          do N=1,NML
	    IC = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+NC(J)-ISADR(ITYP)
	    IH = ISADDR(ITYP)+(N-1)*NSITS(ITYP)+NH(J)-ISADR(ITYP)
            DX = SX(IH)-SX(IC)
            DY = SY(IH)-SY(IC)
            DZ = SZ(IH)-SZ(IC)
            call PBC(DX,DY,DZ)
            RR = sqrt(DX**2+DY**2+DZ**2)
            CORD = DZ/RR                  ! Cos(theta)
            ORD = 0.5*(3*CORD**2-1.)      ! NMR order parameter
            CSHL(J)=CSHL(J)+RR            ! average bond length
            SORD(J)=SORD(J)+ORD
            DCP(J)= DCP(J)+ORD/RR**3       ! NMR dipole coupling
            if(IOD.gt.0)then
              DZ=SZ(IC)-ZRF
              call PBC(DX,DY,DZ)
              if(DZ.le.0.and.LREFL)CORD=-CORD
              IN1=(CORD+1.)*0.5*IOD+1
              ORDR(IN1,J)=ORDR(IN1,J)+1.
            end if
          end do
        end if
        go to 5
 50     continue
      end do     ! while IEND
      open(unit=26,file=FORD,status='unknown')
      write(26,*)' Order parameter calculations'
      write(26,*)' Configurations sampled: ',IAN
      write(26,*)'  Bond                 order       length ',
     &'  DC(1/A**3)'
      do J=1,NEQH
        ITYP=ITS(NC(J))
        NML=NSPEC(ITYP)
        FAC=1./(IAN*NML)
        SO=SORD(J)*FAC
        BL=CSHL(J)*FAC
        DC=DCP(J)*FAC
        write(26,'(a4,a1,a4,2i5,3f12.4)')
     +  NM(NC(J)),'-',NM(NH(J)),NC(J),NH(J),SO,BL,DC
      end do
      write(26,*)
      write(26,*)' Distribution of Cos(theta)          '
      write(26,*)
      do J=1,NEQH
        ITYP=ITS(NC(J))
        NML=NSPEC(ITYP)
        FAC=IOD*1./(IAN*NML)
        write(26,*)
        write(26,'(a,a4,a1,a4,2i5)')
     &  ' vector ',NM(NC(J)),'-',NM(NH(J)),NC(J),NH(J)
        do I=1,IOD
          COT=-1.+(I-0.5)*2./IOD
          write(26,'(2f12.5,5x,a4,a1,a4)')
     &    COT,ORDR(I,J)*FAC,NM(NC(J)),':',NM(NH(J))
        end do
        write(26,*)'-------------------------------------------'
      end do
      stop
 100  write(*,*)' Input file ',FORDIN,' for order parameter ',
     & 'calculations not found '
      stop
 20   write(*,*)' Input error in the order parameter file ',FORDIN
      write(*,*)STR
      stop
      end
