*=============== INTERNAL FORCES =============================================
*
*=============== GETBND ===========================================
*
      SUBROUTINE GETBND(ITYP)
*
      include "mdee.h"
* 
      character*3 NAMN
      DIMENSION BR(NB),BE(NB)
*                    
      NRBS      = NRB(ITYP)
      IF(NRBS.LE.0) RETURN
      timeb0	= cputime(0.d0)
*
      ISBEG     = ISADDR(ITYP)+1
      ISEND     = ISADDR(ITYP +1)
*     
      VIR1	= 0.
      NAMN=NAME(ITYP)(1:3)
      NSS       = NSITS (ITYP)
      IF(NAMN.eq.'H2O'.or.NAMN.eq.'h2o'.or.NAMN.eq.'HOH'
     +.or.NAMN.eq.'hoh')THEN
        CALL STRBND(ITYP)
        timeb	= timeb+cputime(timeb0)
        RETURN
      END IF
*
      NBBEG     = IADB(ITYP)
      NBEND     = NBBEG+NRBS-1
*
      DO M      = NBBEG,NBEND
        BR(M)     = 0.D0
        BE(M)     = 0.D0
      END DO
      NSP       = NSPEC (ITYP)
      FNSPI     = 1.D0/DFLOAT(NSP)
*
      ISHF      = ISADDR(ITYP)
      DO M      = NBBEG,NBEND
        IBB       = IB(M)+ISHF
        JBB       = JB(M)+ISHF
        REQ       = RB(M)
        FCN       = FB(M)
*
        DO K      = 1,NSP
          I         = (K-1)*NSS+IBB
          J         = (K-1)*NSS+JBB
          BX        = SX(J)-SX(I)
          BY        = SY(J)-SY(I)
          BZ        = SZ(J)-SZ(I)
          call PBC(BX,BY,BZ)
          BSQ       = BX*BX+BY*BY+BZ*BZ
          BD        = DSQRT(BSQ)
*     if(ityp.eq.1.and.m.eq.nbbeg) print *,bd
          R1I       = 1.D0/BD
          BBB       =  BD-REQ
          DBND      = BBB*FCN
          EBDD	    = DBND*BBB
          FORB      =  -2.0D0*DBND*R1I
          BR(M)     = BR(M)+BD
          BE(M)     = BE(M)+EBDD
          PINT	    = PINT+EBDD
*          write(*,*)I,J,IS,JS,NM(IS),NM(JS),BD,EBDD*0.001*ENERF
          HX(I)     = HX(I)-BX*FORB
          HY(I)     = HY(I)-BY*FORB
          HZ(I)     = HZ(I)-BZ*FORB
          HX(J)     = HX(J)+BX*FORB
          HY(J)     = HY(J)+BY*FORB
          HZ(J)     = HZ(J)+BZ*FORB
          VIR1	    = VIR1+FORB*BSQ
        END DO! OF K
*
      END DO! OF M
      VIRS	= VIRS+VIR1
      VIRB	= VIRB+VIR1
*
      DO M      = NBBEG,NBEND
        BB(M)     = BB(M)+BR(M)*FNSPI
        EB(M)     = EB(M)+BE(M)*FNSPI
      END DO
*
      timeb	= timeb+cputime(timeb0)
      RETURN
      END
*
*=============== GETANG ===========================================
*
      SUBROUTINE GETANG(ITYP)
*
      include "mdee.h"
*
      DIMENSION AR(NA),AE(NA)
*
      NRAS      = NRA(ITYP)
      IF(NRAS.LE.0) RETURN
*
      timea0	= cputime(0.d0)
      NABEG     = IADA(ITYP)
      NAEND     = NABEG+NRAS-1
*
      DO M      = NABEG,NAEND
        AR(M)     = 0.D0
        AE(M)     = 0.D0
      END DO
*
      ISHF      = ISADDR(ITYP)
      IBEG      = IADDR (ITYP)+1
      IEND      = IADDR (ITYP+1)
      NSP       = NSPEC (ITYP)
      NSS       = NSITS (ITYP)
      FNSPI     = 1.D0/DFLOAT(NSP)
*
      DO M      = NABEG,NAEND
        II        = IA(M)+ISHF
        JJ        = JA(M)+ISHF
        KK        = KA(M)+ISHF
        FCN       = FA(M)
        AEQ       = RA(M)
        DO L      = 1,NSP
          I         = (L-1)*NSS+II
          J         = (L-1)*NSS+JJ
          K         = (L-1)*NSS+KK
*
          AX        = SX(I)-SX(J)
          AY        = SY(I)-SY(J)
          AZ        = SZ(I)-SZ(J)
          call PBC(AX,AY,AZ)
          ASQ       = AX*AX+AY*AY+AZ*AZ
          ASQI      = 1.D0/ASQ
*
          BX        = SX(K)-SX(J)
          BY        = SY(K)-SY(J)
          BZ        = SZ(K)-SZ(J)
          call PBC(BX,BY,BZ)
          BSQ       = BX*BX+BY*BY+BZ*BZ
          BSQI      = 1.D0/BSQ
*
          AB        = DSQRT(ASQ*BSQ)
          ABI       = 1.D0/AB
*
          COSA      = AX*BX+AY*BY+AZ*BZ
          COSA      = COSA*ABI
          COSA      = DMIN1(COSA, 1.D0)
          COSA      = DMAX1(COSA,-1.D0)
          ALF       = DACOS(COSA)
          AVA       = AVA+ALF
          SINA      = DSIN(ALF)
          SINAI     =  1.D0/SINA
          DALF      =  ALF-AEQ
          ABC       = DALF*FCN
          EANG	    = ABC*DALF	
*
          AR(M)     = AR(M)+ALF
          AE(M)     = AE(M)+EANG
          PINT	    = PINT+EANG 
*
*      if(nstep.eq.2)
*     x print  "(2i5,5x,3(1x,a4),f8.2,3x,f18.8)"
*     X ,m,l,nam(i),nam(j),nam(k)
*     x ,alf*180./3.1412,ABC*DALF*enerf*1.d-3
          SAB       =-2.D0*ABC*SINAI
          CAB       = SAB*COSA
          FAB       = SAB*ABI
          FAA       = CAB*ASQI
          FBB       = CAB*BSQI
*
          FAX       = FAB*BX-FAA*AX
          FAY       = FAB*BY-FAA*AY
          FAZ       = FAB*BZ-FAA*AZ
*
          FBX       = FAB*AX-FBB*BX
          FBY       = FAB*AY-FBB*BY
          FBZ       = FAB*AZ-FBB*BZ
*
          HX(I)     = HX(I)-FAX
          HY(I)     = HY(I)-FAY
          HZ(I)     = HZ(I)-FAZ
          HX(J)     = HX(J)+FAX+FBX
          HY(J)     = HY(J)+FAY+FBY
          HZ(J)     = HZ(J)+FAZ+FBZ
          HX(K)     = HX(K)-FBX
          HY(K)     = HY(K)-FBY
          HZ(K)     = HZ(K)-FBZ
        END DO! OF NN
      END DO! OF M
*
      DO M      = NABEG,NAEND
      AA(M)     = AA(M)+AR(M)*FNSPI
      EA(M)     = EA(M)+AE(M)*FNSPI
      END DO
*
      timea	= timea+cputime(timea0)
      RETURN
      END
*
*=============== GETTRS ============================================== 
*
      SUBROUTINE GETTRS(ITYP)
*
      include "mdee.h"
*
      DIMENSION TR(NT),TE(NT)
*
      NRTS      = NRT(ITYP)
      IF(NRTS.LE.0) RETURN
      timet0	= cputime(0.d0)
      IF(ITORS(ITYP).NE.0) THEN
      CALL GETTOR(ITYP)
      timet	= timet+cputime(timet0)
      RETURN
      END IF
*
      NTBEG     = IADT(ITYP)
      NTEND     = NTBEG+NRTS-1
      VIR1	= 0.
*
      DO M      = NTBEG,NTEND
        TR(M)     = 0.D0
        TE(M)     = 0.D0
      END DO
*
      ISHF      = ISADDR(ITYP)
      IBEG      = IADDR (ITYP)+1
      IEND      = IADDR (ITYP+1)
      NSP       = NSPEC (ITYP)
      NSS       = NSITS (ITYP)
      FNSPI     = 1.D0/DFLOAT(NSP)
*
      DO M      = NTBEG,NTEND
        II        = IT(M)+ISHF
        JJ        = JT(M)+ISHF
        KK        = KT(M)+ISHF
        LL        = LT(M)+ISHF
        FCN       = FT(M)
        TEQ       = RT(M)
        MUL       = NMUL(M)
        DO N      = 1,NSP
          I         = (N-1)*NSS+II
          J         = (N-1)*NSS+JJ
          K         = (N-1)*NSS+KK
          L         = (N-1)*NSS+LL
*      if(n.eq.1) print *,'torsio  ',nam(i),nam(j),nam(k),nam(l)
*
          ax        = sx(j)-sx(i)
          ay        = sy(j)-sy(i)
          az        = sz(j)-sz(i)
          bx        = sx(k)-sx(j)
          by        = sy(k)-sy(j)
          bz        = sz(k)-sz(j)
          cx        = sx(l)-sx(k)
          cy        = sy(l)-sy(k)
          cz        = sz(l)-sz(k)
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
*
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
            si        = dsin(arg)
            if( dabs(si).lt.1.0d-12 ) si = dsign(1.0d-12,si)
*
            EARG      = MUL*ARG-TEQ
*
*...Energies:
*
            POTT      = FCN*(1.D0+DCOS(EARG))
            TR(M)     = TR(M)+ARG
            TE(M)     = TE(M)+POTT
            PINT	= PINT+POTT
*
*...Forces:
*
            HELP      = -FCN*MUL*DSIN(EARG)
            de1       = 1.0d0/den/si*HELP
            axb       = axb/den*co
            bxc       = bxc/den*co
*X
      dnum      =  cx*bt - bx*bc
      dden      = (ab*bx - ax*bt)*bxc
      FFI1      = (dnum  - dden) * de1
      dnum      = ((bx-ax)*bc - ab*cx ) + (2.0*ac*bx - cx*bt)
      dden      = axb*(bc*cx-bx*ct) + (ax*bt-at*bx-ab*(bx-ax))*bxc
      FFJ1      = (dnum - dden) * de1
      dnum      = ab*bx - ax*bt
      dden      = axb*( bt*cx - bc*bx )
      FFL1      = (dnum - dden) * de1
      FFK1      = -(ffi1+ffj1+ffl1)
*
      HX(I)     = HX(I)+ffi1
      HX(J)     = HX(J)+ffj1
      HX(K)     = HX(K)+ffk1
      HX(L)     = HX(L)+ffl1
*Y
      dnum      =  cy*bt - by*bc
      dden      = (ab*by - ay*bt)*bxc
      FFI2      = (dnum  - dden) * de1
      dnum      = ((by-ay)*bc - ab*cy ) + (2.0*ac*by - cy*bt)
      dden      = axb*(bc*cy-by*ct) + (ay*bt-at*by-ab*(by-ay))*bxc
      FFJ2      = (dnum - dden) * de1
      dnum      = ab*by - ay*bt
      dden      = axb*( bt*cy - bc*by )
      FFL2      = (dnum - dden) * de1
      FFK2      = -(ffi2+ffj2+ffl2)
      HY(I)     = HY(I)+ffi2
      HY(J)     = HY(J)+ffj2
      HY(K)     = HY(K)+ffk2
      HY(L)     = HY(L)+ffl2
*Z
      dnum      =  cz*bt - bz*bc
      dden      = (ab*bz - az*bt)*bxc
      FFI3      = (dnum  - dden) * de1
      dnum      = ((bz-az)*bc - ab*cz ) + (2.0*ac*bz - cz*bt)
      dden      = axb*(bc*cz-bz*ct) + (az*bt-at*bz-ab*(bz-az))*bxc
      FFJ3      = (dnum - dden) * de1
      dnum      = ab*bz - az*bt
      dden      = axb*( bt*cz - bc*bz )
      FFL3      = (dnum - dden) * de1
      FFK3      = -(ffi3+ffj3+ffl3)
      HZ(I)     = HZ(I)+ffi3
      HZ(J)     = HZ(J)+ffj3
      HZ(K)     = HZ(K)+ffk3
      HZ(L)     = HZ(L)+ffl3
*
          endif ! den>0
        END DO! OF N
      END DO! OF M
*
      DO M      = NTBEG,NTEND
        TT(M)     = TT(M)+TR(M)*FNSPI
        ET(M)     = ET(M)+TE(M)*FNSPI
      END DO
*
*
      timet	= timet+cputime(timet0)
      RETURN
      END
*
*=============== GETTOR ============================================
*
      SUBROUTINE GETTOR(ITYP)
*
      include "mdee.h"
*
      DIMENSION TR(NT),TE(NT)
*
      NRTS  = NRT(ITYP)
      IF(NRTS.LE.0) RETURN
*
      NTBEG = IADT(ITYP)
      NTEND = NTBEG+NRTS-1
      VIR1	= 0.
*
      DO M  = NTBEG,NTEND
      TR(M) = 0.D0
      TE(M) = 0.D0
      END DO
*
      ISHF  = ISADDR(ITYP)
      IBEG  = IADDR (ITYP)+1
      IEND  = IADDR (ITYP+1)
      NSP   =  NSPEC (ITYP)
      NSS   = NSITS (ITYP)
      FNSPI = 1.D0/DFLOAT(NSP)
*
      DO M  = NTBEG,NTEND
      II    = IT(M)+ISHF
      JJ    = JT(M)+ISHF
      KK    = KT(M)+ISHF
      LL    = LT(M)+ISHF
      FN1   = 0.5D0*FT1(M)
      FN2   = 0.5D0*FT2(M)
      FN3   = 0.5D0*FT3(M)
      DO N  = 1,NSP
      I     = (N-1)*NSS+II
      J     = (N-1)*NSS+JJ
      K     = (N-1)*NSS+KK
      L     = (N-1)*NSS+LL
*
      AX    = SX(J)-SX(I)
      AY    = SY(J)-SY(I)
      AZ    = SZ(J)-SZ(I)
*
      BX    = SX(K)-SX(J)
      BY    = SY(K)-SY(J)
      BZ    = SZ(K)-SZ(J)
*
      CX    = SX(L)-SX(K)
      CY    = SY(L)-SY(K)
      CZ    = SZ(L)-SZ(K)
      call PBC(AX,AY,AZ)
      call PBC(BX,BY,BZ)
      call PBC(CX,CY,CZ)
*
       AXAX = AX*AX 
       BXBX = BX*BX 
       CXCX = CX*CX
       AYAY = AY*AY 
       BYBY = BY*BY 
       CYCY = CY*CY
       AZAZ = AZ*AZ 
       BZBZ = BZ*BZ 
       CZCZ = CZ*CZ
*
       AXBX = AX*BX 
       AXCX = AX*CX 
       BXCX = BX*CX
       AXBY = AX*BY 
       AXCY = AX*CY 
       BXCY = BX*CY
       AXBZ = AX*BZ 
       AXCZ = AX*CZ 
       BXCZ = BX*CZ
*
       AYBX = AY*BX 
       AYCX = AY*CX 
       BYCX = BY*CX
       AYBY = AY*BY 
       AYCY = AY*CY 
       BYCY = BY*CY
       AYBZ = AY*BZ 
       AYCZ = AY*CZ 
       BYCZ = BY*CZ
*
       AZBX = AZ*BX 
       AZCX = AZ*CX 
       BZCX = BZ*CX
       AZBY = AZ*BY 
       AZCY = AZ*CY 
       BZCY = BZ*CY
       AZBZ = AZ*BZ 
       AZCZ = AZ*CZ 
       BZCZ = BZ*CZ
*
      DX    = AYBZ-AZBY
      DY    = AZBX-AXBZ
      DZ    = AXBY-AYBX
*
      EX    = BYCZ-BZCY
      EY    = BZCX-BXCZ
      EZ    = BXCY-BYCX
*
      D2    = DX*DX+DY*DY+DZ*DZ
      D1    = DSQRT(D2)
*
      E2    = EX*EX+EY*EY+EZ*EZ
      E1    = DSQRT(E2)
*
      DE    = DSQRT(D2*E2)
*
      DDOTE = DX*EX+DY*EY+DZ*EZ
      DII   = 1.D0/D1
      EII   = 1.D0/E1
      DEI   = DII*EII
*
      COS1F = DDOTE*DEI
      COS1F = DMIN1(COS1F, 1.D0)
      COS1F = DMAX1(COS1F,-1.D0)
      COS2F =  2.D0*COS1F*COS1F-1.D0
      COS3F = COS1F*(2.D0*COS2F-1.D0)
      COS1F2= COS1F*COS1F
      COS2F2= COS2F*COS2F
      COS3F2= COS3F*COS3F
      ARG   = DACOS(COS1F)
*
      FX1   = BZ*EY-BY*EZ
      FY1   = BX*EZ-BZ*EX
      FZ1   = BY*EX-BX*EY
*
      FX4   = BZ*DY-BY*DZ
      FY4   = BX*DZ-BZ*DX
      FZ4   = BY*DX-BX*DY
*
      FX3   =-2.D0*BX*(AYCY+AZCZ)+AX*(BYCY+BZCZ)+CX*(AYBY+AZBZ)
      FY3   =-2.D0*BY*(AXCX+AZCZ)+AY*(BXCX+BZCZ)+CY*(AXBX+AZBZ)
      FZ3   =-2.D0*BZ*(AXCX+AYCY)+AZ*(BXCX+BYCY)+CZ*(AXBX+AYBY)
*
      FX2   =-FX1-FX3 
      FY2   =-FY1-FY3 
      FZ2   =-FZ1-FZ3 
*
      GX1   = 2.D0*FX4
      GY1   = 2.D0*FY4
      GZ1   = 2.D0*FZ4
*
      GX3   = 2.D0*(AZ*DY-DZ*AY)
      GY3   = 2.D0*(AX*DZ-DX*AZ)
      GZ3   = 2.D0*(AY*DX-DY*AX)
*
      GX2   = -GX1-GX3
      GY2   = -GY1-GY3
      GZ2   = -GZ1-GZ3
*
      HX2   = 2.D0*(CZ*EY-CZ*EY)
      HY2   = 2.D0*(CX*EZ-CX*EZ)
      HZ2   = 2.D0*(CY*EX-CY*EX)
*
      HX4   = 2.D0*FX1
      HY4   = 2.D0*FY1
      HZ4   = 2.D0*FZ1
*
      UTORS = FN1 * (1.D0 + COS1F)
     X       +FN2 * (1.D0 - COS2F)
     X       +FN3 * (1.D0 + COS3F)
*
      TR(M)     = TR(M)+ARG
      TE(M)     = TE(M)+UTORS
      PINT	= PINT+UTORS
*
      FOR   = FN1 - 4.D0*FN2*COS1F + 3.D0*FN3*(2.D0*COS2F+1.D0)
      FOR   = -DEI*FOR
      FORK  = 0.5D0*COS1F*E1*DII
      FORL  = 0.5D0*COS1F*D1*EII
*
      F1X   = FOR*(FX1-FORK*GX1)
      F1Y   = FOR*(FY1-FORK*GY1)
      F1Z   = FOR*(FZ1-FORK*GZ1)
*
      F2X   = FOR*(FX2-FORK*GX2-FORL*HX2)
      F2Y   = FOR*(FY2-FORK*GY2-FORL*HY2)      
      F2Z   = FOR*(FZ2-FORK*GZ2-FORL*HZ2)
*
      F4X   = FOR*(FX4-FORL*HX4)
      F4Y   = FOR*(FY4-FORL*HY4)
      F4Z   = FOR*(FZ4-FORL*HZ4)
*
      HX(I) = HX(I)+F1X
      HY(I) = HY(I)+F1Y
      HZ(I) = HZ(I)+F1Z
*
      HX(J) = HX(J)+F2X
      HY(J) = HY(J)+F2Y
      HZ(J) = HZ(J)+F2Z
*
      HX(K) = HX(K)-F1X-F2X-F4X
      HY(K) = HY(K)-F1Y-F2Y-F4Y
      HZ(K) = HZ(K)-F1Z-F2Z-F4Z
*
      HX(L) = HX(L)+F4X
      HY(L) = HY(L)+F4Y
      HZ(L) = HZ(L)+F4Z
*
      END DO! OF N
      END DO! OF M
*
      DO M      = NTBEG,NTEND
      TT(M)     = TT(M)+TR(M)*FNSPI
      ET(M)     = ET(M)+TE(M)*FNSPI
      END DO
*
*
*
      RETURN
      END
*
*=============== GETIMP ===============================================
*
      SUBROUTINE GETIMP(ITYP)
*
      include "mdee.h"
*
      DIMENSION TR(NT),TE(NT)
*
      NRIS      = NRI(ITYP)
      IF(NRIS.LE.0) RETURN
      timei0	= cputime(0.d0)
*
      NTBEG     = IADI(ITYP)
      NTEND     = NTBEG+NRIS-1
*
      DO M      = NTBEG,NTEND
      TR(M)     = 0.D0
      TE(M)     = 0.D0
      END DO
      VIR1	= 0.
*
      ISHF      = ISADDR(ITYP)
      IBEG      = IADDR (ITYP)+1
      IEND      = IADDR (ITYP+1)
      NSP       = NSPEC (ITYP)
      NSS       = NSITS (ITYP)
      FNSPI     = 1.D0/DFLOAT(NSP)
*
      DO M      = NTBEG,NTEND
      II        = IM(M)+ISHF
      JJ        = JM(M)+ISHF
      KK        = KM(M)+ISHF
      LL        = LM(M)+ISHF
      FCN       = FI(M)
      TEQ       = RI(M)
      DO N      = 1,NSP
      I         = (N-1)*NSS+II
      J         = (N-1)*NSS+JJ
      K         = (N-1)*NSS+KK
      L         = (N-1)*NSS+LL
*
      AX        = SX(I)-SX(J)
      AY        = SY(I)-SY(J)
      AZ        = SZ(I)-SZ(J)
*
      BX        = SX(K)-SX(J)
      BY        = SY(K)-SY(J)
      BZ        = SZ(K)-SZ(J)
*
      CX        = SX(K)-SX(L)
      CY        = SY(K)-SY(L)
      CZ        = SZ(K)-SZ(L)
*
      DX        = SX(I)-SX(K)
      DY        = SY(I)-SY(K)
      DZ        = SZ(I)-SZ(K)
*
      EX        = SX(J)-SX(L)
      EY        = SY(J)-SY(L)
      EZ        = SZ(J)-SZ(L)
*
      AKBX      = AY*BZ-AZ*BY
      AKBY      = AZ*BX-AX*BZ
      AKBZ      = AX*BY-AY*BX
*
      BKCX      = BY*CZ-BZ*CY
      BKCY      = BZ*CX-BX*CZ
      BKCZ      = BX*CY-BY*CX
*
      AKB2      = AKBX*AKBX+AKBY*AKBY+AKBZ*AKBZ
      AKB       = DSQRT(AKB2)
*
      BKC2      = BKCX*BKCX+BKCY*BKCY+BKCZ*BKCZ
      BKC       = DSQRT(BKC2)
*
      AKBKC     = SQRT(AKB2*BKC2)
*
      HH        = AKBX*BKCX+AKBY*BKCY+AKBZ*BKCZ
      AI        = 1.D0/AKB
      BI        = 1.D0/BKC
      ABI       = AI*BI
*
      COSF      = HH*ABI
      COSF      = DMIN1(COSF, 1.D0)
      COSF      = DMAX1(COSF,-1.D0)
      PHI       = DCOS (COSF)
      DFI       = PHI-TEQ
      ENER	= FCN*DFI**2
      POTI      = POTI+ENER
      TR(M)     = TR(M)+PHI
      TE(M)     = TE(M)+ENER
      PINT	= PINT+ENER
*
      P         = 2.D0*FCN*DFI
      PN        = DSIN(PHI)
*!    PN        =  DMAX1(PN,1.D-12)
*!    PN        =  DMIN1(PN,1.D+12)
      P         =-P/PN
      PH        = P*COSF
      PKL       = P/AKBKC
      PKK       = PH/AKB2
      PLL       = PH/BKC2
*
      XK        = PKL*BKCX-PKK*AKBX
      YK        = PKL*BKCY-PKK*AKBY
      ZK        = PKL*BKCZ-PKK*AKBZ
*
      XL        = PKL*AKBX-PLL*BKCX
      YL        = PKL*AKBY-PLL*BKCY
      ZL        = PKL*AKBZ-PLL*BKCZ
      P         =-P/PN
      PH        = P*COSF
      PKL       = P/AKBKC
      PKK       = PH/AKB2
      PLL       = PH/BKC2
*
      XAA       = BY*ZK-YK*BZ
      YAA       = BZ*XK-ZK*BX
      ZAA       = BX*YK-XK*BY
*
      XBB       = DY*ZK-YK*DZ-CY*ZL+YL*CZ
      YBB       = DZ*XK-ZK*DX-CZ*XL+ZL*CX
      ZBB       = DX*YK-XK*DY-CX*YL+XL*CY
*
      XCC       = EY*ZL-YL*EZ-AY*ZK+YK*AZ
      YCC       = EZ*XL-ZL*EX-AZ*XK+ZK*AX
      ZCC       = EX*YL-XL*EY-AX*YK+XK*AY
*
      XDD       = BY*ZL-YL*BZ
      YDD       = BZ*XL-ZL*BX
      ZDD       = BX*YL-XL*BY
*
      HX(I)     = HX(I)-XAA
      HY(I)     = HY(I)-YAA
      HZ(I)     = HZ(I)-ZAA
*
      HX(J)     = HX(J)-XBB
      HY(J)     = HY(J)-YBB
      HZ(J)     = HZ(J)-ZBB
*
      HX(K)     = HX(K)-XCC
      HY(K)     = HY(K)-YCC
      HZ(K)     = HZ(K)-ZCC
*
      HX(L)     = HX(L)-XDD
      HY(L)     = HY(L)-YDD
      HZ(L)     = HZ(L)-ZDD
      VIR1	= VIR1-XAA*SX(I)-YAA*SY(I)-ZAA*SZ(I)-
     -                 XBB*SX(J)-YBB*SY(J)-ZBB*SZ(J)-
     -                 XCC*SX(K)-YCC*SY(K)-ZCC*SZ(K)-
     -                 XDD*SX(L)-YDD*SY(L)-ZDD*SZ(L)
*
      END DO! OF N
      END DO! OF M
*
      DO M      = NTBEG,NTEND
      TI(M)     = TI(M)+TR(M)*FNSPI
      EI(M)     = EI(M)+TE(M)*FNSPI
      END DO
      VIRS	= VIRS+VIR1
	VIRT	= VIRT+VIR1
      timei	= timei+cputime(timei0)
*
      RETURN
      END
*
*=============== STRBND ===========================================
*
      SUBROUTINE STRBND(ITYP)
*      
* This routine calculates intramolecular forces for 
* triatomic molecules (like water)
*
      include "mdee.h"
*
      DIMENSION BS(NB),BE(NB)
*
*.... ANHARMONIC FLEXIBLE WATER MODEL:
*
      PARAMETER (AAP  = 2809.56D0 , BBP = 2809.56D0)
      PARAMETER (CCP  =  687.41D0 , DDP =   2.566D0)
      PARAMETER (EEP  = -884.63D0 , FFP = -884.63D0)
      PARAMETER (GGP  =  467.31D0 ,BBMA=1.3,BBSA=1.4  )
*
      IF(NSITS(ITYP).NE.3) STOP '!!! NOT A WATER MOLECULE'
      IF( IPOT(ITYP).GT.2) STOP '!!! IPOT.GT.2  - NOT SUPPORTED! '
*
      NBBEG    = IADB(ITYP)
      NBEND    = NBBEG+NRB(ITYP)-1
*
      I1       = NBBEG
      I2       = I1+1
      I3       = I2+1
*
      VIR1	= 0.
      IBEG     = IADDR(ITYP)+1
      IEND     = IADDR(ITYP+1)
      NSP      = NSPEC(ITYP)
      NSS      = NSITS(ITYP)
      FNSPI    = 1.D0/DFLOAT(NSP)
      FACT     = 1.D3/AVSNO/UNITE
      FACTOR   = FACT
*
      ROH1     = RB(I1)
      ROH2     = RB(I2)
      RH1H2    = RB(I3)
*
      AAA      = AAP*FACTOR
      BBB      = BBP*FACTOR
      DDD      = DDP
      IPT      = IPOT(ITYP)
      IF(IPT.EQ.2) THEN
        AAA      = AAA/DDD**2
        BBB      = BBB/DDD**2
      END IF
*
      CCC      = CCP*FACTOR
      EEE      = EEP*FACTOR
      FFF      = FFP*FACTOR
      GGG      = GGP*FACTOR
      BBM	= BBMA
      BBS	= BBSA
*
      EFAC     =  FKCALJ*1.D+3/AVSNO/UNITE
*
      DO M     = NBBEG,NBEND
      BS(M)    = 0.D0
      BE(M)    = 0.D0
      END DO
*
      AADD     = AAA*DDD
      BBDD     = BBB*DDD
*
      ISHF     = ISADDR(ITYP)
      IBB      = ISHF+1
      JBB      = ISHF+2
      KBB      = ISHF+3
*
      IF(IPOT(ITYP).EQ.2) THEN
*
      DO N     = 1,NSP
      I        = (N-1)*NSS+IBB
      J        = (N-1)*NSS+JBB
      K        = (N-1)*NSS+KBB 
*
      AX       = SX(J)-SX(I)
      AY       = SY(J)-SY(I)
      AZ       = SZ(J)-SZ(I)
*
      BX       = SX(K)-SX(I)
      BY       = SY(K)-SY(I)
      BZ       = SZ(K)-SZ(I)
*
      CX       = SX(K)-SX(J)
      CY       = SY(K)-SY(J)
      CZ       = SZ(K)-SZ(J)
*
      call PBC(AX,AY,AZ)
      call PBC(BX,BY,BZ)
      call PBC(CX,CY,CZ)
      RRA      = AX*AX+AY*AY+AZ*AZ
      RRB      = BX*BX+BY*BY+BZ*BZ
      RRC      = CX*CX+CY*CY+CZ*CZ
      ASS      = DSQRT(RRA)
      BSS      = DSQRT(RRB)
      if(ASS.gt.BBS)write(*,*)' !!! bond length ',ASS,' for ',I,J
      if(BSS.gt.BBS)write(*,*)' !!! bond length ',BSS,' for ',I,K
      CSS      = DSQRT(RRC)
      COSA     = (AX*BX+AY*BY+AZ*BZ)/ASS/BSS
      ARG      = DACOS(COSA)
*
      SA       = ASS-ROH1
      SB       = BSS-ROH2
      SH       = CSS-RH1H2
*
      EXPA     = DEXP(-DDD*SA)
      EXPAM    = 1.D0-EXPA
      EXPB     = DEXP(-DDD*SB)
      EXPBM    = 1.D0-EXPB
*
      WB1       = AAA*EXPAM**2+EEE*SA*SH
      WB2       = BBB*EXPBM**2+FFF*SB*SH
      WB3       = CCC*SH*SH   +GGG*SA*SB
      WB        = WB1+WB2+WB3
*
      BS(I1)    = BS(I1)+ASS
      BS(I2)    = BS(I2)+BSS
      BS(I3)    = BS(I3)+CSS
*!    BS(I3)    = BS(I3)+ARG
*
      BE(I1)    = BE(I1)+WB1
      BE(I2)    = BE(I2)+WB2
      BE(I3)    = BE(I3)+WB3 
	EWB	= WB1+WB2+WB3
      PINT	= PINT+EWB
*	if(IPRINT.ge.7)write(*,'(I5,f13.5,3f10.5)')
*     +N,0.001*EWB*ENERF,ASS,BSS,CSS
*
      if(ASS.gt.BBM)then
        FB3A	= -2.*SA*AADD*DDD
      else
        FB3A	= -2.*AADD*EXPA*EXPAM
      end if
      FB3A     = FB3A   - EEE*SH-GGG*SB
      if(BSS.gt.BBM)then
        FB3B	= -2.*SB*BBDD*DDD
      else
        FB3B	= -2.*BBDD*EXPB*EXPBM
      end if
      FB3B     = FB3B   - FFF*SH-GGG*SA
      FB3C     = -2.0*CCC*SH            - EEE*SA-FFF*SB
*
      FB3A     = FB3A/ASS
      FB3B     = FB3B/BSS
      FB3C     = FB3C/CSS
*
      AXX       = FB3A*AX
      AYY       = FB3A*AY
      AZZ       = FB3A*AZ
*
      BXX       = FB3B*BX
      BYY       = FB3B*BY
      BZZ       = FB3B*BZ
*
      CXX       = FB3C*CX
      CYY       = FB3C*CY
      CZZ       = FB3C*CZ
*
      HX(I)    = HX(I)-AXX-BXX
      HY(I)    = HY(I)-AYY-BYY
      HZ(I)    = HZ(I)-AZZ-BZZ
*
      HX(J)    = HX(J)+AXX   -CXX
      HY(J)    = HY(J)+AYY   -CYY
      HZ(J)    = HZ(J)+AZZ   -CZZ
*
      HX(K)    = HX(K)   +BXX+CXX
      HY(K)    = HY(K)   +BYY+CYY
      HZ(K)    = HZ(K)   +BZZ+CZZ
*
      VIR1	= VIR1+AXX*AX+BXX*BX+CXX*CX+
     +                 AYY*AY+BYY*BY+CYY*CY+
     +                 AZZ*AZ+BZZ*BZ+CZZ*CZ
      END DO! OF N
*
      ELSE
*
      DO N     = 1,NSP
      I        = (N-1)*NSS+IBB
      J        = (N-1)*NSS+JBB
      K        = (N-1)*NSS+KBB
*
      AX       = SX(J)-SX(I)
      AY       = SY(J)-SY(I)
      AZ       = SZ(J)-SZ(I)
*
      BX       = SX(K)-SX(I)
      BY       = SY(K)-SY(I)
      BZ       = SZ(K)-SZ(I)
*
      CX       = SX(K)-SX(J)
      CY       = SY(K)-SY(J)
      CZ       = SZ(K)-SZ(J)
      call PBC(AX,AY,AZ)
      call PBC(BX,BY,BZ)
      call PBC(CX,CY,CZ)
*
      RRA      = AX*AX+AY*AY+AZ*AZ
      RRB      = BX*BX+BY*BY+BZ*BZ
      RRC      = CX*CX+CY*CY+CZ*CZ
*
      ASS      = DSQRT(RRA)
      BSS      = DSQRT(RRB)
      CSS      = DSQRT(RRC)
      COSA     = (AX*BX+AY*BY+AZ*BZ)/ASS/BSS
      ARG      = DACOS(COSA)
*
      SA       = ASS-ROH1
      SB       = BSS-ROH2
      SH       = CSS-RH1H2
*
      WB1      = AAA*SA*SA+EEE*SA*SH
      WB2      = BBB*SB*SB+FFF*SB*SH 
      WB3      = CCC*SH*SH+GGG*SA*SB
      WB       = WB1+WB2+WB3
*
      BS(I1)    = BS(I1)+ASS
      BS(I2)    = BS(I2)+BSS
      BS(I3)    = BS(I3)+CSS
*!    BS(I3)    = BS(I3)+ARG
*
      BE(I1)    = BE(I1)+WB1
      BE(I2)    = BE(I2)+WB2
      BE(I3)    = BE(I3)+WB3
      PINT	= PINT+WB1+WB2+WB3
*
      FB3A     = -2.0*AAA*SA - EEE*SH-GGG*SB
      FB3B     = -2.0*BBB*SB - FFF*SH-GGG*SA
      FB3C     = -2.0*CCC*SH - EEE*SA-FFF*SB
*
      FB3A     = FB3A/ASS
      FB3B     = FB3B/BSS
      FB3C     = FB3C/CSS
*
      AXX       = FB3A*AX
      AYY       = FB3A*AY
      AZZ       = FB3A*AZ
*
      BXX       = FB3B*BX
      BYY       = FB3B*BY
      BZZ       = FB3B*BZ
*
      CXX       = FB3C*CX
      CYY       = FB3C*CY
      CZZ       = FB3C*CZ
*
      HX(I)    = HX(I)-AXX-BXX
      HY(I)    = HY(I)-AYY-BYY
      HZ(I)    = HZ(I)-AZZ-BZZ
*
      HX(J)    = HX(J)+AXX   -CXX
      HY(J)    = HY(J)+AYY   -CYY
      HZ(J)    = HZ(J)+AZZ   -CZZ
*
      HX(K)    = HX(K)   +BXX+CXX
      HY(K)    = HY(K)   +BYY+CYY
      HZ(K)    = HZ(K)   +BZZ+CZZ
*
      VIR1	= VIR1+AXX*AX+BXX*BX+CXX*CX+
     +                 AYY*AY+BYY*BY+CYY*CY+
     +                 AZZ*AZ+BZZ*BZ+CZZ*CZ
      END DO! OF N
*
      END IF
*
      DO M      = NBBEG,NBEND
      BB(M)     = BB(M)+BS(M)*FNSPI
      EB(M)     = EB(M)+BE(M)*FNSPI
      END DO! OF M
*      print *,bs(1)*fnspi
*     X       ,bs(2)*fnspi
*     X       ,bs(3)*fnspi
*
      VIRS	= VIRS+VIR1
      VIRB	= VIRB+VIR1
	RETURN
      END
*
