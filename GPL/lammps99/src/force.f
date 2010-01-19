c LAMMPS 99 - Molecular Dynamics Simulator
c Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
c Steve Plimpton, sjplimp@sandia.gov
c
c Copyright (2003) Sandia Corporation.  Under the terms of Contract
c DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
c certain rights in this software.  This software is distributed under 
c the GNU General Public License.
c
c See the README file in the top-level LAMMPS directory.

c -------------------------------------------------------------------------
c non-bonded forces
c  force components stored for atoms I own
c  variable force = (-1/r * d(phi)/dr)
c  variable virial = (-r * d(phi)/dr)

c cut LJ and no Coulombic

      subroutine lj_cut(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq
            r6inv = r2inv*r2inv*r2inv
            forcelj = r6inv*(lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            force = spfactor*forcelj*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and cut Coulombic

      subroutine lj_cut_coul_cut(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else
              forcecoul = 0.0
            endif

            if (rsq.lt.cutljsq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = spfactor * (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and long-range Coulombic

      subroutine lj_cut_coul_long(iflag)
      include "lammps.h"
      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich) - 2.0
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulsq) then
              r = sqrt(rsq)
              grij = gewsr*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)/r
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor.lt.1.0)
     $             forcecoul = forcecoul - (1.0-spfactor)*prefactor
            else
              forcecoul = 0.0
            endif

            if (rsq.lt.cutljsq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = spfactor *
     $             r6inv*(lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and no Coulombic

      subroutine lj_smooth(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = fskin*r
            endif

            force = spfactor*forcelj*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
          
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and cut Coulombic

      subroutine lj_smooth_coul_cut(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else
              forcecoul = 0.0
            endif

            if (rsq.lt.cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else if (rsq.lt.cutljsq(itype,jtype)) then
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = fskin*r
            else
              forcelj = 0.0
            endif

            force = spfactor*(forcecoul+forcelj)*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
          
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and long-range Coulombic

      subroutine lj_smooth_coul_long(iflag)
      include "lammps.h"
      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich) - 2.0
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulsq) then
              r = sqrt(rsq)
              grij = gewsr*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)/r
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor.lt.1.0)
     $             forcecoul = forcecoul - (1.0-spfactor)*prefactor
            else
              forcecoul = 0.0
            endif

            if (rsq.lt.cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = spfactor * r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else if (rsq.lt.cutljsq(itype,jtype)) then
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = spfactor * fskin*r
            else
              forcelj = 0.0
            endif

            force = (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif

        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c shifted LJ and no Coulombic

      subroutine lj_shift(iflag)
      include "lammps.h"
      
      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r = sqrt(rsq)
            rshift = r - lj5(itype,jtype)
            rshiftsq = rshift*rshift
            r2inv = 1.0/rshiftsq
            r6inv = r2inv*r2inv*r2inv
            forcelj = r6inv*(lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            force = spfactor*forcelj/rshift/r

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c soft cosine potential and no Coulombic
c set current values of ramped pre-factors before computing forces

      subroutine soft(iflag)
      include "lammps.h"
      data pi /3.1415926/

      do i = 1,ntypes
        do j = 1,ntypes
          lj4(i,j) = lj1(i,j) +
     $         float(itime)/nsteps*(lj2(i,j)-lj1(i,j))
        enddo
      enddo
      
      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r = sqrt(rsq)
            arg = pi*r/lj3(itype,jtype)
            force = spfactor * lj4(itype,jtype) *
     $           pi/lj3(itype,jtype)/r * sin(arg)

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c cut class II LJ and no Coulombic 

      subroutine ljclass2_cut(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            rinv = 1.0/sqrt(rsq)
            r2inv = rinv*rinv
            r3inv = r2inv*rinv
            r6inv = r3inv*r3inv
            forcelj = r6inv*(lj1(itype,jtype)*r3inv-lj2(itype,jtype))
            force = spfactor*forcelj*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut class II LJ and cut Coulombic

      subroutine ljclass2_cut_coul_cut(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            rinv = 1.0/sqrt(rsq)
            r2inv = rinv*rinv

            if (rsq.lt.cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*rinv
            else
              forcecoul = 0.0
            endif

            if (rsq.lt.cutljsq(itype,jtype)) then
              r3inv = r2inv*rinv
              r6inv = r3inv*r3inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r3inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = spfactor * (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut class II LJ and long-range Coulombic

      subroutine ljclass2_cut_coul_long(iflag)
      include "lammps.h"
      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich) - 2.0
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r = sqrt(rsq)
            rinv = 1.0/r
            r2inv = rinv*rinv

            if (rsq.lt.cutcoulsq) then
              grij = gewsr*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)*rinv
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor.lt.1.0)
     $             forcecoul = forcecoul - (1.0-spfactor)*prefactor
            else
              forcecoul = 0.0
            endif

            if (rsq.lt.cutljsq(itype,jtype)) then
              r3inv = r2inv*rinv
              r6inv = r3inv*r3inv
              forcelj = spfactor *
     $             r6inv*(lj1(itype,jtype)*r3inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and smooth Coulombic

      subroutine lj_cut_coul_smooth(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulintsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else if (rsq.lt.cutcoulsq) then
              r = sqrt(rsq)
              dd1 = (cutcoul - cutcoulint)**3
              switch_r = (cutcoul - r)**2 *
     $             (cutcoul + (2.0D0*r) - (3.0D0*cutcoulint)) / dd1
              deriv_sw_r = 6.0D0*(r - cutcoul)*(r - cutcoulint) / dd1
              forcecoul = (coulpre *qtmp*q(j)*sqrt(r2inv)*switch_r) -
     $             (coulpre*qtmp*q(j)*deriv_sw_r)
            else
              forcecoul = 0.0D0
            endif

            if (rsq.lt.cutljsq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = spfactor * (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and smooth Coulombic

      subroutine lj_smooth_coul_smooth(iflag)
      include "lammps.h"

      i = atompnt
      do ii = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nliststart(i),nliststop(i)

          j = nlist(k)

          if (j.le.maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          if (abs(delx).gt.xmc.and.perflagx.eq.0) then
            if (delx.lt.0.0) then
              delx = delx + xprd
            else
              delx = delx - xprd
            endif
          endif
          if (abs(dely).gt.ymc.and.perflagy.eq.0) then
            if (dely.lt.0.0) then
              dely = dely + yprd
            else
              dely = dely - yprd
            endif
          endif
          if (abs(delz).gt.zmc.and.perflagz.eq.0) then
            if (delz.lt.0.0) then
              delz = delz + zprd
            else
              delz = delz - zprd
            endif
          endif
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq.lt.cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq.lt.cutcoulintsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else if (rsq.lt.cutcoulsq) then
              r = sqrt(rsq)
              dd1 = (cutcoul - cutcoulint)**3
              switch_r = (cutcoul - r)**2 *
     $             (cutcoul + (2.0D0*r) - (3.0D0*cutcoulint)) / dd1
              deriv_sw_r = 6.0D0*(r - cutcoul)*(r - cutcoulint) / dd1
              forcecoul = (coulpre *qtmp*q(j)*sqrt(r2inv)*switch_r) -
     $             (coulpre*qtmp*q(j)*deriv_sw_r)
            else
              forcecoul = 0.0D0
            endif

            if (rsq.lt.cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else if (rsq.lt.cutljsq(itype,jtype)) then
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = fskin*r
            else
              forcelj = 0.0
            endif

            force = spfactor*(forcecoul+forcelj)*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond.eq.1.or.j.le.maxown) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag.eq.1) then
              if (newton_nonbond.eq.1.or.j.le.maxown) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
          
        enddo
        i = list(i)
      enddo

      return
      end
