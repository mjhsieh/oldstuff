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
c bond forces
c  force components stored for atoms I own
c  variable force = (-1/r * d(phi)/dr)
c  variable virial = (-r * d(phi)/dr)

c harmonic spring

      subroutine bond_harmonic(iflag)
      include "lammps.h"

      e_bond = 0.0

      do m = 1,nbondlocal

        i1 = bondlist(1,m)
        i2 = bondlist(2,m)
        itype = bondlist(3,m)

        delx = x(1,i1) - x(1,i2)
        dely = x(2,i1) - x(2,i2)
        delz = x(3,i1) - x(3,i2)
        call minimg(delx,dely,delz)

        rsq = delx*delx + dely*dely + delz*delz
        r = sqrt(rsq)

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 2
        endif

        dr = r - bondcoeff(2,itype)
        rk = bondcoeff(1,itype) * dr
        e_bond = e_bond + ifactor/2.0 * rk*dr

        if (r.gt.0.0) then
          force = -2.0*rk/r
        else
          force = 0.0
        endif

        if (newton_bond.eq.1.or.i1.le.maxown) then
          f(1,i1) = f(1,i1) + delx*force
          f(2,i1) = f(2,i1) + dely*force
          f(3,i1) = f(3,i1) + delz*force
        endif

        if (newton_bond.eq.1.or.i2.le.maxown) then
          f(1,i2) = f(1,i2) - delx*force
          f(2,i2) = f(2,i2) - dely*force
          f(3,i2) = f(3,i2) - delz*force
        endif

        if (iflag.eq.1) then
          virial(1) = virial(1) + ifactor/2.0*delx*delx*force
          virial(2) = virial(2) + ifactor/2.0*dely*dely*force
          virial(3) = virial(3) + ifactor/2.0*delz*delz*force
          virial(4) = virial(4) + ifactor/2.0*delx*dely*force
          virial(5) = virial(5) + ifactor/2.0*delx*delz*force
          virial(6) = virial(6) + ifactor/2.0*dely*delz*force
        endif

      enddo
      return
      end


c -------------------------------------------------------------------------
c FENE with repulsive LJ
c   if r -> r0, then rlogarg < 0.0 which is an error,
c     so issue a warning and reset rlogarg = epsilon
c   if r > 2*r0 something serious is wrong, just stop

      subroutine bond_fene_standard(iflag)
      include "lammps.h"

      e_bond = 0.0

      do m = 1,nbondlocal

        i1 = bondlist(1,m)
        i2 = bondlist(2,m)
        itype = bondlist(3,m)

        delx = x(1,i1) - x(1,i2)
        dely = x(2,i1) - x(2,i2)
        delz = x(3,i1) - x(3,i2)
        call minimg(delx,dely,delz)

        rsq = delx*delx + dely*dely + delz*delz

        r0sq = bondcoeff(2,itype) * bondcoeff(2,itype)
        rlogarg = 1.0 - rsq/r0sq

        if (rlogarg.lt.0.1) then
          write (6,*) 'WARNING: FENE bond too long',
     $         node,ntimestep,tag(i1),tag(i2),sqrt(rsq)
          if (rlogarg.le.-3.0) call exit(0)
          rlogarg = 0.1
        endif

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 2
        endif

        e_bond = e_bond - ifactor/4.0 *
     $       bondcoeff(1,itype)*r0sq*log(rlogarg)

        force = -bondcoeff(1,itype)/rlogarg

        if (rsq.lt.
     $       two_1_3*bondcoeff(4,itype)*bondcoeff(4,itype)) then
          sr2 = bondcoeff(4,itype)*bondcoeff(4,itype)/rsq
          sr6 = sr2*sr2*sr2
          e_bond = e_bond + ifactor/2.0 *
     $         (4.0*bondcoeff(3,itype)*sr6*(sr6-1.0) +
     $         bondcoeff(3,itype))
          force = force +
     $         48.0*bondcoeff(3,itype)*sr6*(sr6-0.5)/rsq
        endif

        if (newton_bond.eq.1.or.i1.le.maxown) then
          f(1,i1) = f(1,i1) + delx*force
          f(2,i1) = f(2,i1) + dely*force
          f(3,i1) = f(3,i1) + delz*force
        endif

        if (newton_bond.eq.1.or.i2.le.maxown) then
          f(1,i2) = f(1,i2) - delx*force
          f(2,i2) = f(2,i2) - dely*force
          f(3,i2) = f(3,i2) - delz*force
        endif

        if (iflag.eq.1) then
          virial(1) = virial(1) + ifactor/2.0*delx*delx*force
          virial(2) = virial(2) + ifactor/2.0*dely*dely*force
          virial(3) = virial(3) + ifactor/2.0*delz*delz*force
          virial(4) = virial(4) + ifactor/2.0*delx*dely*force
          virial(5) = virial(5) + ifactor/2.0*delx*delz*force
          virial(6) = virial(6) + ifactor/2.0*dely*delz*force
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c shifted FENE with LJ
c   if r -> r0, then rlogarg < 0.0 which is an error,
c     so issue a warning and resets rlogarg = epsilon
c   if r > 2*r0 something serious is wrong, just stop

      subroutine bond_fene_shift(iflag)
      include "lammps.h"

      e_bond = 0.0

      do m = 1,nbondlocal

        i1 = bondlist(1,m)
        i2 = bondlist(2,m)
        itype = bondlist(3,m)
        
        delx = x(1,i1) - x(1,i2)
        dely = x(2,i1) - x(2,i2)
        delz = x(3,i1) - x(3,i2)
        call minimg(delx,dely,delz)

        rsq = delx*delx + dely*dely + delz*delz

        r = sqrt(rsq)
        rshift = r - bondcoeff(5,itype)
        rshiftsq = rshift*rshift

        r0sq = bondcoeff(2,itype) * bondcoeff(2,itype)
        rlogarg = 1.0 - rshiftsq/r0sq

        if (rlogarg.lt.0.1) then
          write (6,*) 'WARNING: Bond too long',
     $         node,ntimestep,tag(i1),tag(i2),sqrt(rsq)
          if (rlogarg.le.-3.0) call exit(0)
          rlogarg = 0.1
        endif
        
        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 2
        endif

        e_bond = e_bond - ifactor/4.0 *
     $       bondcoeff(1,itype)*r0sq*log(rlogarg)

        force = -bondcoeff(1,itype)*rshift/rlogarg/r
        
        if (rshiftsq.lt.
     $       two_1_3*bondcoeff(4,itype)*bondcoeff(4,itype)) then
          sr2 = bondcoeff(4,itype)*bondcoeff(4,itype)/rshiftsq
          sr6 = sr2*sr2*sr2
          e_bond = e_bond + ifactor/2.0 *
     $         (4.0*bondcoeff(3,itype)*sr6*(sr6-1.0) +
     $         bondcoeff(3,itype))
          force = force +
     $         48.0*bondcoeff(3,itype)*sr6*(sr6-0.5)/rshift/r
        endif

        if (newton_bond.eq.1.or.i1.le.maxown) then
          f(1,i1) = f(1,i1) + delx*force
          f(2,i1) = f(2,i1) + dely*force
          f(3,i1) = f(3,i1) + delz*force
        endif

        if (newton_bond.eq.1.or.i2.le.maxown) then
          f(1,i2) = f(1,i2) - delx*force
          f(2,i2) = f(2,i2) - dely*force
          f(3,i2) = f(3,i2) - delz*force
        endif

        if (iflag.eq.1) then
          virial(1) = virial(1) + ifactor/2.0*delx*delx*force
          virial(2) = virial(2) + ifactor/2.0*dely*dely*force
          virial(3) = virial(3) + ifactor/2.0*delz*delz*force
          virial(4) = virial(4) + ifactor/2.0*delx*dely*force
          virial(5) = virial(5) + ifactor/2.0*delx*delz*force
          virial(6) = virial(6) + ifactor/2.0*dely*delz*force
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c non-linear spring (van Swol)

      subroutine bond_nonlinear(iflag)
      include "lammps.h"

      e_bond = 0.0

      do m = 1,nbondlocal

        i1 = bondlist(1,m)
        i2 = bondlist(2,m)
        itype = bondlist(3,m)

        delx = x(1,i1) - x(1,i2)
        dely = x(2,i1) - x(2,i2)
        delz = x(3,i1) - x(3,i2)
        call minimg(delx,dely,delz)

        rsq = delx*delx + dely*dely + delz*delz

        r = sqrt(rsq)
        delr = r - bondcoeff(2,itype)
        delrsq = delr*delr
        denom = bondcoeff(3,itype) - delrsq

        if (denom.le.0.0) then
          write (6,*) 'WARNING: Bond too short/long',
     $         node,ntimestep,tag(i1),tag(i2),r
          call exit(0)
        endif

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 2
        endif

        e_bond = e_bond + ifactor/2.0 * bondcoeff(1,itype)*delrsq /
     $       (bondcoeff(3,itype) - delrsq)

        denomsq = denom*denom
        force = -bondcoeff(1,itype)/r *
     $       2.0*delr*bondcoeff(3,itype)/denomsq

        if (newton_bond.eq.1.or.i1.le.maxown) then
          f(1,i1) = f(1,i1) + delx*force
          f(2,i1) = f(2,i1) + dely*force
          f(3,i1) = f(3,i1) + delz*force
        endif

        if (newton_bond.eq.1.or.i2.le.maxown) then
          f(1,i2) = f(1,i2) - delx*force
          f(2,i2) = f(2,i2) - dely*force
          f(3,i2) = f(3,i2) - delz*force
        endif

        if (iflag.eq.1) then
          virial(1) = virial(1) + ifactor/2.0*delx*delx*force
          virial(2) = virial(2) + ifactor/2.0*dely*dely*force
          virial(3) = virial(3) + ifactor/2.0*delz*delz*force
          virial(4) = virial(4) + ifactor/2.0*delx*dely*force
          virial(5) = virial(5) + ifactor/2.0*delx*delz*force
          virial(6) = virial(6) + ifactor/2.0*dely*delz*force
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c class 2 quartic bond
c class II force field
c written by Eric J. Simon, Cray Research
c assumes the array bondcoeff contains:
c       1       b0  (equilibrium bond length)
c       2       Kb2 (coefficient for quadratic term)
c       3       Kb3 (coefficient for cubic term)
c       4       Kb4 (coefficient for quartic term)

      subroutine bond_class2(iflag)
      include "lammps.h"
      parameter (zero=0.0,two=2.0,three=3.0,four=4.0)

      e_bond = 0.0
      
      do m = 1,nbondlocal

        i1 = bondlist(1,m)
        i2 = bondlist(2,m)
        itype = bondlist(3,m)

        delx = x(1,i1) - x(1,i2)
        dely = x(2,i1) - x(2,i2)
        delz = x(3,i1) - x(3,i2)
        call minimg(delx,dely,delz)

        rsq = delx*delx + dely*dely + delz*delz
        r = sqrt(rsq)

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.maxown) ifactor = ifactor + 1
          if (i2.le.maxown) ifactor = ifactor + 1
        else
          ifactor = 2
        endif

        if (r.gt.zero) then
          dr  = r - bondcoeff(1,itype)
          dr2 = dr  * dr
          dr3 = dr2 * dr
          dr4 = dr3 * dr

          e_bond = e_bond + (ifactor/2.0)*(bondcoeff(2,itype)*dr2 
     1         + bondcoeff(3,itype)*dr3
     2         + bondcoeff(4,itype)*dr4)

          de_bond    = two  *bondcoeff(2,itype)*dr  +
     1         three*bondcoeff(3,itype)*dr2 +
     2         four *bondcoeff(4,itype)*dr3

          force = -de_bond

          fx = force*delx/r
          fy = force*dely/r
          fz = force*delz/r

          if (newton_bond.eq.1.or.i1.le.maxown) then
            f(1,i1) = f(1,i1) + fx
            f(2,i1) = f(2,i1) + fy
            f(3,i1) = f(3,i1) + fz
          endif

          if (newton_bond.eq.1.or.i2.le.maxown) then
            f(1,i2) = f(1,i2) - fx
            f(2,i2) = f(2,i2) - fy
            f(3,i2) = f(3,i2) - fz
          endif

          if (iflag.eq.1) then
            virial(1) = virial(1) + (ifactor/2.0)*delx*fx
            virial(2) = virial(2) + (ifactor/2.0)*dely*fy
            virial(3) = virial(3) + (ifactor/2.0)*delz*fz
            virial(4) = virial(4) + (ifactor/2.0)*delx*fy
            virial(5) = virial(5) + (ifactor/2.0)*delx*fz
            virial(6) = virial(6) + (ifactor/2.0)*dely*fz
          endif

        endif

      enddo

      return
      end
