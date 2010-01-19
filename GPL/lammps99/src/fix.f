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
c assign group of atoms to a constraint
c   fixgroup = 1 -> single atom
c   fixgroup = 2 -> molecule
c   fixgroup = 3 -> atom type
c   fixgroup = 4 -> spatial region
c   fixgroup = 5 -> all remaining unconstrained atoms

      subroutine set_fix
      include "lammps.h"

      if (node.eq.0) write (6,*) 'Assigning atoms to a fix ...'

      if (readflag.eq.0)
     $     call error('Data file must be read before assigning fix')

      if (fixgroup.eq.3) then
        if (fixtype.gt.ntypes)
     $       call error('Fixing too large an atom type')
      endif

      if (fixgroup.eq.1) then
        i = localptr(fixatom)
        if (i.ne.0.and.i.le.maxown) fix(i) = maxfix*fix(i) + fixwhich
      else if (fixgroup.eq.2) then
        i = atompnt
        do ii = 1,nlocal
          if (molecule(i).eq.fixatom)
     $         fix(i) = maxfix*fix(i) + fixwhich
          i = list(i)
        enddo
      else if (fixgroup.eq.3) then
        i = atompnt
        do ii = 1,nlocal
          if (type(i).eq.fixtype) fix(i) = maxfix*fix(i) + fixwhich
          i = list(i)
        enddo
      else if (fixgroup.eq.4) then
        i = atompnt
        do ii = 1,nlocal
          if (x(1,i).ge.fixregion(1).and.x(1,i).le.fixregion(2).and.
     $         x(2,i).ge.fixregion(3).and.x(2,i).le.fixregion(4).and.
     $         x(3,i).ge.fixregion(5).and.x(3,i).le.fixregion(6))
     $         fix(i) = maxfix*fix(i) + fixwhich
          i = list(i)
        enddo
      else if (fixgroup.eq.5) then
        i = atompnt
        do ii = 1,nlocal
          if (fix(i).eq.0) fix(i) = fixwhich
          i = list(i)
        enddo
      endif

      return
      end


c figure out how many global quantities each fix will need to accumulate
c called from start routine before a run
c fixptr(k) = 1st loc in fixstore where quantities associated with fix k
c             will be accumulated
c fixnum = total # of values to be accumulated

      subroutine fix_setup
      include "lammps.h"

      fixnum = 0
      do k = 1,nfixes
        fixptr(k) = fixnum + 1
        istyle = fixstyle(k)
        if (istyle.eq.1) then
          fixnum = fixnum + 5
        else if (istyle.eq.3) then
          fixnum = fixnum + 5
        else if (istyle.eq.4) then
          fixnum = fixnum + 2
        else if (istyle.eq.7) then
          fixnum = fixnum + 5
        endif
      enddo

      return
      end


c fix forces and velocities for each constrained atom at end of timestep
c  applies to all fix styles

      subroutine fix_apply
      include "lammps.h"
      include "mpif.h"

c zero out storage before accumulating info about each fix

      do n = 1,fixnum
        fixstore(n) = 0.0
      enddo

c accumulate each processor's contribution to quantities
c  associated with each fix

      i = atompnt
      do ii = 1,nlocal
        itmp = fix(i)
 10     if (itmp.gt.0) then
          iwhich = mod(itmp,maxfix)
          n = fixptr(iwhich)
          istyle = fixstyle(iwhich)
          if (istyle.eq.1) then
            fixstore(n) = fixstore(n) + 1.0
            fixstore(n+1) = fixstore(n+1) + f(1,i)
            fixstore(n+2) = fixstore(n+2) + f(2,i)
            fixstore(n+3) = fixstore(n+3) + f(3,i)
            fixstore(n+4) = fixstore(n+4) + mass(type(i))
          else if (istyle.eq.3) then
            fixstore(n) = fixstore(n) + 1.0
            fixstore(n+1) = fixstore(n+1) + f(1,i)
            fixstore(n+2) = fixstore(n+2) + f(2,i)
            fixstore(n+3) = fixstore(n+3) + f(3,i)
            fixstore(n+4) = fixstore(n+4) + mass(type(i))
          else if (istyle.eq.4) then
            if (mod(itime,nint(fixcoeff(3,iwhich))).eq.0) then
              fixstore(n) = fixstore(n) + 1.0
              fixstore(n+1) = fixstore(n+1) +
     $             (v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)) *
     $             mass(type(i))
            endif
          else if (istyle.eq.7) then
            fixstore(n) = fixstore(n) + 1.0
            iboxx = mod(true(i),1000) - 500
            iboxy = mod(true(i)/1000,1000) - 500
            iboxz = true(i)/1000000 - 500
            rmass = mass(type(i))
            fixstore(n+1) = fixstore(n+1) + (x(1,i)+iboxx*xprd)*rmass
            fixstore(n+2) = fixstore(n+2) + (x(2,i)+iboxy*yprd)*rmass
            fixstore(n+3) = fixstore(n+3) + (x(3,i)+iboxz*zprd)*rmass
            fixstore(n+4) = fixstore(n+4) + rmass
          endif
          itmp = itmp/maxfix
          goto 10
        endif
        i = list(i)
      enddo

c do a global sum of these local quantities across all processors
c  use buf4 as work array

      if (fixnum.gt.0) then
        call mpi_allreduce(fixstore,buf4,fixnum,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,fixnum
          fixstore(i) = buf4(i)
        enddo
      endif

c normalize global quantities as needed and setup fix coeffs

      do k = 1,nfixes
        n = fixptr(k)
        istyle = fixstyle(k)
        if (istyle.eq.1) then
          ncount = nint(fixstore(n))
          fixstore(n+1) = fixstore(n+1)/ncount
          fixstore(n+2) = fixstore(n+2)/ncount
          fixstore(n+3) = fixstore(n+3)/ncount
          fixstore(n+4) = fixstore(n+4)/ncount
        else if (istyle.eq.3) then
          ncount = nint(fixstore(n))
          fixstore(n+1) = fixstore(n+1)/ncount
          fixstore(n+2) = fixstore(n+2)/ncount
          fixstore(n+3) = fixstore(n+3)/ncount
          fixstore(n+4) = fixstore(n+4)/ncount
          fixcoeff(4,k) = fixstore(n+1) + fixcoeff(1,k)
          fixcoeff(5,k) = fixstore(n+2) + fixcoeff(2,k)
          fixcoeff(6,k) = fixstore(n+3) + fixcoeff(3,k)
        else if (istyle.eq.4) then
          fixcoeff(5,k) = -1.0
          if (mod(itime,nint(fixcoeff(3,k))).eq.0) then
            ncount = nint(fixstore(n))
            fixstore(n+1) = fixstore(n+1) / (boltz*idimension*ncount)
            target = fixcoeff(1,k) +
     $           float(itime)/nsteps*(fixcoeff(2,k)-fixcoeff(1,k))
            if (abs(target-fixstore(n+1)).gt.fixcoeff(4,k).and.
     $           fixstore(n+1).gt.0.0)
     $           fixcoeff(5,k) = sqrt(target/fixstore(n+1))
          endif
        else if (istyle.eq.5) then
          target = fixcoeff(1,k) +
     $         float(itime)/nsteps*(fixcoeff(2,k)-fixcoeff(1,k))
          fixcoeff(4,k) = sqrt(24.0*boltz*target*fixcoeff(3,k)/dt)
        else if (istyle.eq.7) then
          fixstore(n+1) = fixstore(n+1)/fixstore(n+4)
          fixstore(n+2) = fixstore(n+2)/fixstore(n+4)
          fixstore(n+3) = fixstore(n+3)/fixstore(n+4)
          itrue = 500500500
          call remap(fixstore(n+1),fixstore(n+2),fixstore(n+3),itrue)
          dx = fixstore(n+1) - fixcoeff(1,k)
          dy = fixstore(n+2) - fixcoeff(2,k)
          dz = fixstore(n+3) - fixcoeff(3,k)
          if (fixflag(1,k).eq.1) dx = 0.0
          if (fixflag(2,k).eq.1) dy = 0.0
          if (fixflag(3,k).eq.1) dz = 0.0
          call minimg2(dx,dy,dz)
          fixcoeff(5,k) = dx*fixcoeff(4,k)
          fixcoeff(6,k) = dy*fixcoeff(4,k)
          fixcoeff(7,k) = dz*fixcoeff(4,k)
          fixcoeff(8,k) = fixstore(n+4)
        endif
      enddo

c apply constraints to each atom by adjusting force or velocity
c looping over all fixes on that atom in reverse of assignment order

      i = atompnt
      do ii = 1,nlocal
        itmp = fix(i)
 20     if (itmp.gt.0) then
          iwhich = mod(itmp,maxfix)
          istyle = fixstyle(iwhich)
          if (istyle.eq.1) then
            if (fixflag(1,iwhich).eq.0) f(1,i) = fixcoeff(1,iwhich)
            if (fixflag(2,iwhich).eq.0) f(2,i) = fixcoeff(2,iwhich)
            if (fixflag(3,iwhich).eq.0) f(3,i) = fixcoeff(3,iwhich)
          else if (istyle.eq.2) then
            f(1,i) = f(1,i) + fixcoeff(1,iwhich)
            f(2,i) = f(2,i) + fixcoeff(2,iwhich)
            f(3,i) = f(3,i) + fixcoeff(3,iwhich)
          else if (istyle.eq.3) then
            f(1,i) = fixcoeff(4,iwhich)
            f(2,i) = fixcoeff(5,iwhich)
            f(3,i) = fixcoeff(6,iwhich)
          else if (istyle.eq.4) then
            if (fixcoeff(5,iwhich).ge.0.0) then
              v(1,i) = v(1,i)*fixcoeff(5,iwhich)
              v(2,i) = v(2,i)*fixcoeff(5,iwhich)
              v(3,i) = v(3,i)*fixcoeff(5,iwhich)
            endif
          else if (istyle.eq.5) then
            rmass = mass(type(i))
            gamma1 = -fixcoeff(3,iwhich)*rmass
            gamma2 = fixcoeff(4,iwhich)*sqrt(rmass)
            if (fixflag(1,iwhich).eq.1)
     $           f(1,i) = f(1,i) + gamma1*v(1,i) +
     $           gamma2*(ranmars(0)-0.5)
            if (fixflag(2,iwhich).eq.1)
     $           f(2,i) = f(2,i) + gamma1*v(2,i) +
     $           gamma2*(ranmars(0)-0.5)
            if (fixflag(3,iwhich).eq.1)
     $           f(3,i) = f(3,i) + gamma1*v(3,i) +
     $           gamma2*(ranmars(0)-0.5)
          else if (istyle.eq.6) then
            continue
          else if (istyle.eq.7) then
            massfrac = mass(type(i))/fixcoeff(8,iwhich)
            if (fixflag(1,iwhich).eq.0)
     $           f(1,i) = f(1,i) - massfrac*fixcoeff(5,iwhich)
            if (fixflag(2,iwhich).eq.0)
     $           f(2,i) = f(2,i) - massfrac*fixcoeff(6,iwhich)
            if (fixflag(3,iwhich).eq.0)
     $           f(3,i) = f(3,i) - massfrac*fixcoeff(7,iwhich)
          else if (istyle.eq.8) then
            dx = x(1,i) - fixcoeff(1,iwhich)
            dy = x(2,i) - fixcoeff(2,iwhich)
            dz = x(3,i) - fixcoeff(3,iwhich)
            if (fixflag(1,iwhich).eq.1) dx = 0.0
            if (fixflag(2,iwhich).eq.1) dy = 0.0
            if (fixflag(3,iwhich).eq.1) dz = 0.0
            call minimg2(dx,dy,dz)
            r = sqrt(dx*dx + dy*dy + dz*dz)
            if (r.gt.fixcoeff(5,iwhich)) then
              pre = fixcoeff(4,iwhich)/r
              f(1,i) = f(1,i) - pre*dx
              f(2,i) = f(2,i) - pre*dy
              f(3,i) = f(3,i) - pre*dz
            endif
          endif
          itmp = itmp/maxfix
          goto 20
        endif
        i = list(i)
      enddo

      return
      end
