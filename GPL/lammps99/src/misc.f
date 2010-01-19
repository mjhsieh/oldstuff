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
c copy force and virial arrays

      subroutine copy_force(f_target,vir_target)
      include "lammps.h"
      dimension f_target(3,*),vir_target(6)

      do i = 1,6
        vir_target(i) = virial(i)
      enddo

      i = atompnt
      do ii = 1,nlocal
        f_target(1,i) = f(1,i)
        f_target(2,i) = f(2,i)
        f_target(3,i) = f(3,i)
        i = list(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c zero out force and virial arrays

      subroutine zero_force
      include "lammps.h"

      do i = 1,6
        virial(i) = 0.0
      enddo

      if (newton.eq.0) then
        i = atompnt
        do ii = 1,nlocal
          f(1,i) = 0.0
          f(2,i) = 0.0
          f(3,i) = 0.0
          i = list(i)
        enddo
      else
        i = atompnt
        do ii = 1,nlocal+nother
          f(1,i) = 0.0
          f(2,i) = 0.0
          f(3,i) = 0.0
          if (ii.le.nlocal) then
            i = list(i)
          else
            i = i + 1
          endif
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c minimum image convention with neighbor cutoff

      subroutine minimg(delx,dely,delz)
      include "lammps.h"

      if (abs(delx).gt.xms.and.perflagx.eq.0) then
        if (delx.lt.0.0) then
          delx = delx + xprd
        else
          delx = delx - xprd
        endif
      endif
      if (abs(dely).gt.yms.and.perflagy.eq.0) then
        if (dely.lt.0.0) then
          dely = dely + yprd
        else
          dely = dely - yprd
        endif
      endif
      if (abs(delz).gt.zms.and.perflagz.eq.0) then
        if (delz.lt.0.0) then
          delz = delz + zprd
        else
          delz = delz - zprd
        endif
      endif
      
      return
      end


c -------------------------------------------------------------------------
c minimum image convention with 1/2 the box length

      subroutine minimg2(delx,dely,delz)
      include "lammps.h"

      if (abs(delx).gt.0.5*xprd.and.perflagx.eq.0) then
        if (delx.lt.0.0) then
          delx = delx + xprd
        else
          delx = delx - xprd
        endif
      endif
      if (abs(dely).gt.0.5*yprd.and.perflagy.eq.0) then
        if (dely.lt.0.0) then
          dely = dely + yprd
        else
          dely = dely - yprd
        endif
      endif
      if (abs(delz).gt.0.5*zprd.and.perflagz.eq.0) then
        if (delz.lt.0.0) then
          delz = delz + zprd
        else
          delz = delz - zprd
        endif
      endif
      
      return
      end


c -------------------------------------------------------------------------
c remap the point (xx,yy,zz) into the periodic box
c  no matter how far away it is
c only do it if periodicity is on

      subroutine remap(xx,yy,zz,itrue)
      include "lammps.h"

      if (perflagx.eq.0) then
 10     if (xx.ge.xboundhi) then
          xx = xx - xprd
          itrue = itrue + 1
          goto 10
        endif
 20     if (xx.lt.xboundlo) then
          xx = xx + xprd
          itrue = itrue - 1
          goto 20
        endif
      endif

      if (perflagy.eq.0) then
 30     if (yy.ge.yboundhi) then
          yy = yy - yprd
          itrue = itrue + 1000
          goto 30
        endif
 40     if (yy.lt.yboundlo) then
          yy = yy + yprd
          itrue = itrue - 1000
          goto 40
        endif
      endif

      if (perflagz.eq.0) then
 50     if (zz.ge.zboundhi) then
          zz = zz - zprd
          itrue = itrue + 1000000
          goto 50
        endif
 60     if (zz.lt.zboundlo) then
          zz = zz + zprd
          itrue = itrue - 1000000
          goto 60
        endif
      endif

      return
      end


c -------------------------------------------------------------------------
c error handler routine - close log file when exit

      subroutine error(str)
      include "lammps.h"
      include "mpif.h"
      character*(*) str

c close all files

      if (dumpatomfileflag.gt.0) call dumpatomclose
      if (dumpvelfileflag.gt.0) call dumpvelclose
      if (dumpforcefileflag.gt.0) call dumpforceclose
      do idiag = 1,numdiag
        if (diagfileflag(idiag).gt.0) call diag_close(idiag)
      enddo
      if (node.eq.0) close (1)

c final timer

      call mpi_barrier(mpi_comm_world,ierror)
      time = mpi_wtime()
      time_total = time-time_total

      if (str(1:8).eq.'All Done') then
        if (node.eq.0) write (6,*) 'Total time =',time_total
      else
        if (node.eq.0) write (6,*) 'Error: ',str
      endif

c close down MPI

      call mpi_finalize(ierror)

      call exit(0)

      return
      end
