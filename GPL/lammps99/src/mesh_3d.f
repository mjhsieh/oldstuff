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
c create a 3-d mesh of processors for physical volume (xprd,yprd,zprd)

      subroutine mesh_3d(node,nprocs,idimension,xprd,yprd,zprd,
     $     pgrid,me,mpart)
      include "mpif.h"
      integer node,nprocs,idimension
      real*8 xprd,yprd,zprd
      integer pgrid(3),me(3),mpart(2,3)

      logical periods(3),reorder
      integer icommgrid,ierror
      integer mpi_grid(3),mpi_me(3)

c map procs into 3-d grid that most closely matches box volume

      if (pgrid(1).eq.0) call proc2box3d(nprocs,idimension,
     $     xprd,yprd,zprd,pgrid(1),pgrid(2),pgrid(3))
      
c determine where I am and who are my neighboring procs in 3-d proc grid
c order axes so as to agree with old version of code
c this ordering puts procs along x-axis first, then y, then z

      mpi_grid(1) = pgrid(3)
      mpi_grid(2) = pgrid(2)
      mpi_grid(3) = pgrid(1)

      reorder = .FALSE.
      periods(1) = .TRUE.
      periods(2) = .TRUE.
      periods(3) = .TRUE.
      
      call mpi_cart_create(mpi_comm_world,3,mpi_grid,periods,reorder,
     $     icommgrid,ierror)
      call mpi_cart_get(icommgrid,3,mpi_grid,periods,mpi_me,ierror)
      call mpi_cart_shift(icommgrid,2,1,mpart(1,1),mpart(2,1),ierror)
      call mpi_cart_shift(icommgrid,1,1,mpart(1,2),mpart(2,2),ierror)
      call mpi_cart_shift(icommgrid,0,1,mpart(1,3),mpart(2,3),ierror)
      call mpi_comm_free(icommgrid,ierror)

      me(1) = mpi_me(3)
      me(2) = mpi_me(2)
      me(3) = mpi_me(1)

      return
      end


c ----------------------------------------------------------------------
c assign nprocs to 3-d xprd,yprd,zprd box so as to minimize surface area
c if idimension = 2, enforce pz = 1
c return 3-d px,py,pz grid of procs

      subroutine proc2box3d(nprocs,idimension,xprd,yprd,zprd,px,py,pz)
      integer nprocs,idimension
      real*8 xprd,yprd,zprd
      integer px,py,pz

      real*8 surf,bestsurf,boxx,boxy,boxz
      integer ipx,ipy,ipz

      bestsurf = 2.0 * (xprd*yprd + yprd*zprd + zprd*xprd)

c loop thru all possible factorizations of nprocs
c surf = surface area of a proc sub-domain

      ipx = 1
      do while (ipx.le.nprocs)
        if (mod(nprocs,ipx).eq.0) then
          nremain = nprocs/ipx
          ipy = 1
          do while (ipy.le.nremain)
            if (mod(nremain,ipy).eq.0) then
              ipz = nremain/ipy
              boxx = xprd/ipx
              boxy = yprd/ipy
              boxz = zprd/ipz
              surf = boxx*boxy + boxy*boxz + boxz*boxx
              if ((surf.lt.bestsurf).and.
     $             (idimension.eq.3.or.
     $             (idimension.eq.2.and.ipz.eq.1))) then
                bestsurf = surf
                px = ipx
                py = ipy
                pz = ipz
              endif
            endif
            ipy = ipy + 1
          enddo
        endif
        ipx = ipx + 1
      enddo

      if (px*py*pz.ne.nprocs) call error('Bad result in proc2box3d')

      return
      end


c ----------------------------------------------------------------------
c assign nprocs to 2d nx,ny grid so as to minimize surface area
c return 2d px,py grid of procs

      subroutine proc2grid2d(nprocs,nx,ny,px,py)
      integer nprocs
      integer nx,ny
      integer px,py

      integer surf,boxx,boxy
      integer bestsurf,bestboxx,bestboxy
      integer ipx,ipy

      bestsurf = 2 * (nx + ny)
      bestboxx = 0
      bestboxy = 0

c loop thru all possible factorizations of nprocs
c surf = surface area of largest proc sub-domain
c innermost if test minimizes surface area and surface/volume ratio

      ipx = 1
      do while (ipx.le.nprocs)
        if (mod(nprocs,ipx).eq.0) then
          ipy = nprocs/ipx
          boxx = nx/ipx
          if (mod(nx,ipx).ne.0) boxx = boxx + 1
          boxy = ny/ipy
          if (mod(ny,ipy).ne.0) boxy = boxy + 1
          surf = boxx + boxy
          if (surf.lt.bestsurf.or.(surf.eq.bestsurf.and.
     $         boxx*boxy.gt.bestboxx*bestboxy)) then
            bestsurf = surf
            bestboxx = boxx
            bestboxy = boxy
            px = ipx
            py = ipy
          endif
        endif
        ipx = ipx + 1
      enddo

      if (px*py.ne.nprocs) call error('Bad result in proc2grid2d')

      return
      end
