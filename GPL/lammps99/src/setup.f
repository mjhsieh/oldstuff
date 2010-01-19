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
c find box-bounds with non-PBC

      subroutine setup_bounds
      include "lammps.h"
      include "mpif.h"

      parameter (small=1.0E-4)
      real*8 tmp1(3),tmp2(3),tmp1all(3),tmp2all(3)
      
      tmp1(1) = 1.0E20
      tmp1(2) = 1.0E20
      tmp1(3) = 1.0E20
      tmp2(1) = -1.0E20
      tmp2(2) = -1.0E20
      tmp2(3) = -1.0E20

      i = atompnt
      do ii = 1,nlocal
        if (x(1,i).lt.tmp1(1)) tmp1(1) = x(1,i)
        if (x(1,i).gt.tmp2(1)) tmp2(1) = x(1,i)
        if (x(2,i).lt.tmp1(2)) tmp1(2) = x(2,i)
        if (x(2,i).gt.tmp2(2)) tmp2(2) = x(2,i)
        if (x(3,i).lt.tmp1(3)) tmp1(3) = x(3,i)
        if (x(3,i).gt.tmp2(3)) tmp2(3) = x(3,i)
        i = list(i)
      enddo

      call mpi_allreduce(tmp1,tmp1all,3,mpi_double_precision,
     $     mpi_min,mpi_comm_world,ierror)
      call mpi_allreduce(tmp2,tmp2all,3,mpi_double_precision,
     $     mpi_max,mpi_comm_world,ierror)

c reset box bounds in each non-periodic dimension
c  adding small -> prevents exchange routine from trying to move an atom
c   exactly on a sub-domain boundary (xboundhi) to a non-existent processor

      if (perflagx.eq.1) then
        xboundlo = tmp1all(1) - small
        xboundhi = tmp2all(1) + small
      endif
      if (perflagy.eq.1) then
        yboundlo = tmp1all(2) - small
        yboundhi = tmp2all(2) + small
      endif
      if (perflagz.eq.1) then
        zboundlo = tmp1all(3) - small
        zboundhi = tmp2all(3) + small
      endif

      return
      end


c define box-dependent values and sub-domain owned by each processor
      
      subroutine setup_box
      include "lammps.h"
      
c find atom bounds
c  only necessary for non-periodic BC (not for constant P)

      if (perflagx+perflagy+perflagz.gt.0) call setup_bounds

c define box

      xprd = xboundhi - xboundlo
      yprd = yboundhi - yboundlo
      zprd = zboundhi - zboundlo

c set PBC test distances

      xmc = xprd - cutforce
      xms = xprd - cutneigh
      ymc = yprd - cutforce
      yms = yprd - cutneigh
      if (idimension.eq.2) then
        zmc = 1.0
        zms = 1.0
      else
        zmc = zprd - cutforce
        zms = zprd - cutneigh
      endif

c error check on box size

      if (perflagx.eq.0.and.cutneigh.ge.xprd/2.0)
     $     call error('Neighbor list cutoff >= 1/2 box size')
      if (perflagy.eq.0.and.cutneigh.ge.yprd/2.0)
     $     call error('Neighbor list cutoff >= 1/2 box size')
      if (idimension.eq.3.and.perflagz.eq.0.and.cutneigh.ge.zprd/2.0)
     $     call error('Neighbor list cutoff >= 1/2 box size')

c set sub-domain boundaries

      border(1,1) = xboundlo + float(me(1))/pgrid(1) * xprd
      border(2,1) = xboundlo + float(me(1)+1)/pgrid(1) * xprd
      border(1,2) = yboundlo + float(me(2))/pgrid(2) * yprd
      border(2,2) = yboundlo + float(me(2)+1)/pgrid(2) * yprd
      border(1,3) = zboundlo + float(me(3))/pgrid(3) * zprd
      border(2,3) = zboundlo + float(me(3)+1)/pgrid(3) * zprd

c this forces setup_comm and setup_neigh to be (re)called

      boxflag = 1

      return
      end


c setup inter-processor communication

      subroutine setup_comm
      include "lammps.h"

      real*8 prd(3),prdlo(3),prdhi(3)

c find how many boxes I need in each dimension

      need(1) = cutneigh / (xprd/pgrid(1)) + 1
      need(2) = cutneigh / (yprd/pgrid(2)) + 1
      need(3) = cutneigh / (zprd/pgrid(3)) + 1

c don't exchange more than 1/2 way around (e.g. 2 boxes away when pgrid = 5)

      need(1) = min(need(1),pgrid(1)/2)
      need(2) = min(need(2),pgrid(2)/2)
      need(3) = min(need(3),pgrid(3)/2)
      
c setup 4 parameters for each exchange: (spart,rpart,boundlo,boundhi)
c  (1,2) nodes to swap with
c  (3,4) slab boundaries (in correct dimension) of atoms that will be sent
c 1st part of if is sending to the west (south,down)
c 2nd part of if is sending to the east (north,up)
c nbox = box (in this dimension) who originally owned the atoms 
c   I will be sending in this swap
c kk/2+1 vs pgrid/2 test is to make sure a node doesn't get more than 1/2
c   of a box from both directions (else would get some atoms twice)

      prd(1) = xprd
      prd(2) = yprd
      prd(3) = zprd
      prdlo(1) = xboundlo
      prdlo(2) = yboundlo
      prdlo(3) = zboundlo
      prdhi(1) = xboundhi
      prdhi(2) = yboundhi
      prdhi(3) = zboundhi
      nswap = 0
      
      do k = 1,3
        do kk = 0,2*need(k)-1
          nswap = nswap + 1
          if (nswap.le.maxswap) then
            if (mod(kk,2).eq.0) then
              spart(nswap) = mpart(1,k)
              rpart(nswap) = mpart(2,k)
              nbox = me(k) + kk/2
              if (nbox.ge.pgrid(k)) nbox = nbox - pgrid(k)
              blo = prdlo(k) + float(nbox)/pgrid(k) * prd(k)
              bhi = prdlo(k) + float(nbox+1)/pgrid(k) * prd(k)
              if (kk/2+1.eq.need(k)) then
                if (cutneigh.lt.need(k)*(bhi-blo)) then
                  bhi = border(1,k) + cutneigh
                  if (bhi.ge.prdhi(k)) bhi = bhi - prd(k)
                endif
                if (kk/2+1.eq.(pgrid(k)+1)/2) then
                  btmp = prdlo(k) + float(nbox+1)/pgrid(k) * prd(k)
                  bmid = (blo+btmp) / 2.0
                  bhi = min(bhi,bmid)
                endif
              endif
            else
              spart(nswap) = mpart(2,k)
              rpart(nswap) = mpart(1,k)
              nbox = me(k) - kk/2
              if (nbox.lt.0) nbox = nbox + pgrid(k)
              blo = prdlo(k) + float(nbox)/pgrid(k) * prd(k)
              bhi = prdlo(k) + float(nbox+1)/pgrid(k) * prd(k)
              if (kk/2+1.eq.need(k)) then
                if (cutneigh.lt.need(k)*(bhi-blo)) then
                  blo = border(2,k) - cutneigh
                  if (blo.lt.prdlo(k)) blo = blo + prd(k)
                endif
                if (kk/2+1.eq.(pgrid(k)+1)/2) then
                  btmp = prdlo(k) + float(nbox)/pgrid(k) * prd(k)
                  bmid = (btmp+bhi) / 2.0
                  blo = max(blo,bmid)
                endif
              endif
            endif
            boundlo(nswap) = blo
            boundhi(nswap) = bhi
          endif
        enddo
      enddo
      
      if (nswap.gt.maxswap)
     $     call error('Too many swaps - boost maxswap')

      return
      end


c setup neighbor binning parameters in 2-d/3-d box owned by each processor
c  adding small -> bins slightly larger
c  prevents round-off error so that bin # always < nbin_xyz when atoms binned
      
      subroutine setup_neigh
      include "lammps.h"
      include "mpif.h"

      parameter (small=1.0E-4)

c number and size of bins in each dimension

      nbinx = xprd / cutneigh
      binsizex = (xprd+small*xprd) / nbinx
      nbiny = yprd / cutneigh
      binsizey = (yprd+small*yprd) / nbiny

      if (idimension.eq.2) then
        nbinz = 1
        binsizez = 1.0
      else
        nbinz = zprd / cutneigh
        binsizez = (zprd+small*zprd) / nbinz
      endif

c error check for too small a domain

      if (nbinx.le.2.or.nbiny.le.2)
     $     call error('Two or less neighbor bins in a dimension')

      if (idimension.eq.3.and.nbinz.le.2)
     $     call error('Two or less neighbor bins in a dimension')

c compute local bin extent
c  subtract/add 1 to include nearby atoms
      
      mbinxlo = int((border(1,1) - xboundlo) / binsizex) - 1
      ix = int((border(2,1) - xboundlo) / binsizex) + 1
      if (ix.gt.nbinx) ix = nbinx
      mbinx = ix - mbinxlo + 1
      mbinx = min(mbinx,nbinx)
      if (mbinxlo.lt.0) mbinxlo = nbinx - 1

      mbinylo = int((border(1,2) - yboundlo) / binsizey) - 1
      iy = int((border(2,2) - yboundlo) / binsizey) + 1
      if (iy.gt.nbiny) iy = nbiny
      mbiny = iy - mbinylo + 1
      mbiny = min(mbiny,nbiny)
      if (mbinylo.lt.0) mbinylo = nbiny - 1

      mbinzlo = int((border(1,3) - zboundlo) / binsizez) - 1
      iz = int((border(2,3) - zboundlo) / binsizez) + 1
      if (iz.gt.nbinz) iz = nbinz
      mbinz = iz - mbinzlo + 1
      mbinz = min(mbinz,nbinz)
      if (mbinzlo.lt.0) mbinzlo = nbinz - 1

c error check

      iflag = mbinx*mbiny*mbinz
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (jflag.gt.maxbin)
     $     call error('Too many local bins - boost maxbin')
      
      return
      end
