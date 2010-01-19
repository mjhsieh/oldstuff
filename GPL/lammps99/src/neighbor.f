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
c check each timestep to see if neighbor lists should be constructed
c  return nflag=0 for no, nflag=1 for yes
c  for neightrigger=1, test max distance moved against skin distance
c  if this is step where restart file will be written, force a reneighboring
c   to insure perfect restarts

      subroutine check_neighbor(nflag)
      include "lammps.h"
      include "mpif.h"
      integer nflag

      nflag = 0

      if (nrestart_next.eq.ntimestep) then
        nflag = 1
        return
      endif

      neighago = neighago + 1
      if (neighago.lt.neighdelay) return
      if (mod(neighago,neighfreq).ne.0) return

      if (neightrigger.eq.0) then
        nflag = 1
      else
        i = atompnt
        do ii = 1,nlocal
          delx = x(1,i) - xhold(1,i)
          dely = x(2,i) - xhold(2,i)
          delz = x(3,i) - xhold(3,i)
          call minimg(delx,dely,delz)
          rsq = delx*delx + dely*dely + delz*delz
          if (rsq.gt.triggersq) nflag = 1
          i = list(i)
        enddo
        call mpi_allreduce(nflag,nflagall,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (nflagall.gt.0) then
          nflag = 1
          if (neighago.eq.max(neighfreq,neighdelay))
     $         ndanger = ndanger + 1
        endif
      endif

      return
      end


c driver for neighbor-list construction routines
c  store current positions (neightrigger=1)
c  call appropriate routine (neighstyle=0,1)
c  store new maximums

      subroutine neighbor
      include "lammps.h"

      numneigh = numneigh + 1
      neighago = 0

      if (neightrigger.eq.1) then
        i = atompnt
        do ii = 1,nlocal
          xhold(1,i) = x(1,i)
          xhold(2,i) = x(2,i)
          xhold(3,i) = x(3,i)
          i = list(i)
        enddo
      endif
      
      if (neighstyle.eq.0) then
        call neighbor_nsq
      else
        if (boxflag.eq.1) call setup_neigh
        call neighbor_bin
      endif
      
      boxflag = 0
      max_nlocal = max(max_nlocal,nlocal)
      max_nother = max(max_nother,nother)
      max_slist = max(max_slist,nslist(nswap+1)-1)
      max_neigh = max(max_neigh,neightop)

      return
      end
      
      
c 2-d/3-d neighbor list construction with no binning, Newton's 3rd law
c  N^2 / 2 search for neighbor pairs in my box
c  pair stored once in neighbor list if atoms i,j are both in my box
c   (if atom i comes before j in list vector)
c  pair stored by me if j is NOT in my box and tag(j)-tag(i) < N/2 in
c   periodic sense (see goto 30 statements)
      
      subroutine neighbor_nsq
      include "lammps.h"
      
c      nhalf = natoms/2
      npnt = 0
      i = atompnt
      do ii = 1,nlocal
        nliststart(i) = npnt + 1
        itag = tag(i)
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        num1tmp = num1bond(i)
        num2tmp = num2bond(i)
        num3tmp = num3bond(i)
        j = list(i)
        do jj = ii+1,nlocal+nother
          if (newton_nonbond.eq.1.and.j.gt.maxown) then
            jtag = tag(j)
            if (itag.gt.jtag) then
              if (mod(itag+jtag,2).eq.0) goto 30
            else
              if (mod(itag+jtag,2).eq.1) goto 30
            endif
          endif
c          if (newton_nonbond.eq.1.and.j.gt.maxown) then
c            jtag = tag(j)
c            if (jtag.gt.itag) then
c              if (jtag-itag.gt.nhalf) goto 30
c            else
c              if (itag-jtag.le.nhalf) goto 30
c            endif
c          endif
          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
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
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)
          if (rsq.le.cutneighsq(itype,jtype)) then
            jtag = tag(j)
            iflag = 0
            jpnt = 1
 20         if (iflag.eq.0.and.jpnt.le.num3tmp) then
              if (specbond(jpnt,i).eq.jtag) then
                iflag = 1
                if (jpnt.le.num1tmp) then
                  iwhich = 1
                else if (jpnt.le.num2tmp) then
                  iwhich = 2
                else
                  iwhich = 3
                endif
              else
                jpnt = jpnt + 1
              endif
              goto 20
            else if (iflag.eq.0) then
              npnt = npnt + 1
              nlist(npnt) = j
            else if (special(iwhich).eq.0.0) then
              continue
            else if (special(iwhich).eq.1.0) then
              npnt = npnt + 1
              nlist(npnt) = j
            else
              npnt = npnt + 1
              nlist(npnt) = iwhich*maxatomp1 + j-1
            endif
          endif
 30       if (jj.le.nlocal) then
            j = list(j)
          else
            j = j + 1
          endif
        enddo
        if (npnt.gt.maxown*maxneigh) then
          write (6,*) 'Neighbor list too big:',node,npnt,ii,nlocal
          call exit(0)
        endif
        nliststop(i) = npnt
        i = list(i)
      enddo
      
      neightop = npnt

      return
      end
      
      

c neighbor list construction with binning, Newton's 3rd law
c  all of mine and nearby atoms binned once
c  each owned atom i checks 27 (3d) or 9 (2d) surrounding boxes
c  pair stored once in list if atoms i,j are both owned and tag(i) < tag(j)
c  pair stored by me if j is NOT in my box and tag(j)-tag(i) < N/2 in
c   periodic sense (see goto 30 statements)
      
      subroutine neighbor_bin
      include "lammps.h"
      
      do i = 1,mbinx*mbiny*mbinz
        binpnt(i) = 0
      enddo

      i = atompnt
      do ii = 1,nlocal+nother
        ix = (x(1,i) - xboundlo) / binsizex
        iy = (x(2,i) - yboundlo) / binsizey
        iz = (x(3,i) - zboundlo) / binsizez
        ix = ix - mbinxlo
        if (ix.lt.0) ix = ix + nbinx
        iy = iy - mbinylo
        if (iy.lt.0) iy = iy + nbiny
        iz = iz - mbinzlo
        if (iz.lt.0) iz = iz + nbinz
        if (ix.le.mbinx.and.iy.le.mbiny.and.iz.le.mbinz) then
          ib = iz*mbiny*mbinx + iy*mbinx + ix + 1
          bin(i) = binpnt(ib)
          binpnt(ib) = i
        endif
        if (ii.le.nlocal) then
          i = list(i)
        else
          i = i + 1
        endif
      enddo
      
      iboxneigh = 26
      if (idimension.eq.2) iboxneigh = 8
c      nhalf = natoms/2
      npnt = 0
      i = atompnt
      do ii = 1,nlocal
        nliststart(i) = npnt + 1
        itag = tag(i)
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        ixx = (xtmp - xboundlo) / binsizex
        iyy = (ytmp - yboundlo) / binsizey
        izz = (ztmp - zboundlo) / binsizez
        num1tmp = num1bond(i)
        num2tmp = num2bond(i)
        num3tmp = num3bond(i)
        do k = 0,iboxneigh
          ix = ixx + mod(k,3) - 1
          iy = iyy + mod(k/3,3) - 1
          iz = izz + k/9 - 1
          ix = ix - mbinxlo
          if (ix.lt.0) ix = ix + nbinx
          if (ix.eq.nbinx) ix = 0
          iy = iy - mbinylo
          if (iy.lt.0) iy = iy + nbiny
          if (iy.eq.nbiny) iy = 0
          iz = iz - mbinzlo
          if (iz.lt.0) iz = iz + nbinz
          if (iz.eq.nbinz) iz = 0
          ib = iz*mbiny*mbinx + iy*mbinx + ix + 1
          j = binpnt(ib)
 10       if (j.ne.0) then
c            if (j.le.i) goto 30
            if (j.le.maxown) then
              if (tag(j).le.itag) goto 30
            else if (newton_nonbond.eq.1) then
              jtag = tag(j)
              if (itag.gt.jtag) then
                if (mod(itag+jtag,2).eq.0) goto 30
              else
                if (mod(itag+jtag,2).eq.1) goto 30
              endif
            endif
c            if (newton_nonbond.eq.1.and.j.gt.maxown) then
c              jtag = tag(j)
c              if (jtag.gt.itag) then
c                if (jtag-itag.gt.nhalf) goto 30
c              else
c                if (itag-jtag.le.nhalf) goto 30
c              endif
c            endif
            delx = xtmp - x(1,j)
            dely = ytmp - x(2,j)
            delz = ztmp - x(3,j)
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
            rsq = delx*delx + dely*dely + delz*delz
            jtype = type(j)
            if (rsq.le.cutneighsq(itype,jtype)) then
              jtag = tag(j)
              iflag = 0
              jpnt = 1
 20           if (iflag.eq.0.and.jpnt.le.num3tmp) then
                if (specbond(jpnt,i).eq.jtag) then
                  iflag = 1
                  if (jpnt.le.num1tmp) then
                    iwhich = 1
                  else if (jpnt.le.num2tmp) then
                    iwhich = 2
                  else
                    iwhich = 3
                  endif
                else
                  jpnt = jpnt + 1
                endif
                goto 20
              else if (iflag.eq.0) then
                npnt = npnt + 1
                nlist(npnt) = j
              else if (special(iwhich).eq.0.0) then
                continue
              else if (special(iwhich).eq.1.0) then
                npnt = npnt + 1
                nlist(npnt) = j
              else
                npnt = npnt + 1
                nlist(npnt) = iwhich*maxatomp1 + j-1
              endif
            endif
 30         j = bin(j)
            goto 10
          endif
        enddo
        if (npnt.gt.maxown*maxneigh) then
          write (6,*) 'Neighbor list too big:',node,npnt,ii,nlocal
          call exit(0)
        endif
        nliststop(i) = npnt
        i = list(i)
      enddo
      
      neightop = npnt

      return
      end

      
c neighbor list construction of bonds
c  for each bond store 2 local atoms in bond and ptr to bond parameters 
c  j12 = 0 check is to see if atom in bond is not in own or other list
c  if so is an error

      subroutine neighbor_bond
      include "lammps.h"

      nbondlocal = 0
      i = atompnt
      do ii = 1,nlocal
        do k = 1,numbond(i)
          j1 = localptr(bondatom1(k,i))
          j2 = localptr(bondatom2(k,i))
          if (j1.eq.0.or.j2.eq.0) then
            write (6,*) 'Bond partner missing:',node,
     $           bondatom1(k,i),bondatom2(k,i),numbond(i),ntimestep
            call exit(0)
          endif
          if (newton_bond.eq.1.or.(i.le.j1.and.i.le.j2)) then
            nbondlocal = nbondlocal + 1
            if (nbondlocal.le.maxbond) then
              bondlist(1,nbondlocal) = j1
              bondlist(2,nbondlocal) = j2
              bondlist(3,nbondlocal) = bondtype(k,i)
            endif
          endif
        enddo
        i = list(i)
      enddo

c check for too many bonds

      if (nbondlocal.gt.maxbond) then
        write (6,*) 'Bond list too big:',node,nbondlocal
        call exit(0)
      endif

      max_bond = max(max_bond,nbondlocal)

      return
      end


c neighbor list construction of angles
c  for each angle store 3 local atoms in angle and ptr to angle parameters 
c  j123 = 0 check is to see if atom in angle is not in own or other list
c  if so is an error

      subroutine neighbor_angle
      include "lammps.h"

      nanglelocal = 0
      i = atompnt
      do ii = 1,nlocal
        do k = 1,numangle(i)
          j1 = localptr(angleatom1(k,i))
          j2 = localptr(angleatom2(k,i))
          j3 = localptr(angleatom3(k,i))
          if (j1.eq.0.or.j2.eq.0.or.j3.eq.0) then
            write (6,*) 'Angle partner missing:',node,
     $           angleatom1(k,i),angleatom2(k,i),angleatom3(k,i),
     $           numangle(i),ntimestep
            call exit(0)
          endif
          if (newton_bond.eq.1.or.
     $         (i.le.j1.and.i.le.j2.and.i.le.j3)) then
            nanglelocal = nanglelocal + 1
            if (nanglelocal.le.maxangle) then
              anglelist(1,nanglelocal) = j1
              anglelist(2,nanglelocal) = j2
              anglelist(3,nanglelocal) = j3
              anglelist(4,nanglelocal) = angletype(k,i)
            endif
          endif
        enddo
        i = list(i)
      enddo

c check for too many angles

      if (nanglelocal.gt.maxangle) then
        write (6,*) 'Angle list too big:',node,nanglelocal
        call exit(0)
      endif

      max_angle = max(max_angle,nanglelocal)

      return
      end


c neighbor list construction of dihedrals
c  for each dihed store 4 local atoms in dihed and ptr to dihed parameters 
c  j1234 = 0 check is to see if atom in dihed is not in own or other list
c  if so is an error

      subroutine neighbor_dihedral
      include "lammps.h"

      ndihedlocal = 0
      i = atompnt
      do ii = 1,nlocal
        do k = 1,numdihed(i)
          j1 = localptr(dihedatom1(k,i))
          j2 = localptr(dihedatom2(k,i))
          j3 = localptr(dihedatom3(k,i))
          j4 = localptr(dihedatom4(k,i))
          if (j1.eq.0.or.j2.eq.0.or.j3.eq.0.or.j4.eq.0) then
            write (6,*) 'Dihedral partner missing:',node,
     $           dihedatom1(k,i),dihedatom2(k,i),dihedatom3(k,i),
     $           dihedatom4(k,i),numdihed(i),ntimestep
            call exit(0)
          endif
          if (newton_bond.eq.1.or.
     $         (i.le.j1.and.i.le.j2.and.i.le.j3.and.i.le.j4)) then
            ndihedlocal = ndihedlocal + 1
            if (ndihedlocal.le.maxdihed) then
              dihedlist(1,ndihedlocal) = j1
              dihedlist(2,ndihedlocal) = j2
              dihedlist(3,ndihedlocal) = j3
              dihedlist(4,ndihedlocal) = j4
              dihedlist(5,ndihedlocal) = dihedtype(k,i)
            endif
          endif
        enddo
        i = list(i)
      enddo

c check for too many diheds

      if (ndihedlocal.gt.maxdihed) then
        write (6,*) 'Dihedral list too big:',node,ndihedlocal
        call exit(0)
      endif

      max_dihed = max(max_dihed,ndihedlocal)

      return
      end


c neighbor list construction of impropers
c  for each impro store 4 local atoms in impro and ptr to impro parameters 
c  j1234 = 0 check is to see if atom in impro is not in own or other list
c  if so is an error

      subroutine neighbor_improper
      include "lammps.h"

      nimprolocal = 0
      i = atompnt
      do ii = 1,nlocal
        do k = 1,numimpro(i)
          j1 = localptr(improatom1(k,i))
          j2 = localptr(improatom2(k,i))
          j3 = localptr(improatom3(k,i))
          j4 = localptr(improatom4(k,i))
          if (j1.eq.0.or.j2.eq.0.or.j3.eq.0.or.j4.eq.0) then
            write (6,*) 'Improper partner missing:',node,
     $           improatom1(k,i),improatom2(k,i),improatom3(k,i),
     $           improatom4(k,i),numimpro(i),ntimestep
            call exit(0)
          endif
          if (newton_bond.eq.1.or.
     $         (i.le.j1.and.i.le.j2.and.i.le.j3.and.i.le.j4)) then
            nimprolocal = nimprolocal + 1
            if (nimprolocal.le.maximpro) then
              improlist(1,nimprolocal) = j1
              improlist(2,nimprolocal) = j2
              improlist(3,nimprolocal) = j3
              improlist(4,nimprolocal) = j4
              improlist(5,nimprolocal) = improtype(k,i)
            endif
          endif
        enddo
        i = list(i)
      enddo

c check for too many impros

      if (nimprolocal.gt.maximpro) then
        write (6,*) 'Improper list too big:',node,nimprolocal
        call exit(0)
      endif

      max_impro = max(max_impro,nimprolocal)

      return
      end
