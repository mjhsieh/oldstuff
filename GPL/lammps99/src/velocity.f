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
c create velocities for atoms according to createstyle and creategroup
c   createstyle = 1 -> init atoms by uniform RNG to t_create
c   createstyle = 2 -> init atoms by gaussian RNG to t_create
c   createstyle = 3 -> init atoms to velocity of createvec(3)
c   creategroup = 0 -> apply to all atoms
c   creategroup = 1 -> only apply to atoms with 
c                        createtypelo <= type <= createtypehi
c   creategroup = 2 -> apply to atoms within a spatial region
c   creategroup = 3 -> apply to all unassigned atoms

      subroutine temp_create
      include "lammps.h"
      include "mpif.h"

      if (node.eq.0) write (6,*) 'Creating temperature ...'

      if (readflag.eq.0)
     $     call error('Data file must be read before creating temp')

      if (creategroup.eq.1) then
        if (createtypelo.gt.ntypes.or.createtypehi.gt.ntypes)
     $       call error('Creating temp for too large a type')
      endif

c mark atoms whose velocities will be initialized and count them

      ncount = 0
      countmass = 0.0
      i = atompnt
      do ii = 1,nlocal
        iflag = 0
        if (creategroup.eq.0) then
          iflag = 1
        else if (creategroup.eq.1) then
          if (createtypelo.le.type(i).and.
     $         type(i).le.createtypehi) iflag = 1
        else if (creategroup.eq.2) then
          if (x(1,i).ge.createregion(1).and.
     $         x(1,i).le.createregion(2).and.
     $         x(2,i).ge.createregion(3).and.
     $         x(2,i).le.createregion(4).and.
     $         x(3,i).ge.createregion(5).and.
     $         x(3,i).le.createregion(6)) iflag = 1
        else if (creategroup.eq.3) then
          if (velflag(i).eq.0) iflag = 1
        endif
        if (iflag.eq.1) then
          velflag(i) = -1
          ncount = ncount + 1
          countmass = countmass + mass(type(i))
        endif
        i = list(i)
      enddo

      itmp = ncount
      call mpi_allreduce(itmp,ncount,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      if (ncount.eq.0) return

      tmp = countmass
      call mpi_allreduce(tmp,countmass,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      if (createstyle.le.2) then

c for uniform and Gaussian initializers, loop over all atoms to insure
c  RNG usage assigns velocities the same independent of P

        do j = 1,natoms
          if (createstyle.eq.1) then
            vx = ranpark(iseed)
            vy = ranpark(iseed)
            vz = ranpark(iseed)
          else
            vx = gausspark(iseed)
            vy = gausspark(iseed)
            vz = gausspark(iseed)
          endif
          i = localptr(j)
          if (i.ne.0) then
            if (velflag(i).eq.-1) then
              v(1,i) = vx
              v(2,i) = vy
              v(3,i) = vz
              if (idimension.eq.2) v(3,i) = 0.0
            endif
          endif
        enddo

c zero out momentum of only initialized atoms (unless just one atom)

        if (ncount.gt.1) call momentum_zero(ncount,countmass)

c compute temp of only initialized atoms

        t_current = 0.0
        i = atompnt
        do ii = 1,nlocal
          if (velflag(i).eq.-1) t_current = t_current +
     $         (v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)) *
     $         mass(type(i))
          i = list(i)
        enddo
        call mpi_allreduce(t_current,tmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        t_current = tmp / (boltz*idimension*ncount)

c rescale only initialized atoms to temp t_create

        if (t_current.ne.0.0) then
          factor = sqrt(t_create/t_current)
          i = atompnt
          do ii = 1,nlocal
            if (velflag(i).eq.-1) then
              v(1,i) = v(1,i)*factor
              v(2,i) = v(2,i)*factor
              v(3,i) = v(3,i)*factor
            endif
            i = list(i)
          enddo
        endif

      else if (createstyle.eq.3) then

        i = atompnt
        do ii = 1,nlocal
          if (velflag(i).eq.-1) then
            v(1,i) = createvec(1)
            v(2,i) = createvec(2)
            v(3,i) = createvec(3)
          endif
          i = list(i)
        enddo

      endif

c mark atoms as initialized

      i = atompnt
      do ii = 1,nlocal
        if (velflag(i).eq.-1) velflag(i) = 1
        i = list(i)
      enddo

c reset creategroup to all atoms

      creategroup = 0

      return
      end


c rescale temperature to target temp

      subroutine temp_rescale(target)
      include "lammps.h"

      factor = sqrt(target/t_current)
      i = atompnt
      do ii = 1,nlocal
        v(1,i) = v(1,i)*factor
        v(2,i) = v(2,i)*factor
        v(3,i) = v(3,i)*factor
        i = list(i)
      enddo
      t_current = target

      return
      end


c Gaussian replacement of all velocities at target temp

      subroutine temp_replace(target)
      include "lammps.h"
      
      i = atompnt
      do ii = 1,nlocal
        v(1,i) = gaussmars()
        v(2,i) = gaussmars()
        v(3,i) = gaussmars()
        if (idimension.eq.2) v(3,i) = 0.0
        i = list(i)
      enddo
      call momentum_zero(0,masssum)
      call temperature
      if (t_current.ne.0.0) call temp_rescale(target)
      
      return
      end


c subtract center-of-mass momentum from velocities
c  ncount = 0 -> include all atoms
c  ncount > 1 -> include only atoms tagged with velflag(i) = -1

      subroutine momentum_zero(ncount,countmass)
      include "lammps.h"
      include "mpif.h"
      integer ncount
      real*8 countmass

      real*8 ptot(3),ptotall(3)

      ptot(1) = 0.0
      ptot(2) = 0.0
      ptot(3) = 0.0

      if (ncount.eq.0) then
        i = atompnt
        do ii = 1,nlocal
          rmass = mass(type(i))
          ptot(1) = ptot(1) + v(1,i)*rmass
          ptot(2) = ptot(2) + v(2,i)*rmass
          ptot(3) = ptot(3) + v(3,i)*rmass
          i = list(i)
        enddo
      else
        i = atompnt
        do ii = 1,nlocal
          if (velflag(i).eq.-1) then
            rmass = mass(type(i))
            ptot(1) = ptot(1) + v(1,i)*rmass
            ptot(2) = ptot(2) + v(2,i)*rmass
            ptot(3) = ptot(3) + v(3,i)*rmass
          endif
          i = list(i)
        enddo
      endif

      call mpi_allreduce(ptot,ptotall,3,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      if (ncount.eq.0) then
        i = atompnt
        do ii = 1,nlocal
          v(1,i) = v(1,i) - ptotall(1)/countmass
          v(2,i) = v(2,i) - ptotall(2)/countmass
          v(3,i) = v(3,i) - ptotall(3)/countmass
          i = list(i)
        enddo
      else
        i = atompnt
        do ii = 1,nlocal
          if (velflag(i).eq.-1) then
            v(1,i) = v(1,i) - ptotall(1)/countmass
            v(2,i) = v(2,i) - ptotall(2)/countmass
            v(3,i) = v(3,i) - ptotall(3)/countmass
          endif
          i = list(i)
        enddo
      endif

      return
      end
