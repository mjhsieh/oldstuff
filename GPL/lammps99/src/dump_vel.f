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
c dump atom velocities to a file

      subroutine dump_vel
      include "lammps.h"
      include "mpif.h"

      real*8 time1,time2
      integer istatus(mpi_status_size)

 900  format (i6,i3,3f10.5)

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()

c do not write this timestep if previously written, just update counter

      if (ntimestep.eq.ndumpvel_prev) then
        ndumpvel_next = min(ntimestep+ndumpvel,ntime_last)
        return
      endif

c 1st time writing to new dump file, write overall file header

      if (dumpvelfileflag.eq.1) then
        dumpvelfileflag = 2
        if (node.eq.0) then
          write (12,*) 'ITEM: NUMBER OF ATOMS'
          write (12,*) natoms
        endif
      endif

c write timestep header

      if (node.eq.0) then
        write (12,*) 'ITEM: TIMESTEP'
        write (12,*) ntimestep
        write (12,*) 'ITEM: VELOCITIES'
      endif

c number of datums per atom

      nper = 5

c buffer up my atom velocities (scaled by dtfactor)

      rinverse = 1.0/dtfactor
      j = 0
      i = atompnt
      do ii = 1,nlocal
        buf4(j+1) = tag(i)
        buf4(j+2) = type(i)
        buf4(j+3) = v(1,i)*rinverse
        buf4(j+4) = v(2,i)*rinverse
        buf4(j+5) = v(3,i)*rinverse
        j = j + nper
        i = list(i)
      enddo
      nlocal_tmp = nlocal

c node 0 pings each node, receives their buffer, writes to file
c  all other nodes wait for ping, send buffer to node 0

      if (node.eq.0) then

        do inode = 0,nprocs-1
          if (inode.ne.0) then
            call mpi_irecv(buf4,nper*maxown,mpi_double_precision,
     $           inode,0,mpi_comm_world,irequest,ierror)
            call mpi_send(itmp,0,mpi_integer,inode,0,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_double_precision,
     $           nlocal_tmp,ierror)
            nlocal_tmp = nlocal_tmp/nper
          endif

          j = 0
          do i = 1,nlocal_tmp
            write (12,900) nint(buf4(j+1)),nint(buf4(j+2)),
     $           buf4(j+3),buf4(j+4),buf4(j+5)
            j = j + nper
          enddo
        enddo

      else

        call mpi_recv(itmp,0,mpi_integer,0,0,
     $       mpi_comm_world,istatus,ierror)
        call mpi_rsend(buf4,nper*nlocal_tmp,mpi_double_precision,
     $       0,0,mpi_comm_world,ierror)

      endif

c update timestep counters

      ndumpvel_next = min(ntimestep+ndumpvel,ntime_last)
      ndumpvel_prev = ntimestep

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1

      return
      end


c -------------------------------------------------------------------------
c open a velocity dump file

      subroutine dumpvelopen
      include "lammps.h"

      if (node.eq.0) then
        open (12,file=dumpvelfile,status='unknown')
        close (12,status='delete')
        open (12,file=dumpvelfile,status='new')
      endif

c no velocity dumps have been written into new file
      
      ndumpvel_prev = -1

      return
      end


c -------------------------------------------------------------------------
c close a velocity dump file

      subroutine dumpvelclose
      include "lammps.h"

      if (node.eq.0) close (12)

      return
      end
