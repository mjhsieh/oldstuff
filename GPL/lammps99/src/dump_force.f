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
c dump atom forces to a file

      subroutine dump_force
      include "lammps.h"
      include "mpif.h"

      real*8 time1,time2
      integer istatus(mpi_status_size)

 900  format (i6,i3,3f9.4)
 901  format (i6,i3,3f11.6)

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      
c do not write this timestep if previously written, just update counter

      if (ntimestep.eq.ndumpforce_prev) then
        ndumpforce_next = min(ntimestep+ndumpforce,ntime_last)
        return
      endif

c 1st time writing to new dump file, write overall file header

      if (dumpforcefileflag.eq.1) then
        dumpforcefileflag = 2
        if (node.eq.0) then
          write (13,*) 'ITEM: NUMBER OF ATOMS'
          write (13,*) natoms
        endif
      endif

c write timestep header

      if (node.eq.0) then
        write (13,*) 'ITEM: TIMESTEP'
        write (13,*) ntimestep
        write (13,*) 'ITEM: FORCES'
      endif

c number of datums per atom

      nper = 5

c buffer up my atoms
c if respa, then sum 4 partial force arrays

      if (nrespa.eq.0) then
        j = 0
        i = atompnt
        do ii = 1,nlocal
          buf4(j+1) = tag(i)
          buf4(j+2) = type(i)
          buf4(j+3) = f(1,i)
          buf4(j+4) = f(2,i)
          buf4(j+5) = f(3,i)
          j = j + nper
          i = list(i)
        enddo
      else
        j = 0
        i = atompnt
        do ii = 1,nlocal
          buf4(j+1) = tag(i)
          buf4(j+2) = type(i)
          buf4(j+3) = f_stretch(1,i) + f_intra(1,i) +
     $         f_short(1,i) + f_long(1,i)
          buf4(j+4) = f_stretch(2,i) + f_intra(2,i) +
     $         f_short(2,i) + f_long(2,i)
          buf4(j+5) = f_stretch(3,i) + f_intra(3,i) +
     $         f_short(3,i) + f_long(3,i)
          j = j + nper
          i = list(i)
        enddo
      endif
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
            write (13,900) nint(buf4(j+1)),nint(buf4(j+2)),
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

      ndumpforce_next = min(ntimestep+ndumpforce,ntime_last)
      ndumpforce_prev = ntimestep

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1

      return
      end


c -------------------------------------------------------------------------
c open a force dump file

      subroutine dumpforceopen
      include "lammps.h"

      if (node.eq.0) then
        open (13,file=dumpforcefile,status='unknown')
        close (13,status='delete')
        open (13,file=dumpforcefile,status='new')
      endif

c no force dumps have been written into new file
      
      ndumpforce_prev = -1

      return
      end


c -------------------------------------------------------------------------
c close a force dump file

      subroutine dumpforceclose
      include "lammps.h"

      if (node.eq.0) close (13)

      return
      end
