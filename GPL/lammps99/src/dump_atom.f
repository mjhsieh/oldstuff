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
c dump atom positions to a file

      subroutine dump_atom
      include "lammps.h"
      include "mpif.h"

      real*8 time1,time2
      integer istatus(mpi_status_size)

 900  format (i6,i3,3f8.5,3i5)

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()

c do not write this timestep if previously written, just update counter

      if (ntimestep.eq.ndumpatom_prev) then
        ndumpatom_next = min(ntimestep+ndumpatom,ntime_last)
        return
      endif

c 1st time writing to new dump file, write overall file header
c  only write box bounds if periodic and constant volume

      if (dumpatomfileflag.eq.1) then
        dumpatomfileflag = 2
        if (node.eq.0) then
          write (11,*) 'ITEM: NUMBER OF ATOMS'
          write (11,*) natoms
          if (perflagx+perflagy+perflagz.eq.0.and.ensemble.le.2) then
            write (11,*) 'ITEM: BOX BOUNDS'
            write (11,*) xboundlo,xboundhi
            write (11,*) yboundlo,yboundhi
            if (idimension.eq.3) write (11,*) zboundlo,zboundhi
          endif
        endif
      endif

c write timestep header

      if (node.eq.0) then
        write (11,*) 'ITEM: TIMESTEP'
        write (11,*) ntimestep
      endif

c output current box bounds if non-PBC or constant P
c make sure bounds are current if non-PBC
c  doesn't matter that bounds are changed without calling setup_box
c  nothing in code depends on bounds until setup_box is called

      if (perflagx+perflagy+perflagz.gt.0.or.ensemble.eq.3) then
        if (perflagx+perflagy+perflagz.gt.0) call setup_bounds
        if (node.eq.0) then
          write (11,*) 'ITEM: BOX BOUNDS'
          write (11,*) xboundlo,xboundhi
          write (11,*) yboundlo,yboundhi
          if (idimension.eq.3) write (11,*) zboundlo,zboundhi
        endif
      endif

c write timestep header

      if (node.eq.0) write (11,*) 'ITEM: ATOMS'

c number of datums per atom

      nper = 5
      if (mod(trueflag,2).eq.1) nper = 8

c buffer up my atoms

      j = 0
      i = atompnt
      do ii = 1,nlocal
        buf4(j+1) = tag(i)
        buf4(j+2) = type(i)
        buf4(j+3) = (x(1,i)-xboundlo)/xprd
        buf4(j+4) = (x(2,i)-yboundlo)/yprd
        buf4(j+5) = (x(3,i)-zboundlo)/zprd
        if (mod(trueflag,2).eq.1) then
          buf4(j+6) = mod(true(i),1000) - 500
          buf4(j+7) = mod(true(i)/1000,1000) - 500
          buf4(j+8) = true(i)/1000000 - 500
        endif
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
          if (nper.eq.5) then
            do i = 1,nlocal_tmp
              write (11,900) nint(buf4(j+1)),nint(buf4(j+2)),
     $             buf4(j+3),buf4(j+4),buf4(j+5)
              j = j + nper
            enddo
          else
            do i = 1,nlocal_tmp
              write (11,900) nint(buf4(j+1)),nint(buf4(j+2)),
     $             buf4(j+3),buf4(j+4),buf4(j+5),
     $             nint(buf4(j+6)),nint(buf4(j+7)),nint(buf4(j+8))
              j = j + nper
            enddo
          endif
        enddo

      else

        call mpi_recv(itmp,0,mpi_integer,0,0,
     $       mpi_comm_world,istatus,ierror)
        call mpi_rsend(buf4,nper*nlocal_tmp,mpi_double_precision,
     $       0,0,mpi_comm_world,ierror)

      endif

c update timestep counters

      ndumpatom_next = min(ntimestep+ndumpatom,ntime_last)
      ndumpatom_prev = ntimestep

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1

      return
      end


c -------------------------------------------------------------------------
c open a atom dump file

      subroutine dumpatomopen
      include "lammps.h"

      if (node.eq.0) then
        open (11,file=dumpatomfile,status='unknown')
        close (11,status='delete')
        open (11,file=dumpatomfile,status='new')
      endif

c no atom dumps have been written into new file
      
      ndumpatom_prev = -1

      return
      end


c -------------------------------------------------------------------------
c close a atom dump file

      subroutine dumpatomclose
      include "lammps.h"

      if (node.eq.0) close (11)

      return
      end
