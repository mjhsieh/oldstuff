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
c write a restart file in native binary format

      subroutine write_restart
      include "lammps.h"
      include "mpif.h"

      integer istatus(mpi_status_size)
      
      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()

c make sure bounds are current if non-PBC
c doesn't matter that bounds are changed without calling setup_box
c nothing in code depends on bounds until setup_box is called

      if (perflagx+perflagy+perflagz.gt.0) call setup_bounds

c toggle restart files

      if (node.eq.0) then
        if (restartfileflag.eq.0) then
          open (unit=2,file=restart_out1,form='unformatted',
     $         status='unknown')
          close (2,status='delete')
          open (unit=2,file=restart_out1,form='unformatted',
     $         status='new')
          restartfileflag = 1
        else
          open (unit=2,file=restart_out2,form='unformatted',
     $         status='unknown')
          close (2,status='delete')
          open (unit=2,file=restart_out2,form='unformatted',
     $         status='new')
          restartfileflag = 0
        endif
      endif

c write header
c include all parameters that must be set before "read restart" command
c   so can check for match when restart is done
c include all values that insure "perfect" restart capability
c write dummy nbyte,ibyte,eta (version 5)

      if (node.eq.0) then
        nbyte = 0
        ibyte = 0
        write (2) iversion,nbyte,ibyte
        write (2) ntimestep,nprocs,pgrid(1),pgrid(2),pgrid(3)
        write (2) idimension,perflagx,perflagy,perflagz
        write (2) units,newton
        write (2) natoms,nbonds,nangles,ndihedrals,nimpropers
        write (2) ntypes
        if (nbonds.gt.0) write (2) nbondtypes
        if (nangles.gt.0) write (2) nangletypes
        if (ndihedrals.gt.0) write (2) ndihedtypes
        if (nimpropers.gt.0) write (2) nimprotypes
        write (2) max_bondper,max_angleper,max_dihedper,max_improper
        write (2) xboundlo,xboundhi,yboundlo,yboundhi,
     $       zboundlo,zboundhi
        eta = 0.0
        write (2) eta,eta_dot,omega_dot(1),omega_dot(2),omega_dot(3)
      endif

c write all coefficients

      if (node.eq.0) then
        
        write (2) (mass(i),i=1,ntypes)

        write (2) nonstyle
        if (nonstyle.eq.1) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.2) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.3) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.4) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.5) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        endif

        if (nbonds.gt.0) then
          write (2) bondstyle
          if (bondstyle.eq.1) then
            write (2) ((bondcoeff(i,j),i=1,2),j=1,nbondtypes)
          else if (bondstyle.eq.2) then
            write (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          else if (bondstyle.eq.3) then
            write (2) ((bondcoeff(i,j),i=1,5),j=1,nbondtypes)
          else if (bondstyle.eq.4) then
            write (2) ((bondcoeff(i,j),i=1,3),j=1,nbondtypes)
          else if (bondstyle.eq.5) then
            write (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          endif
        endif

        if (nangles.gt.0) then
          write (2) anglestyle
          if (anglestyle.eq.1) then
            write (2) ((anglecoeff(i,j),i=1,2),j=1,nangletypes)
          else if (anglestyle.eq.2) then
            write (2) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
            write (2) ((bondbondcoeff(i,j),i=1,3),j=1,nangletypes)
            write (2) ((bondanglecoeff(i,j),i=1,4),j=1,nangletypes)
          endif
        endif

        if (ndihedrals.gt.0) then
          write (2) dihedstyle
          if (dihedstyle.eq.1) then
            write (2) ((dihedcoeff(i,j),i=1,3),j=1,ndihedtypes)
          else if (dihedstyle.eq.2) then
            write (2) ((dihedcoeff(i,j),i=1,6),j=1,ndihedtypes)
            write (2) ((midbondtorsioncoeff(i,j),i=1,4),
     $           j=1,ndihedtypes)
            write (2) ((endbondtorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            write (2) ((angletorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            write (2) ((angleangletorsioncoeff(i,j),i=1,3),
     $           j=1,ndihedtypes)
          endif
        endif

        if (nimpropers.gt.0) then
          write (2) improstyle
          if (improstyle.eq.1) then
            write (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
          else if (improstyle.eq.2) then
            write (2) ((improcoeff(i,j),i=1,3),j=1,nimprotypes)
          else if (improstyle.eq.3) then
            write (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
            write (2) ((angleanglecoeff(i,j),i=1,6),j=1,nimprotypes)
          endif
        endif

      endif

c create reverse list of atom pointers
c this is so that read-in restart will not be in reversed order
c this should enable exact restarts,
c  so long as number of processors stays the same

      i = atompnt
      do ii = nlocal,1,-1
        buf4(ii) = i
        i = list(i)
      enddo

c load buffer with all of my owned atom information
c use reverse list to load them in buffer in reverse order

      j = 1
      buf5(j) = nlocal
      do ii = 1,nlocal
        i = nint(buf4(ii))
        buf5(j+1) = x(1,i)
        buf5(j+2) = x(2,i)
        buf5(j+3) = x(3,i)
        buf5(j+4) = tag(i)
        buf5(j+5) = molecule(i)
        buf5(j+6) = type(i)
        buf5(j+7) = q(i)
        buf5(j+8) = v(1,i)
        buf5(j+9) = v(2,i)
        buf5(j+10) = v(3,i)
        buf5(j+11) = true(i)
        buf5(j+12) = fix(i)
        buf5(j+13) = numbond(i)
        buf5(j+14) = numangle(i)
        buf5(j+15) = numdihed(i)
        buf5(j+16) = numimpro(i)
        j = j + peratom_io
        do jj = 1,numbond(i)
          buf5(j+1) = bondatom1(jj,i)
          buf5(j+2) = bondatom2(jj,i)
          buf5(j+3) = bondtype(jj,i)
          j = j + 3
        enddo
        do jj = 1,numangle(i)
          buf5(j+1) = angleatom1(jj,i)
          buf5(j+2) = angleatom2(jj,i)
          buf5(j+3) = angleatom3(jj,i)
          buf5(j+4) = angletype(jj,i)
          j = j + 4
        enddo
        do jj = 1,numdihed(i)
          buf5(j+1) = dihedatom1(jj,i)
          buf5(j+2) = dihedatom2(jj,i)
          buf5(j+3) = dihedatom3(jj,i)
          buf5(j+4) = dihedatom4(jj,i)
          buf5(j+5) = dihedtype(jj,i)
          j = j + 5
        enddo
        do jj = 1,numimpro(i)
          buf5(j+1) = improatom1(jj,i)
          buf5(j+2) = improatom2(jj,i)
          buf5(j+3) = improatom3(jj,i)
          buf5(j+4) = improatom4(jj,i)
          buf5(j+5) = improtype(jj,i)
          j = j + 5
        enddo
      enddo
      icnt = j

c node 0 pings each node, receives their buffer, writes to file atom by atom
c  all other nodes wait for ping, send buffer to node 0

      if (node.eq.0) then

        do inode = 0,nprocs-1
          if (inode.ne.0) then
            call mpi_irecv(buf5,maxrestot,mpi_double_precision,
     $           inode,0,mpi_comm_world,irequest,ierror)
            call mpi_send(itmp,0,mpi_integer,inode,0,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
          endif

          j = 1
          nlocal_tmp = nint(buf5(j))
          do i = 1,nlocal_tmp
            write (2) (buf5(j+k),k=1,peratom_io)
            j = j + peratom_io
            numbond_tmp = nint(buf5(j-3))
            numangle_tmp = nint(buf5(j-2))
            numdihed_tmp = nint(buf5(j-1))
            numimpro_tmp = nint(buf5(j))
            do jj = 1,numbond_tmp
              write (2) buf5(j+1),buf5(j+2),buf5(j+3)
              j = j + 3
            enddo
            do jj = 1,numangle_tmp
              write (2) buf5(j+1),buf5(j+2),buf5(j+3),buf5(j+4)
              j = j + 4
            enddo
            do jj = 1,numdihed_tmp
              write (2) buf5(j+1),buf5(j+2),buf5(j+3),
     $             buf5(j+4),buf5(j+5)
              j = j + 5
            enddo
            do jj = 1,numimpro_tmp
              write (2) buf5(j+1),buf5(j+2),buf5(j+3),
     $             buf5(j+4),buf5(j+5)
              j = j + 5
            enddo
          enddo

        enddo

      else

        call mpi_recv(itmp,0,mpi_integer,0,0,
     $       mpi_comm_world,istatus,ierror)
        call mpi_rsend(buf5,icnt,mpi_double_precision,
     $       0,0,mpi_comm_world,ierror)

      endif

c close restart file

      if (node.eq.0) close (2)

c update timestep counter

      nrestart_next = min(ntimestep+nrestart,ntime_last)

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1

      return
      end
