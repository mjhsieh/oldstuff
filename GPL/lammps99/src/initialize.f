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
c initialize run parameters

      subroutine initialize
      include "lammps.h"
      include "mpif.h"

c identify self

      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,node,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)

c print out version and determine starting time
      
      iversion = 6
      if (node.eq.0) then
        time_total = mpi_wtime()
        write (6,*) 'LAMMPS 99 (June 1999)'
      endif
      call mpi_bcast(time_total,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

c open log file

      if (node.eq.0) then
        open (unit=1,file='log.lammps',status='unknown')
        close (1,status='delete')
        open (unit=1,file='log.lammps',status='new')
        write (1,*) 'LAMMPS 99 (June 1999)'
      endif

c constant

      two_1_3 = 2.0**(1.0/3.0)

c no data file has been read in yet

      readflag = 0

c number of quantities per atom

      peratom_comm = 19
      peratom_io = 16

c default settings

      units = 0
      boltz = 0.001987191
      dtfactor = 48.88821
      pfactor = 68589.796

      efactor = 332.054

      idimension = 3
      pgrid(1) = 0
      pgrid(2) = 0
      pgrid(3) = 0
      perflagx = 0
      perflagy = 0
      perflagz = 0
      newton = 3
      newton_bond = 1
      newton_nonbond = 1

      dt = 1.0
      nrespa = 0

      skin = 2.0
      neighstyle = 0
      neighfreq = 1
      neighdelay = 10
      neightrigger = 1

      special(1) = 0.0
      special(2) = 0.0
      special(3) = 0.4

      nthermo = 0
      thermostyle = 0
      trueflag = 0

      ndumpatom = 0
      dumpatomfileflag = 0
      ndumpvel = 0
      dumpvelfileflag = 0
      ndumpforce = 0
      dumpforcefileflag = 0

      nrestart = 0
      restartfileflag = 0

      call diag_init

      tempflag = 0
      pressflag = 0

      nonstyle = 1
      cutlj = 10.0
      offsetflag = 0
      do i = 1,maxtype
        do j = i,maxtype
          nontypeflag(i,j) = 0
        enddo
      enddo
      mixflag = 0
      mixstyle = 1

      coulstyle = 1
      cutcoul = 10.0
      meshflag = 0
      orderflag = 5
      dielectric = 1.0

      bondstyle = 1
      do i = 1,maxbondtype
        bondtypeflag(i) = 0
      enddo

      anglestyle = 1
      do i = 1,maxangletype
        angletypeflag(i) = 0
      enddo

      dihedstyle = 1
      do i = 1,maxdihedtype
        dihedtypeflag(i) = 0
      enddo

      improstyle = 1
      do i = 1,maximprotype
        improtypeflag(i) = 0
      enddo

      creategroup = 0
      nfixes = 0

c Ewald setting

      kmax = 0
      
c optimization settings

      optstyle = 1
      optfileflag = 0

      return
      end
