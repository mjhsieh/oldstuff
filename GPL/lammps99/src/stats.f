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

c -----------------------------------------------------------------------
c statistics on data

      subroutine stats(data,n,ave,xmax,xmin,histo,histotmp,nhisto)
      implicit real*8 (a-h,o-z)
      include "mpif.h"
      real*8 data(*),ave,xmax,xmin
      integer histo(*),histotmp(*)
      
      xmin = 1.0E20
      xmax = -1.0E20
      ave = 0.0
      do i = 1,n
        ave = ave + data(i)
        if (data(i).lt.xmin) xmin = data(i)
        if (data(i).gt.xmax) xmax = data(i)
      enddo

      ntot = n
      call mpi_allreduce(ntot,ntot1,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      ave = ave/ntot1

      call mpi_allreduce(ave,tmp,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      ave = tmp
      call mpi_allreduce(xmax,tmp,1,mpi_double_precision,mpi_max,
     $     mpi_comm_world,ierror)
      xmax = tmp
      call mpi_allreduce(xmin,tmp,1,mpi_double_precision,mpi_min,
     $     mpi_comm_world,ierror)
      xmin = tmp

      do i = 1,nhisto
        histotmp(i) = 0
      enddo
      
      del = xmax-xmin
      do i = 1,n
        if (del.eq.0.0) then
          j = 1
        else
          j = (data(i)-xmin)/del * nhisto + 1
          if (j.gt.nhisto) j = nhisto
        endif
        histotmp(j) = histotmp(j) + 1
      enddo
      
      call mpi_allreduce(histotmp,histo,nhisto,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

      return
      end
