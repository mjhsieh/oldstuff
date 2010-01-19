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
c map density in 3-d brick (with ghosts) to density in FFT decomposition
c density_in = 3-d brick decomp
c density_out = FFT decomp

      subroutine brick2fft(density_in,density_out)
      include "lammps.h"
      include "mpif.h"
      real*8 density_in(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 density_out(*)

      integer istatus(mpi_status_size)

c pack my ghosts for +x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in+1,nxhi_out
            icount = icount + 1
            buf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

c pass data to self or +x processor

      if (mpart(2,1).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(1,1),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(2,1),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_in,nxlo_in+nxlo_ghost-1
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + buf2(icount)
          enddo
        enddo
      enddo

c pack my ghosts for -x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_out,nxlo_in-1
            icount = icount + 1
            buf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

c pass data to self or -x processor

      if (mpart(1,1).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(2,1),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(1,1),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in-nxhi_ghost+1,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + buf2(icount)
          enddo
        enddo
      enddo

c pack my ghosts for +y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in+1,nyhi_out
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            buf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

c pass data to self or +y processor

      if (mpart(2,2).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(1,2),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(2,2),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_in,nylo_in+nylo_ghost-1
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + buf2(icount)
          enddo
        enddo
      enddo

c pack my ghosts for -y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nylo_in-1
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            buf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

c pass data to self or -y processor

      if (mpart(1,2).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(2,2),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(1,2),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in-nyhi_ghost+1,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + buf2(icount)
          enddo
        enddo
      enddo

c pack my ghosts for +z processor

      icount = 0
      do iz = nzhi_in+1,nzhi_out
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            buf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

c pass data to self or +z processor

      if (mpart(2,3).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(1,3),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(2,3),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my real cells

      icount = 0
      do iz = nzlo_in,nzlo_in+nzlo_ghost-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + buf2(icount)
          enddo
        enddo
      enddo

c pack my ghosts for -z processor

      icount = 0
      do iz = nzlo_out,nzlo_in-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            buf1(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

c pass data to self or -z processor

      if (mpart(1,3).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(2,3),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(1,3),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my real cells

      icount = 0
      do iz = nzhi_in-nzhi_ghost+1,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_in(ix,iy,iz) = density_in(ix,iy,iz) + buf2(icount)
          enddo
        enddo
      enddo

c remap from 3-d brick decomposition to FFT decomposition
c copy done first to grab only inner portion of density from 3-d brick
c remap could be done as pre-stage of FFT,
c  but this one works optimally on only real*8 values, not complex values

      icount = 0
      do iz = nzlo_in,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            density_out(icount) = density_in(ix,iy,iz)
          enddo
        enddo
      enddo

      call remap_3d(density_out,density_out,workvec1,plan_remap)

      return
      end


c ----------------------------------------------------------------------
c fill-in ghost values of potential gradients on each 3-d brick

      subroutine fillbrick(vdx,vdy,vdz)
      include "lammps.h"
      include "mpif.h"
      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 vdy(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 vdz(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)

      integer istatus(mpi_status_size)

c pack my real cells for +z processor

      icount = 0
      do iz = nzhi_in-nzhi_ghost+1,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            buf1(icount+1) = vdx(ix,iy,iz)
            buf1(icount+2) = vdy(ix,iy,iz)
            buf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

c pass data to self or +z processor

      if (mpart(2,3).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(1,3),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(2,3),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzlo_in-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = buf2(icount+1)
            vdy(ix,iy,iz) = buf2(icount+2)
            vdz(ix,iy,iz) = buf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

c pack my real cells for -z processor

      icount = 0
      do iz = nzlo_in,nzlo_in+nzlo_ghost-1
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            buf1(icount+1) = vdx(ix,iy,iz)
            buf1(icount+2) = vdy(ix,iy,iz)
            buf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

c pass data to self or +z processor

      if (mpart(1,3).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(2,3),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(1,3),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my ghost cells

      icount = 0
      do iz = nzhi_in+1,nzhi_out
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = buf2(icount+1)
            vdy(ix,iy,iz) = buf2(icount+2)
            vdz(ix,iy,iz) = buf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

c pack my real cells for +y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in-nyhi_ghost+1,nyhi_in
          do ix = nxlo_in,nxhi_in
            buf1(icount+1) = vdx(ix,iy,iz)
            buf1(icount+2) = vdy(ix,iy,iz)
            buf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

c pass data to self or +y processor

      if (mpart(2,2).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(1,2),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(2,2),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nylo_in-1
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = buf2(icount+1)
            vdy(ix,iy,iz) = buf2(icount+2)
            vdz(ix,iy,iz) = buf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

c pack my real cells for -y processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_in,nylo_in+nylo_ghost-1
          do ix = nxlo_in,nxhi_in
            buf1(icount+1) = vdx(ix,iy,iz)
            buf1(icount+2) = vdy(ix,iy,iz)
            buf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

c pass data to self or -y processor

      if (mpart(1,2).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(2,2),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(1,2),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nyhi_in+1,nyhi_out
          do ix = nxlo_in,nxhi_in
            vdx(ix,iy,iz) = buf2(icount+1)
            vdy(ix,iy,iz) = buf2(icount+2)
            vdz(ix,iy,iz) = buf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

c pack my real cells for +x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in-nxhi_ghost+1,nxhi_in
            buf1(icount+1) = vdx(ix,iy,iz)
            buf1(icount+2) = vdy(ix,iy,iz)
            buf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

c pass data to self or +x processor

      if (mpart(2,1).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(1,1),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(2,1),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_out,nxlo_in-1
            vdx(ix,iy,iz) = buf2(icount+1)
            vdy(ix,iy,iz) = buf2(icount+2)
            vdz(ix,iy,iz) = buf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

c pack my real cells for -x processor

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxlo_in,nxlo_in+nxlo_ghost-1
            buf1(icount+1) = vdx(ix,iy,iz)
            buf1(icount+2) = vdy(ix,iy,iz)
            buf1(icount+3) = vdz(ix,iy,iz)
            icount = icount + 3
          enddo
        enddo
      enddo

c pass data to self or +x processor

      if (mpart(1,1).eq.node) then
        do i = 1,icount
          buf2(i) = buf1(i)
        enddo
      else
        call mpi_irecv(buf2,maxexchtot,mpi_double_precision,
     $       mpart(2,1),0,mpi_comm_world,irequest,ierror)
        call mpi_send(buf1,icount,mpi_double_precision,
     $       mpart(1,1),0,mpi_comm_world,ierror)
        call mpi_wait(irequest,istatus,ierror)
      endif

c unpack and sum recv data into my ghost cells

      icount = 0
      do iz = nzlo_out,nzhi_out
        do iy = nylo_out,nyhi_out
          do ix = nxhi_in+1,nxhi_out
            vdx(ix,iy,iz) = buf2(icount+1)
            vdy(ix,iy,iz) = buf2(icount+2)
            vdz(ix,iy,iz) = buf2(icount+3)
            icount = icount + 3
          enddo
        enddo
      enddo

      return
      end
