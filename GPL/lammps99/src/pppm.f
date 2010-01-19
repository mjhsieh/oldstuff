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
c PPPM routines authored by Roy Pollock (modified for LAMMPS by SJP)
c (510) 422-4088, pollockr@llnl.gov
c L-299, LLNL, PO Box 808, Livermore, CA  94550
c -------------------------------------------------------------------------

c -------------------------------------------------------------------------
c PPPM driver

      subroutine pppm(iflag)
      include "lammps.h"
      include "mpif.h"

c insure all particles are within my 3-d brick

      call particle_unmap

c map my particles onto my local 3-d density grid

      time1 = mpi_wtime()
      call make_rho(density_brick)
      time2 = mpi_wtime()
      time_rho = time_rho + time2-time1

c all procs communicate density from their 3-d brick (with ghosts)
c  to others FFT columns - result is filled/summed grid in FFT decomposition

      call brick2fft(density_brick,density_fft)

c compute potential gradient on my FFT grid and
c   portion of e_long on this proc's FFT grid
c return gradients in 3-d brick decomposition

      time1 = mpi_wtime()
      call poisson(vdx_brick,vdy_brick,vdz_brick)
      time2 = mpi_wtime()
      time_poiss = time_poiss + time2-time1

c all procs communicate electric fields from their FFT columns 
c  others 3-d bricks (with ghosts) - result is filled grid in 3-d decomp

      call fillbrick(vdx_brick,vdy_brick,vdz_brick)

c calculate the e-fields on my particles

      time1 = mpi_wtime()
      call electric_field(vdx_brick,vdy_brick,vdz_brick)
      time2 = mpi_wtime()
      time_field = time_field + time2-time1

c convert e-field to force

      i = atompnt
      do ii = 1,nlocal
        f(1,i) = f(1,i) + q(i)*coulpre*ek(1,ii)
        f(2,i) = f(2,i) + q(i)*coulpre*ek(2,ii)
        f(3,i) = f(3,i) + q(i)*coulpre*ek(3,ii)
        i = list(i)
      enddo

c sum total long-range component of energy

      tmp = e_long
      call mpi_allreduce(tmp,e_long,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)

      pi = 4.0*atan(1.0)
      e_long = e_long - ewald_g*qsqsum/1.77245385/sqrt(2.0) -
     $     0.5*pi*qsum*qsum/(ewald_g*ewald_g/2.0*xprd*yprd*zprd)
      e_long = coulpre*e_long

c compute long-range contribution to virial

      if (iflag.eq.1) then
        i = atompnt
        do ii = 1,nlocal
          virial(1) = virial(1) + q(i)*coulpre*ek(1,ii)*x(1,i)
          virial(2) = virial(2) + q(i)*coulpre*ek(2,ii)*x(2,i)
          virial(3) = virial(3) + q(i)*coulpre*ek(3,ii)*x(3,i)
          virial(4) = virial(4) + q(i)*coulpre*ek(2,ii)*x(1,i)
          virial(5) = virial(5) + q(i)*coulpre*ek(3,ii)*x(1,i)
          virial(6) = virial(6) + q(i)*coulpre*ek(3,ii)*x(2,i)
          i = list(i)
        enddo
      endif
 
      return
      end


c -------------------------------------------------------------------------
c check that all my particles and their stencils will fit in my 3-d brick
c  may require unwrapping of PBC
c store adjusted coords and charge in particle(4,nlocal) array
c store grid pt that is center of stencil in partgrid(3,nlocal) array

      subroutine particle_unmap
      include "lammps.h"
      include "mpif.h"

      delxinv = nx_pppm/xprd
      delyinv = ny_pppm/yprd
      delzinv = nz_pppm/zprd

      if (mod(orderflag,2).ne.0) then
        pshift = 0.5
      else
        pshift = 0.0
      endif

      icount = 0
      i = atompnt
      do ii = 1,nlocal

        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)

c (nx,ny,nz) = global coords of grid pt to "lower left" of charge
c original x(*,i) are assumed to be inside global box (values > {xyz}boundlo)
c this will have to change if don't enforce PBC at every timestep

        nx = int((xtmp-xboundlo)*delxinv+pshift)
        ny = int((ytmp-yboundlo)*delyinv+pshift)
        nz = int((ztmp-zboundlo)*delzinv+pshift)

c undo PBC if grid pt is outside my 3-d brick
c adjusted nx,ny,nz and {xyz}tmp may be outside global box

        if (nx.lt.nxlo_out) then
          nx = nx + nx_pppm
          xtmp = xtmp + xprd
        else if (nx.gt.nxhi_out) then
          nx = nx - nx_pppm
          xtmp = xtmp - xprd
        endif
        if (ny.lt.nylo_out) then
          ny = ny + ny_pppm
          ytmp = ytmp + yprd
        else if (ny.gt.nyhi_out) then
          ny = ny - ny_pppm
          ytmp = ytmp - yprd
        endif
        if (nz.lt.nzlo_out) then
          nz = nz + nz_pppm
          ztmp = ztmp + zprd
        else if (nz.gt.nzhi_out) then
          nz = nz - nz_pppm
          ztmp = ztmp - zprd
        endif

c insure that entire stencil around new nx,ny,nz will fit in my 3-d brick

        if (nx+nlower.lt.nxlo_out.or.nx+nupper.gt.nxhi_out.or.
     $       ny+nlower.lt.nylo_out.or.ny+nupper.gt.nyhi_out.or.
     $       nz+nlower.lt.nzlo_out.or.nz+nupper.gt.nzhi_out) then
          write (6,*) 'Out-of-range Atom:',node,x(1,i),x(2,i),x(3,i),
     $         nx,nxlo_out,nxhi_out,ny,nylo_out,nyhi_out,
     $         nz,nzlo_out,nzhi_out
          icount = icount + 1
        endif

c store adjusted particle grid pt and coords in (1:nlocal) arrays

        partgrid(1,ii) = nx
        partgrid(2,ii) = ny
        partgrid(3,ii) = nz

        particle(1,ii) = xtmp
        particle(2,ii) = ytmp
        particle(3,ii) = ztmp
        particle(4,ii) = q(i)

        i = list(i)
      enddo

      call mpi_allreduce(icount,jcount,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jcount.gt.0) call error('Cannot compute PPPM')

      return
      end


c -------------------------------------------------------------------------
c create discretized "density" on section of global grid due to my particles
c density(x,y,z) = charge "density" at grid points of my 3-d brick
c (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
c  in global grid

      subroutine make_rho(density)
      include "lammps.h"
      real*8 density(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)

      parameter (maxorder=16)
      dimension rho1d(-maxorder/2:maxorder/2,3)
      
      delxinv = nx_pppm/xprd
      delyinv = ny_pppm/yprd
      delzinv = nz_pppm/zprd
      delvolinv = delxinv*delyinv*delzinv

      if (mod(orderflag,2).ne.0) then
        qshift = 0.0
      else
        qshift = 0.5
      endif

c clear density array

      do k = nzlo_out,nzhi_out
        do j = nylo_out,nyhi_out
          do i = nxlo_out,nxhi_out
            density(i,j,k) = 0.0
          enddo
        enddo
      enddo

c loop over my charges, add their contribution to nearby grid points
c (nx,ny,nz) = global coords of grid pt to "lower left" of charge
c (dx,dy,dz) = distance to "lower left" grid pt
c (mx,my,mz) = global coords of moving stencil pt

      do i = 1,nlocal

        nx = partgrid(1,i)
        ny = partgrid(2,i)
        nz = partgrid(3,i)
        dx = nx+qshift - (particle(1,i)-xboundlo)*delxinv
        dy = ny+qshift - (particle(2,i)-yboundlo)*delyinv
        dz = nz+qshift - (particle(3,i)-zboundlo)*delzinv

        call calc_rho1d(dx,dy,dz,rho1d(-maxorder/2,1),orderflag)

        a0 = delvolinv * particle(4,i)
        do n = nlower,nupper
          mz = n+nz
          aa = a0*rho1d(n,3)
          do m = nlower,nupper
            my = m+ny
            aaa = aa*rho1d(m,2)
            do l = nlower,nupper
              mx = l+nx
              density(mx,my,mz) = density(mx,my,mz) + aaa*rho1d(l,1)
            enddo
          enddo
        enddo

      enddo

      return
      end


c -------------------------------------------------------------------------
c Poisson solver - FFT based

      subroutine poisson(vdx,vdy,vdz)
      include "lammps.h"
      include "mpif.h"
      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 vdy(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 vdz(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)

      complex*16 zim
      real*8 zero,minus_one

      zero = 0.0
      minus_one = -1.0
      zim = cmplx(zero,minus_one)

c transform charge density (r -> k) 

      do indx = 1,nfft
        workvec1(indx) = cmplx(density_fft(indx),zero)
      enddo

      call fft_3d(workvec1,workvec1,1,plan1_fft)

c scale by 1/total-grid-pts to get rho(k)
c multiply by Green's function to get V(k)
c e_long = portion of k-space potential energy on this proc's FFT grid

      scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm)
      e_long = 0.0
      do indx = 1,nfft
        workvec2(indx) = scaleinv*workvec1(indx)
        workvec1(indx) = greensfn(indx)*workvec2(indx)
        e_long = e_long + real(workvec1(indx)*conjg(workvec2(indx)))
      enddo
      e_long = 0.5*xprd*yprd*zprd*e_long

c compute gradients of V(r) in each of 3 dims by transformimg -ik*V(k)
c FFT leaves data in 3-d brick decomposition
c copy it into inner portion of vd_xyz arrays

c x direction gradient

      icount = 0
      do k = nzlo_fft,nzhi_fft
        do j = nylo_fft,nyhi_fft
          do i = nxlo_fft,nxhi_fft
            icount = icount + 1
            workvec2(icount) = zim*fkvecs(1,i)*workvec1(icount)
          enddo
        enddo
      enddo

      call fft_3d(workvec2,workvec2,-1,plan2_fft)

      icount = 0
      do iz = nzlo_in,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            vdx(ix,iy,iz) = real(workvec2(icount))
          enddo
        enddo
      enddo

c y direction gradient

      icount = 0
      do k = nzlo_fft,nzhi_fft
        do j = nylo_fft,nyhi_fft
          do i = nxlo_fft,nxhi_fft
            icount = icount + 1
            workvec2(icount) = zim*fkvecs(2,j)*workvec1(icount)
          enddo
        enddo
      enddo

      call fft_3d(workvec2,workvec2,-1,plan2_fft)

      icount = 0
      do iz = nzlo_in,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            vdy(ix,iy,iz) = real(workvec2(icount))
          enddo
        enddo
      enddo

c z direction gradient

      icount = 0
      do k = nzlo_fft,nzhi_fft
        do j = nylo_fft,nyhi_fft
          do i = nxlo_fft,nxhi_fft
            icount = icount + 1
            workvec2(icount) = zim*fkvecs(3,k)*workvec1(icount)
          enddo
        enddo
      enddo

      call fft_3d(workvec2,workvec2,-1,plan2_fft)

      icount = 0
      do iz = nzlo_in,nzhi_in
        do iy = nylo_in,nyhi_in
          do ix = nxlo_in,nxhi_in
            icount = icount + 1
            vdz(ix,iy,iz) = real(workvec2(icount))
          enddo
        enddo
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c interpolate grid values to get electric field acting on my particles

      subroutine electric_field(vdx,vdy,vdz)
      include "lammps.h"
      real*8 vdx(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 vdy(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      real*8 vdz(nxlo_out:nxhi_out,nylo_out:nyhi_out,
     $     nzlo_out:nzhi_out)
      parameter (maxorder=16)
      dimension rho1d(-maxorder/2:maxorder/2,3)

      delxinv = nx_pppm/xprd
      delyinv = ny_pppm/yprd
      delzinv = nz_pppm/zprd

      if (mod(orderflag,2).ne.0) then
        qshift = 0.0
      else
        qshift = 0.5
      endif
      
c clear e-fields

      do i = 1,nlocal
        ek(1,i) = 0.0
        ek(2,i) = 0.0
        ek(3,i) = 0.0
      enddo

c loop over my charges, interpolate electric field from nearby grid points
c (nx,ny,nz) = global coords of grid pt to "lower left" of charge
c (dx,dy,dz) = distance to "lower left" grid pt
c (mx,my,mz) = global coords of moving stencil pt

      do i = 1,nlocal

        nx = partgrid(1,i)
        ny = partgrid(2,i)
        nz = partgrid(3,i)
        dx = nx+qshift - (particle(1,i)-xboundlo)*delxinv
        dy = ny+qshift - (particle(2,i)-yboundlo)*delyinv
        dz = nz+qshift - (particle(3,i)-zboundlo)*delzinv
        call calc_rho1d(dx,dy,dz,rho1d(-maxorder/2,1),orderflag)

        do n = nlower,nupper
          mz = n+nz
          z0 = rho1d(n,3)
          do m = nlower,nupper
            my = m+ny
            y0 = z0*rho1d(m,2)
            do l = nlower,nupper
              mx = l+nx
              x0 = y0*rho1d(l,1)
              ek(1,i) = ek(1,i) - x0*vdx(mx,my,mz)
              ek(2,i) = ek(2,i) - x0*vdy(mx,my,mz)
              ek(3,i) = ek(3,i) - x0*vdz(mx,my,mz)
            enddo
          enddo
        enddo

      enddo

      return
      end


c -------------------------------------------------------------------------
c cacl_rho1d - calculate 1d assignment functions
c   INPUT:
c         dx,dy,dz = distance of particle  from "lower, left" grid point
c         n ======== order of assignment scheme
c   OUTPUT:
c         w = assignment functions for particle (3*n values)

      subroutine calc_rho1d(dx,dy,dz,w,n)
      implicit real*8 (a-h,o-z)
      include 'param.h'
      parameter (maxorder=16)
      dimension a(0:maxorder-1,-maxorder:maxorder),
     $     acoef(0:maxorder-1,(1-maxorder)/2:maxorder/2),
     $     w(-maxorder/2:maxorder/2,3)
      save acoef
      data n_old /0/
      save n_old

c generate the coeffients for the weight function of order n
c              (n-1)
c  Wn(x) =     Sum    wn(k,x) , Sum is over every other integer
c           k=-(n-1)
c  For k=-(n-1),-(n-1)+2, ....., (n-1)-2,n-1
c  Note k is odd integers if n is even
c      and  even integers if n is old
c              ---
c             | n-1
c             | Sum a(l,j)*(x-k/2)**l   if abs(x-k/2) < 1/2
c  wn(k,x) = <  l=0
c             |
c             |  0                       otherwise
c              ---
c  The a coeffients are packed into the array acoef to eliminate zeros
c  acoef(l,((k+mod(n+1,2))/2) = a(l,k)

      if (n.ne.n_old) then
        do k = -n,n
          do l = 0,n-1
            a(l,k) = 0.0d0
          enddo
        enddo
        
        a(0,0) = 1.0d0
        do n = 1,n-1
          do k = -n,n,2
            s = 0.0d0
            do l = 0,n-1
              a(l+1,k) = (a(l,k+1)-a(l,k-1))/dfloat(l+1)
              s = s+0.5d0**(l+1)*(a(l,k-1)+(-1.0d0)**l*a(l,k+1))/
     $             dfloat(l+1)
            enddo
            a(0,k) = s
          enddo
        enddo
        
        kk = (1-n)/2
        do k = -(n-1),n-1,2
          do l = 0,n-1
            acoef(l,kk) = a(l,k)
          enddo
          kk = kk + 1
        enddo
        n_old = n
      endif

c if acoef is already calculated then just form the {w}s

      do k = (1-n)/2,n/2
        w(k,1)  =   0.0
        w(k,2)  =   0.0
        w(k,3)  =   0.0
        do l = n-1,0,-1
          w(k,1) = acoef(l,k) + w(k,1)*dx
          w(k,2) = acoef(l,k) + w(k,2)*dy
          w(k,3) = acoef(l,k) + w(k,3)*dz
        enddo
      enddo
      
      return
      end
