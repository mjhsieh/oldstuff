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
c setup of 2 decompositions (3-d brick and 2-d FFT) and various coeffs
c  for PPPM calculation
c iflag = 0 section only need be computed once since those quantities
c  depend only on grid size
c NOTE: if box volume changes enough that a different resolution grid 
c  is needed then iflag = 0 section should be recomputed
c iflag = 1 section must be recomputed whenever box size changes

      subroutine pppm_coeff(iflag)
      include "lammps.h"
      include "mpif.h"
      integer iflag

      integer istatus(mpi_status_size)
      logical factorable
      integer nfacs,factors(3)
c comment one of next 2 lines out
c 1st is for power-of-2 FFTs
c 2nd is for non-power-of-2 FFTs
c      data nfacs /1/, factors /2,0,0/
      data nfacs /3/, factors /2,3,5/

      if (iflag.eq.0) then

c PPPM grid size and cutoff parameter from estimation routine

        call get_p3m_parms(orderflag,long_prec,cutcoul,
     $       xprd,yprd,zprd,ewald_g,nx_pppm,ny_pppm,nz_pppm)

        gewsr = ewald_g/sqrt(2.0)

        if (node.eq.0) then
          write (6,*) 'PPPM g =',ewald_g
          write (6,*) 'Needed PPPM nx,ny,nz =',
     $         nx_pppm,ny_pppm,nz_pppm
          write (1,*) 'PPPM g =',ewald_g
          write (1,*) 'Needed PPPM nx,ny,nz =',
     $         nx_pppm,ny_pppm,nz_pppm
        endif

c convert PPPM grid size into "good" numbers - factorable by 2 or 2,3,5

        do while (.not.factorable(nx_pppm,nfacs,factors))
          nx_pppm = nx_pppm + 1
        enddo
        do while (.not.factorable(ny_pppm,nfacs,factors))
          ny_pppm = ny_pppm + 1
        enddo
        do while (.not.factorable(nz_pppm,nfacs,factors))
          nz_pppm = nz_pppm + 1
        enddo

c overwrite PPPM grid size with input values

        if (meshflag.eq.1) then
          nx_pppm = nx_pppm_input
          ny_pppm = ny_pppm_input
          nz_pppm = nz_pppm_input
        endif
        
c print-out actual PPPM grid to be used
        
        if (node.eq.0) then
          write (6,*) 'Actual PPPM nx,ny,nz = ',
     $         nx_pppm,ny_pppm,nz_pppm
          write (1,*) 'Actual PPPM nx,ny,nz = ',
     $         nx_pppm,ny_pppm,nz_pppm
        endif

c nlo_in,nhi_in = lower/upper limits of the 3-d sub-brick of
c   global PPPM grid that I own in each dimension WITHOUT ghost cells
c indices range from 0 to N-1

        nxlo_in = me(1)*nx_pppm/pgrid(1)
        nxhi_in = (me(1)+1)*nx_pppm/pgrid(1) - 1
        nylo_in = me(2)*ny_pppm/pgrid(2)
        nyhi_in = (me(2)+1)*ny_pppm/pgrid(2) - 1
        nzlo_in = me(3)*nz_pppm/pgrid(3)
        nzhi_in = (me(3)+1)*nz_pppm/pgrid(3) - 1

c nlower,nupper = stencil size for mapping particles to PPPM grid

        nlower = -(orderflag-1)/2
        nupper = orderflag/2

c nlo_out,nhi_out = lower/upper limits of the 3-d sub-brick of
c   global PPPM grid that I own in each dimension WITH ghost cells
c ghost 3-d brick boundaries = inner + stencil + particles moving skin/2.0
c nlo/nhi = global coords of grid pt to "lower left" of smallest/largest
c           position a particle in my box can be at
c add/subtract 4096 to avoid problem of int(-0.75) = 0 when it needs to be -1

        delxinv = nx_pppm/xprd
        delyinv = ny_pppm/yprd
        delzinv = nz_pppm/zprd

        if (mod(orderflag,2).ne.0) then
          pshift = 4096.5
        else
          pshift = 4096.0
        endif

        cuthalf = skin/2.0

        nlo = int((border(1,1)-cuthalf-xboundlo)*delxinv+pshift) - 4096
        nhi = int((border(2,1)+cuthalf-xboundlo)*delxinv+pshift) - 4096
        nxlo_out = nlo + nlower
        nxhi_out = nhi + nupper

        nlo = int((border(1,2)-cuthalf-yboundlo)*delyinv+pshift) - 4096
        nhi = int((border(2,2)+cuthalf-yboundlo)*delyinv+pshift) - 4096
        nylo_out = nlo + nlower
        nyhi_out = nhi + nupper

        nlo = int((border(1,3)-cuthalf-zboundlo)*delzinv+pshift) - 4096
        nhi = int((border(2,3)+cuthalf-zboundlo)*delzinv+pshift) - 4096
        nzlo_out = nlo + nlower
        nzhi_out = nhi + nupper

c check if any 3-d brick w/ ghosts has too many mesh points > maxgrid
        
        itotal = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
     $       (nzhi_out-nzlo_out+1)
        call mpi_allreduce(itotal,ntotal,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)
        if (ntotal.gt.maxgrid)
     $       call error('PPPM full mesh too large - boost maxgrid')

c check if any 3-d brick w/out ghosts has too many mesh points > maxfft
        
        itotal = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
     $       (nzhi_in-nzlo_in+1)
        call mpi_allreduce(itotal,ntotal,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)
        if (ntotal.gt.maxfft)
     $       call error('PPPM mesh too large - boost maxfft')

c exchange messages to determine how many planes I will send,recv in each dir
c if no neighbor proc exists, values comes from self
c   since I have ghosts regardless

        nplanes = nxlo_in - nxlo_out
        if (mpart(1,1).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,1),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nxhi_ghost,1,mpi_integer,mpart(2,1),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nxhi_ghost = nplanes
        endif

        nplanes = nxhi_out - nxhi_in
        if (mpart(2,1).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,1),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nxlo_ghost,1,mpi_integer,mpart(1,1),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nxlo_ghost = nplanes
        endif

        nplanes = nylo_in - nylo_out
        if (mpart(1,2).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,2),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nyhi_ghost,1,mpi_integer,mpart(2,2),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nyhi_ghost = nplanes
        endif

        nplanes = nyhi_out - nyhi_in
        if (mpart(2,2).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,2),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nylo_ghost,1,mpi_integer,mpart(1,2),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nylo_ghost = nplanes
        endif

        nplanes = nzlo_in - nzlo_out
        if (mpart(1,3).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,3),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nzhi_ghost,1,mpi_integer,mpart(2,3),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nzhi_ghost = nplanes
        endif

        nplanes = nzhi_out - nzhi_in
        if (mpart(2,3).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,3),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nzlo_ghost,1,mpi_integer,mpart(1,3),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nzlo_ghost = nplanes
        endif

c test that ghost overlap is not bigger than my sub-domain

        jflag = 0
        if (nxlo_ghost.gt.nxhi_in-nxlo_in+1) jflag = 1
        if (nxhi_ghost.gt.nxhi_in-nxlo_in+1) jflag = 1
        if (nylo_ghost.gt.nyhi_in-nylo_in+1) jflag = 1
        if (nyhi_ghost.gt.nyhi_in-nylo_in+1) jflag = 1
        if (nzlo_ghost.gt.nzhi_in-nzlo_in+1) jflag = 1
        if (nzhi_ghost.gt.nzhi_in-nzlo_in+1) jflag = 1

        call mpi_allreduce(jflag,kflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)

        if (kflag.gt.0)
     $       call error('PPPM stencil extends beyond neighbor proc'//
     $       ' - reduce PPPM order')

c bounds check on buffer usage in brick2fft and fillbrick
c idel = max # of ghost planes to send or recv in +/- direction of each dim
c 3 factor in itotal is for 3 components of vd_xyz in fillbrick

        nx = nxhi_out - nxlo_out + 1
        ny = nyhi_out - nylo_out + 1
        nz = nzhi_out - nzlo_out + 1

        idelx = max(nxlo_ghost,nxhi_ghost)
        idelx = max(idelx,nxhi_out-nxhi_in)
        idelx = max(idelx,nxlo_in-nxlo_out)

        idely = max(nylo_ghost,nyhi_ghost)
        idely = max(idely,nyhi_out-nyhi_in)
        idely = max(idely,nylo_in-nylo_out)

        idelz = max(nzlo_ghost,nzhi_ghost)
        idelz = max(idelz,nzhi_out-nzhi_in)
        idelz = max(idelz,nzlo_in-nzlo_out)

        nxx = idelx * ny * nz
        nyy = idely * nx * nz
        nzz = idelz * nx * ny

        itotal = max(nxx,nyy)
        itotal = max(itotal,nzz)
        itotal = 3*itotal

        call mpi_allreduce(itotal,ntotal,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)

        if (ntotal.gt.maxexchtot)
     $       call error('PPPM remap needs buffer space'//
     $       ' - boost maxexch')

c decomposition of FFT mesh
c proc owns entire x-dimension, clump of columns in y,z dimensions
c npey_fft,npez_fft = # of procs in y,z dims
c if nprocs is small enough, proc can own 1 or more entire xy planes,
c   else proc owns 2-d sub-blocks of yz plane
c me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
c nlo_fft,nhi_fft = lower/upper limit of the section
c   of the global FFT mesh that I own in each dimension
c indices range from 0 to N-1

        if (nz_pppm.ge.nprocs) then
          npey_fft = 1
          npez_fft = nprocs
        else
          call proc2grid2d(nprocs,ny_pppm,nz_pppm,npey_fft,npez_fft)
        endif

        me_y = mod(node,npey_fft)
        me_z = node/npey_fft
        
        nxlo_fft = 0
        nxhi_fft = nx_pppm - 1
        nylo_fft = me_y*ny_pppm/npey_fft
        nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1
        nzlo_fft = me_z*nz_pppm/npez_fft
        nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1

c check if total # of FFT points on any proc is too big

        nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
     $       (nzhi_fft-nzlo_fft+1)
        call mpi_allreduce(nfft,ntotal,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)
        if (ntotal.gt.maxfft)
     $       call error('PPPM FFT mesh too large - boost maxfft')
        
c create 2 FFT plans and remap plan
c 1st FFT plan keeps data in FFT decompostion
c 2nd FFT plan returns data in 3-d brick decomposition
c remap plan takes data from 3-d brick to FFT decomposition
c all indices must be converted to range from 1 to N

        call fft_3d_create_plan(mpi_comm_world,
     $       nx_pppm,ny_pppm,nz_pppm,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       0,0,itmp,plan1_fft)

        call fft_3d_create_plan(mpi_comm_world,
     $       nx_pppm,ny_pppm,nz_pppm,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       nxlo_in+1,nxhi_in+1,nylo_in+1,nyhi_in+1,
     $       nzlo_in+1,nzhi_in+1,
     $       0,0,itmp,plan2_fft)

        call remap_3d_create_plan(mpi_comm_world,
     $       nxlo_in+1,nxhi_in+1,nylo_in+1,nyhi_in+1,
     $       nzlo_in+1,nzhi_in+1,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       1,0,0,2,plan_remap)

      endif

c pre-compute fkvecs(3,*) for my FFT grid pts

      unitk = (2.0*3.141592654/xprd)
      do k = nxlo_fft,nxhi_fft
        kper = k - nx_pppm*int(2*k/nx_pppm)
        fkvecs(1,k) = unitk*kper
      enddo
      
      unitk = (2.0*3.141592654/yprd)
      do l = nylo_fft,nyhi_fft
        lper = l - ny_pppm*int(2*l/ny_pppm)
        fkvecs(2,l) = unitk*lper
      enddo
      
      unitk = (2.0*3.141592654/zprd)
      do m = nzlo_fft,nzhi_fft
        mper = m - nz_pppm*int(2*m/nz_pppm)
        fkvecs(3,m) = unitk*mper
      enddo

c pre-compute Green's function coeffs

      call make_greensfn(3,orderflag,ewald_g,
     $     nx_pppm,ny_pppm,nz_pppm,xprd,yprd,zprd,
     $     nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
     $     greensfn)

      return
      end


c -------------------------------------------------------------------------
c get p3m parameters
c   INPUT:
c           na_ord---------------------p3m assignment scheme order
c           error----------------------desired relative error in forces
c           rcut-----------------------real space cutoff 
c           cellx,celly,cellz----------edge lengths of periodic cell
c   OUTPUT:
c           gew-------------Ewald parameter 
c           nx_pppm,ny_pppm,nz_pppm-----P3m grid sizes
c
c   gew_rc returns gew*rc for a desired real space forces error in P3m
c   gew_delta returns required gew*delta for given error in k-space
c             part of forces (where delta = grid spacing)

      subroutine get_p3m_parms(na_ord,rerror,rcut,
     $     cellx,celly,cellz,gew,nx_pppm,ny_pppm,nz_pppm)
      implicit real*8 (a-h,o-z)
      parameter (maxorder=16)
      dimension a(maxorder),b(maxorder)
      data a /5.0,5.78,6.93,7.76,8.55,9.10,9.66,10.17,10.60,
     $     11.1,14,15,16,17,18,19/
      data b /2.0,2.25,3.41,4.62,5.85,7.05,8.20, 9.28,10.34,
     $     11.3,12,13,14,15,16,17/
      data aek /2.5d0/

c implicit functions

      gew_delta(narg,arg) = (arg/(aek*exp(-a(narg))))**(1.0/b(narg))
      gew_rc(arg) = sqrt(2.0)*(1.35-.15*log(arg))

c error check

      if (na_ord.gt.maxorder)
     $     call error('PPPM order cannot be larger than maxorder')

c compute necessary gew
c based on desired error and real space cutoff
	
      gewrc = gew_rc(rerror)
      gew = gewrc/rcut

c compute optimal nx_pppm,ny_pppm,nz_pppm
c based on assignment order and desired error

      gewdelta = gew_delta(na_ord,rerror)
      nx_pppm = gew*cellx/gewdelta
      ny_pppm = gew*celly/gewdelta
      nz_pppm = gew*cellz/gewdelta

      return
      end


c -------------------------------------------------------------------------
c set up modified (Hockney-Eastwood) Coulomb Green's Fn
c        (replaces 4*pi/k**2) for k-vectors on this PE
c     INPUT:
c         nper_dim=3 (fully periodic), 1 (slab), 0 (cluster)
c         na_ord=order of assignment scheme
c         gew=(G-Ewald) width parameter for Gaussians
c         nx_pppm,ny_pppm,nz_pppm=grid array dimensions
c         cellz,celly,cellz=grid sizes
c         nx1,nx2 = lower, upper x limit of k-space arrays for this PE
c         ny1,ny2 = lower, upper y limit of k-space arrays for this PE
c         nz1,nz2 = lower, upper z limit of k-space arrays for this PE
c     OUTPUT:
c         greensfn=Coulomb Green's Fn for k-vectors on this PE
c     DATA:
c         epshoc=convergence parameter for sums 
c     CALLS: gf_denom

      subroutine make_greensfn(nper_dim,na_ord,gew,nx_pppm,ny_pppm,
     $     nz_pppm,cellx,celly,cellz,nx1,nx2,ny1,ny2,nz1,nz2,greensfn)
      implicit real*8 (a-h,o-z)
      dimension greensfn(*)
      data epshoc /1.e-7/

c implicit function

      shpfn(arg) = exp(-.25*(arg/gew)**2)

      unitkx = (2.*3.141592654/cellx)
      unitky = (2.*3.141592654/celly)
      unitkz = (2.*3.141592654/cellz)

      nbx = (gew*cellx/(3.141592654*nx_pppm))*((-log(epshoc))**.25)
      nby = (gew*cellx/(3.141592654*ny_pppm))*((-log(epshoc))**.25)
      nbz = (gew*cellz/(3.141592654*nz_pppm))*((-log(epshoc))**.25)

      form = 1.0
      indx = 1
      do m = nz1,nz2
        mper = m-nz_pppm*int(2*m/nz_pppm)
        snz2 = (sin(.5*unitkz*mper*cellz/nz_pppm))**2
        do l = ny1,ny2
          lper = l-ny_pppm*int(2*l/ny_pppm)
          sny2 = (sin(.5*unitky*lper*celly/ny_pppm))**2
          do k = nx1,nx2

            kper = k-nx_pppm*int(2*k/nx_pppm)
            snx2 = (sin(.5*unitkx*kper*cellx/nx_pppm))**2
            sqk = ((unitkx*kper)**2+(unitky*lper)**2+(unitkz*mper)**2)
            sqk_perp = ((unitkx*kper)**2+(unitky*lper)**2)

            if (sqk.ne.0.0) then
c              if (nper_dim.eq.0) form = 1.0-cos(sqrt(sqk)*rcut)
c              if (nper_dim.eq.2) then
c                perpk = sqrt(sqk_perp)
c                expcut = exp(-perpk*zcut)
c                argz = unitkz*mper*zcut
c                f1 = cos(argz)
c                f2 = 0.0
c                if (sqk_perp.ne.0.0)
c     $               f2 = unitkz*mper*sin(argz)/perpk
c                form = 1.0d0-expcut*(f1-f2)
c              endif

              greensfn(indx) = form*12.5663706/sqk
              sum1 = 0.0
              do nx = -nbx,nbx
                qx = unitkx*(kper+nx_pppm*nx)
                sx = shpfn(qx)
                wx = 1.0
                argx = .5*qx*cellx/nx_pppm
                if (argx.ne.0.0) wx = (sin(argx)/argx)**na_ord
                do ny = -nby,nby
                  qy = unitky*(lper+ny_pppm*ny)
                  sy = shpfn(qy)
                  wy = 1.0
                  argy = .5*qy*celly/ny_pppm
                  if (argy.ne.0.0) wy = (sin(argy)/argy)**na_ord
                  do nz = -nbz,nbz
                    qz = unitkz*(mper+nz_pppm*nz)
                    sz = shpfn(qz)
                    wz = 1.0
                    argz = .5*qz*cellz/nz_pppm
                    if (argz.ne.0.0) wz = (sin(argz)/argz)**na_ord

                    dot1 = unitkx*kper*qx +
     $                   unitky*lper*qy+unitkz*mper*qz
                    dot2 = qx*qx+qy*qy+qz*qz
                    sum1 = sum1+(dot1/dot2)*(wx*sx*wy*sy*wz*sz)**2
                  enddo
                enddo
              enddo
              call gf_denom(na_ord,snx2,sny2,snz2,sum2sq)
              greensfn(indx) = greensfn(indx)*sum1/sum2sq
            else
              greensfn(indx) = 0.0
            endif

            indx = indx+1
          enddo
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c returns denominator for Hockney-Eastwood Green's function
c               inf                 n-1
c      S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
c              j=-inf               l=0
c
c       = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z=sin(x)
c     INPUT:
c           n ====== order of assignment scheme
c           x,y,z=== sin(kx*deltax/2), etc.
c           b ====== denominator expansion coefficients /Gamma(2n)
c     OUTPUT:
c           s ====== denominator of Hockney-Eastwood expression

      subroutine gf_denom(n,x,y,z,s)
      implicit real*8 (a-h,o-z)
      parameter (maxorder=16)
      dimension b(0:maxorder-1)
      save b
      data nlast /0/
      save nlast

      if (n.ne.nlast) then
        ifact = 1
        do k = 1,2*n-1
          ifact = ifact*k
        enddo
        gaminv = 1./float(ifact)
        do l = 0,n-1
          b(l) = 0.0
        enddo
        b(0) = 1.0
        
        do m = 1,n-1
          do l = m,1,-1
            b(l) = 4.0*(b(l)*(l-m)*(l-m-.5)-b(l-1)*(l-1-m)**2)
          enddo
          b(0) = 4.0*(b(0)*(l-m)*(l-m-.5))
        enddo
        do l = 0,n-1
          b(l) = b(l)*gaminv
        enddo
        nlast = n
      endif
      
      sx = 0.0
      sy = 0.0
      sz = 0.0
      do l = n-1,0,-1
        sx = b(l) + sx*x
        sy = b(l) + sy*y
        sz = b(l) + sz*z
      enddo
      s = sx*sy*sz
      s = s*s
      
      return
      end


c -------------------------------------------------------------------------
c check if all the factors of n are in factors(1:nfacs)
c return TRUE if they are
c return FALSE if they are not

      logical function factorable(n,nfacs,factors)
      integer n,nfacs,factors(nfacs)

      integer nlist,list(20)

c build list of factors of n
c nlist = # of factors
c after inner while loop, i is a factor of n
c nremain = remaining portion of n that has not been factored

      nlist = 0
      nremain = n

      do while (nremain.gt.1)
        i = 2
        do while (mod(nremain,i).ne.0)
          i = i + 1
        enddo
        nlist = nlist + 1
        list(nlist) = i
        nremain = nremain/i
      enddo

c check if every factor of n in list is also is in input factors
c as soon as one list item is not in factors, return FALSE

      do i = 1,nlist
        iflag = 0
        do j = 1,nfacs
          if (list(i).eq.factors(j)) iflag = 1
        enddo
        if (iflag.eq.0) then
          factorable = .FALSE.
          return
        endif
      enddo

      factorable = .TRUE.

      return
      end
