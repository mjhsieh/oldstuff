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
c Ewald routines authored by Roy Pollock (modified for LAMMPS by SJP)
c   (510) 422-4088, pollockr@llnl.gov
c   L-299, LLNL, PO Box 808, Livermore, CA  94550

c -------------------------------------------------------------------------
c compute k-space part of Ewald forces and potential

      subroutine ewald(iflag)
      include "lammps.h"
      include "mpif.h"
      integer iflag

c compute partial structure factors by summing over particles on processor

      call eikdotr

c compute total structure factor by summing over all processors

      call mpi_allreduce(sfacrl,sfacrl_all,kcount,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(sfacim,sfacim_all,kcount,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

c compute k-space part of electrostatic potential energy

      pi = 4.0*atan(1.0)
      gewsq = ewald_g*ewald_g

      e_long = 0.0
      do kk = 1,kcount
        e_long = e_long + ug(kk)*(sfacrl_all(kk)**2+sfacim_all(kk)**2)
      enddo
      e_long = e_long - ewald_g*qsqsum/1.772453851 -
     $     0.5*pi*qsum*qsum/(gewsq*xprd*yprd*zprd)
      e_long = coulpre*e_long

c compute k-space part of electric field for particles on processor

      i = atompnt
      do ii = 1,nlocal
        ek(1,i) = 0.0
        ek(2,i) = 0.0
        ek(3,i) = 0.0
        i = list(i)
      enddo

      do kk = 1,kcount
        kx = kxvecs(kk)
        ky = kyvecs(kk)
        kz = kzvecs(kk)
        i = atompnt
        do ii = 1,nlocal
          cypz = cs(ky,2,i)*cs(kz,3,i) - sn(ky,2,i)*sn(kz,3,i)
          sypz = sn(ky,2,i)*cs(kz,3,i) + cs(ky,2,i)*sn(kz,3,i)
          exprl = cs(kx,1,i)*cypz - sn(kx,1,i)*sypz
          expim = sn(kx,1,i)*cypz + cs(kx,1,i)*sypz
          foo = expim*sfacrl_all(kk) - exprl*sfacim_all(kk)
          ek(1,i) = ek(1,i) + eg(1,kk)*foo
          ek(2,i) = ek(2,i) + eg(2,kk)*foo
          ek(3,i) = ek(3,i) + eg(3,kk)*foo
          i = list(i)
        enddo
      enddo

c convert e-field to force

      i = atompnt
      do ii = 1,nlocal
        f(1,i) = f(1,i) + q(i)*coulpre*ek(1,i)
        f(2,i) = f(2,i) + q(i)*coulpre*ek(2,i)
        f(3,i) = f(3,i) + q(i)*coulpre*ek(3,i)
        i = list(i)
      enddo

c compute virial contribution from Ewald

      if (iflag.eq.1) then
        i = atompnt
        do ii = 1,nlocal
          virial(1) = virial(1) + q(i)*coulpre*ek(1,i)*x(1,i)
          virial(2) = virial(2) + q(i)*coulpre*ek(2,i)*x(2,i)
          virial(3) = virial(3) + q(i)*coulpre*ek(3,i)*x(3,i)
          virial(4) = virial(4) + q(i)*coulpre*ek(2,i)*x(1,i)
          virial(5) = virial(5) + q(i)*coulpre*ek(3,i)*x(1,i)
          virial(6) = virial(6) + q(i)*coulpre*ek(3,i)*x(2,i)
          i = list(i)
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c calculate structure factors and exp(k*rr(i)) 

      subroutine eikdotr
      include "lammps.h"

      real*8 unitk(3)

      pi = 4.0*atan(1.0)
      unitk(1) = 2.0*pi/xprd
      unitk(2) = 2.0*pi/yprd
      unitk(3) = 2.0*pi/zprd
      gewsq = ewald_g*ewald_g
      gsqmx = -4.0*gewsq*log(long_prec)

      kcount = 0
c
c    (k,0,0),(0.l.0).(0,0,m)
c
      do ic = 1,3
        sqk = unitk(ic)**2
        if (sqk.le.gsqmx) then
          cstr1 = 0.0
          sstr1 = 0.0
          i = atompnt
          do ii = 1,nlocal
            cs(0,ic,i) = 1.0
            sn(0,ic,i) = 0.0
            cs(1,ic,i) = cos(unitk(ic)*x(ic,i))
            sn(1,ic,i) = sin(unitk(ic)*x(ic,i))
            cs(-1,ic,i) = cs(1,ic,i)
            sn(-1,ic,i) = -sn(1,ic,i)
            cstr1 = cstr1 + q(i)*cs(1,ic,i)
            sstr1 = sstr1 + q(i)*sn(1,ic,i)
            i = list(i)
          enddo
          kcount = kcount + 1
          sfacrl(kcount) = cstr1
          sfacim(kcount) = sstr1
        endif
      enddo

      do m = 2,kmax
        do ic = 1,3
          sqk = (m*unitk(ic))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            i = atompnt
            do ii = 1,nlocal
              cs(m,ic,i) = cs(m-1,ic,i)*cs(1,ic,i) -
     $             sn(m-1,ic,i)*sn(1,ic,i)
              sn(m,ic,i) = sn(m-1,ic,i)*cs(1,ic,i) +
     $             cs(m-1,ic,i)*sn(1,ic,i)
              cs(-m,ic,i) = cs(m,ic,i)
              sn(-m,ic,i) = -sn(m,ic,i)
              cstr1 = cstr1 + q(i)*cs(m,ic,i)
              sstr1 = sstr1 + q(i)*sn(m,ic,i)
              i = list(i)
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1
          endif
        enddo
      enddo
c
c   1 = (k,l,0), 2 = (k,-l,0)
c
      do k = 1,kmax
        do l = 1,kmax
          sqk = (k*unitk(1))**2 + (l*unitk(2))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            cstr2 = 0.0
            sstr2 = 0.0
            i = atompnt
            do ii = 1,nlocal
              cstr1 = cstr1 + q(i)*(cs(k,1,i)*cs(l,2,i) -
     $             sn(k,1,i)*sn(l,2,i))
              sstr1 = sstr1 + q(i)*(sn(k,1,i)*cs(l,2,i) +
     $             cs(k,1,i)*sn(l,2,i))
              cstr2 = cstr2 + q(i)*(cs(k,1,i)*cs(l,2,i) +
     $             sn(k,1,i)*sn(l,2,i))
              sstr2 = sstr2 + q(i)*(sn(k,1,i)*cs(l,2,i) -
     $             cs(k,1,i)*sn(l,2,i))
              i = list(i)
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1

            kcount = kcount + 1
            sfacrl(kcount) = cstr2
            sfacim(kcount) = sstr2
          endif
        enddo
      enddo
c
c   1 = (0,l,m), 2 = (0,l,-m)
c
      do l = 1,kmax
        do m = 1,kmax
          sqk = (l*unitk(2))**2 + (m*unitk(3))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            cstr2 = 0.0
            sstr2 = 0.0
            i = atompnt
            do ii = 1,nlocal
              cstr1 = cstr1 + q(i)*(cs(l,2,i)*cs(m,3,i) -
     $             sn(l,2,i)*sn(m,3,i))
              sstr1 = sstr1 + q(i)*(sn(l,2,i)*cs(m,3,i) +
     $             cs(l,2,i)*sn(m,3,i))
              cstr2 = cstr2 + q(i)*(cs(l,2,i)*cs(m,3,i) +
     $             sn(l,2,i)*sn(m,3,i))
              sstr2 = sstr2 + q(i)*(sn(l,2,i)*cs(m,3,i) -
     $             cs(l,2,i)*sn(m,3,i))
              i = list(i)
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1

            kcount = kcount + 1
            sfacrl(kcount) = cstr2
            sfacim(kcount) = sstr2
          endif
        enddo
      enddo
c
c   1 = (k,0,m), 2 = (k,0,-m)
c
      do k = 1,kmax
        do m = 1,kmax
          sqk = (k*unitk(1))**2 + (m*unitk(3))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            cstr2 = 0.0
            sstr2 = 0.0
            i = atompnt
            do ii = 1,nlocal
              cstr1 = cstr1 + q(i)*(cs(k,1,i)*cs(m,3,i) -
     $             sn(k,1,i)*sn(m,3,i))
              sstr1 = sstr1 + q(i)*(sn(k,1,i)*cs(m,3,i) +
     $             cs(k,1,i)*sn(m,3,i))
              cstr2 = cstr2 + q(i)*(cs(k,1,i)*cs(m,3,i) +
     $             sn(k,1,i)*sn(m,3,i))
              sstr2 = sstr2 + q(i)*(sn(k,1,i)*cs(m,3,i) -
     $             cs(k,1,i)*sn(m,3,i))
              i = list(i)
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1

            kcount = kcount + 1
            sfacrl(kcount) = cstr2
            sfacim(kcount) = sstr2
          endif
        enddo
      enddo
c
c   1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
c
      do k = 1,kmax
        do l = 1,kmax
          do m = 1,kmax
            sqk = (k*unitk(1))**2 + (l*unitk(2))**2 + (m*unitk(3))**2
            if (sqk.le.gsqmx) then
              cstr1 = 0.0
              sstr1 = 0.0
              cstr2 = 0.0
              sstr2 = 0.0
              cstr3 = 0.0
              sstr3 = 0.0
              cstr4 = 0.0
              sstr4 = 0.0
              i = atompnt
              do ii = 1,nlocal
                clpm = cs(l,2,i)*cs(m,3,i) - sn(l,2,i)*sn(m,3,i)
                slpm = sn(l,2,i)*cs(m,3,i) + cs(l,2,i)*sn(m,3,i)
                cstr1 = cstr1 + q(i)*(cs(k,1,i)*clpm -
     $               sn(k,1,i)*slpm)
                sstr1 = sstr1 + q(i)*(sn(k,1,i)*clpm +
     $               cs(k,1,i)*slpm)

                clpm = cs(l,2,i)*cs(m,3,i) + sn(l,2,i)*sn(m,3,i)
                slpm = -sn(l,2,i)*cs(m,3,i) + cs(l,2,i)*sn(m,3,i)
                cstr2 = cstr2 + q(i)*(cs(k,1,i)*clpm -
     $               sn(k,1,i)*slpm)
                sstr2 = sstr2 + q(i)*(sn(k,1,i)*clpm +
     $               cs(k,1,i)*slpm)

                clpm = cs(l,2,i)*cs(m,3,i) + sn(l,2,i)*sn(m,3,i)
                slpm = sn(l,2,i)*cs(m,3,i) - cs(l,2,i)*sn(m,3,i)
                cstr3 = cstr3 + q(i)*(cs(k,1,i)*clpm -
     $               sn(k,1,i)*slpm)
                sstr3 = sstr3 + q(i)*(sn(k,1,i)*clpm +
     $               cs(k,1,i)*slpm)

                clpm = cs(l,2,i)*cs(m,3,i) - sn(l,2,i)*sn(m,3,i)
                slpm = -sn(l,2,i)*cs(m,3,i) - cs(l,2,i)*sn(m,3,i)
                cstr4 = cstr4 + q(i)*(cs(k,1,i)*clpm -
     $               sn(k,1,i)*slpm)
                sstr4 = sstr4 + q(i)*(sn(k,1,i)*clpm +
     $               cs(k,1,i)*slpm)
                i = list(i)
              enddo
              kcount = kcount + 1
              sfacrl(kcount) = cstr1
              sfacim(kcount) = sstr1

              kcount = kcount + 1
              sfacrl(kcount) = cstr2
              sfacim(kcount) = sstr2

              kcount = kcount + 1
              sfacrl(kcount) = cstr3
              sfacim(kcount) = sstr3

              kcount = kcount + 1
              sfacrl(kcount) = cstr4
              sfacim(kcount) = sstr4
            endif
          enddo
        enddo
      enddo

      return
      end
