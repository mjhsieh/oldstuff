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
c

c pre-compute Ewald coefficients for potential energy and fields

      subroutine ewald_coeff
      include "lammps.h"

      pi = 4.0*atan(1.0)
      ewald_g = (1.35 - 0.15*log(long_prec))/cutcoul
      gewsr = ewald_g

      kmaxold = kmax
      nkxmx = (ewald_g*xprd/pi)*sqrt(-log(long_prec))
      nkymx = (ewald_g*yprd/pi)*sqrt(-log(long_prec))
      nkzmx = (ewald_g*zprd/pi)*sqrt(-log(long_prec))
      kmax = max(nkxmx,nkymx,nkzmx)

      if (kmax.ne.kmaxold) then
        if (node.eq.0) then
          write (6,*) 'Ewald g,kmax = ',ewald_g,kmax
          write (1,*) 'Ewald g,kmax = ',ewald_g,kmax
        endif
      endif

      if (kmax.gt.nkmax)
     $     call error('Not enough k vectors for Ewald - boost nkmax')

      gewsq = ewald_g*ewald_g
      gsqmx = -4.0*gewsq*log(long_prec)

      unitkx = 2.0*pi/xprd
      unitky = 2.0*pi/yprd
      unitkz = 2.0*pi/zprd
      preu = 4.0*pi/(xprd*yprd*zprd) 

      kcount = 0
c
c (k,0,0),(0,l,0),(0,0,m)
c
      do m = 1,kmax
        sqk = (m*unitkx)**2
        if (sqk.le.gsqmx) then
          kcount = kcount + 1
          kxvecs(kcount) = m
          kyvecs(kcount) = 0
          kzvecs(kcount) = 0
          ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
          eg(1,kcount) = 2.0*unitkx*m*ug(kcount)
          eg(2,kcount) = 0.0
          eg(3,kcount) = 0.0
        endif
        sqk = (m*unitky)**2
        if (sqk.le.gsqmx) then
          kcount = kcount + 1
          kxvecs(kcount) = 0
          kyvecs(kcount) = m
          kzvecs(kcount) = 0
          ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
          eg(1,kcount) = 0.0
          eg(2,kcount) = 2.0*unitky*m*ug(kcount)
          eg(3,kcount) = 0.0
        endif
        sqk = (m*unitkz)**2
        if (sqk.le.gsqmx) then
          kcount = kcount + 1
          kxvecs(kcount) = 0
          kyvecs(kcount) = 0
          kzvecs(kcount) = m
          ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
          eg(1,kcount) = 0.0
          eg(2,kcount) = 0.0
          eg(3,kcount) = 2.0*unitkz*m*ug(kcount)
        endif
      enddo
c
c 1 = (k,l,0), 2 = (k,-l,0)
c
      do k = 1,kmax
        do l = 1,kmax      
          sqk = (unitkx*k)**2 + (unitky*l)**2
          if (sqk.le.gsqmx) then
            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = l
            kzvecs(kcount) = 0
            ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
            eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
            eg(2,kcount) = 2.0*unitky*l*ug(kcount)
            eg(3,kcount) = 0.0

            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = -l
            kzvecs(kcount) = 0
            ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
            eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
            eg(2,kcount) = -2.0*unitky*l*ug(kcount)
            eg(3,kcount) = 0.0
          endif
        enddo
      enddo
c
c 1 = (0,l,m), 2 = (0,l,-m)
c
      do l = 1,kmax
        do m = 1,kmax      
          sqk = (unitky*l)**2+(unitkz*m)**2
          if (sqk.le.gsqmx) then
            kcount = kcount + 1
            kxvecs(kcount) = 0
            kyvecs(kcount) = l
            kzvecs(kcount) = m
            sqk = (unitky*l)**2+(unitkz*m)**2
            ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
            eg(1,kcount) =  0.0
            eg(2,kcount) =  2.0*unitky*l*ug(kcount)
            eg(3,kcount) =  2.0*unitkz*m*ug(kcount)

            kcount = kcount + 1
            kxvecs(kcount) = 0
            kyvecs(kcount) = l
            kzvecs(kcount) = -m
            ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
            eg(1,kcount) =  0.0
            eg(2,kcount) =  2.0*unitky*l*ug(kcount)
            eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
          endif
        enddo
      enddo
c
c 1 = (k,0,m), 2 = (k,0,-m)
c
      do k = 1,kmax
        do m = 1,kmax      
          sqk = (unitkx*k)**2 + (unitkz*m)**2
          if (sqk.le.gsqmx) then
            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = 0
            kzvecs(kcount) = m
            ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
            eg(1,kcount) =  2.0*unitkx*k*ug(kcount)
            eg(2,kcount) =  0.0
            eg(3,kcount) =  2.0*unitkz*m*ug(kcount)

            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = 0
            kzvecs(kcount) = -m
            ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
            eg(1,kcount) =  2.0*unitkx*k*ug(kcount)
            eg(2,kcount) =  0.0
            eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
          endif
        enddo
      enddo
c
c 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
c
      do k = 1,kmax
        do l = 1,kmax
          do m = 1,kmax      
            sqk = (unitkx*k)**2 + (unitky*l)**2 + (unitkz*m)**2
            if (sqk.le.gsqmx) then
              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = l
              kzvecs(kcount) = m
              ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = 2.0*unitky*l*ug(kcount)
              eg(3,kcount) = 2.0*unitkz*m*ug(kcount)

              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = -l
              kzvecs(kcount) = m
              ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = -2.0*unitky*l*ug(kcount)
              eg(3,kcount) = 2.0*unitkz*m*ug(kcount)

              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = l
              kzvecs(kcount) = -m
              ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = 2.0*unitky*l*ug(kcount)
              eg(3,kcount) = -2.0*unitkz*m*ug(kcount)

              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = -l
              kzvecs(kcount) = -m
              ug(kcount) = preu*exp(-0.25*sqk/gewsq)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = -2.0*unitky*l*ug(kcount)
              eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
            endif
          enddo
        enddo
      enddo

      return
      end
