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
c Park and Miller RNG
c  compute in double precision, return single or double
c  no initialization needed - do NOT call with iseed = 0
      
      real*8 function ranpark(iseed)
      double precision aa,mm,sseed
      parameter (aa=16807.0D0,mm=2147483647.0D0)
      
      sseed = iseed
      sseed = mod(aa*sseed,mm)
      iseed = sseed
      ranpark = sseed/mm

      return
      end


c Marsaglia RNG
c  compute in double precision, return single or double
c  MUST be initialized by calling with iseed > 0
c  call to return a random # is by calling with iseed = 0

      real*8 function ranmars(iseed)
      double precision u(97),s,t,c,cd,cm,uni
      data iflag /0/
      save
      
      if (iseed.gt.0.or.iflag.eq.0) then
        if (iseed.eq.0) call error('Uninitialized Marsaglia RNG')
        ij = (iseed-1)/30082
        kl = (iseed-1) - 30082*ij
        i = mod(ij/177,177) + 2
        j = mod(ij,177) + 2
        k = mod(kl/169,178) + 1
        l = mod(kl,169)
        do ii = 1,97
          s = 0.0D0
          t = 0.5D0
          do jj = 1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if (mod(l*m,64).ge.32) s = s + t
            t = 0.5*t
          enddo
          u(ii) = s
        enddo
        c = 362436.0D0/16777216.0D0
        cd = 7654321.0D0/16777216.0D0
        cm = 16777213.0D0/16777216.0D0
        i97 = 97
        j97 = 33
        iflag = 1
      endif
      
      uni = u(i97) - u(j97)
      if (uni.lt.0.0) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      if (i97.eq.0) i97 = 97
      j97 = j97 - 1
      if (j97.eq.0) j97 = 97
      c = c - cd
      if (c.lt.0.0) c = c + cm
      uni = uni - c
      if (uni.lt.0.0) uni = uni + 1.0
      ranmars = uni
      
      return
      end


c Gaussian RNG using Park/Miller RNG
c  compute in double precision, return single or double

      real*8 function gausspark(iseed)
      double precision fac,gsave,rsq,v1,v2
      real*8 ranpark
      data isave /0/
      save

      if (isave.eq.0) then
 10     v1 = 2.0*ranpark(iseed)-1.0
        v2 = 2.0*ranpark(iseed)-1.0
        rsq = v1**2 + v2**2
        if (rsq.ge.1.0.or.rsq.eq.0.0) goto 10
        fac = sqrt(-2.0*log(rsq)/rsq)
        gsave = v1*fac
        gausspark = v2*fac
        isave = 1
      else
        gausspark = gsave
        isave = 0
      endif

      return
      end


c Gaussian RNG using Marsaglia RNG
c  compute in double precision, return single or double
c  ranmars is called with a 0 because we assume ranmars has
c   already been initialized

      real*8 function gaussmars()
      double precision fac,gsave,rsq,v1,v2
      real*8 ranmars
      data isave /0/
      save

      if (isave.eq.0) then
 10     v1 = 2.0*ranmars(0)-1.0
        v2 = 2.0*ranmars(0)-1.0
        rsq = v1**2 + v2**2
        if (rsq.ge.1.0.or.rsq.eq.0.0) goto 10
        fac = sqrt(-2.0*log(rsq)/rsq)
        gsave = v1*fac
        gaussmars = v2*fac
        isave = 1
      else
        gaussmars = gsave
        isave = 0
      endif

      return
      end
