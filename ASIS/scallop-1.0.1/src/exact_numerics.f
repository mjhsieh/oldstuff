
      real*8 function cosint(r, r0, m)
      implicit none
      real*8 r
      real*8 r0
      real*8 m
      real*8 answer, mstar, pi, next, nextdiv

      pi = 4.d0*atan(1.d0)
      mstar = 2.d0*pi*m/r0

      nextdiv = 2.d0
      next = - ((2.d0*mstar*r)**2)/4.d0
      answer = next
      do while (abs(next).gt.abs((1.e-14)*answer)) 
         next = - next * nextdiv * (2.d0*mstar*r)**2
         nextdiv = nextdiv + 2.d0
         next = next/(nextdiv*nextdiv*(nextdiv - 1.d0))
         answer = answer + next
      enddo
     
      cosint = answer
      
      return
      end

      real*8 function sinint(r, r0, m)
      implicit none
      real*8 r
      real*8 r0
      real*8 m
      real*8 answer, mstar, pi, next, nextdiv

      pi = 4.d0*atan(1.d0)
      mstar = 2.d0*pi*m/r0

      nextdiv = 1.d0
      next = 2.d0*mstar*r
      answer = next
      do while (abs(next).gt.abs((1.e-14)*answer)) 
         next = - next * nextdiv * (2.d0*mstar*r)**2
         nextdiv = nextdiv + 2.d0
         next = next/(nextdiv*nextdiv*(nextdiv - 1.d0))
         answer = answer + next
      enddo
     
      sinint = answer
      
      return
      end


c     This calculates the exact solution, phi, for a rhs as found in
c     prhsmaker.

      subroutine exact(phi, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, r0, center, h, m)

      implicit none
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 r0, h, center(3), pi, mstar, mstar2, m
      real*8 phi(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      integer i, j, k
      real*8 phi_inc, a, r
      real*8 cosint, sinint

      pi = 4.d0*atan(1.d0)

      if (m.eq.(0.d0)) then
         do k = nlow3, nhigh3
            do j = nlow2, nhigh2
               do i = nlow1, nhigh1
                  r = sqrt( (i*h - center(1))**2
     &                 + (j*h - center(2))**2
     &                 + (k*h - center(3))**2 )
                  if (r.lt.r0) then
                     a = r/r0
                     phi_inc = (((((1.d0/110.d0)*a
     &                    - 2.d0/45.d0)*a
     &                    + 1.d0/12.d0)*a
     &                    - 1.d0/14.d0)*a
     &                    + 1.d0/42.d0)*(a**4)*(r**2)
                     phi(i, j, k) = phi(i, j, k) + phi_inc
     &                    - (r0**2)/(1260.d0)
                     
                  else
                     phi(i, j, k) = phi(i, j, k) - (r0**3)/(2310.d0*r)
                  endif
               enddo
            enddo
         enddo
      endif

      return
      end

c     This makes a test soln in order to test the accuracy of the
c     interpolation subroutine

      subroutine makecsln(phi, nlow1, nlow2, nhigh1, nhigh2, h)

      implicit none
      integer nlow1, nlow2, nhigh1, nhigh2
      real*8 phi(nlow1:nhigh1, nlow2:nhigh2)
      real*8 center(2), h, r
      integer i, j

      center(1) = 0.5d0
      center(2) = 0.5d0

      do i = nlow1, nhigh1
         do j = nlow2, nhigh2
            r = sqrt( (i*h - center(1))**2
     &           + (j*h - center(2))**2 )
            phi(i, j) = -log(r)
         enddo
      enddo

      return
      end

      subroutine makefsln(phi, nlow1, nlow2, nhigh1, nhigh2, h)

      implicit none
      integer nlow1, nlow2, nhigh1, nhigh2
      real*8 phi(nlow1:nhigh1, nlow2:nhigh2)
      real*8 center(2), h, r
      integer i, j

      center(1) = 0.5d0
      center(2) = 0.5d0

      do i = nlow1, nhigh1
         do j = nlow2, nhigh2
            r = sqrt( (i*h - center(1))**2
     &           + (j*h - center(2))**2 )
            phi(i, j) = -log(r)
         enddo
      enddo

      return
      end

