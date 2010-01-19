c     This is the fortran code for making a smooth Rhs.

c     These are defined for our convenience.

c     This makes a cell-centered rhs, nice and round, ideal for testing
c     because we can make a corresponding exact solution.
      subroutine prhsmaker(rhs, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3,
     &     plow1, plow2, plow3, phigh1, phigh2, phigh3,
     &     r0, center, h, m)

      implicit none
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      integer plow1, plow2, plow3, phigh1, phigh2, phigh3
      real*8 r0, h, center(3), m, pi, p
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      integer i, j, k
      real*8 rtemp

      pi = 4.d0*atan(1.d0)
      p = 2.d0*pi*m/r0

      do k = plow3, phigh3
         do j = plow2, phigh2
            do i = plow1, phigh1
c     The following is for a smooth Rhs.
               rtemp = sqrt( (i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2 )

               if (rtemp.gt.r0) then
c     Not sure what the h**7 factor was ever for, exactly.
c                  rhs(i,j) = rhs(i, j) + h**7
c This was a no-op anyway
c                  rhs(i,j, k) = rhs(i, j, k)
               else
                  if (m.eq.(0.d0)) then
                     rhs(i, j, k) = rhs(i, j, k)
     &                    + ( (1.d0 - rtemp/r0)*(rtemp/r0) )**4
                  else
                     rhs(i, j, k) = rhs(i, j, k) + ( sin(p*rtemp)*
     &                    (1.d0 - rtemp/r0)*(rtemp/r0) )**2
                  endif
               endif

            enddo
         enddo
      enddo

      return
      end

c     make a rectangular (solid) patch for a right-hand side.
      subroutine rectrhs(rhs, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, lo, hi, charge, h)

      implicit none
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 lo(3), hi(3), charge, h
      integer rlow1, rlow2, rlow3, rhigh1, rhigh2, rhigh3
      integer i, j, k

      if (lo(1).lt.(0.d0)) then
         rlow1 = int(lo(1)/h)
      else
         rlow1 = int(lo(1)/h) + 1
      endif
      if (lo(2).lt.(0.d0)) then
         rlow2 = int(lo(2)/h)
      else
         rlow2 = int(lo(2)/h) + 1
      endif
      if (lo(3).lt.(0.d0)) then
         rlow3 = int(lo(3)/h)
      else
         rlow3 = int(lo(3)/h) + 1
      endif
      if (hi(1).lt.(0.d0)) then
         rhigh1 = int(hi(1)/h) - 1
      else
         rhigh1 = int(hi(1)/h)
      endif
      if (hi(2).lt.(0.d0)) then
         rhigh2 = int(hi(2)/h) - 1
      else
         rhigh2 = int(hi(2)/h)
      endif
      if (hi(3).lt.(0.d0)) then
         rhigh3 = int(hi(3)/h) - 1
      else
         rhigh3 = int(hi(3)/h)
      endif

      rlow1 = max(nlow1, rlow1)
      rlow2 = max(nlow2, rlow2)
      rlow3 = max(nlow3, rlow3)
      rhigh1 = min(nhigh1, rhigh1)
      rhigh2 = min(nhigh2, rhigh2)
      rhigh3 = min(nhigh3, rhigh3)

      do k = rlow3, rhigh3
         do j = rlow2, rhigh2
            do i = rlow1, rhigh1
               rhs(i, j, k) = rhs(i, j, k) + charge
            enddo
         enddo
      enddo

      return
      end

      
c     This corrects the right hand side for node-centering.
      subroutine fixrhs(rhs, nlow1, nlow2, nhigh1, nhigh2)

      integer nlow1, nlow2, nhigh1, nhigh2
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2)
      integer i, j

      do j = nlow2, nhigh2
         rhs(nlow1, j) = rhs(nlow1, j)/(2.d0)
         rhs(nhigh1, j) = rhs(nhigh1, j)/(2.d0)
      enddo

      do i = nlow1, nhigh1
         rhs(i, nlow2) = rhs(i, nlow2)/(2.d0)
         rhs(i, nhigh2) = rhs(i, nhigh2)/(2.d0)
      enddo

      return
      end


c     This makes a node-centered rhs.
      subroutine nrhsmaker (rhs, nlow1, nlow2, nhigh1, nhigh2,
     &     plow1, plow2, phigh1, phigh2, r0, center, h)

      implicit none
      integer nlow1, nlow2, nhigh1, nhigh2
      integer plow1, plow2, phigh1, phigh2
      real*8 r0, h, center(2)
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2)
      integer i, j
      real*8 rtemp

      do j = plow2, phigh2
         do i = plow1, phigh1            
            rtemp = sqrt( (i*h - center(1))
     &           *(i*h - center(1))
     &           + (j*h - center(2))
     &           *(j*h - center(2)) )

            if (rtemp.gt.r0) then
               rhs(i,j) = 0.d0
            else
               rhs(i,j) = ( (1.d0 - rtemp/r0)*(rtemp/r0) )**4
            endif

         enddo
      enddo

      return
      end



