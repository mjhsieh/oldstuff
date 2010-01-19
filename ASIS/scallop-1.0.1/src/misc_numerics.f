c     for debugging
      subroutine avestencil0(s, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, dir)
      implicit none
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3, dir
      real*8 s(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)

      return
      end

c     Average edges (and corners) of a stencil
      subroutine avestencil(s, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, dir)

      implicit none
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3, dir
      real*8 s(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      integer i, j, k

c     These definitions should match those in scallop_macros.h.  I
c     should probably just use a preprocessor to include those
c     definitions here, but I don't have enough faith in preprocessors.
      integer LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP
      parameter (LEFT = 0, RIGHT = 1, FRONT = 2, BACK = 3,
     &     BOTTOM = 4, TOP = 5)

c     The stencil passed in should actually be a plane of data, so one
c     of the dimensions will be only one unit thick.  We need to go
c     around the boundary in the other two dimensions, dividing by 2
c     (and by 4 in the corners, which happens automatically if we divide
c     by 2 twice).

c     not LEFT or RIGHT: cycle through x values
      if ((dir.ne.LEFT).and.(dir.ne.RIGHT)) then
         do i = nlow1, nhigh1
            s(i, nlow2, nlow3) = s(i, nlow2, nlow3)/(2.d0)
            s(i, nhigh2, nhigh3) = s(i, nhigh2, nhigh3)/(2.d0)
         enddo
      endif

c     not FRONT or BACK: cycle through y values
      if ((dir.ne.FRONT).and.(dir.ne.BACK)) then
         do j = nlow2, nhigh2
            s(nlow1, j, nlow3) = s(nlow1, j, nlow3)/(2.d0)
            s(nhigh1, j, nhigh3) = s(nhigh1, j, nhigh3)/(2.d0)
         enddo
      endif

c     not BOTTOM or TOP: cycle through y values
      if ((dir.ne.BOTTOM).and.(dir.ne.TOP)) then
         do k = nlow3, nhigh3
            s(nlow1, nlow2, k) = s(nlow1, nlow2, k)/(2.d0)
            s(nhigh1, nhigh2, k) = s(nhigh1, nhigh2, k)/(2.d0)
         enddo
      endif

      return
      end



      subroutine interpstencil0(f, nflow1, nflow2, nflow3,
     &     nfhigh1, nfhigh2, nfhigh3, c, nclow1, nclow2, nclow3,
     &     nchigh1, nchigh2, nchigh3, dir, nref, h)

      implicit none
      integer nflow1, nflow2, nflow3, nfhigh1, nfhigh2, nfhigh3
      integer nclow1, nclow2, nclow3, nchigh1, nchigh2, nchigh3
      real*8 f(nflow1:nfhigh1,nflow2:nfhigh2,nflow3:nfhigh3)
      real*8 c(nclow1:nchigh1,nclow2:nchigh2,nclow3:nchigh3)
      integer dir, nref
      real*8 h
      integer i, j, k
      
      do k = nflow3, nfhigh3
         do j = nflow2, nfhigh2
            do i = nflow1, nfhigh1
               f(i, j, k) = 0.d0
            enddo
         enddo
      enddo

      return
      end


      subroutine interpstencil(f, nflow1, nflow2, nflow3,
     &     nfhigh1, nfhigh2, nfhigh3, c, nclow1, nclow2, nclow3,
     &     nchigh1, nchigh2, nchigh3, dir, nref, h)

      implicit none
      integer nflow1, nflow2, nflow3, nfhigh1, nfhigh2, nfhigh3
      integer nclow1, nclow2, nclow3, nchigh1, nchigh2, nchigh3
      real*8 f(nflow1:nfhigh1,nflow2:nfhigh2,nflow3:nfhigh3)
      real*8 c(nclow1:nchigh1,nclow2:nchigh2,nclow3:nchigh3)
      integer dir, nref
      real*8 h
      integer i, j, k, mini, maxi, minj, maxj, mink, maxk
      real*8 phix, phiy, phixx, phiyy, phixxx, phiyyy, phixxxx, phiyyyy
      real*8 x, y
      real*8 result

c     These definitions should match those in scallop_macros.h.
      integer LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP
      parameter (LEFT = 0, RIGHT = 1, FRONT = 2, BACK = 3,
     &     BOTTOM = 4, TOP = 5)

      i = nclow1 + 2
      j = nclow2 + 2
      k = nclow3 + 2
c     For left and right faces
      if ((dir.eq.LEFT).or.(dir.eq.RIGHT)) then
         phixxxx = c(i, j-2, k) - 4.d0*c(i, j-1, k) + 6.d0*c(i, j, k)
     &        - 4.d0*c(i, j+1, k) + c(i, j+2, k)
         phiyyyy = c(i, j, k-2) - 4.d0*c(i, j, k-1) + 6.d0*c(i, j, k)
     &        - 4.d0*c(i, j, k+1) + c(i, j, k+2)
         phixxx = - 0.5d0*c(i, j-2, k) + c(i, j-1, k) - c(i, j+1, k)
     &        + 0.5d0*c(i, j+2, k)
         phiyyy = -0.5d0*c(i, j, k-2) + c(i, j, k-1) - c(i, j, k+1)
     &        + 0.5d0*c(i, j, k+2)
         phixx = - (1.d0/12.d0)*c(i, j-2, k) + (4.d0/3.d0)*c(i, j-1, k)
     &        - 2.5d0*c(i, j, k) + (4.d0/3.d0)*c(i, j+1, k)
     &        - (1.d0/12.d0)*c(i, j+2, k)
         phiyy = - (1.d0/12.d0)*c(i, j, k-2) + (4.d0/3.d0)*c(i, j, k-1)
     &        - 2.5d0*c(i, j, k) + (4.d0/3.d0)*c(i, j, k+1)
     &        - (1.d0/12.d0)*c(i, j, k+2)
         phix = (1.d0/12.d0)*c(i, j-2, k) - (2.d0/3.d0)*c(i, j-1, k)
     &        + (2.d0/3.d0)*c(i, j+1, k) - (1.d0/12.d0)*c(i, j+2, k)
         phiy = (1.d0/12.d0)*c(i, j, k-2) - (2.d0/3.d0)*c(i, j, k-1)
     &        + (2.d0/3.d0)*c(i, j, k+1) - (1.d0/12.d0)*c(i, j, k+2)
         minj = max((nclow2+2)*nref - nref/2, nflow2)
         maxj = min((nclow2+2)*nref + nref/2, nfhigh2)
         mink = max((nclow3+2)*nref - nref/2, nflow3)
         maxk = min((nclow3+2)*nref + nref/2, nfhigh3)
         if (dir.eq.LEFT) then
            i = nflow1
         else
            i = nfhigh1
         endif
         do j = minj, maxj
            do k = mink, maxk
               x = h*((nclow2+2)*nref - j)
               y = h*((nclow3+2)*nref - k)
               f(i, j, k) = f(i, j, k) + x*phix + 0.5d0*(x**2)*phixx
     &              + (x**3)*phixxx/(6.d0) + (x**4)*phixxxx/(24.d0)
     &              + y*phiy + 0.5d0*(y**2)*phiyy + (y**3)*phiyyy/(6.d0)
     &              + (y**4)*phiyyyy/(24.d0)
            enddo
         enddo
c     For front and back faces
      elseif ((dir.eq.FRONT).or.(dir.eq.BACK)) then
         phixxxx = c(i-2, j, k) - 4.d0*c(i-1, j, k) + 6.d0*c(i, j, k)
     &        - 4.d0*c(i+1, j, k) + c(i+2, j, k)
         phiyyyy = c(i, j, k-2) - 4.d0*c(i, j, k-1) + 6.d0*c(i, j, k)
     &        - 4.d0*c(i, j, k+1) + c(i, j, k+2)
         phixxx = - 0.5d0*c(i-2, j, k) + c(i-1, j, k) - c(i+1, j, k)
     &        + 0.5d0*c(i+2, j, k)
         phiyyy = -0.5d0*c(i, j, k-2) + c(i, j, k-1) - c(i, j, k+1)
     &        + 0.5d0*c(i, j, k+2)
         phixx = - (1.d0/12.d0)*c(i-2, j, k) + (4.d0/3.d0)*c(i-1, j, k)
     &        - 2.5d0*c(i, j, k) + (4.d0/3.d0)*c(i+1, j, k)
     &        - (1.d0/12.d0)*c(i+2, j, k)
         phiyy = - (1.d0/12.d0)*c(i, j, k-2) + (4.d0/3.d0)*c(i, j, k-1)
     &        - 2.5d0*c(i, j, k) + (4.d0/3.d0)*c(i, j, k+1)
     &        - (1.d0/12.d0)*c(i, j, k+2)
         phix = (1.d0/12.d0)*c(i-2, j, k) - (2.d0/3.d0)*c(i-1, j, k)
     &        + (2.d0/3.d0)*c(i+1, j, k) - (1.d0/12.d0)*c(i+2, j, k)
         phiy = (1.d0/12.d0)*c(i, j, k-2) - (2.d0/3.d0)*c(i, j, k-1)
     &        + (2.d0/3.d0)*c(i, j, k+1) - (1.d0/12.d0)*c(i, j, k+2)
         mini = max((nclow1+2)*nref - nref/2, nflow1)
         maxi = min((nclow1+2)*nref + nref/2, nfhigh1)
         mink = max((nclow3+2)*nref - nref/2, nflow3)
         maxk = min((nclow3+2)*nref + nref/2, nfhigh3)
         if (dir.eq.FRONT) then
            j = nflow2
         else
            j = nfhigh2
         endif
         do i = mini, maxi
            do k = mink, maxk
               x = h*((nclow1+2)*nref - i)
               y = h*((nclow3+2)*nref - k)
               f(i, j, k) = f(i, j, k) + x*phix + 0.5d0*(x**2)*phixx
     &              + (x**3)*phixxx/(6.d0) + (x**4)*phixxxx/(24.d0)
     &              + y*phiy + 0.5d0*(y**2)*phiyy + (y**3)*phiyyy/(6.d0)
     &              + (y**4)*phiyyyy/(24.d0)
            enddo
         enddo
c     For bottom and top faces
      else
         phixxxx = c(i-2, j, k) - 4.d0*c(i-1, j, k) + 6.d0*c(i, j, k)
     &        - 4.d0*c(i+1, j, k) + c(i+2, j, k)
         phiyyyy = c(i, j-2, k) - 4.d0*c(i, j-1, k) + 6.d0*c(i, j, k)
     &        - 4.d0*c(i, j, k+1) + c(i, j, k+2)
         phixxx = - 0.5d0*c(i-2, j, k) + c(i-1, j, k) - c(i+1, j, k)
     &        + 0.5d0*c(i+2, j, k)
         phiyyy = -0.5d0*c(i, j-2, k) + c(i, j-1, k) - c(i, j+1, k)
     &        + 0.5d0*c(i, j+2, k)
         phixx = - (1.d0/12.d0)*c(i-2, j, k) + (4.d0/3.d0)*c(i-1, j, k)
     &        - 2.5d0*c(i, j, k) + (4.d0/3.d0)*c(i+1, j, k)
     &        - (1.d0/12.d0)*c(i+2, j, k)
         phiyy = - (1.d0/12.d0)*c(i, j-2, k) + (4.d0/3.d0)*c(i, j-1, k)
     &        - 2.5d0*c(i, j, k) + (4.d0/3.d0)*c(i, j+1, k)
     &        - (1.d0/12.d0)*c(i, j+2, k)
         phix = (1.d0/12.d0)*c(i-2, j, k) - (2.d0/3.d0)*c(i-1, j, k)
     &        + (2.d0/3.d0)*c(i+1, j, k) - (1.d0/12.d0)*c(i+2, j, k)
         phiy = (1.d0/12.d0)*c(i, j-2, k) - (2.d0/3.d0)*c(i, j-1, k)
     &        + (2.d0/3.d0)*c(i, j+1, k) - (1.d0/12.d0)*c(i, j+2, k)
         mini = max((nclow1+2)*nref - nref/2, nflow1)
         maxi = min((nclow1+2)*nref + nref/2, nfhigh1)
         minj = max((nclow2+2)*nref - nref/2, nflow2)
         maxj = min((nclow2+2)*nref + nref/2, nfhigh2)
         if (dir.eq.BOTTOM) then
            k = nflow3
         else
            k = nfhigh3
         endif
         do i = mini, maxi
            do j = minj, maxj
               x = h*((nclow1+2)*nref - i)
               y = h*((nclow3+2)*nref - k)
               f(i, j, k) = f(i, j, k) + x*phix + 0.5d0*(x**2)*phixx
     &              + (x**3)*phixxx/(6.d0) + (x**4)*phixxxx/(24.d0)
     &              + y*phiy + 0.5d0*(y**2)*phiyy + (y**3)*phiyyy/(6.d0)
     &              + (y**4)*phiyyyy/(24.d0)
            enddo
         enddo
      endif

      return
      end


      subroutine sample(f, nflow1, nflow2, nflow3, nfhigh1, nfhigh2,
     &     nfhigh3, c, nclow1, nclow2, nclow3, nchigh1, nchigh2,
     &     nchigh3, nslow1, nslow2, nslow3, nshigh1, nshigh2, nshigh3,
     &     nref)

      implicit none
      integer nref
      integer nflow1, nflow2, nflow3, nfhigh1, nfhigh2, nfhigh3
      integer nclow1, nclow2, nclow3, nchigh1, nchigh2, nchigh3
      integer nslow1, nslow2, nslow3, nshigh1, nshigh2, nshigh3
      real*8 f(nflow1:nfhigh1,nflow2:nfhigh2,nflow3:nfhigh3)
      real*8 c(nclow1:nchigh1,nclow2:nchigh2,nclow3:nchigh3)
      integer i, j, k

      do k = nslow3, nshigh3
         do j = nslow2, nshigh2
            do i = nslow1, nshigh1
               c(i, j, k) = f(nref*i, nref*j, nref*k)
            enddo
         enddo
      enddo
      
      return
      end

