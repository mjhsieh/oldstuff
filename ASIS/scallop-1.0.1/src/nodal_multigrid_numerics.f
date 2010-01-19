      subroutine lap7(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3,
     &     l, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3,h)
      integer nlowp1, nlowp2, nlowp3, nhighp1, nhighp2, nhighp3
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 p(nlowp1:nhighp1, nlowp2:nhighp2, nlowp3:nhighp3)
      real*8 l(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k

      do k = nlowp3 + 1, nhighp3 - 1
         do j = nlowp2 + 1, nhighp2 - 1
            do i = nlowp1 + 1, nhighp1 - 1
               l(i,j,k) = (p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1)
     &              - 6.d0*p(i,j,k))/(h*h)
            enddo
         enddo
      enddo

      return
      end


      subroutine gsrb7(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3,
     &     rhs, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, h)
      
      integer nlowp1, nlowp2, nlowp3, nhighp1, nhighp2, nhighp3
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 p(nlowp1:nhighp1, nlowp2:nhighp2, nlowp3:nhighp3)
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k, ioff, iinc

c     red
      call n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)
      ioff = 1
      do k = nlowp3 + 1, nhighp3 - 1
         do j = nlowp2 + 1, nhighp2 - 1
            iinc = abs(mod(j + k + ioff, 2))
            do i = nlowp1 + iinc + 1, nhighp1 - 1, 2
               p(i,j,k) = (p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1)
     &              - rhs(i,j,k)*h*h)/6.d0
            enddo
         enddo
      enddo

c     black
      call n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)
      ioff = 0
      do k = nlowp3 + 1, nhighp3 - 1
         do j = nlowp2 + 1, nhighp2 - 1
            iinc = abs(mod(j + k + ioff, 2))
            do i = nlowp1 + iinc + 1, nhighp1 - 1, 2
               p(i,j,k) = (p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1)
     &              - rhs(i,j,k)*h*h)/6.d0
            enddo
         enddo
      enddo

      return
      end


      subroutine resid7(p, rhs, res, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, h)

      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 p(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 res(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k
      
      call n_boundary(p, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3)

      do k = nlow3 + 1, nhigh3 - 1
         do j = nlow2 + 1, nhigh2 - 1
            do i = nlow1 + 1, nhigh1 - 1
               res(i, j, k) = rhs(i, j, k)
     &              - (p(i + 1, j, k) + p(i - 1, j, k)
     &              + p(i, j + 1, k) + p(i, j - 1, k)
     &              + p(i, j, k + 1) + p(i, j, k - 1)
     &              - 6.d0*p(i, j, k))/(h*h)
            enddo
         enddo
      enddo
      
      return
      end


      subroutine lap19(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3,
     &     l, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3,h)
      implicit none
      integer nlowp1, nlowp2, nlowp3, nhighp1, nhighp2, nhighp3
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      integer low1, low2, low3, high1, high2, high3
      real*8 p(nlowp1:nhighp1, nlowp2:nhighp2, nlowp3:nhighp3)
      real*8 l(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k

      low1 = max(nlowp1+1, nlow1)
      low2 = max(nlowp2+1, nlow2)
      low3 = max(nlowp3+1, nlow3)
      high1 = min(nhighp1-1, nhigh1)
      high2 = min(nhighp2-1, nhigh2)
      high3 = min(nhighp3-1, nhigh3)
      do k = low3, high3
         do j = low2, high2
            do i = low1, high1
               l(i,j,k) = (p(i+1, j+1, k) + p(i+1, j-1, k)
     &              + p(i-1, j+1, k) + p(i-1, j-1, k)
     &              + p(i+1, j, k+1) + p(i+1, j, k-1)
     &              + p(i-1, j, k+1) + p(i-1, j, k-1)
     &              + p(i, j+1, k+1) + p(i, j+1, k-1)
     &              + p(i, j-1, k+1) + p(i, j-1, k-1)
     &              + 2.d0*(p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1))
     &              - 24.d0*p(i,j,k))/(6.d0*h*h)
            enddo
         enddo
      enddo

      return
      end


      subroutine gsfc19(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3,
     &     rhs, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, h)
      
      integer nlowp1, nlowp2, nlowp3, nhighp1, nhighp2, nhighp3
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 p(nlowp1:nhighp1, nlowp2:nhighp2, nlowp3:nhighp3)
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k, ioff, joff, iinc, jinc

c     For this larger stencil, two colors are not enough: we need four
c     (thence gsfc).

c     orange stars
      call n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)
      ioff = 1
      joff = 1
      do k = nlowp3 + 1, nhighp3 - 1
         jinc = abs(mod(k + joff, 2))
         do j = nlowp2 + jinc + 1, nhighp2 - 1, 2
            iinc = abs(mod(k + ioff, 2))
            do i = nlowp1 + iinc + 1, nhighp1 - 1, 2
               p(i,j,k) = ( p(i+1, j+1, k) + p(i+1, j-1, k)
     &              + p(i-1, j+1, k) + p(i-1, j-1, k)
     &              + p(i+1, j, k+1) + p(i+1, j, k-1)
     &              + p(i-1, j, k+1) + p(i-1, j, k-1)
     &              + p(i, j+1, k+1) + p(i, j+1, k-1)
     &              + p(i, j-1, k+1) + p(i, j-1, k-1)
     &              + 2.d0*(p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1))
     &              - 6.d0*h*h*rhs(i,j,k) )/24.d0
            enddo
         enddo
      enddo

c     yellow moons
      call n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)
      ioff = 0
      joff = 0
      do k = nlowp3 + 1, nhighp3 - 1
         jinc = abs(mod(k + joff, 2))
         do j = nlowp2 + jinc + 1, nhighp2 - 1, 2
            iinc = abs(mod(k + ioff, 2))
            do i = nlowp1 + iinc + 1, nhighp1 - 1, 2
               p(i,j,k) = ( p(i+1, j+1, k) + p(i+1, j-1, k)
     &              + p(i-1, j+1, k) + p(i-1, j-1, k)
     &              + p(i+1, j, k+1) + p(i+1, j, k-1)
     &              + p(i-1, j, k+1) + p(i-1, j, k-1)
     &              + p(i, j+1, k+1) + p(i, j+1, k-1)
     &              + p(i, j-1, k+1) + p(i, j-1, k-1)
     &              + 2.d0*(p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1))
     &              - 6.d0*h*h*rhs(i,j,k) )/24.d0
            enddo
         enddo
      enddo

c     blue diamonds
      call n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)
      ioff = 0
      joff = 1
      do k = nlowp3 + 1, nhighp3 - 1
         jinc = abs(mod(k + joff, 2))
         do j = nlowp2 + jinc + 1, nhighp2 - 1, 2
            iinc = abs(mod(k + ioff, 2))
            do i = nlowp1 + iinc + 1, nhighp1 - 1, 2
               p(i,j,k) = ( p(i+1, j+1, k) + p(i+1, j-1, k)
     &              + p(i-1, j+1, k) + p(i-1, j-1, k)
     &              + p(i+1, j, k+1) + p(i+1, j, k-1)
     &              + p(i-1, j, k+1) + p(i-1, j, k-1)
     &              + p(i, j+1, k+1) + p(i, j+1, k-1)
     &              + p(i, j-1, k+1) + p(i, j-1, k-1)
     &              + 2.d0*(p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1))
     &              - 6.d0*h*h*rhs(i,j,k) )/24.d0
            enddo
         enddo
      enddo

c     green clovers
      call n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)
      ioff = 1
      joff = 0
      do k = nlowp3 + 1, nhighp3 - 1
         jinc = abs(mod(k + joff, 2))
         do j = nlowp2 + jinc + 1, nhighp2 - 1, 2
            iinc = abs(mod(k + ioff, 2))
            do i = nlowp1 + iinc + 1, nhighp1 - 1, 2
               p(i,j,k) = ( p(i+1, j+1, k) + p(i+1, j-1, k)
     &              + p(i-1, j+1, k) + p(i-1, j-1, k)
     &              + p(i+1, j, k+1) + p(i+1, j, k-1)
     &              + p(i-1, j, k+1) + p(i-1, j, k-1)
     &              + p(i, j+1, k+1) + p(i, j+1, k-1)
     &              + p(i, j-1, k+1) + p(i, j-1, k-1)
     &              + 2.d0*(p(i+1, j, k) + p(i-1, j, k)
     &              + p(i, j+1, k) + p(i, j-1, k)
     &              + p(i, j, k+1) + p(i, j, k-1))
     &              - 6.d0*h*h*rhs(i,j,k) )/24.d0
            enddo
         enddo
      enddo

      return
      end


      subroutine resid19(p, rhs, res, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, h)

      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 p(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 res(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k

      call n_boundary(p, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3)

      do k = nlow3 + 1, nhigh3 - 1
         do j = nlow2 + 1, nhigh2 - 1
            do i = nlow1 + 1, nhigh1 - 1
               res(i, j, k) = rhs(i, j, k)
     &              - (p(i + 1, j + 1, k) + p(i + 1, j - 1, k)
     &              + p(i - 1, j + 1, k) + p(i - 1, j - 1, k)
     &              + p(i + 1, j, k + 1) + p(i + 1, j, k - 1)
     &              + p(i - 1, j, k + 1) + p(i - 1, j, k - 1)
     &              + p(i, j + 1, k + 1) + p(i, j + 1, k - 1)
     &              + p(i, j - 1, k + 1) + p(i, j - 1, k - 1)
     &              + 2.d0*(p(i + 1, j, k) + p(i - 1, j, k)
     &              + p(i, j + 1, k) + p(i, j - 1, k)
     &              + p(i, j, k + 1) + p(i, j, k - 1))
     &              - 24.d0*p(i, j, k))/(6.d0*h*h)
            enddo
         enddo
      enddo
      
      return
      end


      subroutine naverage(fine, nflow1, nflow2, nflow3,
     &     nfhigh1, nfhigh2, nfhigh3,
     &     coarse, nclow1, nclow2, nclow3,
     &     nchigh1, nchigh2, nchigh3)
      

      integer nflow1, nflow2, nflow3, nfhigh1, nfhigh2, nfhigh3
      integer nclow1, nclow2, nclow3, nchigh1, nchigh2, nchigh3
      real*8 fine(nflow1:nfhigh1, nflow2:nfhigh2, nflow3:nfhigh3)
      real*8 coarse(nclow1:nchigh1, nclow2:nchigh2, nclow3:nchigh3)
      
      integer i, j, k, ifine, jfine, kfine

      do k = nclow3 + 1, nchigh3 - 1
         do j = nclow2 + 1, nchigh2 - 1
            do i = nclow1 + 1, nchigh1 - 1
               
               ifine = 2*i
               jfine = 2*j
               kfine = 2*k

               coarse(i, j, k) = fine(ifine, jfine, kfine)

            enddo
         enddo
      enddo

      return
      end


c     I don't think this version of averaging is necessary for
c     nodal-point schemes, but I've made it for testing purposes, 
c     just in case.
      subroutine nhaverage(fine, nflow1, nflow2, nflow3,
     &     nfhigh1, nfhigh2, nfhigh3,
     &     coarse, nclow1, nclow2, nclow3,
     &     nchigh1, nchigh2, nchigh3)
      

      integer nflow1, nflow2, nflow3, nfhigh1, nfhigh2, nfhigh3
      integer nclow1, nclow2, nclow3, nchigh1, nchigh2, nchigh3
      real*8 fine(nflow1:nfhigh1, nflow2:nfhigh2, nflow3:nfhigh3)
      real*8 coarse(nclow1:nchigh1, nclow2:nchigh2, nclow3:nchigh3)
      
      integer i, j, k, ifine, jfine, kfine

      do k = nclow3 + 1, nchigh3 - 1
         do j = nclow2 + 1, nchigh2 - 1
            do i = nclow1 + 1, nchigh1 - 1
               
               ifine = 2*i
               jfine = 2*j
               kfine = 2*k

c     This is ugly, but we have 27 fine points (nodes) to average into
c     a single coarse value.  Imagine the points laid out on a cube.
c     The formula below is just
c     coarse = (8*(1 center value) + 4*(6 face-center values)
c               + 2*(12 edge-center values) + (8 corner values))/64

               coarse(i,j,k) = 0.015625d0*
     &              (8.d0*fine(ifine, jfine, kfine)
     &              + 4.d0*(fine(ifine + 1, jfine, kfine)
     &              + fine(ifine - 1, jfine, kfine)
     &              + fine(ifine, jfine + 1, kfine)
     &              + fine(ifine, jfine - 1, kfine)
     &              + fine(ifine, jfine, kfine + 1)
     &              + fine(ifine, jfine, kfine - 1))
     &              + 2.d0*(fine(ifine + 1, jfine + 1, kfine)
     &              + fine(ifine + 1, jfine - 1, kfine)
     &              + fine(ifine - 1, jfine + 1, kfine)
     &              + fine(ifine - 1, jfine - 1, kfine)
     &              + fine(ifine + 1, jfine, kfine + 1)
     &              + fine(ifine + 1, jfine, kfine - 1)
     &              + fine(ifine - 1, jfine, kfine + 1)
     &              + fine(ifine - 1, jfine, kfine - 1)
     &              + fine(ifine, jfine + 1, kfine + 1)
     &              + fine(ifine, jfine + 1, kfine - 1)
     &              + fine(ifine, jfine - 1, kfine + 1)
     &              + fine(ifine, jfine - 1, kfine - 1))
     &              + fine(ifine + 1, jfine + 1, kfine + 1)
     &              + fine(ifine + 1, jfine + 1, kfine - 1)
     &              + fine(ifine + 1, jfine - 1, kfine + 1)
     &              + fine(ifine + 1, jfine - 1, kfine - 1)
     &              + fine(ifine - 1, jfine + 1, kfine + 1)
     &              + fine(ifine - 1, jfine + 1, kfine - 1)
     &              + fine(ifine - 1, jfine - 1, kfine + 1)
     &              + fine(ifine - 1, jfine - 1, kfine - 1))
               
            enddo
         enddo
      enddo
      
      return
      end
      

      subroutine ninterp(fine, nflow1, nflow2, nflow3,
     &     nfhigh1, nfhigh2, nfhigh3,
     &     coarse, nclow1, nclow2, nclow3,
     &     nchigh1, nchigh2, nchigh3)
      
      integer nflow1, nflow2, nflow3, nfhigh1, nfhigh2, nfhigh3
      integer nclow1, nclow2, nclow3, nchigh1, nchigh2, nchigh3
      real*8 fine(nflow1:nfhigh1, nflow2:nfhigh2, nflow3:nfhigh3)
      real*8 coarse(nclow1:nchigh1, nclow2:nchigh2, nclow3:nchigh3)
      integer i, j, k

      do k = nclow3 + 1, nchigh3
         do j = nclow2 + 1, nchigh2
            do i = nclow1 + 1, nchigh1
               fine(2*i - 1, 2*j - 1, 2*k - 1) =
     &              0.125d0*(coarse(i, j, k) + coarse(i, j, k - 1)
     &              + coarse(i, j - 1, k) + coarse(i, j - 1, k - 1)
     &              + coarse(i - 1, j, k) + coarse(i - 1, j, k - 1)
     &              + coarse(i - 1, j - 1, k)
     &              + coarse(i - 1, j - 1, k - 1))
     &              + fine(2*i - 1, 2*j - 1, 2*k - 1)
               fine(2*i - 1, 2*j - 1, 2*k) =
     &              0.25d0*(coarse(i, j, k) + coarse(i - 1, j, k)
     &              + coarse(i, j - 1, k) + coarse(i - 1, j - 1, k))
     &              + fine(2*i - 1, 2*j - 1, 2*k)
               fine(2*i - 1, 2*j, 2*k - 1) =
     &              0.25d0*(coarse(i, j, k) + coarse(i - 1, j, k)
     &              + coarse(i, j, k - 1) + coarse(i - 1, j, k - 1))
     &              + fine(2*i - 1, 2*j, 2*k - 1)
               fine(2*i, 2*j - 1, 2*k - 1) =
     &              0.25d0*(coarse(i, j, k) + coarse(i, j - 1, k)
     &              + coarse(i, j, k - 1) + coarse(i, j - 1, k - 1))
     &              + fine(2*i, 2*j - 1, 2*k - 1)
               fine(2*i - 1, 2*j, 2*k) =
     &              0.5d0*(coarse(i, j, k) + coarse(i - 1, j, k))
     &              + fine(2*i - 1, 2*j, 2*k)
               fine(2*i, 2*j - 1, 2*k) =
     &              0.5d0*(coarse(i, j, k) + coarse(i, j - 1, k))
     &              + fine(2*i, 2*j - 1, 2*k)
               fine(2*i, 2*j, 2*k - 1) =
     &              0.5d0*(coarse(i, j, k) + coarse(i, j, k - 1))
     &              + fine(2*i, 2*j, 2*k - 1)
               fine(2*i, 2*j, 2*k) =
     &              coarse(i, j, k) + fine(2*i, 2*j, 2*k)
            enddo
         enddo
      enddo
      
      call n_boundary(fine, nflow1, nflow2, nflow3,
     &     nfhigh1, nfhigh2, nfhigh3)
      
      return
      end


      subroutine n_boundary(p, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3)

      integer nlowp1, nlowp2, nlowp3, nhighp1, nhighp2, nhighp3
      real*8 p(nlowp1:nhighp1, nlowp2:nhighp2, nlowp3:nhighp3)
      integer i, j, k

      do k = nlowp3, nhighp3
         do j = nlowp2, nhighp2
            p(nlowp1, j, k) = 0.d0
            p(nhighp1, j, k) = 0.d0
         enddo
      enddo

      do j = nlowp2, nhighp2
         do i = nlowp1, nhighp1
            p(i, j, nlowp3) = 0
            p(i, j, nhighp3) = 0
         enddo
      enddo

      do k = nlowp3, nhighp3
         do i = nlowp1, nhighp1
            p(i, nlowp2, k) = 0
            p(i, nhighp2, k) = 0
         enddo
      enddo

      return
      end


      subroutine norm(p, nlow1, nlow2, nlow3,
     &     nhigh1, nhigh2, nhigh3, result)
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 p(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 result
      integer i, j, k

      result = 0.
      do k = nlow3, nhigh3
         do j = nlow2, nhigh2
            do i = nlow1, nhigh1
               result = result + p(i, j, k)*p(i, j, k)
            enddo
         enddo
      enddo  
      
      result = sqrt(result)

      return 
      end


      subroutine nidbound(phi, nlowp1, nlowp2, nlowp3,
     &     nhighp1, nhighp2, nhighp3,
     &     rhs, nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3, h)
      
      implicit none
      integer nlowp1, nlowp2, nlowp3, nhighp1, nhighp2, nhighp3
      integer nlow1, nlow2, nlow3, nhigh1, nhigh2, nhigh3
      real*8 phi(nlowp1:nhighp1, nlowp2:nhighp2, nlowp3:nhighp3)
      real*8 rhs(nlow1:nhigh1, nlow2:nhigh2, nlow3:nhigh3)
      real*8 h
      integer i, j, k
      real*8 center(3), pi, dist, q
      
c Find the sum of all the charges.
      
      q = 0.d0
      do k = nlow3, nhigh3
         do j = nlow2, nhigh2
            do i = nlow1, nhigh1
               q = q + rhs(i, j, k)
            enddo
         enddo
      enddo
      q = q*(h**3)
      
      if (q.ne.(0.d0)) then

c     Find the center of charge.
         
         center(1) = 0.d0
         center(2) = 0.d0
         center(3) = 0.d0
         do k = nlow3, nhigh3
            do j = nlow2, nhigh2
               do i = nlow1, nhigh1
                  center(1) = center(1) + i*h*rhs(i, j, k)
                  center(2) = center(2) + j*h*rhs(i, j, k)
                  center(3) = center(3) + k*h*rhs(i, j, k)
               enddo
            enddo
         enddo
         center(1) = (center(1)*(h**3))/q
         center(2) = (center(2)*(h**3))/q
         center(3) = (center(3)*(h**3))/q

         pi = 4.d0*atan(1.d0)
         
c     The border is set based based on the following approximation.
c     All the charge is viewed as a point charge.  The potential at any
c     point x is then p(x) = (1/(2*pi))*q*log(|x-xc|), where |x-xc| is 
c     the distance between the point x and the charge.  For convenience, 
c     I simply leave the distance squared and divide the log by another 
c     factor of 2.  This result is added onto whatever is sitting in 
c     phi, because there are other contributions that we do not want to
c     erase.

         
c     set top boundary
         k = nhighp3
         do j = nlowp2, nhighp2
            do i = nlowp1, nhighp1
               dist = sqrt((i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2)
               phi(i, j, k) = phi(i, j, k) - q/(4.d0*pi*dist)
            enddo
         enddo
c     set bottom boundary
         k = nlowp3
         do j = nlowp2, nhighp2
            do i = nlowp1, nhighp1
               dist = sqrt((i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2)
               phi(i, j, k) = phi(i, j, k) - q/(4.d0*pi*dist)
            enddo
         enddo
c     set back boundary
         j = nhighp2
         do k = nlowp3 + 1, nhighp3 - 1
            do i = nlowp1, nhighp1
               dist = sqrt((i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2)
               phi(i, j, k) = phi(i, j, k) - q/(4.d0*pi*dist)
            enddo
         enddo
c     set front boundary
         j = nlowp2
         do k = nlowp3 + 1, nhighp3 - 1
            do i = nlowp1, nhighp1
               dist = sqrt((i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2)
               phi(i, j, k) = phi(i, j, k) - q/(4.d0*pi*dist)
            enddo
         enddo
c     set right boundary
         i = nhighp1
         do k = nlowp3 + 1, nhighp3 - 1
            do j = nlowp2 + 1, nhighp2 - 1
               dist = sqrt((i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2)
               phi(i, j, k) = phi(i, j, k) - q/(4.d0*pi*dist)
            enddo
         enddo
c     set left boundary
         i = nlowp1
         do k = nlowp3 + 1, nhighp3 - 1
            do j = nlowp2 + 1, nhighp2 - 1
               dist = sqrt((i*h - center(1))**2
     &              + (j*h - center(2))**2
     &              + (k*h - center(3))**2)
               phi(i, j, k) = phi(i, j, k) - q/(4.d0*pi*dist)
            enddo
         enddo

      endif

      return
      end


      subroutine npotcalc(phih, nlhp1, nlhp2, nlhp3, nhhp1, nhhp2,
     &     nhhp3, phi, nlp1, nlp2, nlp3, nhp1, nhp2, nhp3, extra, rhs,
     &     nlr1, nlr2, nlr3, nhr1, nhr2, nhr3, 
     &     nlcp1, nlcp2, nlcp3, nhcp1, nhcp2, nhcp3,
     &     lefth, righth, bottomh, toph, backh, fronth,
     &     left, right, bottom, top, back, front, 
     &     leftc, rightc, bottomc, topc, backc, frontc,
     &     nref, h)

      implicit none
      integer nlhp1, nlhp2, nlhp3, nhhp1, nhhp2, nhhp3
      integer nlp1, nlp2, nlp3, nhp1, nhp2, nhp3
      integer nlcp1, nlcp2, nlcp3, nhcp1, nhcp2, nhcp3
      integer nlr1, nlr2, nlr3, nhr1, nhr2, nhr3
      real*8 phih(nlhp1:nhhp1, nlhp2:nhhp2, nlhp3:nhhp3)
      real*8 phi(nlp1:nhp1, nlp2:nhp2, nlp3:nhp3)
      real*8 extra(nlp1:nhp1, nlp2:nhp2, nlp3:nhp3)
      real*8 rhs(nlr1:nhr1, nlr2:nhr2, nlr3:nhr3)

      real*8 lefth(nlhp2:nhhp2, nlhp3:nhhp3)
      real*8 righth(nlhp2:nhhp2, nlhp3:nhhp3)
      real*8 fronth(nlhp1:nhhp1, nlhp3:nhhp3)
      real*8 backh(nlhp1:nhhp1, nlhp3:nhhp3)
      real*8 bottomh(nlhp1:nhhp1, nlhp2:nhhp2)
      real*8 toph(nlhp1:nhhp1, nlhp2:nhhp2)
      real*8 left(nlp2:nhp2, nlp3:nhp3)
      real*8 right(nlp2:nhp2, nlp3:nhp3)
      real*8 front(nlp1:nhp1, nlp3:nhp3)
      real*8 back(nlp1:nhp1, nlp3:nhp3)
      real*8 bottom(nlp1:nhp1, nlp2:nhp2)
      real*8 top(nlp1:nhp1, nlp2:nhp2)
      real*8 leftc(nlcp2:nhcp2, nlcp3:nhcp3)
      real*8 rightc(nlcp2:nhcp2, nlcp3:nhcp3)
      real*8 frontc(nlcp1:nhcp1, nlcp3:nhcp3)
      real*8 backc(nlcp1:nhcp1, nlcp3:nhcp3)
      real*8 bottomc(nlcp1:nhcp1, nlcp2:nhcp2)
      real*8 topc(nlcp1:nhcp1, nlcp2:nhcp2)

      integer nref
      real*8 h

      real*8 dist, pi, q, center(3)
      integer i, j, k, ih, jh, kh
      real*8 xi, xhl, xhr, xl, xr, xih
      real*8 yj, yhf, yhb, yf, yb, yjh
      real*8 zk, zhb, zht, zb, zt, zkh

c     we need temporary storage for coarsened planes from which
c     we will interpolate boundary values

      real*8 tcxy(nlp1/nref:nhp1/nref, nlp2/nref:nhp2/nref,
     &     nlp3/nref:nhp3/nref)
c      real*8 tcyz()
c      real*8 tcxz()

c Find the sum of all the charges.
      
c      write(*,*) nlhp1, nlp1, nhhp1, nhp1, nlr1, nhr1
      q = 0.d0
      do k = nlr3, nhr3
         do j = nlr2, nhr2
            do i = nlr1, nhr1
               q = q + rhs(i, j, k)
            enddo
         enddo
      enddo
      q = q*(h**3)
      
      if (q.ne.(0.d0)) then

c     Find the center of charge.

         
         center(1) = 0.d0
         center(2) = 0.d0
         center(3) = 0.d0
         do k = nlr3, nhr3
            do j = nlr2, nhr2
               do i = nlr1, nhr1
                  center(1) = center(1) + i*h*rhs(i, j, k)
                  center(2) = center(2) + j*h*rhs(i, j, k)
                  center(3) = center(3) + k*h*rhs(i, j, k)
               enddo
            enddo
         enddo

         center(1) = (center(1)*(h**3))/q
         center(2) = (center(2)*(h**3))/q
         center(3) = (center(3)*(h**3))/q

         pi = 4.d0*atan(1.d0)

c Compute locations of the outer edges of both the phih and phi grids 
c (note that edges are located at integer multiples of h).

         xhl = nlhp1*h
         xhr = nhhp1*h
         xl = nlp1*h
         xr = nhp1*h
         yhf = nlhp2*h
         yhb = nhhp2*h
         yf = nlp2*h
         yb = nhp2*h
         zhb = nlhp3*h
         zht = nhhp3*h
         zb = nlp3*h
         zt = nhp3*h

c Compute field induced by boundary charge at outer edge of large grid:

c We're going to waste some memory here, and worry about it later....

c     right and left planes
         do k = nlhp3, nhhp3
            do j = nlhp2, nhhp2
               lefth(j, k) = -25.d0*phih(nlhp1, j, k)/12.d0
     &              + 4.d0*phih(nlhp1 + 1, j, k)
     &              - 3.d0*phih(nlhp1 + 2, j, k)
     &              + 4.d0*phih(nlhp1 + 3, j, k)/3.d0
     &              - 0.25d0*phih(nlhp1 + 4, j, k)
               lefth(j, k) = h*lefth(j, k)
     &              - h*h*q*(xhl - center(1))/(4.d0*pi
     &              *( (xhl - center(1))*(xhl - center(1))
     &              + (j*h - center(2))*(j*h - center(2))
     &              + (k*h - center(3))*(k*h - center(3)) )**1.5d0)
               righth(j, k) = -25.d0*phih(nhhp1, j, k)/12.d0
     &              + 4.d0*phih(nhhp1 - 1, j, k)
     &              - 3.d0*phih(nhhp1 - 2, j, k)
     &              + 4.d0*phih(nhhp1 - 3, j, k)/3.d0
     &              - 0.25d0*phih(nhhp1 - 4, j, k)
               righth(j, k) = h*righth(j, k)
     &              + h*h*q*(xhr - center(1))/(4.d0*pi
     &              *( (xhr - center(1))*(xhr - center(1))
     &              + (j*h - center(2))*(j*h - center(2))
     &              + (k*h - center(3))*(k*h - center(3)) )**1.5d0)
            enddo
         enddo

c     front and back planes
         do k = nlhp3, nhhp3
            do i = nlhp1 + 1, nhhp1 - 1
               fronth(i, k) = -25.d0*phih(i, nlhp2, k)/12.d0
     &              + 4.d0*phih(i, nlhp2 + 1, k)
     &              - 3.d0*phih(i, nlhp2 + 2, k)
     &              + 4.d0*phih(i, nlhp2 + 3, k)/3.d0
     &              - 0.25d0*phih(i, nlhp2 + 4, k)
               fronth(i, k) = h*fronth(i, k)
     &              - h*h*q*(yhf - center(2))/(4.d0*pi
     &              *( (i*h - center(1))*(i*h - center(1))
     &              + (yhf - center(2))*(yhf - center(2))
     &              + (k*h - center(3))*(k*h - center(3)) )**1.5d0)
               backh(i, k) = -25.d0*phih(i, nhhp2, k)/12.d0
     &              + 4.d0*phih(i, nhhp2 - 1, k)
     &              - 3.d0*phih(i, nhhp2 - 2, k)
     &              + 4.d0*phih(i, nhhp2 - 3, k)/3.d0
     &              - 0.25d0*phih(i, nhhp2 - 4, k)
               backh(i, k) = h*backh(i, k)
     &              + h*h*q*(yhb - center(1))/(4.d0*pi
     &              *( (i*h - center(1))*(i*h - center(1))
     &              + (yhb - center(2))*(yhb - center(2))
     &              + (k*h - center(3))*(k*h - center(3)) )**1.5d0)
            enddo
         enddo

c     top and bottom planes
         do j = nlhp2 + 1, nhhp2 - 1
            do i = nlhp1 + 1, nhhp1 - 1
               bottomh(i, j) = -25.d0*phih(i, j, nlhp3)/12.d0
     &              + 4.d0*phih(i, j, nlhp3 + 1)
     &              - 3.d0*phih(i, j, nlhp3 + 2)
     &              + 4.d0*phih(i, j, nlhp3 + 3)/3.d0
     &              - 0.25d0*phih(i, j, nlhp3 + 4)
               bottom(i, j) = h*bottomh(i, j)
     &              - h*h*q*(zhb - center(3))/(4.d0*pi
     &              *( (i*h - center(1))*(i*h - center(1))
     &              + (j*h - center(2))*(j*h - center(2))
     &              + (zhb - center(3))*(zhb - center(3)) )**1.5d0)
               toph(i, j) = -25.d0*phih(i, j, nhhp3)/12.d0
     &              + 4.d0*phih(i, j, nhhp3 - 1)
     &              - 3.d0*phih(i, j, nhhp3 - 2)
     &              + 4.d0*phih(i, j, nhhp3 - 3)/3.d0
     &              - 0.25d0*phih(i, j, nhhp3 - 4)
               toph(i, j) = h*toph(i, j)
     &              + h*h*q*(zht - center(3))/(4.d0*pi
     &              *( (i*h - center(1))*(i*h - center(1))
     &              + (j*h - center(2))*(j*h - center(2))
     &              + (zht - center(3))*(zht - center(3)) )**1.5d0)
            enddo
         enddo

      endif

c     remove all previous data from phi
      do k = nlp3, nhp3
         do j = nlp2, nhp2
            do i = nlp1, nhp1
               phi(i, j, k) = 0.d0
            enddo
         enddo
      enddo

      if (q.ne.(0.d0)) then
c     In the 2d version we used the trapezoidal rule to provide more
c     accurate integration.  I'll add that back if necessary.

c     In the 2d version, the outer planes were the inner loops (because they
c     were the larger arrays).  In the 3d version, since we'll calculate
c     the field on sampled outer planes, the outer planes will be smaller.
c     So the inner loops are the the inner arrays in the 3d version.

c     compute field for left and right outer planes
         do k = nlcp3, nhcp3
            zk = k*h
            do j = nlcp2, nhcp2
               yj = j*h
c     due to left and right inner planes
               do kh = nlhp3, nhhp3
                  zkh = kh*h
                  do jh = nlhp2, nhhp2
                     yjh = jh*h

                     dist = sqrt( (xl - xhl)**2 + (yj - yjh)**2
     &                    + (zk - zkh)**2 )
                     leftc(j, k) = leftc(j, k)
     &                    - lefth(jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xl - xhr)**2 + (yj - yjh)**2
     &                    + (zk - zkh)**2 )
                     leftc(j, k) = leftc(j, k)
     &                    - righth(jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xr - xhl)**2 + (yj - yjh)**2
     &                    + (zk - zkh)**2 )
                     rightc(j, k) = rightc(j, k)
     &                    - lefth(jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xr - xhr)**2 + (yj - yjh)**2
     &                    + (zk - zkh)**2 )
                     rightc(j, k) = rightc(j, k)
     &                    - righth(jh, kh)/(4.d0*pi*dist)
                  enddo
               enddo
c     due to front and back inner planes
               do kh = nlhp3, nhhp3
                  zkh = kh*h
                  do ih = nlhp1 + 1, nhhp1 - 1
                     xih = ih*h

                     dist = sqrt( (xl - xih)**2 + (yj - yhf)**2
     &                    + (zk - zkh)**2 )
                     leftc(j, k) = leftc(j, k)
     &                    - fronth(ih, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xl - xih)**2 + (yj - yhb)**2
     &                    + (zk - zkh)**2 )
                     leftc(j, k) = leftc(j, k)
     &                    - backh(ih, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xr - xih)**2 + (yj - yhf)**2
     &                    + (zk - zkh)**2 )
                     rightc(j, k) = rightc(j, k)
     &                    - fronth(ih, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xr - xih)**2 + (yj - yhb)**2
     &                    + (zk - zkh)**2 )
                     rightc(j, k) = rightc(j, k)
     &                    - backh(ih, kh)/(4.d0*pi*dist)
                  enddo
               enddo
c     due to bottom and top inner planes
               do jh = nlhp2 + 1, nhhp2 - 1
                  yjh = jh*h
                  do ih = nlhp1 + 1, nhhp1 - 1
                     xih = ih*h

                     dist = sqrt( (xl - xih)**2 + (yj - yjh)**2
     &                    + (zk - zhb)**2 )
                     leftc(j, k) = leftc(j, k)
     &                    - bottomh(ih, jh)/(4.d0*pi*dist)

                     dist = sqrt( (xl - xih)**2 + (yj - yjh)**2
     &                    + (zk - zht)**2 )
                     leftc(j, k) = leftc(j, k)
     &                    - toph(ih, jh)/(4.d0*pi*dist)

                     dist = sqrt( (xr - xih)**2 + (yj - yjh)**2
     &                    + (zk - zhb)**2 )
                     rightc(j, k) = rightc(j, k)
     &                    - bottomh(ih, jh)/(4.d0*pi*dist)

                     dist = sqrt( (xr - xih)**2 + (yj - yjh)**2
     &                    + (zk - zht)**2 )
                     rightc(j, k) = rightc(j, k)
     &                    - toph(ih, jh)/(4.d0*pi*dist)
                  enddo
               enddo

            enddo
         enddo

c     compute field for front and back outer planes
         do k = nlcp3, nhcp3
            zk = k*h
            do i = nlcp1, nhcp1
               xi = i*h
c     due to left and right inner planes
               do kh = nlhp3, nhhp3
                  zkh = kh*h
                  do jh = nlhp2, nhhp2
                     yjh = jh*h

                     dist = sqrt( (xi - xhl)**2 + (yf - yjh)**2
     &                    + (zk - zkh)**2 )
                     frontc(i, k) = frontc(i, k)
c                     phi(i, nlp2, k) = phi(i, nlp2, k)
     &                    - phih(nlhp1, jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xhr)**2 + (yf - yjh)**2
     &                    + (zk - zkh)**2 )
                     frontc(i, k) = frontc(i, k)
c                     phi(i, nlp2, k) = phi(i, nlp2, k)
     &                    - phih(nhhp1, jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xhl)**2 + (yb - yjh)**2
     &                    + (zk - zkh)**2 )
                     backc(i, k) = backc(i, k)
c                     phi(i, nhp2, k) = phi(i, nhp2, k)
     &                    - phih(nlhp1, jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xhr)**2 + (yb - yjh)**2
     &                    + (zk - zkh)**2 )
                     backc(i, k) = backc(i, k)
c                     phi(i, nhp2, k) = phi(i, nhp2, k)
     &                    - phih(nhhp1, jh, kh)/(4.d0*pi*dist)
                  enddo
               enddo
c     due to front and back inner planes
               do kh = nlhp3, nhhp3, nref
                  zkh = kh*h
                  do ih = nlhp1 + 1, nhhp1 - 1, nref
                     xih = ih*h

                     dist = sqrt( (xi - xih)**2 + (yf - yhf)**2
     &                    + (zk - zkh)**2 )
                     frontc(i, k) = frontc(i, k)
c                     phi(i, nlp2, k) = phi(i, nlp2, k)
     &                    - phih(ih, nlhp2, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yf - yhb)**2
     &                    + (zk - zkh)**2 )
                     frontc(i, k) = frontc(i, k)
c                     phi(i, nlp2, k) = phi(i, nlp2, k)
     &                    - phih(ih, nhhp2, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yb - yhf)**2
     &                    + (zk - zkh)**2 )
                     backc(i, k) = backc(i, k)
c                     phi(i, nhp2, k) = phi(i, nhp2, k)
     &                    - phih(ih, nlhp2, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yb - yhb)**2
     &                    + (zk - zkh)**2 )
                     backc(i, k) = backc(i, k)
c                     phi(i, nhp2, k) = phi(i, nhp2, k)
     &                    - phih(ih, nhhp2, kh)/(4.d0*pi*dist)
                  enddo
               enddo
c     due to bottom and top inner planes
               do jh = nlhp2 + 1, nhhp2 - 1, nref
                  yjh = jh*h
                  do ih = nlhp1 + 1, nhhp1 - 1, nref
                     xih = ih*h

                     dist = sqrt( (xi - xih)**2 + (yf - yjh)**2
     &                    + (zk - zhb)**2 )
                     frontc(i, k) = frontc(i, k)
c                     phi(i, nlp2, k) = phi(i, nlp2, k)
     &                    - phih(ih, jh, nlhp3)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yf - yjh)**2
     &                    + (zk - zht)**2 )
                     frontc(i, k) = frontc(i, k)
c                     phi(i, nlp2, k) = phi(i, nlp2, k)
     &                    - phih(ih, jh, nhhp3)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yb - yjh)**2
     &                    + (zk - zhb)**2 )
                     backc(i, k) = backc(i, k)
c                     phi(i, nhp2, k) = phi(i, nhp2, k)
     &                    - phih(ih, jh, nlhp3)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yb - yjh)**2
     &                    + (zk - zht)**2 )
                     backc(i, k) = backc(i, k)
c                     phi(i, nhp2, k) = phi(i, nhp2, k)
     &                    - phih(ih, jh, nhhp3)/(4.d0*pi*dist)
                  enddo
               enddo

            enddo
         enddo

c     compute field for bottom and top outer planes
         do j = nlcp2, nhcp2
            yj = j*h
            do i = nlcp1, nhcp1
               xi = i*h
c     due to left and right inner planes
               do kh = nlhp3, nhhp3
                  zkh = kh*h
                  do jh = nlhp2, nhhp2
                     yjh = jh*h

                     dist = sqrt( (xi - xhl)**2 + (yj - yjh)**2
     &                    + (zb - zkh)**2 )
                     bottomc(i, j) = bottomc(i, j)
c                     phi(i, j, nlp3) = phi(i, j, nlp3)
     &                    - phih(nlhp1, jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xhr)**2 + (yj - yjh)**2
     &                    + (zb - zkh)**2 )
                     bottomc(i, j) = bottomc(i, j)
c                     phi(i, j, nlp3) = phi(i, j, nlp3)
     &                    - phih(nhhp1, jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xhl)**2 + (yj - yjh)**2
     &                    + (zt - zkh)**2 )
                     topc(i, j) = topc(i, j)
c                     phi(i, j, nhp3) = phi(i, j, nlp3)
     &                    - phih(nlhp1, jh, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xhr)**2 + (yj - yjh)**2
     &                    + (zt - zkh)**2 )
                     topc(i, j) = topc(i, j)
c                     phi(i, j, nhp3) = phi(i, j, nlp3)
     &                    - phih(nhhp1, jh, kh)/(4.d0*pi*dist)
                  enddo
               enddo
c     due to front and back inner planes
               do kh = nlhp3, nhhp3
                  zkh = kh*h
                  do ih = nlhp1 + 1, nhhp1 - 1
                     xih = ih*h

                     dist = sqrt( (xi - xih)**2 + (yj - yhf)**2
     &                    + (zb - zkh)**2 )
                     bottomc(i, j) = bottomc(i, j)
c                     phi(i, j, nlp3) = phi(i, j, nlp3)
     &                    - phih(ih, nlhp2, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yj - yhb)**2
     &                    + (zb - zkh)**2 )
                     bottomc(i, j) = bottomc(i, j)
c                     phi(i, j, nlp3) = phi(i, j, nlp3)
     &                    - phih(ih, nhhp2, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yj - yhf)**2
     &                    + (zt - zkh)**2 )
                     topc(i, j) = topc(i ,j)
c                     phi(i, j, nhp3) = phi(i, j, nhp3)
     &                    - phih(ih, nlhp2, kh)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yj - yhb)**2
     &                    + (zt - zkh)**2 )
                     topc(i ,j) = topc(i, j)
c                     phi(i, j, nhp3) = phi(i, j, nhp3)
     &                    - phih(ih, nhhp2, kh)/(4.d0*pi*dist)
                  enddo
               enddo
c     due to bottom and top inner planes
               do jh = nlhp2 + 1, nhhp2 - 1
                  yjh = jh*h
                  do ih = nlhp1 + 1, nhhp1 - 1
                     xih = ih*h

                     dist = sqrt( (xi - xih)**2 + (yj - yjh)**2
     &                    + (zb - zhb)**2 )
                     bottomc(i, j) = bottomc(i, j)
c                     phi(i, j, nlp3) = phi(i, j, nlp3)
     &                    - phih(ih, jh, nlhp3)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yj - yjh)**2
     &                    + (zb - zht)**2 )
                     bottomc(i, j) = bottomc(i, j)
c                     phi(i, j, nlp3) = phi(i, j, nlp3)
     &                    - phih(ih, jh, nhhp3)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yj - yjh)**2
     &                    + (zt - zhb)**2 )
                     topc(i, j) = topc(i, j)
c                     phi(i, j, nhp3) = phi(i, j, nhp3)
     &                    - phih(ih, jh, nlhp3)/(4.d0*pi*dist)

                     dist = sqrt( (xi - xih)**2 + (yj - yjh)**2
     &                    + (zt - zht)**2 )
                     topc(i, j) = topc(i, j)
c                     phi(i, j, nhp3) = phi(i, j, nhp3)
     &                    - phih(ih, jh, nhhp3)/(4.d0*pi*dist)
                  enddo
               enddo

            enddo
         enddo
      endif

c     Now interpolate from coarse to fine.
      
      call interp8(leftc, nlcp2, nlcp3, nhcp2, nhcp3, 
     &     left, nlp2, nlp3, nhp2, nhp3, h, nref)
      
      call interp8(rightc, nlcp2, nlcp3, nhcp2, nhcp3, 
     &     right, nlp2, nlp3, nhp2, nhp3, h, nref)

      call interp8(frontc, nlcp1, nlcp3, nhcp1, nhcp3,
     &     front, nlp1, nlp3, nhp1, nhp3, h, nref)

      call interp8(backc, nlcp1, nlcp3, nhcp1, nhcp3,
     &     back, nlp1, nlp3, nhp1, nhp3, h, nref)

      call interp8(bottomc, nlcp1, nlcp2, nhcp1, nhcp2,
     &     bottom, nlp1, nlp2, nhp1, nhp2, h, nref)

      call interp8(topc, nlcp1, nlcp2, nhcp1, nhcp2,
     &     top, nlp1, nlp2, nhp1, nhp2, h, nref)

c     And copy into the appropriate location

      i = nlp1
      do j = nlp2, nhp2
         do k = nlp3, nhp3
            phi(i, j, k) = left(j, k)
         enddo
      enddo
      i = nhp1
      do j = nlp2, nhp2
         do k = nlp3, nhp3
            phi(i, j, k) = right(j, k)
         enddo
      enddo
      j = nlp2
      do i = nlp1, nhp1
         do k = nlp3, nhp3
            phi(i, j, k) = front(i, k)
         enddo
      enddo
      j = nhp2
      do i = nlp1, nhp1
         do k = nlp3, nhp3
            phi(i, j, k) = back(i, k)
         enddo
      enddo
      k = nlp3
      do i = nlp1, nhp1
         do j = nlp2, nhp2
            phi(i, j, k) = bottom(i, j)
         enddo
      enddo
      k = nhp3
      do i = nlp1, nhp1
         do j = nlp2, nhp2
            phi(i, j, k) = top(i, j)
         enddo
      enddo

      return
      end



c     Do eighth-order interpolation from a coarse plane to a plane
c     (face) of a fine volume.

      subroutine interp8(coarse, clow1, clow2, chigh1, chigh2,
     &     fine, flow1, flow2, fhigh1, fhigh2, fine_h, nref)

      implicit none
      integer clow1, clow2, chigh1, chigh2
      integer flow1, flow2, fhigh1, fhigh2
      integer nref
      real*8 coarse(clow1:chigh1, clow2:chigh2)
      real*8 fine(flow1:fhigh1, flow2:fhigh2)
c     temporary array to hold results from interpolating in one dimension
      real*8 temp(0:nref, -3:4)
      real*8 fine_h, coarse_h
      real*8 f
      real*8 f1, f2, f3, f4, f5, f6, f7
      real*8 f1_2, f2_3, f3_4, f4_5, f5_6, f6_7
      real*8 f1_2_3, f2_3_4, f3_4_5, f4_5_6, f5_6_7
      real*8 f1_2_3_4, f2_3_4_5, f3_4_5_6, f4_5_6_7
      real*8 f1_2_3_4_5, f2_3_4_5_6, f3_4_5_6_7
      real*8 f1_2_3_4_5_6, f2_3_4_5_6_7
      real*8 f1_2_3_4_5_6_7
      real*8 x, xc0, xc1, xc2, xc3, xc4, xc5, xc6
      integer i, j, ii, jj

c     I think that all this could be viewed as a big matrix-matrix
c     operation that could certainly be computed more efficiently by
c     the blas library.

c     There is probably useless expansion of code here, but I don't try
c     for finesse in Fortran.

      do j = clow2 + 3, chigh2 - 4
         do i = clow1 + 3, chigh1 - 4

c     Intepolate in the current i direction
            do jj = -3, 4
               f = coarse(i, j + jj)
               
               f1 = (coarse(i+1, j + jj) - coarse(i, j + jj))
     &              /coarse_h
               f2 = - (coarse(i-1, j + jj)
     &              - coarse(i+1, j + jj))/(2.0*coarse_h)
               f3 = (coarse(i+2, j + jj)
     &              - coarse(i-1, j + jj))/(3.0*coarse_h)
               f4 = - (coarse(i-2, j + jj)
     &              - coarse(i+2, j + jj))/(4.0*coarse_h)
               f5 = (coarse(i+3, j + jj)
     &              - coarse(i-2, j + jj))/(5.0*coarse_h)
               f6 = - (coarse(i-3, j + jj)
     &              - coarse(i+3, j + jj))/(6.0*coarse_h)
               f7 = (coarse(i+4, j + jj)
     &              - coarse(i-3, j + jj))/(7.0*coarse_h)

               f1_2 = - (f2 - f1)/coarse_h
               f2_3 = (f3 - f2)/coarse_h
               f3_4 = - (f4 - f3)/coarse_h
               f4_5 = (f5 - f4)/coarse_h
               f5_6 = - (f6 - f5)/coarse_h
               f6_7 = (f7 - f6)/coarse_h
               
               f1_2_3 = (f2_3 - f1_2)/(2.0*coarse_h)
               f2_3_4 = - (f3_4 - f2_3)/(3.0*coarse_h)
               f3_4_5 = (f4_5 - f3_4)/(4.0*coarse_h)
               f4_5_6 = - (f5_6 - f4_5)/(5.0*coarse_h)
               f5_6_7 = (f6_7 - f5_6)/(6.0*coarse_h)
               
               f1_2_3_4 = - (f2_3_4 - f1_2_3)/(2.0*coarse_h)
               f2_3_4_5 = (f3_4_5 - f2_3_4)/(2.0*coarse_h)
               f3_4_5_6 = - (f4_5_6 - f3_4_5)/(2.0*coarse_h)
               f4_5_6_7 = (f5_6_7 - f4_5_6)/(2.0*coarse_h)
               
               f1_2_3_4_5 = (f2_3_4_5 - f1_2_3_4)/(3.0*coarse_h)
               f2_3_4_5_6 = - (f3_4_5_6 - f2_3_4_5)/(4.0*coarse_h)
               f3_4_5_6_7 = (f4_5_6_7 - f3_4_5_6)/(5.0*coarse_h)

               f1_2_3_4_5_6 = - (f2_3_4_5_6 - f1_2_3_4_5)
     &              /(3.0*coarse_h)
               f2_3_4_5_6_7 = (f3_4_5_6_7 - f2_3_4_5_6)
     &              /(3.0*coarse_h)

               f1_2_3_4_5_6_7 = (f2_3_4_5_6_7 - f1_2_3_4_5_6)
     &              /(4.0*coarse_h)            

               xc0 = 0.0
               xc1 = coarse_h
               xc2 = - coarse_h
               xc3 = 2.0*coarse_h
               xc4 = - 2.0*coarse_h
               xc5 = 3.0*coarse_h
               xc6 = - 3.0*coarse_h
               
               do ii = 0, nref
                  x = ii*fine_h
                  temp(ii, jj) = coarse(i, j + jj)
     &                 + (x-xc0)*f1 + (x-xc0)*(x-xc1)*f1_2
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*f1_2_3
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*f1_2_3_4
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*(x-xc4)
     &                 *f1_2_3_4_5
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*(x-xc4)
     &                 *(x-xc5)*f1_2_3_4_5_6
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*(x-xc4)
     &                 *(x-xc5)*(x-xc6)*f1_2_3_4_5_6_7
               enddo
            enddo

c     Now we have a set of (nref + 1)*8 interpolated coarse-grid values.
c     We need to interpolate these in the current j direction.
            do ii = 0, nref
               f = temp(ii, 0)
               f1 = (temp(ii, 1) - temp(ii, 0))
     &              /coarse_h
               f2 = - (temp(ii, -1) - temp(ii, 1))/(2.0*coarse_h)
               f3 = (temp(ii, 2) - temp(ii, -1))/(3.0*coarse_h)
               f4 = - (temp(ii, -2) - temp(ii, 2))/(4.0*coarse_h)
               f5 = (temp(ii, 3) - temp(ii, -2))/(5.0*coarse_h)
               f6 = - (temp(ii, -3) - temp(ii, 3))/(6.0*coarse_h)
               f7 = (temp(ii, 4) - temp(ii, -3))/(7.0*coarse_h)

               f1_2 = - (f2 - f1)/coarse_h
               f2_3 = (f3 - f2)/coarse_h
               f3_4 = - (f4 - f3)/coarse_h
               f4_5 = (f5 - f4)/coarse_h
               f5_6 = - (f6 - f5)/coarse_h
               f6_7 = (f7 - f6)/coarse_h
               
               f1_2_3 = (f2_3 - f1_2)/(2.0*coarse_h)
               f2_3_4 = - (f3_4 - f2_3)/(3.0*coarse_h)
               f3_4_5 = (f4_5 - f3_4)/(4.0*coarse_h)
               f4_5_6 = - (f5_6 - f4_5)/(5.0*coarse_h)
               f5_6_7 = (f6_7 - f5_6)/(6.0*coarse_h)
               
               f1_2_3_4 = - (f2_3_4 - f1_2_3)/(2.0*coarse_h)
               f2_3_4_5 = (f3_4_5 - f2_3_4)/(2.0*coarse_h)
               f3_4_5_6 = - (f4_5_6 - f3_4_5)/(2.0*coarse_h)
               f4_5_6_7 = (f5_6_7 - f4_5_6)/(2.0*coarse_h)
               
               f1_2_3_4_5 = (f2_3_4_5 - f1_2_3_4)/(3.0*coarse_h)
               f2_3_4_5_6 = - (f3_4_5_6 - f2_3_4_5)/(4.0*coarse_h)
               f3_4_5_6_7 = (f4_5_6_7 - f3_4_5_6)/(5.0*coarse_h)

               f1_2_3_4_5_6 = - (f2_3_4_5_6 - f1_2_3_4_5)
     &              /(3.0*coarse_h)
               f2_3_4_5_6_7 = (f3_4_5_6_7 - f2_3_4_5_6)
     &              /(3.0*coarse_h)

               f1_2_3_4_5_6_7 = (f2_3_4_5_6_7 - f1_2_3_4_5_6)
     &              /(4.0*coarse_h)            

               xc0 = 0.0
               xc1 = coarse_h
               xc2 = - coarse_h
               xc3 = 2.0*coarse_h
               xc4 = - 2.0*coarse_h
               xc5 = 3.0*coarse_h
               xc6 = - 3.0*coarse_h

               do jj = 0, nref
                  x = jj*fine_h
                  fine(i*nref + ii, j*nref + jj) = temp(ii, 0)
     &                 + (x-xc0)*f1 + (x-xc0)*(x-xc1)*f1_2
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*f1_2_3
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*f1_2_3_4
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*(x-xc4)
     &                 *f1_2_3_4_5
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*(x-xc4)
     &                 *(x-xc5)*f1_2_3_4_5_6
     &                 + (x-xc0)*(x-xc1)*(x-xc2)*(x-xc3)*(x-xc4)
     &                 *(x-xc5)*(x-xc6)*f1_2_3_4_5_6_7
               enddo
            enddo
         enddo
      enddo

      return
      end

