c      subroutine f_add2d(u, ul0, ul1, uh0, uh1,
c     &     e, el0, el1, eh0, eh1, wl0, wl1, wh0, wh1)
c
c     add e to u on region w

      subroutine f_add2d(u, ul0, ul1, uh0, uh1,
     &     e, el0, el1, eh0, eh1, wl0, wl1, wh0, wh1)
      integer ul0, ul1, uh0, uh1, el0, el1, eh0, eh1
      integer wl0, wl1, wh0, wh1
      integer i, j
      real*8 u(ul0:uh0, ul1:uh1)
      real*8 e(el0:eh0, el1:eh1)

      do j = wl1, wh1
         do i = wl0, wh0
            u(i,j) = u(i,j) + e(i,j)
         end do
      end do

      return
      end


c      subroutine f_add3d(u, ul0, ul1, ul2, uh0, uh1, uh2,
c     &     e, el0, el1, el2, eh0, eh1, eh2,
c     &     wl0, wl1, wl2, wh0, wh1, wh2)
c
c     add e to u on region w

      subroutine f_add3d(u, ul0, ul1, ul2, uh0, uh1, uh2,
     &     e, el0, el1, el2, eh0, eh1, eh2,
     &     wl0, wl1, wl2, wh0, wh1, wh2)
      integer ul0, ul1, ul2, uh0, uh1, uh2
      integer el0, el1, el2, eh0, eh1, eh2
      integer wl0, wl1, wl2, wh0, wh1, wh2
      integer i, j, k 
      real*8 u(ul0:uh0,ul1:uh1,ul2:uh2)
      real*8 e(el0:eh0,el1:eh1,el2:eh2)

      do k = wl2, wh2
         do j = wl1, wh1
            do i = wl0, wh0
               u(i,j,k) = u(i,j,k) + e(i,j,k)
            end do
         end do
      end do

      return
      end
