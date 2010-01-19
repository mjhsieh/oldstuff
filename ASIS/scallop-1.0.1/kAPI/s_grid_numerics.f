c      subroutine f_fill2d(u, ul0, ul1, uh0, uh1,
c     &     rl0,rl1,rh0,rh1, v)
c
c     fill a subregion of u with value v

      subroutine f_fill2d(u, ul0, ul1, uh0, uh1,
     &     rl0,rl1,rh0,rh1, v)
      integer ul0, ul1, uh0, uh1, rl0, rl1, rh0, rh1
      integer i, j
      real*8 u(ul0:uh0, ul1:uh1)
      real*8 v
      do j = rl1, rh1
         do i = rl0, rh0
            u(i,j) = v
         enddo
      enddo

      return
      end


c      subroutine f_fill3d(u,ul0, ul1, ul2, uh0, uh1, uh2,
c     &     rl0, rl1, rl2, rh0, rh1, rh2, v)	
c
c     fill a subregion of u with value v

      subroutine f_fill3d(u,ul0, ul1, ul2, uh0, uh1, uh2,
     &     rl0, rl1, rl2, rh0, rh1, rh2, v)	
      integer ul0, ul1, ul2, uh0, uh1, uh2
      integer rl0, rl1, rl2, rh0, rh1, rh2
      integer i, j, k
      real*8 u(ul0:uh0, ul1:uh1, ul2:uh2)
      real*8 v

      do k = rl2, rh2
         do j = rl1, rh1
            do i = rl0, rh0
               u(i,j,k) = v
            enddo
         enddo
      enddo
      
      end


c      subroutine f_negate2d(u, ul0, ul1, uh0, uh1)
c
c     negate u

      subroutine f_negate2d(u, ul0, ul1, uh0, uh1)
      integer ul0, ul1, uh0, uh1
      integer i, j
      real*8 u(ul0:uh0, ul1:uh1)

      do j = ul1, uh1
         do i = ul0, uh0
            u(i,j) = -u(i,j)
         enddo
      enddo

      end


c      subroutine f_negate3d(u, ul0, ul1, ul2, uh0, uh1, uh2)
c
c     negate u

      subroutine f_negate3d(u, ul0, ul1, ul2, uh0, uh1, uh2)
      integer ul0, ul1, ul2, uh0, uh1, uh2
      integer i, j, k
      real*8 u(ul0:uh0, ul1:uh1, ul2:uh2)

      do k = ul2, uh2
         do j = ul1, uh1
            do i = ul0, uh0
               u(i,j,k) = -u(i,j,k)
            enddo
         enddo
      enddo

      end


c      subroutine f_mult2d(u, ul0, ul1, uh0, uh1, s)
c
c     multiply array u by scalar s

      subroutine f_mult2d(u, ul0, ul1, uh0, uh1, s)
      integer ul0, ul1, uh0, uh1
      integer i, j
      real*8 u(ul0:uh0, ul1:uh1)
      real*8 s

      do j = ul1, uh1
         do i = ul0, uh0
            u(i,j) = u(i,j)*s
         enddo
      enddo

      end


c      subroutine f_mult3d(u, ul0, ul1, ul2, uh0, uh1, uh2, s)
c
c     multiply array u by scalar s

      subroutine f_mult3d(u, ul0, ul1, ul2, uh0, uh1, uh2, s)
      integer ul0, ul1, ul2, uh0, uh1, uh2
      integer i, j, k
      real*8 u(ul0:uh0, ul1:uh1, ul2:uh2)
      real*8 s

      do k = ul2, uh2
         do j = ul1, uh1
            do i = ul0, uh0
               u(i,j,k) = u(i,j,k)*s
            enddo
         enddo
      enddo

      end


c      subroutine f_norm2d(result, p, u, ul0, ul1, uh0, uh1)
c
c     return the p-norm of u

      subroutine f_norm2d(result, p, u, ul0, ul1, uh0, uh1)
      integer ul0, ul1, uh0, uh1
      integer i, j, p
      real*8 u(ul0:uh0, ul1:uh1)
      real*8 d, result

      result = 0.d0

      if (p .eq. 0) then
         do j = ul1, uh1
            do i = ul0, uh0
               d = dabs(u(i,j))
               if (d .gt. result) then
                  result = d
               endif
            enddo
         enddo
      else if (p .eq. 1) then
         do j = ul1, uh1
            do i = ul0, uh0
               d = dabs(u(i,j))
               result = result + d
            enddo
         enddo
      else if (p .eq. 2) then
         do j = ul1, uh1
            do i = ul0, uh0
               d = (u(i,j))**2
               result = result + d
            enddo
         enddo
         result = sqrt(result)
      else 
         do j = ul1, uh1
            do i = ul0, uh0
               d = (dabs(u(i, j)))**p
               result = result + d
            enddo
         enddo
         result = result**(1.d0/p)
      endif

      return
      end


c      subroutine f_norm3d(result, p, u, ul0, ul1, ul2, uh0, uh1, uh2)
c
c     return the p-norm of u

      subroutine f_norm3d(result, p, u, ul0, ul1, ul2, uh0, uh1, uh2)
      integer ul0, ul1, ul2, uh0, uh1, uh2
      integer i, j, k, p
      real*8 u(ul0:uh0, ul1:uh1, ul2:uh2)
      real*8 d, result

      result = 0.d0

      if (p .eq. 0) then
         do k = ul2, uh2
            do j = ul1, uh1
               do i = ul0, uh0
                  d = dabs(u(i,j,k))
                  if (d .gt. result) then
                     result = d
                  endif
               enddo
            enddo
         enddo
      else if (p .eq. 1) then
         do k = ul2, uh2
            do j = ul1, uh1
               do i = ul0, uh0
                  d = dabs(u(i,j,k))
                  result = result + d
               enddo
            enddo
         enddo
      else if (p .eq. 2) then
         do k = ul2, uh2
            do j = ul1, uh1
               do i = ul0, uh0
                  d = (u(i,j,k))**2
                  result = result + d
               enddo
            enddo
         enddo
         result = sqrt(result)
      else
         do k = ul2, uh2
            do j = ul1, uh1
               do i = ul0, uh0
                  d = (dabs(u(i,j,k)))**p
                  result = result + d
               enddo
            enddo
         enddo
         result = result**(1.d0/p)
      endif

      return
      end
