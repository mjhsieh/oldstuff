program print3DHyperplane
   implicit none
   integer xm,     ym,   zm
   integer xmymzm, xmym, nhpp
   integer ierr,   i
   integer, allocatable :: hppseq(:), totprevhpp(:)
   xm      = 4
   ym      = 4
   zm      = 4
   xmymzm  = xm*ym*zm
   xmym    = xm*ym
   nhpp    = xm+ym+zm-2
   ! The number of hyperplanes in a cube mesh
   ! +-+-+-+  In this 2d example xm = 4, ym = 3, total number of
   ! | | | |  hyperplanes in this example is 6, i.e. xm + ym - 1
   ! +-+-+-+
   ! | | | |
   ! +-+-+-+

   allocate( hppseq(xmymzm),  stat = ierr )
   !   All numbered grids sorted with hyperplane ordering
   allocate( totprevhpp(nhpp+1),stat = ierr )
   !   The total grids of all previous hyperplanes

   call hppidx(xm,ym,zm,xmymzm,xmym,nhpp,hppseq,totprevhpp)
   do i = 1, nhpp
      print *,hppseq(totprevhpp(i)+1:totprevhpp(i+1))
   enddo

   deallocate( hppseq,  stat = ierr )
   deallocate( totprevhpp,stat = ierr )
end program

subroutine hppidx(xm,ym,zm,xmymzm,xmym,nhpp,hppseq,totprevhpp)
!  This provides a 3-d hyperplane ordering sequence in hppseq.
!
!  To use it, refer to this example:
!       do i = 1, nhpp
!          print *,hppseq(totprevhpp(i)+1:totprevhpp(i+1))
!       enddo
!
!  Inspired by Ushiro Yasunori's code in http://netlib.org/ddsv
!  Author: Mengjuei Hsieh, University of California Irvine
   implicit none

   integer xm,ym,zm,xmymzm,xmym,nhpp
   integer i,j,k,mygrid,ierr
   integer hppseq(*),totprevhpp(*)
   integer, allocatable :: queryhpp(:)
   integer, allocatable :: lastinhpp(:)
   allocate( queryhpp(xmymzm),stat = ierr )
   ! Which hyperplane does this grid belong to.
   allocate( lastinhpp(nhpp), stat = ierr )
   ! The last grid's rank in the hp ordering

   totprevhpp(1) = 0
   lastinhpp     = 0

   ! H_m is the hyperplane m start from 1
   ! m = i + j + k - 2, 1 <= m <= xm+ym+zm-2
   do k =  0, zm-1
      do j = 0, ym-1
         mygrid=k*xmym+j*xm
         do i = 1, xm
            queryhpp(mygrid+i)=i+j+k
         enddo 
      enddo
   enddo
   do i = 1, xmymzm
      lastinhpp(queryhpp(i))=lastinhpp(queryhpp(i))+1
   enddo
   totprevhpp(2)=totprevhpp(1)+lastinhpp(1)
   lastinhpp(1) =totprevhpp(1)
   do i = 2, nhpp
      totprevhpp(i+1)=totprevhpp(i)+lastinhpp(i)
      lastinhpp(i)   =totprevhpp(i)
   enddo
   do i = 1, xmymzm
      lastinhpp(queryhpp(i))=lastinhpp(queryhpp(i))+1
      hppseq(lastinhpp(queryhpp(i)))=i
   enddo

   deallocate( queryhpp, stat = ierr )
   deallocate( lastinhpp,stat = ierr )
   return
end subroutine hppidx
