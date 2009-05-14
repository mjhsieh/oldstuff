program print3DHyperplane
   implicit none
   integer xm,     ym,   zm
   integer xmymzm, xmym, nhpp
   integer ierr, i,j
   integer, allocatable :: hpporder(:), totprevhpp(:)
   xm      = 4
   ym      = 4
   zm      = 4
   xmymzm  = xm*ym*zm
   xmym    = xm*ym
   nhpp = xm+ym+zm-2 ! The number of hyperplanes in a cube mesh

   ! +-+-+-+  In this 2d example xm = 4, ym = 3, total number of
   ! | | | |  hyperplanes in this example is 6, i.e. xm + ym - 1
   ! +-+-+-+
   ! | | | |
   ! +-+-+-+

   !   hpporder(xmymzm),whichhpp(xmymzm),totprevhpp(nhpp+1),lastinhpp(nhpp)
   allocate( hpporder(xmymzm),  stat = ierr )
   !   All numbered grids sorted with hyperplane ordering
   allocate( totprevhpp(nhpp+1),stat = ierr )
   !   The total grids of all previous hyperplanes

   call hppidx(xm,ym,zm,xmymzm,xmym,nhpp,hpporder,totprevhpp)
   do i = 1, nhpp
      print *,hpporder(totprevhpp(i)+1:totprevhpp(i+1))
   enddo

   deallocate( hpporder,  stat = ierr )
   deallocate( totprevhpp, stat = ierr )
end program

subroutine hppidx(xm,ym,zm,xmymzm,xmym,nhpp,hpporder,totprevhpp)
!  This provides a 3-d hyperplane ordering sequence in hpporder.
!
!  To use it, refer to this example:
!       do i = 1, nhpp
!          print *,hpporder(totprevhpp(i)+1:totprevhpp(i+1))
!       enddo
!
!  Inspired by Ushiro Yasunori from http://www.netlib.org/ddsv
!  Author: Mengjuei Hsieh, University of California Irvine
   implicit none

   integer xm,ym,zm,xmymzm,xmym,nhpp
   integer i,j,k,mygrid,ierr
   integer hpporder(*),totprevhpp(*)
   integer, allocatable :: whichhpp(:)
   integer, allocatable :: lastinhpp(:)
   allocate( whichhpp(xmymzm),stat = ierr )
   ! Which hyperplane does this grid belong to.
   allocate( lastinhpp(nhpp), stat = ierr )
   ! The last grid's rank in the hp ordering

   totprevhpp(1) = 0
   lastinhpp     = 0

   do k=0, zm-1
      do j=0, ym-1
         mygrid=k*xmym+j*xm
         do i=1, xm
            whichhpp(mygrid+i)=i+j+k
         enddo 
      enddo
   enddo
   do i=1, xmymzm
      lastinhpp(whichhpp(i))=lastinhpp(whichhpp(i))+1
   enddo
   totprevhpp(2)=totprevhpp(1)+lastinhpp(1)
   lastinhpp(1)=totprevhpp(1)
   do i=2, nhpp
      totprevhpp(i+1)=totprevhpp(i)+lastinhpp(i)
      lastinhpp(i)=totprevhpp(i)
   enddo
   do i=1, xmymzm
      lastinhpp(whichhpp(i))=lastinhpp(whichhpp(i))+1
      hpporder(lastinhpp(whichhpp(i)))=j
   enddo

   deallocate( whichhpp, stat = ierr )
   deallocate( lastinhpp,stat = ierr )
   return
end subroutine hppidx
