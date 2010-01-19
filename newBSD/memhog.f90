!
! Author: Mengjuei Hsieh, University of California, Irvine.
!
! Not yet complete

data CHUNK /128 * 1024/

program main
   real*8, dimension(:), allocatable :: a, b
   integer ier
   allocate ( a(200000000), stat = ier )
   if ( ier /= 0 ) then
      print *,"allocation error"
      deallocate ( a, stat=ier )
      if ( ier /= 0 ) print *,"deallocation error"
      stop
   else
      print *,"allocation success"
   endif
   allocate ( b(100000000), stat = ier )
   if ( ier /= 0 ) then
      print *,"allocation error"
      deallocate ( b, stat=ier )
      if ( ier /= 0 ) print *,"deallocation error"
      stop
   else
      print *,"allocation success"
   endif
   deallocate ( a, stat=ier )
   if ( ier /= 0 ) then
      print *,"deallocation error"
      stop
   else
      print *,"deallocation success"
   endif
   deallocate ( b, stat=ier )
   if ( ier /= 0 ) then
      print *,"deallocation error"
      stop
   else
      print *,"deallocation success"
   endif
end program

! If ifc (ver 7) does not compile, use -lPEPCF90 in the compiler command.
subroutine argumentchk
       implicit none
       integer   :: IArgC,ArgC,i
       character :: ArgV(4)*20

       ArgC = IArgC()
       if (ArgC==4) then
          do i=1,ArgC
             call GetArg(i,ArgV(i))
             write(6,*) ArgV(i)
          enddo
       else
          write(6,*) "Usage: mytest arg1 arg2 arg3 arg4"
          stop
       endif
end subroutine
