C
C   This subroutine substitute 2,3,4 columns of the firts file
C   with 2,3,4 columns of the second file
C   (useful to rewtite optimized coordinates for .smol fiel)
C
      character*4 CH,CH1,NM
      character*32 f1,f2,f3
 5    write(*,*)' first file (fragment of .smol)'
      read(*,*)f1
      open(unit=9,file=f1,status='old',err=10)
      go to 15
 10   write(*,*)' file ',f1,' not found'
      go to 5
 15   write(*,*)' second file (.xmol fragment)'
      read(*,*)f2
      open(unit=10,file=f2,status='old',err=20)
      go to 25
 20   write(*,*)' file ',f1,' not found'
      go to 15
 25   write(*,*)' new file '
      read(*,*)f3
      open(unit=11,file=f3,status='unknown')
 30   read(9,*,err=50,end=50)CH,X,Y,Z,Q,NM
      read(10,*,err=40,end=40)CH1,X1,Y1,Z1
      write(11,'(a4,3f9.4,f9.5,2x,a4)')CH,X1,Y1,Z1,Q,NM
      write(*,'(a4,3f9.4,f9.5,2x,a4)')CH,X1,Y1,Z1,Q,NM
      go to 30
 40   write(*,*)' unexpected end of ',f1,' file'
 50   stop
      end
