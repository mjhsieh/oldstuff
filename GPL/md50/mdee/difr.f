      dimension Y(100)
      I=1
 10   read(*,*,end=20)i1,d2,d3,d4,Y(I)
      I=I+1
      go to 10
 20   I=I-1
      do K=2,I
         write(*,*)K-1,Y(K)-Y(K-1)
      end do
      end
