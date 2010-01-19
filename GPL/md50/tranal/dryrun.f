C    Sample trajectory analysis program for tranal v.5.0
C    Just read trajectories and report the time stamp 
C    and coordinates of the first atom  
      program DRYRUN
      include "tranal.h"
      call SETTRAJ
      IEND=0
      do while(IEND.eq.0)
        call READCONF(IEND)
        write(*,'(a,f12.3,a,3f10.3)')
     +  ' time',FULTIM,' 1 at:',SX(1),SY(1),SZ(1)
      end do
      write(*,*)IAN,' configurations analysed'
      stop
      end
