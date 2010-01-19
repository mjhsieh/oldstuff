       real function cputime( t )
       external getcpu
	call getcpu(t0)
        write(*,*)t0
        cputime=t0-t
        end
