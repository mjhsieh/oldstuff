       real*8 function cputime( t )
        real*8 t 
        real*4 t0,r
        t0=t
        call getcpu( t0, r )
	cputime = r
        end
