	real*8 function cputime( t )
	real*8 t, r
	real*4 et(2)
	r = etime_(et) - t
	cputime = r
	end


