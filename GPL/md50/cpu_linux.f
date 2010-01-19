       real*8 function cputime( t )
        real*8 t 
        real*4 r
        call getcpu( r )
        cputime = r-t
        end
