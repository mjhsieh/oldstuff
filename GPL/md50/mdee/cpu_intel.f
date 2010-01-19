       real*8 function cputime( t )
       real*4 r
       real*8 t
        call CPU_TIME( r )
        cputime = r-t
        end
