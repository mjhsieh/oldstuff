*   This procedure get elapsed cpu time using "etime" function
*   (e.g. Portland Group compiler)
*      
       real*8 function cputime( t )
       real*8 t
       real*4 etime,u,s
*        call getcpu( r )
        cputime = etime(u,s)-t
        end
