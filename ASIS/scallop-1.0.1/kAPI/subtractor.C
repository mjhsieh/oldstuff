#include "subtractorNDIM.h"

#define f_subtract F77_FUNC_(f_subtractNDIMd, F_SUBTRACTNDIMD)
#define FORTRAN_ARGS FORTRAN_ARGSNDIM
#define FortranRegion FortranRegionNDIM

extern "C" {
    void f_subtract(double *, FORTRAN_ARGS, const double *const, FORTRAN_ARGS,
		    FORTRAN_ARGS) ;
}


/*
void subtractorNDIM::unpack(s_gridNDIM& G, const RegionNDIM& R, double *buf)

subtract buf from s_grid g when called
*/
void subtractorNDIM::unpack(s_gridNDIM& G, const RegionNDIM& R, double *buf)
{
   FortranRegion g_extents(G.region());
   FortranRegion b_extents(R);
   f_subtract(FORTRAN_DATA(G), FORTRAN_REGIONNDIM(g_extents),
	      buf, FORTRAN_REGIONNDIM(b_extents), 
	      FORTRAN_REGIONNDIM(b_extents));
}


/*
void subtractorNDIM::copy(const s_gridNDIM& S, const RegionNDIM& F, s_gridNDIM& D, 
		      const RegionNDIM& T)

subtract s_grid D from s_grid S when called
*/
void subtractorNDIM::copy(const s_gridNDIM& S, const RegionNDIM& F, s_gridNDIM& D, 
		      const RegionNDIM& T)
{
   FortranRegion s_extents(S.region());
   FortranRegion d_extents(D.region());
   FortranRegion t_extents(T);

   f_subtract(FORTRAN_DATA(D), FORTRAN_REGIONNDIM(d_extents),
	      FORTRAN_DATA(S), FORTRAN_REGIONNDIM(s_extents),
	      FORTRAN_REGIONNDIM(t_extents));
}
