#include "adderNDIM.h"

#define f_add F77_FUNC_(f_addNDIMd, F_ADDNDIMD)
#define FORTRAN_ARGS FORTRAN_ARGSNDIM
#define FortranRegion FortranRegionNDIM

extern "C" {
    void f_add(double *, FORTRAN_ARGS, const double *const, FORTRAN_ARGS,
	       FORTRAN_ARGS);
}


/*
void adderNDIM::unpack(s_gridNDIM& G, const RegionNDIM& R, double *buf)

add buf to s_grid g when called
*/
void adderNDIM::unpack(s_gridNDIM& G, const RegionNDIM& R, double *buf)
{
   FortranRegion g_extents(G.region());
   FortranRegion b_extents(R);
   f_add(FORTRAN_DATA(G), FORTRAN_REGIONNDIM(g_extents),
	 buf, FORTRAN_REGIONNDIM(b_extents), FORTRAN_REGIONNDIM(b_extents));
}


/*
void adderNDIM::copy(const s_gridNDIM& S, const RegionNDIM& F, s_gridNDIM& D, const RegionNDIM& T)

add s_grid D to s_grid S when called
*/
void adderNDIM::copy(const s_gridNDIM& S, const RegionNDIM& F, s_gridNDIM& D, const RegionNDIM& T)
{
   FortranRegion s_extents(S.region());
   FortranRegion d_extents(D.region());
   FortranRegion t_extents(T);

   f_add(FORTRAN_DATA(D), FORTRAN_REGIONNDIM(d_extents),
	 FORTRAN_DATA(S), FORTRAN_REGIONNDIM(s_extents),
         FORTRAN_REGIONNDIM(t_extents));
}
