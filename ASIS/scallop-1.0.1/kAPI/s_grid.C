#include "s_gridNDIM.h"

#define f_fill F77_FUNC_(f_fillNDIMd, F_FILLNDIMD)
#define FORTRAN_ARGS FORTRAN_ARGSNDIM
#define FortranRegion FortranRegionNDIM

extern "C" {
    void f_fill(const double* const, FORTRAN_ARGS,
		FORTRAN_ARGS, const double* const);
}

void s_gridNDIM::fill(const double val, const RegionNDIM& R)
{
   FortranRegion FR(R);
   FortranRegion Fthis(this->region());
   f_fill(FORTRAN_DATA(*this),FORTRAN_REGIONNDIM(Fthis),FORTRAN_REGIONNDIM(FR),&val);
};

/* 
void s_gridNDIM::negate()

negate all the values in the grid 
*/
#define f_negate F77_FUNC_(f_negateNDIMd, F_NEGATENDIMD)

extern "C" {
    void f_negate(const double* const u, FORTRAN_ARGS);
}

void s_gridNDIM::negate()
{
   FortranRegion Fthis(this->region());
   f_negate(FORTRAN_DATA(*this),FORTRAN_REGIONNDIM(Fthis));
};


/* 
void s_gridNDIM::mult(double s)

multiply all the values in the grid by s
*/
#define f_mult F77_FUNC_(f_multNDIMd, F_MULTNDIMD)


extern "C" {
    void f_mult(const double* const u, FORTRAN_ARGS, double* s);
}

void s_gridNDIM::mult(double s)
{
   FortranRegion Fthis(this->region());
   f_mult(FORTRAN_DATA(*this), FORTRAN_REGIONNDIM(Fthis), &s);
};


/*
double s_gridNDIM::norm(const int p)

returns the p-norm of this grid.
*/
#define f_norm F77_FUNC_(f_normNDIMd, F_NORMNDIMD)

extern "C" {
    double f_norm(double* result, const int* p,
		  const double* const u, FORTRAN_ARGS);
}

double s_gridNDIM::norm(const int p)
{
   double result;

   FortranRegion Fthis(this->region());

   f_norm(&result, &p, FORTRAN_DATA(*this), FORTRAN_REGIONNDIM(Fthis));

   return(result);
};
