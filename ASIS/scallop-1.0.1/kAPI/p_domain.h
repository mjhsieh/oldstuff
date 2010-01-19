#ifndef _p_domainNDIM_h
#define _p_domainNDIM_h

/*

This is just a floorplan with a couple extra useful constructs for things I
do a lot.

*/

#include "kelp.h"
#include "dock.h"
#include "GridNDIM.h"
#include "XArrayNDIM.h"
#include "FloorPlanNDIM.h"

class p_domainNDIM: public FloorPlanNDIM
{
 public:
    p_domainNDIM() {}
    p_domainNDIM(const int n) : FloorPlanNDIM(n) {}
    p_domainNDIM(const FloorPlanNDIM& f) : FloorPlanNDIM(f) {}
    p_domainNDIM(const p_domainNDIM& b) : FloorPlanNDIM(b) {}

    void grow(const int index, const int dim) {
	PointNDIM p(dim);
	(*this).FloorPlanNDIM::grow(index, p);
    }
    void grow(const int dim) {
	for_1(i, *this) 
	    (*this).grow(i, i);
	end_for;
    }
    void coarsen(const int nref) {
	PointNDIM lo;
	PointNDIM hi;
	for_1(i, *this) 
	    lo = (*this)(i).lower();
	    hi = (*this)(i).upper();
	    for (int dim = 0; dim < NDIM; dim++) {
	      lo(dim) = ( (lo(dim) < 0) ? ( (lo(dim) + 1)/nref - 1 ) 
			  : lo(dim)/nref );
	      hi(dim) = ( (hi(dim) < 0) ? ( (hi(0) + 1)/nref - 1 ) 
			  : hi(dim)/nref );
	    }
	    setlower(i, lo);
	    setupper(i, hi);
	end_for;
    }
    void refine(const int nref) {
	PointNDIM lo;
	PointNDIM hi;
	for_1(i, *this) 
	    lo = (*this)(i).lower();
	    hi = (*this)(i).upper();
	    lo *= nref;
	    hi = (hi+1)*nref - 1;
	    setlower(i, lo);
	    setupper(i, hi);
	end_for;
    }
    void nrefine(const int nref) {
	PointNDIM lo;
	PointNDIM hi;
	for_1(i, *this) 
	    lo = (*this)(i).lower();
	    hi = (*this)(i).upper();
	    lo *= nref;
	    hi *= nref;
	    setlower(i, lo);
	    setupper(i, hi);
	end_for;
    }
    p_domainNDIM& operator=(const p_domainNDIM& b) {
	resize(b.size());
	for_1(i, b) 
	    setowner(i, b(i).owner());
	    setregion(i, b(i));
	end_for;
	return *this;
    }
};

p_domainNDIM grow(const p_domainNDIM& b, int dim);
p_domainNDIM* grow(const p_domainNDIM* b, int dim);
#endif
