#ifndef _subtractorNDIM_h
#define _subtractorNDIM_h

/*

The subtractor class is an extension of the mover class, replacing copying
with subtraction

*/
#include "kelp.h"
#include "s_gridNDIM.h"
#include "VectorMoverNDIM.h"

class subtractorNDIM: public VectorMoverNDIM<s_gridNDIM,double>
{
   void unpack(s_gridNDIM& G, const RegionNDIM& R, double *buf);
   void copy(const s_gridNDIM& S, const RegionNDIM& F, s_gridNDIM& D, 
	     const RegionNDIM& T);
public:
   subtractorNDIM(XArrayNDIM<s_gridNDIM>& send, 
		  XArrayNDIM<s_gridNDIM>& recv, MotionPlanNDIM& S): 
      VectorMoverNDIM<s_gridNDIM,double>(send, recv, S) {}
};

#endif
