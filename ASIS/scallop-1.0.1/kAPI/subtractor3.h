/* *** This file is automatically generated by dimensions.sh from a ***
   *** dimension-independent file.  Edit at your peril.             *** */


#ifndef _subtractor3_h
#define _subtractor3_h

/*

The subtractor class is an extension of the mover class, replacing copying
with subtraction

*/
#include "kelp.h"
#include "s_grid3.h"
#include "VectorMover3.h"

class subtractor3: public VectorMover3<s_grid3,double>
{
   void unpack(s_grid3& G, const Region3& R, double *buf);
   void copy(const s_grid3& S, const Region3& F, s_grid3& D, 
	     const Region3& T);
public:
   subtractor3(XArray3<s_grid3>& send, 
		  XArray3<s_grid3>& recv, MotionPlan3& S): 
      VectorMover3<s_grid3,double>(send, recv, S) {}
};

#endif
