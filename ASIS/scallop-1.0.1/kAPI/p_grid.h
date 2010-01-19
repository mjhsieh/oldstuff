#ifndef _p_gridNDIM_h
#define _p_gridNDIM_h

/*

p_grid implements a parallel grid: basically a more useful XArray

*/

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "kelp.h"
#include "s_gridNDIM.h"
#include "XArrayNDIM.h"
#include "VectorMoverNDIM.h"
#include "GhostPlanNDIM.h"
#include "GhostIteratorNDIM.h"

class p_gridNDIM: public XArrayNDIM<s_gridNDIM>
{
   PointNDIM _ghost;	/* ghost cell vector			*/
   static int debug_counter;
public:
   MotionPlanNDIM fg_plan;	/* fill ghost region motionplan		*/
public:
   p_gridNDIM(); 
   p_gridNDIM(const FloorPlanNDIM& F, const int build_fgp = TRUE);
   p_gridNDIM(const GhostPlanNDIM& F, const int build_fgp = TRUE);
   p_gridNDIM(const p_gridNDIM& IG, const int build_fgp = TRUE);
   ~p_gridNDIM();

   virtual void instantiate(const FloorPlanNDIM& F, const int build_fgp = TRUE);
   virtual void instantiate(const GhostPlanNDIM& F, const int build_fgp = TRUE);
   virtual void instantiate(const p_gridNDIM& IG, const int build_fgp = TRUE);
   virtual void build_fg_plan();

   const PointNDIM& ghost() const {return _ghost;			}
   int ghost(const int dim) const {return _ghost(dim);		}

   RegionNDIM interior(const int i)
   { return(grow((*this)(i).region(),-ghost()));		}

   void negate();
   void fill(const double v);
   void assign_ghost(const double v);
   void fill_ghost();
   void copy_on_intersect(p_gridNDIM &G);
   void CopyOnIntersection(p_gridNDIM &G);
   void copy(p_gridNDIM& G);
   void copy(p_gridNDIM& G, const FloorPlanNDIM& F);
   void copy_grid(p_gridNDIM& from_grid, int from_index, int to_index);
   void add_on_intersect(p_gridNDIM &G);
   void add_grid(p_gridNDIM& from_grid, int from_index, int to_index);
   void plus(p_gridNDIM& G);
   void plus_all(p_gridNDIM& G);
   void minus(p_gridNDIM& G);
   void mult(double s);
   double norm(int p);

   void debug_print(char *name);
};

#endif
