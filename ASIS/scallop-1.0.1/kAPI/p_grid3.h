/* *** This file is automatically generated by dimensions.sh from a ***
   *** dimension-independent file.  Edit at your peril.             *** */


#ifndef _p_grid3_h
#define _p_grid3_h

/*

p_grid implements a parallel grid: basically a more useful XArray

*/

#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "kelp.h"
#include "s_grid3.h"
#include "XArray3.h"
#include "VectorMover3.h"
#include "GhostPlan3.h"
#include "GhostIterator3.h"

class p_grid3: public XArray3<s_grid3>
{
   Point3 _ghost;	/* ghost cell vector			*/
   static int debug_counter;
public:
   MotionPlan3 fg_plan;	/* fill ghost region motionplan		*/
public:
   p_grid3(); 
   p_grid3(const FloorPlan3& F, const int build_fgp = TRUE);
   p_grid3(const GhostPlan3& F, const int build_fgp = TRUE);
   p_grid3(const p_grid3& IG, const int build_fgp = TRUE);
   ~p_grid3();

   virtual void instantiate(const FloorPlan3& F, const int build_fgp = TRUE);
   virtual void instantiate(const GhostPlan3& F, const int build_fgp = TRUE);
   virtual void instantiate(const p_grid3& IG, const int build_fgp = TRUE);
   virtual void build_fg_plan();

   const Point3& ghost() const {return _ghost;			}
   int ghost(const int dim) const {return _ghost(dim);		}

   Region3 interior(const int i)
   { return(grow((*this)(i).region(),-ghost()));		}

   void negate();
   void fill(const double v);
   void assign_ghost(const double v);
   void fill_ghost();
   void copy_on_intersect(p_grid3 &G);
   void CopyOnIntersection(p_grid3 &G);
   void copy(p_grid3& G);
   void copy(p_grid3& G, const FloorPlan3& F);
   void copy_grid(p_grid3& from_grid, int from_index, int to_index);
   void add_on_intersect(p_grid3 &G);
   void add_grid(p_grid3& from_grid, int from_index, int to_index);
   void plus(p_grid3& G);
   void plus_all(p_grid3& G);
   void minus(p_grid3& G);
   void mult(double s);
   double norm(int p);

   void debug_print(char *name);
};

#endif
