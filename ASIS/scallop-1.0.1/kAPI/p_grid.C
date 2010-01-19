#include "p_gridNDIM.h"
#include "adderNDIM.h"
#include "subtractorNDIM.h"
#include "ioNDIM.h"


/* the compiler wouldn't accept these in the header file	*/
p_gridNDIM::p_gridNDIM(): fg_plan() {}
p_gridNDIM::~p_gridNDIM() {}


/*
p_gridNDIM::p_gridNDIM(const FloorPlanNDIM& F, const int build_fgp)

construct from a floorplan
*/
p_gridNDIM::p_gridNDIM(const FloorPlanNDIM& F, const int build_fgp)
{
   _ghost = 0;
   instantiate(F, build_fgp);
}


/*
p_gridNDIM::p_gridNDIM(const GhostPlanNDIM& F, const int build_fgp)

construct from a ghostplan
*/
p_gridNDIM::p_gridNDIM(const GhostPlanNDIM& F, const int build_fgp)
{
   instantiate(F,build_fgp);
}


/*
p_gridNDIM::p_gridNDIM(const p_gridNDIM& IG, const int build_fgp)

construct from another p_grid
*/
p_gridNDIM::p_gridNDIM(const p_gridNDIM& IG, const int build_fgp)
{
   instantiate(IG,build_fgp);
}


/*
void p_gridNDIM::instantiate(const FloorPlanNDIM& F,const int build_fgp)

instantiate from a floorplan
*/
void p_gridNDIM::instantiate(const FloorPlanNDIM& F,const int build_fgp)
{
   _ghost = 0;
   debug_counter = 0;
   XArrayNDIM<s_gridNDIM>::instantiate(F);
   fill(0.);
   if (build_fgp)
      build_fg_plan();
}


/*
void p_gridNDIM::instantiate(const GhostPlanNDIM& F,const int build_fgp)

instantiate from a ghostplan
*/
void p_gridNDIM::instantiate(const GhostPlanNDIM& F,const int build_fgp)
{
   _ghost = F.ghost();
   debug_counter = 0;
   XArrayNDIM<s_gridNDIM>::instantiate(F);
   fill(0.);
   if (build_fgp)
      build_fg_plan();
}


/*
void p_gridNDIM::instantiate(const p_gridNDIM& IG,const int build_fgp)

instantiate from another p_grid
*/
void p_gridNDIM::instantiate(const p_gridNDIM& IG,const int build_fgp)
{
   _ghost = IG.ghost();
   XArrayNDIM<s_gridNDIM>::instantiate(IG);
   fill(0.);
   if (build_fgp)
      build_fg_plan();
}


/*
void p_gridNDIM::build_fg_plan()

setup for filling ghost cells
*/
void p_gridNDIM::build_fg_plan()
{
   if (ghost() !=  PointNDIM(0)) {
    for_1(i,*this)
     RegionNDIM inside = interior(i);
     for_1(j,*this)
	if (i != j) {
	  fg_plan.CopyOnIntersection(*this,i,*this,j,inside);
	}
     end_for
    end_for
   }
}


/*
void p_gridNDIM::negate()

negate this whole p_grid
 */
void p_gridNDIM::negate()
{
    for_all(i, *this)
	(*this)(i).negate();
    end_for_all;
}


/*
void p_gridNDIM::fill(const double v)

fill the whole p_grid with v
*/
void p_gridNDIM::fill(const double v)
{
   for_all(i,*this)
      (*this)(i).fill(v);
   end_for_all
}


/*
void p_gridNDIM::assign_ghost(const double v)

fill ghost cells with assigned value v
*/
void p_gridNDIM::assign_ghost(const double v)
{
   for_all(i,*this)
      RegionNDIM R = (*this)(i).region();
      PointNDIM G = ghost();
      for (GhostIteratorNDIM Boundary(R,G); Boundary; ++Boundary) {
	const RegionNDIM B = Boundary();
	(*this)(i).fill(v,B);
      }
   end_for_all
}


/*
void p_gridNDIM::fillGhost()

fill ghost cells by copying
*/
void p_gridNDIM::fill_ghost()
{
   VectorMoverNDIM<s_gridNDIM,double> Ex(*this, *this, fg_plan);
   Ex.execute();
}


/*
void p_gridNDIM::copy_on_intersect(p_gridNDIM &G)

copy another p_grid on intersection
*/
void p_gridNDIM::copy_on_intersect(p_gridNDIM &G)
{
   MotionPlanNDIM S;

   for_1(i,G)
     for_1(j,*this)
       S.CopyOnIntersection(G, i, *this, j);
     end_for
   end_for

   VectorMoverNDIM<s_gridNDIM,double> Ex(G,*this,S);
   Ex.execute();
}


/*
void p_gridNDIM::CopyOnIntersection(p_gridNDIM &G)

copy another p_grid on intersection (an alias for copy_on_intersect)

This method will be removed from the API once I'm sure it is no longer 
necessary (i.e. purged from my main source code).
*/

void p_gridNDIM::CopyOnIntersection(p_gridNDIM &G)
{
    copy_on_intersect(G);
}


/*
void p_gridNDIM::copy(p_gridNDIM& G)

for each grid in a given p_grid, copy its contents to this p_grid
on intersection
*/
void p_gridNDIM::copy(p_gridNDIM& G)
{
    MotionPlanNDIM M;
    for_1(i, G) 
	M.CopyOnIntersection(G, i, *this, i);
    end_for;
    VectorMoverNDIM<s_gridNDIM,double> Ex(G, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::copy(p_gridNDIM& G, const FloorPlanNDIM& F)

for each grid in a given p_grid, copy its contents to this p_grid
on the intersection with a given floorplan
*/
void p_gridNDIM::copy(p_gridNDIM& G, const FloorPlanNDIM& F)
{
    MotionPlanNDIM M;
    for_1(i, G) {
	M.CopyOnIntersection(G, i, *this, i, F(i));
    } end_for;
    VectorMoverNDIM<s_gridNDIM,double> Ex(G, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::copy_grid(p_gridNDIM& from_grid, 
		       int from_index, int to_index)

copy one grid from a given p_grid to a specific grid in this p_grid,
on intersection
*/
void p_gridNDIM::copy_grid(p_gridNDIM& from_grid, 
		       int from_index, int to_index)
{
    MotionPlanNDIM M;
    M.CopyOnIntersection(from_grid, from_index, *this, to_index);

    VectorMoverNDIM<s_gridNDIM,double> Ex(from_grid, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::add_on_intersect(p_gridNDIM &G)

add values from a given p_grid to this p_grid on intersection
*/
void p_gridNDIM::add_on_intersect(p_gridNDIM &G)
{
   MotionPlanNDIM S;

   for_1(i,G)
     for_1(j,*this)
       S.CopyOnIntersection(G,i,*this,j);
     end_for
   end_for

   adderNDIM Ex(G,*this,S);
   Ex.execute();
}


/*
void p_gridNDIM::add_grid(p_gridNDIM& from_grid, 
		     int from_index, int to_index)

add values from a single grid in a given p_grid to a specific grid
in this p_grid, on intersection
*/
void p_gridNDIM::add_grid(p_gridNDIM& from_grid, 
		     int from_index, int to_index)
{
    MotionPlanNDIM M;
    M.CopyOnIntersection(from_grid, from_index, *this, to_index);

    adderNDIM Ex(from_grid, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::plus(p_gridNDIM& G)

add values from each grid of a given p_grid to corresponding grids in
this p_grid, on intersection
*/
void p_gridNDIM::plus(p_gridNDIM& G)
{
    MotionPlanNDIM M;
    for_1(i, G) {
	M.CopyOnIntersection(G, i, *this, i);
    } end_for;
    adderNDIM Ex(G, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::plus_all(p_gridNDIM& G)

add on intersection
*/
void p_gridNDIM::plus_all(p_gridNDIM& G)
{
    MotionPlanNDIM M;
    for_1(i, G) {
	for_1(j, *this) {
	    M.CopyOnIntersection(G, i, *this, j);
	} end_for;
    } end_for;
    adderNDIM Ex(G, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::minus(p_gridNDIM& G)

subtract values from each grid of a given p_grid from corresponding grids in
this p_grid, on intersection
*/
void p_gridNDIM::minus(p_gridNDIM& G)
{
    MotionPlanNDIM M;
    for_1(i, G) {
	M.CopyOnIntersection(G, i, *this, i);
    } end_for;
    subtractorNDIM Ex(G, *this, M);
    Ex.execute();
}


/*
void p_gridNDIM::mult(double s)

multiply all values in this p_grid by s
*/
void p_gridNDIM::mult(double s)
{
    for_all(i, *this) {
	(*this)(i).mult(s);
    } end_for_all;
}


/* 
double p_gridNDIM::norm(int p)

take the p-norm of this p_grid
*/
double p_gridNDIM::norm(int p)
{
    double result = 0.0;
    double tmpresult = 0.0;

    if (p == 0) {
	for_all(i, *this) {
	    result = MAX((*this)(i).norm(p), result);
	} end_for_all;
	mpReduceMax(&result, 1);
    }
    else if (p == 1) {
	for_all(i, *this) {
	    result += (*this)(i).norm(p);
	} end_for_all;
	mpReduceAdd(&result, 1);
    }
    else if (p == 2) {
	for_all(i, *this) {
	    tmpresult = (*this)(i).norm(p);
	    tmpresult *= tmpresult;
	    result += tmpresult;
	} end_for_all;
	mpReduceAdd(&result, 1);
	result = sqrt(result);
    }

    return result;
}


/*
void p_gridNDIM::debug_print(char *name)

write out the p_grid values for debugging purposes
*/
void p_gridNDIM::debug_print(char *name)
{
    char *file_name = new char[64];
    char *node_name =new char[8];
    ofstream debug_file;
    for (nodeIterator ni(this->floorplan()); ni; ++ni) {
	file_name = strncpy(file_name, name, 55);
	sprintf(node_name, "%03d.%03d", ni()%1000, debug_counter++);
	debug_file.open(strncat(file_name, node_name, 7), ios::out);
	debug_file << operator()(ni());
	debug_file.close();
    }
    delete file_name;
    delete node_name;
}
    
