#ifndef _infdmn_h_
#define _infdmn_h_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "kelp.h"
#include "dock.h"
#include "def3.h"
#include "s_grid3.h"
#include "io3.h"
#include "p_grid3.h"
#include "nodal_fftwSolver.h"
#include "helpers.h"
#include "scallop_macros.h"

#ifdef DEBUG
#include "util.h"
#endif

/*
Header file for infinite domain solver for multigrid. Inputs are 
rhs : FAB defining right-hand side. Includes no border points.
soln : FAB on which the solution is required . soln has allocated one
	cell wide of border points. The interior of the index set that
	soln is defined on must be strictly larger that that of rhs.
*/

// The following is similar to infdmnSolve, but it is uses a 
// fourth-order method to get a more accurate solution on the coarse
// domain.
void nodal_id_solve(P_GRID *rhs, P_GRID *soln, int nref, double h, 
		    int debug = 0);

#endif
