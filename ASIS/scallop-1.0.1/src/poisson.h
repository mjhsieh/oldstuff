#ifndef _poisson_h_
#define _poisson_h_

#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "scallop_macros.h"
#ifdef DIM2
#include "def2.h"
#include "p_domain2.h"
#include "p_grid2.h"
#include "adder2.h"
#include "subtractor2.h"
#include "io2.h"
#else
#include "def3.h"
#include "p_domain3.h"
#include "p_grid3.h"
#include "adder3.h"
#include "subtractor3.h"
#include "io3.h"
#endif
#include "stencil.h"
#include "nodal_infdmn.h"
#include "helpers.h"
#include "kelp.h"

struct poisson {

    // User input
    REGION fine_domain;		// Where the solution will live.
    P_DOMAIN fine_interior;	// The size of the interior.
    double fine_h;		// Grid spacing for the fine domain.
    int nref;			// The refinement ratio.
    int debug_level;		// Controls how much output to spew.

    P_DOMAIN soln_domain;		// The domain on which a solution
				// is wanted.
    P_DOMAIN coarse_domain;		// The size of Coarse->Soln.
    P_DOMAIN coarse_interior;		// The size of the interior of Coarse.
    
    P_GRID *input_rhs;	// The entire rhs.
    P_GRID *fine_rhs;	// The rhs for the initial solve.
    P_GRID *fine_soln;	// The P_GRID to hold the initial solve.
    P_GRID *coarse_rhs;	// ... to hold the coarse rhs from each patch.
    P_GRID *coarse_soln;	// ... to hold the coarse data on each patch.
    P_GRID *coarse_domain_rhs;  // ... to hold the entire coarse rhs.
    P_GRID *coarse_domain_soln; // ... to hold the entire coarse solution.

    stencil *s;			// Our lovely stencil.
    P_DOMAIN *sf_domain;	// P_DOMAIN describing stencil's fine domain.
    P_DOMAIN *sc_domain;	// ... describing stencil's coarse domain.
    int *stencil_patch_index;	// To assign patches to stencils.
    int *stencil_dir_array;	// To assign directions to stencils.
    int** lc_remote_reference;
    int** lf_remote_reference;
    P_DOMAIN *lc_coarse_domain;
    P_GRID *lc_coarse_soln;
    P_DOMAIN *lf_fine_domain;
    P_GRID *lf_fine_soln;

    nodal_fftwSolver *final_soln;	// For our perform the final solve.

    int c;			// The "correction radius."
    double coarse_h;               // Grid spacing for the coarse domain.


    // Default constructor.
    poisson();

    // General purpose constructor.
    poisson(const REGION& input_domain, const P_DOMAIN& input_patches, 
	double input_h, int input_nref, int debug = 0);

    // This sets the rhs to whatever the user says.
    int set_rhs(P_GRID& rhs);

    void fine_boundary_communication();
    // This sets the boundary of a fine patch for the final solve.
    void set_fine_boundary();

    // I want to be able to just say solve, and have it done.
    void solve();

    // This sets up the stencil data object for later use.
    void stencil_setup();
    

};

#endif
