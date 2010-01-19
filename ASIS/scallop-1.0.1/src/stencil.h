#ifndef _stencil_h_
#define _stencil_h_

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "scallop_macros.h"
#include "helpers.h"
#include "def3.h"
#include "p_grid3.h"

struct stencil {
    // A coarse box defining where the stencil lies.
    P_DOMAIN* coarse_domain;

    // A fine box defining where the stencil lies.
    P_DOMAIN* fine_domain;
    
    // A 3x3 coarse FAB to hold interpolation data.
    P_GRID* coarse_data;

    // A FAB for the fine interpolated data.
    P_GRID* fine_data;

    // Default constructor.
    stencil() {
	coarse_domain = NULL;
	fine_domain = NULL;
	coarse_data = NULL;
	fine_data = NULL;
    }

    // General purpose constructor.
    //stencil(const POINT& Point, int Nref, int Direction, double h);
    stencil(P_DOMAIN *input_coarse_domain, P_DOMAIN *input_fine_domain);

    // Destructor.
    ~stencil() {
	delete coarse_domain;
	delete fine_domain;
	delete coarse_data;
	delete fine_data;
    }

    // For debugging, takes the laplacian of the coarse data.
    double laplacian();

    // Re-instantiate coarse_data and set it to zero.
    void reset();
};

#endif
