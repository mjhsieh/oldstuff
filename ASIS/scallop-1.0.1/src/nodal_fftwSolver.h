#ifndef _nodal_fftwsolver_h_
#define _nodal_fftwsolver_h_

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
#include "adder3.h"
#include "helpers.h"
#include "scallop_macros.h"

class nodal_fftwSolver {

protected:
    int debug_level;
    P_DOMAIN interior;
    P_GRID *phi;
    P_GRID *rhs;
    P_GRID *res;
    P_GRID *restmp;
    P_DOMAIN bcoarse;
    double h;

public:
    // Default constructor.
    nodal_fftwSolver();

    // Sets up storage for multigrid routine.
    nodal_fftwSolver(const P_DOMAIN& domain, double H, int debug = 0);

    // Acces routines.
    void set_rhs(P_GRID& r);
    void set_soln(P_GRID& s);
    P_GRID& get_rhs();
    P_GRID& get_soln();
    void set_domain(const P_DOMAIN& domain, double H);
    double& get_dx();

    // 19 point fftw solver.
    void solve19(int i);

    // 7 point fftw solver.
    void solve7(int i);

    ~nodal_fftwSolver();

};
#endif

