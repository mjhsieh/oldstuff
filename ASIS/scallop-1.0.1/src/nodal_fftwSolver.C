#include "nodal_fftwSolver.h"

#define FORT_POISSON3D F77_FUNC( poisson3d , POISSON3D )

extern "C" {
    void FORT_POISSON3D(double*,int*,double*,int*);
};

// Default constructor. Initializes pointers to NULL. This is for situations
// (e.g. arrays of multigrid objects) in which you want to declare the 
// variables separately from setting up the storage. To set up the
// storage, use the member function multigrid::setDomain defined below.
nodal_fftwSolver::nodal_fftwSolver() {
    phi = NULL;
    rhs = NULL;
    res = NULL;
    restmp = NULL;
    h = 0.;
}

// Constructor.  Sets up memory-managed storage for multigrid object.
nodal_fftwSolver::nodal_fftwSolver(const P_DOMAIN& given_interior, double H,
				 int debug) {
    debug_level = debug;
    phi = NULL;
    rhs = NULL;
    res = NULL;
    restmp = NULL;
    set_domain(given_interior, H);
}


// Destructor. We need one, since we have pointers to memory as 
// part of the data.
nodal_fftwSolver::~nodal_fftwSolver() {
    delete phi;
    delete rhs;
    delete res;
    delete restmp;
}

// SetDomain sets up memory-managed storage for multigrid object
// previously defined using the default constructor.
void nodal_fftwSolver::set_domain(const P_DOMAIN& given_interior, double H) {
    if (phi)
	delete phi;
    if (rhs)
	delete rhs;
    if (res)
	delete res;
    if (restmp)
	delete restmp;
    interior.resize(given_interior.size());
    for_1(i, interior) {
	interior.setowner(i, given_interior(i).owner());
	interior.setregion(i, given_interior(i));
    } end_for;

    phi = new P_GRID(interior);
    rhs = new P_GRID(interior);
    res = new P_GRID(interior);

    h = H;
}

void nodal_fftwSolver::solve19(int i) {
    // Call fftsolver for 19-point stencil.
    int opType=1;
    int nx = phi->ptr(i)->domain.u0 - phi->ptr(i)->domain.l0;
    int ny = phi->ptr(i)->domain.u1 - phi->ptr(i)->domain.l1;
    int nz = phi->ptr(i)->domain.u2 - phi->ptr(i)->domain.l2;
    assert ((nx == ny) && (nx  == nz));
    for (int kk = 0;kk < (nx+1)*(ny+1)*(nz+1);kk++)
    {   
        phi->ptr(i)->data()[kk] = rhs->ptr(i)->data()[kk];
    }
    FORT_POISSON3D(phi->ptr(i)->data(),&nx,&h,&opType);
    
}

void nodal_fftwSolver::solve7(int i) {
    // Call FFT solver for 7-point stencil.
    int opType=2;
    int nx = phi->ptr(i)->domain.u0 - phi->ptr(i)->domain.l0;
    int ny = phi->ptr(i)->domain.u1 - phi->ptr(i)->domain.l1;
    int nz = phi->ptr(i)->domain.u2 - phi->ptr(i)->domain.l2;
    assert ((nx == ny) && (nx  == nz));
    for (int kk = 0;kk < (nx+1)*(ny+1)*(nz+1);kk++)
    {   
        phi->ptr(i)->data()[kk] = rhs->ptr(i)->data()[kk];
    }
    FORT_POISSON3D(phi->ptr(i)->data(),&nx,&h,&opType);
    
}



void nodal_fftwSolver::set_rhs(P_GRID& r) {
    rhs->copy(r);
}

void nodal_fftwSolver::set_soln(P_GRID& s) {
    phi->copy(s);
}

double& nodal_fftwSolver::get_dx() {return h;}

P_GRID& nodal_fftwSolver::get_rhs() {return *rhs;}

P_GRID& nodal_fftwSolver::get_soln() {return *phi;}
