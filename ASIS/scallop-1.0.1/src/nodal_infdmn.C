#include <time.h>
#include "nodal_fftwSolver.h"
#include "nodal_infdmn.h"

#define FORT_NIDBOUND F77_FUNC(nidbound, NIDBOUND)
#define FORT_LAPLACIAN F77_FUNC(lap19, LAP19)
#define FORT_NODALPOTCALC F77_FUNC(npotcalc, NPOTCALC)
#define FORT_ZEROBC F77_FUNC_(n_boundary, N_BOUNDARY)
#define FORT_MAKE_RHS F77_FUNC(prhsmaker, PRHSMAKER)
#define FORT_EXACT F77_FUNC(exact, EXACT)

extern "C" {    
    void FORT_NIDBOUND(FAA_TYPE, FAA_TYPE, double*);
    void FORT_LAPLACIAN(FAA_TYPE, FAA_TYPE, double*);
  void FORT_NODALPOTCALC(FAA_TYPE, FAA_TYPE, double*, FAA_TYPE, FR_TYPE,
			   double*, double*, double*, double*, 
			   double*, double*, double*, double*, 
			   double*, double*, double*, double*,
			   double*, double*, double*, double*,
			   double*, double*,
			   int*, double*);
    void FORT_ZEROBC(FAA_TYPE);
    void FORT_MAKE_RHS(FAA_TYPE, FR_TYPE, double*, double*, double*, double*);
    void FORT_EXACT(FAA_TYPE, double*, double*, double*, double*);
};


void nodal_id_solve(P_GRID *rhs, P_GRID *soln, int nref, double h, 
		    int debug) {
    int debug_level = debug;
    int ii;

    P_GRID* pc_rhs;
    pc_rhs = new P_GRID((FLOORPLAN)(*rhs));
   
    pc_rhs->fill(0.);

    if (debug_level > 0) {
	cout << "starting inner grid solve..." << endl;
    }
    pc_rhs->copy(*rhs);
    if (debug_level > 5) {
	rhs->debug_print("pcrhs");
    }

    // Set up intermediate-sized domain for inner grid solve.  Its
    // boundaries are halfway between those of the rhs and those of
    // the solution.

    POINT low_h, high_h;
    POINT e1(1, 1, 1);
    POINT e2(2, 2, 2);
    P_DOMAIN d_i((FLOORPLAN)(*soln));
    for_1(i, *soln) {
	low_h = coarsen(((*soln)(i).lower() + (*rhs)(i).lower()), 2);
	high_h = coarsen(((*soln)(i).upper() + (*rhs)(i).upper()), 2);

	// We need to have low_h and high_h be divisible by 2 to do the 
	// Richardson-Romberg integration.  (Actually, it's simpler 
	// integration than that.)
	low_h = (coarsen(low_h, 2))*2;
	high_h = (coarsen(high_h - e1, 2) + e1)*2;
	// unfortunately comparison is not yet defined for a point
	//assert ((*soln)(i).lower() < low_h);
	//assert ((*soln)(i).upper() > high_h);
	d_i.setlower(i, low_h);
	d_i.setupper(i, high_h);
    } end_for;
    nodal_fftwSolver a_i(d_i, h);

    P_DOMAIN dF = (FLOORPLAN)(*soln);
    nodal_fftwSolver a(dF, h);

    for_all(i, *soln) {
	// Set the first moment charge at the boundary.
	FORT_NIDBOUND(FORT_ARRAY_ARG(a_i.get_soln(), i),
		      FORT_ARRAY_ARG_PTR(rhs, i),
		      &h);

	// Put in residual correction form.
	FORT_LAPLACIAN(FORT_ARRAY_ARG(a_i.get_soln(), i),
			FORT_ARRAY_ARG(a_i.get_rhs(), i),
			&h);
    } end_for_all;

    mpBarrier();

    if (debug_level > 5) {
	a_i.get_soln().debug_print("soln_after_first_solve");
    }
    a_i.get_soln().fill(0.);
    a_i.get_rhs().negate();
    a_i.get_rhs().plus(*pc_rhs);

    for_all(i, *soln) {

 // call FFTW solver for 19 point stencil.
 
    a_i.solve19(i);
	FORT_ZEROBC(FORT_ARRAY_ARG(a_i.get_soln(), i));
	FORT_NIDBOUND(FORT_ARRAY_ARG(a_i.get_soln(), i),
		      FORT_ARRAY_ARG_PTR(rhs, i),
		      &h);

    } end_for_all;

    if (debug_level > 5) {
	a_i.get_soln().debug_print("post_init");
	rhs->debug_print("post_init_rhs");
    }
    // We're going to waste some memory here to do the Richardson-
    // Romberg integration.  We may want to come up with a better way 
    // to do this later.
    P_GRID *extra;
    extra = new P_GRID((FLOORPLAN)(*soln));
    
    if (debug_level > 0) {
	cout << "starting potential calculation..." << endl;
    }

    /** These are all the arrays necessary for interpolating in the 
	far field.

	l?h and h?h  define planes corresponding to the faces 
	             of the inner grid
	l? and h?    define planes corresponding to the faces of
	             the outer grid
	l?c and h?c  define planes of coarse data from which the
	             outer grid values are interpolated
	the double * arrays hold the actual data for the planes
    */

    int lih, ljh, lkh, hih, hjh, hkh;
    int li, lj, lk, hi, hj, hk;
    int lic, ljc, lkc, hic, hjc, hkc;
    double *lefth, *righth, *downh, *uph, *backh, *fronth;
    double *left, *right, *down, *up, *back, *front;
    double *leftc, *rightc, *downc, *upc, *backc, *frontc;

    for_all(i, *soln) {
      lih = (d_i)(i).lower(0);
      ljh = (d_i)(i).lower(1);
      lkh = (d_i)(i).lower(2);
      hih = (d_i)(i).upper(0);
      hjh = (d_i)(i).upper(1);
      hkh = (d_i)(i).upper(2);
      li = (*soln)(i).lower(0);
      lj = (*soln)(i).lower(1);
      lk = (*soln)(i).lower(2);
      hi = (*soln)(i).upper(0);
      hj = (*soln)(i).upper(1);
      hk = (*soln)(i).upper(2);
      lic = (*soln)(i).lower(0)/nref - 3;
      ljc = (*soln)(i).lower(1)/nref - 3;
      lkc = (*soln)(i).lower(2)/nref - 3;
      hic = (*soln)(i).upper(0)/nref + 3;
      hjc = (*soln)(i).upper(1)/nref + 3;
      hkc = (*soln)(i).upper(2)/nref + 3;

      lefth = new double[(hjh - ljh + 1)*(hkh - lkh + 1)];
      righth = new double[(hjh - ljh + 1)*(hkh - lkh + 1)];
      downh = new double[(hih - lih + 1)*(hkh - lkh + 1)];
      uph = new double[(hih - lih + 1)*(hkh - lkh + 1)];
      backh = new double[(hih - lih + 1)*(hjh - ljh + 1)];
      fronth = new double[(hih - lih + 1)*(hjh - ljh + 1)];
      left = new double[(hj - lj + 1)*(hk - lk + 1)];
      right = new double[(hj - lj + 1)*(hk - lk + 1)];
      down = new double[(hi - li + 1)*(hk - lk + 1)];
      up = new double[(hi - li + 1)*(hk - lk + 1)];
      back = new double[(hi - li + 1)*(hj - lj + 1)];
      front = new double[(hi - li + 1)*(hj - lj + 1)];
      leftc = new double[(hjc - ljc + 1)*(hkc - lkc + 1)];
      rightc = new double[(hjc - ljc + 1)*(hkc - lkc + 1)];
      downc = new double[(hic - lic + 1)*(hkc - lkc + 1)];
      upc = new double[(hic - lic + 1)*(hkc - lkc + 1)];
      backc = new double[(hic - lic + 1)*(hjc - ljc + 1)];
      frontc = new double[(hic - lic + 1)*(hjc - ljc + 1)];

      FORT_NODALPOTCALC(FORT_ARRAY_ARG(a_i.get_soln(), i),
			  FORT_ARRAY_ARG_PTR(soln, i),
			  (*extra)(i).data(),
			  FORT_ARRAY_ARG_PTR(rhs, i),
			&lic, &ljc, &lkc, &hic, &hjc, &hkc,
			  lefth, righth, downh, uph, backh, fronth,
			  left, right, down, up, back, front,
			  leftc, rightc, downc, upc, backc, frontc,
			  &nref, &h);

      delete lefth;
      delete righth;
      delete downh;
      delete uph;
      delete backh;
      delete fronth;
      delete left;
      delete right;
      delete down;
      delete up;
      delete back;
      delete front;
      delete leftc;
      delete rightc;
      delete downc;
      delete upc;
      delete backc;
      delete frontc;
    } end_for_all;

    if (debug_level > 0) {
	cout << "starting outer grid solve..." << endl;
    }
    if (debug_level > 5) {
	soln->debug_print("npc");
	a_i.get_soln().mult(1/(h*h));
	a_i.get_soln().debug_print("charge");
    }
    for_all(i, *soln) {
  	FORT_NIDBOUND(FORT_ARRAY_ARG_PTR(soln, i),
  		      FORT_ARRAY_ARG_PTR(rhs, i),
  		      &h);		  

	FORT_LAPLACIAN(FORT_ARRAY_ARG_PTR(soln, i),
			FORT_ARRAY_ARG(a.get_rhs(), i),
			&h);
    } end_for_all;

    a.get_soln().fill(0.);
    a.get_rhs().negate();
    a.get_rhs().plus(*pc_rhs);

    
    for_all(i, *soln) {
	    a.solve19(i);
    } end_for_all;

    // We need to copy the soln, but only the interior, since the boundary
    // is already correct.
    soln->copy(a.get_soln(), grow((P_DOMAIN)a.get_soln(), -1));

    delete pc_rhs;
    delete extra;
}
    
