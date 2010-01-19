#include "poisson.h"
#include "timer.h"


#define FORT_STENCIL_INTERP F77_FUNC(interpstencil, INTERPSTENCIL)
#define FORT_STENCIL_AVERAGE F77_FUNC(avestencil, AVESTENCIL)
#define FORT_SAMPLE F77_FUNC(sample, SAMPLE)
#define FORT_LAPLACIAN F77_FUNC(lap19, LAP19)
#define FORT_LAPLACIAN7 F77_FUNC(lap7, LAP7)

extern "C" {
    void FORT_STENCIL_INTERP(FAA_TYPE, FAA_TYPE, const int*, const
			     int*, double*);
    void FORT_STENCIL_AVERAGE(FAA_TYPE, const int*);
    void FORT_SAMPLE(FAA_TYPE, FAA_TYPE, FR_TYPE, const int*);
    void FORT_LAPLACIAN(FAA_TYPE, FAA_TYPE, double*);
    void FORT_LAPLACIAN7(FAA_TYPE, FAA_TYPE, double*);
};

poisson::poisson() {
    input_rhs = NULL;
    fine_rhs = NULL;
    fine_soln = NULL;
    coarse_rhs = NULL;
    coarse_domain_rhs = NULL;
    coarse_domain_soln = NULL;
    final_soln = NULL;
}

poisson::poisson(const REGION& input_domain, const P_DOMAIN& input_patches,
		 double input_h, int input_nref, int debug) {

    // Set up some defaults first.  

    // In my experience, I have never seen a need for a correction
    // radius greater than 0.  I'm leaving it in the code in case it
    // needs to be turned on later, at which point I will make another
    // constructor and add a method to set the value.
    // Also, in my experience, id_border of 10 works just fine.  I'll
    // add hooks later to change it, if necessary.  
    // Both these changes could be a bit painful: you can't change
    // them after instantiation without completely changing the
    // grids.  Obviously, that's not something you want to do.  Hm. 
    c = 0;
    int id_border = 10;

    // The value of d is stencil-dependent, and shouldn't be changed
    // unless the stencil is changed.  There shouldn't be a user
    // option to change this.  Change c instead: it will have the same
    // effect.
    int d = 2;

    fine_domain = input_domain;
    fine_interior = input_patches;
    fine_h = input_h;
    nref = input_nref;  // Consider choosing a default for this somehow?
    debug_level = debug;

    coarse_h = input_h*nref;

    int border = MAX((c+d)*nref, id_border*(MAX(MAX(fine_domain.extents(0),
						    fine_domain.extents(1)),
						fine_domain.extents(2)))
		     / (input_patches.size()*100));

    int lcm_ref = ((nref - 1)/4 + 1)*4;
    P_DOMAIN all_fine_domains =
	nrefine(ncoarsen(grow(fine_interior, border), lcm_ref), lcm_ref);
    REGION phi_region;
    POINT low_h, high_h;
    POINT e1(1, 1, 1);
    POINT e2(2, 2, 2);

    for_1(i, input_patches) {
	phi_region = all_fine_domains(i);
	low_h = coarsen((phi_region.lower() + input_patches(i).lower()), 2);
	high_h = (phi_region.upper() + input_patches(i).upper())/e2;
	low_h = (coarsen(low_h, 2))*2;
	high_h = (coarsen(high_h - e1, 2) + e1)*2;
	if (!(grow(phi_region, -1).inside(low_h)
	      && grow(phi_region, -1).inside(high_h))) {
	    all_fine_domains.grow(i, lcm_ref);
	}
    } end_for;

    coarse_interior = ncoarsen(input_patches, nref);
    coarse_domain = grow(ncoarsen(input_patches, nref), c+d);

    coarse_rhs = new P_GRID(grow(coarse_interior, c));
    coarse_soln = new P_GRID(coarse_domain);
    fine_rhs = new P_GRID(fine_interior);
    fine_soln = new P_GRID(all_fine_domains);
    // By the way, input_rhs is not set here, but must be set by the
    // set_rhs method.

    REGION coarse_rhs_region;
    REGION coarse_soln_region;
    coarse_rhs_region = grow(coarsen(fine_domain, nref), c);
    border = id_border*(MAX(MAX(coarse_rhs_region.extents(0),
				coarse_rhs_region.extents(1)),
			    coarse_rhs_region.extents(2)))/100;
    coarse_soln_region = grow(coarse_rhs_region, border);
    coarse_soln_region = nrefine(ncoarsen(coarse_soln_region, 4), 4);

    low_h = coarsen((coarse_soln_region.lower() + coarse_rhs_region.lower()), 2);
    high_h = (coarse_soln_region.upper() + coarse_rhs_region.upper())/e2;
    low_h = (coarsen(low_h, 2))*2;
    high_h = (coarsen(high_h - e1, 2) + e1)*2;
    // This, i.e. the hardcoded 4 below, may need to be adjusted.  It makes
    // sense for the kinds of problems I'm working on now, though.
    if (!(grow(coarse_soln_region, -1).inside(low_h)
	  && grow(coarse_soln_region, -1).inside(high_h))) {
	coarse_soln_region.grow(4);
    }

    P_DOMAIN coarse_rhs_boxes(mpNodes());
    P_DOMAIN coarse_soln_boxes(mpNodes());
    // Assigning one box per processor.
    for_1(i, coarse_rhs_boxes) {
	coarse_rhs_boxes.setregion(i, coarse_rhs_region);
	coarse_soln_boxes.setregion(i, coarse_soln_region);
	coarse_rhs_boxes.setowner(i, i);
	coarse_soln_boxes.setowner(i, i);
    } end_for;

    coarse_domain_rhs = new P_GRID(coarse_rhs_boxes);
    coarse_domain_soln = new P_GRID(coarse_soln_boxes);

    final_soln = new nodal_fftwSolver(fine_interior, fine_h);
}

// Set the right hand side.  I'm using a user supplied P_GRID, which
// should match the P_DOMAIN that instance of the solver was created
// with.  Eventually, I should make it easy for users to use a more
// general array structure (or one of their own), rather than
// expecting them to use the kAPI.
int poisson::set_rhs(P_GRID& rhs) {
    // TODO
    // Do some sanity checking here: does the given rhs match the
    // domains I'm expecting?  If not, return an error code.
    input_rhs = &rhs;
    fine_rhs->fill(0.);
    fine_rhs->copy(*input_rhs);
    return 1;
}

void poisson::fine_boundary_communication() {
    mpBarrier();
    stencil_setup();
    mpBarrier();


    // Make coarse-grain remote copy of necessary local coarse data.

    // This is terribly ugly.
    // For each grid of coarse_domain that I own, make grids corresponding
    // to the intersection of that coarse_domain grid with every other 
    // coarse_domain grid.  Nice idea, but _everyone_ needs a copy of the 
    // whole P_DOMAIN.

    // Count number of grids required.
    REGION intersect;
    int intersect_count = 0;
    for (int current_processor = 0; current_processor < mpNodes();
	 current_processor++) {
      for_1 (i, coarse_domain) {
	if (coarse_domain.owner(i) == current_processor) {
	  for_1 (j, coarse_domain) {
	    intersect = coarse_domain(i) * coarse_domain(j);
	    if (!intersect.empty()) {
	      intersect_count += 1;
	    }
	  } end_for;
	}
      } end_for;
    }


    lc_coarse_domain = new P_DOMAIN(intersect_count);

    int current = 0;
    lc_remote_reference = new int*[coarse_domain.size()];
    for (int counter = 0; counter < coarse_domain.size(); counter++) {
      lc_remote_reference[counter] = new int[coarse_domain.size()];
    }

    // Build floorplan of local coarse domains.
    for (int current_processor = 0; current_processor < mpNodes();
	 current_processor++) {
      for_1 (i, coarse_domain) {
	if (coarse_domain.owner(i) == current_processor) {
	  for_1 (j, coarse_domain) {
	    intersect = coarse_domain(i) * coarse_domain(j);
	    if (!intersect.empty()) {
	      lc_coarse_domain->setregion(current, intersect);
	      lc_coarse_domain->setowner(current, current_processor);
	      // if processor i wants to get a local copy of the data
	      // from processor j, it can look it up through this array
	      lc_remote_reference[i][j] = current;
	      current += current;
	    }
	  } end_for;
	}
      } end_for;
    }

    lc_coarse_soln = new P_GRID(*lc_coarse_domain);

    // Build motion plan to fill in lc_coarse_soln.
    MOTIONPLAN LC;
    current = 0;
    for (int current_processor = 0; current_processor < mpNodes();
	 current_processor++) {
      for_1 (i, coarse_domain) {
	if (coarse_domain.owner(i) == current_processor) {
	  for_1 (j, coarse_domain) {
	    intersect = coarse_domain(i) * coarse_domain(j);
	    if (!intersect.empty()) {
	      LC.CopyOnIntersection(*coarse_soln, j, *lc_coarse_soln, current);
	    }
	  } end_for;
	}
      } end_for;
    }

    MOVER local_coarse_copy(*coarse_soln, *lc_coarse_soln, LC);
    local_coarse_copy.execute();

    /// what the heck???
    ////#if 0
    

    // Same kind of thing, but now for fine grid.
    intersect_count = 0;
    for (int current_processor = 0; current_processor < mpNodes();
	 current_processor++) {
      for_1 (i, fine_interior) {
	if (fine_interior.owner(i) == current_processor) {
	  for_1 (j, fine_interior) {
	    intersect = fine_interior(i) * fine_interior(j);
	    if (!intersect.empty()) {
	      intersect_count += 1;
	    }
	  } end_for;
	}
      } end_for;
    }

    lf_fine_domain = new P_DOMAIN(intersect_count);

    current = 0;
    lf_remote_reference = new int*[fine_interior.size()];
    for (int counter = 0; counter < fine_interior.size(); counter++) {
      lf_remote_reference[counter] = new int[fine_interior.size()];
    }
    // Build floorplan of local coarse domains.
    for (int current_processor = 0; current_processor < mpNodes();
	 current_processor++) {
      for_1 (i, fine_interior) {
	if (fine_interior.owner(i) == current_processor) {
	  for_1 (j, fine_interior) {
	    intersect = fine_interior(i) * fine_interior(j);
	    if (!intersect.empty()) {
	      lf_fine_domain->setregion(current, intersect);
	      lf_fine_domain->setowner(current, current_processor);
	      // if processor i wants to get a local copy of the data
	      // from processor j, it can look it up through this array
	      lf_remote_reference[i][j] = current;
	      current += 1;
	    }
	  } end_for;
	}
      } end_for;
    }
    
    lf_fine_soln = new P_GRID(*lf_fine_domain);

    // Build motion plan to fill in lc_coarse_soln.
    MOTIONPLAN LF;
    current = 0;
    for (int current_processor = 0; current_processor < mpNodes();
	 current_processor++) {
      for_1 (i, fine_interior) {
	if (fine_interior.owner(i) == current_processor) {
	  for_1 (j, fine_interior) {
	    intersect = fine_interior(i) * fine_interior(j);
	    if (!intersect.empty()) {
	      LF.CopyOnIntersection(*fine_soln, j, *lf_fine_soln, current);
	    }
	  } end_for;
	}
      } end_for;
    }

    MOVER local_fine_copy(*fine_soln, *lf_fine_soln, LF);
    local_fine_copy.execute();
    /// end of what the heck???
    ////#endif
}
// Set the boundary for all the fine patches, before the final solve.
void poisson::set_fine_boundary() {

    // Copy from local coarse-grid solution onto local stencils.
    MOTIONPLAN M;
    for_1 (i, *s->coarse_data) {
	M.CopyOnIntersection(*coarse_domain_soln,
			     s->coarse_data->owner(i),
			     *s->coarse_data, i);
    } end_for;
    //    mpBarrier();
    MOVER move(*coarse_domain_soln, *s->coarse_data, M);
    move.execute();
    //    mpBarrier();
    //    OUTPUT("done mover..." << endl);

    // Subtract off local coarse grid data (using local copies).
    MOTIONPLAN S;
    for_1 (i, *s->coarse_domain) {
	for_1(j, coarse_domain) {
	    if (coarse_domain(j).inside((*s->coarse_domain)(i).lower())
		&& coarse_domain(j).inside((*s->coarse_domain)(i).upper())) {
		S.CopyOnIntersection(*lc_coarse_soln, 
				     lc_remote_reference[stencil_patch_index[i]][j], 
				     *s->coarse_data, i);
	    }
	} end_for;
    } end_for;
    // mpBarrier();
    SUBTRACTOR sub(*lc_coarse_soln, *s->coarse_data, S);
    sub.execute();
    //mpBarrier();
    //OUTPUT("done subtractor..." << endl);

    // Interpolate
    for_all(i, *s->coarse_data) {
	FORT_STENCIL_INTERP(FORT_ARRAY_ARG_PTR(s->fine_data, i),
			    FORT_ARRAY_ARG_PTR(s->coarse_data, i),
			    &stencil_dir_array[i], &nref, &coarse_h);
    } end_for_all;
    //mpBarrier();
    //OUTPUT("done interpolate..." << endl);

    // Make coarse-grain remote copy of necessary local fine data.
//     MOTIONPLAN LF;
//     for_1 (i, lf_fine_domain) {
//       for_1 (j, coarse_domain) {
// 	LF.CopyOnIntersection(*fine_soln, j, lf_fine_domain, i);
//       } end_for;
//     } end_for;

    // Add back local fine grid data (using local copies).
    MOTIONPLAN A;
    for_1(i, *s->coarse_domain) {
	for_1(j, coarse_domain) {
	    if (coarse_domain(j).inside((*s->coarse_domain)(i).lower())
		&& coarse_domain(j).inside((*s->coarse_domain)(i).upper())) {
		A.CopyOnIntersection(*lf_fine_soln, 
				     lf_remote_reference[stencil_patch_index[i]][j], 
				     *s->fine_data, i);
	    }
	} end_for;
    } end_for;
    //mpBarrier();

    //OUTPUT("before adder..." << endl);
    ADDER add(*lf_fine_soln, *s->fine_data, A);
    add.execute();
    //OUTPUT("done adder..." << endl);

    for_all(i, *s->fine_data) {
	FORT_STENCIL_AVERAGE(FORT_ARRAY_ARG_PTR(s->fine_data, i),
			     &stencil_dir_array[i]);
    } end_for_all;

    // Copy results to boundary.
    MOTIONPLAN C;
    for_1(i, *s->fine_data) {
	C.CopyOnIntersection(*s->fine_data, i,
			     final_soln->get_soln(), stencil_patch_index[i]);
    } end_for;
    ADDER really_copy(*s->fine_data, final_soln->get_soln(), C);
    really_copy.execute();
    //OUTPUT("done really_copy..." << endl);
}

void poisson::stencil_setup() {
    int stencil_number = 0;
    int ii, jj, kk;
    POINT start_point;
    POINT end_point;
    POINT tmp_point;
    POINT ex(1, 0, 0);
    POINT ey(0, 1, 0);
    POINT ez(0, 0, 1);
    POINT e1(1, 1, 1);  

    // I need a stencil for every point on the face of a
    // coarse_interior.
    for_1(i, *fine_soln) {
	stencil_number += 2*(coarse_interior(i).extents(0)
			     *coarse_interior(i).extents(1)
			     + coarse_interior(i).extents(0)
			     *coarse_interior(i).extents(2)
			     + coarse_interior(i).extents(1)
			     *coarse_interior(i).extents(2));
    } end_for;

    sf_domain = new P_DOMAIN(stencil_number);
    sc_domain = new P_DOMAIN(stencil_number);
    stencil_patch_index = new int[stencil_number];
    stencil_dir_array = new int[stencil_number];

    // I'm pretty sure there is an iterator in KeLP1.4 that would make
    // my life easier here, but such conveniences don't seem to exist
    // in 1.3, and I'm not interested in dealing with the 1.3->1.4
    // switchover yet.
    stencil_number = 0;
    for_1(i, *fine_soln) {
	for (int direction = 0; direction < 6; direction++) {
	    // The start and end points are the lower and upper
	    // extents of the coarse interior, shifted appropriately
	    // to get only the chosen faces.
	    start_point = coarse_interior(i).lower();
	    end_point = coarse_interior(i).upper();
	    // set the left face
	    switch (direction) {
	    case LEFT:
		end_point(0) = coarse_interior(i).lower(0);
		break;
	    case RIGHT:
		start_point(0) = coarse_interior(i).upper(0);
		break;
	    case FRONT:
		end_point(1) = coarse_interior(i).lower(1);
		break;
	    case BACK:
		start_point(1) = coarse_interior(i).upper(1);
		break;
	    case BOTTOM:
		end_point(2) = coarse_interior(i).lower(2);
		break;
	    case TOP:
		start_point(2) = coarse_interior(i).upper(2);
		break;
	    default:
		// perhaps add some catch here, but I'm fairly
		// confident the default never happens
		break;
	    }

	    // Again, the comment about this being easier with an
	    // iterator.
	    switch(direction) {
	    case LEFT:
	    case RIGHT:
		tmp_point(0) = start_point(0);
		for (jj = start_point(1); jj <= end_point(1); jj++) {
		    tmp_point(1) = jj;
		    for (kk = start_point(2); kk <= end_point(2);
			 kk++) {
			tmp_point(2) = kk;
			sc_domain->setregion(stencil_number,
					     grow(REGION(tmp_point,
							 tmp_point),
						  2));
			sc_domain->setowner(stencil_number,
					    fine_soln->owner(i));
			sf_domain->setregion(stencil_number,
					     REGION(tmp_point*nref 
						    - (e1 - ex)*nref/2, 
						    tmp_point*nref 
						    + (e1 - ex)*nref/2));

			sf_domain->setowner(stencil_number,
					    fine_soln->owner(i));
			stencil_patch_index[stencil_number] = i;
			stencil_dir_array[stencil_number] = direction;

			stencil_number++;
		    }
		}
		break;
	    case FRONT:
	    case BACK:
		tmp_point(1) = start_point(1);
		for (ii = start_point(0); ii <= end_point(0); ii++) {
		    tmp_point(0) = ii;
		    for (kk = start_point(2); kk <= end_point(2);
			 kk++) {
			tmp_point(2) = kk;
			sc_domain->setregion(stencil_number,
					     grow(REGION(tmp_point,
							 tmp_point),
						  2));
			sc_domain->setowner(stencil_number,
					    fine_soln->owner(i));
			sf_domain->setregion(stencil_number,
					     REGION(tmp_point*nref 
						    - (e1 - ey)*nref/2, 
						    tmp_point*nref 
						    + (e1 - ey)*nref/2));

			sf_domain->setowner(stencil_number,
					    fine_soln->owner(i));
			stencil_patch_index[stencil_number] = i;
			stencil_dir_array[stencil_number] = direction;

			stencil_number++;
		    }
		}
		break;
	    case TOP:
	    case BOTTOM:
		tmp_point(2) = start_point(2);
		for (ii = start_point(0); ii <= end_point(0); ii++) {
		    tmp_point(0) = ii;
		    for (jj = start_point(1); jj <= end_point(1);
			 jj++) {
			tmp_point(1) = jj;
			sc_domain->setregion(stencil_number,
					     grow(REGION(tmp_point,
							 tmp_point),
						  2));
			sc_domain->setowner(stencil_number,
					    fine_soln->owner(i));
			sf_domain->setregion(stencil_number,
					     REGION(tmp_point*nref 
						    - (e1 - ez)*nref/2, 
						    tmp_point*nref 
						    + (e1 - ez)*nref/2));

			sf_domain->setowner(stencil_number,
					    fine_soln->owner(i));
			stencil_patch_index[stencil_number] = i;
			stencil_dir_array[stencil_number] = direction;

			stencil_number++;
		    }
		}
		break;
	    default:
		// perhaps add some catch here, but I'm fairly
		// confident the default never happens
		break;
	    }
	}
    } end_for;

    s = new stencil(sc_domain, sf_domain);
}

void poisson::solve() {
    Statistics StatMaster(StatsNames, NUM_STATS);

    OUTPUT("First loop over patches..." << endl);

    STATS_START(STATS_TOTAL);
    STATS_START(STATS_FIRSTPASS);

    // Do infinite domain solve on each patch.
    nodal_id_solve(fine_rhs, fine_soln, nref, fine_h, debug_level);

    // Average onto coarse grids.
    for_all(i, *fine_soln) {
	FortranRegion3 sample_region = 
	    FortranRegion3((*coarse_soln)(i).region()
			   *grow(ncoarsen(grow((*fine_soln)(i).region(), 1), 
					  nref), -1));
	FORT_SAMPLE(FORT_ARRAY_ARG_PTR(fine_soln, i),
		    FORT_ARRAY_ARG_PTR(coarse_soln, i),
		    FORT_REGION(sample_region),
		    &nref);
    } end_for_all;

    // Take the laplacian on the coarse grids.
    for_all(i, *coarse_soln) {
  	FORT_LAPLACIAN(FORT_ARRAY_ARG_PTR(coarse_soln, i),
		       FORT_ARRAY_ARG_PTR(coarse_rhs, i),
		       &coarse_h);
    } end_for_all;
    
    mpBarrier();
    STATS_STOP(STATS_FIRSTPASS);
    STATS_START(STATS_REDUCTION);

    // Gather coarse data onto one large grid, coarse_domain.  This
    // means adding up ALL contributions, even if they overlap.
    // This becomes our new rhs.
    coarse_domain_rhs->plus_all(*coarse_rhs);

    STATS_STOP(STATS_REDUCTION);

    // mpBarrier();

    STATS_START(STATS_COARSEGRID);

    OUTPUT("Doing infinite domain solve on the coarse grid..." << endl);
    OUTPUT((*coarse_domain_soln)(0).region() << endl);

    // Do infinite domain solve on coarse grid.
    nodal_id_solve(coarse_domain_rhs, coarse_domain_soln, nref, coarse_h, 
		   debug_level);



    // Again looping over all patches...
    OUTPUT("Second loop over patches..." << endl);

    STATS_STOP(STATS_COARSEGRID);
    mpBarrier();

    // Interpolate from coarse grid to fix fine grid boundary
    // conditions.
    STATS_START(STATS_SETBOUNDARY);
    fine_boundary_communication();
    STATS_STOP(STATS_SETBOUNDARY);
    OUTPUT("Done communicating boundary" << endl);
    STATS_START(STATS_FINALPASS);
    set_fine_boundary();

    OUTPUT("Done setting boundary" << endl);

    mpBarrier();

    P_GRID dummy_soln(final_soln->get_soln());
    P_GRID dummy_rhs(final_soln->get_rhs());
    dummy_soln.fill(0.);
    dummy_rhs.fill(0.);


    dummy_soln.copy(final_soln->get_soln());
    dummy_rhs.copy(final_soln->get_rhs());	

	
    for_all(i, dummy_soln) {
    	FORT_LAPLACIAN7(FORT_ARRAY_ARG(dummy_soln, i),
			 FORT_ARRAY_ARG(final_soln->get_rhs(), i),
			 &fine_h);
    } end_for_all;
    final_soln->get_rhs().negate();
    final_soln->get_rhs().plus(dummy_rhs);
	
    final_soln->get_soln().fill(0.);

    // Poisson solve on each fine patch again, using the boundary
    // conditions we just set.
    
    for_all(i, final_soln->get_soln()) {

// Call FFTW solver for final 7-point solution to Poisson's equation. 
      
        final_soln->solve7(i);

    } end_for_all;

    mpBarrier();

    STATS_STOP(STATS_FINALPASS);
    STATS_STOP(STATS_TOTAL);

    mpBarrier();

    KELP_STATS_REPORT(STATS_AVG);
    STATS_REPORT(STATS_AVG);
}
