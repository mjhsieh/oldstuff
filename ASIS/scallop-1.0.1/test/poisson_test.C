#include "poisson_test.h"
#define FORT_PRHSMAKER F77_FUNC(prhsmaker, PRHSMAKER)
extern "C" {
    void FORT_PRHSMAKER(FAA_TYPE, FR_TYPE, double*, double*, double*, double*);
};

struct inputs {
    int c;
    int nref;
    int id_border;
    int long_dim;
    int num_patches;
    int input;
    double m;
    int debug_level;
};

/****************************************************************
*	int ParseCommandLine (int argc,char **argv)		*
*	returns N, the first command line argument		*
****************************************************************/

inputs ParseCommandLine (int argc, char **argv) //, 
{
    /**************************************************/
    /* Parse input values - make sure the value is OK */
    /**************************************************/
    double m = 0.0;
    int c = 0;
    int nref = 4;
    int id_border = 10;
    int long_dim = 64;
    int num_patches = 1;
    int input = CENTERED;
    int debug_level = 0;
	
    int i;
    int argcount = argc;

#ifdef P4
    argcount -= 2;       // P4 passes extra arguments
#endif

    // Read command line arguments to set values. 
    if (mpMyID() == 0) {
	for (i = 1; (i < argcount && *argv[i] == '-'); i++) {
	    if (!strncmp("-c", argv[i], 2)) {
		c = atoi(argv[++i]);
	    }
	    else if (!strncmp("-nref", argv[i], 3)) {
		nref = atoi(argv[++i]);
	    }
	    else if (!strncmp("-idborder", argv[i], 3)) {
		id_border = atoi(argv[++i]);
	    }
	    else if (!strncmp("-longdim", argv[i], 2)) {
		long_dim = atoi(argv[++i]);
	    }
	    else if (!strncmp("-npatches", argv[i], 3)) {
		num_patches = atoi(argv[++i]);
	    }
	    else if (!strncmp("-wave", argv[i], 2)) {
		m = double(atoi(argv[++i]));
	    }
	    else if (!strncmp("-debug", argv[i], 2)) {
		debug_level = atoi(argv[++i]);
	    }
	    else if (!strncmp("-input", argv[i], 3)) {
		if (!strncmp("file", argv[i+1], 1)) {
		    input = FROMFILE;
		}
		else if (!strncmp("centered", argv[i+1], 1)) {
		    input = CENTERED;
		}
		else if (!strncmp("periodic", argv[i+1], 1)) {
		    input = PERIODIC;
		}
		else if (!strncmp("mperiodic", argv[i+1], 1)) {
		    input = MPERIODIC;
		}
		else if (!strncmp("random", argv[i+1], 1)) {
		    input = RANDOMISH;
		}
		else if (!strncmp("one", argv[i+1], 2)) {
		    input = ONE;
		}
		else if (!strncmp("offset", argv[i+1], 2)) {
		    input = OFFSET;
		}
		i++;
	    }
	}

	if ((num_patches < 1) 
	    || (mpNodes() >  num_patches*num_patches*num_patches))  {
	    OUTPUT(argv[0] << ": Invalid Npatches\nprocessors");
	    exit (ERROR);		
	}
    }

//#ifdef P4

    // MPI problem: only root process gets argv correctly
    // must broadcast ring size to other processes

     MPI_Bcast(&c,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&nref,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&id_border,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&long_dim,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&num_patches,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&input,1,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(&m,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     MPI_Bcast(&debug_level,1,MPI_INT,0,MPI_COMM_WORLD);
    
//#endif

     inputs return_values;
     return_values.c = c;
     return_values.nref = nref;
     return_values.id_border = id_border;
     return_values.long_dim = long_dim;
     return_values.num_patches = num_patches;
     return_values.input = input;
     return_values.m = m;
     return_values.debug_level = debug_level;

     return return_values;
}

/****************************************************************
* main()							*
*								*
* main() takes one argument: N: the number of points		*
* along each axis						*
*****************************************************************/
int main(int argc,char **argv)
{
    double h, m;
    int i, j, k;
    int c, nref, id_border, long_dim, num_patches, patch_length, debug_level;
    int input;

    MPI_Init(&argc, &argv);

    InitKeLP(argc,argv);

// I ought to test which of these is faster sometime...
    KelpConfig(CONTIG_MSG_IN_PLACE,TRUE);
//    KelpConfig(CONTIG_MSG_IN_PLACE,FALSE);

    inputs setup_values;

    setup_values = ParseCommandLine(argc, argv);

    c = setup_values.c;
    nref = setup_values.nref;
    id_border = setup_values.id_border;
    long_dim = setup_values.long_dim;
    num_patches = setup_values.num_patches;
    input = setup_values.input;
    m = setup_values.m;
    debug_level = setup_values.debug_level;

    /* Print header information*/
    OUTPUT("Multigrid computation run on " << mpNodes() << " processors with " 
	   << num_patches << "^3 patches" << endl);	


    REGION interior(0, 0, 0, long_dim, long_dim, long_dim);
    REGION patch;

    h = 1.0/long_dim;

    patch_length = long_dim/num_patches;

    P_DOMAIN patches(num_patches*num_patches*num_patches);

    int n_procs = mpNodes();
    for (i = 0; i < num_patches; i++) {
	for (j = 0; j < num_patches; j++) {
	    for (k = 0; k < num_patches; k++) {
		patch.setlower(i*patch_length, j*patch_length, k*patch_length);
		patch.setupper((i+1)*patch_length, (j+1)*patch_length,
			   (k+1)*patch_length);
		patches.setregion(i*num_patches*num_patches +
				  j*num_patches + k, patch);
		patches.setowner(i*num_patches*num_patches +
				 j*num_patches + k, 
				 (i*num_patches*num_patches +
				  j*num_patches + k)%n_procs);
	    }
	}
    }

    // Make a test rhs.
    double center[3], radius;
    radius = 0.425;
    center[0] = 0.475;
    center[1] = 0.475;
    center[2] = 0.475;

    P_GRID rhs(patches);
    for_all(i, rhs) {
	FORT_PRHSMAKER(FORT_ARRAY_ARG(rhs, i),
		       FORT_ARRAY_BOUND(rhs, i),
		       &radius, center, &h, &m);
    } end_for_all;

    // Call PMG constructor.
    poisson my_problem(interior, patches, h, nref, debug_level);

    my_problem.set_rhs(rhs);

    // solve
    my_problem.solve();

    MPI_Finalize();
    return(0);
}
