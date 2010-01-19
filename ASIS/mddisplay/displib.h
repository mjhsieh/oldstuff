/*
 *  displib.h  core #defines and support library routines module for MD Display
 *     See also displib.c
 *
 *  These functions layer error checking on some ANSI C I/O routines
 *  and provide a platform independent getopt() implementation for
 *  command line parsing.
 *
 *  Copyright (C) 2002 by Chris Moth and Terry Lybrand
 *
 *  Do not distribute this code without explicit permission.
 *  Do not incorporate into other codes without explicit permission.
 *
 *  Current contact information:
 *
 *    Chris Moth
 *    chris.moth@vanderbilt.edu
 *    Phone: 615-936-3569
 *    http://www.structbio.vanderbilt.edu/~cmoth
 *--------------------------------------------------------------------------
 *    Derived from MD Display Version 2.0
 *
 *  Copyright (c) 1990  Timothy J. Callahan, Eric Swanson, and Terry Lybrand
 *  
 *  Please cite:       
 *  
 *    Callahan, T.J., Swanson, E., Lybrand, T.P.: 
 *       MD Display: AN interactive graphics program for visualization
 *       of molecular dynamics trajectories, 
 *    Journal of Molecular Graphics, 1996, 14:39-41
 *
*/

// scaling factor -- resolution is 1/KAPPA A 
// In order to save space in teh binary movie file, all coordinates are
// represented as 16 bit shorts.  Thus, all coordinates must be in 
// original range of +/- 32767/KAPPA
#define KAPPA 200.0       

// I am slowly but surely eliminating dependence on global
// fixed MAX values - and moving to dynamic dimensioning of all
// arrays

#define MAXATM 120000      /* Max # of Atoms */

#define MAXATP 1000		  /* Max # of Atom Types */

#define MAXDRAW (MAXATM + MAXATM/2)     /* for bonds (~1.5*MAXATM) */
#define MAXCHK 3.0
#define MAXBPA 6			  /* Max bonds per atom */

#define MAXSRF 200000	  /* Max Surfaces? */

#define MAXFRAMES 2500

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef FALSE
#undef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#undef TRUE
#define TRUE 1
#endif


/* This type is the "Attribute File Record" type created */

/* by preproc and read by display */

typedef struct {
	char label[5];   /* Atom "name" */

	short moln;      /*MAS*/
	unsigned short res;
	short color;
	short showlabl;
} atomattrib;

extern int debugFlag; /* Turn on with -d command line option */

// Open a file, ala fopen.  However, in case of error, give a robust
// error messag and exit to the OS.
FILE* fopen_or_exit_if_error(const char* FileName,const char* openMode);

// fclose a FILE* - but exit if there is a problem.  This is great
// for double-checking runtime integrity of file handles.
void fclose_or_exit_if_error(const char* FileName,FILE* stream);

// fread and fwrite routines that exit if unable to read/write all bytes.
void fread_or_exit_if_error(void *ptr, size_t size, size_t n, FILE *stream);
void fwrite_or_exit_if_error(const void *ptr, size_t size, size_t n, FILE *stream);

// malloc which exits and prints error in case of failure 
void* malloc_or_exit_if_error(size_t size,const char* context);

// platform-independent "getopt" implementation.
int getopt(int argc,char* const* argv,const char* optstr);

