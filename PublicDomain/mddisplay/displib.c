/*
 *  displib.c  core support library routines module for MD Display
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "displib.h" /* prototypes for these functions */


int debugFlag = 0; /* Turn on with -d command line option */

// Open a file, ala fopen.  However, in case of error, give a robust
// error messag and exit to the OS.

FILE* fopen_or_exit_if_error(const char* FileName,const char* openMode)
{
   FILE* f =fopen(FileName,openMode);

   if (! f)
      {
      fprintf(stderr,"Unable to fopen(\"%s\",\"%s\"): %s\n",FileName,openMode,strerror(errno));
      exit(EXIT_FAILURE);
      }

   return f;
}

// fclose a FILE* - but exit if there is a problem.  This is great
// for double-checking runtime integrity of file handles.

void fclose_or_exit_if_error(const char* FileName,FILE* stream)
{
   if (fclose(stream) != 0)
      {
      fprintf(stderr,"Unable to fclose(\"%s\"): %s\n",FileName,strerror(errno));
      exit(EXIT_FAILURE);
      }
}

// fread, but fail if the data is not there.

void fread_or_exit_if_error(void *ptr, size_t size, size_t n, FILE *stream)
{
   size_t nread = fread(ptr,size,n,stream);

   if ( nread != n )
      {
      fprintf(stderr,"Unable to fread(%d elements of size %d).  fread() returned %d  Error: %s\n",n,size,nread,strerror(errno));
      exit(EXIT_FAILURE);
      }
}

// fwrite, but fail if data can not be written

void fwrite_or_exit_if_error(const void *ptr, size_t size, size_t n, FILE *stream)
{
   size_t nwritten  = fwrite(ptr,size,n,stream);

   if ( nwritten != n )
      {
      fprintf(stderr,"Unable to fwrite(%d elements of size %d): %s\n",n,size,strerror(errno));
      exit(EXIT_FAILURE);
      }
}

// like malloc - but fails with a clear message in case of failure.

void* malloc_or_exit_if_error(size_t size,const char* context)
{
   void* retval = malloc(size);

   if ( retval == 0 )
      {
      fprintf(stderr,"Unable to malloc(%ld): Out of Dynamic Memory in %s\n",(long)size,context);
      exit(EXIT_FAILURE);
      }

   if (debugFlag)
      printf("Successful allocation of %ld bytes in %s\n",(long)size,context);

   return retval;
}


// The code below gives a "getopt" routine which parses out the command line.  Not all 
// C libraries have a "getopt" so this fills in the gap with similar, 
// though perhaps not as robust, semantics.  I vaguely recall finding this on a website
// somewhere - but I've modified it a bit at least.  Apologies to original author.

char* optarg;
int optind = 1;
int opterr,optopt;
void getoptreset(void)
{
   optind = 1;
}

static int getopt_missing_argument(const char* optstr)
{
   /* We come here if we were missing an argument */
   /* after a valid option letter. */


   /* The semantics are to return a colon ':' if the */
   /* first character of optstr is a : or '?' otherwise */
   optarg = "";

   return (*optstr == ':') ? ':' : '?';
}
int getopt(int argc,char* const* argv,const char* optstr)
{
   optarg = ""; /* Point optarg to an empty string for safety. */


   if (optind >= argc) /* No more command lien args - return -1 */

      return -1;

   if (argv[optind][0] == '-') /* Looks promising, we have a hyphen */

      {
      int testchar = argv[optind][1];
      const char* argstr_ptr = strchr(optstr,testchar);
      if (argstr_ptr == 0) /* letter after hyphen not acceptable? */

         return '?'; /* Option character not in optstring... */


      /* OK - we have a valid letter...  Does it need an argument? */

      if (argstr_ptr[1] == ':') /* Yes, it needs an argument... */

         {
         /* go ahead and advance the global - because */
         /* we _might_ return from this block... */

         optind++;

         /* If we are clearly missing an argument then return */

         if (optind >= argc)
            return getopt_missing_argument(optstr);

         /* Well - we seem to have something following the option */

         optarg = argv[optind];

         /* But - if that appears to be a command character option */
         /* We need to return "missing argument" as well. */

         if ((optarg[0] == '-') &&
             (strlen(optarg) == 2) &&
             (strchr(optstr,optarg[1]) != 0))
                  return getopt_missing_argument(optstr);
         }

      /* all is well */

      optind++;
      return testchar;
      }

   return getopt_missing_argument(optstr);
}
