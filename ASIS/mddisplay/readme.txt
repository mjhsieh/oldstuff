=======================
Building MD Display 3.0
=======================

If you find MD Display helpful in your research,
please cite:

    Callahan, T.J., Swanson, E., Lybrand, T.P.:
       MD Display: An interactive graphics program for visualization
       of molecular dynamics trajectories,
    Journal of Molecular Graphics, 1996, 14:39-41
    
If you modify MD Display, or wish to drop in your own
movie frame rendering code, you might find this trade press
article helpful:

   Moth, Chris:
      A 3D Stereo Movie Viewer with GLUT
      C/C++ Users Journal, Sept 2003 pg 16
=============================================

mddisplay and preproc have no special requirements beyond 
the GLUT library and a complete ANSI C implementation.

If you already have built GLUT example programs, then you
should be able to immediately select a makefile.yourplatform and
immediately tailor it to your environment.

In general, you'll type

    make -f makefile.yourplatform clean
    make -f makefile.yourplatform

For information on program options, be sure to see the .pdf manual at

  http://www.structbio.vanderbilt.edu/~cmoth/mddisplay

  
    
If you are new to glut then the platform specific makefiles I've included
may help you get going more quickly.  Specific instructions for preparing
glut and your compilation environment are in each makefile.

The platforms/compilers I have tested are:

makefile.bc:   Borland compiler v 5.x for Windows, including free ver 5.5

makefile.msc:  Microsoft Visual C++ Version 6.0

makefile.lnx:  gcc for Linux, all processors

makefile.sgi:  cc for SGI IRIX ver 6.x.  See notes on legacy stereo.

makefile.osx:  Mac OS/X - note version specifics and issues in the makefile

makefile.osf:  Compaq TruUnix 64 (alpha chip)

I hope you enjoy working with MD Display 3.0

Please send suggestions or bug reports to me at:

chris.moth@vanderbilt.edu
http://www.structbio.vanderbilt.edu/~cmoth

Known bugs:

- There are a couple of issues with the GLUT implementation on Apple OS/X.
  These are documented in makefile.osx. 

- The movieTimerKey strategy documented in the CUJ article has an
  oversight.  If you drag another window over the movie, the environment
  generates additional glutPostRedisplay events.  This causes the
  doubling of movie speed that the algorithm is attempting to avoid.
  The solution will involve moving the call to MovieTimerFunc() out
  of MovieDisplayFunc().  I only discovered this a few days before
  going to press - after two years of running this code in our
  lab.  Apologies.
