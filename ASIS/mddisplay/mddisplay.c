/* ------------------------------------------------------------------------
 * MD Display Version 3.0
 * Copyright (C) 2002 Chris Moth and Terry Lybrand
 * Vanderbilt University
 *
 *    Derived from MD Display Version 2.0
 *
 *  Copyright (c) 1990  Timothy J. Callahan, Eric Swanson, and Terry Lybrand
 *
 *  Please cite:
 *
 *    Callahan, T.J., Swanson, E., Lybrand, T.P.:
 *       MD Display: An interactive graphics program for visualization
 *       of molecular dynamics trajectories,
 *    Journal of Molecular Graphics, 1996, 14:39-41
 *
 * Do not distribute this code without explicit permission.
 * Do not incorporate into other codes without explicit permission.
 *
 * MDDisplay.c manages rendering and coloring, as well as
 * all MD specific features.  Functions in this module are called
 * by dispui.c - the generic 3D Stereo Movie user interface module.  For more
 * information on dispui.c see
 *
 *    Moth, C., C++/C Users Journal, Sept 2003
 *    "A 3D Stereo Movie Viewer with GLUT"
 *
 *  Contact information as of August 1, 2003
 *
 *    Chris Moth
 *    chris.moth@vanderbilt.edu
 *    Phone: 615-936-3569
 *    http://www.structbio.vanderbilt.edu/~cmoth
 *--------------------------------------------------------------------------
 */

// Ansi C headers  
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <time.h>

#ifdef __APPLE__ // "Think different" for OS/X :) 
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifdef _MSC_VER
#include <direct.h>
#endif

#ifdef __TURBOC__
#include <dir.h>
#endif

// Set to 1 if user does nice job of exiting via "Quit"
// vs killing program through OS
int displayGracefulExit = 0;
int displaySaveCurrentSession = 0;

#include "mddisplay.h"
#include "dispui.h"
#include "displib.h"

typedef double Matrix[16];

#define VERSION "3.0"           /* program version */

// Not supported in MD Display 3.0
#define  SPACEBALL 0            /* Spaceball: 1 if present, 0 if not */

#define MAXHB  6000  // Max number of hydrogen bonds
#define MAXCAT 8000  /* max number of cations */
#define MAXFILT 50   // Max number of frames over which filtering is permitted

#define CATSIZE 75              /* size of cation cross */
#define VPSEP 40                /* used in SPLIT mode */

                    
// The original MD Display had a rudimentary dot surface display option
// I've removed that for MD Display 3.0 - but the left the code in for
// a future project.
// #define MAXSURF 200000          /* max number of surface points */
// static short surface[3 * MAXSURF]; // , surf = 0,
// static short surfonoff = 0;
// static short surf;
// static short surfatm[MAXATM];

// The MAX* values are defined in displib.h, and shared with preproc.c
static int bondtype[MAXDRAW];
static int masterbonds, masterBondedAtomPairs[MAXDRAW][2], bondedAtomPairsToDraw[MAXDRAW][2];
static int masterhbonds, masterhblist[MAXHB][2], hblist[MAXHB][2];

// The "drawlist" global goes back in time to move/draw days a bit - but it still works....
// I (CWM) have brought it forward to OpenGL for time being - but added #defines
// to avoid insanity

/* drawlist[i][0] == 0 if it's a move command */
/* drawlist[i][0] == 1 if it's a draw command */
/* drawlist[i][1] == atom to move or draw to */

static int drawindx = 0, drawlist[MAXDRAW][2];  /* array of moves & draws */
#define drawlist_WhichAtom(i) (drawlist[i][1])
#define drawlist_IsMove(i)    (drawlist[i][0] == 0)
#define drawlist_IsDraw(i)    (drawlist[i][0] == 1)

static int catindx, catlist[MAXCAT];   /* array of cations */

static char bondFileName[FILENAME_MAX],
            coordFileName[FILENAME_MAX],
            attribFileName[FILENAME_MAX],
            // binFileName[FILENAME_MAX],
            surfaceFileName[FILENAME_MAX];
static char sesFileName[FILENAME_MAX], setFileName[FILENAME_MAX];

/* # frames, first, last frame requested */
// You can override these with the -f option at startup
// as in -f start,end
static int frames = 9999, startFrame = 0, endFrame = 0;


// I use the OpenGL display list feature to speed up display
// of grey "static" structures dramatically.  I had experimented with converting
// the entire movie trajectory to a series of display lists - but this
// overwhelmed OpenGL on all but the simplest movies.
// This variable is Open GL "base" into "compiled" display lists of frames.
static int openGLDisplayListBaseIndex;

// Global flag to let rendering routines know that a new static structure
// display has been requested.
int staticNeedsRebuild = 1;

static int nat;      // # of atoms in the movie
static int nres;     // # of resideues (AAs or whatever) in the movie
static int bondsToDraw; // Index into actual bonds which will be drawn.  Maintained by movedraw()
static int hbonds;   /* # atoms,res,bonds,hb */
static int colour;              /* used in readparms and drawmol */

// We might want the ability to toggle depth cuing in future - not sure
// why we need it just now.  Removed in V3.0
// static int dcue = 1;            /* depth cueing */

// The showlist tracks atoms which are to be drawn with labeling
// I need to get rid of dependence on MAXATM
static int showlist[MAXATM][2], activeLabelCount;

// Globals
static int filterFlag=0; // True if we are rendering a filtered trajectory
static int filterSize=0; // sort of the # of frames over which we are filtering

// Globals which report on period box information in the input files
static int boxatm = 0;  // The subscript of the coordinate for the box in the
                        // frame of coordinates.
static int oneIfBox=0;  // True if either a periodic box or Truncated Octahedron
                        // is in use

// Note (oneIfBox && (! isTruncatedOctahedron)) implies that we have
// a periodic rectangular boundary                        
static int isTruncatedOctahedron = 0;  // True iff the periodic boundary
                                       // is a truncated octahedron
                                       
// The ps_measures store the computed distances, angles, and dihedrals
// which are dipslayed in the top right corner.  They are only used
// during the postscript output routines - which are not wired in now
// anyway.
// static float ps_measures[nMaxMeasurements];       

// The next group of variables is for maintenance of a 3D cursor
// This feature was removed from MD Display 3.0 - but I left the variables
// here in case we want it again later.

// static int curs3D = 0;          /* 3D cursor switch */
// static Coord curs3Dx, curs3Dy, curs3Dz = 0.; /* Coord of 3D cursor */
// static Matrix cursormat;        /* Transform of 3D cursor position */

// static int curs3Dvx, curs3Dvy;  /* Remembered position */

static int delbond = 2, delbondflag = 0;  /* flag for drawing H-H bonds */
// static int newsurf = 0; // We might re-add surface support later.

// Each rlab[n] contains the text description of
// that residue. Example: rlab[42] might contain "ARG\0"
typedef char rlab_type[4];
rlab_type* rlab;

// I had messed about with using display lists for movie display - maybe
// another time.  This was a command line flag to toggle this on and off.
// static int useDisplayLists = 0;  // Override with 'x' command line option

// array of attributes (colors, labeling status, etc) for each atom
// plus residue information.  If the user saves their session,
// we write all this out to disk so that we can resume a session later
atomattrib* attribs; // see displib.h
int attribytesread;

// I really ought to put all the colors in a central place so that
// you can tweak them easily.  Oh well...

// Color for Dihedral angle display text in measurement display at
// top right of screen.
static float colorDihedral[3] = {0, 150.0/255.0, 0};

// Color for Bond angle display in measurement display at
// top right of screen.
static float colorAngle[3] = {180.0/255.0, 150.0/255.0, 0};

// Color used for some points in Ramachandron plot.
static float colorMyBlue[3] = {75.0/255.0, 75.0/255.0, 1.0};

// mapcolor(MYDARKBLUE, 50, 50, 100);
float colorMyDarkBlue[3] = {50.0/255.0, 50.0/255.0, 100.0/255.0};

typedef float vector[3];        /* used in dihedral angle routine */

// Colors used for postscript output processing
// Not wired in currently
#ifdef POSTSCRIPT_DUMP_REINTEGRATED
static float psrgb[10][3] = {
{1.0, 1.0, 1.0},  /* white */
{0.0, 0.0, 0.0},  /* black */
{1.0, 0.0, 0.0},  /* red */
{0.0, 0.0, 1.0},  /* dark blue */
{0.0, 1.0, 1.0},  /* light blue */
{0.0, 1.0, 0.0},  /* green */
{1.0, 0.0, 1.0},  /* magenta */
{1.0, 1.0, 0.0},  /* yellow (orange) */
{1.0, 1.0, 0.0},  /* yellow (orange) */
{0.5, 0.5, 0.5}
};                              /* gray */
    /* (used in PostScript dump) */
#endif

// Use this to highlight whatever atom we last clicked on
int lastPickedAtom = 0;

// Use this to tell program not to reprocess highlighted atom.
int lastPickedNeedsProcessing = 0;

// ================================================
// Here are some very basic linear algebra routines
// ================================================

void VectorDiff(const vector a, const vector b, vector c)
{
// Let c = a-b;
   int m;
   for (m = 0; m < 3; m++)
      c[m] = a[m] - b[m];
}

void VectorCrossProduct(const vector a, const vector b, vector c)
{
// Let c = a x b;
   c[0] = a[1] * b[2] - a[2] * b[1];
   c[1] = a[2] * b[0] - a[0] * b[2];
   c[2] = a[0] * b[1] - a[1] * b[0];
}

float VectorDotProduct(const vector a, const vector b)
{
   // return a.b
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

float VectorNorm(vector a)
{
   // return || a ||;
   
   return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

// The parameters here are the three points a,b, and c.
// VectorAngle is a bit of a misnomer
float VectorAngle(const vector a,const vector b,const vector c)
{
   // return arcCosine((a-b).(c-b)/(norm (a-b) norm(c-b)) in degrees
   float cosine, mp, ang, dot();
   vector d1, d2;

   VectorDiff(a, b, d1);
   VectorDiff(c, b, d2);

   mp = VectorNorm(d1) * VectorNorm(d2);
   
   // Watch out for numberical problems with coincident points
   if (mp < .00001)
      mp = .00001;
   cosine = VectorDotProduct(d1, d2) / mp;
   if (cosine < (-1.0))
      cosine = (-1.0);
   if (cosine > 1.0)
      cosine = 1.0;

   // Return a value in the interval [0,180] degrees
   ang = 180.0 / M_PI * acos(cosine);
   return (ang);
}

// Let's avoid numerical troubles by flagging very short vectors
// (which are most likely due to bad structures anyway
int DistanceCheck(const vector a)
{
   float r;
   r = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];

   return (r < 0.1) ? 1 : 0;
}

/*==================================================================*/
// Return the Dihedral angle between two planes which consist of:
// plane 1: points as,bs,ds
// plane 2: points bs,cs,ds
/*==================================================================*/

float tors(const short as[3],
           const short bs[3],
           const short cs[3],
           const short ds[3])
{
   vector a, b, c, d;
   vector p, q, x, y;
   static vector z = { 0.0, 0.0, 0.0 };
   float tors_retval;
   int m;

   // Convert short coordinates to full floating point coordinates

   for (m = 0; m < 3; m++)
      {
      a[m] = (float) as[m] / KAPPA;
      b[m] = (float) bs[m] / KAPPA;
      c[m] = (float) cs[m] / KAPPA;
      d[m] = (float) ds[m] / KAPPA;
      }
   VectorDiff(a, b, p);
   VectorDiff(c, b, q);
   if (DistanceCheck(q) || DistanceCheck(p))
      return 999.0; // Things make no sense for these points
   VectorCrossProduct(p, q, x);
   VectorDiff(b, c, p);
   VectorDiff(d, c, q);
   if (DistanceCheck(q))
      return 999.0; // Things make no sense for these points
   VectorCrossProduct(p, q, y);
   tors_retval = VectorAngle(x, z, y);
   VectorCrossProduct(x, y, q);
   if (VectorDotProduct(p, q) > 0.0)
      tors_retval = (-tors_retval);
   return tors_retval;
}

// At the opening of main() we ask the user a few questions in text mode.
// This process is handled, with buffer overrun protection, here.
void SafeInputString(const char* promptString,char* instr,size_t stringSize)
{
   instr[0] = 0;
   printf("%s",promptString);
   fgets(instr,stringSize-1,stdin);
   if (instr[0] == '\n')
      instr[0] = 0;
}

// My idea with the Frame Manager code is to someday expand mddisplay so that
// it can // dynamically load movie frames, and/or display multiple
// simultaneous movies. For now this code is just an idea - not fully developed.
// In the current code, all frames are loaded in as part of FrameManagerConstruct

typedef struct {
   int frames;
   int bytesPerFrame;
   short* xx;
   FILE* f;
   int startFrame;
   int headerSize;
}  FrameManager;

FrameManager* FrameManagerConstruct(FILE* _f,int _headerSize,int _startFrame,int _frames,int _bytesPerFrame)
{
   FrameManager* fm = malloc_or_exit_if_error(sizeof(FrameManager),"FrameManagerConstruct - FrameManager");
   fm->f = _f;
   fm->startFrame = _startFrame;
   fm->headerSize = _headerSize;
   fm->frames = _frames;
   fm->bytesPerFrame = _bytesPerFrame;

   fm->xx = malloc_or_exit_if_error(_frames * _bytesPerFrame,"FrameManagerConstruct - xx");

   if (fseek(fm->f,fm->headerSize + fm->startFrame*fm->bytesPerFrame,SEEK_SET))
      {
      printf("FrameManager: Unable to fseek to position %d  errno=%d: %s ",fm->headerSize + fm->startFrame*fm->bytesPerFrame,errno,strerror(errno));
      exit(EXIT_FAILURE);
      }

   /* read frames n1-n2 */
   _frames = fread(fm->xx, fm->bytesPerFrame,fm->frames,fm->f);
   if (_frames != fm->frames)
      {
      printf("FrameManager: Unable to read %d frames",fm->frames);
      exit(EXIT_FAILURE);
      }

   return fm;
}

void FrameManagerDelete(FrameManager* fm)
{
   fclose(fm->f);
   fm->f = 0;

   free(fm->xx);
   fm->xx =0 ;

   free(fm);
}

// In future versions, this might cause a disk access to read a particular frame
// if it is not in a most recently used buffer.
const short* FrameManagerGetFrame(const FrameManager* fm,int whichFrame)
{
   assert(whichFrame >=0 && whichFrame < fm->frames);
   return &fm->xx[whichFrame*fm->bytesPerFrame/sizeof(short)];
}


// When computing filtered coordinates, or wrapping atoms into
// a periodic box or octahedron, the input movie coordinates must be modified.
// For now, until we have a better frame manager, Simply copying
// the entire movie into xxs allows xxs to subsequently be safely manipulated
/*------------------------------------------------------------*/
static void cpcoords(short* xxs, const FrameManager* fm)
{
   /* COPY COORDS DIRECTLY */
   int frame;
   for (frame=0;frame<frames;frame++)
      {
      memcpy(xxs + 3*(nat+oneIfBox)*frame,
         FrameManagerGetFrame(fm,frame),
         3* (nat + oneIfBox) * sizeof(xxs[0]));
      }
}

// Given a set of coordinates in xxs, we may need wrap them back into a rectangular
// periodic box so that they do not appear to have run away from the simulation.
// The box dimensions for each frame will vary -
// so wrapping must be done on a frame-by-frame basis.
// Another compliation is that we wrap on a residue by residue basis.  This prevents
// us from wrapping half a residue - which would look worse than having a part
// of a residue extend outside the box.

// Wrapping into the truncated octahedron is handled in a separate routine
// wrap_to() below.

void putinbox(short* xxs, int filtsize)
{
   int res = 0, ico, co, fr, i, ncord;
   short box[3], halfbox[3];

   ncord = 3 * (nat + oneIfBox);
   for (fr = 0; fr < frames - filtsize; fr++)
      {
      int frame_baseco = fr * ncord;      // Let's work on coordinates at this point - and forward
      int done_co = ncord - 3 * oneIfBox; // Don't modify any past this

      // printf("put in box %d",fr);
      /* get box dimensions for this frames */
      for (i = 0; i < 3; i++)
         box[i] = xxs[frame_baseco + boxatm * 3 - 3 + i];

      for (i = 0; i < 3; i++)
         halfbox[i] = (short)(box[i] / 2);

      // co will reference x,y,z for atom 1, x,y,z for atom 2, etc.
      for (co = 0; co < done_co; co++)
         {
         // If we are looking at a new x component
         // but it belongs to a residue we've already looked at
         // then skip ahead
         if (!(co % 3))  
            {
            while (attribs[co / 3 + 1].res == res)
               co += 3;  // Advance to next residue
               
            if (co >= done_co)
               continue;
            }                    
               
         
         // Now we are looking at the first atom of the next residue
         res = attribs[co / 3 + 1].res;
         // If the coordinate of its first atom is outside (>) of the mainbox
         if (xxs[frame_baseco + co] > halfbox[co % 3])
            {
            // Back up to the first atom in this residue
            ico = co;
            do
               {
               ico -= 3;
               }
            while (ico >= 0 && attribs[ico / 3 + 1].res == res);
            ico += 3;
            // Adjust the out-of-bounds component for
            // all atoms in this residue.
            do
               {
               xxs[frame_baseco + ico] -= box[co % 3];
               ico += 3;
               }
            while (attribs[ico / 3 + 1].res == res);
            }
         // Same algorithm as above - but if we are outside in the < direction
         else if (xxs[frame_baseco + co] < (-halfbox[co % 3]))
            {
            ico = co;
            do
               {
               ico -= 3;
               }
            while (ico >= 0 && attribs[ico / 3 + 1].res == res);
            ico += 3;
            do
               {
               xxs[frame_baseco + ico] += box[co % 3];
               ico += 3;
               }
            while (attribs[ico / 3 + 1].res == res);
            }
         }
      }
}

// Goal: -3.7 returns -4.0 and 3.7 returns 4.0
// see "man dnint" on an SGI system - or a FORTRAN manual
// This is not in most ANSI C implementations - so I had to
// recode it for version 3.0
double dnint(double x)
{
   if (x >= 0)
      return floor(x + .5);
   else
      return ceil(x - .5);
}

// The code below wraps residues back into a truncated octahedron period
// boundary.
//
// I took the code here from Amber 7's ew_box.f source code and translated to C
//
// AMBER uses a triclinic cell abstraction of all the box types which is not
// relevant to this post-AMBER code.  I got the job done by finding some simpler
// algorithms in:  I ended up bringing in
// Allen, M.P.; Tildesley, D.J.; Computer Simulation of Liquids, Oxford University Press, 1987
//   Specifically, see appendix F program F.1 (found at various web sites)
//
// van Gunsteren - lecture notes on MD simulation:
// www.chem-edu.ethz.ch/rosm/docs/rosm01-6.pdf 

void wrap_to(short* xxs,int filtsize /*,int nspm,int nsp*/)
{
      double cx,cy,cz,x0,y0,z0,dx,dy,dz,xt,yt,zt;

/*     The trunf. oct. has:
c        * center at (0,0,0) where the corner of the triclinic cell is. 
c        * one hex face with the x axis for a normal
c                       face is perpendicular to the x axis and 
c                       1/2 box(1) away from the origin.
c        * one hex face perp to xy plane, second edge vector of the
c                       triclinic cell is its normal, 109 degrees from
c                       the x axis in xy plane (-x,+y quadrant)
c
c      Approach to reconstruct the t.o. is to rotate the coordinates
c         to put the hex faces in the (+-1,+-1,+-1) normal directions, and
c         the square (diamond) faces perp to xyz axes. This is 3 rotations:
c         We did +45 around z, +(90-tetra/2) around y, +90 around x to get
c            the to oriented for the triclinic cell 
c         Now we do the opposite: -90(x), -(90-tetra/2)(y), -45(z)
c            to reproduce the original orientation, map coords into 
c            a t.o. centered at origin, then do the rotations again
c            to put it back in the orientation that matches the restrt.
c */

   int res = 0, co, fr;

   int ncord = 3 * (nat + oneIfBox);

   // Each frame can have varying box sizes - so we must take care
   // to wrap into the box for each frame.
   for (fr = 0; fr < frames - filtsize; fr++)
      {
      // The distance given in the AMBER coordinate file frames is the distance
      // between opposite HEXAGONS.  So, we adjust this using formulae in references

      // This piece of code gets box[0] from the end of the frame
      // Note that box[1] and box[2] are the same - we only deal with
      // "regular" truncated octahedrons in MDDisplay (so far!)
      double inter_hexagon_distance = xxs[fr * ncord + boxatm * 3 - 3];

      // All our algorithms work with the larger distance between opposite
      // diamonds
      double inter_diamond_distance = inter_hexagon_distance *2.0/sqrt(3.0);
      // Save a few cycles later by computing a half value, and a reciprocal now
      double idd_over_two = inter_hexagon_distance/sqrt(3.0);
      double reciprocal_idd = 1.0/inter_diamond_distance;

      // This code is the rotation code taken DIRECTLY from the AMBER ew_box.f code
      // I tried to roll my own code - and ran into trouble.  So, I'll assume this
      // is correct until someone yelps.  After this rotation, the 6 diamonds
      // are speared by the x,y,z axes.  This is the orientation required by the 
      // Allen and Tildesley, and van Gunsteren, algorithms
      double phi=M_PI/4.0;
      double cos1=cos(phi);
      double sin1=sin(phi);
      double cos2=sqrt(2.0)/sqrt(3.0);
      double sin2=1.0/sqrt(3.0);
      
      // Here is our rotation matrix   
      double t11=  cos2*cos1;
      double t12= -cos2*sin1;
      double t13= -sin2;
      double t21= -sin2*cos1;
      double t22=  sin2*sin1;
      double t23= -cos2;
      double t31=  sin1;
      double t32=  cos1;
      double t33=  0.0;
      
      // co will start by deferencing x,y,z components of mol 0, atom 1
      co = 0;
      while (co < ncord - 3 * oneIfBox)
         {
         // Amber goes to great lengths in its routine
         // to compute a geometric center and see if that is outside the periodic truncated
         // octahedron.  Here, we'll only worry about the first atom of a residue!
         // This is what MDDisplay has always done for periodic cubes anyway.  By the time
         // we get to visualization, these out-lying molecules are not hyper-critcal anyway.
      
         cx=xxs[fr * ncord + co];
         cy=xxs[fr * ncord + co+1];
         cz=xxs[fr * ncord + co+2];
//     ------Rotate---------------------------------------------        
         x0=cx*t11+cy*t21+cz*t31;
         y0=cx*t12+cy*t22+cz*t32;
         z0=cx*t13+cy*t23+cz*t33;

         // After rotation, we can now apply the van Gunsteren/Allen and Tildesley algorithms

         // This code is taken from van Gunsteren web site lecture notes
         // See his .pdf file:
         // dnint Coded above to work like Fortran
         xt=x0 - inter_diamond_distance * dnint(x0*reciprocal_idd); 
         yt=y0 - inter_diamond_distance * dnint(y0*reciprocal_idd); 
         zt=z0 - inter_diamond_distance * dnint(z0*reciprocal_idd); 
//     ------ wrap molecules external to diag faces -----
         if ((fabs(xt) + fabs(yt) + fabs(zt)) > 
             (0.75*inter_diamond_distance) )
            {
            if (xt > 0.0)
              xt -= idd_over_two;
            else
              xt += idd_over_two;
              
            if (yt > 0.0 )
              yt -= idd_over_two;
            else
              yt += idd_over_two;
            
            if (zt > 0.0 )
              zt -= idd_over_two;
            else
              zt += idd_over_two;
            }
         
//     ------get the translation in the rotated space for this
//     ------   molecules c-o-geom ----------------------------- 
//     i.e., we have a "deltax, deltay, deltaz" to apply to bring this
//     atom back into the (rotated) box.
         dx=xt-x0;
         dy=yt-y0;
         dz=zt-z0;
         
//     ------Rotate this translation back into our original space
         cx=dx*t11+dy*t12+dz*t13;
         cy=dx*t21+dy*t22+dz*t23;
         cz=dx*t31+dy*t32+dz*t33;
         
//     ------Now move the entire molecule------------------------        
         res = attribs[co / 3 + 1].res;
         while (res == attribs[co / 3 + 1].res)
            {
            // attribs[co/3+1].color = 2;
            xxs[fr * ncord + co]   += cx;
            xxs[fr * ncord + co+1] += cy;
            xxs[fr * ncord + co+2] += cz;
            co +=3;
            }
         } // end for (co /* coordinate loop)
         
      } // end for (fr...)
}      



// This routine is called by dispui.c when the user presses '?'
// or clicks on the Help button.  It pains help text on a
// white background.
void displayRenderHelpBox(void)
{
   int win_height = glutGet(GLUT_WINDOW_HEIGHT);
   int win_width = glutGet(GLUT_WINDOW_WIDTH);

   glColor3ub(0,0,255);
   DrawString(5,win_height-18,GLUT_BITMAP_HELVETICA_18,"MD Display Version ");
   dispuiDrawStringAtCurrentRasterPosition(GLUT_BITMAP_HELVETICA_18,VERSION);

   glColor3ub(0,0,0);
   DrawString(5,win_height - 32,GLUT_BITMAP_HELVETICA_10,
   "Chris Moth, Tim Callahan, Eric Swanson, and Terry Lybrand");

   glColor3ub(255,0,0);
   DrawString(5,win_height-48,GLUT_BITMAP_HELVETICA_12,"Copyright (C) 2002 by C. Moth and T. Lybrand");

   glColor3ub(0,0,0);
   DrawString(win_width-250,win_height-20,GLUT_BITMAP_HELVETICA_12," Vanderbilt University, Nashville TN");
   DrawString(win_width-264,win_height-48,GLUT_BITMAP_HELVETICA_10,"http://www.structbio.vanderbilt.edu/~cmoth/mddisplay");

   glBegin(GL_LINE_STRIP);
      glVertex3i(0,win_height-52,0);
      glVertex3i(win_width-1,win_height-52,0);
   glEnd();

   // Display keyboard examples help
   {
   static char* keyboardExamples[][2] = 
   {  
      {"#0",         "all of molecule zero (the first molecule)"},
      {":2@C32",     "the C32 atom in residue two"},
      {":3@H1@O",    "the H1 and O atoms in residue three"},
      {":1-10",      "all of residues 1 through 10"},
      {":1-5:7-10",  "all of residues 1 through 5, and 7 through 10"},
      {"@C*",        "all carbon atoms"},
      {":1-10@C*",   "all carbon atoms in residues 1 through 10"},
      {":ALA",       "all alanine residues"},
      {":GUA@CA",    "the CA atoms in all guanine residues"},
      {"#1:GUA@CA",  "the CA atoms in all guanine residues in molecule #1"},
      {":1-100@bb",  "the backbone atoms in residues 1 through 100"},
      {":3-10@sc",   "the side chains of residues 3-10"},
      {0,0} // <--- this 0 ends display loop
   };

   int base_offset = 80;
   int example;

   glColor3ub(0,0,255);
   DrawString(2,win_height-base_offset,GLUT_BITMAP_HELVETICA_12,"Keyboard Input Examples");
   glColor3ub(0,0,0);
   glBegin(GL_LINE_STRIP);
      glVertex3i(0,win_height-base_offset-2,0);
      glVertex3i(300,win_height-base_offset-2,0);
   glEnd();

   base_offset += 20;

   for (example=0;keyboardExamples[example][0];example++)
      {
      glColor3ub(0,0,0);
      DrawString(4,win_height-base_offset-example*15,GLUT_BITMAP_9_BY_15,keyboardExamples[example][0]);
      glColor3ub(25,25,25);
      DrawString(4+100,win_height-base_offset-example*15,GLUT_BITMAP_TIMES_ROMAN_10,keyboardExamples[example][1]);
      }
   } // End Keyboard Examples
   

   

   // Display mouse examples help
   {
   static char* mouseExamples[][2] = 
   {  
      {"Left","Pick Atoms or Controls"},
      {"Middle","Translate Scene"},
      {"Right","Rotate Scene"},
      {"Left + Right","Scale Scene"},
      {"",""},
      {"SHIFT+Left","Translate Scene"},

      {0,0} // <--- this 0 ends display loop
   };

   int example;
   int base_offset = 80;
   int right_offset = 275;

   glColor3ub(0,0,255);
   DrawString(win_width-right_offset,win_height-base_offset,GLUT_BITMAP_HELVETICA_12,"Mouse Button Functions");
   glColor3ub(0,0,0);
   glBegin(GL_LINE_STRIP);
      glVertex3i(win_width-right_offset,win_height-base_offset-2,0);
      glVertex3i(win_width-1,win_height-base_offset-2,0);
   glEnd();

   base_offset += 20;

   for (example=0;mouseExamples[example][0];example++)
      {
      glColor3ub(0,0,0);
      DrawString(win_width-right_offset,win_height-base_offset-example*15,GLUT_BITMAP_HELVETICA_12,mouseExamples[example][0]);
      glColor3ub(25,25,25);
      DrawString(win_width-right_offset+100,win_height-base_offset-example*15,GLUT_BITMAP_TIMES_ROMAN_10,mouseExamples[example][1]);
      }
   } // End Mouse Examples

   // Key Character
   {
   static char* keyChars[][2] =
   {
      {"s","superimpose (smear)"},
      {"S","cancel smearing"},
      {"F","Filter coordinates"},
      {"f","toggle filtered movie"},
      {"H","change H-bond cutoff"},

      {"d","next display mode"},
      {"F1","Mono"},
      {"F2","Stereo Pairs"},
      {"F3","TRI"},
      {"F4","True Stereo (if supported)"},

      {"c",  "color palette"},
      {"0-9","solid colors"},
      {"b","color by atom"},
      {"r or 0","remove (color black)"},
      {"",""},


      {0,0} // <--- this 0 ends display loop
   };

   int base_offset = 300;
   int row;

   glColor3ub(0,0,255);
   DrawString(2,win_height-base_offset,GLUT_BITMAP_HELVETICA_12,"Additional Keyboard Command Key Information");
   glColor3ub(0,0,0);
   glBegin(GL_LINE_STRIP);
      glVertex3i(0,win_height-base_offset-2,0);
      glVertex3i(win_width-1,win_height-base_offset-2,0);
   glEnd();

   glBegin(GL_LINE_STRIP);
      glVertex3i(4+200-12,win_height-base_offset-2,0);
      glVertex3i(4+200-12,0,0);
   glEnd();

   glBegin(GL_LINE_STRIP);
      glVertex3i(4+400-12,win_height-base_offset-2,0);
      glVertex3i(4+400-12,0,0);
   glEnd();

   base_offset += 20;

   for (row=0;keyChars[row*3][0];row++)
      {
      int column;
      for (column=0;column<3;column++)
         {
         int example = column*5 + row;
         static int offset[3] = {20,25,45};
         glColor3ub(0,0,0);
         DrawString(4+column*200,win_height-base_offset-row*15,GLUT_BITMAP_HELVETICA_10,keyChars[example][0]);
         glColor3ub(25,25,25);
         DrawString(column*200+offset[column],win_height-base_offset-row*15,GLUT_BITMAP_HELVETICA_10,keyChars[example][1]);
         }

      }
   }
}


/*------------------------------------------------------------*/

/* MAS 3-9-90 parsing expects the following order of commands MAS # > : > @,
 * where '#' refers to molecule number (0-9999), MAS ':' refers to residue
 * number (1-9999), and '@' refers to MAS atom type.  yes, molecule zero is
 * the first molecule.
 *
 * TJC 8-1-90 TJC You can now enter multiple ranges, e.g. #1:20-30:35-40.
 * TJC
 *
 * MAS MAS mas_parse() MAS INPUT: str: a string of the form descibed above
 * for applying MAS some action to a specified subset of the atoms MAS
 * displayed. MAS OUTPUT: arry: an array of resatm structures which specify
 * the MAS subset of atoms to be modified.  see definition MAS of resatm
 * structure at top of program. MAS PURPOSE: to translate the code (of the
 * form #n-m:o-p@w) into MAS "bounds" specified in the relevant fields of MAS
 * a resatm structure. */

// Version 3 note from Chris Moth:
// This format is used for all text input during specification of residues
// See the online help '?', and the .pdf manual, for more information about how
// codes like @SC or :ARG etc are interpreted by MDDisplay.
// The goal of this code is to figure out what range of molecules, residues
// and atoms that the user means by their text input.
// Each arry[] element represents data about ranges of data.  Users
// can specify multiple groups by separating them with commas.

// resatm below is the structure of a unit of parsed input text data

typedef struct {
   int lowmol; // i.e. starting molecule number
   int himol;  // i.e. ending molecule number
   int lores;  // i.e. starting residue number
   int hires;  // i.e. ending residue number
   char reslab[4];   // i.e. "ARG" or other residue label
   char atmlab[5];   // "CA" or other residue label
} resatm;

// Parse up to 20 resatm records form the input text.
void mas_parse(const char* str1, resatm arry[20])
{
   int n, k, i, j = 0, m, curri;
   int hotres = 0; // TRUE if a residue is explcitly specified
   int hotmol = 0; // TRUE if a molecule has been explicity specified
   int prevcmd0;
   int lowres = 1, highres;
   int lowmol = 0, himol = 99999;   /* MAS: note that molecule numbers start
                                     * at 0. */
   char reslab[10], cmdstr[20];
   
   int len = strlen(str1);

   // I wanted to make the parsed string const, so duplicating it
   // here seems like the right thing to do.  I free it at the end
   char* str = strdup(str1);
   char* s= str;
   if (str == 0)
      {
      printf("Unable to execute strdup");
      exit(EXIT_FAILURE);
      }

   // Capitalize this string.
   // Still not sure how this is going to work out with Calcium counter
   // ions and C-alpha atoms.  Oh well.
   while (*s)
      {
      if (islower(*s))
         *s = toupper(*s);
      s++;
      }

   highres = nres + 1;          /* TJC saves a little computation */

   i = 0;

   /* strip leading garbage */
   while (! strchr("#:@",str[i]) && (i < len))
      i++;

   // By default, we want the user input to apply to all residue names
   reslab[0] = '*';
   reslab[1] = 0;
   
   prevcmd0 = 0;

   while (i < len)
   {
      /* MAS 3-4-90 parse() MAS parsing for the molecule number, after '#' */
      memset(cmdstr,0,sizeof(cmdstr)); // Zero out all of cmdstr
      k = 0;
      // build up cmdstr with the non-punctuation characters.
      while ((i< len) && (k < sizeof(cmdstr)-2) && (k == 0 || (! strchr("#:@-,",str[i]))))
         {
         cmdstr[k++] = toupper(str[i++]);
         }
      cmdstr[k] = 0;

      if (cmdstr[0] == ',')
           {
           cmdstr[0] = (char)prevcmd0;
           }

      // User can enter #2 to mean "molecule 2"
      if (cmdstr[0] == '#')
      {
         prevcmd0 = cmdstr[0];
         // If we had earlier parse data, save it now
         // and increment j.  We'll get the next piece now
         if (hotres || hotmol)
            {                      /* TJC */
            arry[j].lowmol = lowmol;
            arry[j].himol = himol;
            arry[j].lores = lowres;
            arry[j].hires = highres;
            for (m = 0; m < 3; m++)
               arry[j].reslab[m] = reslab[m];
            arry[j].reslab[3] = 0;
            arry[j++].atmlab[0] = 0;
            }

         // Let's go get input of molecule # following the '#'
         lowmol = 0;
         himol = 99999;
         lowres = 1;
         highres = nres + 1;
         reslab[0] = '*';
         reslab[1] = 0;

         if (cmdstr[1] >= '0' && cmdstr[1] <= '9')
            {
            sscanf(cmdstr + 1, "%d", &n); /* low mol entered */
            lowmol = n;
            himol = n;
            }
         /* see if there's been a range specified */
         /* skip over blanks */
         curri = i;
         k = 0;
         while ((! strchr("#:@,",str[i])) && i < len)
            {
            cmdstr[k++] = toupper(str[i++]);
            }
         cmdstr[k] = 0;
         if (cmdstr[0] == '-')
            {                   /* high molecule specified */
            sscanf(cmdstr + 1, "%d", &n);
            himol = n;
            }
         else
            {
            i = curri;
            }
         // For next input parse, note that a molecule or range
         // has been explcitly entered!
         hotmol = 1;
      }

      // Did the user specify a specific residue
      if (cmdstr[0] == ':')
      {
         prevcmd0 = cmdstr[0];

         // If input is going on for an earlier residue, be sure
         // to specify that one and increment j
         if (hotres)
            {/* TJC */
            arry[j].lowmol = lowmol;
            arry[j].himol = himol;
            arry[j].lores = lowres;
            arry[j].hires = highres;
            for (m = 0; m < 3; m++)
               arry[j].reslab[m] = reslab[m];
            arry[j].reslab[3] = 0;
            arry[j++].atmlab[0] = 0;
            }

         lowres = 1;
         highres = nres + 1;
         reslab[0] = '*';
         reslab[1] = 0;

         if (isupper(cmdstr[1]))
         {                      /* residue label */
            for (m = 0; m < 3; m++)
               reslab[m] = cmdstr[m + 1];
            reslab[3] = 0;
         }
         else if (isdigit(cmdstr[1]))
         {
            sscanf(cmdstr + 1, "%d", &n); /* (low) residue */
            lowres = n;
            highres = n;
            /* see if there's been a range specified */
            /* skip over the current command */
            curri = i;
            k = 0;
            while ((! strchr("#:@,",str[i])) && (i < len))
               {
               cmdstr[k++] = toupper(str[i++]);
               }
            if (cmdstr[0] == '-')
            {                   /* high residue specified */
               sscanf(cmdstr + 1, "%d", &n);
               highres = n;
            }
            else
            {
               i = curri;
            }
         }
         hotres = 1;
      }

      // If the user has specified an atom, that certainly
      // completes an input sequence.
      if (cmdstr[0] == '@')
      {                         /* atom specified */
         prevcmd0 = cmdstr[0];
         k = 0;
         arry[j].lowmol = lowmol;
         arry[j].himol = himol;
         arry[j].lores = lowres;
         arry[j].hires = highres;
         while (cmdstr[k + 1] && k < 4)
         {
            arry[j].atmlab[k] = cmdstr[k+ 1];
            k++;
         }
         while (k < 4)
            {
            arry[j].atmlab[k++] = ' ';
            }
         for (m = 0; m < 3; m++)
            arry[j].reslab[m] = reslab[m];
         arry[j].reslab[3] = 0;
         // printf("--->%x %x %x %x<----\n",arry[j].atmlab[0],arry[j].atmlab[1],arry[j].atmlab[2],arry[j].atmlab[3]);
         arry[j++].atmlab[4] = 0;
         
         hotres = 0;
         hotmol = 0;
      }
   }

   // To get here, we did not explicity enter an atom.  That's fine
   // Save the specified molecules and residues
   if (hotres || hotmol)
   {
      arry[j].lowmol = lowmol;
      arry[j].himol = himol;
      arry[j].lores = lowres;
      arry[j].hires = highres;
      for (m = 0; m < 3; m++)
         arry[j].reslab[m] = reslab[m];
      arry[j].reslab[3] = 0;
      arry[j++].atmlab[0] = 0;
   }
   
   // This 0 marks the end of parsing
   arry[j].lores = 0;
   free(str);
}                               

// Return 0 if *a == *b, 1 otherwise.
// taking into account windcalds and the BB and SC
// atom specification strings.
int stcmp(const char* a, const char* b)
{
   int i = 0;

   // If we are comparing to "BB" then match on all backbone atoms
   if (b[0] == 'B' && b[1] == 'B')
      {
      // Backbone atoms are N, C, and CA
      if ((a[0] == 'C' && a[1] <= 'A') ||
          (a[0] == 'N' && a[1] < 'A'))
         return 0;
      else
         return 1;
      }
   else if (b[0] == 'S' && b[1] == 'C') // If we are matching to Sidechain
      {
      // Backbone atoms are N, C, and CA
      if ((a[0] == 'C' && a[1] <= 'A') ||
          (a[0] == 'N' && a[1] < 'A'))
         return 1;
      else
         return 0;
      }
   else
      {
      // Match ABC to A*
      // Match ABbb to ABbb
      while (a[i] == b[i] && i < 4)
         i++;
      if (i == 4 || b[i] == '*')
         return 0;
      else
         return 1;
      }
}

#ifdef OLD_STUFF

/*----------------------------------------------------------*/
// The static rendering is handled elwewhere now in V3.0
// I've left some code here because it refers to surface rendering.
// Perhaps we can bring this back in future.

void dostatic(void)
{
   int i, j, n=0;
   static short bit, active;
   long surfmenu;
   int low, high;
   resatm resatmarry[20];

   if (butt)
   {
      if (mx > 0 && mx < 100 && my > 780 && my < 810)
      {
         if (frames > 1)
            newstatic = 2;
         if (surf && frames < 2)
         {
#ifdef GL_STUFF
            if (surfonoff)
               surfmenu =
                  defpup("Dot Surface %t|Off|Remove Some|Restore Some");
            else
               surfmenu =
                  defpup("Dot Surface %t|On|Remove Some|Restore Some");
            n = dopup(surfmenu);
#endif
            if (n == 1)
               surfonoff = (short)(!surfonoff);
            if (n == 2)
            {
               bit = 0;
               active = 1;
            }
            if (n == 3)
            {
               bit = 1;
               active = 1;
               surfonoff = 1;
            }
            newsurf = 1;
         }
      }
   }

   if (active)
   {
#ifdef GL_STUFF
      if (bit)
         getinput("Enter #mol:res@atoms to restore surface to:");
      else
         getinput("Enter #mol:res@atoms to remove surface from:");
      if (answerready)
      {
         mas_parse(getinput(""), resatmarry);
         i = 0;
         while (resatmarry[i].lores)
         {
            low = resatmarry[i].lores;
            high = resatmarry[i].hires;
            if (resatmarry[i].atmlab[0])
            {                   /* atoms specified */
               for (j = 1; j <= nat; j++)
               {
                  if (attribs[j].res >= low && attribs[j].res <= high
                      && !stcmp(attribs[j].label, resatmarry[i].atmlab)
                      && !stcmp(rlab[attribs[j].res], resatmarry[i].reslab))
                     surfatm[j] = bit;
               }
            }
            else
            {                   /* no atoms specified, just residues */
               for (j = 1; j <= nat; j++)
               {
                  if (attribs[j].res >= low && attribs[j].res <= high
                      && !stcmp(rlab[attribs[j].res], resatmarry[i].reslab))
                     surfatm[j] = bit;
               }
            }
            i++;
         }
         active = 0;
      }
#endif
   }
}
#endif

/*------------------------------------------------------------*/



#ifdef POSTSCRIPT_DUMP_REINTEGRATED
/*===================================================================*/
// Used by makeps - this expects each element to be an array
// of 6 floats...
int FloatCompareFunction(const void* ap_asvoid, const void* bp_asvoid)
{
   int ans;
   float* ap = (float*)ap_asvoid;
   float* bp = (float*)bp_asvoid;

   float diff = bp[4] - ap[4];

   if (diff > 0.0)
      ans = 1;
   else if (diff < 0.0)
      ans = -1;
   else
      ans = 0;

   return ans;

#ifdef OLD_CODE
   // CWM thinks this is to MSB/LSB specific!

   if (ap[4] < bp[4])
      ans = 1;
   else
      ans = (-1);

   return (ans);
#endif
}
#endif

/*===================================================================*/
/* This wil be useful if we ever add back the ability to dump postscript
   output */
/*===================================================================*/
#ifdef POSTSCRIPT_DUMP_REINTEGRATED
void getpscolor(int dc,
                int bwcol,
                int wob,
                float z,
                int clr,
                float* rop,
                float* gop,
                float* bop)
{
   float ro, go, bo;

   if (dc)
   {
      if (bwcol && wob)
      {
         ro = psrgb[clr][0] * (-z / 2.0 + 0.5);
         go = psrgb[clr][1] * (-z / 2.0 + 0.5); /* depthcued color */
         bo = psrgb[clr][2] * (-z / 2.0 + 0.5); /* black background */
      }
      if (bwcol && !wob)
      {
         ro = 1.0 - (1.0 - psrgb[clr][0]) * (-z / 2.0 + 0.5);
         go = 1.0 - (1.0 - psrgb[clr][1]) * (-z / 2.0 + 0.5);  /* dc color */
         bo = 1.0 - (1.0 - psrgb[clr][2]) * (-z / 2.0 + 0.5);  /* on white */
      }
      if (!bwcol)
      {
         ro = (z / 2.0 + 0.5);
         go = ro;               /* black & white depthcueing */
         bo = ro;
      }
   }
   else
   {
      if (bwcol)
      {
         ro = psrgb[clr][0];
         go = psrgb[clr][1];    /* color but no depth-cueing */
         bo = psrgb[clr][2];
      }
      else
      {
         ro = 2.0 * (-z / 2.0 + 0.5);
         go = ro;
         bo = ro;
      }
   }

   *rop = ro;
   *gop = go;
   *bop = bo;
}
#endif


/*------------------------------------------------------------*/
/* CREATE A POSTSCRIPT FILE SHOWING IMAGE ON SCREEN */
// More postscript code.  This is just simply not wire in at present
/*------------------------------------------------------------*/

#ifdef POSTSCRIPT_DUMP_REINTEGRATED
typedef float doublevector[6];
void display_makeps(dispinfo* vp, const short* xxx, Matrix mat1, Matrix mat2)
{
   vector* tx;
   float old_coord[4], new_coord[4];   /* holds transformed atom
                                           * coordinates */
   doublevector* segx; // Coordinates of line segments
   int res, atm, nseg, ster = 0, iter;
   FILE *ps;
   float ro = 0.0, go = 0.0, bo = 0.0, roprev=0.0, goprev=0.0, boprev = 0.0, z, dis;
   int a1, a2, a3, a4, n, clr, r, i, j, ii, kc;
   float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
   float slope, nor, xc, yc, xp, yp, sep = 0.02;
   int bwcol = 0, dc = 1, ss = 0, wob = 0, psborder = 0;
   // int bgn, end;
   char inst[20], psname[50];
   Matrix mat;

   tx = malloc_or_exit_if_error(MAXATM*sizeof(vector),"display_makeps:tx");
   segx = malloc_or_exit_if_error(MAXDRAW * sizeof(doublevector),"makeps:segx");

   if (mat2)
      ster = 1;

#ifdef GL_STUFF
   clearall();
   tpon();
#endif

   /* ----- GET INFO ----- */

   printf("Enter name of PostScript file to create: ");
   gets(psname);
   if (psname[0])
   {

      printf("Color or black & white? [c/b] ");
      gets(inst);
      if (inst[0] == 'c' || inst[0] == 'C')
         bwcol = 1;

      if (bwcol)
      {
         printf("On white or black background? [w/b] ");
         gets(inst);
         if (inst[0] == 'b' || inst[0] == 'B')
         {
            wob = 1;
            /* reverse black and white */
            for (i = 0; i < 3; i++)
               psrgb[0][i] = 0.0;
            for (i = 0; i < 3; i++)
               psrgb[1][i] = 1.0;
         }
      }

      printf("Depth-cueing? [y/n] ");
      gets(inst);
      if (inst[0] == 'y' || inst[0] == 'Y')
      {
         if (!bwcol)
         {
            printf("Depth cue by line width or grayscale? [w/g] ");
            gets(inst);
            if (inst[0] == 'g' || inst[0] == 'G')
               dc = 2;
            else
               dc = 1;
         }
      }
      else
         dc = 0;

      if (ster)
      {
         printf("Small stereo (60 mm separation) ? [y/n] ");
         gets(inst);
         if (inst[0] == 'y' || inst[0] == 'Y')
            ss = 1;
      }
      else
      {
         printf("Border ? [y/n] ");
         gets(inst);
         if (inst[0] == 'y' || inst[0] == 'Y')
            psborder = 1;
      }

      printf
         ("\n\nPlease wait while the PostScript file is being created...\n");

      /* ----- PostScript HEADER ----- */

      ps = fopen(psname, "wb");
      fprintf(ps, "%%!\n");
      fprintf(ps, "/seg { newpath 1.5 setlinewidth\n");
      fprintf(ps, "       setrgbcolor moveto lineto stroke} def\n");
      fprintf(ps, "/vseg { newpath setlinewidth\n");
      fprintf(ps, "        moveto lineto stroke} def\n");
      fprintf(ps, "/pnt { newpath setrgbcolor 0.5 setlinewidth\n");
      fprintf(ps, "        moveto 0.5 0 rlineto stroke} def\n");
      fprintf(ps, "/pointsize {\n");
      fprintf(ps, "     /Helvetica findfont\n");
      fprintf(ps, "     exch scalefont setfont\n");
      fprintf(ps, "} def\n");

      if (wob)
      {
         fprintf(ps, "0 0 moveto\n");
         fprintf(ps, "612 0 lineto\n");
         fprintf(ps, "612 792 lineto\n");
         fprintf(ps, "0 792 lineto\n");
         fprintf(ps, "0.0 setgray closepath fill\n"); /* make background
                                                       * black */
      }
      /* if depth cueing, set finer grain gray scale (improves color too) */
      if (dc)
         fprintf(ps, "currentscreen 3 -1 roll 2 mul 3 1 roll setscreen\n");
      fprintf(ps, "8 pointsize\n");
      fprintf(ps, "90 rotate\n");
      fprintf(ps, "0 setgray\n");
      fprintf(ps, "1 setlinecap\n");   /* rounded line ends */
      if (psborder)
      {
         fprintf(ps, "100 -540 translate\n");
         fprintf(ps, "1.2 0.92 scale\n"); /* aspect ratio */
      }
      else
      {
         fprintf(ps, " 36 -582 translate\n");
         fprintf(ps, "1.44 1.104 scale\n");
      }
      if (ss)
      {
         fprintf(ps, "0.55 0.55 scale\n");   /* shrink down for 60 mm
                                              * separation */
         fprintf(ps, "200 250 translate\n"); /* re-center images */
      }

      /* BORDER & CREDITS */

      /* (don't print border for small stereo since it is distracting) */

      if (psborder && !ss)
      {
         fprintf(ps, "0.0 setgray\n");
         fprintf(ps, "0 0 moveto\n");
         fprintf(ps, "500 0 lineto\n");
         fprintf(ps, "500 500 lineto\n");
         fprintf(ps, "0 500 lineto\n");
         fprintf(ps, "0 0 lineto\n");
         fprintf(ps, "0 -10 moveto\n");
         fprintf(ps, "( MD Display) show\n");
         fprintf(ps, "stroke\n");

         fprintf(ps, "0 0 moveto\n");
         fprintf(ps, "500 0 lineto\n");
         fprintf(ps, "500 500 lineto\n");
         fprintf(ps, "0 500 lineto\n");
         fprintf(ps, "closepath clip newpath\n");  /* keep it in the box */
      }

      if (ster)
         fprintf(ps, "0.5 1.0 scale\n");  /* half size for stereo pair */

      for (iter = 0; iter <= ster; iter++)
      {

         nseg = 0;

         if (iter)
         {                      /* mat=mat2; get correct matrix */
         memcpy(mat,mat2,sizeof(mat));
         }

         else
         {                      /* mat=mat1 */
         memcpy(mat,mat1,sizeof(mat));
         }

         if (iter)
            fprintf(ps, "500 0 translate\n");   /* shift right image */

         /* ------ DO DOT SURFACE IF PRESENT ----- */
#ifdef SURFACE_FILE_BROUGHT_BACK
         if (surfonoff)
         {
            printf("Adding dot surface to file...\n");
            old_coord[3] = 1.0;
            for (i = 1; i <= nat; i++)
            {
               clr = attribs[i].color;
               if (clr && surfatm[i])
               {
                  bgn = (int) (surface[2 * i] + 10000 * surface[2 * i + 1]);
                  end =
                     (int) (surface[2 * i + 2] + 10000 * surface[2 * i + 3]);
                  for (j = bgn; j < end; j += 3)
                  {

                     for (r = 0; r < 3; r++)
                        old_coord[r] = (float) surface[j + r];
                     for (r = 0; r < 4; r++)
                     {
                        new_coord[r] = 0.0;
                        // REVIEW - Careful... Note CWM's new Matrix deference code
                        for (c = 0; c < 4; c++)
                           new_coord[r] += old_coord[c] * mat[c*4+r];
                     }
                     for (r = 0; r < 3; r++)
                        tx[0][r] = new_coord[r] / new_coord[3];   /* NORMALIZE */
                     x1 = tx[0][0];
                     y1 = tx[0][1];
                     z1 = tx[0][2]; /* z1 varies between -1.0 and 1.0 */
                     if (z1 > (-1.0) && z1 < 1.0)
                     {
                        segx[nseg][0] = x1;
                        segx[nseg][1] = y1;
                        segx[nseg][2] = x1 + (0.0 / 250.0);
                        segx[nseg][3] = y1;
                        segx[nseg][4] = z1;
                        segx[nseg][5] = (float) clr;
                        nseg++;
                     }
                  }
               }
            }
         }
#endif
         /* -----TRANSFORM EACH ATOM'S COORDS----- */

         old_coord[3] = 1.0;
         for (i = 1; i <= nat; i++)
         {
            for (r = 0; r < 3; r++)
               old_coord[r] = (float) xxx[i * 3 - 3 + r];
            for (r = 0; r < 4; r++)
            {
               new_coord[r] = 0.0;
               new_coord[r] = old_coord[r];
               // for (c = 0; c < 4; c++)
                  // new_coord[r] += old_coord[c] * mat[c*4+r];
            }
            for (r = 0; r < 3; r++)
               tx[i][r] = new_coord[r] / new_coord[3];   /* NORMALIZE */
         }

         /* -----DIVIDE EACH BOND INTO MANY SMALL SEGMENTS----- */

         a1 = bondedAtomPairsToDraw[0][0];
         a2 = bondedAtomPairsToDraw[0][1];
         dis = (tx[a1][0] - tx[a2][0]) * (tx[a1][0] - tx[a2][0]) +
            (tx[a1][1] - tx[a2][1]) * (tx[a1][1] - tx[a2][1]) +
            (tx[a1][2] - tx[a2][2]) * (tx[a1][2] - tx[a2][2]);
         dis = sqrt(dis);
         n = (int) (100.0 * dis);   /* n is number of segments to divide each
                                     * bond into */
         if (n < 2)
            n = 2;
         printf("n=%d\n", n);
         for (i = 0; i < bondsToDraw; i++)
         {
            a1 = bondedAtomPairsToDraw[i][0];
            a2 = bondedAtomPairsToDraw[i][1];
            for (j = 0; j < n; j++)
            {
               x1 =
                  tx[a1][0] + ((float) j / (float) n) * (tx[a2][0] -
                                                         tx[a1][0]);
               x2 =
                  tx[a1][0] + ((float) (1 + j) / (float) n) * (tx[a2][0] -
                                                               tx[a1][0]);
               y1 =
                  tx[a1][1] + ((float) j / (float) n) * (tx[a2][1] -
                                                         tx[a1][1]);
               y2 =
                  tx[a1][1] + ((float) (1 + j) / (float) n) * (tx[a2][1] -
                                                               tx[a1][1]);
               z1 =
                  tx[a1][2] + ((float) j / (float) n) * (tx[a2][2] -
                                                         tx[a1][2]);
               z2 =
                  tx[a1][2] + ((float) (1 + j) / (float) n) * (tx[a2][2] -
                                                               tx[a1][2]);
               if (fabs(x1) <= 1.0 && fabs(y1) <= 1.0 && fabs(z1) <= 1.0 &&
                   fabs(x2) <= 1.0 && fabs(y2) <= 1.0 && fabs(z2) <= 1.0)
               {
                  segx[nseg][0] = x1;
                  if (j == 0)
                     segx[nseg][0] += 0.0001;   /* forces new line */
                  segx[nseg][1] = y1;
                  segx[nseg][2] = x2;
                  segx[nseg][3] = y2;
                  segx[nseg][4] = (z2 + z1) / 2.0; /* average z coord */
                  if (j < n / 2)
                     segx[nseg][5] = (float) attribs[a1].color;
                  else
                     segx[nseg][5] = (float) attribs[a2].color;
                  nseg++;
               }
            }
         }

         /* ---- SORT (if necessary) ---- */

         if (bwcol || dc == 2)
         {
            printf("Sorting segments...\n");

            // This assert makes sure that we have no accidental maintenace
            // problem.
            assert(sizeof(segx[0]) == sizeof(float)*6);
            qsort(segx, nseg, sizeof (segx[0]), FloatCompareFunction);
            printf("Done.\n");
         }

         /* ---- WRITE TO FILE ---- */
         /* An attempt is made to condense segments where possible, */
         /* otherwise the postscript files may be voluminous */

         printf("Writing file.\n");
         ii = 0;
         for (i = 0; i <= nseg; i++)
         {
            if (i < nseg)
            {
               clr = (int) segx[i][5];
               z = segx[i][4];  /* z varies between -1.0 and 1.0 */
               getpscolor(dc, bwcol, wob, z, clr, &ro, &go, &bo);
            }
            if (bwcol || dc == 2)
            {
               kc = 20 * ro + .5;
               ro = (float) kc / 20.;
               kc = 20 * go + .5;
               go = (float) kc / 20.;
               kc = 20 * bo + .5;
               bo = (float) kc / 20.;
               if (i > 0)
               {
                  if (segx[i][0] != segx[i - 1][2] ||
                      segx[i][1] != segx[i - 1][3] || ro != roprev ||
                      go != goprev || bo != boprev || i == nseg)
                  {
                     fprintf(ps, "%3.1f %3.1f %3.1f %3.1f %3.2f %3.2f %3.2f seg\n", (1.0 + segx[ii][0]) * 250.0, (1.0 + segx[ii][1]) * 250.0, (1.0 + segx[i - 1][2]) * 250.0, /* coords
                                                                                                                                                                               * in
                                                                                                                                                                               * range
                                                                                                                                                                               * 0<x<500
                                                                                                                                                                               */
                             (1.0 + segx[i - 1][3]) * 250.0,
                             roprev, goprev, boprev);
                     ii = i;
                  }
               }
               roprev = ro;
               goprev = go;
               boprev = bo;
            }
            else
            {
               if (dc)
               {
                  kc = 10 * ro + .5;
                  ro = 1.0 - (float) kc / 10.;
               }
               else
                  ro = 1.0;
               if (i > 0)
               {
                  if (segx[i][0] != segx[i - 1][2] ||
                      segx[i][1] != segx[i - 1][3] || ro != roprev ||
                      i == nseg)
                  {
                     fprintf(ps, "%3.1f %3.1f %3.1f %3.1f %3.1f vseg\n", (1.0 + segx[ii][0]) * 250.0, (1.0 + segx[ii][1]) * 250.0, (1.0 + segx[i - 1][2]) * 250.0,   /* coords
                                                                                                                                                                      * in
                                                                                                                                                                      * range
                                                                                                                                                                      * 0<x<500
                                                                                                                                                                      */
                             (1.0 + segx[i - 1][3]) * 250.0, roprev);
                     ii = i;
                  }
               }
               roprev = ro;
            }
         }

     /*--- DISTANCES & DIHEDRALS ---*/

         for (i = 0; i < vp->nmeas; i++)
         {

            a1 = vp->measarry[i][0];
            a2 = vp->measarry[i][1];
            a3 = vp->measarry[i][2];
            a4 = vp->measarry[i][3];
            x1 = tx[a1][0];
            y1 = tx[a1][1];
            z1 = tx[a1][2];
            x2 = tx[a2][0];
            y2 = tx[a2][1];
            z2 = tx[a2][2];
            x3 = tx[a3][0];
            y3 = tx[a3][1];
            z3 = tx[a3][2];
            x4 = tx[a4][0];
            y4 = tx[a4][1];
            z4 = tx[a4][2];

            if (!a3)
            {                   /* DISTANCE */
               fprintf(ps, "%d setgray\n", wob);
               fprintf(ps, "[6 4] 0 setdash\n");
               if (fabs(z1) < 1.0 && fabs(z1) < 1.0)
               {
                  xc = (x1 + x2) / 2.0;
                  yc = (y1 + y2) / 2.0;
                  slope = (y1 - y2) / (x1 - x2);
                  nor = (-1.0) / slope;
                  if (fabs(nor) > 1.0)
                  {
                     if (nor > 0.0)
                        yp = yc + sep;
                     else
                        yp = yc - sep;
                     xp = xc + sep / fabs(nor);
                  }
                  else
                  {
                     xp = xc + sep;
                     yp = yc + sep * nor;
                  }
                  fprintf(ps, "%3.2f %3.2f %3.2f %3.2f %3.2f vseg\n",
                          (1.0 + x1) * 250.0, (1.0 + y1) * 250.0,
                          (1.0 + x2) * 250.0, (1.0 + y2) * 250.0, 1.5);
                  fprintf(ps, "%3.2f %3.2f moveto\n",
                          (1.0 + xp) * 250.0, (1.0 + yp) * 250.0);
                  fprintf(ps, "(%3.2f A) show\n", ps_measures[i]);
                  fprintf(ps, "[] 0 setdash\n");
               }

            }
            else if (a4)
            {                   /* DIHEDRAL */

               fprintf(ps, "%d setgray\n", wob);
               fprintf(ps, "[6 4] 0 setdash\n");
               if (fabs(z1) < 1.0 && fabs(z2) < 1.0 && fabs(z3) < 1.0 &&
                   fabs(z4) < 1.0)
               {
                  xc = (x2 + x3) / 2.0;
                  yc = (y2 + y3) / 2.0;
                  slope = (y2 - y3) / (x2 - x3);
                  nor = (-1.0) / slope;
                  if (fabs(nor) > 1.0)
                  {
                     if (nor > 0.0)
                        yp = yc + sep;
                     else
                        yp = yc - sep;
                     xp = xc + sep / fabs(nor);
                  }
                  else
                  {
                     xp = xc + sep;
                     yp = yc + sep * nor;
                  }
                  fprintf(ps, "%3.2f %3.2f %3.2f %3.2f %3.2f vseg\n",
                          (1.0 + x1) * 250.0, (1.0 + y1) * 250.0,
                          (1.0 + x2) * 250.0, (1.0 + y2) * 250.0, 1.5);
                  fprintf(ps, "%3.2f %3.2f %3.2f %3.2f %3.2f vseg\n",
                          (1.0 + x2) * 250.0, (1.0 + y2) * 250.0,
                          (1.0 + x3) * 250.0, (1.0 + y3) * 250.0, 1.5);
                  fprintf(ps, "%3.2f %3.2f %3.2f %3.2f %3.2f vseg\n",
                          (1.0 + x3) * 250.0, (1.0 + y3) * 250.0,
                          (1.0 + x4) * 250.0, (1.0 + y4) * 250.0, 1.5);
                  fprintf(ps, "%3.2f %3.2f moveto\n",
                          (1.0 + xp) * 250.0, (1.0 + yp) * 250.0);
                  fprintf(ps, "(%3.2f) show\n", ps_measures[i]);
                  fprintf(ps, "[] 0 setdash\n");
               }

            }
            else
            {                   /* ANGLE */

               fprintf(ps, "%d setgray\n", wob);
               fprintf(ps, "[6 4] 0 setdash\n");
               if (fabs(z1) < 1.0 && fabs(z2) < 1.0 && fabs(z3) < 1.0)
               {
                  xp = (x1 + x2 + x3) / 3.0;
                  yp = (y1 + y2 + y3) / 3.0;
                  fprintf(ps, "%3.2f %3.2f %3.2f %3.2f %3.2f vseg\n",
                          (1.0 + x1) * 250.0, (1.0 + y1) * 250.0,
                          (1.0 + x2) * 250.0, (1.0 + y2) * 250.0, 1.5);
                  fprintf(ps, "%3.2f %3.2f %3.2f %3.2f %3.2f vseg\n",
                          (1.0 + x2) * 250.0, (1.0 + y2) * 250.0,
                          (1.0 + x3) * 250.0, (1.0 + y3) * 250.0, 1.5);
                  fprintf(ps, "%3.2f %3.2f moveto\n",
                          (1.0 + xp) * 250.0, (1.0 + yp) * 250.0);
                  fprintf(ps, "(%3.2f) show\n", ps_measures[i]);
                  fprintf(ps, "[] 0 setdash\n");
               }
            }
         }

     /*---LABELING---*/

         if (activeLabelCount)
         {
            fprintf(ps, "%d setgray\n", wob);
            for (i = 0; i < activeLabelCount; i++)
            {
               atm = showlist[i][0];
               if ((showlist[i][1] & 2) && tx[atm][2] < 1.0 &&
                   tx[atm][2] > (-1.0))
               {
                  fprintf(ps, "%3.2f %3.2f moveto\n",
                          (1.0 + tx[atm][0]) * 250.0,
                          (1.0 + tx[atm][1]) * 250.0 - 15.0);
                  res = attribs[atm].res;
                  fprintf(ps, "((%s %d)) show\n", rlab[res], res);
               }
               if ((showlist[i][1] & 1) && tx[atm][2] < 1.0 &&
                   tx[atm][2] > (-1.0))
               {
                  fprintf(ps, "%3.2f %3.2f moveto\n",
                          (1.0 + tx[atm][0]) * 250.0,
                          (1.0 + tx[atm][1]) * 250.0);
                  fprintf(ps, "(%s) show\n", attribs[atm].label);
               }
            }
         }

      }

      /* ---- PostScript EPILOG ---- */

      fprintf(ps, "stroke\n");
      fprintf(ps, "showpage\n");
      fclose(ps);
   }

   /* --- BACK TO GRAPHICS --- */
#ifdef GL_STUFF
   tpoff();
#endif
   free(segx);
   free(tx);
   // redotitles(vp);
}
#endif // POSTSCRIPT_DUMP_REINTEGRATED


// These globals are shared by the drawimage_* routines
dispinfo view;               /* view information */
int drawimage_hasmultiframes;
float movieViewBoxSize; // coordinates of Movie Frame
// long xsiz;

// Pointer to a global area where a copy of the movie
// can be stored if needed - and subsequently modified through
// filtering or periodic boundary imaging
short* xxs;

// This is the global read-only frame manager that will give
// us our main, unfiltered, movie.
FrameManager* primaryFM;


/*------------------------------------------------------------*/
/*-------   Hydrogen Bond Routines   -------------------------*/
/*------------------------------------------------------------*/

void checkhb(void) /* cutoff distance (angstroms) */
{
   int i,ii, jj, ihbcut2;
   register int dx, dy, dz, dist2;
   const short* xx;

#ifdef LATER
   if (filterSize && filterFlag || oneIfBox)
      xxx = xxs + i * 3 * (nat + oneIfBox);
   else
#endif
      xx = FrameManagerGetFrame(primaryFM,dispuiMovieFrame()); // ,xx + i * 3 * (nat + oneIfBox);


   ihbcut2 = (int) (KAPPA * KAPPA * view.hbcut * view.hbcut);
   hbonds = 0;

   for (i = 0; i < masterhbonds; i++)
   {
      ii = masterhblist[i][0];
      jj = masterhblist[i][1];
      if (attribs[ii].color && attribs[jj].color)
      {
         dx = (int) (xx[ii * 3 - 3] - xx[jj * 3 - 3]);
         dy = (int) (xx[ii * 3 - 2] - xx[jj * 3 - 2]);
         dz = (int) (xx[ii * 3 - 1] - xx[jj * 3 - 1]);
         dist2 = dx * dx + dy * dy + dz * dz;
         if (dist2 < ihbcut2)
         {
            hblist[hbonds][0] = ii;
            hblist[hbonds][1] = jj;
            hbonds++;
         }
      }
   }

}

// These variables are optionally input by the user at startup
// and display in the bottom right corner of the movie as it
// progresses
static float picoSecondsInitialOffset, picoSecondsPerFrame;

// This function is called one time - by main() to load
// the atom coordinates at each frame and set
// global variables about the dimension, etc, of the movie.
void drawimage_prepare(void)
{
   float minx = 99999.9, maxx = (-99999.9),
     miny = 99999.9;
   float maxy = (-99999.9),
     minz = 99999.9,
     maxz = (-99999.9);
   int i, ii;
   char instr[82];
   FILE* fdes;
   int bytesperframe, nf, cv, nnat, bpf;
   int header[5];
   size_t nbytes;

#ifdef SURFACE_FILE_BROUGHT_BACK

   fdes = fopen(surfaceFileName, "rb");
   if (fdes != 0)
   {
      fread_or_exit_if_error(surface, (unsigned int) (3 * MAXSURF),1,fdes);
      fclose_or_exit_if_error(surfaceFileName,fdes);
      // surf = 1;
         /* surf is flag which says that surface data
                                 * is available */
      surfonoff = 1;            /* surfonoff is flag for whether display on
                                 * or off */
      for (i = 0; i < MAXATM; i++)
         surfatm[i] = 1;        /* initialize all atoms on */
   }
#endif

   if (frames == 9999)
      printf("--loading all frames--\n");
   else
      printf("--loading %d frames--\n", frames);

   // Open the binary movie frame which preproc created
   fdes = fopen_or_exit_if_error(coordFileName, "rb");
   bytesperframe = 3 * (nat + oneIfBox) * sizeof (short);
   /* read header record */

   assert(sizeof(header) == 20); // <--- A portability check...

   fread_or_exit_if_error(header, 20,1,fdes);
   cv = header[0];
   nnat = header[1];
   bpf = header[2];
   nf = header[3];
   /* Check file version, if not 101 assume old format */
   if (cv == 101 && header[4] == -101)
   {
      printf("File version %d:  nat=%d  bpf=%d  #frames=%d\n", cv, nnat, bpf,
             nf);
      if (bpf != bytesperframe)
      {
         printf("Mismatch in number of bytes per frame\n\n");
         exit(2);
      }
      if (nnat < 0 || nnat > MAXATM)
      {
         printf("Unable to handle number of atoms (%d).  \nMAXATM set to %d in dispui.c\n\n", nnat,MAXATM);
         exit(2);
      }
      if (endFrame > nf - 1)
      {
         printf("Only %d frames available in file\n", nf);
         endFrame = nf - 1;
      }
      if (frames == 9999)
         frames = nf;
   }
   else
      {                            /* old file version */
      fclose_or_exit_if_error(coordFileName,fdes);
      fdes = fopen_or_exit_if_error(coordFileName, "rb");
      if (frames == 9999)
      {
#define MALDEF 10               /* default Mbytes of memory to allocate, */
      /* when proper amount cannot be calculated. */
         frames = (10000 * MALDEF / bytesperframe) * 100 + 1;
         if (frames > 2000)
            frames = 2000;
      }
   }

   nbytes = frames * bytesperframe;

// printf("Allocating 2*%d bytes of memory\n", nbytes);

   // This will read all the movie frames in
   primaryFM = FrameManagerConstruct(fdes,sizeof(header),startFrame,frames,bytesperframe);

   // This creates an area that accept a copy of all the atom coordinates of all frames
   // which can in turn be modified as needed (filtering, and box reimaging, for example).
   xxs = malloc_or_exit_if_error(nbytes,"drawimage - xxs");

   printf("%d bytes (%d frames) read.\n\n", bytesperframe*frames, frames);

   if (frames > 1)
      {
      SafeInputString("\nTime interval per frame (picoseconds) :",instr,sizeof(instr));
      sscanf(instr, "%f", &picoSecondsPerFrame);

      SafeInputString("Initial offset (picoseconds) :",instr,sizeof(instr));
      sscanf(instr, "%f", &picoSecondsInitialOffset);
      drawimage_hasmultiframes = 1;
      }
   else
      drawimage_hasmultiframes = 0;

   /* find extreme coordinates for sizing and scaling */

   {
   const short* xx = FrameManagerGetFrame(primaryFM,0);
   for (i = 0; i < nat; i++)
      {
      ii = i * 3;
      if ((float) xx[ii] < minx)
         minx = (float) xx[ii];
      if ((float) xx[ii] > maxx)
         maxx = (float) xx[ii];
      if ((float) xx[ii + 1] < miny)
         miny = (float) xx[ii + 1];
      if ((float) xx[ii + 1] > maxy)
         maxy = (float) xx[ii + 1];
      if ((float) xx[ii + 2] < minz)
         minz = (float) xx[ii + 2];
      if ((float) xx[ii + 2] > maxz)
         maxz = (float) xx[ii + 2];
      }
   }

   // The movieViewBoxSize will be passed to the
   // dispui.c init routines.  It is the cube which
   // dispui.c will manage in cartesian space.

   {
   double size = maxx - minx;
   if ((maxy - miny) > size)
      size = maxy - miny;
   if ((maxz - minz) > size)
      size = maxz - minz;
   movieViewBoxSize = 1.2 * size;
   }

   /* IF PERIODIC BOUNDS, PUT ATOMS IN BOX */
   if (oneIfBox)
   {
      cpcoords(xxs,primaryFM);
      if (isTruncatedOctahedron)
         wrap_to(xxs,0 /* filtsize */);
      else
         putinbox(xxs, 0 /*filtsize*/);
   }
}


void drawimage_initview(int recoverFlag)
{
   // Read in information about previous session
   // Or, set nice defaults for various OpenGL scene transforms
   
   FILE* fdes;
   int i;

   // If the user is attempting to recover a "crashed" session
   // then the session file will have a .xxx extension.
   if (recoverFlag)
      {
      i = strlen(sesFileName);
      strcpy(sesFileName+i,".xxx");
      }

   fdes = fopen(sesFileName, "rb");
   if (fdes != 0)
   {
      printf("\nReading prior session file %s",sesFileName);
      if (recoverFlag)
         printf(" from ungraceful shutdown");
      printf(".\n");
      fread_or_exit_if_error(&movie, (unsigned int)sizeof (movie),1,fdes);
      fread_or_exit_if_error(&view,  (unsigned int)sizeof (view),1,fdes);
      fclose(fdes);

      if (view.hbcut > 0.0 && frames == 1)
         checkhb();
   }
   else // Start an all new session - new clip planes, scaling, no translation
   {
      int i;
      movie.frontClip = 1.0;  // A Clipping from -1 x movieSize to 2 x movieSize
      movie.backClip =  2.0;
      movie.scaleFactor = 1.0; // All movie should be visible on screen at 1.0
      movie.xTranslation = movie.yTranslation = 0.0;
      movie.smear = 0;  // No frame smearing
      movie.imageView = nSingleImage; // No stereo or TRI mode for starters
      movie.hideLeftAndBottom = 0;    // Make UI components visible
      view.statictype = nStaticNone;  // No static image overlay

      view.nmeas = 0;                 // No measurements
      view.hbcut = 0.0;               // No H Bond display for starters

      // Set the rotation matrix to Identity - i.e. NO rotation.  
      // Yes - this is very tedious - but LoadIdentity does not seem
      // to be working until we do more openGL stuff downstream...
      for (i=0;i<16;i++)
         movie.rotationMatrix[i] = 0.0;

      movie.rotationMatrix[0] =
         movie.rotationMatrix[5] =
            movie.rotationMatrix[10] =
               movie.rotationMatrix[15] = 1.0;
   }

   if (recoverFlag) // Restore session file name to "correct" name
      sesFileName[i] = 0;
}

// Free memory when we shut down.
void drawimage_close(void)
{
   free(xxs);
}


/****************************************************************
 *
 * The code below is starting to deal with actually drawing the
 * movie frames - lines from atom to atom, and so forth
 *

 *------------------------------------------------------------*/
/* FAST MOLECULE DRAW - for clipping aid */
// We don't draw any colors - and just blast through the
// movedraw list.

void qkdrawmol(const short* xx)
{
   int i, atom_vertex_subscript;

   glColor3ub(100,100,100); // Simple grey
   glBegin(GL_LINE_STRIP);

   for (i = 0; i < drawindx; i++)
   {
      atom_vertex_subscript = drawlist_WhichAtom(i)*3-3;
      if (drawlist_IsMove(i) || ((i % 200) == 50))
      {
         glEnd();
         glBegin(GL_LINE_STRIP);
      }
      glVertex3sv(&xx[atom_vertex_subscript]);
   }

   glEnd();
}

/*---------------------------------------*/
/* DRAWS 3D CROSS FOR individual CATIONS */
/*---------------------------------------*/
void DrawCationCross(const short* v)
{
   short x, y, z;

   x = v[0];
   y = v[1];
   z = v[2];
   glBegin(GL_LINES); // Each pair of vertices will be connected by a line
   glVertex3i(x - CATSIZE, y, z);
   glVertex3i(x + CATSIZE, y, z);
   glVertex3i(x, y - CATSIZE, z);
   glVertex3i(x, y + CATSIZE, z);
   glVertex3i(x, y, z - CATSIZE);
   glVertex3i(x, y, z + CATSIZE);
   glEnd();
}

// Display the "static" backdrop (first or average) selected by the user!
// To speed this up, I use OpenGL display lists.

void PrepStaticMolecule(void)
{
   register int i, j, atm; // , ii;
   const short *x = 0;
   short* aves = 0;
   /* drawlist[i][0] == 0 if it's a move command */
   /* drawlist[i][0] == 1 if it's a draw command */
   /* drawlist[i][1] == atom to move or draw to */

   if (openGLDisplayListBaseIndex == 0) 
      openGLDisplayListBaseIndex = glGenLists(1);
   if (openGLDisplayListBaseIndex == 0) 
      {
      printf("displayRenderFrame: Unable to create Display Lists\n");
      exit(EXIT_FAILURE);
      }

   // Draw the molecule - but add to list as we go!
   glNewList(openGLDisplayListBaseIndex,GL_COMPILE);

   if (view.statictype == nStaticFirst)
      x = FrameManagerGetFrame(primaryFM,0);
   else if (view.statictype == nStaticAverage)
      {
      int* ave; 
      // aves is freed in another block - sorry.
      aves = malloc_or_exit_if_error(nat * 3 * sizeof(short),"DrawStaticMolecule- aves");

      // This is a temp var, freed in this block
      ave = malloc_or_exit_if_error(nat* 3 * sizeof(int),"DrawStaticMolecule - ave");
      memset(ave,0,nat*3*sizeof(int));

      for (i = 0; i < frames; i++)
         {
         const short* xx = FrameManagerGetFrame(primaryFM,i);
//       ii = i * (nat + oneIfBox) * 3;
         for (j = 0; j < nat * 3; j++)
            {
            ave[j] += xx[j];
            }
         }

      for (j = 0; j < nat * 3; j++)
         aves[j] = (short)(ave[j] / frames);
      x = aves;
      free(ave);
      }

   // Draw the static background molecule
   glColor3ub(180,180,180); // Simple grey
   glBegin(GL_LINE_STRIP);

   for (i = 0; i < drawindx; i++)
      {
      int atom_vertex_subscript = drawlist_WhichAtom(i)*3-3;
      if (drawlist_IsMove(i))
         {
         glEnd();
         glBegin(GL_LINE_STRIP);
         }
      glVertex3sv(&x[atom_vertex_subscript]);
      }
   glEnd();

   if (catindx)                 /*---CATIONS---*/
      for (i = 0; i < catindx; i++)
         {
         atm = catlist[i];
         DrawCationCross(&x[atm * 3 - 3]);
         }

   if (aves != 0)
      free(aves);

   glEndList();
}

void DrawStaticMolecule(void)
{
   // We only create a new display list if there is fundamental
   // change in the movie rendering strategy - or if this is a "first time"
   if (staticNeedsRebuild)
      {
      PrepStaticMolecule();
      staticNeedsRebuild = 0;
      }

   // Draw fast from display list created earlier.
   glCallList(openGLDisplayListBaseIndex);
}




/*-----------------------------------------------------------------*/
/* DRAWS MOLECULE WITH HALF-BOND COLORING, LABELING, ETC.  */
// We'll use the movedraw list information created in the most
// recent call to movedraw() and draw it all out.

void DrawMoleculeAndCations(const short xx[MAXATM * 3])
{
   extern int drawindx;
   int i,j;
   // register int j, bgn, end;
   int current_atom,current_atom_vertex_subscript;
   int previous_atom,previous_atom_vertex_subscript;

   // After a little experimenting with different platforms,
   // I found that the OpenGL controls here gave reasonable
   // anti-aliasing without a performance hit.  You are welcome
   // to tune these parameters to taste!

   glEnable( GL_LINE_SMOOTH);
   glEnable( GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

   // glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   glLineWidth(1.0); // ?REVIEW? DOCUMENT VAREFULLY - and #define AT TOP

   glShadeModel(GL_FLAT);
   glColor3b(255,255,255);

   glBegin(GL_LINE_STRIP);


   for (i = 0; i < drawindx; i++)
      {
      current_atom = drawlist_WhichAtom(i);
      current_atom_vertex_subscript = current_atom * 3 - 3;

      // If we have a move instruction in the drawlist, we
      // pick up the pen, start a new line at the atom coordinate,
      // and set the color for our new element.  Note that
      // we do not drawn anything in this little block. 
      if (drawlist_IsMove(i))
      {
      glEnd(); // Render any lines we have been drawing
      // Start a new line at the current atom position
      glBegin(GL_LINE_STRIP);
      glVertex3sv(&xx[current_atom_vertex_subscript]);
      // Set the color to the current atom color
      dispuiColorByPalette(attribs[current_atom].color);
      }
      else // It must be a "draw"
      {
      // Two different colors for these bonded atoms?
      if (attribs[current_atom].color != attribs[previous_atom].color)
         {
         short middle[3];
         
         // Use prior color to connect previous atom to midpoint
         for (j = 0; j < 3; j++)
            {
            middle[j] = (xx[previous_atom_vertex_subscript + j] + xx[current_atom_vertex_subscript + j]) / 2;
            }
         glVertex3sv(middle);

         // Use our atom color to finish the connection
         dispuiColorByPalette(attribs[current_atom].color);
         glVertex3sv(&xx[current_atom_vertex_subscript]);
         }
      else
         { // Both atoms have the same color - so simple draw from one to the other is fine!
         glVertex3sv(&xx[current_atom_vertex_subscript]);
         }
      } // End of IsMove or IsDraw..  We always come here
      previous_atom = current_atom;
      previous_atom_vertex_subscript = current_atom_vertex_subscript;
      } // End of for loop

   glEnd();


   /*---CATIONS---*/
   // If catindx is zero (no cations) nothing happens here...
   for (i = 0; i < catindx; i++)
      {
         int color;

         current_atom = catlist[i];
         color = attribs[current_atom].color;
         if (color)
            {
            dispuiColorByPalette(color);
            DrawCationCross(&xx[current_atom * 3 - 3]);
            }
      }


   glDisable(GL_BLEND);
   glDisable(GL_LINE_SMOOTH);


#ifdef SURFACE_FILE_BROUGHT_BACK
   /* ---- SURFACE ---- */
   if (surfonoff)
      for (i = 1; i <= nat; i++)
      {
         color = attribs[i].color;
         if (color && surfatm[i])
         {
            bgn = (int) (surface[2 * i] + 10000 * surface[2 * i + 1]);
            end = (int) (surface[2 * i + 2] + 10000 * surface[2 * i + 3]);
            basecolor = 256 * color + 255;
            colorl(basecolor - colorrange / 8, basecolor, 0, 0x7fff);
            bgnpoint();
            if (end - bgn > 250)
               for (j = bgn; j < end; j += 3)
               {
                  if ((j % 300) < 3)
                  {             /* make sure not too many points */
                     endpoint();
                     bgnpoint();
                  }
                  glVertex3sv(&surface[j]);
               }
            else
               for (j = bgn; j < end; j += 3)
                  glVertex3sv(&surface[j]);
            endpoint();
         }
      }

#endif
}

// When an atom is picked, we turn it into a cross - seems to look OK.
// Scaling is a bit odd sometimes - oh well.
void DrawPickedAtom(const short xx[MAXATM * 3])
{
   if (lastPickedAtom > 0)
      {                            /* HIGHLIGHT PICKED ATOM */
      glLineWidth(3.0);
      glColor3f(.7,.7,.7);          // Strongish Grey
      DrawCationCross(&xx[lastPickedAtom * 3 - 3]);
      glLineWidth(1.0);
      }
}

// You can't miss the pink lines for possible Hbonds... :)
void DrawHydrogenBonds(const short xx[MAXATM * 3])
{
   /* ----- HYDROGEN BONDS ----- */

   if (hbonds)
      {
      int i;

      glLineWidth(1);
      glColor3ub(255,128,255);   // Pink
      glEnable(GL_LINE_SMOOTH);
      glBegin(GL_LINES);
      for (i = 0; i < hbonds; i++)
         {
         int atom1 = hblist[i][0] * 3 - 3;
         int atom2 = hblist[i][1] * 3 - 3;
         glVertex3sv(&xx[atom1]);
         glVertex3sv(&xx[atom2]);
         }
      glEnd();
      glDisable(GL_LINE_SMOOTH);
      }
}


// Draw the text atom labels that are active
void DrawLabels(const short xx[MAXATM * 3])
{
   int i;
   char tempstr[80];

   if (activeLabelCount)
   {
      void* font = GLUT_BITMAP_HELVETICA_12;

      // If the window is small, go to smaller font
      if (glutGet(GLUT_WINDOW_WIDTH) < 300) // This aesthetic decision was from the old display.. sort of
         font = GLUT_BITMAP_HELVETICA_10;

      glColor3ub(255,255,0); // Yellow Letters for now...

      for (i = 0; i < activeLabelCount; i++)
      {
         const char* sptr  = 0;

         int atm = showlist[i][0];

         // If we are flagged to show a residue label, show that
         if (showlist[i][1] & 2)
            {
            int res = attribs[atm].res;
            sprintf(tempstr, "(%s %d)", rlab[res], res);
            sptr = tempstr;
            }

         // If we are doing an atom only label, so be it.
         if (showlist[i][1] & 1)
            {
            sptr = attribs[atm].label;
            }

         // Position the graphics "cursor" on the 3D screen and
         // write out text at the right location
         if (sptr)
            {
            glRasterPos3s(xx[atm * 3 - 3], xx[atm * 3 - 2], xx[atm * 3 - 1]);
            dispuiDrawStringAtCurrentRasterPosition(font,sptr);
            }
      }
   }
}

// Render Measurements (distances, angles, dihedrals) with a dashed line
// and the measurement distance, or theta, next to the dashed line.
void DrawMeasurementLines(const short xx[],const dispinfo* vp)
{
   int i;
   int a1, a2, a3, a4;
   void* font = GLUT_BITMAP_HELVETICA_12;

   if (glutGet(GLUT_WINDOW_WIDTH) < 300) // This aesthetic decision was from the old display.. sort of
      font = GLUT_BITMAP_HELVETICA_10;

   glEnable( GL_LINE_SMOOTH);
   glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);

   for (i = 0; i < vp->nmeas; i++)
   {                            /* INDICATE MEASUREMENTS */
      a1 = vp->measarry[i][0] * 3 - 3;
      a2 = vp->measarry[i][1] * 3 - 3;
      a3 = vp->measarry[i][2] * 3 - 3;
      a4 = vp->measarry[i][3] * 3 - 3;
      if (a3 < 0)  // Then this is simply a distance
         {
         float af[3];
         char distance_str[20];

         /* DOTTED LINE FOR DISTANCE */
         glEnable(GL_LINE_STIPPLE);
         glLineStipple(1,0x4444);
         // Make it a bit fatter than bonds
         glLineWidth(1.3);

         glColor3ub(255,255,0); //  Yellow

         glBegin(GL_LINES);
            glVertex3sv(&xx[a1]);
            glVertex3sv(&xx[a2]);
         glEnd();

         glDisable(GL_LINE_STIPPLE);

         // Now add text of measured length
         glColor3ub(255,255,150); //  Yellow Letters too
         // Find mid point of dotted line
         glRasterPos3i(
            (xx[a1]+xx[a2])/2,
            (xx[a1+1]+xx[a2+1])/2,
            (xx[a1+2]+xx[a2+2])/2);

         // Now write the distance.
         af[0] = (float) (xx[a1] - xx[a2]) / KAPPA;
         af[1] = (float) (xx[a1 + 1] - xx[a2 + 1]) / KAPPA;
         af[2] = (float) (xx[a1 + 2] - xx[a2 + 2]) / KAPPA;

         sprintf(distance_str," %1.3f",VectorNorm(af));
         dispuiDrawStringAtCurrentRasterPosition(font,distance_str);
         }
      else
         {                         /* OUTLINE ANGLES & DIHEDRALS */
         glLineWidth(4.0);
         glColor3ub(255,255,0); // Yellow Letters for now...
         glMatrixMode(GL_MODELVIEW);
         /* DRAW LINE TWICE FOR VISIBILITY */
         glPushMatrix();
         // First go back into display
         glTranslatef(0.0, 0.0, 2.5);

         glBegin(GL_LINE_STRIP);
            glVertex3sv(&xx[a1]);
            glVertex3sv(&xx[a2]);
            glVertex3sv(&xx[a3]);
            if (a4 >= 0)               // Add 4th point if dihedral
               glVertex3sv(&xx[a4]);
         glEnd();
         // Now come out of the display
         glTranslatef(0.0, 0.0, -5.0);
         glBegin(GL_LINE_STRIP);
            glVertex3sv(&xx[a1]);
            glVertex3sv(&xx[a2]);
            glVertex3sv(&xx[a3]);
            if (a4 >= 0)   // Add 4th point if dihedral
               glVertex3sv(&xx[a4]);
         glEnd();
         // Restore transformation matrices
         glPopMatrix();
         }
      glLineWidth(1.0);
   }

   glDisable( GL_LINE_SMOOTH);
}

// This is global because when the mouse is clicked we need to know
// if it was clicked in the measureBox
static int measureBoxWidth = 350;
static const int textLineHeight = 18;

int MeasurementTextYpos(int measureSubscript,int movieWindowHeight)
{
   return movieWindowHeight - 5 - textLineHeight - textLineHeight * measureSubscript;
}

// Draw text in the upper right corner of the screen.
void DrawMeasurementText(
            const short xx[MAXATM * 3],
            const dispinfo* vp,
            int movieWindowWidth,
            int movieWindowHeight)
{
   static const int distanceWidth = 250;
   static const int angleWidth = 300;
   static const int dihedralWidth = 350;

   void* textFont = GLUT_BITMAP_HELVETICA_12;

   int i,y;
   int m,j;

   float af[3],bf[3],cf[3];
   float dist,theta;

   char outstr[50];

   measureBoxWidth = distanceWidth; // Try to keep the black box small...
   for (i=0;i<vp->nmeas;i++)
      {
      int a3 = 3 * (vp->measarry[i][2]) - 3;
      int a4 = 3 * (vp->measarry[i][3]) - 3;
      if (a3 < 0) // It's another distance
         continue;
      else if (a4 < 0) // It's an angle
         {
         if (measureBoxWidth <angleWidth)
            measureBoxWidth = angleWidth;
         }
      else // It's a dihedral
         {
         measureBoxWidth = dihedralWidth; // This is as wide as it gets
         break;
         }
      }

   // Don't display text if none to display - or not enough window width!
   if ((! vp->nmeas) || (measureBoxWidth > movieWindowWidth))
      return;

// Originally, I blacked out the top right of the screen
// But users wanted atoms to display there - so I said "OK"      
// Let's have a kind of "Black Box" at the top right of the screen.
//        if (vp->nmeas)
//        {
//                glColor3ub(0,0,0); // Black Box...
//                glRectf(movieWindowWidth-measureBoxWidth, movieWindowHeight - 10 - textLineHeight * vp->nmeas, movieWindowWidth-1, movieWindowHeight-1);
//        }

   // Iterate over all active measurements and display
   // the text for each in the top right of the screen.
   for (i = 0; i < vp->nmeas; i++)
      {
      int a1 = 3 * (vp->measarry[i][0]) - 3;
      int a2 = 3 * (vp->measarry[i][1]) - 3;
      int a3 = 3 * (vp->measarry[i][2]) - 3;
      int a4 = 3 * (vp->measarry[i][3]) - 3;

      y = MeasurementTextYpos(i,movieWindowHeight);

      if (a3 < 0) /* it's got 2 components - so it's a DISTANCE */
         {
         glColor3ub(255,0,0); // color(RED);
         af[0] = (float) (xx[a1] - xx[a2]) / KAPPA;
         af[1] = (float) (xx[a1 + 1] - xx[a2 + 1]) / KAPPA;
         af[2] = (float) (xx[a1 + 2] - xx[a2 + 2]) / KAPPA;
         dist = VectorNorm(af);
         // ps_measures[i] = dist; // used for ps output - not wired in
         a1 = a1 / 3 + 1;
         a2 = a2 / 3 + 1;
         sprintf(outstr, "(Dist) :%d@%s & :%d@%s", attribs[a1].res,
                 attribs[a1].label, attribs[a2].res, attribs[a2].label);

         DrawString(movieWindowWidth-distanceWidth,y,textFont,outstr);

         sprintf(outstr, "%7.3f", dist);
         DrawString(movieWindowWidth-StringWidth(textFont,outstr)-5,y,textFont,outstr);
         }
      else if (a4 < 0) // it's got 3 components - so it's a bond ANGLE
         {
         int res;
         glColor3fv(colorAngle);
         for (m = 0; m < 3; m++)
            {
            af[m] = (float) xx[a1 + m] / KAPPA;
            bf[m] = (float) xx[a2 + m] / KAPPA;
            cf[m] = (float) xx[a3 + m] / KAPPA;
            }
         theta = VectorAngle(af, bf, cf);
         // ps_measures[i] = theta; // used for ps output - not wired in
         DrawString(movieWindowWidth-angleWidth,y,textFont,"(Angle)");
         res = 0;
         for (j = 0; j < 3; j++)
            {
            int atm = vp->measarry[i][j];
            if (res != attribs[atm].res)
               {
               res = attribs[atm].res;
               sprintf(outstr, " :%d", res);
               dispuiDrawStringAtCurrentRasterPosition(textFont,outstr);
               }
            sprintf(outstr, "@%s", attribs[atm].label);
            for (m = 0; outstr[m] > 32; m++) ;
               outstr[m] = 0;
            dispuiDrawStringAtCurrentRasterPosition(textFont,outstr);
            }
         sprintf(outstr, "%7.2f", theta);
         DrawString(movieWindowWidth-StringWidth(textFont,outstr)-5,y,textFont,outstr);
         }
      else // This measurement has 4 components - must be a DIHEDRAL
         {                         
         const short* a = xx + a1;
         const short* b = xx + a2;
         const short* c = xx + a3;
         const short* d = xx + a4;
         int res = 0;

         theta = tors(a, b, c, d);
         // ps_measures[i] = theta;   // used for ps output - not wired in

         glColor3fv(colorDihedral);
         DrawString(movieWindowWidth-dihedralWidth,y,textFont,"(Dihed)");
         for (j = 0; j < 4; j++)
            {
            int atm = vp->measarry[i][j];
            if (res != attribs[atm].res)
               {
               res = attribs[atm].res;
               sprintf(outstr, " :%d", res);
               dispuiDrawStringAtCurrentRasterPosition(textFont,outstr);
               }
            sprintf(outstr, "@%s", attribs[atm].label);
            for (m = 0; outstr[m] > 32; m++) ;
            outstr[m] = 0;
            dispuiDrawStringAtCurrentRasterPosition(textFont,outstr);
            }
         sprintf(outstr, "%7.2f", theta);
         DrawString(movieWindowWidth-StringWidth(textFont,outstr)-5,y,textFont,outstr);
         }
      } // end of for (i = 0; i < vp->nmeas; i++)
}

// When the user clicks int eh movie Window, dispui.c calls this function
// We simply determine if the user has clicked on a measurement in the top
// right of the screen, and delete it if so.
int displayProcessRawMovieClick(int mouseX0,int mouseY0,int movieWindowWidth,int movieWindowHeight)
{
   if ( ( view.nmeas == 0 ) || (measureBoxWidth > movieWindowWidth))
      return 0; // Nothing to do...  Click nowehere near measurement text

   if ( mouseX0 > movieWindowWidth - measureBoxWidth )
      {
      int i;

      // This looks like "dumb code" in that I could have done a little
      // division to immediately determine which line is being deleted.
      // Oh well - this is very fast.
      for ( i=0;i<view.nmeas;i++ )
         {
         int ymin = MeasurementTextYpos(i,movieWindowHeight);

         if ( (mouseY0 >= ymin) && (mouseY0 <= ymin + textLineHeight) )
            {
            int j;
            // Delete the measurement
            view.nmeas--;
            for ( j=i;j<view.nmeas;j++ )
               {
               memcpy(view.measarry[j],
                      view.measarry[j+1],
                      sizeof(view.measarry[0]));

               }
            memset(view.measarry[view.nmeas],0,sizeof(view.measarry[0]));
            dispuiStatusMessage[0] = 0;
            MovieAdvanceTimerKeyAndPostRedisplay();
            return 1; // Tell dispui.c calling func that we handled the Mouse Click!
            }
         }
      }

   return 0; // Tell dispui.c that click did not match up with anything..
}


// Called by dispui.c to render the clipping tool in bottom
// left of screen (if activated by user).
// This just calls qkdrawmol for now - but we might
// use DisplayLists or something else later.
void displayRenderClipToolFrame(int frame)
{
   const short* xxx = FrameManagerGetFrame(primaryFM,frame);

#ifdef NOT_SURE_WHY_THIS_FAILS_BUG_MAYBE_IMPLEMENT_LATER
   // If we have a display list, then use that for the clipping tool!
   if (openGLDisplayListBaseIndex != 0)
      {
      if (openGLFrameReady[frame])
         {
         glCallList(openGLDisplayListBaseIndex+frame);
         }
      else
         qkdrawmol(xxx);
      }
   else
#endif

   qkdrawmol(xxx);
}

// dispui.c calls this function when it needs us to draw the molecule
// frame tells us which frame to draw.  smearFlag lets us know
// If we are drawing a smeared image - and in this case - we won't
// draw labels (which just look bad when smeared).
void displayRenderFrameInPerspectiveBox(int frame,int smearFlag)
{
   const short* xxx;
   int i;


   // If we have an active filter, or we have periodic boundary conditions
   // then we want to get coordinates from the xxs global copy.
   // If NOT, then use the unmodified preproc coordinates directly.
   // In any case xxx[atoms] will be the coordinates to draw.
   i = (frame + 2 * (frames - filterSize)) % (frames - filterSize);
   if (filterSize && filterFlag || oneIfBox)
      xxx = xxs + i * 3 * (nat + oneIfBox);
   else
      xxx = FrameManagerGetFrame(primaryFM,i);

   // printf("*** frame=%d, i=%d %d %d %d\n",frame,i,drawimage_hasmultiframes,nat,oneIfBox);
   glEnable(GL_DEPTH_TEST);

   // If user wants a static image to compare current frame to, then
   // draw that first.
   if (view.statictype == nStaticFirst || view.statictype == nStaticAverage)
      DrawStaticMolecule();

#ifdef DISPLAY_LISTS_WORKING_AND_IN_USE
      if (openGLDisplayListBaseIndex != 0) // Then we have allocated a display list!
         {
         if (openGLFrameReady[frame])
            {

            // Draw lightning fast from display list created earlier.
            glCallList(openGLDisplayListBaseIndex+frame);
            }
         else
            {
            // Draw the molecule - but add to list as we go!
            glNewList(openGLDisplayListBaseIndex+frame,GL_COMPILE_AND_EXECUTE);
               DrawMoleculeAndCations(xxx);
            glEndList();

            openGLFrameReady[frame] = 1;
            }
         }
      else // No Display Lists - so just draw the molecule
#endif

   // Draw the entire molecule and cations.  This is where
   // we do the most drawing
   DrawMoleculeAndCations(xxx/*,&view,smearFlag*/);

   glDisable(GL_DEPTH_TEST);

   // Only draw the h bonds, measurements, picked atoms, and labels
   // if this is the first frame.  I.e., don't draw them on frames
   // 2.. n if the user is smearing.
   if (! smearFlag)
      {
      DrawHydrogenBonds(xxx);
      DrawMeasurementLines(xxx,&view);
      DrawPickedAtom(xxx);
      DrawLabels(xxx);
      }
}

// dispui.c calls this to add some 2D things.
// We use this call as an opportunity to draw text for measurements in the
// top right of the display.
void displayRenderFrameInOrthoBox(int frame,int movieWindowWidth,int movieWindowHeight)
{
   const short* xxx;

   // If we have an active filter, or we have periodic boundary conditions
   // then we want to get coordinates from the xxs global copy.
   // If NOT, then use the unmodified preproc coordinates directly.
   // In any case xxx[atoms] will be the coordinates to draw.
   int i = (frame + 2 * (frames - filterSize)) % (frames - filterSize);
   if (filterSize && filterFlag || oneIfBox)
      xxx = xxs + i * 3 * (nat + oneIfBox);
   else
      xxx = FrameManagerGetFrame(primaryFM,i);
   DrawMeasurementText(xxx,&view,movieWindowWidth,movieWindowHeight);
}


// For picking with the mode, simply draw a virtual cube around
// all the atoms that are not black.
// This code is _not_ displayed to the screen since this code
// runs in GL_SELECT mode.  But, it is used with the OpenGL "picking" code...
void displayRenderFrameInSelectMode(int frame)
{
   int atom;
   const short* xxx;

   // If we have an active filter, or we have periodic boundary conditions
   // then we want to get coordinates from the xxs global copy.
   // If NOT, then use the unmodified preproc coordinates directly.
   // In any case xxx[atoms] will be the coordinates to draw.
   int i = (frame + 2 * (frames - filterSize)) % (frames - filterSize);
   if (filterSize && filterFlag || oneIfBox)
      xxx = xxs + i * 3 * (nat + oneIfBox);
   else
      xxx = FrameManagerGetFrame(primaryFM,i);

   // Set color to white - We won't see it - it is irrelevant.
   // Note that we are not drawing bond lines - just these cubes.
   glColor3ub(255,255,255);
   for (atom=1;atom<= nat;atom++)
      {
      if (attribs[atom].color != 0) // If it is not "black" then draw the "rectangle" around the atom
         {
         GLint x = xxx[atom*3-3];
         GLint y = xxx[atom*3-2];
         GLint z = xxx[atom*3-1];

         glLoadName(atom);          // So - our Open GL "hit records" will consist of atom numbers...
      {
      int err = glGetError();
      if (err != GL_NO_ERROR)
         {
         printf("Terminated at atom %d (%d,%d,%d) due to OpenGL error %d %s\n",atom,x,y,z,err,gluErrorString(err));
         exit(EXIT_FAILURE);
         }
      }
         glBegin(GL_LINE_LOOP);
            glVertex3i(x-1,y-1,z);
            glVertex3i(x-1,y+1,z);
            glVertex3i(x+1,y+1,z);
            glVertex3i(x+1,y-1,z);
         glEnd();
         }
      }
}

void (* AtomPickedFunc)(void);

void AtomPickedFuncStart(void (* _AtomPickedFunc)(void))
{
   AtomPickedFunc = _AtomPickedFunc;

   // If we already had a selected atom - then go ahead and do whatever with it...
   if (lastPickedAtom > 0 && lastPickedNeedsProcessing)
      AtomPickedFunc();
}


// called by dispui.c after we have rendered the atoms in GL_SELECT
// mode.
void displayProcessHits(GLint hits,const GLuint* selectBuf)
{
   // This global lets us know the last atom we clicked on
   lastPickedAtom = 0;
   // Use this to tell program not to reprocess highlighted atom.
   lastPickedNeedsProcessing = 0;

   // We must sift through all atoms which OpenGL thinks are "near" the mouse
   // pointer... and we'll pick the one closest to the user as the "right" one
   // for furter processing as the selected atom.
   // If hits == 0, then the user clicked in the movie but OpenGL
   // thought nothing was nearby.
   
   if (hits > 0)
      {
      int i;
      const GLuint* hitRecord = &selectBuf[0];
      unsigned int bestZ = hitRecord[1]; // minZValue

      lastPickedAtom = hitRecord[3];   // Atom

      for (i=1;i<hits;i++)
         {
         const GLuint* hitRecord = &selectBuf[i*4];

         if (hitRecord[1] < bestZ) // If this hit is "closer" to us then it is better
            {
            bestZ = hitRecord[1];  // minZValue
            lastPickedAtom = hitRecord[3];   // Atom
            }

         }

      lastPickedNeedsProcessing = 1;
      MovieAdvanceTimerKeyAndPostRedisplay();

      // If we have an active "handler" function for dealing with atoms as they are picked, then call away!
      if (AtomPickedFunc != 0)
         AtomPickedFunc();


#ifdef TESTMODE_BUG_ETC
      // When debugging mouse clicks in the movie window, this code is
      // extremely helpful
      printf("Total hits %d\n",hits);
      for (i=0;i<hits;i++)
         {
         const GLuint* hitRecord = &selectBuf[i*4];
         unsigned int numNames = hitRecord[0];
         unsigned int minZValue = hitRecord[1];
         unsigned int maxZValue = hitRecord[2];
         unsigned int atom = hitRecord[3];

         printf("%3d On Stack: %u  Atom: %4u  Residue %d Z:%u-%u Label: %-.5s\n",
                i,
                           numNames,
                                   atom,
                                          attribs[atom].res,
                                             minZValue,maxZValue,

                                                   attribs[atom].label);
         }
#endif
      } // end hits > 0
}


/***************************************
 * Code to support Ramachandron Plot Window
 *
 * Taken from MD Display 2.0, and converted
   as directly as possible to OpenGL
 *
 */

// dofisi() calculates all phi and psi angles
// for all residies in the movie frame.

void dofisi(const short* xxx, // Coordinates this frame
            int c[],          // Backbone C atom subscripts
            int ca[],         // Backbone C-alpha atom subscripts
            int n[],          // Backbone N atom subscripts
            float phi[],      // Returned phi angles
            float psi[])      // Returnes psi angles
{
   int i;
   int c1, n1, ca1;
   int c2, n2, ca2;

   for (i = 1; i < nres; i++)
   {
      ca1 = ca[i];
      n1 = n[i];
      c1 = c[i];
      ca2 = ca[i + 1];
      n2 = n[i + 1];
      c2 = c[i + 1];
      
      // Note that we can't compute at ends
      // due to lack of data.
      if (xxx[c1 * 3 - 3] != 0 || xxx[c1 * 3 - 2] != 0 ||
          xxx[c1 * 3 - 1] != 0)
         {
         // psi = torsion from N-Ca-C-next N
         psi[i] =
            tors(&xxx[n1 * 3 - 3], &xxx[ca1 * 3 - 3], &xxx[c1 * 3 - 3],
                 &xxx[n2 * 3 - 3]);
                 
         // phi of next residue = torsion from C-nextN-nextCa-nextC
         phi[i + 1] =
            tors(&xxx[c1 * 3 - 3], &xxx[n2 * 3 - 3], &xxx[ca2 * 3 - 3],
                 &xxx[c2 * 3 - 3]);
         }
      else
         {
         psi[i] = 0.0;
         phi[i + 1] = 0.0;
         psi[i - 1] = 0.0;
         phi[i] = 0.0;
         }
   }
}

// Given phi or psi, determine whether the residue is likely part of
// an Alpha Helix (A) Beta Sheet (B) or in an illegal conformation (*)
// Exempt Glycine from this :)

void findab(const float phi[], const float psi[], char* code)
{
   int gly, i;
   float ph, ps;

   code[0] = '-';
   code[nres + 1] = '-';
   code[nres + 2] = '-';
   code[nres + 3] = '-';
   for (i = 1; i < nres + 1; i++)
   {
      code[i] = ' ';
      ph = phi[i];
      ps = psi[i];
      // If this is a glycine, we have lots of flexibility - don't
      // show "disallowed" unless we are way off.
      if (rlab[i][0] == 'G' && rlab[i][2] == 'Y')
         gly = 1;
      else
         gly = 0;
      // Beta sheet?
      if (ph < (-45) && ps > 100)
         code[i] = 'B';
      // Alpha helix?
      if (ph < (-40) && ph > (-75) && ps > (-65) && ps < (-30))
         code[i] = 'A';
      // Wider range is still OK - but not too far outside of this
      // Intervales which are not allowed are:
      // ps in [-155,-70]
      // ph in [80,180]
      // ps in [-180,0] and ph in [20,180]
      // ps in [115,180] and ph in [-20,180]
      // ph in [-20,20]
      if ((ps < (-70.0) && ps > (-155.0) || ph > 80.0 || ps < 0.0 &&
           ph > (-20.0) || ps > 115.0 && ph > (-20) || ph > (-20.0) &&
           ph < 20.0) && !gly)
         code[i] = '*';         /* in disallowed region */
      if (ph == 0.0 && ps == 0.0)   /* non-existant */
         code[i] = '-';
   }
}

/*=================== */
/* Draw THE RAMA PLOT */
/*=================== */

// The following coordinate system allows axes to be displayed nicely - assuming
// a reasonable size Rama Window

#define nRamaLeft (-240.0)
#define nRamaRight 200.0
#define nRamaBottom (-225.0)
#define nRamaTop 190.0

void DisplayDrawRamaPlot(
   const float phi[], // computed in findab
   const float psi[], // computed in findab
   const char code[], // computed in findab
   int ca[]) // C-alpha subscript - Determines whether res is colored 
{
   int gly, angle;
   int residue;
   static char str[10];
   long width;
   float dotWidth;

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   // Impose our desired coordiante system on the Rama Window, whatever
   // it's size
   gluOrtho2D(nRamaLeft,nRamaRight,nRamaBottom,nRamaTop);

   width = glutGet(GLUT_WINDOW_WIDTH);

   dotWidth = (float)width / (nRamaRight - (nRamaLeft));

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   {
   GLfloat saveLineWidth;
   glGetFloatv(GL_LINE_WIDTH,&saveLineWidth);

   glLineWidth(3.0);

   // Frame the entire RAMA Window in a square
   glColor3f(1.0,1.0,1.0); // REVIEW - colors should be in a .h file
   glBegin(GL_LINE_STRIP);
      glVertex2i(-180,180);
      glVertex2i(-180, -180);
      glVertex2i(180, -180);
      glVertex2i(180, 180);
      glVertex2i(-180, 180);
   glEnd();

   glLineWidth(saveLineWidth);

   // REVIEW - colors should be in a .h file
   glColor3f(1.0,1.0,0.0);

   // Draw the axes legends
   // Draw out the number of degrees of the angles.  Ticks first and then number legend
   for (angle = -180; angle < 181; angle += 45)
      {
      glBegin(GL_LINES);
         glVertex2i(angle, -180);
         glVertex2i(angle, -195);

         glVertex2i(-180, angle);
         glVertex2i(-190, angle);
      glEnd();

      sprintf(str, "%d", angle);
      DrawString(angle - (1.0/dotWidth)*(StringWidth(GLUT_BITMAP_HELVETICA_10,str)/2),
                  -206,
                 GLUT_BITMAP_HELVETICA_10,
                 str);

      sprintf(str, "%d", angle);
      DrawString(-190 - (1.0/dotWidth)*StringWidth(GLUT_BITMAP_HELVETICA_10,str) - 3,
                  angle-4,
                 GLUT_BITMAP_HELVETICA_10,
                 str);
      }

   DrawString(-234,
               15,
              GLUT_BITMAP_HELVETICA_10,
              "PSI");

   DrawString(-15,
               -219,
              GLUT_BITMAP_HELVETICA_10,
              "PHI");


   /* draw approximate guide-lines in gray */
   // Marks outside the guideliness will
   // will clue the user to residue conformations
   // that may be in trouble.

   // REVIEW - colors should be in the .h file
   glColor3f(0.5,0.5,0.5); 

   glBegin(GL_LINE_STRIP);
      glVertex2i(-40, 180);
      glVertex2i(-35, 150);
      glVertex2i(-35, 90);
      glVertex2i(-75, 40);
      glVertex2i(-35, -30);
      glVertex2i(-35, -70);
      glVertex2i(-170, -70);
      glVertex2i(-170, -25);
      glVertex2i(-145, -35);
      glVertex2i(-120, -20);
      glVertex2i(-170, 40);
      glVertex2i(-170, 180);
   glEnd();

   glBegin(GL_LINE_STRIP);
      glVertex2i(-170, -180);
      glVertex2i(-170, -170);
      glVertex2i(-45, -170);
      glVertex2i(-40, -180);
   glEnd();

   glBegin(GL_LINE_STRIP);
      glVertex2i(65, 5);
      glVertex2i(65, 100);
      glVertex2i(35, 70);
      glVertex2i(35, 35);
      glVertex2i(65, 5);
   glEnd();

   glBegin(GL_LINE_STRIP);
      /* draw alpha - beta guidelines */
      glVertex2i(-70, 180);
      glVertex2i(-45, 140);
      glVertex2i(-45, 100);
      glVertex2i(-150, 100);
      glVertex2i(-150, 180);
   glEnd();

   glBegin(GL_LINE_STRIP);
      glVertex2i(-45, -55);
      glVertex2i(-45, -35);
      glVertex2i(-110, -35);
      glVertex2i(-140, -45);
      glVertex2i(-140, -55);
      glVertex2i(-45, -55);
   glEnd();

   }

   // We can and will now draw marks for the conformation of each residue
   for (residue = 1; residue <= nres; residue++)
   {
      // If a residue had both phi and psi angles AND is colored (vs black)
      // then we want to mark it on the RAMA plot window.
      if (code[residue] != '-' && code[residue - 1] != '-' && attribs[ca[residue]].color)
      {
         float x = (float) phi[residue];
         float y = (float) psi[residue];
         if (code[residue] == 'A')
            glColor3f(1.0,0.0,1.0); // Strong RED fo res in alpha helix
         else if (code[residue] == 'B') // Blue for res in Beta Sheet
            glColor3fv(colorMyBlue);
         else glColor3f(1.0,1.0,1.0); // White for others

         if (rlab[residue][0] == 'G' && rlab[residue][2] == 'Y')
            gly = 1;
         else
            gly = 0;

         /// Biggish boxes for glycines
         if (gly)
            {
            glBegin(GL_LINES);
               glVertex2f(x + 2.0, y + 2.0);
               glVertex2f(x - 2.0, y - 2.0);
               glVertex2f(x - 2.0, y + 2.0);
               glVertex2f(x + 2.0, y - 2.0);
            glEnd();
            }
         else // Little boxes for other residues
            glRectf(x - 0.5, y - 0.5, x + 0.5, y + 0.5);

         /* LABEL THE BAD GUYS - i.e. residues with very suspect conformations*/
         if (code[residue] == '*')
            {
            glColor3f(1.0,1.0,1.0);
            sprintf(str, "%d", residue);
            DrawString(x + 0.8, y + 0.8, GLUT_BITMAP_HELVETICA_10,str);
            }
      }
   }
}

// When the user clicks in the Rama window, we want to know which residue
// corresponds to the mouse click.  Other code will take it from there and
// do labelling or whatever the user is looking for.
int DisplayComputeClosestRamaResidue(
              const float phi[],
              const float psi[],
              const char code[],
              int ca[],
              int mouseX,int mouseY)
{
   long width = glutGet(GLUT_WINDOW_WIDTH);
   long height = glutGet(GLUT_WINDOW_HEIGHT);

   // scale the Raw GLUT mouse clicks into our coordinate system

      float scaledMouseX = (float)mouseX/(float) width*(nRamaRight-nRamaLeft)+nRamaLeft;
   float scaledMouseY = (float)(height-mouseY-1)/(float) height*(nRamaTop-nRamaBottom)+nRamaBottom;

   // Seems to me that we ought to be able to pick within 15...
   float curdistance = 15.0;
   int curres = 0;

   int residue;
   // printf("%d %d %f %f\n",mouseX,mouseY,scaledMouseX,scaledMouseY);

   for (residue = 1; residue <= nres; residue++)
      {
      if (code[residue] != '-' && code[residue - 1] != '-' && attribs[ca[residue]].color)
         {
         float x = (float) phi[residue];
         float y = (float) psi[residue];

         float squareDistToMouse = (x-scaledMouseX)*(x-scaledMouseX)+(y-scaledMouseY)*(y-scaledMouseY);

         if (squareDistToMouse < curdistance)
            {
            curres = residue;
            curdistance= squareDistToMouse;
            }
         }
      }

   return curres;
}


// Called by code below - to draw all the hash marks or dots for each
// residue in the Rama window.
void DisplayRenderOrSelectRamaFrame(int frame,int selectMode,int mouseX,int mouseY)
{
   int maxres = nres + 10;
   int* c =  malloc_or_exit_if_error(maxres*sizeof(int),"displayRenderRama - c");
   int* ca = malloc_or_exit_if_error(maxres*sizeof(int),"displayRenderRama - ca");
   int* nn = malloc_or_exit_if_error(maxres*sizeof(int),"displayRenderRama - nn");

   float* phi = malloc_or_exit_if_error(maxres*sizeof(float),"displayRenderRama - phi");
   float* psi = malloc_or_exit_if_error(maxres*sizeof(float),"displayRenderRama - psi");

   char* code = malloc_or_exit_if_error(maxres*sizeof(char),"displayRenderRama - code");

   const short* xxx;
   int i; 

   memset(c,0,maxres*sizeof(*c));
   memset(ca,0,maxres*sizeof(*ca));
   memset(nn,0,maxres*sizeof(*nn));

   // REVIEW THIS_USED_TO_BE_LESS_HEAVY_HANDED_NEED_EVERY_TIME???
   /* CONDENSE BACKBONE COORDS (FOR RAMA ROUTINES) */
   for (i = 1; i <= nat; i++)
      {
      if (attribs[i].label[0] == 'C' && (attribs[i].label[1] == ' ' ||
                                         attribs[i].label[1] == 0))
         c[attribs[i].res] = i;
      if (attribs[i].label[0] == 'C' && attribs[i].label[1] == 'A' &&
          (attribs[i].label[2] == ' ' || attribs[i].label[2] == 0))
         ca[attribs[i].res] = i;
      if (attribs[i].label[0] == 'N' && (attribs[i].label[1] == ' ' ||
                                         attribs[i].label[1] == 0))
         nn[attribs[i].res] = i;
      }


   // REVIEW - Let's get this "Frame" and "i" thing settled once and for all!
   if (filterSize && filterFlag || oneIfBox)
      xxx = xxs + frame * 3 * (nat + oneIfBox);
   else
      xxx = FrameManagerGetFrame(primaryFM,frame);

   dofisi(xxx,c,ca,nn,phi,psi);
   findab(phi,psi,code);

   if (selectMode) // Then what residue is closet to Mouse Click?
      {
      int pickedRamaResidue =
         DisplayComputeClosestRamaResidue(phi,psi,code,ca,mouseX,mouseY);

      if (pickedRamaResidue)
         {
         char newWindowTitle[100];
         char keyString[100];
         sprintf(newWindowTitle,"Res %s%d   Press <R> to close.",rlab[pickedRamaResidue],pickedRamaResidue);
         glutSetWindowTitle(newWindowTitle);

         // If the keyboard handler is active - then fake entry of :nnn
         // should be very cool indeed!
         sprintf(keyString,":%d\r",pickedRamaResidue);
         KeyboardStuff(keyString); // This will in turn call code to highlight the residue or whatever...
         }
      else
         glutSetWindowTitle("No residue selected.   Press <R> to close.");
      }
   else // Just draw the Rama plot
      {
      DisplayDrawRamaPlot(phi,psi,code,ca);
      }

   free(code);
   free(psi);
   free(phi);
   free(nn);
   free(ca);
   free(c);
}

// This is a funciton called by dispui.c when it needs
// to render the Rama window
void displayRenderRamaFrame(int frame)
{
   DisplayRenderOrSelectRamaFrame(frame,0,-1,-1);
}

// This is a funciton called by dispui.c when there
// has been a mouse clieck in the Rama Window
void displaySelectRamaResidue(int frame,int mouseX,int mouseY)
{
   DisplayRenderOrSelectRamaFrame(frame,1,mouseX,mouseY);
}

// Banner for startup and exit
void DisplayStdoutCredits(void)
{
   printf("\n");
   printf("___________________________ MD Display Version %s ___________________________\n\n",VERSION);

   printf("           Chris Moth, Tim Callahan, Eric Swanson, and Terry Lybrand\n");
   printf("                  Copyright (C) 2002 by C. Moth and T. Lybrand\n\n");
   printf("                       Vanderbilt Universtiy, Nashville TN \n");
   printf("               http://www.structbio.vanderbilt.edu/~cmoth/mddisplay \n");
   printf("______________________________________________________________________________\n");
}

/*------------------------------------------------------------*/
/* SAVE SESSION (COLORING, VIEWING TRANSFORMATION, ETC.) */

void saveit(void) // so the user can resume later, where they left off.
{
   FILE* f;
   // int stat;

   if (! displayGracefulExit)
      {
      strcat(attribFileName,".xxx");
      strcat(sesFileName,".xxx");
      }

   // Write the attributes, which have the colors and labels of the atoms
   printf("\nWriting attribute file %s (%d bytes)\n", attribFileName,attribytesread);
   f = fopen_or_exit_if_error(attribFileName,"wb");
   fwrite_or_exit_if_error(attribs,(unsigned) attribytesread,1,f);
   fclose_or_exit_if_error(attribFileName,f);

   // Write the movie and view structs, which contain transformation matrices,
   // as well as other info such as scaling
   printf("   Writing session file %s (%d bytes).\n", sesFileName,sizeof(movie) + sizeof (view));
   f = fopen_or_exit_if_error(sesFileName,"wb");
   fwrite_or_exit_if_error(&movie,sizeof (movie),1,f);
   fwrite_or_exit_if_error(&view,sizeof (view),1,f);
   fclose_or_exit_if_error(attribFileName,f);
}


// Unfortunately, because of the way glut shuts down, the text displayed by "bye" does not appear until
// after it returns - when the glut graphics window is finally closed...

void bye(void) // Called via ANSI C "atexit" system.
{
   if (displayGracefulExit) // For graceful exit, clear console pretty well
      {
      int i;
      for (i=0;i<15;i++)
         printf("\n");
      DisplayStdoutCredits();
      
      printf("\n\nMD Display exiting....\n\n");
      // Did the user already answer 'y' to save session?
      if (displaySaveCurrentSession)
         saveit();
      }
   else
      {
      printf("\n\n");
      printf("---------------------------------------------------------------------------\n");
      printf("You have ungracefully exited MD Display.  In future, use the <Q>uit option.\n");
      printf("Your current session will be saved in *.xxx files and may be recoved with\n");
      printf("the -r option next time you run mddisplay.\n");
      printf("---------------------------------------------------------------------------\n");
      saveit();
      }
}


#if SPACEBALL == 1
/* Spaceball stuff - not supported in MDDisplay version 3.0 */

spcbal_do(rotvec, rotscale, transvec, transcale)
                    Coord rotvec[3];
                    float rotscale;
                    Coord transvec[3];
                    float transcale;
{

   sb_trx = transvec[0] * transcale / 3.0;
   sb_try = transvec[1] * transcale / 3.0;
   sb_trz = 0;
   sb_rtx = rotvec[0] * rotscale * 10.0;
   sb_rty = rotvec[1] * rotscale * 10.0;
   sb_rtz = -rotvec[2] * rotscale * 10.0;

}

spchome()
{
   return;
}

void SpaceBallRotate(dispinfo* vp)
{

   pushmatrix();
   loadmatrix(idmat);
   rot(sb_rtz, 'z');
   rot(sb_rty, 'y');
   rot(sb_rtx, 'x');
   multmatrix(vp->rotmat);
   getmatrix(vp->rotmat);
   popmatrix();
   sb_rtx = 0;
   sb_rty = 0;
   sb_rtz = 0;
}
#endif

// FUNCTION TO RETURN CORRECT COLOR FOR TYPE OF ATOM 
// We've checked and those little Oxygen atoms are
// indeed red :)
int clrbyatm(const char* str)
{
   int clr;

   switch (*str)
      {
      case 'C':   // Carbon
         clr = nGreen;
         break;
      case 'O':   // Oxygen
         clr = nRed;
         break;
      case 'N':   // Nitrogen
         clr = nBlue;
         break;
      case 'S':   // Sulfur
         clr = nOrange;
         break;
      default:    // Alles anders
         clr = nWhite;
         break;
      }

   return clr;
}

/*------------------------------------------------------------*/
/* ROUTINES USED IN MAKING THE DRAW LIST
   movedraw, 
   drawchain, 
   moveto

   The drawlist global is documented at the top of this
   module.
------------------------------------------------------------*/
void moveto(int whichAtom)
{
   drawlist[drawindx][1] = whichAtom;
   drawlist[drawindx++][0] = 0;  // 0 means "move to"
}

/*------------------------------------------------*/
void drawto(int whichAtom)
{
   drawlist[drawindx][1] = whichAtom;
   drawlist[drawindx++][0] = 1; // 1 means "draw to"
}


/*------------------------------------------------*/

typedef int maxBpaArray[MAXBPA];

// This routine worls hard to maximize chain lengths - and
// thereby cut down on the number of glBegin/glEnd pairs with
// line segment drawing.
void drawchain(int start, maxBpaArray* bondsFromAtom, int* numBondsFromAtom)
{
   int   j,
         pres = start,
         next,
         endofchain = 0,
         npoint = 0;

   do
      {
      j = 0;
      /* first try to find a bonded atom of the same color */
      while (j < MAXBPA && (attribs[bondsFromAtom[pres][j]].color != colour ||
                            bondsFromAtom[pres][j] == 0 || numBondsFromAtom[bondsFromAtom[pres][j]] % 2 == 1))
         {
         j++;
         }
      if (j == MAXBPA) // No luck - so look for even # left to draw
         {
         j = 0;
         /* now try to find one w/ even # of bonds left to draw */
         while (j < MAXBPA && (bondsFromAtom[pres][j] == 0 || numBondsFromAtom[bondsFromAtom[pres][j]] % 2 == 1))
            {
            j++;
            }
         if (j < MAXBPA)
            colour = attribs[bondsFromAtom[pres][j]].color;
         }
      if (j == MAXBPA)
         {
         j = 0;
         /* now take any */
         while (j < MAXBPA && bondsFromAtom[pres][j] == 0)
            {
            j++;
            }
         colour = attribs[bondsFromAtom[pres][j]].color;
         }
         
      if (j == MAXBPA)
         {
/*------------------------------------------------*/
         printf("fatal error: maximum bonds per atom exceeded. MAXBPA=%d in dispui.h\n",MAXBPA);
         exit(EXIT_FAILURE);
         }
      next = bondsFromAtom[pres][j];
      bondsFromAtom[pres][j] = 0;
      for (j = 0; j < MAXBPA; j++)
         {
         if (bondsFromAtom[next][j] == pres)
            bondsFromAtom[next][j] = 0;
         }
      numBondsFromAtom[pres]--;
      numBondsFromAtom[next]--;
      if (numBondsFromAtom[next] == 0)
         endofchain = 1; // End the while loop afte rthis one
      pres = next;
      drawto(pres);
      npoint++;
      
      // If we get a very long chain, let's start anew.
      if (npoint > 250)
         {
         moveto(pres);
         npoint = 0;
         }
      }
   while (!endofchain);
}

// Whenever the colors are changed, a new movedraw list must be computed.
// This allows for more efficient drawing, since sequences of same colored
// atoms which are bonded together can subsequently be rendered in a single
// glBegin/glEnd bracket.  This was probably especially important to efficiency
// in the original GL based version of MD Display - and it does not hurt now.
//
// Concept check: This routine does no rendering - it just prepares a new
// drawlist for when we next do a rendering.
void movedraw(void)
{
   register int i;
   int pres, btm, oddbtm, doneodd;
   int atm1, atm2;

   int *numBondsFromAtom = malloc_or_exit_if_error(MAXATM*sizeof(int),"movedraw - n");
   maxBpaArray *bondsFromAtom = malloc_or_exit_if_error(MAXATM*sizeof(maxBpaArray),"movedraw - n");

   /*-- INITIALIZE THE BOND ARRAYS --*/
   drawindx = 0;

   memset(numBondsFromAtom,0,MAXATM*sizeof(int));
   memset(bondsFromAtom,0,MAXATM*sizeof(maxBpaArray));

   /*
    * n[atm] contains number of bonds involving atm;
    * bnd[atm][*] contains the
    * atoms that atm is bonded to. */

   bondsToDraw = 0;
   // masterbonds is the global count o fbonds.
   /* only copy needed bonds */
   /* to working bond list */
   for (i = 0; i < masterbonds; i++)
      {
      atm1 = masterBondedAtomPairs[i][0];
      atm2 = masterBondedAtomPairs[i][1];
      if (((attribs[atm1].color != nBlack) && (attribs[atm2].color != nBlack))
          && (!delbondflag || bondtype[i] != delbond) || (catindx == (-1)))
         {
         bondedAtomPairsToDraw[bondsToDraw][0] = masterBondedAtomPairs[i][0];
         bondedAtomPairsToDraw[bondsToDraw][1] = masterBondedAtomPairs[i][1];
         bondsToDraw++;
         }
      }


   for (i = 0; i < bondsToDraw; i++)
      {
      atm1 = bondedAtomPairsToDraw[i][0];
      atm2 = bondedAtomPairsToDraw[i][1];

      bondsFromAtom[atm1][numBondsFromAtom[atm1]++] = atm2;
      bondsFromAtom[atm2][numBondsFromAtom[atm2]++] = atm1;
      }

   /* find atoms that aren't bonded to any others */
   // These are cations!
   // But, only do this if "cations" are passed in as -1 - Don't ask me why!
   // I guess we onl want to iterate over all this once.
   if (catindx == (-1)) // This is a first time flag.
      {
      catindx = 0;
      // We only do this the first time we run movedraw
      for (i = 1; i <= nat && catindx < MAXCAT; i++)
         if (!numBondsFromAtom[i])
            catlist[catindx++] = i;
      }

   btm = 1;
   oddbtm = 1;
   doneodd = 0;

/*-----------MAIN LOOP-----------------------*/
   do
      {
   /*-----find start of chain-----*/
      i = btm;
      // Skip over atoms with 0 bonds (i.e. cations)
      while (numBondsFromAtom[i] == 0 && i < nat)
         i++;

      // Let's start from here and try to build a chain
      btm = i;
      if (i < nat)
         {                         /* if there's anything left */
         if (!doneodd)
            {
            if (i > oddbtm)
               oddbtm = i;
            i = oddbtm;
            while ((attribs[i].color != colour || numBondsFromAtom[i] % 2 == 0) && i < nat)
               {
               i++;
               }

            // REVIEW  CWM changed from (if i= nat) to (if i== nat)
            if (i == nat)
               {
               i = oddbtm;
               while (numBondsFromAtom[i] % 2 == 0 && i < nat)
                  {
                  i++;
                  }
               }
            oddbtm = i;
            if (i < nat)
               {
               pres = i;
               colour = attribs[pres].color;
               moveto(pres);
               drawchain(pres, bondsFromAtom, numBondsFromAtom);
               }
            else
               {
               doneodd = 1;
               i = 0;
               }
            }
         else // (dooneodd == TRUE)
            {
            pres = i;
            while ((attribs[i].color != colour || numBondsFromAtom[i] == 0) && i < nat)
               i++;

            if (i < nat)
               {
               pres = i;
               }
            else
               {
               colour = attribs[pres].color;
               i = pres;
               }
            moveto(pres);
            drawchain(pres, bondsFromAtom, numBondsFromAtom);
            }
         }
      }
   while (i < nat);

   // printf("movedraw() complete for %d atoms... drawindx = %d\n",
      // nat,drawindx);

   free(bondsFromAtom);
   free(numBondsFromAtom);

}

// Forward declaration -
// By default, if user clicks we just want to get a label
void LabelOrUNLabelAtomPicked(void);
static int labelFlag; // Set to 1 if "labelling" 0 if "UNlabelling"

/*------------------------------------------------------------*/
/* COLORING & REMOVING */

static int newColorJustRequested;

// dispui.c calls here on completion of text input of
// a residue to color.
void NewColorInputComplete(int lastKey,const char* inputString)
{
   int res, low, high;
   resatm resatmarry[20];
   int lowmol, himol, moln;     /* MAS */
   if (lastKey == 13) // User pressed <ENTER> to complete input
      {
      int i = 0;
      // Parse user input into an array of atom ranges
      mas_parse(inputString, resatmarry);
      
      // whie user input remains...
      while (resatmarry[i].lores)
         {
         int j;

         // Get the range of molecules and residues
         // from this MAS molrange:resrage:atom segment
         low = resatmarry[i].lores;
         high = resatmarry[i].hires;
         lowmol = resatmarry[i].lowmol;
         himol = resatmarry[i].himol;  
         for (j = 1; j <= nat; j++)
            {
            res = attribs[j].res;
            moln = attribs[j].moln;

            // Is atom j in the user input?
            if (moln >= lowmol && moln <= himol
                && res >= low && res <= high
                && !stcmp(rlab[res], resatmarry[i].reslab)
                && (!stcmp(attribs[j].label, resatmarry[i].atmlab) ||
                    !resatmarry[i].atmlab[0]))
               {
               // Yes, we need to color atom j.  Recall 99 means "by atom"
               attribs[j].color =
                     (newColorJustRequested == 99 ? clrbyatm(attribs[j].label) : newColorJustRequested);

               // If we go to black be sure to turn off label!
               if (!newColorJustRequested)
                  attribs[j].showlabl = 0;
               }
            }
            i++;
         }

#ifdef DISPLAY_LISTS_IN_USE // Maybe we'll wire in later
      // We have to blow away the existing display lists, unfortunately.
      if (openGLDisplayListBaseIndex != 0)

         {
         glDeleteLists(openGLDisplayListBaseIndex,frames);
         openGLDisplayListBaseIndex = 0;
         }
#endif

      movedraw(); // We also must create a new "movedraw" list...
      checkhb();  // Also, make a list of h bonds in range

      // Rebuild the static molecule display list if needed
      staticNeedsRebuild = 1;
      // Update the movie
      MovieAdvanceTimerKeyAndPostRedisplay();
      }
   labelFlag = 1;
   AtomPickedFunc = LabelOrUNLabelAtomPicked;
}

// This function is called if the color input box is active and the user
// clicks on an atom
void NewColorAtomPicked(void)
{
   if (lastPickedAtom > 0 && lastPickedNeedsProcessing)
      {
      int newColorForThisAtom =
          (newColorJustRequested == 99 ? // Recall 99 means "color by atom"
               clrbyatm(attribs[lastPickedAtom].label) :
               newColorJustRequested);

      // If the color of the atom we selected is already "correct" then
      // color entire residue...
      // I.e. - clicking twice will always color the entire residue!
      if (attribs[lastPickedAtom].color == newColorForThisAtom)   /* COLOR ENTIRE RESIDUE */
         {
         int res = attribs[lastPickedAtom].res;
         int j;
         for (j = 1; j <= nat; j++)
            {
            if (attribs[j].res == res && attribs[j].color)
               {
               attribs[j].color = // Recall 99 means color by atom
                  (newColorJustRequested == 99 ? clrbyatm(attribs[j].label) : newColorJustRequested);
               }
            }
         }
      else
         {                      /* JUST COLOR THE ATOM */
         attribs[lastPickedAtom].color = newColorForThisAtom;

         // If we are going black, be sure to turn off label!
         if (!newColorJustRequested)
            attribs[lastPickedAtom].showlabl = 0;
         }


      lastPickedNeedsProcessing = 0;
      // We have to blow away the existing display lists, unfortunately.

#ifdef DISPLAY_LISTS_IN_USE // Maybe we'll wire in later
      if (openGLDisplayListBaseIndex != 0)
         {
         glDeleteLists(openGLDisplayListBaseIndex,frames);
         openGLDisplayListBaseIndex = 0;
         }
#endif

      movedraw(); // We also must create a new "movedraw" list...
      checkhb();
      staticNeedsRebuild = 1;
      MovieAdvanceTimerKeyAndPostRedisplay();
      }
}



// Called from dispui.c when color palette selection made or keystroke selected.
void displayNewColorRequested(int _newColor)
{
   static char prompt_string[100];

   newColorJustRequested = _newColor;
   sprintf(prompt_string,"Enter #mol:res@atoms to color %s:",
               (_newColor == 99 ? "by atom" : dispuiPaletteColorString(_newColor))
            );

   KeyboardInputStart(prompt_string,
                     NewColorInputComplete);

   AtomPickedFuncStart(NewColorAtomPicked);
}

// Remove is really just a matter of selecting "black" for residues
void displayRemoveRequested(void)
{
   newColorJustRequested = 0; // Now that will be black
   KeyboardInputStart("Enter #mol:res@atoms to remove:",
                     NewColorInputComplete);
}

// When the user has selected <Origin> from the menu and finishes
// the text box, this routine is called.
void CenterOfRotationComplete(int lastKey,const char* inputString)
{
   // static int active;
   int m,i,j,n,res,low,high;
   float com[3];
   resatm resatmarry[20];
   int lowmol,himol,moln;

  if (lastKey == 13) // Did user end input with a CR
      {
      movie.xTranslation = movie.yTranslation = 0;
      if (inputString[0] == 0) // Then no input entered...
         {
         // With no input, we remove any tranlation in the movie.
         movie.centerOfRotation[0] =
         movie.centerOfRotation[1] =
         movie.centerOfRotation[2] = 0.0;
         }
      else
         {
         const short* xx = FrameManagerGetFrame(primaryFM,0);

         // Use rmight have specified many molecules, residues,
         // or atoms.  We need to compute a center of maass, and go
         // with that!
         mas_parse(inputString,resatmarry);
         i=0;
         n=0;
         for (m=0;m<3;m++)
            com[m]=0.0;

         // Iterate over all user input fragments
         while (resatmarry[i].lores)
            {
            low=resatmarry[i].lores;
            high=resatmarry[i].hires;
            lowmol=resatmarry[i].lowmol;
            himol=resatmarry[i].himol;
            // Consider all atoms.  Are they in this input fragment?
            for (j=1;j<=nat;j++)
               {
               res=attribs[j].res;
               moln=attribs[j].moln;
               if (moln>=lowmol && moln<=himol
                   && res>=low && res<=high
                   && !stcmp(rlab[res],resatmarry[i].reslab)
                   && (!stcmp(attribs[j].label,resatmarry[i].atmlab) ||
                   !resatmarry[i].atmlab[0]))
                  {
                  // User, user wanted this atom to be included in
                  // com calculation
                   n+=1;
                   for (m=0;m<3;m++)
                     com[m] += (float)xx[3*j-3+m];
                  }
               }
            i++;
            }
            if (n>0)
               {
               for (m=0;m<3;m++)
                  movie.centerOfRotation[m] = com[m]/n;
               movie.xTranslation = movie.yTranslation = 0.0;
               }

// The old SGI code would move the mouse to the center of the screen
// I don't think a modern user interface would do this....  I don't
// even think GLUT can do this for us if we wanted to :(
#ifdef BUG_NEED_THIS_FEATURE_LATER_MAYBE
            cursor2origin(vp->fr,vp->bk);
#endif
         } // End else
      MovieAdvanceTimerKeyAndPostRedisplay();
      } // end if (lastKey == 13)
   labelFlag = 1;
   AtomPickedFunc = LabelOrUNLabelAtomPicked;
}

// If the mouse is clicked on an atom while the Center
// of rotation box is active, then make that atom the new
// center of rotation.
void CenterOfRotationAtomPicked(void)
{
   if (lastPickedAtom > 0 && lastPickedNeedsProcessing)
      {
      int m;
      const short* xx = FrameManagerGetFrame(primaryFM,0);
      for (m=0;m<3;m++)
         movie.centerOfRotation[m] = (float)xx[3*lastPickedAtom-3+m];

      movie.xTranslation = movie.yTranslation = 0.0;
      MovieAdvanceTimerKeyAndPostRedisplay();
      lastPickedNeedsProcessing = 0;
      }

}

// Begin keyboard input when the user clicks on the
// Origin menu box.
void displayCenterOfRotationRequested(void)
{
   KeyboardInputStart("Enter #mol:res@atom to make center of rotation:",
                     CenterOfRotationComplete);

   AtomPickedFuncStart(CenterOfRotationAtomPicked);
}

// We come here in the event the user has clicked on the
// label or unlabel box, and the user clicks on an atom.
void LabelOrUNLabelAtomPicked(void)
{
   if (lastPickedAtom > 0 && lastPickedNeedsProcessing)
      {
      int atom;

      if (labelFlag)
         attribs[lastPickedAtom].showlabl = (1 + attribs[lastPickedAtom].showlabl) & 3;
      else
         attribs[lastPickedAtom].showlabl = 0;
      activeLabelCount = 0;
      for (atom = 1; atom <= nat; atom++)
         {                      /* RE-DO SHOWLIST */
         if (attribs[atom].showlabl)
            {
            showlist[activeLabelCount][0] = atom;
            showlist[activeLabelCount++][1] = attribs[atom].showlabl;
            }
         }
      }
   MovieAdvanceTimerKeyAndPostRedisplay();
}

// When the user completes input during labeling or unlabeling,
// the thread of execution is tranferred here.
void LabelOrUNLabelComplete(int lastKey,const char* inputString)
{
   int i,j,low,high; // ,res,low,high;
   resatm resatmarry[20];
   int lowmol,himol;

   if (lastKey == 13)
      {
      // parse the user's text input into resatmarry segments
      mas_parse(inputString,resatmarry);
      i=0;
      // look for matching atoms in each fragment
      while (resatmarry[i].lores)
         {
         low=resatmarry[i].lores;
         high=resatmarry[i].hires;
         lowmol=resatmarry[i].lowmol;   /* MAS dolabel(), same as docolor() */
         himol=resatmarry[i].himol;     /* MAS */
         if (resatmarry[i].atmlab[0])
            {    /* if atoms specified, label them */
            for (j=1;j<=nat;j++)
               {
               if (attribs[j].res>=low && attribs[j].res<=high
                    && attribs[j].moln>=lowmol && attribs[j].moln<=himol /*MAS*/
                    && !stcmp(attribs[j].label, resatmarry[i].atmlab)
                    && !stcmp(rlab[attribs[j].res],resatmarry[i].reslab)
          && attribs[j].color)
                    attribs[j].showlabl=(attribs[j].showlabl&2)|labelFlag;
               }
            }
         else
            {                      /* no atoms specified, just label res */
            for (j=1;j<=nat;j++)
               {
               if ( attribs[j].moln>=lowmol && attribs[j].moln<=himol
                 && attribs[j].res>=low && attribs[j].res<=high
                  && attribs[j].color)
                  {
                  /* LABEL 2ND ATOM IN RES */
                  if (attribs[j+1].res==attribs[j].res && labelFlag
                        && attribs[j+1].color) j++;
                  if ( !stcmp(rlab[attribs[j].res],resatmarry[i].reslab) )
                       attribs[j].showlabl=(attribs[j].showlabl&1)|(2*labelFlag);
             /* SKIP OVER REST OF ATOMS IN RESIDUE IF LABELING */
                 if (labelFlag)
                    while (attribs[j].res==attribs[j+1].res) j++;
                 }
               }
            }
         i++;
         }
      activeLabelCount=0;
      for (i=1;i<=nat;i++)
         {
         if (attribs[i].showlabl)
            {
            showlist[activeLabelCount][0]=i;
            showlist[activeLabelCount++][1]=attribs[i].showlabl;
            }
         }

      MovieAdvanceTimerKeyAndPostRedisplay();
      } // End if lastKey == 13

   labelFlag = 1; // Return to labelling if clicks made
}

// Code called by dispui.c when Label option taken by user
void displayLabelRequested(void)
{
   labelFlag = 1;
   KeyboardInputStart("Enter #mol:res@atoms to label:",
      LabelOrUNLabelComplete);

   AtomPickedFuncStart(LabelOrUNLabelAtomPicked);
}

// Code called by dispui.c when UnLabel option taken by user
void displayUNLabelRequested(void)
{
   labelFlag = 0;
   KeyboardInputStart("Enter #mol:res@atoms to UN-label:",
      LabelOrUNLabelComplete);
   AtomPickedFuncStart(LabelOrUNLabelAtomPicked);
}



/*-----------------------------------------------------------*/
/* ROUTINES FOR COLLECTING/Displayinh REQUESTED MEASUREMENTS */
/*-----------------------------------------------------------*/

// This is just an error message in case we maxed out on measurements.
int displayMaxMeasurementsReached(void)
{  
   int retval = FALSE;
   if (view.nmeas >= nMaxMeasurements)
      {
      retval = TRUE;
      // this message will appear on top of screen for 5 secons - if movie
      // is running.
      sprintf(dispuiStatusMessage,"Maximum Measurements Reached.  Remove by clicking at right.");
      dispuiStatusMessageExpiration = time(NULL)+5;
      MovieAdvanceTimerKeyAndPostRedisplay();
      }
   return retval;

}

// Here are the three types of measurements MD Display can
// make
enum {nDistance,nAngle,nDihedral} measurementOfInterest;

// Number of atom inputs wanted
static int nwanted;

// Number of atoms (0..nwanted-1) picked so far.
static int npicked;

// When an atom is clicked on after user selection of -, > or Z
// we come here.
void MeasurementAtomPicked(void)
{
   if (  lastPickedAtom > 0 &&
         lastPickedNeedsProcessing &&
         view.nmeas < nMaxMeasurements)
      {
      int i;
      int duplicate = FALSE;
      // In case user reclicks on same atom, ignore that.
      for (i=0;i<npicked;i++)
         if (lastPickedAtom == view.measarry[view.nmeas][i])
            duplicate = TRUE;
      if (! duplicate)
         {
         // Add this atom to the list.
         view.measarry[view.nmeas][npicked++]=lastPickedAtom;
         // If we have picked all the atoms that we need for this
         // type of measurement, then update the display.
         if (npicked==nwanted)
            {
            for ( ; npicked<4; npicked++)
               view.measarry[view.nmeas][npicked] = 0;
            view.nmeas++;
            if (displayMaxMeasurementsReached())
               KeyboardStuff("\033"); // This should close keyboard input!
            MovieAdvanceTimerKeyAndPostRedisplay();
            npicked = 0;
            }
         }
      lastPickedNeedsProcessing = 0;
      }
}

// The user can specify measreuments between atoms using text input
// We come here when the user types <ENTER> or <ESCAPE>
void MeasurementComplete(int lastKey,const char* inputString)
{
   int gres, atm = 0, i;
   char lbres[5];
   resatm raarry[10];

   if ((lastKey == 13) && strlen(inputString) && view.nmeas < nMaxMeasurements)
      { // Keyboard Input complete
      mas_parse(inputString, raarry);

      do
         {
         gres = raarry[atm].lores;
         for (i = 0; i < 4; i++)
            {
            lbres[i] = raarry[atm].atmlab[i];
            }
         lbres[4] = 0;
         i = 0;
         while (i <= nat &&
                (attribs[i].res != gres || stcmp(attribs[i].label, lbres)))
            {
            i++;
            }
         view.measarry[view.nmeas][atm++] = i;
         }
      while (gres && i <= nat && atm < nwanted);

      if (atm == nwanted && gres && i <= nat)
         {
         for (; atm < 4; atm++)
            view.measarry[view.nmeas][atm] = 0;
         view.nmeas++;
         }
      nwanted = 0;
      } // End if lastKey == 13

   labelFlag = 1;
   AtomPickedFunc = LabelOrUNLabelAtomPicked;
}

// Start text input and mouse selection when
// user clicks on '-'
void displayDistanceMeasurementRequested(void)
{
   if (! displayMaxMeasurementsReached())
      {
      measurementOfInterest = nDistance;
      nwanted = 2;
      npicked = 0;
      KeyboardInputStart("Enter #mol:res@atoms (or pick them) for new DISTANCE:",
         MeasurementComplete);
      AtomPickedFuncStart(MeasurementAtomPicked);
      }
}

// Start text input and mouse selection when
// user clicks on '>'
void displayAngleMeasurementRequested(void)
{
   if (! displayMaxMeasurementsReached())
      {
      measurementOfInterest = nAngle;
      nwanted = 3;
      npicked = 0;
      KeyboardInputStart("Enter #mol:res@atoms (or pick them) for new ANGLE:",
         MeasurementComplete);
      AtomPickedFuncStart(MeasurementAtomPicked);
      }
}

// Start text input and mouse selection when
// user clicks on 'z'
void displayDihedralMeasurementRequested(void)
{
   if (! displayMaxMeasurementsReached())
      {
      nwanted = 4;
      npicked = 0;
      measurementOfInterest = nDihedral;
      KeyboardInputStart("Enter #mol:res@atoms (or pick them) for new DIHEDRAL:",
         MeasurementComplete);
      AtomPickedFuncStart(MeasurementAtomPicked);
      }
}

// When user finished entry of Hbond cutoff, we come here
// to parse that value, and redo the hbond list.
void HBondComplete(int lastKey,const char* inputString)
{
      if (lastKey == 13)
         {
         sscanf(inputString, "%f", &view.hbcut);

         checkhb();
         }
}

// Input an new hbond cutoff when user clicks on H-bond 
void displayHBondRequested(void)
{
   KeyboardInputStart("Enter cutoff length for H-bonds (A):",
      HBondComplete);
}

// As the user whether or not they really want to quit.
// If YES, then exit... which will close down GLUT and call
// functions registered for atexit() processing.
void displayQuitSaveCurrentSessionComplete(int lastKey,const char* inputString)
{
   if (lastKey == 13)
      {
      if (inputString[0] == 'y' || inputString[0] == 'Y')
         {
         displaySaveCurrentSession = 1;
         }
      displayGracefulExit = 1;
      exit(EXIT_SUCCESS);     // Good bye - see atexit() references
      }

}

// If the user really did mean to quite, then prompt them about whether
// they want to save the session state.
void displayQuitAreYouSureComplete(int lastKey,const char* inputString)
{
   if (lastKey == 13)
      {
      if (inputString[0] == 'y' || inputString[0] == 'Y')
         KeyboardInputStart("Save Current Session?",
            displayQuitSaveCurrentSessionComplete);
      }

}

// If the user accidentlly clicks on Quit or types Q, then
// let them say "no - I did not mean it"
void displayQuitRequested(void)
{
   KeyboardInputStart("\001Quit Requested: Are you sure?",
      displayQuitAreYouSureComplete);
}

/*-----------------------------------------------------------------*/
/* CREATE A PDB FILE CONTAINING AN IMAGE OF THE CURRENT FRAME */
void displayPDBComplete(int lastKey,const char* inputString);
static int PDBframe;

void dopdb(const char* PDBfilename,int frame)
{
   FILE *pdbFILE;
   int res, j;


   const short* xx = FrameManagerGetFrame(primaryFM,frame); // ,xx + i * 3 * (nat + oneIfBox);

#ifdef BAD_IDEA_TESTING_PS
   Matrix dummy1,dummy2;
   display_makeps(&view,xx,dummy1,dummy2);
   return;
#endif

   if (PDBfilename[0] != ' ')
   {
   pdbFILE = fopen(PDBfilename, "wb");
   if (pdbFILE)
      {
      fprintf(pdbFILE, "REMARK Created by MD Display\n");
      fprintf(pdbFILE, "REMARK Frame %d of data set %s\n", frame, setFileName);
      for (j = 1; j <= nat; j++)
         {
         // Note the pdb conversion code below has been adapted from
         // the ambpdb utility supplied with Amber 7.0
         char lbres[4];
         char atnam[5];

         res = attribs[j].res;

         memcpy(lbres,rlab[res],3);
         lbres[3] = 0;
#ifdef RESIDUE_TRANSLATION_DESIRED
// Note from Chris on version 3.0...  We had a little bit of debate
// inside the Lybrand lab on this - and decided we did not want this
// translation scheme after all for generated .pdbs
         {
         int i;
         static char* xlate_list[][2] =
           {{"HID","HIS"},
            {"HIE","HIS"},
            {"HIP","HIS"},
            {"HIC","HIS"},
            {"CYX","CYS"},
            {"CYM","CYS"},
            {"MEM","MET"},
            {"ASH","ASP"},
            {"GLH","GLU"},
            {0,0}
            };

         for (i=0;xlate_list[i][0];i++)
            if (! memcmp(lbres,xlate_list[i][0],3))
               {
               memcpy(lbres,xlate_list[i][1],3);
               break;
               }
         }

         { // Change to G, A, C, U, T if RG DG RA DA etc..
         int i;
         static char* nucleic_list[] =
           {"RG","DG","RA","DA","RC","DC","RU","DT",0};

         for (i=0;nucleic_list[i];i++)
            if (! memcmp(lbres,nucleic_list[i],3))
               {
               sprintf(lbres,"  %c",nucleic_list[i][1]);
               break;
               }
         }
#endif


// This code converts " XXX" atom names to "XXX "
         sprintf(atnam+1,"%-3.3s",attribs[j].label);
         atnam[0] = ' ';
         atnam[4] = 0;

         if (attribs[j].label[3] && (attribs[j].label[3] != ' '))
            {
            atnam[0] = attribs[j].label[3];
            }

#ifdef ATOM_NAME_TRANSLATION_DESIRED
//       --- here are some more Brookhaven wraparounds:
//           This gives files that look very much like Brookhaven, EXCEPT
//           that Brookhaven uses "1" and "2" for beta-protons (for example)
//           whereas the standard Amber database (along with many in
//           the NMR field) use "2" and "3", i.e. we have 2HB and 3HB,
//           whereas Brookhaven files use 1HB and 2HB.
//
// Note from Chris on version 3.0...  We had a little bit of debate
// inside the Lybrand lab on this - and decided we did not want this
// translation scheme after all for generated .pdbs
//
            if (atnam[0] == ' ' && atnam[1] == 'H' &&
               (atnam[3] == '1' || atnam[3] == '2' || 
                atnam[3] == '3') &&
               ((atnam[1] == 'H' && atnam[2] == 'B') || 
                     (memcmp(lbres,"PHE",3) &&
                      memcmp(lbres,"TYR",3) &&
                      memcmp(lbres,"TRP",3) &&
                      memcmp(lbres,"HI",2) )))
               {
               atnam[0] = atnam[3];
               atnam[3] = ' ';
               }
//
//       --- Convert nucleic acid primed names into asterisk: these
//           are always in the fourth column of the atom name:
//
  
            if (atnam[3] == '\'')
               atnam[3] = '*';    
//
//       --- Now, special case out the two-character element names:
//

         { // Change " CL" to "Cl " etc
         int i;
         static char* ion_list[] = 
           {"NA+","FE","CL","ZN","LI+","CA+","MG+","BR-",0};
            
         if (atnam[0] == ' ')
            for (i=0;ion_list[i];i++)
               if (atnam[1] == ion_list[i][1] && 
                   toupper(atnam[2]) == ion_list[i][2] &&
                  (ion_list[i][3] == 0 || atnam[3] == ion_list[i][3]))
            {
             atnam[0] = atnam[1];
             atnam[1] = atnam[2];
             atnam[2] = atnam[3];
             atnam[3] = ' ';
            }
         }
#endif

         fprintf(pdbFILE, "ATOM  %5d %-4.4s %-3.3s %5d    %8.3f%8.3f%8.3f\n",
                 j, atnam, lbres, res,
                 (float) xx[j * 3 - 3] / KAPPA, (float) xx[j * 3 - 2] / KAPPA,
                 (float) xx[j * 3 - 1] / KAPPA);
         }
      sprintf(dispuiStatusMessage,"PDB file %s of frame %d successfully created.",PDBfilename,frame);
      fclose(pdbFILE);
      }
   else
      sprintf(dispuiStatusMessage,"Unable to open %s: PDB not created!",PDBfilename);
   }
   else
      sprintf(dispuiStatusMessage,"Invalid filename entry: [%s] PDB not created!",PDBfilename);

   // Show the dispui Status Message for 4 seconds.
   dispuiStatusMessageExpiration = time(NULL) + 4; // Display status for 3 seconds
   MovieAdvanceTimerKeyAndPostRedisplay();
}

// This function is called by dispui.c when the user finishes
// the input text for the name of the pdb file.
void displayPDBComplete(int lastKey,const char* inputString)
{
   if ((lastKey == 13) && strlen(inputString))
      dopdb(inputString,PDBframe);
}

// Called by dispui.c when the user clicks on PDB
// Called frmo below when the user presses 'P'
void displayPDBRequested(void)
{
   char outname[200];
   char prompt[50];

   PDBframe = dispuiMovieFrame();

   sprintf(prompt,"Enter PDB Filename for dump of frame %d ",PDBframe);
   KeyboardInputStart(prompt,displayPDBComplete);
   sprintf(outname, "%s.%d", setFileName, PDBframe);
   KeyboardInputSetDefault(outname);
}

// If dispui.c does not act on a pressed key, it calls
// displayKeyboardFunc and we have a chance here to act on
// the key at our level.
void displayKeyboardFunc(int key)
{
   switch (key)
      {
      case 'b':   // Keystroke for "by atom"
      case 'B':   // Keystroke for "by atom"
         displayNewColorRequested(99); // Call display to do the real work...
         break;

      case 'F':   // Prompt for filtering parameters and build filtered coords
         displayInputNewFilter();
         break;

      case 'f':   // Toggle filtered display if filtering performed before
         if (filterSize)
            {
            filterFlag = (! filterFlag);
            MovieAdvanceTimerKeyAndPostRedisplay();
            }
         break;

      case 'l':   // Label atoms and/or residues
      case 'L':
         displayLabelRequested();
         break;

      case 'u':   // Remove labels - "UnLabel"
      case 'U':
         displayUNLabelRequested();
         break;

      case '-':
         displayDistanceMeasurementRequested();
         break;

      case '>':
         displayAngleMeasurementRequested();
         break;

      case 'Z':
      case 'z':
         displayDihedralMeasurementRequested();
         break;

      case 'H':
         displayHBondRequested();
         break;

      case 'P':
         displayPDBRequested();
         break;

      case 'Q':
      case 'q':
         displayQuitRequested();
         break;

      default:
         break;
      }
}


//------------------------------------------------------------
// ROUTINES USED IN FILTERING THE ATOMIC TRAJECTORIES
// to attenuate high frequency motions from the trajectory
// dofilt, cpcoords, putinbox, sa, safilt, simpfilt
//
// Be sure to see the 1996 paper in J Molecular Graphics
// which has several paragraphs explaining the algorithm.
//
// The 1996 paper gives this reference for more information on
// FIR filtering:
// 8 0ppenheim, A.V., and Schafer, R.W. Discrete-Time
//   Signal Processing. Prentice Hall, Englewood Cliffs,
//   New Jersey, 1989, Ch, 7, pp. 444-452
//------------------------------------------------------------

static float filterSamplingWaveNumber;

float sa(float x)
{
   float ans;

   if (x == 0.0)
      ans = 1.0;
   else
      ans = sin(x) / x;

   return ans;
}

/*------------------------------------------------------------*/

int safilt(short* xxs, float cof)
{
   int ord = 4, ncord, co, fr, i;
   int imp[MAXFILT*2];
   register float tmp;
   float sum;

   const short* xx = FrameManagerGetFrame(primaryFM,0);

   ncord = 3 * (nat + oneIfBox);
   for (i = 0; i <= ord; i++)
   {                            /* windowed Sa function */
      // REVIEW - This seems very prone to overflow to CWM.. Especially look at SA
      imp[ord - i] =
         (short)((ord - i + 1) * (short) (1000.0 * sa(2.0 * M_PI * cof * (float) i)));
      imp[ord + i] = imp[ord - i];
   }

   sum = 0.0;
   for (i = 0; i <= 2 * ord; i++)
      sum += (float) imp[i];

   /* FILTER COORDS & BOX DIMENSIONS */
   for (co = 0; co < ncord; co++)
   {
      for (fr = 0; fr < frames - 2 * ord; fr++)
      {
         tmp = 0.0;
         /* CONVOLUTION */
         for (i = 0; i < 2 * ord + 1; i++)
            tmp += (float) (imp[i]) * (float) (xx[(fr + i) * ncord + co]);
         xxs[fr * ncord + co] = (short) (tmp / sum);
      }
   }

   return (9);
}

/*------------------------------------------------------------*/
// Returns filtsize
int simpfilt(short* xxs, int ord)
{
   int filtsize, n, ncord, co, fr, i;
   int imp[MAXFILT*2];
   register float tmp;
   float sum;

   const short* xx = FrameManagerGetFrame(primaryFM,0);

   ncord = 3 * (nat + oneIfBox);

   n = (ord > (MAXFILT / 2) ? (MAXFILT / 2) : ord);

   filtsize = n * 2;

   for (i = 0; i <= n; i++)
   {                            /* /\ window function */
      imp[n - i] = (n - i + 1);
      imp[n + i] = imp[n - i];
   }

   sum = 0.0;
   for (i = 0; i <= 2 * n; i++)
      sum += imp[i];

   /* FILTER COORDS , BOX DIMENSIONS */
   for (co = 0; co < ncord; co++)
   {
      for (fr = 0; fr < frames - filtsize; fr++)
      {
         tmp = 0.0;
         for (i = 0; i <= filtsize; i++)
            tmp += (float) (imp[i]) * (float) (xx[(fr + i) * ncord + co]);
         xxs[fr * ncord + co] = (short) (tmp / sum);
      }
   }

   return filtsize;

}

// After the user completes input (CR or ESC) of
// a cutoff freq or Number of frames to average over...
// we need to parse that out and add filtering
void FilterComplete(int lastKey,const char* inputString)
{
   int n;

  if (lastKey == 13) // Did user end input with a CR
      {
      sscanf(inputString, "%d", &n);
      if (n > 1)
         {
         // The computations can take time - so show
         // an hour glass in all subwindows
         dispuiSetWaitCursors();

         filterFlag = 1;
         if (picoSecondsPerFrame > 0.0)
            {
            /* cutoff freq expressed as fraction of sampling freq */
            float df = (float) n / filterSamplingWaveNumber;
            if (df > 0.1 && df < 0.5)
               filterSize = safilt(xxs, df);
            if (df <= 0.1)
               filterSize = simpfilt(xxs, (int) (0.33 / df + 1.0));
            }
         else
            if (n >= frames)
               n = frames - 1;

         filterSize = simpfilt(xxs, n / 2);
         if (oneIfBox)
            {
            if (isTruncatedOctahedron)
               wrap_to(xxs,filterSize);
            else
               putinbox(xxs, filterSize);
            }

         // Return to the cursor styles we had before in each subwindow
         dispuiPopCursors();
         }
         else
         {
            filterSize = 0;
            if (oneIfBox)
            {
               cpcoords(xxs,primaryFM);
               if (isTruncatedOctahedron)
                  wrap_to(xxs,filterSize);
               else
                  putinbox(xxs, filterSize);
               
            }
         }


#ifdef DISPLAY_LISTS_IN_USE
      // We have to blow away the existing display lists, unfortunately.
      if (openGLDisplayListBaseIndex != 0)

         {
         glDeleteLists(openGLDisplayListBaseIndex,frames);
         openGLDisplayListBaseIndex = 0;
         }
#endif

      MovieAdvanceTimerKeyAndPostRedisplay();
      }
}

// When the user presses 'F' or clicks on the filter menu box, we come here
// and ask a cutoff freq question, the flavor of which depends on the time
// interval (or lack thereof) between frames.
void displayInputNewFilter(void)
{
   const char* promptString;

   if (picoSecondsPerFrame > 0.0)
      {
      filterSamplingWaveNumber =
         (1.0e+12) / (picoSecondsPerFrame * (3.0e+10));   /* sampling wave number */

      promptString = "Enter desired cutoff frequency (1/cm) for low-pass filter:";
     }
   else
      {
      promptString = "Enter number of frames to average over:";
      }

   KeyboardInputStart(promptString,
                     FilterComplete);
}

// When we first crank up the program at the OS prompt, there are
// a number of questions which are asked in text mode.
// Here we give the user the option of setting some colors before
// the movie is displayed.  It seems very clunky do be prompting
// for things in this way in 2003.  But, I'll leave it in case someone
// out there is running the program with scripts.

void dorescolor(void)
{
   int high, low, i, j, res, clr = (-1);
   char instr[60];
   resatm resatmarry[20];
   int lowmol, himol, moln;     /* MAS */

   printf("\n________ Initial coloring of residues & atoms ________\n\n");
   do
   {
      int palette_color;

      for (palette_color = 0;palette_color < 10;palette_color+=2 )
         printf("         %2d) %-20.20s    %2d) %-25.25s\n",
                 palette_color,dispuiPaletteColorString(palette_color),
                 palette_color+1,dispuiPaletteColorString(palette_color+1)
                 );

      printf("         98) each molecule a different color\n");
      printf("         99) or <b>: automatic color by atom type\n\n");
      SafeInputString("Enter optional initial coloring (0-99 or <b>) or press <ENTER> :",
         instr,sizeof(instr)  );
      clr = -1;
      if (instr[0] == 'b')
         clr = 99;
      else
         sscanf(instr, "%d", &clr);
      /* MAS 3-22-90 dorescolor() MAS color each molecule a different color */
      if (clr == 98)
      {
         for (j = 1; j <= nat; j++)
         {
            /* colors 1 through 7 for each molecule */
            attribs[j].color = (short)((attribs[j].moln % 6) + 1);
         }
      }

      if ((clr >= 0 && clr <= 9 || clr == 99) && instr[0])
      {
         SafeInputString("Enter #mol:res@atoms --> ",
            instr,sizeof(instr)  );

         mas_parse(instr, resatmarry);
         i = 0;
         while (resatmarry[i].lores)
         {
            low = resatmarry[i].lores;
            high = resatmarry[i].hires;
            lowmol = resatmarry[i].lowmol;   /* MAS */
            himol = resatmarry[i].himol;  /* MAS */
            for (j = 1; j <= nat; j++)
            {
               res = attribs[j].res;
               moln = attribs[j].moln; /* MAS */
               if (moln >= lowmol && moln <= himol /* MAS */
                   && res >= low && res <= high
                   && !stcmp(rlab[res], resatmarry[i].reslab)
                   && (!stcmp(attribs[j].label, resatmarry[i].atmlab) ||
                       !resatmarry[i].atmlab[0]))
               {
                  attribs[j].color =
                     (short)(clr == 99 ? clrbyatm(attribs[j].label) : clr);
                  if (!clr)
                     attribs[j].showlabl = 0;
               }
            }
            i++;
         }
      }
   } while (clr >= 0 && clr != 98 && instr[0]);
}

/*----------------------------------------------------------*/

#ifdef CUbIC_3D_CURSOR_OPTION
// The original MDDisplay had a 3D cubic cursor option.
// This may be something we want to bring back later - so I've
// left here in the original source code for a reminder.
void draw3Dcursor(int i, float size, float front, float back)
{
   float x, y, fw, z;

   if (curs3D && !getbutton(MOUSE1))
   {
      pushmatrix();
      x = (curs3Dx / XMAXSCREEN - 0.5) * 2 * size;
      y = (curs3Dy / YMAXSCREEN - 0.5) * size;
      fw = (100. + curs3Dz) / 200.;
      z = size * (fw * (front + back) - back);
      translate(x, y, z);
      if (i)
         getmatrix(cursormat);
      color(RED);
      recti(-0.01 * size, -0.01 * size, 0.01 * size, 0.01 * size);
      popmatrix();
      color(YELLOW);
   }

}

void dozcurs(const dispinfo* vp)

{
   if (keypress == 'z' || butt && my > 330 && my < 370 && mx > 0 && mx < 90)
   {
      curs3D = !curs3D;
      if (curs3D)
         cursor2origin(vp->fr, vp->bk);
   }
}
#endif






// On program startup, the first thing we need to do is figre out
// what atoms and residues we have, what colors they are, and what
// is bonded to what.  We also need to allocate some globals.
// We read preproc's output file to accomplish all this.

void readparms_and_allocate_globals(int recoverFlag)
{
   int* inpt;
   FILE *f;
   int bytesread;
   int i, totbonds, bndtype;

   // REVIEW - I think we mean long here!
   inpt = (int*) malloc_or_exit_if_error(3*MAXDRAW * sizeof(int),"readparms - inpt");

   // 1981 style platform independent clear screen function :)
   printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
   DisplayStdoutCredits();

   printf("Loading Connectivity Information from %s:  Standby..... \n",bondFileName);

   f = fopen_or_exit_if_error(bondFileName, "rb");

   bytesread = fread(inpt, 1, 3 * MAXDRAW * sizeof(*inpt), f);

   printf("%d bytes read: ", bytesread);

   if (bytesread < (3 * sizeof(*inpt)))
      {
      printf("Insufficient Input from %s\n",bondFileName);
      exit(EXIT_FAILURE);
      }

   // Important globals are in the header of the binary file we just read in
   nat = inpt[0];
   nres = inpt[1];
   totbonds = inpt[2];

#ifdef REVIEW_CWM_ELIMINATED_THIS_ANCIENT_BACKWARDS_COMPATABILITY
   if (nat > 1000000 || nres > 100000)
   {                            /* probably old (short) bond data */
      fclose_or_exit_if_error(f,bondFileName);
      printf("Assuming old (short) bond data based on header\n");

      f = fopen_or_exit_if_error(bondFileName, "rb");

      bytes = sizeof (inpt_s);
      bytesread = read(fdes, (char *) inpt_s, (unsigned int) bytes);
      printf("%d bytes read.\n", bytesread);
      nat = inpt_s[0];
      nres = inpt_s[1];
      totbonds = inpt_s[2];
      for (i = 3; i < bytesread / 2; i++)
         inpt[i] = inpt_s[i];
   }
#endif

   printf("  %d atoms, %d residues, %d bonds.\n", nat, nres, totbonds);

   if (nat > MAXATM)
      {
      printf("Exceeded limit of #define MAXATM=%d atoms in displib.h.\n\n(Apologies: I need to fix obsolete hard dimensions.)\n\n", MAXATM);
      exit(EXIT_FAILURE);
      }

   fclose_or_exit_if_error(bondFileName,f); // New CWM code

   // If -r keyed on command line, then we are resuming a crashed movie
   if (recoverFlag) // apend .xxx temporarily to this file name
      {
      i = strlen(attribFileName);
      strcpy(attribFileName+i,".xxx");
      printf("\nRecovering atom attributes from ungraceful shutdown: %s\n",attribFileName);
      }

   printf("Loading Connectivity Information from %s:  Standby..... \n",attribFileName);
   attribs = malloc_or_exit_if_error(sizeof(atomattrib)*(nat+nres+2),"readparms");
   f = fopen_or_exit_if_error(attribFileName, "rb");

   fread_or_exit_if_error(attribs,sizeof(atomattrib),nat + nres + 1,f);
   fclose_or_exit_if_error(attribFileName,f); // New CWM code

   if (recoverFlag)
      attribFileName[i] = 0;

   attribytesread = sizeof(atomattrib) * (nat + nres + 1);

   // Some attribute fails have an addition reacord if "oneIfBox" was
   // turned on during preproc.

   attribytesread +=
      sizeof(atomattrib) *
            fread(&attribs[nat+nres+1],sizeof(atomattrib),1,f);

   printf("%d attribute bytes read (%d attribute records).\n",
      attribytesread,attribytesread/sizeof(atomattrib));

   /* CHECK FOR BOX COORDS */
   if (attribs[nat].moln == 999)
   {
      boxatm = nat;
      nat--;
      oneIfBox = 1;
      // if (getenv("KRISTINA"))
         // oneIfBox = 1;
   }
   
   /* CHECK FOR BOX COORDS */
   if (attribs[nat].moln == 1000)
   {
      boxatm = nat;
      nat--;
      oneIfBox = 1;
      isTruncatedOctahedron = 1; // Bounding box is Truncated octahedron
      // if (getenv("KRISTINA"))
         // oneIfBox = 1;
   }

   /* COMPILE label SHOWLIST FROM the ATTRIBUTES JUST READ IN */
   activeLabelCount = 0;
   for (i = 1; i <= nat; i++)
   {
      if (attribs[i].showlabl)
      {
         showlist[activeLabelCount][0] = i;
         showlist[activeLabelCount++][1] = attribs[i].showlabl;
      }
   }

   /* ---- COPY RESIDUE LABELS TO THEIR OWN ARRAY ---- */
   rlab = malloc_or_exit_if_error(sizeof(rlab_type) * (nres+1),"readparms: rlab");
   for (i = 1; i <= nres; i++)
   {
      if (attribs[i + nat + oneIfBox].res != i ||
          attribs[i + nat + oneIfBox].color != (-1))
         rlab[i][0] = 0;
      else
      {
         int j;
         for (j = 0; j < 3; j++)
            rlab[i][j] = attribs[i + nat + oneIfBox].label[j];
         rlab[i][3] = 0;
      }
   }

   if (attribs[0].label[0] == (char)0xff) // CWM notes - was 255
   {
      delbondflag = 1;
      delbond = attribs[0].res;
   }

   // Get list of what is bonded to what
   // And what might be a hydrogen bond
   masterbonds = 0;
   masterhbonds = 0;
   for (i = 0; i < totbonds; i++)
   {
      bndtype = inpt[3 * i + 5];
      if (bndtype < 0)
      {                         /* HYDROGEN BOND */
         masterhblist[masterhbonds][0] = inpt[3 * i + 3] / 3 + 1;
         masterhblist[masterhbonds][1] = inpt[3 * i + 4] / 3 + 1;
         masterhbonds++;
      }
      else
      {                         /* NORMAL BOND */
         masterBondedAtomPairs[masterbonds][0] = inpt[3 * i + 3] / 3 + 1;
         masterBondedAtomPairs[masterbonds][1] = inpt[3 * i + 4] / 3 + 1;
         bondtype[masterbonds] = bndtype;
         masterbonds++;
      }
#ifdef LATER_REVIEW
NOTE NEED FOR bondtype overwrite protection somewhere
#endif
   }
   printf("%d bonds, %d hydrogen bonds.\n", masterbonds, masterhbonds);

   // prompt user for initial coloring
   dorescolor();

   /* MAKE THE MOVE/DRAW LIST -- FIRST FOR CATIONS */
   catindx = (-1);
   movedraw();
   // set the global startNeedsRebuild flag so we know to
   // make a new display list when the user clicks on Static
   staticNeedsRebuild = 1;

   /* Now call movedraw "for real" */
   // Actually, I have no idea why we do this twice :(
   movedraw();

   free(inpt);
}

// Set the view.statictype based on user's submenu
// selection in the Static Windows/menu box
void staticMenuCallback(int whichOption)
{
   if (view.statictype != whichOption) // nStaticFirst, nStaticLast, nStaticNone
      {
      view.statictype = whichOption;
      staticNeedsRebuild = 1;

      MovieAdvanceTimerKeyAndPostRedisplay();
      }
}


// dispui.c needs to know how the maximum frame # it might
// have to show.  We index ours from 0 - so here is the code:
int displayMovieMaxFrame(void)
{
   if (filterSize && filterFlag)
      return frames-filterSize-1;
   else
      return frames-1;
}

// dispui.c needs to know the minimum frame number when it is time to
// roll past start or end of movie.
int displayMovieMinFrame(void)
{
   return startFrame;
}

// If users command line makes no sense, or user wants help, here
// it is:
void mddisplay_usage_message(void)
{
   printf("\nUsage: mddisplay -n dataset [-f #frames]\n\n");
   printf("Additional mddisplay options and information:\n");
   printf("--------------------------------------------\n");
   printf("         -r resume from interrupted session\n");
   printf("         -q use quad buffered stereo\n");
   printf("         -f can specify a range Ex: -f 50,100\n");
   printf("\n\n");
   DisplayStdoutCredits();
}

int main(int argc, char** argv)
{
   char ch;
   int i, c;
   extern char *optarg;
   int recoverFlag = 0; // set to 1 if user requests -r in command line
   short quad_stereo_flag = 0; // True if user starts with -q option
   static char windowTitle[500];
   sprintf(windowTitle,"MDDisplay");

   /* read in command line arguments */

   while ((c = getopt(argc, argv, "qa:b:c:f:h:n:r")) != -1)
      switch (c)
      {
         // In the old days, users might take trouble to specify
         // names for attributes, coordinates, and bond files separately
         case 'a':
            sprintf(attribFileName, "%s", optarg);
            break;
         case 'c':
            sprintf(coordFileName, "%s", optarg);
            break;
         case 'b':
            sprintf(bondFileName, "%s", optarg);
            break;

         // Without this flag, F4 gives you legacy stereo on SGI
         // and nothing on non-SGI platforms.
         case 'q':
            quad_stereo_flag = 1;
            break;

         // You can optionally specify a range of frames at startup
         // You don't have to load in all the binary data created by preproc!
         case 'f':
            i = sscanf(optarg, "%d%c%d", &startFrame, &ch, &endFrame);
            if (i < 3)
            {
               frames = startFrame;
               startFrame = 0;
               endFrame = frames - 1;
            }
            else
               frames = endFrame - startFrame + 1;
               
            if (frames <= 0)
               {
               printf("\n  -f:  Invalid frame counts.\n\n");
               exit(2);
               }
            break;
            /*
             * case 'h': sscanf(optarg,"%f",&hbcut); break; */

         // -n moveFile is usually what we use to start MD Display 3.0.
         // I add .bnd, .cor, .attr etc to make all the component file names
         case 'n':
            {
            char cwd_str[300];
            char pathsep = '/';
#ifdef _MSC_VER
            _getcwd(cwd_str,sizeof(cwd_str)-1);
#else
            getcwd(cwd_str,sizeof(cwd_str)-1);
#endif
            // printf("%s\n", optarg);
            i = strlen(optarg);
            if (i > 5 && strncmp(&optarg[i - 4], ".attr", 5) == 0)
               {
               optarg[i - 5] = 0;
               }

            if (strchr(cwd_str,'\\') || strchr(optarg,'\\'))
               pathsep = '\\';
            sprintf(windowTitle,"MDDisplay: %s%c%s",cwd_str,pathsep,optarg);
            sprintf(bondFileName, "%s.bnd", optarg);
            sprintf(coordFileName, "%s.cor", optarg);
            sprintf(attribFileName, "%s.attr", optarg);
            sprintf(surfaceFileName, "%s.srf", optarg);
            sprintf(sesFileName, "%s.se3", optarg);
            sprintf(setFileName, "%s", optarg);
            }
            break;

         case 'r':
            recoverFlag = 1;
            break;

#ifdef REVIEW_IN_CASE_WE_DO_DISPLAY_LISTS_AGAIN
         case 'x':
            useDisplayLists = 1;
            break;
#endif

         case '?':
            mddisplay_usage_message();
            exit(EXIT_SUCCESS);
      }

   // Without all three files on the command line (through -n or
   // the long way, we must stop.
   if (coordFileName[0] == 0 || bondFileName[0] == 0 || frames == 0)
      {
      mddisplay_usage_message();

      exit(EXIT_FAILURE);
      }

   // Open the files, read the data.
   readparms_and_allocate_globals(recoverFlag);

   // Read in the coordinates into memory
   drawimage_prepare();

   // Display where we left off.  Or, if -r specified, where
   // we crashed.
   drawimage_initview(recoverFlag);

   // See above...
   // Let's Give ourselves a chance to report things 
   // in cases where the user does not exit with 'Q' or we
   // have a crash of some kind.
   atexit(bye);

#ifdef OLD_GL_METHOD
   drawimage();
#endif

   // Crank up the GLUT 3D library
   glutInit(&argc,argv);

   // Prepare dispui.c
   dispuiInit(
         movieViewBoxSize,
         0,
         picoSecondsInitialOffset,
         picoSecondsPerFrame,
         quad_stereo_flag);

   // When an atom is picked, simply turn its label on or off
   // as the default action
   labelFlag = 1;
   AtomPickedFunc = LabelOrUNLabelAtomPicked;

   // Let's show the movie.   
   dispuiMain(windowTitle);

   // The following code are nice ideas - but GLUT never returns to here...
   // We 'exit' above when QUIT is wrapped up - and that exit()s immediately.
   drawimage_close();

   return EXIT_SUCCESS;
}

