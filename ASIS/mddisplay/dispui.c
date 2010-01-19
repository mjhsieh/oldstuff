/*
 *  MD Display Version 3.0
 *
 *  dispui.c  user interface management
 *
 *  This module, dispui.c, manges the majority of the user interface components
 *  seen in MD Display.  
 *
 *  When compiled with the -DDISPUI_TESTMODE flag, this module can be used
 *  as a generic movie viewer, as shown in the C/C++ Users Journal
 *  article: A 3D Stereo Movie Viewer with GLUT, Sept 2003
 *
 *  dispui.c calls GLUT library routines to manage windows as well as
 *  keystroke and mouse events.  dispui.c provides complete control of
 *  scene scaling, translation, rotation, and frame update.
 *  This module knows very little about the molecular movie that it displays,
 *  and depends on display.c functions to do all the movie frame rendering.
 *
 *  As shown in the CUJ article, it is quite easy to drop in an alternate
 *  movie rendering application, and take advantage of the user interface
 *  requirements.
 *
 *  Copyright (C) 2002 by Chris Moth and Terry Lybrand
 *
 * Please cite:
 *    A 3D Stereo Movie Viewer with GLUT
 *    C/C++ Users Journal, Sept 2003
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
 *  If this code is used for Molecular Dynamics visualization,
 *  Please cite:       
 *  
 *    Callahan, T.J., Swanson, E., Lybrand, T.P.: 
 *       MD Display: An interactive graphics program for visualization
 *       of molecular dynamics trajectories, 
 *    Journal of Molecular Graphics, 1996, 14:39-41
 *
*/

// ANSI C standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

// Microsoft Visual C++ 6.0 lacks M_PI
// Remove these lines if you are using a later compiler
// or they give you trouble
#ifdef _MSC_VER
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif
#include <assert.h>
#include <ctype.h>

// GLUT window management library header
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// Function prototypes for this and other MDDisplay-specific modules.
#ifdef DISPUI_TESTMODE
#include "cujmovie.h"
#else
#include "mddisplay.h"
#include "displib.h"
#endif
#include "dispui.h"


#ifdef __TURBOC__
// Without this pragma, Borland's bcc32 compiler
// will complain about required parameters in GLUT callback functions.
#pragma warn -par
#endif

// Globals where a text status message can be stored until time out.
char   dispuiStatusMessage[80];
time_t dispuiStatusMessageExpiration;

// The selection of legacy split-buffer stereo on SGI platforms
// is a bit convoluted - and aided by some odd variables.
int SGIStereoRequested = 0;
int SGIStereoOffToFullScreenRequested = 0;
int quadStereoAvailable; // True if passed in via dispuiInit()

// A global flag that answers the question "are we in full screen mode"
// This variable is set when the user presses <u> or clicks 
// in the [Full Screen] menu box
int dispuiFullScreenFlag;


// SGI_STEREO is defined on the cc command line 
// in the MD Display makefile for SGI
// It is the only major platform-dependent aspect of this code, and allows
// MD Display to display stereo on legacy SGI systems with 
// the original SGI split-buffer stereo implementation.

#ifdef SGI_STEREO
// With SGI legacy split-buffer stereo, 
// special calls are required to select 
// the two eye buffers.  Here they are:
#include "fullscreen_stereo.h"

// This value is the number of rows lost in the split buffer mode of operation
const int VPSEP = 40; 
#endif


// As the program evolves, it would be nice to have all colors of user
// interface objects come from a central table in the source code
// for easy tailoring of the color schemes.  Here is a bad start towards
// that goal.

GLdouble menuBoxColor[] = 
   {100.0/255.0,60.0/255.0,0.0}; // brownish
GLdouble scalingBarColor[] = 
   {75.0/255.0,75.0/255.0,150.0/255.0,0.0}; // Dark Blue

// Cursor style that is used to select items in the movie
int movieCursor = GLUT_CURSOR_CROSSHAIR;

// Text font used for keyboard input (needs to be fixed width for
// cursor to look good (I think)
#define keyboardTextFont GLUT_BITMAP_9_BY_15



// When the user initiates rotation operations with the mouse, 
// a blue disk appears on the screen to help cue the user.  
// That disk is stored in the global movieTrackBallCircleQuadric
GLUquadricObj* movieTrackBallCircleQuadric;

void SetupMovieTrackBallCircleQuadric(void)
{
   assert(movieTrackBallCircleQuadric == 0); // Only called once!
	glClearColor(0.0f,0.0f,0.0f,0.0f);


	movieTrackBallCircleQuadric = gluNewQuadric();
	assert(movieTrackBallCircleQuadric != 0);
	gluQuadricDrawStyle(movieTrackBallCircleQuadric,GLU_SILHOUETTE);
	gluQuadricNormals(movieTrackBallCircleQuadric,GLU_NONE);
	gluQuadricTexture(movieTrackBallCircleQuadric,GL_FALSE);
	assert(glGetError() == GL_NO_ERROR);
}

void DeleteMovieTrackBallCircleQuadric(void)
{
	assert(movieTrackBallCircleQuadric != 0);
	gluDeleteQuadric(movieTrackBallCircleQuadric);
	movieTrackBallCircleQuadric = 0;
}

// In GLUT, windows and subwindows are assigned a unique integer
// "handle" when they are created.  It is easiest to store these 
// handles in global variables.

// The Master GLUT window is just a backdrop for the subwindows.  
// No drawing is done into the Master Window.
int masterWindow;	
// Forward declaration - see below
void MasterCreateOrChangeAllSubwindows(void); 

// handle to GLUT help Window
int helpWindow = -1;

// handle to GLUT subWindow displayed at screen bottom - 
// We draw the the scale, speed, and hyper controls on this
// subwindow
int bottomWindow = -1;
int bottomWindowHeight;


// Subwindow containing menu bottons, and misc graphics
int leftWindow = -1;		

// It is cheesy to have global variables for the dimensioning of a GLUT
// window.  Many other sizing decisions are made from the m x n pixel
// dimensions of the leftWindow.  Apologies, for these globals.
int leftWindowHeight;
int leftWindowWidth;

// One of the menu buttons, "static" is actually overlayed by a subwindow - so
// that a GLUT popup menu ("First", "Average", "None") can appear 
// when the button is clicked on.
int staticMenuButtonWindow = -1; 
int staticMenuIdentifier; // Menu identifier for the little statuc popup window

// Subwindow where the "movie" is displayed
int movieWindow = -1;


// This struct keeps up with user selected
// origin, rotation, view, etc of the running movie.
// It is declared in dispui.h
TMovieWindowControl movie;

void dispuiTranslation(float xTranslation,float yTranslation)
{
   movie.xTranslation = xTranslation;
   movie.yTranslation = yTranslation;
}

void dispuiCenterOfRotation(float x,float y,float z)
{
   movie.centerOfRotation[0] = x;
   movie.centerOfRotation[1] = y;
   movie.centerOfRotation[2] = z;
}

void dispuiInitRotationMatrix(void)
{
   int i;
   for (i=0;i<16;i++)
         {
 			movie.rotationMatrix[i] = 0.0;
         }
      // Yes - this is very tedious - but LoadIdentity does not seem
      // to be working until we do more openGL stuff downstream...
      movie.rotationMatrix[0] =
         movie.rotationMatrix[5] =
            movie.rotationMatrix[10] =
               movie.rotationMatrix[15] = 1.0;
}               



// Again, apologies for making these variables global - but they are used
// in many routines.  These variables are set in MovieChangeSizeFunc()
int movieWindowWidth;  
int movieWindowHeight;
double movieWindowAspectRatio; // w/h computed by ChangeSize

// MD Display 2.x was hard coded to the 1280 x 1024 resolution
// of earlier SGI hardware.  This "FixupFactor" allowed more of the old code
// to port directly over - and I scale drawing into the new size by this factor.
// 
double movieWindowFixupFactor; 

// Bottom left corner SubWindow where Clipping control is displayed
int clipWindow = -1;				
// Next frame to display in Clip Window... 
// This global keeps the clip image it in sync with main window.
int clipWindowFrameToDraw = -1; 

// The Palette subwindow allows selecting of atom coloring from color squres.
// It appears when the <color> menu button is clicked on or 'c' is typed.
void PaletteWindowCreateOrChange(void);
int paletteWindow = -1;		// Palette selection window, cranked up after menu option

// Cursor style when mouse pointer enters palette window
int paletteCursor= GLUT_CURSOR_LEFT_ARROW;   

#ifndef DISPUI_TESTMODE
// GLUT Window (NOT a subWindow) where Ramachandron plot is displayed
int ramaWindow = -1;
int ramaWindowLastFrameDisplayed;
int ramaWindowFrameToDraw = -1; // Next frame to display in Rama Window... Keeps in synch with main window.
void RamaWindowToggle(void);
#endif

// Forward declaration of 
// Function which is called when <?> key is pressed or help menu clicked 
void HelpWindowToggle(void);

// When scaling with the '[' and ']' keys or the bar, these constants set a reasonable range
const float movieMaxScale = 20.0f;
const float movieMinScale =  0.2f;


// Size in cubic pixels of the 3D movie image.  
// This is passed in at dispuiInit - and is
// constant throughout the program session
static int  movieImageSize; 

#ifndef DISPUI_TESTMODE                
// These two variables probably should not even be in dispui.c as they 
// are movie dependent.
// They are read to display the current frame number - 
// in "picoseconds" of simulation time
static float moviePicoSecondsInitialOffset;
#endif
static float moviePicoSecondsPerFrame;


// dispui.c has an abstracts notion of "frame number" but knows nothing about
// the contents of the movie window.
static int movieFrame;

// Sometimes, external modules need to know the frame in view.  They call
// this function.
int dispuiMovieFrame(void)
{
	return movieFrame;
}


// The speed of the display can be regulated using GLUT Timer
// functions.  These variables control the maximum slowness
// as well as default delay between frames.

const unsigned int movieMaxMsecs = 1000;
// On APPLE, it seems like the minimum value we can take here is 1
#ifdef __APPLE__
#define movieMinMsecs 1
#else
#define movieMinMsecs 0
#endif
// Time between frames using TimerFunc
unsigned int movieMsecsBetweenFrames = movieMinMsecs; 

// There is only one function call hook in GLUT for timer functions.  
// So, when a new step interval (i.e. movie frame rate) is selected, 
// the timerKey is incremented.  Pending timer Callbacks associated
// with other keys are thereby ignored when they are trapped.

int movieTimerKey = 0; 
const int maxMovieTimerKey = 10000; // Roll over if key gets this large

// There is a separate timer key for "flip text cursor" in keyboard input box.
// Each time the user keys a letter, a new key must be issued as the cursor
// should be active again for the right look and feel
#define minKeyboardTimerKey (10000+100)
#define maxKeyboardTimerKey (minKeyboardTimerKey + 1000)

int keyboardTimerKey = minKeyboardTimerKey;
void (* keyboardInputCompleteCallback)(int lastKey,const char* inputString);

// Default is to run movie Forwards.  Only other possibility is -1
int movieDirection = 1;	
int movieStepMode = 0;	// Default is to run continuously, via glutTimerFunc

// Add this to the frame each time we display a new one!  When the "Hyper"
// option is clicked, this varible is increased
int movieFrameIncrement = 1;	

// Set true during user initiated rotation operations.
// The user sees the blue disc while rotating
int movieTrackBallCircleVisible = 0;

// When the user presses '[' or ']' the scale bar must be refreshed to
// show new position.  Hence the forward declaration.
void RefreshScale(void);

// One enum for each (brown) menu box seen in the leftWindow
// The menu options are tied to the molecular movie
// display.c capabilities.
// Perhaps this should be decoupled in future.
enum {nHelp,
		nStatic,		
		nHide,
#ifndef DISPUI_TESTMODE
		nRemove,
#endif
		nOrigin,
#ifndef DISPUI_TESTMODE
		nDistanceOrAngleOrTorsion,
		nRama,
#endif
		nColors,
#ifndef DISPUI_TESTMODE
		nHBond,
		nUnLabelOrLabel,
#endif
		nFullScreen,
		nDisplay,
#ifndef DISPUI_TESTMODE
		nFiltering,
		nPDB,
#endif
		nClipLock,
		nClipping,
#ifndef DISPUI_TESTMODE
		nQuit,
#endif
		nMenuItemCount};

// Here are the descriptive strings that are displayed in each menu box.
// The \001 (ctrl-a) tells my display code to underline the next letter
// and thus remind the user of the keyboard equivalent.
const char* menuStrings[] = {
	"\001? Help",
	"Static",
	"\001hide",
#ifndef DISPUI_TESTMODE
	"\001remove",
#endif
	"\001Origin",
#ifndef DISPUI_TESTMODE
	"- | > | Z",	// Note 3 submenus
	"\001Rama",
#endif
	"\001colors",
#ifndef DISPUI_TESTMODE
	"\001H-bond",
	"\001Un|\001Label",		// Note 2 submenus
#endif
	"F\001ull Screen",
	"\001disp \001F\0011-\001F\0014",
#ifndef DISPUI_TESTMODE
	"\001Filter",
	"\001PDB Dump",
#endif
	"Clip Lock",
	"Clipping",
#ifndef DISPUI_TESTMODE
	"\001Quit"
#endif
};


// Note that some menu boxes actually have two virtual "submenu" boxes
// These are separated by a vertical bar.  So, there is a bit of management
// to be done to figure out what menu the user wants, given the mouse click.
// The width, in pixels of the '|' character in the current font is core to 
// this:
int menuAsciiBarWidth;
#define nMaxSubMenus 3
int subMenus[nMenuItemCount];
int subMenuXTextpos[nMenuItemCount][nMaxSubMenus];


// These menu* variables hold font, placement, etc, data
// of the menuItems in the LeftWindow.
// The variables are set whenever the leftWindow is
// created - in the MasterChangeSize routine.

float menuFontHeight = 12;
void* menuFont;

int menuTop = 60;			// Y position of first (topmost) menu box.
int menuBoxHeight = 28;	// Height of each brown menu box
int menuBoxSpace = 35;	// Spacing from a box - to box below
int menuBoxWidth = 88;	// Width of the box
int menuBoxLeft = 2;		// Starting X coordinate of Menu box in leftWindow.

// When the user picks the (c)olor menu option, or types 'c', we present
// a pallete of colors.  Here are the blends of colors.  0=None 255=Max

#define PALETTE_COLORS 10

const GLubyte palette[PALETTE_COLORS][3] =
{
	{  0,  0,  0},	// 0: Black
	{255,255,255}, // 1: White
	{255,  0,  0}, // 2: Red
	{ 64, 64,255}, // 3: Blue
	{  0,255,255}, // 4: Aqua
	{  0,255,  0}, // 5: Green
	{255,  0,255}, // 6: Magenta
	{255,128,  0}, // 7: Orange
	{255,255,043}, // 8: Yellow
	{155,155,155}	// 9: Gray
};

// Onto each color square in the pallette, we'll write the
// name of the color.
const char* dispuiPaletteColorString(int paletteColor)
{
	static const char* color_descriptions[PALETTE_COLORS] =
	{	"invisible (remove)",
	   "white",
	   "red",
		"blue",
		"aqua",
		"green",
		"magenta",
		"orange",
		"yellow",
		"gray"};

	assert(paletteColor >= 0 && paletteColor < PALETTE_COLORS);
	return color_descriptions[paletteColor];
}

// Given input color # from 0 to PALETTE_COLORS-1
// set the color from the table above.
void dispuiColorByPalette(int paletteColor)
{
	assert(paletteColor >= 0 && paletteColor < PALETTE_COLORS);
	glColor3ubv(palette[paletteColor]);
}

// In SGI legacy stereo mode, we must set width=(2 x height) for the user
// to see nice "square" palette color blocks
int paletteElementHeight;
int paletteElementWidth;
int paletteVerticalCenter;


// When we output a text string in GLUT, we write to the current
// "rastor" position.  This is where the \001 CTRL-A underline
// indicators are taken into account.
void dispuiDrawStringAtCurrentRasterPosition(void* font,const char* s)
{
	GLfloat underline_position[4];
	int underline_next = 0;

   while (*s)
      {
		if (*s == '\001') // "Underlined next character"
			{
			glGetFloatv(GL_CURRENT_RASTER_POSITION,underline_position);
			underline_position[1] -= 2; // Drop the y value down a bit for underline
			underline_next = 1;
			}
		else
			{
			if (underline_next) // Then we draw a simple "underline" under this letter first
				{
				glRectf(underline_position[0],
						underline_position[1],
					underline_position[0]+ glutBitmapWidth(font,*s),
						underline_position[1]-1.0);
				underline_next = 0;
				}
	      glutBitmapCharacter(font, *s);
			}
		s++;
      }
}

// Position the raster cursor to x,y and call above.
void DrawString(float x,float y,void* font,const char* s)
{
   glRasterPos2f(x, y);

	dispuiDrawStringAtCurrentRasterPosition(font,s);
}

// For centering of strings, etc, we need to know the width
// of a stnrg in pixels (excluding the \001 underline indicators)
int StringWidth(void* font, const char* s)
{
	int retval = 0;
	while (*s)
		{
		if (*s != '\001') // Don't count underline markers
			retval += glutBitmapWidth(font,*s);
		s++;
		}

	return retval;
}


// GLUT has no provision for keyboard input
// So, I had to write a routine.  EVERY window and subwindow
// MUST set it's keyboard handler to the same routine,
//   MasterKeyboardFunc()
// which users this code.
// by some of this code.

int keyboardWindow = -1;
int keyboardInputActive;

int keyboardCursorOn=1;

char* keyboardPromptString;
char keyboardInputString[40];

// I wanted a nice input cursor.  The idea is that when the Timer
// expires, the cursor will flip.  That said, if a new keystroke is pressed
// we want the cursor to appear again, which is handled elsewhere.
void KeyboardTimerFunc(int keyValue)
{
	if ((keyboardWindow != -1) && (keyValue == keyboardTimerKey))
		{
		int currentWindow = glutGetWindow();
		// FlipCursor
		keyboardCursorOn = ! keyboardCursorOn;
		glutSetWindow(keyboardWindow);
		glutPostRedisplay();
		glutSetWindow(currentWindow);
		}

}

// Return 2 if the user has pressed F4 to go into stereo mode.
// Return 1 otherwise.
int MovieEyeCount(void)
{
	return (movie.imageView == nSGILegacyHardwareStereo ||
            movie.imageView == nQuadStereo) ? 2 : 1;
}

// Return 2 only if in SGI Stereo mode.  In SGI Legacy stereo mode, every subWindow needs
// to have left and right eyes drawn.  Newer quad stereo displays are smarter about this.
int SGIEyeCount(void)
{
	return (movie.imageView == nSGILegacyHardwareStereo ? 2 : 1);
}

// When drawing in stereo, you have to draw separately in the left and right
// eye buffers.  This routine gets the job done, whether SGI Legacy or Quad-buffer
// stereo is being used.
void dispuiSelectEye(int eye)
{
#ifdef SGI_STEREO
	if (movie.imageView == nSGILegacyHardwareStereo)
		{
		if (eye)
		   stereo_right_buffer();
		else
			stereo_left_buffer();
		}
	else 
#endif   
   if (movie.imageView == nQuadStereo)
      {
		if (eye)
			glDrawBuffer(GL_BACK_RIGHT);
		else
			glDrawBuffer(GL_BACK_LEFT);
      }
}

// Display the user's text input and the cursor (if keyboardCursorOn)
// Yes, everytime the cursor flips we draw all of this stuff!
void KeyboardDisplayFunc(void)
{
	int i;

	int input_strlen = strlen(keyboardInputString);

	int w = glutGet(GLUT_WINDOW_WIDTH);
	int h = glutGet(GLUT_WINDOW_HEIGHT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

// We want the keyboard input to be visible in both eyes 
// hence we draw in both buffers
for (i=0;i < SGIEyeCount();i++)
	{
	dispuiSelectEye(i);

	// Reset coordinate System to a simple orthographic 2D system
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);

   // Let's have a blue box
	glClearColor(0.0,0.0,1.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT);

   // With a nice white border around it.
	glColor3f(1.0,1.0,1.0);

	glBegin(GL_LINE_LOOP);
		glVertex2i(1,1);
		glVertex2i(w-1,1);
		glVertex2i(w-1,h-1);
		glVertex2i(1,h-1);
	glEnd();

#define nKeyboardWindowOffset 5

   // Draw the "prompt" text
	DrawString(nKeyboardWindowOffset,h-16,keyboardTextFont,keyboardPromptString);

   // Draw any input that the user has entered so far.
	if (input_strlen)
		DrawString(nKeyboardWindowOffset,h-37,keyboardTextFont,keyboardInputString);

   // Draw the underline cursor to the right of the keyboardInputString
	if (keyboardCursorOn)
		{
		int inputWidth = StringWidth(keyboardTextFont,keyboardInputString);
		glRecti(nKeyboardWindowOffset + inputWidth,h-39,inputWidth+15,h-38);
		}

	}
	glutSwapBuffers();

	// Set up cursor flip...
   // 400 msecs seems to feel right, relative to other graphics cursors we see all the time
	keyboardTimerKey++;
   if (keyboardTimerKey >= maxKeyboardTimerKey)
   	keyboardTimerKey = minKeyboardTimerKey;
	glutTimerFunc(400,KeyboardTimerFunc,keyboardTimerKey);
}

// When we set full screen, then (newFlag != 0)
// When we want to restore windowed screen, newFlag = 0
// This function can toggle eitherway. to be toggleable as much as possible
void MasterSetFullScreen(int newFlag)
{
	int current_window = glutGetWindow();
	static int lastMasterWindowX = -1;
	static int lastMasterWindowY = -1;
	static int lastMasterWindowWidth = -1;
	static int lastMasterWindowHeight = -1;

	glutSetWindow(masterWindow);

	if (newFlag) // Prior to setting full screen mode, save old window placement and size.
		{
		lastMasterWindowX = glutGet(GLUT_WINDOW_X);;
		lastMasterWindowY = glutGet(GLUT_WINDOW_Y);;
		lastMasterWindowWidth = glutGet(GLUT_WINDOW_WIDTH);
		lastMasterWindowHeight = glutGet(GLUT_WINDOW_HEIGHT);
#ifdef __APPLE_STILL_CANT_FULLSCREEN__
// glutFullScreen does many many things on OS/X.  So, I've included
// a workaround here - although desktop controls are not covered
      glutPositionWindow(0,0);
      glutReshapeWindow(glutGet(GLUT_SCREEN_WIDTH),glutGet(GLUT_SCREEN_HEIGHT));
#else
		glutFullScreen();
#endif
		}
	else // We are returning from full screen mode
		{ // i.e. ('u' pressed 2nd time) restore if possible
		int screen_h = glutGet(GLUT_SCREEN_HEIGHT);
		int screen_w = glutGet(GLUT_SCREEN_WIDTH);

		// If we don't know how to get back to "normal" go to 90% of full screen...
		if ((lastMasterWindowWidth == -1) ||
			 ((lastMasterWindowHeight > (int)(screen_h * .95)) && (lastMasterWindowWidth > (int)(screen_w * 0.95))))
			{
			glutPositionWindow((int)((.1 * screen_w)/2.0),		// Center left to right
								  (int)((.05 * screen_w)/2.0));	// Close to top - but not completely there
			glutReshapeWindow((int)(screen_w * .9),(int)(screen_h* .9));
			}
		else // We can simply restore position and size we had were before full screen change...
			{
			glutPositionWindow(lastMasterWindowX,lastMasterWindowY);
			glutReshapeWindow(lastMasterWindowWidth,lastMasterWindowHeight);
			}
		}

	glutSetWindow(current_window);
	dispuiFullScreenFlag = newFlag;
}

// Sometimes, such as in case of resizing, we want a redisplay of the 
// current movie frame without advancing the current frame number.
// We achieve this by advancing the new movieTimerKey and telling
// the moveie subWindow to refresh itself
// When the delayed timer settings come in, they will be simply ignored.
void MovieAdvanceTimerKeyAndPostRedisplay(void)
{
	int current_window = glutGetWindow();
	glutSetWindow(movieWindow);
	movieTimerKey++;	// All prior timer events will be ignored!
	if (movieTimerKey > maxMovieTimerKey)
			movieTimerKey = 0;
	glutPostRedisplay();    // This will cause refresh of the main movie subWindow
	glutSetWindow(current_window);
}

// Forward declarations of some needed functions.
void dispuiSetMovieFunctions(void); // set GLUT callbacks for the Movie Window
void MasterKeyboardFunc(unsigned char key,int x, int y);
void MasterSpecialFunc(int key, int x,int y);

// When the user wants to shift into or out of stereo, there is substantial
// subwindow management which must take place.
void MasterChangeImageViewMode(TImageView oldImageView,TImageView newImageView)
{
   if (oldImageView != newImageView) // If the user repeatedly presses Fn, then do nothing.
  		{
      // Otherwise, a new Fn key was pressed so action may be needed.
		int current_window = glutGetWindow();
      
      // If we were in full screen mode before SGI legacy stereo, then we want to revert
      // to full screen on exit from stereo.  This variable allows this to happen.
#ifdef SGI_STEREO
   	static int full_screen_flag_before_stereo;
#endif
      
      // If we are going into full Quad Stereo, or leaving Quad Stereo.
      // then we need to destroy the current movie Window, and remake it
      // in QuadStereo mode.
      if ((newImageView == nQuadStereo) || (oldImageView == nQuadStereo))
         {
         int x,y,w,h;
         int kx,ky,kw,kh;
         unsigned int quad_stereo_bits;
         int oldMovieWindow = movieWindow;
         glutSetWindow(movieWindow);
         
         x = glutGet(GLUT_WINDOW_X);
         y = glutGet(GLUT_WINDOW_Y);
         w = glutGet(GLUT_WINDOW_WIDTH);
         h = glutGet(GLUT_WINDOW_HEIGHT); 
         
         // The paletteWindow and keyboardWindow (subWindows on the movieWindow)
         // might be active.
         // So take great care to destroy and reinstate them properly!
         if (paletteWindow != -1)
            {
            if (current_window == paletteWindow)
               current_window = -1;
            glutDestroyWindow(paletteWindow);
            }
         if (keyboardWindow != -1)
            {
            glutSetWindow(keyboardWindow);
            kx = glutGet(GLUT_WINDOW_X);
            ky = glutGet(GLUT_WINDOW_Y);
            kw = glutGet(GLUT_WINDOW_WIDTH);
            kh = glutGet(GLUT_WINDOW_HEIGHT); 

            if (current_window == keyboardWindow)
               current_window = -1;
            glutDestroyWindow(keyboardWindow);
            }

         // Finally, destroy the movieWindow         
         if (current_window == movieWindow)
            current_window = -1;
         
			DeleteMovieTrackBallCircleQuadric();
         glutDestroyWindow(movieWindow);
         
         // GLUT_STEREO is part of the flag then if the -q flag was on the
         // command line AND user pressed F4.  Otherwise, we're going back to 
         // regular graphics mode (quad_stereo_bits = 0)
   		quad_stereo_bits = quadStereoAvailable && (newImageView == nQuadStereo) ? GLUT_STEREO : 0;
	   	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | quad_stereo_bits);
         movieWindow = glutCreateSubWindow(
                           masterWindow,
                           leftWindowWidth,
                           0,
                           w,
                           h);
	SetupMovieTrackBallCircleQuadric(); 

         
         if (oldMovieWindow == current_window)
            current_window = movieWindow;
         dispuiSetMovieFunctions();
         
         // Reinstate the palette subWindow if we had one before
         if (paletteWindow != -1)
            {
            paletteWindow = -1;
            PaletteWindowCreateOrChange();
            }
            
         // Reinstate the keyboard subWindow if we had one before
         if (keyboardWindow != -1)
            {
            keyboardWindow = glutCreateSubWindow(movieWindow,kx-x,ky-y,kw,kh);
          	glutKeyboardFunc(MasterKeyboardFunc);
         	glutSpecialFunc(MasterSpecialFunc);
   	      glutDisplayFunc(KeyboardDisplayFunc);
            }
         }
#ifdef SGI_STEREO
// With SGI Legacy stereo, we have a fundamental shift in 
// the video hardware to cope with.  Full screen is required
// in stereo mode.
  		else if (oldImageView == nSGILegacyHardwareStereo)
			{
         // We are coming out of stereo - back to normal
			if (full_screen_flag_before_stereo)
				{
				SGIStereoOffToFullScreenRequested = 1;
				glutSetWindow(masterWindow);
				glutReshapeWindow(1,1); // Trapped by change size func
				glutPositionWindow(19,32);
				}
			else
				{
  				stop_fullscreen_stereo();
				MasterSetFullScreen(0);
				}
         // In the process of the above, the Rama and help
         // windows will have been covered.  Restore them
         // now if they were active.
#ifndef DISPUI_TESTMODE
			if (ramaWindow != -1)
				{
				RamaWindowToggle();
				RamaWindowToggle();
				}
#endif
			if (helpWindow != -1)
				{
				HelpWindowToggle();
				HelpWindowToggle();
				}
			}
  		else
  		if (newImageView == nSGILegacyHardwareStereo)
			{
         // We are going into SGI stereo
         // Note whether we are currently in full screen mode.
		   full_screen_flag_before_stereo = dispuiFullScreenFlag;
			glutSetWindow(masterWindow);
			// Generate a call to ChangeSize
			MasterSetFullScreen(0);
			// Crzy stuff to force resize
			glutReshapeWindow(1,1); // Trapped by change size func
			glutPositionWindow(19,32);
			// During resize - this guy will go right into stereo
			SGIStereoRequested = 1;
         // In the process of the above, the Rama and help
         // windows will have been covered.  Restore them
         // now if they were active.
#ifndef DISPUI_TESTMODE
			if (ramaWindow != -1)
				{
				RamaWindowToggle();
				RamaWindowToggle();
				}
#endif
			if (helpWindow != -1)
				{
				HelpWindowToggle();
				HelpWindowToggle();
				}
		  	}
#endif         
      // If our current window was not deleted, refocus it now
      if (current_window != -1)
   		glutSetWindow(current_window);
      // We certainly need to redraw the movie frame.         
		MovieAdvanceTimerKeyAndPostRedisplay();
		}
}

// In case of an Fkey press, GLUT has a callback function, which
// we point to this routine below:
void MasterSpecialFunc(int key, int x,int y)
{
	TImageView oldImageView = movie.imageView;

	switch ( key )
		{
		case GLUT_KEY_F1:
			movie.imageView = nSingleImage;
			break;
		case GLUT_KEY_F2:
			movie.imageView = nPairOfImages;
			break;
		case GLUT_KEY_F3:
			movie.imageView = nThreeImages;
			break;
		case GLUT_KEY_F4:
#ifdef SGI_STEREO
         // Some SGI systems DO support quad_buffered stereo
         // We'll use quad buffer if -q was specified on the 
         // command line.
			movie.imageView = 
               quadStereoAvailable ? 
                  nQuadStereo : nSGILegacyHardwareStereo;
#else                  
			if (quadStereoAvailable) // i.e. was -q on the command line?
				movie.imageView = nQuadStereo;
			else
				{
            // Rather than letting system crash, we tell user about -q option and need
            // for hardware support.
				sprintf(dispuiStatusMessage,"Quadbuffered stereo requires the -q option and hardware support.");
				dispuiStatusMessageExpiration = time(NULL)+4;
				MovieAdvanceTimerKeyAndPostRedisplay();
				}
#endif                  

			break;

		case GLUT_KEY_F6:
#ifdef BUG_FOR_LATER_TOGGLES_DEPTH_CUE_IN_ORIGINAL
			movie.imageView = nSGILegacyHardwareStereo;
#endif
			break;

		default:
				break;

		}

	MasterChangeImageViewMode(oldImageView,movie.imageView);
}

// When the user clicks on the display menu button, 
// or presses 'd' we advance the view to split, tri, stereo, 
// and back to normal again.
void MovieAdvanceImageView(void)
{
	TImageView oldImageView = movie.imageView;
	switch (movie.imageView)
		{
		case nSingleImage:
			movie.imageView = nPairOfImages;
			break;
		case nPairOfImages:
			movie.imageView = nThreeImages;
			break;
		case nThreeImages:
#ifdef SGI_STEREO      
			movie.imageView = 
               quadStereoAvailable ? 
                  nQuadStereo : nSGILegacyHardwareStereo;
#else
			movie.imageView = 
               quadStereoAvailable ? 
                  nQuadStereo : nSingleImage;
#endif      
                  
			break;
         
		case nQuadStereo:
		case nSGILegacyHardwareStereo:
		default:	
			movie.imageView = nSingleImage;
			break;
		}

//	MovieAdvanceTimerKeyAndPostRedisplay();
	MasterChangeImageViewMode(oldImageView,movie.imageView);
}

void ClipWindowToggle(void);

// This function deals with user pressing 'h' to hide menus (maximize movie display)
void MovieHideLeftAndBottom(int hideFlag)
{
	static int hadClipWindow = 0;
   
   // Set the global to indicate the new hidden (or not) state
	movie.hideLeftAndBottom = hideFlag;
   
   // If we are hiding, record whether or not we have a clip tool
	if (hideFlag)
		hadClipWindow = (clipWindow != -1);
   
   // With the flag
	MasterCreateOrChangeAllSubwindows();

   // If we are unhiding (second press of 'h' or whatever)
   // then restore the clip tool if we had it before hiding.
	if ((! hideFlag) && (hadClipWindow))
	   ClipWindowToggle();	
}

// You can stuff key strokes into the input buffer to simulate user input.
// This can be a convenient way of organizing your software- so that
// clicking of the mouse results in keyboard input to the text box.

void KeyboardStuff(const char* stuffString)
{
if (keyboardWindow != -1)
	{
	char* s = keyboardInputString;
	*s = 0;

	while (*stuffString)
		{
		char key = *stuffString;
		if (key == 13) // Then it is a CR
			{
			// This is "your" function that we'll call when user press <ENTER> or <ESC>ape
			keyboardInputCompleteCallback(key,keyboardInputString);
			// Keep stuffing
			s = keyboardInputString;
			}
		else if (key == 27) // Then an ESCAPE is being stuffed! CANCEL this input!
			{
			glutDestroyWindow(keyboardWindow);
			keyboardWindow = -1;
			// We must do the callback so it can update
			// other things!
			keyboardInputCompleteCallback(key,keyboardInputString);
			free(keyboardPromptString);
			return;
			}
		else
			{
			*s++ = key;
			}
		stuffString++;
		*s = 0;
		}
	}
}

// All windows and subwindows in this application must point their
// Keyboard handler to this function with glutKeyboardFunc(MasterKeyboardFucn).
// When keystrokes are received (regardless of which subwindow is active), they
// are processed as text input or as menu equivalents.
void MasterKeyboardFunc(unsigned char key,int x, int y)
{
	// Always toggle help if ? is pressed!
   // Right: '?' can not be input by the user in a text string.
	if (key == '?')
		HelpWindowToggle();
	else
	if (keyboardWindow == -1) // No active input, keystrokes as UI controls 
		{
		if (isdigit(key))
			displayNewColorRequested(key-'0');
      else
      
		switch (key)
			{
			case 'h': // Toggle visibility of user interface components
                   // Clearing these components with the staticButton is
                   // is being used, however, is not doable.
				if (glutGetWindow() != staticMenuButtonWindow)
					MovieHideLeftAndBottom(! movie.hideLeftAndBottom);
				break;

			case 'u': // Toggle Full Screen display
				if (movie.imageView != nSGILegacyHardwareStereo)
					MasterSetFullScreen(! dispuiFullScreenFlag);
				break;

			case 'd': // Shift from F1,F2,F3,F4 and F1 again display modes.
				MovieAdvanceImageView();
				break;

			case 's': // Cause smearing of movie images.
				movie.smear++;
				break;

			case 'S': // Turn off smearing
				movie.smear = 0;
				break;

			case '[': // Zoom out
				movie.scaleFactor /= 1.05;
				RefreshScale();
				MovieAdvanceTimerKeyAndPostRedisplay();
			break;
			case ']': // Zoom in
				movie.scaleFactor *= 1.05;
				RefreshScale();
				MovieAdvanceTimerKeyAndPostRedisplay();
				break;

#ifndef DISPUI_TESTMODE
			case 'R': // Toggle Ramachandron Window
				RamaWindowToggle();
				break;
   		case 'r':
	   		displayRemoveRequested();
		   	break;
#endif

			case 'c': // Toggle Color input and Palette Window
			case 'C':
				if (paletteWindow != -1)
					{
					glutDestroyWindow(paletteWindow);
					paletteWindow = -1;
					}
				else
					{
					PaletteWindowCreateOrChange();
					}
				break;
            

   		case 'o':
   		case 'O':
   			displayCenterOfRotationRequested();
   			break;

			default:
				// If dispui.c does not process the keystroke, hand it off to
            // display.c  Maybe display.c will care about it and do something
            // with it.
				displayKeyboardFunc(key); // Let "display" know we have an uninteresting key
				break;

			}
		}
	else // if we get here, we are in the middle of active text string input 
		{
		int inputChanged = 0;

		int input_strlen = strlen(keyboardInputString);

      // If ASCII text entered, add it to the keyboard string
		if (isprint(key) && input_strlen < sizeof(keyboardInputString)-1)
			{
//			if (islower(key))
//				key = toupper(key);

			keyboardInputString[input_strlen++] = key;
			keyboardInputString[input_strlen] = 0;
			inputChanged = 1;
			}

		// If Backspace entered, wipe out last keystroke.
		if ((key == 8) && // ASCII Backspace
			 (input_strlen > 0))
			{
			input_strlen--;
			keyboardInputString[input_strlen] = 0;
			inputChanged = 1;
			}

      // We are "done" if the user presses <ENTER> or <ESC>
		if ((key == 13) || // ASCII Carriage Return
			 (key == 27))   // ASCII Escape
			{
			glutDestroyWindow(keyboardWindow);
			keyboardWindow = -1;
			free(keyboardPromptString);
			// This is "your" function that we'll call when user press <ENTER> or <ESC>ape
			keyboardInputCompleteCallback(key,keyboardInputString);
			}
		else
		if (inputChanged)
			{
         // We need to redraw the entire text input box when the input
         // has changed througha keystroke or Backspace key.
			int currentWindow = glutGetWindow();

			// Always Turn cursor On when we get a new key or input change!
         // Otherwise, it just does not look right when compared with
         // other GUI text entry routines.
			keyboardCursorOn = 1;
			keyboardTimerKey++; // Ignore pending cursor flip timers
		   if (keyboardTimerKey >= maxKeyboardTimerKey)
		   	keyboardTimerKey = minKeyboardTimerKey;
			glutSetWindow(keyboardWindow);
			glutPostRedisplay();
			glutSetWindow(currentWindow);
			}
		}
}

void MovieDisplayFunc(void);
void MovieChangeSizeFunc(int newWidth,int newHeight);
void MovieMouseFunc(int button, int state,int x, int y);
void MovieMotionFunc(int x,int y);

// Setup GLUT callbacks for the Movie Window
// This has to be done a couple of different places.  So,
// we just do it here to save code.
void dispuiSetMovieFunctions(void)
{
		// This cursor is nice for picking things in the movie
		glutSetCursor(movieCursor = GLUT_CURSOR_CROSSHAIR);

      // As with all other windows and subwindows, any keystrokes pressed
      // should be processed by the Master window
		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);

      // Specify callbacks in case of movie window refresh or size change
		glutDisplayFunc(MovieDisplayFunc);
		glutReshapeFunc(MovieChangeSizeFunc);

      // Specify callbacks in case of Mouse clicks or moves
		glutMouseFunc(MovieMouseFunc);
		glutMotionFunc(MovieMotionFunc);
}


void KeyboardInputStart(
	// Text string prompt for user input
	const char* _promptString,
   // "your" Function which is to be caled when user input is complete (13=ENTER)
   // or cancelled (27 = ESC)
   void (* _InputCompleteCallback)(int lastKey,const char* inputString))
{
	if (keyboardWindow != -1) // Then, we have an active input going on from another caller... Clean up!
		{
		// If the user has (accidentally) reselected the same UI component - don't even bother to destroy...
      // Just let user proceed with current input
		if ((keyboardInputCompleteCallback == _InputCompleteCallback) &&
			 (! strcmp(keyboardPromptString,_promptString)))
			return;

      // Otherwise, we need to act like the user pressed <ESC> to cancel
      // the current input - because now he wants to input something else.
		glutDestroyWindow(keyboardWindow);
		keyboardWindow = -1;
		free(keyboardPromptString);

      // Pass <ESC> (27) to the Callback function - which will then likely
      // treat the input string as being cancelled.
		keyboardInputCompleteCallback(27,"");
		}

   // Start with an empty input string and the cursor in the on state
	keyboardInputString[0] = 0;
	keyboardCursorOn = 1;

	// This is "your" function that we'll call when the user finally
   // presses <ENTER> or <ESC>ape to terminate string input
	keyboardInputCompleteCallback = _InputCompleteCallback;
	keyboardPromptString = strdup(_promptString);

   // Compute the least obstrusive size for the text input window, given
   // the input string length, and display that. 
	{
	int prompt_string_width = StringWidth(keyboardTextFont,keyboardPromptString)+40;
	int input_string_width = StringWidth(keyboardTextFont," ")*sizeof(keyboardInputString)+40;
	int keyboard_window_width = (prompt_string_width > input_string_width) ?
					prompt_string_width : input_string_width;

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	keyboardWindow = glutCreateSubWindow(
												movieWindow,
												20,
												10,
												keyboard_window_width,
												50);
	}

	glutKeyboardFunc(MasterKeyboardFunc);
	glutSpecialFunc(MasterSpecialFunc);

	glutDisplayFunc(KeyboardDisplayFunc);
}

// Setting a default input string can sometimes save end user keystrokes
// (such as when a .pdb is output and the filename default is constructed
// from the current frame number
void KeyboardInputSetDefault(const char* defaultInputString)
{
	// Copy over the supplied default, taking care to not exceed
   // the sizeof the static keyboardInputString
	int slen = strlen(defaultInputString);
	if (slen >= sizeof(keyboardInputString))
		slen = sizeof(keyboardInputString)-1;
	if (slen)
		memcpy(keyboardInputString,defaultInputString,slen);
	keyboardInputString[slen] = 0;
}

// These routines manage the speed and scaling bars.  The vertical line should be
// positioned to give the user a feel for having made a selection along a
// range of possibilities.

// Given a range of values, where should we draw the vertical line?
int PercentageBarXpos(float value,float minValue,float maxValue,int barLeftX,int barWidth)
{
	float pctg = (value - minValue) / (maxValue-minValue);

	return barLeftX + 2 + pctg*(barWidth - 6);
}

// Draw the vertical line
void DrawPercentageBar(float value,float minValue,float maxValue,
							int barLeftX,int barBottomY,
							int barRightX,int barTopY)
{
	// Get the correct location for the vertical bar
	int barXpos = PercentageBarXpos(value,minValue,maxValue,barLeftX,barRightX-barLeftX+1);

   // Now draw it 3 "pixels (units)" wide
	glRecti(barXpos-1,barBottomY,
					barXpos+1,barTopY);
}

// When the user clicks in a speed or scaling control bar, we
// need to know how far
// into the bar the new mouse press lies - as a value from 0.0 (extreme left)
// to 1.0 (extreme right)
// The '2' and '4' factor that the left and right edges of the bar are "off limits"
// and always give the user a value of 0.0 or 1.0.

float BarPercentage(int mouseX,int barLeftX,int barRightX)
{

	float bar_pctg = (float)(mouseX - barLeftX - 2)/
		 (float)(barRightX - barLeftX - 4);

	if (bar_pctg < 0.0)
			bar_pctg = 0.0;
	if (bar_pctg > 1.0)
			bar_pctg = 1.0;

	return bar_pctg;
}

// Convert 0.0 through 1.0 of the above function into a value ranging
// from minValue to maxValue
float BarValue(int mouseX,float minValue,float maxValue,int barLeftX,int barRightX)
{
	return minValue + BarPercentage(mouseX,barLeftX,barRightX) * (maxValue-minValue);
}

#ifndef DISPUI_TESTMODE
void RamaDisplayFunc(void)
{
	int eye;

	for (eye = 0;eye < SGIEyeCount();eye++)
		{
		dispuiSelectEye(eye);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	displayRenderRamaFrame(ramaWindowFrameToDraw);
		}

	glutSwapBuffers();
}
#endif

#ifndef DISPUI_TESTMODE
void RamaMouseFunc(int button,int state,int mouseX,int mouseY)
{
	if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN))
		{
		int current_window = glutGetWindow();

		glutSetWindow(ramaWindow);

		displaySelectRamaResidue(ramaWindowFrameToDraw,mouseX,mouseY);

		glutSetWindow(current_window);
		}
}
#endif

#ifndef DISPUI_TESTMODE
void RamaChangeSizeFunc(int w,int h)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);
	glClearColor(0.0,0.0,0.0,1.0);
}
#endif

// If the Palette subWindow should change size, simply reset
// the coordinate system to the new size. 
void PaletteChangeSizeFunc(int w,int h)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);
}

void PaletteRectangle(int color,float x1,float y1,float x2,float y2)
{
		void* font = 0;
		int digit_char = color + '0';
		int font_height = 0;  // If window is very small, no digit_char display
		float letter_color = 1.0f; // Default is white letters on palette

		// But we need black letters on bright colors
		if (strchr("1458",digit_char))
				letter_color = 0.0f;

      // Set the OpenGL RGB color from our 0-9 integer (see above)
		dispuiColorByPalette(color);

      // Draw the palette rectangle
		glRectf(x1,y1,x2,y2);

      // If there is room, display the digit_char
		if (y1 - y2 > 14)
			{
			font = GLUT_BITMAP_HELVETICA_10;
			font_height = 10;
			}
      // and pick a large font if there is room for that.
		if (y1 - y2 > 18)
			{
			font = GLUT_BITMAP_HELVETICA_12;
			font_height = 12;
			}

		if (font_height)
			{
			char digit_str[4];
			glColor3f(letter_color,letter_color,letter_color);

         // The \001 is our internal signal to DrawString (below)
         // that the following digit_char should be underlined
			sprintf(digit_str,"\001%c",digit_char);

         // Draw the digit_char nicely centered in the color square
			DrawString(x2 + (x1-x2-glutBitmapWidth(font,digit_char))/2,y2 + (y1-y2-font_height)/2,
				font,digit_str);
			}
}

void PaletteDisplayFunc(void)
{
	int color;
	int eye;

for (eye = 0;eye < SGIEyeCount();eye++)
	{
	dispuiSelectEye(eye);

	glClearColor(1.0,1.0,1.0,0.0);
	glClear(GL_COLOR_BUFFER_BIT); // Should have White everywhere to start...

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(0.0,(float)paletteElementHeight,0.0); // Go up by one element!

   // Draw the color squares on the palette
	for (color=0;color<PALETTE_COLORS;color+=2)
		{
		PaletteRectangle(color,
					1.0,paletteElementHeight*5-1-paletteElementHeight*(color/2),
					paletteElementWidth+1,paletteElementHeight*5-1-paletteElementHeight*(color/2)-paletteElementHeight);

		PaletteRectangle(color+1,
				paletteElementWidth+1,paletteElementHeight*5-1-paletteElementHeight*(color/2),
				paletteElementWidth*2-1,paletteElementHeight*5-1-paletteElementHeight*(color/2)-paletteElementHeight);
		}

	glPopMatrix(); // undo the translate above

	// Create teh bottom area where hte "by atom" option can be selected
	glColor3ub(0,0,0); // Black Box
	glRectf(1.0,1.0,
			2*paletteElementWidth-1,paletteElementHeight-2);
#ifndef DISPUI_TESTMODE
	{
	static const char* paletteByAtomString = "\001by atom";
   int byAtomStringWidth;

	byAtomStringWidth = StringWidth(GLUT_BITMAP_HELVETICA_12,paletteByAtomString);
	glColor3ub(255,255,255); // White text
	DrawString(paletteElementWidth-byAtomStringWidth/2,1+paletteElementHeight/2 - 12/2,GLUT_BITMAP_HELVETICA_12,paletteByAtomString); // Center up string
   }
#endif

	}
	glutSwapBuffers();
}

// When the user clicks in palette square, inform the user program
// that a new color has been requested.
void PaletteMouseFunc(int button,int state,int x,int y)
{
	int w = glutGet(GLUT_WINDOW_WIDTH)-2;
	int line0to4;
	int side0to1;

	if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_UP))
		{
      // The palette is arranged in 2 columns.  Figure out the
      // row (line0to4) and column (side0to1) and report back the color
      // selected (0-9)

		line0to4 = y/paletteElementHeight;
		side0to1 = (x > w/2) ? 1: 0;

		glutDestroyWindow(paletteWindow);
		paletteWindow = -1;

#ifndef DISPUI_TESTMODE
		if (line0to4 == 5) // User selected "by Atom"
			displayNewColorRequested(99); // Call display to do the real work...
		else
#endif
		// Let the display code know that the user has requested a color change      
		displayNewColorRequested(line0to4*2+side0to1);
		}
}

void PaletteWindowCreateOrChange(void)
{
	// Routine to size color palette subWindow of movie Window
	// The Idea here is that vertical midpoint of palette window is aligned with bottom of COLOR menu button
	// And, each palette item is a SQUARE with height 1.5xmenuBoxSpace
	// If nothing else - this looks good and is like original display
	// This code fails if we change leftWindow to _not_ reach to top of screen
	paletteElementHeight  = paletteElementWidth = menuBoxSpace*1.5;
	if (movie.imageView == nSGILegacyHardwareStereo)
		paletteElementWidth *= 2;
	paletteVerticalCenter = menuTop + nColors*menuBoxSpace + menuBoxHeight;

	if (paletteWindow == -1) // Then Create the new palette subWindow
		{
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		paletteWindow = glutCreateSubWindow(
												movieWindow,
												0,
												paletteVerticalCenter - 3*paletteElementHeight,
												paletteElementWidth*2,
												paletteElementHeight*6);

		// Tell GLUT what functions to call to respond to events                                    
		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);

		glutDisplayFunc(PaletteDisplayFunc);
		glutReshapeFunc(PaletteChangeSizeFunc);
		glutMouseFunc(PaletteMouseFunc);
		glutSetCursor(paletteCursor = GLUT_CURSOR_LEFT_ARROW);
		}
	else
		{
      // If the subWindow was already on the screen, just resize it.
      // This makes for a smoother U/I look.
		glutSetWindow(paletteWindow);
		glutPositionWindow(0,paletteVerticalCenter - 3*paletteElementHeight);
		glutReshapeWindow(paletteElementWidth*2,paletteElementHeight*6);
		}

}

// Simply increment or decrement the current movieFrame number
// Multiplying by movieDirection takes into account the "Hyper" control -
// which allows frames to be skipped.
// The Frame number wraps around past either end of the movie
void MovieNextFrame(void)
{
	movieFrame+= (movieFrameIncrement * movieDirection);

	if (movieDirection == 1) // If we are going Forward
		{
		if (movieFrame > displayMovieMaxFrame()) // And were past end
			movieFrame = displayMovieMinFrame();  // Then back to beginning
		}

	if (movieDirection == -1)	// If we are going backwards
		{
		if (movieFrame < displayMovieMinFrame()) // And we are before start
			movieFrame = displayMovieMaxFrame();  // wrap to end frame
		}
}

// When the Timer Event expires and this function is called,
// we advance the movieFrame number and generate a Redisplay event
// This is the core of the animation process.

// static int timerPostRedisplayInitiated = 1;
void MovieTimerFunc(int keyValue)
{
	// If the keyValue is a timer event to which we should
   // respond, and we are not in step mode then...
   if ((keyValue == movieTimerKey) && (! movieStepMode))
		{
      // get current focus
      int currentWindow = glutGetWindow();
      // advance the movieFrame
		MovieNextFrame();
      // Post a Redisplay event to trigger redrawing of the movie
      // window.  It will be redrawn with at the new movieFrame
		glutSetWindow(movieWindow);
		glutPostRedisplay();
		// timerPostRedisplayInitiated = 1;
		glutSetWindow(currentWindow);
		}
}

// This routine sets an hour glass (or whatever the GLUT implementation
// provides - beachball on Macintosh) to tell the user that the system
// is tied up in a calculation.
// The hour class is seen in all active subwindows.
void dispuiSetWaitCursors(void)
{
	int current_window = glutGetWindow();

	glutSetWindow(movieWindow);
	glutSetCursor(GLUT_CURSOR_WAIT);
	glutSetWindow(masterWindow);
	glutSetCursor(GLUT_CURSOR_WAIT);
	if (leftWindow != -1)
		{
		glutSetWindow(leftWindow);
		glutSetCursor(GLUT_CURSOR_WAIT);
		}

	if (clipWindow != -1)
		{
		glutSetWindow(clipWindow);
		glutSetCursor(GLUT_CURSOR_WAIT);
		}

	if (bottomWindow != -1)
		{
		glutSetWindow(bottomWindow);
		glutSetCursor(GLUT_CURSOR_WAIT);
		}
	glutSetWindow(current_window);
}

// When calculations are complete, and the hourglass is no longer needed
// call this function to restore the "normal" cursor.
void dispuiPopCursors(void)
{
	int current_window = glutGetWindow();

	glutSetWindow(movieWindow);
	glutSetCursor(movieCursor);
	glutSetWindow(masterWindow);
	glutSetCursor(GLUT_CURSOR_INHERIT);

	if (leftWindow != -1)
		{
		glutSetWindow(leftWindow);
		glutSetCursor(GLUT_CURSOR_INHERIT);
		}

	if (clipWindow != -1)
		{
		glutSetWindow(clipWindow);
		glutSetCursor(GLUT_CURSOR_INHERIT);
		}

	if (bottomWindow != -1)
		{
		glutSetWindow(bottomWindow);
		glutSetCursor(GLUT_CURSOR_INHERIT);
		}

	glutSetWindow(current_window);
}

// This code simply calls your supplied display* routines
// which make calls to OpenGL.  For each frame display,
// this routine is called once for "normal" mode,
// twice for stereo modes, and thrice for Tri mode.
// The transformation matrices must not be altered in
// your display* code.
void MovieDisplayCurrentFrame(GLenum renderMode)
{
	if (renderMode == GL_SELECT) // Then we are "rendering" for benefit of mouse click only...
		{
		displayRenderFrameInSelectMode(movieFrame);
		}
	else
	if (movie.smear == 0)
		{
		displayRenderFrameInPerspectiveBox(movieFrame,0);
		}
	else
		{
		int end_frame = movieFrame+movie.smear;
		int f;
		if (end_frame > displayMovieMaxFrame())
			end_frame = displayMovieMaxFrame();
		for (f=movieFrame;f<=end_frame;f++)
			{
			int smearFlag = f-movieFrame;
			glPushMatrix();
			displayRenderFrameInPerspectiveBox(f,smearFlag);
			glPopMatrix();
			}
		}
}


// If we are "Selecting" the main frame, then this function is called to set
// the sensitivity of the mouse to selecting our atom or object...
static void InitPickMatrix(
							 GLint mouseX,
							 GLint mouseY,
							 GLdouble selectionWidth,
							 GLdouble selectionHeight)
{
	GLint viewport[4];

	glGetIntegerv(GL_VIEWPORT,viewport);
   
   // See OpenGL references for more info on the important and complex
   // functions which setup picking objects in 3D.
	gluPickMatrix(mouseX,mouseY,selectionWidth,selectionHeight,viewport);

	assert(glGetError() == GL_NO_ERROR);
}


/*------------------------------------------------------------
*
* The following was taken from Silicon Graphics 4DGifts and
* modified slightly:
* Copyright (c) 1991 Silicon Graphics, Inc.
*
* Permission to use, copy, modify, distribute, and sell this software and
* its documentation for any purpose is hereby granted without fee, provided
* that (i) the above copyright notices and this permission notice appear in
* all copies of the software and related documentation, and (ii) the name of
* Silicon Graphics may not be used in any advertising or
* publicity relating to the software without the specific, prior written
* permission of Silicon Graphics.
*/

// Also thanks to Eric Swanson.  Some inspiration was gained from
// his psshow module stereopersp.c

// There are many web sites which explain the offset required to place
// left and right eye images correctly.  Apologies for lack of detailed
// comments here.  

void StereoPerspective(float fovy,
						float aspect,
						float near_,
						float far_,
						float conv,
						float eye)
{
	float left, right, top, bottom;
	float gltan;
	float fovy_radians = fovy * M_PI/180.0;

	if (conv == 0.0)
		conv = (near_ + far_)/2.0;

	top = tan(fovy_radians/2.0) * near_;
	bottom = -top;

	gltan = tan(fovy_radians/2.0 * aspect);

	left = (-gltan - eye / conv) * near_;
	right = (gltan - eye / conv) * near_;

	glFrustum(left, right, bottom, top, near_, far_);

	glTranslatef(-eye, 0.0, 0.0);
}

// In Tri (F3) mode where we display Front, side and top views, we need
// a central routine to give us viewport boundaries.  Here it is:
enum {nNormalImage=0,nTopImage=1,nSideImage=2};

void MovieQuadrantViewport(int whichImage,int* x, int* y,int* w, int* h)
{
	*w = movieWindowWidth/2-5;
	*h = movieWindowHeight/2-5;

	if (whichImage == nNormalImage)
		{*x = movieWindowWidth/2+5; *y = 0;}
	else if (whichImage == nTopImage)
		{*x = movieWindowWidth/2+5,*y = movieWindowHeight/2+5;}
	else if (whichImage == nSideImage)
		{*x = *y = 0;}
}

// When we are showing stereo pairs (F2) mode, we need to
// split the screen in half, less a little margin in the middle.
// This gives us teh window width.
const  int stereoPairsMarginPixels = 4;
static int StereoPairsViewPortWidth(int movieWindowWidth)
{
	return movieWindowWidth/2-stereoPairsMarginPixels/2;
}

// Whether we are reponding to a PostRedisplay event in our Movie Window,
// or trying to determine which object was selected by the mouse, it is helpful
// to go through this shared function which manages the setup of the viewports
// and transformation matrices, and then calls MovieDisplayCurrentFrame
// above as needed.
// **** There are a number of aesthetic decisions which are made in this
// code regarding depth cueing (fog) etc.  You are welcome and encouraged
// to make changes to these to suit your needs and tastes. 
void MovieRenderOrSelectMainScene
								(GLenum renderMode,
								 GLint mouseX,	// Only used if renderMode == GL_SELECT
								 GLint mouseY, // Only used if renderMode == GL_SELECT
								 GLdouble selectionWidth,  // Only used if renderMode == GL_SELECT
								 GLdouble selectionHeight // Only used if renderMode == GL_SELECT
								 )
{
	glClear(GL_COLOR_BUFFER_BIT |  GL_DEPTH_BUFFER_BIT);

   // Don't bother with FOG if we are in select mode
	if (renderMode != GL_SELECT)
		{
		static GLfloat fogColor[4] = {0.0, 0.0, 0.0, 1.0};
		float z_ClippingStartOfFog;
		float z_MovieImageStartOfFog;
		float z_fog_start;

		glEnable(GL_FOG);
      // Image will fade "linearly" into the background.  See
      // Open GL references.
		glFogi(GL_FOG_MODE, GL_LINEAR);

		glHint(GL_FOG_HINT, GL_FASTEST); /* per pixel */

		// I (personally) want the brightest image possible in the foreground...
      // So, I make the front of the fog area either the molecule
      // size or the clipping pane
      // whichever makes more sense.

		// First compute "near" if we start fogging at clipping outset...
		z_ClippingStartOfFog = (3.0-movie.frontClip)*movieImageSize;

		// But, where is the nearest point that we could actually find part of the movie???
		z_MovieImageStartOfFog = (3.0*movieImageSize - (movieImageSize/2.0) * movie.scaleFactor);

		if (z_MovieImageStartOfFog > z_ClippingStartOfFog)
			z_fog_start = z_MovieImageStartOfFog;
		else
			z_fog_start = z_ClippingStartOfFog;

		// I have personally set this rate of fade out.  It seems to give
      // the right effect on most platforms.
		glFogf(GL_FOG_START, z_fog_start);
		glFogf(GL_FOG_END, z_fog_start + 1.2*movieImageSize*movie.scaleFactor);

		glFogfv(GL_FOG_COLOR, fogColor);
		glClearColor(0.0, 0.0, 0.0, 1.0);
		glDepthFunc(GL_LESS);
		}

	glEnable(GL_DEPTH_TEST);

	assert(glGetError() == GL_NO_ERROR);

   // Whether we are going to draw the movie frame in an F1 (nSingleImage)
   // F2 (nPairOfImages), F3 (TRI) or F4 (hardware stereo) we can compute
   // the clipping box for any of them.
	{
		float final_back_clip = movie.backClip;
		float final_front_clip = movie.frontClip;

		float znear = (3.0-final_front_clip)*movieImageSize; // At program startup, 3-1 = 2
		float zfar = (3.0+final_back_clip)*movieImageSize; // At program startup Starts at 3+2 = 5

		float fov = 22.5/movie.scaleFactor;	// From original MD Display

// Let's prepare the window and call the central routine to
// draw the image
glClearDepth(1.0);

glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

switch ( movie.imageView )
	{
	case nSingleImage:
   	// For this simplest case - F1 single image - the
      // viewport is the entire movieWindow.
		glViewport(0,0,movieWindowWidth,movieWindowHeight);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

      // If we are selecting with the mouse, activate the pick matrix.
		if (renderMode == GL_SELECT)
			InitPickMatrix(mouseX,mouseY,selectionWidth,selectionHeight);

      // Make a perspective transformation (distant objects will be
      // smaller).
		gluPerspective(fov,
							movieWindowAspectRatio,
							znear,
							zfar);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

      // The following calls which apply transformations to
      // the OpenGL scene transformation matrix should be
      // understood as being last-in-first-out.  I.e., OpenGL will
      // translate to center of rotation FIRST, then rotate, etc.

      // Move the origin back into the screen
      // OpenGL does this LAST
		glTranslatef(0,0,(-3.0)*movieImageSize);
      glScalef(movie.scaleFactor,movie.scaleFactor,movie.scaleFactor);

      // Translate the movie left-right and/or up-down
      // based on latest user request
		glTranslatef(movie.xTranslation, movie.yTranslation, 0.0);

      // Rotate based on latest user request
		glMultMatrixf(movie.rotationMatrix);

      // Place center of translation in the middle of things.
      // OpenGL will do this FIRST
		glTranslatef(-movie.centerOfRotation[0],
						 -movie.centerOfRotation[1],
						 -movie.centerOfRotation[2]);

		{
		int err = glGetError();
		if (err != GL_NO_ERROR)
			{
			printf("Terminated point 1 due to OpenGL error %d\n",err);
			exit(EXIT_FAILURE);
			}
		}

		// With the transformation matrix prepared, we are now ready to
      // do the rendering.  See how your code is called from this function
      // above.
		MovieDisplayCurrentFrame(renderMode);

		{
		int err = glGetError();
		if (err != GL_NO_ERROR)
			{
			printf("Terminated point 2 due to OpenGL error %d\n",err);
			exit(EXIT_FAILURE);
			}
		}

		// If the user is rotating with the mouse, the blue disc is a helpful
      // cue.  Draw it now.
		if (movieTrackBallCircleVisible)
			{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

			gluOrtho2D(-(movieWindowWidth*movieWindowFixupFactor)/2,(movieWindowWidth*movieWindowFixupFactor)/2,
				-(movieWindowHeight*movieWindowFixupFactor)/2,(movieWindowHeight*movieWindowFixupFactor)/2);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			glColor3f(0.0,0.0,0.8);
			gluDisk(movieTrackBallCircleQuadric,300.0,300.0,100,1);
			}

		break;


	case nPairOfImages: // User pressed F2 for stereo pairs
		{
		int whichImage;
		enum {nLeftImage=0,nRightImage=1};

      // loop over images for each eye
		for ( whichImage = nLeftImage;
				whichImage <= nRightImage;
				whichImage++ )
			{
         // The width of the view port is half the movieWindow width
         // less the small margin
         int viewPortWidth = StereoPairsViewPortWidth(movieWindowWidth);

         // Place left eye image on left half of window, right image
         // on right half.
			glViewport((whichImage == nLeftImage) ? 0 :
         			movieWindowWidth/2+stereoPairsMarginPixels/2,
						  0, // y starts at 0 for both Renderings...
						  viewPortWidth,
						  movieWindowHeight);

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

         // Prepare the matrix for picking if we are not rendering
			if (renderMode == GL_SELECT)
				InitPickMatrix(mouseX,mouseY,selectionWidth,selectionHeight);

			// Setup a perspective transformation.
			gluPerspective(fov,
								movieWindowAspectRatio,
								znear,
								zfar);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

         // Perform typical translations and rotations - see above for
         // more comments.

         // Lastly, translate movie back into z axis.
			glTranslatef(0,0,(-3.0)*movieImageSize);
			glScalef(movie.scaleFactor * 2.0,movie.scaleFactor,movie.scaleFactor);
			glTranslatef(movie.xTranslation, movie.yTranslation, 0.0);

			// Rotate left image 4 degrees about y axis...
         // This is from old display software (Lybrand, et. al.)
			if ( whichImage == nLeftImage )
				glRotatef(4.0,0.0,1.0,0.0);

			glMultMatrixf(movie.rotationMatrix);

         // First, center coordinate system on center of rotation
			glTranslatef(-movie.centerOfRotation[0],
							 -movie.centerOfRotation[1],
							 -movie.centerOfRotation[2]);

			// With the transformation matrix prepared, we are now ready to
         // do the rendering.  See how your code is called from this function
         // above.
			MovieDisplayCurrentFrame(renderMode);

			// If the user is rotating with the mouse, the blue disc is a helpful
	      // cue.  Draw it now.
			if (movieTrackBallCircleVisible)
				{
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();

				gluOrtho2D(-(viewPortWidth*movieWindowFixupFactor)/2.0,(viewPortWidth*movieWindowFixupFactor)/2,

					-(movieWindowHeight*movieWindowFixupFactor)/2.0,(movieWindowHeight*movieWindowFixupFactor)/2.0);

				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();

				glColor3f(0.0,0.0,0.8);
				gluDisk(movieTrackBallCircleQuadric,300.0,300.0,100,1);
				}


			} // End for loop over left and right images
		} // End Pair of images...
		break;

   // If we are in F3 (TRI) mode, we need to render
   // the image three times - in 3 quadrants of the screen.
	case nThreeImages:
		{
		int whichImage;

      // Loop from 0 to 2
		for ( whichImage = nNormalImage;
				whichImage <= nSideImage;
				whichImage++ )
			{
			int x,y,w,h;
			MovieQuadrantViewport(whichImage,&x,&y,&w,&h);
			glViewport(x,y,w,h);


			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

			if (renderMode == GL_SELECT)
				InitPickMatrix(mouseX,mouseY,selectionWidth,selectionHeight);

			gluPerspective(fov,
							movieWindowAspectRatio,
							znear,
							zfar);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			glTranslatef(0,0,(-3.0)*movieImageSize);
			glScalef(movie.scaleFactor,movie.scaleFactor,movie.scaleFactor);

			glTranslatef(movie.xTranslation, movie.yTranslation, 0.0);

			if (whichImage == nTopImage)
				glRotatef(90.0,1.0,0.0,0.0); // Rotate about X axis to show top
			else if (whichImage == nSideImage)
				glRotatef(90.0,0.0,1.0,0.0); // Rotate about Y axis to show side

			glMultMatrixf(movie.rotationMatrix);

         // First, center coordinate system on center of rotation
			glTranslatef(-movie.centerOfRotation[0],
							 -movie.centerOfRotation[1],
							 -movie.centerOfRotation[2]);

			// With the transformation matrix prepared, we are now ready to
         // do the rendering.  See how your code is called from this function
         // above.
			MovieDisplayCurrentFrame(renderMode);
			} // End loop over 3 images

#ifdef TRACKBALL_CIRCLE_USED_ON_TRI_MODE
		// If the user is rotating with the mouse, the blue disc is a helpful
  	   // cue.  Draw it now - over the middle of the screen.
		if (movieTrackBallCircleVisible)
			{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

			gluOrtho2D(-(movieWindowWidth*movieWindowFixupFactor)/2,(movieWindowWidth*movieWindowFixupFactor)/2,
				-(movieWindowHeight*movieWindowFixupFactor)/2,(movieWindowHeight*movieWindowFixupFactor)/2);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			glColor3f(0.0,0.0,0.8);
			gluDisk(movieTrackBallCircleQuadric,300.0,300.0,100,1);
			}
#endif         
       }
		break;

	// In F4 - hardware stereo mode - we need to superimpose
   // two images in the full movie Window.  In SGILegacy stereo mode,
   // the pixels are rectangluar (double height).  The code here (in
   // conjunction with other code) manages this nuisance fine.
   case nSGILegacyHardwareStereo:
   case nQuadStereo:
		{
		int whichImage;
		enum {nLeftEye=0,nRightEye=1};
		int lastEye = nRightEye;

		glViewport(0,0,movieWindowWidth,movieWindowHeight);

		for ( whichImage = nLeftEye;
				whichImage <= lastEye;
				whichImage++ )
			{
			float eyeFactor;

			if (whichImage == nLeftEye)
				{
				eyeFactor = -0.07;
            // Call GLUT or SGI code to select the left eye buffer
				dispuiSelectEye(0);
				glClearDepth(1.0);
				//glViewport(0,(1024 + 20) /2,1279,(1024-20)/2);
				glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
				}
			else
				{
				eyeFactor = 0.07;
            // Call GLUT or SGI code to select the right eye buffer
				dispuiSelectEye(1);
				//glViewport(0,20,1279,(1024-20)/2);
				// stereo_right_buffer();
				glClearDepth(1.0);
				glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
				}

         // The viewport is the entire window.
			glViewport(0,
						  0, // y starts at 0 for both Renderings...
						  movieWindowWidth,
						  movieWindowHeight);

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

			if (renderMode == GL_SELECT)
				InitPickMatrix(mouseX,mouseY,selectionWidth,selectionHeight);

         // Setup the correct stereo perspective transformations.
         StereoPerspective(fov, // Parameter fovy
							1.29, // Parameter aspect
	                     znear,zfar,
      	            3.0*movieImageSize, // Parameter conv
							eyeFactor*movieImageSize);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

         // See comments above for reading OpenGL transformations
         // bottom up.
         
	      // Lastly, Move the origin back into the screen
			glTranslatef(0,0,(-3.0)*movieImageSize);

         // Scale, translate, and rotate under user control.
         // controls.
			glScalef(movie.scaleFactor,movie.scaleFactor,movie.scaleFactor);
			glTranslatef(movie.xTranslation, movie.yTranslation, 0.0);
			glMultMatrixf(movie.rotationMatrix);

         // First, center coordinate system on center of rotation
			glTranslatef(-movie.centerOfRotation[0],
							 -movie.centerOfRotation[1],
							 -movie.centerOfRotation[2]);

			// Draw the image.
			MovieDisplayCurrentFrame(renderMode);

			// If the user is rotating with the mouse, the blue disc is a helpful
	  	   // cue.  Draw it now - over the middle of the screen.
			if (movieTrackBallCircleVisible)
				{
				float heightDivider = 2.0;
	         if (movie.imageView == nSGILegacyHardwareStereo)
					heightDivider = 1.0;

				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();

				gluOrtho2D(-(movieWindowWidth*movieWindowFixupFactor)/2.0,(movieWindowWidth*movieWindowFixupFactor)/2.0,
					-(movieWindowHeight*movieWindowFixupFactor)/heightDivider,(movieWindowHeight*movieWindowFixupFactor)/heightDivider);

				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();

				glColor3f(0.0,0.0,0.8);
				gluDisk(movieTrackBallCircleQuadric,600.0/heightDivider,600.0/heightDivider,100,1);
				}

			glFlush();
			} // End for loop over left and right eyes
         // For good measure, reselect left eye buffer.  Otherwise, we might
         // have foulups with other subwindow processing on legacy platforms.
			dispuiSelectEye(0);
		} // End Hardware stereo
		break;


#ifdef DEBUG_MODE_OR_BUG_ETC
	// This option lets you see the little rectangles drawn for atoms
	// in pick mode...
	case nTestPickDraw:
		glViewport(0,0,movieWindowWidth,movieWindowHeight);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		if (renderMode == GL_SELECT)
			InitPickMatrix(mouseX,mouseY,selectionWidth,selectionHeight);

		gluPerspective(fov,
							movieWindowAspectRatio,
							znear,
							zfar);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glTranslatef(0,0,(-3.0)*movieImageSize);
      glScalef(movie.scaleFactor,movie.scaleFactor,movie.scaleFactor);

		glTranslatef(movie.xTranslation, movie.yTranslation, 0.0);

		glMultMatrixf(movie.rotationMatrix);

		glTranslatef(-movie.centerOfRotation[0],
						 -movie.centerOfRotation[1],
						 -movie.centerOfRotation[2]);

		{
		int err = glGetError();
		if (err != GL_NO_ERROR)
			{
			printf("Terminated point 1 due to OpenGL error %d\n",err);
			exit(EXIT_FAILURE);
			}
		}

		displayRenderFrameInSelectMode(movieFrame);

		{
		int err = glGetError();
		if (err != GL_NO_ERROR)
			{
			printf("Terminated point 2 due to OpenGL error %d\n",err);
			exit(EXIT_FAILURE);
			}
		}


		if (movieTrackBallCircleVisible)
			{
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

			gluOrtho2D(-(movieWindowWidth*movieWindowFixupFactor)/2,(movieWindowWidth*movieWindowFixupFactor)/2,
				-(movieWindowHeight*movieWindowFixupFactor)/2,(movieWindowHeight*movieWindowFixupFactor)/2);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			glColor3f(0.0,0.0,0.8);
			gluDisk(movieTrackBallCircleQuadric,300.0,300.0,100,1);
			}

		break;
#endif
	} // End switch


   // If OpenGL is in trouble, this is a great point to crash out
   // with an error message.
		{
		int err = glGetError();
		if (err != GL_NO_ERROR)
			{
			printf("Terminated due to OpenGL error %d\n",err);
			exit(EXIT_FAILURE);
			}
		}

		glDisable(GL_DEPTH_TEST);
	}
}

void MovieDisplayFunc(void)
{
	char framestr[20];
	int i;

	glRenderMode(GL_RENDER);	// GL_RENDER or GL_SELECT.....

	// Display the main picture of interest
   // This function above is the big kahuna that sets up the
   // necessary transformations, and calls display* code to
   // render the movie
	MovieRenderOrSelectMainScene(GL_RENDER,0,0,0,0);

   // Finally, we'll give your code the chance to render some graphics
   // in an ortho box.  This could be text or whatever.
for (i=0;i < MovieEyeCount();i++)
	{
	dispuiSelectEye(i);
	// Now let "display" do some final labelling in Ortho mode...
	glViewport(0,0,movieWindowWidth,movieWindowHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluOrtho2D(0,movieWindowWidth,0,movieWindowHeight);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

   // You provide this function.
	displayRenderFrameInOrthoBox(
			movieFrame,
			movieWindowWidth,
			movieWindowHeight);

	// Now update the "Frame Number" in the bottom right corner.
	glColor3f(1.0,1.0,1.0);

	sprintf(framestr,"Frame #%d   ",movieFrame + displayMovieMinFrame());
	DrawString(movieWindowWidth - StringWidth(GLUT_BITMAP_HELVETICA_12,framestr)-4,
					(moviePicoSecondsPerFrame > 0.0) ? 16 : 4,
					GLUT_BITMAP_HELVETICA_12,
					framestr);
#ifndef DISPUI_TESTMODE

	if (moviePicoSecondsPerFrame > 0.0)
		{
		char timestr[20];
		sprintf(timestr,"Time:%3.3f ps",
			moviePicoSecondsPerFrame * movieFrame +
			moviePicoSecondsInitialOffset);
	DrawString(movieWindowWidth - StringWidth(GLUT_BITMAP_HELVETICA_12,timestr)-4,
					3,
					GLUT_BITMAP_HELVETICA_12,
					timestr);
		}
#endif
      
   // If there is a warning/status message to be displayed,
   // and the timer has not expired, put that on the top of the movie
   // window in large white text.
	if (dispuiStatusMessage[0] && (time(NULL) < dispuiStatusMessageExpiration))
		{
		if (movieWindowHeight > 20) // If smaller than this, who cares anyway
			{
			void* textFont = GLUT_BITMAP_HELVETICA_18;
			glColor3f(1.0,1.0,1.0); // White Text
			DrawString(0,movieWindowHeight-20,textFont,dispuiStatusMessage);
			}
		}
	else
		dispuiStatusMessage[0] = 0;
	}
   
	// Let the user see all our work...
	glutSwapBuffers();

	clipWindowFrameToDraw = movieFrame;

#ifndef DISPUI_TESTMODE
   ramaWindowFrameToDraw = movieFrame;

	// If we have a ramaWindow active, only update the rama display if the
	// frame has changed - I.e. don't update it if we are just rotating or
	// scaling or coloring a frame.
	if (ramaWindow != -1)
		{
		glutSetWindow(ramaWindow);
		glutPopWindow();

		if (ramaWindowLastFrameDisplayed != movieFrame)
			{
			glutPostRedisplay();
			ramaWindowLastFrameDisplayed = movieFrame;
			}
		}
#endif      

	if (clipWindow != -1)
		{
		glutSetWindow(clipWindow);
		glutPostRedisplay();
		}
	// If we're showing this in continuous movie mode,
   // (which is the case if the user has not placed the speed bar
   // into step mode)
   // tell glut to call us again, a little later...
	if (/* timerPostRedisplayInitiated &&*/ (! movieStepMode))
		{	
		// timerPostRedisplayInitiated = 0;		
		// if (movieMsecsBetweenFrames)
			glutTimerFunc(movieMsecsBetweenFrames,MovieTimerFunc,movieTimerKey);
		// else
			// MovieTimerFunc(movieTimerKey); // Immediately draw next frame
		}
}

// the left window contains the menu buttons.
// In the case of a window resize, just updatthe viewport
// and move on.
void leftChangeSizeFunc(int w, int h)
{
	if (h == 0)
		h= 1;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);
}


/**********************************************************************************
The bottom window is a long strip which has three user interface components:
	1) A graphical sliding bar which shows the "scale" of display, defaulting
		to 1.0.  This is managed by the bottomScaleXXX variables.

   2) A graphical sliding bar which shows the speed of the animation.  This bar
			has green (fast) and yellow (slow) areas in forwards and reverse directions.

   3) A "Hyper" Box which tells the number of frames (0 if OFF) to skip with
       each timestep.
*************************************************************************************/


int bottomTopTextY;			// Line for "top" text above graphics, -1 if text not shown
int bottomGraphicsY;			// Line at which all three "graphics" appear.
int bottomGraphicsHeight; 	// Height for rectangular graphics boxes
int bottomBottomTextY;		// Line for "bottom" text under graphics, -1 if text not shown

int bottomScaleGraphicWidth;	// Width, in pixels, of Graphics scale box.
const char* bottomScaleString = "<[>  --Scale--  <]>";
int bottomScaleStringX;			// X coordinate of where "--Scale--" is drawn
int bottomScaleLeftX;			// X coordinate of left end of graphics display.

// Breaking this into fragments lets us match graphics pieces to it more easily.
const char* bottomFastSlowString = "Fast-----Slow---";
const char* bottomStepString = "--Step--";
const char* bottomSlowFastString = "---Slow-----Fast";

const char* bottomReverseForwardString = "<-- Reverse  Forward -->"; // Below Graphic


// These variables contain the widths of the pieces of the strings.
// If the fonts are "turned off" due to small window size then these widths are 
// just pixel counts and the strings do not matter.
int bottomFastSlowStringWidth;
int bottomStepStringWidth;
int bottomSlowFastStringWidth;
int bottomReverseForwardStringWidth;

int bottomFastSlowY = 40;
int bottomFastSlowHeight = 20;

int bottomFastReverseX;  // Left Edge of Fast Reverse Area
int bottomSlowReverseX;  // Left Edge of Slow Reverse Aaea
int bottomStepReverseX;	 // Left Edge of Step (Reverse) Area
int bottomStepForwardX;	 // Left Edge of Step (Forward) Area
int bottomSlowForwardX;  // Left Edge of Slow Forward Area
int bottomFastForwardX;  // Left Edge of Fast Forward Area
int bottomFastForwardX2; // Right Edge of Fast Forward Area

int bottomSpeedStringHeight;
void* bottomFont = 0;
int bottomFontHeight = 0;
int bottomSpeedStringX;


int bottomReverseForwardStringX;

const char* bottomHyperString = "Hyper";
const char* bottomHyperOFFString = "OFF";
int bottomHyperStringWidth;
int bottomHyperStringX;	// Also X for "Hyper Box" Graphic
int bottomHyperOFFStringX;
int bottomHyperOFFStringWidth;


int hyperModeFlag = 0; // If true - then we go "full throttle" using the IdleFunc


int BottomMovieScaleBarXpos(void)
{
	return PercentageBarXpos(log(movie.scaleFactor*10),log(movieMinScale*10),log(movieMaxScale*10),
										bottomScaleLeftX,bottomScaleGraphicWidth);

}

// Check font from most to least preferable, and return best preferance

/*****************************************************
 This code sorts out placement of the scaling and speed bars, as well
 as the Hyper box.
 
 As the user resizes the main window, the bottom window may need to 
 grow or shrink.  Since GLUT's bitmapped fonts are being
 used, we have to be smart about shrinking fonts aesthetically.
 
 If font is 0, this routine assumes we are making a 
 last-ditch effort (the window is small) to load the GUI parts
 without their text descriptive strings.
******************************************************/

int BottomPlaceGraphics(int w,int h,void* font,int fontHeight)
{
	int success = 0;

	bottomFont = 0;
	bottomFontHeight = 0;

   // Job one is to determine some element positions within the
   // bottom window.
	// Can we fit a nice scaling bar, Fast Strings, etc in there?
	if ( font != 0 )  // Then, we want to be able to display the text!
		{
		if (movie.imageView == nSGILegacyHardwareStereo)
			bottomGraphicsHeight = (int)((double)fontHeight * 2.0);
		else
			bottomGraphicsHeight = (int)((double)fontHeight * 2.8);

		// Allow height for everything - and some border area.
		bottomWindowHeight = bottomGraphicsHeight + 2*fontHeight + 12;
		bottomTopTextY = bottomWindowHeight - fontHeight;
		bottomGraphicsY = bottomTopTextY-3-bottomGraphicsHeight-3;
		bottomBottomTextY = bottomGraphicsY-fontHeight-2;

		if (movie.imageView == nSGILegacyHardwareStereo)
			{
			bottomTopTextY-=2;
			bottomBottomTextY--;
			}
		}
	else // If no text - this is easier
		{
		// 20 is smallish - but should be enough, by my eye...
		bottomGraphicsHeight = 20;
		bottomWindowHeight = bottomGraphicsHeight + 6;
		bottomTopTextY = -1;
		bottomGraphicsY = bottomWindowHeight-3-bottomGraphicsHeight;
		bottomBottomTextY = -1;
		}

	// *4 is aesthetic.  If room for height is here, then press on...
	// else we'll return success=0
	if ( h > bottomWindowHeight * 4)
		{
		// Now that "height" is OK - let's see about width.
		// If we have a font, and the long "Reverse/Forward etc"
		// text can fit, everything else is a "go".
		bottomHyperStringWidth = StringWidth((font == 0) ? GLUT_BITMAP_HELVETICA_10 : font,bottomHyperString);

		bottomHyperStringX = w -
									bottomHyperStringWidth -
									bottomHyperStringWidth / 2;
                           
      // In the Forward and Reverse sides of the speed bar are different colored areas for 
      // going fast, slow or stepping.  The widths of these areas must be computed to allow
      // for decoding of mouse events in these areas.       

		if ( font != 0 )
			{
			// Break up pieces of the GUI slider for speed...
			int fastSlowWidth = StringWidth(font,bottomFastSlowString);
			int stepWidth = StringWidth(font,bottomStepString);
			int slowFastWidth = StringWidth(font,bottomSlowFastString);
			int fastOnlyWidth = StringWidth(font,"Fast---");

			int bottomSpeedStringWidth =
					fastSlowWidth + stepWidth + slowFastWidth;

			bottomFont = font;
			bottomFontHeight = fontHeight;

			// Make sure that the string can fit in the space of bigger than 1/4 or so.
			if (bottomSpeedStringWidth <= ((double)(w)/2.5)) // We have a good one!
				{
				bottomFastReverseX = bottomHyperStringX - bottomSpeedStringWidth - bottomHyperStringWidth/2;
				bottomSlowReverseX = bottomFastReverseX + fastOnlyWidth;
				bottomStepReverseX = bottomFastReverseX + fastSlowWidth;
				bottomStepForwardX = bottomStepReverseX + stepWidth/2;

				bottomSlowForwardX = bottomStepReverseX + stepWidth;
				bottomFastForwardX = bottomSlowForwardX + fastSlowWidth - fastOnlyWidth;
				bottomFastForwardX2 = bottomFastForwardX + fastOnlyWidth-1;

				// Compute Center of Speed String, and then center ReverseForward about that.
				bottomReverseForwardStringX = bottomFastReverseX + bottomSpeedStringWidth/2
															- StringWidth(bottomFont,bottomReverseForwardString)/2;

				success = 1;



				bottomHyperOFFStringWidth = StringWidth(font,bottomHyperOFFString);
				bottomHyperOFFStringX = bottomHyperStringX +
													bottomHyperStringWidth/2 -
													bottomHyperOFFStringWidth/2;
				}
			}
		else // We'll just try the graphics - font is turned off..
			{
			int bottomSpeedStringWidth = w/2.5; // Arbitrary aesthetic minimum!
			// Divide into 5 panes of equal length - since no strings to worry with.
			int paneWidth = bottomSpeedStringWidth / 5;

			// "40" gives us headroom for other elements along the bar
			if (bottomSpeedStringWidth >= 40) // We have a good one!
				{
				bottomFastReverseX = bottomHyperStringX - bottomSpeedStringWidth - bottomHyperStringWidth/2;
				bottomSlowReverseX = bottomFastReverseX + paneWidth;

				bottomStepReverseX = bottomSlowReverseX + paneWidth;
				bottomStepForwardX = bottomStepReverseX + paneWidth/2;

				bottomSlowForwardX = bottomStepReverseX + paneWidth;
				bottomFastForwardX = bottomSlowForwardX + paneWidth;
				bottomFastForwardX2 = bottomFastForwardX + paneWidth-1;

				bottomReverseForwardStringX = 0;

				success = 1;
				}
			}
		}

	if ( success ) // Looking good - so stick the scaling graphic to the left side...
		{
		bottomScaleGraphicWidth = 4*w/9; // Was bottomHyperStringX/2;
		bottomScaleLeftX = 10; // Push all the way over left in the bottom Box to the left a good bit;

		bottomScaleStringX = 0;

		// Cetner the "--Scale--" string over the graphic
		if (bottomFont != 0)
			bottomScaleStringX = bottomScaleLeftX +
											(bottomScaleGraphicWidth -
												StringWidth(bottomFont,bottomScaleString))/2;
		}
	else // Zero some variables - though we should not be drawing (referring to these) anyway.
		{
		bottomScaleGraphicWidth = 	bottomScaleLeftX = bottomScaleStringX = 0;

		bottomWindowHeight = 0;
		bottomFont = 0;
		}

	return success;
}

// called by GLUT when the bottom subWindow chanes size.
void BottomChangeSizeFunc(int w, int h)
{
	if (h == 0)
		h= 1;

	glViewport(0,0,w,h);

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluOrtho2D(0,w,0,h);
}

// In the left subWindow we display Menu Boxes.  This function draws the menu box for a specific 
// menu item - which may have submenues delineated by | characters.  The parameters are simply
// the menuItem number, text to display, and upper left/bottom right menu box positions.

void leftDrawAndInitMenuBox(int menuItem,const char* menuBoxText,float x1,float y1,float x2, float y2)
{
	int menuButtonTextWidth;

	int subMenuCount = 0;	// Submenus are the options separated by | marks
	const char* subMenuText[nMaxSubMenus];
	int subMenu;

	char* tempstr = strdup(menuBoxText);
	char* s = tempstr;
	int menuFontYpos,menuBarYpos;

   // The | character delineates active areas on a single menu button
	menuAsciiBarWidth = StringWidth(menuFont,"|");

	assert(tempstr != 0);
   
	// Center text in button
	menuButtonTextWidth = StringWidth(menuFont,menuBoxText);

   // Locate the ASCII | characters in the string, and zero them out.
	while (*s)
		{
		char* nextBar = strchr(s,'|');
		subMenuText[subMenuCount++] = s;

		if (nextBar != 0)
			{
			s = nextBar;
			*s++ = 0;
			}
		else
			{
			s+= strlen(s); // Will point us at final null
			}
		}

	subMenus[menuItem] = subMenuCount;


	glColor3f(1.0,1.0,1.0); // White border "around" menu box
	glRectf(x1,y1,x2,y2);

	glColor3dv(menuBoxColor); // Brownish box inset in white border
	glRectf(x1+1.0,
			  y1-1.0,
			  x2-1.0,
			  y2+1.0);

	glColor3f(1.0,1.0,1.0); // White for Menu Button text

	subMenuXTextpos[menuItem][0] = x1 + (x2-x1)/2 - (menuButtonTextWidth/2);

	menuFontYpos = y1 - menuFontHeight - (menuBoxHeight-menuFontHeight)/2.0 + 2;

	menuBarYpos = y2;

	DrawString( subMenuXTextpos[menuItem][0], // X
					menuFontYpos, // Y
					menuFont,
					subMenuText[0]);

   // For each subMenu output a full vertical line and draw the next subMenu text.
	for (subMenu=1;subMenu<subMenuCount;subMenu++)
		{
		// Draw a full height bar using glRectf.  I.e. - we don't output a | character
      // It's just not as attractive.
		int lastSubMenuEndPos = subMenuXTextpos[menuItem][subMenu-1] + 
                   StringWidth(menuFont,subMenuText[subMenu-1]);

      glRectf(lastSubMenuEndPos+menuAsciiBarWidth/2,
               menuBarYpos,
               1+lastSubMenuEndPos+menuAsciiBarWidth/2,
            menuBarYpos+menuBoxHeight);

      // Advace current output position
		subMenuXTextpos[menuItem][subMenu] = lastSubMenuEndPos + menuAsciiBarWidth;

      // Write text for this subMenu option.
		DrawString(  subMenuXTextpos[menuItem][subMenu], // X
						menuFontYpos, // Y
					menuFont,
					subMenuText[subMenu]);
		}

	free(tempstr);
}

// The Clipping tool lets the user section through the 3D image.  When activated from the
// menu the clip Window shows a tiny version of the movie, birdseye fashion, in the bottom
// left corner of the screen.

int clipWindowHeight,clipWindowWidth;
int clipMouseY0;              // position where mouse button first went down.
int clipMouseDownFlag = 0;    // TRUE if mouse button is currently down in the clipping window
int clipLockFlag = 0;		   // TRUE if User clicks "Clip Locked"
float clipLockDistance = 0;   // Distance between locked clipping panes in tool.


void ClipScaleMouse(int* x,int* y)
{
	// The clipping window tool is always a 60 x 250 _scaled_ window.
	// But, GLUT mouse clicks are always in physical window coordinates.

	// Let's fix this up.

	// First, y needs to range from 0 _up_ to height - and not down
	(*y)=clipWindowHeight - (*y);

	// And, we need to "scale" the mouse coordinates to the 60 x 250 window
	// that we work with elsewhere.  So, if the window is physically smaller than 250,
	// we scale y UP - and vice versa.

	(*y) =  (float)(*y) * 250.0/(float)clipWindowHeight;

	// Same idea for x
	(*x) =  (float)(*x) * 60.0/(float)clipWindowWidth;
}

enum {nClipFrontControl,nClipBackControl,nClipUnknown};  // info about where was button pressed.
int clipMouseDownLocation; // a value which will record the above state at mouse down time

int ClipMouseLocation (int scaledMouseX,int scaledMouseY)
{
	int retval = nClipUnknown;

	if (scaledMouseX < 29) // We are in the "Front Clip Control Area"
		retval = nClipFrontControl;

	if (scaledMouseX > 31) // We are in the "Back Clip Control Area"
		retval = nClipBackControl;

	return retval;
}


// A minimal Change SizeFunction which will be called by GLUT when needed.
void ClipChangeSizeFunc(int w, int h)
{
	clipWindowHeight = h;
	clipWindowWidth = w;

	glViewport(0,0,clipWindowWidth,clipWindowHeight);

	glClearColor(0.0,0.0,0.0,1.0);
}


void ClipDisplayFunc(void)
{
	int frpos,bkpos;
	int eye;

	// printf("Here with clip frame %d\n",clipWindowFrameToDraw);

// If we are on a legacy SGI platform, take time to draw in both eyes if in stereo
for (eye = 0;eye < SGIEyeCount();eye++)
	{
	dispuiSelectEye(eye);

	glClear(GL_COLOR_BUFFER_BIT);

	// Rese4t coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// This is a "hard coded" 60 x 250 drawing area.
	// In other Clip* functions, the mouse clicks, etc. are scaled into this space.
	glOrtho(0,60.0,
			  0,250.0,
			  0,movieImageSize);


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(30.0f,100.0f,-movieImageSize/2);



   // Rotate 90 degrees around X axis - so Y becoms Z and we have a bird's
   // eye view down onto the movie frame in the clipping window
	glRotatef(90.0,1.0,0.0,0.0);	

	glMultMatrixf(movie.rotationMatrix);

	glScalef(50.0/(movieImageSize/* *movieWindowAspectRatio */),
				50.0/movieImageSize,
				50.0/movieImageSize);
   glScalef(movie.scaleFactor,movie.scaleFactor,movie.scaleFactor);

   // Your code is welcome to draw the entire movie frame.  Or, it might elect to draw
   // a simpler wireframe image to enjoy a speedup.
	displayRenderClipToolFrame(clipWindowFrameToDraw);

   // We're now ready to draw the clipping tool bars (thin rectangles) on your image.
   glLoadIdentity();
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
	glOrtho(0,60,
			  0,250,
			  0,1);

   // Compute positions of front and back bars
   frpos=(100.0-movie.frontClip*50.0);
   bkpos=(100.0+movie.backClip*50.0);

   glColor3f(.9,.9,.9);		// Nice strong grey

	// First, a vertical line up through center
   glRectf(29,0,31,250);

   if (clipLockFlag)
		{
	   glColor3f(1.0,0.0,0.0);
	   glRectf(29,frpos-1,31,bkpos+1); // Red Vertical "lock bar"
		}
   else
	  glColor3f(1.0,1.0,0.0); // Yellow

   // The two clip bars.. first front - then back
   glRectf(0,frpos-1,35,frpos+1);
   glRectf(25,bkpos-1,60,bkpos+1);
	}

  glutSwapBuffers();
}

// Whether the user is dragging the mouse or clicking in a new location, the 
// software needs to respond immediately.
void ClipHandleMouse(int mouseLocation,int glutMouseState,int y)
{
	int clip_refresh_needed = 0; // If the mouse is pressed in an odd area do nothing.
   
	switch (mouseLocation) {
		case nClipFrontControl: // User is adjusting the front clipping pane
			movie.frontClip = (100-y)/50.0;
			//printf("clip LockDistance = %f,movie.frontClip1 = %f",clipLockDistance,movie.frontClip);
			if (clipLockFlag) // If locked, move back pane simultaneously.
				{
				if (movie.frontClip < clipLockDistance-3.0)
					movie.frontClip  = clipLockDistance-3.0;
				movie.backClip=clipLockDistance-movie.frontClip;
				//printf("movie.frontClip2 = %f\n",movie.frontClip);
				}
			else
				{
				if (-movie.frontClip > movie.backClip)
					movie.frontClip = -movie.backClip;
				}
			clip_refresh_needed = 1;
		break;

		case nClipBackControl: // User is adjusting the back clipping pane
			movie.backClip = (float)((y-100)/50.0);
			// printf("clip LockDistance = %f,movie.backClip1 = %f",clipLockDistance,movie.backClip);
			if (clipLockFlag) // If lock is on, move front pane simultaneously
				{
				if (movie.backClip < clipLockDistance-2.0)
					movie.backClip = clipLockDistance-2.0;
				movie.frontClip=clipLockDistance-movie.backClip;
				//printf("movie.backClip2 = %f\n",movie.backClip);
				}
			else
				{
				if (-movie.frontClip > movie.backClip)
					movie.backClip = -movie.frontClip;
				}
			clip_refresh_needed = 1;
		break;

		case nClipUnknown:
		default:  
         // Not sure if this is ever taken - mouse coverage in this window is broad!
         // But, we'll not generated a refresh by leaving clip_refresh_needed = 0;
			break;

		} // End Switch

	if (clip_refresh_needed)
		{
		int curwin = glutGetWindow();
		glutSetWindow(clipWindow);
		glutPostRedisplay();
		glutSetWindow(curwin);
      // Show Movie Frame with new clip area..
      // The new redisplay will trigger a redisplay of the clipping window
		MovieAdvanceTimerKeyAndPostRedisplay(); 
		}
}

void ClipMouseFunc(int button,int state,int x,int y)
{
	// Convert "mouse" GLUT x,y coordinates to OpenGL coordinate!
	ClipScaleMouse(&x,&y);

	if (button == GLUT_LEFT_BUTTON)
		{
		if (state == GLUT_DOWN)  // User has pressed mouse
			{
			clipMouseY0 = y;

         // Determine where the press was
			clipMouseDownLocation = ClipMouseLocation(x,y);

         // Move the clip bar to the mouse
			ClipHandleMouse(clipMouseDownLocation,state,y);
         
         // So now the mouse is being dragged.
			clipMouseDownFlag = 1;
			}
		else
			{
         // Release the mouse
			clipMouseDownFlag = 0;

			if (state == GLUT_UP) // On release of mouse we want to..
				{
            // move the clip bar to release point - but ONLY if release
            // point is in same general area of the mouse down location.
				if (clipMouseDownLocation == ClipMouseLocation(x,y))
					ClipHandleMouse(clipMouseDownLocation,state,y);
				else // If we release in a "funny" place, go back to original
					{
					ClipHandleMouse(clipMouseDownLocation,state,clipMouseY0);
					}

				}
			}

		}
}

// This function is called by GLUT when the user is moving the mouse
// in the clip window.
void ClipMotionFunc(int x,int y)
{
	// Convert "mouse" x,y coordinates to OpenGL coordinate!
	ClipScaleMouse(&x,&y);

   // If we are dragging the mouse in the same clip tool area (front or back)
   // as where it went down, then move the clip tool bar right along.
	if (clipMouseDownFlag && (ClipMouseLocation(x,y) == clipMouseDownLocation))
		ClipHandleMouse(clipMouseDownLocation,-1000,y);
}

// this function is called when the user clicks on the clipping tool menu botton, or 
// when the master window is resized, and other code needs to resize the clipping tool window.
// This code requires that the leftWindow, just above the clipping window, have been created first
void ClipWindowCreateOrChange(void)
{
	int save_current_window = glutGetWindow();
	int h;
	int clip_window_height;
	int clip_window_width;

	// For compatability with old Display, we want this thing to be 60x250 -
	// and we are very happy to scale to get there!
	glutSetWindow(masterWindow);
	h = glutGet(GLUT_WINDOW_HEIGHT);

#ifdef SGI_STEREO
	if (movie.imageView == nSGILegacyHardwareStereo)
		h = (h - VPSEP) / 2;
#endif

   // Use all remaining space below the leftWindow
	clip_window_height = h-leftWindowHeight;
	clip_window_width = clip_window_height * (60.0/250.0);

   // Scale the clipping Widow so that it makes sense.
	if (clip_window_width > leftWindowWidth)
		{
		clip_window_width = leftWindowWidth;
		clip_window_height = clip_window_width * 250.0/60.0;
		}

	if (clipWindow == -1) // Then we need to create it
		{
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		clipWindow = glutCreateSubWindow(masterWindow,
													0,							// x
													leftWindowHeight,		// y
													clip_window_width,		// width
													clip_window_height);	// height

		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);

		glutDisplayFunc(ClipDisplayFunc);
		glutReshapeFunc(ClipChangeSizeFunc);
		glutMotionFunc(ClipMotionFunc);
		glutMouseFunc(ClipMouseFunc);
		}
	else // We just need to resize it and reposition it
		{
		glutSetWindow(clipWindow);
		glutPositionWindow(0,leftWindowHeight);
		glutReshapeWindow(clip_window_width,clip_window_height);
		}


	glutSetWindow(save_current_window);
}

// This is called by GLUT in response to creation of the help Window when the
// user presses '?'
void HelpDisplayFunc(void)
{
	int eye;

	for (eye = 0;eye < SGIEyeCount();eye++)
		{
		dispuiSelectEye(eye);

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		displayRenderHelpBox();
		}

	glutSwapBuffers();
}

// GLUT calls this function when the help window is first created and when it is subsequently resized.
void HelpChangeSizeFunc(int w, int h)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);
	glClearColor(1.0,1.0,1.0,1.0);
}


// Make the helpwindow visible/invisible after user presses '?' or clicks on the 
// help window.  Depending on your particular application, you may want to
// move the help window default position, or change the size.

void HelpWindowToggle(void)
{
	if (helpWindow == -1)
		{
		int currentWindow = glutGetWindow();

		// Let's make start window about 90% of screen real estate.
		int screen_w = glutGet(GLUT_SCREEN_WIDTH);
		int screen_h = glutGet(GLUT_SCREEN_HEIGHT);

		if (movie.imageView == nSGILegacyHardwareStereo)
			screen_h /= 2;
		glutInitWindowSize(600,400);

		glutInitWindowPosition(screen_w-600-20,screen_h-400-50);

		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		helpWindow = glutCreateWindow("Help:   Press <?> again to close.");
		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);

		glutPopWindow();
		glutDisplayFunc(HelpDisplayFunc);
		glutReshapeFunc(HelpChangeSizeFunc);
		glutSetWindow(currentWindow);
		}
	else
		{
		glutDestroyWindow(helpWindow);
		helpWindow = -1;
		}
}

void ClipWindowToggle(void)
{
	if (clipWindow == -1)
		{
      // only show a clip window if left sub window is displayed...
      // Otherwise, we're likely in "hide" mode or the master window has
      // been shrunken dramatically.
		if (leftWindow != -1) 
			{
			ClipWindowCreateOrChange();
			clipWindowFrameToDraw = movieFrame;
			}
		}
	else
		{
		glutDestroyWindow(clipWindow);
		clipWindow = -1;
		}
}

// Toggle the appearance of the color palette
void PaletteWindowToggle(void)
{
	if (paletteWindow != -1)
		{
		glutDestroyWindow(paletteWindow);
		paletteWindow = -1;
		}
	else
		{
		PaletteWindowCreateOrChange();
		}
}

#ifndef DISPUI_TESTMODE

void RamaWindowToggle(void)
{
	if (ramaWindow == -1)
		{
		int currentWindow = glutGetWindow();
		if (movie.imageView == nSGILegacyHardwareStereo)
			glutInitWindowSize(300,150);
		else
			glutInitWindowSize(300,300);
		glutInitWindowPosition(leftWindowWidth,60);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		ramaWindow = glutCreateWindow("Ramachandron Plot   Press <R> to close.");
		glutMouseFunc(RamaMouseFunc);
		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);

		glutPopWindow();
		glutDisplayFunc(RamaDisplayFunc);
		glutReshapeFunc(RamaChangeSizeFunc);
		glutSetWindow(currentWindow);
		ramaWindowFrameToDraw = movieFrame;
		}
	else
		{
		glutDestroyWindow(ramaWindow);
		ramaWindow = -1;
		}

	ramaWindowLastFrameDisplayed = displayMovieMinFrame()-1000; // Presumably, a safe way to say "not displayed"
}

#endif

// Given a menu button with several subMenus, compute the "submenu"
int LeftSubMenuCompute(int menuItem,int mouseX)
{
	int last = subMenus[menuItem]-1;
	int i;

	if (last == 0)	// If we only have one submenu... return it..
		return 0;

	// Check out first submenu - from left to just left of first bar
	if (mouseX < subMenuXTextpos[menuItem][1]-menuAsciiBarWidth/2-2)
		return 0;

	// Check out last submenu - from right of last bar to end of menu box...
	if (mouseX > subMenuXTextpos[menuItem][last]-menuAsciiBarWidth/2)
		return last;


	for (i=1;i<last;i++)
		{
		if (mouseX > subMenuXTextpos[menuItem][i]-menuAsciiBarWidth/2 &&
			 mouseX < subMenuXTextpos[menuItem][i+1]-menuAsciiBarWidth/2-2)
		return i;
		}

	return -1;	// Must have clicked right on a bar - so do nothing...
}

// When the mouse if clicked in the left window (where the menu boxes are)
// we need to figure out which menu button is in play.
void LeftMouseFunc(int button,int state,int x,int y)
{
	static int lastDownCandidate=-1;

	// Note: Here y = 0 means mouse click was at top of window

	// printf("Mouse func (%d,%d)\n",x,y);
	//Are we in the ballpart of a button?
	if (y > menuTop && y < menuTop + menuBoxSpace*nMenuItemCount &&
		 x >= menuBoxLeft && x < menuBoxLeft + menuBoxWidth)
		{
		int candidateItem = (y-menuTop) / menuBoxSpace;
		int candidateLine = (y-menuTop) % menuBoxSpace;

		if ((candidateLine >= 0) && (candidateLine < menuBoxHeight))
		{
		if (state == GLUT_DOWN)
			{
			lastDownCandidate = candidateItem;
			}
		else
		if (state == GLUT_UP) // Did we release where we pressed?
			{
			if (candidateItem == lastDownCandidate)
				{
				// printf("MouseFunc on candidateItem=%d\n",candidateItem);
				switch (candidateItem)
				{
				case nHelp:
					HelpWindowToggle();
					break;

				case nColors:
					PaletteWindowToggle();
					break;
				case nHide:
					MovieHideLeftAndBottom(1);
					break;

				case nOrigin:
					displayCenterOfRotationRequested();
					break;

#ifndef DISPUI_TESTMODE               
				case nRemove:
					displayRemoveRequested();
					break;

				case nHBond:
					displayHBondRequested();
					break;

				case nDistanceOrAngleOrTorsion:
					{
					int subMenu = LeftSubMenuCompute(nDistanceOrAngleOrTorsion,x);
					if (subMenu == 0)
						displayDistanceMeasurementRequested();
					else if (subMenu == 1)
						displayAngleMeasurementRequested();
					else if (subMenu == 2)
						displayDihedralMeasurementRequested();
					}
					break;

				case nRama:
					RamaWindowToggle();
					break;

				case nUnLabelOrLabel:
					{
					int subMenu = LeftSubMenuCompute(nUnLabelOrLabel,x);
					if (subMenu == 0)
						displayUNLabelRequested();
					else if (subMenu == 1)
						displayLabelRequested();
					}
					break;
#endif

				case nFullScreen:
					if (movie.imageView != nSGILegacyHardwareStereo)
						MasterSetFullScreen(! dispuiFullScreenFlag);
					break;

				case nDisplay:
					MovieAdvanceImageView();
					break;

				case nClipping:
					ClipWindowToggle();
					break;


				case nClipLock:
					clipLockDistance = movie.frontClip+movie.backClip;
					if (clipWindow == -1)
						{
						clipLockFlag = 1;

						ClipWindowToggle();
						}
					else
						{
						clipLockFlag = ! clipLockFlag;
						glutSetWindow(clipWindow);
						glutPostRedisplay();
						}
					break;

#ifndef DISPUI_TESTMODE

				case nFiltering:
					displayInputNewFilter();
					break;

				case nPDB:
					displayPDBRequested();
					break;

				case nQuit:
					displayQuitRequested();
					break;
#endif


				default:
					// KeyboardInputStart("Enter \"wire\" or \"solid\" :",InputCompleteCallback);
					break;
				}
				}
			lastDownCandidate = -1;
			}
		else
			{
			printf("Bad mouse state %d",state);
			exit(EXIT_FAILURE);
			}
		}
		}
	else
		lastDownCandidate = -1;
}


// GLUT calls this function to display the left sub window - which
// contains the menu buttons and a coordinate system axes graphic
void leftDisplayFunc(void)
{
	int menu;

   int eye;
	int h = glutGet(GLUT_WINDOW_HEIGHT);
	int w = glutGet(GLUT_WINDOW_WIDTH);

	static float xAxisColor[] = {1.0,1.0,1.0};
	static float yAxisColor[] = {1.0,1.0,0.4};
	static float zAxisColor[] = {0.4,1.0,1.0};

	int axis_length = w/3;

	if (axis_length > menuTop/3)
		axis_length = menuTop/3;

for (eye = 0;eye < SGIEyeCount();eye++)
	{
	dispuiSelectEye(eye);

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	glOrtho(0,w,0,h,-30.0,30.0);

   // Black out the entire window first
	glClearColor(0.0,0.0,0.0,0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// First, draw the axes at the top of the menu
	glTranslatef(w/2,h-menuTop/2,0.0);

	glMultMatrixf(movie.rotationMatrix);

   // Draw the three axis lines
	glBegin(GL_LINES);
		glColor3fv(xAxisColor);
		glVertex3i(0,0,0);
		glVertex3i(axis_length,0,0);

		glColor3fv(yAxisColor);
		glVertex3i(0,0,0);
		glVertex3i(0,axis_length,0);

		glColor3fv(zAxisColor);
		glVertex3i(0,0,0);
		glVertex3i(0,0,axis_length);
	glEnd();

	// Draw arrows at end of axes
	glBegin(GL_LINE_STRIP);
		glColor3fv(xAxisColor);
		glVertex3s(15,5,0);
		glVertex3s(20,0,0);
		glVertex3s(15,0,5);
	glEnd();

	glBegin(GL_LINE_STRIP);
		glColor3fv(yAxisColor);
		glVertex3s(5.0,15.0,0.0);
		glVertex3s(0.0,20.0,0.0);
		glVertex3s(0.0,15.0,5.0);
	glEnd();

	glBegin(GL_LINE_STRIP);
		glColor3fv(zAxisColor);
		glVertex3s(5.0,0.0,15.0);
		glVertex3s(0.0,0.0,20.0);
		glVertex3s(0.0,5.0,15.0);
	glEnd();


   // Now label the three axes
	glColor3fv(xAxisColor);
	glRasterPos3i(axis_length+3,0,0);
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'x');

	glColor3fv(yAxisColor);
	glRasterPos3i(0,axis_length+3,0);
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'y');

	glColor3fv(zAxisColor);
	glRasterPos3i(0,0,axis_length+3);
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,'z');


	// End Axes Draw
	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);

	// And draw some menu buttons!

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	for (menu=0;menu<nMenuItemCount;menu++)
		{
		float x1 = menuBoxLeft;
		float y1=  leftWindowHeight - menu*menuBoxSpace - menuTop;
		float x2=  x1 + menuBoxWidth;
		float y2 = y1 - menuBoxHeight;

		leftDrawAndInitMenuBox(menu,menuStrings[menu],x1,y1,x2,y2);
		}

//	glFlush();
	}
	glutSwapBuffers();
//	glutTimerFunc(2000,TimerFunc,2);
}

// The static button allows you to superimpose your movie on a first, or average
// view.  The menu button is actually another subwindow to which GLUT menu items are associated.
// So, when the static "window" is displayed, we actually display on the leftWindow
void staticDisplayFunc(void)
{
	int eye;

for (eye = 0;eye < MovieEyeCount();eye++)
	{
	dispuiSelectEye(eye);
	glClear(GL_COLOR_BUFFER_BIT);

	leftDrawAndInitMenuBox(nStatic,menuStrings[nStatic],0,menuBoxHeight,menuBoxWidth,0);
	glFlush();
	}
	glutSwapBuffers();
}

// GLUT calls this function when the static subWindow changes size.
void staticChangeSizeFunc(int w, int h)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,w,h);

	gluOrtho2D(0,w,0,h);
	// Clear to Brown - same color as menu box we are overwriting
	glClearColor(
			menuBoxColor[0],
			menuBoxColor[1],
			menuBoxColor[2],1.0);
}

// GLUT calls this function to display the subWindow at the bottom - which
// includes scale bar, speed bar, and HyperControl.
void BottomDisplayFunc(void)
{
	int eye;

// Draw in both eyes if we are in SGI legacy stereo mode.
for (eye = 0;eye < SGIEyeCount();eye++)
	{
	dispuiSelectEye(eye);

	glClearColor(0.0,0.0,0.0,0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	if (bottomFont != 0) // Then we were able to fit some text in here.
		{
		char tempstr[10];
		char *decPoint;
		int width_before_decimal_point;

      // First draw the scale bar - and the decimal representation of scaling
		sprintf(tempstr,"%02.2f",movie.scaleFactor);
		decPoint = strchr(tempstr,'.');
		if (decPoint != 0) // Should always be found
				*decPoint = 0; // Nuke it temporarily

		width_before_decimal_point = StringWidth(bottomFont,tempstr);

		// Put back the '.' for display
		if (decPoint != 0)
				*decPoint = '.';

		glColor3f(1.0,1.0,1.0);

		// Draw "scale" above the graphic
		DrawString(bottomScaleStringX,bottomTopTextY,
					  		bottomFont,bottomScaleString);

		// Draw the value of scale aligned with the black bar
		DrawString(BottomMovieScaleBarXpos()-width_before_decimal_point-StringWidth(bottomFont,".")/2,
							bottomBottomTextY,bottomFont,tempstr);

		// Draw the top "Fast---Slow---" string
		DrawString(bottomFastReverseX,bottomBottomTextY,
					  bottomFont,bottomFastSlowString);

		// Draw the top "--Step--" string
		DrawString(bottomStepReverseX,bottomBottomTextY,
					  bottomFont,bottomStepString);

		// Draw the top "---Slow---Fast" string
		DrawString(bottomSlowForwardX,bottomBottomTextY,
					  bottomFont,bottomSlowFastString);

		// Under the graphic, draw the <-- Reverse -- Forward --> string
		DrawString(bottomReverseForwardStringX,bottomTopTextY,
					  bottomFont,bottomReverseForwardString);
                 
                 
      // *****************************
      // Now render the Hyper Box control                 
      // *****************************

		// Draw "Hyper" out to right
		DrawString(bottomHyperStringX,bottomTopTextY,
					  bottomFont,bottomHyperString);

		if (movieFrameIncrement > 1)
			{
			// Draw "OFF" below Hyper
			DrawString(bottomHyperOFFStringX,bottomBottomTextY+1,
					  bottomFont,bottomHyperOFFString);


         glColor3f(1.0,0.0,0.0); // Red Box
		   glRecti(bottomHyperStringX,bottomGraphicsY,
						bottomHyperStringX+bottomHyperStringWidth-1,bottomGraphicsY+bottomGraphicsHeight);

			glColor3f(1.0,1.0,1.0); // White box around "OFF" Text...
			glBegin(GL_LINE_LOOP);
				glVertex2f(bottomHyperStringX-1,bottomBottomTextY-2);
				glVertex2f(bottomHyperStringX+bottomHyperStringWidth,bottomBottomTextY-2);
				glVertex2f(bottomHyperStringX+bottomHyperStringWidth,bottomBottomTextY+bottomFontHeight+1);
				glVertex2f(bottomHyperStringX-1,bottomBottomTextY+bottomFontHeight+1);
			glEnd();
			}
		}


	// White box around scaling bar
	glColor3f(1.0,1.0,1.0);

	glRecti(bottomScaleLeftX-1,bottomGraphicsY-1,
			bottomScaleLeftX+bottomScaleGraphicWidth-1+1,bottomGraphicsY+bottomGraphicsHeight+1);

	// Draw the long scaling bar - color defined at top of file, typically dark Blue
	glColor3dv(scalingBarColor);

	glRecti(bottomScaleLeftX,bottomGraphicsY,
			bottomScaleLeftX+bottomScaleGraphicWidth-1,bottomGraphicsY+bottomGraphicsHeight);

	// Draw the black line in the bar to indicate relative scale
	glColor3f(0.0,0.0,0.0);
	DrawPercentageBar(
				log(movie.scaleFactor*10),log(movieMinScale*10),log(movieMaxScale*10),
				//movie.scaleFactor,movieMinScale,movieMaxScale,
				bottomScaleLeftX,bottomGraphicsY,
				bottomScaleLeftX+bottomScaleGraphicWidth-1,bottomGraphicsY+bottomGraphicsHeight);

   // *****************************
   // Now draw the speed control bar
   // *****************************

	// White box around Reverse/Forward Area
	glColor3f(1.0,1.0,1.0);
	glRecti(bottomFastReverseX-1,bottomGraphicsY-1,
					bottomFastForwardX2+1,bottomGraphicsY+bottomGraphicsHeight+1);

	// Draw the color areas.
	glColor3f(0.0,1.0,0.0); // Green - (the Go Fast reverse and forward ares)

	glRecti(bottomFastReverseX,bottomGraphicsY,
					bottomSlowReverseX,bottomGraphicsY+bottomGraphicsHeight);
	glRecti(bottomFastForwardX,bottomGraphicsY,
					bottomFastForwardX2,bottomGraphicsY+bottomGraphicsHeight);

	glColor3f(1.0,1.0,0.0); // Yellow - Go Slow Areas
	glRecti(bottomSlowReverseX,bottomGraphicsY,
					bottomStepReverseX,bottomGraphicsY+bottomGraphicsHeight);
	glRecti(bottomSlowForwardX,bottomGraphicsY,
					bottomFastForwardX,bottomGraphicsY+bottomGraphicsHeight);

	glColor3f(1.0,1.0,1.0); // White - Step area
	glRecti(bottomStepReverseX,bottomGraphicsY,
					bottomSlowForwardX,bottomGraphicsY+bottomGraphicsHeight);

	glColor3f(0.0,0.0,0.0);
	if (movieStepMode) // If we are stepping we get a black rectangle in the white area
		{
		int inset = (bottomSlowForwardX-bottomStepForwardX) / 4;
		if (inset < 1)
				inset = 1;
		if (inset > 3)
				inset = 3;

		if (movieDirection == 1)
			{
			glRecti(bottomStepForwardX + inset,bottomGraphicsY+inset,
				bottomSlowForwardX-inset,bottomGraphicsY+bottomGraphicsHeight-inset);
			}
		else
			{
			glRecti(bottomStepReverseX + inset,bottomGraphicsY+inset,
				bottomStepForwardX-inset,bottomGraphicsY+bottomGraphicsHeight-inset);
			}
		}
	else // we are showing the movie continuously - draw a black line (tihn rectangle) in green
		{ // or yellow area to convey relative speed to user.
		if (movieDirection == 1)
			{
			DrawPercentageBar(0.0 + movieMaxMsecs - movieMsecsBetweenFrames,
							movieMinMsecs,movieMaxMsecs,
							bottomSlowForwardX,bottomGraphicsY,
							bottomFastForwardX2,bottomGraphicsY+bottomGraphicsHeight);
			}
		else
			{
			DrawPercentageBar(movieMsecsBetweenFrames,
							movieMinMsecs,movieMaxMsecs,
							bottomFastReverseX,bottomGraphicsY,
							bottomStepReverseX,bottomGraphicsY+bottomGraphicsHeight);
			}
		}

	glColor3f(1.0,1.0,1.0); // White box around HyperBox area
	glBegin(GL_LINE_LOOP);
		glVertex2f(bottomHyperStringX-1,bottomGraphicsY-1);
		glVertex2f(bottomHyperStringX+bottomHyperStringWidth,bottomGraphicsY-1);
		glVertex2f(bottomHyperStringX+bottomHyperStringWidth,bottomGraphicsY+bottomGraphicsHeight+1);
		glVertex2f(bottomHyperStringX-1,bottomGraphicsY+bottomGraphicsHeight+1);
	glEnd();

	// If frames are being skipped, note that here.
  if (movieFrameIncrement > 1)
		{
		char tempstr[10];
		int  tempstrWidth;
		void*  tempstrFont = GLUT_BITMAP_8_BY_13;
		int  tempstrFontHeight = 13;

 		glColor3f(1.0,0.0,0.0); // Red Box
		glRecti(bottomHyperStringX,bottomGraphicsY,
						bottomHyperStringX+bottomHyperStringWidth-1,bottomGraphicsY+bottomGraphicsHeight);

		glColor3f(1.0,1.0,1.0); // White Text
		sprintf(tempstr,"%d",movieFrameIncrement);

		tempstrWidth = StringWidth(tempstrFont,tempstr);

		DrawString(bottomHyperStringX+(bottomHyperStringWidth-tempstrWidth)/2,
					  bottomGraphicsY+2+bottomGraphicsHeight/2-tempstrFontHeight/2,
					  tempstrFont,tempstr);
		}

	glFlush();
	} // End left and right eyes
	glutSwapBuffers();
}



// After the globalScale has been changed, we use this function to make sure it is within
// bounds and to update the movie window.  Recall that scale can be changed with '[' and ']'
// keys in addition to mouse.

void RefreshScale(void)
{
	if (movie.scaleFactor < movieMinScale)
		movie.scaleFactor = movieMinScale;

	if (movie.scaleFactor > movieMaxScale)
		movie.scaleFactor = movieMaxScale;

	if (bottomWindow != -1) // Update scaling bar if we have one
		{
		int current_window = glutGetWindow();

		glutSetWindow(bottomWindow);
		glutPostRedisplay();

		glutSetWindow(current_window);
		}
}

// Cause the bottom window to reflect the new speed bar position.
void RefreshSpeed(void)
{
#if movieMinMsecs > 0
	if (movieMsecsBetweenFrames < movieMinMsecs)
		movieMsecsBetweenFrames = movieMinMsecs;
#endif
	if (movieMsecsBetweenFrames > movieMaxMsecs)
		movieMsecsBetweenFrames = movieMaxMsecs;

	if (bottomWindow != -1) // Update speed control if it is on screen
		{
		int current_window = glutGetWindow();

		glutSetWindow(bottomWindow);
		glutPostRedisplay();

		glutSetWindow(current_window);
		}
}

// In the bottom Window are three UI components.  These routines figure out
// where the mouse was clicked in the bottom subWindow.

int bottomMouseX0;
int bottomMouseY0;
int bottomMouseDownFlag = 0;

// Where was button pressed

enum {nNowhere,nInScale,nInSpeed,nInHyper,nInHyperOFF};  
int bottomMouseDownLocation;

int BottomMouseLocation(int x,int y)
{
	int retval = nNowhere; // Default to nothing (user clicked outside of a control)

	if (y >= bottomGraphicsY && y < bottomGraphicsY + bottomGraphicsHeight)
		{
      // Figure out if user clicked in scale bar, speed bar, or hyper control
      // (or missed all of these)
		if (x >= bottomScaleLeftX && x < bottomScaleLeftX+bottomScaleGraphicWidth)
			retval = nInScale;
		else if (x >= bottomFastReverseX && x <= bottomFastForwardX2)
			retval = nInSpeed;
		else if (x >= bottomHyperStringX && x <= bottomHyperStringX+bottomHyperStringWidth )
			retval = nInHyper;
		}
	else // Was click on OFF of Hyper control?
      if ((y > bottomBottomTextY) && (y < bottomGraphicsY))
			{ 
			if (x >= bottomHyperStringX && x <= bottomHyperStringX+bottomHyperStringWidth )
				{
				retval = nInHyperOFF;
				}
			}

	return retval;

}

// The mouse can be dragged around or clicked - and we need to update
// the user interface components either way.
void BottomHandleMouse(int mouseLocation,int glutMouseState,int x)
{
   static int movieStepIgnoreUp; // See below
   
   // Depending on where the mouse was clicked, handle the UI component
	switch (mouseLocation) 
      {
      // For a mouse event in the scale bar, an exponential scaling of 
      // the click location feels "right"
		case nInScale:
			movie.scaleFactor = (float)
				exp(BarValue(x,log(movieMinScale*10.0),log(movieMaxScale*10.0),
								bottomScaleLeftX,bottomScaleLeftX+bottomScaleGraphicWidth-1))/10.0f;

			// If user is trying for 1.0 help him out a bit.
			if (movie.scaleFactor >= .97 && movie.scaleFactor <= 1.03)
					movie.scaleFactor = 1.0;
			RefreshScale();
			MovieAdvanceTimerKeyAndPostRedisplay();
			break;

      // In the speed bar, we need to figure out if it is stepwise motion that
      // has been requested.
		case nInSpeed:
			// Was the click in the Slow or Fast Forward Area?
			if (x >= bottomSlowForwardX)
				{
				movieDirection = 1;

				movieMsecsBetweenFrames = movieMaxMsecs + 0.0 -
					BarValue(x,movieMinMsecs,movieMaxMsecs,
									bottomSlowForwardX,bottomFastForwardX2);

				movieStepMode = 0;
				MovieAdvanceTimerKeyAndPostRedisplay();
				}
			// Or, was the click in the Slow or Fast Reverse Area?
			else if (x < bottomStepReverseX)
				{
				movieDirection = -1;

				movieMsecsBetweenFrames =
					BarValue(x,movieMinMsecs,movieMaxMsecs,
									bottomFastReverseX,bottomStepReverseX-1);

				movieStepMode = 0;
				MovieAdvanceTimerKeyAndPostRedisplay();
				}

			// Or, Are we in the Slow Step Forward Area
			else if (x >= bottomStepForwardX)
				{
				// In the step area, only "advance a frame" if we click and release
				// in this box - and we were already in stepMode
				if ((glutMouseState == GLUT_UP) && (movieStepMode))
					{
					if (movieStepIgnoreUp)
						movieStepIgnoreUp = 0; // Next time is a keeper...
					else
						{
						// Just make sure our down click was in this box!
					   if (bottomMouseX0 >= bottomStepReverseX &&
							 bottomMouseX0 < bottomSlowForwardX)
							{
							MovieNextFrame();
							glutSetWindow(movieWindow);
							glutPostRedisplay();
							glutSetWindow(bottomWindow);
							}
						}
					}

				// But this will bring us to a stop for sure when we drag or click in here...
				if ((movieDirection != 1) || (! movieStepMode))
					{
					movieDirection = 1;
					if (! movieStepMode)
						movieStepIgnoreUp = 1;
					movieStepMode = 1;
					glutSetWindow(movieWindow);
					glutPostRedisplay();
					glutSetWindow(bottomWindow);
					}
				}
			// Are we in the Slow Step Reverse Area
			else if (x >= bottomStepReverseX)
				{
				// In the step area, only "advance a frame" if we click and release
				// in this box - and we are already in step mode....
				if ((glutMouseState == GLUT_UP) && (movieStepMode))
					{
					if (movieStepIgnoreUp)
						movieStepIgnoreUp = 0; // Next time is a keeper...
					else
						{
						// Just make sure our down click was in a step box!
					   if (bottomMouseX0 >= bottomStepReverseX &&
							 bottomMouseX0 < bottomSlowForwardX)
							{
							MovieNextFrame();
							glutSetWindow(movieWindow);
							glutPostRedisplay();
							glutSetWindow(bottomWindow);
							}
						}
					}

				// But this will bring us to a stop for sure.
				// as we drag mouse through the bar...
				if ((movieDirection != -1) || (! movieStepMode))
					{
					if (! movieStepMode)
						movieStepIgnoreUp = 1;
					movieDirection = -1;
					movieStepMode = 1;
					glutSetWindow(movieWindow);
					glutPostRedisplay();
					glutSetWindow(bottomWindow);
					}
				}

			RefreshSpeed();
			break;

      // In the red HyperBox, means we increment the movieFrameIncrement. 
		case nInHyper:
			if (glutMouseState == GLUT_UP &&
				 movieFrameIncrement < (displayMovieMaxFrame()-displayMovieMinFrame()-1))
				{
				movieFrameIncrement++;
				glutPostRedisplay();
				MovieAdvanceTimerKeyAndPostRedisplay();
				}
			break;

		case nInHyperOFF:
         // In the "OFF" Box below the Hyper Box means we go back to a movieFrame Increment
         // of 1 - which means every frame will be displayed.
			if (glutMouseState == GLUT_UP &&
				 movieFrameIncrement > 1)
				{
				movieFrameIncrement = 1;
				glutPostRedisplay();
				MovieAdvanceTimerKeyAndPostRedisplay();
				}
			break;

		case nNowhere:
		default:
			break;

		} // End Switch
}

// This function is called by GLUT when there is a mouse
// clock in the bottom window.
void BottomMouseFunc(int button,int state,int x,int y)
{
	// Convert "mouse" y coordinate to OpenGL coordinate!
	y=bottomWindowHeight - y - 1;

	if (button == GLUT_LEFT_BUTTON)
		{
		if (state == GLUT_DOWN)
			{
			bottomMouseX0 = x;
			bottomMouseY0 = y;

			bottomMouseDownLocation = BottomMouseLocation(x,y);

         BottomHandleMouse(bottomMouseDownLocation,state,x);
			bottomMouseDownFlag = 1;
			}
		else
			{
			bottomMouseDownFlag = 0;

			if (state == GLUT_UP)
				{
				if (bottomMouseDownLocation == BottomMouseLocation(x,y))
               BottomHandleMouse(bottomMouseDownLocation,state,x);
				else // If we release in a "funny" place, go back to original
					{ // Ignore stray mouse if Hyper
					if (bottomMouseDownLocation != nInHyper)
                  BottomHandleMouse(bottomMouseDownLocation,state,bottomMouseX0);
					}
				}
			}

		}
}

// If the mouse is dragged in the bottom area, then update the current UI component
// with it.

void BottomMotionFunc(int x,int y)
{
	if (bottomMouseDownFlag && (BottomMouseLocation(x,y) == bottomMouseDownLocation))
      BottomHandleMouse(bottomMouseDownLocation,-1000,x);
}

// GLUT calls this routine when the Movie subwindow is resized.
// There are some global variables setup to make other code
// compatible with earlier GLUT versions.
void MovieChangeSizeFunc(int w,int h)
{
	movieWindowHeight  = h;
	movieWindowWidth  = w;
	movieWindowAspectRatio = ((float)w)/((float)h);
	{
	// For compatability with earlier display, we want a range of SGI-like xy values
	// of at least 1280 x 1024.  But, we also want a lot more flexibility with
	// window size than the rigid "fixedaspect" of earlier MD Display version.

	float x_fixup_required = 1280.0/movieWindowWidth;
	float y_fixup_required = 1024.0/movieWindowHeight;

	// Pick the max of the two... and this will ensure that our drawing space
	// is good for the range we are interested in...
	movieWindowFixupFactor = (x_fixup_required > y_fixup_required) ?
					x_fixup_required : y_fixup_required;
	}

	// Reset coordinate System
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
}

// *****************************************************************
// Movie Window Mouse Management
//
// Rotation, translation, scaling, and picking of elements are
// all possible in the movie Subwindow.
// The variables and functions below manage this
// *****************************************************************

int movieY0,movieX0; // Coordinate of last MOUSE_DOWN event

// One challenge is that if the user presses two buttons simultaneously, we get
// two separate events - so we have to be a little clever!

#define NO_BUTTON_DOWN (-2098) // A Zany number clearly distinct
// Set at mouse click in movie window - and then used when mouse is dragged
// for rotation, translation, or scaling
int movieButtonDown = NO_BUTTON_DOWN; 
int movieMouseScaleMode = 0;

// Some Rotation modes
enum {nXYRot,nZRot1,nZRot2,nZRot3}; 

// When the user requests a new rotation (or translation), we
// need to add this new reauest to the existing rotation 
// (or translation)
float movieOldXRotation,movieOldYRotation,movieOldZRotation;
float movieOldXTranslation,movieOldYTranslation;
float movieOldScaleFactor;

// The current state of the scene rotation
double movieXRotation; // radians about X axis
double movieYRotation; // "" about Y axis
double movieZRotation; // "" about Z axis
int movieRotationMode = -1;

void MovieMouseButtonRelease(void)
{
	// Note that we don't use "button" parameter here...
	// Basically - we just go back to the "no buttons down" world.
	if (movieButtonDown != NO_BUTTON_DOWN)
		{
		glutSetCursor(movieCursor = GLUT_CURSOR_CROSSHAIR);
		movieRotationMode = -1;
		movieX0 = movieY0 = 0;
		MovieAdvanceTimerKeyAndPostRedisplay();
		movieButtonDown = NO_BUTTON_DOWN;
		movieMouseScaleMode = 0;
		}
}

// See Open GL Reference Manual "selectObjects()" routine..
#define nSelectBufsize 40	

void MoviePickClosestElement(int mouseX,int mouseY)
{
	GLuint selectBuf[nSelectBufsize]; // A buffer of objects which are near the mouse click.
	GLint hits;
	int i;

	// Intialize the selection buffer to invalid object ids
	for (i=0;i< nSelectBufsize;i++)
		selectBuf[i] = 999999;

   // Prepare teh selectBuffer
	glRenderMode(GL_RENDER);
	glSelectBuffer(nSelectBufsize,selectBuf);

	{
		int err = glGetError();
		if (err != GL_NO_ERROR)
			{
			printf("Terminated inPick ClosestElement due to OpenGL error %d %s\n",err,gluErrorString(err));
			exit(EXIT_FAILURE);
			}
	}

	glRenderMode(GL_SELECT);
	assert(glGetError() == GL_NO_ERROR);


	glInitNames();
	assert(glGetError() == GL_NO_ERROR);
	glPushName(0);

	assert(glGetError() == GL_NO_ERROR);

   // Render the scene - but in GL_SELECT mode.  Nothing happens
   // on the screen.  Rather, we get a list of elements which are near
   // the mouse click, within the 6.0 tolerance
	MovieRenderOrSelectMainScene(GL_SELECT,mouseX,mouseY,6.0,6.0);
	assert(glGetError() == GL_NO_ERROR);

	glPopName();

   // Now we may have some "hits" (elements near the mouse)
	hits = glRenderMode(GL_SELECT);
	glRenderMode(GL_RENDER);
   // Let display* code decide what to do with them.
	displayProcessHits(hits,selectBuf);
}

#undef nSelectBufsize

// The program determines whether rotation is about Z or XY based
// location of the mouse at the time that the right button is pressed
// So, where is the center?  If we have a pair of images the rotation center
// point wil be in the middle of th  the left or right half of the movie Window.

void MovieCenterPixel(int* centerX, int* centerY,int mouseX)
{
	if (movie.imageView == nPairOfImages)
      {
		if (mouseX < movieWindowWidth/2)
			{
			*centerX = 
      		StereoPairsViewPortWidth(movieWindowWidth)/2;
			}
      else 
			{
         *centerX = (movieWindowWidth - 
                    	StereoPairsViewPortWidth(movieWindowWidth)/2);
			}
		}
	else
		{
		*centerX = movieWindowWidth/2;
		}
	*centerY = movieWindowHeight/2;
}

// How far are we from center point?
void MovieMouseClickPixelsFromCenter(
	int* x0_from_center,
	int* y0_from_center,
	int mouseX,int mouseY,
   int mouseXForCenter)
{
	int centerX,centerY;
   MovieCenterPixel(&centerX, &centerY,mouseXForCenter);

	*x0_from_center = mouseX - centerX;
	*y0_from_center = mouseY - centerY;
	
   if (movie.imageView == nSGILegacyHardwareStereo)
           *x0_from_center /= 2;
	*y0_from_center = mouseY - movieWindowHeight/2;
}

// When the user does a rotation, we have to, in a sense, append the new
// user requested rotation to the one already in play.  The saving
// of the old matrix happens when the user goes into rotation mode.
// The appending happens when the user decides on a final rotation.
// Putting this intermediate data in a global saves some headaches
float movieOldRotationMatrix[16];

// GLUT calls this function when the mouse is pressed or released
// in the movie window.
void MovieMouseFunc(int button, int state,int x, int y)
{
if (state == GLUT_DOWN)
	{
	// If there is already a button down, then user is _adding_ a button...
	// So, we need to close down current mode, and crank up "scaling" mode.
   // I.e. - any combination of 2 buttons causes scaling mode
	if (movieButtonDown != NO_BUTTON_DOWN)
		{
		MovieMouseButtonRelease();
		movieMouseScaleMode = 1;
		}

	movieX0 = x;

	if (movie.imageView == nSGILegacyHardwareStereo)
		{
		//printf("Stereo Click(%d,%d)\n",x,y);
		movieY0 = movieWindowHeight - y - 1;
		if (movieY0 < 0)
			movieY0 = 0;
		//printf("Translated Stereo Click(%d,%d)\n",movieX0,movieY0);
		}
	else
		movieY0 = movieWindowHeight - y - 1; // Convert to "OpenGL" coordinate


	movieButtonDown = button;

	// Seems that pressing the button is a good time to turn off mouse
	// But, if this is the right button, flag is reversed below
	movieTrackBallCircleVisible = 0;

   // Turn off the cursor if we are in scaling mode.
	if (movieMouseScaleMode)
		{
		movieOldScaleFactor = movie.scaleFactor;
		glutSetCursor(movieCursor = GLUT_CURSOR_NONE);
		}
	else
		{
		// On a two button mouse system, pressing SHIFT with left mouse gives middle
		if (button == GLUT_LEFT_BUTTON)
			{
			if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
				movieButtonDown = GLUT_MIDDLE_BUTTON;
			else
				{
            // We have pressed the left mouse button
				if (state == GLUT_DOWN)
					{
					// The idea here is that if display can just "handle" the click
					// right away (as in removal of distance metric) then we are
					// done with this click!
					if ( ! displayProcessRawMovieClick(movieX0,movieY0,movieWindowWidth,movieWindowHeight) )
                  {
                  // Otherwise, assume a picking operation is going on.
						MoviePickClosestElement(movieX0,movieY0);
                  }
					if (movie.hideLeftAndBottom)
						MovieHideLeftAndBottom(! movie.hideLeftAndBottom);
					}
				}
			}
         
      // If the right button is pressed, then we have
      // initiated a rotation of some kind.  The specific
      // kind depends on how far from the center that the
      // right button was pressed.
		if (movieButtonDown == GLUT_RIGHT_BUTTON)
			{
			movieOldXRotation = movieXRotation;
			movieOldYRotation = movieYRotation;
			movieOldZRotation = movieZRotation;

			// Save the "old" rotation matrix
			memcpy(movieOldRotationMatrix,movie.rotationMatrix,sizeof(movieOldRotationMatrix));

         // Default is a rotation around X and Y axes - when user clicks near
         // center of image.
			movieRotationMode = nXYRot;

			{
			int x0_from_center;
			int y0_from_center;
			// Compute mouse click offset from "center"
			MovieMouseClickPixelsFromCenter(
					&x0_from_center,
					&y0_from_center,
					movieX0,
               movieY0,movieX0);

			if (movie.imageView == nThreeImages)  // Are we in F3 (TRI) mode?
				{
				// No need to do anything - center is center of main window!
            // In this mode, rotation is always the nXYRot flavir
				}
			else
				{
				float dist_from_center;

				dist_from_center = sqrt(x0_from_center*x0_from_center + y0_from_center*y0_from_center);

            // If the user clicks outside the blue circle, then the rotation is going to be
            // around the Z axis
				if (dist_from_center > 300.0 / movieWindowFixupFactor)
					movieRotationMode = nZRot1;
				}
			}

         // Changing the cursor to reflect the rotation type has a very nice look
			if (movieRotationMode == nXYRot)
				glutSetCursor(movieCursor = GLUT_CURSOR_CROSSHAIR);
			else
				glutSetCursor(movieCursor = GLUT_CURSOR_CYCLE);

         // Display the blue disk for user reference.
			if (movie.imageView != nThreeImages)
				{
				movieTrackBallCircleVisible = 1;
				glutPostRedisplay();
				}
			}
		else // If the middle button is pressed, then the user is translating
         if (movieButtonDown == GLUT_MIDDLE_BUTTON)
   			{
   			movieOldXTranslation = movie.xTranslation;
   			movieOldYTranslation = movie.yTranslation;
   			glutSetCursor(movieCursor = GLUT_CURSOR_INFO); // A hand for translation mode... Nice..
   			}
		}
	} // End if (state == GLUT_DOWN)
else
if (state == GLUT_UP)
	{
	MovieMouseButtonRelease();
	}
}


// GLUT calls this function when the user is dragging the mouse
// The mode of scaling, translation, or rotation will have already been set
// in the above code when the mouse button(s) were clicked down.
void MovieMotionFunc(int x,int y) 
{
   // Distances mouse has moved since original mouse down
	int dmx = x-movieX0; // x distance
	int dmy;             // y distance
	GLint current_matrix_mode;

	y = movieWindowHeight-y-1; // Reverse y coordinate for compatability with old code... 0 is not at bottom.

	dmy = y-movieY0;
   
  	if (movie.imageView == nSGILegacyHardwareStereo)
   	dmy *= 2.0;

if (movieMouseScaleMode) // Then we are scaling
	{
   // Rescale the movie
	movie.scaleFactor = movieOldScaleFactor * (1.0 + 3.0 * (float)(y-movieY0)/(movieWindowHeight));
   // Refresh the scale bar on the bottom of the screen
	RefreshScale();
   // Force a redisplay of the current movie frame
	MovieAdvanceTimerKeyAndPostRedisplay();
	}
else if (movieButtonDown == GLUT_RIGHT_BUTTON) // Then we are rotating!
	{
	glGetIntegerv(GL_MATRIX_MODE,&current_matrix_mode);

	// Here we closely copy original display code which results in a
	// "post multiply" matrix instead of pre multiply...
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if (movieRotationMode >= nZRot1)
		{
		int x1,y1,x0,y0;

		if (movieRotationMode == nZRot1)
			{
         MovieMouseClickPixelsFromCenter(&x1,&y1,movieX0,movieY0,movieX0);
         MovieMouseClickPixelsFromCenter(&x0,&y0,x,y,movieX0);
			}

      // Do some trig to figure angle that charactizes to change in position
		if (fabs(x0)>10. || fabs(y0)>10.)
			{
		   double ang1 = atan2(y1,x1);
			double ang0 = atan2(y0,x0);

			movieZRotation = (360.0/(2.0*M_PI))*(ang0-ang1) + movieOldZRotation;
			glRotatef((360.0/(2.0*M_PI))*(ang0-ang1),0.0,0.0,1.0);
			}
		}
	else // XYRot mode
		{
		// We need to scale our physical movie window width and height
		// as we translate mouse moves to old rotations!
		// printf("%d %d\n",dmx,dmy);
		float xrot,yrot;
		movieXRotation = dmy + movieOldXRotation;
		movieYRotation = dmx + movieOldYRotation;
		xrot = .3*(double)movieWindowWidth/1279.0;
		yrot = .3*(double)movieWindowHeight/1023.0;

		if (movie.imageView == nThreeImages)
			{
         glRotatef(dmx*xrot/2.0,0.0,1.0,0.0);
   		glRotatef(-dmy*yrot/2.0,1.0,0.0,0.0);
			}
		else
			{
			glRotatef(dmx*xrot,0.0,1.0,0.0);
			glRotatef(-dmy*yrot,1.0,0.0,0.0);
			}
		}

	glMultMatrixf(movieOldRotationMatrix);
   
// Now replace old movie.rotationMatrix with this one.
	glGetFloatv(GL_MODELVIEW_MATRIX,movie.rotationMatrix);

	glMatrixMode(current_matrix_mode);

   // Post a redisplay event to th emovie Window
	MovieAdvanceTimerKeyAndPostRedisplay();

	// Update the axes display
	if (leftWindow != -1)
		{
		glutSetWindow(leftWindow);
		glutPostRedisplay();
		}
	} // End Right button down (rotation)
else if (movieButtonDown == GLUT_MIDDLE_BUTTON) // Then we are translating.
	{
   // Just update the translation x and y - and then post a redisplay to
   // the movie window
	movie.xTranslation = movieOldXTranslation + (float)(dmx)*movieImageSize/(1000.0*movie.scaleFactor*movie.scaleFactor);
	movie.yTranslation = movieOldYTranslation + (float)(dmy)*movieImageSize/(1000.0*movie.scaleFactor*movie.scaleFactor);
	MovieAdvanceTimerKeyAndPostRedisplay();
	}
}

// Set global variables for a 12 point Helvetica font
// in the menu boxes
static void MenuSettingsHelvetica12_First(void)
{
	// First, try a 12 point font...  Set the globals and see if they fit!
	menuFontHeight = 12;
	menuFont = GLUT_BITMAP_HELVETICA_12;
	menuTop = 85;
	menuBoxHeight = 27;
	menuBoxSpace = 30;		// Spacing from a box - to box below
	menuBoxWidth = 74;		// Width of the box
	menuBoxLeft = 3;		// Starting X coordinate of Menu box in leftWindow
	leftWindowWidth = menuBoxWidth + menuBoxLeft * 2;
	leftWindowHeight = (nMenuItemCount * menuBoxSpace)  // The menu buttons
							  + menuTop
							  + 3;
							   // Allow for 4x1 aspect clip viewer
}

// Set global variables for a 10 point Helvetica font
// in the menu boxes
static void MenuSettingsHelvetica10_Second(void)
{
	menuFontHeight = 10;
	menuFont = GLUT_BITMAP_HELVETICA_10;
	menuTop = 75;
	menuBoxHeight = 25;
	menuBoxSpace = 27;		// Spacing from a box - to box below
	menuBoxWidth = 64;		// Width of the box
	menuBoxLeft = 2;		// Starting X coordinate of Menu box in leftWindow
	leftWindowWidth = menuBoxWidth + menuBoxLeft * 2;
	leftWindowHeight = (nMenuItemCount * menuBoxSpace) + menuTop + 4;
}

// Set global variables for a 10 point Helvetica font
// but tightest possible fit into the smallest of menu boxes
void MenuSettingsHelvetica10_Third(void)
{
	menuFontHeight = 10;
	menuFont = GLUT_BITMAP_HELVETICA_10;
	menuTop = 60;
	menuBoxHeight = 20;
	menuBoxSpace = 22;		// Spacing from a box - to box below
	menuBoxWidth = 60;		// Width of the box
	menuBoxLeft = 1;		// Starting X coordinate of Menu box in leftWindow
	leftWindowWidth = menuBoxWidth + menuBoxLeft * 2;
	leftWindowHeight = (nMenuItemCount * menuBoxSpace) + menuTop + 4;
}

// This spacing is good for SGI Legacy stereo
void MenuSettingsForStereo(void)
{
	menuFontHeight = 10;
	menuFont = GLUT_BITMAP_HELVETICA_10;
	menuTop = 37;
	menuBoxHeight = 14;
	menuBoxSpace = 15;		// Spacing from a box - to box below
	menuBoxWidth = 66;		// Width of the box
	menuBoxLeft = 2;		// Starting X coordinate of Menu box in leftWindow
	leftWindowWidth = menuBoxWidth + menuBoxLeft * 2;
	leftWindowHeight = (nMenuItemCount * menuBoxSpace) + menuTop + 4;
}

// This large function manges all the gory details of placing subwindows on the Master window
// It is called, for example, when the user chooses the 'hide' option to increase
// screen display.
void MasterCreateOrChangeAllSubwindows(void)
{
	int w,h;
   
   // If our current_window gets deleted, we will take care
   // to not attemp to "restore it"
	int save_current_window = glutGetWindow();
   int clipWindowActive = (clipWindow != -1);
   
   // Lets start with an array of pointers to subWindows.

   static int* window_addresses[] =
		{&leftWindow,&staticMenuButtonWindow,&movieWindow,&paletteWindow,&bottomWindow,&clipWindow,NULL};

	int* current_window_ptr;

   // Goal here is to record the address of a GLUT window handle which 
   // we might be destroying below - and make that is current after OK 
   // recreation! 	

	current_window_ptr = window_addresses[0];
	while ((current_window_ptr != 0) && 
			 (*current_window_ptr != save_current_window))
		{
		current_window_ptr++;
		}


   // Let's start by getting some info about the Master Window
   // which is the big black unseen window underlying all 
   // subWindows
	glutSetWindow(masterWindow);
	w = glutGet(GLUT_WINDOW_WIDTH);
	h = glutGet(GLUT_WINDOW_HEIGHT);

#ifdef SGI_STEREO
	if (movie.imageView == nSGILegacyHardwareStereo)
		h = (h - VPSEP) / 2;
#endif


	// We'll want a leftWindow for Menu buttons on left side of screen.
   // Destroy it if it is already active.
   // Also, destroy the staticMenuButton subWindow
	if (leftWindow != -1)
			{
			glutSetWindow(staticMenuButtonWindow);
			glutDetachMenu(GLUT_RIGHT_BUTTON);
			glutDetachMenu(GLUT_MIDDLE_BUTTON);
			glutDetachMenu(GLUT_LEFT_BUTTON);
			glutDestroyMenu(staticMenuIdentifier);
			glutDestroyWindow(staticMenuButtonWindow);

			glutDestroyWindow(leftWindow);

			leftWindow = staticMenuButtonWindow = -1;
			}

#ifdef APPLE_SEEMS_TO_NEED_THIS_BLOCK
	if (movieWindow != -1)
		{
		DeleteMovieTrackBallCircleQuadric();
		glutDestroyWindow(movieWindow);
		movieWindow = -1;
		}
#endif

	// If a clipWindow is displayed (always at left below leftWindow)
   // that needs to go too.

	if (clipWindowActive)
			{
			glutDestroyWindow(clipWindow);
			clipWindow = -1;
			}

   
   // If active, kill off the bottom window where speed, scale, and hyper
   // controls are found
	if (bottomWindow != -1)
		{
		glutDestroyWindow(bottomWindow);
		bottomWindow = -1;
		}

   // If the movie.hide... flag is set at this point, then no 
   // components will be created on the masterWindow - other than the movieWindow
	if (! movie.hideLeftAndBottom)
	{

/*********************************************************************************
  If room, we'll allocate a subwindow on the left side for the Menu Buttons
  If the Window is tall and wide, use a nice 12 point font.  Otherwise, scale down to 10
  point... or get rid of altogether.
  
  SGI Legacy stereo demands its own "funny" font and box heights
*********************************************************************************/


	if (movie.imageView == nSGILegacyHardwareStereo)
		MenuSettingsForStereo();
	else
	{
	MenuSettingsHelvetica12_First();

	// If we don't have enought room in this window - scale down font.
	// and try again.  Idea: We want 33% of window height for Clipping
	if ((h < leftWindowHeight + (h/3)) ||
		 (w < leftWindowWidth*3)) // *3 is for aesthtics
		{
		MenuSettingsHelvetica10_Second();
		}

	// If we don't have enought room in this window - scale down menu buttons.
	// and try again.  We can stand 20% if that is as good as we can do.
	if ((h < leftWindowHeight + (h/5)) ||
		 (w < leftWindowWidth*3)) // *3 is for aesthtics
		{
		MenuSettingsHelvetica10_Third();
		}
	}


   // If we can fit the menu items - then let's have a left subWindow
   // with a static subWindow in the first position.
	if ((h >= leftWindowHeight + 10) && 
		 (w >= leftWindowWidth*3))
		{
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
		leftWindow = glutCreateSubWindow(masterWindow,0,0,leftWindowWidth,leftWindowHeight);
		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);

		glutDisplayFunc(leftDisplayFunc);
		glutReshapeFunc(leftChangeSizeFunc);
		glutMouseFunc(LeftMouseFunc);

		staticMenuButtonWindow = glutCreateSubWindow(
			leftWindow,
			menuBoxLeft,
			nStatic*menuBoxSpace + menuTop,
			menuBoxWidth,
			menuBoxHeight);
				
		glutKeyboardFunc(MasterKeyboardFunc);
		glutSpecialFunc(MasterSpecialFunc);
		glutDisplayFunc(staticDisplayFunc);
		glutReshapeFunc(staticChangeSizeFunc);
      
      // For the static subWindow, we'll attach a GLUT menu system.
		staticMenuIdentifier = glutCreateMenu(staticMenuCallback);
		glutAddMenuEntry("First",nStaticFirst);
		glutAddMenuEntry("Average",nStaticAverage);
		glutAddMenuEntry("None",nStaticNone);
		glutAttachMenu(GLUT_LEFT_BUTTON);
		glutAttachMenu(GLUT_RIGHT_BUTTON);
		glutAttachMenu(GLUT_MIDDLE_BUTTON);
		}
	else // Could not fit menus...
		{
		leftWindowWidth = leftWindowHeight = 0;
		}

/*********************************************************************************
  If room, we'll allocate a subwindow on the bottom side for scaling, speed, and Hyper
  UI elements.  This subwindow lies to the right of the menu window.

  As above, if the Window is too small, we'll axe these altogether.
*********************************************************************************/


	if ((w > 400) && (h > 100))
		{
		int success = 0;

		if (movie.imageView == nSGILegacyHardwareStereo)
			success = BottomPlaceGraphics(w-leftWindowWidth,h,GLUT_BITMAP_HELVETICA_10,6);

		if (! success)
         success  = BottomPlaceGraphics(w-leftWindowWidth,h,GLUT_BITMAP_HELVETICA_12,12);

		if ( ! success )
			success = BottomPlaceGraphics(w-leftWindowWidth,h,GLUT_BITMAP_HELVETICA_10,10);

      // If we can't fit the bottom graphics with text, then just put the UI elements
      // only (no text)
		if ( ! success )
			success = BottomPlaceGraphics(w-leftWindowWidth,h,0,0);

		if ( success ) // great let's have a bottom subWindow
			{
			glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

			bottomWindow =
				glutCreateSubWindow(
					masterWindow,
					leftWindowWidth,     // Start subwin to right of left menu win
					h-bottomWindowHeight,
					w-leftWindowWidth,   // Finish on right of screen
					bottomWindowHeight);

         // connect Keyboard Handler to our global handler
			glutKeyboardFunc(MasterKeyboardFunc);
			glutSpecialFunc(MasterSpecialFunc);

         // provide GLUT callbacks for usual stuff
			glutMouseFunc(BottomMouseFunc);
			glutMotionFunc(BottomMotionFunc);
			glutDisplayFunc(BottomDisplayFunc);
			glutReshapeFunc(BottomChangeSizeFunc);
			}
		}
	} // end if (! masterHideLeftAndBottom)
	else
	{
   // the user pressed 'h' - so no left or bottom subWindows!
	leftWindowWidth = leftWindowHeight = bottomWindowHeight = 0;
	}

   // If we already have a movieWindow, then we need only reposition it.
	if (movieWindow != -1)
		{
		glutSetWindow(movieWindow);
		glutPositionWindow(leftWindowWidth,0);
		glutReshapeWindow(w-leftWindowWidth,h-bottomWindowHeight);
		if ((paletteWindow != -1) && (leftWindow == -1))
			{
			glutDestroyWindow(paletteWindow);
			paletteWindow = -1;
			}
		}
	else // Create Movie subWindow on Master Winow for First Time
		{
		unsigned int quad_stereo_bits = quadStereoAvailable && 
         (movie.imageView == nQuadStereo) ? GLUT_STEREO : 0;
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | quad_stereo_bits);
		movieWindow = glutCreateSubWindow(
                  masterWindow,
                  leftWindowWidth,
                  0,
                  w-leftWindowWidth,
                  h-bottomWindowHeight);

		SetupMovieTrackBallCircleQuadric(); 
      // This function specifies GLUT callbacks for the movie Window.
      dispuiSetMovieFunctions();
		}

	// If the clipping tool was in use - we need to resize it
	if ((clipWindowActive) && (leftWindow != -1))
		{
		// Then we need the clip window back again.
		ClipWindowCreateOrChange();
		}

	// If paletteWindow was in use - then we need to resize it
	if (paletteWindow != -1)
		{
		if (leftWindow == -1)  // If there is no room or desire for a "left Window" then
			PaletteWindowToggle(); // This will destroy the clipping window - which is what we need.
		PaletteWindowCreateOrChange();
		}

   // Posting a redisplay to the master should cause all subwindows to get a redisplay.
	glutSetWindow(masterWindow);
	glutPostRedisplay();


   // If we were able to get the full current window information
   // at outset, and if that window has not been destroyed,
   // then make that window the current one now!
	if (current_window_ptr) // If we could identify window up top
		{
		if	(*current_window_ptr != -1) // And not destroyed?
			glutSetWindow(*current_window_ptr);
		}
	else // No problem - must have been a window that we would not have destroyed
		{
		glutSetWindow(save_current_window);
		}
}



// GLUT calls this function when the user resizes the entire active window
// Also it is called in the switch to and from SGI legacy stereo mode.
void MasterChangeSizeFunc(int w, int h)
{
	if (h == 0)
		h= 1;

	// This odd code ensures that we _were_ not full screen
	// before we do the switch to stereo!
	if (SGIStereoRequested)
		{
#ifdef SGI_STEREO   
		glutSetWindow(masterWindow);
		glutFullScreen();
      if (movie.imageView == nSGILegacyHardwareStereo)
   		start_fullscreen_stereo();
#endif         
		SGIStereoRequested = 0;
		}

	if (SGIStereoOffToFullScreenRequested)
		{
#ifdef SGI_STEREO   
      if (movie.imageView == nSGILegacyHardwareStereo)
   		stop_fullscreen_stereo();
		MasterSetFullScreen(1);
#endif         
		SGIStereoOffToFullScreenRequested = 0;
		}

	// This means we don't do a bunch of work on the
	// "Fake" stereo resize
	if (w > 1 && h > 1)
		MasterCreateOrChangeAllSubwindows();
}

// GLUT will rarely call the master display function - since all action
// is in the subWindows.  But, keep it black for good measure.
void MasterDisplayFunc(void)
{
	int eye;
	for (eye = 0;eye < SGIEyeCount(); eye++)
		{
		dispuiSelectEye(eye);
	// printf("MasterDisplayFunc called\n");
	// Don't do anything!
	// All the subwindows do the work!
		glClearColor(0.0,0.0,0.0,0.0); // BLACK
		glClear(GL_COLOR_BUFFER_BIT);
		glFlush();
		}
	glutSwapBuffers();
}

// Shutdown the dispui system.  Restores the
// screen if in SGI stereo mode.  This is called
// atexit() - which is set below.
void dispuiShutDown(void)
{
#ifdef SGI_STEREO
	if (movie.imageView == nSGILegacyHardwareStereo)
      stop_fullscreen_stereo();
#endif      
}

static int dispui_InitCalledFlag = 0;

// dispuiMain is called from yourdisplay.c after it parses out arguments.
// From here, glutMainLoop is called and everything else is a result of
// callbacks from the glutLibrary - a completely event driven architecture
void dispuiMain(const char* winTitle)
{
	int screen_w,screen_h;

	assert(dispui_InitCalledFlag);

	// Let's make start window about 90% of screen real estate.
	screen_w = glutGet(GLUT_SCREEN_WIDTH);
	screen_h = glutGet(GLUT_SCREEN_HEIGHT);
	glutInitWindowSize((int)(screen_w * .9),
							 (int)(screen_h * .9));

	glutInitWindowPosition((int)((.1 * screen_w)/2.0),		// Center left to right
								  (int)((.05 * screen_w)/2.0));	// Close to top - but not completely there


	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

	masterWindow = glutCreateWindow(winTitle);
   
   // Setup GLUT callbacks for the Master window.
	glutDisplayFunc(MasterDisplayFunc);
	glutReshapeFunc(MasterChangeSizeFunc);
	glutKeyboardFunc(MasterKeyboardFunc);
	glutSpecialFunc(MasterSpecialFunc);

   // Setup the blue rotation-cue circle for subsequent
  	// display during user rotations.

#ifdef SGI_STEREO
   if (movie.imageView == nSGILegacyHardwareStereo)
		SGIStereoRequested = 1;
#endif      
	glutMainLoop();
}

void dispuiInit(float _movieImageSize,
					 int _movieFirstFrame,
#ifndef DISPUI_TESTMODE                
					 float picoSecondsInitialOffset,
					 float picoSecondsPerFrame,
#endif                
					 int quad_stereo_flag)
{
	assert(! dispui_InitCalledFlag);
	movieImageSize = _movieImageSize;
	movieFrame = _movieFirstFrame;
   quadStereoAvailable = quad_stereo_flag;

	if (movieFrame < displayMovieMinFrame())
		movieFrame = displayMovieMinFrame();

#ifndef DISPUI_TESTMODE                
	moviePicoSecondsInitialOffset = picoSecondsInitialOffset;
	moviePicoSecondsPerFrame = picoSecondsPerFrame;
#endif

	dispui_InitCalledFlag = 1;
	atexit(dispuiShutDown);
}



