/*
 *  MD Display Version 3.0
 *
 *  function and typedefs for dispui.c:  user interface management
 *
 *  See comments in dispui.c
 */


/* Before dispuiMain is called, you must call dispuiInit 
   with some preliminary info */

void dispuiInit(float _movieImageSize,
					 int _movieFirstFrame,
#ifndef DISPUI_TESTMODE
					 float picoSecondsInitialOffset,
					 float picoSecondsPerFrame,
#endif
					 int _quad_stereo_flag);
 
// Start GLUT, Creates Windows, goes to work on the "UI" 
// Once your application calls this, it only gets callbacks
// to the display* functions as dispui.c needs to make them.
void dispuiMain(const char* winTitle);     


// A global message that will be displayed prominently in 
// the movie window in case of user input errors.
extern char   dispuiStatusMessage[80];
extern time_t dispuiStatusMessageExpiration;

// Set the rotation matrix to identity - needs to be done at startup
// but could be done at other times as well
void dispuiInitRotationMatrix(void);

// Set the scene translation, usually called in response to
// middle mouse button
void dispuiTranslation(float xTranslation,float yTranslation);

// Set a new center of rotation, initially 0.0,0.0,0.0 - but can
// be moved.
void dispuiCenterOfRotation(float x,float y,float z);

// Do our best to get back to original Text Screen 
// On SGI, this can involve some low level hardware controls
void dispuiShutDown(void); 

// The object coloring scheme is very simple.  Here are the colors 0 to 9
enum {nBlack,nWhite,nRed,nBlue,nAqua,nGreen,nMagenta,nOrange,nYellow,nGray};

// Your display* code may sometimes need to set the color
// by the above enums, prior to drawing an object.  Use this function
void dispuiColorByPalette(int paletteColor);

// Your display* code may want to output a string like "Orange" for
// user input or other purposes.  Call this function:
const char* dispuiPaletteColorString(int paletteColor);

// If your application needs time to do a calculation, then you
// will want to show an hour glass (or beach ball on Apple)
// Call this function:
void dispuiSetWaitCursors(void);

// Call this function when calculations are complete:
void dispuiPopCursors(void);


// Commands for Drawing Strings.  A "\001" embedded in a string
// causes subsequent characters to be underlined.
// The (void* font) parameter is a GLUT font.
void dispuiDrawStringAtCurrentRasterPosition(void* font,const char* s);
void DrawString(float x,float y,void* font,const char* s);
// Compute width of a string in pixels, based on font and string contents
int StringWidth(void* font, const char* s);

/* The state information for the running movie is stored 
   in the "movie" variable */

/* Grouping it into a structure allows reading and writing it to disk. */
typedef enum { nSingleImage,		/* MONO in old Display */
               nPairOfImages,		/* PAIR in old Display */
               nSGILegacyHardwareStereo,	/* SPLIT Stereo for old SGI machines	 */
               nThreeImages,     //TRI view in older Display
               nQuadStereo       // Quad Buffered Stereo
            } TImageView;		/* TRI in old Display */

typedef struct  {
   // The rotation to be applied to the entire scene.
	float rotationMatrix[16];
   
   // Scaling: 1.0 is the usualy start point but it can be varied with
   // the '['/']' keys and scaling bar contrl
	float scaleFactor;	

   // Clipping panes are parallel to the monitor screen.
   // Typical start values are 1.0 for front, and 2.0 for back.
   // These are multiplied by movieImageSize as needed during display
	float frontClip;		
	float backClip;		

   // Translation of scene selected by middle mouse button
	float xTranslation,yTranslation;
   
   // Current center of scene.  0.0,0.0,0.0 is "normal"
   // But, you can change this with the <o>rigin command at runtime.
	float centerOfRotation[3];

   // Are we showing single image, paris, stereo, or TRI?
   // See above for options
	TImageView imageView;

   /* True if we should eliminate Menu and Scale bar etc for max screen */
	short hideLeftAndBottom;	

   // Count of frames which should be smeared together
	short smear;	
} TMovieWindowControl;

// This global in dispui.c can be directly read by other code.
extern TMovieWindowControl movie;

// Call this function to being keyboard input.  See examples in
// display.c or cujmovie.c
void KeyboardInputStart(const char* _promptString,
								void (* _InputCompleteCallback)(int lastKey,const char* inputString));

// Sometimes, it is helpful to have a default in the input box.
// See display.c for example of how this is used.
void KeyboardInputSetDefault(const char* defaultString);

// You can "fake" keyboard input by passing a sequence to 
// this function.  This allows a great deal of flexibility in handling
// some kinds of events when keyboard string input is active.
void KeyboardStuff(const char* stuffString); 

// return the current movie frame #
int dispuiMovieFrame(void);

// This function tells dispui to redisplay the current frame
// and should be called whenever changes in any component of the
// scene have been made.  Apologies for the lengthy name of this
// often called function.
void MovieAdvanceTimerKeyAndPostRedisplay(void);

// When the static menu botton is clicked, dispui calls the function below
// to alert the client application.
enum {nStaticFirst,nStaticAverage,nStaticNone};

void staticMenuCallback(int whichStaticOption); /* Function called when Static Menu choice selected. */

