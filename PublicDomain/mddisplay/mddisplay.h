// display.h - functions and variables in display.c which are
// called from and read from dispui.c

// dispinfo

#ifndef DISPUI_TESTMODE
#define nMaxMeasurements 15

typedef struct {
	int statictype;
	int nmeas, measarry[nMaxMeasurements][4];  /* measurements */
	float hbcut;
} dispinfo;
extern dispinfo view;
#endif


// dispui.c will prepare a perspective box for the "molecule" and
// call this function in display.c to render the molecule.
void displayRenderFrameInPerspectiveBox(int frame,int smearFlag); // Called from dispui.c to draw a complete "frame"

// Called from dispui to render any final text - in this case we display "MEasurement" info
void displayRenderFrameInOrthoBox(int frame,int movieWindowWidth,int movieWindowHeight);

// What are min and max frame numbers to ask for!
int displayMovieMinFrame(void);
int displayMovieMaxFrame(void);

// Called frmo dispui to aid in selection by the mouse!
void displayRenderFrameInSelectMode(int frame);
void displayProcessHits(GLint hits,const GLuint* selectBuf);
int  displayProcessRawMovieClick(int movieX0,int movieY0,int movieWindowWidth,int movieWindowHeight);


void displayRenderClipToolFrame(int frame); // Called from dispui.c to draw a complete "frame" in clip tool
#ifndef DISPUI_TESTMODE
void displayRenderRamaFrame(int frame); // Called by UI code when a new Rama Frame is needed.
void displaySelectRamaResidue(int frame,int mouseX,int mouseY);
#endif
void displayInputNewFilter(void); // Called by UI when user clicks on Filter

void displayKeyboardFunc(int key);

// Called when user interface detects a press on the palette or other stuff...
void displayNewColorRequested(int newColor);
void displayRemoveRequested(void);
void displayCenterOfRotationRequested(void);
void displayLabelRequested(void);
void displayUNLabelRequested(void);

#ifndef DISPUI_TESTMODE
void displayDistanceMeasurementRequested(void);
void displayAngleMeasurementRequested(void);
void displayDihedralMeasurementRequested(void);


void displayHBondRequested(void);
void displayPDBRequested(void);
#endif
void displayQuitRequested(void);
void displayRenderHelpBox(void);
