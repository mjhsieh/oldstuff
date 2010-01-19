/* fullscreen_stereo.c  --  GLUT support for full screen stereo mode  on SGI
   workstations. */

/* 24-Oct-95 Mike Blackwell  mkb@cs.cmu.edu */

#include <stdlib.h>
#include <stdio.h>
#include <get.h>
#include <X11/Xlib.h>
#include <X11/extensions/XSGIvcstr.h>
#include <X11/extensions/XSGIvc.h>
#include <X11/extensions/SGIStereo.h>

/* We need to access some GLUT internal variables - this include file  is
   found in the GLUT source code distribution. */

/* XXX I do not normally encourage programs to use GLUT internals.  Programs
   that do (like this one) are inherently unportable GLUT programs.  In the
   case of SGI's low-end stereo there was enough demand to warrant supplying
   an example, and the low-end stereo is not clean enough to be supported
   directly in GLUT. -mjk */

#include "glutint.h"

#include "fullscreen_stereo.h"

/* XXX Video display modes for stereo are selected by running
   /usr/gfx/setmon; in IRIX 6.2 and later releases, the XSGIvc API supplies
   the functionality of setmon and more. */

static char current_display_string[200];
void get_current_display_string(void)
{
	// Let's shell outto the system - and read the output of this
	// program to figure out what kind of display to switch back to!
	// This is much sloppier than calling the API - but I just don't 
	// see how else to do it... API is poorly documented.  Sorry.

	FILE* f = popen("/usr/gfx/gfxinfo","r");
	int display_info_found = 0;
	char buf[200];

   if (f != 0)
		{
   	while ((! display_info_found) && fgets(buf,sizeof(buf)-1,f))
			{
			char* s;
			buf[sizeof(buf)-1] = 0;
			s = strstr(buf,"Video Output:");
			if (s != 0)
				{
				char *leftParen = strstr(s," (");
				if (leftParen != 0)
					{
					char *displayString = leftParen+2;
					char *rightParen = strchr(displayString,')');
					if (rightParen != 0)
						{
						*rightParen = 0;
						if (strlen(displayString) < sizeof(current_display_string))
							{
							strcpy(current_display_string,displayString);
							printf("Current Display: %s\n",current_display_string);
							display_info_found = 1;
							}
						}
					}
				}
			}	
		pclose(f);
		}

	if (! display_info_found)
		strcpy(current_display_string,"60hz"); // Sucks - but works
}
   

void
start_fullscreen_stereo(void)
{
  int event, error;

  if (!XSGIStereoQueryExtension(__glutDisplay, &event, &error)) {
    fprintf(stderr, "Stereo not supported on this display!\n");
    exit(0);
  }
  if (XSGIQueryStereoMode(__glutDisplay, __glutCurrentWindow->win) < 0) {
    fprintf(stderr, "Stereo not supported on this window!\n");
    exit(0);
  }

  get_current_display_string(); // See above

 #ifdef BAD_IDEA 
  if (! XSGIvcQueryVideoScreenInfo(__glutDisplay,__glutCurrentWindow->win,&sinfo)) {
    fprintf(stderr, "Unable to get VideoScreenInfo\n");
    exit(0);
  }
#endif
  
#ifdef NOT_WORKING
  if (XSGISetStereoMode(__glutDisplay,__glutCurrentWindow->win,
		492,532,STEREO_BOTTOM) < 0)
	{
    fprintf(stderr, "SetStereoMode not supported on this window!\n");
    exit(0);
	}
#endif

//	if (system("/usr/gfx/setmon -n 1280x1024_100s") != 0)
   if (system("/usr/gfx/setmon -n STR_TOP") != 0) 
		{
		fprintf(stderr,"setmon to STR_TOP failed\n");
		
	    stop_fullscreen_stereo();
    	exit(0);
		}
}

void
stop_fullscreen_stereo(void)
{
	char command[400];
	sprintf(command,"/usr/gfx/setmon -n %s",current_display_string);
  if (system(command) != 0) 
		{
    	fprintf(stderr, "setmon to %s attempt failed!\n",current_display_string);
		// Try some others in desparation
  	 	system("/usr/gfx/setmon -n 60hz");
  	 	system("/usr/gfx/setmon -n 72hz");
		exit(0);
		}
}

void
stereo_left_buffer(void)
{
//   glViewport(0,(1024+YSTEREO)/2,XMAXSCREEN,(1024-YSTEREO)/2);
  XSGISetStereoBuffer(__glutDisplay, __glutCurrentWindow->win, STEREO_BUFFER_LEFT);
  XSync(__glutDisplay, False);
  // glViewport(0, 0, XMAXSCREEN, YSTEREO);
}

void
stereo_right_buffer(void)
{
  XSGISetStereoBuffer(__glutDisplay, __glutCurrentWindow->win, STEREO_BUFFER_RIGHT);
  XSync(__glutDisplay, False);
  // glViewport(0, 0, XMAXSCREEN, (1024-YSTEREO)/2);
}
