#include <stdio.h>
#include <stdlib.h>

#include <X11/IntrinsicP.h>
#ifdef sun
#include <X11/ObjectP.h>	/* why don't they just use X from mit!?! */
#include <X11/RectObjP.h>
#endif

#include <X11/StringDefs.h>

#include <X11/Xraw/XawAll.h>




#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef char* user_data[6];

static user_data title_user_data ={
  "login" , "uid", "gid", "name in life", "home directory", "shell"
};

static user_data all_user_data[] = {
#include "data_to_look"
};

static char stuff_font[] =
  "-*-terminal-medium-*-normal-*-14-*-*-*-*-*-iso8859-*";

static char title_font[] =
  "-*-courier-bold-o-*-*-18-*-*-*-*-*-iso8859-*";

static char number_font[] =
  "-*-courier-bold-r-normal-*-18-*-*-*-*-*-iso8859-*";

static char command_font[] =
  "-*-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-*";

static char stuff_color[]  = "#aea292";
static char number_color[] = "#5a7e86";
static char title_color[]  = "#5a7e86";


static void ChangeColor();
static void test_callback();
static void Quit();
static void ChangedColumnWidth();
static void ChangedRowHeight();
static void AddRow();
static int  error_handler();

static Widget app_shell     = NULL;
static Widget main_widget   = NULL;
static Widget paned_widget  = NULL; 
static Widget command_panel = NULL;
static Widget stuff_widget  = NULL;
static Widget title_widget  = NULL;
static Widget frame_widget  = NULL;
static Widget view_widget   = NULL;
static Widget number_widget = NULL;
static Widget logo_widget   = NULL;
static Widget file_widget   = NULL;
static Widget edit_widget   = NULL;
static Widget menu_button_widget   = NULL;
  
/* Command Widget */

static Widget insert_row  = NULL;
static Widget append_row  = NULL;
static Widget prepend_row = NULL;
static Widget delete_row  = NULL;


void main( argc, argv )
    int   argc;
    char* argv[];
{   
  int i = 0;
  int j = 0;
  int column_width;
  
  int row_height;
  int num_all_user_data = sizeof(all_user_data) / sizeof(all_user_data[0]);



  app_shell = XtInitialize ("list_user", "List_user",
			   (XrmOptionDescRec*)NULL, (Cardinal)0,
			   &argc, argv );


  main_widget = 
    XtCreateManagedWidget ("main",
			   panedWidgetClass, app_shell,
			   NULL, (Cardinal)0);
  
  
  paned_widget = 
    XtVaCreateManagedWidget ("menu_paned",
			     panedWidgetClass, main_widget,
			     XtNorientation, XtorientHorizontal,
			     XtNshadowWidth, 2,
  			     XtNshowGrip,    False, 
			     NULL);

/*********************** Command Widget ************************/

  command_panel = 
    XtVaCreateManagedWidget ("command_panel",
			     panedWidgetClass, main_widget,
			     XtNorientation, XtorientHorizontal,
  			     XtNshowGrip,    True,  
			     NULL);

  insert_row = 
    XtVaCreateManagedWidget ("insert_row",
			     commandWidgetClass, command_panel,

			     XtNlabel,    "Insert Row",
			     XtNshowGrip, False,
			     
			     XtVaTypedArg, XtNfont,
			     XtRString,    command_font,
			     strlen(command_font) + 1,
			     
			     NULL);

  append_row = 
    XtVaCreateManagedWidget ("append_row",
			     commandWidgetClass,  command_panel,

			     XtNlabel,    "Append Row",
			     XtNshowGrip, False,
			     
			     XtVaTypedArg, XtNfont,
			     XtRString,    command_font,
			     strlen(command_font) + 1,
			     
			     NULL);
  
  prepend_row = 
    XtVaCreateManagedWidget ("prepend_row",
			     commandWidgetClass, command_panel,

			     XtNlabel,    "Prepend Row",
			     XtNshowGrip, False,
			     
			     XtVaTypedArg, XtNfont,
			     XtRString,    command_font,
			     strlen(command_font) + 1,
			     
			     NULL);

  delete_row = 
    XtVaCreateManagedWidget ("delete_row",
			     commandWidgetClass, command_panel,

			     XtNlabel,    "Delete Row",
			     XtNshowGrip, False,
			     
			     
			     XtVaTypedArg, XtNfont,
			     XtRString,    command_font,
			     strlen(command_font) + 1,
			     
			     NULL);
  
  (void)XtVaCreateManagedWidget ("simple_row",
				 simpleWidgetClass, command_panel,

				 XtNwidth,    1,
				 XtNheight,   1,
				 XtNshowGrip, False,

				 NULL);


/*********************** Menu Widget ************************/

  (void)XtVaCreateManagedWidget ("button_file",
				 menuButtonWidgetClass, paned_widget,

				 XtNlabel ,   "File",
				 XtNshowGrip, False,
				 XtNmenuName, "file",

				 NULL);
  
  (void)XtVaCreateManagedWidget ("button_edit",
				 menuButtonWidgetClass, paned_widget,

				 XtNlabel ,   "Edit",
				 XtNmenuName, "edit",
				 XtNshowGrip, False,

				 NULL);
  
  (void)XtVaCreateManagedWidget ("empty",
				 simpleWidgetClass, paned_widget,

				 XtNshowGrip,    False,
				 XtNwidth,       1,
				 XtNheight,      1,

				 NULL);
  

  file_widget = XtVaCreateWidget ("file", 
				  simpleMenuWidgetClass, main_widget,
				  XtNshowGrip,    False,
				  NULL);

  menu_button_widget = 
    XtVaCreateManagedWidget ("save",
			     smeBSBObjectClass, file_widget,
			     XtNlabel, "Empty yet..",
			     NULL);
  
  menu_button_widget = 
    XtVaCreateManagedWidget ("quit",
			     smeBSBObjectClass, file_widget,
			     XtNlabel, "Quit",
			     NULL);
  
  XtAddCallback (menu_button_widget, XtNcallback, Quit, NULL);
  
  edit_widget = XtVaCreateWidget ("edit", 
				simpleMenuWidgetClass, main_widget,
				XtNshowGrip,    False,
				NULL);

  menu_button_widget = 
    XtVaCreateManagedWidget ("add_user",
			     smeBSBObjectClass, edit_widget,
			     XtNlabel, "Add row to the end..",
			     NULL);
  
  XtAddCallback (menu_button_widget, XtNcallback, AddRow, NULL);
  
  
  (void)XtCreateManagedWidget ("add_column",
				 smeLineObjectClass, edit_widget, NULL, 0);
  
  (void)XtVaCreateManagedWidget ("delete_user",
				 smeBSBObjectClass, edit_widget,
				 XtNlabel, "Empty yet..",
				 NULL);

  
  /*******************************************************************
   *
   *
   *                       ScrolledTable
   *
   *
   *******************************************************************/

  
  paned_widget = 
    XtVaCreateManagedWidget ("table_paned",
			     scrolledTableWidgetClass, main_widget,

			     XtNallowHoriz,  True,
			     XtNshadowWidth, 2,
			     XtNshowGrip,    True,

			     NULL);


  logo_widget = 
    XtCreateManagedWidget ("logo",
			   logoWidgetClass, paned_widget,
			   NULL, 0);
  


  stuff_widget = 
    XtVaCreateManagedWidget("table",
			    tableWidgetClass, paned_widget,
			    XtNliteral,      True,
			    XtNcolumns,      5,
			    XtNrows,         num_all_user_data,
			    XtNjustify,      XtJustifyLeft,

			    XtVaTypedArg, XtNfont,
			    XtRString,    stuff_font,
			    strlen(stuff_font) + 1,
			    
			    XtVaTypedArg, XtNbackground,
			    XtRString, stuff_color,
			    strlen (stuff_color) + 1,
			    
			    NULL);

  XawTableSetColumnJustify (stuff_widget, 0, XtJustifyCenter);
  XawTableSetColumnJustify (stuff_widget, 1, XtJustifyCenter);

  
  title_widget = 
    XtVaCreateManagedWidget ("title",
			     tableWidgetClass, paned_widget,
			     XtNliteral,      False,
			     XtNcolumns,      5,
			     XtNrows,         1,
			     XtNeditable,     False,
			     
			     XtVaTypedArg, XtNfont,
			     XtRString,    title_font,
			     strlen(title_font) + 1,
			     
			     XtVaTypedArg, XtNbackground,
			     XtRString, title_color,
			     strlen(title_color) + 1,
			     
			     NULL);

  XtAddCallback(stuff_widget, XtNchangedColumnWidth, 
		ChangedColumnWidth, title_widget);

  for (j = 0; j < 5; j++) 
  {
    column_width = 2;
    
    for (i = 0; i < num_all_user_data; i++)
      column_width = MAX(column_width, strlen(all_user_data[i][j+1]));
    
    XawTableSetColumnWidth (stuff_widget, j, column_width);
  }
  
  XtVaGetValues (stuff_widget,
		 XtNrowHeight, &row_height,
		 NULL);


  number_widget = 
    XtVaCreateManagedWidget ("number",
			     tableWidgetClass, paned_widget,
			     XtNliteral,      True,
			     XtNcolumns,      2,
			     XtNeditable,     False,
			     XtNrowHeight,    row_height,
			     XtNrows,         num_all_user_data,
			     XtNjustify,      XtJustifyLeft,
			     
			     XtVaTypedArg, XtNfont,
			     XtRString,    number_font,
			     strlen(number_font) + 1,
			     
			     XtVaTypedArg, XtNbackground,
			     XtRString, number_color,
			     strlen(number_color) + 1,

					   NULL);
  XtAddCallback(stuff_widget, XtNchangedRowHeight,
		ChangedRowHeight, number_widget);

  
  XawTableSetColumnWidth (number_widget, 0, 3);
  XawTableSetColumnJustify (number_widget, 0, XtJustifyCenter);

  column_width = 2;
  
  for (i = 0; i < num_all_user_data; i++)
    column_width = MAX(column_width, strlen(all_user_data[i][0]));
  XawTableSetColumnWidth (number_widget, 1, column_width);
  

  
  for (i = 0; i < num_all_user_data; i++)
  {
    char number[10];
    sprintf(number,"%-d\0",i+1);
    XawTableSetLabel (number_widget, i, 0, number);
    XawTableSetLabel (number_widget, i, 1, all_user_data[i][0]);
  }

  for (j = 0; j < 5; j++)
  {
    for (i = 0; i < num_all_user_data; i++)
    {
      XawTableSetLabel (stuff_widget, i, j, all_user_data[i][j+1]);
    }
    XawTableSetLabel (title_widget, 0, j, title_user_data[j+1]);
  }

  XtAddCallback(stuff_widget, XtNwhatCell, ChangeColor, stuff_widget);
  XtAddCallback(stuff_widget, XtNchangedCell, test_callback, stuff_widget);
  

  XtVaSetValues(paned_widget,
		XtNrowWidget,    title_widget,
		XtNcolumnWidget, number_widget,
		XtNstuffWidget,  stuff_widget,
		XtNsignWidget,   logo_widget,
		NULL);
  

  XtVaSetValues(app_shell,
		XtNwidth, 600,
		XtNheight, 600,
		NULL);


#if defined(EBUG) || defined(DEBUG)  
  XSynchronize( XtDisplay(app_shell), True );
  (void)XSetErrorHandler(error_handler);
#endif

  XtRealizeWidget (app_shell);


  
  XtMainLoop ();
  
}

static void test_callback(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data;
    XtPointer call_data;
{
  XawTableCallbackStruct* table_struct = (XawTableCallbackStruct*) call_data;
  Pixel fore, back;

}


/**************************************************************
 *
 *                          Callbacks
 *
 **************************************************************/
/* ARGSUSED */
static void ChangeColor(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data;
    XtPointer call_data;
{
  XawTableCallbackStruct* table_struct = (XawTableCallbackStruct*) call_data;
  Pixel fore, back;
  Pixel table_fore, table_back;


  XtVaGetValues(widget, 
		XtNforeground, &table_fore,
		XtNbackground, &table_back,
		NULL);

  XawTableGetCellColoursByCell(widget, table_struct->old_cell, &fore, &back);

  if (table_fore == fore)
  {
    if (!FetchPixel ((Widget)client_data, "brown", &fore))
      fore = BlackPixelOfScreen(XtScreen((Widget)client_data));


    if (!FetchPixel ((Widget)client_data, "gray", &back))
      back = WhitePixelOfScreen(XtScreen((Widget)client_data));
  


    XawTableSetCellColours ((Widget)client_data,
		     table_struct->row,
		     table_struct->column,
		     fore,
		     back);
  }
  else
  {
    XawTableSetCellDefaultColours ((Widget)client_data,
		     table_struct->row,
		     table_struct->column);
  }
}

/* ARGSUSED */
static void Quit(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data;
    XtPointer call_data;
{
  exit(0);
}

/* ARGSUSED */
static void AddRow(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data;
    XtPointer call_data;
{
  char number[10];
  int rows;
  Widget scroll = NULL;
  float shown;

  XawTableAppendRow(number_widget);
  XawTableAppendRow(stuff_widget);
  XtVaGetValues(number_widget, XtNrows, &rows, NULL);

  sprintf(number,"%-d\0",rows-1);
  XawTableSetLabel (number_widget, rows-1, 0, number);
  scroll = XtNameToWidget(paned_widget, "vertical");

  if (scroll)
  {
    XtVaGetValues(scroll, XtNshown, &shown, NULL);
    XawScrollbarSetThumb(scroll, 
#if NeedWidePrototypes			 
			 (float)(1.0-shown), (float)shown);
#else    
			 (double)(1.0-shown), (double)shown);
#endif
  }
}


/* ARGSUSED */
static void ChangedColumnWidth (widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data; /* unused */
    XtPointer call_data;
{
  Widget title = (Widget)client_data;
  XawTableCallbackStruct* table_struct = (XawTableCallbackStruct*) call_data;
  int column_width;

  column_width = XawTableGetColumnPixelWidth (stuff_widget,
					      table_struct->column);

  XawTableSetColumnWidth (title, table_struct->column, column_width);
}


/* ARGSUSED */
static void ChangedRowHeight(widget, client_data, call_data)
    Widget    widget;
    XtPointer client_data; /* unused */
    XtPointer call_data;
{
  Widget number = (Widget)client_data;
  int row_height;

  XtVaGetValues (widget, XtNrowHeight, &row_height, NULL);
  XtVaSetValues (number, XtNrowHeight, row_height, NULL);
}



static int error_handler(dpy, err)
          Display *dpy;
          XErrorEvent *err;
{
  char *p = NULL;
  /* Sigmantation fault is needed to get coredump */
  *p = ' ';
  return 1;
}
