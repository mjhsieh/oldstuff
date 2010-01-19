/* 
  File of little C output routines - callable from Fortran
  Purpose:  to do concise, precision-limited output,
            without extra spaces - from Fortran
  Scope:    works to either file or screen
  Overview: 
            2 routines to open/close file
	    8 routines for file output (2 int, 2 real, 2 double, 2 string)
	    8 routines for screen output (2 int, 2 real, 2 double, 2 string)
  
File open/close: (uses global variable FILE *fp)

  call fopen('file.name') - open a file of given name for writing "w"
  call fclose - close the file

File write: (uses global variable FILE *fp)

  call fprinti(ii) - prints an integer with default format string '%d '
  call fprintsi('Answer = %d ',ii) - prints an integer with format string
  
  call fprintf(rr) - prints a real*4 with default format string '%g '
  call fprintsf('Answer = %g ',rr) - prints a real*4 with format string
  
  call fprintd(rr) - prints a real*8 with default format string '%g '
  call fprintsd('Answer = %g',rr) - prints a real*8 with format string
  
  call fprintn - prints a newline '\n'
  call fprints('This is the end\n') - prints any string

Screen write: (same as preceeding 8 routines without leading "f")

------------------------------------------------------------------------ */

#include <stdio.h>

FILE *fp;

void fopen_(s,length)
     char *s;
     int length;
{
  int i;
  char t[256];
  FILE *fopen();
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  fp = fopen(t,"w");
}

void fclose_()
{
  fclose(fp);
}


void fprinti_(variable)
     int *variable;
{
  fprintf(fp,"%d ",*variable);
}


void fprintsi_(s,variable,length)
     char *s;
     int *variable;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  fprintf(fp,t,*variable);
}


void fprintf_(variable)
     float *variable;
{
  fprintf(fp,"%g ",*variable);
}


void fprint_sf_(s,variable,length)
     char *s;
     float *variable;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  fprintf(fp,t,*variable);
}


void fprintd_(variable)
     double *variable;
{
  fprintf(fp,"%g ",*variable);
}


void fprint_sd_(s,variable,length)
     char *s;
     double *variable;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  fprintf(fp,t,*variable);
}


void fprintn_()
{
  fprintf(fp,"\n");
}


void fprints_(s,length)
     char *s;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  fprintf(fp,t);
}


void printi_(variable)
     int *variable;
{
  printf("%d ",*variable);
}


void printsi_(s,variable,length)
     char *s;
     int *variable;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  printf(t,*variable);
}


void printf_(variable)
     float *variable;
{
  printf("%g ",*variable);
}


void print_sf_(s,variable,length)
     char *s;
     float *variable;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  printf(t,*variable);
}


void printd_(variable)
     double *variable;
{
  printf("%g ",*variable);
}


void print_sd_(s,variable,length)
     char *s;
     double *variable;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  printf(t,*variable);
}


void printn_()
{
  printf("\n");
}


void prints_(s,length)
     char *s;
     int length;
{
  int i;
  char t[100];
  
  for (i = 0; i < length; i++) t[i] = s[i];
  t[length] = 0;

  printf(t);
}

