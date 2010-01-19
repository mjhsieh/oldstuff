/*  Author: Ray Luo, Phd., UC Irvine
 *          Mengjuei Hsieh, UC Irvine
 *  Licensing: GPLv2
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#define DIMEN 1000000

/* small code to find the running average of an array */
int main(int argc, char **argv)
{
  FILE  *ifp;
  int   i, j;
  int   count, window;
  float *buffer,sum,avg,dev;
  float ravg;
  float tmpread;
  int   tmpdim;

  tmpdim=DIMEN;
  buffer=(float *)malloc(sizeof(float)*DIMEN);
  if(buffer==NULL){
    exit(1);//Insufficient memory
  }

  if(argc!=3) {
    printf("Usage: %s <window> <file>\n",argv[0]);
    exit(0);
  }

  sscanf(argv[1], "%d", &window);
  if (window <= 0 ) {
    printf("Error: nonpositive window %d\n", window);
    exit(0);
  }
  
  if(!(ifp=fopen(argv[2],"r"))) {
    printf("Error: Unable to open file %s\n",argv[2]);
    exit(0);
  }

  sum=0.0;
  dev=0.0;
  
  i=0;
  while(fscanf(ifp,"%f\n", &tmpread)!=EOF) {
    sum=sum+tmpread;
    if(tmpdim <= i) {
      tmpdim+=1000;
      buffer = (float *)realloc(buffer, sizeof(float)*tmpdim);
      if(buffer==NULL){
        exit(1);//Insufficient memory
      }
    }
    buffer[i]=tmpread;
    i++;
  }

  avg=sum/i;
  
  for(j=0;j<=i-1;j++) {
    dev=dev+(avg-buffer[j])*(avg-buffer[j]);
  }
  
  dev=sqrt(dev/(i-1));
  //printf("points %5d, average %8.4f, std dev %8.4f\n",i,avg,dev);
  
  sum = 0.0;
  count = 0;
  for(j=0;j<=i-1;j++) {
    count++;
    sum = sum + buffer[j];
    if ( count == window ) {
      ravg = sum/count;
      printf("%f\n", ravg);
      count--;
      sum = sum - buffer[j-(window-1)];
    }
  }
}
