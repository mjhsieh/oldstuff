//
// avg: a small code to find the average of an array 
// Author: Mengjuei Hsieh, University of California Irvine
// Original author: Prof. Ray Luo
//
#include <stdio.h>
#include <math.h>

#define DIMEN 1000000

int main(int argc, char **argv)
{
  FILE *ifp;
  
  int i,j,tmpdim;
  float tmpread,sum,avg,dev;
  float *buffer;
  tmpread=0;
  tmpdim=DIMEN;
  buffer=(float *)malloc(sizeof(float)*DIMEN);
  if(buffer==NULL){
    exit(1);//Insufficient memory
  }
  
  if(argc!=2){
    printf("Mengjuei's averaging tool\n");
    printf("Usage: %s <file>\n",argv[0]);
    exit(0);
  }
  
  if(!(ifp=fopen(argv[1],"r"))){
    printf("Error: Unable to open file %s\n",argv[1]);
    exit(1);
  }
  
  sum=0.0;
  dev=0.0;
  
  i=0;
  while(fscanf(ifp,"%f\n",&tmpread)!=EOF){
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
  
  for(j=0;j<=i-1;j++){
    dev=dev+(avg-buffer[j])*(avg-buffer[j]);
  }
  
  dev=sqrt(dev/(i-1));
  printf("points: %5d average: %8.4f std_dev: %8.4f\n",i,avg,dev);

  free(buffer);
}
