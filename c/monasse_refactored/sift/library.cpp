// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// David Lowe  "Method and apparatus for identifying scale invariant 
// features in an image and use of same for locating an object in an 
// image",  U.S. Patent 6,711,293.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.

#define _USE_MATH_DEFINES // For Windows
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "library.h"

void copy(float *u,float *v,int size) {
    for(int i=0; i<size ;i++)
        v[i]=u[i];
}

void combine(float *u,float a,float *v,float b, float *w,  int size)  {
    for(int i=0;i<size ;i++)
        w[i]= a*u[i] + b*v[i];
}

static int normalize(float *u,int size) {
  float s = 0.0;  
  for(int i=0;i<size;i++) s+=u[i];
  if (s != 0.0) {	for(int i=0;i<size;i++) u[i]/=s;  return 1;}  
  else return 0;
}
	
float* gauss(int sflag,float std,int *size) {
   float *u,prec = 4.0,shift;
   double v;
   int n,i;

   if (sflag) n=*size;
   else
	n = 1+2*(int)ceil((double)std*sqrt(prec*2.*log(10.)));   
   
   u =new float[n];

   if (n==1) 
    u[0]=1.0;
   else{

      shift = 0.5f*(float)(n-1);

      for (i=(n+1)/2;i--;) {
         v = ((double)i - (double) shift)/(double)std;
         u[i] = u[n-1-i] = (float) exp(-0.5*v*v); 
      }
   }	

   if (normalize(u,n)) {
       *size=n;
   } else {
	printf("ERROR: _gauss: _normalize: normalization equals zero.\n");
   }
   return u;
}

void compute_gradient_orientation(float* igray,float *grad, int w, int h) {
    float* out=grad;
    for(int r=0; r<h; r++)
      for(int c=0; c<w; c++) {
        float xgrad, ygrad;
        if (c == 0)
          xgrad = 2.0f * (igray[r*w+c+1] - igray[r*w+c]);
        else if (c == w-1)
          xgrad = 2.0f * (igray[r*w+c] - igray[r*w+c-1]);
        else
          xgrad = igray[r*w+c+1] - igray[r*w+c-1];
        if (r == 0)
          ygrad = 2.0f * (igray[r*w+c] - igray[(r+1)*w+c]);
        else if (r == h-1)
          ygrad = 2.0f * (igray[(r-1)*w+c] - igray[r*w+c]);
        else
          ygrad = igray[(r-1)*w+c] - igray[(r+1)*w+c];

        *out++ = (float)sqrt((double)(xgrad * xgrad + ygrad * ygrad));
        *out++ = (float)atan2 (-(double)ygrad,(double)xgrad);
      }
}

void sample(float *igray,float *ogray, float factor, int width, int height) {
	int swidth = (int)((float) width / factor);
	int sheight = (int)((float) height / factor);

	for(int j=0; j < sheight; j++)
        for(int i=0; i < swidth; i++)
            ogray[j*swidth+i] = igray[int((float)j*factor)*width +
                                      int((float)i*factor)];
}

void draw_line(float *igray, int a0, int b0, int a1, int b1, float value,
               int w, int h) {
  int sx,sy,dx,dy,x,y,z,l;

  if (a0 < 0) a0=0; 
  else if (a0>=w) a0=w-1;
   
  if (a1<0)  a1=0; 
  else  if (a1>=w)   a1=w-1;
	  
  if (b0<0) b0=0; 
  else if (b0>=h) b0=h-1;
   
  if (b1<0) 	b1=0; 
  else if (b1>=h) b1=h-1; 

  if (a0<a1) { sx = 1; dx = a1-a0; } else { sx = -1; dx = a0-a1; }
  if (b0<b1) { sy = 1; dy = b1-b0; } else { sy = -1; dy = b0-b1; }
  x=0; y=0;
  
  if (dx>=dy) 
    {
      z = (-dx) / 2;
      while (abs(x) <= dx) 
	{

	  l =  (y+b0)*w+x+a0; 
	
	  igray[l] = value;
	  
	  x+=sx;
	  z+=dy;
	  if (z>0) { y+=sy; z-=dx; }

	} 

    }
  else 
    {
      z = (-dy) / 2;
      while (abs(y) <= dy) {

	l = (y+b0)*w+x+a0;
  	igray[l] = value;
 
	y+=sy;
	z+=dx;
	if (z>0) { x+=sx; z-=dy; }
      }
    }
}
