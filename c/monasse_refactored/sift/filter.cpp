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

#define _USE_MATH_DEFINES
#include <cmath>
#include "filter.h"
#include "library.h"
#include <cfloat>
#include <cstdlib>
#include <algorithm>

static void buffer_convolution(float *buffer,float *kernel,int size,int ksize)
{

    for (int i = 0; i < size; i++) {

      float sum = 0.0;
      float *bp = &buffer[i];
      float *kp = &kernel[0];
      

      /* Loop unrolling: do 5 multiplications at a time. */
      int k=0;
      for(;k + 4 < ksize;  bp += 5, kp += 5, k += 5) 
	        sum += bp[0] * kp[0] +  bp[1] * kp[1] + bp[2] * kp[2] +
        	       bp[3] * kp[3] +  bp[4] * kp[4];

      /* Do multiplications at a time on remaining items. */
      for(; k < ksize; bp++ , kp++, k++)  sum += *bp * (*kp);

      buffer[i] = sum;
    }
}

/* Convolve image with the 1-D kernel vector along image rows.  This
   is designed to be as efficient as possible.  
*/
static void horizontal_convolution(float *u, float *v, int width, int height,
                                   float *kernel, int ksize, int boundary)
{

    int halfsize = ksize / 2;
    int buffersize = width + ksize;
    float *buffer = new float[buffersize];

    for (int r = 0; r < height; r++) {

	/// symmetry
	int l = r*width;
	if (boundary == 1)
        	for (int i = 0; i < halfsize; i++)
            		buffer[i] = u[l + halfsize - 1 - i ];
	else
		for (int i = 0; i < halfsize; i++)
            		buffer[i] = 0.0;


        for (int i = 0; i < width; i++)
            buffer[halfsize + i] = u[l + i];


	if (boundary == 1)
	        for (int i = 0; i <  halfsize; i++)
        	    buffer[i + width + halfsize] = u[l + width - 1 - i];
	else 
		for (int i = 0; i <  halfsize; i++)
        	    buffer[i + width + halfsize] = 0.0;

        buffer_convolution(buffer, kernel, width, ksize);
        for (int c = 0; c < width; c++)
          v[r*width+c] = buffer[c];
    }
    delete [] buffer;
}

static void vertical_convolution(float *u, float *v, int width, int height,
                                 float *kernel,int ksize, int boundary)
{
    int halfsize = ksize / 2;
    int buffersize = height + ksize;
    float *buffer = new float[buffersize];

    for (int c = 0; c < width; c++) {

	if (boundary == 1)
		for (int i = 0; i < halfsize; i++)
				buffer[i] = u[(halfsize-i-1)*width + c];
	else 
		for (int i = 0; i < halfsize; i++)
				buffer[i] = 0.0f;

      for (int i = 0; i < height; i++)
		buffer[halfsize + i] = u[i*width + c];

	if (boundary == 1)
		for (int i = 0; i < halfsize; i++)
				buffer[halfsize + height + i] = u[(height - i - 1)*width+c];
	else
		for (int i = 0; i < halfsize; i++)
			buffer[halfsize + height + i] = 0.0f;

      buffer_convolution(buffer, kernel, height, ksize);

      for (int r = 0; r < height; r++)
          v[r*width+c] = buffer[r];

    }
    delete [] buffer;
}

/* Convolution with a kernel */
/* No padding applied to the image */
void convol(float *u,float *v,int width,int height,float *kernel,int kwidth,int kheight)
{
 
  int K2 = kwidth / 2;
  int L2 = kheight / 2;

  for(int y=0 ; y < height; y++) 
    for (int x=0 ; x < width; x++) {

	float S = 0.0;
	
	for (int l = -L2; l <= L2; l++) 
		for (int k = -K2 ; k<= K2; k++)
		{ 
			int px=x+k;
			int py=y+l;

			if (px>=0 && px < width && py>=0 && py<height)
			  S += u[width*py + px] * kernel[kwidth*(l+L2) + k+K2];
		}

      v[y*width+x] = (float) S;

    }
}

void gaussian_convolution(float *u, float *v, int width, int height, float sigma)
{
	int ksize;	
	float * kernel;

    ksize = (int)(2.0 * 4.0 * sigma + 1.0);
    if(ksize % 2 == 0)
        ++ksize;
	kernel = gauss(1,sigma,&ksize);

	int boundary = 1;

	copy(u,v,width*height);
	horizontal_convolution(v, v, width, height, kernel, ksize, boundary);
    vertical_convolution(v, v, width, height,  kernel,  ksize, boundary);
    delete [] kernel;
}
