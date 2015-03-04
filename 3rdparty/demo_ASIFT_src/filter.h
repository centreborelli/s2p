// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef _FILTER_H_
#define _FILTER_H_


#include "library.h"


float * directional_gauss_filter(float xsigma, float ysigma, float angle, int *kwidth, int *kheight);

	
void median(float *u,float *v, float radius, int niter, int width,int height);
void remove_outliers(float *igray,float *ogray,int width, int height);

/// Convolution with a  separable kernel,  boundary condition: 0=zero, 1=symmetry
void separable_convolution(float *u, float *v, int width, int height, float *xkernel, int xsize, float *ykernel, int ysize,int boundary);     

void buffer_convolution(float *buffer,float *kernel,int size,int ksize);
void horizontal_convolution(float *u, float *v, int width, int height, float *kernel, int ksize, int boundary);
void vertical_convolution(float *u, float *v, int width, int height, float *kernel,int ksize, int boundary);

void fast_separable_convolution(float *u, float *v, int width, int height,float * xkernel, int xsize,float *ykernel,int ysize,int boundary);

/// Can be called with u=v
void gaussian_convolution(float *u, float *v, int width, int height, float sigma);
void gaussian_convolution(float *u, float *v, int width, int height, float sigma, int ksize);

void convol(float *u, float *v, int width, int height, float *kernel, int kwidth, int kheight);    /// Convolution with a kernel,  No padding applied to the image

void heat(float *u, float *v, float step, int niter, float sigma, int width, int height);


#endif // _FILTER_H_

