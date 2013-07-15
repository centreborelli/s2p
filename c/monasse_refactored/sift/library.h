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

#ifndef _LIBRARY_H_
#define _LIBRARY_H_

void copy(float *u,float *v,int size);	/// v = u

void combine(float *u,float a,float *v,float b,float *w, int size); ///w=a*u+b*v

float* gauss(int sflag,float std,int *size); /// 1D Gauss kernel

void compute_gradient_orientation(float* igray, float *grad, int w, int h);

void sample(float *igray, float *ogray, float factor, int w, int h);
 
void draw_line(float *igray, int a0, int b0, int a1, int b1,
               float value, int width, int height);

#endif
