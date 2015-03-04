// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#ifndef _SPLINES_H_
#define _SPLINES_H_

#include "numerics1.h"
#include "library.h"
#include "string.h"
#include <vector>
using namespace std;


//float v(float *in,int x,int y,float bg, int width, int height);

// Guoshen Yu, 2010.09.21, Windows version
float v(vector<float>& in,int x,int y,float bg, int width, int height);
//float v(float *in, int x,int y,float bg, int width, int height);

void keys(float *c,float t,float a);
void spline3(float *c,float t);
void init_splinen(float *a,int n);
void splinen(float *c,float t,float *a,int n);

//void finvspline(float *in,int order,float *out, int width, int height);

// Guoshen Yu, 2010.09.22, Windows versions
void finvspline(vector<float>& in,int order,vector<float>& out, int width, int height);
// void finvspline(float *in,int order,float *out, int width, int height);


#endif

