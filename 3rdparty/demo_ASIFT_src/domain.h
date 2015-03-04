// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef _DOMAIN_H_
#define _DOMAIN_H_


#include "numerics1.h"
#include "library.h"
#include "splines.h"


/// Compute homography from n points using svd
void compute_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, float **H);


/// Compute homography using svd + Ransac
void compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, float **H);

/// Compute homography by using Lionel Code
void compute_moisan_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float &epsilon, float **H, int &counter, int recursivity);


/// Apply planar homography to image
void compute_planar_homography_bounding_box(int width, int height, float **H, float *x0, float *y0,  int *nwidth, int *nheight);


void apply_planar_homography(float *input, int width, int height, float **H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight);


/// Apply zoom of factor z
void apply_zoom(float *input, float *out, float zoom, int order, int width, int height);

/// Apply general transformation
void apply_general_transformation(float *input, float *transformx, float* transformy, float *out, float bg, int order, int width, int height, int nwidth,int nheight);

#endif

