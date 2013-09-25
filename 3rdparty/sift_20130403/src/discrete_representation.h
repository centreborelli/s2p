/*
IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan 
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20130318 (March 18, 2013)

== Patent Warning and Licence =================================================

The SIFT method is pattented 

    [3] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89
  
 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

*/
/**
 * @file discrete_representation.h
 * @brief simple image transformations 
 *
 * This is a front-end to libpng, with routines to:
 *      @li Separable discrete convolutions
 *      @li Gaussian blur via discrete convolution with truncated kernel
 *      @li 2D Gradient
 *      @li Subsampling by integer factor
 *      @li bilinear interpolation
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#ifndef _DISCRETE_REPRESENTATION_H_
#define _DISCRETE_REPRESENTATION_H_

// TEMP for printf
#include <stdio.h>



#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define MAX(i,j) ( (i)<(j) ? (j):(i) )

/** @brief Discrete Gaussian convolution on image
 * 
 *   Applies in each direction monodimensional sampled Gaussian kernel and
 *   truncated at radius 4\sigma.
 * 
 * \param in input image of size  width X height.
 * \param out output image of size  width X height.
 * \param sigma standard deviation of the Gaussian kernel
 */
void add_gaussian_blur(double* in, double* out, int width, int height, double sigma);



/** @brief Image Gradient
 * 
 *  Computes the gradient via the centered finite difference scheme
 *  [-1/2,0,+1/2]
 * 
 * \param im input image of size  width X height.
 * \param im_x gradient component along x (|  top-bottom)  ( width X height samples).
 * \param im_y gradient component along y (-> left-right)  ( width X height samples).
 * \param sigma standard deviation of the Gaussian kernel 
 * 
 */
void compute_gradient(double* im, double* im_x, double* im_y, int width, int height);




/** @brief Image subsampling by a factor 2
 * 
 * \param in [i,j] , (0 <= i <= height-1) (0 <= j <= width-1)  . 
 * \param out [i,j]= in[2*i,2*j] ,  (0 <= i <= int(height/2)-1) (0 <= j <= int(width/2)-1)
 * 
 */
void subsample_by2(double* in, double* out, int in_width, int in_height);



/** @brief Image subsampling by an integer factor k
 * 
 *  Computes the gradient via the centered finite difference scheme
 *  [-1/2,0,+1/2]
 * 
 * \param in [i,j] , (0 <= i <= height-1) (0 <= j <= width-1)  . 
 * \param out [i,j]= in[k*i,k*j] ,  (0 <= i <= int(height/k)-1) (0 <= j <= int(width/k)-1)
 * 
 */
void subsample_by_intfactor(double* in, double* out, int in_width, int in_height, int factor);







void oversample_by2_bilin(double* in, double* out, int in_width, int in_height);

void oversample_bilin(double* in , int in_width,  int in_height,
                      double* out, int out_width, int out_height,
                      double delta_out);




/** @brief Separable convolution
 * 
 * 
 * 
 * 
 */
void convolve(double* in, double* out,
              int width, int height,     /* input image and output image have same image */
              double* xker, int r_xker,
              double* yker, int r_yker);


#endif
