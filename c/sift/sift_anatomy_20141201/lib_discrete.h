/*
IPOL SIFT
Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan 
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20140911 (September 11th, 2014)

This C ANSI source code is related to the IPOL publication

    [1] "Anatomy of the SIFT Method." 
        I. Rey Otero  and  M. Delbracio
        Image Processing Online, 2013.
        http://www.ipol.im/pub/algo/rd_anatomy_sift/

An IPOL demo is available at
        http://www.ipol.im/pub/demo/rd_anatomy_sift/





== Patent Warning and License =================================================

The SIFT method is patented 

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
 * @file lib_discrete.h
 * @brief simple image transformations 
 *
 * This is a front-end to libpng, with routines to:
 *      @li Separable discrete convolutions
 *      @li Gaussian blur via discrete convolution with truncated kernel
 *      @li 2D Gradient
 *      @li Subsampling by integer factor
 *      @li bilinear interpolation
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#ifndef _LIB_DISCRETE_H_
#define _LIB_DISCRETE_H_





/** @brief Discrete Gaussian convolution on image
 * 
 *   Applies in each direction monodimensional sampled Gaussian kernel and
 *   truncated at radius 4\sigma.
 * 
 * \param in input image of size  w X h.
 * \param out output image of size  w X h.
 * \param sigma standard deviation of the Gaussian kernel
 */
void sift_add_gaussian_blur(const float* in, float* out, int w, int h, float sigma);



/** @brief Image Gradient
 * 
 *  Computes the gradient via the centered finite difference scheme
 *  [-1/2,0,+1/2]
 * 
 * \param im input image of size  w X h.
 * \param im_x gradient component along x (|  top-bottom)  ( w X h samples).
 * \param im_y gradient component along y (-> left-right)  ( w X h samples).
 * \param sigma standard deviation of the Gaussian kernel 
 * 
 */
void sift_compute_gradient(const float* im, float* im_x, float* im_y, int w, int h);




/** @brief Image subsampling by a factor 2
 * 
 * \param in [i,j] , (0 <= i <= h-1) (0 <= j <= w-1)  . 
 * \param out [i,j]= in[2*i,2*j] ,  (0 <= i <= int(h/2)-1) (0 <= j <= int(w/2)-1)
 * 
 */
void sift_subsample_by2(const float* in, float* out, int wi, int hi);



void oversample_by2_bilin(const float* in, float* out, int wi, int hi);

/** @brief Interpolate the image with a bilinear model
 * 
 *  the inter-pixel distance in the output image is delta_out
 * 
 *  in  : input digital image with (wi X hi) samples.
 *  out : output digital image with (wo X ho) samples,
 *        with wo  = \lfloor wi  / delta_out \rfloor
 *         and ho = \lfloor hi / delta_out \rfloor
 * 
 * 
 */
void sift_oversample_bilin(const float* in , int wi,  int hi,
                      float* out, int wo, int ho,
                      float delta_out);



#endif  // _LIB_DISCRETE_H_
