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
 * @file lib_discrete.c
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




#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "lib_discrete.h"
#include "lib_util.h"







/** @brief Compute image gradient via symmetric finite difference schemes
 *
 *   image extension :  symmetrization at x=-1/2 
 * 
 */
void sift_compute_gradient(const float* im, float* im_x, float* im_y, int w, int h){

    const float* p_in;
    const float* p_in_p;
    const float* p_in_m;
    float* p_out;
    
    /** Computing im_y      */
    /* pixels in the center */
    p_in   = &im[1];
    p_out  = &im_y[1];
    p_in_p = p_in+1;
    p_in_m = p_in-1;
    for(int i=1;i<h*w-1;i++){   /* produces false value on borders but it's ok */
        *p_out = (*p_in_p-*p_in_m)*0.5;
        p_out++;
        p_in_m++;
        p_in_p++;
    }
    /*    pixels on borders - symmetrization at y=-1/2  */
    for(int i=0;i<h;i++){
        im_y[i*w]           = im[i*w+1]- im[i*w];
        im_y[i*w+(w-1)] = im[i*w+(w-1)] - im[i*w+(w-2)];
    }


    /** Computing im_x      */
    /* pixels in the center */
    p_in   = &im[w];
    p_out  = &im_x[w];
    p_in_p = p_in+w;
    p_in_m = p_in-w;
    for(int i=w;i<h*w-w;i++){   /* produces false value on borders */
        *p_out = (*p_in_p-*p_in_m)*0.5;
        p_out++;
        p_in_m++;
        p_in_p++;
    }
    /*    pixels on borders - symmetrization at x=-1/2  */
    for(int i=0;i<w;i++){
        im_x[i]                   = im[i+w]-im[i];
        im_x[(h-1)*w+i]  = im[(h-1)*w+i] - im[(h-1-1)*w+i];
    }
}




/**
 * @brief Builds mono-dimensional sampled Gaussian kernel.
 *
 *        Only considers the right side of the kernel 
 *        (including the center
 *
 *        Returns a (rad+1) long vector
 *
 *
 */
static float* malloc_gaussian_symm_kernel(float sigma, int rad)
{
    assert(sigma>=0);
    float* gker = xmalloc((rad+1)*sizeof(float));
    gker[0] = 1.;
    if(sigma>0){
        float sum = gker[0];
        for(int i = 1; i <= rad; i++){
            gker[i] = exp(-0.5*(float)i*(float)i/sigma/sigma);
            sum += 2*gker[i];
        }
        for(int i = 0; i <= rad; i++)
            gker[i] /= sum;
    }else{
       for(int i = 1; i <= rad; i++){
           gker[i] = 0.0;
       }
    }
    return gker;
}


static void convolve_symm(const float* in, float* out, int w, int h,
                          const float* xker, int r_xker,
                          const float* yker, int r_yker);


void sift_add_gaussian_blur(const float* in, float* out, int w, int h, float sigma){
    int r_gker = (int)ceil(4*sigma);  // ceil(4*sigma)
    float* gker = malloc_gaussian_symm_kernel(sigma, r_gker);
    convolve_symm(in,out,w,h,gker,r_gker,gker,r_gker);
    xfree(gker);
}

/** @brief sub sampling by factor 2, keeping sample (0,0) */
void sift_subsample_by2(const float* in, float* out, int wi, int hi){
    int wo = wi/2;
    int ho = hi/2;
    int i, j, i_p, j_p;
    for(i=0;i<ho;i++){
        i_p = 2*i;
        for(j=0;j<wo;j++){
            j_p=2*j;
            out[i*wo+j] = in[i_p*wi+j_p];
        }
    }
}




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
                          float delta_out)
{
    assert(delta_out<=1);

    for(int i = 0; i < ho; i++){
        for(int j = 0; j < wo; j++){

            float x = i*delta_out;
            float y= j*delta_out;
            int im = (int)x;
            int jm = (int)y;
            int ip = im +1;
            int jp = jm +1;

            //image extension by symmetrization
            if(ip >= hi){ip = 2*hi-1-ip;}
            if(im >= hi){im = 2*hi-1-im;}
            if(jp >= wi){jp = 2*wi-1-jp;}
            if(jm >= wi){jm = 2*wi-1-jm;}

            const float fractional_x = x - floor(x);
            const float fractional_y = y - floor(y);
            out[i*wo+j] = fractional_x  * ( fractional_y  * in[ip*wi+jp]
                                   + (1 - fractional_y) * in[ip*wi+jm] )
                     + (1-fractional_x) * ( fractional_y  * in[im*wi+jp]
                                   + (1 - fractional_y) * in[im*wi+jm] );
        }
    }
}



// Returns the corresponding coordinate 0<=i<l under symmetrized signal
// extension 
static inline int symmetrized_coordinates(int i, int l)
{
    int ll = 2*l;
    i = (i+ll)%(ll);
    if(i>l-1){i = ll-1-i;}
    return i;
}



/** @brief Apply a convolution with a separable kernel
 * and signal extension by symmetrization
 *
 * @param in             Input image of w X h samples
 * @param out            Output image (same dimension)
 *
 * @param xker           Kernel applied along direction x
 *                          radius = r_xker
 *                          w = 2*r_xker+1
 *
 * @param yker           Kernel applied along direction y
 *                          radius = r_yker
 *                          w = 2*r_yker+1
 *
 *  Compute:
 *
 *   out[i,j] = \sum_{-r_xker \leq k \leq +r_xker} 
 *                    xker[k]
 *                    . \sum_{-r_xker \leq l \leq +r_xker}
 *                              yker[l].in[i-k,j-l]
 *
 *  Border condition:  symmetrization at border
 *                     (at x = -1./2 and x = w-1./2)
 *
 */
static void convolve_symm(const float* in, float* out, int w, int h,
                          const float* xker, int r_xker,
                          const float* yker, int r_yker)
{
    float* im_tmp = xmalloc(w*h*sizeof(float));
    /* convolution along x coordinates */
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++){
            float sum = in[i*w+j] * xker[0] ;
            for(int k = 1; k <= r_xker; k++){
                int i_p_left = symmetrized_coordinates(i-k, h);
                int i_p_right = symmetrized_coordinates(i+k, h);
                sum += xker[k]*(in[i_p_left*w+j] + in[i_p_right*w+j]);
            }
            im_tmp[i*w+j] = sum;
        }
    }
    /* convolution along y coordinates */
    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++){
            float sum = im_tmp[i*w+j] * xker[0];
            for(int k = 1; k <= r_yker; k++){
                int j_p_left = symmetrized_coordinates(j-k, w);
                int j_p_right = symmetrized_coordinates(j+k, w);
                sum += yker[k]*(im_tmp[i*w+j_p_left] + im_tmp[i*w+j_p_right]);
            }
            out[i*w+j] = sum;
        }
    }
    xfree(im_tmp);
}
