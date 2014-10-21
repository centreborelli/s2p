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

    [2] "Method and apparatus for identifying scale invariant features
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
 * @file sift_description.c
 * @brief Computation the SIFT feature vector
 *
 * @li Attribution of a principal orientation
 * @li Computation of the SIFT feature vector
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */





#include <math.h>
#include "lib_description.h"


#include "lib_util.h"

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )
#define ABS(x) ((x)<0?-(x):(x))


#ifndef EPSILON
    #define EPSILON 0.0001
#endif

static inline int ori_to_bin(float ori, int nbins)
{
    if (ori < 0)
        ori +=2 *M_PI;
    int bin = (int)(ori/(2*M_PI)*nbins+0.5)%nbins;
    return bin;
}


static inline float bin_to_ori(float bin, int nbins)
{
    float ori = (bin+0.5)*2*M_PI/(float)nbins;
    if (ori > M_PI)
        ori -= 2*M_PI;
    return ori;
}


/** @brief Accumulate gradient orientation histogram around a keypoint
 *
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key      keypoint coordinates.
 *
 * @param imX
 * @param imX          precomputed gradient relative to the nearest scale.
 *
 * ----PARAM----
 * @param lambda_ori     (= 1.5)
 *                       - The patch P^ori is ( 6 X lambda_ori X sigma_key )
 *                       - The Gaussian window has a standard deviation of lambda_ori X sigma_key
 *
 * @param nbins          (= 36)  number of bins covering the range [0,2pi]
 *
 *
 * ----OUTPUT----
 * @output hist          the output is vec of bins value
 *
 * ----RETURN----
 * @return count         number of pixels contributing to the histogram.
 *
 */
void sift_accumulate_orientation_histogram(float x_key,
                                           float y_key,
                                           float sigma_key,
                                           const float* imX,
                                           const float* imY,
                                           int w,
                                           int h,
                                           int nbins,
                                           float lambda_ori,
                                           float* hist)
{
    /// Initialize output vector
    for(int i= 0;i<nbins;i++){hist[i] = 0.0;}

    /// Contributing pixels are inside a patch [siMin;siMax] X [sjMin;sjMax]
    // of width w 6*lambda_ori*sigma_key (=9*sigma_key)
    float R = 3*lambda_ori*sigma_key;
    int siMin = MAX(0, (int)(x_key-R+0.5));
    int sjMin = MAX(0, (int)(y_key-R+0.5));
    int siMax = MIN((int)(x_key+R+0.5), h-1);
    int sjMax = MIN((int)(y_key+R+0.5), w-1);

    /// For each pixel inside the patch.
    for(int si = siMin; si <= siMax; si++){
        for(int sj = sjMin; sj <= sjMax; sj++){

            /// Compute pixel coordinates (sX,sY) on keypoint's invariant
            //referential.
            float sX = (si-x_key)/sigma_key;
            float sY = (sj-y_key)/sigma_key;

            // gradient orientation (theta)
            float dx = imX[si*w+sj];
            float dy = imY[si*w+sj];
            float ori = modulus(atan2(dy, dx), 2*M_PI);

            // gradient magnitude with Gaussian weighing
            float r2 = sX*sX+sY*sY;
            float M = hypot(dx,dy) * exp(-r2/(2*lambda_ori*lambda_ori));

            /// Determine the bin index in the circular histogram
            int gamma = ori_to_bin(ori, nbins);

            /// Add the contribution to the orientation histogram
            hist[gamma] += M;
        }
    }
}


static float interpolate_peak(float h1, float h2, float h3)
{
    float offset = (h1-h3)/(2*(h1+h3-2*h2));
    return offset;
}


static void smooth_circular_histogram(int niter, float* hist, int nbins);



/** @brief Extract principal orientations from gradient orientation histogram
 *
 * ----INPUT----
 * @param hist           histogram of gradient orientation.
 * @param nbins          (=36) number of bins covering range [0,2pi].
 *
 * ----PARAM----
 * @param threshold      (0.8) local maxima are considered secondary
 *                       principal orientation if they exceed
 *                       threshold times the absolute maximum.
 *                       expressed in percentage of the max value
 *
 * ----OUTPUT----
 * @output principal_orientations      pointer to tab storing principal orientations.
 *
 * ----RETURN----
 * @return nori          number of principal orientation (>=1)
 */
int sift_extract_principal_orientations(float* hist,
                                        int nbins,
                                        float threshold,
                                        float* principal_orientations)
{
    // number of principal orientations ( the return value).
    int o = 0;

    // Smooth histogram : 6 iterated box filters
    smooth_circular_histogram(6, hist, nbins);
    // What is the value of the global maximum
    float max_value = array_max(hist, nbins);
    // Search for local extrema in the histogram
    for(int i = 0; i < nbins; i++){
        int i_prev = (i-1+nbins)%nbins;
        int i_next = (i+1)%nbins;
        if ( (hist[i] > threshold*max_value) && (hist[i]>hist[i_prev]) && (hist[i]>hist[i_next])){
            // Quadratic interpolation of the position of each local maximum
            float offset = interpolate_peak(hist[i_prev], hist[i], hist[i_next]);
            // Add to vector of principal orientations (expressed in [0,2pi]
            principal_orientations[o] = bin_to_ori((float)i + offset, nbins);
            o++;
        }
    }
    // return the number of principal orientations
    return o;
}



float sift_extract_one_orientation(float* hist, int nbins)
{
    float ori; // return value

    int i,i_prev,i_next;   // bin indices
    float offset;          // for normalization and interpolation

    // Smooth histogram : 6 iterated box filters
    smooth_circular_histogram(6, hist, nbins);

    // Find the histogram global extrema
    float t = find_array_max(hist, nbins, &i);
    (void)t;
    i_prev=(i-1+nbins)%nbins;
    i_next=(i+1)%nbins;
    // Quadratic interpolation of the position of each local maximum
    offset = interpolate_peak(hist[i_prev], hist[i], hist[i_next]);
    ori = bin_to_ori((float)i + offset, nbins);
    return ori;
}


/** @brief Extract keypoint feature vector
 *
 * ----INPUT----
 * @param x_key
 * @param y_key
 * @param sigma_key
 * @param theta_key      keypoint coordinates.
 *
 * @param imX
 * @param imX          precomputed gradient relative to the nearest scale.
 *
 * ----PARAM----
 * @param lambda_descr   (=6)
 *                       - The gaussian window has a standard deviation of lambda_descr X=* sigma_key
 *                       - The patch P^descr is ( 2 * lambda_descr * sigma_key * (1+1/Nhist) wide
 *
 * @param Nhist          (=4) number of histograms in each of the two directions,
 * @param Nbins          (=8) number of bins covering the range [0,2pi]
 *
 * ----OUTPUT----
 * @output descr
 *
 * ----RETURN----
 * @return count         number of sample contributing to the descriptor
 */
void sift_extract_feature_vector(float x_key, float y_key, float sigma_key,
                                 float theta_key,
                                 const float* imX, const float* imY,
                                 int w, int h,
                                 int Nhist,
                                 int Nbins,
                                 float lambda_descr,
                                 float* descr)
{
    // Initialize descr tab
    for(int i = 0; i < Nhist*Nhist*Nbins; i++){descr[i] = 0.0;}
    // Contributing pixels are inside a patch [siMin;siMax]X[sjMin;sjMax] of
    // width 2*lambda_descr*sigma_key*(nhist+1)/nhist
    float R =  (1+1/(float)Nhist)*lambda_descr*sigma_key;
    float Rp =  M_SQRT2*R;
    int siMin = MAX(0, (int)(x_key - Rp +0.5));
    int sjMin = MAX(0, (int)(y_key - Rp +0.5));
    int siMax = MIN((int)(x_key + Rp +0.5), h-1);
    int sjMax = MIN((int)(y_key + Rp +0.5), w-1);
    /// For each pixel inside the patch.
    for(int si = siMin; si < siMax; si++){
        for(int sj = sjMin; sj < sjMax; sj++){
            // Compute pixel coordinates (sX,sY) on keypoint's invariant referential.
            float X = si - x_key;
            float Y = sj - y_key;
            apply_rotation(X, Y, &X, &Y, -theta_key);
            // Does this sample fall inside the descriptor area ?
            if (MAX(ABS(X),ABS(Y)) < R) {
                // Compute the gradient orientation (theta) on keypoint referential.
                double dx = imX[si*w+sj];
                double dy = imY[si*w+sj];
                float ori = atan2(dy, dx) - theta_key;
                ori = modulus(ori, 2*M_PI);
                // Compute the gradient magnitude and apply a Gaussian weighing to give less emphasis to distant sample
                double t = lambda_descr*sigma_key;
                double M = hypot(dx, dy) * exp(-(X*X+Y*Y)/(2*t*t));

                // bin indices, Compute the (tri)linear weightings ...
                float alpha = X/(2*lambda_descr*sigma_key/Nhist) + (Nhist-1.0)/2.0;
                float beta  = Y/(2*lambda_descr*sigma_key/Nhist) + (Nhist-1.0)/2.0;
                float gamma = ori/(2*M_PI)*Nbins;
                //    ...and add contributions to respective bins in different histograms.
                // a loop with 1 or two elements
                int i0 = floor(alpha);
                int j0 = floor(beta);
                for(int i = MAX(0,i0);i<=MIN(i0+1,Nhist-1);i++){
                    for(int j = MAX(0,j0);j<=MIN(j0+1,Nhist-1);j++){ // looping through all surrounding histograms.

                        int k;
                        // Contribution to left bin.
                        k = ((int)gamma+Nbins)%Nbins;
                        descr[i*Nhist*Nbins+j*Nbins+k] += (1.-(gamma-floor(gamma)))
                                                         *(1.0-ABS((float)i-alpha))
                                                         *(1.0-ABS((float)j-beta))
                                                         *M;

                        // Contribution to right bin.
                        k = ((int)gamma+1+Nbins)%Nbins;
                        descr[i*Nhist*Nbins+j*Nbins+k] += (1.0-(floor(gamma)+1-gamma))
                                                         *(1.0-ABS((float)i-alpha))
                                                         *(1.0-ABS((float)j-beta))
                                                         *M;
                        

                    }
                }
            }
        }
    }
}





/** @brief Modify descriptor for robustness and speed
 *
 *  - Threshold bins exceeding the value of threshold times the l2 norm
 *  (to increase robustness to saturation effects)
 *
 *  - Normalized the vector (l2 norm = 512)
 *
 *  - Quantize vector values
 *  (to increase the speed of distance computations)
 *
 *
 * ----INPUT OUTPUT----
 * @param descr                 float histogram of length Nhist*Nhist*Nbins (=128)
 *
 * ----PARAM----
 * @param threshold             threshold applied to the l2-norm of the descriptor
 *
 * ----OUTPUT----
 * @param qtdescr               int histogram of length Nhist*Nhist*Nbins and values in [0..255]
 *
 */
void sift_threshold_and_quantize_feature_vector(float* descr, int n, float threshold)
{
    // Normalize
    float l2norm = array_l2norm(descr, n);
    // Threshold bins
    for(int i = 0; i < n; i++){
        descr[i] = MIN(descr[i],threshold*l2norm);
    }
    // Renormalize
    l2norm = array_l2norm(descr, n);
    // Quantization
    for(int i = 0; i < n; i++){
        descr[i] = (int)(descr[i]*512.0/l2norm);
        descr[i] = MIN(descr[i], 255);
    }
}




/** @brief Iterative box filter of w 3 bins
 *
 * ----INPUT OUTPUT----
 * @param hist           : histogram
 *
 * ----INPUT----
 * @param nbins          : number of bins
 *
 * ----PARAM----
 * @param niter          : number of iteration
 *
 */
static void smooth_circular_histogram(int niter, float* hist, int nbins)
{
    int i,i_prev,i_next;
    float tmp[nbins];
    /// Initialization
    for(i = 0; i < nbins; i++)
        tmp[i] = hist[i];
    /// Convolution with box filters
    for(;niter>0;niter--){
        for(i=0;i<nbins;i++)
            tmp[i] = hist[i];
        for(i=0;i<nbins;i++){
            i_prev = (i-1+nbins)%nbins;
            i_next = (i+1)%nbins;
            hist[i] = (tmp[i_prev]+tmp[i]+tmp[i_next])/3.;
        }
    }
}
