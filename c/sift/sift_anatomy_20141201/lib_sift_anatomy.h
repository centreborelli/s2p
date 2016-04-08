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



*/
/**
 * @file lib_sift_anatomy.h
 * @brief SIFT anatomy interface
 *
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#ifndef _LIB_SIFT_ANATOMY_H_
#define _LIB_SIFT_ANATOMY_H_


#include "lib_scalespace.h"
#include "lib_keypoint.h"


struct sift_parameters
{
    int n_oct, n_spo, n_hist, n_bins, n_ori, itermax;
    float sigma_min, delta_min, sigma_in, C_DoG, C_edge, lambda_ori, t, lambda_descr;
};

struct sift_parameters* sift_assign_default_parameters();

struct sift_keypoints* sift_anatomy(const float* x, int w, int h, const struct sift_parameters* p,
                                    struct sift_scalespace* ss[4],
                                    struct sift_keypoints* kk[6]);

struct sift_keypoints* sift_anatomy_without_description(const float* x, int w, int h, const struct sift_parameters* p,
                                    struct sift_scalespace* ss[2],
                                    struct sift_keypoints* kk[5]);

void sift_anatomy_only_description(const float* x, int w, int h, const struct sift_parameters* p, struct sift_keypoints* k);

void sift_anatomy_orientation_and_description(const float* x, int w, int h, const struct sift_parameters* p, struct sift_keypoints* k);


void scalespace_compute(struct sift_scalespace* ss,
                               const float* image,
                               int im_w,
                               int im_h,
                               float sigma_in);

void scalespace_compute_gradient(const struct sift_scalespace* scalespace,
                                 struct sift_scalespace* sx,
                                 struct sift_scalespace* sy);




#endif // _LIB_SIFT_ANATOMY_H_
