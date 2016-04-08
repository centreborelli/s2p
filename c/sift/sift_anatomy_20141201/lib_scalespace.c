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
 * @file sift_scalespace.c
 * @brief data structures to store the scale-space
 *
 * @li struct keypoint      keypoint data structure.
 * @li struct sift_keypoints  list of keypoints with variable length.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "lib_discrete.h"
#include "lib_scalespace.h"
#include "lib_util.h"



/** ************************************ ALLOCATION *******************************************/


/** @brief Octave structure memory allocation
 *
 *  An octave is a stack of images sharing same w, h and intersample distance.
 *  Each image of the stack has a level of blur.
 * 
 * @param nSca   : number of images in the stack.
 * @param delta  : intersample distance
 * @param w 
 * @param h
 * @param sigma  : levels of blur for each image in the stack
 *
 *
 */
static struct octa *malloc_octa(float delta, int w, int h, int nSca, const float* sigmas){

    /* Memory allocation */
    struct octa* octave = xmalloc(sizeof(struct octa));   /* pointer to a (structure octave) */

    /* Modifying structure members */
    octave->delta  = delta;
    octave->w  = w;
    octave->h = h;
    octave->nSca   = nSca;
    octave->sigmas  = xmalloc(nSca*sizeof(float));                /* pointer to the array of level of simulated blurs */
    octave->imStack = xmalloc(nSca*w*h*sizeof(float));   /* pointer to the array of pixels  */
    for(int i=0;i<nSca;i++){
        octave->sigmas[i] = sigmas[i];
    }
    return octave;
}

/** @brief allocation and attributes modification
 *
 *       Allocate memory for a scalespace structure 
 *       and modify attributes that will be used to actually build the scalespace
 * 
 * @param nOct     : number of octaves
 * @param deltas   : intersample distance in each octave of the scalespace
 * @param ws   : image w in each octave of the scalespace
 * @param hs  :
 * @param nScas    : number of images in each octave of the scalespace
 * @param sigmas   : array of blurs for each image in each octave of the scalespace
 * 
 */
struct sift_scalespace* sift_malloc_scalespace(int nOct,
                                               const float* deltas,
                                               const int* ws,
                                               const int* hs,
                                               const int* nScas,
                                               float** sigmas)
{
    struct sift_scalespace* scalespace = xmalloc(sizeof(struct sift_scalespace));
    scalespace->nOct = nOct;
    for(int o=0;o<nOct;o++){
        scalespace->octaves[o] = malloc_octa(deltas[o], ws[o], hs[o], nScas[o], sigmas[o]);
    }
    return scalespace;
}



static void free_octa(struct octa *octave){
    xfree(octave->sigmas);    /* the memory allocated to store simulated values */
    xfree(octave->imStack);   /* the memory allocated to store the octave samples */
    xfree(octave);
}


void sift_free_scalespace(struct sift_scalespace *scalespace){
    int nOct = scalespace->nOct;
    for(int o=0;o<nOct;o++)
        free_octa(scalespace->octaves[o]);
    xfree(scalespace);
}



/** **********************************  ALLOCATION WRAPPERS ****************************************/


/** @brief Allocation via copy a scalespace structure
 * 
 * - Copy structure (number of octaves, dimensions, sampling rates, levels of blur)
 * - Doesn't copy the image stacks
 * 
 */
struct sift_scalespace* sift_malloc_scalespace_from_model(struct sift_scalespace* model_sift_scalespace){

    struct sift_scalespace * copy_sift_scalespace;

    int nOct = model_sift_scalespace->nOct;

    int* nScas = xmalloc(nOct*sizeof(int));
    int* ws = xmalloc(nOct*sizeof(int));
    int* hs = xmalloc(nOct*sizeof(int));
    float* deltas  = xmalloc(nOct*sizeof(float));
    float** sigmas = xmalloc(nOct*sizeof(float*));

    for(int o=0;o<nOct;o++){
        nScas[o]   = model_sift_scalespace->octaves[o]->nSca;
        ws[o] = model_sift_scalespace->octaves[o]->w;
        hs[o] = model_sift_scalespace->octaves[o]->h;
        deltas[o]  = model_sift_scalespace->octaves[o]->delta;

        sigmas[o]  = xmalloc(nScas[o]*sizeof(float));
        for(int i=0;i<nScas[o];i++)
            sigmas[o][i] = model_sift_scalespace->octaves[o]->sigmas[i];
    }

    copy_sift_scalespace = sift_malloc_scalespace(nOct, deltas, ws, hs, nScas, sigmas);

    xfree(deltas);
    xfree(ws);
    xfree(hs);
    xfree(nScas);
    for(int o=0;o<nOct;o++)
        xfree(sigmas[o]);
    xfree(sigmas);

    return copy_sift_scalespace;
}









/** @brief Allocation of a scalespace with David Lowe's original structure.
 *
 *
 * The set of simulated level of blur
 *    sigmas[o][s] = deltas[o]*sigma_min*pow(2.0,(float)s/(float)nSca);
 *
 *
 */
struct sift_scalespace* sift_malloc_scalespace_lowe(int nOct,                 /* # of octaves            */
                                                    int nSca,           /* # of scales of detection (excluding the 3 auxiliary scales */
                                                    int im_w, int im_h,  /* # input image dimension */
                                                    float delta_min,             /* minimal inter distance sample     */
                                                    float sigma_min)            /* minimal scale in each octave (relatively to the sampling rate) */
{
    int* ws  = xmalloc(nOct*sizeof(int));
    int* hs  = xmalloc(nOct*sizeof(int));
    int* nScas      = xmalloc(nOct*sizeof(int));
    float* deltas  = xmalloc(nOct*sizeof(float));
    float** sigmas = xmalloc(nOct*sizeof(float*));    /* nSca+3 = nSca + 2 extra scales to search 3d extrema + 1 extra scale to compute DoG */
    assert(delta_min <=1);
    deltas[0] = delta_min;
    hs[0] = (int)(im_h/delta_min);
    ws[0] = (int)(im_w/delta_min);
    for(int o=1;o<nOct;o++){
        ws[o] = ws[o-1]/2;     /*integer division*/
        hs[o] = hs[o-1]/2;
        deltas[o] = deltas[o-1]*2.0;
    }
    for(int o=0;o<nOct;o++){
        nScas[o] = nSca+3;  /* 3 extra images in the stack, 1 for dog computation and 2 for 3d discrete extrema definition*/
        sigmas[o] = xmalloc(nScas[o]*sizeof(float));
        for(int s=0;s<nSca+3;s++){ /* nSca images + 3 auxiliary images*/
            sigmas[o][s] = deltas[o]/deltas[0]*sigma_min*pow(2.0,(float)s/(float)nSca);
        }
    }
    struct sift_scalespace* scalespace = sift_malloc_scalespace(nOct, deltas, ws, hs, nScas, sigmas);   
    xfree(deltas);
    xfree(ws);
    xfree(hs);
    xfree(nScas);
    for(int o=0;o<nOct;o++)
        xfree(sigmas[o]);
    xfree(sigmas);    
    return scalespace;
}






/** @brief allocation for Difference of Gaussian scalespace.
 * 
 *   Characteristics of a DoG scalespace
 *      - one image less in each octave than the scalespace it is based on
 *      - same level of blurs than the scalespace it is based on (this is an approximation)
 * 
 */
struct sift_scalespace* sift_malloc_scalespace_dog_from_scalespace(struct sift_scalespace* scalespace)
{
    /* minimal scale in each octave (relatively to the sampling rate) */

    /* copy and allocate using a model */
    struct sift_scalespace* dog = sift_malloc_scalespace_from_model(scalespace);

    /* modify some allocations */
    int nOct = scalespace->nOct;
    int nSca;
    int w;
    int h;

    for(int o=0;o<nOct;o++){

        xfree(dog->octaves[o]->sigmas);
        xfree(dog->octaves[o]->imStack);

        nSca   = scalespace->octaves[o]->nSca-1;
        w  = scalespace->octaves[o]->w;
        h = scalespace->octaves[o]->h;

        dog->octaves[o]->nSca = nSca;
        dog->octaves[o]->sigmas  = xmalloc(nSca*sizeof(float));
        dog->octaves[o]->imStack = xmalloc(nSca*w*h*sizeof(float));

        for(int i=0;i<nSca;i++){
            dog->octaves[o]->sigmas[i] = scalespace->octaves[o]->sigmas[i];
        }
    }
    return dog;
}
