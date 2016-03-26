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



#ifndef _LIB_SCALESPACE_H_
#define _LIB_SCALESPACE_H_



/** ***********************  STRUCTURES *************************/

struct octa{
    float delta;       /* sampling rate in this octave */
    int w;             /* the same for all images in the stack */
    int h;
    int nSca;          /* number of images in the stack */
    float* sigmas;     /* intrinsic levels of blur */
    float* imStack;    /* stack of nSca images of w*h samples (indexed from fine to coarse)*/
};

struct sift_scalespace{
    int nOct;                   /* number of octaves */
    struct octa* octaves[100];  /* array of pointer to octave structure */
};


/** ************************ ALLOCATION **************************/

struct sift_scalespace* sift_malloc_scalespace(int nOct, const float* deltas, const int* ws, const int* hs, const int* nScas, float** sigmas);
void sift_free_scalespace(struct sift_scalespace *sift_scalespace); /*frees memory for scalespace structure */


/** ********************* ALLOCATION WRAPPERS *********************************/

struct sift_scalespace* sift_malloc_scalespace_from_model(struct sift_scalespace * model_sift_scalespace);
struct sift_scalespace* sift_malloc_scalespace_dog_from_scalespace(struct sift_scalespace* scalespace);


/** @brief  Lowe's original structure.
 */
struct sift_scalespace* sift_malloc_scalespace_lowe(int nOct,  /* # of octaves  */
                                         int nSca,            /* # of scales of detection (excluding the 3 auxiliary scales) */
                                         int im_w, int im_h,  /* # input image dimension */
                                         float delta_min,     /* minimal inter-sample distance */
                                         float sigma_min);    /* minimal scale in each octave (relatively to the sampling rate) */


#endif
