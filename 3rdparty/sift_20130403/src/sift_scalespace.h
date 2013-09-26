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



#ifndef _SIFT_SCALESPACE_H_
#define _SIFT_SCALESPACE_H_

#include <stdio.h>
#include <stdbool.h>

#include "discrete_representation.h"
#include "std_image_lib.h"

#include <assert.h>


/** ***********************  STRUCTURES *************************/

struct octa{
    double delta;       /* sampling rate in this octave */
    int width;          /* the same for all images in the stack */
    int height;         
    int nSca;           /* number of images in the stack */
    double* sigmas;     /* intrinsec levels of blur */
    double* imStack;    /* stack of nSca images of width*height samples (indexed from fine to coarse)*/
};

struct scsp{
    int nOct;                   /* number of octaves */
    struct octa* octaves[100];  /* array of pointer to octave structure */
};


/** ************************ ALLOCATION **************************/

struct octa* malloc_octa(double delta, int width, int height, int nSca, double* sigmas);
struct scsp* malloc_scsp(int nOct, double* deltas, int* widths, int* heights, int* nScas, double** sigmas);

void free_octa(struct octa *octa); /*frees memory of octave culture*/
void free_scsp(struct scsp *scsp); /*frees memory for scalespace structure*/


/** ********************* ALLOCATION WRAPPERS *********************************/

struct scsp* malloc_scsp_from_model(struct scsp * model_scsp);

struct scsp* malloc_scsp_lowe_scalespace(int nOct,                    
                                         int nSca,                    /* # of scales of detection (excluding the 3 auxiliary scales)*/
                                         int width, int height,       /* relative to input image*/
                                         bool zoomin,                 /* flag for prior zoom     */
                                         double sigma_min);           /* minimal scale in each octave (relatively to the sampling rate) */

struct scsp* malloc_scsp_dog_from_scalespace(struct scsp* scalespace);

/** @brief Allocation of a scalespace with David Lowe's original structure.
 *
 */
struct scsp* malloc_scsp_lowe_scalespace_bis(int nOct,                 /* # of octaves            */
                                         int nSca,                     /* # of scales of detection (excluding the 3 auxiliary scales) */
                                         int im_width, int im_height,  /* # input image dimension */
                                         double delta_min,             /* minimal inter-sample distance */
                                         double sigma_min);            /* minimal scale in each octave (relatively to the sampling rate) */


void build_octave(double* in,                         /* input image */
                  int width, int height,              /* image dimension */
                  double sigma0, double delta,        /* image blur  / sampling rate */
                  int nSca, double* sigmas,           /* # of scales / array of simulated scales */ 
                  double* oct);                       /* array of nSca*width*height samples */


void print_scsp_gray(struct scsp* scalespace, char* basename);
void print_scsp_gray_nearestneighor(struct scsp* scalespace, char* basename);
void print_scsp_rgb(struct scsp* scalespace, char* basename);


#endif
