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
 * @file sift_keypoint.h
 * @brief data structures to store information relative to keypoint
 *
 * @li struct kypt      keypoint data structure.
 * @li struct kypt_lst  list of keypoints with variable length.
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */




#ifndef _SIFT_KEYPOINT_H_
#define _SIFT_KEYPOINT_H_

#include <stdio.h>
#include <stdlib.h>

/** 4 structures are declared here :
 *     - kypt         (keypoint)
 *     - kypt_lst     (list of keypoints)
 *     - intTrio      (3d coordinate)
 *     - intTrio_lst  (list of 3D coordinates)*/


/** @brief kypt structure, related to a keypoint
 *
 * stores SIFT output (position, scale, orientation, feature vector...)
 * stores intermediary results (orientation histogram, summarizes...)
 *
 */
struct kypt{

    /* absolute coordinates, orientation and feature vector  */
    double x;
    double y;
    double sigma;   /* level of blur (includes the initial image blur), not the added blur */
    double theta;   
    double* descr;   

    /* discrete coordinates  */
    int o;
    int s;
    int i;
    int j;
    
    double val;         /* stored value of the normalized analysis operator */
                        /* independant of the sampling of the scalespace    */
    
    double hsEdgeResp; /* Harris and Stephen edge response applied on */
                       /* on the normalized analysis operator.        */
             
    /* feature vector dimension */
    int n_ori;          /* number of bins per histogram */ 
    int n_hist;         /* number of histogram in each direction */
    int nsamples_descr; /* number of samples in the normalized patch P_descr */
    
    /* extrema interpolation */
    int nMaxIntrp;        /* number of max iteration for extrema interpolation */
    int nIntrp;           /* number of iteration of the extrema interpolation */
    int* pos_intrp;       /* successive position where the interpolation is computed */
    double* offset_intrp; /* successive offsets */
     
    /* computation of the reference orientation */ 
    int n_bins   ;       /* number of bins in the orientation histogram */
    double* hist_prOri;  /* histogram of gradient orientation */
    int n_prOri;         /* number of principal orientations  */
    double* ori_prOri;   /* list of principal orientations    */
    int nsamples_ori;    /* number of samples in normalized patch P_ori*/
};

struct kypt* malloc_kypt(int n_ori, int n_hist);
struct kypt* malloc_kypt_extra(int n_ori, int n_hist, int n_bins, int nMaxIntrp);
struct kypt* malloc_kypt_from_model_and_copy(struct kypt* key);



/** @brief list of pointers to kypt structures (variable size list)
 *      "vector" of kypt struct (list of variable size)
 */
struct kypt_lst{
    int card;             /* number of elements in the list */
    int msize;            /* current max number os elements in th list (card<=msize)*/
    struct kypt ** lst;   /* array of pointers to kypt structure */
};


/** @brief allocate list of keypoints */
struct kypt_lst* malloc_kypt_lst();
/** @brief increase the size of list of keypoints */
void realloc_kypt_lst(struct kypt_lst* keys);
/** @brief add element in list of keypoints */
void add_kypt_to_lst(struct kypt* key, struct kypt_lst* keys);

/** @brief free keypoint */
void free_kypt(struct kypt* keypoint);
void free_kypt_extra(struct kypt* keypoint);
/** @brief free list of keypoints */
void free_kypt_lst(struct kypt_lst *keys);
// void free_kypt_extra_lst(struct kypt_lst* keys);

/** @brief Save list of keys on a txt file (3 formats) */
void save_kypt_lst_detection(struct kypt_lst * keys, char *name);
void save_kypt_lst_description(struct kypt_lst* keys, int dim, char* name);
void save_kypt_lst_extra(struct kypt_lst * keys, int dim, int n_bins, char* name);

/** @brief Print in the standard output (2 formats) */
void print_kypt_lst(struct kypt_lst* keys, int dim);
void print_kypt_lst_extra(struct kypt_lst* keys, int dim, int n_bins);

/** @brief Read list of keys from a txt file  */
void read_kypt_lst(struct kypt_lst* keys, char* basename, int n_hist, int n_ori);
void read_kypt_lst_extra(struct kypt_lst* keys, char* basename, int n_hist, int n_ori, int n_bins);




/** @brief integer triplets */
struct intTrio{
    int i;
    int j;
    int s;
};


/** @brief list of pointers to triplets of integers
 *    "vector" of triplets (list of variable size)
 */
struct intTrio_lst{
    struct intTrio ** lst;  /* array o pointer to triplets */
    int card;               /* size of list of keypoint */
    int msize;              /* size of allocated memory - ( need card<=msize)  */
};

/** @brief Allocate list of keypoints */
struct intTrio_lst* malloc_intTrio_lst();
/** @brief Increase the maximum size of the list */
void realloc_intTrio_lst(struct intTrio_lst* coords);
/** @brief Add element to the list (push_back)*/
void add_intTrio_to_lst(struct intTrio* coord, struct intTrio_lst * coords);
/** @brief Free allocated memory */
void free_intTrio_lst(struct intTrio_lst * coords);





#endif