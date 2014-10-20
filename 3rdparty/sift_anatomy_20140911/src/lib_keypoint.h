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
 * @file lib_keypoint.h
 * @brief data structures to store information relative to keypoint
 *
 * @li struct keypoint      keypoint data structure.
 * @li struct sift_keypoints  list of keypoints with variable length.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */




#ifndef _LIB_KEYPOINT_H_
#define _LIB_KEYPOINT_H_


/** @brief keypoint structure, related to a keypoint
 *
 * stores SIFT output (keypoint position, scale, orientation, feature vector...)
 * stores intermediary results (orientation histogram, summarizes...)
 *
 *
 */
struct keypoint
{
    float x; // coordinates
    float y;
    float sigma; // level of blur (it includes the assumed image blur)
    float theta; // orientation

    int o; // discrete coordinates
    int s;
    int i;
    int j;

    float val; // normalized operator value (independant of the scalespace sampling)
    float edgeResp; // edge response
    int n_hist;     // number of histograms in each direction
    int n_ori;      // number of bins per histogram
    float* descr;

    int n_bins; // number of bins in the orientation histogram
    float* orihist; // gradient orientation histogram
};

struct keypoint* sift_malloc_keypoint(int n_ori, int n_hist, int n_bins);

struct keypoint* sift_malloc_keypoint_from_model_and_copy(const struct keypoint* key);


/** @brief list of pointers to keypoint structures (variable size list)
 *      "vector" of keypoint struct (list of variable size)
 */
struct sift_keypoints{
    int size;             /* number of elements in the list */
    int capacity;            /* current max number os elements in th list (size<=capacity)*/
    struct keypoint ** list;   /* array of pointers to keypoint structure */
};

/** @brief allocate list of keypoints */
struct sift_keypoints* sift_malloc_keypoints();

/** @brief add element in list of keypoints */
void sift_add_keypoint_to_list(struct keypoint* key, struct sift_keypoints* keys);

struct keypoint* sift_add_keypoint_to_listPASC(struct sift_keypoints* keypoints, int n_ori, int n_hist, int n_bins);

/** @brief free list of keypoints */
void sift_free_keypoints(struct sift_keypoints *keys);

    
/** @brief Save list of keys on a txt filea (3 formats) 
 *  flag values
 *       flag >= 0  ->  prints coordinates (x,y,sigma,theta)
 *       flag >= 1  ->  prints descriptor
 *       flag >= 2  ->  prints scalespace sample coordinates
 *                             + orientation histogram.
 */                       
void sift_save_keypoints(const struct sift_keypoints* keys, const char* name, int flag);

/** @brief Save list of keys on the standard output
 */
void sift_print_keypoints(const struct sift_keypoints* keys, int flag);


/** @brief Read list of keys from a txt file  */
void sift_read_keypoints(struct sift_keypoints* keys,
                         const char* name,
                         int n_hist,
                         int n_ori,
                         int n_bins, 
                         int flag);

// Print a single keypoint
void fprintf_one_keypoint(FILE* f, const struct keypoint *k, int n_descr, int n_bins, int flag);

#endif // _LIB_KEYPOINT_H_
