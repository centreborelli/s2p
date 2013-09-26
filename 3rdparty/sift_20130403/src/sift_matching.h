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
 * @file sift_matching.h
 * @brief data structures to store information relative to a pair of keypoints
 *
 * @li struct kyptPr     : Pair of keypoint data structure.
 * @li struct kyptPr_lst : List of pairs.
 * @li print,save, read for lists of pairs.
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */




#ifndef _SIFT_MATCHING_H_
#define _SIFT_MATCHING_H_

#include <math.h> // for the sqrt of the Euclidian norm
#include "sift_keypoint.h"


/** @brief data structure : keypoint pair.
 * 
 *   - points to three keypoints
 *   - stores distance to first and second neighbor
 */
struct kyptPr{
                            /* pointers to the keypoints  */
    struct kypt *key1;      /* key in image 1 */
    struct kypt *key2a;     /* first neighbor in image 2  */
    struct kypt *key2b;     /* second neighbor in image 2 */
    
    double distance_1_2a;   /* distance to first neighbor  */
    double distance_1_2b;   /* distance to second neighbor */
};

/** @brief list of keypoint pairs
 * 
 */
struct kyptPr_lst{
    int card;
    int msize;
    struct kyptPr** lst;
};
struct kyptPr_lst* malloc_kyptPr_lst();
void add_kyptPr_to_lst(struct kyptPr* pair, struct kyptPr_lst *pairs);
void realloc_kyptPr_lst(struct kyptPr_lst *pairs);
void free_kyptPr_lst(struct kyptPr_lst *pairs);

/** @brief Pairing two sets of keypoints */
void pairing_keyslist_to_keyslist(struct kyptPr_lst* pairs,
                                  struct kypt_lst* keys1,
                                  struct kypt_lst* keys2,
                                  int dim);


void only_matching_kyptPr(double thresh,
                          int flag,
                          struct kyptPr_lst* pairs,
                          struct kyptPr_lst* matches,
                          struct kypt_lst* keys1Matching,
                          struct kypt_lst* keys2Matching,
                          struct kypt_lst* keys1NotMatching);

void save_kyptPr_extra(struct kyptPr_lst* pairs, int dim, int n_bins, char* name);
void print_pairs(struct kyptPr_lst* pairs);


#endif
