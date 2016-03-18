/*
IPOL SIFT
Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20140911 (September 11th, 2014)

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
 * @file sift_matching.c
 * @brief data structures to store information relative to a pair of keypoints
 *
 * @li struct keypointPr     : Pair of keypoint data structure.
 * @li struct keypointPr_list : List of pairs.
 * @li print,save, read for lists of pairs.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lib_keypoint.h"
#include "lib_matching.h"
#include "lib_util.h"

static void compute_keypoints_distance(float* dist,
                                       const struct sift_keypoints *k1,
                                       const struct sift_keypoints *k2)
{
    int n_hist = k1->list[0]->n_hist;
    int n_ori  = k1->list[0]->n_ori;
    int dim = n_hist*n_hist*n_ori;
    int n1 = k1->size;
    int n2 = k2->size;
    for(int i = 0; i < n1; i++){
        const float * v1 =  k1->list[i]->descr;
        for(int j = 0; j < n2; j++){
            const float * v2 =  k2->list[j]->descr;
            float d = euclidean_distance(v1, v2, dim);
            dist[i*n2+j] = d;
        }
    }
}


static void find_the_two_nearest_keys(const float* dist, int n1, int n2,
                                     int* indexA, int* indexB,
                                     float* distA, float* distB)
{
    for(int i = 0; i < n1; i++){
        int iA, iB;
        float dA, dB;
        find_array_two_min(&dist[i*n2], n2, &dA, &dB, &iA, &iB);
        indexA[i] = iA;
        indexB[i] = iB;
        distA[i] = dA;
        distB[i] = dB;
    }
}

void matching(struct sift_keypoints *k1,
              struct sift_keypoints *k2,
              struct sift_keypoints *out_k1,
              struct sift_keypoints *out_k2A,
              struct sift_keypoints *out_k2B,
              float thresh,
              int flag)
{
    int n1 = k1->size;
    int n2 = k2->size;

    float* dist  = (float*)xmalloc(n1*n2*sizeof(float));
    float* distA = (float*)xmalloc(n1*sizeof(float));
    float* distB = (float*)xmalloc(n1*sizeof(float));
    int* indexA  = (int*)xmalloc(n1*sizeof(int));
    int* indexB  = (int*)xmalloc(n1*sizeof(int));

    compute_keypoints_distance(dist, k1, k2);
    find_the_two_nearest_keys(dist, n1, n2, indexA, indexB, distA, distB);

    int j = 0;
    for(int i = 0; i < n1; i++){
        float val;
        val = (flag == 1 ? distA[i]/distB[i] : distA[i]);
        if (val < thresh){
            int iA = indexA[i];
            int iB = indexB[i];
            struct keypoint* k;
            k = sift_malloc_keypoint_from_model_and_copy(k1->list[i]);
            sift_add_keypoint_to_list(k, out_k1);
            k = sift_malloc_keypoint_from_model_and_copy(k2->list[iA]);
            sift_add_keypoint_to_list(k, out_k2A);
            k = sift_malloc_keypoint_from_model_and_copy(k2->list[iB]);
            sift_add_keypoint_to_list(k, out_k2B);
            j++;
        }
    }

    free(dist);
    free(indexA);
    free(indexB);
    free(distA);
    free(distB);
}

void print_pairs(const struct sift_keypoints *k1,
                 const struct sift_keypoints *k2)
{
    if (k1->size > 0){
        int n = k1->size;
        for(int i = 0; i < n ;i++){
            fprintf_one_keypoint(stdout, k1->list[i], 0, 0, 0);
            fprintf_one_keypoint(stdout, k2->list[i], 0, 0, 0);
            fprintf(stdout, "\n");
        }
    }
}

void save_pairs_extra(const char* name,
                      const struct sift_keypoints *k1,
                      const struct sift_keypoints *k2A,
                      const struct sift_keypoints *k2B)
{
    FILE* f = fopen(name,"w");
    
    if (k1->size > 0){

        int n_hist = k1->list[0]->n_hist;
        int n_ori = k1->list[0]->n_ori;
        int dim = n_hist*n_hist*n_ori;
        int n_bins  = k1->list[0]->n_bins;
        int n = k1->size;
        for(int i = 0; i < n; i++){
            fprintf_one_keypoint(f, k1->list[i], dim, n_bins, 2);
            fprintf_one_keypoint(f, k2A->list[i], dim, n_bins, 2);
            fprintf_one_keypoint(f, k2B->list[i], dim, n_bins, 2);
            fprintf(f, "\n");
        }
    }
    fclose(f);
}
