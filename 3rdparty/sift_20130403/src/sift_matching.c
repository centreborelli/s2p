/*
IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20130318 (March 18, 2013)

== Patent Warning and Licence =================================================

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
 * @li struct kyptPr     : Pair of keypoint data structure.
 * @li struct kyptPr_lst : List of pairs.
 * @li print,save, read for lists of pairs.
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#include "sift_matching.h"


/**  @brief Compute the L2 distance between the feature vectors of two keypoints
 *  
 *  
 * outputs the 
 * 
 */
double compute_keys_distance(struct kypt* k1, struct kypt* k2, int dim){
    double d = 0.0 ;
    double* v1 = k1->descr;
    double* v2 = k2->descr;
    for(int i=0;i<dim;i++){
        d +=(v1[i]-v2[i])*(v1[i]-v2[i]);
    }
    return sqrt(d);
}


/** @brief Allocate memory for list of keypoints pairs */
struct kyptPr_lst* malloc_kyptPr_lst(){
    struct kyptPr_lst* pairs = (struct kyptPr_lst*)malloc(sizeof(struct kyptPr_lst));
    pairs->lst = (struct kyptPr**)malloc(200*sizeof(struct kyptPr*));
    pairs->card=0;
    pairs->msize=200;
    return pairs;
}
/** @brief Increase the size of the list of keypoints pairs */
void realloc_kyptPr_lst(struct kyptPr_lst *pairs){
    pairs->lst = (struct kyptPr**)realloc(pairs->lst, (2*pairs->msize)*sizeof(struct kyptPr*));
    pairs->msize = 2*pairs->msize;
}
/** @brief Add one element to the list of keypoint pairs */
void add_kyptPr_to_lst(struct kyptPr *pair, struct kyptPr_lst *pairs){
    if (pairs->card > pairs->msize-5){
        realloc_kyptPr_lst(pairs);
    }
    pairs->lst[pairs->card] = pair;
    pairs->card +=1;
}

/** @brief free memory */
void free_kyptPr_lst(struct kyptPr_lst *pairs){
    for(int k=0;k<pairs->card;k++){
        free(pairs->lst[k]);
    }
    free(pairs->lst);
    free(pairs);
}



/** @brief Find the nearest neighbor o of keypoints to another
 * 
 *  in  @param key1      one keypoints in image 1
 *  in  @param keys2     list of all keypoints in image 2
 * 
 *  out @param  pair     keypoint 1 paired with its nearest neighbor 
 * 
 * 
 *  @param dim dimension of the feature vector
 * 
 */
void pairing_key_to_keyslist(struct kyptPr* pair,
                             struct kypt* key1,
                             struct kypt_lst* keys2,
                             int dim){
    
    /* compute distance of key1 to all keys in list keys2 */
    double* dists = (double*)malloc(keys2->card*sizeof(double));
    for(int i=0;i<keys2->card;i++){
        dists[i] = compute_keys_distance(key1, keys2->lst[i], dim);
    }

    /* Where are the first and second nearest neighbours in the list of keys ?*/
    /* Loop through all keypoints in second list */
    int idx_first;
    int idx_second;
    double dist_min;
    
    if (dists[0]<dists[1]) {
        idx_first  = 0;
        idx_second = 1;
        dist_min = dists[0];
    }
    else {
        idx_first  = 1;
        idx_second = 0;
        dist_min = dists[1];
    }
    for(int i=2;i<keys2->card; i++) {
        if (dists[i] < dist_min) {
            idx_second = idx_first;
            idx_first = i;
            dist_min = dists[i];
        }
    }
    
    /* Saving in kyptPr structure */
    pair->key1  = key1;
    pair->key2a = keys2->lst[idx_first];
    pair->key2b = keys2->lst[idx_second];
    pair->distance_1_2a = dists[idx_first];
    pair->distance_1_2b = dists[idx_second];
    free(dists);
}


/** @brief Pairing a list of keypoints to another
 * 
 *  out @param  pairs    list of keypoint pairs 
 *  in  @param keys1     list of keypoints in image 1
 *  in  @param keys2     list of keypoints in image 2
 * 
 *   @param dim dimension of the feature vector
 * 
 * Each keypoint in keys1 is paired with one keypoint in keys2
 *                             its nearest neighbor
 * 
 */
void pairing_keyslist_to_keyslist(struct kyptPr_lst* pairs,
                                  struct kypt_lst* keys1,
                                  struct kypt_lst* keys2,
                                  int dim){
    for(int i=0;i<keys1->card;i++){
        struct kyptPr* pair = (struct kyptPr*)malloc(1*sizeof(struct kyptPr));
        pairing_key_to_keyslist(pair, keys1->lst[i], keys2, dim);
        add_kyptPr_to_lst(pair, pairs);
    }
}



/** @brief Selecting the matching pairs 
 * 
 *  @param flag : method flag defining the matching criteria
 *  
 *  @param thresh (depending on the value of flag)
 *         if flag == 0 thresh over the distance to the nearest neighbor.
 *         if flag == 1 thresh over the distance ratio to the two nearest neighbors.
 *
 *  in  @param pairs  : all pairs of keypoints
 *  out @param matches : only the matching pairs.
 * 
 *  out @param keys1Matching :    Copy of the list of the matching keypoints in image 1.
 *  out @param keys2Matching :    Copy of the list of the matching keypoints in image 2.
 *  out @param keys1NotMatching : Copy of the list of keypoints not matching in image 1.
 *
 *
 * 
 */
void only_matching_kyptPr(double thresh, int flag,
                          struct kyptPr_lst* pairs,
                          struct kyptPr_lst* matches,
                          struct kypt_lst* keys1Matching,
                          struct kypt_lst* keys2Matching,
                          struct kypt_lst* keys1NotMatching){

    int isMatch;
    
    struct kyptPr *pair;
    struct kypt *k1_copy;
    struct kypt *k2_copy;
    
    for(int i=0;i<pairs->card;i++){     
        pair = pairs->lst[i];
        if(flag==0){
            /** Threshold on distance to first neighbor  */
            isMatch = (pair->distance_1_2a < thresh);
        }
        else{
            /** Threshold on distance ratio between first and second neighbors */
            isMatch = (pair->distance_1_2a/pair->distance_1_2b)<thresh;
        }
        
        if (isMatch){    
            add_kyptPr_to_lst(pair,matches);
            k1_copy = malloc_kypt_from_model_and_copy(pair->key1);
            k2_copy = malloc_kypt_from_model_and_copy(pair->key2a);
            add_kypt_to_lst(k1_copy, keys1Matching);
            add_kypt_to_lst(k2_copy, keys2Matching);
        }
        else {
            k1_copy = malloc_kypt_from_model_and_copy(pair->key1);
            add_kypt_to_lst(k1_copy, keys1NotMatching);
        }
    }
}



/** @brief Save exhaustive information on the output of the matching 
 * 
 *     The format of the output is
 *               keydata1   keydata2a   keydata2b
 *                # keypoint in image 1
 *                           # keydata 2a (nearest neighbor in image 2)
 *                                       # keydata 2b (second nearest image in image 2) 
 *     with 
 *           keydata = x y sigma theta o s  fvec[] ori_hist[]
 *                         #scale
 *                               #principal orientation
 *                                     # octave index
 *                                       # scale index  
 *                                          # feature vector.
 *                                                 # gradient orientation vector.
 * 
 *  @param pairs : data structure (list of kypt pairs)
 *  @param dim : dimension of the feature vector
 *  @param n_bins : number of bins in the orientation vector
 * 
 */
void save_kyptPr_extra(struct kyptPr_lst* pairs, int dim, int n_bins, char* name){
    FILE* file = fopen(name,"w");
    for(int i=0; i<pairs->card;i++){
        
        /* x1 , y1 , sigma1, theta1 */        
        fprintf(file,"%20.5f%20.5f%20.5f%20.5f", pairs->lst[i]->key1->x,
                                                 pairs->lst[i]->key1->y,
                                                 pairs->lst[i]->key1->sigma,
                                                 pairs->lst[i]->key1->theta);
        /* feature vector 1 */
        for(int n=0;n<dim;n++){fprintf(file," %10.3i", (int)pairs->lst[i]->key1->descr[n]);} 
        /* octave index and scale index */
        fprintf(file,"     %10i%10i ",pairs->lst[i]->key1->o, pairs->lst[i]->key1->s);
        /* orientation histo 1 */
        for(int n=0;n<n_bins;n++){fprintf(file," %10.3f", pairs->lst[i]->key1->hist_prOri[n]);} 
            
        /* x2a , y2a , sigma2a, theta2a fvec o s orihist*/
        fprintf(file,"%20.5f%20.5f%20.5f%20.5f", pairs->lst[i]->key2a->x,
                                                 pairs->lst[i]->key2a->y,
                                                 pairs->lst[i]->key2a->sigma,
                                                 pairs->lst[i]->key2a->theta);
        for(int n=0;n<dim;n++){fprintf(file," %10.3i", (int)pairs->lst[i]->key2a->descr[n]);}       
        fprintf(file,"     %10i%10i ",pairs->lst[i]->key2a->o, pairs->lst[i]->key2a->s);
        for(int n=0;n<n_bins;n++){fprintf(file," %10.3f", pairs->lst[i]->key2a->hist_prOri[n]);}
        
        /* x2b , y2b , sigma2b, theta2b */
        fprintf(file,"%20.5f%20.5f%20.5f%20.5f", pairs->lst[i]->key2b->x,
                                                 pairs->lst[i]->key2b->y,
                                                 pairs->lst[i]->key2b->sigma,
                                                 pairs->lst[i]->key2b->theta);
        for(int n=0;n<dim;n++){fprintf(file," %10.3i", (int)pairs->lst[i]->key2b->descr[n]);}       
        fprintf(file,"  %10i %10i ",pairs->lst[i]->key2b->o, pairs->lst[i]->key2b->s);
        for(int n=0;n<n_bins;n++){fprintf(file," %10.3f", pairs->lst[i]->key2b->hist_prOri[n]);}
        
        fprintf(file,"\n");
    }
    fclose(file);
}


void print_pairs(struct kyptPr_lst* pairs){
    for(int i=0; i<pairs->card;i++){   
         printf("%20.5f%20.5f%20.5f%20.5f%20.5f%20.5f%20.5f%20.5f\n", pairs->lst[i]->key1->x,
                                                                      pairs->lst[i]->key1->y,
                                                                      pairs->lst[i]->key1->sigma,
                                                                      pairs->lst[i]->key1->theta,
                                                                      pairs->lst[i]->key2a->x,
                                                                      pairs->lst[i]->key2a->y,
                                                                      pairs->lst[i]->key2a->sigma,
                                                                      pairs->lst[i]->key2a->theta);
    }
}
