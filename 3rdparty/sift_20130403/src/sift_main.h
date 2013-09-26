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




*/
/**
 * @file sift_main.h
 * @brief The SIFT algorithmic chain
 *
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#ifndef _SIFT_MAIN_H_
#define _SIFT_MAIN_H_

#include <stdio.h>
#include "sift_keypoint.h"
#include "sift_scalespace.h"
#include "sift_description.h"
// #include <iostream>
#include <time.h>
#include <stdbool.h>
#include <assert.h>

#include <omp.h>



void sift_transform(double* image, int width, int height,
                    struct kypt_lst* keysOUT,
                    char* basename,
                    int flag,
                    int n_oct,
                    int n_spo,
                    double sigma_min,
                    double delta_min,
                    double sigma_in,
                    double C_DoG,
                    double C_edge,
                    int n_bins,
                    double lambda_ori,
                    double t,
                    double n_hist,
                    double n_ori,
                    double lambda_descr);

void build_scalespace(double* im,
                      int im_width, int im_height, double im_sigma,
                      struct scsp* scalespace);

void build_lowe_scalespace(struct scsp* scalespace,
                           double* image, int im_width, int im_height,
                           double sigma0);

void build_lowe_scalespace_ter(struct scsp* scalespace,
                               double* image,
                               int im_width,
                               int im_height,
                               double sigma_in);




/** Build Lowe's scale space with a trick (using subsampling of already computed images */
void build_lowe_scalespace_bis(struct scsp* scalespace,
                               double* image, int im_width, int im_height,
                               double sigma0);

void build_dog(struct scsp *scalespace,
               struct scsp *dog);

void find_candidate_keypoints(struct scsp *dog,
                              struct kypt_lst*keys,
                              int n_ori,   /* required for the allocation off the keypoint structures */
                              int n_hist,
                              int n_bins);

void interpolate_position_of_candidate_keypoints(struct scsp *dog,
                                                 struct kypt_lst *keys,
                                                 struct kypt_lst *keysInterpol,
                                                 struct kypt_lst *keysReject);

void build_scalespace_gradient(struct scsp *scalespace,
                               struct scsp *dx_scalespace,
                               struct scsp *dy_scalespace);



void discard_detections_duetonoise(struct kypt_lst *keys_in,
                                   struct kypt_lst *keys_accept,
                                   struct kypt_lst *keys_reject,
                                   double thresh);

void compute_HarrisStephen_edge_response(struct scsp *dog,
                                         struct kypt_lst *keys);

void discard_detections_onedge(struct kypt_lst *keysIn,
                               struct kypt_lst *keysAccept,
                               struct kypt_lst *keysReject,
                               double thresh);
                                     
void attribute_orientation_to_keypoints(struct scsp *dx_scalespace,
                                        struct scsp *dy_scalespace, /*gradient scalespaces*/
                                        struct kypt_lst *keysIn,
                                        struct kypt_lst *keysOut,
                                        struct kypt_lst *cp_keys_multiori, /* a copy of keys with multiple detections */
                                        int n_bins,
                                        double lambda_ori,
                                        double t); 

void attribute_featurevector_to_keypoints(struct scsp *dx_scalespace,
                                          struct scsp *dy_scalespace, /*gradient scalespaces*/
                                          struct kypt_lst *keys,
                                          int n_hist,
                                          int n_ori,
                                          double lambda_descr);


/** 10 complex parameters !!
 * maybe we should SPLIT this in half , first "prepare" the scalespace structure,
 * by prepare I mean :
 * 1) allocate the memory using the parameters  
 * 2) change attributes using again the parameters
 *        (maybe 1 and 2 is already covered by malloc..)
 * and there is no other computation otherwise 
we should avoid too many parameters, or only for very simple things such as allocation
EDIT this is done via malloc_scsp or associated wrapper EDIT
 */




#endif
