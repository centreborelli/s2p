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
 * @file matching.c
 * @brief [[MAIN]] Matching of keypoints
 *
 * @li Method 1 : threshold on the ratio of distances to the two nearest keypoints
 * @li Method 2 : threshold on the distance to the nearest keypoint
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lib_keypoint.h"
#include "lib_matching.h"
#include "lib_util.h"


void print_usage()
{
    fprintf(stderr, "Anatomy of the SIFT method (www.ipol.im/pub/pre/82/)  ver 20140911         \n");
    fprintf(stderr, "Usage:  match_cli keys1 keys2 [options...]                                 \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "    -ori_nbins    (36)  number of bins in the orientation histogram        \n");
    fprintf(stderr, "                        (used only for keypoints input/output)             \n");
    fprintf(stderr, "    -descr_nhist   (4)  number of histograms per dimension                 \n");
    fprintf(stderr, "    -descr_nori    (8)  number of bins in each histogram                   \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "    -absolute thresh (250) threshold applied on the euclidean distance     \n");
    fprintf(stderr, "    -relative thresh (0.6) threshold applied on the ratio of  distance     \n");
    fprintf(stderr, "                                                                           \n");
    fprintf(stderr, "    -verb         label  flag for output                                   \n");
}


/**
 *
 * Output 
 *   -1 : malformed argument
 *    0 : option not found  
 *    1 : option found
 */
static int pick_option(int* c, char*** v, char* opt, char* val)
{
    int output = 0;
    int argc = *c;
    char **argv = *v;
    // scan the command line for '-opt'
    for(int i = 0; i < argc; i++){
        if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, opt))
        {
            // check for a corresponding value
            if (i == argc-1){
                output = -1;
            }
            else{
                if (argv[i+1][0] == '-'){
                    output  = -1;
                }
                // if the option call is well formed
                else{
                    // copy the option value ...
                    strcpy(val, argv[i+1]);
                    // ... and remove from the command line
                    for (int j = i; j < argc - 2; j++){
                        (*v)[j] = (*v)[j+2];
                    }
                    *c -= 2;
                    output = 1;
                }
            }
            // print an error if not succes
            if(output == -1){
                fprintf(stderr, "Fatal error: option %s requires an argument.\n", opt);
            }
        }
    }
    return output;
}

static int parse_options(int argc, char** argv,
                         int *n_bins,
                         int *n_hist,
                         int *n_ori,
                         int *meth_flag,
                         float *thresh,
                         int *verb_flag,
                         char* label)
{
    int isfound;
    char val[128];

    isfound = pick_option(&argc, &argv, "ori_nbins", val);
    if (isfound ==  1)    *n_bins = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "descr_nhist", val);
    if (isfound ==  1)    *n_hist = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "descr_nori", val);
    if (isfound ==  1)    *n_ori = atoi(val);
    if (isfound == -1)    return EXIT_FAILURE;


    isfound = pick_option(&argc, &argv, "absolute", val);
    if (isfound ==  1){
        *meth_flag = 0;
        *thresh = atof(val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "relative", val);
    if (isfound ==  1){
        *meth_flag = 1;
        *thresh = atof(val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    isfound = pick_option(&argc, &argv, "verb", val);
    if (isfound ==  1){
        *verb_flag = 1;
        strcpy(label, val);
    }
    if (isfound == -1)    return EXIT_FAILURE;

    // check for unknown option call
    for(int i = 0; i < argc; i++){
        if (argv[i][0] == '-'){
            fprintf(stderr, "Fatal error: option \"-%s\" is unknown.\n", argv[i]+1);
            print_usage();
            return EXIT_FAILURE;
        }
    }
    // check for input image
    if (argc != 3){
        print_usage();
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}



int main(int argc, char **argv)
{
    // Setting default parameters
    int n_hist = 4;
    int n_ori = 8;
    int n_bins = 36;
    int meth_flag = 1;
    float thresh = 0.6;
    int verb_flag = 0;
    char label[256];
    strcpy(label, "extra");

    // Parsing command line
    int res = parse_options(argc, argv, &n_bins, &n_hist, &n_ori,
                            &meth_flag, &thresh, &verb_flag, label);
    if (res == EXIT_FAILURE)
        return EXIT_FAILURE;

    // Memory allocation
    struct sift_keypoints* k1 = sift_malloc_keypoints();
    struct sift_keypoints* k2 = sift_malloc_keypoints();
    struct sift_keypoints* out_k1 = sift_malloc_keypoints();
    struct sift_keypoints* out_k2A = sift_malloc_keypoints();
    struct sift_keypoints* out_k2B = sift_malloc_keypoints();

    // Read input keypoint ASCII files
    int readflag = verb_flag + 1;
    sift_read_keypoints(k1, argv[1], n_hist, n_ori, n_bins, readflag);
    sift_read_keypoints(k2, argv[2], n_hist, n_ori, n_bins, readflag);

    // Matching
    matching(k1, k2, out_k1, out_k2A, out_k2B, thresh, meth_flag);

    // Print
    print_pairs(out_k1, out_k2A);
    char name[FILENAME_MAX];
    if(verb_flag == 1){
        save_pairs_extra("OUTmatches.txt", out_k1, out_k2A, out_k2B);
        sprintf(name, "%s_im0.txt", label);
        sift_save_keypoints(out_k1, name, 1);
        sprintf(name, "%s_im1.txt", label);
        sift_save_keypoints(out_k2A, name, 1);
    }

    // Free memory
    sift_free_keypoints(k1);
    sift_free_keypoints(k2);
    sift_free_keypoints(out_k1);
    sift_free_keypoints(out_k2A);
    sift_free_keypoints(out_k2B);

    return EXIT_SUCCESS;
}
