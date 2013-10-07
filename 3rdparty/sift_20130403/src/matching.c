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
 * @file matching.c
 * @brief [[MAIN]] Matching of keypoints
 *
 * @li Method 1 : threshold on the ratio of distances to the two nearest keypoints
 * @li Method 2 : threshold on the distance to the nearest keypoint
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */





#include "sift_main.h"
#include "sift_matching.h"
#include "io_png.h"




int main(int argc, char **argv){
    
    
    int n_hist;
    int n_ori;
    int n_bins;
    int flag;
    double C;
    int verbose;
    
    switch(argc){
        case 3:
            verbose = 0;
            flag   = 1;
            C      = 0.6; 
            n_hist = 4;
            n_ori  = 8;
            n_bins = 36;
            break;
        case 8:
            verbose = 1;
            flag    = atoi(argv[3]);
            C       = atof(argv[4]);
            n_hist  = atoi(argv[5]);
            n_ori   = atoi(argv[6]);
            n_bins  = atoi(argv[7]);
            break;
        default:
            printf(" Usage  matching keys1.txt keys2.txt methodflag(absolute=0,relative=1) Cte n_hist n_ori n_bins\n");
            printf("        matching keys1.txt keys2.txt 1 0.6 4 8 36 \n");
            return -1;
    }
        
    

    /* Memory allocation */
    struct kypt_lst* keys1 = malloc_kypt_lst();
    struct kypt_lst* keys2 = malloc_kypt_lst();
    struct kypt_lst* keys1Matching = malloc_kypt_lst();
    struct kypt_lst* keys2Matching = malloc_kypt_lst();
    struct kypt_lst* keys1NotMatch = malloc_kypt_lst();
    struct kyptPr_lst* pairs   = malloc_kyptPr_lst();
    struct kyptPr_lst* matches = malloc_kyptPr_lst();
    
    
    
        
    /** Read input keypoint ASCII files */
    if(verbose==1){
        /* File expected format
         * [ x y sigma theta fv[1]...fv[d] octa sca orihist[1]...orihist[n_bins] */
        read_kypt_lst_extra(keys1,argv[1],n_hist,n_ori,n_bins);
        read_kypt_lst_extra(keys2,argv[2],n_hist,n_ori,n_bins);
    }else{
        /* File expected format
         * [ x y sigma theta fv[1]...fv[128] */
        read_kypt_lst(keys1,argv[1],n_hist,n_ori);
        read_kypt_lst(keys2,argv[2],n_hist,n_ori);
    }
    
    /** Pairing and matching *************************************************/
    pairing_keyslist_to_keyslist(pairs, keys1, keys2, n_hist*n_hist*n_ori);
    
    only_matching_kyptPr(C,
                         flag,
                         pairs,
                         matches,
                         keys1Matching,
                         keys2Matching,
                         keys1NotMatch);
    
    /** OUTPUT ***************************************************************/
    print_pairs(matches);
    
//    /** EXTRA OUTPUT */
//    if(verbose==1){
//        save_kyptPr_extra(matches, n_hist*n_hist*n_ori, n_bins, "OUTmatches.txt");
//        save_kyptPr_extra(pairs, n_hist*n_hist*n_ori, n_bins, "OUTpairsAll.txt"); 
//        save_kypt_lst_extra(keys1Matching, n_hist*n_hist*n_ori,n_bins, "OUTkeys1Matching.txt");
//        save_kypt_lst_extra(keys2Matching, n_hist*n_hist*n_ori,n_bins, "OUTkeys2Matching.txt");
//        save_kypt_lst_extra(keys1NotMatch, n_hist*n_hist*n_ori,n_bins, "OUTkeys1NotMatching.txt");
//        // TODO -- AVOID duplicates
//        save_kypt_lst_description(keys1Matching, n_hist*n_hist*n_ori, "matching_keys_im0.txt");
//        save_kypt_lst_description(keys2Matching, n_hist*n_hist*n_ori, "matching_keys_im1.txt");        
//        
//    }else{
//        save_kypt_lst_description(keys1Matching, n_hist*n_hist*n_ori, "matching_keys_im0.txt");
//        save_kypt_lst_description(keys2Matching, n_hist*n_hist*n_ori, "matching_keys_im1.txt");
//    }
    
    
    
    /* Freeing memory */
    free_kypt_lst(keys1Matching);
    free_kypt_lst(keys2Matching);
    free_kypt_lst(keys1NotMatch);
    free_kypt_lst(keys1);
    free_kypt_lst(keys2);
    
    free_kyptPr_lst(pairs);  // it also includes the elements of match
    //free_kyptPr_lst(matches);
    
}
