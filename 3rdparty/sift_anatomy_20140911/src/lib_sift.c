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


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.



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

*/
/**
 * @file lib_sift.c 
 * @brief A simplified interface to the anatomy with standard parameters.
 *
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lib_sift.h"
#include "lib_sift_anatomy.h"
#include "lib_keypoint.h"
#include "lib_util.h"


static struct sift_keypoints* sift_translate_standard_into_anatomy(const struct sift_keypoint_std* k, int n)
{
    // load the default parameters are required
    struct sift_parameters* p = sift_assign_default_parameters();
    float sigma_min = p->sigma_min;  // 0.8
    float delta_min = p->delta_min;  // 0.5
    int n_spo = p->n_spo;   // 3
    int n_ori = p->n_ori;   // 8
    int n_hist = p->n_hist; // 4
    int n_bins = p->n_bins; // 36

    struct sift_keypoints* keys = sift_malloc_keypoints();
    for(int i = 0; i < n; i++){
        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
        /* reading the extremum continuous coordinates */
        key->x = k[i].x;
        key->y = k[i].y;
        key->sigma = k[i].scale;
        key->theta = k[i].orientation;
        /*  inferring the discrete coordinates in the scale-space grid */
        // We look for the pair of integer (o,s) such that nspo*o+s is the nearest
        // to alpha = nspo * log( k[i].scale / sigma_min) /M_LN2
        // with the constraint s=1,2,..,nspo.
        int o,s;
        int a = (int)(round( n_spo * log( k[i].scale / sigma_min) /M_LN2  ));
        o = (a-1)/n_spo;
        if (o > -1){
            s = (a-1)%n_spo + 1;
        }
        else{
            o = 0;
            s = 0;
        }
        key->o = o;
        key->s = s;
        key->i = (int)( key->x / ( delta_min * exp( key->o * M_LN2)) + 0.5 );
        key->j = (int)( key->y / ( delta_min * exp( key->o * M_LN2)) + 0.5 );
        sift_add_keypoint_to_list(key,keys);
    }
    return keys;
}


//static void sift_translate_anatomy_into_standard(const struct sift_keypoints *keys, struct sift_keypoint_std* k, int *n)
void sift_translate_anatomy_into_standard(const struct sift_keypoints *keys, struct sift_keypoint_std* k, int *n)
{
    *n = keys->size;
    k = (struct sift_keypoint_std*)xrealloc(k, (*n)*sizeof(*k));
    for(int i = 0; i < *n; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;
        for(int j = 0; j < 128; j++){
            k[i].descriptor[j] = keys->list[i]->descr[j];
        }
    }
}





/** @brief Extracts oriented keypoints (without description
 *
 *
 */
struct sift_keypoint_std* sift_compute_features(const float* x, int w, int h, int *n)
{

    /** assign default parameters **/
    struct sift_parameters* p = sift_assign_default_parameters();

    /** Memory dynamic allocation */
    // WARNING 6 lists of keypoints containing intermediary states of the algorithm
    struct sift_keypoints **kk = xmalloc(6*sizeof(struct sift_keypoints*));
    for(int i = 0; i < 6; i++){
        kk[i] = sift_malloc_keypoints();
    }
    // WARNING 4 scalespace structures
    struct sift_scalespace **ss = xmalloc(4*sizeof(struct sift_scalespace*));

    /** Algorithm */
    struct sift_keypoints* keys = sift_anatomy(x, w, h, p, ss, kk);

    /* Copy to a list of keypoints */
    *n = keys->size;
    struct sift_keypoint_std* k = xmalloc((*n)*sizeof(*k));
    for(int i = 0; i < keys->size; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;
        for(int j = 0; j < 128; j++){
            k[i].descriptor[j] = keys->list[i]->descr[j];
        }
    }

    /* memory deallocation */
    xfree(p);
    sift_free_keypoints(keys);
    for(int i = 0; i < 6; i++){
        sift_free_keypoints(kk[i]);
    }
    xfree(kk);
    for(int i = 0; i < 4; i++){
        sift_free_scalespace(ss[i]);
    }
    xfree(ss);

    return k;
}



/** @brief Extracts oriented keypoints (without description
 *
 *
 */
struct sift_keypoint_std* sift_compute_points(const float* x, int w, int h, int *n)
{

    /** assign default parameters **/
    struct sift_parameters* p = sift_assign_default_parameters();

    /** Memory dynamic allocation */
    // WARNING 5 lists of keypoints containing intermediary states of the algorithm
    struct sift_keypoints **kk = xmalloc(5*sizeof(struct sift_keypoints*));
    for(int i = 0; i < 5; i++){
        kk[i] = sift_malloc_keypoints();
    }
    // WARNING 4 scalespace structure containing the DoG and Gaussian scalespace and the gradient
    struct sift_scalespace **ss = xmalloc(2*sizeof(struct sift_scalespace*));

    /** Algorithm */
    struct sift_keypoints* keys = sift_anatomy_without_description(x, w, h, p, ss, kk);

    /* Copy to a list of keypoints */
    *n = keys->size;
    struct sift_keypoint_std* k = xmalloc((*n)*sizeof(*k));
    for(int i = 0; i < keys->size; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;
        for(int j = 0; j < 128; j++){
            k[i].descriptor[j] = 0;
        }
    }

    /* memory deallocation */
    xfree(p);
    sift_free_keypoints(keys);
    for(int i = 0; i < 5; i++){
        sift_free_keypoints(kk[i]);
    }
    xfree(kk);
    for(int i = 0; i < 2; i++){
        sift_free_scalespace(ss[i]);
    }
    xfree(ss);

    return k;
}

/** @brief Computes SIFT descriptors on given oriented keypoints
 *
 */
void sift_fill_descriptors(const float *x, int w, int h, struct sift_keypoint_std *k, int n)
{

    struct sift_keypoints* keys = sift_translate_standard_into_anatomy(k, n);

    /** assign default parameters **/
    struct sift_parameters* p = sift_assign_default_parameters();

    /** algorithm */
    sift_anatomy_only_description(x, w, h, p, keys);

    /* Copy back to the input flat list of keypoints */
    for(int i = 0; i < keys->size; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;
        for(int j=0; j<128; j++){
            k[i].descriptor[j] = (unsigned char)(keys->list[i]->descr[j]);
        }
    }
}


void sift_find_ori_and_fill_descriptors(const float *x, int w, int h, struct sift_keypoint_std *k, int n)
{
    struct sift_keypoints* keys = sift_translate_standard_into_anatomy(k, n);

    /** assign default parameters **/
    struct sift_parameters* p = sift_assign_default_parameters();

    /** algorithm */
    sift_anatomy_orientation_and_description(x, w, h, p, keys);

    /* Copy back to the input flat list of keypoints */
    for(int i = 0; i < keys->size; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;
        for(int j=0; j<128; j++){
            k[i].descriptor[j] = (unsigned char)(keys->list[i]->descr[j]);
        }
    }
}

void fprintf_keypoint_std(FILE* f, const struct sift_keypoint_std* k, int n)
{
    for(int i=0; i<n; i++){
        fprintf(f, "%f %f %f %f ", k[i].x, k[i].y, k[i].scale, k[i].orientation);
        for(int j=0; j<128; j++){
            fprintf(f, "%u ", k[i].descriptor[j]);
        }
        fprintf(f, "\n");
    }
}


void sift_write_to_file(const char *filename, const struct sift_keypoint_std *k, int n)
{
    FILE* file = fopen(filename,"w");
    fprintf_keypoint_std(file, k, n);
    fclose(file);
}



struct sift_keypoint_std * sift_read_from_file(const char *filename, int *n)
{
    struct sift_parameters* p = sift_assign_default_parameters();
    int n_ori = p->n_ori;   // 8
    int n_hist = p->n_hist; // 4
    int l = n_hist * n_hist * n_ori;
    int n_bins = p->n_bins;

    struct sift_keypoints* keys = sift_malloc_keypoints();
    int flag = 1; // read coordinates + descriptor
    sift_read_keypoints(keys, filename, n_hist, n_ori, n_bins, flag);

    // translating into a flat list
    *n = keys->size;
    struct sift_keypoint_std* k = xmalloc((*n)*sizeof(*k));
    for(int i = 0; i < keys->size; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;

        for(int j = 0; j < l; j++){
            k[i].descriptor[j] = keys->list[i]->descr[j];
        }
    }
    sift_free_keypoints(keys);
    return k;
}


struct sift_keypoint_std *sift_read_keyslocation_from_file(char *filename, int *n)
{
    struct sift_parameters* p = sift_assign_default_parameters();
    int n_ori = p->n_ori;   // 8
    int n_hist = p->n_hist; // 4
    int l = n_hist*n_hist*n_ori;
    int n_bins = p->n_bins;

    // read keypoints locations from a file and
    // save them into a sift_keypoints structure
    struct sift_keypoints* keys = sift_malloc_keypoints();
    int flag = 0; // read coordinates 
    sift_read_keypoints(keys, filename, n_hist, n_ori, n_bins, flag);

    // translate the sift_keypoints structure into a flat list
    *n = keys->size;
    struct sift_keypoint_std* k = xmalloc((*n)*sizeof(struct sift_keypoint_std));
    for(int i = 0; i < keys->size; i++){
        k[i].x = keys->list[i]->x;
        k[i].y = keys->list[i]->y;
        k[i].scale = keys->list[i]->sigma;
        k[i].orientation = keys->list[i]->theta;  // 0
        for(int j=0; j<l; j++){
            k[i].descriptor[j] = keys->list[i]->descr[j]; // 0
        }
    }
    sift_free_keypoints(keys);

    return k;
}
