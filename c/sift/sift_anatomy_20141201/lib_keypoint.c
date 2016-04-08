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
 @file sift_keypoint.c
 * @brief data structures to store information relative to keypoint
 *
 * @li struct keypoint      keypoint data structure.
 * @li struct sift_keypoints  list of keypoints with variable length.
 * @li print,save, read  for lists of keypoints.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "lib_keypoint.h"
#include "lib_util.h"


struct keypoint* sift_malloc_keypoint(int n_ori, int n_hist, int n_bins)
{
    struct keypoint* key =  xmalloc(sizeof(struct keypoint));
    key->n_ori  = n_ori;
    key->n_hist = n_hist;
    key->n_bins = n_bins;
    key->descr = xmalloc(n_ori*n_hist*n_hist*sizeof(float));
    key->orihist = xmalloc(n_bins*sizeof(float));
    // set default value
    key->x = 0.;
    key->y = 0.;
    key->sigma = 0.;
    key->theta = 0.;
    key->val = 0;
    key->o = 0;
    key->s = 0;
    key->i = 0;
    key->j = 0;
    for(int n=0;n<n_ori*n_hist*n_hist;n++){
        key->descr[n] = 0.;
    }
    for(int i = 0; i < n_bins; i++){
        key->orihist[i] = 0.;
    }
    // return pointer
    return key;
}


static void copy_keypoint(const struct keypoint* kA, struct keypoint* kB)
{
    int l = kB->n_hist * kB->n_hist * kB->n_ori; // length of the feature vector
    assert(l == (kA->n_hist * kA->n_hist * kA->n_ori));
    // copy struct
    kB->i        = kA->i;
    kB->j        = kA->j;
    kB->s        = kA->s;
    kB->o        = kA->o;
    kB->x        = kA->x;
    kB->y        = kA->y;
    kB->sigma    = kA->sigma;
    kB->val      = kA->val;
    kB->theta    = kA->theta;
    kB->edgeResp = kA->edgeResp;
    for(int n = 0; n < l; n++){
        kB->descr[n] = kA->descr[n];
    }
    for(int n = 0; n < kB->n_bins; n++){
        kB->orihist[n] = kA->orihist[n];
    }
}

struct keypoint* sift_malloc_keypoint_from_model_and_copy(const struct keypoint* key)
{
    int n_hist  = key->n_hist;   /* number of histograms per direction in the feature vector */
    int n_ori   = key->n_ori;    /* number of bins in each orientation histogram */
    int n_bins  = key->n_bins;   /* number of bins in the orientation histogram for orientation attribution */
    struct keypoint* keycp = sift_malloc_keypoint(n_ori,n_hist,n_bins);
    copy_keypoint(key, keycp);
    return keycp;
}


struct sift_keypoints* sift_malloc_keypoints()
{
    struct sift_keypoints *keys = xmalloc(sizeof(struct sift_keypoints));
    keys->list = xmalloc(100 * sizeof(struct keypoint*));
    keys->capacity = 100;
    keys->size = 0;
    return keys;
}


static void realloc_sift_keypoints(struct sift_keypoints* keys)
{
    keys->list = (struct keypoint **) xrealloc(keys->list,
                                               2 * keys->capacity * sizeof(struct keypoint*));
    keys->capacity *= 2;
}

void sift_add_keypoint_to_list(struct keypoint* key,
                               struct sift_keypoints* keys)
{
    if (keys->size > keys->capacity - 1) // checks the memory limit
        realloc_sift_keypoints(keys);
    keys->list[keys->size] = key;
    keys->size += 1;
}


void sift_free_keypoint(struct keypoint* key)
{
    xfree(key->orihist);
    xfree(key->descr);
    xfree(key);
}


void sift_free_keypoints(struct sift_keypoints* keys)
{
    for (int k = 0; k < keys->size; k++)
        sift_free_keypoint(keys->list[k]);
    xfree(keys->list);
    xfree(keys);
}

void fprintf_one_keypoint(FILE* f, const struct keypoint* k, int n_descr, int n_bins, int flag)
{
    // coordinates
    fprintf(f,"%f %f ", k->y, k->x);

    // scale and orientation
    if (flag > -1)
        fprintf(f,"%f %f ", k->sigma, k->theta);

    // descriptor
    if (flag > 0)
        for (int n = 0; n < n_descr; n++)
            fprintf(f,"%i ", (int) k->descr[n]);

    // orientation histogram
    if (flag > 1)
        for (int n = 0; n < n_bins; n++)
            fprintf(f,"%f ", k->orihist[n]);
}


// dump coordinates, scale, orientation and descriptor to a binary file
void raw_dump_keypoint(FILE* f, const struct keypoint* k, int n_descr)
{
    char buf[4 * sizeof(float) + n_descr * sizeof(uint8_t)];
    float *keypoint = (float *) buf;
    uint8_t *descriptor = (uint8_t *) (buf + 4 * sizeof(float)) ;

    // copy the data to a buffer
    keypoint[0] = k->y;
    keypoint[1] = k->x;
    keypoint[2] = k->sigma;
    keypoint[3] = k->theta;
    for (int j = 0; j < n_descr; j++)
        descriptor[j] = (k->descr)[j];

    // write the buffer to the output file
    fwrite(buf, sizeof(char), 4 * sizeof(float) + n_descr, f);
}


void fprintf_keypoints(FILE* f, const struct sift_keypoints* keys, int nb,
                       bool binary, int flag)
{
    // use at most nb keypoints unless nb is zero
    int n = nb ? min_2(keys->size, nb) : keys->size;
    if (n > 0) {
        int n_hist = keys->list[0]->n_hist;
        int n_ori  = keys->list[0]->n_ori;
        int n_descr = n_hist * n_hist * n_ori;
        int n_bins = keys->list[0]->n_bins;
        for (int i = 0; i < n; i++) {
            if (binary)
                raw_dump_keypoint(f, keys->list[i], n_descr);
            else {
                fprintf_one_keypoint(f, keys->list[i], n_descr, n_bins, flag);
                fprintf(f,"\n");
            }
        }
    }
}


void sift_save_keypoints(const struct sift_keypoints* keys, const char* name, int flag)
{
    FILE* f = fopen(name,"w");
    if (!f){
        fatal_error("Failed to open %s for writing\n", name);
    }
    fprintf_keypoints(f, keys, 0, false, flag);
    fclose(f);
}


void sift_print_keypoints(const struct sift_keypoints* keys, int flag)
{
    fprintf_keypoints(stdout, keys, 0, false, flag);
}






/** @brief read list of oriented keypoints from txt file
 *
 * @param keys    =  output list of keypoints
 * @param name    =  input filename
 *
 * @param n_hist  the descriptor has (n_hist*n_hist) weighted coefficients.
 * @param n_ori   the descriptor has (n_hist*n_hist*n_ori) coefficients
 * @param n_bins  the size of the orientation histogram
 *
 */
void sift_read_keypoints(struct sift_keypoints* keys,
                         const char* name,
                         int n_hist,
                         int n_ori,
                         int n_bins,
                         int flag)
{
    size_t buffer_size = 1024 * 1024;  // 1MB buffer for long lines.
    char* buffer = xmalloc(buffer_size);
    FILE* stream = fopen(name,"r");
    if ( !stream)
        fatal_error("File \"%s\" not found.", name);
    while(fgets(buffer, buffer_size, stream) != NULL){
        int pos = 0;
        int read = 0;
        struct keypoint* key = sift_malloc_keypoint(n_ori, n_hist, n_bins);
        // read coordinates
        sscanf(buffer+pos,"%f  %f  %f  %f %n", &key->y
                                             , &key->x
                                             , &key->sigma
                                             , &key->theta
                                             , &read);
        pos+=read;
        if (flag > 0){
            // read descriptor
            for(int i = 0; i < n_hist*n_hist*n_ori; i++){
                sscanf(buffer+pos, "%f %n",&(key->descr[i]),&read);
                pos +=read;
            }
        }
        if (flag > 1){
            // read orientation histogram
            for(int i=0;i<n_bins;i++){
                sscanf(buffer+pos, "%f %n",&(key->orihist[i]),&read);
                pos += read;
            }
        }
        sift_add_keypoint_to_list(key,keys);
    }
    xfree(buffer);
}
