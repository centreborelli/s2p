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


*/
/**
 * @file lib_sift.h
 * @brief The SIFT algorithmic chain
 *
 *
 * @author Ives Rey-Otero <ivesreyotero@gmail.com>
 */
/**
 * @file lib_sift.c 
 * @brief A simplified interface to the anatomy with standard parameters.
 *
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */



#ifndef _LIB_SIFT_H_
#define _LIB_SIFT_H_
#include <stdio.h>



/** @brief the SIFT keypoint structure
 *
 * contains the position and orientation of the keypoint and the 128 descriptor from the original SIFT method.
 *
 */
struct sift_keypoint_std
{
    float x;
    float y;
    float scale;
    float orientation;
    unsigned char descriptor[128];
};


/** @brief: compute the keypoints of an image (position, scale and orientation)
 *    Does not extract the SIFT descriptor
 *
 *  x input image of dimensions w \times h
 *
 *  n points to the integer that stores the number of extracted keypoints.
 *
 */
struct sift_keypoint_std* sift_compute_points(const float* x, int w, int h, int* n);


/** @brief: compute the descriptors of a given set of oriented keypoints
 *
 * The oriented keypoints are provided by the user as a flat list of keypoint structure.
 *
 */
void sift_fill_descriptors(const float *x, int w, int h, struct sift_keypoint_std *k, int n);

/** @brief: compute the descriptors of a given set of oriented keypoints
 *
 * The keypoints are provided by the user as a flat list of keypoint structure.
 * The routine computes one principal orientation and one descriptor per keypoint.
 *
 */
void sift_find_ori_and_fill_descriptors(const float *x, int w, int h, struct sift_keypoint_std *k, int n);

/** @brief: compute keypoints and their descriptors (The standard SIFT method)
 *
 */
struct sift_keypoint_std *sift_compute_features(const float *x, int w, int h, int *n);

// input and output of SIFT features
struct sift_keypoint_std *sift_read_from_file(const char *filename, int *n);


// Read SIFT keypoints location (x,y,sigma) from a file
struct sift_keypoint_std *sift_read_keyslocation_from_file(char *filename, int *n);


void sift_write_to_file(const char *filename, const struct sift_keypoint_std *k, int n);
void fprintf_keypoint_std(FILE* f, const struct sift_keypoint_std* k, int n);

#endif
