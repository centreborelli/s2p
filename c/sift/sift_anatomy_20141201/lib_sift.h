/**
 * @file lib_sift.h
 * @brief The SIFT algorithmic chain
 *
 *
 * @author Ives Rey-Otero <ivesreyotero@gmail.com>
 */
#ifndef _LIB_SIFT_H_
#define _LIB_SIFT_H_
#include <stdio.h>
#include <stdbool.h>



/** @brief the SIFT keypoint structure
 *
 * contains the position and orientation of the keypoint and the 128 descriptor
 * from the original SIFT method.
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
 *  n points to the integer that stores the number of extracted keypoints.
 */
struct sift_keypoint_std* sift_compute_points(const float* x, int w, int h,
        int* n);


/** @brief: compute the descriptors of a given set of oriented keypoints
 *
 * The oriented keypoints are provided by the user as a flat list of keypoint
 * structure.
 */
void sift_fill_descriptors(const float *x, int w, int h,
        struct sift_keypoint_std *k, int n);


/** @brief: compute the descriptors of a given set of oriented keypoints
 *
 * The keypoints are provided by the user as a flat list of keypoint structure.
 * The routine computes one principal orientation and one descriptor per keypoint.
 */
void sift_find_ori_and_fill_descriptors(const float *x, int w, int h,
        struct sift_keypoint_std *k, int n);


/** @brief: compute keypoints and their descriptors (The standard SIFT method)
 */
struct sift_keypoint_std *sift_compute_features(const float *x, int w, int h,
        int *n);

// input and output of SIFT features
struct sift_keypoint_std *sift_read_from_file(const char *filename, int *n);

// Read SIFT keypoints location (x, y, sigma) from a file
struct sift_keypoint_std *sift_read_keyslocation_from_file(char *filename, int *n);

void sift_write_to_file(const char *filename, const struct sift_keypoint_std *k, int n, bool binary);
void fprintf_keypoint_std(FILE* f, const struct sift_keypoint_std* k, int n);
#endif
