/**
 * @file lib_matching.h
 * @brief data structures to store information relative to a pair of keypoints
 *
 * @li struct keypointPr     : Pair of keypoint data structure.
 * @li struct keypointPr_list : List of pairs.
 * @li print,save, read for lists of pairs.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */
#ifndef _LIB_MATCHING_H_
#define _LIB_MATCHING_H_

void matching(struct sift_keypoints *k1, struct sift_keypoints *k2,
              struct sift_keypoints *out_k1, struct sift_keypoints *out_k2,
              float sift_thresh, int flag,
              double fund_mat[5], float epi_thresh);

void fprintf_pairs(const char *filename,
                   const struct sift_keypoints *k1,
                   const struct sift_keypoints *k2);

void print_pairs(const struct sift_keypoints *k1,
                 const struct sift_keypoints *k2);

void save_pairs_extra(const char* name,
                      const struct sift_keypoints *k1,
                      const struct sift_keypoints *k2A,
                      const struct sift_keypoints *k2B);

#endif
