/**
 * @file sift_matching.c
 * @brief data structures to store information relative to a pair of keypoints
 *
 * @li struct keypointPr     : Pair of keypoint data structure.
 * @li struct keypointPr_list : List of pairs.
 * @li print,save, read for lists of pairs.
 *
 * @author (original) Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 * @author (modified) Carlo de Franchis <carlodef@gmail.com>
 * @author (modified) David Youssefi <david.youssefi@c-s.fr>
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../linalg.h"
#include "lib_keypoint.h"
#include "lib_matching.h"
#include "lib_util.h"

static float keypoints_distance_epipolar(struct keypoint* k1_i,
					 struct keypoint* k2_j,
					 float s1[3],
					 float s2[3],
					 float epi_thresh,
					 int n)
{
    float dist = INFINITY;
    float a = s1[0];
    float b = s1[1];
    float c = s1[2];
    float d = s2[0];
    float e = s2[1];
    float f = s2[2]; 

    // keypoints x, y coordinates
    float x1 = k1_i->x;
    float y1 = k1_i->y;
    float x2 = k2_j->x;
    float y2 = k2_j->y;

    // /!\ in Ives' conventions x is the row index
    // rectified x coordinates
    float xx1 = b * y1 + a * x1 + c; // scalar product of (b, a, c) and (x1, y1, 1)
    float xx2 = e * y2 + d * x2 + f; // scalar product of (e, d, f) and (x2, y2, 1)

    // points satisfy the epipolar constraint when the rectified x are equal
    if (fabs(xx1 - xx2) < epi_thresh){
      dist = euclidean_distance(k1_i->descr,
				k2_j->descr, n);
    }
    return dist;
}

static float keypoints_distance(struct keypoint* k1_i,
				struct keypoint* k2_j,
				float s1[3],
				float s2[3],
				float epi_thresh,
				int n)
{
    // parameters used in keypoints_distance_epipolar
    (void) s1;
    (void) s2;
    (void) epi_thresh;

    return euclidean_distance(k1_i->descr,
			      k2_j->descr, n);
}



void matching(struct sift_keypoints *k1, struct sift_keypoints *k2,
	      struct sift_keypoints *out_k1, struct sift_keypoints *out_k2,
	      float sift_thresh, int flag,
	      double fund_mat[5], float epi_thresh)
{
    int n, n1, n2, n_hist, n_ori;
    float dist, distA, distB;
    float s1[3]; float s2[3];
    int indexA;
    float val;

    float (*keypoints_distance_overloaded)(struct keypoint*,
					   struct keypoint*,
					   float*, float*,
					   float, int);

    if (fund_mat) {
      rectifying_similarities_from_affine_fundamental_matrix(s1, s2, fund_mat);
      keypoints_distance_overloaded = keypoints_distance_epipolar;
    }
    else {
      keypoints_distance_overloaded = keypoints_distance;
    }

    n1 = k1->size;
    n2 = k2->size;  
    n_hist = k1->list[0]->n_hist;
    n_ori = k1->list[0]->n_ori;
    n = n_hist * n_hist * n_ori;

    for (int i = 0; i < n1; i++) {
      distA = INFINITY;
      distB = INFINITY;
      indexA = -1;

      for (int j = 0; j < n2; j++) {
	dist = keypoints_distance_overloaded(k1->list[i],
					     k2->list[j],
					     s1,s2,
					     epi_thresh,n);

	// find_the_two_nearest_keys
	if (dist < distA) {
	  distB = distA;
	  distA = dist;
	  indexA = j;
	}
	else if (dist < distB) {
	  distB = dist;
	}
      }

      val = (flag == 1 ? distA / distB : distA);
      if (val < sift_thresh) {
	sift_add_keypoint_to_list(k1->list[i], out_k1);
	sift_add_keypoint_to_list(k2->list[indexA], out_k2);
      }
    }
}

void fprintf_pairs(const char *filename, const struct sift_keypoints *k1,
                   const struct sift_keypoints *k2)
{
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < k1->size; i++) {
        fprintf_one_keypoint(f, k1->list[i], 0, 0, -1);
        fprintf_one_keypoint(f, k2->list[i], 0, 0, -1);
        fprintf(f, "\n");
    }
    fclose(f);
}


void print_pairs(const struct sift_keypoints *k1,
                 const struct sift_keypoints *k2)
{
    fprintf_pairs("stdout", k1, k2);
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
