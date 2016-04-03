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
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../timing.h"
#include "../linalg.h"
#include "lib_keypoint.h"
#include "lib_matching.h"
#include "lib_util.h"

static void compute_keypoints_distance(float* dist,
                                       const struct sift_keypoints *k1,
                                       const struct sift_keypoints *k2)
{
    int n_hist = k1->list[0]->n_hist;
    int n_ori  = k1->list[0]->n_ori;
    int n1 = k1->size;
    int n2 = k2->size;
    for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
        dist[i*n2+j] = euclidean_distance(k1->list[i]->descr,
                                          k2->list[j]->descr,
                                          n_hist * n_hist * n_ori);
}


static void compute_keypoints_distance_epipolar(float* dist,
                                                const struct sift_keypoints *k1,
                                                const struct sift_keypoints *k2,
                                                const float s1[3],
                                                const float s2[3],
                                                float epi_thresh)
{
    int n = k1->list[0]->n_hist * k1->list[0]->n_hist * k1->list[0]->n_ori;
    int n1 = k1->size;
    int n2 = k2->size;
    float a = s1[0];
    float b = s1[1];
    float c = s1[2];
    float d = s2[0];
    float e = s2[1];
    float f = s2[2];
    for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
        // keypoints x, y coordinates
        float x1 = k1->list[i]->x;
        float y1 = k1->list[i]->y;
        float x2 = k2->list[j]->x;
        float y2 = k2->list[j]->y;

        // rectified y coordinates
        float yy1 = b * x1 + a * y1 + c; // scalar product of (b, a, c) and (x1, y1, 1)
        float yy2 = e * x2 + d * y2 + f; // scalar product of (e, d, f) and (x2, y2, 1)

        // points satisfy the epipolar constraint when the rectified y are equal
        if (fabs(yy1 - yy2) < epi_thresh)
            dist[i * n2 + j] = euclidean_distance(k1->list[i]->descr,
                                                  k2->list[j]->descr, n);
    }
}


static void find_the_two_nearest_keys(const float* dist, int n1, int n2,
                                      int* indexA, int* indexB,
                                      float* distA, float* distB)
{
    for(int i = 0; i < n1; i++){
        int iA, iB;
        float dA, dB;
        find_array_two_min(&dist[i*n2], n2, &dA, &dB, &iA, &iB);
        indexA[i] = iA;
        indexB[i] = iB;
        distA[i] = dA;
        distB[i] = dB;
    }
}


void matching(struct sift_keypoints *k1, struct sift_keypoints *k2,
              struct sift_keypoints *out_k1, struct sift_keypoints *out_k2,
              float sift_thresh, int flag,
              double fund_mat[5], float epi_thresh, bool verbose)
{
    int n1 = k1->size;
    int n2 = k2->size;

    struct timespec ts; portable_gettime(&ts);
    float *dist  = (float *) xmalloc(n1 * n2 * sizeof(float));
    float *distA = (float *) xmalloc(n1 * sizeof(float));
    float *distB = (float *) xmalloc(n1 * sizeof(float));
    int *indexA  = (int *) xmalloc(n1 * sizeof(int));
    int *indexB  = (int *) xmalloc(n1 * sizeof(int));

    if (fund_mat) {
        float s1[3]; float s2[3];
        rectifying_similarities_from_affine_fundamental_matrix(s1, s2, fund_mat);
        //fprintf(stderr, "s1: %f %f %f\n", s1[0], s1[1], s1[2]);
        //fprintf(stderr, "s2: %f %f %f\n", s2[0], s2[1], s2[2]);
        for (int i = 0; i < n1 * n2; i++)
            dist[i] = INFINITY;
        if (verbose) print_elapsed_time(&ts, " - init distances:", 33);
        compute_keypoints_distance_epipolar(dist, k1, k2, s1, s2, epi_thresh);
    } else
        compute_keypoints_distance(dist, k1, k2);
    if (verbose) print_elapsed_time(&ts, " - compute distances:", 33);

    find_the_two_nearest_keys(dist, n1, n2, indexA, indexB, distA, distB);
    if (verbose) print_elapsed_time(&ts, " - find two nearest:", 33);

    for (int i = 0; i < n1; i++) {
        float val = (flag == 1 ? distA[i] / distB[i] : distA[i]);
        if (val < sift_thresh) {
            sift_add_keypoint_to_list(k1->list[i], out_k1);
            sift_add_keypoint_to_list(k2->list[indexA[i]], out_k2);
        }
    }
    if (verbose) print_elapsed_time(&ts, " - select matches:", 33);

    free(dist);
    free(indexA);
    free(indexB);
    free(distA);
    free(distB);
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
