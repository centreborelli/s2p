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

#include "timing.h"
#include "linalg.h"
#include "pickopt.h"
#include "sift_anatomy_20141201/lib_keypoint.h"
#include "sift_anatomy_20141201/lib_matching.h"


void print_usage(char *v[])
{
    fprintf(stderr, "Anatomy of the SIFT method (www.ipol.im/pub/pre/82/)  ver 20140911         \n");
    fprintf(stderr, "Usage:  %s keys1 keys2 [options...]                                  \n", v[0]);
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


int main(int c, char *v[])
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

    // parse arguments
    const char *output_file = pick_option(&c, &v, "o", "/dev/stdout");
    bool verbose = pick_option(&c, &v, "-verbose", NULL);

    // read the (optional) fundamental matrix
    const char *fund_mat_str = pick_option(&c, &v, "f", "");
    double *fund_mat = NULL;
    if (strcmp(fund_mat_str, "")) {
        // the string is not empty
        int n_fund;
        fund_mat = alloc_parse_doubles(5, fund_mat_str, &n_fund);
        if (n_fund != 5) {
            fprintf(stderr, "can't read affine fundamental matrix from \"%s\"\n", fund_mat_str);
            return EXIT_FAILURE;
        }
    }

    // initialise time counter
    struct timespec ts; portable_gettime(&ts);

    // Memory allocation
    struct sift_keypoints* k1 = sift_malloc_keypoints();
    struct sift_keypoints* k2 = sift_malloc_keypoints();
    struct sift_keypoints* out_k1 = sift_malloc_keypoints();
    struct sift_keypoints* out_k2 = sift_malloc_keypoints();
    struct sift_keypoints* out_k2B = sift_malloc_keypoints();

    // Read input keypoint ASCII files
    int readflag = verb_flag + 1;
    sift_read_keypoints(k1, v[1], n_hist, n_ori, n_bins, readflag);
    sift_read_keypoints(k2, v[2], n_hist, n_ori, n_bins, readflag);
    print_elapsed_time(&ts, "read input keypoints");

    // rectify keypoints coordinates
    if (fund_mat) {
        double s1[9] = {0};
        double s2[9] = {0};
        rectifying_similarities_from_affine_fundamental_matrix(s1, s2, fund_mat);
        //fprintf(stderr, "  : %f %f %f\n", s1[0], s1[1], s1[2]);
        //fprintf(stderr, "s1: %f %f %f\n", s1[3], s1[4], s1[5]);
        //fprintf(stderr, "  : %f %f %f\n", s1[6], s1[7], s1[8]);
        //fprintf(stderr, "  : %f %f %f\n", s2[0], s2[1], s2[2]);
        //fprintf(stderr, "s2: %f %f %f\n", s2[3], s2[4], s2[5]);
        //fprintf(stderr, "  : %f %f %f\n", s2[6], s2[7], s2[8]);
        //apply_homography_to_keypoints(k1, h1);
        //apply_homography_to_keypoints(k2, h1);
    }

    // Matching
    matching(k1, k2, out_k1, out_k2, thresh, meth_flag);
    print_elapsed_time(&ts, "compute matches");

    // Print
    fprintf_pairs(output_file, out_k1, out_k2);
    print_elapsed_time(&ts, "print output");
    fprintf(stderr, "%d matches\n", out_k1->size);

//    char name[FILENAME_MAX];
//    if(verb_flag == 1){
//        save_pairs_extra("OUTmatches.txt", out_k1, out_k2A, out_k2B);
//        sprintf(name, "%s_im0.txt", label);
//        sift_save_keypoints(out_k1, name, 1);
//        sprintf(name, "%s_im1.txt", label);
//        sift_save_keypoints(out_k2A, name, 1);
//    }
//
//    // Free memory
    sift_free_keypoints(k1);
    sift_free_keypoints(k2);
    free(out_k1);
    free(out_k2);

    return EXIT_SUCCESS;
}
