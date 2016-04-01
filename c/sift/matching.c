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
    // set default parameters
    int n_hist = 4;
    int n_ori = 8;
    int n_bins = 36;
    int meth_flag = 1;
    float sift_thresh = 0.6;
    int verb_flag = 0;
    char label[256];
    strcpy(label, "extra");

    // parse arguments
    const char *output_file = pick_option(&c, &v, "o", "/dev/stdout");
    bool verbose = pick_option(&c, &v, "-verbose", NULL);
    float epi_thresh = atof(pick_option(&c, &v, "-epipolar-threshold", "10"));

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

    // initialise timer
    struct timespec ts; portable_gettime(&ts);

    // memory allocation
    struct sift_keypoints* k1 = sift_malloc_keypoints();
    struct sift_keypoints* k2 = sift_malloc_keypoints();
    struct sift_keypoints* out_k1 = sift_malloc_keypoints();
    struct sift_keypoints* out_k2 = sift_malloc_keypoints();

    // Read input keypoint ASCII files
    int readflag = verb_flag + 1;
    sift_read_keypoints(k1, v[1], n_hist, n_ori, n_bins, readflag);
    sift_read_keypoints(k2, v[2], n_hist, n_ori, n_bins, readflag);
    print_elapsed_time(&ts, "read input keypoints:", 35);

    // matching
    matching(k1, k2, out_k1, out_k2, sift_thresh, meth_flag, fund_mat, epi_thresh);
    print_elapsed_time(&ts, "compute matches:", 35);

    // print
    fprintf_pairs(output_file, out_k1, out_k2);
    print_elapsed_time(&ts, "print output:", 35);
    fprintf(stderr, "%d matches\n", out_k1->size);

    // cleanup
    sift_free_keypoints(k1);
    sift_free_keypoints(k2);
    free(out_k1);
    free(out_k2);

    return EXIT_SUCCESS;
}
