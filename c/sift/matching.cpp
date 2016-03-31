/**
 * @file matching.cpp
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

#include "Time.h"
extern "C" {
    #include "pickopt.h"
    #include "sift_anatomy_20141201/lib_keypoint.h"
    #include "sift_anatomy_20141201/lib_matching.h"
}


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

    // Parsing command line
    char *output_file = pick_option(&c, &v, "o", "/dev/stdout");
    bool verbose = pick_option(&c, &v, "-verbose", NULL);

    // read the fundamental matrix
    char *fund_mat_string = pick_option(&c, &v, "f", "");
    int n_fund;
    double *fund_mat = alloc_parse_doubles(9, fund_mat_string, &n_fund);
    if (n_fund != 9) {
        fprintf(stderr, "can not read 3x3 matrix from \"%s\"", fund_mat_string);
        return EXIT_FAILURE;
    }

    // initialise time counter
    Time time;

    // Memory allocation
    struct sift_keypoints* k1 = sift_malloc_keypoints();
    struct sift_keypoints* k2 = sift_malloc_keypoints();
    struct sift_keypoints* out_k1 = sift_malloc_keypoints();
    struct sift_keypoints* out_k2A = sift_malloc_keypoints();
    struct sift_keypoints* out_k2B = sift_malloc_keypoints();

    // Read input keypoint ASCII files
    int readflag = verb_flag + 1;
    sift_read_keypoints(k1, v[1], n_hist, n_ori, n_bins, readflag);
    sift_read_keypoints(k2, v[2], n_hist, n_ori, n_bins, readflag);
    time.get_time("read input keypoints");

    // Matching
    matching(k1, k2, out_k1, out_k2A, out_k2B, thresh, meth_flag);
    time.get_time("compute matches");

    // Print
    fprintf_pairs(output_file, out_k1, out_k2A);
    time.get_time("print output");
    fprintf(stderr, "%d matches\n", out_k1->size);
    char name[FILENAME_MAX];
    if(verb_flag == 1){
        save_pairs_extra("OUTmatches.txt", out_k1, out_k2A, out_k2B);
        sprintf(name, "%s_im0.txt", label);
        sift_save_keypoints(out_k1, name, 1);
        sprintf(name, "%s_im1.txt", label);
        sift_save_keypoints(out_k2A, name, 1);
    }

    // Free memory
    sift_free_keypoints(k1);
    sift_free_keypoints(k2);
    sift_free_keypoints(out_k1);
    sift_free_keypoints(out_k2A);
    sift_free_keypoints(out_k2B);

    return EXIT_SUCCESS;
}
