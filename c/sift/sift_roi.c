#include <stdlib.h>

#include "timing.h"
#include "pickopt.h"
#include "sift_anatomy_20141201/lib_sift_anatomy.h"
#include "gdal.h"
#include "cpl_conv.h"


void print_help(char *v[])
{
    fprintf(stderr, "usage:\n\t%s file.tif x y w h [-o file] [--max-nb-pts n]"
    //                          0 1        2 3 4 5
            " [-b] [--verbose] [--thresh-dog t (0.0133)]"
            " [--scale-space-noct n (8)] [--scale-space-nspo n (3)]\n", *v);
}


int main(int c, char *v[]) {
    // process input arguments
    if (c < 2) {
        print_help(v);
        return 1;
    }

    // optional arguments
    const char *output_file = pick_option(&c, &v, "o", "/dev/stdout");
    bool binary = (bool) pick_option(&c, &v, "b", NULL);
    bool verbose = (bool) pick_option(&c, &v, "-verbose", NULL);
    int max_nb_pts = atoi(pick_option(&c, &v, "-max-nb-pts", "INT_MAX"));
    float thresh_dog = (float) atof(pick_option(&c, &v, "-thresh-dog", "0.0133"));
    int ss_noct = atoi(pick_option(&c, &v, "-scale-space-noct", "8"));
    int ss_nspo = atoi(pick_option(&c, &v, "-scale-space-nspo", "3"));

    // initialise time
    struct timespec ts;
    portable_gettime(&ts);

    // define the rectangular region of interest (roi)
    int x, y, w, h;
    if (c == 6) {
        x = atoi(v[2]);
        y = atoi(v[3]);
        w = atoi(v[4]);
        h = atoi(v[5]);
    } else {
        print_help(v);
        return 1;
    }

    // open the input image
    GDALAllRegister();
    GDALDatasetH hDataset = GDALOpen(v[1], GA_ReadOnly);
    if (hDataset == NULL) {
        fprintf(stderr, "ERROR: can't open %s\n", v[1]);
        return 1;
    }

    // clip roi to stay inside the image boundaries
    if (x < 0) {
        w += x;
        x = 0;
    }
    if (y < 0) {
        h += y;
        y = 0;
    }
    int size_x = GDALGetRasterXSize(hDataset);
    int size_y = GDALGetRasterYSize(hDataset);
    if (x + w > size_x)
        w = size_x - x;
    if (y + h > size_y)
        h = size_y - y;
    if (w <= 0 || h <= 0) {
        fprintf(stderr, "WARNING: empty roi\n");
        return 0;
    }

    // read roi
    GDALRasterBandH hBand = GDALGetRasterBand(hDataset, 1);
    float *roi = (float *) CPLMalloc(sizeof(float)*w*h);
    GDALRasterIO(hBand, GF_Read, x, y, w, h, roi, w, h, GDT_Float32, 0, 0);
    GDALClose(hBand);
    GDALClose(hDataset);
    if (verbose) print_elapsed_time(&ts, "read ROI", 35);

    // prepare sift parameters
    struct sift_parameters* p = sift_assign_default_parameters();
    p->C_DoG = thresh_dog;
    p->n_oct = ss_noct;
    p->n_spo = ss_nspo;

    // compute sift keypoints
    struct sift_scalespace **ss = (struct sift_scalespace**) malloc(4 * sizeof(struct sift_scalespace*));
    struct sift_keypoints **kk = (struct sift_keypoints**) malloc(6 * sizeof(struct sift_keypoints*));
    for (int i = 0; i < 6; i++)
        kk[i] = sift_malloc_keypoints();
    struct sift_keypoints* kpts = sift_anatomy(roi, w, h, p, ss, kk);
    if (verbose) print_elapsed_time(&ts, "run SIFT", 35);

    // add (x, y) offset to keypoints coordinates
    for (int i = 0; i < kpts->size; i++) {
        kpts->list[i]->x += y;  // in Ives' conventions x is the row index
        kpts->list[i]->y += x;
    }
    if (verbose) print_elapsed_time(&ts, "add offset", 35);

    // write to standard output
    FILE *f = fopen(output_file, "w");
    fprintf_keypoints(f, kpts, max_nb_pts, binary, 1);
    fclose(f);
    if (verbose) print_elapsed_time(&ts, "write output", 35);

    // cleanup
    CPLFree(roi);
    sift_free_keypoints(kpts);
    for (int i = 0; i < 6; i++)
        sift_free_keypoints(kk[i]);
    free(kk);
    for (int i = 0; i < 4; i++)
        sift_free_scalespace(ss[i]);
    free(ss);
    free(p);
    return 0;
}
