#include <ctime>
#include <stdlib.h>

#include <gdal.h>
#include <cpl_conv.h>

#include "Utilities/Time.h"
#include "Utilities/Parameters.h"
#include "LibImages/LibImages.h"
#include "LibSift/LibSift.h"

#include "pickopt.c"


static void print_help(char *v[])
{
    fprintf(stderr, "usage:\n\t%s file.tif [x y w h] [-o file]"
    //                          0 1        2 3 4 5
    //        " [-b] [--verbose] [--thresh-dog t (0.0133)]"
            " [--verbose] [--thresh-dog t (0.0133)]"
            " [--scale-space-noct n (8)] [--scale-space-nspo n (3)]\n", *v);
}


int main(int c, char *v[]) {
    if (c < 2) {
        print_help(v);
        return 1;
    }

    // initialise time
    Time time;

    // optional arguments
    const char *output_file = pick_option(&c, &v, "o", "/dev/stdout");
    //bool binary = (bool) pick_option(&c, &v, "b", NULL);
    bool verbose = (bool) pick_option(&c, &v, "-verbose", NULL);
    //int max_nb_pts = atoi(pick_option(&c, &v, "-max-nb-pts", "INT_MAX"));
    float thresh_dog = (float) atof(pick_option(&c, &v, "-thresh-dog", "0.0133"));
    int ss_noct = atoi(pick_option(&c, &v, "-scale-space-noct", "8"));
    int ss_nspo = atoi(pick_option(&c, &v, "-scale-space-nspo", "3"));

    // open the input image
    GDALAllRegister();
    GDALDatasetH hDataset = GDALOpen(v[1], GA_ReadOnly);
    if (hDataset == NULL) {
        fprintf(stderr, "ERROR: can't open %s\n", v[1]);
        return 1;
    }
    int size_x = GDALGetRasterXSize(hDataset);
    int size_y = GDALGetRasterYSize(hDataset);

    // define the rectangular region of interest (roi)
    int x, y, w, h;
    if (c == 6) {
        x = atoi(v[2]);
        y = atoi(v[3]);
        w = atoi(v[4]);
        h = atoi(v[5]);
    } else {
        x = 0;
        y = 0;
        w = size_x;
        h = size_y;
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
    if (x + w > size_x)
        w = size_x - x;
    if (y + h > size_y)
        h = size_y - y;
    if (w <= 0 || h <= 0) {
        fprintf(stderr, "ERROR: empty roi\n");
        return 1;
    }

    // read roi
    GDALRasterBandH hBand = GDALGetRasterBand(hDataset, 1);
    float *roi = (float *) CPLMalloc(sizeof(float)*w*h);
    CPLErr ee = GDALRasterIO(hBand, GF_Read, x, y, w, h, roi, w, h, GDT_Float32, 0, 0);
    GDALClose(hDataset);
    if (verbose) time.get_time("read roi", 35);
    Image im(roi, (const size_t) w, (const size_t) h, 1);
    if (verbose) time.get_time("copy roi into Marc's image object", 35);
    CPLFree(roi);

    // prepare params object
    Parameters params;
    params.setDefaultValues();
    params.set_thresh_dog(thresh_dog);
    params.set_noct(ss_noct);
    params.set_nspo(ss_nspo);

    // run sift
    Sift sift(params);
    if (verbose) time.get_time("initialization", 35);
    sift.computeKeyPoints(im);
    if (verbose) time.get_time("compute keypoints", 35);

    // add (x, y) offset to keypoints coordinates
    std::list<KeyPoint*>::iterator key = sift.m_keyPoints->begin();
    for (; key != sift.m_keyPoints->end(); key++) {
        (*key)->setX((*key)->getX() + y); // in Ives' conventions x is the row index
        (*key)->setY((*key)->getY() + x);
    }
    if (verbose) time.get_time("add offset", 35);

    // write output
    sift.writeKeyPoints(output_file);
    if (verbose) time.get_time("write output", 35);
    return 0;
}
