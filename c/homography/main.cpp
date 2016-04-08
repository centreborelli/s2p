/**
 * @file main.cpp
 *
 * @brief Main function to resample an image according to a homography.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com> (original version)
 * @author Carlo de Franchis <carlodef@gmail.com> (modified version)
 **/

//! Global includes
#include <ctime>
#include <stdio.h>


//! Local includes
#include "LibImages/LibImages.h"
#include "Utilities/Time.h"
#include "Utilities/Utilities.h"
#include "Utilities/Parameters.h"
#include "LibHomography/Homography.h"

extern "C" {
    #include "linalg.h"
    #include "pickopt.h"
    #include "fancy_image.h"
}


static void int_bounding_box(int output[4], double input[4][2])
{
    double x[4] = {input[0][0], input[1][0], input[2][0], input[3][0]};
    double y[4] = {input[0][1], input[1][1], input[2][1], input[3][1]};
    output[0] = (int) floor(min_n(x, 4));
    output[1] = (int) floor(min_n(y, 4));
    output[2] = (int) ceil(max_n(x, 4) - output[0]);
    output[3] = (int) ceil(max_n(y, 4) - output[1]);
}


static void compute_needed_roi(int *out, double *hom, int w, int h)
{
    double hom_inv[9];
    matrix_33_inverse(hom_inv, hom);
    double roi_after_hom[4][2] = {{(double) 0, (double) 0},
                                  {(double) w, (double) 0},
                                  {(double) w, (double) h},
                                  {(double) 0, (double) h}};
    double roi_before_hom[4][2];
    for (int i = 0; i < 4; i++)
        apply_homography(roi_before_hom[i], hom_inv, roi_after_hom[i]);

    // as crop uses integer coordinates, be careful to round off
    // (x0, y0) before modifying the homography. We want the crop and the
    // translation representing it to do exactly the same thing.
    int_bounding_box(out, roi_before_hom);
}


static void print_help(char *s)
{
    fprintf(stderr, "\t usage: %s input.tif [-h \"h1 ... h9\"] output.tif width"
            " height [--verbose]\n", s);
}


int main(int c, char* v[])
{
    // read the homography. Default value is identity
    const char *hom_string = pick_option(&c, &v, "h", "1 0 0 0 1 0 0 0 1");
    int n_hom;
    double *hom = alloc_parse_doubles(9, hom_string, &n_hom);
    if (n_hom != 9) {
        fprintf(stderr, "can not read 3x3 matrix from \"%s\"", hom_string);
        return EXIT_FAILURE;
    }

    // verbosity
    bool verbose = pick_option(&c, &v, "-verbose", NULL);

    // parse the remaining arguments
    if (c != 5) {
        print_help(*v);
        return EXIT_FAILURE;
    }

    char *fname_input = v[1];
    char *fname_output = v[2];
    int out_w = atoi(v[3]);
    int out_h = atoi(v[4]);

    // initialization of time
    Time time;

    // compute coordinates of the needed ROI
    int roi_coords[4];
    compute_needed_roi(roi_coords, hom, out_w, out_h);
    int x = roi_coords[0];
    int y = roi_coords[1];
    int w = roi_coords[2];
    int h = roi_coords[3];
    //fprintf(stderr, "roi %d %d %d %d\n", x, y, w, h);

    // compensate the homography for the translation due to the crop
    double translation[9] = {1, 0, (double) x, 0, 1, (double) y, 0, 0, 1};
    double hom_compensated[9];
    matrix_33_product(hom_compensated, hom, translation);
    if (verbose) time.get_time("Compute needed ROI");

    // read the needed ROI in the input image
    struct fancy_image *fimg = fancy_image_open(fname_input, (char *) "");
    float *roi = (float*) malloc(w*h*sizeof(float));
    fancy_image_fill_rectangle_float_split(roi, w, h, fimg, 0, x, y);
    if (verbose) time.get_time("Read needed ROI");
    Image imI(roi, (const size_t) w, (const size_t) h, 1);
    if (verbose) time.get_time("Copy ROI to an Image instance");

    // call the mapping function
    Image imO;
    Parameters params(0, (const size_t) out_w, (const size_t) out_h, true);
    runHomography(imI, hom_compensated, imO, params);
    if (verbose) time.get_time("Apply the homography");

    // write the output image
    imO.write(fname_output);
    if (verbose) time.get_time("Write image");

    return EXIT_SUCCESS;
}
