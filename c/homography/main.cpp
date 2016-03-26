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


//! Local includes
#include "LibImages/LibImages.h"
#include "Utilities/Time.h"
#include "Utilities/Utilities.h"
#include "Utilities/Parameters.h"
#include "LibHomography/Homography.h"

extern "C" {
    #include "pickopt.h"
    #include "fancy_image.h"
}


static void apply_homography(double y[2], double h[9], double x[2])
{
    //                    h[0] h[1] h[2]
    // The convention is: h[3] h[4] h[5]
    //                    h[6] h[7] h[8]
    double z = h[6]*x[0] + h[7]*x[1] + h[8];
    double tmp = x[0];  // to enable calls like 'apply_homography(x, h, x)'
    y[0] = (h[0]*x[0] + h[1]*x[1] + h[2]) / z;
    y[1] = (h[3]*tmp  + h[4]*x[1] + h[5]) / z;
}


static double invert_homography(double o[9], double i[9])
{
    double det = i[0]*i[4]*i[8] + i[2]*i[3]*i[7] + i[1]*i[5]*i[6]
               - i[2]*i[4]*i[6] - i[1]*i[3]*i[8] - i[0]*i[5]*i[7];
    o[0] = (i[4]*i[8] - i[5]*i[7]) / det;
    o[1] = (i[2]*i[7] - i[1]*i[8]) / det;
    o[2] = (i[1]*i[5] - i[2]*i[4]) / det;
    o[3] = (i[5]*i[6] - i[3]*i[8]) / det;
    o[4] = (i[0]*i[8] - i[2]*i[6]) / det;
    o[5] = (i[2]*i[3] - i[0]*i[5]) / det;
    o[6] = (i[3]*i[7] - i[4]*i[6]) / det;
    o[7] = (i[1]*i[6] - i[0]*i[7]) / det;
    o[8] = (i[0]*i[4] - i[1]*i[3]) / det;
    return det;
}


static void matrix_product_3x3(double ab[9], double a[9], double b[9])
{
    ab[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
    ab[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
    ab[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
    ab[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
    ab[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
    ab[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
    ab[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
    ab[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
    ab[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}


static double min(double *x, int n)
{
    double out = *x;
    for (int i = 1; i < n; i++)
        if (*(++x) < out)
            out = *x;
    return out;
}


static double max(double *x, int n)
{
    double out = *x;
    for (int i = 1; i < n; i++)
        if (*(++x) > out)
            out = *x;
    return out;
}


static void int_bounding_box(int output[4], double input[4][2])
{
    double x[4] = {input[0][0], input[1][0], input[2][0], input[3][0]};
    double y[4] = {input[0][1], input[1][1], input[2][1], input[3][1]};
    output[0] = (int) floor(min(x, 4));
    output[1] = (int) floor(min(y, 4));
    output[2] = (int) ceil(max(x, 4) - output[0]);
    output[3] = (int) ceil(max(y, 4) - output[1]);
}


static double *alloc_parse_doubles(int nmax, const char *ss, int *n)
{
    // add a space to s, so that gabriele's expression works as intended
    int ns = strlen(ss);
    char t[2+ns], *s = t;
    s[0] = ' ';
    for (int i = 0; i <= ns; i++)
        s[i+1] = ss[i];

    double *r = (double *) malloc(nmax * sizeof*r);
    int i = 0, w;
    while (i < nmax && 1 == sscanf(s, "%*[][ \n\t,:;]%lf%n", r + i, &w)) {
        i += 1;
        s += w;
    }
    *n = i;
    return r;
}


static void compute_needed_roi(int *out, double *hom, int w, int h)
{
    double hom_inv[9];
    double det = invert_homography(hom_inv, hom);
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
    char *hom_string = pick_option(&c, &v, (char *) "h", (char *) "1 0 0 0 1 0 0 0 1");
    int n_hom;
    double *hom = alloc_parse_doubles(9, hom_string, &n_hom);
    if (n_hom != 9) {
        fprintf(stderr, "can not read 3x3 matrix from \"%s\"", hom_string);
        return EXIT_FAILURE;
    }

    // verbosity
    bool verbose = pick_option(&c, &v, (char *) "-verbose", NULL);

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
    double translation[9] = {1, 0, x, 0, 1, y, 0, 0, 1};
    double hom_compensated[9];
    matrix_product_3x3(hom_compensated, hom, translation);

    // read the needed ROI in the input image
    struct fancy_image *fimg = fancy_image_open(fname_input, (char *) "");
    float *roi = (float*) malloc(w*h*sizeof(float));
    fancy_image_fill_rectangle_float_split(roi, w, h, fimg, 0, x, y);
    Image imI(roi, (const size_t) w, (const size_t) h, 1);  // TODO: avoid this useless copy
    if (verbose) time.getTime("Read image");

    // call the mapping function
    Image imO;
    Parameters params(0, (const size_t) out_w, (const size_t) out_h, true);
    runHomography(imI, hom_compensated, imO, params);
    if (verbose) time.getTime("Apply the homography");

    // write the output image
    imO.write(fname_output);
    if (verbose) time.getTime("Write image");

    return EXIT_SUCCESS;
}
