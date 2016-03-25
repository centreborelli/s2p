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


static void print_help(char *s)
{
    fprintf(stderr, "\t usage: %s input.tif [-h \"h1 ... h9\"] output.tif width"
            " height [--verbose]\n", s);
}


int main(int c, char* v[])
{
    // read the homography. Default value is identity
    char *hom_string = pick_option(&c, &v, "h", "1 0 0 0 1 0 0 0 1");
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
    int x, y, w, h;
    x = 0;
    y = 0;
    w = 256;
    h = 256;

    // read the needed ROI in the input image
    struct fancy_image *fimg = fancy_image_open(fname_input, (char *) "");
    float *roi = (float*) malloc(w*h*sizeof(float));
    fancy_image_fill_rectangle_float_split(roi, w, h, fimg, 0, x, y);
    Image imI(roi, (const size_t) w, (const size_t) h, 1);  // TODO: avoid this useless copy
    if (verbose) time.getTime("Read image");

    // call the mapping function
    Image imO;
    Parameters params(0, (const size_t) w, (const size_t) h, true);
    runHomography(imI, hom, imO, params);
    if (verbose) time.getTime("Apply the homography");

    // write the output image
    imO.write(fname_output);
    if (verbose) time.getTime("Write image");

    return EXIT_SUCCESS;
}
