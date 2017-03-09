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

//! Gdal includes
#include <gdal_priv.h>
#include <cpl_conv.h>

//! Local includes
#include "LibImages/LibImages.h"
#include "Utilities/Time.h"
#include "Utilities/Utilities.h"
#include "Utilities/Parameters.h"
#include "LibHomography/Homography.h"

#include "linalg.c"
#include "pickopt.c"

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

    // open the input image
    GDALAllRegister();
    GDALDataset *poDataset = (GDALDataset *) GDALOpen(fname_input, GA_ReadOnly);
    if (poDataset == NULL) {
        fprintf(stderr, "ERROR: can't open %s\n", fname_input);
        return EXIT_FAILURE;
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
    int size_x = poDataset->GetRasterXSize();
    int size_y = poDataset->GetRasterYSize();
    int pixel_dimension = poDataset->GetRasterCount();
    if (x + w > size_x)
        w = size_x - x;
    if (y + h > size_y)
        h = size_y - y;
    if (w <= 0 || h <= 0) {
        fprintf(stderr, "ERROR: empty roi\n");
        return EXIT_FAILURE;
    }

    // compensate the homography for the translation due to the crop
    double translation[9] = {1, 0, (double) x, 0, 1, (double) y, 0, 0, 1};
    double hom_compensated[9];
    matrix_33_product(hom_compensated, hom, translation);
    if (verbose) time.get_time("Compute needed ROI");

	// create the output image
	Image out(out_w, out_h, pixel_dimension);

	// crop and warp each input dimension separately
	for (int l = 0; l < pixel_dimension; l++)
	{
		// read the needed ROI in the input image
		// note: GDAL bands are 1-based
		GDALRasterBand *poBand = poDataset->GetRasterBand(1 + l);
		float *roi_data = (float *) CPLMalloc(sizeof(float)*w*h);
		int e = poBand->RasterIO(GF_Read, x, y, w, h, roi_data,
				w, h, GDT_Float32, 0, 0);
		if (e != CPLE_None)
			fprintf(stderr, "errorRasterIO(%d) = %d\n", l, e);

		// copy the ROI data to marc's image struct
		Image roi(roi_data, w, h, 1);
		if (verbose) time.get_time("Read needed ROI");

		// call the mapping function
		Image out_l(out_w, out_h, 1);
		Parameters params(0, out_w, out_h, true);
		runHomography(roi, hom_compensated, out_l, params);
		if (verbose) time.get_time("Apply homography");

		// copy the data into the final output image
		float *from = out_l.getPtr(0);
		float *to = out.getPtr(l);
		int n = out_w * out_h;
		memcpy(to, from, n*sizeof*to);

	}
	out.write(fname_output);

	// NOTE:
	// The GDAL objects used above are not freed automatically.
	// (it would be pointless to free them manually just before
	// exiting the program)

    return EXIT_SUCCESS;
}
