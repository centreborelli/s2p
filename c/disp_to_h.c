#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "iio.h"
#include "rpc.h"
#include "fail.c"
#include "parsenumbers.c"
#include "pickopt.c"


void utm_alt_zone(double *out, double lat, double lon, int zone);

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


void disp_to_xyz(float *xyz, float *err,  // outputs
                 float *dispx, float *dispy, float *msk, int nx, int ny,  // inputs
                 double ha[9], double hb[9],
                 struct rpc *rpca, struct rpc *rpcb, int zone,
                 float orig_img_bounding_box[4])
{
    // invert homographies
    double ha_inv[9];
    double hb_inv[9];
    invert_homography(ha_inv, ha);
    invert_homography(hb_inv, hb);

    // read image bounding box
    float col_min = orig_img_bounding_box[0];
    float col_max = orig_img_bounding_box[1];
    float row_min = orig_img_bounding_box[2];
    float row_max = orig_img_bounding_box[3];

    // initialize output images to nan
    for (int row = 0; row < ny; row++)
    for (int col = 0; col < nx; col++) {
        int pix = col + nx*row;
        err[pix] = NAN;
        for (int k = 0; k < 3; k++)
            xyz[3 * pix + k] = NAN;
    }

    // intermediate buffers
    double p[2], q[2], lonlat[2], utm[2];
    double e, z;

    // loop over all the pixels of the input disp map
    // a 3D point is produced for each non-masked disparity
    for (int row = 0; row < ny; row++)
    for (int col = 0; col < nx; col++) {
        int pix = col + nx*row;
        if (!msk[pix])
            continue;

        // compute coordinates of pix in the full reference image
        double a[2] = {col, row};
        apply_homography(p, ha_inv, a);

        // check that it lies in the image domain bounding box
        if (round(p[0]) < col_min || round(p[0]) > col_max ||
            round(p[1]) < row_min || round(p[1]) > row_max)
            continue;

        // compute (lon, lat, alt) of the 3D point
        double dx = dispx[pix];
        double dy = dispy[pix];
        double b[2] = {col + dx, row + dy};
        apply_homography(q, hb_inv, b);
        z = rpc_height(rpca, rpcb, p[0], p[1], q[0], q[1], &e);
        eval_rpc(lonlat, rpca, p[0], p[1], z);

        // convert (lon, lat, alt) to utm
        utm_alt_zone(utm, lonlat[1], lonlat[0], zone);

        // store the output values
        xyz[3 * pix + 0] = utm[0];
        xyz[3 * pix + 1] = utm[1];
        xyz[3 * pix + 2] = z;
        err[pix] = e;
    }
}


float squared_distance_between_3d_points(float a[3], float b[3])
{
    float x = (a[0] - b[0]);
    float y = (a[1] - b[1]);
    float z = (a[2] - b[2]);
    return x*x + y*y + z*z;
}


void count_3d_neighbors(int *count, float *xyz, int nx, int ny, float r, int p)
{
    // count the 3d neighbors of each point
    for (int y = 0; y < ny; y++)
    for (int x = 0; x < nx; x++) {
        int pos = x + nx * y;
        float *v = xyz + pos * 3;
        int c = 0;
        int i0 = y > p ? -p : -y;
        int i1 = y < ny - p ? p : ny - y - 1;
        int j0 = x > p ? -p : -x;
        int j1 = x < nx - p ? p : nx - x - 1;
        for (int i = i0; i <= i1; i++)
        for (int j = j0; j <= j1; j++) {
            float *u = xyz + (x + j + nx * (y + i)) * 3;
            float d = squared_distance_between_3d_points(u, v);
            if (d < r*r) {
                c++;
            }
        }
        count[pos] = c;
    }
}


void remove_isolated_3d_points(
    float* xyz,  // input (and output) image, shape = (h, w, 3)
    int nx,      // width w
    int ny,      // height h
    float r,     // filtering radius, in meters
    int p,       // filtering window (square of width is 2p+1 pixels)
    int n,       // minimal number of neighbors to be an inlier
    int q)       // neighborhood for the saving step (square of width 2q+1)
{
    int *count = (int*) malloc(nx * ny * sizeof(int));
    bool *rejected = (bool*) malloc(nx * ny * sizeof(bool));

    // count the 3d neighbors of each point
    count_3d_neighbors(count, xyz, nx, ny, r, p);

    // brutally reject any point with less than n neighbors
    for (int i = 0; i < ny * nx; i++)
        rejected[i] = count[i] < n;

    // show mercy; save points with at least one close and non-rejected neighbor
    bool need_more_iterations = true;
    while (need_more_iterations) {
        need_more_iterations = false;
        // scan the grid and stop at rejected points
        for (int y = 0; y < ny; y++)
        for (int x = 0; x < nx; x++)
        if (rejected[x + y * nx])
        // explore the neighborhood (square of width 2q+1)
        for (int yy = y - q; yy < y + q + 1; yy++) {
            if (yy < 0) continue; else if (yy > ny-1) break;
            for (int xx = x - q; xx < x + q + 1; xx++) {
                if (xx < 0) continue; else if (xx > nx-1) break;
                // is the current rejected point's neighbor non-rejected?
                if (!rejected[xx + yy * nx])
                // is this connected neighbor close (in 3d)?
                if (squared_distance_between_3d_points(
                        xyz + (x + y * nx)*3, xyz + (xx + yy * nx)*3) < r*r) {
                    rejected[x + y * nx] = false; // save the point
                    yy = xx = ny + nx + 2*q + 2;  // break loops over yy and xx
                    need_more_iterations = true;  // this point may save others!
                }
            }
        }
    }

    // set to NAN the rejected pixels
    for (int i = 0; i < ny * nx; i++)
        if (rejected[i])
            for (int c = 0; c < 3; c++)
                xyz[c + i * 3] = NAN;

    free(rejected);
    free(count);
}


static void help(char *s)
{
    fprintf(stderr, "usage:\n\t"
            "%s rpc_ref.xml rpc_sec.xml disp.tif heights.tif err.tif "
            "[--mask mask.png] "
            "[-href \"h1 ... h9\"] [-hsec \"h1 ... h9\"] "
            "[--utm-zone ZONE] [--mask-orig msk.png] "
            "[--col-m x0] [--col-M xf] [--row-m y0] [--row-M yf]\n", s);
}


int main_disp_to_h(int c, char *v[])
{
    if (c < 6 || c > 48) {
        help(v[0]);
        return EXIT_FAILURE;
    }

    // read input mask
    char *mask_path = pick_option(&c, &v, "-mask", "");

    // rectifying homographies
    double ha[9], hb[9];
    int n_hom;
    const char *hom_string_ref = pick_option(&c, &v, "href", "");
    if (*hom_string_ref) {
        double *ha = alloc_parse_doubles(9, hom_string_ref, &n_hom);
        if (n_hom != 9)
            fail("can not read 3x3 matrix from \"%s\"", hom_string_ref);
    }
    const char *hom_string_sec = pick_option(&c, &v, "hsec", "");
    if (*hom_string_sec) {
        double *hb = alloc_parse_doubles(9, hom_string_sec, &n_hom);
        if (n_hom != 9)
            fail("can not read 3x3 matrix from \"%s\"", hom_string_sec);
    }

    // utm zone
    int zone = atoi(pick_option(&c, &v, "-utm-zone", ""));

    // x-y bounding box
    double col_m = atof(pick_option(&c, &v, "-col-m", "-inf"));
    double col_M = atof(pick_option(&c, &v, "-col-M", "inf"));
    double row_m = atof(pick_option(&c, &v, "-row-m", "-inf"));
    double row_M = atof(pick_option(&c, &v, "-row-M", "inf"));

    // remaining positional arguments: rpcs, disparity map, output files
    struct rpc rpca[1], rpcb[1];
    read_rpc_file_xml(rpca, v[1]);
    read_rpc_file_xml(rpcb, v[2]);
    char *disp_path = v[3];
    char *fout_heights  = v[4];
    char *fout_err = v[5];

    int nx, ny, nch;
    float *dispy;
    float *dispx = iio_read_image_float_split(disp_path, &nx, &ny, &nch);
    if (nch > 1) dispy = dispx + nx*ny;
    else dispy = calloc(nx*ny, sizeof(*dispy));
    float *msk  = iio_read_image_float_split(mask_path, &nx, &ny, &nch);

    // triangulation
    float *xyz_map = calloc(nx*ny*3, sizeof(*xyz_map));
    float *err_map = calloc(nx*ny, sizeof(*err_map));
    float img_bbx[4] = {col_m, col_M, row_m, row_M};
    disp_to_xyz(xyz_map, err_map, dispx, dispy, msk, nx, ny, ha, hb,
                rpca, rpcb, zone, img_bbx);

    // save the height map and error map
    iio_write_image_float_vec(fout_heights, xyz_map, nx, ny, 3);
    iio_write_image_float_vec(fout_err, err_map, nx, ny, 1);
    return 0;
}


int main_count_3d_neighbors(int c, char *v[])
{
    if (c != 5) {
        fprintf(stderr, "usage:\n\t"
                "%s xyz.tif r p out.tif"
              // 0   1      2 3 4
                "\n", *v);
        return EXIT_FAILURE;
    }

    // read input data
    int nx, ny, nch;
    float *xyz = iio_read_image_float_vec(v[1], &nx, &ny, &nch);
    if (nch != 3) fprintf(stderr, "xyz image must have 3 channels\n");
    float r = atof(v[2]);
    int p = atoi(v[3]);
    char *output_filename = v[4];

    // allocate output data
    int *out = calloc(nx*ny, sizeof(*out));

    // do the job
    count_3d_neighbors(out, xyz, nx, ny, r, p);

    // save output
    iio_write_image_int(output_filename, out, nx, ny);
    return 0;
}

int main(int c, char *v[])
{
    return main_disp_to_h(c, v);
//    return main_count_3d_neighbors(c, v);
}
