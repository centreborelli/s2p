#include <stdio.h>
#include <stdint.h>

#include "iio.h"
#include "rpc.h"
#include "fail.c"
#include "parsenumbers.c"
#include "pickopt.c"
#include "read_matrix.c"


void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);


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


unsigned char test_little_endian(void)
{
    int x = 1;
    return (*(char*) & (x) == 1);
}


void write_ply_header(FILE* f, bool ascii, int npoints, int zone, bool hem,
        bool colors, bool normals, bool extra)
{
    if (!ascii)
        if (!test_little_endian())
              fail("BINARY PLY NOT SUPPORTED ON BIG ENDIAN SYSTEMS!\n");

    fprintf(f, "ply\n");
    if (ascii)
        fprintf(f, "format ascii 1.0\n");
    else
        fprintf(f, "format binary_little_endian 1.0\n");

    fprintf(f, "comment created by S2P\n");
    if (zone >= 0)
        fprintf(f, "comment projection: UTM %02d%s\n", zone, (hem ? "N" : "S"));
    fprintf(f, "element vertex %d\n", npoints);
    fprintf(f, "property double x\n");
    fprintf(f, "property double y\n");
    fprintf(f, "property double z\n");
    if (normals) {
        fprintf(f, "property double nx\n");
        fprintf(f, "property double ny\n");
        fprintf(f, "property double nz\n");
    }
    if (colors) {
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
    }
    if (extra) {
        fprintf(f, "property double extr\n");
    }
    fprintf(f, "end_header\n");
}


static void parse_utm_string(int *zone, bool *hem, const char *s)
{
    if (0 == strcmp(s, "no_utm_zone")) {
        *zone = -1;
        return;
    }
    char hem_string[FILENAME_MAX];
    if (2 == sscanf(s, "%02d%s", zone, hem_string)) {
        // hem_string must be equal to "N" or "S"
        if (hem_string[1] == '\0') {
            if (hem_string[0] == 'N' || hem_string[0] == 'S') {
                *hem = (hem_string[0] == 'N');
                return;
            }
        }
    }
    fprintf(stderr, "zone: %d\themisphere: %s\n", *zone, hem_string);
    fprintf(stderr, "incorrect value for --utm-zone."
            " It must look like '27N'\n");
    *zone = -1;
}


void intersect_rays(double out[3], double p[2], double q[2], struct rpc *r1,
        struct rpc *r2)
{
    // compute height
    double err;
    out[2] = rpc_height(r1, r2, p[0], p[1], q[0], q[1], &err);

    // compute lon, lat
    eval_rpc(out, r1, p[0], p[1], out[2]);
}


static void help(char *s)
{
    fprintf(stderr, "\t usage: %s out.ply disp.tif msk.png "
            "rpc_ref.xml rpc_sec.xml [colors.png] [extra.tif]"
            "[-href \"h1 ... h9\"] [-hsec \"h1 ... h9\"] "
            "[--utm-zone ZONE] [--ascii] [--mask-orig msk.png] "
            "[--lon-m l0] [--lon-M lf] [--lat-m l0] [--lat-M lf] "
            "[--col-m x0] [--col-M xf] [--row-m y0] [--row-M yf]\n", s);
}


int main(int c, char *v[])
{
    if (c < 6 || c > 48) {
        help(v[0]);
        return 1;
    }

    // utm zone and hemisphere: true for 'N' and false for 'S'
    int zone;
    bool hem;
    const char *utm_string = pick_option(&c, &v, "-utm-zone", "no_utm_zone");
    parse_utm_string(&zone, &hem, utm_string);

    // ascii flag
    bool ascii = pick_option(&c, &v, "-ascii", NULL);

    // longitude-latitude bounding box
    double lon_m = atof(pick_option(&c, &v, "-lon-m", "-inf"));
    double lon_M = atof(pick_option(&c, &v, "-lon-M", "inf"));
    double lat_m = atof(pick_option(&c, &v, "-lat-m", "-inf"));
    double lat_M = atof(pick_option(&c, &v, "-lat-M", "inf"));

    // x-y bounding box
    double col_m = atof(pick_option(&c, &v, "-col-m", "-inf"));
    double col_M = atof(pick_option(&c, &v, "-col-M", "inf"));
    double row_m = atof(pick_option(&c, &v, "-row-m", "-inf"));
    double row_M = atof(pick_option(&c, &v, "-row-M", "inf"));

    // mask on the unrectified image grid
    const char *msk_orig_fname = pick_option(&c, &v, "-mask-orig", "");
    int msk_orig_w, msk_orig_h;
    float *msk_orig = iio_read_image_float(msk_orig_fname, &msk_orig_w,
                                           &msk_orig_h);

    // rectifying homography
    double href_inv[9], hsec_inv[9];
    int n_hom;
    const char *hom_string_ref = pick_option(&c, &v, "href", "");
    if (*hom_string_ref) {
        double *hom = alloc_parse_doubles(9, hom_string_ref, &n_hom);
        if (n_hom != 9)
            fail("can not read 3x3 matrix from \"%s\"", hom_string_ref);
        invert_homography(href_inv, hom);
    }
    const char *hom_string_sec = pick_option(&c, &v, "hsec", "");
    if (*hom_string_sec) {
        double *hom = alloc_parse_doubles(9, hom_string_sec, &n_hom);
        if (n_hom != 9)
            fail("can not read 3x3 matrix from \"%s\"", hom_string_sec);
        invert_homography(hsec_inv, hom);
    }

    // open disp and mask input images
    int w, h, nch, ww, hh, pd;
    float *dispy;
    float *dispx = iio_read_image_float_split(v[2], &w, &h, &nch);
    if (nch > 1) dispy = dispx + w*h;
    else dispy = calloc(w*h, sizeof(*dispy));

    float *mask = iio_read_image_float(v[3], &ww, &hh);
    if (w != ww || h != hh) fail("disp and mask image size mismatch\n");

    // open color images if provided
    uint8_t *clr = NULL;
    if (c > 6) {
        clr = iio_read_image_uint8_vec(v[6], &ww, &hh, &pd);
        if (w != ww || h != hh) fail("disp and color image size mismatch\n");
    }

    // open extra channel image if provided
    uint8_t *extr = NULL;
    if (c > 7) {
        extr = iio_read_image_uint8_vec(v[7], &ww, &hh, &pd);
        if (w != ww || h != hh) fail("disp and extra image size mismatch\n");
    }

    // read input rpc models
    struct rpc rpc_ref[1], rpc_sec[1];
    read_rpc_file_xml(rpc_ref, v[4]);
    read_rpc_file_xml(rpc_sec, v[5]);

    // outputs
    double p[2], q[2], X[3];

    // count number of valid pixels, and determine utm zone
    int npoints = 0;
    for (int row=0; row<h; row++)
    for (int col=0; col<w; col++) {
        int pix = row*w + col;
        if (!mask[pix]) continue;

        // compute coordinates of pix in the full reference image
        double a[2] = {col, row};
        apply_homography(p, href_inv, a);

        // check that it lies in the image domain bounding box
        if (round(p[0]) < col_m || round(p[0]) > col_M ||
            round(p[1]) < row_m || round(p[1]) > row_M)
            continue;

        // check that it passes the image domain mask
        int x = (int) round(p[0]) - col_m;
        int y = (int) round(p[1]) - row_m;
        if ((x < msk_orig_w) && (y < msk_orig_h))
            if (!msk_orig[y * msk_orig_w + x])
                continue;

        // compute (lon, lat, alt) of the 3D point
        float dx = dispx[pix];
        float dy = dispy[pix];
        double b[2] = {col + dx, row + dy};
        apply_homography(q, hsec_inv, b);
        intersect_rays(X, p, q, rpc_ref, rpc_sec);

        // check with lon/lat bounding box
        if (X[0] < lon_m || X[0] > lon_M || X[1] < lat_m || X[1] > lat_M)
            continue;

        // if it passed all these tests then it's a valid point
        npoints++;

        // if not defined, utm zone is that of the first point
        if (zone < 0)
            utm_zone(&zone, &hem, X[1], X[0]);
    }

    // print header for ply file
    FILE *ply_file = fopen(v[1], "w");
    write_ply_header(ply_file, ascii, npoints, zone, hem, (bool) clr, false, (bool) extr);

    // loop over all the pixels of the input disp map
    // a 3D point is produced for each non-masked disparity
    for (int row=0; row<h; row++)
    for (int col=0; col<w; col++) {
        int pix = row*w + col;
        if (!mask[pix]) continue;

        // compute coordinates of pix in the full reference image
        double a[2] = {col, row};
        apply_homography(p, href_inv, a);

        // check that it lies in the image domain bounding box
        if (round(p[0]) < col_m || round(p[0]) > col_M ||
            round(p[1]) < row_m || round(p[1]) > row_M)
            continue;

        // check that it passes the image domain mask
        int x = (int) round(p[0]) - col_m;
        int y = (int) round(p[1]) - row_m;
        if ((x < msk_orig_w) && (y < msk_orig_h))
            if (!msk_orig[y * msk_orig_w + x])
                continue;

        // compute (lon, lat, alt) of the 3D point
        float dx = dispx[pix];
        float dy = dispy[pix];
        double b[2] = {col + dx, row + dy};
        apply_homography(q, hsec_inv, b);
        intersect_rays(X, p, q, rpc_ref, rpc_sec);

        // check with lon/lat bounding box
        if (X[0] < lon_m || X[0] > lon_M || X[1] < lat_m || X[1] > lat_M)
            continue;

        // convert (lon, lat, alt) to utm
        utm_alt_zone(X, X[1], X[0], zone);

        // colorization: if greyscale, copy the grey level on each channel
        uint8_t rgb[3];
        if (clr) {
            for (int k = 0; k < pd; k++) rgb[k] = clr[k + pd*pix];
            for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
        }
        double extra[1];
        if (extr) {
            extra[0] = extr[pix];
        }

        // write to ply
        if (ascii) {
            fprintf(ply_file, "%0.17g %0.17g %0.17g ", X[0], X[1], X[2]);
            if (clr)
                fprintf(ply_file, "%d %d %d", rgb[0], rgb[1], rgb[2]);
            fprintf(ply_file, "\n");
        } else {
            double XX[3] = {X[0], X[1], X[2]};
            fwrite(XX, sizeof(double), 3, ply_file);
            if (clr) {
                fwrite(rgb, sizeof(unsigned char), 3, ply_file);
            }
            if (extr) {
                fwrite(extra, sizeof(double), 1, ply_file);
            }
        }
    }

    fclose(ply_file);
    return 0;
}
