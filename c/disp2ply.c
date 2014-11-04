#include <stdio.h>
#include <stdint.h>

#include "iio.h"
#include "rpc.h"
#include "fail.c"
#include "pickopt.c"
#include "read_matrix.c"


void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);

static void apply_homography(double y[2], double h[3][3], double x[2])
{
    double z = h[2][0]*x[0] + h[2][1]*x[1] + h[2][2];
    y[0] = (h[0][0]*x[0] + h[0][1]*x[1] + h[0][2]) / z;
    y[1] = (h[1][0]*x[0] + h[1][1]*x[1] + h[1][2]) / z;
}


static double invert_homography(double invh[3][3], double h[3][3])
{
    double *a = h[0], *r = invh[0];
    double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
               - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
    r[0] = (a[4]*a[8] - a[5]*a[7]) / det;
    r[1] = (a[2]*a[7] - a[1]*a[8]) / det;
    r[2] = (a[1]*a[5] - a[2]*a[4]) / det;
    r[3] = (a[5]*a[6] - a[3]*a[8]) / det;
    r[4] = (a[0]*a[8] - a[2]*a[6]) / det;
    r[5] = (a[2]*a[3] - a[0]*a[5]) / det;
    r[6] = (a[3]*a[7] - a[4]*a[6]) / det;
    r[7] = (a[1]*a[6] - a[0]*a[7]) / det;
    r[8] = (a[0]*a[4] - a[1]*a[3]) / det;
    return det;
}


unsigned char test_little_endian(void)
{
    int x = 1;
    return (*(char*) & (x) == 1);
}


void write_ply_header(FILE* f, bool ascii, int npoints, int zone, bool hem,
        bool colors, bool normals)
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
        fprintf(f, "comment projection: UTM %i%s\n", zone, (hem ? "N" : "S"));
    fprintf(f, "element vertex %d\n", npoints);
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    if (normals) {
        fprintf(f, "property float nx\n");
        fprintf(f, "property float ny\n");
        fprintf(f, "property float nz\n");
    }
    if (colors) {
        fprintf(f, "property uchar red\n");
        fprintf(f, "property uchar green\n");
        fprintf(f, "property uchar blue\n");
    }
    fprintf(f, "end_header\n");
}


static void parse_utm_string(int *zone, bool *hem, char *s)
{
    if (s == "no_utm_zone") {
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
    fprintf(stderr, "\t usage: %s out.ply img_disp img_msk"
            " hom_ref.txt hom_sec.txt rpc_ref.xml rpc_sec.xml"
            " [img_ref] [--utm-zone ZONE] [--ascii]\n", s);
}


int main(int c, char *v[])
{
    if (c != 8 && c != 9 && c != 10 && c != 11) {
        help(v[0]);
        return 1;
    }

    // utm zone and hemisphere: true for 'N' and false for 'S'
    int zone;
    bool hem;
    char *utm_string = pick_option(&c, &v, "-utm-zone", "no_utm_zone");
    parse_utm_string(&zone, &hem, utm_string);

    // ascii flag
    bool ascii   = (pick_option(&c, &v, "-ascii", NULL) != NULL);

    // open disp and mask input images
    int w, h, ww, hh, pd;
    float *disp = iio_read_image_float(v[2], &w, &h);
    float *mask = iio_read_image_float(v[3], &ww, &hh);
    if (w != ww || h != hh) fail("disp and mask image size mismatch\n");

    // open color images if provided
    uint8_t *clr = NULL;
    if (c > 8) {
        clr = iio_read_image_uint8_vec(v[8], &ww, &hh, &pd);
        if (w != ww || h != hh) fail("disp and color image size mismatch\n");
    }

    // read and invert input homographies
    double hom_ref[3][3], hom_sec[3][3], hom_ref_inv[3][3], hom_sec_inv[3][3];
    read_matrix(hom_ref, v[4]);
    read_matrix(hom_sec, v[5]);
    invert_homography(hom_ref_inv, hom_ref);
    invert_homography(hom_sec_inv, hom_sec);

    // read input rpc models
    struct rpc rpc_ref[1], rpc_sec[1];
    read_rpc_file_xml(rpc_ref, v[6]);
    read_rpc_file_xml(rpc_sec, v[7]);

    // outputs
    double p[2], q[2], X[3];

    // count number of valid pixels, and determine utm zone
    int npoints = 0;
    for (int row=0; row<h; row++)
    for (int col=0; col<w; col++) {
        int pix = row*w + col;
        if (mask[pix]) {
            npoints++;
            // if not defined, utm zone is that of the first point
            if (zone < 0) {
                float d = disp[pix];

                // compute coordinates of pix in the big image
                double a[2] = {col, row};
                double b[2] = {col + d, row};
                apply_homography(p, hom_ref_inv, a);
                apply_homography(q, hom_sec_inv, b);

                // compute (lon, lat, alt) of the 3D point
                intersect_rays(X, p, q, rpc_ref, rpc_sec);
                utm_zone(&zone, &hem, X[1], X[0]);
            }
        }
    }

    // print header for ply file
    FILE *ply_file = fopen(v[1], "w");
    write_ply_header(ply_file, ascii, npoints, zone, hem, (bool) clr, false);

    // loop over all the pixels of the input disp map
    // a 3D point is produced for each non-masked disparity
    for (int row=0; row<h; row++)
    for (int col=0; col<w; col++) {
        int pix = row*w + col;
        if (mask[pix]) {
            float d = disp[pix];

            // compute coordinates of pix in the big image
            double a[2] = {col, row};
            double b[2] = {col + d, row};
            apply_homography(p, hom_ref_inv, a);
            apply_homography(q, hom_sec_inv, b);

            // compute (lon, lat, alt) of the 3D point
            intersect_rays(X, p, q, rpc_ref, rpc_sec);

            // convert (lon, lat, alt) to utm
            utm_alt_zone(X, X[1], X[0], zone);

            // colorization: if greyscale, copy the grey level on each channel
            uint8_t rgb[3];
            if (clr) {
                for (int k = 0; k < pd; k++) rgb[k] = clr[k + pd*pix];
                for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
            }

            // write to ply
            if (ascii) {
                fprintf(ply_file, "%.3f %.3f %.3f ", X[0], X[1], X[2]);
                if (clr)
                    fprintf(ply_file, "%d %d %d", rgb[0], rgb[1], rgb[2]);
                fprintf(ply_file, "\n");
            } else {
                float XX[3] = {X[0], X[1], X[2]};
                fwrite(XX, sizeof(float), 3, ply_file);
                if (clr) {
                    unsigned char C[3] = {rgb[0], rgb[1], rgb[2]};
                    fwrite(rgb, sizeof(unsigned char), 3, ply_file);
                }
            }
        }
    }

    fclose(ply_file);
    return 0;
}
