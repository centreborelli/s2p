#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fail.c"
#include "xmalloc.c"
#include "parsenumbers.c"
#include "rpc.c"
#include "srtm4.c"


// homographic transform y=H(x)
static void apply_homography(double y[2], double H[9], double x[2])
{
    double z[3];
    z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
    z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
    z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
    y[0] = z[0]/z[2];
    y[1] = z[1]/z[2];
}

static double invert_homography(double out[9], double h[9])
{
    double det = h[0]*h[4]*h[8] + h[2]*h[3]*h[7] + h[1]*h[5]*h[6]
        - h[2]*h[4]*h[6] - h[1]*h[3]*h[8] - h[0]*h[5]*h[7];

    out[0] = (h[4]*h[8] - h[5]*h[7]) / det;
    out[1] = (h[2]*h[7] - h[1]*h[8]) / det;
    out[2] = (h[1]*h[5] - h[2]*h[4]) / det;
    out[3] = (h[5]*h[6] - h[3]*h[8]) / det;
    out[4] = (h[0]*h[8] - h[2]*h[6]) / det;
    out[5] = (h[2]*h[3] - h[0]*h[5]) / det;
    out[6] = (h[3]*h[7] - h[4]*h[6]) / det;
    out[7] = (h[1]*h[6] - h[0]*h[7]) / det;
    out[8] = (h[0]*h[4] - h[1]*h[3]) / det;
    return det;
}

void water_mask_fill(int *x, int w, int h, double H[9], struct rpc *r)
{
    double inv_H[9];
    invert_homography(inv_H, H);
    for (int row = 0; row < h; row++) {
        for (int col = 0; col < w; col++) {
            // apply inverse homography H to (col, row) to get the coordinates
            // in the full image
            double p[2] = {(double) col, (double) row};
            double q[2];
            apply_homography(q, inv_H, p);
            // TODO: this should be the geoid height with respect to the ellipsoid
            double geoid_alt = 0.0;
            double geo[2];
            eval_rpc(geo, r, q[0], q[1], geoid_alt);
            double alt = srtm4(geo[0], geo[1]);
            if (alt + 32768.0 < 1) // -32768 is the flag for water
                x[w*row + col] = 0;
        }
    }
}

void print_help(char *v)
{
    fprintf(stderr, "usage:\n\t%s "
        "width height -h \"h1 ... h9\" [rpc.xml [out.png]]\n", v);
        //   1 2                          3           4
}

#define WATERMASK_MAIN
#ifdef WATERMASK_MAIN

#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
    // read input arguments
    char *Hstring = pick_option(&c, &v, "h", "");
    if (c != 5 && c!= 4 && c != 3) {
        print_help(v[0]);
        return 1;
    }
    int w = atoi(v[1]);
    int h = atoi(v[2]);
    char *filename_rpc = c > 3 ? v[3] : "-";
    char *filename_out = c > 4 ? v[4] : "PNG:-";

    // acquire space for output image
    int *x = xmalloc(w*h*sizeof*x);
    for (int i = 0; i < w*h; i++)
        x[i] = 255;

    // read input rpc file
    struct rpc r[1]; read_rpc_file_xml(r, filename_rpc);

    // read the rectifying homography
    if (!*Hstring) {
        print_help(v[0]);
        return 1;
    }
    int nH;
    double *H = alloc_parse_doubles(9, Hstring, &nH);
    if (nH != 9)
        fail("can not read 3x3 matrix from \"%s\"", Hstring);

    // draw mask over output image
    water_mask_fill(x, w, h, H, r);

    // save output image
    iio_save_image_int(filename_out, x, w, h);

    //cleanup
    free(x);
    return 0;
}

#endif//WATERMASK_MAIN
