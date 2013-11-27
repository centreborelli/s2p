#include <stdint.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


#include "iio.h"

#include "fail.c"

#include "rpc.h"
#include "read_matrix.c"

static void mercator(double m[2], double x[2])
{
    double R = 6378100;
    double deg = M_PI/180;
    m[0] = R * x[0] * deg;
    m[1] = R * log( ( 1 + sin(x[1]*deg) ) / cos(x[1]*deg) );
}

static void getxyz(double xyz[3], struct rpc *r, double i, double j, double h)
{
    double tmp[2];
    eval_rpc(tmp, r, i, j, h);
    mercator(xyz, tmp);
    xyz[2] = h;
}

static void apply_homography(double y[2], double H[3][3], double x[2])
{
    double z = H[2][0]*x[0] + H[2][1]*x[1] + H[2][2];
    y[0] = (H[0][0]*x[0] + H[0][1]*x[1] + H[0][2]) / z;
    y[1] = (H[1][0]*x[0] + H[1][1]*x[1] + H[1][2]) / z;
}

static double invert_homography(double invH[3][3], double H[3][3])
{
    double *a = H[0], *r = invH[0];
    double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
               - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
    r[0] = (a[4]*a[8]-a[5]*a[7])/det;
    r[1] = (a[2]*a[7]-a[1]*a[8])/det;
    r[2] = (a[1]*a[5]-a[2]*a[4])/det;
    r[3] = (a[5]*a[6]-a[3]*a[8])/det;
    r[4] = (a[0]*a[8]-a[2]*a[6])/det;
    r[5] = (a[2]*a[3]-a[0]*a[5])/det;
    r[6] = (a[3]*a[7]-a[4]*a[6])/det;
    r[7] = (a[1]*a[6]-a[0]*a[7])/det;
    r[8] = (a[0]*a[4]-a[1]*a[3])/det;
    return det;
}

// normalize in place a 3d vector
static void normalize_vector_3d(double vec[3])
{
    const int dim = 3;
    double norm = 0;
    for(int i = 0; i < dim; i++)
        norm += vec[i]*vec[i];
    norm = sqrt(norm);
    for(int i = 0; i < dim; i++)
        vec[i] /= norm;
}


void write_ply_header(FILE* f, int npoints) {
    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
    fprintf(f, "comment created by S2P\n");
    fprintf(f, "element vertex %d\n", npoints);
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "property float nx\n");
    fprintf(f, "property float ny\n");
    fprintf(f, "property float nz\n");
    fprintf(f, "property uchar red\n");
    fprintf(f, "property uchar green\n");
    fprintf(f, "property uchar blue\n");
    fprintf(f, "end_header\n");
}

#include "smapa.h"
SMART_PARAMETER_SILENT(IJMESH, 0)
SMART_PARAMETER_SILENT(IJMESHFAC, 2)

int main(int c, char *v[])
{
    if (c != 6) {
        fprintf(stderr, "usage:\n\t"
                "%s colors heights rpc Hfile.txt out.ply\n", *v);
        //       0    1      2       3   4        5
        return 1;
    }
    char *fname_colors = v[1];
    char *fname_heights = v[2];
    char *fname_rpc = v[3];
    double H[3][3], invH[3][3];
    read_matrix(H, v[4]);
    FILE *out = fopen(v[5], "w");
    invert_homography(invH, H);

    int w, h, pd, ww, hh;
    uint8_t *colors = iio_read_image_uint8_vec(fname_colors, &w, &h, &pd);
    float *heights = iio_read_image_float(fname_heights, &ww, &hh);
    if (w != ww || h != hh) fail("color and height image size mismatch");
    if (pd != 1 && pd != 3) fail("expecting a gray or color image");

    struct rpc r[1]; read_rpc_file_xml(r, fname_rpc);

    uint8_t (*color)[w][pd] = (void*)colors;
    float (*height)[w] = (void*)heights;

    // count number of valid pixels
    int npoints = 0;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
            if (!isnan(height[j][i]))
                npoints++;

    // print header for ply file
    write_ply_header(out, npoints);

    // print points coordinates and values
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
        if (!isnan(height[j][i])) {
            // if it is a greyscale image, copy the grey level on each one
            // of the rgb channels
            uint8_t rgb[3];
            for (int k = 0; k < pd; k++) rgb[k] = color[j][i][k];
            for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];

            // convert the pixel local coordinates (ie in the crop) to the
            // full image coordinates, then to mercator coordinates through
            // the rpc_eval direct estimation function normals (ie unit 3D
            // vector pointing in the direction of the camera) are saved
            // together with the coordinates of the 3D points
            double xy[2] = {i, j}, pq[2];
            apply_homography(pq, invH, xy);
            double xyz[3] = {pq[1], pq[0], IJMESHFAC() * height[j][i]};
            double nrm[3] = {0, 0, 1};
            if (!IJMESH()) {
                double tmp[3];
                getxyz(xyz, r, pq[0], pq[1], height[j][i]);
                getxyz(tmp, r, pq[0], pq[1], height[j][i] + 10);
                nrm[0] = tmp[0] - xyz[0];
                nrm[1] = tmp[1] - xyz[1];
                nrm[2] = tmp[2] - xyz[2];
                normalize_vector_3d(nrm);
            }

            // print the voxel in the ply output file
            fprintf(out, "%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d %d %d\n",
                    xyz[0], xyz[1], xyz[2], nrm[0], nrm[1], nrm[2], rgb[0],
                    rgb[1], rgb[2]);
        }

    return 0;
}
