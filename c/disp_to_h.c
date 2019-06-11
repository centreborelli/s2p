#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "vvector.h"
#include "iio.h"
#include "rpc.h"
#include "read_matrix.c"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

void utm_alt_zone(double *out, double lat, double lon, int zone);

void applyHom(double outv[3], double M[3][3], double v[3]) {
    MAT_DOT_VEC_3X3(outv, M, v);
    outv[0] /= outv[2];
    outv[1] /= outv[2];
    outv[2] /= outv[2];
}


void disp_to_xyz(float *xyz, float *err,  // outputs
                 float *dispx, float *dispy, float *msk, int nx, int ny,  // inputs
                 double Ha[3][3], double Hb[3][3],
                 struct rpc *rpca, struct rpc *rpcb, int zone)
{
    // invert homographies
    double det;
    double invHa[3][3];
    double invHb[3][3];
    INVERT_3X3(invHa, det, Ha);
    INVERT_3X3(invHb, det, Hb);

    // handle UTM zone
    bool hem;

    int npoints = 0;
    for (int y = 0; y < ny; y++)
    for (int x = 0; x < nx; x++) {
        int pos = x + nx*y;
        if (msk[pos] <= 0) {
            for (int k = 0; k < 3; k++)
                xyz[3 * pos + k] = NAN;
            err[pos] = NAN;
            continue;
        }
        double q0[3], q1[3];
        double lonlat[2];
        double utm[2];
        double e, z;
        double dx = dispx[pos];
        double dy = dispy[pos];
        double p0[3] = {x, y, 1};
        double p1[3] = {x+dx, y+dy, 1};
        applyHom(q0, invHa, p0);
        applyHom(q1, invHb, p1);

        // compute the coordinates
        z = rpc_height(rpca, rpcb, q0[0], q0[1], q1[0], q1[1], &e);
        eval_rpc(lonlat, rpca, q0[0], q0[1], z);

        // convert (lon, lat, alt) to utm
        utm_alt_zone(utm, lonlat[1], lonlat[0], zone);

        // store the output values
        xyz[3 * pos + 0] = utm[0];
        xyz[3 * pos + 1] = utm[1];
        xyz[3 * pos + 2] = z;
        err[pos] = e;
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
    for (int y = p; y < ny - p; y++)
    for (int x = p; x < nx - p; x++) {
        int pos = x + nx * y;
        float *v = xyz + pos * 3;
        int c = 0;
        for (int i = -p; i <= p; i++)
        for (int j = -p; j <= p; j++) {
            float *u = xyz + (x + j + nx * (y + i)) * 3;
            float d = squared_distance_between_3d_points(u, v);
            if (d < r*r) {
                c++;
            }
        }
        count[pos] = c;
    }
}


int main_disp_to_h(int c, char *v[])
{
    if (c != 9) {
        fprintf(stderr, "usage:\n\t"
                "%s rpca rpcb Ha Hb dispAB mskAB out_heights RPCerr utm_zone"
              // 0   1   2    3   4   5      6       7         8      9
                "\n", *v);
        return EXIT_FAILURE;
    }

    // read input data
    struct rpc rpca[1], rpcb[1];
    read_rpc_file_xml(rpca, v[1]);
    read_rpc_file_xml(rpcb, v[2]);
    double Ha[3][3], Hb[3][3];
    read_matrix(Ha, v[3]);
    read_matrix(Hb, v[4]);
    char *fout_heights  = v[7];
    char *fout_err = v[8];
    int zone = atoi(v[9]);

    int nx, ny, nch;
    float *dispy;
    float *dispx = iio_read_image_float_split(v[5], &nx, &ny, &nch);
    if (nch > 1) dispy = dispx + nx*ny;
    else dispy = calloc(nx*ny, sizeof(*dispy));
    float *msk  = iio_read_image_float_split(v[6], &nx, &ny, &nch);

    float *xyz_map = calloc(nx*ny*3, sizeof(*xyz_map));
    float *errMap = calloc(nx*ny, sizeof(*errMap));
    disp_to_xyz(xyz_map, errMap, dispx, dispy, msk, nx, ny, Ha, Hb, rpca, rpcb, zone);

    // save the height map and error map
    iio_save_image_float_vec(fout_heights, xyz_map, nx, ny, 3);
    iio_save_image_float_vec(fout_err, errMap, nx, ny, 1);
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
