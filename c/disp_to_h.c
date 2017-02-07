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


void applyHom(double outv[3], double M[3][3], double v[3]) {
    MAT_DOT_VEC_3X3(outv, M, v);
    outv[0] /= outv[2];
    outv[1] /= outv[2];
    outv[2] /= outv[2];
}


/*// convert geodetic coordinates to mercator using a reference longitude
static void convert_geodetic_to_mercator(double mercator[2], double
        geodetic[2], double reference_longitude) {
    double lon0 = reference_longitude;
    double R = 6378.1*1000;
    double cte = 1/360.*2*M_PI;
    mercator[0] = R * (geodetic[0] - lon0) * cte;
    mercator[1] = R * log((1 + sin(geodetic[1] * cte)) / cos(geodetic[1] * cte));
}


// normalize in place a 3d vector
static void normalize_vector_3d(double vec[3]) {
    const int dim = 3;
    double norm = 0;
    for (int i = 0; i < dim ; i++)
        norm += vec[i] * vec[i];
    norm = sqrt(norm);
    for (int i = 0; i < dim ; i++)
        vec[i] /= norm;
}


// Structure that contains all the possible outputs we could wish for points of
// the point cloud
typedef struct world_point world_point;
struct world_point {
    // geodetic coordinates of the point
    float lon;
    float lat;
    float h;
    // surface normal for the geodetic coordinates
    float normal_lon;
    float normal_lat;
    float normal_h;
    // coordinates on a mercator projection
    float x_mercator;
    float y_mercator;
    float exagerated_h;
    // surface normal for the mercator projection
    float normal_x_mercator;
    float normal_y_mercator;
    float normal_h_mercator;
    // projection error
    float rpc_error;
    // image coordinates
    float x;
    float y;
    // color of the pixel
    float r;
    float g;
    float b;
};*/


int main_disp_to_h(int c, char *v[])
{
    if (c != 9) {
        fprintf(stderr, "usage:\n\t"
                "%s rpca rpcb Ha Hb dispAB mskAB out_heights RPCerr"
              // 0   1   2    3   4   5      6       7         8
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

    int nx, ny, nch;
    float *dispy;
    float *dispx = iio_read_image_float_split(v[5], &nx, &ny, &nch);
    if (nch > 1) dispy = dispx + nx*ny;
    else dispy = calloc(nx*ny, sizeof(*dispy));

    float *msk  = iio_read_image_float_split(v[6], &nx, &ny, &nch);
    char *fout_heights  = v[7];
    char *fout_err = v[8];
    float *heightMap = calloc(nx*ny, sizeof(*heightMap));
    float *errMap = calloc(nx*ny, sizeof(*errMap));

    // invert homographies
    double det;
    double invHa[3][3];
    double invHb[3][3];
    INVERT_3X3(invHa, det, Ha);
    INVERT_3X3(invHb, det, Hb);

    // allocate structure for the output data
//    struct world_point *outbuf = malloc(nx * ny * sizeof(*outbuf));

    int npoints = 0;
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int pos = x + nx*y;
            if (msk[pos] <= 0) {
                heightMap[pos] = NAN;
                errMap[pos] = NAN;
            } else {
                double q0[3], q1[3];
//                double groundCoords[2], groundCoordsPlus10[2],
//                       groundCoordsNorm[3];
                double err, h;
                double dx = dispx[pos];
                double dy = dispy[pos];
                double p0[3] = {x, y, 1};
                double p1[3] = {x+dx, y+dy, 1};
                applyHom(q0, invHa, p0);
                applyHom(q1, invHb, p1);

                // compute the coordinates
                h = rpc_height(rpca, rpcb, q0[0], q0[1], q1[0], q1[1], &err);
                heightMap[pos] = h;
                errMap[pos] = err;

//                eval_rpc(groundCoords, rpca, q0[0], q0[1], h);
//                // compute normal
//                eval_rpc(groundCoordsPlus10, rpca, q0[0], q0[1], h+10);
//                groundCoordsNorm[0] = groundCoordsPlus10[0] - groundCoords[0];
//                groundCoordsNorm[1] = groundCoordsPlus10[1] - groundCoords[1];
//                groundCoordsNorm[2] = h+10 - h;
//                normalize_vector_3d(groundCoordsNorm);
//
//                // mercator conversion
//                double mercator[2], mercatorPlus10[2], mercatorNorm[3],
//                       lon0 = 0.0; // reference longitude
//                convert_geodetic_to_mercator(mercator, groundCoords, lon0);
//                // compute normal
//                convert_geodetic_to_mercator(mercatorPlus10, groundCoordsPlus10, lon0);
//                mercatorNorm[0] = mercator[0] - mercatorPlus10[0];
//                mercatorNorm[1] = mercator[1] - mercatorPlus10[1];
//                mercatorNorm[2] = h+10 -h;
//                normalize_vector_3d(mercatorNorm);
//
//                // relief exageration:
//                // 1   --> no exageration
//                // 0.1 --> x10 factor
//                double reliefExagerationFactor = 1;
//
//                outbuf[npoints].lon          = groundCoords[0];
//                outbuf[npoints].lat          = groundCoords[1];
//                outbuf[npoints].h            = h;
//                outbuf[npoints].normal_lon   = groundCoordsNorm[0];
//                outbuf[npoints].normal_lat   = groundCoordsNorm[1];
//                outbuf[npoints].normal_h     = groundCoordsNorm[2];
//
//                outbuf[npoints].x_mercator   = mercator[0];
//                outbuf[npoints].y_mercator   = mercator[1];
//                outbuf[npoints].exagerated_h = h/reliefExagerationFactor;
//                outbuf[npoints].normal_x_mercator = mercatorNorm[0];
//                outbuf[npoints].normal_y_mercator = mercatorNorm[1];
//                outbuf[npoints].normal_h_mercator = mercatorNorm[2];
//
//                outbuf[npoints].rpc_error    = err;
//                outbuf[npoints].x            = x;
//                outbuf[npoints].y            = y;

            }
        }
    }
    // save the height map and error map
    iio_save_image_float_vec(fout_heights, heightMap, nx, ny, 1);
    iio_save_image_float_vec(fout_err, errMap, nx, ny, 1);
    return 0;
}

int main(int c, char *v[])
{
    return main_disp_to_h(c, v);
}
