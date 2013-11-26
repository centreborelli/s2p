#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "vvector.h"
#include "iio.h"
#include "read_matrix.c"

#include "svd.c"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


void applyHom(double outv[3], double M[3][3], double v[3]) {
    MAT_DOT_VEC_3X3(outv, M, v);
    outv[0] /= outv[2];
    outv[1] /= outv[2];
    outv[2] /= outv[2];
}




// normalize in place a 3d vector
static void normalize_vector_3d(double vec[3]){
    const int dim=3;
    double norm = 0;
    for(int i=0;i<dim;i++)
        norm += vec[i]*vec[i];
    norm = sqrt(norm);
    for(int i=0;i<dim;i++)
        vec[i] /= norm;
}

static void apply_3x4(double y[3], double H[3][4], double x[4])
{
	y[0] = H[0][0]*x[0] + H[0][1]*x[1] + H[0][2]*x[2] + H[0][3]*x[3];
	y[1] = H[1][0]*x[0] + H[1][1]*x[1] + H[1][2]*x[2] + H[1][3]*x[3];
	y[2] = H[2][0]*x[0] + H[2][1]*x[1] + H[2][2]*x[2] + H[2][3]*x[3];
}


void triangulate_DLT(double P1[3][4], double P2[3][4], double X1[3], double X2[3], double X3D[4], double *error) {
   // Truiangulate the points X1,X2 using the DLT method (Hartley p312)
   X1[0] /= X1[2]; X1[1] /= X1[2]; X1[2] /= X1[2];
   X2[0] /= X2[2]; X2[1] /= X2[2]; X2[2] /= X2[2];

   // build homogeneous system
   double A[4][4];
   for(int i=0;i<4;i++) {
      A[0][i] = X1[0] * P1[2][i] - P1[0][i];
      A[1][i] = X1[1] * P1[2][i] - P1[1][i];
      A[2][i] = X2[0] * P2[2][i] - P2[0][i];
      A[3][i] = X2[1] * P2[2][i] - P2[1][i];
   }

   //solve
   double d[4], u[4][4], v[4][4];
   // compute the decomposition: a = u d v^T
   // d: output size n (singular values)
   // a: input size m*n (m>=n)
   // u: output size m*m
   // v: output size n*n
   // the returned v is not transposed 
   svd_double(d, A[0], u[0], 4, v[0], 4);

   double vmin = INFINITY;
   for(int i=0;i<4;i++){
      if(d[i]<vmin) {
         vmin = d[i];
         X3D[0] = v[0][i]/v[3][i];
         X3D[1] = v[1][i]/v[3][i];
         X3D[2] = v[2][i]/v[3][i];
         X3D[3] = v[3][i]/v[3][i];
      }
   }
   // compute the distance of the projected point 
   // can do the same for P2
   double err[2],pp[3];
   apply_3x4(pp,P1,X3D);
   err[0] = pp[0]/pp[2] - X1[0];
   err[1] = pp[1]/pp[2] - X1[1];
   *error = hypot(err[0],err[1]);


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
};


int main_disp_to_h(int c, char *v[])
{
    if (c != 9) {
        fprintf(stderr, "usage:\n\t"
                "%s Pa  Pb Ha Hb dispAB mskAB out_heights     err"
              // 0   1   2    3   4   5      6       7         8
                "\n", *v);
        return EXIT_FAILURE;
    }

    // read input data
   double Pa[3][4], Pb[3][4];
	read_matrix_3x4(Pa, v[1]);
	read_matrix_3x4(Pb, v[2]);
    double Ha[3][3], Hb[3][3];
    read_matrix(Ha,v[3]);
    read_matrix(Hb,v[4]);

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
    struct world_point *outbuf = malloc(nx*ny*sizeof(*outbuf));

    int npoints = 0;
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int pos = x + nx*y;
            if (msk[pos] <= 0) {
                heightMap[pos] = NAN;
                errMap[pos] = NAN;
            } else {
                double q0[3], q1[3],
                       groundCoords[2], groundCoordsPlus10[2],
                       groundCoordsNorm[3];
                double err, h;
                double dx = dispx[pos];
                double dy = dispy[pos];
                double p0[3] = {x, y, 1};
                double p1[3] = {x+dx, y+dy, 1};
                applyHom(q0, invHa, p0);
                applyHom(q1, invHb, p1);

                // compute the coordinates
                double X3D[4];
                triangulate_DLT(Pa, Pb, q0, q1, X3D, &err);

                // distance to the projective plane. Equation 6.15 (Hartley)
                double M[3][3], detM; 
                for (int i=0;i<3;i++) for (int j=0;j<3;j++) M[i][j] = Pa[i][j];
                DETERMINANT_3X3(detM,M);
                h = ( detM>0 ? 1:-1 ) * ( (Pa[2][0]*X3D[0] + Pa[2][1]*X3D[1] + Pa[2][2]*X3D[2] + Pa[2][3]*X3D[3]) / X3D[3] ) ;

                groundCoords[0] = X3D[0] / X3D[3];
                groundCoords[1] = X3D[1] / X3D[3];

                //  TODO NORMAL?  
//                eval_rpc(groundCoords, rpca, q0[0], q0[1], h);
//                // compute normal
//                eval_rpc(groundCoordsPlus10, rpca, q0[0], q0[1], h+10);
//                groundCoordsNorm[0] = groundCoordsPlus10[0] - groundCoords[0];
//                groundCoordsNorm[1] = groundCoordsPlus10[1] - groundCoords[1];
//                groundCoordsNorm[2] = h+10 - h;
//                normalize_vector_3d(groundCoordsNorm);


                // relief exageration:
                // 1   --> no exageration
                // 0.1 --> x10 factor
                double reliefExagerationFactor = 1;

                outbuf[npoints].lon          = groundCoords[0];
                outbuf[npoints].lat          = groundCoords[1];
                outbuf[npoints].h            = h;
                outbuf[npoints].normal_lon   = groundCoordsNorm[0];
                outbuf[npoints].normal_lat   = groundCoordsNorm[1];
                outbuf[npoints].normal_h     = groundCoordsNorm[2];


                outbuf[npoints].rpc_error    = err;
                outbuf[npoints].x            = x;
                outbuf[npoints].y            = y;

                    heightMap[pos] = h;
                    errMap[pos] = err;
//                if (err < 5) {
//                    heightMap[pos] = h;
//                    errMap[pos] = err;
//                } else {
//                    heightMap[pos] = NAN;
//                    errMap[pos] = NAN;
//                }

            }
        }
    }
    // save the height map and error map
    iio_save_image_float_vec(fout_heights, heightMap, nx,ny, 1);
    iio_save_image_float_vec(fout_err, errMap, nx,ny, 1);
    return 0;
}

int main(int c, char *v[])
{
    return main_disp_to_h(c, v);
}
