#include <unistd.h>
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#ifndef m_pi
#define m_pi 3.14159265358979323846264338328
#endif

#include "iio.h"
#include "fail.c"
#include "rpc.h"
#include "read_matrix.c"
#include "pickopt.c"

// static void mercator(double m[2], double x[2])
// {
//     double r = 6378100;
//     double deg = m_pi/180;
//     m[0] = r * x[0] * deg;
//     m[1] = r * log( ( 1 + sin(x[1]*deg) ) / cos(x[1]*deg) );
// }
//
// static void getxyz(double xyz[3], struct rpc *r, double i, double j, double h)
// {
//     double tmp[2];
//     eval_rpc(tmp, r, i, j, h);
//     mercator(xyz, tmp);
//     xyz[2] = h;
// }

void utm(double *out, double lat, double lon);
void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);


static void getxyz(double xyz[3], struct rpc *r, double i, double j, double h,
        int zone)
{
    double lon_lat[2];
    eval_rpc(lon_lat, r, i, j, h);
    if (zone >= 0)
        utm_alt_zone(xyz, lon_lat[1], lon_lat[0], zone);
    else
        utm(xyz, lon_lat[1], lon_lat[0]);
    xyz[2] = h;
    //printf("%f %f\n", xyz[0], xyz[1]);
}


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


void write_ply_header(FILE* f, int npoints, int zone, bool hem, bool normals)
{
    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
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
    fprintf(f, "property uchar red\n");
    fprintf(f, "property uchar green\n");
    fprintf(f, "property uchar blue\n");
    fprintf(f, "end_header\n");
}


unsigned char test_little_endian( void )
{
      int x=1;   return (*(char*)&(x)==1);
}

void write_ply_header_binary(FILE* f, int npoints, int zone, bool hem, bool
        normals)
{
    if (!test_little_endian())
       for (int i = 1; i < 100; i++)
          printf("BINARY PLY NOT SUPPORTED ON BIG ENDIAN SYSTEMS!!\n");
    fprintf(f, "ply\n");
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
    fprintf(f, "property uchar red\n");
    fprintf(f, "property uchar green\n");
    fprintf(f, "property uchar blue\n");
    fprintf(f, "end_header\n");
}

#include "smapa.h"
SMART_PARAMETER_SILENT(IJMESH, 0)
SMART_PARAMETER_SILENT(IJMESHFAC, 2)


// void print_help(char *bin_name) {
//     fprintf(stderr, "usage:\n\t"
//             "%s [-a] colors heights rpc Hfile.txt out.ply [x0 y0]\n", bin_name);
//     //       0         1      2       3   4        5       6  7
//     fprintf(stderr, "\t use -a flag to write an ascii ply (default: binary)\n");
// }

void print_help(char *bin_name)
{
    fprintf(stderr, "usage:\n\t"
            "%s colors heights rpc Hfile.txt out.ply [x0 y0]"
    //       0    1      2       3   4        5       6  7
            " [--with-normals]\n", bin_name);
}


int main(int c, char *v[])
{
    if (c != 6 & c != 7 & c != 8 & c!= 9) {
        print_help(*v);
        return 1;
    }

    // with_normals flag
    bool normals = (pick_option(&c, &v, "-with-normals", NULL) != NULL);

    // binary | ascii flag
    int binary = 1;
//    int k;
//
//    // parse args with getopt to search for the ascii '-a' flag
//    while ((k = getopt (c, v, "a")) != -1)
//        switch (k)
//        {
//            case 'a':
//                binary = 0;
//                printf("ascii\n");
//                break;
//            case '?':
//                if (isprint (optopt))
//                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
//                else
//                    fprintf (stderr, "Unknown option character `\\x%x'.\n",
//                            optopt);
//                print_help(*v);
//                return 1;
//            default:
//                fprintf(stderr, "default case\n");
//                print_help(*v);
//                abort ();
//        }

    // read the other arguments
    // x0 and y0 are meant to allow the user to choose the origin in the
    // mercator coordinates system, in order to avoid visualisation problems
    // due to huge values of the coordinates (for which float precision is
    // often not enough)
//    int i = optind;
    int i = 1;
    char *fname_colors = v[i++];
    char *fname_heights = v[i++];
    char *fname_rpc = v[i++];
    double H[3][3], invH[3][3];
    read_matrix(H, v[i++]);
    FILE *out = fopen(v[i++], "w");
    invert_homography(invH, H);

    int x0 = c > i ? atoi(v[i++]) : 0;
    int y0 = c > i ? atoi(v[i]) : 0;

    int w, h, pd, ww, hh;
    uint8_t *colors = iio_read_image_uint8_vec(fname_colors, &w, &h, &pd);
    float *heights = iio_read_image_float(fname_heights, &ww, &hh);
    if (w != ww || h != hh) fail("color and height image size mismatch");
    if (pd != 1 && pd != 3) fail("expecting a gray or color image");

    struct rpc r[1]; read_rpc_file_xml(r, fname_rpc);

    uint8_t (*color)[w][pd] = (void*)colors;
    float (*height)[w] = (void*)heights;

    int zone = -1;
    bool hem = true;

    // count number of valid pixels
    int npoints = 0;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
            if (!isnan(height[j][i]))
            {
                npoints++;

                // UTM Zone will be the zone of first not nan point
                if (zone < 0)
                {
                    double xy[2] = {i, j}, pq[2];
                    apply_homography(pq, invH, xy);
                    double lon_lat[2];
                    eval_rpc(lon_lat, r, pq[0], pq[1], height[j][i]);
                    utm_zone(&zone, &hem, lon_lat[1], lon_lat[0]);
                }
            }

    // print header for ply file
    if (binary)
      write_ply_header_binary(out, npoints, zone, hem, normals);
    else
      write_ply_header(out, npoints, zone, hem, normals);

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
            // the rpc_eval direct estimation function.
            // Normals (ie unit 3D
            // vector pointing in the direction of the camera) are saved
            // together with the coordinates of the 3D points.
            double xy[2] = {i, j}, pq[2];
            apply_homography(pq, invH, xy);
            double xyz[3] = {pq[1], pq[0], IJMESHFAC() * height[j][i]};
            double nrm[3] = {0, 0, 1};
            if (!IJMESH()) {
                double tmp[3];
                getxyz(xyz, r, pq[0], pq[1], height[j][i], zone);
                if (normals) {
                    getxyz(tmp, r, pq[0], pq[1], height[j][i] + 10, zone);
                    nrm[0] = tmp[0] - xyz[0];
                    nrm[1] = tmp[1] - xyz[1];
                    nrm[2] = tmp[2] - xyz[2];
                    normalize_vector_3d(nrm);
                }
                // translate the x, y coordinates, to avoid
                // visualisation problems
                xyz[0] -= x0;
                xyz[1] -= y0;
            }

            // print the voxel in the ply output file
            if (binary) {
               float X[3] = {xyz[0], xyz[1], xyz[2]};
               fwrite(X, sizeof(float), 3, out);
               if (normals) {
                   float N[3] = {nrm[0], nrm[1], nrm[2]};
                   fwrite(N, sizeof(float), 3, out);
               }
               unsigned char C[3] = {rgb[0], rgb[1], rgb[2]};
               fwrite(C, sizeof(unsigned char), 3, out);
            } else {
               fprintf(out, "%.16f %.16f %.16f ", xyz[0], xyz[1], xyz[2]);
               if (normals)
                   fprintf(out, "%.1f %.1f %.1f ", nrm[0], nrm[1], nrm[2]);
               fprintf(out, "%d %d %d\n", rgb[0], rgb[1], rgb[2]);
            }
        }

    return 0;
}
