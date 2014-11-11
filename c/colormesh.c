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
#include "parsenumbers.c"
#include "pickopt.c"


void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);


static void getxyz(double xyz[3], struct rpc *r, double i, double j, double h,
        int zone)
{
    double lon_lat[2];
    eval_rpc(lon_lat, r, i, j, h);
    utm_alt_zone(xyz, lon_lat[1], lon_lat[0], zone);
    xyz[2] = h;
    //printf("%f %f\n", xyz[0], xyz[1]);
}


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


static void help(char *s)
{
    fprintf(stderr, "\t usage: %s out.ply heights.tif rpc.xml "
            "[colors.png] [-h \"h1 ... h9\"] [--utm-zone ZONE] "
            "[--offset_x x0] [--offset_y y0] [--with-normals] [--ascii]\n", s);

    // offset allows the user to choose the origin of the coordinates system,
    // in order to avoid visualisation problems due to huge values of the
    // coordinates (for which float precision is often not enough)
}


int main(int c, char *v[])
{
    if (c < 4 || c > 15) {
        help(*v);
        return 1;
    }

    // with_normals and ascii flags
    bool normals = pick_option(&c, &v, "-with-normals", NULL);
    bool ascii   = pick_option(&c, &v, "-ascii", NULL);

    // offset
    char *offset_x = pick_option(&c, &v, "-offset_x", "0");
    char *offset_y = pick_option(&c, &v, "-offset_y", "0");
    int x0 = atoi(offset_x);
    int y0 = atoi(offset_y);
    bool there_is_an_offset = x0 != 0 || y0 != 0;

    // utm zone and hemisphere: true for 'N' and false for 'S'
    int zone;
    bool hem;
    char *utm_string = pick_option(&c, &v, "-utm-zone", "no_utm_zone");
    parse_utm_string(&zone, &hem, utm_string);

    // rectifying homography. If not provided, it is identity (ie full images)
    char *hom_string = pick_option(&c, &v, "h", "");
    bool there_is_a_homography = *hom_string;
    double inv_hom[9];
    if (there_is_a_homography) {
        int n_hom;
        double *hom = alloc_parse_doubles(9, hom_string, &n_hom);
        if (n_hom != 9)
            fail("can not read 3x3 matrix from \"%s\"", hom_string);
        invert_homography(inv_hom, hom);
    }

    // parse the remaining arguments
    char *fname_ply = v[1];
    char *fname_heights = v[2];
    char *fname_rpc = v[3];
    char *fname_colors = NULL;
    bool there_is_color = c > 4;
    if (there_is_color)
        fname_colors = v[4];

    // read input images
    int w, h, pd;
    float *height = iio_read_image_float(fname_heights, &w, &h);
    uint8_t *color = NULL;
    if (there_is_color) {
        int ww, hh;
        color = iio_read_image_uint8_vec(fname_colors, &ww, &hh, &pd);
        if (w != ww || h != hh) fail("color and height image size mismatch");
        if (pd != 1 && pd != 3) fail("expecting a gray or color image");
    }

    // read rpc
    struct rpc r[1];
    read_rpc_file_xml(r, fname_rpc);

    // count number of valid pixels
    int npoints = 0;
    printf("counting valid points...\r");
    for (int row = 0; row < h; row++)
    for (int col = 0; col < w; col++) {
        uint64_t pix = (uint64_t) row * w + col;
        if (!isnan(height[pix])) {
            npoints++;

            // UTM Zone will be the zone of first 'not NaN' point
            if (zone < 0) {
                double xy[2] = {col, row};
                if (there_is_a_homography)
                    apply_homography(xy, inv_hom, xy);
                double lon_lat[2];
                eval_rpc(lon_lat, r, xy[0], xy[1], height[pix]);
                utm_zone(&zone, &hem, lon_lat[1], lon_lat[0]);
            }
        }
    }
    printf("found %06d valid points\n", npoints);

    // print header for ply file
    FILE *ply_file = fopen(fname_ply, "w");
    write_ply_header(ply_file, ascii, npoints, zone, hem, there_is_color,
            normals);

    // loop over all the pixels of the input height map
    // a 3D point is produced for each 'non Nan' height
    for (int row = 0; row < h; row++)
    for (int col = 0; col < w; col++) {
        if (row % 1000 == 0)
            printf("processing row %06d...\r", row);
        uint64_t pix = (uint64_t) row * w + col;
        if (!isnan(height[pix])) {

            // compute coordinates of pix in the big image
            double xy[2] = {col, row};
            if (there_is_a_homography)
                apply_homography(xy, inv_hom, xy);

            // compute utm coordinates
            double xyz[3], nrm[3], tmp[3];
            getxyz(xyz, r, xy[0], xy[1], height[pix], zone);

            // compute the normal (unit 3D vector with direction of the camera)
            if (normals) {
                getxyz(tmp, r, xy[0], xy[1], height[pix] + 10, zone);
                nrm[0] = tmp[0] - xyz[0];
                nrm[1] = tmp[1] - xyz[1];
                nrm[2] = tmp[2] - xyz[2];
                normalize_vector_3d(nrm);
            }
            if (there_is_an_offset) {
                xyz[0] -= x0;
                xyz[1] -= y0;
            }

            // colorization: if greyscale, copy the grey level on each channel
            uint8_t rgb[3];
            if (there_is_color) {
                for (int k = 0; k < pd; k++) rgb[k] = color[k + pd*pix];
                for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
            }

            // write to ply
            if (ascii) {
               fprintf(ply_file, "%.16f %.16f %.16f ", xyz[0], xyz[1], xyz[2]);
               if (normals)
                   fprintf(ply_file, "%.1f %.1f %.1f ", nrm[0], nrm[1], nrm[2]);
               fprintf(ply_file, "%d %d %d\n", rgb[0], rgb[1], rgb[2]);
            } else {
               float X[3] = {xyz[0], xyz[1], xyz[2]};
               fwrite(X, sizeof(float), 3, ply_file);
               if (normals) {
                   float N[3] = {nrm[0], nrm[1], nrm[2]};
                   fwrite(N, sizeof(float), 3, ply_file);
               }
               if (there_is_color)
                   fwrite(rgb, sizeof(uint8_t), 3, ply_file);
            }
        }
    }

    fclose(ply_file);
    return 0;
}
