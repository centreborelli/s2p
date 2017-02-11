#include <unistd.h>
#include <ctype.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef m_pi
#define m_pi 3.14159265358979323846264338328
#endif

#include "iio.h"
#include "fail.c"
#include "rpc.h"
#include "parsenumbers.c"
#include "pickopt.c"

#undef USE_TIMING
#include "timing.h"


void utm_alt_zone(double *out, double lat, double lon, int zone);
void utm_zone(int *zone, bool *northp, double lat, double lon);

static void lonlat_from_ijh(double ll[2],
        struct rpc *r, double i, double j, double h)
{
    eval_rpc(ll, r, i, j, h);
}

void utm_from_lonlat_and_zone(double xy[2], double lon, double lat, int zone)
{
    utm_alt_zone(xy, lon, lat, zone);
}

static void getxyz(double xyz[3], struct rpc *r, double i, double j, double h,
        int zone)
{
    double ll[3];
    lonlat_from_ijh(ll, r, i, j, h);
    utm_from_lonlat_and_zone(xyz, ll[1], ll[0], zone);
    xyz[2] = h;
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


void write_ply_header(FILE* f, uint64_t npoints, int zone,
        bool hem, bool colors, bool normals)
{
    fprintf(f, "ply\n");
    fprintf(f, "format binary_little_endian 1.0\n");
    fprintf(f, "comment created by S2P\n");
    if (zone >= 0)
        fprintf(f, "comment projection: UTM %02d%s\n", zone, (hem ? "N" : "S"));
    fprintf(f, "element vertex %" PRIu64 "\n", npoints);
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
    fprintf(f, "end_header\n");
}

static void parse_utm_string(int *zone, bool *hem, char *s)
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

static void help(char *s)
{
    fprintf(stderr, "\t usage: %s out.ply heights.tif rpc.xml "
            "[colors.png] [-h \"h1 ... h9\"] [--utm-zone ZONE] "
            "[--offset_x x0] [--offset_y y0] [--with-normals] "
        "[--lon-m l0] [--lon-M lf] [--lat-m l0] [--lat-M lf] "
        "[--lonlat-clip polygon.kml] [--max-memory MB]\n", s);

    // offset allows the user to choose the origin of the coordinates system,
    // in order to avoid visualisation problems due to huge values of the
    // coordinates (for which float precision is often not enough)
    // NOTE: now we use double, so the "offset" is unnecessary
}


int main(int c, char *v[])
{
    if (c < 4 || c > 17) {
        help(*v);
        return 1;
    }

    // with_normals flag
    bool normals = pick_option(&c, &v, "-with-normals", NULL);

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

    // longitude-latitude bounding box
    double lon_m = atof(pick_option(&c, &v, "-lon-m", "-inf"));
    double lon_M = atof(pick_option(&c, &v, "-lon-M", "inf"));
    double lat_m = atof(pick_option(&c, &v, "-lat-m", "-inf"));
    double lat_M = atof(pick_option(&c, &v, "-lat-M", "inf"));

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

    // max memory
    int nthreads = 1;
    #ifdef _OPENMP
    nthreads = omp_get_max_threads();
    #endif
    char *max_mem = pick_option(&c, &v, "-max-memory", "1024");
    uint64_t buf_size = 1024*1024*atoi(max_mem) / nthreads;

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
    uint64_t npoints = 0;
    printf("counting valid points...\n");
    TIMING_CPUCLOCK_START(0);
    for (uint64_t pix = 0; pix < (uint64_t) w*h; pix++) {
        if (!isnan(height[pix])) {
            // compute coordinates of pix in the big image
            int col = pix % w;
            int row = pix / w;
            double xy[2] = {col, row};
            if (there_is_a_homography)
                apply_homography(xy, inv_hom, xy);

            // reject points outside lonlat bbx
            double ll[2];
            lonlat_from_ijh(ll, r, xy[0], xy[1], height[pix]);
            if (ll[0]<lon_m || ll[0]>lon_M || ll[1]<lat_m || ll[1]>lat_M)
                continue;

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
    TIMING_CPUCLOCK_TOGGLE(0);
    TIMING_PRINTF("CPU time spent counting points: %0.6fs\n", TIMING_CPUCLOCK_S(0));
    printf("found %" PRIu64 " valid points\n", npoints);

    // print header for ply file
    FILE *ply_file = fopen(fname_ply, "w");
    write_ply_header(ply_file, npoints, zone, hem, there_is_color,
            normals);

    // prepare one buffer per thread
    char *buf[nthreads], *ptr[nthreads];
    for (int i = 0; i < nthreads; i++) {
        buf[i] = malloc(buf_size);
        ptr[i] = buf[i];
    }

    // loop over all the pixels of the input height map
    // a 3D point is produced for each 'non Nan' height
    TIMING_WALLCLOCK_RESET(0);
    TIMING_WALLCLOCK_TOGGLE(0);
    size_t point_size = 3*sizeof(double);
    if (normals)
        point_size += 3*sizeof(double);
    if (there_is_color)
        point_size += 3*sizeof(uint8_t);

#ifdef _OPENMP
    # pragma omp parallel for
#endif
    for (uint64_t pix = 0; pix < (uint64_t) w*h; pix++) {
        if (!isnan(height[pix])) {

            // if the buffer is full, write to file and reset
            int i = 0;
            #ifdef _OPENMP
            i = omp_get_thread_num();
            #endif
            if (buf[i] + buf_size - ptr[i] < point_size) {
                int nbytes = ptr[i] - buf[i];
#ifdef _OPENMP
                # pragma omp critical
#endif
                    fwrite(buf[i], sizeof(char), nbytes, ply_file);
                ptr[i] = buf[i];
            }

            // compute coordinates of pix in the big image
            int col = pix % w;
            int row = pix / w;
            double xy[2] = {col, row};
            if (there_is_a_homography)
                apply_homography(xy, inv_hom, xy);

            // compute utm coordinates (rejecting points outside lonlat bbx)
            double xyz[3], nrm[3], tmp[3], ll[2];
            //getxyz(xyz, r, xy[0], xy[1], height[pix], zone);
            lonlat_from_ijh(ll, r, xy[0], xy[1], height[pix]);
            if (ll[0]<lon_m || ll[0]>lon_M || ll[1]<lat_m || ll[1]>lat_M)
                continue;
            utm_from_lonlat_and_zone(xyz, ll[1], ll[0], zone);
            xyz[2] = height[pix];

            if (there_is_an_offset) {
                xyz[0] -= x0;
                xyz[1] -= y0;
            }

            // write to memory
            double *ptr_double = (double *) ptr[i];
            ptr_double[0] = xyz[0];
            ptr_double[1] = xyz[1];
            ptr_double[2] = xyz[2];

            // normals (unit 3D vector with direction of the camera)
            if (normals) {
                getxyz(tmp, r, col, row, height[pix] + 10, zone);
                nrm[0] = tmp[0] - xyz[0];
                nrm[1] = tmp[1] - xyz[1];
                nrm[2] = tmp[2] - xyz[2];
                normalize_vector_3d(nrm);
                ptr_double[3] = nrm[0];
                ptr_double[4] = nrm[1];
                ptr_double[5] = nrm[2];
            }

            // colorization: if greyscale, copy the grey on each channel
            uint8_t rgb[3];
            if (there_is_color) {
                for (int k = 0; k < pd; k++) rgb[k] = color[k + pd*pix];
                for (int k = pd; k < 3; k++) rgb[k] = rgb[k-1];
                char *ptr_char = ptr[i] + 3*sizeof(double);
                ptr_char[0] = rgb[0];
                ptr_char[1] = rgb[1];
                ptr_char[2] = rgb[2];
            }

            ptr[i] += point_size;
        }
    }
    TIMING_WALLCLOCK_TOGGLE(0);
    TIMING_PRINTF("WALL time spent computing the points: %0.6fs\n",
            TIMING_WALLCLOCK_S(0));

    // write the last buffers to file
    for (int i = 0; i < nthreads; i++)
        fwrite(buf[i], sizeof(char), ptr[i] - buf[i], ply_file);

    fclose(ply_file);
    return 0;
}
