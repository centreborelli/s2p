// take a series of ply files and produce a digital elevation map

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "lists.c"
#include "fail.c"
#include "xmalloc.c"

#define USE_GDAL // TODO: add an alternative interface (e.g. geotiff)

#ifdef USE_GDAL
#include "gdal.h"
#include "ogr_api.h"
#include "ogr_srs_api.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#endif


struct ply_property {
	enum {UCHAR,FLOAT,DOUBLE} type;
	char name[0x100];
	size_t len;
};

static bool parse_property_line(struct ply_property *t, char *buf)
{
	char typename[0x100];
	bool r = 2 == sscanf(buf, "property %s %s\n", typename, t->name);
	if (!r) return r;
	else if (!strcmp(typename, "uchar")) { t->type = UCHAR;  t->len = 1;}
	else if (!strcmp(typename, "float")) { t->type = FLOAT;  t->len = 4;}
	else if (!strcmp(typename, "double")){ t->type = DOUBLE; t->len = 8;}
	else fail("Unknown property type: %s\n", buf);
	return r;
}



// fast forward "f" until "end_header" is found
// returns the number of 'properties'
// the array of structures *t, contains the names and sizes
// the properties in bytes, isbin is set if binary encoded
// and reads the utm zone
static size_t header_get_record_length_and_utm_zone(FILE *f_in, char *utm,
		int *isbin, struct ply_property *t)
{
	size_t n = 0;
	*isbin = 0;

	char buf[FILENAME_MAX] = {0};
	while (fgets(buf, FILENAME_MAX, f_in)) {
		if (0 == strcmp(buf, "format binary_little_endian 1.0\n"))
			*isbin=1;
		else if (0 == strcmp(buf, "format ascii 1.0\n")) *isbin=0;
		else {
			if (parse_property_line(t+n, buf))
				n += 1;
			else if (0 == strncmp(buf, "comment projection:", 19)) {
				sscanf(buf,"comment projection: UTM %3s\n",utm);
			}
		}
		if (0 == strcmp(buf, "end_header\n"))
			break;
	}
	return n;
}

static void update_min_max(double *min, double *max, double x)
{
	if (x < *min) *min = x;
	if (x > *max) *max = x;
}

// rescale a double:
// min: start of interval
// resolution: spacing between values
static int rescale_double_to_int(double x, double min, double resolution)
{
	int r = floor( (x - min) / resolution);
	return r;
}

static float recenter_double(int i, double xmin, double resolution)
{
	return xmin + resolution * (0.5 + i);
}

struct accumulator_image {
	float *min;
	float *max;
	float *cnt;
	float *avg;
	int w, h;
};

// update the output images with a new height
static void accumulate_height(
		struct accumulator_image *x,
		int i, int j,         // position of the new height
		double v,             // new height
		float weight,         // relative weight
		bool updateminmax     // whether to update min,max fields
		                      // (only makes sense when radius=1)
		)
{
	uint64_t k = (uint64_t) x->w * j + i;
	if (updateminmax) {
		x->min[k] = fmin(x->min[k], v);
		x->max[k] = fmax(x->max[k], v);
	}
	x->avg[k] = (v * weight + x->cnt[k] * x->avg[k]) / (weight + x->cnt[k]);
	x->cnt[k] += weight;
}

size_t get_record(FILE *f_in, int isbin, struct ply_property *t, int n, double *data){
	size_t rec = 0;
	if(isbin) {
		for (int i = 0; i < n; i++) {
			switch(t[i].type) {
			case UCHAR: {
				unsigned char X;
				size_t r = fread(&X, 1, 1, f_in);
				if (r != 1)
					return rec;
				rec += r;
				data[i] = X;
				break; }
			case FLOAT: {
				float X;
				size_t r = fread(&X, sizeof(float), 1, f_in);
				if (r != 1)
					return rec;
				rec += r;
				data[i] = X;
				break; }
			case DOUBLE: {
				 double X;
				 size_t r = fread(&X, sizeof(double), 1, f_in);
				 if (r != 1)
					 return rec;
				 rec += r;
				 data[i] = X;
				 break; }
			}
		}
	} else {
		int i=0;
		while (i < n && !feof(f_in)) {
			rec += fscanf(f_in,"%lf", &data[i]);  i++;
		}
	}
	return rec;
}

// open a ply file, read utm zone in the header, and update the known extrema
static void parse_ply_points_for_extrema(double *xmin, double *xmax,
		double *ymin, double *ymax,
		char *utm, char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return;
	}

	int isbin=0;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	//fprintf(stderr, "%d\n", n);
	//fprintf(stderr, "%s\n", utm);

	double data[n];
	while ( n == get_record(f, isbin, t, n, data) ) {
		update_min_max(xmin, xmax, data[0]);
		update_min_max(ymin, ymax, data[1]);
	}
	fclose(f);
}

// check whether a point is inside the image domain
static int insideP(int w, int h, int i, int j)
{
	return i>=0 && j>=0 && i<w && j<h;
}

static float distance_weight(float sigma, float d)
{
	if (isinf(sigma))
		return 1;
	else
		return exp(-d*d/(2*sigma*sigma));
}

// open a ply file, and accumulate its points to the grand image
static void accumulate_ply_points_to_images(
		struct accumulator_image *x,
		double xmin,  double ymax,  // origin of output grid
		double R,                   // resolution of output grid
		char utm_zone[3],
		char *fname,
		int col_idx,
		int radius,
		double sigma)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return;
	}

	// check that the utm zone is the same as the provided one
	char utm[5];
	int isbin=1;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	if (0 != strncmp(utm_zone, utm, 3))
		fprintf(stderr, "error: different UTM zones among ply files\n");

	if (col_idx < 2 || col_idx > n-1)
		exit(fprintf(stderr, "error: bad col_idx %d\n", col_idx));

	double data[n];
	double sigma2mult2 = 2*sigma*sigma;
	bool updatemmx = radius == 0;

	while ( n == get_record(f, isbin, t, n, data) )
	{
		int i = rescale_double_to_int(data[0], xmin, R);
		int j = rescale_double_to_int(-data[1], -ymax, R);

		for (int k1 = -radius; k1 <= radius; k1++)
		for (int k2 = -radius; k2 <= radius; k2++) {
			int ii = i + k1;
			int jj = j + k2;
			float dist_x = data[0] - recenter_double(ii, xmin, R);
			float dist_y = data[1] - recenter_double(jj, ymax, -R);
			float dist = hypot(dist_x, dist_y);
			float weight = distance_weight(sigma, dist);

			if (insideP(x->w, x->h, ii, jj)) {
				if (col_idx == 2) {
					accumulate_height(x, ii, jj,
						data[2], weight, updatemmx);
					assert(isfinite(data[2]));
				}
				else {
					unsigned int rgb = data[col_idx];
					accumulate_height(x, ii, jj,
						rgb, weight, updatemmx);
				}
			}
		}
	}
	fclose(f);
}

static void save_output_image_with_utm_georeferencing(
		char *filename,
		float *x, int w, int h,                    // data to save
		double g_xoff, double g_xdx, double g_xdy, // geo-transform
		double g_yoff, double g_ydx, double g_ydy,
		char *utm                                  // utm zone string
		)
{
#ifdef USE_GDAL
	GDALAllRegister();
	char **papszOptions = NULL;
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALDatasetH hDstDS = GDALCreate( hDriver, filename,
					  w, h, 1, GDT_Float32,
					  papszOptions );

	double adfGeoTransform[6] = {g_xoff,g_xdx,g_xdy, g_yoff,g_ydx,g_ydy};
	OGRSpatialReferenceH hSRS;
	char *pszSRS_WKT = NULL;
	GDALRasterBandH hBand;
	GDALSetGeoTransform( hDstDS, adfGeoTransform );
	hSRS = OSRNewSpatialReference( NULL );
	char utmNumber[3] = {utm[0], utm[1], '\0'};
	int nZone = atoi(utmNumber);
	int bNorth = (utm[2] == 'N');
	OSRSetUTM( hSRS, nZone, bNorth );
	OSRSetWellKnownGeogCS( hSRS, "WGS84" );
	OSRExportToWkt( hSRS, &pszSRS_WKT );
	OSRDestroySpatialReference( hSRS );
	GDALSetProjection( hDstDS, pszSRS_WKT );
	CPLFree( pszSRS_WKT );
	hBand = GDALGetRasterBand( hDstDS, 1 );
	int r = GDALRasterIO( hBand, GF_Write, 0, 0, w, h,
			  x, w, h, GDT_Float32,
			  0, 0 );
	if (r != 0)
		fprintf(stderr, "ERROR: cannot write %s\n", filename);
	GDALClose( hDstDS );
#endif//USE_GDAL
}


void help(char *s)
{
	fprintf(stderr, "usage:\n\t"
			"ls files | %s [-c column] "
			"[-srcwin \"xoff yoff xsize ysize\"] "
			"[-radius 0] [-sigma resolution] "
			"resolution out.tif\n", s);
	fprintf(stderr, "\t the resolution is in meters per pixel\n");
}

#include "pickopt.c"

int main(int c, char *v[])
{
	int col_idx = atoi(pick_option(&c, &v, "c", "2"));
	char *srcwin = pick_option(&c, &v, "srcwin", "");
	int radius = atoi(pick_option(&c, &v, "radius", "0"));
	float sigma = atof(pick_option(&c, &v, "sigma", "inf"));

	// process input arguments
	if (c != 3) {
		help(*v);
		return 1;
	}

	double resolution = atof(v[1]);
	char *filename_out = v[2];

	// initialize x, y extrema values
	double xmin = INFINITY;
	double xmax = -INFINITY;
	double ymin = INFINITY;
	double ymax = -INFINITY;

	// Subwindow from the source
	double xoff, yoff;
	int xsize, ysize;

	// process each filename from stdin to determine x, y extrema and store
	// the filenames in a list of strings, to be able to open the files
	// again
	char fname[FILENAME_MAX], utm[5];
	struct list *l = NULL;
	while (fgets(fname, FILENAME_MAX, stdin))
	{
		strtok(fname, "\n");
		l = push(l, fname);
		parse_ply_points_for_extrema(&xmin, &xmax, &ymin, &ymax,
				utm, fname);
	}

	if (0 != strcmp(srcwin, "") ) {
		sscanf(srcwin, "%lf %lf %d %d", &xoff, &yoff, &xsize, &ysize);
	}
	else {
		xsize = 1 + floor((xmax - xmin) / resolution);
		ysize = 1 + floor((ymax - ymin) / resolution);
		xoff = (xmax + xmin - resolution * xsize) / 2;
		yoff = (ymax + ymin + resolution * ysize) / 2;
	}
	fprintf(stderr, "xmin: %lf, xmax: %lf, ymin: %lf, ymax: %lf\n", xmin, xmax, ymin, ymax);
	fprintf(stderr, "xoff: %lf, yoff: %lf, xsize: %d, ysize: %d\n", xoff, yoff, xsize, ysize);

	// allocate and initialize accumulator
	struct accumulator_image x[1];
	x->w = xsize;
	x->h = ysize;
	x->min = xmalloc(xsize*ysize*sizeof(float));
	x->max = xmalloc(xsize*ysize*sizeof(float));
	x->cnt = xmalloc(xsize*ysize*sizeof(float));
	x->avg = xmalloc(xsize*ysize*sizeof(float));
	for (uint64_t i = 0; i < (uint64_t) xsize*ysize; i++)
	{
		x->min[i] = INFINITY;
		x->max[i] = -INFINITY;
		x->cnt[i] = 0;
		x->avg[i] = 0;
	}

	// process each filename to accumulate points in the dem
	while (l != NULL)
	{
		// printf("FILENAME: \"%s\"\n", l->current);
		accumulate_ply_points_to_images(x, xoff, yoff, resolution, utm,
				l->current, col_idx, radius, sigma);
		l = l->next;
	}

	// set unknown values to NAN
	for (uint64_t i = 0; i < (uint64_t) xsize*ysize; i++)
		if (!x->cnt[i])
			x->avg[i] = NAN;

	// save output image
	save_output_image_with_utm_georeferencing(
			filename_out,
			x->avg, x->w, x->h,
			xoff, resolution, 0, yoff, 0, -resolution, utm
			);

	// cleanup and exit
	free(x->min);
	free(x->max);
	free(x->cnt);
	free(x->avg);
	return 0;
}
