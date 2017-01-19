// take a series of ply files and produce a digital elevation map

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <geotiff/xtiffio.h>
#include <geotiff/geotiffio.h>
#include <geotiff/geo_tiffp.h>

#include "../3rdparty/iio/iio.h"
#include "lists.c"
#include "fail.c"
#include "xmalloc.c"


// convert string like '28N' into a number like 32628, according to:
// WGS84 / UTM northern hemisphere: 326zz where zz is UTM zone number
// WGS84 / UTM southern hemisphere: 327zz where zz is UTM zone number
// http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.3.3.1
static int get_utm_zone_index_for_geotiff(char *utm_zone)
{
	int out = 32000;
	if (utm_zone[2] == 'N')
		out += 600;
	else if (utm_zone[2] == 'S')
		out += 700;
	else
		fprintf(stderr, "error: bad utm zone value: %s\n", utm_zone);
	utm_zone[2] = '\0';
	out += atoi(utm_zone);
	return out;
}

void set_geotif_header(char *tiff_fname, char *utm_zone, float xoff,
		float yoff, float scale)
{
	// open tiff file
	TIFF *tif = XTIFFOpen(tiff_fname, "r+");
	if (!tif)
		fail("failed in XTIFFOpen\n");

	GTIF *gtif = GTIFNew(tif);
	if (!gtif)
		fail("failed in GTIFNew\n");

	// write TIFF tags
	double pixsize[3] = {scale, scale, 0.0};
	TIFFSetField(tif, GTIFF_PIXELSCALE, 3, pixsize);

	double tiepoint[6] = {0.0, 0.0, 0.0, xoff, yoff, 0.0};
	TIFFSetField(tif, GTIFF_TIEPOINTS, 6, tiepoint);

	// write GEOTIFF keys
	int utm_ind = get_utm_zone_index_for_geotiff(utm_zone);
	GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, utm_ind);
	GTIFWriteKeys(gtif);

	// free and close
	GTIFFree(gtif);
	XTIFFClose(tif);
}



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
	else if (0 == strcmp(typename, "uchar")) { t->type = UCHAR;  t->len = 1;}
	else if (0 == strcmp(typename, "float")) { t->type = FLOAT;  t->len = 4;}
	else if (0 == strcmp(typename, "double")){ t->type = DOUBLE; t->len = 8;}
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
		if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) *isbin=1;
		else if (0 == strcmp(buf, "format ascii 1.0\n")) *isbin=0;
		else {
			if (parse_property_line(t+n, buf))
				n += 1;
			else if (0 == strncmp(buf, "comment projection:", 19)) {
				sscanf(buf, "comment projection: UTM %s", utm);
			}
		}
		if (0 == strcmp(buf, "end_header\n"))
			break;
	}
	return n;
}

static void update_min_max(float *min, float *max, float x)
{
	if (x < *min) *min = x;
	if (x > *max) *max = x;
}

// re-scale a float between 0 and w
static int rescale_float_to_int(double x, double min, double max, int w)
{
	int r = w * (x - min)/(max - min);
	if (r < 0) r = 0;
	if (r >= w) r = w -1;
	return r;
}

struct images {
	float *min;
	float *max;
	float *cnt;
	float *avg;
	int w, h;
};

// update the output images with a new height
static void add_height_to_images(struct images *x, int i, int j, float v)
{
	uint64_t k = (uint64_t) x->w * j + i;
	x->min[k] = fmin(x->min[k], v);
	x->max[k] = fmax(x->max[k], v);
	x->avg[k] = (v + x->cnt[k] * x->avg[k]) / (1 + x->cnt[k]);
	x->cnt[k] += 1;
}

size_t get_record(FILE *f_in, int isbin, struct ply_property *t, int n, double *data){
	size_t rec = 0;
	if(isbin) {
		for (int i = 0; i < n; i++) {
			switch(t[i].type) {
				case UCHAR: {
						    unsigned char X;
						    rec += fread(&X, 1, 1, f_in);
						    data[i] = X;
						    break; }
				case FLOAT: {
						    float X;
						    rec += fread(&X, sizeof(float), 1, f_in);
						    data[i] = X;
						    break; }
				case DOUBLE: {
						     double X;
						     rec += fread(&X, sizeof(double), 1, f_in);
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
static void parse_ply_points_for_extrema(float *xmin, float *xmax, float *ymin,
		float *ymax, char *utm, char *fname)
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

// open a ply file, and accumulate its points to the image
static void add_ply_points_to_images(struct images *x,
		float xmin, float xmax, float ymin, float ymax,
		char utm_zone[3], char *fname, int col_idx)
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

	if (col_idx < 2 || col_idx > 5)
		exit(fprintf(stderr, "error: bad col_idx %d\n", col_idx));


	double data[n];
	while ( n == get_record(f, isbin, t, n, data) ) {

		int i = rescale_float_to_int(data[0], xmin, xmax, x->w);
		int j = rescale_float_to_int(-data[1], -ymax, -ymin, x->h);
		if (col_idx == 2) {
			add_height_to_images(x, i, j, data[2]);
			assert(isfinite(data[2]));
		}
		else
		{
			unsigned int rgb = data[col_idx];
			add_height_to_images(x, i, j, rgb);
		}

	}

	fclose(f);
}


void help(char *s)
{
	fprintf(stderr, "usage:\n\t"
			"ls files | %s [-c column] [-bb \"xmin xmax ymin ymax\"] resolution out.tif\n", s);
	fprintf(stderr, "\t the resolution is in meters per pixel\n");
}

#include "pickopt.c"

int main(int c, char *v[])
{
	int col_idx = atoi(pick_option(&c, &v, "c", "2"));
	const char *bbminmax = pick_option(&c, &v, "bb", "");

	// process input arguments
	if (c != 3) {
		help(*v);
		return 1;
	}
	float resolution = atof(v[1]);
	char *filename_out = v[2];

	// initialize x, y extrema values
	float xmin = INFINITY;
	float xmax = -INFINITY;
	float ymin = INFINITY;
	float ymax = -INFINITY;

	// process each filename from stdin to determine x, y extremas and store the
	// filenames in a list of strings, to be able to open the files again
	char fname[FILENAME_MAX], utm[5];
	struct list *l = NULL;
	while (fgets(fname, FILENAME_MAX, stdin))
	{
		strtok(fname, "\n");
		l = push(l, fname);
		parse_ply_points_for_extrema(&xmin, &xmax, &ymin, &ymax, utm, fname);
	}
	if (0 != strcmp(bbminmax, "") ) {
		sscanf(bbminmax, "%f %f %f %f", &xmin, &xmax, &ymin, &ymax);
	}
	fprintf(stderr, "xmin: %20f, xmax: %20f, ymin: %20f, ymax: %20f\n", xmin, xmax, ymin, ymax);

	// compute output image dimensions
	int w = 1 + (xmax - xmin) / resolution;
	int h = 1 + (ymax - ymin) / resolution;

	// allocate and initialize output images
	struct images x;
	x.w = w;
	x.h = h;
	x.min = xmalloc(w*h*sizeof(float));
	x.max = xmalloc(w*h*sizeof(float));
	x.cnt = xmalloc(w*h*sizeof(float));
	x.avg = xmalloc(w*h*sizeof(float));
	for (uint64_t i = 0; i < (uint64_t) w*h; i++)
	{
		x.min[i] = INFINITY;
		x.max[i] = -INFINITY;
		x.cnt[i] = 0;
		x.avg[i] = 0;
	}

	// process each filename to accumulate points in the dem
	while (l != NULL)
	{
		// printf("FILENAME: \"%s\"\n", l->current);
		add_ply_points_to_images(&x, xmin, xmax, ymin, ymax, utm, l->current, col_idx);
		l = l->next;
	}

	// set unknown values to NAN
	for (uint64_t i = 0; i < (uint64_t) w*h; i++)
		if (!x.cnt[i])
			x.avg[i] = NAN;

	// save output image
	iio_save_image_float(filename_out, x.avg, w, h);
	set_geotif_header(filename_out, utm, xmin, ymax, resolution);

	// cleanup and exit
	free(x.min);
	free(x.max);
	free(x.cnt);
	free(x.avg);
	return 0;
}
