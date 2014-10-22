// take a series of ply files and produce a digital elevation map


#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xmalloc.c"
#include "smapa.h"
SMART_PARAMETER(PLY_RECORD_LENGTH, 27)


// fast forward "f" until "lin" is found
static void eat_until_this_line(FILE *f, char *lin)
{
	char buf[FILENAME_MAX] = {0};
	while (fgets(buf, FILENAME_MAX, f))
		if (0 == strcmp(buf, lin))
			return;
}

// re-scale a float between 0 and w
static int rescale_float_to_int(float x, float min, float max, int w)
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
	int k = x->w * j + i;
	x->min[k] = fmin(x->min[k], v);
	x->max[k] = fmax(x->max[k], v);
	x->avg[k] = (v + x->cnt[k] * x->avg[k]) / (1 + x->cnt[k]);
	x->cnt[k] += 1;
}

// open a ply file, and accumulata its points to the image
static void add_ply_points_to_images(struct images *x,
		float xmin, float xmax, float ymin, float ymax,
		char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return;
	}

	eat_until_this_line(f, "end_header\n");

	size_t n = PLY_RECORD_LENGTH();
	char cbuf[n];
	float *fbuf = (void*)cbuf;
	while (n == fread(cbuf, 1, n, f))
	{
		int i = rescale_float_to_int(fbuf[0], xmin, xmax, x->w);
		int j = rescale_float_to_int(fbuf[1], ymin, ymax, x->h);
		//fprintf(stderr, "\t%8.8lf %8.8lf %8.8lf %d %d\n",
		//		fbuf[0], fbuf[1], fbuf[2], i, j);
		add_height_to_images(x, i, j, fbuf[2]);
	}

	fclose(f);
}


#include "iio.h"
int main(int c, char *v[])
{
	// process input arguments
	if (c != 8) {
		fprintf(stderr, "usage:\n\t"
				"ls files|%s x0 xf y0 yf w h out.tiff\n", *v);
		//                         0 1  2  3  4  5 6 7
		return 1;
	}
	float xmin = atof(v[1]);
	float xmax = atof(v[2]);
	float ymin = atof(v[3]);
	float ymax = atof(v[4]);
	int w = atoi(v[5]);
	int h = atoi(v[6]);
	char *filename_out = v[7];

	// allocate and initialize output images
	struct images x;
	x.w = w;
	x.h = h;
	x.min = xmalloc(w*h*sizeof(float));
	x.max = xmalloc(w*h*sizeof(float));
	x.cnt = xmalloc(w*h*sizeof(float));
	x.avg = xmalloc(w*h*sizeof(float));
	for (int i = 0; i < w*h; i++)
	{
		x.min[i] = INFINITY;
		x.max[i] = -INFINITY;
		x.cnt[i] = 0;
		x.avg[i] = 0;
	}

	// process each filename from stdin
	char fname[FILENAME_MAX];
	while (fgets(fname, FILENAME_MAX, stdin))
	{
		strtok(fname, "\n");
		printf("FILENAME: \"%s\"\n", fname);
		add_ply_points_to_images(&x, xmin, xmax, ymin, ymax, fname);
	}

	// set unknown values to NAN
	for (int i = 0; i < w*h; i++)
		if (!x.cnt[i])
			x.avg[i] = NAN;

	// save output image
	iio_save_image_float(filename_out, x.avg, w, h);
	//iio_save_image_float("/tmp/flattened_min.tiff", x.min, w, h);
	//iio_save_image_float("/tmp/flattened_max.tiff", x.max, w, h);
	//iio_save_image_float("/tmp/flattened_cnt.tiff", x.cnt, w, h);
	//iio_save_image_float("/tmp/flattened_avg.tiff", x.avg, w, h);

	// cleanup and exit
	free(x.min);
	free(x.max);
	free(x.cnt);
	free(x.avg);
	return 0;
}
