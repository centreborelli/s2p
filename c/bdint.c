// lowest neighbor interpolation
// optionally: highest-neighbor and average-neighbor
// todo: IDW/shepard of all the neighbors



#include <assert.h>
#include <stdlib.h>
#include <math.h>


#define xmalloc malloc
#include "abstract_dsf.c"

// check wether a pixel is inside the image domain
static int insideP(int w, int h, int i, int j)
{
	return i >= 0 && j >= 0 && i < w && j < h;
}

// type of an accumulator function (like fmin, fmax or sum)
typedef float (accumulator_t)(float, float);

// the "sum" function (gives 0 for INFINITY-INFINITY)
static float sumf(float x, float y)
{
	if (!isfinite(x) && !isfinite(y))
		return 0;
	return x + y;
}

// build a DSF of the NAN holes from image x
static void fill_nan_reps(int *rep, float *x, int w, int h)
{
	adsf_begin(rep,w*h);

	// remove from dsf pixels with known values in input
	for (int i = 0; i < w*h; i++)
		if (!isnan(x[i]))
			rep[i] = -1;

	// join neighboring NANs
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		int p0 = j*w + i;
		int p1 = j*w + i+1;
		int p2  = (j+1)*w + i;
		if (isnan(x[p0]) && isnan(x[p1]))
			adsf_union(rep, w*h, p0, p1);
		if (isnan(x[p0]) && isnan(x[p2]))
			adsf_union(rep, w*h, p0, p2);
	}

	// canonicalize dsf (after this, the DSF is not changed anymore)
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			rep[i] = adsf_find(rep, w*h, i);
}

// fill-in holes by accumulating values at the boundary of each hole
void bdint_gen(float *x, int w, int h, accumulator_t *a)
{
	// create the dsf
	int *rep = xmalloc(w*h*sizeof*rep);
	fill_nan_reps(rep, x, w, h);

	// initialize the accumulators
	float *acc_v = xmalloc(w*h*sizeof*acc_v); // accumulator for values
	float *acc_n = xmalloc(w*h*sizeof*acc_n); // value count
	for (int i = 0; i < w*h; i++)
	if (rep[i] == i)
	{
		acc_v[i] = -a(INFINITY,-INFINITY);
		acc_n[i] = 0;
	}

	// for each value that has a NAN neighbor, update the accumulators
	int n[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int k = 0; k < 4; k++)
	{
		int ii = i + n[k][0];
		int jj = j + n[k][1];
		int ij = j  * w + i;  // point just outside the hole
		int IJ = jj * w + ii; // point just inside the hole
		if (insideP(w, h, i, j) && insideP(w, h, ii, jj) &&
				rep[ij] == -1 && rep[IJ] != -1)
		{
			int ridx = rep[IJ];   // representative of the hole
			acc_v[ridx] = a(acc_v[ridx], x[ij]);
			acc_n[ridx] = acc_n[ridx] + 1;
		}
	}

	// fill-in using the computed value
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			x[i] = acc_v[rep[i]] / (a==sumf?acc_n[rep[i]]:1);

	//cleanup
	free(rep);
	free(acc_v);
	free(acc_n);
}

void bdint_gen_split(float *x, int w, int h, int pd, accumulator_t *a)
{
	for (int l = 0; l < pd; l++)
		bdint_gen(x + w*h*l, w, h, a);
}


#ifndef OMIT_BDINT_MAIN
#define USE_BDINT_MAIN
#endif

#ifdef USE_BDINT_MAIN
#include <stdio.h>
#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	char *opt_a = pick_option(&c, &v, "a", "min");
	char *filename_mask = pick_option(&c, &v, "m", "");
	int help_argument = (int)pick_option(&c, &v, "h", 0);
	if (help_argument || (c != 1 && c != 2 && c != 3)) {
		fprintf(stderr, "usage:\n\t"
			"%s [-a {min|max|avg}] [in.tiff [out.tiff]]\n", *v);
		//        0                     1        2
		return 1;
	}
	char *filename_in   = c > 1 ? v[1] : "-";
	char *filename_out  = c > 2 ? v[2] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);

	if (filename_mask && *filename_mask) {
		int mw, mh;
		float *m = iio_read_image_float(filename_mask, &mw, &mh);
		for (int l = 0; l < pd; l++)
		for (int i = 0; i < mw*mh; i++)
			if (i < w*h && m[i])
				x[l*w*h+i] = NAN;
		free(m);
	}

	accumulator_t *a = fminf;
	if (strstr(opt_a, "ma")) a = fmaxf;
	if (strstr(opt_a, "me") || strstr(opt_a, "av") ) a = sumf;

	bdint_gen_split(x, w, h, pd, a);

	iio_save_image_float_split(filename_out, x, w, h, pd);

	return 0;
}
#endif//USE_BDINT_MAIN
