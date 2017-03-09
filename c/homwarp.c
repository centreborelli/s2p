
// y = H(x)
static void apply_homography(double y[2], double H[9], double x[2])
{
	double X = H[0]*x[0] + H[1]*x[1] + H[2];
	double Y = H[3]*x[0] + H[4]*x[1] + H[5];
	double Z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = X / Z;
	y[1] = Y / Z;
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

#include <math.h>
#include "bilinear_interpolation.c"
#include "marching_interpolation.c"
#include "bicubic_gray.c"
#include "spline.c"

static
float nearest_neighbor_interpolator(float *x, int w, int h, float p, float q)
{
	int ip = round(p);
	int iq = round(q);
	if (ip < 0) ip = 0;
	if (iq < 0) iq = 0;
	if (ip >= w) ip = w - 1;
	if (iq >= h) iq = h - 1;
	return x[w*iq+ip];
}

typedef float (*gray_interpolator_t)(float *,int,int,float,float);


int homwarp(float *X, int W, int H, double M[9], float *x,
		int w, int h, int o)
{
	gray_interpolator_t u = bicubic_interpolation_gray;
	if (o == 0) u = nearest_neighbor_interpolator;
	if (o == 1) u = marching_interpolation_at;
	if (o == 2) u = bilinear_interpolation_at;

	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		double p[2] = {i, j}, Mp[2];
		apply_homography(p, M, p);
		X[j*W+i] = u(x, w, h, p[0], p[1]);
	}

	return 0;
}

int shomwarp(float *X, int W, int H, double M[9], float *x,
		int w, int h, int o)
{
	// if low order-interpolation, evaluate right away
	if (o == 0 || o == 1 || o == 2 || o == -3)
		return homwarp(X, W, H, M, x, w, h, o);

	// otherwise, pre-filtering is required
	bool r = prepare_spline(x, w, h, 1, o);
	if (!r) return 2;

	// warp the points
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	{
		double p[2] = {i, j};
		apply_homography(p, M, p);
		p[0] += 0.5; // solve a mis-alignement convention
		p[1] += 0.5;
		float *out = X + (j*W + i);
		evaluate_spline_at(out, x, w, h, 1, o, p[0], p[1]);
	}

	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include "xmalloc.c"
#include "parsenumbers.c"
#include "pickopt.c"
int main(int c, char *v[])
{
	bool do_invert = pick_option(&c, &v, "i", NULL);
	int order = atoi(pick_option(&c, &v, "o", "-3"));
	if (c < 4 || c > 6)
		return fprintf(stderr, "usage:\n\t"
		"%s [-i] [-o {0|1|2|-3|3|5|7}] hom w h [in [out]]\n"
		//0                            1   2 3  4   5
		"\t-i\tinvert input homography\n"
		"\t-o\tchose interpolation order (default -3 = bicubic)\n"
		, *v);
	double H_direct[9], H_inv[9];
	read_n_doubles_from_string(H_direct, v[1], 9);
	invert_homography(H_inv, H_direct);
	double *H = do_invert ? H_inv : H_direct;
	int ow = atoi(v[2]);
	int oh = atoi(v[3]);
	char *filename_in  = c > 4 ? v[4] : "-";
	char *filename_out = c > 5 ? v[5] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_split(filename_in, &w, &h, &pd);
	float *y = xmalloc(ow * oh * pd * sizeof*y);

	int r = 0;
	for (int i = 0; i < pd; i++)
		r += shomwarp(y + i*ow*oh, ow, oh, H, x + i*w*h, w, h, order);

	iio_save_image_float_split(filename_out, y, ow, oh, pd);
	return r;
}
