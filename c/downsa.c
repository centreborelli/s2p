#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"


#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef ODDP
#define ODDP(x) ((x)&1)
#endif
#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

#include "fail.c"
#include "xmalloc.c"

struct statistics_float {
	float min, max, median, average, sample, variance, middle, laverage;
};

//SMART_PARAMETER_INT(STATISTIC_MEDIAN_BIAS,0)
//SMART_PARAMETER_SILENT(STATISTIC_MIDDLE_BIAS,-1)

#define STATISTIC_MEDIAN_BIAS 0
#define STATISTIC_MIDDLE_BIAS -1

int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

int compare_floats(const void *a, const void *b)
{
	const float *da = (const float *) a;
	const float *db = (const float *) b;
	return (*da > *db) - (*da < *db);
}

static void statistics_getf_spoilable(struct statistics_float *s, float *f,
		int n)
{
	s->middle = f[n/2-1];
	int mt = STATISTIC_MEDIAN_BIAS;
	int mi = STATISTIC_MIDDLE_BIAS;
	switch(mi)
	{
		case -1: break;
		case 0: s->middle += f[n/2]; s->middle /=2; break;
		case 1: s->middle = f[n/2]; break;
		default: fail("bad STATISTIC_MEDIAN_BIAS %d", mt);
	}
	//
	qsort(f, n, sizeof*f, compare_floats);
	s->min = f[0];
	s->max = f[n-1];
	s->median = f[n/2-1];
	if (EVENP(n))
	{
		int mtype = STATISTIC_MEDIAN_BIAS;
		switch(mtype)
		{
			case -1: break;
			case 0: s->median += f[n/2]; s->median /=2; break;
			case 1: s->median = f[n/2]; break;
			default: fail("bad STATISTIC_MEDIAN_BIAS %d", mtype);
		}
	}
	s->average = 0;
	for (int i = 0; i < n; i++)
		s->average += f[i];
	s->average /= n;
	s->laverage = 0;
	for (int i = 0; i < n; i++)
		s->laverage += exp(f[i]/255);
	s->laverage = log(s->laverage/n)*255;
	s->variance = 0;
	for (int i = 0; i < n; i++)
		s->variance = hypot(s->variance, s->average - f[i]);
	s->sample = f[randombounds(0, n-1)];
}

static void statistics_getf(struct statistics_float *s, float *fin, int n)
{
	if (n > 1000) {
		float *f = xmalloc(n*sizeof*f);
		memcpy(f, fin, n*sizeof*f);
		statistics_getf_spoilable(s, f, n);
		xfree(f);
	} else {
		float f[n];
		memcpy(f, fin, n*sizeof*f);
		statistics_getf_spoilable(s, f, n);
	}
}

static bool innerP(int w, int h, int i, int j)
{
	if (i < 0) return false;
	if (j < 0) return false;
	if (i >= w) return false;
	if (j >= h) return false;
	return true;
}

static void downsa2d(float *oy, float *ox, int w, int h, int pd, int n, int ty)
{
	int W = w/n;
	int H = h/n;
	float (*x)[w][pd] = (void*)ox;
	float (*y)[W][pd] = (void*)oy;
	for (int j = 0; j < H; j++)
	for (int i = 0; i < W; i++)
	for (int l = 0; l < pd; l++)
	{
		int nv = 0;
		float vv[n*n];
		for (int jj = 0; jj < n; jj++)
		for (int ii = 0; ii < n; ii++)
		if (innerP(w, h, i*n+ii, j*n+jj))
			vv[nv++] = x[j*n+jj][i*n+ii][l];
		struct statistics_float s;
		statistics_getf(&s, vv, nv);
		float g;
		switch (ty)
		{
		case 'i': g = s.min;      break;
		case 'e': g = s.median;   break;
		case 'a': g = s.max;      break;
		case 'v': g = s.average;  break;
		case 'V': g = s.laverage; break;
		case 'r': g = s.sample;   break;
		case 'f': g = vv[0];      break;
		default:  fail("downsa type %c not implemented", ty);
		}
		y[j][i][l] = g;
	}
}

int main_downsa(int c, char *v[])
{
	if (c != 5 && c != 4 && c != 3) {
		fprintf(stderr,"usage:\n\t%s {i|e|a|v|r|f} n [in [out]]\n", *v);
		//                         0        1      2  3   4
		return EXIT_FAILURE;
	}
	int n = atoi(v[2]);
	char *in = c > 3 ? v[3] : "-";
	char *out = c > 4 ? v[4] : "-";

	int w, h, pd;
	float *x = iio_read_image_float_vec(in, &w, &h, &pd);
	int W = w/n;
	int H = h/n;
	float *y = xmalloc(W*H*pd*sizeof*y);
	downsa2d(y, x, w, h, pd, n, v[1][0]);
	iio_save_image_float_vec(out, y, W, H, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	return main_downsa(c, v);
}
#endif//OMIT_MAIN
