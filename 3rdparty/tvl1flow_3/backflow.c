#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

#include "fragments.c"

static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static float getsample(float *fx, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	float (*x)[w][pd] = (void*)fx;
	return x[j][i][l];
	//return x[(i+j*w)*pd + l];
}

static float getsamplen(float *fx, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return NAN;
	float (*x)[w][pd] = (void*)fx;
	return x[j][i][l];
	//return x[(i+j*w)*pd + l];
}

static void bilinear_interpolation_at(float *result,
		float *x, int w, int h, int pd,
		float p, float q)
{
	int ip = p;
	int iq = q;
	FORL(pd) {
		float a = getsamplen(x, w, h, pd, ip  , iq  , l);
		float b = getsamplen(x, w, h, pd, ip+1, iq  , l);
		float c = getsamplen(x, w, h, pd, ip  , iq+1, l);
		float d = getsamplen(x, w, h, pd, ip+1, iq+1, l);
		float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
		result[l] = r;
	}
}


static void invflow(float *ou, float *flo, float *pin, int w, int h, int pd)
{
	float (*out)[w][pd] = (void*)ou;
	float (*in)[w][pd] = (void*)pin;
	float (*flow)[w][2] = (void*)flo;

	FORJ(h) FORI(w) {
		float p[2] = {i + flow[j][i][0], j + flow[j][i][1]};
		float result[pd];
		bilinear_interpolation_at(result, pin, w, h, pd, p[0], p[1]);
		FORL(pd)
			out[j][i][l] = result[l];
	}
}

int main_backflow(int c, char *v[])
{
	if (c != 2 && c != 4 && c != 3) {
		fprintf(stderr, "usage:\n\t%s flow [in [out]]\n", *v);
		//                          0 1     2   3
		return EXIT_FAILURE;
	}
	char *inname = c > 2 ? v[2] : "-";
	char *outname = c > 3 ? v[3] : "-";

	int w, h, pd;
	float *flow = iio_read_image_float_vec(v[1], &w, &h, &pd);
	fprintf(stderr, "w h pd P = %d %d %d %d\n", w, h, pd, w*h*pd);

	if (pd != 2)
		error("flow must have two-dimensional pixels\n");

	int iw, ih;
	float *in = iio_read_image_float_vec(inname, &iw, &ih, &pd);
	fprintf(stderr, "w h pd P = %d %d %d %d\n", iw, ih, pd, iw*ih*pd);
	if (iw != w || ih != h)
		error("flow and image size mismatch\n");
	float *out = xmalloc(w*h*pd*sizeof*out);
	fprintf(stderr, "p = %p\n", (void*)out);
	invflow(out, flow, in, w, h, pd);
	iio_save_image_float_vec(outname, out, w, h, pd);
	return EXIT_SUCCESS;
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	return main_backflow(c, v);
}
#endif//OMIT_MAIN
