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

#include "fail.c"
#include "xmalloc.c"
#include "getpixel.c"

#include "smapa.h"

SMART_PARAMETER_SILENT(NEAREST,0)
SMART_PARAMETER_SILENT(BILINEAR,0)

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

static void interpolate_nearest(float *result,
		float *x, int w, int h, int pd,
		float p, float q)
{
	int ip = round(p);
	int iq = round(q);
	FORL(pd) {
		float r = getsamplen(x, w, h, pd, ip  , iq  , l);
		result[l] = r;
	}
}


SMART_PARAMETER_SILENT(BACKDIV,0)
SMART_PARAMETER_SILENT(BACKDET,0)
SMART_PARAMETER_SILENT(BFBOUND,0)

static void compute_flow_div(float *d, float *u, int w, int h)
{
	getsample_operator p = getsample_1;
	FORJ(h) FORI(w) {
		float ux = 0.5*(p(u,w,h,2, i+1, j, 0) - p(u,w,h,2, i-1, j, 0));
		float vy = 0.5*(p(u,w,h,2, i, j+1, 1) - p(u,w,h,2, i, j-1, 1));
		d[j*w+i] = ux + vy;
	}
}

static void compute_flow_det(float *d, float *u, int w, int h)
{
	getsample_operator p = getsample_1;
	FORJ(h) FORI(w) {
		float ux = 0.5*(p(u,w,h,2, i+1, j, 0) - p(u,w,h,2, i-1, j, 0));
		float vy = 0.5*(p(u,w,h,2, i, j+1, 1) - p(u,w,h,2, i, j-1, 1));
		float uy = 0.5*(p(u,w,h,2, i, j+1, 0) - p(u,w,h,2, i, j-1, 0));
		float vx = 0.5*(p(u,w,h,2, i+1, j, 1) - p(u,w,h,2, i-1, j, 1));
		d[j*w+i] = (1+ux)*(1+vy) - uy*vx;;
	}
}



static void env_interpolate_at(float *out,
		float *x, int w, int h, int pd,
		float p, float q)
{
	if (BILINEAR())
		bilinear_interpolation_at(out, x, w, h, pd, p, q);
	else
		interpolate_nearest(out, x, w,h, pd, p, q);
}

static void invflow(float *ou, float *flo, float *pin, int w, int h, int pd, int win, int hin)
{
	float (*out)[w][pd] = (void*)ou;
	float (*in)[win][pd] = (void*)pin;
	float (*flow)[w][2] = (void*)flo;
	float *flowdiv = NULL;
	float *flowdet = NULL;

	if (BACKDIV() > 0) {
		flowdiv = xmalloc(w * h * sizeof*flowdiv);
		compute_flow_div(flowdiv, flo, w, h);
	}

	if (BACKDET() > 0) {
		flowdet = xmalloc(w * h * sizeof*flowdiv);
		compute_flow_det(flowdet, flo, w, h);
	}

	FORJ(h) FORI(w) {
		float p[2] = {i + flow[j][i][0], j + flow[j][i][1]};
		float result[pd];

		env_interpolate_at(result, pin, win, hin, pd, p[0], p[1]);

		float factor = 1;
		if (flowdiv)
			factor = exp(BACKDIV() * flowdiv[j*w+i]);
		FORL(pd)
			out[j][i][l] = factor * result[l];
			//out[j][i][l] = 100*log(factor * exp(result[l]/100));
		if (flowdet) {
			float bd = BACKDET();
			float det = flowdet[j*w+i];
			if (det > bd) det = bd;
			//if (det < 1/bd) det = 1/bd;
			FORL(pd)
				out[j][i][l] = det*result[l];
				//out[j][i][l] = bd*log(det*(exp((result[l])/bd)));
		}
	}

	if (flowdiv)
		free(flowdiv);
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
	//fprintf(stderr, "w h pd P = %d %d %d %d\n", w, h, pd, w*h*pd);

	if (pd != 2)
		fail("flow must have two-dimensional pixels\n");

	int iw, ih;
	float *in = iio_read_image_float_vec(inname, &iw, &ih, &pd);
	//fprintf(stderr, "w h pd P = %d %d %d %d\n", iw, ih, pd, iw*ih*pd);
	float *out = xmalloc(w*h*pd*sizeof*out);
	//fprintf(stderr, "p = %p\n", (void*)out);
	invflow(out, flow, in, w, h, pd, iw, ih);
	iio_save_image_float_vec(outname, out, w, h, pd);
	return EXIT_SUCCESS;
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	return main_backflow(c, v);
}
#endif//OMIT_MAIN
