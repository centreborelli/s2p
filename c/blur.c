// blur of a 2D image, using various kernels

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>



#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
#endif

#include "fail.c"
#include "xmalloc.c"

static void *fftwf_xmalloc(size_t n)
{
	float *r = fftwf_malloc(n);
	if (!r)
		fail("coult not fftwf_malloc %zu bytes\n", n);
	return r;
}

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORK(n) for(int k=0;k<(n);k++)
#define FORL(n) for(int l=0;l<(n);l++)



#ifndef USE_WISDOM
void evoke_wisdom(void) {}
void bequeath_wisdom(void) {}
#else//USE_WISDOM
#include "fftwisdom.c"
#endif//USE_WISDOM



// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued image
static void fft_2dfloat(fftwf_complex *fx, float *x, int w, int h)
{
	fftwf_complex *a = fftwf_xmalloc(w*h*sizeof*a);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_dft_2d(h, w, a, fx,
						FFTW_FORWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h) a[i] = x[i]; // complex assignment!
	fftwf_execute(p);

	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_cleanup();
}

#include "smapa.h"
SMART_PARAMETER_SILENT(FIWARN,1)

// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2dfloat(float *ifx,  fftwf_complex *fx, int w, int h)
{
	fftwf_complex *a = fftwf_xmalloc(w*h*sizeof*a);
	fftwf_complex *b = fftwf_xmalloc(w*h*sizeof*b);

	//fprintf(stderr, "planning...\n");
	evoke_wisdom();
	fftwf_plan p = fftwf_plan_dft_2d(h, w, a, b,
						FFTW_BACKWARD, FFTW_ESTIMATE);
	bequeath_wisdom();
	//fprintf(stderr, "...planned!\n");

	FORI(w*h) a[i] = fx[i];
	fftwf_execute(p);
	float scale = 1.0/(w*h);
	FORI(w*h) {
		fftwf_complex z = b[i] * scale;
		ifx[i] = crealf(z);
		if (FIWARN() > 0)
		{
			if (cimagf(z) > 0.001)
				fail("z is not real {cimagf(z)=%g} (set FIWARN=0 to run anyway)", cimagf(z));
			//assert(cimagf(z) < 0.001);
		}
	}
	fftwf_destroy_plan(p);
	fftwf_free(a);
	fftwf_free(b);
	fftwf_cleanup();
}

SMART_PARAMETER_SILENT(BLUR_INVERSE,0)
SMART_PARAMETER_SILENT(BLUR_INVERSE_WIENER,0)
#define UGLY_HACK_FOR_WIENER_FILTERING 1

static void pointwise_complex_multiplication(fftwf_complex *w,
		fftwf_complex *z, fftwf_complex *x, int n)
{
#ifdef UGLY_HACK_FOR_WIENER_FILTERING
	if (BLUR_INVERSE() > 0 || BLUR_INVERSE_WIENER() > 0) {
		if (BLUR_INVERSE_WIENER() > 0) {
			float t = BLUR_INVERSE_WIENER();
			FORI(n)
				w[i] = z[i] * x[i]/(cabs(x[i])*cabs(x[i])+t);
		} else
			FORI(n)
				w[i] = z[i] / x[i];
	}
	else
#endif//UGLY_HACK_FOR_WIENER_FILTERING
		FORI(n)
			w[i] = z[i] * x[i];
}

static void fill_2d_gaussian_image(float *gg, int w, int h, float inv_s)
{
	float (*g)[w] = (void *)gg;
	float alpha = inv_s*inv_s/(M_PI);

	FORJ(h) FORI(w) {
		float x = i < w/2 ? i : i - w;
		float y = j < h/2 ? j : j - h;
		float r = hypot(x, y);
		g[j][i] = alpha * exp(-r*r*inv_s*inv_s);
	}

	// if the kernel is too large, it escapes the domain, so the
	// normalization above must be corrected
	double m = 0;
	FORJ(h) FORI(w) m += g[j][i];
	FORJ(h) FORI(w) g[j][i] /= m;
}


// gaussian blur of a gray 2D image
void gblur_gray(float *y, float *x, int w, int h, float s)
{
	s = 1/s;

	fftwf_complex *fx = fftwf_xmalloc(w*h*sizeof*fx);
	fft_2dfloat(fx, x, w, h);

	float *g = xmalloc(w*h*sizeof*g);
	fill_2d_gaussian_image(g, w, h, s);

	fftwf_complex *fg = fftwf_xmalloc(w*h*sizeof*fg);
	fft_2dfloat(fg, g, w, h);

	pointwise_complex_multiplication(fx, fx, fg, w*h);
	ifft_2dfloat(y, fx, w, h);

	fftwf_free(fx);
	fftwf_free(fg);
	free(g);
}


// gausian blur of a 2D image with pd-dimensional pixels
// (the blurring is performed independently for each co-ordinate)
void gblur(float *y, float *x, int w, int h, int pd, float s)
{
	float *c = xmalloc(w*h*sizeof*c);
	float *gc = xmalloc(w*h*sizeof*gc);
	FORL(pd) {
		FORI(w*h) {
			float tmp = x[i*pd + l];
			if (!isfinite(tmp))
				tmp = 0;
			c[i] = tmp;//x[i*pd + l];
		}
		if (s)
			gblur_gray(gc, c, w, h, s);
		else FORI(w*h) gc[i] = c[i];
		FORI(w*h)
			y[i*pd + l] = gc[i];
	}
	free(c);
	free(gc);
}


static float kernel_2d_square(float x, float y, float *p)
{
	int nx, ny;
	if (p[0] == 1)
		nx = ny = p[1];
	else if (p[0] == 2) {
		nx = p[1];
		ny = p[2];
	} else
		fail("square kernel needs 1 or 2 parameters (got %g)",p[0]);

	float r = 0;
	if (2*fabs(x) < nx && 2*fabs(y) < ny)
	       r = 1;

	return r;
}

static float kernel_2d_disk(float x, float y, float *p)
{
	float radius = p[1];

	float r = 0;
	if (hypot(x, y) < radius)
	       r = 1;

	return r;
}

static float kernel_2d_gaussian(float x, float y, float *p)
{
	float sigma = p[1];

	float a = x*x + y*y;
	float r = exp(-a/(2*sigma*sigma));
	return r;
}

static float kernel_2d_cauchy(float x, float y, float *p)
{
	float sigma = p[1];

	float a = (x*x + y*y)/(sigma*sigma);
	float r = 1/(1+a);
	return r;
}

static float kernel_2d_powerlaw2(float x, float y, float *p)
{
	float sigma = p[1];
	float power = 2;

	float a = (x*x + y*y)/(sigma*sigma);
	float r = 1.0/(1.0 + a*a);
	return r;
}

static float kernel_2d_laplace(float x, float y, float *p)
{
	float sigma = p[1];

	float r = exp(-M_SQRT2*hypot(x,y)/sigma);
	return r;
}

static void fill_kernel_image(float *kk, int w, int h,
		float (*f)(float,float,float*), float *p)
{

	float (*k)[w] = (void *)kk;

	FORJ(h) FORI(w) {
		float x = i < w/2 ? i : i - w;
		float y = j < h/2 ? j : j - h;
		k[j][i] = f(x, y, p);
	}

	// if the kernel is too large, it escapes the domain, so the
	// normalization above must be corrected
	double m = 0;
	FORJ(h) FORI(w) m += k[j][i];
	FORJ(h) FORI(w) k[j][i] /= m;

}

static void substract_from_identity(float *k, int w, int h)
{
	for (int i = 0; i < w*h; i++)
		k[i] *= -1;
	k[0] += 1;
}

static void gray_fconvolution_2d(float *y, float *x, fftwf_complex *fk,
		int w, int h)
{
	fftwf_complex *fx = fftwf_xmalloc(w*h*sizeof*fx);
	fft_2dfloat(fx, x, w, h);

	pointwise_complex_multiplication(fx, fx, fk, w*h);
	ifft_2dfloat(y, fx, w, h);

	fftwf_free(fx);
}

static void color_fconvolution_2d(float *y, float *x, fftwf_complex *fk,
		int w, int h, int pd)
{
	float *c = xmalloc(w*h*sizeof*c);
	float *kc = xmalloc(w*h*sizeof*kc);
	FORL(pd) {
		FORI(w*h) {
			float tmp = x[i*pd + l];
			if (!isfinite(tmp))
				tmp = 0;
			c[i] = tmp;//x[i*pd + l];
		}
		gray_fconvolution_2d(kc, c, fk, w, h);
		FORI(w*h)
			y[i*pd + l] = kc[i];
	}
	free(c);
	free(kc);
}

void blur_2d(float *y, float *x, int w, int h, int pd,
		char *kernel_id, float *param, int nparams)
{
	float p[1+nparams];
	for (int i = 0; i < nparams; i++)
		p[1+i] = param[i];
	p[0] = nparams;
	//FORI(nparams+1) fprintf(stderr, "p[%d:%d] = %g\n", i,nparams,p[i]);

	if (nparams == 1 && param[0] == 0) {
		for (int i = 0; i < w*h*pd; i++)
			y[i] = x[i];
		return;
	}

	float (*f)(float,float,float*);
	switch(tolower(kernel_id[0])) {
	case 'g': f = kernel_2d_gaussian; break;
	case 'l': f = kernel_2d_laplace;  break;
	case 'c': f = kernel_2d_cauchy;   break;
	case 'd': f = kernel_2d_disk;     break;
	case 's': f = kernel_2d_square;   break;
	case 'p': f = kernel_2d_powerlaw2;   break;
	default: fail("unrecognized kernel name \"%s\"", kernel_id);
	}

	float *k = xmalloc(w*h*sizeof*k);
	fill_kernel_image(k, w, h, f, p);
	if (isupper(kernel_id[0]))
		substract_from_identity(k, w, h);
	//void iio_save_image_float(char*,float*,int,int);
	//iio_save_image_float("/tmp/blurk.tiff", k, w, h);

	fftwf_complex *fk = fftwf_xmalloc(w*h*sizeof*fk);
	fft_2dfloat(fk, k, w, h);
	free(k);

	color_fconvolution_2d(y, x, fk, w, h, pd);

	fftwf_free(fk);
}

#ifndef OMIT_BLUR_MAIN

static int print_version(void)
{
	printf("blur 1.0\n\n"
	"Written by Enric Meinhardt-Llopis\n");
	return 0;
}

static int print_help(void)
{
	printf(
	"Smooth a 2D image using the selected kernel\n"
	"\n"
	"Usage: blur KERNEL SIZE [input [output]]\n"
	"\n"
	"Options:\n"
	"\tKERNEL\tthe name of the kernel.\n\tPossible values: "
	"gaussian, laplace, cauchy, disk, square\n\n"
	"\tSIZE\tthe width of the kernel\n"
	"\n"
	"Examples:\n"
	"blur gaussian 1.6 in.png out.png\n"
	"\n"
	"Report bugs to <enric.meinhardt@cmla.ens-cachan.fr>.\n"
	);
	return 0;
}

#include "iio.h"
#include "parsenumbers.c"
int main(int c, char *v[])
{
	if (c == 2 && 0 == strcmp(v[1], "--version")) return print_version();
	if (c == 2 && 0 == strcmp(v[1], "--help")) return print_help();
	if (c != 5 && c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t"
				"%s kernel \"params\" [in [out]]\n", *v);
		//                0 1        2         3   4
		return EXIT_FAILURE;
	}
	char *kernel_id = v[1];
	char *kernel_params = v[2];
	char *filename_in = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "-";

	int maxparam = 10;
	float param[maxparam];
	int nparams = parse_floats(param, maxparam, kernel_params);
	if (nparams < 1) fail("please, give at least one parameter");
	//fprintf(stderr, "nparams = %d\n", nparams);
	//FORI(nparams)
	//	fprintf(stderr, "param[%d] = %g\n", i, param[i]);

	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *y = xmalloc(w*h*pd*sizeof*y);

	blur_2d(y, x, w, h, pd, kernel_id, param, nparams);

	iio_save_image_float_vec(filename_out, y, w, h, pd);
	free(x);
	free(y);
	return EXIT_SUCCESS;
}
#endif//OMIT_BLUR_MAIN
