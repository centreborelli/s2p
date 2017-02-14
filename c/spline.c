/* Spline interpolation.
    Copyright (C) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
    Copyright (C) 2010 Pascal Monasse <monasse@imagine.enpc.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// NOTE: minor changes by Enric Meinhardt-Llopis (2015)
//
// These changes only affect the types of the external interface, the
// algorithmic part is the same as written by Moisan and Monasse.
//
// The changes are the following:
// 	1. Make the code portable C (remove C++ headers and casts)
// 	2. Remove the dependence on an external image struct
// 	3. Change spaces to tabs
// 	4. Remove Doxygen-style markup in comments
// 	5. Remove coefficients of even-order splines
// 	6. Added an optional "main" function
//
// Please, look at file "spline_orig_mm.c" for the original file.


#include <assert.h>
#include <stdbool.h>
#include <math.h>

static double initcausal(float *c, int step, int n, double z)
{
	double zk,z2k,iz,sum;

	zk = z; iz = 1/z;
	z2k = pow(z,n-1);
	sum = c[0] + z2k * c[step*(n-1)];
	z2k = z2k*z2k*iz;
	for (int k = 1; k <= n-2; k++) {
		sum += (zk + z2k) * c[step*k];
		zk *= z;
		z2k *= iz;
	}
	return (sum/(1-zk*zk));
}

static double initanticausal(float *c, int step, int n, double z)
{
	return (z/(z*z-1)) * (z * c[step*(n-2)] + c[step*(n-1)]);
}

static void invspline1D(float *c, int step, int size, double *z, int npoles)
{
	/* normalization */
	double lambda=1;
	for (int k = npoles-1; k >= 0; k--)
		lambda *= (1-z[k])*(1-1/z[k]);
	for (int n = size-1; n >= 0; n--)
		c[step*n] *= lambda;

	for (int k = 0 ; k < npoles; k++) { // Loop on poles
		/* forward recursion */
		c[0] = initcausal(c, step, size, z[k]);
		for (int n = 1; n < size; n++)
			c[step*n] += z[k]*c[step*(n-1)];
		/* backward recursion */
		c[step*(size-1)] = initanticausal(c, step, size, z[k]);
		for (int n = size-2; n >= 0; n--)
			c[step*n] = z[k]*(c[step*(n+1)]-c[step*n]);
	}
}

// Put in array z the poles of the spline of given order.
static bool fill_poles(double* z, int order)
{
	switch(order) {
	case 1: break;
	case 3: z[0] = sqrt(3)-2;  /* sqrt(3)-2 */
		break;
	case 5: z[0] = -4.305753470999738e-1;
		z[1] = -4.309628820326465e-2;
		break;
	case 7: z[0] = -5.352804307964382e-1;
		z[1] = -1.225546151923267e-1;
		z[2] = -9.148694809608277e-3;
		break;
	case 9: z[0] = -6.079973891686259e-1;
		z[1] = -2.017505201931532e-1;
		z[2] = -4.322260854048175e-2;
		z[3] = -2.121306903180818e-3;
		break;
	case 11:z[0] = -6.612660689007345e-1;
		z[1] = -2.721803492947859e-1;
		z[2] = -8.975959979371331e-2;
		z[3] = -1.666962736623466e-2;
		z[4] = -5.105575344465021e-4;
		break;
	default: return false;
	}
	return true;
}

// Prepare image (in-place) for cardinal spline interpolation.
bool prepare_spline(float *img, int w, int h, int pd, int order)
{
	if(order < 3)
		return true;

	// Replace nans and infinities with 0
	for (int i = 0; i < w * h * pd; i++)
		if (!isfinite(img[i]))
			img[i] = 0;

	// Init poles of associated z-filter
	double z[5];
	if (! fill_poles(z, order))
		return false;
	int npoles = order / 2;

	for (int k = 0; k < pd; k++) { // Loop on image components
		for (int y = 0; y < h; y++) // Filter on lines
			invspline1D(img + (y*w+0)*pd + k, pd * 1, w, z, npoles);
		for (int x = 0; x < w; x++) // Filter on columns
			invspline1D(img + (0*w+x)*pd + k, pd * w, h, z, npoles);
	}

	return true;
}

/* c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... */

/* coefficients for cubic interpolant (Keys' function) */
static void keys(float *c, float t, float a)
{
	float t2 = t * t;
	float at = a * t;
	c[0] = a * t2 * (1 - t);
	c[1] = ( 2*a + 3 - (a + 2)*t ) * t2 - at;
	c[2] = ((a + 2)*t - a - 3) * t2 + 1;
	c[3] = a * (t - 2) * t2 + at;
}

/* coefficients for cubic spline */
static void spline3(float *c, float t)
{
	float T = 1 - t;
	c[0] = t * t * t / 6;
	c[1] = ( 4 - 3 * (1 + t) * T * T ) / 6;
	c[2] = ( 4 - 3 * (2 - t) * t * t ) / 6;
	c[3] = T * T * T / 6;
}

/* pre-computation for spline of order >3 */
static void init_splinen(float *a, int n)
{
	a[0] = 1;
	for ( int k = 2; k <= n; k++)
	          a[0] /= k;
	for ( int k = 1; k <= n+1; k++)
		a[k] = - a[k-1] *(n+2-k)/k;
}

/* fast integral power function */
static float ipow(float x, int n)
{
	float res = 1;
	for(; n; n>>=1) {
		if(n&1)
			res *= x;
		x *= x;
	}
	return res;
}

/* coefficients for spline of order >3 */
static void splinen(float* c, float t, float* a, int n)
{
	for (int i = 0; i < n+1; i++)
		c[i] = 0;
	for(int k=0; k <= n+1; k++) {
		float xn = ipow(t+k, n);
		for(int i=k; i <= n; i++)
			c[i] += a[i-k]*xn;
	}
}

static bool insideP(int w, int h, int i, int j)
{
	return (i >= 0) && (j >= 0) && (i < w) && (j < h);
}

typedef float (*getsample_t)(float*,int,int,int,int,int,int);

static float getsample_ass(float *x, int w, int h, int pd, int i, int j, int l)
{
	assert(i >= 0);
	assert(j >= 0);
	assert(l >= 0);
	assert(i < w);
	assert(j < h);
	assert(l < pd);
	return x[(i+j*w)*pd + l];
}

// like n%p, but works for all numbers
static int good_modulus(int n, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(n, -p);

	int r = n % p;
	r = r < 0 ? r + p : r;

	assert(r >= 0);
	assert(r < p);
	return r;
}

// symmetrized and periodized index
static int positive_reflex(int n, int p)
{
	int r = good_modulus(n, 2*p);
	if (r == p)
		r -= 1;
	if (r > p)
		r = 2*p - r;
	assert(r >= 0);
	assert(r < p);
	return r;
}

// extrapolate by periodized reflection
static float getsample_2(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = positive_reflex(i, w);
	j = positive_reflex(j, h);
	return getsample_ass(x, w, h, pd, i, j, l);
}

// Spline interpolation of given order of image im at point (x,y).
// out must be an array of size the number of components.
// Supported orders: 0(nn), 1(bilinear), -3(Keys's bicubic), 3, 5, 7, 9, 11.
// Success means a valid order and pixel in image.
bool evaluate_spline_at(float *out,
		float *img, int w, int h, int pd,
		int order, float x, float y)
{
	float  cx[12],cy[12];

	/* CHECK ORDER */
	if (order != 0 && order != 1 && order != -3 && order != 3 &&
			order != 5 && order != 7 && order != 9 && order != 11)
		return false;

	float ak[13];
	if (order > 3)
		init_splinen(ak, order);

	/* INTERPOLATION */
	if(order == 0) { /* zero order interpolation (pixel replication) */
		int xi = x;
		int yi = y;
		if (!insideP(w, h, xi, yi))
			return false;
		float* p = img + (w*yi+xi) * pd;
		for(int i = 0; i < pd; i++)
			out[i] = p[i];
	} else { /* higher order interpolations */
		if (!(x>=0 && x<=w && y>=0 && y<=h))
			return false;
		x -= 0.5; y -= 0.5;
		int xi = (x<0)? -1: x;
		int yi = (y<0)? -1: y;
		float ux = x - xi;
		float uy = y - yi;
		float paramKeys = -0.5;
		switch(order)  {
		case 1: /* first order interpolation (bilinear) */
			cx[0] = ux; cx[1] = 1-ux;
			cy[0] = uy; cy[1] = 1-uy;
			break;
		case -3: /* third order interpolation (bicubic Keys) */
			keys(cx, ux, paramKeys);
			keys(cy, uy, paramKeys);
			break;
		case 3: /* spline of order 3 */
			spline3(cx, ux);
			spline3(cy, uy);
			break;
		default: /* spline of order >3 */
			splinen(cx, ux, ak, order);
			splinen(cy, uy, ak, order);
			break;
		}
		int n2 = (order == -3) ? 2 : (order+1)/2;
		int n1 = 1 - n2;
		/* this test saves computation time */
		if (insideP(w, h, xi+n1, yi+n1) && insideP(w, h, xi+n1, yi+n2))
		{
			for (int k = 0; k < pd; k++) {
				out[k] = 0;
				for (int dy = n1; dy <= n2; dy++) {
					int ii = xi + n1;
					int jj = yi + dy;
					float* v = img + (jj*w+ii)*pd+k;
					for (int dx = n1; dx <= n2; dx++) {
						float f = *v;
						out[k]+=cy[n2-dy]*cx[n2-dx] * f;
						v += pd;
					}
				}
			}
		} else
			for (int k = 0; k < pd; k++)
			{
				out[k] = 0;
				for (int dy = n1; dy <= n2; dy++)
				for (int dx = n1; dx <= n2; dx++) {
					int ii = xi + dx;
					int jj = yi + dy;
					getsample_t S = getsample_2;
					float v = S(img, w, h, pd, ii, jj, k);
					out[k] += cy[n2-dy] * cx[n2-dx] * v;
				}
			}
	}
	return true;
}

// this main serves as a unit test of the code above
#ifdef MAIN_SPLINE
#include <stdio.h>
#include <stdlib.h>
#include "iio.h"

static double invert_affinity(double invA[6], double A[6])
{
	double a, b, c, d, p, q;
	a=A[0]; b=A[1]; p=A[2];
	c=A[3]; d=A[4]; q=A[5];
	double det = a*d - b*c;
	invA[0] =  d; invA[1] = -b;
	invA[3] = -c; invA[4] =  a;
	invA[2] = b*q - d*p;
	invA[5] = c*p - a*q;
	for (int i = 0; i < 6; i++)
		invA[i] /= det;
	return det;
}

static void apply_affinity(double y[2], double A[6], double x[2])
{
	double p = A[0]*x[0] + A[1]*x[1] + A[2];
	double q = A[3]*x[0] + A[4]*x[1] + A[5];
	y[0] = p;
	y[1] = q;
}

static void naive_affine_map_using_spline(float *y, int out_w, int out_h,
		float *x, int w, int h, int pd, double A[6], int order)
{
	float *fx = malloc(w * h * pd * sizeof*fx);
	for (int i = 0; i < w * h * pd; i++)
		fx[i] = x[i];
	bool r = prepare_spline(fx, w, h, pd, order);
	if (!r) exit(fprintf(stderr,"prepare spline failed (ord=%d)\n",order));
	//if (1) {
	//	char buf[FILENAME_MAX];
	//	snprintf(buf,FILENAME_MAX,"/tmp/prepared_spline_%d.tiff",order);
	//	iio_save_image_float_vec(buf, fx, w, h, pd);
	//}
	double invA[6]; invert_affinity(invA, A);
	for (int j = 0; j < out_h; j++)
	for (int i = 0; i < out_w; i++)
	{
		double p[2] = {i, j};
		apply_affinity(p, A, p);
		float *out = y + (j * out_w + i) * pd;
		evaluate_spline_at(out, fx, w, h, pd, order, p[0], p[1]);
	}
	free(fx);
}

int main(int c, char *v[])
{
	// process input arguments
	if (c < 10 || c > 13 ) {
		fprintf(stderr, "usage:\n\t"
			"%s a b p c d q W H ord [in.png [out.png]]\n", *v);
		//       0  1 2 3 4 5 6 7 8  9   10      11
		return 1;
	}
	double A[6]; for (int i = 0; i < 6; i++) A[i] = atof(v[1+i]);
	int out_w = atoi(v[7]);
	int out_h = atoi(v[8]);
	int order = atoi(v[9]);
	char *filename_in  = c > 10 ? v[10] : "-";
	char *filename_out = c > 11 ? v[11] : "-";

	// read input image
	int w, h, pd;
	float *x = iio_read_image_float_vec(filename_in, &w, &h, &pd);

	// allocate space for output image
	float *y = malloc(out_w * out_h * pd * sizeof*y);
	if (!y) return 2;

	// perform the computation
	naive_affine_map_using_spline(y, out_w, out_h, x, w, h, pd, A, order);

	// save output result
	iio_save_image_float_vec(filename_out, y, out_w, out_h, pd);

	// cleanup and exit
	free(x);
	free(y);
	return 0;
}
#endif//MAIN_SPLINE
