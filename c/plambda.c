// PLAMBDA(1)                  imscript                  PLAMBDA(1)        {{{1
//
// NAME                                                                    {{{2
//	plambda - a RPN calculator for image pixels
//
// SYNOPSIS                                                                {{{2
//	plambda a.pnb b.png c.png ... "lambda expression" > output
//	plambda a.pnb b.png c.png ... "lambda expression" -o output.png
//	plambda -c num1 num2 ... "lambda expression"
//
// DESCRIPTION                                                             {{{2
//	Plambda applies an expression to all the pixels of a collection of
//	images, and produces a single output image.  Each input image
//	corresponds to one of the variables of the expression (in alphabetical
//	order).  There are modifiers to the variables that allow access to the
//	values of neighboring pixels, or to particular components of a pixel.
//
// LANGUAGE                                                                {{{2
//	A "plambda" program is a sequence of tokens.  Tokens may be constants,
//	variables, or operators.  Constants and variables get their value
//	computed and pushed to the stack.  Operators pop values from the stack,
//	apply a function to them, and push back the results.
//
//	CONSTANTS: numeric constants written in scientific notation, and "pi"
//	OPERATORS: +, -, *, ^, /, and all the functions from math.h
//	VARIABLES: anything not recognized as a constant or operator.  There
//	must be as many variables as input images, and they are assigned to
//	images in alphabetical order.
//
//	All operators (unary, binary and ternary) are vectorializable.
//
//	Some "sugar" is added to the language:
//
//	Predefined variables (always preceeded by a colon):
//
//		TOKEN	MEANING
//
//		:i	horizontal coordinate of the pixel
//		:j	vertical coordinate of the pixel
//		:w	width of the image
//		:h	heigth of the image
//		:n	number of pixels in the image
//		:x	relative horizontal coordinate of the pixel
//		:y	relative horizontal coordinate of the pixel
//		:r	relative distance to the center of the image
//		:t	relative angle from the center of the image
//		:I	horizontal coordinate of the pixel (centered)
//		:J	vertical coordinate of the pixel (centered)
//		:W	width of the image divided by 2*pi
//		:H	height of the image divided by 2*pi
//
//
//	Variable modifiers acting on regular variables:
//
//		TOKEN	MEANING
//
//		x	value of pixel (i,j)
//		x(0,0)	value of pixel (i,j)
//		x(1,0)	value of pixel (i+1,j)
//		x(0,-1)	value of pixel (i,j-1)
//		...
//
//		x	value of pixel (i,j)
//		x[0]	value of first component of pixel (i,j)
//		x[1]	value of second component of pixel (i,j)
//
//		x(1,-1)[2] value of third component of pixel (i+1,j-1)
//
//
//	Stack operators (allow direct manipulation of the stack):
//
//		TOKEN	MEANING
//
//		del	remove the value at the top of the stack (ATTOTS)
//		dup	duplicate the value ATTOTS
//		rot	swap the two values ATTOTS
//		split	split the vector ATTOTS into scalar components
//		join	join the components of two vectors ATTOTS
//		join3	join the components of three vectors ATTOTS
//		n njoin	join the components of n vectors ATTOTS
//
//
//	Magic modifiers (which are not computed locally for each image):
//
//		x%i	value of the smallest sample of image "x"
//		x%a	value of the largest sample
//		x%v	average sample value
//		x%m	median sample value
//		x%I	value of the smallest pixel
//		x%A	value of the largest pixel
//		x%V	average pixel value
//		x%M	median pixel value (not implemented)
//		x%qn	nth sample percentile
//		x%Qn	nth pixel percentile (not implemented)
//		x%r	random sample of the image (not implemented)
//		x%R	random pixel of the image (not implemented)
//
//	Notice that the scalar "magic" modifiers may act upon individual
//	components, e.g. x[2]%i is the minimum value of the blue component.
//
//
//	Vectorial operators (non-vectorializable):
//
//		TOKEN	MEANING
//		randu	push a random number U(0,1)
//		randn	push a random number N(0,1)
//		randg	push a random number N(0,1)
//		randc	push a random number C(0,1)
//		randl	push a random number L(0,1)
//		rande	push a random number E(1)
//		randp	push a random number P(1)
//		rand	push a random integer returned from "rand()"
//		topolar	convert a 2-vector from cartesian to polar
//		frompolar
//		cprod	multiply two 2-vectors as complex numbers
//		mprod	multiply two vectors as matrices (4vector=2x2 matrix...)
//		vprod	vector product of two 3-vectors
//		sprod	scalar product of two n-vectors
//		mdet	determinant of a n-matrix (a n*n-vector)
//		mtrans	transpose of a matrix
//		mtrace	trace of a matrix
//		minv	inverse of a matrix
//		vavg	average value of a vector
//		vsum	sum of the components of a vector
//		vmul	product of the components of a vector
//		vmax	max component of a vector
//		vmin	min component of a vector
//		vnorm	euclidean norm of a vector
//		vdim	length of a vector
//
//
// OPTIONS                                                                 {{{2
//	-o file		save output to named file
//	-c		act as a symbolic calculator
//	-h		print short help message
//	--help		print longer help message
//	--man		print manpage (requires help2man)
//	--version	print program version
//
// EXAMPLES                                                                {{{2
// 	Sum two images:
//
// 		plambda a.png b.png "a b +" > aplusb.png
//
//	Add a gaussian to half of lena:
//
//		plambda /tmp/lena.png "x 2 / :r :r * -1 * 40 * exp 200 * +"
//
//	Forward differences to compute the derivative in horizontal direction:
//
//		plambda lena.png "x(1,0) x -"
//
//	Sobel edge detector:
//		plambda lena.png "x(1,0) 2 * x(1,1) x(1,-1) + + x(-1,0) 2 * x(-1,1) x(-1,-1) + + - x(0,1) 2 * x(1,1) x(-1,1) + + x(0,-1) 2 * x(1,-1) x(-1,-1) + + - hypot"
//
//	Color to gray:
//		plambda lena.png "x[0] x[1] x[2] + + 3 /"
//
//	Pick the blue channel of a RGB image:
//		plambda lena.png "x[2]"
//
//	Swap the blue an green channels of a RGB image (6 equivalent ways):
//		plambda lena.png "x[0] x[2] x[1] join3"
//		plambda lena.png "x[0] x[2] x[1] join join"
//		plambda lena.png "x[0] x[1] x[2] rot join3"
//		plambda lena.png "x[0] x[1] x[2] rot join join"
//		plambda lena.png "x split rot join join"
//		plambda lena.png "x split rot join3"
//
//	Merge the two components of a vector field into a single file
//		plambda x.tiff y.tiff "x y join" > xy.tiff
//
//	Set to 0 the green component of a RGB image
//		plambda lena.png "x[0] 0 x[2] join3"
//
//	Naive Canny filter:
//		cat lena.png | gblur 2 | plambda - "x(1,0) 2 * x(1,1) x(1,-1) + + x(-1,0) 2 * x(-1,1) x(-1,-1) + + - >1 x(0,1) 2 * x(1,1) x(-1,1) + + x(0,-1) 2 * x(1,-1) x(-1,-1) + + - >2 <1 <2 hypot <2 <1 atan2 join" | plambda - "x[0] 4 > >1 x[1] fabs pi 4 / > x[1] fabs pi 4 / 3 * < * >2 x[1] fabs pi 4 / < x[1] fabs pi 4 / 3 * > + >3 x[0] x[0](0,1) > x[0] x[0](0,-1) > * >4 x[0] x[0](1,0) > x[0] x[0](-1,0) > * >5 <1 <3 <5 * * <1 <2 <4 * * + x[0] *" | qauto | display
//
//	Anti-Lalpacian (solve Poisson equation):
//		cat lena.png | fft 1 | plambda - "x  :I :I * :J :J * + / -1 *" | fft -1 | qauto | display
//
//	Wiener Filter (for real kernels):
//		PREC=0.01
//		plambda kernel.fft image.fft "h[0] dup dup * $PREC + / y *"
//
//	Deconvolution using max frequency cut (for real kernels):
//		FCUT=80
//		plambda kernel.fft image.fft ":I :J hypot $FCUT < y h[0] / 0 if"
//
//	Generate a U(-1,1) scalar field with gaussian grain
//		GRAINSIZE=7
//		plambda zero:WxH "randn"|blur g $GRAINSIZE|plambda - "x $GRAINSIZE * pi sqrt * 2 * 2 sqrt / erf"
//
//	Generate a N(0,1) scalar field with gaussian grain
//		plambda zero:WxH "randn"|blur g $GRAINSIZE|plambda - "x $GRAINSIZE * pi sqrt * 2 *"
//
//	Generate a L(0,sigma=1) scalar field with gaussian grain
//		plambda zero:WxH "randn randn randn randn  4 njoin $GRAINSIZE * pi sqrt * 2 *"|blur g $GRAINSIZE|plambda - "x[0] x[1] * x[2] x[3] * - 2 sqrt /"
//
//	Periodic component of an image
//		  cat image|fftsym|fft|plambda - "x :I :I * :J :J * + *"|ifft|crop 0 0 `imprintf "%w %h" image`|fft|plambda - "x :I :I * :J :J * + / 4 /"|ifft >pcomponent
//
//
//
// AUTHOR                                                                  {{{2
//
//	Written by Enric Meinhardt-Llopis
//
//
// REPORTING BUGS                                                          {{{2
//
//	Download last version from https://github.com/mnhrdt
//	Report bugs to <enric.meinhardt@cmla.ens-cachan.fr>
//
// TODO LIST                                                               {{{2
//
//	* implement shunting-yard algorithm to admit infix notation
//	* handle 3D and nD images
//	* merge colonvars and magicvars (the only difficulty lies in naming)


// #includes {{{1

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//#define __STDC_IEC_559_COMPLEX__ 1
#ifdef __STDC_IEC_559_COMPLEX__
#include <complex.h>
#endif

#include "smapa.h"

#include "fail.c"
#include "xmalloc.c"
#include "random.c"
#include "parsenumbers.c"
#include "colorcoordsf.c"


// #defines {{{1

#define PLAMBDA_MAX_TOKENS 2049
#define PLAMBDA_MAX_VARLEN 0x100
#define PLAMBDA_MAX_PIXELDIM 0x100
#define PLAMBDA_MAX_MAGIC 42


#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif
#ifndef FORJ
#define FORJ(n) for(int j=0;j<(n);j++)
#endif
#ifndef FORK
#define FORK(n) for(int k=0;k<(n);k++)
#endif
#ifndef FORL
#define FORL(n) for(int l=0;l<(n);l++)
#endif


#define PLAMBDA_CONSTANT 0 // numeric constant
#define PLAMBDA_SCALAR 1   // pixel component
#define PLAMBDA_VECTOR 2   // whole pixel
#define PLAMBDA_OPERATOR 3 // function
#define PLAMBDA_COLONVAR 4 // colon-type variable
#define PLAMBDA_STACKOP 5  // stack operator
#define PLAMBDA_VARDEF 6   // register variable definition (hacky)
#define PLAMBDA_MAGIC 7    // "magic" modifier (requiring cached global data)
#define PLAMBDA_IMAGEOP 8    // comma-modified variable

#define IMAGEOP_IDENTITY 0
#define IMAGEOP_X 1
#define IMAGEOP_Y 2
#define IMAGEOP_XX 3
#define IMAGEOP_XY 4
#define IMAGEOP_YX 5
#define IMAGEOP_YY 6
#define IMAGEOP_NGRAD 8
#define IMAGEOP_LAP 9
#define IMAGEOP_CURV 10
#define IMAGEOP_ILAP 12
#define IMAGEOP_HESS 1001
#define IMAGEOP_GRAD 1002
#define IMAGEOP_DIV 1003
#define IMAGEOP_SHADOW 1004

#define SCHEME_FORWARD 0
#define SCHEME_BACKWARD 1
#define SCHEME_CENTERED 2
#define SCHEME_CONSISTENT 3
#define SCHEME_PRATT 4
#define SCHEME_SOBEL 5
#define SCHEME_PREWITT 6
#define SCHEME_SCHARR 7
#define SCHEME_MORPHO5_FORWARD 8
#define SCHEME_MORPHO5_BACKWARD 9
#define SCHEME_MORPHO5_CENTERED 10
#define SCHEME_MORPHO9_FORWARD 11
#define SCHEME_MORPHO9_BACKWARD 12
#define SCHEME_MORPHO9_CENTERED 13


// local functions {{{1

static double sum_two_doubles      (double a, double b) { return a + b; }
static double substract_two_doubles(double a, double b) { return a - b; }
static double multiply_two_doubles (double a, double b) { return a * b; }
static double divide_two_doubles   (double a, double b) {
	if (!b && !a) return 0;
	return a / b;
}
static double logic_g      (double a, double b) { return a > b; }
static double logic_l      (double a, double b) { return a < b; }
static double logic_e      (double a, double b) { return a == b; }
static double logic_ge     (double a, double b) { return a >= b; }
static double logic_le     (double a, double b) { return a <= b; }
static double logic_ne     (double a, double b) { return a != b; }
static double logic_if (double a, double b, double c) { return a ? b : c; }
static double logic_or (double a, double b) { return a || b; }
static double logic_and (double a, double b) { return a && b; }
static double logic_not (double a) { return !a; }

static double function_isfinite  (double x) { return isfinite(x); }
static double function_isinf     (double x) { return isinf(x);    }
static double function_isnan     (double x) { return isnan(x);    }
static double function_isnormal  (double x) { return isnormal(x); }
static double function_signbit   (double x) { return signbit(x);  }

static double inftozero (double x) { return isinf(x) ? 0 : x; }
static double nantozero (double x) { return isnan(x) ? 0 : x; }
static double notfintozero (double x) { return isfinite(x) ? x : 0; }
static double force_finite (double x) { return isfinite(x) ? x : 0; }
static double force_normal (double x) { return isnormal(x) ? x : 0; }

static double quantize_255 (double x)
{
	int ix = x;
	if (ix < 0) return 0;
	if (ix > 255) return 255;
	return ix;
}

static double quantize_easy(double x, double a, double b)
{
	return quantize_255(255.0*(x-a)/(b-a));
}

static double range(double x, double a, double b)
{
	return (x-a)/(b-a);
}

static double bound_double(double x, double a, double b)
{
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static void from_cartesian_to_polar(float *y, float *x)
{
	y[0] = hypot(x[0], x[1]);
	y[1] = atan2(x[1], x[0]);
}

static void from_polar_to_cartesian(float *y, float *x)
{
	y[0] = x[0] * cos(x[1]);
	y[1] = x[0] * sin(x[1]);
}

static void complex_product(float *xy, float *x, float *y)
{
	xy[0] = x[0]*y[0] - x[1]*y[1];
	xy[1] = x[0]*y[1] + x[1]*y[0];
}

static void complex_exp(float *y, float *x)
{
#ifdef __STDC_IEC_559_COMPLEX__
	*(complex float *)y = cexp(*(complex float *)x);
#else
	y[0] = exp(x[0]) * cos(x[1]);
	y[1] = exp(x[0]) * sin(x[1]);
#endif
}

#ifdef __STDC_IEC_559_COMPLEX__
#define REGISTERC(f) static void complex_ ## f(float *y, float *x) {\
	*(complex float *)y = f(*(complex float *)x); }
REGISTERC(cacos)
REGISTERC(cacosh)
REGISTERC(casin)
REGISTERC(casinh)
REGISTERC(catan)
REGISTERC(catanh)
REGISTERC(ccos)
REGISTERC(ccosh)
REGISTERC(cexp)
REGISTERC(clog)
REGISTERC(conj)
REGISTERC(cproj)
REGISTERC(csin)
REGISTERC(csinh)
REGISTERC(csqrt)
REGISTERC(ctan)
REGISTERC(ctanh)
#endif


static void matrix_product_clean(
		float *ab, int *ab_nrows, int *ab_ncols,
		float *a, int a_nrows, int a_ncols,
		float *b, int b_nrows, int b_ncols)
{
	if (a_ncols != b_nrows)
		fail("impossible matrix product (%d %d) * (%d %d)",
				a_nrows, a_ncols, b_nrows, b_ncols);
	*ab_nrows = a_nrows;
	*ab_ncols = b_ncols;
	float (*A)[a_ncols] = (void*)a;
	float (*B)[b_ncols] = (void*)b;
	float (*AB)[*ab_ncols] = (void*)ab;
	for (int i = 0; i < *ab_nrows; i++)
	for (int j = 0; j < *ab_ncols; j++)
	{
		AB[i][j] = 0;
		for (int k = 0; k < a_ncols; k++)
			AB[i][j] += A[i][k] * B[k][j];
	}
}

// instance of "bivector_function"
static int matrix_product(float *ab, float *a, float *b, int na, int nb)
{
	int a_nrows, a_ncols, b_nrows, b_ncols;
	if (na == 4 && nb == 4) {
		a_nrows = a_ncols = b_nrows = b_ncols = 2;
	} else if (na == 9 && nb == 9) {
		a_nrows = a_ncols = b_nrows = b_ncols = 3;
	} else if (na == 16 && nb == 16) {
		a_nrows = a_ncols = b_nrows = b_ncols = 4;
	} else if (na == 9 && nb == 3) {
		a_nrows = a_ncols = b_nrows = 3;
		b_ncols = 1;
	} else if (na == 4 && nb == 2) {
		a_nrows = a_ncols = b_nrows = 2;
		b_ncols = 1;
	} else if (na == 1 && nb == 1) {
		a_nrows = a_ncols = b_nrows = b_ncols = 1;
	} else if (na == 6 && nb == 2) {
		ab[0] = a[0]*b[0] + a[1]*b[1] + a[2];
		ab[1] = a[3]*b[0] + a[4]*b[1] + a[5];
		return 2;
	} else fail("bad matrix product (%d %d)", na, nb);
	assert(a_ncols == b_nrows);

	int ab_nrows, ab_ncols;
	matrix_product_clean(ab, &ab_nrows, &ab_ncols,
			a, a_nrows, a_ncols, b, b_nrows, b_ncols);
	return ab_nrows * ab_ncols;
}

// instance of "bivector_function"
static int vector_product(float *ab, float *a, float *b, int na, int nb)
{
	if (na != 3 || nb != 3)
		fail("bad vector product (%d %d)", na, nb);
	ab[0] = a[1]*b[2] - a[2]*b[1];
	ab[1] = a[2]*b[0] - a[0]*b[2];
	ab[2] = a[0]*b[1] - a[1]*b[0];
	return 3;
}

// instance of "bivector_function"
static int scalar_product(float *ab, float *a, float *b, int na, int nb)
{
	if (na != nb)
		fail("bad scalar product (%d %d)", na, nb);
	*ab = 0;
	for (int i = 0 ; i < na; i++)
		*ab += a[i] * b[i];
	return 1;
}

// instance of "univector_function"
static int matrix_determinant(float *r, float *a, int n)
{
	switch(n) {
	case 1: *r = *a; break;
	case 4: *r = a[0]*a[3] - a[1]*a[2]; break;
	case 6: *r = a[0]*a[4] - a[1]*a[3]; break;
	case 9: *r = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7]; break;
	default: fail("can not compute determinant of object of size %d", n);
	}
	return 1;
}

// instance of "univector_function"
static int matrix_inverse(float *r, float *a, int n)
{
	float det;
	matrix_determinant(&det, a, n);
	switch(n) {
	case 1:
		r[0] = 1/a[0];
		return 1;
	case 4:
		r[0] = a[3]/det;
		r[1] = -a[1]/det;
		r[2] = -a[2]/det;
		r[3] = a[0]/det;
		return 4;
	case 6:
		r[0] = a[4]/det;
		r[1] = -a[1]/det;
		r[2] = (a[1]*a[5]-a[2]*a[4])/det;
		r[3] = -a[3]/det;
		r[4] = a[0]/det;
		r[5] = (a[2]*a[3]-a[0]*a[5])/det;
		return 6;
	case 9:
		r[0] = (a[4]*a[8]-a[5]*a[7])/det;
		r[1] = (a[2]*a[7]-a[1]*a[8])/det;
		r[2] = (a[1]*a[5]-a[2]*a[4])/det;
		r[3] = (a[5]*a[6]-a[3]*a[8])/det;
		r[4] = (a[0]*a[8]-a[2]*a[6])/det;
		r[5] = (a[2]*a[3]-a[0]*a[5])/det;
		r[6] = (a[3]*a[7]-a[4]*a[6])/det;
		r[7] = (a[1]*a[6]-a[0]*a[7])/det;
		r[8] = (a[0]*a[4]-a[1]*a[3])/det;
		return 9;
	default: fail("can not compute inversion of object of size %d", n);
	}
}

// instance of "univector_function"
static int matrix_trace(float *r, float *a, int nn)
{
	int n;
	switch(nn) {
	case 1: n = 1; break;
	case 4: n = 2; break;
	case 9: n = 3; break;
	default: fail("can not compute trace of object of size %d", nn);
	}
	assert(n*n == nn);
	*r = 0;
	for (int i = 0; i < n; i++)
		*r += a[i*n+i];
	return 1;
}

// instance of "univector_function"
static int matrix_3x3to12x9(float *r, float *a, int n)
{
	if (n != 9)
		fail("mto12x9 needs a 3x3 matrix");
	if (108 > PLAMBDA_MAX_PIXELDIM)
		fail("MAX_PIXELDIM needs to be at least 108");
	for (int i = 0; i < 108; i++)
		r[i] = 0;
	r[0] = r[39] = r[78] = a[0];
	r[1] = r[40] = r[79] = a[1];
	r[2] = r[41] = r[80] = a[2];
	r[12] = r[51] = r[90] = a[3];
	r[13] = r[52] = r[91] = a[4];
	r[14] = r[53] = r[92] = a[5];
	r[24] = r[63] = r[102] = a[6];
	r[25] = r[64] = r[103] = a[7];
	r[26] = r[65] = r[104] = a[8];
	r[21] = r[34] = r[71] = 1;
	r[45] = r[82] = r[95] = -1;
	return 108;
}

// instance of "univector_function"
static int matrix_transpose(float *r, float *a, int nn)
{
	int n;
	switch(nn) {
	case 1: n = 1; break;
	case 4: n = 2; break;
	case 9: n = 3; break;
	default: fail("can not transpose object of size %d", nn);
	}
	assert(n*n == nn);
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		r[i*n+j] = a[j*n+i];
	return nn;
}

// instances of "univector_function"
static int xyz2rgb(float *y, float *x, int n)
	{ assert(n == 3); xyz_to_rgb_floats(y, x); return 3; }
static int rgb2xyz(float *y, float *x, int n)
	{ assert(n == 3); rgb_to_xyz_floats(y, x); return 3; }
static int hsv2rgb(float *y, float *x, int n)
	{ assert(n == 3); hsv_to_rgb_floats(y, x); return 3; }
static int rgb2hsv(float *y, float *x, int n)
	{ assert(n == 3); rgb_to_hsv_floats(y, x); return 3; }

// instance of "univector_function"
static int vector_avg(float *r, float *a, int n)
{
	*r = 0;
	for (int i = 0; i < n; i++)
		*r += a[i]/n;
	return 1;
}

// instance of "univector_function"
static int vector_sum(float *r, float *a, int n)
{
	*r = 0;
	for (int i = 0; i < n; i++)
		*r += a[i];
	return 1;
}

// instance of "univector_function"
static int vector_min(float *r, float *a, int n)
{
	*r = *a;
	for (int i = 1; i < n; i++)
		*r = fmin(*r, a[i]);
	return 1;
}

// instance of "univector_function"
static int vector_max(float *r, float *a, int n)
{
	*r = *a;
	for (int i = 1; i < n; i++)
		*r = fmax(*r, a[i]);
	return 1;
}

// instance of "univector_function"
static int vector_mul(float *r, float *a, int n)
{
	*r = 1;
	for (int i = 0; i < n; i++)
		*r *= a[i];
	return 1;
}

// instance of "univector_function"
static int vector_norm(float *r, float *a, int n)
{
	*r = 0;
	for (int i = 0; i < n; i++)
		*r = hypot(*r, a[i]);
	return 1;
}

// instance of "univector_function"
static int vector_dimension(float *r, float *a, int n)
{
	(void)a;
	*r = n;
	return 1;
}

// instance of "univector_function"
static int vector_rgb2gray(float *r, float *a, int n)
{
	if (n != 3) fail("vgray needs a 3-vector");
	*r = 0.299 * a[0] + 0.587 * a[1] + 0.114 * a[2];
	return 1;
}

// instance of "univector_function"
static int vector_colorsign(float *r, float *a, int n)
{
	if (n != 1) fail("vcolorsign needs a scalar");
	r[0] = fabs(*a) * (*a < 0);
	r[1] = fabs(*a) * (*a > 0);
	r[2] = 0;
	return 3;
}

// table of all functions (local and from math.h) {{{1
static struct predefined_function {
	void (*f)(void);
	char *name;
	int nargs;
	float value;
} global_table_of_predefined_functions[] = {
#define REGISTER_FUNCTION(x,n) {(void(*)(void))x, #x, n, 0}
#define REGISTER_FUNCTIONN(x,xn,n) {(void(*)(void))x, xn, n, 0}
//#define REGISTER_FUNCTIONC(x,n) {(void(*)(void))x, "complex#x, n, 0}
	REGISTER_FUNCTION(acos,1),
	REGISTER_FUNCTION(acosh,1),
	REGISTER_FUNCTION(asin,1),
	REGISTER_FUNCTION(asinh,1),
	REGISTER_FUNCTION(atan,1),
	REGISTER_FUNCTION(atanh,1),
	REGISTER_FUNCTION(cbrt,1),
	REGISTER_FUNCTION(ceil,1),
	REGISTER_FUNCTION(cos,1),
	REGISTER_FUNCTION(cosh,1),
	REGISTER_FUNCTION(erf,1),
	REGISTER_FUNCTION(erfc,1),
	REGISTER_FUNCTION(exp,1),
	REGISTER_FUNCTION(exp2,1),
	REGISTER_FUNCTION(expm1,1),
	REGISTER_FUNCTION(fabs,1),
	REGISTER_FUNCTION(floor,1),
	REGISTER_FUNCTION(lgamma,1),
	REGISTER_FUNCTION(log,1),
	REGISTER_FUNCTION(log10,1),
	REGISTER_FUNCTION(log1p,1),
	REGISTER_FUNCTION(log2,1),
	REGISTER_FUNCTION(logb,1),
	REGISTER_FUNCTION(nearbyint,1),
	REGISTER_FUNCTION(rint,1),
	REGISTER_FUNCTION(round,1),
	REGISTER_FUNCTION(sin,1),
	REGISTER_FUNCTION(sinh,1),
	REGISTER_FUNCTION(sqrt,1),
	REGISTER_FUNCTION(tan,1),
	REGISTER_FUNCTION(tanh,1),
	REGISTER_FUNCTION(tgamma,1),
	REGISTER_FUNCTION(trunc,1),
	REGISTER_FUNCTION(inftozero,1),
	REGISTER_FUNCTION(nantozero,1),
	REGISTER_FUNCTION(notfintozero,1),
	REGISTER_FUNCTION(force_finite,1),
	REGISTER_FUNCTION(force_normal,1),
	REGISTER_FUNCTION(atan2,2),
	REGISTER_FUNCTION(copysign,2),
	REGISTER_FUNCTION(fdim,2),
	REGISTER_FUNCTION(fma,2),
	REGISTER_FUNCTION(fmax,2),
	REGISTER_FUNCTION(fmin,2),
	REGISTER_FUNCTION(fmod,2),
	REGISTER_FUNCTION(hypot,2),
	REGISTER_FUNCTION(ldexp,2),
	REGISTER_FUNCTION(nextafter,2),
	REGISTER_FUNCTION(nexttoward,2),
	REGISTER_FUNCTION(pow,2),
	REGISTER_FUNCTION(remainder,2),
	REGISTER_FUNCTIONN(quantize_255,"q255",1),
	REGISTER_FUNCTIONN(quantize_easy,"qe",3),
	REGISTER_FUNCTIONN(range,"range",3),
	REGISTER_FUNCTIONN(bound_double,"bound",3),
	REGISTER_FUNCTIONN(pow,"^",2),
	REGISTER_FUNCTIONN(sum_two_doubles,"+",2),
	REGISTER_FUNCTIONN(logic_g,">",2),
	REGISTER_FUNCTIONN(logic_l,"<",2),
	REGISTER_FUNCTIONN(logic_e,"=",2),
	REGISTER_FUNCTIONN(logic_ge,">=",2),
	REGISTER_FUNCTIONN(logic_le,"<=",2),
	REGISTER_FUNCTIONN(logic_ne,"!=",2),
	REGISTER_FUNCTIONN(logic_if,"if",3),
	REGISTER_FUNCTIONN(logic_and,"and",2),
	REGISTER_FUNCTIONN(logic_or,"or",2),
	REGISTER_FUNCTIONN(logic_not,"not",1),
	REGISTER_FUNCTIONN(function_isfinite,"isfinite",1),
	REGISTER_FUNCTIONN(function_isinf,"isinf",1),
	REGISTER_FUNCTIONN(function_isnan,"isnan",1),
	REGISTER_FUNCTIONN(function_isnormal,"isnormal",1),
	REGISTER_FUNCTIONN(function_signbit,"signbit",1),
	REGISTER_FUNCTIONN(divide_two_doubles,"/",2),
	REGISTER_FUNCTIONN(multiply_two_doubles,"*",2),
	REGISTER_FUNCTIONN(substract_two_doubles,"-",2),
	REGISTER_FUNCTIONN(random_uniform,"randu",-1),
	REGISTER_FUNCTIONN(random_normal,"randn",-1),
	REGISTER_FUNCTIONN(random_normal,"randg",-1),
	REGISTER_FUNCTIONN(random_cauchy,"randc",-1),
	REGISTER_FUNCTIONN(random_laplace,"randl",-1),
	REGISTER_FUNCTIONN(random_exponential,"rande",-1),
	REGISTER_FUNCTIONN(random_pareto,"randp",-1),
	REGISTER_FUNCTIONN(random_raw,"rand",-1),
	REGISTER_FUNCTIONN(random_stable,"rands",2),
	REGISTER_FUNCTIONN(from_cartesian_to_polar,"topolar", -2),
	REGISTER_FUNCTIONN(from_polar_to_cartesian,"frompolar", -2),
	REGISTER_FUNCTIONN(complex_exp,"cexp", -2),
#ifdef __STDC_IEC_559_COMPLEX__
	REGISTER_FUNCTIONN(complex_cacos , "cacos", -2),
	REGISTER_FUNCTIONN(complex_cacosh, "cacosh", -2),
	REGISTER_FUNCTIONN(complex_casin , "casin", -2),
	REGISTER_FUNCTIONN(complex_casinh, "casinh", -2),
	REGISTER_FUNCTIONN(complex_catan , "catan", -2),
	REGISTER_FUNCTIONN(complex_catanh, "catanh", -2),
	REGISTER_FUNCTIONN(complex_ccos  , "ccos", -2),
	REGISTER_FUNCTIONN(complex_ccosh , "ccosh", -2),
	REGISTER_FUNCTIONN(complex_cexp  , "ccexp", -2),
	REGISTER_FUNCTIONN(complex_clog  , "clog", -2),
	REGISTER_FUNCTIONN(complex_conj  , "conj", -2),
	REGISTER_FUNCTIONN(complex_cproj , "cproj", -2),
	REGISTER_FUNCTIONN(complex_csin  , "csin", -2),
	REGISTER_FUNCTIONN(complex_csinh , "csinh", -2),
	REGISTER_FUNCTIONN(complex_csqrt , "csqrt", -2),
	REGISTER_FUNCTIONN(complex_ctan  , "ctan", -2),
	REGISTER_FUNCTIONN(complex_ctanh , "ctanh", -2),
#endif
	REGISTER_FUNCTIONN(complex_exp,"cexp", -2),
	REGISTER_FUNCTIONN(complex_product,"cprod", -3),
	REGISTER_FUNCTIONN(matrix_product,"mprod",-5),
	REGISTER_FUNCTIONN(vector_product,"vprod",-5),
	REGISTER_FUNCTIONN(scalar_product,"sprod",-5),
	REGISTER_FUNCTIONN(matrix_determinant,"mdet",-6),
	REGISTER_FUNCTIONN(matrix_transpose,"mtrans",-6),
	REGISTER_FUNCTIONN(matrix_inverse,"minv",-6),
	REGISTER_FUNCTIONN(matrix_trace,"mtrace",-6),
	REGISTER_FUNCTIONN(matrix_3x3to12x9,"mto12x9",-6),
	REGISTER_FUNCTIONN(vector_avg,"vavg",-6),
	REGISTER_FUNCTIONN(vector_sum,"vsum",-6),
	REGISTER_FUNCTIONN(vector_min,"vmin",-6),
	REGISTER_FUNCTIONN(vector_max,"vmax",-6),
	REGISTER_FUNCTIONN(vector_mul,"vmul",-6),
	REGISTER_FUNCTIONN(vector_norm,"vnorm",-6),
	REGISTER_FUNCTIONN(vector_dimension,"vdim",-6),
	REGISTER_FUNCTIONN(vector_rgb2gray,"vgray",-6),
	REGISTER_FUNCTIONN(vector_colorsign,"vcsign",-6),
	REGISTER_FUNCTIONN(hsv2rgb,"hsv2rgb",-6),
	REGISTER_FUNCTIONN(rgb2hsv,"rgb2hsv",-6),
	REGISTER_FUNCTIONN(xyz2rgb,"xyz2rgb",-6),
	REGISTER_FUNCTIONN(rgb2xyz,"rgb2xyz",-6),
#undef REGISTER_FUNCTION
#undef REGISTER_FUNCTIONN
	{NULL, "pi", 0, M_PI},
#ifdef M_E
#define REGISTER_CONSTANT(x) {NULL, #x, 0, x}
	REGISTER_CONSTANT(M_E),
	REGISTER_CONSTANT(M_LOG2E),
	REGISTER_CONSTANT(M_LOG10E),
	REGISTER_CONSTANT(M_LN2),
	REGISTER_CONSTANT(M_LN10),
	//REGISTER_CONSTANT(M_PI),
	REGISTER_CONSTANT(M_PI_2),
	REGISTER_CONSTANT(M_PI_4),
	REGISTER_CONSTANT(M_1_PI),
	REGISTER_CONSTANT(M_2_PI),
	REGISTER_CONSTANT(M_2_SQRTPI),
	REGISTER_CONSTANT(M_SQRT2),
	REGISTER_CONSTANT(M_SQRT1_2),
#undef REGISTER_CONSTANT
#endif
};


static float apply_function(struct predefined_function *f, float *v)
{
	switch(f->nargs) {
	case 0: return f->value;
	case 1: return ((double(*)(double))(f->f))(v[0]);
	case 2: return ((double(*)(double,double))f->f)(v[1], v[0]);
	case 3: return ((double(*)(double,double,double))f->f)(v[2],v[1],v[0]);
	case -1: return ((double(*)(void))(f->f))();
	default: fail("bizarre{%d}", f->nargs);
	}
	//return 0;
}

static int symmetrize_index_inside(int i, int m)
{
	assert( i >= 0 && i < m);
	int r = 0;
	if (i >= m/2) r = i-m;
	if (i < m/2) r = i;
	return r;
}

// the value of colon variables depends on the position within the image
static float eval_colonvar(int w, int h, int i, int j, int c)
{
	double x, y;
	switch(c) {
	case 'i': return i;
	case 'j': return j;
	case 'w': return w;
	case 'h': return h;
	case 'n': return w*h;
	case 'x': return (2.0/(w-1))*i - 1;
	case 'y': return (2.0/(h-1))*j - 1;
	case 'r': return hypot((2.0/(h-1))*j-1,(2.0/(w-1))*i-1);
	case 't': return atan2((2.0/(h-1))*j-1,(2.0/(w-1))*i-1);
	case 'I': return symmetrize_index_inside(i,w);
	case 'J': return symmetrize_index_inside(j,h);
	case 'P': return symmetrize_index_inside(i,w)*2*M_PI/w;
	case 'Q': return symmetrize_index_inside(j,h)*2*M_PI/h;
	case 'L': x = symmetrize_index_inside(i,w);
		  y = symmetrize_index_inside(j,h);
		  return -(x*x+y*y);
	case 'R': x = symmetrize_index_inside(i,w);
		  y = symmetrize_index_inside(j,h);
		  return hypot(x, y);
	case 'W': return w/(2*M_PI);
	case 'H': return h/(2*M_PI);
	default: fail("unrecognized colonvar \":%c\"", c);
	}
}


// struct plambda_program {{{1

struct plambda_token {
	int type;
	float value;         // if type==constant, value
	int index;           // if type==variable, its index
	                     // if type==operator, its index
	int component;       // if type==variable, index of selected component
	int displacement[2]; // if type==variable, relative displacement
	int colonvar;        // if type==colon, the letter

	char *tmphack;       // temporary place for storing the unsorted index

	// only for imageop hack:
	int imageop_operator;
	int imageop_scheme;
};

struct collection_of_varnames {
	int n;
	char *t[PLAMBDA_MAX_TOKENS];
};

struct plambda_program {
	int n;
	struct plambda_token t[PLAMBDA_MAX_TOKENS];
	struct collection_of_varnames var[1];

};


// image statistics {{{1

struct image_stats {
	bool init_simple, init_vsimple, init_ordered, init_vordered;
	float scalar_min, scalar_max, scalar_avg, scalar_med, scalar_sum;
	float scalar_std;
	//float vector_n1min[PLAMBDA_MAX_PIXELDIM]; // exemplar with min L1
	float vector_n2min[PLAMBDA_MAX_PIXELDIM]; // exemplar with min L2
	//float vector_nimin[PLAMBDA_MAX_PIXELDIM]; // exemplar with min Linf
	//float vector_n1max[PLAMBDA_MAX_PIXELDIM];
	float vector_n2max[PLAMBDA_MAX_PIXELDIM];
	//float vector_nimax[PLAMBDA_MAX_PIXELDIM];
	float vector_avg[PLAMBDA_MAX_PIXELDIM];
	float vector_med[PLAMBDA_MAX_PIXELDIM];
	float vector_sum[PLAMBDA_MAX_PIXELDIM];
	bool init_csimple, init_cordered;
	float component_min[PLAMBDA_MAX_PIXELDIM];
	float component_max[PLAMBDA_MAX_PIXELDIM];
	float component_avg[PLAMBDA_MAX_PIXELDIM];
	float component_med[PLAMBDA_MAX_PIXELDIM];
	float component_sum[PLAMBDA_MAX_PIXELDIM];
	float component_std[PLAMBDA_MAX_PIXELDIM];
	float *sorted_samples, *sorted_components[PLAMBDA_MAX_PIXELDIM];

	// original image data, for debugging purposes
	int w, h, pd;
	float *x;
};

// silly function to print the components of a vector to stderr
//static int devec(float *x, int n)
//{
//	for (int i = 0; i < n; i++)
//		fprintf(stderr, "%g%s", x[i], i==n-1?"":", ");
//}
//
//TODO: make a reasonable way to call this debugging function
//static void print_image_stats_to_stderr(struct image_stats *s)
//{
//	int n = s->pd;
//	fprintf(stderr, "image [%dx%d,%d] (%p) stats %p\n", s->w, s->h, n,
//			(void*)s->x, (void*)s);
//	if (s->init_simple) {
//		fprintf(stderr, "\tscalar_min = %g\n", s->scalar_min);
//		fprintf(stderr, "\tscalar_max = %g\n", s->scalar_min);
//		fprintf(stderr, "\tscalar_avg = %g\n", s->scalar_min);
//		fprintf(stderr, "\tscalar_sum = %g\n", s->scalar_min);
//		fprintf(stderr, "\tscalar_std = %g\n", s->scalar_min);
//	}
//	if (s->init_ordered) {
//		fprintf(stderr, "\tscalar_med = %g\n", s->scalar_med);
//	}
//	if (s->init_vsimple) {
//		fprintf(stderr,"\tvector_sum = "); devec(s->vector_sum,n);
//		fprintf(stderr,"\tvector_avg = "); devec(s->vector_avg,n);
//		fprintf(stderr,"\tvector_n2min = "); devec(s->vector_n2min,n);
//		fprintf(stderr,"\tvector_n2max = "); devec(s->vector_n2max,n);
//	}
//	if (s->init_csimple) {
//		fprintf(stderr,"\tcomponent_min = "); devec(s->component_min,n);
//		fprintf(stderr,"\tcomponent_max = "); devec(s->component_max,n);
//		fprintf(stderr,"\tcomponent_avg = "); devec(s->component_avg,n);
//		fprintf(stderr,"\tcomponent_sum = "); devec(s->component_sum,n);
//		fprintf(stderr,"\tcomponent_std = "); devec(s->component_std,n);
//	}
//	if (s->init_cordered) {
//		fprintf(stderr,"\tcomponent_med = "); devec(s->component_med,n);
//	}
//}

struct linear_statistics {
	float min, max, avg, avgnz, sum, std;
	int n, rns, rnz, nnan, ninf;
};

static void compute_linstats(struct linear_statistics *s,
		float *x, int n, int stride, int offset)
{
	int rns = 0, rnz = 0, nnan = 0, ninf = 0;
	float min = INFINITY, max = -INFINITY;
	long double avg = 0, avgnz = 0, sumsq = 0;
	for (int i = 0; i < n; i++) {
		float y = x[i*stride + offset];
		if (isnan(y)) {
			nnan += 1;
			continue;
		}
		if (!isfinite(y)) ninf += 1;
		if (y < min) min = y;
		if (y > max) max = y;
		avg += y;
		rns += 1;
		if (y) {
			avgnz += y;
			rnz += 1;
		}
	}
	s->sum = avg;
	avg /= rns; avgnz /= rnz;
	for (int i = 0; i < n; i++) {
		float y = x[i*stride + offset];
		if (!isnan(y))
			sumsq += (y - avg) * (y - avg);
	}
	s->std = sqrt(sumsq / rns);
	s->min=min; s->max=max; s->avg=avg; s->avgnz=avgnz;
	s->n=n; s->rns=rns; s->rnz=rnz; s->nnan=nnan; s->ninf=ninf;
}

static void compute_simple_sample_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_simple) return;
	if (w*h > 1) s->init_simple = true;
	struct linear_statistics ls[1];
	compute_linstats(ls, x, w*h*pd, 1, 0);
	s->scalar_min = ls->min;
	s->scalar_max = ls->max;
	s->scalar_avg = ls->avg;
	s->scalar_sum = ls->sum;
	s->scalar_std = ls->std;
}

static void compute_simple_component_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_csimple) return;
	if (w*h > 1) s->init_csimple = true;
	for (int l = 0; l < pd; l++)
	{
		struct linear_statistics ls[1];
		compute_linstats(ls, x, w*h, pd, l);
		s->component_min[l] = ls->min;
		s->component_max[l] = ls->max;
		s->component_avg[l] = ls->avg;
		s->component_std[l] = ls->std;
	}
}

static float euclidean_norm_of_float_vector(const float *x, int n)
{
	if (n == 1) return fabs(x[0]);
	else {
		float r = 0;
		for (int i = 0; i < n; i++)
			r = hypot(r, x[i]);
		return r;
	}
}

static void compute_simple_vector_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_vsimple) return;
	if (w*h > 1) s->init_vsimple = true;
	int np = w * h, rnp = 0;
	float minpixel = INFINITY, maxpixel = -INFINITY;
	long double avgpixel[PLAMBDA_MAX_PIXELDIM];
	int minidx=-1, maxidx=-1;
	for (int j = 0; j < pd; j++)
		avgpixel[j] = 0;
	for (int i = 0; i < np; i++)
	{
		float xnorm = euclidean_norm_of_float_vector(x + pd*i, pd);
		if (isnan(xnorm)) continue;
		if (xnorm < minpixel) { minidx = i; minpixel = xnorm; }
		if (xnorm > maxpixel) { maxidx = i; maxpixel = xnorm; }
		for (int j = 0; j < pd; j++)
			avgpixel[j] += x[pd*i+j];
		rnp += 1;
	}
	//assert(rnp);
	FORI(pd) s->vector_sum[i] = avgpixel[i];
	long double mipi[pd], mapi[pd];
	for (int j = 0; j < pd; j++) {
		mipi[j] = x[minidx*pd+j];
		mapi[j] = x[maxidx*pd+j];
		avgpixel[j] /= rnp;
	}
	FORI(pd) s->vector_n2min[i] = mipi[i];
	FORI(pd) s->vector_n2max[i] = mapi[i];
	FORI(pd) s->vector_avg[i] = avgpixel[i];
}

static int compare_floats(const void *aa, const void *bb)
{
	const float *a = (const float *)aa;
	const float *b = (const float *)bb;
	return (*a > *b) - (*a < *b);
}

static void compute_ordered_sample_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_ordered) return;
	if (w*h > 1) s->init_ordered = true;
	int ns = w * h * pd;
	s->sorted_samples = xmalloc(ns*sizeof(float));
	FORI(ns) s->sorted_samples[i] = x[i];
	qsort(s->sorted_samples, ns, sizeof(float), compare_floats);
	s->scalar_med = s->sorted_samples[ns/2];
}

static void compute_ordered_component_stats(struct image_stats *s,
		float *x, int w, int h, int pd)
{
	if (s->init_cordered) return;
	if (w*h > 1) s->init_cordered = true;
	int ns = w * h;
	float *t = xmalloc(pd*ns*sizeof(float));
	for (int l = 0; l < pd; l++)
	{
		s->sorted_components[l] = t + l*ns;
		FORI(ns) s->sorted_components[l][i] = x[i*pd+l];
		qsort(s->sorted_components[l],ns,sizeof(float),compare_floats);
		s->component_med[l] = s->sorted_components[l][ns/2];
	}
}

//static void compute_ordered_vector_stats(struct image_stats *s,
//		float *x, int w, int h, int pd)
//{
//	(void)x;
//	(void)w;
//	(void)h;
//	(void)pd;
//	fail("ordered vector stats not implemented");
//	// there is some bizarre trickery waiting to be coded in here
//	//s->init_vordered = true;
//}

static int bound(int a, int x, int b)
{
	if (b < a) return bound(b, x, a);
	if (x < a) return a;
	if (x > b) return b;
	return x;
}


// the value of magic variables depends on some globally cached data
static int eval_magicvar(float *out, int magic, int img_index, int comp, int qq,
		float *x, int w, int h, int pd) // only needed on the first run
{
	// XXX WARNING : global variables here (leading to non-re-entrant code)
	static bool initt = false;
	static struct image_stats *t = 0;
	//static struct image_stats t[PLAMBDA_MAX_MAGIC];
	if (!initt) {
		t = xmalloc(PLAMBDA_MAX_MAGIC * sizeof*t);
		for (int i = 0; i < PLAMBDA_MAX_MAGIC; i++) {
			t[i].init_simple = false;
			t[i].init_ordered = false;
			t[i].init_vsimple = false;
			t[i].init_vordered = false;
			t[i].init_csimple = false;
			t[i].init_cordered = false;
			t[i].w=w;t[i].h=h;t[i].pd=pd;t[i].x=x; // for debug only
		}
		initt = true;
	}
	//fprintf(stderr, "magic=%c index=%d comp=%d\n",magic,img_index,comp);

	if (img_index >= PLAMBDA_MAX_MAGIC)
		fail("%d magic images is too much for me!", PLAMBDA_MAX_MAGIC);

	struct image_stats *ti = t + img_index;

	if(magic=='i' || magic=='a' || magic=='v' || magic=='s' || magic=='r') {
		if (comp < 0) { // use all samples
			compute_simple_sample_stats(ti, x, w, h, pd);
			switch(magic) {
				case 'i': *out = ti->scalar_min; break;
				case 'a': *out = ti->scalar_max; break;
				case 'v': *out = ti->scalar_avg; break;
				case 's': *out = ti->scalar_sum; break;
				case 'r': *out = ti->scalar_std; break;
				default: fail("this can not happen");
			}
			return 1;
		} else { // use samples from the specified component
			compute_simple_component_stats(ti, x, w, h, pd);
			switch(magic) {
				case 'i': *out = ti->component_min[comp]; break;
				case 'a': *out = ti->component_max[comp]; break;
				case 'v': *out = ti->component_avg[comp]; break;
				case 's': *out = ti->component_sum[comp]; break;
				case 'r': *out = ti->component_std[comp]; break;
				default: fail("this can not happen");
			}
			return 1;
		}
	} else if (magic=='I' || magic=='A' || magic=='V' || magic=='S') {
		compute_simple_vector_stats(ti, x, w, h, pd);
		switch(magic) {
		case 'I': FORI(pd) out[i] = ti->vector_n2min[i]; break;
		case 'A': FORI(pd) out[i] = ti->vector_n2max[i]; break;
		case 'V': FORI(pd) out[i] = ti->vector_avg[i]; break;
		case 'S': FORI(pd) out[i] = ti->vector_sum[i]; break;
		default: fail("this can not happen");
		}
		return pd;
	} else if (magic=='Y' || magic=='E' || magic=='R') {
		compute_simple_component_stats(ti, x, w, h, pd);
		switch(magic) {
		case 'Y': FORI(pd) out[i] = ti->component_min[i]; break;
		case 'E': FORI(pd) out[i] = ti->component_max[i]; break;
		case 'R': FORI(pd) out[i] = ti->component_std[i]; break;
		default: fail("this can not happen");
		}
		return pd;
	} else if (magic == 'm' || magic == 'q') {
		if (comp < 0) { // use all samples
			compute_ordered_sample_stats(ti, x, w, h, pd);
			if (magic == 'm') {
				*out = ti->scalar_med;
				return 1;
			}
			if (magic == 'q') {
				int qpos = round(qq*w*h*pd/100.0);
				qpos = bound(0, qpos, w*h*pd-1);
				*out = ti->sorted_samples[qpos];
				return 1;
			}
		} else {
			compute_ordered_component_stats(ti, x, w, h, pd);
			if (magic == 'm') {
				*out = ti->component_med[comp];
				return 1;
			}
			if (magic == 'q') {
				int qpos = round(qq*w*h/100.0);
				qpos = bound(0, qpos, w*h-1);
				*out = ti->sorted_components[comp][qpos];
				return 1;
			}
		}
	} else if (magic == 'O') {
		compute_ordered_component_stats(ti, x, w, h, pd);
		FORI(pd) {
			int qposi = round(qq*w*h/100.0);
			qposi = bound(0, qposi, w*h-1);
			out[i] = ti->sorted_components[i][qposi];
		}
		return pd;
	} else if (magic == 'W') {
		compute_ordered_component_stats(ti, x, w, h, pd);
		FORI(pd) {
			int qposi = round(qq*(w*h/1000000.0));
			qposi = bound(0, qposi, w*h-1);
			out[i] = ti->sorted_components[i][qposi];
		}
		return pd;
	} else if (magic == '0') {
		compute_ordered_component_stats(ti, x, w, h, pd);
		FORI(pd) {
			int qposi = qq;//round(qq*w*h/1000000.0);
			qposi = bound(0, qposi, w*h-1);
			out[i] = ti->sorted_components[i][qposi];
		}
		return pd;
	} else if (magic == '9') {
		compute_ordered_component_stats(ti, x, w, h, pd);
		FORI(pd) {
			int qposi = w*h-1-qq;//round(qq*w*h/1000000.0);
			qposi = bound(0, qposi, w*h-1);
			out[i] = ti->sorted_components[i][qposi];
		}
		return pd;
	} else
		fail("magic of kind '%c' is not yet implemented", magic);

	return 0;
}

// lexing and parsing {{{1

// if the token resolves to a numeric constant, store it in *x and return true
// otherwise, return false
// if trailing characters are ignored, print a warning message
static bool token_is_number(float *x, const char *t)
{
	char *endptr;
	*x = strtof(t, &endptr);
	if (endptr == t) return false;
	if (*endptr != '\0')
		fprintf(stderr, "TOKEN "
				"WARNING: trailing characters (\"%s\") "
				"ignored "
			       	"in numeric constant\n", endptr);
	return true;
}

// if token is colonvar, return the id
// otherwise, return zero
static int token_is_colonvar(const char *t)
{
	if (t[0] != ':') return 0;
	if (isalpha(t[1]) && t[2]=='\0') return t[1];
	return 0;
}

// if token is a variable definition, return the index
// otherwise, return zero
static int token_is_vardef(const char *t)
{
	if (t[0]=='>' && isdigit(t[1]) && t[1]>'0' && t[2]=='\0')
		return t[1] - '0';
	if (t[0]=='<' && isdigit(t[1]) && t[1]>'0' && t[2]=='\0')
		return -(t[1] - '0');
	return 0;
}


#define PLAMBDA_STACKOP_NO 0
#define PLAMBDA_STACKOP_DEL 1
#define PLAMBDA_STACKOP_DUP 2
#define PLAMBDA_STACKOP_VSPLIT 3
#define PLAMBDA_STACKOP_VMERGE 4
#define PLAMBDA_STACKOP_ROT 5
#define PLAMBDA_STACKOP_VMERGE3 6
#define PLAMBDA_STACKOP_VMERGEALL 7
#define PLAMBDA_STACKOP_NMERGE 10
#define PLAMBDA_STACKOP_INTERLEAVE 11
#define PLAMBDA_STACKOP_DEINTERLEAVE 12
#define PLAMBDA_STACKOP_HALVE 13
#define PLAMBDA_STACKOP_NSPLIT 14

// if token is a stack operation, return its id
// otherwise, return zero
static int token_is_stackop(const char *t)
{
	if (0 == strcmp(t, "del")) return PLAMBDA_STACKOP_DEL;
	if (0 == strcmp(t, "dup")) return PLAMBDA_STACKOP_DUP;
	if (0 == strcmp(t, "rot")) return PLAMBDA_STACKOP_ROT;
	if (0 == strcmp(t, "split")) return PLAMBDA_STACKOP_VSPLIT;
	if (0 == strcmp(t, "merge")) return PLAMBDA_STACKOP_VMERGE;
	if (0 == strcmp(t, "join")) return PLAMBDA_STACKOP_VMERGE;
	if (0 == strcmp(t, "merge3")) return PLAMBDA_STACKOP_VMERGE3;
	if (0 == strcmp(t, "join3")) return PLAMBDA_STACKOP_VMERGE3;
	if (0 == strcmp(t, "mergeall")) return PLAMBDA_STACKOP_VMERGEALL;
	if (0 == strcmp(t, "joinall")) return PLAMBDA_STACKOP_VMERGEALL;
	if (0 == strcmp(t, "njoin")) return PLAMBDA_STACKOP_NMERGE;
	if (0 == strcmp(t, "nmerge")) return PLAMBDA_STACKOP_NMERGE;
	if (0 == strcmp(t, "interleave")) return PLAMBDA_STACKOP_INTERLEAVE;
	if (0 == strcmp(t, "deinterleave")) return PLAMBDA_STACKOP_DEINTERLEAVE;
	if (0 == strcmp(t, "halve")) return PLAMBDA_STACKOP_HALVE;
	if (0 == strcmp(t, "nsplit")) return PLAMBDA_STACKOP_NSPLIT;
	return 0;
}

// if the token is a valid word, return its length
//         and if the token is followed by modifiers, fill *endptr
// otherwise, return zero
static int token_is_word(const char *t, const char **endptr)
{
	*endptr = NULL;
	if ((*t=='+'||*t=='-'||*t=='/'||*t=='^'||*t=='*'||*t=='>'||*t=='<'||*t=='=')&&t[1]=='\0')
		return 1;
	if (!isalpha(t[0])) {
		return 0;
	}
	int n = 1;
	while (t[n]) {
		if  (!(isalnum(t[n])||t[n]=='_')) {
			*endptr = t+n;
			//return (t[n]=='('||t[n]=='['||t[n]=='%'||t[n]==';')?n:0;
			return ispunct(t[n]) ? n : 0;
		}
		n += 1;
	}
	return n;
}

static int word_is_predefined(const char *id)
{
	int n = sizeof(global_table_of_predefined_functions)/
		sizeof(global_table_of_predefined_functions[0]);
	struct predefined_function *r = global_table_of_predefined_functions;
	FORI(n)
		if (0 == strcmp(r[i].name, id))
			return i;
	return -1;
}

// fills the modifiers with their defined values, otherwise with the default
static void parse_modifiers(const char *mods,
		int *ocomp, int *odx, int *ody, int *omagic)
{
	*ocomp = -1;
	*odx = 0;
	*ody = 0;
	*omagic = 0;
	int comp, dx, dy;
	char magic;
	// NOTE: the order of the following conditions is important for a
	// correct parsing
	if (!mods) {
		return;
	} else if (3 == sscanf(mods, "[%d](%d,%d)", &comp, &dx, &dy)) {
		*odx = dx;
		*ody = dy;
		*ocomp = comp;
	 	return;
	} else if (3 == sscanf(mods, "(%d,%d)[%d]", &dx, &dy, &comp)) {
		*odx = dx;
		*ody = dy;
		*ocomp = comp;
	 	return;
	} else if (2 == sscanf(mods, "(%d,%d)", &dx, &dy)) {
		*odx = dx;
		*ody = dy;
	 	return;
	} else if (2 == sscanf(mods, "%%%c%d", &magic, &dx)) {
		*omagic = magic;
		*odx = dx;
		return;
	} else if (1 == sscanf(mods, "%%%c", &magic)) {
		*omagic = magic;
		return;
	} else if (3 == sscanf(mods, "[%d]%%%c%d", &comp, &magic, &dx)) {
		*omagic = magic;
		*ocomp = comp;
		*odx = dx;
		return;
	} else if (2 == sscanf(mods, "[%d]%%%c", &comp, &magic)) {
		*omagic = magic;
		*ocomp = comp;
		return;
	} else if (1 == sscanf(mods, "[%d]", &comp)) {
		*ocomp = comp;
		return;
	}
}

static bool hasprefix(const char *s, const char *p)
{
	return s == strstr(s, p);
}

static bool hassuffix(const char *s, const char *suf)
{
	int len_s = strlen(s);
	int len_suf = strlen(suf);
	if (len_s < len_suf)
		return false;
	return 0 == strcmp(suf, s + (len_s - len_suf));
}

static void parse_imageop(const char *s, int *op, int *scheme)
{
	*op = IMAGEOP_IDENTITY;
	if (false) ;
	else if (hasprefix(s, "xx")) *op = IMAGEOP_XX;
	else if (hasprefix(s, "yy")) *op = IMAGEOP_YY;
	else if (hasprefix(s, "xy")) *op = IMAGEOP_XY;
	else if (hasprefix(s, "yx")) *op = IMAGEOP_XY;
	else if (hasprefix(s, "l")) *op = IMAGEOP_LAP;
	else if (hasprefix(s, "x")) *op = IMAGEOP_X;
	else if (hasprefix(s, "y")) *op = IMAGEOP_Y;
	else if (hasprefix(s, "n")) *op = IMAGEOP_NGRAD;
	else if (hasprefix(s, "g")) *op = IMAGEOP_GRAD;
	else if (hasprefix(s, "d")) *op = IMAGEOP_DIV;
	else if (hasprefix(s, "S")) *op = IMAGEOP_SHADOW;
	//else if (hasprefix(s, "k")) *op = IMAGEOP_CURV;
	//else fail("unrecognized comma modifier \",%s\"", s);
	*scheme = SCHEME_SOBEL;
	if (*op == IMAGEOP_XY) *scheme = SCHEME_CENTERED;
	if (false) ;
	else if (hassuffix(s, "f")) *scheme = SCHEME_FORWARD;
	else if (hassuffix(s, "b")) *scheme = SCHEME_BACKWARD;
	else if (hassuffix(s, "c")) *scheme = SCHEME_CENTERED;
	else if (hassuffix(s, "s")) *scheme = SCHEME_SOBEL;
	else if (hassuffix(s, "p")) *scheme = SCHEME_PREWITT;
}

static void collection_of_varnames_init(struct collection_of_varnames *x)
{
	x->n = 0;
}

static int collection_of_varnames_find(struct collection_of_varnames *x,
		const char *s)
{
	FORI(x->n)
		if (0 == strcmp(s, x->t[i]))
			return i;
	return -1;
}

static char *collection_of_varnames_add(struct collection_of_varnames *x,
		const char *s)
{
	char *r;
	int i = collection_of_varnames_find(x, s);
	if (i < 0) {
		if (x->n+1 >= PLAMBDA_MAX_TOKENS)
			fail("caca");
		r = xmalloc(1+strlen(s));
		strcpy(r, s);
		x->t[x->n] = r;
		x->n += 1;
	} else {
		r = x->t[i];
	}
	return r;
}

static void collection_of_varnames_end(struct collection_of_varnames *x)
{
	FORI(x->n)
		free(x->t[i]);
	x->n = 0;
}

static int strcmp_for_qsort(const void *aa, const void *bb)
{
	const char **a = (const char **)aa;
	const char **b = (const char **)bb;
	return strcmp(*a, *b);
}

static void collection_of_varnames_sort(struct collection_of_varnames *x)
{
	qsort(x->t, x->n, sizeof*x->t, strcmp_for_qsort);
}

// this function takes a string which contains one token,
// and compiles the corresponding info into p->t[p->n]
//
// TODO (maybe): split token identification from info gathering
// (this will produce longer code but shorter functions)
static void process_token(struct plambda_program *p, const char *tokke)
{
	char tok[1+strlen(tokke)];             // the string of the token
	strcpy(tok, tokke);
	struct plambda_token *t = p->t + p->n; // the compiled token

	int tok_id;
	const char *tok_end;

	float x;
	if (token_is_number(&x, tok)) {
		t->type = PLAMBDA_CONSTANT;
		t->value = x;
		goto endtok;
	}

	if ((tok_id = token_is_colonvar(tok))) {
		t->type = PLAMBDA_COLONVAR;
		t->colonvar = tok_id;
		goto endtok;
	}

	if ((tok_id = token_is_stackop(tok))) {
		t->type = PLAMBDA_STACKOP;
		t->index = tok_id;
		goto endtok;
	}

	if ((tok_id = token_is_vardef(tok))) {
		t->type = PLAMBDA_VARDEF;
		t->index = tok_id;
		goto endtok;
	}

	if ((token_is_word(tok, &tok_end)))
	{
		int idx = word_is_predefined(tok);
		if (idx < 0) {
			char varname[PLAMBDA_MAX_VARLEN+1];
			int varlen = strlen(tok);
			if (tok_end) varlen = tok_end-tok;
			if (varlen >= PLAMBDA_MAX_VARLEN)
				varlen = PLAMBDA_MAX_VARLEN;
			FORI(varlen) varname[i] = tok[i];
			varname[varlen] = '\0';
			int comp, disp[2], magic;
			t->tmphack =collection_of_varnames_add(p->var, varname);
			//fprintf(stderr, "varname, mods = \"%s\" , \"%s\"\n", varname, tok_end);
			parse_modifiers(tok_end, &comp, disp, disp+1, &magic);
			//fprintf(stderr, "comp=%d disp=[%d %d] magic=%d\n",
			//		comp, *disp, disp[1], magic);
			t->type = comp<0 ? PLAMBDA_VECTOR : PLAMBDA_SCALAR;
			if (magic > 0) {
				t->type = PLAMBDA_MAGIC;
				t->colonvar = magic;
			}
			char *iopos = tok_end ? strchr(tok_end, ';') : 0;
			if (!iopos && tok_end) iopos = strchr(tok_end, ',');
			if (iopos) {
				t->type = PLAMBDA_IMAGEOP; // comma operator
				parse_imageop(1+iopos, &t->imageop_operator,
							&t->imageop_scheme);
				//fprintf(stderr, "imageop %d %d\n", t->imageop_operator, t->imageop_scheme);
			}
			t->component = comp;
			t->displacement[0] = disp[0];
			t->displacement[1] = disp[1];
		} else {
			//struct predefined_function *f =
			//	global_table_of_predefined_functions + idx;
			t->type = PLAMBDA_OPERATOR;
			t->index = idx;
		}
		goto endtok;
	}

endtok:
	p->n += 1;
}

// this function updates the indexes of a
// collection of variables which is sorted in alphabetical order
static void update_variable_indices(struct plambda_program *p)
{
	FORI(p->n)
	{
		struct plambda_token *t = p->t + i;
		if (t->type == PLAMBDA_SCALAR || t->type == PLAMBDA_VECTOR
				|| t->type == PLAMBDA_MAGIC
				|| t->type == PLAMBDA_IMAGEOP)
		{
			t->index = collection_of_varnames_find(p->var,
								t->tmphack);
			if (t->index < 0)
				fail("unexpected bad variable \"%s\"",
								t->tmphack);
		}
	}
}

static void plambda_compile_program(struct plambda_program *p, const char *str)
{
	char s[1+strlen(str)], *spacing = " \n\t";
	strcpy(s, str);

	collection_of_varnames_init(p->var);
	p->n = 0;
	char *tok = strtok(s, spacing);
	while (tok) {
		//fprintf(stderr, "TOK \"%s\"\n", tok);
		process_token(p, tok);
		tok = strtok(NULL, spacing);
	}

	collection_of_varnames_sort(p->var);
	update_variable_indices(p);
}

static const char *arity(struct predefined_function *f)
{
	switch(f->nargs) {
	case 0: return "0-ary";
	case 1: return "unary";
	case 2: return "binary";
	case 3: return "ternary";
	case -1: return "strange";
	case -2: return "strange2";
	case -3: return "strange3";
	default: return "unrecognized";
	}
}

inline
static void print_compiled_program(struct plambda_program *p)
{
	fprintf(stderr, "COMPILED PROGRAM OF %d TOKENS:\n", p->n);
	FORI(p->n) {
		struct plambda_token *t = p->t + i;
		fprintf(stderr, "TOKEN[%d]: ", i);
		if (t->type == PLAMBDA_CONSTANT)
			fprintf(stderr, "constant %g", t->value);
		if (t->type == PLAMBDA_COLONVAR)
			fprintf(stderr, "colonvar \"%c\"", t->colonvar);
		if (t->type == PLAMBDA_VECTOR) {
			fprintf(stderr, "variable vector %d \"%s\"",
					t->index, p->var->t[t->index]);
			fprintf(stderr, ", displacement (%d,%d)",
					t->displacement[0], t->displacement[1]);
			if (t->component < -1)
				fprintf(stderr, ", component %d", t->component);
		}
		if (t->type == PLAMBDA_SCALAR) {
			fprintf(stderr, "variable scalar %d \"%s\"",
					t->index, p->var->t[t->index]);
			fprintf(stderr, ", displacement (%d,%d)",
					t->displacement[0], t->displacement[1]);
			fprintf(stderr, ", component %d", t->component);
		}
		if (t->type == PLAMBDA_OPERATOR) {
			struct predefined_function *f =
				global_table_of_predefined_functions+t->index;
			fprintf(stderr, "%s operator %s", arity(f), f->name);
		}
		if (t->type == PLAMBDA_STACKOP)
			fprintf(stderr, "stack manipulation");
		if (t->type == PLAMBDA_VARDEF) {
			fprintf(stderr, "register variable %s %d",
					t->index<0 ? "read" : "write",
					abs(t->index));
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "The program uses %d variables:\n", p->var->n);
	FORI(p->var->n)
		fprintf(stderr, "VARIABLE[%d] = \"%s\"\n", i, p->var->t[i]);
}

// evaluation {{{1

// struct value_vstack {{{2

// stack of vectorial (or possibly scalar) values
//
// This structure actually holds the state of a plambda program as it is being
// executed
struct value_vstack {
	// stack
	int n;
	int d[PLAMBDA_MAX_TOKENS];
	float t[PLAMBDA_MAX_TOKENS][PLAMBDA_MAX_PIXELDIM];

	// registers 0 to 9
	int regn[10];
	float regv[10][PLAMBDA_MAX_PIXELDIM];

	//// stack of "forgotten" values (accessible via lastx, lasty, lastz...)
	//int lasttop, lastn[10];
	//float lastv[10][PLAMBDA_MAX_PIXELDIM];
};

// struct value_vstack manipulation {{{2

static int vstack_pop_vector(float *val, struct value_vstack *s)
{
	if (s->n > 0) {
		s->n -= 1;
		int d = s->d[s->n];
		if (val) FORI(d) val[i] = s->t[s->n][i];
		// TODO: push this value onto a "lastx" register, à la HP
		return d;
	} else fail("popping from empty stack");
}

static void vstack_push_vector(struct value_vstack *s, float *v, int n)
{
	if (s->n+1 < PLAMBDA_MAX_TOKENS) {
		s->d[s->n] = n;
		FORI(n)
			s->t[s->n][i] = v[i];
		s->n += 1;
	} else fail("full stack");
}

static void vstack_push_scalar(struct value_vstack *s, float x)
{
	vstack_push_vector(s, &x, 1);
	//if (s->n+1 < PLAMBDA_MAX_TOKENS) {
	//	s->d[s->n] = 1;
	//	s->t[s->n][0] = x;
	//	s->n += 1;
	//} else fail("full stack");
}

//static void vstack_print(FILE *f, struct value_vstack *s)
//{
//	FORI(s->n) {
//		fprintf(f, "STACK[%d/%d]: {%d}", 1+i, s->n, s->d[i]);
//		FORJ(s->d[i])
//			fprintf(f, " %g", s->t[i][j]);
//		fprintf(f, "\n");
//	}
//}


// function application (strange cases) {{{2
// XXX TODO: refactor the following strange cases into the general setup

// a function that takes no arguments but must be called nonetheless
static void treat_strange_case(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -1);
	float r = apply_function(f, NULL);
	vstack_push_vector(s, &r, 1);
}

// pop a 2-vector x from the stack and push f(x) as a 2-vector
static void treat_strange_case2(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -2);
	float v[PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	int n = vstack_pop_vector(v, s);
	if (n != 2) fail("function \"%s\" requires a 2-vector", f->name);
	((void(*)(float*,float*))(f->f))(r, v);
	vstack_push_vector(s, r, n);
}

#ifndef ODDP
#define ODDP(x) ((x)&1)
#endif
#ifndef EVENP
#define EVENP(x) (!((x)&1))
#endif

// pop two 2-vectors x,y from the stack and push f(x,y) as a 2-vector
// NOTE: works also on vectors of even length, doing the "obvious" thing
static void treat_strange_case3(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -3);
	float a[PLAMBDA_MAX_PIXELDIM];
	float b[PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	int na = vstack_pop_vector(a, s);
	int nb = vstack_pop_vector(b, s);
	if (ODDP(na)) fail("function \"%s\" requires 2-vectors", f->name);
	if (ODDP(nb)) fail("function \"%s\" requires 2-vectors", f->name);
	void (*ff)(float*,float*,float*) = (void(*)(float*,float*,float*))f->f;
	int ca = na / 2;
	int cb = nb / 2;
	int n = 0;
	if (ca == cb) {
		n = na;
		for (int i = 0; i < ca; i++)
			ff(r+2*i, a+2*i, b+2*i);
	} else if (ca == 1) {
		n = nb;
		for (int i = 0; i < cb; i++)
			ff(r+2*i, a, b+2*i);
	} else if (cb == 1) {
		n = na;
		for (int i = 0; i < ca; i++)
			ff(r+2*i, a+2*i, b);
	} else fail("function \"%s\" can not operate on lengths (%d,%d)",
			f->name, na, nb);
	assert(n);
	//if (na != nb) fail("this can not happen");
	//((void(*)(float*,float*,float*))(f->f))(r, a, b);
	vstack_push_vector(s, r, n);
}

// codifies a function R^a x R^b -> R^c
// the dimensions must be smaller than, say, 40
static void treat_bivector_function(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -5);
	float a[PLAMBDA_MAX_PIXELDIM];
	float b[PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	int nb = vstack_pop_vector(b, s);
	int na = vstack_pop_vector(a, s);
	int nr = ((int(*)(float*,float*,float*,int,int))(f->f))(r,a,b,na,nb);
	vstack_push_vector(s, r, nr);
}

// codifies a function R^a -> R^b
// the dimensions must be smaller than, say, 40
static void treat_univector_function(struct value_vstack *s,
		struct predefined_function *f)
{
	assert(f->nargs == -6);
	float a[PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	int na = vstack_pop_vector(a, s);
	int nr = ((int(*)(float*,float*,int))(f->f))(r,a,na);
	vstack_push_vector(s, r, nr);
}

// function application (general case) {{{2

// this function is complicated because it contains the scalar+vector
// semantics, which is complicated
static void vstack_apply_function(struct value_vstack *s,
					struct predefined_function *f)
{
	if (f->nargs == -1) {treat_strange_case(s,f); return;}
	if (f->nargs == -2) {treat_strange_case2(s,f); return;}
	if (f->nargs == -3) {treat_strange_case3(s,f); return;}
	//if (f->nargs < 0) {non_vectorialized_vectorial_operation(s,f);return;}
	if (f->nargs == -5) {treat_bivector_function(s,f); return;}
	if (f->nargs == -6) {treat_univector_function(s,f); return;}
	int d[f->nargs], rd = 1;
	float v[f->nargs][PLAMBDA_MAX_PIXELDIM];
	float r[PLAMBDA_MAX_PIXELDIM];
	FORI(f->nargs)
		d[i] = vstack_pop_vector(v[i], s);
	FORI(f->nargs) // the d[i] which are larger than one must be equal
		if (d[i] > 1) {
			if (rd > 1 && d[i] != rd)
				fail("can not vectorize (%d %d)", rd, d[i]);
			else
				rd = d[i];
		}
	if (rd > 1)
		FORI(f->nargs)
			if (d[i] == 1)
				FORL(rd)
					v[i][l] = v[i][0];
	FORL(rd) {
		float a[f->nargs];
		FORI(f->nargs)
			a[i] = v[i][l];
		r[l] = apply_function(f, a);
	}
	vstack_push_vector(s, r, rd);
}


// vstack_process_op (evaluate a program instruction) {{{2

static void vstack_process_op(struct value_vstack *s, int opid)
{
	switch(opid) {
	case PLAMBDA_STACKOP_DEL:
		vstack_pop_vector(NULL, s);
		break;
	case PLAMBDA_STACKOP_DUP: {
		float x[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		vstack_push_vector(s, x, n);
		vstack_push_vector(s, x, n);
				  }
		break;
	case PLAMBDA_STACKOP_VSPLIT: {
		float x[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		FORI(n)
			vstack_push_scalar(s, x[i]);
				     }
		break;
	case PLAMBDA_STACKOP_VMERGE: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int m = vstack_pop_vector(y, s);
		int n = vstack_pop_vector(x, s);
		if (n+m >= PLAMBDA_MAX_PIXELDIM)
			fail("merging vectors results in large vector");
		FORI(m)
			x[n+i] = y[i];
		vstack_push_vector(s, x, n+m);
				     }
		break;
	case PLAMBDA_STACKOP_VMERGE3: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		float z[PLAMBDA_MAX_PIXELDIM];
		int nz = vstack_pop_vector(z, s);
		int ny = vstack_pop_vector(y, s);
		int nx = vstack_pop_vector(x, s);
		if (nx+ny+nz >= PLAMBDA_MAX_PIXELDIM)
			fail("merging vectors results in large vector");
		FORI(ny) x[nx+i] = y[i];
		FORI(nz) x[nx+ny+i] = z[i];
		vstack_push_vector(s, x, nx+ny+nz);
				     }
		break;
	case PLAMBDA_STACKOP_ROT: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		int m = vstack_pop_vector(y, s);
		vstack_push_vector(s, x, n);
		vstack_push_vector(s, y, m);
		break;
				  }
	case PLAMBDA_STACKOP_VMERGEALL:
		fail("mergeall not implemented");
		break;
	case PLAMBDA_STACKOP_NMERGE: {
		float nn[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(nn, s);
		if (n != 1 || nn[0] < 1 || round(nn[0]) != nn[0]
				|| nn[0] >= PLAMBDA_MAX_PIXELDIM)
			fail("can not nmerge \"%g\" things", nn[0]);
		n = nn[0];
		float x[n][PLAMBDA_MAX_PIXELDIM], y[PLAMBDA_MAX_PIXELDIM];
		int d[n], sdi = 0;
		FORI(n) {
			int j = n-i-1;
			d[j] = vstack_pop_vector(x[j], s);
			sdi += d[j];
		}
		if (sdi >= PLAMBDA_MAX_PIXELDIM)
			fail("merging vectors results in large vector");
		int cx = 0;
		FORI(n) FORJ(d[i]) y[cx++] = x[i][j];
		assert(cx == sdi);
		vstack_push_vector(s, y, sdi);
		break;
				     }
	case PLAMBDA_STACKOP_INTERLEAVE: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (ODDP(n)) fail("can not interleave an odd number %d!", n);
		FORI(n/2) {
			y[2*i] = x[i];
			y[2*i+1] = x[i+n/2];
		}
		vstack_push_vector(s, y, n);
		break;
					 }
	case PLAMBDA_STACKOP_DEINTERLEAVE: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (ODDP(n)) fail("can not deinterleave an odd number %d!", n);
		FORI(n/2) {
			y[i] = x[2*i];
			y[i+n/2] = x[2*i+1];
		}
		vstack_push_vector(s, y, n);
		break;
					 }
	case PLAMBDA_STACKOP_HALVE: {
		float x[PLAMBDA_MAX_PIXELDIM];
		float y[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(x, s);
		if (ODDP(n))
			fail("can not halve a vector of odd length %d!", n);
		FORI(n/2)
			y[i] = x[i+n/2];
		vstack_push_vector(s, x, n/2);
		vstack_push_vector(s, y, n/2);
		break;
					 }
	case PLAMBDA_STACKOP_NSPLIT: {
		float nn[PLAMBDA_MAX_PIXELDIM];
		int n = vstack_pop_vector(nn, s);
		if (n != 1 || nn[0] < 1 || round(nn[0]) != nn[0])
			fail("can not split into \"%g\" things", nn[0]);
		int nparts = nn[0];
		float x[PLAMBDA_MAX_PIXELDIM];
		n = vstack_pop_vector(x, s);
		if (0 != n % nparts)
			fail("can not split \"%d\" in \"%d\" parts", n, nparts);
		int partsize = n/nparts;
		assert(n == nparts*partsize);
		float y[nparts][partsize];
		FORI(nparts) FORJ(partsize)
			y[i][j] = x[i*partsize+j];
		FORI(nparts)
			vstack_push_vector(s, y[i], partsize);
				     }
		break;
	default:
		fail("impossible condition (stackop %d)", opid);
	}
}


// run_program_vectorially_at {{{2
#include "getpixel.c"

SMART_PARAMETER_SILENT(PLAMBDA_GETPIXEL,-1)
static float getsample_cfg(float *x, int w, int h, int pd, int i, int j, int l)
{
	getsample_operator p = get_sample_operator(getsample_1);
	int option = PLAMBDA_GETPIXEL();
	switch (option) {
	case -1: break;
	case 0: p = getsample_0; break;
	case 1: p = getsample_1; break;
	case 2: p = getsample_2; break;
	case 3: p = getsample_per; break;
	case 4: p = getsample_nan; break;
	default: fail("unrecognized PLAMBDA_GETPIXEL value %d", option);
	}
	return p(x, w, h, pd, i, j, l);
}

#define H 0.5
#define Q 0.25
#define O 0.125
static float stencil_3x3_identity[9] =  {0,0,0,  0,1,0, 0,0,0};
static float stencil_3x3_dx_forward[9] =  {0,0,0,  0,-1,1, 0,0,0};
static float stencil_3x3_dx_backward[9] = {0,0,0,  -1,1,0, 0,0,0};
static float stencil_3x3_dx_centered[9] = {0,0,0,  -H,0,H, 0,0,0};
static float stencil_3x3_dy_forward[9] =  {0,0,0,  0,-1,0, 0,1,0};
static float stencil_3x3_dy_backward[9] = {0,-1,0, 0,1,0,  0,0,0};
static float stencil_3x3_dy_centered[9] = {0,-H,0, 0,0,0,  0,H,0};
static float stencil_3x3_dy_sobel[9] = {-O,-2*O,-O,  0,0,0, O,2*O,O};
static float stencil_3x3_dx_sobel[9] = {-O,0,O,  -2*O,0,2*O, -O,0,O};
static float stencil_3x3_dx_prewitt[9] =  {0,0,0,  0,-H,H, 0,-H,H};
static float stencil_3x3_dy_prewitt[9] =  {0,0,0,  0,-H,-H, 0,H,H};
static float stencil_3x3_laplace[9] =  {0,1,0,  1,-4,1, 0,1,0};
static float stencil_3x3_dxx[9] =  {0,0,0,  1,-2,1, 0,0,0};
static float stencil_3x3_dyy[9] =  {0,1,0,  0,-2,0, 0,1,0};
static float stencil_3x3_dxy_4point[9] =  {-Q,0,Q,  0,0,0, Q,0,-Q};
static float stencil_3x3_dxy_7point[9] =  {0,-H,H,  -H,1,-H, H,-H,0};
#undef H
#undef Q
#undef O

static float *get_stencil_3x3(int operator, int scheme)
{
	switch(operator) {
	case IMAGEOP_LAP: return stencil_3x3_laplace;
	case IMAGEOP_XX: return stencil_3x3_dxx;
	case IMAGEOP_YY: return stencil_3x3_dyy;
	case IMAGEOP_XY: { switch(scheme) {
			 case SCHEME_CENTERED: return stencil_3x3_dxy_4point;
			 case SCHEME_SOBEL: return stencil_3x3_dxy_7point;
			 default: fail("unrecognized stencil,xy scheme %d");
			 }
		}
	case IMAGEOP_X: { switch(scheme) {
			case SCHEME_FORWARD: return stencil_3x3_dx_forward;
			case SCHEME_BACKWARD: return stencil_3x3_dx_backward;
			case SCHEME_CENTERED: return stencil_3x3_dx_centered;
			case SCHEME_SOBEL: return stencil_3x3_dx_sobel;
			case SCHEME_PREWITT: return stencil_3x3_dx_prewitt;
			default: fail("unrecognized stencil,x scheme %d");
			}
		}
	case IMAGEOP_Y: { switch(scheme) {
			case SCHEME_FORWARD: return stencil_3x3_dy_forward;
			case SCHEME_BACKWARD: return stencil_3x3_dy_backward;
			case SCHEME_CENTERED: return stencil_3x3_dy_centered;
			case SCHEME_SOBEL: return stencil_3x3_dy_sobel;
			case SCHEME_PREWITT: return stencil_3x3_dy_prewitt;
			default: fail("unrecognized stencil,y scheme %d");
			}
		}
	case IMAGEOP_IDENTITY: return stencil_3x3_identity;
	default: return NULL;//fail("unrecognized stencil operator %d");
	}
}

static float apply_3x3_stencil(float *img, int w, int h, int pd,
		int ai, int aj, int channel, float *s)
{
	assert(s);
	getsample_operator P = getsample_cfg;
	float r = 0;
	for (int i = 0; i < 9; i++)
		r += s[i] * P(img, w, h, pd, ai-1+i%3, aj-1+i/3, channel);
	return r;
}

//static float imageop_scalar_old(float *img, int w, int h, int pd,
//		int ai, int aj, int al, struct plambda_token *t)
//{
//	float *s = get_stencil_3x3(t->imageop_operator, t->imageop_scheme);
//	return apply_3x3_stencil(img, w, h, pd, ai, aj, al, s);
//}

static float imageop_scalar(float *img, int w, int h, int pd,
		int ai, int aj, int al, struct plambda_token *t)
{
	float *s = get_stencil_3x3(t->imageop_operator, t->imageop_scheme);
	if (s)
		return apply_3x3_stencil(img, w, h, pd, ai, aj, al, s);
	else {
		switch(t->imageop_operator) {
		case IMAGEOP_NGRAD:
			{
			float *sx=get_stencil_3x3(IMAGEOP_X,t->imageop_scheme);
			float *sy=get_stencil_3x3(IMAGEOP_Y,t->imageop_scheme);
			float gx = apply_3x3_stencil(img, w,h,pd, ai,aj,al, sx);
			float gy = apply_3x3_stencil(img, w,h,pd, ai,aj,al, sy);
			return hypot(gx, gy);
			}
		default: fail("unrecognized imageop operator %d\n", t->imageop_operator);
		}
	}
	return 0;
}

SMART_PARAMETER_SILENT(SHADOWX,1)
SMART_PARAMETER_SILENT(SHADOWY,1)
SMART_PARAMETER_SILENT(SHADOWZ,1)

static int imageop_vector(float *out, float *img, int w, int h, int pd,
		int ai, int aj, struct plambda_token *t)
{
	float *sx = get_stencil_3x3(IMAGEOP_X, t->imageop_scheme);
	float *sy = get_stencil_3x3(IMAGEOP_Y, t->imageop_scheme);
	switch (t->imageop_operator) {
	case IMAGEOP_GRAD:
		//if (pd != 1) fail("can not yet compute gradient of a vector");
		//out[0] = apply_3x3_stencil(img, w,h,pd, ai,aj,0, sx);
		//out[1] = apply_3x3_stencil(img, w,h,pd, ai,aj,0, sy);
		//return 2;
		for (int l = 0; l < pd; l++) {
			out[2*l+0] = apply_3x3_stencil(img,w,h,pd,ai,aj,l, sx);
			out[2*l+1] = apply_3x3_stencil(img,w,h,pd,ai,aj,l, sy);
		}
		return 2*pd;
	case IMAGEOP_DIV:
		//if (pd!=2)fail("can not compute divergence of a %d-vector",pd);
		//float ax = apply_3x3_stencil(img, w,h,pd, ai,aj,0, sx);
		//float by = apply_3x3_stencil(img, w,h,pd, ai,aj,1, sy);
		//out[0] = ax + by;
		//return 1;
		if (pd%2)fail("can not compute divergence of a %d-vector",pd);
		for (int l = 0; l < pd/2; l++) {
			float ax=apply_3x3_stencil(img,w,h,pd,ai,aj,2*l+0,sx);
			float by=apply_3x3_stencil(img,w,h,pd,ai,aj,2*l+1,sy);
			out[l] = ax + by;
		}
		return pd/2;
	case IMAGEOP_SHADOW: {
		if (pd != 1) fail("can not yet compute shadow of a vector");
		float vdx[3]={1,0,apply_3x3_stencil(img, w,h,pd, ai,aj,0, sx)};
		float vdy[3]={0,1,apply_3x3_stencil(img, w,h,pd, ai,aj,0, sy)};
		//float sun[3] = {-1, -1, 1}, nor[3];
		float sun[3] = {-SHADOWX(), -SHADOWY(), SHADOWZ()}, nor[3];
		vector_product(nor, vdx, vdy, 3, 3);
		return scalar_product(out, nor, sun, 3, 3);
		}
	default: fail("unrecognized imageop %d\n", t->imageop_operator);
	}
}


// compute the requested imageop at the given point
static int imageop(float *out, float *img, int w, int h, int pd,
				int ai, int aj, struct plambda_token *t)
{
	int retval = 1;
	int pi = ai + t->displacement[0];
	int pj = aj + t->displacement[1];
	int channel = t->component;
	if (t->imageop_operator > 1000)
		return imageop_vector(out, img, w, h, pd, pi, pj, t);
	if (channel < 0) { // means the whole of it
		retval = pd;
		FORL(pd)
			out[l] = imageop_scalar(img, w, h, pd, pi, pj, l, t);
	} else
		*out = imageop_scalar(img, w, h, pd, pi, pj, channel, t);
	return retval;
}

// returns the dimension of the output
static int run_program_vectorially_at(float *out, struct plambda_program *p,
		float **val, int *w, int *h, int *pd, int ai, int aj)
{
	getsample_operator P = getsample_cfg;
	struct value_vstack s[1];
	s->n = 0;
	FORI(p->n) {
		struct plambda_token *t = p->t + i;
		switch(t->type) {
		case PLAMBDA_STACKOP:
			vstack_process_op(s, t->index);
			break;
		case PLAMBDA_CONSTANT:
			vstack_push_scalar(s, t->value);
			break;
		case PLAMBDA_COLONVAR: {
			int imw = w ? *w : 1;
			int imh = h ? *h : 1;
			/*hack*/if ('X' == t->colonvar) {
				float v[2] = {ai, aj};
				vstack_push_vector(s, v, 2);
				break;
			}
			float x = eval_colonvar(imw, imh, ai, aj, t->colonvar);
			vstack_push_scalar(s, x);
			break;
				       }
		case PLAMBDA_SCALAR: {
			float *img = val[t->index];
			int imw = w ? w[t->index] : 1;
			int imh = h ? h[t->index] : 1;
			int pdv = pd[t->index];
			int dai = ai + t->displacement[0];
			int daj = aj + t->displacement[1];
			int cmp = t->component;
			float x = P(img, imw, imh, pdv, dai, daj, cmp);
			vstack_push_scalar(s, x);
			break;
				     }
		case PLAMBDA_VECTOR: {
			float *img = val[t->index];
			int imw = w ? w[t->index] : 1;
			int imh = h ? h[t->index] : 1;
			int pdv = pd[t->index];
			int dai = ai + t->displacement[0];
			int daj = aj + t->displacement[1];
			float x[pdv];
			if (t->component == -1) { // regular vector
				FORL(pdv)
				x[l] = P(img, imw, imh, pdv, dai, daj, l);
				vstack_push_vector(s, x, pdv);
			} else if (t->component == -2 && 0==pdv%2) {// 1st half
				FORL(pdv/2)
				x[l] = P(img, imw, imh, pdv, dai, daj, l);
				vstack_push_vector(s, x, pdv/2);
			} else if (t->component == -3 && 0==pdv%2) {// 2nd half
				FORL(pdv/2)
				x[l] = P(img, imw, imh, pdv, dai, daj, pdv/2+l);
				vstack_push_vector(s, x, pdv/2);
			}
				     }
			break;
		case PLAMBDA_IMAGEOP: {
			float *img = val[t->index], lout[PLAMBDA_MAX_PIXELDIM];
			int pdv = pd[t->index];
			int imw = w ? w[t->index] : 1;
			int imh = h ? h[t->index] : 1;
			int rdim = imageop(lout, img, imw, imh, pdv, ai, aj, t);
			vstack_push_vector(s, lout, rdim);
			break;
				      }
		case PLAMBDA_OPERATOR: {
			struct predefined_function *f =
				global_table_of_predefined_functions+t->index;
			vstack_apply_function(s, f);
				       }
			break;
		case PLAMBDA_VARDEF: {
			int n = abs(t->index);
			if (t->index > 0)
				s->regn[n] = vstack_pop_vector(s->regv[n], s);
			if (t->index < 0)
				vstack_push_vector(s, s->regv[n], s->regn[n]);
				     }
			break;
		case PLAMBDA_MAGIC: {
#ifdef _OPENMP
			fail("magic variables are not available in "
					"parallel plambda");
#endif//_OPENMP
			int imw = w ? w[t->index] : 1;
			int imh = h ? h[t->index] : 1;
			int pdv = pd[t->index];
			float *img = val[t->index], x[pdv];
			int rm = eval_magicvar(x, t->colonvar, t->index,
					t->component, t->displacement[0],
					img, imw, imh, pdv);
			vstack_push_vector(s, x, rm);
				    }
			break;
		default:
			fail("unknown tag type %d", t->type);
		}
	}
	return vstack_pop_vector(out, s);
}


// evaluation (higher level) {{{1

static int eval_dim(struct plambda_program *p, float **val, int *pd)
{
	int r = run_program_vectorially_at(NULL, p, val, NULL, NULL, pd, 0, 0);
	return r;
}

// returns the dimension of the output
static int run_program_vectorially(float *out, int pdmax,
		struct plambda_program *p,
		float **val, int *w, int *h, int *pd)
{
	for (int j = 0; j < *h; j++)
	for (int i = 0; i < *w; i++)
	{
		float result[pdmax];
		int r = run_program_vectorially_at(result, p,val, w,h,pd, i,j);
		assert(r == pdmax);
		if (r != pdmax) fail("r != pdmax");
		for (int l = 0; l < r; l++)
			setsample_0(out, *w, *h, pdmax, i, j, l, result[l]);
	}
	return pdmax;
}

// mains {{{1

static void add_hidden_variables(char *out, int maxplen, int newvars, char *in)
{
	int pos = 0;
	for (int i = 0; i < newvars; i++)
		pos += snprintf(out + pos, maxplen - pos, "hidden_%02d ", i);
	snprintf(out + pos, maxplen - pos, "%s", in);
	//fprintf(stderr, "HIVA: %s\n", out);
}

SMART_PARAMETER_SILENT(SRAND,0)

static int main_calc(int c, char **v)
{
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s v1 v2 ... \"plambda\"\n", *v);
		//                          0 1  2        c-1
		return EXIT_FAILURE;
	}

	struct plambda_program p[1];
	plambda_compile_program(p, v[c-1]);

	int n = c - 2, pd[n], pdmax = PLAMBDA_MAX_PIXELDIM;
	if (n > 0 && p->var->n == 0) {
		int maxplen = n*20 + strlen(v[c-1]) + 100;
		char newprogram[maxplen];
		add_hidden_variables(newprogram, maxplen, n, v[c-1]);
		plambda_compile_program(p, newprogram);
	}
	if (n != p->var->n)
		fail("the program expects %d variables but %d vectors "
					"were given", p->var->n, n);

	float *x[n];
	FORI(n) x[i] = alloc_parse_floats(pdmax, v[i+1], pd+i);

	FORI(n) if (!strstr(p->var->t[i], "hidden"))
		fprintf(stderr, "calculator correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	xsrand(100+SRAND());

	float out[pdmax];
	int od = run_program_vectorially_at(out, p, x, NULL, NULL, pd, 0, 0);

	char *fmt = getenv("PLAMBDA_FFMT");
	if (!fmt) fmt = "%.15lf";
	for (int i = 0; i < od; i++)
	{
		printf(fmt, out[i]);
		putchar(i==(od-1)?'\n':' ');
	}

	collection_of_varnames_end(p->var);
	FORI(n) free(x[i]);


	return EXIT_SUCCESS;
}

// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
static char *pick_option(int *c, char ***v, char *o, char *d)
{
	int argc = *c;
	char **argv = *v;
	for (int i = 0; i < argc - 1; i++)
		if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
		{
			char *r = argv[i+1];
			*c -= 2;
			for (int j = i; j < argc - 1; j++)
				(*v)[j] = (*v)[j+2];
			return r;
		}
	return d;
}

#include "iio.h"
static int main_images(int c, char **v)
{
	//fprintf(stderr, "main images c = %d\n", c);
	//for (int i = 0; i < c; i++)
	//	fprintf(stderr, "main images argv[%d] = %s\n", i, v[i]);
	if (c < 2) {
		fprintf(stderr, "usage:\n\t%s in1 in2 ... \"plambda\"\n", *v);
		//                          0 1   2         c-1
		return EXIT_FAILURE;
	}
	char *filename_out = pick_option(&c, &v, "o", "-");

	struct plambda_program p[1];

	plambda_compile_program(p, v[c-1]);

	int n = c - 2;
	//fprintf(stderr, "n = %d\n", n);
	if (n > 0 && p->var->n == 0) {
		//fprintf(stderr, "will add hidden variables! n=%d, vn=%d\n", n, p->var->n);
		int maxplen = n*10 + strlen(v[c-1]) + 100;
		char newprogram[maxplen];
		add_hidden_variables(newprogram, maxplen, n, v[c-1]);
		plambda_compile_program(p, newprogram);
	}
	if (n != p->var->n && !(n == 1 && p->var->n == 0))
		fail("the program expects %d variables but %d images "
					"were given", p->var->n, n);
	int w[n], h[n], pd[n];
	float *x[n];
	FORI(n) x[i] = iio_read_image_float_vec(v[i+1], w + i, h + i, pd + i);
	//FORI(n-1)
	//	if (w[0] != w[i+1] || h[0] != h[i+1])// || pd[0] != pd[i+1])
	//		fail("input images size mismatch");

	if (n>1) FORI(n) if (!strstr(p->var->t[i], "hidden"))
		fprintf(stderr, "plambda correspondence \"%s\" = \"%s\"\n",
				p->var->t[i], v[i+1]);

	xsrand(100+SRAND());

	//print_compiled_program(p);
	int pdreal = eval_dim(p, x, pd);

	float *out = xmalloc(*w * *h * pdreal * sizeof*out);
	int opd = run_program_vectorially(out, pdreal, p, x, w, h, pd);
	assert(opd == pdreal);

	iio_write_image_float_vec(filename_out, out, *w, *h, opd);

	FORI(n) free(x[i]);
	free(out);
	collection_of_varnames_end(p->var);

	return EXIT_SUCCESS;
}

static int print_version(void)
{
	printf("plambda 1.0\n\n"
	"Written by Enric Meinhardt-Llopis\n");
	return 0;
}

static int print_help(char *v, int verbosity)
{
	printf(
"Plambda evaluates an expression with images as variables.\n"
"\n"
"The expression is written in reverse polish notation using common\n"
"operators and functions from `math.h'.  The variables appearing on the\n"
"expression are assigned to each input image in alphabetical order.\n"
//"The resulting image is printed to standard output.  The expression\n"
//"should be written in reverse polish notation using common operators\n"
//"and functions from `math.h'.  The variables appearing on the\n"
//"expression are assigned to each input image in alphabetical order.\n"
"%s"
"\n"
"Usage: %s a.png b.png c.png ... \"EXPRESSION\" > output\n"
"   or: %s a.png b.png c.png ... \"EXPRESSION\" -o output.png\n"
"   or: %s -c num1 num2 num3  ... \"EXPRESSION\"\n"
"\n"
"Options:\n"
" -o file\tsave output to named file\n"
" -c\t\tact as a symbolic calculator\n"
" -h\t\tdisplay short help message\n"
" --help\t\tdisplay longer help message\n"
//" --version\tdisplay version\n"
//" --man\tdisplay manpage\n"
"\n"
"Examples:\n"
" plambda a.tiff b.tiff \"x y +\" > sum.tiff\tCompute the sum of two images.\n"
" plambda -c \"1 atan 4 *\"\t\t\tPrint pi\n"
"%s"
"\n"
"Report bugs to <enric.meinhardt@cmla.ens-cachan.fr>.\n",
verbosity>0?
"\n"
"EXPRESSIONS:\n\n"
"A \"plambda\" expression is a sequence of tokens.\nTokens may be constants,\n"
"variables, or operators.  Constants and variables get their value\n"
"computed and pushed to the stack.  Operators pop values from the stack,\n"
"apply a function to them, and push back the results.\n"
"\n"
"CONSTANTS: numeric constants written in scientific notation, and \"pi\"\n"
"\n"
"OPERATORS: +, -, *, ^, /, <, >, ==, and all the functions from math.h\n"
"\n"
"LOGIC OPS: if, and, or, not\n"
"\n"
"VARIABLES: anything not recognized as a constant or operator.  There\n"
"must be as many variables as input images, and they are assigned to\n"
"images in alphabetical order.  If there are no variables, the input\n"
"images are pushed to the stack.\n"
"\n"
"All operators (unary, binary and ternary) are vectorizable.  Thus, you can\n"
"add a scalar to a vector, divide two vectors of the same size, and so on.\n"
"The semantics of each operation follows the principle of least surprise.\n"
"\n"
"Some \"sugar\" is added to the language:\n"
"\n"
"Predefined variables (always preceeded by a colon):\n"
" :i\thorizontal coordinate of the pixel\n"
" :j\tvertical coordinate of the pixel\n"
" :w\twidth of the image\n"
" :h\theigth of the image\n"
" :n\tnumber of pixels in the image\n"
" :x\trelative horizontal coordinate of the pixel\n"
" :y\trelative horizontal coordinate of the pixel\n"
" :r\trelative distance to the center of the image\n"
" :t\trelative angle from the center of the image\n"
" :I\thorizontal coordinate of the pixel (centered)\n"
" :J\tvertical coordinate of the pixel (centered)\n"
" :P\thorizontal coordinate of the pixel (phased)\n"
" :Q\tvertical coordinate of the pixel (phased)\n"
" :R\tcentered distance to the center\n"
" :L\tminus squared centered distance to the center\n"
" :W\twidth of the image divided by 2*pi\n"
" :H\theight of the image divided by 2*pi\n"
"\n"
"Variable modifiers acting on regular variables:\n"
" x\t\tvalue of pixel (i,j)\n"
" x(0,0)\t\tvalue of pixel (i,j)\n"
" x(1,0)\t\tvalue of pixel (i+1,j)\n"
" x(0,-1)\tvalue of pixel (i,j-1)\n"
" x[0]\t\tvalue of first component of pixel (i,j)\n"
" x[1]\t\tvalue of second component of pixel (i,j)\n"
" x(1,2)[3]\tvalue of fourth component of pixel (i+1,j+2)\n"
"Comma modifiers (pre-defined local operators):\n"
" a,x\tx-derivative of the image a\n"
" a,y\ty-derivative\n"
" a,xx\tsecond x-derivative\n"
" a,yy\tsecond y-derivative\n"
" a,xy\tcrossed second derivative\n"
" a,l\tLaplacian\n"
" a,g\tgradient\n"
" a,n\tgradient norm\n"
" a,d\tdivergence\n"
" a,S\tshadow operator\n"
" a,xf\tx-derivative, forward differences\n"
" a,xb\tx-derivative, backward differences\n"
" a,xc\tx-derivative, centered differences\n"
" a,xs\tx-derivative, sobel\n"
" a,xp\tx-derivative, prewitt\n"
" etc\n"
"\n"
"Stack operators (allow direct manipulation of the stack):\n"
" del\tremove the value at the top of the stack (ATTTOS)\n"
" dup\tduplicate the value ATTTOS\n"
" rot\tswap the two values ATTTOS\n"
" split\tsplit the vector ATTTOS into scalar components\n"
" join\tjoin the components of two vectors ATTOTS\n"
" join3\tjoin the components of three vectors ATTOTS\n"
" njoin\tjoin the components of n vectors\n"
" halve\tsplit an even-sized vector ATTOTS into two equal-sized parts\n"
//" interleave\tinterleave\n"
//" deinterleave\tdeinterleave\n"
//" nsplit\tnsplit\n"
"\n"
//"Magic variable modifiers:\n"
"Magic variable modifiers (global data associated to each input image):\n"
" x%i\tvalue of the smallest sample of image x\n"
" x%a\tvalue of the largest sample\n"
" x%v\taverage sample value\n"
" x%m\tmedian sample value\n"
" x%s\tsum of all samples\n"
" x%I\tvalue of the smallest pixel (in euclidean norm)\n"
" x%A\tvalue of the largest pixel\n"
" x%V\taverage pixel value\n"
" x%S\tsum of all pixels\n"
" x%Y\tcomponent-wise minimum of all pixels\n"
" x%E\tcomponent-wise maximum of all pixels\n"
" x%qn\tnth sample percentile\n"
" x%On\tcomponent-wise nth percentile\n"
" x%Wn\tcomponent-wise nth millionth part\n"
" x%0n\tcomponent-wise nth order statistic\n"
" x%9n\tcomponent-wise nth order statistic (from the right)\n"
//" x[2]%i\tminimum value of the blue channel\n"
//" \n"
//" x%M\tmedian pixel value\n"
"\n"
"Random numbers (seeded by the SRAND environment variable):\n"
" randu\tpush a random number with distribution Uniform(0,1)\n"
" randn\tpush a random number with distribution Normal(0,1)\n"
" randc\tpush a random number with distribution Cauchy(0,1)\n"
" randl\tpush a random number with distribution Laplace(0,1)\n"
" rande\tpush a random number with distribution Exponential(1)\n"
" randp\tpush a random number with distribution Pareto(1)\n"
" rand\tpush a random integer returned from rand(3)\n"
"\n"
"Vectorial operations (acting over vectors of a certain length):\n"
" topolar\tconvert a 2-vector from cartesian to polar\n"
" frompolar\tconvert a 2-vector from polar to cartesian\n"
" hsv2rgb\tconvert a 3-vector from HSV to RGB\n"
" rgb2hsv\tconvert a 3-vector from RGB to HSV\n"
" xyz2rgb\tconvert a 3-vector from XYZ to RGB\n"
" rgb2xyz\tconvert a 3-vector from RGB to XYZ\n"
" cprod\t\tmultiply two 2-vectrs as complex numbers\n"
" mprod\t\tmultiply two 2-vectrs as matrices (4-vector = 2x2 matrix, etc)\n"
" vprod\t\tvector product of two 3-vectors\n"
" sprod\t\tscalar product of two n-vectors\n"
" mdet\t\tdeterminant of a n-matrix (a n*n-vector)\n"
" mtrans\t\ttranspose of a matrix\n"
" mtrace\t\ttrace of a matrix\n"
" minv\t\tinverse of a matrix\n"
"\n"
"Registers (numbered from 1 to 9):\n"
" >7\tcopy to register 7\n"
" <3\tcopy from register 3\n"
"\n"
//"Environment:\n"
//" SRAND\tseed of the random number generator (default=1)\n"
//" CAFMT\tformat of the number printed by the calculator (default=%.15lf)\n"
	:
	"See the manual page for details on the syntax for expressions.\n"
	,
	v, v, v,
	verbosity < 1 ? "" :
	" plambda -c \"355 113 /\"\t\t\t\tPrint an approximation of pi\n"
		);
	return 0;
}

static int do_man(void)
{
#ifdef __OpenBSD__
#define MANPIPE "|mandoc -a"
#else
#define MANPIPE "|man -l -"
#endif
	return system("help2man -N -S imscript -n \"evaluate an expression "
				"with images as variables\" plambda" MANPIPE);
}

int main_plambda(int c, char **v)
{
	if (c == 1) return print_help(*v, 0);
	if (c == 2 && 0 == strcmp(v[1], "-h")) return print_help(*v,0);
	if (c == 2 && 0 == strcmp(v[1], "--help")) return print_help(*v,1);
	if (c == 2 && 0 == strcmp(v[1], "--version")) return print_version();
	if (c == 2 && 0 == strcmp(v[1], "--man")) return do_man();

	int (*f)(int, char**) = (**v=='c' || c==2) ?  main_calc : main_images;
	if (f == main_images && c > 2 && 0 == strcmp(v[1], "-c"))
	{
		for (int i = 1; i <= c; i++)
			v[i] = v[i+1];
		f = main_calc;
		c -= 1;
	}
	return f(c,v);
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_plambda(c, v); }
#endif

// vim:set foldmethod=marker:
