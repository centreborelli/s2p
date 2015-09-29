
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

#include "string.h"
#include "img.h"



// used for census transform
static void pack_bits_into_bytes(unsigned char *out, int *bits, int nbits)
{
	int nbytes = ceil(nbits/8.0);
	for (int i = 0; i < nbytes; i++)
	{
		out[i] = 0;
		for (int j = 0; j < 8; j++)
			out[i] = out[i] * 2 + bits[8*i+j];
	}
}

// used for census transform
// this function doesn't use the color plane convention
static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return NAN;
	return x[(i+j*w)*pd + l];
}


// used for census transform
static void census_at(unsigned char *out, float *x, int w, int h, int pd,
		int winradius, int p, int q)
{
	int side = 2*winradius + 1;
	int nbits = pd * (side * side - 1);
	int bits[nbits];
	int cx = 0;
	for (int l = 0; l < pd; l++)
	for (int j = -winradius; j <= winradius; j++)
	for (int i = -winradius; i <= winradius; i++)
	{
		float a = getsample_nan(x, w, h, pd, p    , q    , l);
		float b = getsample_nan(x, w, h, pd, p + i, q + j, l);
		if (i || j)
			bits[cx++] = a < b;
	}
	assert(cx == nbits);

	pack_bits_into_bytes(out, bits, nbits);
}

// census transform
static void color_census_transform(unsigned char *y, int opd,
		float *x, int w, int h, int pd, int winradius)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		census_at(y + opd * (w * j + i), x, w, h, pd, winradius, i, j);
}

static void pack_bytes_into_floats(float *y, uint8_t *x, int nf, int nb)
{
	int sf = sizeof(*y);
	assert(nb <= nf * sf);

	memcpy(y, x, nb);
}

static float *malloc_floating_census_joint_channels(float *x, int w, int h, int pd,
		int winradius, int *out_pd)
{
	int side = 2 * winradius + 1;
	int nbits = pd * (side * side - 1);
	assert(0 == nbits % 8);
	int nbytes = nbits / 8;
	int nfloats = ceil(nbytes / (float)sizeof(float));

	uint8_t *y = (uint8_t*) malloc(w * h * nbytes);
	color_census_transform(y, nbytes, x, w, h, pd, winradius);

	float *fy = (float*) malloc(w * h * nfloats * sizeof*fy);
	for (int i = 0; i < w*h; i++)
	{
		float *to = fy + nfloats * i;
		uint8_t *from = y + nbytes * i;
		pack_bytes_into_floats(to, from, nfloats, nbytes);
	}
	free(y);

	*out_pd = nfloats;
	return fy;
}





float compute_census_distance_array(uint8_t *a, uint8_t *b, int n)
{
	int count_bits[256] = {
#               define B2(n) n,     n+1,     n+1,     n+2
#               define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#               define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
                B6(0), B6(1), B6(1), B6(2)
	}, r = 0;
	for (int i = 0; i < n; i++)
		r += count_bits[a[i] ^ b[i]];
	return r;
}


float distance_patch_census_packed_in_floats(float *u1, float *u2, int n)
{
	return compute_census_distance_array((uint8_t*)u1, (uint8_t*)u2, n * sizeof(float));
}



// generates the census transformed image from an image bx 
struct Img census_transform(struct Img &bx, int winradius)
{
   int w  = bx.nx;
   int h  = bx.ny;
   int pd = bx.nch;
	float *x = (float*) malloc(w * h * pd * sizeof*x);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < pd; l++)
		x[(w*j + i)*pd + l] = bx[(w*h)*l + j*w + i];

	int opd;
	float *y = malloc_floating_census_joint_channels(x, w, h, pd, winradius, &opd);

   struct Img by(w,h,opd);

	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < opd; l++)
		by[(w*h)*l + j*w + i] = y[(w*j + i)*opd + l];

	free(x);
	free(y);

	return by;
}


