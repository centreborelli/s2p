
#ifndef _BICUBIC_C
#define _BICUBIC_C


typedef float (*getpixel_operator)(float*,int,int,int,int);

// extrapolate by 0
inline
static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i+j*w];
}

// extrapolate by nearest value
inline
static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}


static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));

	float y = 3.0*(v[1] - v[2]) + v[3] - v[0];
	y = x*y + 2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3];
	y = x*y + v[2] - v[0];
	y = 0.5*x*y + v[1];
	return y;
}

static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

float bicubic_interpolation_gray(float *img, int w, int h, float x, float y)
{
	x -= 1;
	y -= 1;

	getpixel_operator p = getpixel_0;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, ix + i, iy + j);
	float r = bicubic_interpolation_cell(c, x - ix, y - iy);
	return r;
}

#endif//_BICUBIC_C
