#ifndef _GETPIXEL_C
#define _GETPIXEL_C

typedef float (*getsample_operator)(float*,int,int,int,int,int,int);
//typedef void (*setsample_operator)(float*,int,int,int,int,int,int,float);

// extrapolate by 0
inline
static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return x[(i+j*w)*pd + l];
}

// extrapolate by nearest value
inline
static float getsample_1(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l >= pd) l = pd-1;
	return x[(i+j*w)*pd + l];
}

#ifdef NAN
// extrapolate by nan
inline
static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return NAN;
	return x[(i+j*w)*pd + l];
}
#endif//NAN

// force a segfault if extrapolation is required
inline
static float getsample_error(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return *(volatile int*)0;
	return x[(i+j*w)*pd + l];
}

// abort the program extrapolation is required
inline
static float getsample_abort(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		abort();
	return x[(i+j*w)*pd + l];
}

// exit the program extrapolation is required
inline
static float getsample_exit(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		exit(42);
	return x[(i+j*w)*pd + l];
}

// like n%p, but works for all numbers
static int good_modulus(int n, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(n, -p);

	int r = n % p;
	r = r < 0 ? r + p : r;

//	assert(r >= 0);
//	assert(r < p);
	return r;
}


static int gmod(int x, int m)
{
	int r = x % m;
	return r < 0 ? r + m : r;
}

static int positive_reflex(int n, int p)
{
	int r = good_modulus(n, 2*p);
	if (r == p) r -= 1;
	if (r > p)
		r = 2*p - r;
	if (n < 0 && p > 1) r += 1;
	//assert(r >= 0);
	//assert(r < p);
	return r;
}

// extrapolate by reflection
inline
static float getsample_2(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = positive_reflex(i, w);
	j = positive_reflex(j, h);
	return getsample_abort(x, w, h, pd, i, j, l);
}

// extrapolate by periodicity
inline
static float getsample_per(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = gmod(i, w);
	j = gmod(j, h);
	return getsample_abort(x, w, h, pd, i, j, l);
}

// extrapolate by constant (set by calling it with zero sizes)
inline static
float getsample_constant(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return value;
	return x[(i+j*w)*pd + l];
}

// test for inclusion of stdlib.h and string.h
#if defined(EXIT_SUCCESS)
#if defined(_STRING_H) || defined(_STRING_H_)
static getsample_operator get_sample_operator(getsample_operator o)
{
	char *option = getenv("GETPIXEL"), *endptr;
	if (!option) return o;
#ifdef NAN
	if (0 == strcmp(option, "nan"      ))  return getsample_nan;
#endif//NAN
	if (0 == strcmp(option, "segfault" ))  return getsample_error;
	if (0 == strcmp(option, "error"    ))  return getsample_error;
	if (0 == strcmp(option, "abort"    ))  return getsample_abort;
	if (0 == strcmp(option, "exit"    ))   return getsample_exit;
	if (0 == strcmp(option, "periodic" ))  return getsample_per;
	if (0 == strcmp(option, "constant" ))  return getsample_0;
	if (0 == strcmp(option, "nearest"   )) return getsample_1;
	if (0 == strcmp(option, "reflex"   ))  return getsample_2;
	if (0 == strcmp(option, "symmetric"))  return getsample_2;
	float value = strtof(option, &endptr);
	if (endptr != option) {
		getsample_constant(&value, 0, 0, 0, 0, 0, 0);
		return getsample_constant;
	}
	return getsample_0;
}
#endif//_STRING_H
#endif//EXIT_SUCCESS


inline
static void setsample_0(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
	x[(i+j*w)*pd + l] = v;
}

inline
static void setsample_segf(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		*(volatile int*)0 = 0;
	x[(i+j*w)*pd + l] = v;
}

typedef float (*getpixel_operator)(float*,int,int,int,int);

static float getpixel_0(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	return x[i + j*w];
}

static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[i+j*w];
}

//
//static void setpixel(float *x, int w, int h, int i, int j, float v)
//{
//	if (i < 0 || i >= w || j < 0 || j >= h)
//		return;
//	x[i + j*w] = v;
//}

#endif//_GETPIXEL_C
