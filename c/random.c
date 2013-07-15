#ifndef _RANDOM_C
#define _RANDOM_C




#ifdef RNG_SPECIFIC_WELL1024
#include "well1024.c"
#else
#include <stdlib.h>
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif


static void xsrand(unsigned int seed)
{
#ifdef RNG_SPECIFIC_WELL1024
	well1024_seed(seed);
#else
	srand(seed);
#endif
}

static int xrand(void)
{
#ifdef RNG_SPECIFIC_WELL1024
	int r;
# ifdef _OPENMP
# pragma omp critical
# endif
	r = RAND_MAX * well1024();
	return r;
#else
	return rand();
#endif
}

static double random_raw(void)
{
	return xrand();
}

static double random_uniform(void)
{
#ifdef RNG_SPECIFIC_WELL1024
	double r;
# ifdef _OPENMP
# pragma omp critical
# endif
	r = well1024();
	return r;
#else
	return rand()/(1.0+RAND_MAX);
#endif
}

static double random_normal(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	//double y2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	return y;
}

int randombounds(int a, int b)
{
	if (b < a)
		return randombounds(b, a);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

static double random_laplace(void)
{
	double x = random_uniform();
	double y = random_uniform();
	double r = log(x/y);
	return isfinite(r)?r:0;
}

static double random_cauchy(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y_1 = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	double y_2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	double r = y_1/y_2;
	return isfinite(r)?r:0;
}

static double random_exponential(void)
{
	//double u = random_uniform();
	//double r = -log(1-u);
	//return r;
	return fabs(random_laplace());
}

static double random_pareto(void)
{
	return exp(random_exponential());
}

#endif//_RANDOM_C
