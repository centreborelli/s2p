#ifndef _RANDOM_C
#define _RANDOM_C

#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>


// Why write a simple random number generator here?
// It is the easiest way to fullfill the following requirements
// 	1. very fast, negligible memory usage
// 	2. generates the same sequence on any architecture
// 	3. do not have annoying warnings (e.g. unsafe rand() in openbsd).
// Strong non-requirements:
// 	1. cryptographicaly secure
// 	2. pass fancy statistical tests

static uint64_t lcg_knuth_seed = 0;

void lcg_knuth_srand(uint32_t x)
{
	lcg_knuth_seed = x;
}

uint32_t lcg_knuth_rand(void)
{
	lcg_knuth_seed *= 6364136223846793005;
	lcg_knuth_seed += 1442695040888963407;
	return lcg_knuth_seed >> 32;
}


static void xsrand(unsigned int seed)
{
	lcg_knuth_srand(seed);
}

static int xrand(void)
{
	return lcg_knuth_rand();
}

static double random_raw(void)
{
	return xrand();
}

static double random_uniform(void)
{
	return lcg_knuth_rand()/(0.0+UINT_MAX);
}

static double random_ramp(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	//double y = sqrt(random_uniform());
	double y = fmax(x1, x2);
	return y;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

static double random_normal(void)
{
	double x1 = random_uniform();
	double x2 = random_uniform();
	double y = sqrt(-2*log(x1)) * cos(2*M_PI*x2);
	//double y2 = sqrt(-2*log(x1)) * sin(2*M_PI*x2);
	return y;
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

// This function samples a stable variable of parameters (alpha,beta)
//
// The algorithm is copied verbatim from formulas (2.3) and (2.4) in
//     Chambers, J. M.; Mallows, C. L.; Stuck, B. W. (1976).
//     "A Method for Simulating Stable Random Variables".
//     Journal of the American Statistical Association 71 (354): 340–344.
//     doi:10.1080/01621459.1976.10480344.
//
// Observation: the algorithm is numerically imprecise when alpha approaches 1.
// TODO: implement appropriate rearrangements as suggested in the article.
static double random_stable(double alpha, double beta)
{
	double U = (random_uniform() - 0.5) * M_PI;
	double W = random_exponential();
	double z = -beta * tan(M_PI * alpha / 2);
	double x = alpha == 1 ? M_PI/2 : atan(-z) / alpha;
	//fprintf(stderr, "U=%g W=%g z=%g x=%g\t\t", U, W, z, x);
	double r = NAN;
	if (alpha == 1) {
		double a = (M_PI/2 + beta * U) * tan(U);
		double b = log(((M_PI/2) * W * cos(U)) / ((M_PI/2) + beta * U));
		//fprintf(stderr, "a=%g b=%g\n", a, b);
		r = (a - beta * b) / x;
	} else {
		double a = pow(1 + z * z, 1 / (2*alpha));
		double b = sin(alpha * (U + x)) / pow(cos(U), 1/alpha);
		double c = pow(cos(U - alpha*(U + x)) / W, (1 - alpha) / alpha);
		//fprintf(stderr, "a=%g b=%g c=%g\n", a, b, c);
		r = a * b * c;
	}
	//fprintf(stderr, "s(%g,%g) = %g\n", alpha, beta, r);
	return r;
}

#endif//_RANDOM_C
