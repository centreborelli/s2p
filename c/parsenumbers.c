#ifndef _PARSENUMBERS_C
#define _PARSENUMBERS_C

#include <string.h>
#include "xmalloc.c"

// utility function: parse floats from a text file
// returns a pointer to a malloc'ed array of the parsed floats
// fills *no with the number of floats
static float *read_ascii_floats(FILE *f, int *no)
{
	int r, n = 0, nt = 0;
	float *t = NULL;
	while(1) {
		if (n >= nt)
		{
			nt = 2 * (nt + 1);
			t = xrealloc(t, nt * sizeof * t);
		}
		r = fscanf(f, "%f", t+n);
		if (r != 1)
			break;
		n += 1;
	}
	*no = n;
	return t;
}

// utility function: parse doubles from a text file
// returns a pointer to a malloc'ed array of the parsed doubles
// fills *no with the number of doubles
static double *read_ascii_doubles(FILE *f, int *no)
{
	int r, n = 0, nt = 0;
	double *t = NULL;
	while(1) {
		if (n >= nt)
		{
			nt = 2 * (nt + 1);
			t = xrealloc(t, nt * sizeof * t);
		}
		r = fscanf(f, "%lf", t+n);
		if (r != 1)
			break;
		n += 1;
	}
	*no = n;
	return t;
}

//// utility function: parse a matrix of doubles from a text file
//// returns a pointer to a malloc'ed array of the parsed doubles
//// fills *nc with the number of cols and *nr with the number of rows
//static double *read_ascii_matrix_doubles(FILE *f, int *nc, int *nl)
//{
//	int r, n = 0, nt = 0;
//	double *t = NULL;
//	while(1) {
//		if (n >= nt)
//		{
//			nt = 2 * (nt + 1);
//			t = xrealloc(t, nt * sizeof * t);
//		}
//		r = fscanf(f, "%lf", t+n);
//		if (r != 1)
//			break;
//		n += 1;
//	}
//	*no = n;
//	return t;
//}

static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

static int parse_floats(float *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}

static float *alloc_parse_floats(int nmax, const char *ss, int *n)
{
	// add a space to s, so that gabriele's expression works as intended
	int ns = strlen(ss);
	char t[2+ns], *s = t;
	s[0] = ' ';
	for (int i = 0; i <= ns; i++)
		s[i+1] = ss[i];

	float *r = xmalloc(nmax * sizeof*r);
	int i = 0, w;
	//while (i < nmax && 1 == sscanf(s, "%g %n", r + i, &w)) {
	while (i < nmax && 1 == sscanf(s, "%*[][ \n\t,:;]%g%n", r + i, &w)) {
		//fprintf(stderr, "read %g\n", r[i]);
		i += 1;
		s += w;
	}
	//fprintf(stderr, "finished, got i=%d\n", i);
	*n = i;
	return r;
}

#endif//_PARSENUMBERS_C
