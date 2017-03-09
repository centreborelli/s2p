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

inline
static double *read_ascii_doubles_fn(const char *fname, int *no)
{
	FILE *f = fopen(fname, "r");
	if (!f) { *no = 0; return NULL; }
	double *r = read_ascii_doubles(f, no);
	fclose(f);
	return r;
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

static float *alloc_parse_floats2(int nmax, const char *s, int *n)
{
	if (!nmax) nmax = strlen(s);
	float *r = xmalloc(nmax*sizeof*r);
	s += strspn(s, "[] \n\t,:;()");
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%g%*[][ \n\t,:;]%n", r + i, &w))
	{
		i += 1;
		s += w;
	}
	*n = i;
	return r;
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

static double *alloc_parse_doubles(int nmax, const char *ss, int *n)
{
	// add a space to s, so that gabriele's expression works as intended
	int ns = strlen(ss);
	char t[2+ns], *s = t;
	s[0] = ' ';
	for (int i = 0; i <= ns; i++)
		s[i+1] = ss[i];

	double *r = xmalloc(nmax * sizeof*r);
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%*[][ \n\t,:;]%lf%n", r + i, &w)) {
		i += 1;
		s += w;
	}
	*n = i;
	return r;
}

// Obtain n numbers from string.
// The string contains a list of ascii numbers
// or the name of a file containing a list of ascii numbers.
// In case of failure, unread numbers are set to zero.
// Returns the number of read numbers.
inline
static int read_n_doubles_from_string(double *out, char *string, int n)
{
	for (int i = 0; i < n; i++)
		out[i] = 0;

	int no;
	double *buf = NULL;
	FILE *f = fopen(string, "r");
	if (f) {
		buf = read_ascii_doubles(f, &no);
		fclose(f);
	} else {
		buf = alloc_parse_doubles(n, string, &no);
	}

	if (no > n) no = n;
	for (int i = 0; i < no; i++)
		out[i] = buf[i];
	free(buf);
	return no;
}

#endif//_PARSENUMBERS_C
