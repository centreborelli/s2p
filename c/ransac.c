#include <assert.h>
#include <stdbool.h>
#include <math.h>

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"

#include "cmphomod.c"

// generic function
// evaluate the error of a datapoint according to a model
// (implementing this function is necessary for each ransac case)
typedef float (ransac_error_evaluation_function)(
		float *model,
		float *datapoint,
		void *usr
		);


// generic function
// compute the model defined from a few data points
// (shall return 0 if no model could be computed)
// (implementing this function is necessary for each ransac case)
typedef int (ransac_model_generating_function)(
		float *out_model,  // parameters of the computed model
		float *data,       // data points
		void *usr
		);

// generic function
// tell whether a given model is good enough (e.g., not severely distorted)
// (this function is optional, and only serves as an optimization hint)
typedef bool (ransac_model_accepting_function)(
		float *model,
		void  *usr);


// API function: evaluate a given model over the data, and fill a mask with the
// inliers (according to the given allowed error).  This function returns the
// number of inliers.
int ransac_trial(
		// output
		bool *out_mask,    // array mask identifying the inliers

		// input data
		float *data,       // array of input data
		float *model,      // parameters of the model
		float max_error,   // maximum allowed error

		// input context
		int datadim,       // dimension of each data point
		int n,             // number of data points
		ransac_error_evaluation_function *mev,

		// decoration
		void *usr
		)
{
	int cx = 0;
	for (int i = 0; i < n; i++)
	{
		float *datai = data + i*datadim;
		float e = mev(model, datai, usr);
		if (!(e >= 0)) fprintf(stderr, "WARNING e = %g\n", e);
		assert(e >= 0);
		out_mask[i] = e < max_error;
		if (out_mask[i])
			cx += 1;
	}
	return cx;
}

// utility function: return a random number in the interval [a, b)
static int random_index(int a, int b)
{
	int r = a + rand()%(b - a);
	assert(r >= a);
	assert(r < b);
	return r;
}

// comparison function for the qsort call below
static int compare_ints(const void *aa, const void *bb)
{
	const int *a = (const int *)aa;
	const int *b = (const int *)bb;
	return (*a > *b) - (*a < *b);
}

// check whether a vector of n ints has different entries
// (sorts the vector inplace)
static bool are_different(int *t, int n)
{
	qsort(t, n, sizeof*t, compare_ints);
	for (int i = 1; i < n; i++)
		if (t[i-1] == t[i])
			return false;
	return true;
}

static int randombounds(int a, int b)
{
	if (b < a)
		fail("the interval [%d, %d] is empty!", a, b);
	if (b == a)
		return b;
	return a + rand()%(b - a + 1);
}

static void swap(void *a, void *b, size_t s)
{
#if 0
#error "memcpy is way slower!"
	char t[s];
	memcpy(t, a, s);
	memcpy(a, b, s);
	memcpy(b, t, s);
#else
	char *x = a;
	char *y = b;
	for (unsigned int i = 0; i < s; i++, x++, y++)
	{
		char t = *x;
		*x = *y;
		*y = t;
	}
#endif
}

// fisher-yates
void shuffle(void *t, int n, size_t s)
{
	char *c = t;

	for (int i = 0; i < n-1; i++)
		swap(c + s*i, c + s*randombounds(i, n-1), s);
}

static void fill_random_shuffle(int *idx, int n, int a, int b)
{
	int m = b - a;
	assert(n > 0);
	assert(m > 0);
	assert(m >= n);
	int t[m];
	for (int i = a; i < b; i++)
		t[i] = i;
	shuffle(t, m, sizeof*t);
	for (int i = 0; i < n; i++)
		idx[i] = t[i];
}

// generate a set of n different ints between a and b
static void fill_random_indices(int *idx, int n, int a, int b)
{
	if (b-a==n) {for(int i=0;i<n;i++)idx[i]=a+i;}
	if (5*n > (b-a)) {fill_random_shuffle(idx, n, a, b);return;}
	// TODO fisher yates shuffle and traverse it by blocks of length nfit
	int safecount = 0;
	do {
		for (int i = 0; i < n; i++)
			idx[i] = random_index(a, b);
		safecount += 1;
	} while (safecount < 100 && !are_different(idx, n));
	if (safecount == 100)
		fail("could not generate any model");
	//fprintf(stderr, "fri");
	//for (int i = 0; i < n; i++)
	//	fprintf(stderr, "\t%d", idx[i]);
	//fprintf(stderr, "\n");
}

#define MAX_MODELS 10


// RANSAC
//
// Given a list of data points, find the parameters of a model that fits to
// those points.  Several models are tried, and the model with the highest
// number of inliers is kept.
//
// A basic idea of this kind of ransac is that a maximum allowed error is fixed
// by hand, and then the inliers of a model are defined as the data points
// which fit the model up to the allowed error.  The RANSAC algorithm randomly
// tries several models and keeps the one with the largest number of inliers.
int ransac(
		// output
		//int *out_ninliers, // number of inliers
		bool *out_mask,    // array mask identifying the inliers
		float *out_model,  // model parameters

		// input data
		float *data,       // array of input data

		// input context
		int datadim,       // dimension of each data point
		int n,             // number of data points
		int modeldim,      // number of model parameters
		ransac_error_evaluation_function *mev,
		ransac_model_generating_function *mgen,
		int nfit,          // data points needed to produce a model

		// input parameters
		int ntrials,       // number of models to try
		int min_inliers,   // minimum allowed number of inliers
		float max_error,   // maximum allowed error

		// decoration
		ransac_model_accepting_function *macc,
		void *usr
		)
{
	fprintf(stderr, "running RANSAC over %d datapoints of dimension %d\n",
			n, datadim);
	fprintf(stderr, "will try to find a model of size %d from %d points\n",
		       	modeldim, nfit);
	fprintf(stderr, "we will make %d trials and keep the best with e<%g\n",
			ntrials, max_error);
	fprintf(stderr, "a model must have more than %d inliers\n",
			min_inliers);

	int best_ninliers = 0;
	float best_model[modeldim];
	bool *best_mask = xmalloc(n * sizeof*best_mask);
	bool *tmp_mask = xmalloc(n * sizeof*best_mask);

	for (int i = 0; i < ntrials; i++)
	{
		int indices[nfit];
		fill_random_indices(indices, nfit, 0, n);

		float x[nfit*datadim];
		for (int j = 0; j < nfit; j++)
		for (int k = 0; k < datadim; k++)
			x[datadim*j + k] = data[datadim*indices[j] + k];

		float model[modeldim*MAX_MODELS];
		int nm = mgen(model, x, usr);
		if (!nm)
			continue;
		if (macc && !macc(model, usr))
			continue;

		// generally, nm=1
		for (int j = 0; j < nm; j++)
		{
			float *modelj = model + j*modeldim;
			int n_inliers = ransac_trial(tmp_mask, data, modelj,
					max_error, datadim, n, mev, usr);

			if (n_inliers > best_ninliers)
			{
				best_ninliers = n_inliers;
				for(int k = 0; k < modeldim; k++)
					best_model[k] = modelj[k];
				for(int k = 0; k < n; k++)
					best_mask[k] = tmp_mask[k];
			}
		}
	}

	fprintf(stderr, "RANSAC found this best model:");
	for (int i = 0; i < modeldim; i++)
		fprintf(stderr, " %g", best_model[i]);
	fprintf(stderr, "\n");
	if (0) {
		FILE *f = xfopen("/tmp/ramo.txt", "w");
		for (int i = 0; i < modeldim; i++)
			fprintf(f,"%lf%c",best_model[i],i==modeldim-1?'\n':' ');
		xfclose(f);
	}
	//fprintf(stderr, "errors of outliers:");
	//for (int i = 0; i < n; i++)
	//	if (!best_mask[i]) {
	//		float e = mev(best_model, data+i*datadim, usr);
	//		fprintf(stderr, " %g", e);
	//	}
	//fprintf(stderr, "\n");
	//fprintf(stderr, "errors of inliers:");
	//for (int i = 0; i < n; i++)
	//	if (best_mask[i]) {
	//		float e = mev(best_model, data+i*datadim, usr);
	//		fprintf(stderr, " %g", e);
	//	}
	//fprintf(stderr, "\n");
	//fprintf(stderr, "errors of data points:\n");
	//for (int i = 0; i < n; i++) {
	//	float e = mev(best_model, data+i*datadim, usr);
	//	fprintf(stderr, "\t%g\t%s\n", e, best_mask[i]?"GOOD":"bad");
	//}

	int return_value = 0;
	if (best_ninliers >= min_inliers)
	{
		return_value =  best_ninliers;
	} else
		return_value = 0;

	for (int j = 0; j < modeldim; j++)
		if (!isfinite(best_model[j]))
			fail("model_%d not finite", j);

	if (out_model)
		for(int j = 0; j < modeldim; j++)
			out_model[j] = best_model[j];
	if (out_mask)
		for(int j = 0; j < n; j++)
			out_mask[j] = best_mask[j];

	free(best_mask);
	free(tmp_mask);

	return return_value;
}

#ifndef OMIT_MAIN

#include <stdio.h>
#include <string.h>

#include "ransac_cases.c" // example functions for RANSAC input
#include "parsenumbers.c" // function "read_ascii_floats"

int main_cases(int c, char *v[])
{
	if (c != 6 && c != 7 && c != 8) {
		fprintf(stderr, "usage:\n\t%s {line,aff,affn,fm} "
		//                         0   1
		"ntrials maxerr minliers omodel [omask [oinliers]] <data\n",*v);
		//2      3      4        5       6      7
		return EXIT_FAILURE;
	}

	// parse input options
	char *model_id = v[1];
	int ntrials = atoi(v[2]);
	float maxerr = atof(v[3]);
	int minliers = atoi(v[4]);
	char *filename_omodel = c > 5 ? v[5] : 0;
	char *filename_omask = c > 6 ? v[6] : 0;
	char *filename_inliers = c > 7 ? v[7] : 0;

	// declare context variables
	int modeldim, datadim, nfit;
	ransac_error_evaluation_function *model_evaluation;
	ransac_model_generating_function *model_generation;
	ransac_model_accepting_function *model_acceptation = NULL;
	void *user_data = NULL;


	// fill context variables according to ransac case
	if (0 == strcmp(model_id, "line")) {
		datadim = 2;
		modeldim = 3;
		nfit = 2;
		model_evaluation = distance_of_point_to_straight_line;
		model_generation = straight_line_through_two_points;

	} else if (0 == strcmp(model_id, "aff")) {
		datadim = 4;
		modeldim = 6;
		nfit = 3;
		model_evaluation = affine_match_error;
		model_generation = affine_map_from_three_pairs;

	} else if (0 == strcmp(model_id, "affn")) {
		datadim = 4;
		modeldim = 6;
		nfit = 3;
		model_evaluation = affine_match_error;
		model_generation = affine_map_from_three_pairs;
		model_acceptation = affine_map_is_reasonable;

	} else if (0 == strcmp(model_id, "hom")) {
		datadim = 4;
		modeldim = 9;
		nfit = 4;
		model_evaluation = homographic_match_error;
		model_generation = homography_from_four;
		//model_acceptation = homography_is_reasonable;
		model_acceptation = NULL;

	} else if (0 == strcmp(model_id, "aff3d")) {
		datadim = 6;
		modeldim = 12;
		nfit = 4;
		model_evaluation = affine3d_match_error;
		model_generation = affine3d_map_from_four_pairs;

	} else if (0 == strcmp(model_id, "fm")) { // fundamental matrix
		datadim = 4;
		modeldim = 9;
		nfit = 7;
		model_evaluation = epipolar_error;
		model_generation = seven_point_algorithm;
		//model_acceptation = fundamental_matrix_is_reasonable;

	} else if (0 == strcmp(model_id, "fmn")) { // fundamental matrix
		int main_hack_fundamental_matrix(int,char*[]);
		return main_hack_fundamental_matrix(c-1, v+1);
	} else if (0 == strcmp(model_id, "fmnt")) { // fundamental matrix
		int main_hack_fundamental_trimatrix(int,char*[]);
		return main_hack_fundamental_trimatrix(c-1, v+1);
	} else {
		printf("unrecognized model \"%s\"\n", model_id);
		return EXIT_FAILURE;
	}

	// read input data
	int n;
	float *data = read_ascii_floats(stdin, &n);
	n /= datadim;

	// call the ransac function to fit a model to data
	float model[modeldim];
	bool *mask = xmalloc(n * sizeof*mask);
	int n_inliers = ransac(mask, model, data, datadim, n, modeldim,
			model_evaluation, model_generation,
			nfit, ntrials, minliers, maxerr,
			model_acceptation, user_data);


	// print a summary of the results
	if (n_inliers > 0) {
		printf("RANSAC found a model with %d inliers\n", n_inliers);
		printf("parameters =");
		for (int i = 0; i < modeldim; i++)
			printf(" %g", model[i]);
		printf("\n");
	} else printf("RANSAC found no model\n");

	// if requested, save the inlying data points
	if (filename_inliers) {
		FILE *f = xfopen(filename_inliers, "w");
		for (int i = 0; i < n; i++)
			if (mask[i]) {
				for(int d = 0; d < datadim; d++)
					fprintf(f,"%g ",data[i*datadim+d]);
				fprintf(f,"\n");
			}
		xfclose(f);
	}

	// if requested, save the inlier mask
	if (filename_omask) {
		FILE *f = xfopen(filename_omask, "w");
		for (int i = 0; i < n; i++)
			fprintf(f, mask[i]?" 1":" 0");
		fprintf(f, "\n");
		xfclose(f);
	}

	// if requested, save the model
	if (filename_omodel) {
		FILE *f = xfopen(filename_omodel, "w");
		for (int i = 0; i < modeldim; i++)
			fprintf(f, "%lf%c", model[i], i==modeldim-1?'\n':' ');
		xfclose(f);
	}


	return EXIT_SUCCESS;
}

int main_hack_fundamental_matrix(int c, char *v[])
{
	if (c < 4) {
		fprintf(stderr, "usage:\n\t%s fmn "
		//                         -1   0
			"ntrials maxerr minliers [inliers] <data\n", *v);
		//       1       2      3        4
		return EXIT_FAILURE;
	}
	int ntrials = atoi(v[1]);
	float maxerr = atof(v[2]);
	int minliers = atoi(v[3]);
	char *inliers_filename = v[4];

	fprintf(stderr, "WARNING: ignoring parameter minliers=%d\n", minliers);

	int datadim = 4;
	int modeldim = 9;
	int n;
	float *data = read_ascii_floats(stdin, &n);
	n /= datadim;

	float model[modeldim];
	bool *mask = xmalloc(n * sizeof*mask);
	for (int i = 0; i < n; i++) mask[i] = false;
	int n_inliers = find_fundamental_matrix_by_ransac(mask, model,
			data, n, ntrials, maxerr);

	// print a summary of the results
	if (n_inliers > 0) {
		printf("RANSAC found a model with %d inliers\n", n_inliers);
		printf("parameters =");
		for (int i = 0; i < modeldim; i++)
			printf(" %g", model[i]);
		printf("\n");
		//fprintf(stderr, "errors of data points:\n");
		//ransac_error_evaluation_function *ep = epipolar_error;
		//for (int i = 0; i < n; i++) {
		//	float e = ep(model, data+i*datadim, NULL);
		//	fprintf(stderr, "\t%g\t%s\n", e, mask[i]?"GOOD":"bad");
		//}
	} else printf("RANSAC found no model\n");

#if 0
	float gamod[9] = {2.8709e-09,  -4.3669e-08,   1.0966e-02,
   4.0071e-08,   3.6767e-10,   3.4892e-03,
  -1.2060e-02,  -4.3969e-03,   1.0000e+00};
	printf("as a comparison, use a fixed model:");
	printf("parameters =");
	for (int i = 0; i < modeldim; i++)
		printf(" %g", gamod[i]);
	printf("\n");
	fprintf(stderr, "errors of data points:\n");
	ransac_error_evaluation_function *ep = epipolar_error;
	for (int i = 0; i < n; i++) {
		float e = ep(gamod, data+i*datadim, NULL);
		fprintf(stderr, "\t%g\t%s\n", e,"unknown");
	}
#endif

	// if needed, print the inliers
	if (inliers_filename) {
		FILE *f = xfopen(inliers_filename, "w");
		for (int i = 0; i < n; i++)
			if (mask[i]) {
				for(int d = 0; d < datadim; d++)
					fprintf(f,"%g ",data[i*datadim+d]);
				fprintf(f,"\n");
			}
		xfclose(f);
	}
	if (false) {
		FILE *f = xfopen("/tmp/omask.txt", "w");
		for (int i = 0; i < n; i++)
			fprintf(f, mask[i]?" 1":" 0");
		fprintf(f, "\n");
		xfclose(f);
	}

	free(mask);
	free(data);

	return EXIT_SUCCESS;
}

int main_hack_fundamental_trimatrix(int c, char *v[])
{
	if (c < 4) {
		fprintf(stderr, "usage:\n\t%s fmnt "
		//                         -1   0
			"ntrials maxerr minliers [inliers] <data\n", *v);
		//       1       2      3        4
		return EXIT_FAILURE;
	}
	int ntrials = atoi(v[1]);
	float maxerr = atof(v[2]);
	int minliers = atoi(v[3]);
	char *inliers_filename = v[4];

	fprintf(stderr, "WARNING: ignoring parameter minliers=%d\n", minliers);

	int datadim = 6;
	int modeldim = 18;
	int n;
	float *data = read_ascii_floats(stdin, &n);
	n /= datadim;

	float model[modeldim];
	bool *mask = xmalloc(n * sizeof*mask);
	for (int i = 0; i < n; i++) mask[i] = false;
	int n_inliers = find_fundamental_pair_by_ransac(mask, model,
			data, n, ntrials, maxerr);

	// print a summary of the results
	if (n_inliers > 0) {
		printf("RANSAC found a model with %d inliers\n", n_inliers);
		printf("parameters =");
		for (int i = 0; i < modeldim; i++)
			printf(" %g", model[i]);
		printf("\n");
		//fprintf(stderr, "errors of data points:\n");
		//ransac_error_evaluation_function *ep = epipolar_error_triplet;
		//for (int i = 0; i < n; i++) {
		//	float e = ep(model, data+i*datadim, NULL);
		//	fprintf(stderr, "\t%g\t%s\n", e, mask[i]?"GOOD":"bad");
		//}
	} else printf("RANSAC found no model\n");

#if 0
	float gamod[9] = {2.8709e-09,  -4.3669e-08,   1.0966e-02,
   4.0071e-08,   3.6767e-10,   3.4892e-03,
  -1.2060e-02,  -4.3969e-03,   1.0000e+00};
	printf("as a comparison, use a fixed model:");
	printf("parameters =");
	for (int i = 0; i < modeldim; i++)
		printf(" %g", gamod[i]);
	printf("\n");
	fprintf(stderr, "errors of data points:\n");
	ransac_error_evaluation_function *ep = epipolar_error;
	for (int i = 0; i < n; i++) {
		float e = ep(gamod, data+i*datadim, NULL);
		fprintf(stderr, "\t%g\t%s\n", e,"unknown");
	}
#endif

	// if needed, print the inliers
	if (inliers_filename) {
		FILE *f = xfopen(inliers_filename, "w");
		for (int i = 0; i < n; i++)
			if (mask[i]) {
				for(int d = 0; d < datadim; d++)
					fprintf(f,"%g ",data[i*datadim+d]);
				fprintf(f,"\n");
			}
		xfclose(f);
	}
	if (false) {
		FILE *f = xfopen("/tmp/omask.txt", "w");
		for (int i = 0; i < n; i++)
			fprintf(f, mask[i]?" 1":" 0");
		fprintf(f, "\n");
		xfclose(f);
	}


	free(mask);
	free(data);

	return EXIT_SUCCESS;
}



int main(int c, char *v[])
{
	return main_cases(c, v);
}
#endif//OMIT_MAIN
