//
// C++ Implementation: stereomatch
//
// Description: eliminate the false matches with epipolar geometry constraint. 
//		See http://www.math-info.univ-paris5.fr/~moisan/epipolar/
//
// Copyright (c) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
// Changelog : 2011 Use Eigen SVD <Pierre Moulon>
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef STEREOMATCH_H
#define STEREOMATCH_H

#include <vector>

#include "libNumerics/numerics.h"
#include "libMatch/match.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include <cstdio>
#include <cmath>
#include <cstdlib>

/*-------------------- GENERAL PURPOSE ROUTINES --------------------*/

/* routines for vectors and matrices */

//float *vector(int nl, int nh);

float **matrix(int nrl, int nrh, int ncl, int nch);

void free_vector(float *v, int nl, int nh);

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);

/* Singular Value Decomposition routine */
void svdcmp(float **a, int m, int n, float *w, float **v);

/* Compute the real roots of a third order polynomial */
/* returns 1 or 3, the number of roots found */
int FindCubicRoots(float coeff[4], float x[3]);

/* logarithm (base 10) of binomial coefficient */
float logcombi(int k, int n);

/* tabulate logcombi(.,n) */
float *makelogcombi_n(int n);


/* tabulate logcombi(k,.) */
float *makelogcombi_k(int k, int nmax);


/* get a (sorted) random 7-uple of 0..n-1 */
void random_p7(int *k, int n);

/*-------------------- END OF GENERAL PURPOSE ROUTINES --------------------*/


/* float comparison for qsort() */
//According to http://www.cplusplus.com/reference/clibrary/cstdlib/qsort/, 
//we should have: void qsort ( void * base, size_t num, size_t size, int ( * comparator ) ( const void *, const void * ) ); that means, for "qsort", the "comparator" has two constant void* type input parameters
int compf(const void *i, const void *j);

void matcherrorn(float **F, const std::vector<float>& p1, const std::vector<float>& p2, float *e);

int epipolar(std::vector<float>& m1, std::vector<float>& m2, int *k, float *z, float **F1, float **F2);

float orsa(int width, int height, std::vector<Match>& match, std::vector<float>& index, int t_value, int verb_value, int n_flag_value, int mode_value, int stop_value);

#endif
