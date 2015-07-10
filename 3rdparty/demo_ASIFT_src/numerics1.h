// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef _NUMERICS_H_
#define _NUMERICS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>


#include "library.h"

#define NRMAX(i,j) ( (i)<(j) ? (j):(i) )
#define NRMIN(i,j) ( (i)<(j) ? (i):(j) )
#define NRTINY 1.0e-10


// **********************************************
//  float **  basic functions
// **********************************************

float ** allocate_float_matrix(int nrows, int ncols);

void desallocate_float_matrix(float **matrix, int nrows, int ncols);

// **********************************************
//  LU based algorithms
// **********************************************

// Solves Ax=b by using lu decomposition 
// a matrix a[1..n][1..n] is replaced by the LU decompositions of a rowwise permutation of itself
// b[1..n] and x[1..n]
int lusolve(float **a, float *x, float *b, int n);


/*-- LU decomposition */
/* Given a matrix a[1..n][1..n] this routine replacess it by the LU decompositions of a rowwise permutation of itself. */
/* a and n are input, a is output, arranged as in equation (2.3.14) above; indx[1..n] in an output vector that records */
/* the row permutation effected by the partial pivoting; d is output as +-1 depending on whether the number of row     */
/* interchanges was even or odd respectively.                                                                          */
int ludcmp(float **a, int n, int *indx, float *d); /* LU decomposition */

/* Solves the set of n linear equations Ax=b. Here a[0..n-1][0..n-1] as input, not as the matrix A but rather as its LU decomposition,*/
/* determined by the routine ludcmp. indx[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is input as the    */
/* right hand side vector and returns with the solution vector x. */
void lubksb(float **a, int n, int *indx, float *b); /* LU linear solution */




#endif

