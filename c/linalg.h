#ifndef LINALG_H
#define LINALG_H

// 3x3 matrices are represented by arrays of length 9.
//                   a[0] a[1] a[2]
// The ordering  is: a[3] a[4] a[5]
//                   a[6] a[7] a[8]


// Apply a homography given by a 3x3 matrix h to a 2D point
// y: output point
// h: homography
// x: input point
void apply_homography(double y[2], double h[9], double x[2]);


// Compute the inverse of a 3x3 matrix
// o: output matrix
// i: input matrix
void matrix_33_inverse(double o[9], double i[9]);


// Compute the product of two 3x3 matrices
// ab: output product
// a, b: input matrices
void matrix_33_product(double ab[9], double a[9], double b[9]);


// Compute the minimum and maximum of on array of doubles
// x: pointer to the first element of the array
// n: length of the array
double min_n(double *x, int n);
double max_n(double *x, int n);


// Parse a string for doubles
// nmax: maximum number of doubles to read from the input string
// ss: input string
// *n: output number of doubles that have actually been read
double *alloc_parse_doubles(int nmax, const char *ss, int *n);

// Compute two rectifying similarities from an affine fundamental matrix
void rectifying_similarities_from_affine_fundamental_matrix(float s1[3],
                                                            float s2[3],
                                                            double f[5]);
#endif //LINALG_H
