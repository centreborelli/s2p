#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// 3x3 matrices are represented by arrays of length 9.
//                   a[0] a[1] a[2]
// The ordering  is: a[3] a[4] a[5]
//                   a[6] a[7] a[8]


void apply_homography(double y[2], double h[9], double x[2])
{
    double z = h[6]*x[0] + h[7]*x[1] + h[8];
    double tmp = x[0];  // to enable calls like apply_homography(x, h, x)
    y[0] = (h[0]*x[0] + h[1]*x[1] + h[2]) / z;
    y[1] = (h[3]*tmp  + h[4]*x[1] + h[5]) / z;
}


void matrix_33_inverse(double o[9], double i[9])
{
    double det = i[0]*i[4]*i[8] + i[2]*i[3]*i[7] + i[1]*i[5]*i[6]
               - i[2]*i[4]*i[6] - i[1]*i[3]*i[8] - i[0]*i[5]*i[7];
    o[0] = (i[4]*i[8] - i[5]*i[7]) / det;
    o[1] = (i[2]*i[7] - i[1]*i[8]) / det;
    o[2] = (i[1]*i[5] - i[2]*i[4]) / det;
    o[3] = (i[5]*i[6] - i[3]*i[8]) / det;
    o[4] = (i[0]*i[8] - i[2]*i[6]) / det;
    o[5] = (i[2]*i[3] - i[0]*i[5]) / det;
    o[6] = (i[3]*i[7] - i[4]*i[6]) / det;
    o[7] = (i[1]*i[6] - i[0]*i[7]) / det;
    o[8] = (i[0]*i[4] - i[1]*i[3]) / det;
}


void matrix_33_product(double ab[9], double a[9], double b[9])
{
    ab[0] = a[0]*b[0] + a[1]*b[3] + a[2]*b[6];
    ab[1] = a[0]*b[1] + a[1]*b[4] + a[2]*b[7];
    ab[2] = a[0]*b[2] + a[1]*b[5] + a[2]*b[8];
    ab[3] = a[3]*b[0] + a[4]*b[3] + a[5]*b[6];
    ab[4] = a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
    ab[5] = a[3]*b[2] + a[4]*b[5] + a[5]*b[8];
    ab[6] = a[6]*b[0] + a[7]*b[3] + a[8]*b[6];
    ab[7] = a[6]*b[1] + a[7]*b[4] + a[8]*b[7];
    ab[8] = a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
}


double min_n(double *x, int n)
{
    double out = *x;
    for (int i = 1; i < n; i++)
        if (*(++x) < out)
            out = *x;
    return out;
}


double max_n(double *x, int n)
{
    double out = *x;
    for (int i = 1; i < n; i++)
        if (*(++x) > out)
            out = *x;
    return out;
}


double *alloc_parse_doubles(int nmax, const char *ss, int *n)
{
    // add a space to s, so that gabriele's expression works as intended
    int ns = strlen(ss);
    char t[2+ns], *s = t;
    s[0] = ' ';
    for (int i = 0; i <= ns; i++)
        s[i+1] = ss[i];

    // parse the string
    double *r = (double *) malloc(nmax * sizeof*r);
    int i = 0, w;
    while (i < nmax && 1 == sscanf(s, "%*[][ \n\t,:;]%lf%n", r + i, &w)) {
        i += 1;
        s += w;
    }
    *n = i;
    return r;
}


// An affine fundamental matrix is a 3x3 matrix with 5 non-zero coefficients.
// It's represented here by an array of doubles of length 5.
//                     0    0  f[0]      0  0  c
// The ordering is:    0    0  f[1]   =  0  0  d
//                   f[2] f[3] f[4]      a  b  e
// with the notations from the reference book of Hartley and Zisserman (p. 345)

// A planar similarity is represented in homogeneous coordinates by a 3x3
// matrix of the form
//                     a  -b   p
//                     b   a   q
//                     0   0   1
// The rectification algorithm implemented by this function produces,
// by design, similarities with coefficient p equal to zero. Each
// output similarity is completely described by (a, b, q), thus is represented
// by an array of 3 floats with this convention: s[3] = {a, b, q}.
void rectifying_similarities_from_affine_fundamental_matrix(float s1[3],
                                                            float s2[3],
                                                            double f[5])
{
    // notations
    double c = f[0];
    double d = f[1];
    double a = f[2];
    double b = f[3];
    double e = f[4];

    // rotations
    double r = sqrt(a*a + b*b);
    double s = sqrt(c*c + d*d);
    //double r1[4] = { b/r, -a/r,   a/r,  b/r};
    //double r2[4] = {-d/s,  c/s,  -c/s, -d/s};

    // zoom and translation
    double z = sqrt(r / s);
    double t = 0.5 * e / sqrt(r * s);

    // output similarities
    s1[0] = (float) (z * b / r);
    s1[1] = (float) (z * a / r);
    s1[2] = (float) t;

    double zz = 1 / z;
    s2[0] = (float) (zz * -d / s);
    s2[1] = (float) (zz * -c / s);
    s2[2] = -t;
}
