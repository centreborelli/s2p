// y = H(x)
static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

// compute the inverse homography (inverse of a 3x3 matrix)
static double invert_homography(double invH[3][3], double H[3][3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = ( a[4] * a[8] - a[5] * a[7] ) / det;
	r[1] = ( a[2] * a[7] - a[1] * a[8] ) / det;
	r[2] = ( a[1] * a[5] - a[2] * a[4] ) / det;
	r[3] = ( a[5] * a[6] - a[3] * a[8] ) / det;
	r[4] = ( a[0] * a[8] - a[2] * a[6] ) / det;
	r[5] = ( a[2] * a[3] - a[0] * a[5] ) / det;
	r[6] = ( a[3] * a[7] - a[4] * a[6] ) / det;
	r[7] = ( a[1] * a[6] - a[0] * a[7] ) / det;
	r[8] = ( a[0] * a[4] - a[1] * a[3] ) / det;
	return det;
}

// C = AoB, composition of two homographies (product of 3x3 matrices)
static void compose_homographies(double C[3][3], double A[3][3], double B[3][3])
{
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 3; k++)
			C[i][j] += A[i][k] * B[k][j];
	}
}

// Find the homography that changes the canonical projective basis
// into the given four points (x, y, z, w)
static void homography_from_four_points(double H[3][3],
		double x[2], double w[2], double z[2], double y[2])
{
	// We have to compute the following 9 coefficients of H:
	//
	//     a b p
	//     c d q
	//     r s t

	// fix the degree of freedom (assuming the four points are finite)
	double t = 1;

	// translation coefficients
	double p = x[0];
	double q = x[1];

	// "core" 2x2 system
	double A = w[0] - z[0];
	double B = y[0] - z[0];
	double C = w[1] - z[1];
	double D = y[1] - z[1];
	double P = z[0] - y[0] - w[0] + p;
	double Q = z[1] - y[1] - w[1] + q;
	double DET = A * D - B * C;
	double r = (D * P - B * Q) / DET;
	double s = (A * Q - C * P) / DET;

	// solve the rest of the diagonal system
	double a = w[0] * ( 1 + r ) - p;
	double b = y[0] * ( 1 + s ) - p;
	double c = w[1] * ( 1 + r ) - q;
	double d = y[1] * ( 1 + s ) - q;

	// fill-in the output
	H[0][0] = a; H[0][1] = b; H[0][2] = p;
	H[1][0] = c; H[1][1] = d; H[1][2] = q;
	H[2][0] = r; H[2][1] = s; H[2][2] = t;
}

// Find the homography that moves the four points (x,y,z,w) to (a,b,c,d)
static void homography_from_eight_points(double H[3][3],
		double x[2], double y[2], double z[2], double w[2],
		double a[2], double b[2], double c[2], double d[2])
{
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}
