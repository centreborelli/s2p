// this is common code for the "synflow" and "paraflow" programs


#include "fail.c"
#include "getpixel.c"
#include "marching_interpolation.c"

static float interpolate_bilinear(float a, float b, float c, float d,
					float x, float y)
{
	float r = 0;
	r += a*(1-x)*(1-y);
	r += b*(1-x)*(y);
	r += c*(x)*(1-y);
	r += d*(x)*(y);
	return r;
}

static float interpolate_nearest(float a, float b, float c, float d,
					float x, float y)
{
	if (x<0.5) return y<0.5 ? a : b;
	else return y<0.5 ? c : d;
}

static float interpolate_cell(float a, float b, float c, float d,
					float x, float y, int method)
{
	switch(method) {
	case 0: return interpolate_nearest(a, b, c, d, x, y);
	case 1: return marchi(a, b, c, d, x, y);
	case 2: return interpolate_bilinear(a, b, c, d, x, y);
	default: fail("caca de vaca");
	}
}

static void general_interpolate(float *result,
		float *x, int w, int h, int pd, float p, float q,
		int m) // method
{
	getsample_operator P = get_sample_operator(getsample_0);
	int ip = floor(p);
	int iq = floor(q);
	FORL(pd) {
		float a = P(x, w, h, pd, ip  , iq  , l);
		float b = P(x, w, h, pd, ip  , iq+1, l);
		float c = P(x, w, h, pd, ip+1, iq  , l);
		float d = P(x, w, h, pd, ip+1, iq+1, l);
		float v = interpolate_cell(a, b, c, d, p-ip, q-iq, m);
		result[l] = v;
	}
}

static void affine_map(double y[2], double A[6], double x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

static void projective_map(double y[2], double H[9], double x[2])
{
	double z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	y[1] = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
}

// (x,y) |-> (x+(a*R^2), y+(a*R^2))
// R^2 = x*x + y*y
// x' = x + x*(a*x*x + a*y*y)
// y' = y + y*(a*x*x + a*y*y)
// X' = X + a*X*X*X

//#include <complex.h>
static double solvecubicspecial(double a, double b)
{
	if (fabs(a) < 1e-29) {
		return b;
	}
	long double x;
	long double r;
	if (a < 0) {
// a*X*X*X + X - X' = 0
// X*X*X + (1/a)*X - X'/a = 0
// p = 1/a, q=b/a;
		long double p = 1/a;
		long double q = -b/a;
		int k = 1;
		long double cosarg = acos((3*q)/(2*p)*sqrt(-3/p))/3-k*2*M_PI/3;
		r = 2*sqrt(-p/3) * cos(cosarg);
		return r;
	} else {
		x = cbrt(sqrt((27*a*b*b+4)/a)/(2*sqrt(27)*a)+b/(2*a));
		r = x-1/(3*a*x);
		return r;
	}
//	double x = sqrt((27*a*b + 4)/a) / (2*sqrt(27)*a) + b/(2*a);
//	double y = pow(x, 1.0/3);
//	//double r = y - 1/(3*a*y);
//	double r = (3*a*y*y - 1)/(3*a*y);
//	return r;
}

// X' = X + a*X*X*X
// a*X*X*X + X - X' = 0
static double invertparabolicdistortion(double a, double xp)
{
	return solvecubicspecial(a, xp);
}

// X' = X + a*X*X*X
static double parabolicdistortion(double a, double x)
{
	return x + a*x*x*x;
}

static void invert_affinity(double invA[6], double A[6])
{
	double a, b, c, d, p, q;
	a=A[0]; b=A[1]; p=A[2];
	c=A[3]; d=A[4]; q=A[5];
	double det = a*d - b*c;
	invA[0] = d;
	invA[1] = -b;
	invA[2] = b*q-d*p;
	invA[3] = -c;
	invA[4] = a;
	invA[5] = c*p-a*q;
	FORI(6) invA[i] /= det;
}

// computes the affinity sending [0,0], [1,0] and [0,1] to x, y and z
static void affinity_from_3pt(double x[2], double y[2], double z[2],
		double A[6])
{
	double a, b, c, d, p, q;
	p = x[0];
	q = x[1];
	a = y[0] - p;
	c = y[1] - q;
	b = z[0] - p;
	d = z[1] - q;
	A[0]=a; A[1]=b; A[2]=p;
	A[3]=c; A[4]=d; A[5]=q;
}

// computes the affinity sending a,b,c to x,y,z
static void affinity_from_3corresp(double a[2], double b[2], double c[2],
		double x[2], double y[2], double z[2], double R[6])
{
	double A[6], B[6], iA[6];
	affinity_from_3pt(a, b, c, A);
	affinity_from_3pt(x, y, z, B);
	invert_affinity(iA, A);
	R[0] = B[0]*iA[0] + B[1]*iA[3];
	R[1] = B[0]*iA[1] + B[1]*iA[4];
	R[2] = B[0]*iA[2] + B[1]*iA[5] + B[2];
	R[3] = B[3]*iA[0] + B[4]*iA[3];
	R[4] = B[3]*iA[1] + B[4]*iA[4];
	R[5] = B[3]*iA[2] + B[4]*iA[5] + B[5];
}

#include "vvector.h"

static void invert_homography9(double invH[9], double H[9])
{
	double h[3][3] = { {H[0], H[1], H[2]},
			{H[3], H[4], H[5]},
			{H[6], H[7], H[8]}};
	double det, ih[3][3];
	INVERT_3X3(ih, det, h);
	FORI(9) invH[i] = ih[i/3][i%3];
}

#include "homographies.c"

static double produce_homography(double H[9], int w, int h,
		char *homtype, double *v)
{
	if (0 == strcmp(homtype, "hom")) { // actual parameters
		FORI(9) H[i] = v[i];
	} else if (0 == strcmp(homtype, "homi")) {
		invert_homography9(H, v);
	} else if (0 == strcmp(homtype, "shomi")) {
		invert_homography9(H, 1+v);
		H[2] *= *v;
		H[5] *= *v;
		H[6] /= *v;
		H[7] /= *v;
	} else if (0 == strcmp(homtype, "hom4p")) {
		// absolute displacement of the image corners
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + v[0], 0 + v[1]},
			{w + v[2], 0 + v[3]},
			{0 + v[4], h + v[5]},
			{w + v[6], h + v[7]}
		};
		double R[3][3];
		homography_from_eight_points(R,
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3]);
		FORI(9) H[i] = R[i/3][i%3];
	} else if (0 == strcmp(homtype, "hom4pr")) {
		// absolute displacement of the image corners
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + w*v[0], 0 + h*v[1]},
			{w + w*v[2], 0 + h*v[3]},
			{0 + w*v[4], h + h*v[5]},
			{w + w*v[6], h + h*v[7]}
		};
		double R[3][3];
		homography_from_eight_points(R,
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3]);
		FORI(9) H[i] = R[i/3][i%3];
	} else if (0 == strcmp(homtype, "hom4prc")) {
		// percentual relative displacement of the image corners
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + w*v[0]/100, 0 + h*v[1]/100},
			{w + w*v[2]/100, 0 + h*v[3]/100},
			{0 + w*v[4]/100, h + h*v[5]/100},
			{w + w*v[6]/100, h + h*v[7]/100}
		};
		double R[3][3];
		homography_from_eight_points(R,
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3]);
		FORI(9) H[i] = R[i/3][i%3];
	} else if (0 == strcmp(homtype, "hom16")) {
		// absolute coordinates of 4 point pairs
		double a[4][2] = {
			{v[0],v[1]},{v[2],v[3]},{v[4],v[5]},{v[6],v[7]}
		};
		double b[4][2] = {
			{v[8],v[9]},{v[10],v[11]},{v[12],v[13]},{v[14],v[15]}
		};
		double R[3][3];
		homography_from_eight_points(R, a[0], a[1], a[2], a[3],
				b[0], b[1], b[2], b[3]);
		FORI(9) H[i] = R[i/3][i%3];
	} else if (0 == strcmp(homtype, "hom16r")) {
		// relative coordinates of 4 point pairs
		double a[4][2] = {
			{w*v[0],h*v[1]},
			{w*v[2],h*v[3]},
			{w*v[4],h*v[5]},
			{w*v[6],h*v[7]}
		};
		double b[4][2] = {
			{w*v[8], h*v[9]},
			{w*v[10],h*v[11]},
			{w*v[12],h*v[13]},
			{w*v[14],h*v[15]}
		};
		double R[3][3];
		homography_from_eight_points(R, a[0], a[1], a[2], a[3],
				b[0], b[1], b[2], b[3]);
		FORI(9) H[i] = R[i/3][i%3];
	} else if (0 == strcmp(homtype, "hom16rc")) {
		// percentual relative coordinates of 4 point pairs
		double a[4][2] = {
			{w*v[0]/100,h*v[1]/100},
			{w*v[2]/100,h*v[3]/100},
			{w*v[4]/100,h*v[5]/100},
			{w*v[6]/100,h*v[7]/100}
		};
		double b[4][2] = {
			{w*v[8] /100,h*v[9] /100},
			{w*v[10]/100,h*v[11]/100},
			{w*v[12]/100,h*v[13]/100},
			{w*v[14]/100,h*v[15]/100}
		};
		double R[3][3];
		homography_from_eight_points(R, a[0], a[1], a[2], a[3],
				b[0], b[1], b[2], b[3]);
		FORI(9) H[i] = R[i/3][i%3];
	} else fail("unrecognized homography type \"%s\"", homtype);
	return 0;
}

static double produce_affinity(double A[6], int w, int h,
		char *afftype, double *v)
{
	int W = w - 1;
	int H = h - 1;
	if (0 == strcmp(afftype, "aff")) { // actual parameters
		FORI(6) A[i] = v[i];
	} else if (0 == strcmp(afftype, "aff3p")) {
		// absolute displacement of the image corners
		double corner[3][2] = {{0,0}, {W,0}, {0,H}};
		double other[3][2] = {
			{0 + v[0], 0 + v[1]},
			{W + v[2], 0 + v[3]},
			{0 + v[4], H + v[5]}
		};
		affinity_from_3corresp(corner[0], corner[1], corner[2],
					other[0], other[1], other[2], A);
	} else if (0 == strcmp(afftype, "aff3pr")) {
		// relative displacement of the image corners
		double corner[3][2] = {{0,0}, {W,0}, {0,H}};
		double other[3][2] = {
			{0 + W*v[0], 0 + H*v[1]},
			{W + W*v[2], 0 + H*v[3]},
			{0 + W*v[4], H + H*v[5]}
		};
		affinity_from_3corresp(corner[0], corner[1], corner[2],
					other[0], other[1], other[2], A);
	} else if (0 == strcmp(afftype, "aff3prc")) {
		// percentual relative displacement of the image corners
		double corner[3][2] = {{0,0}, {w,0}, {0,h}};
		double other[3][2] = {
			{0 + W*v[0]/100, 0 + H*v[1]/100},
			{W + W*v[2]/100, 0 + H*v[3]/100},
			{0 + W*v[4]/100, H + H*v[5]/100}
		};
		affinity_from_3corresp(corner[0], corner[1], corner[2],
					other[0], other[1], other[2], A);
	} else if (0 == strcmp(afftype, "aff12")) {
		// absolute coordinates of 3 point pairs
		double a[3][2] = { {v[0],v[1]},{v[2],v[3]},{v[4],v[5]} };
		double b[4][2] = { {v[6],v[7]},{v[8],v[9]},{v[10],v[11]} };
		affinity_from_3corresp(a[0], a[1], a[2], b[0], b[1], b[2], A);
	} else if (0 == strcmp(afftype, "aff12r")) {
		// relative coordinates of 3 point pairs
		double a[3][2] = {
			{W*v[0],H*v[1]},
			{W*v[2],H*v[3]},
			{W*v[4],H*v[5]}
		};
		double b[3][2] = {
			{W*v[6], H*v[7]},
			{W*v[8], H*v[9]},
			{W*v[10],H*v[11]}
		};
		affinity_from_3corresp(a[0], a[1], a[2], b[0], b[1], b[2], A);
	} else if (0 == strcmp(afftype, "aff12rc")) {
		// percentual relative coordinates of 3 point pairs
		double a[3][2] = {
			{W*v[0]/100,H*v[1]/100},
			{W*v[2]/100,H*v[3]/100},
			{W*v[4]/100,H*v[5]/100}
		};
		double b[3][2] = {
			{W*v[6] /100,H*v[7] /100},
			{W*v[8] /100,H*v[9] /100},
			{W*v[10]/100,H*v[11]/100}
		};
		affinity_from_3corresp(a[0], a[1], a[2], b[0], b[1], b[2], A);
	} else if (0 == strcmp(afftype, "traslation")) {
		double R[6] = {1, 0, v[0], 0, 1, v[1]};
		FORI(6) A[i] = R[i];
	} else if (0 == strcmp(afftype, "traslation_rc")) {
		double R[6] = {1, 0, W*v[0]/100, 0, 1, H*v[1]/100};
		FORI(6) A[i] = R[i];
	} else if (0 == strcmp(afftype, "euclidean")) {
		double theta = v[0];
		double R[6] = {cos(theta), sin(theta), v[1],
			-sin(theta), cos(theta), v[2]};
		FORI(6) A[i] = R[i];
	} else if (0 == strcmp(afftype, "euclidean_rc")) {
		double theta = M_PI*v[0]/180;
		double c[2] = {W/2.0, H/2.0}, rc[2];
		double R[6] = {cos(theta), sin(theta), 0,
			-sin(theta), cos(theta), 0};
		affine_map(rc, R, c);
		R[2] = W*v[1]/100 - rc[0] + c[0];
		R[5] = H*v[2]/100 - rc[1] + c[1];
		FORI(6) A[i] = R[i];
	} else if (0 == strcmp(afftype, "similar")) {
		double theta = v[0];
		double rho = v[1];
		double R[6] = {rho*cos(theta), rho*sin(theta), v[2],
			-rho*sin(theta), rho*cos(theta), v[3]};
		FORI(6) A[i] = R[i];
	} else if (0 == strcmp(afftype, "similar_rc")) {
		double theta = M_PI*v[0]/180;
		double rho = v[1];
		double c[2] = {W/2.0, H/2.0}, rc[2];
		double R[6] = {rho*cos(theta), rho*sin(theta), 0,
			-rho*sin(theta), rho*cos(theta), 0};
		affine_map(rc, R, c);
		R[2] = W*v[2]/100 - rc[0] + c[0];
		R[5] = H*v[3]/100 - rc[1] + c[1];
		FORI(6) A[i] = R[i];
	} else if (0 == strcmp(afftype, "colorwheel")) {
		double disp = v[0];
		double scale = 1 + disp/100.0;
		fprintf(stderr, "scale = %g\n", scale);
		double R[6] = {
			scale, 0, (1-scale)*W/2.0,
			0, scale, (1-scale)*H/2.0
		};
		FORI(6) A[i] = R[i];
	} else fail("unrecognized affinity type \"%s\"", afftype);
	return 0;
}

static double produce_radial_model(double A[3], double iA[3], int w, int h,
		char *pradtype, double *v)
{
	if (0 == strcmp(pradtype, "pararaw")) { // actual parameters
		FORI(3) A[i] = v[i];
	} else if (0 == strcmp(pradtype, "pradial_pixelic")) {
		FORI(2) A[i] = v[i];
		double p = v[2];
		double a = p;
		A[2] = (8*a)/(1.0*w*w*w);
	} else if (0 == strcmp(pradtype, "pradialrc")) {
		A[0] = w/2.0 + w*v[9]/100.0;
		A[1] = h/2.0 + h*v[1]/100.0;
		double p = v[2];
		double a = p*w/100.0;
		A[2] = (8*a)/(1.0*w*w*w);
	} else fail("non-raw radial model not yet supported");
	if (A[2] < 0) {
		// TODO: some more fine-grained warnings about warping
		double a = -A[2];
		fprintf(stderr, "negative a = %g\n", a);
		if (hypot(h,w)/2.0 > 1/sqrt(3*a))
			fprintf(stderr, "WARNING: corner of image warped!\n");
		if (w/2.0 > 1/sqrt(3*a))
			fprintf(stderr, "WARNING: sides of image warped!\n");
		if (h/2.0 > 1/sqrt(3*a))
			fprintf(stderr, "WARNING: tops of image warped!\n");
	}
	return 0;
}

static double produce_iradial_model(double A[3], double iA[3], int w, int h,
		char *pradtype, double *v)
{
	if (0 == strcmp(pradtype, "ipararaw")) { // actual parameters
		FORI(3) iA[i] = A[i] = v[i];
	} else if (0 == strcmp(pradtype, "ipradial_pixelic")) {
		FORI(2) A[i] = v[i];
		double p = v[2];
		double a = p;
		A[2] = (8*a)/(1.0*w*w*w);
	} else if (0 == strcmp(pradtype, "ipradialrc")) {
		A[0] = w/2.0 + w*v[9]/100.0;
		A[1] = h/2.0 + h*v[1]/100.0;
		double p = v[2];
		double a = p*w/100.0;
		A[2] = (8*a)/(1.0*w*w*w);
	} else fail("non-raw radial model not yet supported");
	return 0;
}

static double produce_combi2_model(double H[15], int w, int h,
		char *combitype, double *v)
{
	if (0 == strcmp(combitype, "combi2raw")) {
		FORI(15) H[i] = v[i];
	} else if (0 == strcmp(combitype, "combi2rc")) {
		//error("nice parameters not implemented");
		double c[2][2] = {
			{w/2.0 + w*v[0]/100.0, h/2.0 + h*v[1]/100.0},
			{w/2.0 + w*v[11]/100.0, h/2.0 + h*v[12]/100.0}
		};
		double a[2] = {8*v[2]/(1.0*w*w*w), 8*v[13]/(1.0*w*w*w)};
		double corner[4][2] = {{0,0}, {w,0}, {0,h}, {w,h}};
		double other[4][2] = {
			{0 + w*v[3]/100, 0 + h*v[4]/100},
			{w + w*v[5]/100, 0 + h*v[6]/100},
			{0 + w*v[7]/100, h + h*v[8]/100},
			{w + w*v[9]/100, h + h*v[10]/100}
		};
		double R[3][3];
		homography_from_eight_points(R,
				corner[0], corner[1], corner[2], corner[3],
				other[0], other[1], other[2], other[3]);
		FORI(2) H[i] = c[0][i];
		H[2] = a[0];
		FORI(9) H[i+3] = R[i/3][i%3];
		FORI(2) H[12+i] = c[1][i];
		H[14] = a[1];
	}

	return 0;
}
#define SYNFLOW_MAXPARAM 40 // whatever


#define FLOWMODEL_HIDDEN_AFFINE 1     // 6 parameters
#define FLOWMODEL_HIDDEN_PROJECTIVE 2 // 9 parameters
#define FLOWMODEL_HIDDEN_PRADIAL 3    // 3 parameter
#define FLOWMODEL_HIDDEN_IPRADIAL 4   // 3 parameter
#define FLOWMODEL_HIDDEN_COMBI1 5     // 12 parameters, H(pradial(x))
#define FLOWMODEL_HIDDEN_COMBI2 6     // 12 parameters, ipradial(H(pradial(x)))
#define FLOWMODEL_HIDDEN_COMBI3 7     // 15 parameters, ipradial'(H(pradial(x)))
#define FLOWMODEL_HIDDEN_ICOMBI1 5    // 12 parameters, ipradial(H(x))
#define FLOWMODEL_HIDDEN_ICOMBI2 6    // 12 parameters, pradial(H(ipradial(x)))
#define FLOWMODEL_HIDDEN_ICOMBI3 7    // 15 parameters, pradial'(H(ipradial(x)))


// data structure to store models for parametric (synthethic) movements
// together with their inverses
struct flow_model {
	bool given_by_field;
	char *model_name; // only for reference

	int n; // number of "visible" parameters
	double p[SYNFLOW_MAXPARAM];  // "visible" parameters of forward model

	int hidden_id;
	int nh; // number of "hidden" parameters
	double H[SYNFLOW_MAXPARAM];  // "hidden" parameters of forward model
	double iH[SYNFLOW_MAXPARAM]; // "hidden" parameters of backward model

	int w, h;
	float (*field)[2];
};

// returns the number of hidden parameters
static int parse_flow_name(int *vp, int *hidden_id, char *model_name)
{
	struct {
		char *model_name;
		int visible_params, hidden_params, hidden_id;
	} name_data[] = {
		{"traslation",    2, 6, FLOWMODEL_HIDDEN_AFFINE},
//		{"traslation_r",  2, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"traslation_rc", 2, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"euclidean",     3, 6, FLOWMODEL_HIDDEN_AFFINE},
//		{"euclidean_r",   3, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"euclidean_rc",  3, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"similar",       4, 6, FLOWMODEL_HIDDEN_AFFINE},
//		{"similar_r",     4, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"similar_rc",    4, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"aff",           6, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"aff3p",         6, 6, FLOWMODEL_HIDDEN_AFFINE},
//		{"aff3pr",        6, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"aff3prc",       6, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"aff12",         12, 6, FLOWMODEL_HIDDEN_AFFINE},
//		{"aff12r",        12, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"aff12rc",       12, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"colorwheel",    1, 6, FLOWMODEL_HIDDEN_AFFINE},
		{"hom",           9, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"homi",          9, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"shomi",         10, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"hom4p",         8, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
//		{"hom4pr",        8, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"hom4prc",       8, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"hom16",         16, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
//		{"hom16r",        16, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"hom16rc",       16, 9, FLOWMODEL_HIDDEN_PROJECTIVE},
		{"pararaw",       3, 3, FLOWMODEL_HIDDEN_PRADIAL},
		{"ipararaw",       3, 3, FLOWMODEL_HIDDEN_IPRADIAL},
		{"cpradial",      1, 3, FLOWMODEL_HIDDEN_PRADIAL},
//		{"cpradialr",     1, 3, FLOWMODEL_HIDDEN_PRADIAL},
		{"cpradialrc",    1, 3, FLOWMODEL_HIDDEN_PRADIAL},
		{"pradial_pixelic",       3, 3, FLOWMODEL_HIDDEN_PRADIAL},
		{"pradial",       3, 3, FLOWMODEL_HIDDEN_PRADIAL},
//		{"pradialr",      3, 3, FLOWMODEL_HIDDEN_PRADIAL},
		{"pradialrc",     3, 3, FLOWMODEL_HIDDEN_PRADIAL},
		{"cipradial",     1, 3, FLOWMODEL_HIDDEN_IPRADIAL},
//		{"cipradialr",    1, 3, FLOWMODEL_HIDDEN_IPRADIAL},
		{"cipradialrc",   1, 3, FLOWMODEL_HIDDEN_IPRADIAL},
		{"ipradial_pixelic",      3, 3, FLOWMODEL_HIDDEN_IPRADIAL},
//		{"ipradialr",     3, 3, FLOWMODEL_HIDDEN_IPRADIAL},
		{"ipradialrc",    3, 3, FLOWMODEL_HIDDEN_IPRADIAL},
		{"combi2raw",     15, 15, FLOWMODEL_HIDDEN_COMBI2},
		{"combi2rc",      14, 15, FLOWMODEL_HIDDEN_COMBI2},
		// combined, etc
		{"",0,0,0}
	};
	int i = 0;
	while (name_data[i].visible_params) {
		if (0 == strcmp(model_name, name_data[i].model_name)) {
			*hidden_id = name_data[i].hidden_id;
			*vp = name_data[i].visible_params;
			return name_data[i].hidden_params;
		}
		i = i + 1;
	}
	fail("unrecognized transform name %s", model_name);
	//return 0;
}

// "API"
// fills H, iH and other fields
static void produce_flow_model(struct flow_model *f,
		double *p, int np, char *name, int w, int h)
{
	f->given_by_field = false;
	int vp;
	f->nh = parse_flow_name(&vp, &f->hidden_id, name);
	if (vp != np)
		fail("flow model \"%s\" expects %d parameters but got %d",
				name, vp, np);
	f->model_name = name;
	f->w = w;
	f->h = h;
	f->n = np;
	FORI(np) f->p[i] = p[i];

	if (f->hidden_id == FLOWMODEL_HIDDEN_PROJECTIVE) {
		assert(f->nh == 9);
		double H[9], invH[9];
		produce_homography(H, w, h, f->model_name, f->p);
		invert_homography9(invH, H);
		FORI(9) f->H[i] = H[i];
		FORI(9) f->iH[i] = invH[i];
	} else if (f->hidden_id == FLOWMODEL_HIDDEN_AFFINE) {
		assert(f->nh == 6);
		double H[6], invH[6];
		produce_affinity(H, w, h, f->model_name, f->p);
		invert_affinity(invH, H);
		FORI(6) f->H[i] = H[i];
		FORI(6) f->iH[i] = invH[i];
	} else if (f->hidden_id == FLOWMODEL_HIDDEN_PRADIAL) {
		assert(f->nh == 3);
		double H[3], invH[3];
		produce_radial_model(H, invH, w, h, f->model_name, f->p);
		FORI(3) f->H[i] = H[i];
		FORI(3) f->iH[i] = invH[i];
	} else if (f->hidden_id == FLOWMODEL_HIDDEN_IPRADIAL) {
		assert(f->nh == 3);
		double H[3], invH[3];
		produce_iradial_model(H, invH, w, h, f->model_name, f->p);
		FORI(3) f->H[i] = H[i];
		FORI(3) f->iH[i] = invH[i];
	} else if (f->hidden_id == FLOWMODEL_HIDDEN_COMBI2) {
		assert(f->nh == 15);
		double H[15], invH[9];
		produce_combi2_model(H, w, h, f->model_name, f->p);
		invert_homography9(invH, H+3);
		FORI(15) f->H[i] = H[i];
		FORI(9) f->iH[i] = invH[i];
	} else fail("flow model \"%s\" not yet implemented", name);

	//if (SYNFLOW_VERBOSE()) {
	//	FORI(np) fprintf(stderr, "pfm p[%d] = %g\n", i, f->p[i]);
	//	FORI(f->nh) fprintf(stderr, "H[%d] = %g\n", i, f->H[i]);
	//	FORI(f->nh) fprintf(stderr, "invH[%d] = %g\n", i, f->iH[i]);
	//}
}

static void apply_flowmodel_affine(float y[2], float x[2], double A[6])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

static void apply_flowmodel_projective(float y[2], float x[2], double H[9])
{
	float z[3];
	z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
	z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
	z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = z[0]/z[2];
	y[1] = z[1]/z[2];
}

static void apply_flowmodel_pradial(float y[2], float x[2], double p[3],
		bool inv)
{
	double c[2] = {p[0], p[1]};
	double a = p[2];
	double r = hypot(x[0] - c[0], x[1] - c[1]);
	double R = inv ? invertparabolicdistortion(a, r) :
						parabolicdistortion(a, r);
	if (r > 0.000001)
		FORL(2) y[l] = c[l] + (R/r)*(x[l] - c[l]);
	else
		FORL(2) y[l] = 0;
}

// evaluate the flow vector at one given source point
static void apply_flow(float y[2], struct flow_model *f, float x[2], bool inv)
{
	switch(f->hidden_id) {
	case FLOWMODEL_HIDDEN_AFFINE: {
		assert(f->nh == 6);
		double *p = inv ? f->iH : f->H;
		apply_flowmodel_affine(y, x, p);
		break;
				      }
	case FLOWMODEL_HIDDEN_PROJECTIVE: {
		assert(f->nh == 9);
		double *p = inv ? f->iH : f->H;
		//fprintf(stderr, "H = %g %g %g %g %g %g %g %g %g\n",
		//	       	p[0], p[1], p[2],
		//	       	p[3], p[4], p[5],
		//	       	p[6], p[7], p[8]);
		apply_flowmodel_projective(y, x, p);
		break;
					  }
	case FLOWMODEL_HIDDEN_PRADIAL: {
		assert(f->nh == 3);
		apply_flowmodel_pradial(y, x, f->H, inv);
		break;
				       }
	case FLOWMODEL_HIDDEN_IPRADIAL: {
		assert(f->nh == 3);
		apply_flowmodel_pradial(y, x, f->H, !inv);
		break;
				       }
	case FLOWMODEL_HIDDEN_COMBI2: {
		assert(f->nh == 15);
		//double c[2][2] = {{f->H[0], f->H[1]}, {f->H[12], f->H[13]}};
		//double a[2] = {f->H[2], f->H[14]};
		float tmp[2], tmp2[2];
		double *H = inv ? f->iH : f->H + 3;
		apply_flowmodel_pradial(tmp, x, f->H, inv);
		apply_flowmodel_projective(tmp2, tmp, H);
		apply_flowmodel_pradial(y, tmp2, f->H + 12, !inv);
		break;
				      }
	default: fail("bizarre");
	}
}

// "API"
// fill a image with the vector field of the given flow
static void fill_flow_field(float *xx, struct flow_model *f, int w, int h)
{
	assert(f->w == w);
	assert(f->h == h);
	float (*x)[w][2] = (void*)xx;
	FORJ(h) FORI(w) {
		float p[2] = {i, j}, q[2];
		apply_flow(q, f, p, 0);
		FORL(2) x[j][i][l] = q[l] - p[l];
	}
}

// "API"
// morph an image according to a given flow model
static void transform_back(float *yy, struct flow_model *f, float *xx,
							int w, int h, int pd)
{
	float (*y)[w][pd] = (void *)yy;
	assert(f->w == w);
	assert(f->h == h);
	FORJ(h) FORI(w) {
		float p[2] = {i, j}, q[2];
		apply_flow(q, f, p, 0);
		float val[pd];
		general_interpolate(val, xx, w, h, pd, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}

// "API"
static void transform_forward(float *yy, struct flow_model *f, float *xx,
							int w, int h, int pd)
{
	float (*y)[w][pd] = (void *)yy;
	assert(f->w == w);
	assert(f->h == h);
	FORJ(h) FORI(w) {
		float p[2] = {i, j}, q[2];
		apply_flow(q, f, p, 1);
		float val[pd];
		general_interpolate(val, xx, w, h, pd, q[0], q[1], 2);
		FORL(pd)
			y[j][i][l] = val[l];
	}
}
