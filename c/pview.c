// display points, pairs or triplets
// (optionally transformed, optionally with color masks)
// all programs below produce a color image as output
// the last two arguments "w h" are the size of the desired image
//
// pview points w h < points.txt
// pview pairs h1 ... h9 w h [mask.txt] < pairs.txt
// pview triplets w h < triplets.txt
// pview triplets w h [mask.txt] < triplets.txt
// pview epipolar f1 ... f9 w h < pairs.txt
// pview epipolar f1 ... f9 w h [mask.txt] < pairs.txt
// pview polygons w h < polygons.txt
//
// points.txt   = file with two columns of numbers (list of 2D points)
// pairs.txt    = file with four columns of numbers (list of 2D point pairs)
// triplets.txt = file with six columns of numbers (list of 2D point triplets)
// polygons.txt = file with one polygonal curve per line (list of 2D points)


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"
#include "drawsegment.c"
#include "pickopt.c"

struct rgb_value {
	uint8_t r, g, b;
};
struct rgba_value {
	uint8_t r, g, b, a;
};
#define RGBA_BLACK    (struct rgba_value){.r=0x00,.g=0x00,.b=0x00,.a=0xff}
#define RGBA_RED      (struct rgba_value){.r=0xff,.g=0x00,.b=0x00,.a=0xff}
#define RGBA_GREEN    (struct rgba_value){.r=0x00,.g=0xff,.b=0x00,.a=0xff}
#define RGBA_BLUE     (struct rgba_value){.r=0x00,.g=0x00,.b=0xff,.a=0xff}
#define RGBA_MAGENTA  (struct rgba_value){.r=0xff,.g=0x00,.b=0xff,.a=0xff}
#define RGBA_YELLOW   (struct rgba_value){.r=0xff,.g=0xff,.b=0x00,.a=0xff}
#define RGBA_PHANTOM  (struct rgba_value){.r=0x30,.g=0x20,.b=0x10,.a=0xff}
#define RGBA_BRIGHT   (struct rgba_value){.r= 129,.g=  86,.b=  43,.a=0xff}
#define RGBA_DARKGRAY (struct rgba_value){.r=0x30,.g=0x30,.b=0x30,.a=0xff}
#define RGBA_GRAY10   (struct rgba_value){.r=0x10,.g=0x10,.b=0x10,.a=0xff}
#define RGBA_GRAY20   (struct rgba_value){.r=0x20,.g=0x20,.b=0x20,.a=0xff}
#define RGBA_GRAY30   (struct rgba_value){.r=0x30,.g=0x30,.b=0x30,.a=0xff}
#define RGBA_GRAY40   (struct rgba_value){.r=0x40,.g=0x40,.b=0x40,.a=0xff}
#define RGBA_GRAY50   (struct rgba_value){.r=0x50,.g=0x50,.b=0x50,.a=0xff}
#define RGB_BLACK    (struct rgb_value){.r=0x00,.g=0x00,.b=0x00}
#define RGB_RED      (struct rgb_value){.r=0xff,.g=0x00,.b=0x00}
#define RGB_GREEN    (struct rgb_value){.r=0x00,.g=0xff,.b=0x00}
#define RGB_BLUE     (struct rgb_value){.r=0x00,.g=0x00,.b=0xff}
#define RGB_MAGENTA  (struct rgb_value){.r=0xff,.g=0x00,.b=0xff}
#define RGB_YELLOW   (struct rgb_value){.r=0xff,.g=0xff,.b=0x00}
#define RGB_PHANTOM  (struct rgb_value){.r=0x30,.g=0x20,.b=0x10}
#define RGB_BRIGHT   (struct rgb_value){.r= 129,.g=  86,.b=  43}
#define RGB_DARKGRAY (struct rgb_value){.r=0x30,.g=0x30,.b=0x30}
#define RGB_GRAY10   (struct rgb_value){.r=0x10,.g=0x10,.b=0x10}
#define RGB_GRAY20   (struct rgb_value){.r=0x20,.g=0x20,.b=0x20}
#define RGB_GRAY30   (struct rgb_value){.r=0x30,.g=0x30,.b=0x30}
#define RGB_GRAY40   (struct rgb_value){.r=0x40,.g=0x40,.b=0x40}
#define RGB_GRAY50   (struct rgb_value){.r=0x50,.g=0x50,.b=0x50}



#ifndef FORI
#define FORI(n) for(int i=0;i<(n);i++)
#endif//FORI
#ifndef FORJ
#define FORJ(n) for(int j=0;j<(n);j++)
#endif//FORJ


// decide whether a point falls within the image domain
static bool inner_point(int w, int h, int x, int y)
{
	return (x>=0) && (y>=0) && (x<w) && (y<h);
}

// CLI utility to view a set of planar points
// (produces a transparent PNG image)
static int main_viewp(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s sx sy < pairs.txt\n", *v);
		//                         0  1  2
		return EXIT_FAILURE;
	}
	int n;
	float *t = read_ascii_floats(stdin, &n);
	n /= 2;
	int sizex = atoi(v[1]);
	int sizey = atoi(v[2]);
	struct rgba_value (*x)[sizex] = xmalloc(4*sizex*sizey);
	FORI(sizex*sizey)
		x[0][i] = RGBA_BLACK;
	FORI(n)
	{
		int a = t[2*i+0];
		int b = t[2*i+1];
		if (inner_point(sizex, sizey, a, b))
			x[b][a] = RGBA_GREEN;
	}
	iio_write_image_uint8_vec("-", (uint8_t*)x, sizex, sizey, 4);
	free(t); free(x);
	return EXIT_SUCCESS;
}

// histogram of points
static int main_viewhp(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s sx sy < pairs.txt\n", *v);
		//                         0  1  2
		return EXIT_FAILURE;
	}
	int n;
	float *t = read_ascii_floats(stdin, &n);
	n /= 2;
	int sizex = atoi(v[1]);
	int sizey = atoi(v[2]);
	float (*x)[sizex] = xmalloc(sizex*sizey*sizeof(float));
	FORI(sizex*sizey)
		x[0][i] = 0;;
	FORI(n)
	{
		int a = t[2*i+0];
		int b = t[2*i+1];
		if (inner_point(sizex, sizey, a, b))
			x[b][a] += 1;
	}
	iio_write_image_float("-", x[0], sizex, sizey);
	free(t); free(x);
	return EXIT_SUCCESS;
}

bool identityP(double A[9])
{
	return A[0]==1 && A[1]==0 && A[2]==0 &&
	       A[3]==0 && A[4]==1 && A[5]==0 &&
	       A[6]==0 && A[7]==0 && A[8]==1;
}

// projective map
void projective_map(double y[2], double A[9], double x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
	double y2 = A[6]*x[0] + A[7]*x[1] + A[8];
	y[0] /= y2;
	y[1] /= y2;
}

// argument for the "traverse_segment" function
//static void put_pixel(int a, int b, void *pp)
//{
//	struct {int w, h; struct rgba_value *x, c;} *p = pp;
//	if (inner_point(p->w, p->h, a, b))
//		p->x[b*p->w+a] = p->c;
//}

#include "smapa.h"
SMART_PARAMETER_SILENT(LINN,70)

// linear combination of two intensities in linear intensity space
static double lincombin(double a, double b, double t)
{
	double linn = LINN();
	if (linn < -1) return t<0.5?a:b;
	if (linn < 0) return a*(1-t)+b*t;
	assert(t >= 0 && t <= 1);
	double la = exp(a/linn);
	double lb = exp(b/linn);
	double lr = la*(1-t) + lb*t;
	double r = linn*log(lr);
	//if(!((a<=r && r<=b) || (b<=r && r<=a)))
	//{
	//	fprintf(stderr, "a b t = %g %g %g\n", a, b, t);
	//	fprintf(stderr, "la lb lr = %g %g %g\n", la, lb, lr);
	//	fprintf(stderr, "r = %g\n", r);
	//}
	//if (r <= 0) r = 0;
	//if (r >= 255) r = 255;
	return r;
}

// argument for the "traverse_segment_aa" function
static void put_pixel_color_aa(int a, int b, float f, void *pp)
{
	struct {int w, h; struct rgba_value *x, c;} *p = pp;
	if (inner_point(p->w, p->h, a, b))
	{
		struct rgba_value *g = p->x + b*p->w + a;
		struct rgba_value *k = &p->c;
		g->r = lincombin(g->r, k->r, f);
		g->g = lincombin(g->g, k->g, f);
		g->b = lincombin(g->b, k->b, f);
		//g->r = g->r*(1-f) + k->r*f;
		//g->g = g->g*(1-f) + k->g*f;
		//g->b = g->b*(1-f) + k->b*f;
	}
}

// argument for the "traverse_segment_aa" function
static void put_pixel_gray_aa(int a, int b, float f, void *pp)
{
	struct {int w, h; uint8_t *x, c;} *p = pp;
	if (inner_point(p->w, p->h, a, b))
	{
		uint8_t *g = p->x + b*p->w + a;
		uint8_t *k = &p->c;
		*g = lincombin(*g, *k, f);
	}
}

// draw a color segment over a color image
static void overlay_segment_color(struct rgba_value *x, int w, int h,
		int px, int py, int qx, int qy, struct rgba_value c)
{
	//struct {int w, h; struct rgba_value *x, c; } e = {w, h, x, c};
	struct {int w, h; struct rgba_value *x, c; } e;
	e.w = w;
	e.h = h;
	e.x = x;
	e.c = c;
	traverse_segment_aa(px, py, qx, qy, put_pixel_color_aa, &e);
}

// draw a color segment over a color image
static void overlay_segment_gray(uint8_t *x, int w, int h,
		int px, int py, int qx, int qy)
{
	//struct {int w, h; struct rgba_value *x, c; } e = {w, h, x, c};
	struct {int w, h; uint8_t *x, c; } e;
	e.w = w;
	e.h = h;
	e.x = x;
	e.c = 255;
	traverse_segment_aa(px, py, qx, qy, put_pixel_gray_aa, &e);
}


// draw a color line over a color image
static void overlay_line(double a, double b, double c,
		struct rgba_value *x, int w, int h, struct rgba_value k)
{
	if (b == 0) {
		int f[2] = {-c/a, 0};
		int t[2] = {-c/a, h-1};
		overlay_segment_color(x, w, h, f[0], f[1], t[0], t[1], k);
	} else {
		double alpha = -a/b;
		double beta = -c/b;
		int f[2] = {0, beta};
		int t[2] = {w-1, alpha*(w-1)+beta};
		overlay_segment_color(x, w, h, f[0], f[1], t[0], t[1], k);
	}
}


// CLI utility to view a set of pairs of points,
// the set of "left" points is transformed by the given homography
// (produces a transparent PNG image)
int main_viewpairs(int c, char *v[])
{
	if (c != 12 && c != 13) {
		fprintf(stderr, "usage:\n\t"
			"%s a b r c d s p q 1 sx sy [mask] < pairs.txt\n", *v);
		//       0  1 2 3 4 5 6 7 8 9 10 11 12
		return EXIT_FAILURE;
	}
	double A[9]; FORI(9) A[i] = atof(v[1+i]);
	int sizex = atoi(v[10]);
	int sizey = atoi(v[11]);
	int n;
	double (*p)[4] = (void*)read_ascii_doubles(stdin, &n);
	n /= 4;
	struct rgba_value (*o)[sizex] = xmalloc(sizex*sizey*4);
	bool mask=c>12, bmask[n]; FORI(n) bmask[i] = !mask;
	if (mask) { FILE *f = xfopen(v[12], "r");
		FORI(n) {int t,q=fscanf(f, "%d", &t);bmask[i]=q&&t;}
		xfclose(f);
	}
	FORI(sizex*sizey) o[0][i] = RGBA_BLACK;
	FORI(n) {
		double *pxi = p[i];
		double *pyi = p[i]+2;
		double tt[2]; projective_map(tt, A, pxi);
		int t[2] = {tt[0], tt[1]};
		int z[2] = {pyi[0], pyi[1]};
		if (inner_point(sizex, sizey, t[0], t[1]) &&
				inner_point(sizex, sizey, z[0], z[1]))
			if (!bmask[i]) {
				overlay_segment_color(*o, sizex, sizey,
						t[0], t[1], z[0], z[1],
						RGBA_PHANTOM);
				o[t[1]][t[0]] = RGBA_RED;
				o[z[1]][z[0]] = RGBA_BLUE;
			}
	}
	FORI(n) {
		double tt[2]; projective_map(tt, A, p[i]);
		int t[2] = {tt[0], tt[1]};
		int z[2] = {p[i][2], p[i][3]};
		if (inner_point(sizex, sizey, t[0], t[1]) &&
				inner_point(sizex, sizey, z[0], z[1])
		   ) {
			if (bmask[i]) {
				overlay_segment_color(*o, sizex, sizey,
						t[0], t[1], z[0], z[1],
						RGBA_BRIGHT);
				o[t[1]][t[0]] = RGBA_RED;
				o[z[1]][z[0]] = RGBA_GREEN;
			}
		}
	}
	if (true) { // compute and show statistics
		int n_inliers = 0, n_outliers = 0;
		//stats: min,max,avg
		double instats[3]  = {INFINITY, -INFINITY, 0};
		double outstats[3] = {INFINITY, -INFINITY, 0};
		double allstats[3] = {INFINITY, -INFINITY, 0};
		FORI(n) {
			double tt[2]; projective_map(tt, A, p[i]);
			double e = hypot(tt[0]-p[i][2], tt[1]-p[i][3]);
			allstats[0] = fmin(allstats[0], e);
			allstats[1] = fmax(allstats[1], e);
			allstats[2] += e;
			if (bmask[i]) {
				instats[0] = fmin(instats[0], e);
				instats[1] = fmax(instats[1], e);
				instats[2] += e;
				n_inliers += 1;
			} else {
				outstats[0] = fmin(outstats[0], e);
				outstats[1] = fmax(outstats[1], e);
				outstats[2] += e;
				n_outliers += 1;
			}
		}
		assert(n == n_inliers + n_outliers);
		allstats[2] /= n;
		instats[2] /= n_inliers;
		outstats[2] /= n_outliers;
		fprintf(stderr, "got %d points%c", n, mask?':':'\n');
		if (mask)
			fprintf(stderr, " %d inliers (%g%%) and "
					"%d outliers (%g%%)\n",
				       	n_inliers, n_inliers*100.0/n,
					n_outliers, n_outliers*100.0/n);
		fprintf(stderr, "errors: min=%g max=%g avg=%g\n",
				allstats[0], allstats[1], allstats[2]);
		if (mask & !identityP(A)) {
			fprintf(stderr, "errors (inliers): "
					"min=%g max=%g avg=%g\n",
					instats[0], instats[1], instats[2]);
			fprintf(stderr, "errors (outliers): "
					"min=%g max=%g avg=%g\n",
					outstats[0], outstats[1], outstats[2]);
			if (outstats[0] < instats[2])
				fprintf(stderr, "WARNING: there are outliers "
					"with small error.\nProbably you are "
				"looking at the wrong model or mask file.\n");
		}
	}
	iio_write_image_uint8_vec("-", (uint8_t*)o, sizex, sizey, 4);
	return EXIT_SUCCESS;
}

// CLI utility to view a set of segments,
// (produces a binary image)
int main_viewsegs(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t"
			"%s w h < 4col.txt\n", *v);
		//       0  1 2
		return 1;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);

	int n;
	double (*p)[4] = (void*)read_ascii_doubles(stdin, &n);
	n /= 4;

	if (w == -1) for (int i = 0; i < n; i++) {
			if (w < lrint(p[i][0])) w = lrint(p[i][0]);
			if (w < lrint(p[i][2])) w = lrint(p[i][2]);
		}
	if (h == -1) for (int i = 0; i < n; i++) {
			if (h < lrint(p[i][1])) h = lrint(p[i][1]);
			if (h < lrint(p[i][3])) h = lrint(p[i][3]);
		}

	uint8_t (*o)[w] = xmalloc(w*h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
		o[j][i] = 0;
	for (int i = 0; i < n; i++)
	{
		int a[2] = { lrint(p[i][0]), lrint(p[i][1]) };
		int b[2] = { lrint(p[i][2]), lrint(p[i][3]) };
		if (inner_point(w,h, a[0],a[1]) && inner_point(w,h, b[0],b[1]))
			overlay_segment_gray(*o, w, h, a[0], a[1], b[0], b[1]);
	}

	iio_write_image_uint8_vec("-", (uint8_t*)o, w, h, 1);
	return 0;
}

// read stream until character "stop" is found
// if EOF is reached, return NULL
// otherwise return line
static char *fgets_until(char *line, int n, FILE *f, int stop)
{
	int i = 0;
	while(1) {
		if (i >= n-1) break;
		int c = fgetc(f);
		if (c == EOF) return NULL;
		line[i] = c;
		if (c == stop) break;
		i += 1;
	}
	line[i+1] = '\0';
	//fprintf(stderr, "FGETS UNTIL %d \"%s\"\n", i, line);
	return line;
}

// CLI utility to view a set of polygonal curves
// boolean option "-c" : close the polygon
int main_viewpolygons(int c, char *v[])
{
	bool option_c = pick_option(&c, &v, "c", NULL);
	if (c != 3 && c != 3) {
		fprintf(stderr, "usage:\n\t"
			"%s w h [-c] < polygons.txt\n", *v);
		//       0  1 2
		return EXIT_FAILURE;
	}
	int w = atoi(v[1]);
	int h = atoi(v[2]);
	struct rgba_value (*o)[w] = xmalloc(w * h * 4);
	for (int i = 0; i < w*h; i++)
		o[0][i] = RGBA_BLACK;
	while (1) {
		int n, maxlin = 1000*12*4;
		char line[maxlin], *sl = fgets_until(line, maxlin, stdin, '\n');
		if (!sl) break;
		double *t = alloc_parse_doubles(maxlin, line, &n);
		n /= 2;
		for (int i = 0; i < n - 1 + option_c; i++)
			overlay_segment_color(*o, w, h, t[2*i+0], t[2*i+1],
				t[(2*i+2)%(2*n)], t[(2*i+3)%(2*n)], RGBA_GREEN);
		free(t);
	}
	iio_write_image_uint8_vec("-", (uint8_t*)o, w, h, 4);
	return EXIT_SUCCESS;
}

// CLI utility to display a set of triplets of planar points
// (produces a transparent PNG image)
int main_viewtrips(int c, char *v[])
{
	if (c != 3 && c != 4) {
		fprintf(stderr, "usage:\n\t%s sx sy [mask] < pairs.txt\n", *v);
		//                         0  1  2   3
		return EXIT_FAILURE;
	}
	int s[2], n; FORI(2) s[i] = atoi(v[1+i]);
	double (*p)[6] = (void*)read_ascii_doubles(stdin, &n);
	n /= 6;
	struct rgba_value (*o)[s[0]] = xmalloc(s[0]*s[1]*4);
	bool bmask[n]; FORI(n) bmask[i] = c<4;
	if (c>3) { FILE *f = xfopen(v[3], "r");
		FORI(n) {int t,q=fscanf(f, "%d", &t);bmask[i]=q&&t;}
		xfclose(f);
	}
	FORI(s[0]*s[1]) o[0][i] = RGBA_BLACK;
	FORI(n) {
		int t[2] = {p[i][0], p[i][1]};
		int z[2] = {p[i][2], p[i][3]};
		int Z[2] = {p[i][4], p[i][5]};
		if (inner_point(s[0], s[1], t[0], t[1]) &&
				inner_point(s[0], s[1], z[0], z[1]) &&
				inner_point(s[0], s[1], Z[0], Z[1]))
		{
			//void (*pix)(int,int,void*) = (bmask[i]&&c>3)?
			//	draw_brighter_pixel:draw_phantom_pixel;
			struct rgba_value kk = (bmask[i]&&c>3)?
				RGBA_BRIGHT:RGBA_PHANTOM;
			overlay_segment_color(*o, s[0], s[1],
					t[0], t[1], z[0], z[1], kk);
			overlay_segment_color(*o, s[0], s[1],
					z[0], z[1], Z[0], Z[1], kk);
		}
	}
	FORI(n) {
		int t[2] = {p[i][0], p[i][1]};
		int z[2] = {p[i][2], p[i][3]};
		int Z[2] = {p[i][4], p[i][5]};
		if (inner_point(s[0], s[1], t[0], t[1]) &&
				inner_point(s[0], s[1], z[0], z[1]) &&
				inner_point(s[0], s[1], Z[0], Z[1]))
		{
			o[t[1]][t[0]] = RGBA_YELLOW;
			o[z[1]][z[0]] = RGBA_GREEN;
			o[Z[1]][Z[0]] = RGBA_MAGENTA;
		}
	}
	iio_write_image_uint8_vec("-", (uint8_t*)o, s[0], s[1], 4);
	return EXIT_SUCCESS;
}


// L := x' * F'
// compute the F-epipolar line defined by the point x
static void epipolar_line(double L[3], double F[9], double x[2])
{
	L[0] = x[0]*F[0] + x[1]*F[3] + F[6];
	L[1] = x[0]*F[1] + x[1]*F[4] + F[7];
	L[2] = x[0]*F[2] + x[1]*F[5] + F[8];
}

// project point x to the line L
static void project_point_to_line(double px[2], double L[3], double x[2])
{
	double a = L[0];
	double b = L[1];
	double c = L[2];
	double d = a*a + b*b;
	double e = a*x[1] - b*x[0];
	px[0] = (-a*c - b*e)/d;
	px[1] = (-b*c + a*e)/d;
}

static int innpair(int s[2], double p[4])
{
	return ( 1
			&& p[0] >= 0 && p[0] < s[0]
			&& p[1] >= 0 && p[1] < s[1]
			&& p[2] >= 0 && p[2] < s[0]
			&& p[3] >= 0 && p[3] < s[1]
	       );
}

// CLI utility to display a fundamental matrix and some matched points
// (produces a transparent PNG image)
int main_viewepi(int c, char *v[])
{
	if (c != 12 && c != 13) {
		fprintf(stderr, "usage:\n\t%s f0 ... f8 sx sy [mask] < pairs\n", *v);
		//                          0 1      9  10 11  12
		return EXIT_FAILURE;
	}
	double A[9]; FORI(9) A[i] = atof(v[1+i]);
	int s[2], n; FORI(2) s[i] = atoi(v[10+i]);

	// 0. read input pairs
	double (*p)[4] = (void*)read_ascii_doubles(stdin, &n);
	n /= 4;


	// 1. create output image and plot black background
	struct rgba_value (*o)[s[0]] = xmalloc(s[0]*s[1]*4);
	FORI(s[0]*s[1]) o[0][i] = RGBA_BLACK;

	// 2. if there is a mask file, read it
	bool mask=c>12, bmask[n]; FORI(n) bmask[i] = !mask;
	if (mask) { FILE *f = xfopen(v[12], "r");
		FORI(n) {int t,q=fscanf(f, "%d", &t);bmask[i]=q&&t;}
		xfclose(f);
	}


	// for each point p (=px) in image A
	// 3. draw the epipolar line L of p
	FORI(n) if (!bmask[i] && innpair(s,p[i])) {
			double L[3]; epipolar_line(L, A, p[i]);
			overlay_line(L[0],L[1],L[2],o[0],s[0],s[1],RGBA_GRAY10);
	}
	FORI(n) if (bmask[i] && innpair(s,p[i])) {
		double L[3]; epipolar_line(L, A, p[i]);
			overlay_line(L[0],L[1],L[2],o[0],s[0],s[1],RGBA_GRAY50);
	}
	// 4. draw the projection line of q to L
	FORI(n) {
		if (!innpair(s,p[i])) continue;
		double L[3];
		epipolar_line(L, A, p[i]);
		double Ly[2];
		int y[2] = {lrint(p[i][2]), lrint(p[i][3])};
		project_point_to_line(Ly, L, p[i]+2);
		int iLy[2] = {lrint(Ly[0]), lrint(Ly[1])};
		if (inner_point(s[0], s[1], iLy[0], iLy[1]))
		{
			if (inner_point(s[0], s[1], y[0], y[1]))
				overlay_segment_color(*o, s[0], s[1],
						y[0], y[1], iLy[0], iLy[1],
					bmask[i]?RGBA_BRIGHT:RGBA_PHANTOM);
			if (bmask[i])
				o[iLy[1]][iLy[0]] = RGBA_MAGENTA;
			else
				o[iLy[1]][iLy[0]] = RGBA_RED;
		}
	}
	// 5. draw the corresponding point q (py)
	FORI(n) {
		if (!innpair(s,p[i])) continue;
		int y[2] = {lrint(p[i][2]), lrint(p[i][3])};
		if (inner_point(s[0], s[1], y[0], y[1])) {
			if (bmask[i])
				o[y[1]][y[0]] = RGBA_GREEN;
			else
				o[y[1]][y[0]] = RGBA_BLUE;
		}
	}
	iio_write_image_uint8_vec("-", (uint8_t*)o, s[0], s[1], 4);

	if (true) { // show statistics
		int n_inliers = 0, n_outliers = 0;
		//stats: min,max,avg
		double instats[3]  = {INFINITY, -INFINITY, 0};
		double outstats[3] = {INFINITY, -INFINITY, 0};
		double allstats[3] = {INFINITY, -INFINITY, 0};
		FORI(n) {
			double L[3], Ly[2];
			epipolar_line(L, A, p[i]);
			project_point_to_line(Ly, L, p[i]+2);
			double e = hypot(Ly[0]-p[i][2], Ly[1]-p[i][3]);
			allstats[0] = fmin(allstats[0], e);
			allstats[1] = fmax(allstats[1], e);
			allstats[2] += e;
			if (bmask[i]) {
				instats[0] = fmin(instats[0], e);
				instats[1] = fmax(instats[1], e);
				instats[2] += e;
				n_inliers += 1;
			} else {
				outstats[0] = fmin(outstats[0], e);
				outstats[1] = fmax(outstats[1], e);
				outstats[2] += e;
				n_outliers += 1;
			}
		}
		assert(n == n_inliers + n_outliers);
		allstats[2] /= n;
		instats[2] /= n_inliers;
		outstats[2] /= n_outliers;
		fprintf(stderr, "got %d points%c", n, mask?':':'\n');
		if (mask)
			fprintf(stderr, " %d inliers (%g%%) and "
					"%d outliers (%g%%)\n",
				       	n_inliers, n_inliers*100.0/n,
					n_outliers, n_outliers*100.0/n);
		fprintf(stderr, "errors: min=%g max=%g avg=%g\n",
				allstats[0], allstats[1], allstats[2]);
		if (mask && !identityP(A)) {
			fprintf(stderr, "errors (inliers): "
					"min=%g max=%g avg=%g\n",
					instats[0], instats[1], instats[2]);
			fprintf(stderr, "errors (outliers): "
					"min=%g max=%g avg=%g\n",
					outstats[0], outstats[1], outstats[2]);
			if (outstats[0] < instats[2])
				fprintf(stderr, "WARNING: there are outliers "
					"with small error.\nProbably you are "
				"looking at the wrong model or mask file.\n");
		}
	}
	return EXIT_SUCCESS;
}




static int plot_line(int (*P)[2], int w, int h, double a, double b, double c)
{
	int r = 0;
	if (fabs(a) < fabs(b)) { // slope less than 1
		double p = -a/b;
		double q = -c/b;
		for (int x = 0; x < w; x++) {
			int y = round(p*x + q);
			P[r][0] = x;
			if (y >= 0 && y < h)
				P[r++][1] = y;
		}
	} else {
		double p = -b/a;
		double q = -c/a;
		for (int y = 0; y < h; y++) {
			int x = round(p*y + q);
			P[r][1] = y;
			if (x >= 0 && x < w)
				P[r++][0] = x;
		}
	}
	return r;
}

static int plot_epipolar(int (*p)[2], double fm[9], int w, int h, int i, int j)
{
	double a = i*fm[0] + j*fm[3] + fm[6];
	double b = i*fm[1] + j*fm[4] + fm[7];
	double c = i*fm[2] + j*fm[5] + fm[8];
	return plot_line(p, w, h, a, b, c);
}

#include "random.c"

static int randombounds(int a, int b)
{
	if (b < a)
		fail("the interval [%d, %d] is empty!", a, b);
	if (b == a)
		return b;
	return a + lcg_knuth_rand()%(b - a + 1);
}

// CLI utility to display the epipolar lines defined by a fundamental matrix
int main_viewfmpair(int c, char *v[])
{
	if (c != 13) {
		fprintf(stderr, "usage:\n\t%s f0 ... f8 sx sy [0|1]\n", *v);
		//                          0 1      9  10 11 12
		return EXIT_FAILURE;
	}
	double A[9]; FORI(9) A[i] = atof(v[1+i]);
	int s[2], n; FORI(2) s[i] = atoi(v[10+i]);
	int overlay = atoi(v[12]);

	double Atr[9];
	FORI(3)FORJ(3) Atr[3*i+j]=A[3*j+i];

	struct rgba_value (*o)[2*s[0]] = xmalloc(2*s[0]*s[1]*4);
	FORI(2*s[0]*s[1]) o[0][i] = RGBA_BLACK;

	int maxlpoints = 4 * (s[0]*s[1]);
	int (*p)[2] = xmalloc(maxlpoints*sizeof*p);

	for (int cx = 0; cx < 20; cx += 1)
	{
		int i = randombounds(0, s[0]-1);
		int j = randombounds(0, s[1]-1);
		int np = plot_epipolar(p, A, s[0], s[1], i, j);
		for (int pi = 0; pi < np; pi++)
		{
			int ii = p[pi][0];
			int jj = p[pi][1];
			//fprintf(stderr, "p[%d / %d] = %d %d\n", pi, np, ii, jj);
			assert(ii >= 0);
			assert(ii < s[0]);
			assert(jj >= 0);
			assert(jj < s[1]);
			o[jj][ii] = RGBA_BRIGHT;
		}
		o[j][i] = RGBA_GREEN;

		np = plot_epipolar(p, Atr, s[0], s[1], i, j);
		for (int pi = 0; pi < np; pi++)
			o[p[pi][1]][p[pi][0]+s[0]] = RGBA_BRIGHT;
		o[j][i+s[0]] = RGBA_GREEN;
	}


	iio_write_image_uint8_vec("-", (uint8_t*)o, 2*s[0], s[1], 4);
	return EXIT_SUCCESS;
}

//#define BAD_MAX(x,y) (((x)>(y))?(x):(y))
//
//// CLI utility to display the epipolar lines defined by a fundamental matrix
//int main_viewfmpairi(int c, char *v[])
//{
//	float dimfactor = 0.5;
//
//	if (c != 13) {
//		fprintf(stderr, "usage:\n\t"
//			"%s f0 ... f8 a.png b.png <points >view\n", *v);
//		//        0 1      9  10    11
//		return EXIT_FAILURE;
//	}
//	double A[9]; FORI(9) A[i] = atof(v[1+i]);
//	char *filename_a = v[10];
//	char *filename_b = v[11];
//
//	int wa, ha, wb, hb, pda, pdb;
//	float *a = iio_load_image_uint8_vec(filename_a, &wa, &ha, &pda);
//	float *b = iio_load_image_uint8_vec(filename_b, &wb, &hb, &pdb);
//	int w = wa+wb;
//	int h = BAD_MAX(ha, hb);
//	pd = pda;
//	if (pda != 3 || pdb != 3) fail("works with color images");
//
//	assert(3 == sizeof(struct rgb_value));
//	uint8_t *out = xmalloc(w*h*pd);
//	struct rgb_value (*o)[w] = (void*)out;
//
//	// set background
//	for (int i = 0; i < w*h*pd; i++)
//		out[i] = 0;
//	for (int j = 0; j < ha; j++)
//	for (int i = 0; i < wa; i++)
//		o[j][i] = *(struct rgb_value*)(a+3*(w*j+i));
//	for (int j = 0; j < hb; j++)
//	for (int i = 0; i < wb; i++)
//		o[j][i+wa] = *(struct rgb_value*)(b+3*(w*j+i));
//	for (int i = 0; i < w*h*pd; i++)
//		out[i] *= dimfactor;
//
//	double Atr[9];
//	FORI(3)FORJ(3) Atr[3*i+j]=A[3*j+i];
//
//	int n;
//	float *t = read_ascii_floats(stdin, &n);
//	n /= 4;
//
//	int maxlpoints = 4 * w*h;
//	int (*p)[2] = xmalloc(maxlpoints*sizeof*p);
//
//	for (int cx = 0; cx < 20; cx += 1)
//	{
//		int i = randombounds(0, s[0]-1);
//		int j = randombounds(0, s[1]-1);
//		int np = plot_epipolar(p, A, s[0], s[1], i, j);
//		for (int pi = 0; pi < np; pi++)
//		{
//			int ii = p[pi][0];
//			int jj = p[pi][1];
//			//fprintf(stderr, "p[%d / %d] = %d %d\n", pi, np, ii, jj);
//			assert(ii >= 0);
//			assert(ii < s[0]);
//			assert(jj >= 0);
//			assert(jj < s[1]);
//			o[jj][ii] = RGBA_BRIGHT;
//		}
//		o[j][i] = RGBA_GREEN;
//
//		np = plot_epipolar(p, Atr, s[0], s[1], i, j);
//		for (int pi = 0; pi < np; pi++)
//			o[p[pi][1]][p[pi][0]+s[0]] = RGBA_BRIGHT;
//		o[j][i+s[0]] = RGBA_GREEN;
//	}
//
//
//	iio_write_image_uint8_vec("-", (uint8_t*)o, 2*s[0], s[1], 4);
//	return EXIT_SUCCESS;
//}

// CLI utility to access some visualization programs
// all programs read text file from stdin and write a transparent PNG to stdout
int main_pview(int c, char *v[])
{
	assert(4 == sizeof(struct rgba_value));
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "points")) return main_viewp(c-1, v+1);
	else if (0 == strcmp(v[1], "hpoints")) return main_viewhp(c-1, v+1);
	else if (0 == strcmp(v[1], "pairs")) return main_viewpairs(c-1, v+1);
	else if (0 == strcmp(v[1], "segments")) return main_viewsegs(c-1, v+1);
	else if (0 == strcmp(v[1], "polygons")) return main_viewpolygons(c-1, v+1);
	else if (0 == strcmp(v[1], "triplets")) return main_viewtrips(c-1, v+1);
	else if (0 == strcmp(v[1], "epipolar")) return main_viewepi(c-1, v+1);
	else if (0 == strcmp(v[1], "fmpair")) return main_viewfmpair(c-1, v+1);
	//else if (0 == strcmp(v[1], "fmpairi")) return main_viewfmpairi(c-1,v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s [points|pairs|triplets|epipolar] "
			       "params... < data.txt | display\n", *v);
	       return EXIT_FAILURE;
	}
}

#ifndef HIDE_ALL_MAINS
int main(int c, char **v) { return main_pview(c, v); }
#endif
