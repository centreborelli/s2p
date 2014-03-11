//////////////////
/// INTERFACE ////
//////////////////


struct cloud_polygon {
	int n;                   // number of vertices
	double *v;               // vertex coordinates (array of length 2*n)
};

struct cloud_mask {
	int n;                   // number of polygons
	struct cloud_polygon *t; // array of polygons
	double low[2], up[2];    // "rectangle"
};

int read_cloud_mask_from_gml_file(struct cloud_mask *m, char *filename);

void clouds_mask_fill(int *img, int w, int h, struct cloud_mask *m);





///////////////////////
/// IMPLEMENTATION ////
///////////////////////

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"
#include "drawsegment.c"

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

static void read_until_newline(FILE *f)
{
	while (1) {
		int c = fgetc(f);
		if (c == EOF || c == '\n')
			return;
	}
}

// like strcmp, but finds a needle
static int strhas(char *haystack, char *needle)
{
	char *r = strstr(haystack, needle);
	return r ? 0 : 1;
}

static void cloud_add_polygon(struct cloud_mask *m, struct cloud_polygon p)
{
	m->n += 1;
	m->t = xrealloc(m->t, m->n * sizeof*m->t);
	m->t[m->n-1] = p;
}

// needed only for checking consistency at the end
void free_cloud(struct cloud_mask *m)
{
	for (int i = 0; i < m->n; i++)
		free(m->t[i].v);
	free(m->t);
}

int read_cloud_mask_from_gml_file(struct cloud_mask *m, char *filename)
{
	m->n = 0;
	m->t = NULL;

	FILE *f = xfopen(filename, "r");
	int n = 0x100, nf;
	while(1) {
		char line[n], *sl = fgets_until(line, n, f, '>');
		if (!sl) break;
		if (0 == strhas(line, "lowerCorner")) {
			double *ff = read_ascii_doubles(f, &nf);
			if (nf == 2) for (int i = 0; i < 2; i++)
				m->low[i] = ff[i];
			free(ff);
		}
		if (0 == strhas(line, "upperCorner")) {
			double *ff = read_ascii_doubles(f, &nf);
			if (nf == 2) for (int i = 0; i < 2; i++)
				m->up[i] = ff[i];
			free(ff);
		}
		if (0 == strhas(line, "posList")) {
			struct cloud_polygon p;
			p.v = read_ascii_doubles(f, &p.n);
			p.n /= 2;
			cloud_add_polygon(m, p);
		}
		read_until_newline(f);
	}
	xfclose(f);
	return 0;
}

static void putpixel_0(int *img, int w, int h, float x, float y, int v)
{
	int i = round(x);
	int j = round(y);
	if (i >= 0 && j >= 0 && i < w && j < h)
		img[j*w+i] = v;
}

struct plot_data { int *x; int w, h; int c; };
static void plot_pixel(int i, int j, void *e)
{
	struct plot_data *d = e;
	putpixel_0(d->x, d->w, d->h, i, j, d->c);
}

static void plot_segment_gray(int *img, int w, int h,
		float a[2], float b[2], int c)
{
	int p[2] = {round(a[0]), round(a[1])};
	int q[2] = {round(b[0]), round(b[1])};
	struct plot_data d = { img, w, h, c };
	traverse_segment(p[0], p[1], q[0], q[1], plot_pixel, &d);
}

static int dsf_find(int *t, int a)
{
	if (a != t[a])
		t[a] = dsf_find(t, t[a]);
	return t[a];
}

static int dsf_make_link(int *t, int a, int b)
{
	if (a < b) { // arbitrary choice
		t[b] = a;
		return a;
	} else {
		t[a] = b;
		return b;
	}
}

static int dsf_join(int *t, int a, int b)
{
	a = dsf_find(t, a);
	b = dsf_find(t, b);
	if (a != b)
		b = dsf_make_link(t, a, b);
	return b;
}

// connected components of positive pixels of the image rep
static void positive_connected_component_filter(int *rep, int w, int h)
{
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			rep[i] = i;
	for (int j = 0; j < h - 1; j++)
	for (int i = 0; i < w - 1; i++)
	{
		int p0 = j*w + i;
		int p1 = j*w + i+1;
		int p2  = (j+1)*w + i;
		if (rep[p0] >= 0 && rep[p1] >= 0)
			dsf_join(rep, p0, p1);
		if (rep[p0] >= 0 && rep[p2] >= 0)
			dsf_join(rep, p0, p2);

	}
	for (int i = 0; i < w*h; i++)
		if (rep[i] >= 0)
			rep[i] = dsf_find(rep, i);
}

// area of triangle ABC
static double triangle_area(double *A, double *B, double *C)
{
	double X[2] = {B[0] - A[0], B[1] - A[1]};
	double Y[2] = {C[0] - A[0], C[1] - A[1]};
	return X[0]*Y[1] - X[1]*Y[0];
}

// test wether point X is inside triangle ABC
static int winding_triangle(double *A, double *B, double *C, double *X)
{
	double v1 = triangle_area(A, B, X);
	double v2 = triangle_area(B, C, X);
	double v3 = triangle_area(C, A, X);
	int r = 0;
	if (v1 >= 0 && v2 >= 0 && v3 >= 0) r = 1;
	if (v1 < 0 && v2 < 0 && v3 < 0) r = -1;
	return r;
}

// winding number of point (x,y) with respect to a single polygon
static int winding_number_polygon(struct cloud_polygon *p, int x, int y)
{
	int r = 0;
	for (int i = 1; i < p->n - 2; i++)
	{
		double *A = p->v;
		double *B = p->v + 2*(i);
		double *C = p->v + 2*(i+1);
		double X[2] = {x, y};
		r += abs(winding_triangle(A, B, C, X));
	}
	return r;
}

// winding number of point (x,w) with respect to all polygons of the mask
static int winding_number_clouds(struct cloud_mask *m, int x, int y)
{
	int r = 0;
	for (int i = 0; i < m->n; i++)
		r += winding_number_polygon(m->t + i, x, y);
	return r;
}

// rescale a cloud of points to fit in the given rectangle
static void cloud_mask_rescale(struct cloud_mask *m, int w, int h)
{
	for (int i = 0; i < m->n; i++)
	{
		struct cloud_polygon *p = m->t + i;
		for (int j = 0; j < p->n; j++)
		{
			double *a = p->v + 2*j;
			double A[2] = {
				(a[0] - m->low[0]) / (m->up[0] - m->low[0]) * w,
				(a[1] - m->low[1]) / (m->up[1] - m->low[1]) * h
			};
			a[0] = A[0];
			a[1] = A[1];
		}
	}
}

// homographic transform y=H(x) 
static void apply_homography(double y[2], double H[9], double x[2])
{
	double z[3];
	z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
	z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
	z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = z[0]/z[2];
	y[1] = z[1]/z[2];
}

// transform the coordinates of a cloud_mask by the given homography
static void cloud_mask_homography(struct cloud_mask *m, double *H)
{
	for (int i = 0; i < m->n; i++)
		for (int j = 0; j < m->t[i].n; j++)
			apply_homography(2*j+m->t[i].v, H, 2*j+m->t[i].v);
}

#include "iio.h"
void clouds_mask_fill(int *img, int w, int h, struct cloud_mask *m)
{
	// initialize the image to 0
	for (int i = 0; i < w*h; i++)
		img[i] = 0;

	// plot the pixels at the cloud edges with the value -1
	for (int i = 0; i < m->n; i++)
	{
		struct cloud_polygon *p = m->t + i;
		for (int j = 0; j < p->n - 1; j++)
		{
			float a[2] = {p->v[2*j+0], p->v[2*j+1]};
			float b[2] = {p->v[2*j+2], p->v[2*j+3]};
			plot_segment_gray(img, w, h, a, b, -1);
			putpixel_0(img, w, h, a[0], a[1], -2);
			putpixel_0(img, w, h, b[0], b[1], -2);
		}
	}

	// identify the connected components of positive values
	positive_connected_component_filter(img, w, h);


	// identify the connected components that are background
	//int background_ids[m->n], n_background_ids = 0;
	int background_ids[(int) fmax(1000, m->n)], n_background_ids = 0;
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = j*w + i;
		if (img[idx] == idx)
		{
			int r = winding_number_clouds(m, i, j);
			if (!r) {
				background_ids[n_background_ids] = idx;
				n_background_ids += 1;
//				if (n_background_ids >= m->n)
//					fail("bad heuristic of background"
//						       " components %d %d",
//						       n_background_ids, m->n
//						       );
			}
		}
	}

	// paint the background pixels in black, the clouds in white
	for (int i = 0; i < w*h; i++)
	{
		bool bP = false;
		for (int j = 0; j < n_background_ids; j++)
			if (img[i] == background_ids[j])
				bP = true;
		//int j = 0;
		//while (j < n_background_ids)
		//	if (img[i] == background_ids[j++])
		//		break;
		img[i] = bP ? 0 : 255;
	}
}

#define CLDMASK_MAIN
#ifdef CLDMASK_MAIN

#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	// read input arguments
	char *Hstring = pick_option(&c, &v, "h", "");
	if (c != 5 && c!= 4 && c != 3) {
		return fprintf(stderr, "usage:\n\t%s "
		"width height [-h \"h1 ... h9\"] [clouds.gml [out.png]]\n", *v);
		//   1 2                          3           4
	}
	int out_width = atoi(v[1]);
	int out_height = atoi(v[2]);
	char *filename_clg = c > 3 ? v[3] : "-";
	char *filename_out = c > 4 ? v[4] : "PNG:-";

	// read input cloud file
	struct cloud_mask m[1];
	read_cloud_mask_from_gml_file(m, filename_clg);

	// acquire space for output image
	int w = out_width;
	int h = out_height;
	int *x = xmalloc(w*h*sizeof*x);
	for (int i = 0; i < w*h; i++)
		x[i] = 0;

	// scale the co-ordinates of the cloud
	if (*Hstring) {
		int nH;
		double *H = alloc_parse_doubles(9, Hstring, &nH);
		if (nH != 9)
			fail("can not read 3x3 matrix from \"%s\"", Hstring);
		cloud_mask_homography(m, H);
		free(H);

	} else
		cloud_mask_rescale(m, w, h);

	// draw mask over output image
	clouds_mask_fill(x, w, h, m);

	// save output image
	iio_save_image_int(filename_out, x, w, h);

	//cleanup
	free(x);
	free_cloud(m);
	return 0;
}

#endif//CLDMASK_MAIN
