#include <assert.h>
#include <stdio.h>
#include <math.h>


#include "fail.c"
#include "xmalloc.c"
#include "xfopen.c"
#include "parsenumbers.c"

#include "siftie.h"

#include "smapa.h"

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)



#if USE_KDTREE
#include <kdtree.h>
#endif

void write_raw_sift(FILE *f, struct sift_keypoint *k)
{
	fprintf(f, "%g %g %g %g", k->pos[0], k->pos[1], k->scale,
			k->orientation);
	FORJ(SIFT_LENGTH)
		fprintf(f, " %g", k->sift[j]);
	fprintf(f, "\n");
}

//void write_raw_sift_bin(FILE *f, struct sift_keypoint *k)
//{
//	assert(SIFT_LENGTH==128);
//	size_t n = 4*sizeof(float)+128;
//	char buf[n];
//	float *p = (float*)buf;
//	p[0] = k->pos[0];
//	p[1] = k->pos[1];
//	p[2] = k->scale;
//	p[3] = k->orientation;
//	char *pp = (char *)(4+p);
//	FORI(128) pp[i] = k->sift[i];
//	size_t r = fwrite(buf, n, 1, f);
//	if (r != n) fail("could not write some sift descriptor");
//}
//
//SMART_PARAMETER_SILENT(SIFT_OUTBIN,0)
//SMART_PARAMETER_SILENT(SIFT_INBIN,0)

void write_raw_sifts(FILE *f, struct sift_keypoint *k, int n)
{
	FORI(n) {
//		if (BIN_OSIFT() > 0)
//			write_raw_sift_bin(f, k);
//		else
			write_raw_sift(f, k);
		k += 1;
	}
}

struct sift_keypoint *read_raw_sifts(FILE *f, int *no)
{
	int n, N = 132;
	float *t = read_ascii_floats(f, &n);
	if (n == 0) {*no=0; return NULL;}
	if (0 != n % N) fail("bad raw SIFT keypoints format (%d %d)", n, n%N);
	n /= N;
	struct sift_keypoint *r = xmalloc(n * sizeof * r);
	FORI(n) {
		int off = i*N;
		r->pos[0] = t[0+off];
		r->pos[1] = t[1+off];
		r->scale = t[2+off];
		r->orientation = t[3+off];
		FORJ(SIFT_LENGTH)
		{
			int c = t[4+j+off];
			if (c < 0 || c > 255) fail("bad sift value %d\n", c);
			r->sift[j] = c;
		}
		r += 1;
	}
	xfree(t);
	*no = n;
	return r-n;
}

void write_raw_siftb(FILE *f, struct sift_keypoint *k)
{
	float t[4] = {k->pos[0], k->pos[1], k->scale, k->orientation};
	size_t r = fwrite(t, sizeof(*t), 4, f);
	if (r != sizeof*t)
		fail("could not write SIFT descriptor table (r = %zu)", r);
	char T[SIFT_LENGTH];
	FORI(SIFT_LENGTH)
		if (k->sift[i] < 0 || k->sift[i] > 255)
			fail("bad sift descriptor entry %g", k->sift[i]);
		else
			T[i] = k->sift[i];
	r = fwrite(T, SIFT_LENGTH, 1, f);
	if (r != 1)
		fail("could not write SIFT descriptor table (r=%zu)", r);
}

void write_raw_siftsb(FILE *f, struct sift_keypoint *k, int n)
{
	FORI(n)
		write_raw_siftb(f, k+i);
}

static void *freadwhole_f(FILE *f, long *on)
{
#if 0
	int r = fseek(f, 0, SEEK_END);
	if (r)
		error("can not determine size of file (%d)", r);
	long n = ftell(f);
	if (n < 0)
		error("can not determine size of file (%ld)", n);
	void *ret = xmalloc(n);
	long rr = fread(ret, 1, n, f);
	if (rr != n)
		error("could not read %ld bytes from file (only %ld)", n, rr);
	*on = n;
	return ret;
#else
	int r, n = 0, nt = 0;
	char *t = NULL;
	while(1) {
		if (n >= nt)
		{
			nt = 1000+2 * (nt + 1);
			t = xrealloc(t, nt);
		}
		r = fgetc(f);
		if (r == EOF)
			break;
		t[n] = r;
		n += 1;
	}
	*on = n;
	return t;
#endif
}

static void *freadwhole(const char *fname, long *on)
{
	FILE *f = xfopen(fname, "r");
	void *ret = freadwhole_f(f, on);
	xfclose(f);
	return ret;
}


// sift descriptor saved length
static const int SDSLEN = SIFT_LENGTH*1 + 4*sizeof(float);

struct sift_keypoint *read_raw_siftsb(FILE *f, int *no)
{
	long nn;
	void *p = freadwhole_f(f, &nn);
	long n = nn / SDSLEN;
	if (n*SDSLEN != nn)
		fail("can not read binary sift file (%ld %ld)", n*SDSLEN, nn);
	struct sift_keypoint *ret = xmalloc(n * sizeof * ret);
	FORI(n) {
		struct sift_keypoint *k = ret + i;
		void *pi = SDSLEN * i + (char *)p;
		float *t = pi;
		k->pos[0] = t[0];
		k->pos[1] = t[1];
		k->scale = t[2];
		k->orientation = t[3];
		unsigned char *tt = 4*sizeof(float) + (unsigned char *)pi;
		FORJ(SIFT_LENGTH)
			k->sift[j] = tt[j];

	}
	xfree(p);
	*no = n;
	return ret;
}

SMART_PARAMETER(SIFT_BINARY,0)

struct sift_keypoint *read_raw_sifts_gen(FILE *f, int *no)
{
	return SIFT_BINARY() ? read_raw_siftsb(f, no) : read_raw_sifts(f, no);
}
void write_raw_sifts_gen(FILE *f, struct sift_keypoint *k, int n)
{
	SIFT_BINARY() ? write_raw_siftsb(f, k, n) : write_raw_sifts(f, k, n);
}
void write_raw_sift_gen(FILE *f, struct sift_keypoint *k)
{
	SIFT_BINARY() ? write_raw_siftb(f, k) : write_raw_sift(f, k);
}

static float emvdistppf(float *a, float *b, int n)
{
	//float ac = a[0];
	//float bc = b[0];
	//float r = fabs(ac - bc);
	//FORI(n-1)
	//	r += fabs(ac + a[i+1] - bc - b[i+1]);
	//return r;
	float ac[n]; ac[0] = a[0]; FORI(n-1) ac[i+1] = ac[i] + a[i+1];
	float bc[n]; bc[0] = b[0]; FORI(n-1) bc[i+1] = bc[i] + b[i+1];
	float r = 0;
	FORI(n) r += fabs(ac[i] - bc[i]);
	fprintf(stderr, "%g\n", r);
	return r;
}

static double sqr(double x)
{
	return x*x;
}

static float distppft(float *a, float *b, int n, float t)
{
	double tt = t * t;
	double r = 0;
	for (int i = 0; i < n; i++)
		if (r > tt)
			return t + 1;
		else
			r = hypot(r, b[i] - a[i]);
			//r += sqr(b[i] - a[i]);
	return r;//sqrt(r);
}

static float distppf(float *a, float *b, int n)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += sqr(b[i] - a[i]);
	return sqrt(r);
}

static float distlpf(float *a, float *b, int n, float p)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		r += pow(fabs(b[i] - a[i]), p);
	return pow(r, 1.0/p);
}

SMART_PARAMETER(DIST_DESCS_EMV,0)
SMART_PARAMETER(DIST_DESCS_LP,2)

static double dist_descs(struct sift_keypoint *a, struct sift_keypoint *b)
{
	//float af[SIFT_LENGTH]; FORI(SIFT_LENGTH) af[i] = a->sift[i];
	//float bf[SIFT_LENGTH]; FORI(SIFT_LENGTH) bf[i] = b->sift[i];
	//return distppf(a->sift, b->sift, SIFT_LENGTH);
	if (DIST_DESCS_EMV() > 0.5)
		return emvdistppf(a->sift, b->sift, SIFT_LENGTH);
	else if (DIST_DESCS_LP() != 2)
		return distlpf(a->sift, b->sift, SIFT_LENGTH, DIST_DESCS_LP());
	else
		return distppf(a->sift, b->sift, SIFT_LENGTH);
}

static double dist_descst(struct sift_keypoint *a, struct sift_keypoint *b,
		float t)
{
	if (DIST_DESCS_EMV() > 0.5) {
		float r = emvdistppf(a->sift, b->sift, SIFT_LENGTH);
		return r < t ? r : t;
	} else if (DIST_DESCS_LP() != 2) {
		float r=distlpf(a->sift, b->sift, SIFT_LENGTH,DIST_DESCS_LP());
		return r < t ? r : t;
	} else
		return distppft(a->sift, b->sift, SIFT_LENGTH, t);
}

static float euclidean_distance_topped(float *a, float *b, int n, float top)
{
	double r = 0;
	for (int i = 0; i < n; i++)
		if (r > top)
			return top + 1;
		else
			r = hypot(r, b[i] - a[i]);
	return r;
}

static double sift_distance_topped(
		struct sift_keypoint *a,
		struct sift_keypoint *b,
		float top)
{
	return euclidean_distance_topped(a->sift, b->sift, SIFT_LENGTH, top);
}

static int nearestone(struct sift_keypoint *q, struct sift_keypoint *t, int nt, double *od)
{
	int bi = -1;
	double bd = INFINITY;
	FORI(nt) {
		double nb = dist_descs(q, t+i);
		if (nb < bd) {
			bd = nb;
			bi = i;
		}
	}
	assert(bi >= 0);
	*od = bd;
	return bi;
}

// returns i such that t[i] is as close as possible to q
static int fancynearest(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		double *od, double *opd)
{
	int besti = -1;
	double bestd = INFINITY, bestprev = bestd;
	FORI(nt) {
		double nb = dist_descs(q, t+i);
		if (nb < bestd) {
			bestprev = bestd;
			bestd = nb;
			besti = i;
		}
	}
	assert(besti >= 0);
	*od = bestd;
	*opd = bestprev;
	return besti;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int fancynearestt(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax)
{
	int besti = -1;
	float bestd = INFINITY;
	FORI(nt) {
		// TODO: be bold and change "dmax" to "bestd" in the call below
		float nb = dist_descst(q, t+i, dmax);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int fancynearestt_rad(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax, float dx, float dy)
{
	int besti = -1;
	float bestd = INFINITY;
	FORI(nt) {
		if (fabs(q->pos[0] - t[i].pos[0]) > dx) continue;
		if (fabs(q->pos[1] - t[i].pos[1]) > dy) continue;
		float nb = dist_descst(q, t+i, dmax);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int find_closest_keypoint(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax, float dx, float dy)
{
	int besti = -1;
	float bestd = dmax;
	FORI(nt) {
		if (fabs(q->pos[0] - t[i].pos[0]) > dx) continue;
		if (fabs(q->pos[1] - t[i].pos[1]) > dy) continue;
		float nb = sift_distance_topped(q, t+i, bestd);
		if (nb < bestd) {
			bestd = nb;
			besti = i;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

int (*siftlike_getpairs(
		struct sift_keypoint *ka, int na, 
		struct sift_keypoint *kb, int nb,
		int *np
		))[3]
{
	int (*p)[3] = xmalloc(na * sizeof * p);
	FORI(na) {
		double d;
		p[i][0] = i;
		p[i][1] = nearestone(ka+i, kb, nb, &d);
		p[i][2] = 10*log(d);
	}
	*np = na;
	return p;
}

static int compare_annpairs(const void *a, const void *b)
{
	float x = ((const struct ann_pair *)a)->v[0];
	float y = ((const struct ann_pair *)b)->v[0];
	int r = (x > y) - (x < y);
	//fprintf(stderr, "CMP %g %g = %d\n", x, y, r);
	return r;
}
static void sort_annpairs(struct ann_pair *t, int n)
{
	qsort(t, n, sizeof*t, compare_annpairs);
}

// get two lists of points, and produce a list of pairs
// (nearest match from a to b)
struct ann_pair *siftlike_get_annpairs(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *np
		)
{
	struct ann_pair *p = xmalloc(na * sizeof * p);
	FORI(na) {
		double d, dp;
		p[i].from = i;
		p[i].to = fancynearest(ka+i, kb, nb, &d, &dp);
		p[i].v[0] = d;
		p[i].v[1] = dp;
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, na);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*np = na;
	return p;
}

// get two lists of points, and produce a list of pairs
// (first nearest matches)
struct ann_pair *siftlike_get_accpairs(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *onp,
		float t
		)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float d;
		int j = fancynearestt(ka + i, kb, nb, &d, t);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

// get two lists of points, and produce a list of pairs
// (first nearest matches)
struct ann_pair *siftlike_get_accpairsrad(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *onp,
		float t, float dx, float dy
		)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float d;
		int j = fancynearestt_rad(ka + i, kb, nb, &d, t, dx, dy);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

#include "ok_list.c"
#include "grid.c"

struct ok_grid {
	struct ok_list l[1];
	struct grid g[1];
	int *buf;
};

static void ok_grid_init(struct ok_grid *o, int nr, int np,
		float x0[2], float dx[2], int n[2])
{
	assert(nr == n[0] * n[1]);
	ok_init(o->l, nr, np);
	grid_init(o->g, 2, x0, dx, n);
	o->buf = xmalloc(np*sizeof*o->buf);
}

static void ok_grid_free(struct ok_grid *o)
{
	ok_free(o->l);
	free(o->buf);
}

static int ok_grid_add_point(struct ok_grid *o, int p, float x[2])
{
	int r = grid_locate(o->g, x);
	ok_add_point(o->l, r, p);
	return r;
}

static int ok_neighboring_points(struct ok_grid *o, float x[2])
{
	int r[4], nr = grid_locate_overlapping(r, o->g, x);
	assert(nr <= 4);
	//for (int i=0;i<nr;i++) fprintf(stderr, "\trglos{%g %g} [%d:%d] = %d\n", x[0], x[1], i, nr, r[i]);
	int cx = 0;
	for (int i = 0; i < nr; i++)
	{
		int nri = ok_which_points(o->l, r[i]);
		for (int j = 0; j < nri; j++)
			o->buf[cx++] = o->l->buf[j];
	}
	return cx;
}

// returns i such that t[i] is as close as possible to q
// if the distance is largest than dmax, return -1
static int find_closest_in_grid(struct sift_keypoint *q,
		struct sift_keypoint *t, int nt,
		float *od, float dmax, struct ok_grid *g)
{
	// build the list of neighbors to traverse
	int nn = ok_neighboring_points(g, q->pos);
	int *nbuf = g->buf;
	//fprintf(stderr, "site ( %g , %g ) has %d neighbors\n", q->pos[0], q->pos[1], nn);
	//for(int i = 0; i < nn; i++) fprintf(stderr, "\t%d\n", nbuf[i]);
	
	// compute the closest point on this list
	int besti = -1;
	float bestd = dmax;
	FORI(nn) {
		int idx = nbuf[i];
		if (fabs(q->pos[0] - t[idx].pos[0]) > g->g->dx[0]) continue;
		if (fabs(q->pos[1] - t[idx].pos[1]) > g->g->dx[1]) continue;
		float nb = sift_distance_topped(q, t+idx, bestd);
		if (nb < bestd) {
			bestd = nb;
			besti = idx;
		}
	}
	if (besti < 0) return -1;
	assert(besti >= 0);
	*od = bestd;
	return bestd < dmax ? besti : -1;
}

struct ann_pair *compute_sift_matches_locally(int *onp,
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		float t, float dx, float dy, int w, int h)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }

	// build a grid structure for each list of keypoints
	// NOTE: only the grid of "ka" is actually used
	float x0[2] = {0, 0};
	float dxy[2] = {dx, dy};
	int n[2] = {1+(w-1)/dx, 1+(h-1)/dy};
	int nr = n[0] * n[1];
	//struct ok_grid ga[1]; ok_grid_init(ga, nr, na, x0, dxy, n);
	struct ok_grid gb[1]; ok_grid_init(gb, nr, nb, x0, dxy, n);
	//for (int i = 0; i < na; i++) ok_grid_add_point(ga, i, ka[i].pos);
	for (int i = 0; i < nb; i++) ok_grid_add_point(gb, i, kb[i].pos);
	//ok_display_tables(gb->l);
	//ok_assert_consistency(gb->l);

	// compute the pairs
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	for (int i = 0; i < na; i++) {
		float d;
		int j = find_closest_in_grid(ka + i, kb, nb, &d, t, gb);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	sort_annpairs(p, np);

	// free the grid structures
	//ok_grid_free(ga);
	ok_grid_free(gb);

	// return values
	*onp = np;
	return p;
}


struct ann_pair *compute_sift_matches(int *onp,
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		float t, float dx, float dy)
{
	if (na == 0 || nb == 0) { *onp=0; return NULL; }
	struct ann_pair *p = xmalloc(na * sizeof * p);
	int np = 0;
	for (int i = 0; i < na; i++) {
		float d;
		int j = find_closest_keypoint(ka + i, kb, nb, &d, t, dx, dy);
		if (j >= 0) {
			p[np].from = i;
			p[np].to = j;
			p[np].v[0] = d;
			p[np].v[1] = NAN;
			np += 1;
		}
	}
	sort_annpairs(p, np);
	*onp = np;
	return p;
}


// get three lists of points, and produce a list of matching triplets
// (first nearest matches)
struct ann_trip *siftlike_get_triplets(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		struct sift_keypoint *kc, int nc,
		int *onp,
		float t
		)
{
	if (na == 0 || nb == 0 || nc == 0) { *onp=0; return NULL; }
	struct ann_trip *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float db, dc;
		int jb = fancynearestt(ka + i, kb, nb, &db, t);
		int jc = fancynearestt(ka + i, kc, nc, &dc, t);
		if (jb >= 0 && jc >= 0) {
			p[np].froma = i;
			p[np].tob = jb;
			p[np].toc = jc;
			p[np].v[0] = hypot(db,dc);
			p[np].v[1] = db;
			p[np].v[2] = dc;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	//sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

// get three lists of points, and produce a list of matching triplets
// (first nearest matches)
struct ann_trip *siftlike_get_tripletsrad(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		struct sift_keypoint *kc, int nc,
		int *onp,
		float t, float rx, float ry
		)
{
	if (na == 0 || nb == 0 || nc == 0) { *onp=0; return NULL; }
	struct ann_trip *p = xmalloc(na * sizeof * p);
	int np = 0;
	FORI(na) {
		float db, dc;
		int jb = fancynearestt_rad(ka + i, kb, nb, &db, t, rx, ry);
		int jc = fancynearestt_rad(ka + i, kc, nc, &dc, t, rx, ry);
		if (jb >= 0 && jc >= 0) {
			p[np].froma = i;
			p[np].tob = jb;
			p[np].toc = jc;
			p[np].v[0] = hypot(db,dc);
			p[np].v[1] = db;
			p[np].v[2] = dc;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	//sort_annpairs(p, np);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

// get two lists of points, and produce a list of pairs
// (all nearest matches)
struct ann_pair *siftlike_get_allpairs(
		struct sift_keypoint *ka, int na, 
		struct sift_keypoint *kb, int nb,
		int *onp,
		float t
		)
{
	struct ann_pair *p = NULL;//xmalloc(na * sizeof * p);
	int np = 0, np_top = 0;
	FORI(na) FORJ(nb) {
		//double d = dist_descs(ka+i, kb+j);
		double d = dist_descst(ka+i, kb+j, t);
		if (d < t) {
			if (np >= np_top) {
				np_top = 0x100 + 2 * (np_top + 1);
				p = xrealloc(p, np_top * sizeof * p);
			}
			struct ann_pair *pi = p + np;
			pi->from = i;
			pi->to = j;
			pi->v[0] = d;
			pi->v[1] = NAN;
			np += 1;
		}
	}
	//FORI(na) fprintf(stderr, "BEFORE p[%d].from=%d\n", i, p[i].from);
	fprintf(stderr, "NOW, to sort the %d pairs\n", np);
	sort_annpairs(p, np);
	FILE *g = xfopen("/tmp/siftpairsd.txt", "w");
	FORI(np) fprintf(g, "%g\n", p[i].v[0]);
	xfclose(g);
	//FORI(na) fprintf(stderr, "AFTER p[%d].from=%d\n", i, p[i].from);
	*onp = np;
	return p;
}

#if USE_KDTREE

#define POINTY_HACK 4387
static void *int_to_pointer(int x)
{
	void *p = (void *)(x + POINTY_HACK);
	return p;
}
static int pointer_to_int(void *p)
{
	int x = -POINTY_HACK + (int)p;
	return x;
}

// using kdtrees,
// get two lists of points, and produce a list of pairs
// (nearest match from a to b)
#endif
struct ann_pair *siftlike_get_annpairs_kd(
		struct sift_keypoint *ka, int na,
		struct sift_keypoint *kb, int nb,
		int *np
		)
{
#if USE_KDTREE
	struct ann_pair *p = xmalloc(na * sizeof * p);
	struct kdtree *kd = kd_create(SIFT_LENGTH);
	fprintf(stderr, "inserting points on kdtree\n");
	FORI(nb)
		kd_insertf(kd, kb[i].sift, int_to_pointer(i));
	fprintf(stderr, "\t...done\n");
	fprintf(stderr, "finding nearest neighbors\n");
	FORI(na) {
		float *f = ka[i].sift;
		struct kdres *res = kd_nearestf(kd, f);
		if (1 == kd_res_size(res)) {
			float fn[SIFT_LENGTH];
			void *pi = kd_res_itemf(res, fn);
			if (!pi) fail("epa aquí!");
			p[i].from = i;
			p[i].to = pointer_to_int(pi);
			p[i].v[0] = distppf(fn, f, SIFT_LENGTH);
			p[i].v[1] = INFINITY;
		} else
			fail("no en tenim cap");
	}
	fprintf(stderr, "\t...done\n");
	kd_free(kd);

	sort_annpairs(p, na);
	*np = na;
	return p;
#else
	return siftlike_get_annpairs(ka, na, kb, nb, np);
#endif
}

//int compare_annpairs(const struct ann_pair *a, const struct ann_pair *b)

SMART_PARAMETER(PSIFT_DISTTHRE,-1)
SMART_PARAMETER(PSIFT_LOWERATIO,-1)

// mask pairs by using some thresholds
void siftlike_maskpairs(bool *mask, struct ann_pair *t, int n)
{
	FORI(n) mask[i] = true;
	float dthre = PSIFT_DISTTHRE();
	float lrat = PSIFT_LOWERATIO();
	if (dthre > 0)
		FORI(n)
			if (t[i].v[0] > dthre)
				mask[i] = false;
	if (lrat > 0)
		FORI(n)
			if (t[i].v[0] > lrat * t[i].v[1])
				mask[i] = false;
}


static float sift_d2(struct sift_keypoint *a, struct sift_keypoint *b)
{
	float dx = b->pos[0] - a->pos[0];
	float dy = b->pos[1] - a->pos[1];
	return hypot(dx,dy);
}

static float sift_d128(struct sift_keypoint *a, struct sift_keypoint *b)
{
	return dist_descs(a, b);
	//return distppf(a->sift, b->sift, SIFT_LENGTH);
}

// TODO: make this function run in linear or near-linear time
// (by assuming that the points are uniformly distributed on the plane, and
// using an efficient data structure for localizing them)
// Right now, it is prohibitively slow for the common case of more than 10.000
// points
void sift_remove_redundancy(bool *mask,
		struct sift_keypoint *t, int n,
		float dist_plane, float dist_sift)
{
	FORI(n) mask[i] = true;

	FORI(n) FORJ(i)
		if (mask[i])
			if (sift_d2(t+i, t+j) < dist_plane)
				if (sift_d128(t+i, t+j) < dist_sift)
					mask[j] = false;
}

static
int splitloc(int *t, int w, float rx, float ry, float ox, float oy, float *x)
{
	int r = 0;
	int ix = (int)(x[0]/(rx-ox));
	int iy = (int)(x[1]/(ry-oy));
	fprintf(stderr, "\tix=%d, rx=%g, ox=%g, x[0]=%g\n", ix,rx,ox,x[0]);
	fprintf(stderr, "\tiy=%d, ry=%g, oy=%g, x[1]=%g\n", iy,ry,oy,x[1]);
	assert(ix >= 0);
	assert(ix < w);
	assert(iy >= 0);
	assert(ix * (rx-ox) <= x[0]);
	assert(iy * (ry-oy) <= x[1]);
	assert((ix+1) * (rx-ox) >= x[0]);
	assert((iy+1) * (ry-oy) >= x[1]);
	t[r++] = iy*w + ix;
	bool ovx = ix > 0 && ix * (rx - ox) + ox >= x[0];
	bool ovy = iy > 0 && iy * (ry - oy) + oy >= x[1];
	if (ovx) t[r++] = iy*w + ix - 1;
	if (ovy) t[r++] = (iy-1)*w + ix;
	if (ovx && ovy) t[r++] = (iy-1)*w + ix - 1;
	fprintf(stderr, "\t\tt={%d %d %d %d}\n", t[0], t[1], t[2], t[3]);
	return r;
}

int siftsplit(struct sift_keypoint *p, int n,
		float rx, float ry, float ox, float oy,
		int (*mask)[5])
{
	FORI(n) mask[i][0] = 0;
	FORI(n) FORJ(4) mask[i][1+j] = -1;
	float maxx = -INFINITY; FORI(n) maxx = fmax(maxx,p[i].pos[0]);
	float maxy = -INFINITY; FORI(n) maxy = fmax(maxy,p[i].pos[1]);
	//int w = calcula_amplada(...);
	//int r = calcula_nombre_de_rectangles(rx,ry,ox,oy,maxx,maxy);
	int w = 1+(int)(maxx/(rx-ox));
	int r = (1+w) * (1+(int)(maxy/(ry-oy)));
	fprintf(stderr, "w = %d, r = %d\n", w, r);
	FORI(n) {
		fprintf(stderr, "p[%d] = %g , %g\n",i,p[i].pos[0],p[i].pos[1]);
		mask[i][0] = splitloc(mask[i]+1, w, rx, ry, ox, oy, p[i].pos);
		fprintf(stderr, "belongs to %d region:\n", mask[i][0]);
		FORJ(mask[i][0])
			fprintf(stderr, "\t%d\n", mask[i][1+j]);
		fprintf(stderr, "\n");
	}
	return r;
}

static void affine_mapf(float y[2], float A[6], float x[2])
{
	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
}

static void homographic_mapf(float y[2], float H[9], float x[2])
{
	float z = H[6]*x[0] + H[7]*x[1] + H[8];
	y[0] = (H[0]*x[0] + H[1]*x[1] + H[2])/z;
	y[1] = (H[3]*x[0] + H[4]*x[1] + H[5])/z;
	//y[0] = x[0];
	//y[1] = x[1];
}

void siftaff(struct sift_keypoint *t, int n, float A[9])
{
	float det = A[0]*A[4] - A[1]*A[3];
	fprintf(stderr, "det = %g\n", det);
	FORI(n) {
		struct sift_keypoint *k = t+i;
		float vec[2] = {cos(k->orientation), sin(k->orientation)};
		float x[2], rvec[2];
		affine_mapf(x, A, k->pos);
		rvec[0] = A[0]*vec[0] + A[1]*vec[1];
		rvec[1] = A[3]*vec[0] + A[4]*vec[1];
		FORJ(2) k->pos[j] = x[j];
		k->scale *= det;
		k->orientation = atan2(rvec[1], rvec[0]);
	}
}

void sifthom(struct sift_keypoint *t, int n, float H[9])
{
	// TODO XXX ERROR FIXME : update the scale and orientation accordingly!
	FORI(n) {
		struct sift_keypoint *k = t+i;
		float x[2];
		homographic_mapf(x, H, k->pos);
		FORJ(2) k->pos[j] = x[j];
	}
}
