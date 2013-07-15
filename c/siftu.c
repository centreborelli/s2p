#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "siftie.c"

// compute pairs using sift-nn (non-sym, exhaustive, explicit)
static int main_siftcpairs(int c, char *v[])
{
	if (c != 5) {
		fprintf(stderr,"usage:\n\t%s t k1 k2 pairs.txt\n",*v);
		return EXIT_FAILURE;
	}
	struct sift_keypoint *p[2];
	int n[2];
	FORI(2) {
		FILE *f = xfopen(v[2+i], "r");
		p[i] = read_raw_sifts(f, n+i);
		xfclose(f);
	}
	int npairs;
	struct ann_pair *pairs;
	float t = atof(v[1]);
	pairs = siftlike_get_accpairs(p[0], n[0], p[1], n[1], &npairs, t);
	fprintf(stderr, "SIFTCPAIRS: produced %d pairs "
			"(from %d and %d){%d}[%g%%]\n",
			npairs, n[0], n[1], n[0]*n[1],npairs*100.0/(n[0]*n[1]));
	FILE *f = xfopen(v[4], "w");
	FORI(npairs) {
		struct sift_keypoint *ka = p[0] + pairs[i].from;
		struct sift_keypoint *kb = p[1] + pairs[i].to;
		fprintf(f, "%g %g %g %g\n",
				ka->pos[0], ka->pos[1],
				kb->pos[0], kb->pos[1]);
	}
	xfclose(f);
	FORI(2) if (p[i]) xfree(p[i]);
	if (pairs) xfree(pairs);
	return EXIT_SUCCESS;
}

// compute pairs using sift-nn (non-sym, exhaustive, explicit)
static int main_siftcpairsr(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr,"usage:\n\t%s t k1 k2 rx ry pairs.txt\n",*v);
		//                         0 1 2  3  4  5  6
		return EXIT_FAILURE;
	}
	struct sift_keypoint *p[2];
	int n[2];
	FORI(2) {
		FILE *f = xfopen(v[2+i], "r");
		p[i] = read_raw_sifts(f, n+i);
		xfclose(f);
	}
	int npairs;
	struct ann_pair *pairs;
	float t = atof(v[1]);
	float rx = atof(v[4]);
	float ry = atof(v[5]);
	pairs = siftlike_get_accpairsrad(p[0], n[0], p[1], n[1], &npairs, t, rx, ry);
	fprintf(stderr, "SIFTCPAIRSR: produced %d pairs "
			"(from %d and %d){%d}[%g%%]\n",
			npairs, n[0], n[1], n[0]*n[1],npairs*100.0/(n[0]*n[1]));
	FILE *f = xfopen(v[6], "w");
	FORI(npairs) {
		struct sift_keypoint *ka = p[0] + pairs[i].from;
		struct sift_keypoint *kb = p[1] + pairs[i].to;
		fprintf(f, "%g %g %g %g\n",
				ka->pos[0], ka->pos[1],
				kb->pos[0], kb->pos[1]);
	}
	xfclose(f);
	FORI(2) if (p[i]) xfree(p[i]);
	if (pairs) xfree(pairs);
	return EXIT_SUCCESS;
}

// compute pairs using sift-nn (non-sym, initial segment, explicit)
static int main_siftcpairst(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr,"usage:\n\t%s t k1 k2 top pairs.txt\n",*v);
		//                         0 1 2  3  4   5
		return EXIT_FAILURE;
	}
	struct sift_keypoint *p[2];
	int n[2];
	FORI(2) {
		FILE *f = xfopen(v[2+i], "r");
		p[i] = read_raw_sifts(f, n+i);
		xfclose(f);
	}
	int npairs;
	struct ann_pair *pairs;
	float t = atof(v[1]);
	int top = atoi(v[4]);
	if (top < n[0]) n[0] = top;
	if (top < n[1]) n[1] = top;
	pairs = siftlike_get_accpairs(p[0], n[0], p[1], n[1], &npairs, t);
	fprintf(stderr, "SIFTCPAIRST: produced %d pairs "
			"(from %d and %d){%d}[%g%%]\n",
			npairs, n[0], n[1], n[0]*n[1],npairs*100.0/(n[0]*n[1]));
	FILE *f = xfopen(v[5], "w");
	FORI(npairs) {
		struct sift_keypoint *ka = p[0] + pairs[i].from;
		struct sift_keypoint *kb = p[1] + pairs[i].to;
		fprintf(f, "%g %g %g %g\n",
				ka->pos[0], ka->pos[1],
				kb->pos[0], kb->pos[1]);
	}
	xfclose(f);
	FORI(2) if (p[i]) xfree(p[i]);
	if (pairs) xfree(pairs);
	return EXIT_SUCCESS;
}


// compute triplets using sift-nn
int main_sifttriplets(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr,"usage:\n\t%s t k1 k2 k3 trips.txt\n",*v);
		//                         0 1 2  3  4  5
		return EXIT_FAILURE;
	}
	struct sift_keypoint *p[3];
	int n[3];
	FORI(3) {
		FILE *f = xfopen(v[2+i], "r");
		p[i] = read_raw_sifts(f, n+i);
		xfclose(f);
	}
	int ntrips;
	float t = atof(v[1]);
	struct ann_trip *trips = siftlike_get_triplets(p[0], n[0],
			                               p[1], n[1],
						       p[2], n[2], &ntrips, t);
	fprintf(stderr, "SIFTTRIPS: produced %d triplets "
			"(from %d %d %d points)\n",
			ntrips, n[0], n[1], n[2]);
	FILE *f = xfopen(v[5], "w");
	FORI(ntrips) {
		struct sift_keypoint *ka = p[0] + trips[i].froma;
		struct sift_keypoint *kb = p[1] + trips[i].tob;
		struct sift_keypoint *kc = p[2] + trips[i].toc;
		fprintf(f, "%g %g %g %g %g %g\n",
				ka->pos[0], ka->pos[1],
				kb->pos[0], kb->pos[1],
				kc->pos[0], kc->pos[1]);
	}
	xfclose(f);
	FORI(3) if (p[i]) xfree(p[i]);
	if (trips) xfree(trips);
	return EXIT_SUCCESS;
}

// compute triplets using sift-nn
int main_sifttripletsr(int c, char *v[])
{
	if (c != 8) {
		fprintf(stderr,"usage:\n\t%s t k1 k2 k3 rx ry trips.txt\n",*v);
		//                         0 1 2  3  4  5  6  7
		return EXIT_FAILURE;
	}
	struct sift_keypoint *p[3];
	int n[3];
	FORI(3) {
		FILE *f = xfopen(v[2+i], "r");
		p[i] = read_raw_sifts(f, n+i);
		xfclose(f);
	}
	int ntrips;
	float t = atof(v[1]);
	float rx = atof(v[5]);
	float ry = atof(v[6]);
	struct ann_trip *trips = siftlike_get_tripletsrad(p[0], n[0],
			                               p[1], n[1],
						       p[2], n[2], &ntrips, t,
						       rx, ry);
	fprintf(stderr, "SIFTTRIPS: produced %d triplets "
			"(from %d %d %d points)\n",
			ntrips, n[0], n[1], n[2]);
	FILE *f = xfopen(v[7], "w");
	FORI(ntrips) {
		struct sift_keypoint *ka = p[0] + trips[i].froma;
		struct sift_keypoint *kb = p[1] + trips[i].tob;
		struct sift_keypoint *kc = p[2] + trips[i].toc;
		fprintf(f, "%g %g %g %g %g %g\n",
				ka->pos[0], ka->pos[1],
				kb->pos[0], kb->pos[1],
				kc->pos[0], kc->pos[1]);
	}
	xfclose(f);
	FORI(3) if (p[i]) xfree(p[i]);
	if (trips) xfree(trips);
	return EXIT_SUCCESS;
}

// split a sift list into overlapping rectangular slices
int main_siftsplit(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr,"usage:\n\t%s rx ry ox oy out_fpattern\n",*v);
		return EXIT_FAILURE;
	}
	float rx = atof(v[1]);
	float ry = atof(v[2]);
	float ox = atof(v[3]);
	float oy = atof(v[4]);
	if (rx-ox <= 0 || ry-oy <= 0) fail("negative overlap!");
	struct sift_keypoint *p;
	int n;
	{
		FILE *f = xfopen("-", "r");
		p = read_raw_sifts(f, &n);
		xfclose(f);
	}
	int (*smask)[5] = xmalloc(n*sizeof*smask);
	int ns = siftsplit(p, n, rx, ry, ox, oy, smask);
	if (ns > 400) fail("too many slices (%d), aborting", ns);
	FILE *f[1+ns];
	int nfname = 10 + strlen(v[5]);
	char fname[nfname];
	FORI(ns) {
		snprintf(fname, nfname, v[5], i);
		fprintf(stderr, "i = %d, fname = %s\n", i, fname);
		f[i] = xfopen(fname, "w");
	}
	FORI(n) FORJ(smask[i][0]) {
		int idx = smask[i][1+j];
		if (idx < 0 || idx >= ns) fail("bad mask value %d", idx);
		write_raw_sift(f[idx], p+i);
	}
	xfree(smask);
	FORI(ns) xfclose(f[i]);
	return EXIT_SUCCESS;
}

static void affine_inversion(float invA[6], float A[6])
{
	float a, b, c, d, p, q;
	a=A[0]; b=A[1]; p=A[2];
	c=A[3]; d=A[4]; q=A[5];
	float det = a*d - b*c;
	invA[0] = d;
	invA[1] = -b;
	invA[2] = b*q-d*p;
	invA[3] = -c;
	invA[4] = a;
	invA[5] = c*p-a*q;
	FORI(6) invA[i] /= det;
}
//static void affine_mapf(float y[2], float A[6], float x[2])
//{
//	y[0] = A[0]*x[0] + A[1]*x[1] + A[2];
//	y[1] = A[3]*x[0] + A[4]*x[1] + A[5];
//}

//static void siftaff(struct sift_keypoint *t, int n, float A[6])
//{
//	float det = A[0]*A[4] - A[1]*A[3];
//	fprintf(stderr, "det = %g\n", det);
//	FORI(n) {
//		struct sift_keypoint *k = t+i;
//		float vec[2] = {cos(k->orientation), sin(k->orientation)};
//		float x[2], rvec[2];
//		affine_mapf(x, A, k->pos);
//		rvec[0] = A[0]*vec[0] + A[1]*vec[1];
//		rvec[1] = A[3]*vec[0] + A[4]*vec[1];
//		FORJ(2) k->pos[j] = x[j];
//		k->scale *= det;
//		k->orientation = atan2(rvec[1], rvec[0]);
//	}
//}

SMART_PARAMETER_SILENT(SIFTAFF_INV,0)

// affinely transform sift descriptors
int main_siftaff(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t%s aff <sift.in >sift.out\n", *v);
		return EXIT_FAILURE;
	}
	float A[6], B[6]; FORI(6) B[i] = atof(v[i+1]);
	if (SIFTAFF_INV() > 0.5)
		affine_inversion(A, B);
	else
		FORI(6) A[i] = B[i];
	int n;
	struct sift_keypoint *s = read_raw_sifts(stdin, &n);
	siftaff(s, n, A);
	write_raw_sifts(stdout, s, n);
	return EXIT_SUCCESS;
}

//// mask sift descriptors according to a binary mask
//int main_siftmask(int c, char *v[])
//{
//	if (c != 2) {
//		fprintf(stderr, "usage:\n\t%s m.pgm <sift.in >sift.out\n", *v);
//		return EXIT_FAILURE;
//	}
//	imatge mask[1]; carrega_imatge_general(mask, v[1]);
//	int n;
//	struct sift_keypoint *s = read_raw_sifts(stdin, &n);
//	FORI(n) {
//		int a = s->pos[0];
//		int b = s->pos[1];
//		if (punt_interiorP(mask, a, b, 0) && mask->t[0][b][a] > 0)
//			write_raw_sift(stdout, s);
//		s += 1;
//	}
//	return EXIT_SUCCESS;
//}

// sift file format conversion
int main_siftcon(int c, char *v[])
{
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s [a|b|g] [a|b|g] <in >out\n", *v);
		return EXIT_FAILURE;
	}
	int forma = v[1][0];
	int formb = v[2][0];
	int n;
	struct sift_keypoint *k;
	switch(forma) {
	case 'a': k = read_raw_sifts(stdin, &n); break;
	case 'b': k = read_raw_siftsb(stdin, &n); break;
	case 'g': k = read_raw_sifts_gen(stdin, &n); break;
	default: fail("unrecognized input format '%c'", forma);
	}
	switch(formb) {
	case 'a': write_raw_sifts(stdout, k, n); break;
	case 'b': write_raw_siftsb(stdout, k, n); break;
	case 'g': write_raw_sifts_gen(stdout, k, n); break;
	default: fail("unrecognized output format '%c'", formb);
	}
	xfree(k);
	return EXIT_SUCCESS;
}

// remove redundancy in list of sift descriptors
int main_siftrr(int c, char *v[])
{
	// TODO: make this program run in real-time
	// (by using an appropriate data structure)
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s t2 t128 <sift.i >sift.o\n", *v);
		return EXIT_FAILURE;
	}
	int n;
	struct sift_keypoint *s = read_raw_sifts(stdin, &n);
	bool mask[n];
	sift_remove_redundancy(mask, s, n, atof(v[1]), atof(v[2]));
	FORI(n)
		if (mask[i])
			write_raw_sift(stdout, s+i);
	return EXIT_SUCCESS;
}

int main(int c, char *v[])
{
	if (c < 2) goto usage;
	else if (0 == strcmp(v[1], "pair")) return main_siftcpairs(c-1, v+1);
	else if (0 == strcmp(v[1], "pairr")) return main_siftcpairsr(c-1, v+1);
	else if (0 == strcmp(v[1], "pairt")) return main_siftcpairst(c-1, v+1);
	else if (0 == strcmp(v[1], "trip")) return main_sifttriplets(c-1, v+1);
	else if (0 == strcmp(v[1], "tripr")) return main_sifttripletsr(c-1,v+1);
	else if (0 == strcmp(v[1], "aff")) return main_siftaff(c-1, v+1);
	else if (0 == strcmp(v[1], "split")) return main_siftsplit(c-1, v+1);
	//else if (0 == strcmp(v[1], "mask")) return main_siftmask(c-1, v+1);
	else if (0 == strcmp(v[1], "clean")) return main_siftrr(c-1, v+1);
	else if (0 == strcmp(v[1], "convert")) return main_siftcon(c-1, v+1);
	else {
	usage: fprintf(stderr, "usage:\n\t%s "
				"[pair|trip|aff|split|clean|convert] params\n",
				*v);
		return EXIT_FAILURE;
	}
}
