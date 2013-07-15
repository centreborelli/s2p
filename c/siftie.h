#ifndef _SIFTIE_H
#define _SIFTIE_H

#include <stdbool.h>
#include <stdio.h>


#define SIFT_LENGTH 128

struct sift_keypoint {
	float pos[2], scale, orientation;
	float affinity[6];
	float id;
	float sift[SIFT_LENGTH];
};

struct ann_pair {
	int from, to;
	float v[2];
};

struct ann_trip {
	int froma, tob, toc;
	float v[3];
};

struct sift_keypoint *read_raw_sifts(FILE *f, int *no);
void write_raw_sifts(FILE *f, struct sift_keypoint *, int n);
void write_raw_sift(FILE *f, struct sift_keypoint *);

struct sift_keypoint *read_raw_siftsb(FILE *f, int *no);
void write_raw_siftsb(FILE *f, struct sift_keypoint *, int n);
void write_raw_siftb(FILE *f, struct sift_keypoint *);

struct sift_keypoint *read_raw_sifts_gen(FILE *f, int *no);
void write_raw_sifts_gen(FILE *f, struct sift_keypoint *, int n);
void write_raw_sift_gen(FILE *f, struct sift_keypoint *);


struct sift_keypoint *read_annotated_sifts(FILE *f, int *no);
// a "pair" is actually a triad.  Both indexes, plus a "score".
int (*siftlike_getpairs(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *
		))[3];
struct ann_pair *siftlike_get_annpairs(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *);
struct ann_pair *siftlike_get_annpairs_kd(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *);

struct ann_pair *siftlike_get_allpairs(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *, float);
struct ann_pair *siftlike_get_accpairs(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *, float);
struct ann_pair *siftlike_get_accpairsrad(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *, float, float, float);
struct ann_trip *siftlike_get_triplets(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *, float);
struct ann_trip *siftlike_get_tripletsrad(
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		struct sift_keypoint *, int,
		int *, float, float, float);

void siftlike_maskpairs(bool *mask, struct ann_pair *t, int n);

void sift_remove_redundancy(bool *mask,
		struct sift_keypoint *, int n,
		float dist_plane, float dist_sift);

//void sort_annpairs(struct ann_pair *, int );

int siftsplit(struct sift_keypoint *p, int n,
		float rx, float ry, float ox, float oy,
		int (*mask)[5]);

void sifthom(struct sift_keypoint *t, int n, float h[9]);

#endif//_SIFTIE_H
