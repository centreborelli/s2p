#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>


#include "ok_list.h"

#include "fail.c"
#include "xmalloc.c"


#ifndef FORI
#define FORI(n) for(int i = 0; i < (n); i++)
#endif
#ifndef FORJ
#define FORJ(n) for(int j = 0; j < (n); j++)
#endif

#define GETBIT(x,i) ((x)&(1<<(i)))

#ifndef VERBOSE
#define VERBOSE false
#endif

#ifndef DEBUG
#define DEBUG(...) if(VERBOSE)fprintf(stderr,__VA_ARGS__)
#endif

#define INULL (-42)
#define REMOVED (-4300)

void ok_display_tables(struct ok_list *l)
{
	printf("r:");
	FORI(l->number_of_regions)
		if (l->r[i] != INULL)
			printf(" [%d]%d", i, l->r[i]);
	printf("\np:");
	FORI(l->number_of_points)
		if (l->p[i] != INULL)
			printf(" [%d]%d", i, l->p[i]);
	printf("\nplist:");
	FORI(l->number_of_points)
		printf(" [%d]%d", i, l->plist[i]);
	printf("\n");
}

static void ok_assert_consistency(struct ok_list *l)
{
	//fail("caqueta");
	DEBUG("\nASSERTING CONSISTENCY...\n");
	//ok_display_tables(l);
	FORI(l->number_of_regions)
		if (l->r[i] != INULL)
		{
			assert(l->r[i] >= 0);
			assert(l->r[i] < l->number_of_points);
			if (l->p[l->r[i]] != REMOVED)
				assert(l->p[l->r[i]] == i);
		}
	FORI(l->number_of_points)
		if (l->p[i] != INULL && l->p[i] != REMOVED)
		{
			assert(l->p[i] >= 0);
			assert(l->p[i] < l->number_of_regions);
			if (l->p[l->r[l->p[i]]] != REMOVED)
				assert(l->p[l->r[l->p[i]]] == l->p[i]);
			if (l->p[l->plist[i]] != REMOVED)
				assert(l->p[l->plist[i]] == l->p[i]);
		}
	DEBUG("...ASSERTED!\n");
}
void ok_hack_assert_consistency(struct ok_list *l){ok_assert_consistency(l);
	/*fprintf(stderr,"...");*/}

void ok_init(struct ok_list *l, int nr, int np)
{
	l->number_of_regions = nr;
	l->number_of_points = np;
	if (np < 26) np += 26;
	l->r =     xmalloc(nr * sizeof(*l->r));
	l->p =     xmalloc(np * sizeof(*l->p));
	l->plist = xmalloc(np * sizeof(*l->p));
	l->buf =   xmalloc(np * sizeof(*l->p));
	FORI(l->number_of_regions)
		l->r[i] = INULL;
	FORI(l->number_of_points)
	{
		l->p[i] = l->buf[i] = INULL;
		l->plist[i] = i;
	}
	//ok_assert_consistency(l);
}

void ok_free(struct ok_list *l)
{
	xfree(l->r);
	xfree(l->p);
	xfree(l->plist);
	xfree(l->buf);
	l->number_of_regions = l->number_of_points = 0;
}

// TODO: optimize this so that it skips the removed elements
int ok_which_points(struct ok_list *l, int r)
{
	//DEBUG("going to list the points of region %d\n", r);
	assert(r >= 0);
	if (r >= l->number_of_regions)
		fail("caca");
	assert(r < l->number_of_regions);
	int cx = 0, p = l->r[r];
	if (p != INULL)
	{
		if (l->p[p] != REMOVED)
			l->buf[cx++] = p;
		while (p != l->plist[p])
		{
			p = l->plist[p];
			if (l->p[p] != REMOVED)
				l->buf[cx++] = p;
		}
	}
	//DEBUG("\tthere were %d of them\n", cx);
	return cx;
}

int ok_which_region(struct ok_list *l, int p)
{
	assert(p >= 0);
	assert(p < l->number_of_points);
	int r = l->p[p];
	if (r == INULL)
		fail("el punt %d no pertany a cap regió", p);
	if (r == REMOVED)
	{
		assert(REMOVED < 0);
		return REMOVED;
	}
	return r;
}

//static
void ok_display(struct ok_list *l)
{
	printf("ok_list %p\n", (void *)l);
	FORI(l->number_of_regions)
	{
		int c = ok_which_points(l, i);
		if (c) {
			printf("region %d has %d points:", i, c);
			FORJ(c)
				printf(" %d", l->buf[j]);
			printf("\n");
		}
	}
}

//// returns the next point not removed, or INULL
//static int next_not_removed(ok_list *l, int p)
//{
//	assert(REMOVED == l->p[l->plist[p]]);
//	do {
//		p = l->plist[p];
//	} while (REMOVED == l->p[l->plist[p]] && p != l->plist[p]);
//	return INULL;
//}

//TODO: optimize this function so that it always leaves a non-removed element
//as a representatitve
void ok_remove_point(struct ok_list *l, int p)
{
	assert(p >= 0);
	assert(p < l->number_of_points);
	if (l->p[p] < 0)
		fail("error, point of index %d had a rreg of %d\n",p,l->p[p]);
	assert(l->p[p] >= 0);
	DEBUG("{{{}}} REMOVEing point of index %d (which belang to region %d)\n",
			p, l->p[p]);
	//if (l->r[l->p[p]] == p)
	//{
	//	l->r[l->p[p]] = (l->plist[p] == p) ? INULL : l->plist[p];
	//}
	l->p[p] = REMOVED;
#ifdef USE_MAIN
	ok_display_tables(l);
	ok_assert_consistency(l);
#endif
}

void ok_add_point(struct ok_list *l, int r, int p)
{
	assert(r >= 0);
	assert(r < l->number_of_regions);
	assert(p >= 0);
	assert(p < l->number_of_points);

	if (l->p[p] != INULL)
		fail("point %d is already added to region %d, can not be newly added to region %d\n", p, l->p[p], r);

	DEBUG("{{{}}} ADDing point %d to region %d\n", p, r);

	if (l->r[r] == INULL)
	{
		l->r[r] = p;
		l->p[p] = r;
		l->plist[p] = p;
	} else {
		l->plist[p] = l->r[r];
		l->r[r] = p;
		l->p[p] = r;
	}

#ifdef USE_MAIN
	ok_assert_consistency(l);
#endif
}

//#ifdef USE_IMAGE_STRUCTURES

// input: fill x0, dx and nx
void ok_init_grid(struct ok_list *l, int np)
//float x0[3], float dx[3], int nx[3], int np)
{
	assert(l->dx[0] > 0 && l->dx[1] > 0 && l->dx[2] > 0);
	assert(l->nx[0] > 0 && l->nx[1] > 0 && l->nx[2] > 0);
	ok_init(l, l->nx[0] * l->nx[1] * l->nx[2], np);
	//l->px = NULL;
	//FORI(3) l->x0[i] = x0[i];
	//FORI(3) l->dx[i] = dx[i];
	//FORI(3) l->nx[i] = nx[i];
	fprintf(stderr, "init grid\n\tnx=(%d %d %d)\n"
			"\tdx=(%f %f %f)\n\tx0=(%f %f %f)\n",
			l->nx[0], l->nx[1], l->nx[2],
			l->dx[0], l->dx[1], l->dx[2],
			l->x0[0], l->x0[1], l->x0[2]);
}

// fills buf[0], buf[1], buf[2] with integer coordinates of the region
int ok_regionindex(struct ok_list *l, float x[3])
{

	bool ass = true;
	for (int i = 0; i < 3; i++)
		if (x[i] < l->x0[i] || x[i] > l->x0[i] + l->nx[i] * l->dx[i])
			ass = false;
	if (!ass)
	{
		for (int i = 0; i < 3; i++)
			x[i] = fmax(x[i], l->x0[i]);
		for (int i = 0; i < 3; i++)
			x[i] = fmin(x[i], l->x0[i] + l->nx[i] * l->dx[i]);
		//error("cacona de punt %g %g %g", x[0], x[1], x[2]);
		//goto ridxerr;
	}
	// fer les divisions i tal que toquin:
	int r, ix[3];
	for (int i = 0; i < 3; i++)
	{
		ix[i] = floor((x[i] - l->x0[i])/l->dx[i]);
		if (ix[i] == l->nx[i]) ix[i]--;
		if(1){if(ix[i]<0||ix[i]>=l->nx[i])goto ridxerr;}
		assert(ix[i] >= 0);
		assert(ix[i] < l->nx[i]);
		l->buf[i] = ix[i];
	}
	r = ix[2] * (l->nx[0] * l->nx[1]) + ix[1] * l->nx[0] + ix[0];
	if(1){if(r<0||r>=l->nx[0]*l->nx[1]*l->nx[2])goto ridxerr;}
	assert(r >= 0);
	assert(r < l->nx[0] * l->nx[1] * l->nx[2]);
	return r;
ridxerr:
	fail("x=(%g,%g,%g)\n"
			"x0=(%g,%g,%g)\n"
			"dx=(%g,%g,%g)\n"
			"nx=(%d,%d,%d)\n"
			"\tix=(%d,%d,%d)", 
			x[0], x[1], x[2],
			l->x0[0], l->x0[1], l->x0[2],
			l->dx[0], l->dx[1], l->dx[2],
			l->nx[0], l->nx[1], l->nx[2],
			ix[0], ix[1], ix[2]);
	return -1;
}

static bool ok_innerP(struct ok_list *l, int x[3])
{
	FORI(3) if (x[i] < 0 || x[i] >= l->nx[i]) return false;
	return true;
}

static int ok_idx(struct ok_list *l, int x[3])
{
	FORI(3) assert(x[i] < l->nx[i]);
	int r = x[0] + x[1] * (l->nx[0]) + x[2] * (l->nx[0] * l->nx[1]);
	assert(r >= 0);
	assert(r < l->nx[0] * l->nx[1] * l->nx[2]);
	return r;
}


// return the 8 closest regions to a point (4, in 2D), caring for borders
// (fills buf)
#if 0
int ok_regionindex_neigs(ok_list *l, float x[3])
{
	int ridx = ok_regionindex(l, x);
	l->buf[0] = ridx;
	return 1;
}
#else
//THIS IS INCOMPREHENSIBLE
int ok_regionindex_neigs(struct ok_list *l, float x[3])
{
	DEBUG("OK_RIDX_N: treating point %f,%f,%f...\n", x[0], x[1], x[2]);
	int ridx = ok_regionindex(l, x);
	DEBUG("OK_RIDX_N: ridx = %d\n", ridx);
	int *pcor = l->buf;
	int *desp = l->buf + 3;
	FORI(3)
		desp[i] = x[i] - (l->x0[i] + l->dx[i]*pcor[i]) > 0.5*l->dx[i]
			? 1 : -1;
	//fprintf(stderr, "\tpcor = %d,%d,%d\n", pcor[0], pcor[1], pcor[2]);
	//fprintf(stderr, "\tdesp = %d,%d,%d\n", desp[0], desp[1], desp[2]);

	int p[8][3];
	FORI(8) FORJ(3) {
		int d = GETBIT(i,j) ? desp[j] : 0;
	       	p[i][j] = pcor[j] + d;
		//fprintf(stderr, "\t\td=%d   p[%d][%d]=%d\n", d, i, j, p[i][j]);
	}
	//FORI(8) fprintf(stderr, "\tp[%d] = (%d %d %d)\n", i,
	//		p[i][0], p[i][1], p[i][2]);

	int cx = 0;
	FORI(8) if (ok_innerP(l, p[i])) l->buf[cx++] = ok_idx(l, p[i]);
	//fprintf(stderr, "there are %d voxels surrounding %f,%f,%f\n", cx, x[0], x[1], x[2]);
	//FORI(cx) fprintf(stderr, "buf[%d] = %d\n", i, l->buf[i]);
	FORI(cx) assert(l->buf[i] >= 0 && l->buf[i] < l->number_of_regions);
	assert(cx > 0);
	return cx;
}
#endif

int ok_add_geo_point(struct ok_list *l, float x[3], int p)
{
	{
		DEBUG("OKUPA: agp[%d %d][%g %g %g][%g %g %g] [%d %d %d]"
				"(%g %g %g), %d\n",
				l->number_of_regions, l->number_of_points,
				l->x0[0], l->x0[1], l->x0[2],
				l->dx[0], l->dx[1], l->dx[2],
				l->nx[0], l->nx[1], l->nx[2],
				x[0], x[1], x[2], p);
	}
	int ri = ok_regionindex(l, x);
	ok_add_point(l, ri, p);
	return ri;
}
//#endif /* USE_IMAGE_STRUCTURES */

void ok_svg_layer(void *ff, struct ok_list *l)
{
	FILE *f = ff;
	FORI(l->nx[0]) // segments verticals, en variant la x
	{
		float x1 = l->x0[0] + l->dx[0] * i;
		float x2 = x1;
		float y1 = l->x0[1];
		float y2 = l->x0[1] + l->dx[1] * l->nx[1];
		fprintf(f, "<line "
				"x1=\"%f\" y1=\"%f\" "
				"x2=\"%f\" y2=\"%f\" "
				"stroke=\"green\" "
				"stroke-width=\"0.1\" "
				"/>\n", x1, y1, x2, y2);
	}
	FORI(l->nx[1]) // segments horitzontals, en variant la y
	{
		float x1 = l->x0[0];
		float x2 = l->x0[0] + l->dx[0] * l->nx[0];
		float y1 = l->x0[1] + l->dx[1] * i;
		float y2 = y1;
		fprintf(f, "<line "
				"x1=\"%f\" y1=\"%f\" "
				"x2=\"%f\" y2=\"%f\" "
				"stroke=\"green\" "
				"stroke-width=\"0.1\" "
				"/>\n", x1, y1, x2, y2);
	}
	FORJ(l->nx[1]) FORI(l->nx[0]) {
		float xt[3];
		xt[0] = l->x0[0] + l->dx[0] * (i + 0.1);
		xt[1] = l->x0[1] + l->dx[1] * (j + 0.2);
		xt[2] = 0.1;
		//fprintf(stderr, "i,j,x= %d, %d, (%f %f %f)\n", i, j, xt[0], xt[1], xt[2]);
		int ridx = ok_regionindex(l, xt);
		int rnp = ok_which_points(l, ridx);
		if (rnp)
			fprintf(f, "<text "
					"x=\"%f\" y=\"%f\" "
					"fill=\"green\" "
					"font-size=\"1\">"
					"%d, %d</text>\n",
					xt[0], xt[1], ridx, rnp);
	}
	/*
	if (l->px)
		FORI(l->number_of_points)
		{
			int ri = ok_which_region(l, i);
			if (ri >= 0)
			{
				fprintf(f, "<circle cx=\"%f\" cy=\"%f\" r="
						"\"0.1\" fill=\"green\" />\n",
						l->px[i][0], l->px[i][1]);
				fprintf(f, "<text "
						"x=\"%f\" y=\"%f\" "
						"fill=\"green\" "
						"font-size=\"0.1\">"
						"%d, %d</text>\n",
						l->px[i][0], l->px[i][1], i, ri);
			}
		}
		*/
}

#ifdef USE_MAIN
static void programa_interactiu()
{
	int c, r, p;
	struct ok_list l;
	printf("introduce number of regions and number of points\n");
	scanf("%d %d", &r, &p);
	fflush(stdin);
	fflush(stdout);
	ok_init(&l, r, p);
	printf("using %d regions and %d points.  Available commands:\n"
			"\ta r p (add to region r the point p)\n"
			"\td p (delete point p)\n"
			"\tt (display tables)\n"
			"\tq (quit)\n",
			r, p);
	while ((c = getchar()) != EOF)
		switch(c){
		case 'a':
			scanf(" %d %d", &r, &p);
			ok_add_point(&l, r, p);
			ok_display(&l);
			break;
		case 'd':
			scanf(" %d", &p);
			ok_remove_point(&l, p);
			ok_display(&l);
			break;
		case 't':
			ok_display_tables(&l);
			break;
		case '\n': continue;
		case 'q': goto ent;
		default: printf("unrecognized command %c\n", c);
		}
ent:
	ok_display(&l);
	ok_free(&l);
}
int main()
{
	struct ok_list l;
	ok_init(&l, 3, 10);
	ok_display(&l); printf("\n");
	ok_add_point(&l, 1, 3); ok_display(&l); printf("\n");
	ok_add_point(&l, 1, 5); ok_display(&l); printf("\n");
	ok_add_point(&l, 1, 7); ok_display(&l); printf("\n");
	ok_add_point(&l, 0, 2); ok_display(&l); printf("\n");
	ok_add_point(&l, 1, 0); ok_display(&l); printf("\n");
	//ok_add_point(&l, 2, 5); ok_display(&l); printf("\n");
	ok_free(&l);

	programa_interactiu();
	return 0;
}
#endif /* USE_MAIN */
