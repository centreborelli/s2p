#include <assert.h>

#include "fail.c"


// structure to store a grid of rectangular cells,
// useful to locate points in them
struct grid {
	int dim;     // dimension of the base space
	float x0[4]; // first corner of the grid
	float dx[4]; // cell size
	int n[4];    // number of cells in each dimension
	int nc;      // total number of cells
};


// fill a grid structure with the specified data
void grid_init(struct grid *g, int dim, float *x0, float *dx, int *n)
{
	if (dim < 2 || dim > 4) fail("bad grid dim %d", dim);
	g->dim = dim;
	g->nc = 1;
	for (int i = 0; i < dim; i++) {
		g->x0[i] = x0[i];
		g->dx[i] = dx[i];
		g->n[i] = n[i];
		g->nc *= n[i];
	}
}


// get the index of a cell, from its integer coordinates
int grid_index_of_cell(struct grid *g, int *i)
{
	// TODO write a proper "for" loop here
	int *n = g->n;
	switch(g->dim) {
	case 2: return i[0] + n[0] * i[1];
	case 3: return i[0] + n[0] * (i[1] + n[1] * i[2]);
	case 4: return i[0] + n[0] * (i[1] + n[1] * (i[2] + n[2] * i[3]));
	default: fail("caca");
	}
}


static
int locate_interval_1d(float x0, float dx, int n, float x)
{
	float xf = x0 + dx*n;
	if (x <= x0) return 0;
	if (x >= xf) return n-1;
	int r = n * (x - x0) / (xf - x0);
	assert(r >= 0);
	assert(r < n);
	return r;
}

static
int locate_widened_interval_1d(int ox[2], float x0, float dx, int n, float x)
{
	int r[3];
	r[0] = locate_interval_1d(x0, dx, n, x - dx/2);
	r[1] = locate_interval_1d(x0, dx, n, x       );
	r[2] = locate_interval_1d(x0, dx, n, x + dx/2);
	assert(r[0] <= r[1]);
	assert(r[1] <= r[2]);
	assert(r[1] - r[0] <= 1);
	assert(r[2] - r[1] <= 1);
	assert(r[2] - r[0] <= 1);
	ox[0] = r[0];
	ox[1] = r[2];
	//fprintf(stderr, "\t\tlwi %g = %d %d %d\n", x, r[0], r[1], r[2]);
	return r[0] == r[2] ? 1 : 2;
}

// compute the integer cell coordinates from the float coordinates of a point
void grid_locate_float_point(int *ox, struct grid *g, float *x)
{
	for (int i = 0; i < g->dim; i++)
		ox[i] = locate_interval_1d(g->x0[i], g->dx[i], g->n[i], x[i]);
}

// compute the index of the cell that contains a given point
int grid_locate(struct grid *g, float *x)
{
	int ix[g->dim];
	grid_locate_float_point(ix, g, x);
	return grid_index_of_cell(g, ix);
}

// computes the indexes of the cells that may contain neighbors
// of the given point
int grid_locate_overlapping(int *buf, struct grid *g, float *x)
{
	int ox[g->dim][2];
	int no[g->dim];
	int ret = 1;
	for (int i = 0; i < g->dim; i++) {
		no[i] = locate_widened_interval_1d(ox[i],
				g->x0[i], g->dx[i], g->n[i], x[i]);
		ret *= no[i];
	}
	//fprintf(stderr, "ret = %d\n", ret);
	// TODO write a proper "for" loop here (if possible, which I'm not sure)
	int cx = 0;
	switch (g->dim) {
	case 2:
		for (int j = ox[1][0]; j <= ox[1][1]; j++)
		for (int i = ox[0][0]; i <= ox[0][1]; i++)
		{
			int p[2] = {i, j};
			buf[cx] = grid_index_of_cell(g, p);
			//fprintf(stderr, "\t\tbuf[%d] = grioc(%d %d)=%d\n",cx,i,j,buf[cx]);
			cx += 1;
		}
		break;
	default: fail("undexined grid overlaps in dimension %d", g->dim);
	}
	//fprintf(stderr, "cx ret = %d %d\n", cx, ret);
	assert(cx == ret);
	return ret;
}
