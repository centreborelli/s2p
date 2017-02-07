#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "iio.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)


#include "xmalloc.c"
#include "synflow_core.c"


static int parse_doubles(double *t, int nmax, const char *s)
{
	int i = 0, w;
	while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
		i += 1;
		s += w;
	}
	return i;
}



int main_synflow(int c, char *v[])
{
	if (c != 6) {
		fprintf(stderr, "usage:\n\t%s model \"params\""
				//         0  1       2
				" in out flow\n", *v);
				//3  4   5
		return EXIT_FAILURE;
	}

	int w, h, pd;
	float *x = iio_read_image_float_vec(v[3], &w, &h, &pd);

	float *y = xmalloc(w * h * pd * sizeof*y);
	float *f = xmalloc(w * h * 2 * sizeof*y);

	int maxparam = 40;
	double param[maxparam];
	int nparams = parse_doubles(param, maxparam, v[2]);

	struct flow_model fm[1];
	produce_flow_model(fm, param, nparams, v[1], w, h);
	fill_flow_field(f, fm, w, h);
	transform_forward(y, fm, x, w, h, pd);

	iio_save_image_float_vec(v[4], y, w, h, pd);
	iio_save_image_float_vec(v[5], f, w, h, 2);

	free(x);
	free(y);
	free(f);

	return EXIT_SUCCESS;
}

#ifndef OMIT_MAIN
int main(int c, char *v[])
{
	return main_synflow(c, v);
}
#endif//OMIT_MAIN
