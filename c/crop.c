#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "iio.h"

int main (int argc, char **argv)
{
    int x0 = 0;
    int y0 = 0;
    int w = 1;
    int h = 1;

    /* parameter parsing - parameters*/
    if (argc < 3)
    {
        fprintf(stderr, "too few parameters\n");
        fprintf(stderr, "crop an image of size WxH with X0,Y0 (integers) as\
                upper left corner,\nthe cropped region may fall outside of the\
                image.\n");
        fprintf(stderr, "   usage: %s in out X0 Y0 W H\n", argv[0]);
        return 1;
    }

    if (argc >= 4) x0 = atoi (argv[3]);
    if (argc >= 5) y0 = atoi (argv[4]);
    if (argc >= 6) w  = atoi (argv[5]);
    if (argc >= 7) h  = atoi (argv[6]);
    assert (w>0 && h>0);

    int nc, nr, nch;
    float *in = iio_read_image_float_vec(argv[1], &nc, &nr, &nch);

    float *out = malloc(w*h*nch*sizeof*out);

    // set default
    for (int y=0; y < w*h*nch; y++) out[y] = NAN;

    // copy
    for (int y = y0; y < y0 + h; y++)
        for (int x = x0; x < x0 + w; x++)
            for (int c = 0; c < nch; c++)
            {
                if (x>=0 && y>=0 && x<nc && y<nr)
                    out[c + ( (x-x0) + (y-y0)*w )*nch] = in[c+(x+y*nc)*nch];
            }

    iio_save_image_float_vec(argv[2], out, w, h, nch);

    free(in);
    free(out);

    return 0;
}
