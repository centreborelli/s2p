#include <stdlib.h>
#include <stdio.h>
#include "iio.h"
#include "zoom.h"

int main(int c, char *v[])
{
    if (c < 4)
    {
	printf("Usage: %s <input file> <output file> <output width> <output height>\n", v[0]);
        return EXIT_FAILURE;
    }

    // Parameters loading
    char *file_in = v[1];
    char *file_out = v[2];
    int w_out= atoi(v[3]);
    int h_out= atoi(v[4]);

    // Input image loading
    int w, h, pd;
    float *input = iio_read_image_float_split(file_in, &w, &h, &pd);
    printf("Width : %d, ", w);
    printf("Height : %d, ", h);
    printf("Channels : %d\n", pd);


    // Memory allocations
    int n = w * h;
    int n_out = w_out * h_out;
    float *output = malloc(n_out*pd*sizeof(float));
    int off_in = 0;
    int off_out = 0;

    // Image processing, channel by channel
    for (int channel=0; channel<pd; channel++)
    {
        // Apply the transformation on the current channel
        image_zoom_2d(input + off_in, output + off_out, w, h, w_out, h_out);
        off_in += n;
        off_out += n_out;
    }

    // Save and free memory
    iio_save_image_float_split(file_out, output, w_out, h_out, pd);
    free(input);
    free(output);

    return EXIT_SUCCESS;
}
