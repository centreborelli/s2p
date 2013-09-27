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
    float *input = iio_read_image_float_vec(file_in, &w, &h, &pd);
    printf("Width : %d, ", w);
    printf("Height : %d, ", h);
    printf("Channels : %d\n", pd);


    // Memory allocations
    float *output = malloc(w_out*h_out*pd*sizeof(float));
    float *in_ch = malloc(w*h*sizeof(float));
    float *out_ch = malloc(w_out*h_out*sizeof(float));

    // Image processing, channel by channel
    for (int channel=0; channel<pd; channel++)
    {
        // Copy the current channel
        for (int pix=0; pix<w*h; pix++)
            in_ch[pix] = input[pd*pix+channel];

        // Apply the transformation on the current channel
        image_zoom_2d(in_ch, out_ch, w, h, w_out, h_out);

        // Copy the result in the out image :
        for (int pix=0; pix<w_out*h_out; pix++)
            output[pd*pix+channel] = out_ch[pix];
    }

    iio_save_image_float_vec(file_out, output, w_out, h_out, pd);

    // Free memory
    free(in_ch);
    free(out_ch);
    free(input);
    free(output);

    return EXIT_SUCCESS;
}
