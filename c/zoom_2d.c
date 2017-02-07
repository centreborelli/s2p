#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "iio.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/*
This function modifies the frequencies of the dft vector given as input, in a
very simple way. If the input size is bigger than the output, it cut the highest
frequencies. If the input size is smaller it add zeros at the highest
frequencies.

Parameters:
-fftw_complex *in: address of the input array. It's an array of complex numbers
-fftw_complex *out: address of the output array.
-int n_in: size of the input array
-int n_out: size of the output array

The complex arrays are obtained with the dft_r2c_1d functions of libfftw, so
they contain only one half of the DFT coefficients (because of hermitian symmetry)
*/
void cut_or_add_freq(fftw_complex *in, fftw_complex *out, int n_in, int n_out)
{
    if (n_out <= n_in)
        // In that case we keep the first n_out frequencies
    {
        for (int k=0; k<n_out; k++)
            out[k] = in[k];

        // The last dft coefficient has to be real (because it's the n/2
        // coefficient of the "real" signal
        out[n_out-1] = crealf(out[n_out-1]);
    }
    else
        // In that case we keep all the frequencies and put a zero for the new ones
    {
        //first copy the in vector
        for (int k=0; k<n_in; k++)
            out[k] = in[k];

        //then put zeros in the remaining freq
        for (int k=n_in; k<n_out; k++)
            out[k] = 0;
    }
}


/*
This function applies a 1d-zoom to the mono-channel image stored at
address 'in', and write the result at address 'out'.
The width and height of input image are provided in parameters 'w_in' and
'h'. The image is compressed or stretched along the horizontal axis to get a
width of 'w_out'.

The transformed image is computed line by line. Each line of the input image
is compressed or stretched with zoom factor equal to w_out/w_in. The basic flow
of the algorithm is: for each line, symmetrize the signal by replicating it,
compute its DFT, then filter it with a gaussian kernel (if it's a compression),
and then erase the high frequency coefficients, or add zeros if it's a
stretching. Then compute the inverse DFT of this resized vector and store the
result.
*/

void image_zoom_1d(float *in, float *out, int w_in, int w_out, int h)
{
    // If w_out < w_in, we need to filter the signal
    int filter_needed = (w_out < w_in);
    float s = 0; // Gaussian std deviation
    if (filter_needed)
        s = 0.8*sqrt((w_in/w_out)*(w_in/w_out)-1);

    // Length of the 1d signal
    int n_in = 2*w_in;
    int n_out = 2*w_out;

    // Memory allocation
    double *line_in = fftw_malloc(sizeof(double) * n_in);
    double *line_out = fftw_malloc(sizeof(double) * n_out);
    fftw_complex *dft_in = fftw_malloc(sizeof(fftw_complex) * (w_in+1));
    fftw_complex *dft_out = fftw_malloc(sizeof(fftw_complex) * (w_out+1));

    // Plans computation
    fftw_plan p_fwd = fftw_plan_dft_r2c_1d(n_in, line_in, dft_in, FFTW_ESTIMATE);
    fftw_plan p_back = fftw_plan_dft_c2r_1d(n_out, dft_out, line_out, FFTW_ESTIMATE);

    // Compute the output image row by row :
    for (int row = 0; row < h; row++)
    {
        // Copy and duplicate the current row, with symmetry towards the last term
        for (int i=0; i<w_in; i++)
        {
            const double x = in[row*w_in+i];
            line_in[i] = isfinite(x) ? x : 0;
            line_in[n_in-1-i] = isfinite(x) ? x : 0;
        }

        // Compute the DFT
        fftw_execute(p_fwd);

        if (filter_needed)
        {
            // Filter the signal : each coefficient of the DFT is multiplied by
            // exp(-(1/2)*s*s*omega*omega) where omega is the pulsation
            // This is equivalent to convolve the trigonometric polynomial
            // with a gaussian of parameter s
            for (int k=1; k<=w_in; k++)
                dft_in[k] *= exp(-0.5*s*s*(k*2.0*M_PI/n_in)*(k*2.0*M_PI/n_in));
        }

        // Normalize the DFT
        for (int i=0; i<=w_in; i++)
            dft_in[i] /= n_in;

        // Frequency modification
        cut_or_add_freq(dft_in, dft_out, w_in+1, w_out+1);

        // Compute the iDFT to get the signal zoomed
        fftw_execute(p_back);

        // Copy zoomed row in the output image (just the half that we need)
        for (int i=0; i<w_out; i++)
            out[row*w_out+i] = line_out[i];
    }

    // Free memory
    fftw_destroy_plan(p_fwd);
    fftw_destroy_plan(p_back);
    fftw_free(line_in);
    fftw_free(dft_in);
    fftw_free(dft_out);
    fftw_free(line_out);
}


/*
This function transposes an image
*/
void transpose(float *input, float *output, int width, int height)
{
    for (int row=0; row<height; row++)
        for (int col=0; col<width; col++)
            output[col*height+row] = input[row*width+col];
}


/*
This function applies a 2d-zoom to the mono-channel image stored at address
'in', and write the result at address 'out'.
The width and height of input image are provided in parameters 'w_in' and
'h_in'. The image is compressed or stretched along the horizontal axis to get
a width of 'w_out' and a height of 'h_out'.

The algorithm uses the image_zoom_1d function to do the horizontal zoom, then
transposes the image, applies the same function again, and transpose to get the
final result.
*/
void image_zoom_2d(float *in, float *out, int w_in, int h_in, int w_out, int h_out)
{

    float *tmp_out_1 = malloc(sizeof(float)*w_out*h_in);
    float *tmp_out_2 = malloc(sizeof(float)*w_out*h_in);
    float *tmp_out_3 = malloc(sizeof(float)*w_out*h_out);

    image_zoom_1d(in, tmp_out_1, w_in, w_out, h_in);
    transpose(tmp_out_1, tmp_out_2, w_out, h_in);
    image_zoom_1d(tmp_out_2, tmp_out_3, h_in, h_out, w_out);
    transpose(tmp_out_3, out, h_out, w_out);

    // Free memory
    free(tmp_out_1);
    free(tmp_out_2);
    free(tmp_out_3);
}

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
