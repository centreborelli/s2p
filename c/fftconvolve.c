#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>
#include "iio.h"

#define FORI(n) for (int i=0; i<(n); i++)

// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued image
static void fft_2dfloat(float complex *fx,
      const float *x, const int w, const int h)
{
    float complex *a = fftwf_malloc(w*h*sizeof*a);

    fftwf_plan p = fftwf_plan_dft_2d(h, w, a, fx,
            FFTW_FORWARD, FFTW_ESTIMATE);

    FORI(w*h) a[i] = x[i]; // complex assignment!
    fftwf_execute(p);

    fftwf_destroy_plan(p);
    fftwf_free(a);
    fftwf_cleanup();
}

// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2dfloat(float *ifx,
        const float complex *fx, const int w, const int h)
{
    float complex *a = fftwf_malloc(w*h*sizeof*a);
    float complex *b = fftwf_malloc(w*h*sizeof*b);

    fftwf_plan p = fftwf_plan_dft_2d(h, w, a, b,
    					FFTW_BACKWARD, FFTW_ESTIMATE);

    FORI(w*h) a[i] = fx[i];
    fftwf_execute(p);
    float scale = 1.0/(w*h);
    FORI(w*h) {
        float complex z = b[i] * scale;
        ifx[i] = crealf(z);
        if (cimagf(z) > 0.001) printf("inverse %d %f\n", i, cimagf(z));
        //assert(cimagf(z) < 0.001);
    }
    fftwf_destroy_plan(p);
    fftwf_free(a);
    fftwf_free(b);
    fftwf_cleanup();
}




//fft convolution of input by a filter (MTF ctr at 0)
void fftconvolve(const float* in, const float* filt, float* out, const int w,
        const int h)
{
    // alloc temporary storage for the fft's
    float complex *fin  = fftwf_malloc(w*h*sizeof*fin);
    float complex *fout = fftwf_malloc(w*h*sizeof*fout);

    // fft
    fft_2dfloat(fin,in,w,h);

    // product
    for (int i=0; i<w*h; i++)
       fout[i]= fin[i] * (filt[i] + 0*I);

    // ifft
    ifft_2dfloat(out,fout,w,h);

    fftwf_free(fin);
    fftwf_free(fout);
}


int main (int argc, char **argv)
{
    /* parameter parsing - parameters*/
    if (argc<3)
    {
        fprintf (stderr, "too few parameters\n");
        fprintf (stderr, "fft convolution of input by a filter (MTF ctr at 0)\n");
        fprintf (stderr, "   usage: %s filter [input [output]]\n",argv[0]);
        return 1;
    }

    //read the parameters
    int i = 1;
    char* filt_im = argv[i]; i++;
    char* in_im   = (argc>i)? argv[i]: "-";       i++;
    char* out_im  = (argc>i)? argv[i]: "-";       i++;

    // read input
    int fw,fh,w,h,ch;
    float *filt = iio_read_image_float_split(filt_im, &fw, &fh, &ch);
    float *in = iio_read_image_float_split(in_im, &w, &h, &ch);
    if (h!=fh || w!=fw) {
        fprintf (stderr, "filter and input sizes don't match\n");
        return 1;
    }
    // allocate output
    float *out = malloc(w*h*ch*sizeof*out);

    // call the convolution for each channel
    for (int i=0;i<ch;i++)
        fftconvolve(in + w*h*i, filt, out + w*h*i, w,h);

    // save output
    iio_save_image_float_split(out_im, out, w, h, ch);

    free(filt);
    free(in);
    free(out);

    return 0;
}
