#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>
#include "iio.h"

#define FORI(n) for(int i=0; i<(n); i++)

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


static int idx2freq(const int N, const int i)
{
    int hN=N-((int)N/2);
    if(i<hN) return i;
    else return i-N;
}

static int freq2idx(const int N, const int w)
{
    // TODO consider spectrum folding
    int hN=N-((int)N/2);
    assert(w<hN);
    assert(w>=hN-N);
    if (w>=0) return w;
    else return N+w;
}


/**
 * Performs the zero padding on a fourier transform.
 * The size of the output can be larger (add zeros)
 * or smaller that the input (remove frequencyes).
 * This function is safe because it preserves
 * the conjugated symetry of the spectrum,
 * either by splitting the last frequency (even input),
 * or by adding the conjugates of the last freqency (even output).
 * */
static void zeropadding_fourier_safe(const float complex* in,
        const int w, const int h,
        float complex* out,
        const int ow, const int oh)
{
    // signed min and max freq of the input image (w is hor, h is ver)
    const int mhw=-w/2;
    const int hw = w+mhw;
    const int mhh=-h/2;
    const int hh = h+mhh;

    // signed min frequ of the output image (w is hor, h is ver)
    const int mhow=-ow/2;
    const int mhoh=-oh/2;


    // go over all pixels of the output image
    for(int i=0;i<ow*oh;i++)  out[i] = 0+0*I;
    for(int j=0; j<oh;j++)
        for(int i=0; i<ow;i++)
        {
            // for each point determine the signed frequency
            const int fhor = idx2freq(ow,i);
            const int fver = idx2freq(oh,j);

            // check if the this freq is represented in the input
            if(fhor < hw && fhor >= mhw &&
                    fver < hh && fver >= mhh ){
                // compute the indices on the input image and copy the value
                const int k = freq2idx(w,fhor);
                const int l = freq2idx(h,fver);
                out[i+j*ow] = in[k+l*w];

                // IF output is smaller than input and output is even sized
                // so we collapse (add) the highest freq with its conjugate
                const int kk = freq2idx(w,-fhor);
                const int ll = freq2idx(h,-fver);
                if(fver==mhoh && oh%2==0 && oh < h){
                    out[i+j*ow] += in[k+ll*w];
                }

                if(fhor==mhow && ow%2==0 && ow < w){
                    out[i+j*ow] += in[kk+l*w];
                    if(fver==mhoh && oh%2==0 && oh < h){
                        out[i+j*ow] += in[kk+ll*w];
                    }
                }

                // IF output is larger than the input and the input is even
                // so we split (divide) the frequency at freq: mhh and mhw
                const int ii = freq2idx(ow,-fhor);
                const int jj = freq2idx(oh,-fver);
                if(oh > h && h%2==0 && fver==mhh ){
                    out[i+j*ow]  = creal(in[k+l*w]/2)+0*I;
                    out[i+jj*ow] = creal(in[k+l*w]/2)+0*I;
                }
                if(ow > w && w%2==0 && fhor==mhw ){
                    out[i+j*ow]  = creal(in[k+l*w]/2)+0*I;
                    out[ii+j*ow] = creal(in[k+l*w]/2)+0*I;
                    if(oh > h && h%2==0 && fver==mhh ){
                        out[i+j*ow]  = creal(in[k+l*w]/4)+0*I;
                        out[ii+jj*ow]= creal(in[k+l*w]/4)+0*I;
                        out[i+jj*ow] = creal(in[k+l*w]/4)+0*I;
                        out[ii+j*ow] = creal(in[k+l*w]/4)+0*I;
                    }
                }
            }
        }
}



void show(float* in, int nc, int nr, int nch) {
    for(int c=0;c<nch;c++)
    {
        for(int y=0;y<nr;y++)
        {
            for(int x=0;x<nc;x++)
                printf("%f ", in[c+(x+y*nc)*nch]);
            printf("\n");
        }
        printf("\n");
    }
}


void zoom_zeropadding(const float* in, const int w, const int h,
        float* out, const int ow, const int oh)
{
    // alloc temporary storage for the fft's
    float complex *fin  = fftwf_malloc(w*h*sizeof*fin);
    float complex *fout = fftwf_malloc(ow*oh*sizeof*fout);

    // fft
    fft_2dfloat(fin,in,w,h);
    //   show(fin,w,h,2);

    // zeropadding
    zeropadding_fourier_safe(fin, w, h, fout, ow, oh);
    //   show(fout,ow,oh,2);

    // ifft
    ifft_2dfloat(out,fout,ow,oh);

    // scale values
    double scale= ((double) ow*oh)/((double) w*h);
    for (int i=0; i<ow*oh; i++) out[i]*=scale;

    fftwf_free(fin);
    fftwf_free(fout);
}


int main (int argc, char **argv)
{
    /* parameter parsing - parameters*/
    if (argc<3)
    {
        fprintf (stderr, "too few parameters\n");
        fprintf (stderr, "zoom by zero padding\n");
        fprintf (stderr, "\tusage: %s outsize [input [output]]\n",argv[0]);
        return 1;
    }

    //read the parameters
    int i = 1;
    char* tgt_im = argv[i]; i++;
    char* in_im  = (argc>i)? argv[i]: "-";       i++;
    char* out_im = (argc>i)? argv[i]: "-";       i++;

    // read input
    int outw,outh;
    int w,h,wh;
    float *tmp = iio_read_image_float_split(tgt_im, &outw, &outh, &wh);
    free(tmp);
    float *in = iio_read_image_float_split(in_im, &w, &h, &wh);
    // allocate output
    float *out = malloc(outw*outh*wh*sizeof*out);

    // call the algorithm
    for (int i=0; i<wh; i++)
        zoom_zeropadding(in + w*h*i,w,h, out + outw*outh*i, outw, outh);

    // save output
    iio_save_image_float_split(out_im, out, outw, outh, wh);

    free(in);
    free(out);

    return 0;
}
