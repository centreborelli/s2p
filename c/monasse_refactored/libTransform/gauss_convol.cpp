#include "gauss_convol.h"
#include <cassert>
#include <cmath>

/* Gaussian convolution kernels are truncated at this many sigmas from
   the center.  While it is more efficient to keep this value small,
   experiments show that for consistent scale-space analysis it needs
   a value of about 3.0, at which point the Gaussian has fallen to
   only 1% of its central value.  A value of 2.0 greatly reduces
   keypoint consistency, and a value of 4.0 is better than 3.0.
*/
static const flnum GaussTruncate = 4.0;


/* --------------------------- Blur image --------------------------- */

/* Same as ConvBuffer, but implemented with loop unrolling for increased
   speed.  This is the most time intensive routine in keypoint detection,
   so deserves careful attention to efficiency.  Loop unrolling simply
   sums 5 multiplications at a time to allow the compiler to schedule
   operations better and avoid loop overhead.  This almost triples
   speed of previous version on a Pentium with gcc.
*/
static void ConvBufferFast(flnum *buffer, flnum *kernel, int rsize, int ksize)
{
    int i;
    flnum sum, *bp, *kp, *endkp;

    for (i = 0; i < rsize; i++) {
      sum = 0.0;
      bp = &buffer[i];
      kp = &kernel[0];
      endkp = &kernel[ksize];

      /* Loop unrolling: do 5 multiplications at a time. */
      while (kp + 4 < endkp) {
        sum += bp[0] * kp[0] +  bp[1] * kp[1] + bp[2] * kp[2] +
               bp[3] * kp[3] +  bp[4] * kp[4];
        bp += 5;
        kp += 5;
      }
      /* Do 2 multiplications at a time on remaining items. */
      while (kp + 1 < endkp) {
        sum += bp[0] * kp[0] +  bp[1] * kp[1];
        bp += 2;
        kp += 2;
      }
      /* Finish last one if needed. */
      if (kp < endkp)
        sum += *bp * *kp;

      buffer[i] = sum;
    }
}

/* Convolve image with the 1-D kernel vector along image rows.  This
   is designed to be as efficient as possible.  Pixels outside the
   image are set to the value of the closest image pixel.
*/
static void ConvHorizontal(LWImage<flnum>& image, flnum *kernel, int ksize)
{
    flnum buffer[8000];
    
    const int rows = image.h;
    const int cols = image.w;
    const int halfsize = ksize / 2;
    assert(cols + ksize < 8000); /*TANG: this will give a limit of image size*/

    for(int comp = 0; comp < image.comps; comp++) {
        const int deltaComp = comp*image.stepComp();
        for (int r = 0; r < rows; r++) {
            /* Copy the row into buffer with pixels at ends replicated for
            half the mask size.  This avoids need to check for ends
            within inner loop. */
            for (int i = 0; i < halfsize; i++)
                buffer[i] = image.pixel(0,r)[deltaComp];
            for (int i = 0; i < cols; i++)
                buffer[halfsize + i] = image.pixel(i,r)[deltaComp]; 
            for (int i = 0; i < halfsize; i++)
                buffer[halfsize + cols + i] = image.pixel(cols-1,r)[deltaComp];

            ConvBufferFast(buffer, kernel, cols, ksize);
            for (int c = 0; c < cols; c++)
                image.pixel(c,r)[deltaComp] = buffer[c];
        }
    }
}

/* Same as ConvHorizontal, but apply to vertical columns of image.
*/
static void ConvVertical(LWImage<flnum>& image, flnum *kernel, int ksize)
{
    flnum buffer[8000];

    const int rows = image.h;
    const int cols = image.w;
    const int halfsize = ksize / 2;
    assert(rows + ksize < 8000);  /*TANG: this will give a limit of image size*/

    for(int comp = 0; comp < image.comps; comp++) {
        const int deltaComp = comp*image.stepComp();
        for (int c = 0; c < cols; c++) {
            for (int i = 0; i < halfsize; i++)
                buffer[i] = image.pixel(c,0)[deltaComp];
            for (int i = 0; i < rows; i++)
                buffer[halfsize + i] = image.pixel(c,i)[deltaComp];
            for (int i = 0; i < halfsize; i++)
                buffer[halfsize + rows + i] = image.pixel(c,rows-1)[deltaComp];

            ConvBufferFast(buffer, kernel, rows, ksize);
            for (int r = 0; r < rows; r++)
                image.pixel(c,r)[deltaComp] = buffer[r];
        }
    }
}

/* Convolve image with a Gaussian of width sigma and store result back
   in image.   This routine creates the Gaussian kernel, and then applies
   it sequentially in horizontal and vertical directions.
*/
void gauss_convol(LWImage<flnum>& image, flnum sigma)
{
    /*float x, kernel[100], sum = 0.0;*/
    flnum x, kernel[1000], sum = 0.0; /*TANG*/
    int ksize, i;

    /* The Gaussian kernel is truncated at GaussTruncate sigmas from
       center.  The kernel size should be odd.
    */
    ksize = (int)(2.0 * GaussTruncate * sigma + 1.0);
    /*ksize = MAX(3, ksize);*/    /* Kernel must be at least 3. */
    ksize = ksize > 3 ? ksize : 3;
    if (ksize % 2 == 0)       /* Make kernel size odd. */
      ksize++;
    /*assert(ksize < 100);*/
    assert(ksize < 1000); /*TANG*/
    
    /* Fill in kernel values. */
    for (i = 0; i <= ksize; i++) {
      x = static_cast<flnum>(i - ksize/2);
      kernel[i] = exp(- x * x / (2 * sigma * sigma));
      sum += kernel[i];
    }
    /* Normalize kernel values to sum to 1.0. */
    for (i = 0; i < ksize; i++)
      kernel[i] /= sum;

    ConvHorizontal(image, kernel, ksize);
    ConvVertical(image, kernel, ksize);
} 
