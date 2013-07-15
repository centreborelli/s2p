extern "C" {
#include "libIIO/iio.h"
}

float *read_png_f32_gray(const char *fname, size_t *nxp, size_t *nyp)
{
    int w,h;
    float *out = iio_read_image_float(fname, &w, &h);
    *nxp = w;
    *nyp = h;
    return out;
}

float *read_png_f32(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp) 
{
    int w,h,d;
    float *out = iio_read_image_float_split(fname, &w, &h, &d);
    *nxp = w;
    *nyp = h;
    *ncp = d;
    return out;
}

int write_png_f32(const char *fname, const float *data, size_t nx, size_t ny, size_t nc)
{
    iio_save_image_float_split((char*)fname, (float*)data, (int) nx, (int) ny, (int) nc);
    return 0;
}




float *read_tiff_f32_gray(const char *fname, size_t *nx, size_t *ny)
{
    int w,h;
    float *out = iio_read_image_float(fname, &w, &h);
    *nx = w;
    *ny = h;
    return out;
}

int write_tiff_f32(const char *fname, const float *data, size_t nx, size_t ny, size_t nc)
{
    iio_save_image_float_split((char*)fname, (float*)data, (int) nx, (int) ny, (int) nc);
    return 0;
}
