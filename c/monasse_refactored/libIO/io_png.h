#ifndef _IO_PNG_H
#define _IO_PNG_H

#ifdef __cplusplus
extern "C" {
#endif

#define IO_PNG_VERSION "0.20100917"

#include <stddef.h>

unsigned char *read_png_u8(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp);
unsigned char *read_png_u8_rgb(const char *fname, size_t *nxp, size_t *nyp);
unsigned char *read_png_u8_gray(const char *fname, size_t *nxp, size_t *nyp);
float *read_png_f32(const char *fname, size_t *nxp, size_t *nyp, size_t *ncp);
float *read_png_f32_rgb(const char *fname, size_t *nxp, size_t *nyp);
float *read_png_f32_gray(const char *fname, size_t *nxp, size_t *nyp);
int write_png_u8(const char *fname, const unsigned char *data, size_t nx, size_t ny, size_t nc);
int write_png_f32(const char *fname, const float *data, size_t nx, size_t ny, size_t nc);

#ifdef __cplusplus
}
#endif

#endif /* !_IO_PNG_H */
