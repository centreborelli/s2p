#ifdef __cplusplus
extern "C" {
#endif

#define IO_PNG_VERSION 0.20100727

/* io_png.c */
unsigned char *read_png_u8(const char *fname, size_t *nx, size_t *ny, size_t *nc);
unsigned char *read_png_u8_rgb(const char *fname, size_t *nx, size_t *ny);
unsigned char *read_png_u8_gray(const char *fname, size_t *nx, size_t *ny);
float *read_png_f32(const char *fname, size_t *nx, size_t *ny, size_t *nc);
float *read_png_f32_rgb(const char *fname, size_t *nx, size_t *ny);
float *read_png_f32_gray(const char *fname, size_t *nx, size_t *ny);
int write_png_u8(const char *fname, const unsigned char *data, size_t nx, size_t ny, size_t nc);
int write_png_f32(const char *fname, const float *data, size_t nx, size_t ny, size_t nc);

#ifdef __cplusplus
}
#endif
