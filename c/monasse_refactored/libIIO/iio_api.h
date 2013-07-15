#ifndef _IIO_API_H
#define _IIO_API_H

//
// FLOAT API for 2D images
//

float *iio_read_image_float(const char *fname, int *w, int *h);
// x[i+j*w]

float *iio_read_image_float_vec(const char *fname, int *w, int *h, int *nc);
// x[(i + j*w)*nc + l]

float (*iio_read_image_float_2vec(const char *fname, int *w, int *h)[2];
float (*iio_read_image_float_3vec(const char *fname, int *w, int *h)[3];
float (*iio_read_image_float_4vec(const char *fname, int *w, int *h)[4];
float (*iio_read_image_float_rgb(const char *fname, int *w, int *h)[3];
float (*iio_read_image_float_rgba(const char *fname, int *w, int *h)[4];
// x[i+j*w][l]
// (the "rgb" and "rgba" functions may re-order the channels according to file
// metadata)

float **iio_read_image_float_matrix(const char *fname, int *w, int *h);
// x[j][i]

float (**iio_read_image_float_matrix_2vec(const char *fnam, int *w, int *h))[2];
float (**iio_read_image_float_matrix_3vec(const char *fnam, int *w, int *h))[3];
float (**iio_read_image_float_matrix_4vec(const char *fnam, int *w, int *h))[4];
float (**iio_read_image_float_matrix_rgb(const char *fnam, int *w, int *h))[3];
float (**iio_read_image_float_matrix_rgba(const char *fnam, int *w, int *h))[4];
// x[j][i][channel]
// (the "rgb" and "rgba" functions may re-order the channels according to file
// metadata)



//
// FLOAT API for 3D images
//

float *iio_read_3d_image_float(const char *fname, int *w, int *h, int *d);
// x[i+(j+k*h)*w]

float *iio_read_3d_image_float_vec(char *f, int *w, int *h, int *d, int *nc);
// x[(i + (j+k*h)*w)*nc + l]

float (*iio_read_3d_image_float_2vec(char *fname, int *w, int *h, int *d)[2];
float (*iio_read_3d_image_float_3vec(char *fname, int *w, int *h, int *d)[3];
float (*iio_read_3d_image_float_4vec(char *fname, int *w, int *h, int *d)[4];
float (*iio_read_3d_image_float_rgb (char *fname, int *w, int *h, int *d)[3];
float (*iio_read_3d_image_float_rgba(char *fname, int *w, int *h, int *d)[4];
// x[i+(j+k*h)*w][l]
// (the "rgb" and "rgba" functions may re-order the channels according to file
// metadata)

float ***iio_read_3d_image_float_matrix(char *fname, int *w, int *h, int *d);
// x[k][j][i]

float (***iio_read_3d_image_float_matrix_2vec(char*f, int*w, int*h, int*d))[2];
float (***iio_read_3d_image_float_matrix_3vec(char*f, int*w, int*h, int*d))[3];
float (***iio_read_3d_image_float_matrix_4vec(char*f, int*w, int*h, int*d))[4];
float (***iio_read_3d_image_float_matrix_rgb (char*f, int*w, int*h, int*d))[3];
float (***iio_read_3d_image_float_matrix_rgba(char*f, int*w, int*h, int*d))[4];
// x[k][j][i][channel]
// (the "rgb" and "rgba" functions may re-order the channels according to file
// metadata)



//
// So far, we have seen the functions for reading multi-channel "2D" images,
// whose samples are of type "FLOAT", from a "NAMED FILE".  The rest of the
// API consists in variations of these functions, changing the quoted options
// of the previous sentence in all possible ways.  For instance, there are
// the following variations of the function "iio_read_image_float":
//
// 	// read 3D image into floats, from a named file
// 	float *iio_read_3d_image_float(char *fname, int *w, int *h, int *d);
//
//	// read 2D image into bytes, from a named file
// 	uint8_t *iio_read_image_uint8(char *fname, int *w, int *h);
//
// 	// read 2D image into floats, from an open stream
// 	float *iio_read_image_float_f(FILE *f, int *w, int *h);
//
// 	// read 2D image into floats, from memory
// 	float *iio_read_image_float_m(void *f, size_t s, int *w, int *h);
//
//	// general forms of these functions
//	TYPE *iio_read{,_3d,_4d,_nd}_image_TYPE{,_f,_m}(...);
//	TYPE **iio_read{,_3d,_4d,_nd}_image_TYPE_matrix{,_f,_m}(...);
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_2vec{,_f,_m}(...))[2];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_3vec{,_f,_m}(...))[3];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_4vec{,_f,_m}(...))[4];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_rgb{,_f,_m}(...))[3];
//	TYPE (**iio_read{,_3d,_4d,_nd}_image_TYPE_matrix_rgba{,_f,_m}(...))[4];
//

// All these functions are boring  variations, and they are defined at the
// end of this file.  More interesting are the two following general
// functions:

void *iio_read_image_numbers_as_they_are_stored(char *fname, int *w, int *h,
		int *samples_per_pixel, int *sample_size,
		bool *ieeefp_samples, bool *signed_samples);

void *iio_read_image_numbers_in_requested_format(char *fname, int *w, int *h,
		int *samples_per_pixel, int requested_sample_size,
		bool requested_ieefp, bool requested_sign);

// These two general functions have the usual versions for nD images and
// streams.  There exist also the following truly general functions, that
// read images of any dimension:

void *iio_read_nd_image_as_stored(char *fname,
		int *dimension, int *sizes,
		int *samples_per_pixel, int *sample_size,
		bool *ieeefp_samples, bool *signed_samples);

void *iio_read_nd_image_as_desired(char *fname,
		int *dimension, int *sizes,
		int *samples_per_pixel, int desired_sample_size,
		bool desired_ieeefp_samples, bool desired_signed_samples);







// versions of the API using FILE*
// ...

// rest of the API: chars and ints
// ...

// _4d versions, etc
// ...

#define IIO_MAX_DIMENSION 5

//
void *iio_read_image_raw(const char *fname,
		int *dimension,
		int sizes[IIO_MAX_DIMENSION],
		int *pixel_dimension,
		size_t *sample_integer_size,
		size_t *sample_float_size,
		int *metadata_id);



#endif//_IIO_API_H
