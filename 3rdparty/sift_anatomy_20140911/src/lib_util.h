#ifndef _LIB_UTIL_H_
#define _LIB_UTIL_H_ 

#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT2
    #define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_LN2
    #define M_LN2 0.693147180559945309417
#endif

// Allocate memory or abord on failure
void* xmalloc(size_t size);

// Reallocate memory of abort on failure
void* xrealloc(void* p, size_t size);

// Free memory allocated by xmalloc or xrealloc.
void xfree(void* p);

// Write a formatted message to standard error and abort the program, for
// example, fatal_error("Failed to open file %s for writing", filename);
void fatal_error(const char* format, ...);

// Write formatted message to standard error for debugging.
void debug(const char* format, ...);

// Linearly rescale input such that its values are in the range [0, 250].
void linear_conversion(const float *in, float *out, int l);

// Find the maximum value of an array.
float array_max(const float* array, int length);

// Find the minimum value of an array.
float array_min(const float* array, int length);

// Find the maximum value of an array.
float find_array_max(const float* array, int length, int *position);

// Find the minimum value of an array.
float find_array_min(const float* array, int length, int *position);

// Find the two minimal values in an array.
void find_array_two_min(const float* array, int length, float* minA, float* minB, int* iA, int* iB);

// L2 norm of an array.
float array_l2norm(const float* array, int length);

// Compute the SQUARE of the euclidean distance
float euclidean_distance_square(const float* a1, const float* a2, int length);

// Compute the euclidean distance
float euclidean_distance(const float* a1, const float* a2, int length);

// Compute the x modulus y
float modulus(float x, float y);

// Multiply the rotation matric R_alpha with [x,y]^T
void apply_rotation(float x, float y, float *rx, float *ry, float alpha);

#endif // _LIB_UTIL_H_
