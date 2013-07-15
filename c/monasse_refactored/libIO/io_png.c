/*
 * Copyright (c) 2010, Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under, at your option, the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version, or
 * the terms of the simplified BSD license.
 *
 * You should have received a copy of these licenses along this
 * program. If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * @file io_png.c
 * @brief PNG read/write simplified interface
 *
 * This is a front-end to libpng, with routines to:
 * @li read a PNG file as a deinterlaced 8bit integer or float array
 * @li write a 8bit integer or float array to a PNG file
 *
 * Multi-channel images are handled : grey, grey+alpha, rgb and
 * rgb+alpha, as well as on-the-fly color model conversion.
 *
 * @todo handle lossless 16bit data
 * @todo add a test suite
 * @todo internally handle RGB/gray conversion in read_png_raw()
 * @todo handle deinterlacing as a libpng transform function
 *
 * @author Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* option to use a local version of the libpng */
#ifdef WITH_LOCAL_LIBPNG
#include "png.h"
#else
#include <png.h>
#endif

/* ensure consistency */
#include "io_png.h"

#define PNG_SIG_LEN 4

/* internal only datatype identifiers */
#define IO_PNG_U8  0x0001       /* 8bit unsigned integer */
#define IO_PNG_F32 0x0002       /* 32bit float */

/* string tag inserted into the binary, helps tracking versions */
char _io_png_tag[] = "using io_png " IO_PNG_VERSION;

/*
 * READ
 */

/**
 * @brief internal function used to cleanup the memory when
 * png_read_raw() fails
 *
 * @param fp file pointer to close, ignored if NULL
 * @param png_ptr_p, info_ptr_p, pointers to PNG structure pointers,
 *        ignored if NULL
 * @return NULL
 */
static void *read_png_abort(FILE * fp,
                            png_structp * png_ptr_p, png_infop * info_ptr_p)
{
    png_destroy_read_struct(png_ptr_p, info_ptr_p, NULL);
    if (NULL != fp && stdin != fp)
        (void) fclose(fp);
    return NULL;
}

/**
 * @brief internal function used to read a PNG file into an array
 *
 * @todo don't loose 16bit info
 *
 * @param fname PNG file name, "-" means stdin
 * @param nxp, nyp, ncp pointers to variables to be filled
 *        with the number of columns, lines and channels of the image
 * @param transform a PNG_TRANSFORM to be added to the default read transforms
 * @param dtype identifier for the data type to be used for output
 * @return pointer to an allocated array of pixels,
 *         or NULL if an error happens
 */
static void *read_png_raw(const char *fname,
                          size_t * nxp, size_t * nyp, size_t * ncp,
                          int transform, int dtype)
{
    png_byte png_sig[PNG_SIG_LEN];
    png_structp png_ptr;
    png_infop info_ptr;
    png_bytepp row_pointers;
    png_bytep row_ptr;
    /* volatile : because of setjmp/longjmp */
    FILE *volatile fp = NULL;
    void *data = NULL;
    unsigned char *data_u8 = NULL;
    unsigned char *data_u8_ptr = NULL;
    float *data_f32 = NULL;
    float *data_f32_ptr = NULL;
    size_t size;
    size_t i, j, k;

    /* parameters check */
    if (NULL == fname || NULL == nxp || NULL == nyp || NULL == ncp)
        return NULL;
    if (IO_PNG_U8 != dtype && IO_PNG_F32 != dtype)
        return NULL;

    /* open the PNG input file */
    if (0 == strcmp(fname, "-"))
        fp = stdin;
    else if (NULL == (fp = fopen(fname, "rb")))
        return NULL;

    /* read in some of the signature bytes and check this signature */
    if ((PNG_SIG_LEN != fread(png_sig, 1, PNG_SIG_LEN, fp))
        || 0 != png_sig_cmp(png_sig, (png_size_t) 0, PNG_SIG_LEN))
        return read_png_abort(fp, NULL, NULL);

    /*
     * create and initialize the png_struct
     * with the default stderr and error handling
     */
    if (NULL == (png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                                  NULL, NULL, NULL)))
        return read_png_abort(fp, NULL, NULL);

    /* allocate/initialize the memory for image information */
    if (NULL == (info_ptr = png_create_info_struct(png_ptr)))
        return read_png_abort(fp, &png_ptr, NULL);

    /* set error handling */
    if (0 != setjmp(png_jmpbuf(png_ptr)))
        /* if we get here, we had a problem reading the file */
        /* free all of the memory associated with the png_ptr and info_ptr */
        return read_png_abort(fp, &png_ptr, &info_ptr);

    /* set up the input control using standard C streams */
    png_init_io(png_ptr, fp);

    /* let libpng know that some bytes have been read */
    png_set_sig_bytes(png_ptr, PNG_SIG_LEN);

    /*
     * set the read filter transforms, to get 8bit RGB whatever the
     * original file may contain:
     * PNG_TRANSFORM_STRIP_16      strip 16-bit samples to 8 bits
     * PNG_TRANSFORM_PACKING       expand 1, 2 and 4-bit
     *                             samples to bytes
     */
    transform |= (PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING);

    /* read in the entire image at once */
    png_read_png(png_ptr, info_ptr, transform, NULL);

    /* get image informations */
    *nxp = (size_t) png_get_image_width(png_ptr, info_ptr);
    *nyp = (size_t) png_get_image_height(png_ptr, info_ptr);
    *ncp = (size_t) png_get_channels(png_ptr, info_ptr);
    row_pointers = png_get_rows(png_ptr, info_ptr);

    /*
     * allocate the output data RGB array
     * deinterlace and convert png RGB RGB RGB 8bit to RRR GGG BBB
     * the image is deinterlaced layer after layer
     * this generic loop also works for one single channel
     */
    size = *nxp * *nyp * *ncp;
    switch (dtype)
    {
    case IO_PNG_U8:
        if (NULL == (data_u8 =
                     (unsigned char *) malloc(size * sizeof(unsigned char))))
            return read_png_abort(fp, &png_ptr, &info_ptr);
        data = (void *) data_u8;
        for (k = 0; k < *ncp; k++)
        {
            /* channel loop */
            data_u8_ptr = data_u8 + (size_t) (*nxp * *nyp * k);
            for (j = 0; j < *nyp; j++)
            {
                /* row loop */
                row_ptr = row_pointers[j] + k;
                for (i = 0; i < *nxp; i++)
                {
                    /* pixel loop */
                    *data_u8_ptr++ = (unsigned char) *row_ptr;
                    row_ptr += *ncp;
                }
            }
        }
        break;
    case IO_PNG_F32:
        if (NULL == (data_f32 = (float *) malloc(size * sizeof(float))))
            return read_png_abort(fp, &png_ptr, &info_ptr);
        data = (void *) data_f32;
        for (k = 0; k < *ncp; k++)
        {
            /* channel loop */
            data_f32_ptr = data_f32 + (size_t) (*nxp * *nyp * k);
            for (j = 0; j < *nyp; j++)
            {
                /* row loop */
                row_ptr = row_pointers[j] + k;
                for (i = 0; i < *nxp; i++)
                {
                    /* pixel loop */
                    *data_f32_ptr++ = (float) *row_ptr;
                    row_ptr += *ncp;
                }
            }
        }
        break;
    }

    /* clean up and free any memory allocated, close the file */
    (void) read_png_abort(fp, &png_ptr, &info_ptr);

    return data;
}

/**
 * @brief read a PNG file into a 8bit integer array
 *
 * The array contains the deinterlaced channels.
 * 1, 2 and 4bit images are converted to 8bit.
 * 16bit images are previously downscaled to 8bit.
 *
 * @todo don't downscale 16bit images.
 *
 * @param fname PNG file name
 * @param nxp, nyp, ncp pointers to variables to be filled with the number of
 *        columns, lines and channels of the image
 * @return pointer to an allocated unsigned char array of pixels,
 *         or NULL if an error happens
 */
unsigned char *read_png_u8(const char *fname,
                           size_t * nxp, size_t * nyp, size_t * ncp)
{
    /* read the image as unsigned char */
    return (unsigned char *) read_png_raw(fname, nxp, nyp, ncp,
                                          PNG_TRANSFORM_IDENTITY, IO_PNG_U8);
}

/**
 * @brief read a PNG file into a 8bit integer array, converted to RGB
 *
 * See read_png_u8() for details.
 */
unsigned char *read_png_u8_rgb(const char *fname, size_t * nxp, size_t * nyp)
{
    size_t nc;
    unsigned char *img;

    /* read the image */
    img = (unsigned char *) read_png_raw(fname, nxp, nyp, &nc,
                                         PNG_TRANSFORM_STRIP_ALPHA,
                                         IO_PNG_U8);
    if (NULL == img)
        /* error */
        return NULL;
    if (3 == nc)
        /* already RGB */
        return img;
    else
    {
        /* convert to RGB */
        unsigned char *ptr_r, *ptr_g, *ptr_b, *ptr_end;

        /* resize the image */
        img = realloc(img, 3 * *nxp * *nyp * sizeof(unsigned char));

        /* gray->RGB conversion */
        ptr_r = img;
        ptr_end = ptr_r + *nxp * *nyp;
        ptr_g = img + *nxp * *nyp;
        ptr_b = img + 2 * *nxp * *nyp;
        while (ptr_r < ptr_end)
        {
            *ptr_g++ = *ptr_r;
            *ptr_b++ = *ptr_r++;
        }
        return img;
    }
}

/**
 * @brief read a PNG file into a 8bit integer array, converted to gray
 *
 * See read_png_u8() for details.
 */
unsigned char *read_png_u8_gray(const char *fname, size_t * nxp, size_t * nyp)
{
    size_t nc;
    unsigned char *img;

    /* read the image */
    img = (unsigned char *) read_png_raw(fname, nxp, nyp, &nc,
                                         PNG_TRANSFORM_STRIP_ALPHA,
                                         IO_PNG_U8);
    if (NULL == img)
        /* error */
        return NULL;
    if (1 == nc)
        /* already gray */
        return img;
    else
    {
        /* convert to gray */
        unsigned char *ptr_r, *ptr_g, *ptr_b, *ptr_gray, *ptr_end;

        /*
         * RGB->gray conversion
         * Y = (6969 * R + 23434 * G + 2365 * B)/32768
         * integer approximation of
         * Y = 0.212671 * R + 0.715160 * G + 0.072169 * B
         */
        ptr_r = img;
        ptr_g = img + *nxp * *nyp;
        ptr_b = img + 2 * *nxp * *nyp;
        ptr_gray = img;
        ptr_end = ptr_gray + *nxp * *nyp;
        while (ptr_gray < ptr_end)
            *ptr_gray++ = (unsigned char) (6969 * *ptr_r++
                                           + 23434 * *ptr_g++
                                           + 2365 * *ptr_b++) / 32768;
        /* resize and return the image */
        img = realloc(img, *nxp * *nyp * sizeof(unsigned char));
        return img;
    }
}

/**
 * @brief read a PNG file into a 32bit float array
 *
 * The array contains the deinterlaced channels.
 * 1, 2, 4 and 8bit images are converted to float values
 * between 0. and 1., 3., 15. or 255.
 * 16bit images are also downscaled to 8bit before conversion.
 *
 * @param fname PNG file name
 * @param nxp, nyp, ncp pointers to variables to be filled with the number of
 *        columns, lines and channels of the image
 * @return pointer to an allocated unsigned char array of pixels,
 *         or NULL if an error happens
 */
float *read_png_f32(const char *fname, size_t * nxp, size_t * nyp,
                    size_t * ncp)
{
    /* read the image as float */
    return (float *) read_png_raw(fname, nxp, nyp, ncp,
                                  PNG_TRANSFORM_IDENTITY, IO_PNG_F32);
}

/**
 * @brief read a PNG file into a 32bit float array, converted to RGB
 *
 * See read_png_f32() for details.
 */
float *read_png_f32_rgb(const char *fname, size_t * nxp, size_t * nyp)
{
    size_t nc;
    float *img;

    /* read the image */
    img = (float *) read_png_raw(fname, nxp, nyp, &nc,
                                 PNG_TRANSFORM_STRIP_ALPHA, IO_PNG_F32);
    if (NULL == img)
        /* error */
        return NULL;
    if (3 == nc)
        /* already RGB */
        return img;
    else
    {
        /* convert to RGB */
        float *ptr_r, *ptr_g, *ptr_b, *ptr_end;

        /* resize the image */
        img = realloc(img, 3 * *nxp * *nyp * sizeof(float));

        /* gray->RGB conversion */
        ptr_r = img;
        ptr_end = ptr_r + *nxp * *nyp;
        ptr_g = img + *nxp * *nyp;
        ptr_b = img + 2 * *nxp * *nyp;
        while (ptr_r < ptr_end)
        {
            *ptr_g++ = *ptr_r;
            *ptr_b++ = *ptr_r++;
        }
        return img;
    }
}

/**
 * @brief read a PNG file into a 32bit float array, converted to gray
 *
 * See read_png_f32() for details.
 */
float *read_png_f32_gray(const char *fname, size_t * nxp, size_t * nyp)
{
    size_t nc;
    float *img;

    /* read the image */
    img = (float *) read_png_raw(fname, nxp, nyp, &nc,
                                 PNG_TRANSFORM_STRIP_ALPHA, IO_PNG_F32);
    if (NULL == img)
        /* error */
        return NULL;
    if (1 == nc)
        /* already gray */
        return img;
    else
    {
        /* convert to gray */
        float *ptr_r, *ptr_g, *ptr_b, *ptr_gray, *ptr_end;

        /*
         * RGB->gray conversion
         * Y = (6969 * R + 23434 * G + 2365 * B)/32768
         * integer approximation of
         * Y = 0.212671 * R + 0.715160 * G + 0.072169 * B
         */
        ptr_r = img;
        ptr_g = img + *nxp * *nyp;
        ptr_b = img + 2 * *nxp * *nyp;
        ptr_gray = img;
        ptr_end = ptr_gray + *nxp * *nyp;
        while (ptr_gray < ptr_end)
            *ptr_gray++ = (float) (6969 * *ptr_r++
                                   + 23434 * *ptr_g++
                                   + 2365 * *ptr_b++) / 32768;
        /* resize and return the image */
        img = realloc(img, *nxp * *nyp * sizeof(float));
        return img;
    }
}

/*
 * WRITE
 */

/**
 * @brief internal function used to cleanup the memory when
 * png_write_raw() fails
 *
 * @param fp file pointer to close, ignored if NULL
 * @param idata, row_pointers arrays to free, ignored if NULL
 * @param png_ptr_p, info_ptr_p, pointers to PNG structure pointers,
 *        ignored if NULL
 * @return -1
 */
static int write_png_abort(FILE * fp,
                           png_byte * idata, png_bytep * row_pointers,
                           png_structp * png_ptr_p, png_infop * info_ptr_p)
{
    png_destroy_write_struct(png_ptr_p, info_ptr_p);
    if (NULL != row_pointers)
        free(row_pointers);
    if (NULL != idata)
        free(idata);
    if (NULL != fp && stdout != fp)
        (void) fclose(fp);
    return -1;
}

/**
 * @brief internal function used to write a byte array as a PNG file
 *
 * The PNG file is written as a 8bit image file, interlaced,
 * truecolor. Depending on the number of channels, the color model is
 * gray, gray+alpha, rgb, rgb+alpha.
 *
 * @todo handle 16bit
 *
 * @param fname PNG file name, "-" means stdout
 * @param data deinterlaced (RRR..GGG..BBB..AAA) image byte array
 * @param nx, ny, nc number of columns, lines and channels
 * @param dtype identifier for the data type to be used for output
 * @return 0 if everything OK, -1 if an error occured
 */
static int write_png_raw(const char *fname, const void *data,
                         size_t nx, size_t ny, size_t nc, int dtype)
{
    png_structp png_ptr;
    png_infop info_ptr;
    png_byte *idata = NULL, *idata_ptr = NULL;
    png_bytep *row_pointers = NULL;
    png_byte bit_depth;
    /* volatile : because of setjmp/longjmp */
    FILE *volatile fp;
    const unsigned char *data_u8 = NULL;
    const unsigned char *data_u8_ptr = NULL;
    const float *data_f32 = NULL;
    const float *data_f32_ptr = NULL;
    float tmp;
    int color_type, interlace, compression, filter;
    size_t size;
    size_t i, j, k;

    /* parameters check */
    if (0 >= nx || 0 >= ny || 0 >= nc)
        return -1;
    if (NULL == fname || NULL == data)
        return -1;
    if (IO_PNG_U8 != dtype && IO_PNG_F32 != dtype)
        return -1;

    /* open the PNG output file */
    if (0 == strcmp(fname, "-"))
        fp = stdout;
    else if (NULL == (fp = fopen(fname, "wb")))
        return -1;

    /* allocate the interlaced array and row pointers */
    size = nx * ny * nc;
    if (NULL == (idata = (png_byte *) malloc(size * sizeof(png_byte))))
        return write_png_abort(fp, NULL, NULL, NULL, NULL);

    if (NULL == (row_pointers = (png_bytep *) malloc(ny * sizeof(png_bytep))))
        return write_png_abort(fp, idata, NULL, NULL, NULL);

    /*
     * create and initialize the png_struct
     * with the default stderr and error handling
     */
    if (NULL == (png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                   NULL, NULL, NULL)))
        return write_png_abort(fp, idata, row_pointers, NULL, NULL);

    /* allocate/initialize the memory for image information */
    if (NULL == (info_ptr = png_create_info_struct(png_ptr)))
        return write_png_abort(fp, idata, row_pointers, &png_ptr, NULL);

    /* set error handling */
    if (0 != setjmp(png_jmpbuf(png_ptr)))
        /* if we get here, we had a problem reading the file */
        return write_png_abort(fp, idata, row_pointers, &png_ptr, &info_ptr);

    /* set up the input control using standard C streams */
    png_init_io(png_ptr, fp);

    /* set image informations */
    bit_depth = 8;
    switch (nc)
    {
    case 1:
        color_type = PNG_COLOR_TYPE_GRAY;
        break;
    case 2:
        color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
        break;
    case 3:
        color_type = PNG_COLOR_TYPE_RGB;
        break;
    case 4:
        color_type = PNG_COLOR_TYPE_RGB_ALPHA;
        break;
    default:
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        free(row_pointers);
        free(idata);
        (void) fclose(fp);
        return -1;
    }
    interlace = PNG_INTERLACE_ADAM7;
    compression = PNG_COMPRESSION_TYPE_BASE;
    filter = PNG_FILTER_TYPE_BASE;

    /* set image header */
    png_set_IHDR(png_ptr, info_ptr, (png_uint_32) nx, (png_uint_32) ny,
                 bit_depth, color_type, interlace, compression, filter);
    /* TODO : significant bit (sBIT), gamma (gAMA), comments (text) chunks */
    png_write_info(png_ptr, info_ptr);

    /*
     * interlace and convert RRR GGG BBB to RGB RGB RGB
     * the image is interlaced layer after layer
     * this involves more memory exchange, but allows a generic loop
     */
    switch (dtype)
    {
    case IO_PNG_U8:
        data_u8 = (unsigned char *) data;
        for (k = 0; k < nc; k++)
        {
            /* channel loop */
            data_u8_ptr = data_u8 + (size_t) (nx * ny * k);
            idata_ptr = idata + (size_t) k;
            for (j = 0; j < ny; j++)
            {
                /* row loop */
                for (i = 0; i < nx; i++)
                {
                    /* pixel loop */
                    *idata_ptr = (png_byte) * data_u8_ptr++;
                    idata_ptr += nc;
                }
            }
        }
        break;
    case IO_PNG_F32:
        data_f32 = (float *) data;
        for (k = 0; k < nc; k++)
        {
            /* channel loop */
            data_f32_ptr = data_f32 + (size_t) (nx * ny * k);
            idata_ptr = idata + (size_t) k;
            for (j = 0; j < ny; j++)
            {
                /* row loop */
                for (i = 0; i < nx; i++)
                {
                    /* pixel loop */
                    tmp = floor(*data_f32_ptr++ + .5f);
                    *idata_ptr = (png_byte) (tmp < 0. ? 0. :
                                             (tmp > 255. ? 255. : tmp));
                    idata_ptr += nc;
                }
            }
        }
        break;
    }

    /* set row pointers */
    for (j = 0; j < ny; j++)
        row_pointers[j] = idata + (size_t) (nc * nx * j);

    /* write out the entire image and end it */
    png_write_image(png_ptr, row_pointers);
    png_write_end(png_ptr, info_ptr);

    /* clean up and free any memory allocated, close the file */
    (void) write_png_abort(fp, idata, row_pointers, &png_ptr, &info_ptr);

    return 0;
}

/**
 * @brief write a 8bit unsigned integer array into a PNG file
 *
 * @param fname PNG file name
 * @param data array to write
 * @param nx, ny, nc number of columns, lines and channels of the image
 * @return 0 if everything OK, -1 if an error occured
 */
int write_png_u8(const char *fname, const unsigned char *data,
                 size_t nx, size_t ny, size_t nc)
{
    return write_png_raw(fname, (void *) data,
                         (png_uint_32) nx, (png_uint_32) ny, (png_byte) nc,
                         IO_PNG_U8);
}

/**
 * @brief write a float array into a PNG file
 *
 * The float values are rounded to 8bit integers, and bounded to [0, 255].
 *
 * @todo handle 16bit images and flexible min/max
 *
 * @param fname PNG file name
 * @param data array to write
 * @param nx, ny, nc number of columns, lines and channels of the image
 * @return 0 if everything OK, -1 if an error occured
 */
int write_png_f32(const char *fname, const float *data,
                  size_t nx, size_t ny, size_t nc)
{
    return write_png_raw(fname, (void *) data,
                         (png_uint_32) nx, (png_uint_32) ny, (png_byte) nc,
                         IO_PNG_F32);
}
