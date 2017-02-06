#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fancy_image.h"

#include "iio.h"
#include "tiff_octaves_rw.c"

// the following struct is an implementation detail,
// it is only used on this file
//
struct FI {
	// part of the public interface, "struct fancy_image"
	int w;
	int h;
	int pd;
	int no;

	// options
	double megabytes;
	int max_octaves;
	bool option_read;
	bool option_write;
	bool option_creat;
	int option_verbose;
	int option_w, option_h, option_tw, option_th, option_spp, option_bps;
	int option_fmt; // SAMPLEFORMAT_UINT, etc.
	int option_compressed;

	// implementation details
	bool tiffo;
	struct tiff_octaves t[1];
	float *x, *pyr_x[MAX_OCTAVES];
	int pyr_w[MAX_OCTAVES], pyr_h[MAX_OCTAVES];
	bool x_changed;
	char x_filename[FILENAME_MAX];

};

// Compiler trick to check that "FI" can fit inside a "fancy_image"
// (if this line fails, increase the padding at the struct fancy_image on fi.h)
typedef char check_FI_size[sizeof(struct FI)<=sizeof(struct fancy_image)?1:-1];


// check whether a filename corresponds to a small image or a tiled tiff
static bool filename_corresponds_to_tiffo(char *filename)
{
	if (0 == strcmp(filename, "-"))
		return false;
	struct tiff_info ti[1];
	disable_tiff_warnings_and_errors();
	bool r = get_tiff_info_filename_e(ti, filename);
	if (!r) {
		char buf[FILENAME_MAX];
		snprintf(buf, FILENAME_MAX, filename, 0);
		r = get_tiff_info_filename_e(ti, buf);
		if (!r)
			return false;
		return ti->tiled;
	}
	return ti->tiled;
}

// type of a "zoom-out" function
typedef void (*zoom_out_function_t)(float*,int,int,float*,int,int,int);

static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w) return 0;
	if (j < 0 || j >= h) return 0;
	return x[pd*(j*w+i)+l];
}

static void zoom_out_by_factor_two(float *out, int ow, int oh,
		float *in, int iw, int ih, int pd)
{
	assert(abs(2*ow-iw) < 2);
	assert(abs(2*oh-ih) < 2);
	for (int j = 0; j < oh; j++)
	for (int i = 0; i < ow; i++)
	for (int l = 0; l < pd; l++)
	{
		float a[4];
		a[0] = getsample_0(in, iw, ih, pd, 2*i  , 2*j  , l);
		a[1] = getsample_0(in, iw, ih, pd, 2*i+1, 2*j  , l);
		a[2] = getsample_0(in, iw, ih, pd, 2*i  , 2*j+1, l);
		a[3] = getsample_0(in, iw, ih, pd, 2*i+1, 2*j+1, l);
		out[pd*(ow*j + i)+l] = (a[0] + a[1] + a[2] + a[3])/4;
	}
}

static int build_pyramid(struct FI *f, int max_octaves)
{
	zoom_out_function_t z = zoom_out_by_factor_two;
	f->pyr_w[0] = f->w;
	f->pyr_h[0] = f->h;
	f->pyr_x[0] = f->x;
	int s = 1;
	while(1) {
		if (s > max_octaves) break;
		int      lw   = f->pyr_w[s-1];
		int      lh   = f->pyr_h[s-1];
		float   *lx   = f->pyr_x[s-1];
		int      sw   = ceil(lw / 2.0);
		int      sh   = ceil(lh / 2.0);
		float   *sx = xmalloc(f->pd * sw * sh * sizeof*sx);
		z(sx, sw, sh, lx, lw, lh, f->pd);
		f->pyr_w[s]   = sw;
		f->pyr_h[s]   = sh;
		f->pyr_x[s]   = sx;
		s += 1;
		if (sw + sh <= 2) break;
	}
	if (f->option_verbose)
	{
		fprintf(stderr, "IMG %dx%d,%d (no=%d)\n", f->w, f->h, f->pd, s);
		for (int i = 0; i < s; i++)
			fprintf(stderr, "\tpyr[%d] : %dx%d\n", i,
					f->pyr_w[i], f->pyr_h[i]);
	}
	return s;
}

static void free_pyramid(struct FI *f)
{
	for (int s = 0; s < f->no; s++)
		free(f->pyr_x[s]);
}

static void set_default_options(struct FI *f)
{
	f->megabytes = 100;
	f->option_read = true;
	f->option_write = false;
	f->option_creat = false;
	f->option_verbose = 0;
	f->max_octaves = MAX_OCTAVES;
	f->option_fmt = SAMPLEFORMAT_UINT;
	f->option_w = 256;
	f->option_h = 256;
	f->option_th = 0;
	f->option_tw = 0;
	f->option_compressed = 0;
	f->option_spp = 1;
	f->option_bps = 8;
}

static void interpret_options(struct FI *f, char *options_arg)
{
	set_default_options(f);

	int n = 1 + strlen(options_arg);
	char o[n];
	strncpy(o, options_arg, n);
	char *tok = strtok(o, ",");
	while (tok) {
		//fprintf(stderr, "TOKEN \"%s\"\n", tok);
		if (0 == strcmp(tok, "r"   ))  f->option_read  = true;
		if (0 == strcmp(tok, "rw"  ))  f->option_read  = true;
		if (0 == strcmp(tok, "wr"  ))  f->option_read  = true;
		if (0 == strcmp(tok, "read"))  f->option_read  = true;
		if (0 == strcmp(tok, "w"   ))  f->option_write = true;
		if (0 == strcmp(tok, "rw"  ))  f->option_write = true;
		if (0 == strcmp(tok, "wr"  ))  f->option_write = true;
		if (0 == strcmp(tok, "write")) f->option_write = true;
		if (0 == strcmp(tok, "c"   ))  f->option_creat = true;
		if (0 == strcmp(tok, "creat"   ))  f->option_creat = true;
		double x;
		if (1 == sscanf(tok, "megabytes=%lf",&x)) f->megabytes      = x;
		if (1 == sscanf(tok, "octaves=%lf", &x))  f->max_octaves    = x;
		if (1 == sscanf(tok, "verbose=%lf", &x))  f->option_verbose = x;
		if (1 == sscanf(tok, "width=%lf", &x))    f->option_w       = x;
		if (1 == sscanf(tok, "height=%lf", &x))   f->option_h       = x;
		if (1 == sscanf(tok, "w=%lf", &x))        f->option_w       = x;
		if (1 == sscanf(tok, "h=%lf", &x))        f->option_h       = x;
		if (1 == sscanf(tok, "tw=%lf", &x))       f->option_tw      = x;
		if (1 == sscanf(tok, "th=%lf", &x))       f->option_th      = x;
		if (1 == sscanf(tok, "spp=%lf", &x))      f->option_spp     = x;
		if (1 == sscanf(tok, "pd=%lf", &x))       f->option_spp     = x;
		if (1 == sscanf(tok, "bps=%lf", &x))      f->option_bps     = x;
		if (1 == sscanf(tok, "fmt=%lf", &x))      f->option_fmt     = x;
		if (1 == sscanf(tok, "tilewidth=%lf", &x))f->option_tw      = x;
		if (1 == sscanf(tok, "tileheight=%lf",&x))f->option_tw      = x;
		tok = strtok(NULL, ",");
	}

	// pyramidal writing is complicated, so we disable it
	if (f->option_write)
		f->max_octaves = 1;

	// when c, complete the missing options in a consistent way
	if (f->option_creat)
		f->option_write = 1;
	

}

void create_iio_file(char *filename, int w, int h, int spp)
{
	int n = w * h * spp;
	uint8_t *buf = xmalloc(n);
	memset(buf, 0, n);
	iio_save_image_uint8_vec(filename, buf, w, h, spp);
	free(buf);
}


// API: open a fancy image
struct fancy_image *fancy_image_open(char *filename, char *options)
{
	// create return struct and its alias
	//struct fancy_image r[1];
	struct fancy_image *r = xmalloc(sizeof*r); // I hate this malloc!
	struct FI *f = (void*)r;

	// process options parameter
	interpret_options(f, options);

	// if "c", do create the file
	if (f->option_creat) {
		if (filename_corresponds_to_tiffo(filename) || f->option_tw > 0)
			create_zero_tiff_file(filename,
					f->option_w, f->option_h,
					f->option_tw, f->option_th,
					f->option_spp, f->option_bps,
					f->option_fmt,
					true, f->option_compressed);
		else
			create_iio_file(filename, f->option_w, f->option_h,
					f->option_spp);
	}

	// read the image
	if (filename_corresponds_to_tiffo(filename)) {
		f->tiffo = true;
		tiff_octaves_init0(f->t, filename, f->megabytes,f->max_octaves);
		if (f->option_write) f->t->option_write = true;
		f->w = f->t->i->w;
		f->h = f->t->i->h;
		f->pd = f->t->i->spp;
		f->no = f->t->noctaves;
	} else {
		f->tiffo = false;
		f->x = iio_read_image_float_vec(filename, &f->w, &f->h, &f->pd);
		f->no = build_pyramid(f, f->max_octaves);
		strncpy(f->x_filename, filename, FILENAME_MAX);
		f->x_changed = false;
	}

	if (f->option_verbose) {
		fprintf(stderr, "FANCY IMAGE \"%s\"\n", filename);
		fprintf(stderr, "\tw = %d\n", f->w);
		fprintf(stderr, "\th = %d\n", f->h);
		fprintf(stderr, "\tpd = %d\n", f->pd);
		fprintf(stderr, "\tno = %d\n", f->no);
		fprintf(stderr, "\n");
		fprintf(stderr, "\tmax_octaves= %d\n", f->max_octaves);
		fprintf(stderr, "\ttiffo = %d\n", f->tiffo);
		fprintf(stderr, "\tmegabytes = %g\n", f->megabytes);
	}

	// return image struct
	return r;
}

// call "fancy_image_open" with named options
struct fancy_image *fancy_image_create(char *filename, char *fmt, ...)
{
	int blen = 2000;
	char buf[blen];
	buf[0] = 'c';
	buf[1] = ',';
	va_list argp;
	va_start(argp, fmt);
	vsnprintf(buf + 2, blen - 2, fmt, argp);
	va_end(argp);
	return fancy_image_open(filename, buf);
}

int fancy_image_width_octave(struct fancy_image *fi, int octave)
{
	struct FI *f = (void*)fi;
	if (octave < 0 || octave >= f->no) return 0;
	return f->tiffo ? f->t->i[octave].w : f->pyr_w[octave];
}

int fancy_image_height_octave(struct fancy_image *fi, int octave)
{
	struct FI *f = (void*)fi;
	if (octave < 0 || octave >= f->no) return 0;
	return f->tiffo ? f->t->i[octave].h : f->pyr_h[octave];
}

// API: close a fancy image
void fancy_image_close(struct fancy_image *fi)
{
	struct FI *f = (void*)fi;

	if (f->tiffo)
		tiff_octaves_free(f->t);
	else {
		if ((f->option_write && f->x_changed) || f->option_creat)
			iio_save_image_float_vec(f->x_filename, f->x,
					f->w, f->h, f->pd);
		if (f->no > 1)
			free_pyramid(f);
		else
			free(f->x);
	}
	free(f);
}

// internal conversion function function
static float convert_sample_to_float(struct tiff_info *ti, void *p)
{
	switch(ti->fmt) {
	case SAMPLEFORMAT_UINT:
		if (8 == ti->bps)        return *(uint8_t *)p;
		else if (16 == ti->bps)  return *(uint16_t *)p;
		else if (32 == ti->bps)  return *(uint16_t *)p;
		break;
	case SAMPLEFORMAT_INT:
		if (8 == ti->bps)        return *(int8_t *)p;
		else if (16 == ti->bps)  return *(int16_t *)p;
		else if (32 == ti->bps)  return *(int16_t *)p;
		break;
	case SAMPLEFORMAT_IEEEFP:
		if (32 == ti->bps)       return *(float *)p;
		else if (64 == ti->bps)  return *(double *)p;
		break;
	}
	return NAN;
}

// API: get a sample from an image
float fancy_image_getsample(struct fancy_image *fi, int i, int j, int l)
{
	return fancy_image_getsample_oct(fi, 0, i, j, l);
}

// API: get a sample from an image, at the requested octave
float fancy_image_getsample_oct(struct fancy_image *fi,
		int octave, int i,int j, int l)
{
	struct FI *f = (void*)fi;

	if (octave < 0 || octave >= f->no)
		return NAN;
	if (l < 0) l = 0;
	if (l >= f->pd) l = f->pd - 1;

	if (f->tiffo) {
		uint8_t *p_pixel = tiff_octaves_getpixel(f->t, octave, i, j);
		if (!p_pixel) return NAN;
		uint8_t *p_sample = p_pixel + (l * f->t->i->bps) / 8;
		return convert_sample_to_float(f->t->i, p_sample);
	} else {
		float *x = f->pyr_x[octave];
		int    w = f->pyr_w[octave];
		int    h = f->pyr_h[octave];
		if (i < 0 || j < 0 || i >= w || j >= h)
			return NAN;
		int  idx = (j * w + i) * f->pd + l;
		return x[idx];
	}
}

// API: set a sample of an image
int fancy_image_setsample(struct fancy_image *fi, int i, int j, int l, float v)
{
	struct FI *f = (void*)fi;
	if (!f->option_write) return false;

	if (i < 0 || i >= f->w) return false;
	if (j < 0 || j >= f->h) return false;
	if (l < 0 || l >= f->pd) return false;

	if (f->tiffo) {
		float p[f->pd];
		// TODO: remove this silly loop
		for (int k = 0; k < f->pd; k++)
			p[k] = fancy_image_getsample(fi, i, j, k);
		p[l] = v;
		tiff_octaves_setpixel_float(f->t, i, j, p);
		return true;
	} else {
		int idx = (j * f->w + i) * f->pd + l;
		f->x[idx] = v;
		return f->x_changed = true;
	}
}

int fancy_image_leak_tiff_info(int *tw, int *th, int *fmt, int *bps,
		struct fancy_image *fi)
{
	struct FI *f = (void*)fi;
	if (f->tiffo) {
		if (tw)  *tw = f->t->i->tw;
		if (th)  *th = f->t->i->th;
		if (fmt) *fmt = f->t->i->fmt;
		if (bps) *bps = f->t->i->bps;
		return 1;
	}
	return 0;
}

void fancy_image_fill_rectangle_float_vec(
		float *out, int w, int h,
		struct fancy_image *f, int o, int x0, int y0)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	for (int l = 0; l < f->pd; l++)
	{
		int idx = (j * w + i) * f->pd + l;
		out[idx] = fancy_image_getsample_oct(f, o, x0 + i, y0 + j, l);
	}
}

void fancy_image_fill_rectangle_float_split(
		float *out, int w, int h,
		struct fancy_image *f, int o, int x0, int y0)
{
	for (int l = 0; l < f->pd; l++)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		int idx = l*w*h + j*w + i;
		out[idx] = fancy_image_getsample_oct(f, o, x0 + i, y0 + j, l);
	}
}

#ifdef MAIN_FI
#include <stdio.h>
int main_example(int c, char *v[])
{
	// process input arguments
	if (c != 7) {
		fprintf(stderr, "usage:\n\t%s image opts o i j l\n", *v);
		//                          0 1     2    3 4 5 6
		return 1;
	}
	char *filename = v[1];
	char *opts  = v[2];
	int arg_o = atoi(v[3]);
	int arg_i = atoi(v[4]);
	int arg_j = atoi(v[5]);
	int arg_l = atoi(v[6]);

	// do stuff
	struct fancy_image *f = fancy_image_open(filename, opts);
	printf("image \"%s\"\n", filename);
	printf("\tw  = %d\n", f->w);
	printf("\th  = %d\n", f->h);
	printf("\tpd = %d\n", f->pd);
	printf("\tno = %d\n", f->no);

	float s = fancy_image_getsample_oct(f, arg_o, arg_i, arg_j, arg_l);
	printf("\t (%d)[%d,%d]{%d} = %g\n", arg_o, arg_i, arg_j, arg_l, s);

	//int x = 1;
	//float s;
	//do {
	//	s = fancy_image_getsample(&f, x, x, 0);
	//	printf("\tsample(%d,%d,0) = %g\n", x, x, s);
	//	x *= 10;
	//} while(isfinite(s));
	fancy_image_close(f);

	// exit
	return 0;
}
int main_croparound(int c, char *v[])
{
	// process input arguments
	if (c != 8) {
		fprintf(stderr, "usage:\n\t"
				"%s in.tiff opts o cx cy ww out.png\n",*v);
		//                0 1       2    3 4  5  6  7
		return 1;
	}
	char *filename_in = v[1];
	char *opts   = v[2];
	int octave = atoi(v[3]);
	int cent_x = atoi(v[4]);
	int cent_y = atoi(v[5]);
	int diamet = atoi(v[6]);
	char *filename_out = v[7];

	struct fancy_image *f = fancy_image_open(filename_in, opts);
	if (!f->pd) return 2;
	//if (octave < 0) octave = 0;
	//if (octave >= f.no) octave = f.no - 1;
	float *x = xmalloc(diamet * diamet * f->pd * sizeof*x);
	for (int j = 0; j < diamet; j++)
	for (int i = 0; i < diamet; i++)
	for (int l = 0; l < f->pd; l++)
	{
		int ii = cent_x - diamet/2 + i;
		int jj = cent_y - diamet/2 + j;
		int idx_o = (j * diamet + i) * f->pd + l;
		x[idx_o] = fancy_image_getsample_oct(f, octave, ii, jj, l);
	}
	fancy_image_close(f);
	iio_save_image_float_vec(filename_out, x, diamet, diamet, f->pd);
	free(x);
	return 0;
}
int main_setsample(int c, char *v[])
{
	if (c != 7) {
		fprintf(stderr, "usage:\n\t"
				"%s inout.tiff opts i j l v\n",*v);
		//                0 1          2    3 4 5 6
		return 1;
	}
	char *filename = v[1];
	char *opts = v[2];
	int arg_i = atoi(v[3]);
	int arg_j = atoi(v[4]);
	int arg_l = atoi(v[5]);
	float arg_v = atof(v[6]);

	struct fancy_image *f = fancy_image_open(filename, opts);
	fancy_image_setsample(f, arg_i, arg_j, arg_l, arg_v);
	fancy_image_close(f);

	return 0;
}

int main_times(int c, char *v[])
{
	// process input arguments
	if (c != 3) {
		fprintf(stderr, "usage:\n\t%s file.tiff factor\n", *v);
		//                          0 1         2
		return 1;
	}
	char *filename = v[1];
	double factor = atof(v[2]);

	// open image
	struct fancy_image *f = fancy_image_open(filename, "rw,megabytes=33");

	// process data
	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	for (int l = 0; l < f->pd; l++)
	{
		double x = fancy_image_getsample(f, i, j, l);
		x = x * factor;
		fancy_image_setsample(f, i, j, l, x);
	}

	// close image (and save remaining updated tiles)
	fancy_image_close(f);

	// exit
	return 0;
}





//int main(int c, char *v[]) { return main_croparound(c, v); }
//int main(int c, char *v[]) { return main_setsample(c, v); }
int main(int c, char *v[]) { return main_times(c, v); }
#endif//MAIN_FI
