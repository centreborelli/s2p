#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <tiffio.h>


// structs {{{1

struct tiff_tile {
	int32_t w, h;
	int16_t spp, bps, fmt;
	bool broken;
	uint8_t *data;
};

struct tiff_info {
	int32_t w;   // image width
	int32_t h;   // image height
	int16_t spp; // samples per pixel
	int16_t bps; // bits per sample
	int16_t fmt; // sample format
	bool broken; // whether pixels are contiguous or broken
	bool packed; // whether bps=1,2 or 4
	bool tiled;  // whether data is organized into tiles
	bool compressed;
	int ntiles;

	// only if tiled
	int32_t tw;  // tile width
	int32_t th;  // tile height
	int32_t ta;  // tiles across
	int32_t td;  // tiles down
};

// general utility functions {{{1


// exit the program printing an error message
#ifndef _FAIL_C
#define _FAIL_C
static void fail(const char *fmt, ...)
{
	va_list argp;
	fprintf(stderr, "\nERROR: ");
	va_start(argp, fmt);
	vfprintf(stderr, fmt, argp);
	va_end(argp);
	fprintf(stderr, "\n\n");
	fflush(NULL);
#  ifdef NDEBUG
	exit(-1);
#  else//NDEBUG
	exit(*(int *)0x43);
#  endif//NDEBUG
}
#endif//_FAIL_C

//static double global_accumulated_size = 0;
//static double global_tile_size = 0;


// enable "XMALLOC_STATS" for debugging purposes
//#define XMALLOC_STATS

#ifdef XMALLOC_STATS
#include "xmalloc_stats.c"
#else

#ifndef _XMALLOC_C
#define _XMALLOC_C
static void *xmalloc(size_t size)
{
	if (size == 0)
		fail("xmalloc: zero size");
	void *p = malloc(size);
	if (!p)
	{
		double sm = size / (0x100000 * 1.0);
		fail("xmalloc: out of memory when requesting "
			"%zu bytes (%gMB)",//:\"%s\"",
			size, sm);//, strerror(errno));
	}
	return p;
}
#endif//_XMALLOC_C
#define xfree free
#endif//XMALLOC_STATS



static int my_computetile(struct tiff_info *t, int i, int j)
{
	if (i < 0 || j < 0 || i >= t->w || j >= t->h)
		return -1;
		//fail("got bad pixel (%d %d) [%d %d]", i, j, t->w, t->h);
	int ti = i / t->tw;
	int tj = j / t->th;
	int r = tj * t->ta + ti;
	if (r < 0 || r >= t->ntiles)
		fail("bad tile index %d for point (%d %d)", r, i, j);
	return r;
}


// open a TIFF file, with some magic to access subimages
// (i.e., filename "file.tif,3" refers to the third sub-image)
static TIFF *tiffopen_fancy(char *filename, char *mode)
{
	//fprintf(stderr, "tiffopen fancy \"%s\",\"%s\"\n", filename, mode);
	char *comma = strrchr(filename, ',');
	if (*mode != 'r' || !comma)
	def:	return TIFFOpen(filename, mode);

	int aftercomma = strlen(comma + 1);
	int ndigits = strspn(comma + 1, "0123456789");

	if (aftercomma != ndigits) goto def;

	char buf[FILENAME_MAX];
	strncpy(buf, filename, FILENAME_MAX);
	comma = strrchr(buf, ',');
	*comma = '\0';
	int index = atoi(comma + 1);

	TIFF *tif = TIFFOpen(buf, mode);
	if (!tif) return tif;
	for (int i = 0; i < index; i++)
		TIFFReadDirectory(tif);

	return tif;
}


// tell how many units "u" are needed to cover a length "n"
static int how_many(int n, int u)
{
	assert(n > 0);
	assert(u > 0);
	return n/u + (bool)(n%u);
}



static void get_tiff_info(struct tiff_info *t, TIFF *tif)
{
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH,      &t->w);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH,     &t->h);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &t->spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,   &t->bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT,    &t->fmt);

	uint16_t planarity, compression;
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG,    &planarity);
	TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION,     &compression);
	t->broken = planarity != PLANARCONFIG_CONTIG;
	t->compressed = compression != COMPRESSION_NONE;
	t->packed = 0 != t->bps % 8;
	t->tiled = TIFFIsTiled(tif);

	if (t->tiled) {
		TIFFGetField(tif, TIFFTAG_TILEWIDTH,  &t->tw);
		TIFFGetField(tif, TIFFTAG_TILELENGTH, &t->th);
		t->ta = how_many(t->w, t->tw);
		t->td = how_many(t->h, t->th);
		t->ntiles = TIFFNumberOfTiles(tif);
		assert(t->ta * t->td == t->ntiles);
	} else {
		t->ta = t->td = 1;
		t->tw = t->w;
		t->th = t->h;
		t->ntiles = 1;
	}
}

static bool get_tiff_info_filename_e(struct tiff_info *t, char *fname)
{
	TIFF *tif = tiffopen_fancy(fname, "r");
	if (!tif)
		return false;
	get_tiff_info(t, tif);
	TIFFClose(tif);
	return true;
}

static int tiff_imagewidth(TIFF *tif)
{
	uint32_t w, r = TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	if (r != 1) fail("can not read TIFFTAG_IMAGEWIDTH");
	return w;
}

static int tiff_imagelength(TIFF *tif)
{
	uint32_t h, r = TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	if (r != 1) fail("can not read TIFFTAG_IMAGELENGTH");
	return h;
}

static int tiff_tilewidth(TIFF *tif)
{
	uint32_t R;
	int r = TIFFGetField(tif, TIFFTAG_TILEWIDTH, &R);
	if (r != 1) fail("TIFF tiles of unknown width");
	return R;
}

static int tiff_tilelength(TIFF *tif)
{
	uint32_t R;
	int r = TIFFGetField(tif, TIFFTAG_TILELENGTH, &R);
	if (r != 1) fail("TIFF tiles of unknown width");
	return R;
}

static int tiff_tilesacross(TIFF *tif)
{
	int w = tiff_imagewidth(tif);
	int tw = tiff_tilewidth(tif);
	return how_many(w, tw);
}

static int tiff_tilesdown(TIFF *tif)
{
	int h = tiff_imagelength(tif);
	int th = tiff_tilelength(tif);
	return how_many(h, th);
}

static int tiff_tile_corner(int p[2], TIFF *tif, int tidx)
{
	p[0] = p[1] = -1;

	int tw = tiff_tilewidth(tif);
	int th = tiff_tilelength(tif);
	int ta = tiff_tilesacross(tif);
	int td = tiff_tilesdown(tif);

	int i = tidx % ta;
	int j = tidx / ta;

	if (!(i < ta && j < td))
		return 0;

	p[0] = tw*i;
	p[1] = th*j;
	return 1;
}

static int tiff_samplesperpixel(TIFF *tif)
{
	uint16_t spp;
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	return spp;
}

// divide by 8 to obtain the size per sample (then, 0=packed data)
static int tiff_bitspersample(TIFF *tif)
{
	uint16_t bps;
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bps);
	return bps;
}

//static int tiff_sampleformat(TIFF *tif)
//{
//	uint16_t fmt;
//	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &fmt);
//	return fmt;
//}



static void read_scanlines(struct tiff_tile *tout, TIFF *tif)
{
	// read all file info
	struct tiff_info tinfo[1];
	get_tiff_info(tinfo, tif);

	// fill-in output information
	tout->w = tinfo->w;
	tout->h = tinfo->h;
	tout->bps = tinfo->bps;
	tout->fmt = tinfo->fmt;
	tout->spp = tinfo->spp;
	tout->broken = false;

	// define useful constants
	int pixel_size = tinfo->spp * tinfo->bps/8;
	int output_size = tout->w * tout->h * pixel_size;

	// allocate space for output data
	tout->data = xmalloc(output_size);

	// copy scanlines
	int scanline_size = TIFFScanlineSize(tif);
	assert(scanline_size == tinfo->w * pixel_size);
	for (int j = 0; j < tinfo->h; j++)
	{
		uint8_t *buf = tout->data + scanline_size*j;
		int r = TIFFReadScanline(tif, buf, j, 0);
		if (r < 0) fail("could not read scanline %d", j);
	}
}

static tsize_t my_readtile(TIFF *tif, tdata_t buf,
		uint32 x, uint32 y, uint32 z, tsample_t sample)
{
	tsize_t r = TIFFReadTile(tif, buf, x, y, z, sample);
	if (r == -1) memset(buf, 0, r = TIFFTileSize(tif));
	return r;
}

static void read_tile_from_file(struct tiff_tile *t, char *filename, int tidx)
{
	TIFF *tif = tiffopen_fancy(filename, "r");
	if (!tif) fail("could not open TIFF file \"%s\" for reading", filename);

	uint32_t w, h;
	uint16_t spp, bps, fmt, planarity;
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH,      &w);
	TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH,     &h);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE,   &bps);
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT,    &fmt);
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG,    &planarity);
	t->spp = spp;
	t->bps = bps;
	t->fmt = fmt;
	t->broken = false;

	if (planarity != PLANARCONFIG_CONTIG)
		fail("broken pixels not supported yet");

	if (TIFFIsTiled(tif)) {
		int nt = TIFFNumberOfTiles(tif);
		if (tidx < 0 || tidx >= nt)
			fail("bad tile index %d", tidx);

		t->w = tiff_tilewidth(tif);;
		t->h = tiff_tilelength(tif);

		int ii[2];
		tiff_tile_corner(ii, tif, tidx);
		//int ii = tw*i;
		//int jj = th*j;

		int tbytes = TIFFTileSize(tif);
		t->data = xmalloc(tbytes);
		memset(t->data, 0, tbytes);
		int r = my_readtile(tif, t->data, ii[0], ii[1], 0, 0);
		if (r != tbytes) fail("could not read tile");
	} else { // not tiled, read the whole image into 0th tile
		read_scanlines(t, tif);
	}

	TIFFClose(tif);
}

// overwrite a tile on an existing tiled TIFF image
static void put_tile_into_file(char *filename, struct tiff_tile *t, int tidx)
{
	if (t->broken) fail("can not save broken tiles yet");

	// Note, the mode "r+" is officially undocumented, but its behaviour is
	// described on file tif_open.c from libtiff.
	TIFF *tif = tiffopen_fancy(filename, "r+");
	if (!tif) fail("could not open TIFF file \"%s\" for writing", filename);

	int tw = tiff_tilewidth(tif);
	int th = tiff_tilewidth(tif);
	int spp = tiff_samplesperpixel(tif);
	int bps = tiff_bitspersample(tif);
	//int fmt = tiff_sampleformat(tif);
	if (tw != t->w) fail("tw=%d different to t->w=%d", tw, t->w);
	if (th != t->h) fail("th=%d different to t->h=%d", th, t->h);
	if (spp != t->spp) fail("spp=%d different to t->spp=%d", spp, t->spp);
	if (bps != t->bps) fail("bps=%d different to t->bps=%d", bps, t->bps);
	if (spp != t->spp) fail("spp=%d different to t->spp=%d", spp, t->spp);

	int ii[2];
	int r = tiff_tile_corner(ii, tif, tidx);
	if (!r) fail("bad tile %d", tidx);

	r = TIFFWriteTile(tif, t->data, ii[0], ii[1], 0, 0);

	TIFFClose(tif);
}

//static int fmt_from_string(char *f)
//{
//	if (0 == strcmp(f, "uint")) return SAMPLEFORMAT_UINT;
//	if (0 == strcmp(f, "int")) return SAMPLEFORMAT_INT;
//	if (0 == strcmp(f, "ieeefp")) return SAMPLEFORMAT_IEEEFP;
//	//if (0 == strcmp(f, "void")) return SAMPLEFORMAT_VOID;
//	if (0 == strcmp(f, "complexint")) return SAMPLEFORMAT_COMPLEXINT;
//	if (0 == strcmp(f, "complexieeefp")) return SAMPLEFORMAT_COMPLEXIEEEFP;
//	return SAMPLEFORMAT_VOID;
//}

static void create_zero_tiff_file(char *filename, int w, int h,
		int tw, int th, int spp, int bps, int fmt, bool incomplete,
		bool compressed)
{
	//if (bps != 8 && bps != 16 && bps == 32 && bps != 64
	//		&& bps != 128 && bps != 92)
	//	fail("bad bps=%d", bps);
	//if (spp < 1) fail("bad spp=%d", spp);
	//int fmt_id = fmt_from_string(fmt);
	double gigabytes = (spp/8.0) * w * h * bps / 1024.0 / 1024.0 / 1024.0;
	//fprintf(stderr, "gigabytes = %g\n", gigabytes);
	TIFF *tif = TIFFOpen(filename, gigabytes > 1 ? "w8" : "w");
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, spp);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bps);
	TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, fmt);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_TILEWIDTH, tw);
	TIFFSetField(tif, TIFFTAG_TILELENGTH, th);
	if (compressed)
		TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);


	int tilesize = tw * th * bps/8 * spp;
	uint8_t *buf = xmalloc(tilesize);
	for (int i = 0; i < tilesize; i++)
	{
		float x = ((i/((spp*bps)/8))%tw)/((double)tw)-0.5;
		float y = ((i/((spp*bps)/8))/th)/((double)th)-0.5;
		//buf[i] = 0;
		buf[i] = 127.5+128*cos(90*pow(hypot(x,y+0.3*x),1+(i%spp-1)/1.3));
	}
	TIFFWriteTile(tif, buf, 0, 0, 0, 0);
	TIFFClose(tif);
	if (!incomplete) {
		for (int j = 0; j < h; j += th)
		for (int i = 0; i < w; i += tw)
		{
			tif = TIFFOpen(filename, "r+");
			TIFFWriteTile(tif, buf, i, j, 0, 0);
			TIFFClose(tif);
		}
	}
	free(buf);
}



// getpixel cache with octaves {{{1

#define MAX_OCTAVES 25
struct tiff_octaves {
	// essential data
	//
	int noctaves;
	char filename[MAX_OCTAVES][FILENAME_MAX];
	struct tiff_info i[MAX_OCTAVES];
	void **c[MAX_OCTAVES];        // pointers to cached tiles

	bool option_read;
	bool option_write;
	bool *changed;

	// data only necessary to delete tiles when the memory is full
	//
	unsigned int *a[MAX_OCTAVES]; // access counter for each tile
	unsigned int ax; // global access counter
	int curtiles;    // current number of tiles in memory
	int maxtiles;    // number of tiles allowed to be in memory at once
};

//#include "smapa.h"
//SMART_PARAMETER(FIRST_OCTAVE,0)

//void my_tifferror(const char *module, const char *fmt, va_list ap)
//{
//	(void)ap;
//	if (0 == strcmp(fmt, "%llu: Invalid tile byte count, tile %lu"))
//		fprintf(stderr, "got a zero tile\n");
//	else
//		fprintf(stderr, "TIFF ERROR(%s): \"%s\"\n", module, fmt);
//}
static void disable_tiff_warnings_and_errors(void)
{
	TIFFSetWarningHandler(NULL);//suppress warnings
	TIFFSetErrorHandler(NULL);
}

void tiff_octaves_init0(struct tiff_octaves *t, char *filepattern,
		double megabytes, int max_octaves)
{
	//fprintf(stderr, "tiff octaves init \"%s\"(%gMB)\n", filepattern, megabytes);
	// create filenames until possible
	t->noctaves = 0;
	for (int o = 0; o < max_octaves; o++)
	{
		snprintf(t->filename[o], FILENAME_MAX, filepattern, o);
		if (!get_tiff_info_filename_e(t->i + o, t->filename[o]))
			break;
		if (t->i[o].bps < 8 || t->i[o].packed)
			fail("caching of packed samples is not supported");
		if (o > 0) { // check consistency
			if (0 == strcmp(t->filename[o], t->filename[0])) break;
			if (t->i[o].bps != t->i->bps) fail("inconsistent bps");
			if (t->i[o].spp != t->i->spp) fail("inconsistent spp");
			if (t->i[o].fmt != t->i->fmt) fail("inconsistent fmt");
			if (t->i[o].tw != t->i->tw) fail("inconsistent tw");
			if (t->i[o].th != t->i->th) fail("inconsistent th");
		}
		t->noctaves += 1;

	}
	if (t->noctaves < 1)
		fail("Could not get any file with pattern \"%s\"", filepattern);

	// set up essential data
	for (int o = 0; o < t->noctaves; o++)
	{
		t->c[o] = xmalloc((1 + t->i[o].ntiles) * sizeof*t->c);
		for (int j = 0; j < t->i[o].ntiles; j++)
			t->c[o][j] = 0;
	}

	// set up writing cache for setpixel (just in case)
	t->option_read = true;
	t->option_write = false;
	t->changed = xmalloc(t->i->ntiles * sizeof*t->changed);
	for (int i = 0; i < t->i->ntiles; i++)
		t->changed[i] = false;

	// set up data for old tile deletion
	if (megabytes) {
		for (int o = 0; o < t->noctaves; o++)
			t->a[o] = xmalloc(t->i[o].ntiles * sizeof*t->a[o]);
		t->ax = 0;
		int tilesize = t->i->tw * t->i->th * (t->i->bps/8) * t->i->spp;
		double mbts = tilesize / (1024.0 * 1024);
		t->maxtiles = megabytes / mbts;
		//fprintf(stderr, "maxtiles = %d\n", t->maxtiles);
		t->curtiles = 0;
	} else  {
		// unlimited tile usage
		t->a[0] = NULL;
	}
}

void tiff_octaves_init(struct tiff_octaves *t, char *filepattern,
		double megabytes)
{
	tiff_octaves_init0(t, filepattern, megabytes, MAX_OCTAVES);
}

static void re_write_tile(struct tiff_octaves *t, int tidx)
{
	if (!t->option_write)
		return;
	if (tidx < 0 || tidx >= t->i->ntiles)
		return;

	assert(t->c[0][tidx]);
	struct tiff_tile ti[1];
	//read_tile_from_file(ti, t->filename[0], tidx);
	ti->w = t->i->tw;
	ti->h = t->i->th;
	ti->spp = t->i->spp;
	ti->bps = t->i->bps;
	ti->fmt = t->i->fmt;
	ti->broken = false;
	ti->data = t->c[0][tidx];
	put_tile_into_file(t->filename[0], ti, tidx);
	t->changed[tidx] = false;
}


void tiff_octaves_free(struct tiff_octaves *t)
{
	if (t->option_write)
		for (int i = 0; i < t->i->ntiles; i++)
			if (t->changed[i])
				re_write_tile(t, i);
	for (int i = 0; i < t->noctaves; i++)
	{
		for (int j = 0; j < t->i[i].ntiles; j++)
			if (t->c[i][j])
				xfree(t->c[i][j]);
		xfree(t->c[i]);
		if (t->a[0])
			xfree(t->a[i]);
	}
	xfree(t->changed);
}


static void free_oldest_tile_octave(struct tiff_octaves *t)
{
	// find oldest tile
	int omin = -1, imin = -1;
	for (int o = 0; o < t->noctaves; o++)
		for (int i = 0; i < t->i[o].ntiles; i++)
			if (t->c[o][i]) {
				if (imin >= 0) {
					if (t->a[o][i] < t->a[omin][imin]) {
						imin = i;
						omin = o;
					}
				} else {
					imin = i;
					omin = o;
				}
			}
	assert(imin >= 0);
	assert(t->a[omin][imin] > 0);

	// free it
	//
	//fprintf(stderr, "CACHE: FREEing tile %d of octave %d (%g)\n", imin, omin,global_accumulated_size);
	if (t->option_write && omin == 0 && t->changed[imin])
		re_write_tile(t, imin);
	assert(t->c[omin][imin]);
	xfree(t->c[omin][imin]);
	t->c[omin][imin] = 0;
	t->a[omin][imin] = 0;
	t->curtiles -= 1;
}

//static void free_oldest_half_of_tiles(struct tiff_octaves *t)
//{
//
//}

static void notify_tile_access_octave(struct tiff_octaves *t, int o, int i)
{
	//fprintf(stderr, "notify tile %d\n", i);
	t->a[o][i] = ++t->ax;
}

void *tiff_octaves_gettile(struct tiff_octaves *t, int o, int i, int j)
{
	// sanitize input
	if (o < 0 || o >= t->noctaves) return NULL;
	if (i < 0 || i >= t->i[o].w) return NULL;
	if (j < 0 || j >= t->i[o].h) return NULL;

	// get valid tile index
	int tidx = my_computetile(t->i + o, i, j);
	if (tidx < 0) return NULL;

	// if tile does not exist, read it from file
	if (!t->c[o][tidx])
//#pragma omp critical
	{
		if (t->a[0] && t->curtiles == t->maxtiles)
			free_oldest_tile_octave(t);

		//fprintf(stderr,"CACHE: LOADing tile %d of octave %d (%g)\n",tidx,o, global_accumulated_size);
		struct tiff_tile tmp[1];
		read_tile_from_file(tmp, t->filename[o], tidx);
		t->c[o][tidx] = tmp->data;

		t->curtiles += 1;
	}
	if (t->a[0])
		notify_tile_access_octave(t, o, tidx);

	return t->c[o][tidx];
}

void *tiff_octaves_getpixel(struct tiff_octaves *t, int o, int i, int j)
{
	//fprintf(stderr, "t_o_g(%d, %d, %d)\n", o, i, j);
	void *tile = tiff_octaves_gettile(t, o, i, j);
	if (!tile) return NULL;

	// get pointer to requested pixel
	struct tiff_info *ti = t->i + o;
	int ii = i % ti->tw;
	int jj = j % ti->th;
	int pixel_index = jj * ti->tw + ii;
	int pixel_position = pixel_index * ti->spp * (ti->bps / 8);
	return pixel_position + (char*)tile;
}

void tiff_octaves_setpixel(struct tiff_octaves *t, int i, int j, void *p)
{
	if (i < 0 || i >= t->i->w) return;
	if (j < 0 || j >= t->i->h) return;

	int tidx = my_computetile(t->i, i, j);
	if (tidx < 0) return;
	assert(tidx < t->i->ntiles);
	void *tile = tiff_octaves_gettile(t, 0, i, j);
	if (!tile) return;

	int ii = i % t->i->tw;
	int jj = j % t->i->th;
	int pixel_index = jj * t->i->tw + ii;
	int pixel_size = (t->i->spp * t->i->bps) / 8;
	int pixel_position = pixel_index * pixel_size;
	void *where = pixel_position + (char*)tile;
	memcpy(where, p, pixel_size);
	//{
	//	char *from = (char*)p;
	//	char *to = pixel_position + (char*)tile;
	//	for (int k = 0; k < pixel_size; k++)
	//		to[k] = from[k];
	//}
	t->changed[tidx] = true;
}

static
void convert_floats_to_samples(struct tiff_info *ti, void *s, float *f, int n)
{
#define FIN for(int i=0;i<n;i++)
	switch(ti->fmt) {
	case SAMPLEFORMAT_UINT:
		if (8 == ti->bps)        FIN ((uint8_t *)s)[i] = f[i];
		else if (16 == ti->bps)  FIN ((uint16_t*)s)[i] = f[i];
		else if (32 == ti->bps)  FIN ((uint32_t*)s)[i] = f[i];
		break;
	case SAMPLEFORMAT_INT:
		if (8 == ti->bps)        FIN (( int8_t *)s)[i] = f[i];
		else if (16 == ti->bps)  FIN (( int16_t*)s)[i] = f[i];
		else if (32 == ti->bps)  FIN (( int32_t*)s)[i] = f[i];
		break;
	case SAMPLEFORMAT_IEEEFP:
		if (32 == ti->bps)       FIN (( float  *)s)[i] = f[i];
		else if (64 == ti->bps)  FIN (( double *)s)[i] = f[i];
		break;
	default:
		fail("unrecognized sample fmt(%d) size(%d)", ti->fmt, ti->bps);
	}
#undef FIN
}

void tiff_octaves_setpixel_float(struct tiff_octaves *t, int i, int j, float *p)
{
	int pixel_size = (t->i->spp * t->i->bps) / 8;
	char pix[pixel_size];
	convert_floats_to_samples(t->i, pix, p, t->i->spp);
	tiff_octaves_setpixel(t, i, j, pix);
}

//static int fmt_from_string(char *f)
//{
//	if (0 == strcmp(f, "uint")) return SAMPLEFORMAT_UINT;
//	if (0 == strcmp(f, "int")) return SAMPLEFORMAT_INT;
//	if (0 == strcmp(f, "ieeefp")) return SAMPLEFORMAT_IEEEFP;
//	//if (0 == strcmp(f, "void")) return SAMPLEFORMAT_VOID;
//	if (0 == strcmp(f, "complexint")) return SAMPLEFORMAT_COMPLEXINT;
//	if (0 == strcmp(f, "complexieeefp")) return SAMPLEFORMAT_COMPLEXIEEEFP;
//	return SAMPLEFORMAT_VOID;
//}
