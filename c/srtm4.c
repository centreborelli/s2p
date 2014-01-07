// SRTM4 files are of size 6000x6000 pixels and cover an area of 5x5 degrees

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>



//#define SRTM4_URL "ftp://xftp.jrc.it/pub/srtmV4/arcasci/srtm_%02d_%02d.zip"
#define SRTM4_URL "--http-user=data_public --http-password=GDdci http://data.cgiar-csi.org/srtm/tiles/ASCII/srtm_%02d_%02d.zip"
#define SRTM4_ASC "%s/srtm_%02d_%02d.asc"



// download the contents of an url into a file
static int download(const char *to_file, const char *from_url)
{
	int nbuf = 2*FILENAME_MAX;
	char buf[nbuf];
	snprintf(buf, nbuf, "wget %s -O %s", from_url, to_file);
	int r = system(buf);
	return r;
}

// return the name of the cache directory
static char *cachedir(void)
{
	// XXX TODO FIXME WRONG WARNING : THIS code is not reentrant
	static char output_dirname[FILENAME_MAX];
	char *env_cache = getenv("SRTM4_CACHE");
	if (env_cache) {
		snprintf(output_dirname, FILENAME_MAX, "%s", env_cache);
	} else {
		char *homedir = getenv("HOME");
		if (!homedir)
			homedir = "/tmp";
		snprintf(output_dirname, FILENAME_MAX, "%s/.srtm4", homedir);
	}
	int r = mkdir(output_dirname, 0777);
	if (r != 0 && errno != EEXIST) {
		fprintf(stderr, "SRTM4: cannot create cachedir \"%s\"\n",
				output_dirname);
		exit(2);
	}
	return output_dirname;
}


static void get_tile_index_and_position(
		int *tlon, int *tlat,
		float *xlon, float *xlat,
		double lon, double lat)
{
	if (lat > 60) lat = 60;
	if (lat < -60) lat = -60;

    // tiles longitude indexes go from 1 to 72, covering the range from -180 to
    // +180
    *tlon = fmod(1+floor((lon+180)/5), 72);
    if (*tlon == 0) *tlon = 72;
	lon = fmod(lon + 180, 360);
	*xlon = 1200*fmod(lon, 5);

    // tiles latitude indexes go from 1 to 24, covering the range from 60 to
    // -60
    *tlat = 1+floor((60-lat)/5);
    if (*tlat == 25) *tlat = 24;
	*xlat = 1200*fmod(60-lat, 5);

//	fprintf(stderr, "tlon = %d\n", *tlon);
//	fprintf(stderr, "tlat = %d\n", *tlat);
//	fprintf(stderr, "xlon = %g\n", *xlon);
//	fprintf(stderr, "xlat = %g\n", *xlat);

	assert(*tlon > 0); assert(*tlon <= 72);
	assert(*tlat > 0); assert(*tlat <= 24);
	assert(*xlon >= 0); assert(*xlon < 6000);
	assert(*xlat >= 0); assert(*xlat < 6000);
}

static bool file_exists(const char *fname)
{
	FILE *f = fopen(fname, "r");
	if (f)
	{
		fclose(f);
		return true;
	}
	return false;
}

static void download_tile_file(int tlon, int tlat)
{
	int n = FILENAME_MAX;
	char url[n], zipname[n];
	snprintf(url, n, SRTM4_URL, tlon, tlat);
	snprintf(zipname, n, "%s/tmp.zip", cachedir());
	int rd = download(zipname, url);
	if (0 == rd) {
		char cmd[n];
		snprintf(cmd, n, "unzip -qq -o %s -d %s", zipname, cachedir());
		rd = system(cmd);
		if (rd) {
			fprintf(stderr,"failed unzipping file %s\n", zipname);
		}
		// TODO: do something if unzip fails
	}
}

static char *get_asc_tile_filename(int tlon, int tlat)
{
	// XXX TODO FIXME WRONG WARNING : THIS code is not reentrant
	static char fname[FILENAME_MAX];
	snprintf(fname, FILENAME_MAX, SRTM4_ASC, cachedir(), tlon, tlat);
	return fname;
}

static char *malloc_file_contents(char *filename, int *size)
{
	FILE *f = fopen(filename, "r");
	if (!f)
		return NULL;
	fseek(f, 0, SEEK_END);
	int n = ftell(f);
	fseek(f, 0, SEEK_SET);
	char *r = malloc(n+1);
	if (!r)
		return NULL;
	int rr = fread(r, 1, n, f);
	if (rr != n)
		return NULL;
	fclose(f);
	r[n] = '\0';
	*size = n+1;
	return r;
}

static int my_atoi(char *s)
{
	int r = 0, sign = 1;
	if (*s == '-') {
		sign = -1;
		s += 1;
	}
	while (*s)
	{
		r *= 10;
		r += *s - '0';
		s += 1;
	}
	return sign * r;
}

static const int spacing[0x100] = {
	[' '] = 1, ['\t'] = 1, ['\f'] = 1, ['\n'] = 1, ['\r'] = 1
};

static int my_isspace(unsigned char c)
{
	return spacing[c];
}

// this function is like regular strtok, but faster
// the delimiter list is, implicitly, spaces
static char *my_strtok(char *str)
{
	static char *next; // points to the remaining part of the string
	char *begin = str ? str : next, *s = begin;
	if (!*s)
		return NULL;
	while (*s && !my_isspace(*s))
		s++;
	if (!*s)
		return begin;
	assert(my_isspace(*s));
	while(*s && my_isspace(*s))
		s++;
	if (!*s)
		return begin;
	assert(!my_isspace(*s));
	s--;
	*s = '\0';
	next = s + 1;
	return begin;
}

// parse a tile file into memory
// (this function is ugly due to the error checking)
static float *malloc_tile_data(char *tile_filename)
{
	int fsize;
	char *f = malloc_file_contents(tile_filename, &fsize);
	if (!f)
		return NULL;
	int finc, ncols, nrows, xllcorner, yllcorner;
	double cellsize, nodatavalue;
	cellsize=nodatavalue=finc=ncols=nrows=xllcorner=yllcorner=-42;

	int r = sscanf(f, "ncols %d\nnrows %d\nxllcorner %d\nyllcorner %d\n"
			"cellsize %lf\nNODATA_value %lf%n", &ncols, &nrows,
			&xllcorner, &yllcorner, &cellsize, &nodatavalue, &finc);

	if (r != 6) {
		fprintf(stderr, "strange tile file format on \"%s\" (%d)\n",
				tile_filename, r);
		exit(2);
	}
	if (ncols != 6000 || nrows != 6000) {
		fprintf(stderr,"ncols,nrows = %d,%d != 6000,6000",ncols,nrows);
		exit(2);
	}
	int n = 6000*6000;
	float *t = malloc(n*sizeof*t);
	if (!t) {
		fprintf(stderr, "out of memory!\n");
		exit(2);
	}
	int cx = 0;
	char *fp = f + finc;
	char *tok = my_strtok(fp);
	while (tok && cx < n) {
		int x = my_atoi(tok);
		t[cx++] = x > -100 ? x : NAN;
		tok = my_strtok(NULL);
	}
	free(f);
	if (cx != n) {
		fprintf(stderr, "could only read %d numbers from file \n", cx);
		exit(2);
	}
	return t;
}

static float *global_table_of_tiles[360][180] = {{0}};

static float *produce_tile(int tlon, int tlat)
{
	float *t = global_table_of_tiles[tlon][tlat];
	if (!t) {
		char *fname = get_asc_tile_filename(tlon, tlat);
		if (!file_exists(fname))
			download_tile_file(tlon, tlat);
		if (!file_exists(fname))
			return NULL;
		t = malloc_tile_data(fname);
		global_table_of_tiles[tlon][tlat] = t;
	}
	return t;
}

static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	float r = 0;
	r += a * (1-x) * (1-y);
	r += b * ( x ) * (1-y);
	r += c * (1-x) * ( y );
	r += d * ( x ) * ( y );
	return r;
}

static float getpixel_1(float *x, int w, int h, int i, int j)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	return x[j*w+i];
}

static float bilinear_interpolation_at(float *x, int w, int h, float p, float q)
{
	int ip = p;
	int iq = q;
	float a = getpixel_1(x, w, h, ip  , iq  );
	float b = getpixel_1(x, w, h, ip+1, iq  );
	float c = getpixel_1(x, w, h, ip  , iq+1);
	float d = getpixel_1(x, w, h, ip+1, iq+1);
	float r = evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
	return r;
}

static float nearest_neighbor_interpolation_at(float *x,
        int w, int h, float p, float q)
{
	int ip = rintf(p);
	int iq = rintf(q);
	float r = getpixel_1(x, w, h, ip, iq);
	return r;
}


double srtm4(double lon, double lat)
{
	int tlon, tlat;
	float xlon, xlat;
	get_tile_index_and_position(&tlon, &tlat, &xlon, &xlat, lon, lat);
	float *t = produce_tile(tlon, tlat);
	return bilinear_interpolation_at(t, 6000, 6000, xlon, xlat);
}

double srtm4_nn(double lon, double lat)
{
	int tlon, tlat;
	float xlon, xlat;
	get_tile_index_and_position(&tlon, &tlat, &xlon, &xlat, lon, lat);
	float *t = produce_tile(tlon, tlat);
	return nearest_neighbor_interpolation_at(t, 6000, 6000, xlon, xlat);
}

void srtm4_free_tiles(void)
{
	for (int j = 0; j < 360; j++)
	for (int i = 0; i < 180; i++)
		if (global_table_of_tiles[j][i])
			free(global_table_of_tiles[j][i]);
}

#ifdef MAIN_SRTM4
int main(int c, char *v[])
{
	if (c != 1 && c != 3) {
		fprintf(stderr, "usage:\n\t%s longitude latitude\n", *v);
		return 1;
	}
    if (c == 3) {
	    double lon = atof(v[1]);
	    double lat = atof(v[2]);
	    double r = srtm4(lon, lat);
	    printf("%g\n", r);
	    return 0;
    }
    else {
        double lon, lat, r;
        while(2 == scanf("%lf %lf\n", &lon, &lat)) {
            r = srtm4_nn(lon, lat);
            printf("%g\n", r);
        }
    }
}
#endif//MAIN_SRTM4
