// take a series of ply files and produce a digital elevation map

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <geotiff/xtiffio.h>
#include <geotiff/geotiffio.h>
#include <geotiff/geo_tiffp.h>

#include "iio.h"
#include "lists.c"
#include "fail.c"
#include "xmalloc.c"


static void parse_utm_string(int *zone, bool *hem, char *s)
{
	char hem_string[FILENAME_MAX];
	if (2 == sscanf(s, "%02d%s", zone, hem_string)) {
		// hem_string must be equal to "N" or "S"
		if (hem_string[1] == '\0')
			if (hem_string[0] == 'N' || hem_string[0] == 'S') {
				*hem = (hem_string[0] == 'N');
				fprintf(stderr, "zone: %d\them: %s\n", *zone, hem_string);
				return;
			}
	}
}

// convert string like '28N' into a number like 32628, according to:
// WGS84 / UTM northern hemisphere: 326zz where zz is UTM zone number
// WGS84 / UTM southern hemisphere: 327zz where zz is UTM zone number
// http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.3.3.1
static int get_utm_zone_index_for_geotiff(char *utm_zone)
{
	int out = 32000;
	if (utm_zone[2] == 'N')
		out += 600;
	else if (utm_zone[2] == 'S')
		out += 700;
	else
		fprintf(stderr, "error: bad utm zone value: %s\n", utm_zone);
	utm_zone[2] = '\0';
	out += atoi(utm_zone);
	return out;
}

void set_geotif_header(char *tiff_fname, char *utm_zone, float xoff,
		float yoff, float scale)
{
	// open tiff file
	TIFF *tif = XTIFFOpen(tiff_fname, "r+");
	if (!tif)
		fail("failed in XTIFFOpen\n");

	GTIF *gtif = GTIFNew(tif);
	if (!gtif)
		fail("failed in GTIFNew\n");

	// write TIFF tags
	double pixsize[3] = {scale, scale, 0.0};
	TIFFSetField(tif, GTIFF_PIXELSCALE, 3, pixsize);

	double tiepoint[6] = {0.0, 0.0, 0.0, xoff, yoff, 0.0};
	TIFFSetField(tif, GTIFF_TIEPOINTS, 6, tiepoint);

	// write GEOTIFF keys
	int utm_ind = get_utm_zone_index_for_geotiff(utm_zone);
	GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, utm_ind);
	GTIFWriteKeys(gtif);

	// free and close
	GTIFFree(gtif);
	XTIFFClose(tif);
}



struct ply_property {
	enum {UCHAR,FLOAT,DOUBLE,UNKNOWN} type;
	char name[0x100];
	size_t len;
};

static bool parse_property_line(struct ply_property *t, char *buf)
{
	char typename[0x100];
	bool r = 2 == sscanf(buf, "property %s %s\n", typename, t->name);
	t->type = UNKNOWN;
	if (0 == strcmp(typename, "uchar")) { t->type = UCHAR;  t->len = 1;}
	if (0 == strcmp(typename, "float")) { t->type = FLOAT;  t->len = 4;}
	if (0 == strcmp(typename, "double")){ t->type = DOUBLE; t->len = 8;}
	return r;
}



// fast forward "f" until "end_header" is found
// returns the number of 'properties' 
// the array of structures *t, contains the names and sizes 
// the properties in bytes, isbin is set if binary encoded
// and reads the utm zone
static size_t header_get_record_length_and_utm_zone(FILE *f_in, char *utm, 
		int *isbin, struct ply_property *t)
{
	size_t n = 0;
	*isbin = 0;

	char buf[FILENAME_MAX] = {0};
	while (fgets(buf, FILENAME_MAX, f_in)) {
		if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) *isbin=1;
		else if (0 == strcmp(buf, "format ascii 1.0\n")) *isbin=0;
		else {
			if (parse_property_line(t+n, buf))
				n += 1;
			else if (0 == strncmp(buf, "comment projection:", 19)) {
				sscanf(buf, "comment projection: UTM %s", utm);
			}
		}
		if (0 == strcmp(buf, "end_header\n"))
			break;
	}
	return n;
}

static void update_min_max(float *min, float *max, float x)
{
	if (x < *min) *min = x;
	if (x > *max) *max = x;
}

// re-scale a float between 0 and w
static int rescale_float_to_int(double x, double min, double max, int w, bool *flag)
{
	*flag=1;
	int r = w * (x - min)/(max - min);
	if (r < 0) *flag=0;
	if (r >= w) *flag=0;
	return r;
}


struct images {
	float *cnt;
	float *avg;
	int w, h;
};

// update the output images with a new height
static void add_height_to_images(struct images *x, int i, int j, float v)
{
	uint64_t k = (uint64_t) x->w * j + i;
	x->avg[k] = (v + x->cnt[k] * x->avg[k]) / (1 + x->cnt[k]);
	x->cnt[k] += 1;
}

int get_record(FILE *f_in, int isbin, struct ply_property *t, int n, double *data){
	int rec = 0;
	if(isbin) {
		for (int i = 0; i < n; i++) {
			switch(t[i].type) {
				case UCHAR: {
						    unsigned char X;
						    rec += fread(&X, 1, 1, f_in);
						    data[i] = X;
						    break; }
				case FLOAT: {
						    float X;
						    rec += fread(&X, sizeof(float), 1, f_in);
						    data[i] = X;
						    break; }
				case DOUBLE: {
						     double X;
						     rec += fread(&X, sizeof(double), 1, f_in);
						     data[i] = X;
						     break; }
			}
		}
	} else {
		int i=0;
		while (i < n && !feof(f_in)) {
			rec += fscanf(f_in,"%lf", &data[i]);  i++;
		}
	}
	return rec;
}

// open a ply file, read utm zone in the header, and update the known extrema
static void parse_ply_points_for_extrema(float *xmin, float *xmax, float *ymin,
		float *ymax, char *utm, char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return;
	}

	int isbin=0;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	//fprintf(stderr, "%d\n", n);
	//fprintf(stderr, "%s\n", utm);

	double data[n];
	while ( n == get_record(f, isbin, t, n, data) ) {
		update_min_max(xmin, xmax, data[0]);
		update_min_max(ymin, ymax, data[1]);
	}
	fclose(f);
}

// open a ply file, and accumulate its points to the image
static void add_ply_points_to_images(struct images *x,
		float xmin, float xmax, float ymin, float ymax,
		char utm_zone[3], char *fname, int col_idx)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return;
	}

	// check that the utm zone is the same as the provided one
	char utm[3];
	int isbin=1;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	if (0 != strncmp(utm_zone, utm, 3))
		fprintf(stderr, "error: different UTM zones among ply files\n");

	if (col_idx < 2 || col_idx > 5)
		exit(fprintf(stderr, "error: bad col_idx %d\n", col_idx));


	double data[n];
	while ( n == get_record(f, isbin, t, n, data) ) {
		bool flag1,flag2;
		int i = rescale_float_to_int(data[0], xmin, xmax, x->w, &flag1);
		int j = rescale_float_to_int(-data[1], -ymax, -ymin, x->h, &flag2);
		
		if ( (flag1) && (flag2))
		{
		    if (col_idx == 2) {
			    add_height_to_images(x, i, j, data[2]);
			    assert(isfinite(data[2]));
		    }
		    else
		    {
			    unsigned int rgb = data[col_idx];
			    add_height_to_images(x, i, j, rgb);
		    }
		}
	}

	fclose(f);
}


void help(char *s)
{
	fprintf(stderr, "usage:\n\t"
			"%s [-c column] [-bb \"xmin xmax ymin ymax\"] resolution out_dir cutting_info cloud_dir number_of_tiles\n", s);
	fprintf(stderr, "\t the resolution is in meters per pixel\n");
}

#include "pickopt.c"

int main(int c, char *v[])
{
	int col_idx = atoi(pick_option(&c, &v, "c", "2"));
	char *bbminmax = pick_option(&c, &v, "bb", "");

	// process input arguments
	if (c != 6) {
		help(*v);
		return 1;
	}
	float resolution = atof(v[1]);
	char *out_dir = v[2];
	
	int tw,th,rowmin,rowmax,steprow,colmin,colmax,stepcol;
	FILE* cutting_info = NULL;
	cutting_info = fopen(v[3], "r");
	if (cutting_info != NULL)
	{
	    // On peut lire et écrire dans le fichier
	    fscanf(cutting_info,"%d %d %d %d %d %d %d %d",&tw,&th,&rowmin,&rowmax,&steprow,&colmin,&colmax,&stepcol);
	    fclose(cutting_info);
	}
	else
	{
	    // On affiche un message d'erreur si on veut
	    fprintf(stderr,"ERROR : can't read %s",v[3]);
	    return 1;
	}
	
	// initialize x, y extrema values
	float xmin = INFINITY;
	float xmax = -INFINITY;
	float yymin = INFINITY;
	float yymax = -INFINITY;
	float ymin = INFINITY;
	float ymax = -INFINITY;

	// process each filename from stdin to determine x, y extremas and store the
	// filenames in a list of strings, to be able to open the files again
	
	char utm[3];
	struct list *l = NULL;
	char filename[1000];
	for(int c=colmin; c<=colmax; c+=stepcol)
	    for(int r=rowmin; r<=rowmax; r+=steprow)
	    {
		sprintf(filename,"%s/cloud_%d_%d_row_%d_col_%d.ply",v[4],tw,th,r,c);
		strtok(filename, "\n");
		l = push(l, filename);
		parse_ply_points_for_extrema(&xmin, &xmax, &yymin, &yymax, utm, filename);
	    }
	
	if (0 != strcmp(bbminmax, "") ) {
		sscanf(bbminmax, "%f %f %f %f", &xmin, &xmax, &yymin, &yymax);
	}
	fprintf(stderr, "xmin: %20f, xmax: %20f, ymin: %20f, ymax: %20f\n", xmin, xmax, yymin, yymax);

	struct list *begin=l;
	
	int piece_index=0,n=atoi(v[5]);
	while (yymin+(piece_index+1)*(yymax-yymin)/( (float) n) <= yymax)
	{
	    ymin = yymin+piece_index*(yymax-yymin)/( (float) n);
	    ymax = yymin+(piece_index+1)*(yymax-yymin)/( (float) n);
	
	    // compute output image dimensions
	    int w = 1 + (xmax - xmin) / resolution;
	    int h = 1 + (ymax - ymin) / resolution;

	    // allocate and initialize output images
	    struct images x;
	    x.w = w;
	    x.h = h;
	    x.cnt = xmalloc(w*h*sizeof(float));
	    x.avg = xmalloc(w*h*sizeof(float));
	    for (uint64_t i = 0; i < (uint64_t) w*h; i++)
	    {
		    x.cnt[i] = 0;
		    x.avg[i] = 0;
	    }

	    char looputm[3];
	    strcpy(looputm,utm);

	    // process each filename to accumulate points in the dem
	    l=begin;
	    while (l != NULL)
	    {
		    // printf("FILENAME: \"%s\"\n", l->current);
		    add_ply_points_to_images(&x, xmin, xmax, ymin, ymax, looputm, l->current, col_idx);
		    l = l->next;
	    }

	    // set unknown values to NAN
	    for (uint64_t i = 0; i < (uint64_t) w*h; i++)
		    if (!x.cnt[i])
			    x.avg[i] = NAN;

	    // save output image
	    char out[1000];
	    sprintf(out,"%s/dsm_%d.tif",out_dir,piece_index);
	    iio_save_image_float(out, x.avg, w, h);
	    set_geotif_header(out, looputm, xmin, ymax, resolution);

	    // cleanup and exit
	    free(x.cnt);
	    free(x.avg);
	    
	    piece_index++; 
	}    
	return 0;
}
