// take a series of ply files and produce a digital elevation map

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "iio.h"
#include "lists.c"
#include "fail.c"
#include "xmalloc.c"


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

static void update_min_max(double *min, double *max, double x)
{
	if (x < *min) *min = x;
	if (x > *max) *max = x;
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
				default: break;
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
static bool parse_ply_points_for_extrema(double *xmin, double *xmax, double *ymin,
		double *ymax, char *utm, char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING: can not open file \"%s\"\n", fname);
		return false;
	}

	int isbin=0;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);

    if (n > 0) {
	    double data[n];
	    while ( n == get_record(f, isbin, t, n, data) ) {
	    	update_min_max(xmin, xmax, data[0]);
	    	update_min_max(ymin, ymax, data[1]);
	    }
	    fclose(f);
        return true;
    } else {
	    fclose(f);
        return false;
    }
}


void help(char *s)
{
	fprintf(stderr, "usage:\n\t"
		"%s cloud.ply plyextrema.txt\n", s);
}

#include "pickopt.c"

int main(int c, char *v[])
{
	// process input arguments
	if (c != 3) {
		help(*v);
		return 1;
	}

	// initialize x, y extrema values
	double xmin = INFINITY;
	double xmax = -INFINITY;
	double ymin = INFINITY;
	double ymax = -INFINITY;

	char utm[3];
	
	char *ply_file = v[1];
	char *plyextrema_file = v[2];
	
	if (parse_ply_points_for_extrema(&xmin, &xmax, &ymin, &ymax, utm, ply_file))
	  fprintf(stderr, "xmin: %lf, xmax: %lf, ymin: %lf, ymax: %lf\n",
		  xmin, xmax, ymin, ymax);
	else {
	  fprintf(stderr, "plyextrema: empty input ply file\n");
	  return EXIT_FAILURE;
	}


	FILE* out_extrema = NULL;
	out_extrema = fopen(plyextrema_file, "w");
	if (out_extrema != NULL)
	{
            fprintf(out_extrema, "%lf %lf %lf %lf", xmin, xmax, ymin, ymax);
	    fclose(out_extrema);
	}
	else
	{
	    fprintf(stderr,"ERROR : can't create %s",v[2]);
	    return 1;
	}
	return 0;
}
