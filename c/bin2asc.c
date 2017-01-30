#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "pickopt.c"

void print_help(char *bin_name)
{
    fprintf(stderr, "usage:\n\t %s input_binary.ply [-o output_ascii.ply]\n",
            bin_name);
}

// fast forward "f" until "last_lin" is found. "last_lin" is read.
// checks if there are lines starting with given patterns
static void eat_until_this_line(bool *out, FILE *f, char *last_lin,
        char *pattern[2])
{
    char buf[FILENAME_MAX] = {0};
    while (fgets(buf, FILENAME_MAX, f)) {
        if (0 == strcmp(buf, last_lin))
            break;
        else if (0 == strncmp(buf, pattern[0], strlen(pattern[0])))
            *out = true;
        else if (0 == strncmp(buf, pattern[1], strlen(pattern[1])))
            *(out+1) = true;
    }
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


// prints "f_in" until "last_lin" is found. "last_lin" is printed.
// checks if there are lines starting with given patterns
static int process_header(FILE *f_out, FILE *f_in, struct ply_property *t)
{
	int n = 0;
	char buf[FILENAME_MAX] = {0};
	while (fgets(buf, FILENAME_MAX, f_in)) {
		if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) {
			fprintf(f_out, "format ascii 1.0\n");
		} else {
			if (parse_property_line(t+n, buf))
				n += 1;
			fprintf(f_out, "%s", buf);
		}
		if (0 == strcmp(buf, "end_header\n"))
			break;
	}
	return n;
}

int main(int c, char *v[])
{
	char *output_file = pick_option(&c, &v, "o", "none");
	if (c != 2) {
		print_help(*v);
		return 1;
	}
	const char *input_file = v[1];
	//printf("strip_n %s\n", strip_n ? "true" : "false");
	//printf("strip_c %s\n", strip_c ? "true" : "false");
	//printf("strip_h %s\n", strip_h ? "true" : "false");

	// open the output file
	FILE *f_out = stdout;
	if (strcmp(output_file, "none"))
		f_out = fopen(output_file, "w");

	// open the input binary file, and parse its header
	FILE *f_in = fopen(input_file, "rb");
	struct ply_property t[100];
	int nproperties = process_header(f_out, f_in, t);


	// then read the binary body of the file
	// each 'line' contains a position X = (x, y, z), eventually a normal N =
	// (nx, ny, nz) and color C = (r, g, b)

	while (!feof(f_in)) {
		for (int i = 0; i < nproperties; i++)
		{
			switch(t[i].type) {
			case UCHAR: {
				unsigned char X;
				if( fread(&X, 1, 1, f_in))
					fprintf(f_out, "%d", X);
				break; }
			case FLOAT: {
				float X;
				if( fread(&X, sizeof(float), 1, f_in))
					fprintf(f_out, "%a", X);
				break; }
			case DOUBLE: {
				double X;
				if( fread(&X, sizeof(double), 1, f_in))
					fprintf(f_out, "%a", X);
				break; }
			default: break;
			}
			fprintf(f_out, "%c", i==nproperties-1?'\n':' ');
		}
	}
	fclose(f_in);
	fclose(f_out);
	return 0;
}
