#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "pickopt.c"

void print_help(char *bin_name)
{
    fprintf(stderr, "usage:\n\t %s input_binary.ply [--strip-normals] "
            "[--strip-colors] [--strip-header] [-o output_ascii.ply]\n",
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

// prints "f_in" until "last_lin" is found. "last_lin" is printed.
// checks if there are lines starting with given patterns
static void print_until_this_line(FILE *f_out, bool *out, FILE *f_in,
        char *last_lin, char *pattern[2], bool strip[2])
{
    char buf[FILENAME_MAX] = {0};
    while (fgets(buf, FILENAME_MAX, f_in)) {
        if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) {
            fprintf(f_out, "format ascii 1.0\n");
        } else if (0 == strncmp(buf, pattern[0], strlen(pattern[0]))) {
            *out = true;
            if (!strip[0]) fprintf(f_out, "%s", buf);
        } else if (0 == strncmp(buf, pattern[1], strlen(pattern[1]))) {
            *(out+1) = true;
            if (!strip[1]) fprintf(f_out, "%s", buf);
        } else if (0 == strcmp(buf, last_lin)) {
            fprintf(f_out, "%s", buf);
            break;
        } else
            fprintf(f_out, "%s", buf);
    }
}

int main(int c, char *v[])
{
    // read args
	bool strip_n = (bool) pick_option(&c, &v, "-strip-normals", NULL);
	bool strip_c = (bool) pick_option(&c, &v, "-strip-colors", NULL);
	bool strip_h = (bool) pick_option(&c, &v, "-strip-header", NULL);
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
    char *patterns[2] = {"property uchar", "property float n"};
    bool strip_cn[2] = {strip_c, strip_n};
    bool cn[2] = {false, false};
    if (strip_h)
        eat_until_this_line(cn, f_in, "end_header\n", patterns);
    else
        print_until_this_line(f_out, cn, f_in, "end_header\n", patterns,
                strip_cn);

    bool colors = cn[0];
    bool normals = cn[1];

    // then read the binary body of the file
    // each 'line' contains a position X = (x, y, z), eventually a normal N =
    // (nx, ny, nz) and color C = (r, g, b)
    size_t n;
    float X[3];
    float N[3];
    unsigned char C[3];

    while (!feof(f_in)) {
        n = fread(X, sizeof(float), 3, f_in);
        if (normals)
            n = fread(N, sizeof(float), 3, f_in);
        if (colors)
            n = fread(C, sizeof(unsigned char), 3, f_in);
        fprintf(f_out, "%.10f %.10f %.10f", X[0], X[1], X[2]);
        if (normals & !strip_n)
            fprintf(f_out, " %.1f %.1f %.1f", N[0], N[1], N[2]);
        if (colors & !strip_c)
            fprintf(f_out, " %d %d %d\n", C[0], C[1], C[2]);
        else
            fprintf(f_out, "\n");
    }

    fclose(f_in);
    fclose(f_out);
    return 0;
}
