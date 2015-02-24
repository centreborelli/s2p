#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "pickopt.c"

void print_help(char *bin_name)
{
    fprintf(stderr, "usage:\n\t %s binary_ply [--strip-normals] "
            "[--strip-colors] [--strip-header] > ascii_ply\n", bin_name);
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

// prints "f" until "last_lin" is found. "last_lin" is printed.
// checks if there are lines starting with given patterns
static void print_until_this_line(bool *out, FILE *f, char *last_lin,
        char *pattern[2], bool strip[2])
{
    char buf[FILENAME_MAX] = {0};
    while (fgets(buf, FILENAME_MAX, f)) {
        if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) {
            printf("format ascii 1.0\n");
        } else if (0 == strncmp(buf, pattern[0], strlen(pattern[0]))) {
            *out = true;
            if (!strip[0]) printf("%s", buf);
        } else if (0 == strncmp(buf, pattern[1], strlen(pattern[1]))) {
            *(out+1) = true;
            if (!strip[1]) printf("%s", buf);
        } else if (0 == strcmp(buf, last_lin)) {
            printf("%s", buf);
            break;
        } else
            printf("%s", buf);
    }
}

int main(int c, char *v[])
{
    if (c != 2 && c != 3 && c != 4 && c != 5) {
        print_help(*v);
        return 1;
    }

    // read args
	bool strip_n = (pick_option(&c, &v, "-strip-normals", NULL) != NULL);
	bool strip_c = (pick_option(&c, &v, "-strip-colors", NULL) != NULL);
	bool strip_h = (pick_option(&c, &v, "-strip-header", NULL) != NULL);
    const char *filename = v[1];

    // open the input binary file, and parse its header
    FILE *f = fopen(filename, "rb");
    char *patterns[2] = {"property uchar", "property float n"};
    bool strip_cn[2] = {strip_c, strip_n};
    bool cn[2] = {false, false};
    if (strip_h)
        eat_until_this_line(cn, f, "end_header\n", patterns);
    else
        print_until_this_line(cn, f, "end_header\n", patterns, strip_cn);

    bool colors = cn[0];
    bool normals = cn[1];

    // then read the binary body of the file
    // each 'line' contains a position X = (x, y, z), eventually a normal N =
    // (nx, ny, nz) and color C = (r, g, b)
    size_t n;
    float X[3];
    float N[3];
    unsigned char C[3];

    while (!feof(f)) {
        n = fread(X, sizeof(float), 3, f);
        if (normals)
            n = fread(N, sizeof(float), 3, f);
        if (colors)
            n = fread(C, sizeof(unsigned char), 3, f);
        printf("%.10f %.10f %.10f", X[0], X[1], X[2]);
        if (normals & !strip_n)
            printf(" %.1f %.1f %.1f", N[0], N[1], N[2]);
        if (colors & !strip_c)
            printf(" %d %d %d\n", C[0], C[1], C[2]);
        else
            printf("\n");
    }

    fclose(f);
    return 0;
}
