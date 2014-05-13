#include <stdio.h>
#include <string.h>

#include "pickopt.c"

void print_help(char *bin_name)
{
    fprintf(stderr, "usage:\n\t" "%s binary_ply [--strip-normals] >"
            " ascii_ply\n", bin_name);
}

// fast forward "f" until "lin" is found. "lin" is read.
static void eat_until_this_line(FILE *f, char *lin)
{
    char buf[FILENAME_MAX] = {0};
    while (fgets(buf, FILENAME_MAX, f))
        if (0 == strcmp(buf, lin))
            return;
}

// prints "f" on stdout until "lin" is found. Does not print "lin".
static void print_until_this_line(FILE *f, char *lin)
{
    char buf[FILENAME_MAX] = {0};
    while (fgets(buf, FILENAME_MAX, f))
        if (0 == strcmp(buf, lin))
            return;
        else
            printf("%s", buf);
}

int main(int c, char *v[])
{
    if (c != 2 && c != 3) {
        print_help(*v);
        return 1;
    }

	char *strip_n = pick_option(&c, &v, "-strip-normals", NULL);

    // read args
    const char *filename = v[1];

    // open the input binary file
    FILE *f = fopen(filename, "rb");

    // copy the header, replace the 'format' line, and remove the line about
    // normals if needed
    print_until_this_line(f, "format binary_little_endian 1.0\n");
    printf("format ascii 1.0\n");
    if (strip_n) {
        print_until_this_line(f, "property float nx\n");
        eat_until_this_line(f, "property float nz\n");
    }
    print_until_this_line(f, "end_header\n");
    printf("end_header\n");

    // then read the binary body of the file
    // each 'line' contains position X = (x, y, z), normal N = (nx, ny, nz) and
    // color C = (r, g, b)
    float X[3];
    float N[3];
    unsigned char C[3];
    while (!feof(f)) {
        fread(X, sizeof(float), 3, f);
        fread(N, sizeof(float), 3, f);
        fread(C, sizeof(unsigned char), 3, f);
        printf("%.10f %.10f %.10f ", X[0], X[1], X[2]);
        if (!strip_n)
            printf("%.1f %.1f %.1f ", N[0], N[1], N[2]);
        printf("%d %d %d\n", C[0], C[1], C[2]);
    }

    fclose(f);
    return 0;
}
