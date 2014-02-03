#include <stdio.h>

void print_help(char *bin_name) {
    fprintf(stderr, "usage:\n\t" "%s binary_ply > ascii_ply\n", bin_name);
}

int main(int c, char *v[])
{
    if (c != 2) {
        print_help(*v);
        return 1;
    }

    // read args
    const char *filename = v[1];

    // open the input binary file
    FILE *f = fopen(filename, "rb");

    // copy the header from line 3 to 15 (yes, it's not really generic, but
    // this util is specific for ply files written by the colormesh binary)
    // first skip the first two lines of the file, and write them to output
    char buf[100];
    fgets(buf, 100, f);
    fgets(buf, 100, f);
    printf("ply\n");
    printf("format ascii 1.0\n");
    for (int i=2; i < 15; ++i) {
        fgets(buf, 100, f);
        printf("%s", buf);
    }

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
        printf("%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %d %d %d\n", X[0],
                X[1], X[2], N[0], N[1], N[2], C[0], C[1], C[2]);
    }

    fclose(f);
    return 0;
}
