#include <stdio.h>
#include <stdlib.h>

int utm_from_lonlat(double out_eastnorth[2], double lon, double lat);
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z);

int main(int c, char *v[])
{
	if (c == 3) {
		double lon = atof(v[1]);
		double lat = atof(v[2]);
		double utm[2];
		int z = utm_from_lonlat(utm, lon, lat);
		printf("%lf %lf %d\n", utm[0], utm[1], z);
		return 0;
	}

	if (c == 4) {
		double e = atof(v[1]);
		double n = atof(v[2]);
		int z = atoi(v[3]);
		double lonlat[2];
		lonlat_from_eastnorthzone(lonlat, e, n, z);
		printf("%lf %lf\n", lonlat[0], lonlat[1]);
		return 0;
	}

	return fprintf(stderr, "usage:\n\t%s {lon lat|east north zone}\n", *v);
}
