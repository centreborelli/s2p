#ifndef _COORDCONVERT_H
#define _COORDCONVERT_H

// convert (long,lat,h) to ECEF coord sys. (X,Y,Z)
void geotedic_to_ECEF(double lg, double lt, double h, 
					  double *X, double *Y, double *Z);

// given (X,Y,Z) ECEF coord, computes the alt above the WGS 84 ellipsoid
// (faster than ECEF_to_lgt_lat_alt, but only gives the altitude)
double get_altitude_from_ECEF(double X, double Y, double Z);

// given (X,Y,Z) ECEF coord, computes the lat/long/alt 
// above the WGS 84 ellipsoid
void ECEF_to_lgt_lat_alt(double X, double Y, double Z,
						 double *lgt, double *lat, double *h);

#endif // _COORDCONVERT_H
