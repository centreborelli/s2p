#include "coordconvert.h"
#include "math.h"

// Some numerical values related to WGS 84 ellipsoid
#define A 6378137           // semi major axis
#define E 0.081819190842622 // first eccentricity

// convert (long,lat,h) to ECEF coord sys. (X,Y,Z)
void geotedic_to_ECEF(double lg, double lt, double h, 
					  double *X, double *Y, double *Z)
{
	double pi = acos(-1);
	lg = lg*pi/180.;
	lt = lt*pi/180.;
	
	double N = A/sqrt(1.-E*E*sin(lt)*sin(lt));
	*X=(N+h)*cos(lt)*cos(lg); 
	*Y=(N+h)*cos(lt)*sin(lg); 
	*Z=(N*(1.-E*E)+h)*sin(lt);
}

// given (X,Y,Z) ECEF coord, computes the alt above the WGS 84 ellipsoid
// (faster than ECEF_to_lgt_lat_alt, but only gives the altitude)
double get_altitude_from_ECEF(double X, double Y, double Z)
{
	double p = sqrt(X*X+Y*Y);
	double k0 = 1./(1.-E*E);
	double c0 = pow(p*p+(1.-E*E)*Z*Z*k0*k0,1.5)/(A*E*E);
	double k1 = k0; 
	double c1 = c0;
	double h0=-10000.,h;
	double diff=10000.;
	
	// newton-raphson; most of the time,
	// only takes 2 iterations
	int max_nb_iter=15,nb_iter=0;
	while( (fabs(diff)>0.001) && (nb_iter <= max_nb_iter) ) 
	{
		k1 = 1. + (p*p+(1.-E*E)*Z*Z*k1*k1*k1)/(c1-p*p);
		c1 = pow(p*p+(1.-E*E)*Z*Z*k1*k1,1.5)/(A*E*E);
		h = sqrt(p*p+Z*Z*k1*k1)*(1./k1-1./k0)/(E*E);
		diff=h0-h;
		h0=h;
		
		nb_iter++;
	}
	
	return h;
}

// given (X,Y,Z) ECEF coord, computes the lat/long/alt 
// above the WGS 84 ellipsoid
void ECEF_to_lgt_lat_alt(double X, double Y, double Z,
						 double *lgt, double *lat, double *h)
{         
	double E2=E*E;
	double deg2rad=atan2(1,1)/45.;
	double epsilon_phi = 1e-12;
	double kpi=45./atan2(1,1);

	double R=sqrt(X*X+Y*Y);
	double phi1 = atan2(Z,R*(1-A*E2/sqrt(X*X+Y*Y+Z*Z)));
	double phi0 = phi1 + 10.*epsilon_phi;

	// altitude
	int max_nb_iter=15,nb_iter=0;
	while ( (fabs(phi1-phi0) > epsilon_phi) && (nb_iter <= max_nb_iter) ) 
	{
		phi0 = phi1;
      	phi1 = atan2(Z/R,1-A*E2*cos(phi0)/(R*sqrt(1-E2*sin(phi0)*sin(phi0))));
		*h=R/cos(phi1)-A/sqrt(1-E2*sin(phi1)*sin(phi1));
		nb_iter++;
	}
	
	// longitude
	*lgt = atan2(Y,X)*kpi;
	
	//latitude
	*lat = phi1*kpi;
}
