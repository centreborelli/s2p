#include <string>
#include <cstdlib>
#include <GeographicLib/GeoCoords.hpp>

extern "C" void utm(double *out, double lat, double lon)
{
    GeographicLib::GeoCoords p(lat, lon);
    out[0] = p.Easting();
    out[1] = p.Northing();
}

extern "C" void utm_alt_zone(double *out, double lat, double lon, int zone)
{
    GeographicLib::GeoCoords p(lat, lon);
    p.SetAltZone(zone);
    out[0] = p.AltEasting();
    out[1] = p.AltNorthing();
}

extern "C" void utm_zone(int *zone, bool *northp, double lat, double lon)
{
    GeographicLib::GeoCoords p(lat, lon);
    zone[0] = p.Zone();
    northp[0] = p.Northp();
}

// returns the signed zone
extern "C"
int utm_from_lonlat(double out_eastnorth[2], double lon, double lat)
{
    GeographicLib::GeoCoords p(lat, lon);
    out_eastnorth[0] = p.Easting();
    out_eastnorth[1] = p.Northing();
    return p.Zone() * (p.Northp() ? 1 : -1);
}

// expects a signed zone
extern "C"
void lonlat_from_eastnorthzone(double out_lonlat[2], double e, double n, int z)
{
    GeographicLib::GeoCoords p(std::abs(z), z>0, e, n);
    out_lonlat[0] = p.Longitude();
    out_lonlat[1] = p.Latitude();
}
