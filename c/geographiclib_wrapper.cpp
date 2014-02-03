#include <string>
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
