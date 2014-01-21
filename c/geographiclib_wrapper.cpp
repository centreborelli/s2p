#include <string>
#include <GeographicLib/GeoCoords.hpp>

extern "C" void utm(double *out, double lat, double lon)
{
    GeographicLib::GeoCoords p(lat, lon);
    out[0] = p.Easting();
    out[1] = p.Northing();
}
