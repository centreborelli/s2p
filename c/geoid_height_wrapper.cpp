#include <string>
#include <exception>
#include "GeographicLib/Geoid.hpp"

extern "C" void geoid_height(double *out, double lat, double lon)
{
    GeographicLib::Geoid egm96("egm96-5", "./data");
    *out = egm96(lat, lon);
}
