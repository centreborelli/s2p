// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "homography.h"

namespace libNumerics {

/// Constructor.
Homography::Homography()
: m_H( matrix<double>::eye(3) )
{}

/// Put to identity.
void Homography::setId()
{
    m_H = matrix<double>::eye(3);
}

/// Set to translation.
void Homography::setTrans(double dx, double dy)
{
    setId();
    m_H(0,2) = dx;
    m_H(1,2) = dy;
}

/// Set to zoom.
void Homography::setZoom(double zx, double zy)
{
    setId();
    m_H(0,0) = zx;
    m_H(1,1) = zy;
}

/// Apply homography.
void Homography::operator()(double& x, double& y) const
{
    vector<double> m(3);
    m(0) = x;
    m(1) = y;
    m(2) = 1.0f;
    m = m_H * m;
    double z_1 = 1.0 / m(2);
    x = m(0) * z_1;
    y = m(1) * z_1;
}

/// Compose homographies.
Homography Homography::operator*(const Homography& rhs) const
{
    Homography h;
    h.m_H = m_H * rhs.m_H;
    h.normalize();
    return h;
}

/// Inverse homography.
Homography Homography::inverse() const
{
    Homography h;
    h.m_H = m_H.inv();
    h.normalize();
    return h;
}

/// Put coef(2,2) to 1.
void Homography::normalize()
{
    m_H /= m_H(2,2);
}

} // libNumerics
