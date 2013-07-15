/* Class for planar homography transform.
    Copyright (C) 2010 Pascal Monasse <monasse@imagine.enpc.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
