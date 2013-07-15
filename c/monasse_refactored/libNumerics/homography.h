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

#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY_H

#include "matrix.h"

namespace libNumerics {

/// 2-D homography transform.
class Homography {
public:
    Homography();

    void setId();
    void setTrans(double dx, double dy);
    void setZoom(double zx, double zy);

    matrix<double>& mat() { return m_H; }
    const matrix<double>& mat() const { return m_H; }

    void operator()(double& x, double& y) const;
    Homography operator*(const Homography& rhs) const;
    Homography inverse() const;
private:
    matrix<double> m_H;
    void normalize();
};

} // libNumerics

#endif
