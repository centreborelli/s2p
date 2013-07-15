/* Bounding box of transformed rectangle.
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

#ifndef TRANSFORMSIZE_H
#define TRANSFORMSIZE_H

#include "libNumerics/homography.h"

/// Compute bounding box of transformed rectangle(s).
///
/// This is used to compute the size of an image transformed by an
/// \c Homography, or of a mosaic.
class TransformSize
{
public:
    TransformSize();

    void add(const libNumerics::Homography& map, int w, int h);
    int x() const;
    int y() const;
    int w() const;
    int h() const;
private:
    float x1, y1, x2, y2; ///< Rectangle coordinates
};

/// Left coordinate of bounding box.
inline int TransformSize::x() const
{ return (int)x1; }

/// Top coordinate of bounding box.
inline int TransformSize::y() const
{ return (int)y1; }

/// Width of bounding box.
inline int TransformSize::w() const
{ return (x2 < x1)? 0: int(x2-x1+1); }

/// Height of bounding box.
inline int TransformSize::h() const
{ return (y2 < y1)? 0: int(y2-y1+1); }

#endif
