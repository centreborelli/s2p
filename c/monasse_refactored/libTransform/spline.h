/* Spline interpolation.
    Copyright (C) 2007 Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
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

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "libLWImage/LWImage.h"
#include "libIO/nan.h"

bool prepare_spline(LWImage<float>& im, int order);
bool interpolate_spline(LWImage<float>& im, int order,
                        float x, float y,
                        float* out,
                        float paramKeys=-.5f);

#endif
