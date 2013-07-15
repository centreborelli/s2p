/* Apply homography to an image.
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

#ifndef MAP_IMAGE_H
#define MAP_IMAGE_H

#include "libNumerics/homography.h"
#include "libLWImage/LWImage.h"
#include "libIO/nan.h"
#include <utility>

std::pair<int,int> map_image(LWImage<float> in,
                             libNumerics::Homography map,
                             LWImage<float>& im,
                             int order=1, bool adjustSize=true,
                             bool antiAlias=true,
                             float vOut=NaN);

#endif
