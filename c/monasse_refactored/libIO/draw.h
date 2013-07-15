/* Draw dashed line in image.
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

#ifndef DRAW_H
#define DRAW_H

#ifdef __cplusplus
extern "C" {
#endif

struct color {
    unsigned char r, g, b;
};

void draw_horizontal_dashed_line(unsigned char* data, int w, int h,
                                 int y, int length, int gap, struct color c);
void draw_cross(unsigned char* data, int w, int h,
                int x, int y, int halfLength, struct color c);

#ifdef __cplusplus
}
#endif

#endif
