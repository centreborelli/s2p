/* Matches between points and their I/O.
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

#include "match.h"
#include <fstream>
#include <sstream>

/// Read file.
bool loadMatch(const char* nameFile, std::vector<Match>& match) {
    match.clear();
    std::ifstream f(nameFile);
    while( f.good() ) {
        std::string str;
        std::getline(f, str);
        if( f.good() ) {
            std::istringstream s(str);
            Match m;
            s >> m.x1 >> m.y1 >> m.x2 >> m.y2;
            if(!s.fail() )
                match.push_back(m);
        }
    }
    return !match.empty();
}

/// Write file.
bool saveMatch(const char* nameFile, const std::vector<Match>& match) {
    std::ofstream f(nameFile);
    if( f.is_open() ) {
        std::vector<Match>::const_iterator it = match.begin();
        for(; it != match.end(); ++it)
            f << it->x1 << " " << it->y1 << " "
              << it->x2 << " " << it->y2 << std::endl;
    }
    return f.is_open();
}
