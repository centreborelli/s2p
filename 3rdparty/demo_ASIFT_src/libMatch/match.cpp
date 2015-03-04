// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "match.h"
#include <fstream>
#include <sstream>

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
