CXX="c++ -mavx -march=native -O3"
CXX="c++ -mavx"
CXXFLAGS=`pkg-config gdal --cflags`
CXXLIBS=`pkg-config gdal --libs`
rm -f *.o
$CXX $CXXFLAGS -c main.cpp
$CXX $CXXFLAGS -c LibImages/LibImages.cpp
$CXX $CXXFLAGS -c Utilities/Memory.cpp
$CXX $CXXFLAGS -c Utilities/Parameters.cpp
$CXX $CXXFLAGS -c Utilities/Time.cpp
$CXX $CXXFLAGS -c Utilities/Utilities.cpp
$CXX $CXXFLAGS -c LibHomography/Homography.cpp
$CXX $CXXFLAGS -c LibHomography/Splines.cpp
$CXX *.o -o homography $CXXLIBS
rm -f *.o
