CXXFLAGS ?= -march=native -O3
override CXXFLAGS := $(CXXFLAGS) `gdal-config --cflags` -std=c++11
override CXXLIBS  := $(CXXLIBS) `gdal-config --libs`

OBJ = main.o \
      LibImages/LibImages.o \
      Utilities/Utilities.o \
      Utilities/Memory.o \
      Utilities/Parameters.o \
      Utilities/Time.o \
      LibHomography/Homography.o \
      LibHomography/Splines.o

homography: $(OBJ)
	$(CXX) -o $@ $^ $(CXXLIBS)

clean:
	rm -f homography $(OBJ)
