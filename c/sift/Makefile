CXXFLAGS = -march=native -O3 `gdal-config --cflags`
CXXLIBS = `gdal-config --libs`

CFLAGS = -std=gnu99 -march=native -O3
LDLIBS = -lm

OBJ = sift_roi.o \
      LibImages/LibImages.o \
      LibSSE/LibSSE.o \
      LibSift/KeyPoint.o \
      LibSift/LibSift.o \
      LibSift/Octave.o \
      LibSift/ScaleSpace.o \
      Utilities/Memory.o \
      Utilities/Parameters.o \
      Utilities/Time.o \
      Utilities/Utilities.o \

all: sift_roi matching

sift_roi: $(OBJ)
	$(CXX) -o $@ $^ $(CXXLIBS)

clean:
	rm -f sift_roi matching $(OBJ)
