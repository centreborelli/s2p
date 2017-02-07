CFLAGS   ?= -march=native -O3
CXXFLAGS ?= -march=native -O3

LDLIBS = -lpng -ltiff -ljpeg

mgm: iio.o mgm.cc img.cc point.cc
	$(CXX) $(CXXFLAGS) -I. -Icommon $^ -o $@ $(LDLIBS)

clean:
	$(RM) -f mgm iio.o
