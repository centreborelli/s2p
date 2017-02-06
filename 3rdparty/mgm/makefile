LDFLAGS = -lpng -ltiff -ljpeg

mgm: iio.o mgm.cc img.cc point.cc
	$(CXX) $(CXXFLAGS) -I. -Icommon $^ -o $@ $(LDFLAGS)

clean:
	$(RM) -f mgm iio.o
