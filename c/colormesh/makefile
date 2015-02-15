CC = gcc -std=c99 -g -lprofiler
CFLAGS = -O2 -fopenmp
LDLIBS = -lpng -ltiff -ljpeg -lm -lgeotiff -lstdc++ -lrt
SRCDIR = ..

default:
	$(CC) $(CFLAGS) $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/DMS.o $(SRCDIR)/GeoCoords.o $(SRCDIR)/MGRS.o $(SRCDIR)/PolarStereographic.o $(SRCDIR)/TransverseMercator.o $(SRCDIR)/UTMUPS.o colormesh.c timing.c $(LDLIBS) -o colormesh

clean:
	rm colormesh
