CC = gcc -std=c99 -g -lprofiler
CFLAGS = -O2 -fopenmp
LDLIBS = -lpng -ltiff -ljpeg -lm -lgeotiff -lstdc++ -lrt
SRCDIR = ..
NTH=1

default:
	$(CC) $(CFLAGS) $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/DMS.o $(SRCDIR)/GeoCoords.o $(SRCDIR)/MGRS.o $(SRCDIR)/PolarStereographic.o $(SRCDIR)/TransverseMercator.o $(SRCDIR)/UTMUPS.o $(SRCDIR)/colormesh.c $(SRCDIR)/timing.c $(LDLIBS) -o colormesh

test:
	OMP_NUM_THREADS=$(NTH) ./colormesh out.ply height_map.tif rpc.xml roi_ref.png

test_nocolor:
	OMP_NUM_THREADS=$(NTH) ./colormesh out.ply height_map.tif rpc.xml

test_with_normals:
	OMP_NUM_THREADS=$(NTH) ./colormesh out.ply height_map.tif rpc.xml roi_ref.png --with-normals

valgrind:
	OMP_NUM_THREADS=$(NTH) valgrind --tool=memcheck --leak-check=full ./colormesh out.ply height_map.tif rpc.xml roi_ref.png

clean:
	rm colormesh
