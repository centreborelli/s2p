CC = gcc -std=c99
CFLAGS = -g -O2 -DNDEBUG -DDONT_USE_TEST_MAIN
CXX = g++
CPPFLAGS = -g -O2
LDLIBS = -lpng -ltiff -ljpeg -lm
IIOLIBS = -lpng -ltiff -ljpeg -lm
FFTLIBS = -lfftw3f -lfftw3

BINDIR = bin
SRCDIR = c
GEODIR = 3rdparty/GeographicLib-1.32

default: $(BINDIR) geographiclib monasse sift imscript sgbm

all: $(BINDIR) geographiclib monasse sift imscript msmw tvl1 sgbm

$(BINDIR):
	mkdir -p $(BINDIR)

geographiclib: $(BINDIR)/CartConvert

$(BINDIR)/CartConvert: $(BINDIR) $(GEODIR)/tools/CartConvert.cpp\
	$(GEODIR)/src/DMS.cpp $(GEODIR)/src/Geocentric.cpp\
	$(GEODIR)/src/LocalCartesian.cpp
	$(CXX) $(CPPFLAGS) -I $(GEODIR)/include -I $(GEODIR)/man -o $(BINDIR)/CartConvert\
		$(GEODIR)/tools/CartConvert.cpp $(GEODIR)/src/DMS.cpp\
		$(GEODIR)/src/Geocentric.cpp $(GEODIR)/src/LocalCartesian.cpp

monasse:
	mkdir -p $(BINDIR)/monasse_refactored_build
	cd $(BINDIR)/monasse_refactored_build; cmake ../../c/monasse_refactored; make
	cp $(BINDIR)/monasse_refactored_build/bin/homography $(BINDIR)
	cp $(BINDIR)/monasse_refactored_build/bin/rectify_mindistortion $(BINDIR)
	cp $(BINDIR)/monasse_refactored_build/bin/sift_keypoints $(BINDIR)

sift: $(BINDIR)
	cd 3rdparty/sift_anatomy_20140911; make
	cp 3rdparty/sift_anatomy_20140911/bin/sift_cli $(BINDIR)
	cp 3rdparty/sift_anatomy_20140911/bin/match_cli $(BINDIR)

sgbm:
	cd 3rdparty/sgbm; make
	cp 3rdparty/sgbm/sgbm $(BINDIR)
	cp 3rdparty/sgbm/call_sgbm.sh $(BINDIR)

sgbm_opencv:
	mkdir -p bin/sgbm_build
	cd bin/sgbm_build; cmake -D CMAKE_PREFIX_PATH=~/local ../../3rdparty/stereo_hirschmuller_2008; make
	cp bin/sgbm_build/sgbm2 $(BINDIR)
	cp bin/sgbm_build/SGBM $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM.sh $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_lap.sh $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_cauchy.sh $(BINDIR)

msmw:
	mkdir -p $(BINDIR)/msmw_build
	cd $(BINDIR)/msmw_build; cmake ../../3rdparty/msmw; make
	cp $(BINDIR)/msmw_build/libstereo/iip_stereo_correlation_multi_win2 $(BINDIR)

tvl1:
	cd 3rdparty/tvl1flow_3; make
	cp 3rdparty/tvl1flow_3/tvl1flow $(BINDIR)
	cp 3rdparty/tvl1flow_3/callTVL1.sh $(BINDIR)

# WITHOUT IIO:
# bin2asc siftu srtm4 srtm4_which_tile ransac
#
$(BINDIR)/bin2asc: c/bin2asc.c
	$(CC) $(CFLAGS) c/bin2asc.c -o $(BINDIR)/bin2asc

$(BINDIR)/siftu: c/siftu.c c/siftie.c
	$(CC) $(CFLAGS) c/siftu.c -lm -o $(BINDIR)/siftu

$(BINDIR)/ransac: c/ransac.c c/fail.c c/xmalloc.c c/xfopen.c c/cmphomod.c\
	c/ransac_cases.c c/parsenumbers.c
	$(CC) $(CFLAGS) c/ransac.c -lm -o $(BINDIR)/ransac

$(BINDIR)/srtm4: c/srtm4.c Geoid.o geoid_height_wrapper.o
	$(CC) $(CFLAGS) -DMAIN_SRTM4 c/srtm4.c geoid_height_wrapper.o Geoid.o \
	-ltiff -lm -lstdc++ -o $(BINDIR)/srtm4

$(BINDIR)/srtm4_which_tile: c/srtm4.c Geoid.o geoid_height_wrapper.o
	$(CC) $(CFLAGS) -DMAIN_SRTM4_WHICH_TILE c/srtm4.c geoid_height_wrapper.o \
	Geoid.o -ltiff -lm -lstdc++ -o $(BINDIR)/srtm4_which_tile

# WITH IIO ONLY:
SRCIIO = downsa backflow synflow imprintf iion plambda qauto qeasy crop morsi\
	veco morphoop cldmask disp_to_h_projective colormesh_projective
SRCFFT = gblur blur fftconvolve zoom_zeropadding zoom_2d
SRC = $(SRCIIO) $(SRCFFT)
PROGRAMS = $(addprefix $(BINDIR)/,$(SRC))
imscript: $(BINDIR) $(PROGRAMS) $(addprefix $(BINDIR)/, watermask disp_to_h colormesh bin2asc siftu ransac srtm4 srtm4_which_tile)

$(addprefix $(BINDIR)/,$(SRCIIO)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	    $(CC) $(CFLAGS) $^ -o $@ $(IIOLIBS)

$(addprefix $(BINDIR)/,$(SRCFFT)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	    $(CC) $(CFLAGS) $^ -o $@ $(IIOLIBS) $(FFTLIBS)

$(SRCDIR)/iio.o: $(SRCDIR)/iio.c $(SRCDIR)/iio.h
	$(CC) $(CFLAGS) -c -Wno-deprecated-declarations $< -o $@

$(BINDIR)/watermask: $(SRCDIR)/iio.o Geoid.o geoid_height_wrapper.o c/watermask.c c/fail.c\
	c/xmalloc.c c/pickopt.c c/rpc.c c/srtm4.c c/iio.h c/parsenumbers.c
	$(CC) $(CFLAGS) $(SRCDIR)/iio.o Geoid.o geoid_height_wrapper.o c/watermask.c \
	$(LDLIBS) -lstdc++ -o $(BINDIR)/watermask




# WITH IIO AND  RPC:
#   disp_to_h colormesh
#
rpc.o: c/rpc.c c/xfopen.c
	$(CC) $(CFLAGS) -c c/rpc.c

$(BINDIR)/disp_to_h: $(SRCDIR)/iio.o rpc.o c/disp_to_h.c c/vvector.h c/iio.h c/rpc.h c/read_matrix.c
	$(CC) $(CFLAGS) $(SRCDIR)/iio.o rpc.o c/disp_to_h.c $(LDLIBS) -o $(BINDIR)/disp_to_h

$(BINDIR)/colormesh: $(SRCDIR)/iio.o rpc.o geographiclib_wrapper.o DMS.o GeoCoords.o MGRS.o\
	PolarStereographic.o TransverseMercator.o UTMUPS.o c/colormesh.c c/iio.h\
	c/fail.c c/rpc.h c/read_matrix.c c/smapa.h
	$(CC) $(CFLAGS) $(SRCDIR)/iio.o rpc.o geographiclib_wrapper.o DMS.o GeoCoords.o \
	MGRS.o PolarStereographic.o TransverseMercator.o UTMUPS.o c/colormesh.c \
	$(LDLIBS) -lstdc++ -o $(BINDIR)/colormesh

# GEOGRAPHICLIB STUFF
geographiclib_wrapper.o: c/geographiclib_wrapper.cpp
	$(CXX) $(CPPFLAGS) -c c/geographiclib_wrapper.cpp -I.

geoid_height_wrapper.o: c/geoid_height_wrapper.cpp
	$(CXX) $(CPPFLAGS) -c c/geoid_height_wrapper.cpp -I.

DMS.o: c/DMS.cpp
	$(CXX) $(CPPFLAGS) -c c/DMS.cpp -I.

GeoCoords.o: c/GeoCoords.cpp
	$(CXX) $(CPPFLAGS) -c c/GeoCoords.cpp -I.

MGRS.o: c/MGRS.cpp
	$(CXX) $(CPPFLAGS) -c c/MGRS.cpp -I.

PolarStereographic.o: c/PolarStereographic.cpp
	$(CXX) $(CPPFLAGS) -c c/PolarStereographic.cpp -I.

TransverseMercator.o: c/TransverseMercator.cpp
	$(CXX) $(CPPFLAGS) -c c/TransverseMercator.cpp -I.

UTMUPS.o: c/UTMUPS.cpp
	$(CXX) $(CPPFLAGS) -c c/UTMUPS.cpp -I.

Geoid.o: c/Geoid.cpp
	$(CXX) $(CPPFLAGS) -c c/Geoid.cpp -I.


clean:
	rm *.o c/*.o
	rm -r bin

.PHONY: default all geographiclib monasse sift sgbm sgbm_opencv msmw tvl1\
	imscript clean
