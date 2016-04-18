C99 = $(CC) -std=c99
CFLAGS = -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN
CPPFLAGS = -g -O3
LDLIBS = -lstdc++
IIOLIBS = $(TIFDIR)/lib/libtiff.a -lz -lpng -ljpeg -lm
GEOLIBS = -lgeotiff -ltiff
FFTLIBS = -lfftw3f -lfftw3

ifeq ($(CC), gcc)
	CFLAGS += -fopenmp
endif

OS := $(shell uname -s)
ifeq ($(OS), Linux)
	LDLIBS += -lrt
endif

BINDIR = bin
SRCDIR = c
GEODIR = 3rdparty/GeographicLib-1.32
TIFDIR = 3rdparty/tiff-4.0.4beta

default: $(BINDIR) libtiff geographiclib monasse homography sift asift imscript mgm msmw3 sgbm piio

all: default tvl1


$(BINDIR):
	mkdir -p $(BINDIR)


libtiff: $(TIFDIR)/lib/libtiff.a $(BINDIR)/raw2tiff $(BINDIR)/tiffinfo

$(TIFDIR)/lib/libtiff.a:
	cd $(TIFDIR); ./configure --prefix=`pwd` --disable-lzma --disable-jbig --with-pic; $(MAKE) install -j
	cd $(TIFDIR); $(MAKE) distclean

$(BINDIR)/raw2tiff: $(TIFDIR)/lib/libtiff.a $(BINDIR)
	cp $(TIFDIR)/bin/raw2tiff $(BINDIR)

$(BINDIR)/tiffinfo: $(TIFDIR)/lib/libtiff.a $(BINDIR)
	cp $(TIFDIR)/bin/tiffinfo $(BINDIR)

piio: python/piio/libiio.so
python/piio/libiio.so: python/piio/setup.py python/piio/freemem.c c/iio.c c/iio.h
	cd python/piio; python setup.py build

geographiclib: $(BINDIR) $(BINDIR)/CartConvert $(BINDIR)/GeoConvert

$(BINDIR)/CartConvert: $(GEODIR)/tools/CartConvert.cpp $(GEODIR)/src/DMS.cpp\
	$(GEODIR)/src/Geocentric.cpp $(GEODIR)/src/LocalCartesian.cpp
	$(CXX) $(CPPFLAGS) -I $(GEODIR)/include -I $(GEODIR)/man $^ -o $@

$(BINDIR)/GeoConvert: $(GEODIR)/tools/GeoConvert.cpp $(GEODIR)/src/DMS.cpp\
	$(GEODIR)/src/GeoCoords.cpp $(GEODIR)/src/MGRS.cpp\
	$(GEODIR)/src/PolarStereographic.cpp $(GEODIR)/src/TransverseMercator.cpp\
	$(GEODIR)/src/UTMUPS.cpp
	$(CXX) $(CPPFLAGS) -I $(GEODIR)/include -I $(GEODIR)/man $^ -o $@

asift:
	mkdir -p $(BINDIR)/build_asift
	cd $(BINDIR)/build_asift; cmake ../../3rdparty/demo_ASIFT_src; $(MAKE)
	cp $(BINDIR)/build_asift/demo_ASIFT $(BINDIR)

monasse:
	mkdir -p $(BINDIR)/build_monasse_refactored
	cd $(BINDIR)/build_monasse_refactored; cmake ../../c/monasse_refactored; $(MAKE)
	cp $(BINDIR)/build_monasse_refactored/bin/rectify_mindistortion $(BINDIR)

homography: $(BINDIR)
	mkdir -p $(BINDIR)/build_homography
	cd $(BINDIR)/build_homography; cmake ../../c/homography; $(MAKE)
	cp $(BINDIR)/build_homography/homography $(BINDIR)

sift: $(BINDIR)
	mkdir -p $(BINDIR)/build_sift
	cd $(BINDIR)/build_sift; cmake ../../c/sift; $(MAKE)
	cp $(BINDIR)/build_sift/sift_roi $(BINDIR)
	cp $(BINDIR)/build_sift/matching $(BINDIR)

mgm:
	cd 3rdparty/mgm; $(MAKE)
	cp 3rdparty/mgm/mgm $(BINDIR)

sgbm:
	cd 3rdparty/sgbm; $(MAKE)
	cp 3rdparty/sgbm/sgbm $(BINDIR)
	cp 3rdparty/sgbm/call_sgbm.sh $(BINDIR)

sgbm_opencv:
	mkdir -p bin/sgbm_build
	cd bin/sgbm_build; cmake -D CMAKE_PREFIX_PATH=~/local ../../3rdparty/stereo_hirschmuller_2008; $(MAKE)
	cp bin/sgbm_build/sgbm2 $(BINDIR)
	cp bin/sgbm_build/SGBM $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM.sh $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_lap.sh $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_cauchy.sh $(BINDIR)

msmw:
	mkdir -p $(BINDIR)/msmw_build
	cd $(BINDIR)/msmw_build; cmake ../../3rdparty/msmw; $(MAKE)
	cp $(BINDIR)/msmw_build/libstereo/iip_stereo_correlation_multi_win2 $(BINDIR)

msmw2:
	mkdir -p $(BINDIR)/msmw2_build
	cd $(BINDIR)/msmw2_build; cmake ../../3rdparty/msmw2; $(MAKE)
	cp $(BINDIR)/msmw2_build/libstereo_newversion/iip_stereo_correlation_multi_win2_newversion $(BINDIR)

msmw3:
	mkdir -p $(BINDIR)/build_msmw3
	cd $(BINDIR)/build_msmw3; cmake ../../c/msmw; $(MAKE)
	cp $(BINDIR)/build_msmw3/msmw $(BINDIR)

tvl1:
	cd 3rdparty/tvl1flow_3; $(MAKE)
	cp 3rdparty/tvl1flow_3/tvl1flow $(BINDIR)
	cp 3rdparty/tvl1flow_3/callTVL1.sh $(BINDIR)

PROGRAMS = $(addprefix $(BINDIR)/,$(SRC)) plambda_without_fopenmp
SRC = $(SRCIIO) $(SRCFFT) $(SRCKKK)
SRCIIO = downsa backflow synflow imprintf iion qauto getminmax rescaleintensities qeasy crop morsi\
	morphoop cldmask disp_to_h_projective colormesh_projective tiffu
SRCFFT = gblur blur fftconvolve zoom_zeropadding zoom_2d
SRCKKK = watermask disp_to_h colormesh disp2ply bin2asc siftu ransac srtm4\
	srtm4_which_tile plyflatten

imscript: $(BINDIR) $(TIFDIR)/lib/libtiff.a $(PROGRAMS)

$(addprefix $(BINDIR)/,$(SRCIIO)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(C99) $(CFLAGS) $^ -o $@ $(IIOLIBS)

$(addprefix $(BINDIR)/,$(SRCFFT)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(C99) $(CFLAGS) $^ -o $@ $(IIOLIBS) $(FFTLIBS)

plambda_without_fopenmp:
	$(C99) -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN c/plambda.c c/iio.o -o bin/plambda $(IIOLIBS)

$(SRCDIR)/iio.o: $(SRCDIR)/iio.c $(SRCDIR)/iio.h
	$(C99) $(CFLAGS) -c -DIIO_ABORT_ON_ERROR -Wno-deprecated-declarations $< -o $@

$(SRCDIR)/rpc.o: c/rpc.c c/xfopen.c
	$(C99) $(CFLAGS) -c $< -o $@

$(SRCDIR)/mt/mt.o: c/mt/mt.c c/mt/mt.h
	$(C99) $(CFLAGS) -c $< -o $@

$(BINDIR)/bin2asc: c/bin2asc.c
	$(C99) $(CFLAGS) $^ -o $@

$(BINDIR)/siftu: c/siftu.c c/siftie.c
	$(C99) $(CFLAGS) $< -lm -o $@

$(BINDIR)/ransac: c/ransac.c c/fail.c c/xmalloc.c c/xfopen.c c/cmphomod.c\
	c/ransac_cases.c c/parsenumbers.c c/mt/mt.h c/mt/mt.o
	$(C99) $(CFLAGS) $< $(SRCDIR)/mt/mt.o -lm -o $@

$(BINDIR)/srtm4: c/srtm4.c $(SRCDIR)/Geoid.o $(SRCDIR)/geoid_height_wrapper.o
	$(C99) $(CFLAGS) -DMAIN_SRTM4 $^ $(IIOLIBS) $(LDLIBS) -o $@

$(BINDIR)/srtm4_which_tile: c/srtm4.c $(SRCDIR)/Geoid.o $(SRCDIR)/geoid_height_wrapper.o
	$(C99) $(CFLAGS) -DMAIN_SRTM4_WHICH_TILE $^ $(IIOLIBS) $(LDLIBS) -o $@

$(BINDIR)/watermask: $(SRCDIR)/iio.o $(SRCDIR)/Geoid.o\
	$(SRCDIR)/geoid_height_wrapper.o $(SRCDIR)/watermask.c $(SRCDIR)/fail.c\
	$(SRCDIR)/xmalloc.c $(SRCDIR)/pickopt.c $(SRCDIR)/rpc.c $(SRCDIR)/srtm4.c\
	$(SRCDIR)/iio.h $(SRCDIR)/parsenumbers.c
	$(C99) $(CFLAGS) $(SRCDIR)/iio.o $(SRCDIR)/Geoid.o $(SRCDIR)/geoid_height_wrapper.o $(SRCDIR)/watermask.c $(IIOLIBS) $(LDLIBS) -o $@

$(BINDIR)/disp_to_h: $(SRCDIR)/iio.o $(SRCDIR)/rpc.o c/disp_to_h.c c/vvector.h c/iio.h c/rpc.h c/read_matrix.c
	$(C99) $(CFLAGS) $(SRCDIR)/iio.o $(SRCDIR)/rpc.o c/disp_to_h.c $(IIOLIBS) -o $@

$(BINDIR)/colormesh: $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/DMS.o $(SRCDIR)/GeoCoords.o $(SRCDIR)/MGRS.o\
	$(SRCDIR)/PolarStereographic.o $(SRCDIR)/TransverseMercator.o $(SRCDIR)/UTMUPS.o c/colormesh.c c/iio.h\
	c/fail.c c/rpc.h c/read_matrix.c c/smapa.h c/timing.c c/timing.h
	$(C99) $(CFLAGS) $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/DMS.o $(SRCDIR)/GeoCoords.o $(SRCDIR)/MGRS.o $(SRCDIR)/PolarStereographic.o $(SRCDIR)/TransverseMercator.o $(SRCDIR)/UTMUPS.o c/colormesh.c c/timing.c $(IIOLIBS) $(LDLIBS) -o $@

$(BINDIR)/disp2ply: $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/DMS.o $(SRCDIR)/GeoCoords.o $(SRCDIR)/MGRS.o\
	$(SRCDIR)/PolarStereographic.o $(SRCDIR)/TransverseMercator.o $(SRCDIR)/UTMUPS.o c/disp2ply.c c/iio.h\
	c/fail.c c/rpc.h c/read_matrix.c c/smapa.h
	$(C99) $(CFLAGS) $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/DMS.o $(SRCDIR)/GeoCoords.o $(SRCDIR)/MGRS.o $(SRCDIR)/PolarStereographic.o $(SRCDIR)/TransverseMercator.o $(SRCDIR)/UTMUPS.o c/disp2ply.c $(IIOLIBS) $(LDLIBS) -o $@

$(BINDIR)/plyflatten: $(SRCDIR)/plyflatten.c $(SRCDIR)/iio.o
	$(C99) $(CFLAGS) $^ -o $@ $(IIOLIBS) $(GEOLIBS)


# GEOGRAPHICLIB STUFF
$(SRCDIR)/geographiclib_wrapper.o: c/geographiclib_wrapper.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/geoid_height_wrapper.o: c/geoid_height_wrapper.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@ -DGEOID_DATA_FILE_PATH="\"$(CURDIR)/data\""

$(SRCDIR)/DMS.o: c/DMS.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/GeoCoords.o: c/GeoCoords.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/MGRS.o: c/MGRS.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/PolarStereographic.o: c/PolarStereographic.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/TransverseMercator.o: c/TransverseMercator.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/UTMUPS.o: c/UTMUPS.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

$(SRCDIR)/Geoid.o: c/Geoid.cpp
	$(CXX) $(CPPFLAGS) -c $^ -I. -o $@

test:
	python s2p.py test.json

clean: clean_libtiff clean_geographiclib clean_monasse clean_homography\
	clean_sift clean_imscript clean_msmw clean_msmw2 clean_tvl1 clean_sgbm\
	clean_mgm

clean_libtiff:
	-rm $(TIFDIR)/lib/libtiff.a
	-rm $(TIFDIR)/bin/{raw2tiff,tiffinfo}
	-rm $(BINDIR)/{raw2tiff,tiffinfo}

clean_geographiclib:
	-rm $(BINDIR)/{Cart,Geo}Convert

clean_monasse:
	-rm -r $(BINDIR)/build_monasse_refactored
	-rm $(BINDIR)/{rectify_mindistortion}

clean_homography:
	-rm -r $(BINDIR)/build_homography
	-rm $(BINDIR)/homography

clean_sift:
	-rm -r $(BINDIR)/build_sift
	-rm $(BINDIR)/{sift_roi,match_cli}

clean_asift:
	-rm -r $(BINDIR)/build_asift
	-rm $(BINDIR)/demo_ASIFT

clean_imscript:
	-rm $(PROGRAMS)
	-rm $(SRCDIR)/iio.o
	-rm $(SRCDIR)/rpc.o
	#rm -r $(addsuffix .dSYM, $(PROGRAMS))

clean_msmw:
	-rm -r $(BINDIR)/msmw_build
	-rm $(BINDIR)/iip_stereo_correlation_multi_win2

clean_msmw2:
	-rm -r $(BINDIR)/msmw2_build
	-rm $(BINDIR)/iip_stereo_correlation_multi_win2_newversion

clean_tvl1:
	cd 3rdparty/tvl1flow_3; $(MAKE) clean
	-rm $(BINDIR)/{tvl1flow,callTVL1.sh}

clean_sgbm:
	cd 3rdparty/sgbm; $(MAKE) clean
	-rm $(BINDIR)/sgbm
	-rm $(BINDIR)/call_sgbm.sh

clean_mgm:
	cd 3rdparty/mgm; $(MAKE) clean
	-rm $(BINDIR)/mgm

.PHONY: default all geographiclib monasse sift sgbm sgbm_opencv msmw tvl1\
	imscript clean clean_libtiff clean_geographiclib clean_monasse clean_sift\
	clean_imscript clean_msmw clean_msmw2 clean_tvl1 clean_sgbm test
