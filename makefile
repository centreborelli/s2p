C99 = $(CC) -std=gnu99
CFLAGS = -g -O3 -DNDEBUG
CPPFLAGS = -g -O3 -DNDEBUG
LDLIBS = -lstdc++
IIOLIBS = -lz -ltiff -lpng -ljpeg -lm
GEOLIBS = -lgeotiff -ltiff
FFTLIBS = -lfftw3f -lfftw3

ifeq ($(CC), gcc)
	CFLAGS += -fopenmp
endif

ifeq ($(CXX), g++)
	CPPFLAGS += -fopenmp
endif

#OS := $(shell uname -s)
#ifeq ($(OS), Linux)
#	LDLIBS += -lrt
#endif

BINDIR = bin
SRCDIR = c

default: $(BINDIR) homography sift imscript mgm piio

all: default msmw3 sgbm tvl1


$(BINDIR):
	mkdir -p $(BINDIR)

piio: s2plib/piio/libiio.so

s2plib/piio/libiio.so: s2plib/piio/setup.py s2plib/piio/freemem.c 3rdparty/iio/iio.c 3rdparty/iio/iio.h
	cd s2plib/piio; python setup.py build

asift:
	mkdir -p $(BINDIR)/build_asift
	cd $(BINDIR)/build_asift; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/demo_ASIFT_src; $(MAKE)
	cp $(BINDIR)/build_asift/demo_ASIFT $(BINDIR)

homography: $(BINDIR)
	$(MAKE) -j -C c/homography
	cp c/homography/homography $(BINDIR)

sift: $(BINDIR)
	$(MAKE) -j -C c/sift
	cp c/sift/sift_roi $(BINDIR)
	cp c/sift/matching $(BINDIR)

mgm:
	$(MAKE) -C 3rdparty/mgm
	cp 3rdparty/mgm/mgm $(BINDIR)

sgbm:
	$(MAKE) -C 3rdparty/sgbm
	cp 3rdparty/sgbm/sgbm $(BINDIR)

sgbm_opencv:
	mkdir -p bin/build_sgbm
	cd bin/build_sgbm; cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_PREFIX_PATH=~/local ../../3rdparty/stereo_hirschmuller_2008; $(MAKE)
	cp bin/build_sgbm/sgbm2 $(BINDIR)
	cp bin/build_sgbm/SGBM $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM.sh $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_lap.sh $(BINDIR)
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_cauchy.sh $(BINDIR)

msmw:
	mkdir -p $(BINDIR)/build_msmw
	cd $(BINDIR)/build_msmw; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/msmw; $(MAKE)
	cp $(BINDIR)/build_msmw/libstereo/iip_stereo_correlation_multi_win2 $(BINDIR)

msmw2:
	mkdir -p $(BINDIR)/build_msmw2
	cd $(BINDIR)/build_msmw2; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/msmw2; $(MAKE)
	cp $(BINDIR)/build_msmw2/libstereo_newversion/iip_stereo_correlation_multi_win2_newversion $(BINDIR)

msmw3:
	mkdir -p $(BINDIR)/build_msmw3
	cd $(BINDIR)/build_msmw3; cmake -D CMAKE_BUILD_TYPE=Release ../../c/msmw; $(MAKE)
	cp $(BINDIR)/build_msmw3/msmw $(BINDIR)

tvl1:
	$(MAKE) -C 3rdparty/tvl1flow_3
	cp 3rdparty/tvl1flow_3/tvl1flow $(BINDIR)
	cp 3rdparty/tvl1flow_3/callTVL1.sh $(BINDIR)

PROGRAMS = $(addprefix $(BINDIR)/,$(SRC))
SRC = $(SRCIIO) $(SRCFFT) $(SRCKKK)
SRCIIO = downsa backflow synflow imprintf iion qauto qeasy crop bdint morsi\
	morphoop cldmask disp_to_h_projective colormesh_projective\
	remove_small_cc plambda
SRCFFT = gblur blur fftconvolve zoom_zeropadding zoom_2d
SRCKKK = watermask disp_to_h colormesh disp2ply bin2asc siftu ransac srtm4\
	srtm4_which_tile plyflatten plyextrema plytodsm

imscript: $(BINDIR) $(PROGRAMS)

$(addprefix $(BINDIR)/,$(SRCIIO)) : $(BINDIR)/% : $(SRCDIR)/%.c 3rdparty/iio/iio.o
	$(C99) $(CFLAGS) $^ -o $@ $(IIOLIBS)

$(addprefix $(BINDIR)/,$(SRCFFT)) : $(BINDIR)/% : $(SRCDIR)/%.c 3rdparty/iio/iio.o
	$(C99) $(CFLAGS) $^ -o $@ $(IIOLIBS) $(FFTLIBS)

3rdparty/iio/iio.o: 3rdparty/iio/iio.c 3rdparty/iio/iio.h
	$(C99) $(CFLAGS) -c $< -o $@

$(SRCDIR)/rpc.o: c/rpc.c c/xfopen.c
	$(C99) $(CFLAGS) -c $< -o $@

$(BINDIR)/bin2asc: c/bin2asc.c
	$(C99) $(CFLAGS) $^ -o $@

$(BINDIR)/siftu: c/siftu.c c/siftie.c
	$(C99) $(CFLAGS) $< -lm -o $@

$(BINDIR)/ransac: c/ransac.c c/fail.c c/xmalloc.c c/xfopen.c c/homographies.c\
	c/ransac_cases.c c/parsenumbers.c c/random.c
	$(C99) $(CFLAGS) $< -lm -o $@

$(BINDIR)/srtm4: c/srtm4.c $(SRCDIR)/geoid_height_wrapper.o
	$(C99) $(CFLAGS) -DMAIN_SRTM4 $^ $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(BINDIR)/srtm4_which_tile: c/srtm4.c $(SRCDIR)/geoid_height_wrapper.o
	$(C99) $(CFLAGS) -DMAIN_SRTM4_WHICH_TILE $^ $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(BINDIR)/watermask: 3rdparty/iio/iio.o $(SRCDIR)/geoid_height_wrapper.o $(SRCDIR)/watermask.c $(SRCDIR)/fail.c\
	$(SRCDIR)/xmalloc.c $(SRCDIR)/pickopt.c $(SRCDIR)/rpc.c $(SRCDIR)/srtm4.c 3rdparty/iio/iio.h $(SRCDIR)/parsenumbers.c
	$(C99) $(CFLAGS) 3rdparty/iio/iio.o $(SRCDIR)/geoid_height_wrapper.o $(SRCDIR)/watermask.c $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(BINDIR)/disp_to_h: 3rdparty/iio/iio.o $(SRCDIR)/rpc.o c/disp_to_h.c c/vvector.h 3rdparty/iio/iio.h c/rpc.h c/read_matrix.c
	$(C99) $(CFLAGS) 3rdparty/iio/iio.o $(SRCDIR)/rpc.o c/disp_to_h.c $(IIOLIBS) -o $@

$(BINDIR)/colormesh: 3rdparty/iio/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/colormesh.c 3rdparty/iio/iio.h c/fail.c c/rpc.h c/read_matrix.c c/smapa.h c/timing.c c/timing.h
	$(C99) $(CFLAGS) 3rdparty/iio/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/colormesh.c c/timing.c $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(BINDIR)/disp2ply: 3rdparty/iio/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/disp2ply.c 3rdparty/iio/iio.h c/fail.c c/rpc.h c/read_matrix.c c/smapa.h
	$(C99) $(CFLAGS) 3rdparty/iio/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/disp2ply.c $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(BINDIR)/plyextrema: $(SRCDIR)/plyextrema.c 3rdparty/iio/iio.o
	$(C99) $(CFLAGS)  $^ -o $@ $(IIOLIBS)

$(BINDIR)/plyflatten: $(SRCDIR)/plyflatten.c 3rdparty/iio/iio.o
	$(C99) $(CFLAGS) -I/usr/include/geotiff $^ -o $@ $(IIOLIBS) $(GEOLIBS)

$(BINDIR)/plytodsm: $(SRCDIR)/plytodsm.c 3rdparty/iio/iio.o
	$(C99) $(CFLAGS) -I/usr/include/geotiff $^ -o $@ $(IIOLIBS) $(GEOLIBS)

# Geographiclib wrappers
$(SRCDIR)/geographiclib_wrapper.o: c/geographiclib_wrapper.cpp
	$(CXX) $(CPPFLAGS) -c $^ -o $@

$(SRCDIR)/geoid_height_wrapper.o: c/geoid_height_wrapper.cpp
	$(CXX) $(CPPFLAGS) -c $^ -o $@ -DGEOID_DATA_FILE_PATH="\"$(CURDIR)/c\""

test:
	python -u s2p_test.py

clean: clean_homography clean_asift clean_sift clean_imscript clean_msmw\
	clean_msmw2 clean_msmw3 clean_tvl1 clean_sgbm clean_mgm

clean_homography:
	$(MAKE) -C c/homography clean
	$(RM) $(BINDIR)/homography

clean_sift:
	$(MAKE) -C c/sift clean
	$(RM) $(BINDIR)/sift_roi
	$(RM) $(BINDIR)/matching

clean_asift:
	$(RM) -r $(BINDIR)/build_asift
	$(RM) $(BINDIR)/demo_ASIFT

clean_imscript:
	$(RM) $(PROGRAMS)
	$(RM) 3rdparty/iio/iio.o
	$(RM) $(SRCDIR)/rpc.o

clean_msmw:
	$(RM) -r $(BINDIR)/build_msmw
	$(RM) $(BINDIR)/iip_stereo_correlation_multi_win2

clean_msmw2:
	$(RM) -r $(BINDIR)/build_msmw2
	$(RM) $(BINDIR)/iip_stereo_correlation_multi_win2_newversion

clean_msmw3:
	$(RM) -r $(BINDIR)/build_msmw3
	$(RM) $(BINDIR)/msmw

clean_tvl1:
	$(MAKE) -C 3rdparty/tvl1flow_3 clean
	$(RM) $(BINDIR)/tvl1flow
	$(RM) $(BINDIR)/callTVL1.sh

clean_sgbm:
	$(MAKE) -C 3rdparty/sgbm clean
	$(RM) $(BINDIR)/sgbm

clean_mgm:
	$(MAKE) -C 3rdparty/mgm clean
	$(RM) $(BINDIR)/mgm

.PHONY: default all sift sgbm sgbm_opencv msmw tvl1 imscript clean clean_sift\
	clean_imscript clean_msmw clean_msmw2 clean_tvl1 clean_sgbm test
