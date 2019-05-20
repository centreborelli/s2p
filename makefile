# the following two options are used to control all C and C++ compilations
export CFLAGS =   -march=native -O3
export CXXFLAGS = -march=native -O3

# these options are only used for the programs directly inside "./c/"
LDLIBS = -lstdc++
IIOLIBS = -lz -ltiff -lpng -ljpeg -lm
GEOLIBS = -lgeotiff -ltiff
GDAL_LIBS=`gdal-config --libs`
GDAL_CFLAGS=`gdal-config --cflags`

# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif

# names of source and destination directories
SRCDIR = c
BINDIR = bin
LIBDIR = lib

# default rule builds only the programs necessary for the test
default: $(BINDIR) $(LIBDIR) homography sift imscript mgm mgm_multi tvl1 lsd

# the "all" rule builds four further correlators
all: default msmw3 sgbm mgm_multi

# test for the default configuration
test: default
	pytest tests

# make sure that the destination directory is built
$(BINDIR):
	mkdir -p $(BINDIR)
$(LIBDIR):
	mkdir -p $(LIBDIR)

#
# four standard "modules": homography, sift, mgm, and mgm_multi
#

homography: $(BINDIR)
	$(MAKE) -j -C c/homography
	cp c/homography/homography $(BINDIR)

sift: $(BINDIR)
	$(MAKE) -j -C c/sift
	cp c/sift/libsift4ctypes.so $(LIBDIR)
mgm:
	$(MAKE) -C 3rdparty/mgm
	#cp 3rdparty/mgm/mgm $(BINDIR)

mgm_multi:
	mkdir -p $(BINDIR)/build_mgm_multi
	cd $(BINDIR)/build_mgm_multi; cmake ../../3rdparty/mgm_multi; $(MAKE)
	cp $(BINDIR)/build_mgm_multi/mgm_multi $(BINDIR)
	cp $(BINDIR)/build_mgm_multi/mgm $(BINDIR)

lsd:
	$(MAKE) -C 3rdparty/lsd
	cp 3rdparty/lsd/lsd $(BINDIR)

#
# rules for optional "modules": msmw, asift, sgbm, tvl1, etc
#

asift:
	mkdir -p $(BINDIR)/build_asift
	cd $(BINDIR)/build_asift; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/demo_ASIFT_src; $(MAKE)
	cp $(BINDIR)/build_asift/demo_ASIFT $(BINDIR)

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
	$(MAKE) -C 3rdparty/tvl1flow
	cp 3rdparty/tvl1flow/tvl1flow $(BINDIR)
	cp 3rdparty/tvl1flow/callTVL1.sh $(BINDIR)



#
# rules to build the programs under the source directory
#

PROGRAMS = $(addprefix $(BINDIR)/,$(SRC))
SRC = $(SRCIIO) $(SRCKKK)
SRCIIO = downsa backflow synflow imprintf qauto morsi\
	morphoop cldmask remove_small_cc plambda homwarp pview
SRCKKK = disp_to_h colormesh disp2ply multidisp2ply bin2asc plyflatten plyextrema

imscript: $(BINDIR) $(PROGRAMS)

$(addprefix $(BINDIR)/,$(SRCIIO)) : $(BINDIR)/% : $(SRCDIR)/%.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $^ -o $@ $(IIOLIBS)

$(SRCDIR)/iio.o: c/iio.c c/iio.h
	$(CC) $(CFLAGS) -c $< -o $@

$(SRCDIR)/rpc.o: c/rpc.c c/xfopen.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BINDIR)/bin2asc: c/bin2asc.c
	$(CC) $(CFLAGS) $^ -o $@

$(BINDIR)/disp_to_h: $(SRCDIR)/iio.o $(SRCDIR)/rpc.o c/disp_to_h.c c/vvector.h c/rpc.h c/read_matrix.c
	$(CC) $(CFLAGS) c/iio.o $(SRCDIR)/rpc.o c/disp_to_h.c $(IIOLIBS) -o $@

$(BINDIR)/colormesh: $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/colormesh.c c/fail.c c/rpc.h c/read_matrix.c c/smapa.h
	$(CC) $(CFLAGS) c/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/colormesh.c $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(SRCDIR)/triangulation.o: c/triangulation.c c/triangulation.h
	$(CC) $(CFLAGS) -c $< -lm -o $@

$(SRCDIR)/coordconvert.o: c/coordconvert.c c/coordconvert.h
	$(CC) $(CFLAGS) -c $< -lm -o $@

$(BINDIR)/multidisp2ply: 3rdparty/iio/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/triangulation.o $(SRCDIR)/coordconvert.o c/multidisp2ply.c c/vvector.h 3rdparty/iio/iio.h c/rpc.h c/triangulation.h c/coordconvert.h c/read_matrix.c
	$(CC) $(CFLAGS) 3rdparty/iio/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o $(SRCDIR)/triangulation.o $(SRCDIR)/coordconvert.o c/multidisp2ply.c $(IIOLIBS) -lGeographic -o $@

$(BINDIR)/disp2ply: $(SRCDIR)/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/disp2ply.c c/fail.c c/rpc.h c/read_matrix.c c/smapa.h
	$(CC) $(CFLAGS) c/iio.o $(SRCDIR)/rpc.o $(SRCDIR)/geographiclib_wrapper.o c/disp2ply.c $(IIOLIBS) $(LDLIBS) -lGeographic -o $@

$(BINDIR)/plyextrema: $(SRCDIR)/plyextrema.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS)  $^ -o $@ $(IIOLIBS)

$(BINDIR)/plyflatten: $(SRCDIR)/plyflatten.c $(SRCDIR)/iio.o
	$(CC) $(CFLAGS) $(GDAL_CFLAGS) $^ -o $@ $(IIOLIBS) $(GDAL_LIBS)

# Geographiclib wrappers
$(SRCDIR)/geographiclib_wrapper.o: c/geographiclib_wrapper.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

$(SRCDIR)/geoid_height_wrapper.o: c/geoid_height_wrapper.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@ -DGEOID_DATA_FILE_PATH="\"$(CURDIR)/c\""

# automatic dependency generation
-include makefile.dep
ALL_SOURCES=`ls c/*.c c/*.cc c/*.cpp`
.PHONY:
depend:
	$(CC) -MM $(ALL_SOURCES) | sed '/^[^ ]/s/^/c\//' > makefile.dep


# rules for cleaning, nothing interesting below this point
clean: clean_homography clean_asift clean_sift clean_imscript clean_msmw\
	clean_msmw2 clean_msmw3 clean_tvl1 clean_sgbm clean_mgm \
	clean_depend

clean_depend:
	$(RM) makefile.dep

clean_homography:
	$(MAKE) -C c/homography clean
	$(RM) $(BINDIR)/homography

clean_sift:
	$(MAKE) -C c/sift clean
	$(RM) $(LIBDIR)/libsift4ctypes.so

clean_asift:
	$(RM) -r $(BINDIR)/build_asift
	$(RM) $(BINDIR)/demo_ASIFT

clean_imscript:
	$(RM) $(PROGRAMS)
	$(RM) $(SRCDIR)/iio.o
	$(RM) $(SRCDIR)/rpc.o
	$(RM) $(SRCDIR)/geographiclib_wrapper.o
	$(RM) $(SRCDIR)/geoid_height_wrapper.o

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
	$(MAKE) -C 3rdparty/tvl1flow clean
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
