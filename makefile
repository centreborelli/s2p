# the following two options are used to control all C and C++ compilations
export CFLAGS   = -march=native -O3
export CXXFLAGS = -march=native -O3

# these options are only used for the programs directly inside "./c/"
IIOLIBS     = -lz -ltiff -lpng -ljpeg -lm


# default rule builds only the programs necessary for the test
default: homography sift mgm_multi tvl1 lsd executables libraries

# the "all" rule builds three further correlators
all: default msmw3 sgbm

# test for the default configuration
test: default
	env PYTHONPATH=. pytest tests

#
# four standard "modules": homography, sift, mgm, and mgm_multi
#

homography:
	$(MAKE) -j -C 3rdparty/homography
	cp 3rdparty/homography/homography bin

sift:
	$(MAKE) -j -C 3rdparty/sift/simd libsift4ctypes.so
	cp 3rdparty/sift/simd/libsift4ctypes.so lib


mgm_multi:
	$(MAKE) -C 3rdparty/mgm_multi
	cp 3rdparty/mgm_multi/mgm       bin
	cp 3rdparty/mgm_multi/mgm_multi bin

lsd:
	$(MAKE) -C 3rdparty/lsd
	cp 3rdparty/lsd/lsd bin

tvl1:
	$(MAKE) -C 3rdparty/tvl1flow
	cp 3rdparty/tvl1flow/tvl1flow bin
	cp 3rdparty/tvl1flow/callTVL1.sh bin

# compiled but not used (plain mgm is already included from multi_mgm)
mgm:
	$(MAKE) -C 3rdparty/mgm
	#cp 3rdparty/mgm/mgm bin


#
# rules for optional "modules": msmw, asift, sgbm, tvl1, etc
#

msmw3:
	make -C 3rdparty/msmw3
	cp 3rdparty/msmw3/msmw bin

asift:
	mkdir -p bin/build_asift
	cd bin/build_asift; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/demo_ASIFT_src; $(MAKE)
	cp bin/build_asift/demo_ASIFT bin

sgbm:
	$(MAKE) -C 3rdparty/sgbm
	cp 3rdparty/sgbm/sgbm bin

sgbm_opencv:
	mkdir -p bin/build_sgbm
	cd bin/build_sgbm; cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_PREFIX_PATH=~/local ../../3rdparty/stereo_hirschmuller_2008; $(MAKE)
	cp bin/build_sgbm/sgbm2 bin
	cp bin/build_sgbm/SGBM bin
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM.sh bin
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_lap.sh bin
	cp 3rdparty/stereo_hirschmuller_2008/callSGBM_cauchy.sh bin

msmw:
	mkdir -p bin/build_msmw
	cd bin/build_msmw; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/msmw; $(MAKE)
	cp bin/build_msmw/libstereo/iip_stereo_correlation_multi_win2 bin

msmw2:
	mkdir -p bin/build_msmw2
	cd bin/build_msmw2; cmake -D CMAKE_BUILD_TYPE=Release ../../3rdparty/msmw2; $(MAKE)
	cp bin/build_msmw2/libstereo_newversion/iip_stereo_correlation_multi_win2_newversion bin



#
# rules to build the programs under the source directory
#

SRCIIO   = downsa backflow qauto morsi cldmask remove_small_cc\
           plambda pview morphoop plyextrema colormesh
PROGRAMS = $(addprefix bin/,$(SRCIIO))

executables: $(PROGRAMS)


# generic rule for building binary objects from C sources
c/%.o : c/%.c
	$(CC) -fpic $(CFLAGS) -c $< -o $@

# generic rule for building binary objects from C++ sources
c/%.o: c/%.cpp
	$(CXX) -fpic $(CXXFLAGS) -c $^ -o $@

# generic rule to build most imscript binaries
bin/% : c/%.o c/iio.o
	$(CC) $^ -o $@ $(IIOLIBS)

# this particular object requires a hardcoded filename
c/geoid_height_wrapper.o: c/geoid_height_wrapper.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@ -DGEOID_DATA_FILE_PATH="\"$(CURDIR)/c\""

# this particular program combines different objects in a non-standard way
bin/colormesh: c/colormesh.o c/iio.o c/rpc.o c/geographiclib_wrapper.o
	$(CC) $^ $(IIOLIBS) -lstdc++ -lGeographic -o $@



#
# rules to build the dynamic objects that are used via ctypes
#

libraries: lib/libplyflatten.so lib/disp_to_h.so

lib/disp_to_h.so: c/disp_to_h.o c/geographiclib_wrapper.o c/iio.o c/rpc.o
	$(CC) -shared $^ $(IIOLIBS) -lGeographic -o $@

lib/libplyflatten.so: c/plyflatten.o
	$(CC) -shared $^ -o $@




# automatic dependency generation
-include .deps.mk
.PHONY:
depend:
	$(CC) -MM `ls c/*.c c/*.cpp` | sed '/^[^ ]/s/^/c\//' > .deps.mk


# rules for cleaning, nothing interesting below this point
clean: clean_homography clean_asift clean_sift clean_imscript clean_msmw\
       clean_msmw2 clean_msmw3 clean_tvl1 clean_sgbm clean_mgm clean_mgm_multi\
       clean_lsd clean_s2p
	$(RM) c/*.o bin/* lib/*
	$(RM) -r s2p_tmp

distclean: clean ; $(RM) .deps.mk


# clean targets that use recursive makefiles
clean_homography: ; $(MAKE) clean -C 3rdparty/homography
clean_sift:       ; $(MAKE) clean -C 3rdparty/sift/simd
clean_tvl1:       ; $(MAKE) clean -C 3rdparty/tvl1flow
clean_sgbm:       ; $(MAKE) clean -C 3rdparty/sgbm
clean_mgm:        ; $(MAKE) clean -C 3rdparty/mgm
clean_mgm_multi:  ; $(MAKE) clean -C 3rdparty/mgm_multi
clean_lsd:        ; $(MAKE) clean -C 3rdparty/lsd
clean_msmw3:      ; $(MAKE) clean -C 3rdparty/msmw3

# clean targets that use a build dir
clean_asift:      ; $(RM) -r bin/build_asift
clean_msmw:       ; $(RM) -r bin/build_msmw
clean_msmw2:      ; $(RM) -r bin/build_msmw2


.PHONY: default all sift sgbm sgbm_opencv msmw tvl1 imscript clean clean_sift\
	clean_imscript clean_msmw clean_msmw2 clean_tvl1 clean_sgbm clean_mgm\
	clean_mgm_multi clean_lsd clean_s2p test distclean


# The following conditional statement appends "-std=gnu99" to CFLAGS when the
# compiler does not define __STDC_VERSION__.  The idea is that many older
# compilers are able to compile standard C when given that option.
# This hack seems to work for all versions of gcc, clang and icc.
CVERSION = $(shell $(CC) -dM -E - < /dev/null | grep __STDC_VERSION__)
ifeq ($(CVERSION),)
CFLAGS := $(CFLAGS) -std=gnu99
endif
