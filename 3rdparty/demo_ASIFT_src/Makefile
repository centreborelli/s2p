# C source code
CSRC	= io_png/io_png.c
# C++ source code
CXXSRC	= numerics1.cpp frot.cpp splines.cpp fproj.cpp \
	library.cpp flimage.cpp filter.cpp \
	demo_lib_sift.cpp compute_asift_keypoints.cpp \
	compute_asift_matches.cpp \
	libMatch/match.cpp libNumerics/numerics.cpp \
	orsa.cpp demo_ASIFT.cpp
# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= demo_ASIFT

default	: $(BIN)

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -ansi
# link flags
LDFLAGS	= -lpng -lm

# use local embedded libraries with `make LOCAL_LIBS=1`
ifdef LOCAL_LIBS
# library location
LIBDIR = io_png/libs/build/lib
INCDIR = io_png/libs/build/include
# libpng is required
LIBDEPS += libpng
# compile options to use the local libpng header
CFLAGS 	+= -I$(INCDIR) -D_LOCAL_LIBS
# link options to use the local libraries
LDFLAGS = $(LIBDIR)/libpng.a $(LIBDIR)/libz.a -lm
# io_png.o needs png.h
io_png/io_png.o	:  $(LIBDEPS)
endif

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas
endif

# build the local png library
.PHONY	: libpng
libpng	:
	$(MAKE) -C io_png/libs libpng

CFLAGS 	+= -I./io_png/
CXXFLAGS 	+= -I./

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the opject code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

# housekeeping
.PHONY	: clean distclean
clean	:
	$(RM) $(OBJ)
	$(MAKE) -C ./io_png/libs $@
distclean	: clean
	$(RM) $(BIN)
	$(MAKE) -C ./io_png/libs $@
