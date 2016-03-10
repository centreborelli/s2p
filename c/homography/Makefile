# C source code

# C++ source code
CXX = g++
CXXSRC	= main.cpp \
    ./LibHomography/Homography.cpp \
    ./LibHomography/Splines.cpp \
    ./LibImages/LibImages.cpp \
    ./Utilities/Memory.cpp \
    ./Utilities/Parameters.cpp \
    ./Utilities/Time.cpp \
    ./Utilities/Utilities.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= homography

# C optimization flags
#COPT	= -O3 -ftree-vectorize -funroll-loops -g -msse2 -mavx
COPT	= -O3 -ftree-vectorize -funroll-loops -g -msse2

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -Wno-unused-parameter -ansi
# link flags
LDFLAGS	= -lpng -ltiff -lm -lpthread

# OS detection hacks
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	LDFLAGS += -lrt
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

# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the object code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

test: homography
	./homography -i test_data/input_uint16.tif -t test_data/identity.txt -o test_data/output_id.tif
	./homography -i test_data/input_uint16.tif -t test_data/homography.txt -o test_data/output_hom.tif
	./homography -i test_data/input_uint16_tiled.tif -t test_data/homography.txt -o test_data/output_hom.tif

clean:
	-rm homography
	-rm *.o
	-rm */*.o
