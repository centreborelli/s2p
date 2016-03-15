# source code
CXX = g++
CXXSRC = main.cpp \
    LibHomography/Homography.cpp \
    LibHomography/Splines.cpp \
    LibImages/LibImages.cpp \
    Utilities/Memory.cpp \
    Utilities/Parameters.cpp \
    Utilities/Time.cpp \
    Utilities/Utilities.cpp

# objects
COBJ = $(CSRC:.c=.o)
CXXOBJ = $(CXXSRC:.cpp=.o)
OBJ	= $(COBJ) $(CXXOBJ)

# binary target
BIN	= homography

# optimization flags
COPT = -O3 -ftree-vectorize -funroll-loops -g -msse2 #-mavx
CXXOPT = $(COPT)

# compiler and linker flags
CFLAGS = $(COPT) -Wall -Wextra -Wno-write-strings -ansi
CXXFLAGS = $(CXXOPT) -Wall -Wextra -Wno-write-strings -Wno-deprecated -Wno-unused-parameter -ansi
LDFLAGS	= -lpng -ltiff -lm -lpthread

# OS detection hacks
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	LDFLAGS += -lrt
endif

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS += -fopenmp
CXXFLAGS += -fopenmp
LDFLAGS += -lgomp
else
CFLAGS += -Wno-unknown-pragmas
CXXFLAGS += -Wno-unknown-pragmas
endif

# partial compilation
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

# link all the object code
$(BIN): $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)

test: homography
	./homography -i test_data/input_uint16.tif -t test_data/identity.txt -o test_data/output_id.tif
	./homography -i test_data/input_uint16.tif -t test_data/homography.txt -o test_data/output_hom.tif
	./homography -i test_data/input_uint16_tiled.tif -t test_data/homography.txt -o test_data/output_hom.tif

clean:
	-rm homography
	-rm *.o
	-rm */*.o
