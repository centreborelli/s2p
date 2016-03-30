# C source code

# C++ source code
CXXSRC	= main.cpp \
    ./LibImages/LibImages.cpp \
    ./LibMSMW/ConnectedComponents.cpp \
    ./LibMSMW/LibMSMW.cpp \
    ./LibMSMW/UtilitiesMSMW.cpp \
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
BIN	= msmw

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops -g -msse2 -mavx
#COPT	= -O3 -ftree-vectorize -funroll-loops -g -msse2

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra \
	-Wno-write-strings -Wno-deprecated -ansi
# link flags
LDFLAGS	= -lpng -ltiff -lfftw3f -lm -lrt -lpthread

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
