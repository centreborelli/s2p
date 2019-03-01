CXXFLAGS ?= -march=native -O3
override CXXFLAGS := $(CXXFLAGS) -fPIC
override CXXLIBS  := $(CXXLIBS)
override LDLIBS   := $(LDLIBS) -lm

OBJ = LibImages/LibImages.o \
      LibSSE/LibSSE.o \
      LibSift/KeyPoint.o \
      LibSift/LibSift.o \
      LibSift/Octave.o \
      LibSift/ScaleSpace.o \
      Utilities/Memory.o \
      Utilities/Parameters.o \
      Utilities/Time.o \
      Utilities/Utilities.o \
      sift4ctypes.o \

LIB = libsift4ctypes.so

all: $(LIB)

$(LIB): $(OBJ)
	$(CXX) -fPIC -shared -o $@ $^

clean:
	rm -f $(LIB) $(OBJ)
