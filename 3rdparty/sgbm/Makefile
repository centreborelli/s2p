CC = gcc
CXX = g++
CFLAGS = -std=c99 -O3
CPPFLAGS = -fpic -O3
CV_CORE_OBJS = alloc.o arithm.o array.o convert.o copy.o datastructs.o gpumat.o lapack.o mathfuncs.o matmul.o matop.o matrix.o opengl_interop.o opengl_interop_deprecated.o parallel.o persistence.o stat.o system.o tables_core.o
CV_IMGPROC_OBJS = filter.o smooth.o tables_imgproc.o utils.o

all: sgbm

sgbm: sgbm.cpp stereosgbm.cpp libiio.so libcv_core.so libcv_imgproc.so
		$(CXX) -o sgbm sgbm.cpp stereosgbm.cpp -liio -lcv_core -lcv_imgproc -L.
libiio.so: iio.c include/iio.h
		$(CC) $(CFLAGS) -o libiio.so -shared -fpic iio.c -lpng -ltiff -ljpeg
libcv_core.so: $(CV_CORE_OBJS)
		$(CXX) -o libcv_core.so -shared -fpic $(CV_CORE_OBJS) -lz
libcv_imgproc.so: $(CV_IMGPROC_OBJS) libcv_core.so
		$(CXX) -o libcv_imgproc.so -shared -fpic $(CV_IMGPROC_OBJS) -lcv_core -L.
.cpp.o:
		$(CXX) $(CPPFLAGS) -c $< -I.
clean:
		rm *.o
		rm *.so
		rm sgbm
