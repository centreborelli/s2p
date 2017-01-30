CC = gcc
CXX = c++
CFLAGS = -std=c99 -O3
CPPFLAGS = -O3
CV_CORE_OBJS = alloc.o arithm.o array.o convert.o copy.o datastructs.o gpumat.o lapack.o mathfuncs.o matmul.o matop.o matrix.o opengl_interop.o opengl_interop_deprecated.o parallel.o persistence.o stat.o system.o tables_core.o
CV_IMGPROC_OBJS = smooth.o tables_imgproc.o

all: sgbm

sgbm: sgbm.cpp stereosgbm.cpp libiio.a libcv_core.a libcv_imgproc.a
		$(CXX) -o sgbm sgbm.cpp stereosgbm.cpp -L. -liio -lcv_core -lcv_imgproc -lpng -ltiff -ljpeg -lz
libiio.a: iio.o
		ar rcs libiio.a iio.o
libcv_core.a: $(CV_CORE_OBJS)
		ar rcs libcv_core.a $(CV_CORE_OBJS)
libcv_imgproc.a: $(CV_IMGPROC_OBJS)
		ar rcs libcv_imgproc.a $(CV_IMGPROC_OBJS)
.cpp.o:
		$(CXX) $(CPPFLAGS) -c $< -I.
iio.o:
		$(CC) $(CFLAGS) -c iio.c -o iio.o -Wno-deprecated-declarations
clean:
		$(RM) *.o
		$(RM) *.a
		$(RM) sgbm
