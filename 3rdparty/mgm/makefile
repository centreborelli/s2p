CC = gcc
CXX = g++

default:
	$(CC) -std=c99 -O3 -DNDEBUG -ffast-math -fopenmp -DIIO_ABORT_ON_ERROR -Wno-deprecated-declarations -c iio/iio.c -o iio.o
	$(CXX) -O3 -DNDEBUG -ffast-math -fopenmp -Iiio -Icommon -DTEST_MAIN mgm.cc img.cc point.cc iio.o -o mgm -lpng -ltiff -ljpeg

clean:
	rm -f iio.o mgm
