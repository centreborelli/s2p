mgm: ../iio/iio.o mgm.cc img.cc point.cc
	$(CXX) -O3 -DNDEBUG -ffast-math -I../iio -Icommon -DTEST_MAIN $^ -o mgm -lpng -ltiff -ljpeg

../iio/iio.o: ../iio/iio.c ../iio/iio.h
	$(CC) -std=c99 -O3 -g -fopenmp -DNDEBUG -DDONT_USE_TEST_MAIN -DIIO_ABORT_ON_ERROR -Wno-deprecated-declarations -c $< -o $@

clean:
	$(RM) -f mgm
