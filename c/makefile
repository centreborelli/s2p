all:
	mkdir -p ../bin
	cc -o ../bin/imprintf -std=c99 -g -O3 -DNDEBUG iio.c imprintf.c -lpng -ltiff -ljpeg
	cc -o ../bin/plambda -std=c99 -g -O3 -DNDEBUG iio.c plambda.c -lpng -ltiff -ljpeg
	cc -o ../bin/iion -std=c99 -g -O3 -DNDEBUG iio.c iion.c -lpng -ltiff -ljpeg
	cc -o ../bin/qauto -std=c99 -g -O3 -DNDEBUG iio.c qauto.c -lpng -ltiff -ljpeg
	cc -o ../bin/qeasy -std=c99 -g -O3 -DNDEBUG iio.c qeasy.c -lpng -ltiff -ljpeg
	cc -o ../bin/gblur -std=c99 -g -O3 -DNDEBUG iio.c gblur.c -lpng -ltiff -ljpeg -lfftw3f
	cc -o ../bin/crop -std=c99 -g -O3 -DNDEBUG iio.c crop.c -lpng -ltiff -ljpeg
	cc -o ../bin/fftconvolve -std=c99 -g -O3 -DNDEBUG iio.c fftconvolve.c -lpng -ltiff -ljpeg -lfftw3f
	cc -o ../bin/zoom_zeropadding -std=c99 -g -O3 -DNDEBUG iio.c zoom_zeropadding.c -lpng -ltiff -ljpeg -lfftw3f
	cc -o ../bin/zoom_2d -std=c99 -g -O3 -DNDEBUG iio.c zoom_2d.c zoom_r2c.c -lpng -ltiff -ljpeg -lfftw3
	cc -o ../bin/siftu -std=c99 -g -O3 -DNDEBUG siftu.c -lm
	cc -o ../bin/srtm4 -std=c99 -g -O3 -DNDEBUG -DMAIN_SRTM4 srtm4.c -lm
	cc -o ../bin/disp_to_h -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN iio.c rpc.c disp_to_h.c -lpng -ltiff -ljpeg
	cc -o ../bin/colormesh -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN iio.c rpc.c colormesh.c -lpng -ltiff -ljpeg
	mkdir -p ../bin/monasse_refactored_build
	cd ../bin/monasse_refactored_build; cmake ../../c/monasse_refactored; make
	cp ../bin/monasse_refactored_build/bin/homography ../bin
	cp ../bin/monasse_refactored_build/bin/rectify_mindistortion ../bin
	cp ../bin/monasse_refactored_build/bin/sift_keypoints ../bin
