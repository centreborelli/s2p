BINDIR = ../bin

all: imscript monasse

monasse: build_dir
	cd ../bin/monasse_refactored_build; cmake ../../c/monasse_refactored; make
	cp ../bin/monasse_refactored_build/bin/homography $(BINDIR)
	cp ../bin/monasse_refactored_build/bin/rectify_mindistortion $(BINDIR)
	cp ../bin/monasse_refactored_build/bin/sift_keypoints $(BINDIR)

build_dir:
	mkdir -p ../bin/monasse_refactored_build

imscript: bin_dir bin2asc downsa backflow synflow imprintf iion plambda qauto\
	qeasy crop disp_to_h colormesh_projective disp_to_h_projective\
	height_rpc_move veco morphoop gblur fftconvolve zoom_zeropadding zoom_2d\
	siftu srtm4 srtm4_which_tile ransac colormesh cldmask morsi watermask

bin_dir:
	mkdir -p $(BINDIR)

bin2asc: bin2asc.c
	cc -std=c99 bin2asc.c -o $(BINDIR)/bin2asc

morsi: iio.o morsi.o
	cc iio.o morsi.o -lpng -ltiff -ljpeg -o $(BINDIR)/morsi

watermask: iio.o watermask.o
	cc iio.o watermask.o -lpng -ltiff -ljpeg -o $(BINDIR)/watermask

cldmask: iio.o cldmask.o
	cc iio.o cldmask.o -lpng -ltiff -ljpeg -o $(BINDIR)/cldmask

colormesh: iio.o rpc.o colormesh.o geographiclib_wrapper.o
	cc iio.o rpc.o colormesh.o geographiclib_wrapper.o -lpng -ltiff -ljpeg -lGeographic -lstdc++ -o $(BINDIR)/colormesh -L../3rdparty/GeographicLib-1.32/src/

downsa: iio.o downsa.o
	cc iio.o downsa.o -lpng -ltiff -ljpeg -o $(BINDIR)/downsa

backflow: iio.o backflow.o
	cc iio.o backflow.o -lpng -ltiff -ljpeg -o $(BINDIR)/backflow

synflow: iio.o synflow.o
	cc iio.o synflow.o -lpng -ltiff -ljpeg -o $(BINDIR)/synflow

imprintf: iio.o imprintf.o
	cc iio.o imprintf.o -lpng -ltiff -ljpeg -o $(BINDIR)/imprintf

iion: iio.o iion.c
	cc iio.o iion.c -lpng -ltiff -ljpeg -o $(BINDIR)/iion

plambda: iio.o plambda.o
	cc iio.o plambda.o -lpng -ltiff -ljpeg -o $(BINDIR)/plambda

qauto: iio.o qauto.o
	cc iio.o qauto.o -lpng -ltiff -ljpeg -o $(BINDIR)/qauto

qeasy: iio.o qeasy.o
	cc iio.o qeasy.o -lpng -ltiff -ljpeg -o $(BINDIR)/qeasy

crop: iio.o crop.o
	cc iio.o crop.o -lpng -ltiff -ljpeg -o $(BINDIR)/crop

disp_to_h: iio.o rpc.o disp_to_h.o
	cc iio.o rpc.o disp_to_h.o -lpng -ltiff -ljpeg -o $(BINDIR)/disp_to_h

colormesh_projective: iio.o colormesh_projective.o
	cc iio.o colormesh_projective.o -lpng -ltiff -ljpeg -o $(BINDIR)/colormesh_projective

disp_to_h_projective: iio.o disp_to_h_projective.o
	cc iio.o disp_to_h_projective.o -lpng -ltiff -ljpeg -o $(BINDIR)/disp_to_h_projective

height_rpc_move: iio.o rpc.o height_rpc_move.o
	cc iio.o rpc.o height_rpc_move.o -lpng -ltiff -ljpeg -o $(BINDIR)/height_rpc_move

veco: iio.o veco.o
	cc iio.o veco.o -lpng -ltiff -ljpeg -o $(BINDIR)/veco

morphoop: iio.o morphoop.o
	cc iio.o morphoop.o -lpng -ltiff -ljpeg -o $(BINDIR)/morphoop

gblur: iio.o gblur.o
	cc iio.o gblur.o -lpng -ltiff -ljpeg -lfftw3f -o $(BINDIR)/gblur

fftconvolve: iio.o fftconvolve.o
	cc iio.o fftconvolve.o -lpng -ltiff -ljpeg -lfftw3f -o $(BINDIR)/fftconvolve

zoom_zeropadding: iio.o zoom_zeropadding.o
	cc iio.o zoom_zeropadding.o -lpng -ltiff -ljpeg -lfftw3f -o $(BINDIR)/zoom_zeropadding

zoom_2d: iio.o zoom_2d.o zoom_r2c.o
	cc iio.o zoom_2d.o zoom_r2c.o -lpng -ltiff -ljpeg -lfftw3 -o $(BINDIR)/zoom_2d

siftu: siftu.o
	cc siftu.o -lm -o $(BINDIR)/siftu

srtm4.o: srtm4.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DMAIN_SRTM4 srtm4.c

srtm4_which_tile.o: srtm4.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DMAIN_SRTM4_WHICH_TILE srtm4.c -o srtm4_which_tile.o

srtm4: srtm4.o iio.o
	cc srtm4.o iio.o -lm -lpng -ltiff -ljpeg -o $(BINDIR)/srtm4

srtm4_which_tile: srtm4_which_tile.o iio.o
	cc srtm4_which_tile.o iio.o -lm -lpng -ltiff -ljpeg -o $(BINDIR)/srtm4_which_tile

ransac: ransac.o
	cc ransac.o -lm -o $(BINDIR)/ransac

iio.o: iio.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN iio.c

rpc.o: rpc.c xfopen.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN rpc.c

colormesh.o: colormesh.c iio.h fail.c rpc.h read_matrix.c smapa.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN colormesh.c

geographiclib_wrapper.o: geographiclib_wrapper.cpp
	g++ -c -g -O3 geographiclib_wrapper.cpp

downsa.o: downsa.c fragments.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN downsa.c

backflow.o: backflow.c iio.h fail.c xmalloc.c getpixel.c smapa.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN backflow.c

synflow.o: synflow.c iio.h fragments.c synflow_core.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN synflow.c

imprintf.o: imprintf.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN imprintf.c

plambda.o: plambda.c smapa.h fail.c xmalloc.c random.c parsenumbers.c colorcoords.c getpixel.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN plambda.c

qauto.o: qauto.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN qauto.c

qeasy.o: qeasy.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN qeasy.c

crop.o: crop.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN crop.c

disp_to_h.o: disp_to_h.c vvector.h iio.h rpc.h read_matrix.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN disp_to_h.c

colormesh_projective.o: colormesh_projective.c iio.h fail.c read_matrix.c smapa.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN colormesh_projective.c

disp_to_h_projective.o: disp_to_h_projective.c vvector.h iio.h svd.c read_matrix.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN disp_to_h_projective.c

height_rpc_move.o: height_rpc_move.c vvector.h iio.h rpc.h read_matrix.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN height_rpc_move.c

veco.o: veco.c iio.h fail.c xmalloc.c random.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN veco.c

morphoop.o: morphoop.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN morphoop.c

gblur.o: gblur.c iio.h fail.c xmalloc.c vvector.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN gblur.c

fftconvolve.o: fftconvolve.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN fftconvolve.c

zoom_zeropadding.o: zoom_zeropadding.c iio.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN zoom_zeropadding.c

zoom_2d.o: zoom_2d.c iio.h zoom.h
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN zoom_2d.c

zoom_r2c.o: zoom_r2c.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN zoom_r2c.c

siftu.o: siftu.c siftie.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN siftu.c


ransac.o: ransac.c fail.c xmalloc.c xfopen.c cmphomod.c ransac_cases.c parsenumbers.c
	cc -c -std=c99 -g -O3 -DNDEBUG -DDONT_USE_TEST_MAIN ransac.c

cldmask.o: cldmask.c fail.c xmalloc.c xfopen.c parsenumbers.c drawsegment.c pickopt.c
	cc -c -std=c99 -g -O3 cldmask.c

watermask.o: watermask.c fail.c xmalloc.c pickopt.c rpc.c srtm4.c iio.h parsenumbers.c
	cc -c -std=c99 -g -O3 -DDONT_USE_TEST_MAIN watermask.c

morsi.o: morsi.c iio.h
	cc -c -std=c99 -g -O3 morsi.c

clean:
	rm *.o
