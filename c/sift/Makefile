CFLAGS = -g -O3 -Wall -Wextra
SIFT_DIR = sift_anatomy_20141201

SIFT_OBJa = lib_sift.o lib_sift_anatomy.o lib_scalespace.o lib_keypoint.o lib_description.o lib_discrete.o lib_util.o
SIFT_OBJb = lib_util.o lib_keypoint.o lib_matching.o
OBJa = $(addprefix $(SIFT_DIR)/,$(SIFT_OBJa)) fancy_image.o ../iio.o
OBJb = $(addprefix $(SIFT_DIR)/,$(SIFT_OBJb))

default: sift_roi match_cli

sift_roi: main.c  $(OBJa)
	$(CC) $^ -o $@ $(CFLAGS) -lpng -ltiff -ljpeg

match_cli: $(SIFT_DIR)/match_cli.c $(OBJb)
	$(CC) $^ -o $@ $(CFLAGS)

../iio.o: ../iio.c ../iio.h
	$(CC) -c -o $@ $< $(CFLAGS) -Wno-unused-function -Wno-deprecated-declarations

%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

test: sift_roi
	./sift_roi test_data/img_tiled.tif 100 120 200 178

clean:
	-rm sift_roi match_cli $(OBJa) $(OBJb)
