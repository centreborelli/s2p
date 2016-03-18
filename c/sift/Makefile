CFLAGS = -g -O3 -Wall -Wextra
BIN = sift_roi
SIFT_DIR = sift_anatomy_20141201
SIFT_OBJ = lib_sift.o lib_sift_anatomy.o lib_scalespace.o lib_keypoint.o lib_description.o lib_discrete.o lib_util.o
OBJ = $(addprefix $(SIFT_DIR)/,$(SIFT_OBJ)) fancy_image.o ../iio.o

$(BIN): main.c  $(OBJ)
	gcc $^ -o $@ -lpng -ltiff -ljpeg

../iio.o: ../iio.c ../iio.h
	$(CC) -c -o $@ $< $(CFLAGS) -Wno-unused-function -Wno-deprecated-declarations

%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

test: sift_roi
	./$(BIN) test_data/img_tiled.tif 100 120 200 178

clean:
	-rm $(BIN) $(OBJ)
