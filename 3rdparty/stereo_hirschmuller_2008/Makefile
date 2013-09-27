all:
	gcc -std=c99 -g -o sgbm_simple stereo_match_simplified.cpp iio.c -lopencv_core -lopencv_highgui -lopencv_calib3d -lstdc++

test1:
	sgbm_simple testdata/venus/{im1,im2,d}.tif -30 0
	old_vflip testdata/venus/{d,im1,im2}.tif
test:
	sgbm_simple testdata/teddyF/{im1,im2,d}.tif -214 -1 
	old_vflip testdata/teddyF/{d,im1,im2}.tif
