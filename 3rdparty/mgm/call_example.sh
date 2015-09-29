#!/bin/bash

OMP_NUM_THREADS=4 MEDIAN=1 CENSUS_NCC_WIN=5 TSGM=3 ./mgm -r -22 -R 19 -s vfit -t census -O 8 test_data/rectified_{ref,sec}.tif disp.tif
