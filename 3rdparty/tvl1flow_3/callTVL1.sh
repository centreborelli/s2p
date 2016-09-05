#!/bin/bash

if [ "$4" == "" ]; then
   echo "Usage:"
   echo "  $0 im1 im2 out_disp out_mask"
   echo ""
   echo ""
   echo "  Wrapper to the TVL1 function with Left-Right consistency"
   echo "  Javier Sánchez Pérez, Enric Meinhardt-Llopis, and Gabriele Facciolo,"
   echo "  TV-L1 Optical Flow Estimation, Image Processing On Line, vol. 2013, pp. 137–150."
   echo "  http://dx.doi.org/10.5201/ipol.2013.26"
   exit 1
fi

a=$1
b=$2
disp=$3
mask=$4

# ! make it invariant to illumination changes
gblur 2 $a | plambda  - " x 4 * x(-1,0) -1 * + x(1,0) -1 * + x(0,-1) -1 * + x(0,1) -1 * + " -o $a
gblur 2 $b | plambda  - " x 4 * x(-1,0) -1 * + x(1,0) -1 * + x(0,-1) -1 * + x(0,1) -1 * + " -o $b

# convert input images to range 0:255
thresholds=`qauto $a $a 2>&1 |  cut -f2 -d=`
# use the same threshods for the second image
qeasy $thresholds $b $b

# call
dl=$a.dl.tif
dr=$a.dr.tif
dlr=$a.dlr.tif
tvl1flow $a $b $dl 0 0.3 0.15 0.3 9 0.5 3 0.01 1
tvl1flow $b $a $dr 0 0.3 0.15 0.3 9 0.5 3 0.01 1
backflow $dl $dr $dlr
plambda $dl $dlr "x y + split hypot 1 < 255 *" -o $mask
plambda $dl $mask "y 0 > x[0] nan if" -o $disp

rm $dl $dr $dlr
