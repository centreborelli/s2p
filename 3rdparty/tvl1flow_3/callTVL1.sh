#!/bin/bash

# get real path of the binary file
rel_path_script=`dirname $0`

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

# just make a working copy of the input images
aa=$(basename "$a")
bb=$(basename "$b")
iion $a /tmp/$aa.tif
iion $b /tmp/$bb.tif
a=/tmp/$aa.tif
b=/tmp/$bb.tif


# ! make it invariant to illumination changes
gblur 2 $a   | plambda  - " x 4 *    x(-1,0)  -1 *  +  x(1,0)  -1 * +  x(0,-1)  -1 * +    x(0,1)  -1 * + " | iion - $a
gblur 2 $b   | plambda  - " x 4 *    x(-1,0)  -1 *  +  x(1,0)  -1 * +  x(0,-1)  -1 * +    x(0,1)  -1 * + " | iion - $b


# convert input images to range 0:255
thresholds=`qauto $a $a 2>&1 |  cut -f2 -d=`
# use the same threshods for the second image
qeasy $thresholds $b $b


# call
dl=$a.dl.tif
dr=$a.dr.tif
dlr=$a.dlr.tif
$rel_path_script/tvl1flow $a $b $dl 0 0.3 0.15 0.3 9 0.5 3 0.01 1
$rel_path_script/tvl1flow $b $a $dr 0 0.3 0.15 0.3 9 0.5 3 0.01 1
$rel_path_script/backflow $dl $dr $dlr
plambda $dl $dlr " x y + split hypot 1 < 255 *" | iion - $mask
plambda $dl $mask "y 0 > x[0] nan if" | iion - $disp

rm $dl $dr $a $b $dlr
