#!/bin/bash

# get real path of the jar file
rel_path_script=`dirname $0`

if [ "$3" == "" ]; then
   echo "Usage:"
   $rel_path_script/SGBM
   echo ""
   echo "   Implements: H. Hirschmuller, P.R. Innocent, and J.Garibaldi. "
   echo "   \"Real-Time Correlation-Based Stereo Vision with Reduced Border Errors.\" "
   echo "   Int. J. Comput. Vision 47, 1-3 2002"
   exit 1
fi

a=$1
b=$2


# convert input images to png
aa=$(basename "$a")
a_extension="${aa##*.}"
a_name="${aa%.*}"
bb=$(basename "$b")
b_extension="${bb##*.}"
b_name="${bb%.*}"
if [ $a_extension == "tif" ]; then
    thresholds=`qauto $a /tmp/$a_name.png 2>&1 |  cut -f2 -d=`
    a=/tmp/$a_name.png
    # use the same threshods for the second image
    qeasy $thresholds $b /tmp/$b_name.png
    b=/tmp/$b_name.png

    # HACK! make it invariant to illumination changes?
    gblur 2 $a   | plambda  - " x 4 *    x(-1,0)  -1 *  +  x(1,0)  -1 * +  x(0,-1)  -1 * +    x(0,1)  -1 * + " | qauto - $a
    gblur 2 $b   | plambda  - " x 4 *    x(-1,0)  -1 *  +  x(1,0)  -1 * +  x(0,-1)  -1 * +    x(0,1)  -1 * + " | qauto - $b

fi


$rel_path_script/SGBM $a $b $3 $4 $5 $6 $7 $8 $9
plambda $3 "x isnan 0 255 if" | iion - ${3}.mask.png

rm $a $b
