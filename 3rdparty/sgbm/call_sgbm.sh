#!/bin/bash

if [ "$3" == "" ]; then
   echo "Usage:"
   echo "  $0 im1 im2 out_disp out_cost out_mask mindisp maxdisp win P1 P2 lr"
   echo ""
   echo "  win: matched block size. It must be an odd number >=1."
   echo "  P1: The first parameter controlling the disparity smoothness."
   echo "  P2: The second parameter controlling the disparity smoothness. The"
   echo "    larger the values are, the smoother the disparity is. P1 is the penalty on"
   echo "    the disparity change by plus or minus 1 between neighbor pixels. P2 is the"
   echo "    penalty on the disparity change by more than 1 between neighbor pixels. The"
   echo "    algorithm requires P2 > P1"
   echo "  lr: max allowed difference in the left-right disparity check."
   echo ""
   echo "  Wrapper to opencv SGBM function, which implements a modified version"
   echo "  of Hirschmuller's Semi-Global Matching (SGM):"
   echo "  Hirschmuller, H. \"Stereo Processing by Semiglobal Matching and Mutual Information\""
   echo "  PAMI(30), No. 2, February 2008, pp. 328-34"
   exit 1
fi

a=$1
b=$2
disp=$3
cost=$4
mask=$5
im=$6
iM=$7
SAD_win=$8
P1=$9
P2=${10}
lr=${11}

# run sgbm
echo "sgbm $a $b $disp $cost $im $iM $SAD_win $P1 $P2 $lr"
sgbm $a $b $disp $cost $im $iM $SAD_win $P1 $P2 $lr

# create rejection mask. 0 means the pixel is rejected, 1 means the pixel is accepted.
# the points are either unmatched or not present in the original images $1/$2, we remove these points
plambda $disp "x 0 join" | backflow - $2 | plambda $disp $1 - "x isfinite y isfinite z isfinite or or" -o $mask
