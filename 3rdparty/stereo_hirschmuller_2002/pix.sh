#!/bin/bash

if [ "$3" == "" ]; then
   echo "Useage:"
   echo "   $0 im1 im2 out.pgm [mindisp(0) maxdisp(60) LoG(0) regionRadius(3) maxPerPixelError(20) validateRtoL(1) texture(0.2)]"
   echo ""
   echo "   LoG: Laplacian of Gaussian preprocess 1:enabled 0:disabled"
   echo "   regionRadius: radius of the window"
   echo "   maxPerPixelError: <0 disabled"
   echo "   validateRtoL: 1:enabled 0:disabled"
   echo "   texture: normalized diffrence between best and second best minimum  0:disabled"
   echo ""
   echo "   Implements: H. Hirschmuller, P.R. Innocent, and J.Garibaldi. "
   echo "   \"Real-Time Correlation-Based Stereo Vision with Reduced Border Errors.\" "
   echo "   Int. J. Comput. Vision 47, 1-3 2002"
   exit 1
fi

a=$1
b=$2
out=$3
im=$4
iM=$5

LoG=$6
rad=$7
ppe=$8
lr=$9
txt=${10}

echo "  Params:: LoG=$LoG regionRadius=$rad maxPerPixelError=$ppe validateRtoL=$lr texture=$txt"

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
fi


# HACKS TO FIX THE DISPARITY CONVENTIONS
let m=-iM
let M=-im

#### Logarithm of Gaussians?
if [ "$LoG" == "1" ]; then
   LoG="LoG"
else
   LoG=""
fi


if(( m < 0 )); 
then 
   ## HACK BECAUSE THE ALGORITHM DOES NOT COMPUTE NEGATIVE DISPARITIES
   w=`imprintf %w $a`
   h=`imprintf %h $a`
   let neww=w-m
   crop $b /tmp/b.png 0 0 $neww $h
   crop $a /tmp/a.png $m 0 $neww $h
   let newM=M-m
   let newm=0
   #echo java -classpath build/boofcv_stereo.jar pix${LoG} /tmp/a.png /tmp/b.png /tmp/out.pgm $newm $newM $rad $ppe $lr $txt
   time java -classpath build/boofcv_stereo.jar pix${LoG} /tmp/a.png /tmp/b.png /tmp/out.pgm $newm $newM  $rad $ppe $lr $txt 
   
   # mismatches
   let nomatchtag=newM+1
   plambda /tmp/out.pgm "1 x $nomatchtag = - 255 *" | iion - $out.mask.png

   # restore the order of disparities
   let mm=-m
   plambda /tmp/out.pgm "0 x $mm - -" | iion - $3

   ## also recover the false color disparity 
   crop $3 $3 $mm 0 $w $h
   crop  /tmp/out.pgm.png $3.png $mm 0 $w $h
   crop  $out.mask.png $out.mask.png $mm 0 $w $h

else 
   #echo java -classpath build/boofcv_stereo.jar pix${LoG} $a $b $out $m $M  $rad $ppe $lr $txt
   time java -classpath build/boofcv_stereo.jar pix${LoG} $a $b $out $m $M  $rad $ppe $lr $txt

   # mismatches
   let nomatchtag=M+1
   echo $nomatchtag
   plambda $out "1 x $nomatchtag = - 255 *" | iion - $out.mask.png

   # restore the order of disparities
   plambda $out "0 x -" | iion - $out
fi


