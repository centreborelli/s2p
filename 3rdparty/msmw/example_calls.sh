export PATH=$PATH:./bin/

a=input/aNY2.tif
b=input/bNY2.tif
d=""
nm=NY
disp="-m -50 -M 90"


##### MULTI WINDOWS, NO MULTISCALE 


# 1 scale, 1 window 5x5 
code=S1W1
time ./bin/iip_stereo_correlation_multi_win2 -W 1 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 1 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi

# 1 scale, 3 windows 5x5 (equivalent)
code=S1W3
time ./bin/iip_stereo_correlation_multi_win2 -W 3 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 1 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi
# 1 scale, 5 windows 5x5 (equivalent)
code=S1W5
time ./bin/iip_stereo_correlation_multi_win2 -W 5 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 1 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi
# 1 scale, 9 windows 5x5 (equivalent)
code=S1W9
time ./bin/iip_stereo_correlation_multi_win2 -W 9 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 1 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi





# MULTI WINDOWS, MULTISCALE

# 4 scale, 1 window 5x5 
code=S4W1
time ./bin/iip_stereo_correlation_multi_win2 -W 1 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 4 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi


# 4 scale, 1 window 5x5 
code=Ss4W1
time DONT_USE_SUBPIX_AT_LOW_SCALES=1 ./bin/iip_stereo_correlation_multi_win2 -W 1 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 4 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi



# 4 scale, 1 window 5x5 
code=S4W5
time ./bin/iip_stereo_correlation_multi_win2 -W 5 -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 4 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi



# 4 scales, 5 windows 5x5, multiples ventanas equivalantes a 5x5 multiescala
code=S4W5Par
time ./bin/iip_stereo_correlation_multi_window -d 0 -p 1 -i 1 -r 1 -s 0 -f 0.2 $disp -x 5 -n 4 $a $b output/${nm}_${code}_{out,mask,out2,mask2}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi



# ADDITIONAL FILTERING

# 4 scale, 1 window 5x5 
code=S4W5Md
time ./bin/iip_stereo_correlation_multi_win2 -W 5 -d 1 -p 1 -i 1 -r 1 -s 0 -f 0.2 -n 4 $disp -x 5 -y 5 $a $b output/${nm}_${code}_{out,mask}.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_out.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outmask.tif
if [ $d ]; then  
plambda output/${nm}_${code}_out.tif $d "x y - fabs" | iion - output/${nm}_${code}_outdiffmask.tif
./bin/iip_stereo_apply_mask output/${nm}_${code}_outdiffmask.tif output/${nm}_${code}_mask.tif output/${nm}_${code}_outdiffmask.tif
fi

