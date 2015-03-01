
LR_REVERSE=1 CC_POSTPROCESS=0 time ../bin/iip_stereo_correlation_multi_win2_newversion -i 1 -n 4 -p 4 -W 5 -x 9 -y 9 -r 1 -D 1 -d 1 -t -1 -s 0 -b 1.25 -o -0.25 -O 25 -c 0 -f 0 -m -15 -M 5 -P 32 a.tif b.tif pixell.tif maskl.tif [pixelr.tif maskr.tif]


# in python
#common.run("%s -i 1 -n 4 -p 4 -W 5 -x 9 -y 9 -r 1 -D 1 -d 1 -t -1 -s 0 -b 1.25 -o -0.25 -O 25 -c 0 -f 0 -P 32 -m %d -M %d %s %s %s %s" %(msmw2_binary,
#            disp_min, disp_max, im1, im2, out_disp, out_mask))

