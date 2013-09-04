#!/bin/bash

# arguments:
# $1: absolute path to the folder containing the pleiades data

# check input (ie the pleiades data folder)
if [ ! $1 ] ; then
    printf "\tusage: %s pleiades_folder_path\n" $0
    exit
fi

pleiades_dir=$1
if [ ! -d $pleiades_dir ] ; then
    printf "\tincorrect path to your pleiades folder\n"
    exit
fi

# create the symlinks
for dataset in images/*
    do
        # panchro
        if [ -f "$dataset/images_paths_panchro.txt" ] ; then
            i=0
            for image in `cat $dataset/images_paths_panchro.txt`
                do
                    i=$(($i+1))
                    abs_path=$pleiades_dir/$image
                    #echo $abs_path
                    link_name=`printf im%02d.tif $i`
                    #echo $dataset/$link_name
                    ln -sf $abs_path $dataset/$link_name
            done
        fi

        # multi-spectral
        if [ -f "$dataset/images_paths_ms.txt" ] ; then
            i=0
            for image in `cat $dataset/images_paths_ms.txt`
                do
                    i=$(($i+1))
                    abs_path=$pleiades_dir/$image
                    #echo $abs_path
                    link_name=`printf im%02d_color.tif $i`
                    #echo $dataset/$link_name
                    ln -sf $abs_path $dataset/$link_name
            done
        fi
    done
