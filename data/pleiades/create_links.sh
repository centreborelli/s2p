#!/bin/bash

function create_links ()
# $1: directory where to create the links, and containing the txt file $2
# $2: txt file containing the list of paths to images
# $3: suffix for link names. "" for panchro, "_color" for multispectral images
{
    dataset=$1
    files_list=$2
    suffix=$3
    i=0
    for image in `cat $dataset/$files_list`
        do
            i=$(($i+1))
            # image
            abs_path=$pleiades_dir/$image
            link_name=`printf im%02d$suffix.tif $i`
            ln -sf $abs_path $dataset/$link_name
            # preview
            dir_path=`dirname $image`
            abs_path=$pleiades_dir/$dir_path/PREVIEW_*.JPG
            link_name=`printf prev%02d$suffix.jpg $i`
            ln -sf $abs_path $dataset/$link_name
            # rpc
            abs_path=$pleiades_dir/$dir_path/RPC_*.XML
            link_name=`printf rpc%02d$suffix.xml $i`
            cp $abs_path $dataset/
            ln -sf `basename $abs_path` $dataset/$link_name
    done
}


###############
# Main script #
###############

# arguments:
# $1: absolute path to the folder containing the pleiades data

# check input (ie the pleiades data folder)
if [ ! $1 ] ; then
    printf "\tusage: %s pleiades_data_folder_path\n" $0
    exit
fi

pleiades_dir=$1
if [ ! -d $pleiades_dir ] ; then
    printf "\tincorrect path to pleiades data folder\n"
    exit
fi


# create the symlinks
for dataset in `find * -type d`
    do
        if [ -f "$dataset/images_paths_panchro.txt" ] ; then
            create_links $dataset "images_paths_panchro.txt" ""
        fi

        if [ -f "$dataset/images_paths_ms.txt" ] ; then
            create_links $dataset "images_paths_ms.txt" "_color"
        fi
    done
