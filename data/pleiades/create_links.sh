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
            link_name=`printf im%02d$suffix.tif $i`
            ln -sf $image $dataset/$link_name
            # preview
            abs_path=`dirname $image`/PREVIEW_*.JPG
            link_name=`printf prev%02d$suffix.jpg $i`
            ln -sf $abs_path $dataset/$link_name
            # rpc
            abs_path=`dirname $image`/RPC_*.XML
            link_name=`printf rpc%02d$suffix.xml $i`
            cp $abs_path $dataset/
            ln -sf `basename $abs_path` $dataset/$link_name
            # dim (other xml files with dimensions informations)
            abs_path=`dirname $image`/DIM_*.XML
            link_name=`printf dim%02d$suffix.xml $i`
            cp $abs_path $dataset/
            ln -sf `basename $abs_path` $dataset/$link_name
            # roi mask
            abs_path=`dirname $image`/MASKS/ROI*.GML
            link_name=`printf roi_msk%02d$suffix.gml $i`
            cp $abs_path $dataset/
            ln -sf `basename $abs_path` $dataset/$link_name
            # cloud mask
            abs_path=`dirname $image`/MASKS/CLD*.GML
            link_name=`printf cld_msk%02d$suffix.gml $i`
            if [ -f $abs_path ]; then
                cp $abs_path $dataset/
                ln -sf `basename $abs_path` $dataset/$link_name
            fi
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

# step 1: parse the pleiades data folder to extract the paths to JP2.TIF images
for f in $pleiades_dir/*; do
    if [ -d $f ]; then
        mkdir -p `basename $f`
        if ls $f/dataset_* &> /dev/null; then
            # the dataset has subdatasets (multidate)
            for ff in $f/dataset_*; do
                mkdir -p `basename $f`/`basename $ff`
                find $ff | grep "JP2.TIF" | grep -v "JP2.TIF.ovr" | grep "_P_" >  `basename $f`/`basename $ff`/paths_panchro.txt
                find $ff | grep "JP2.TIF" | grep -v "JP2.TIF.ovr" | grep "_MS_" > `basename $f`/`basename $ff`/paths_ms.txt
            done
        else
            # the dataset has no subdatasets
            find $f | grep "JP2.TIF" | grep -v "JP2.TIF.ovr" | grep "_P_" >  `basename $f`/paths_panchro.txt
            find $f | grep "JP2.TIF" | grep -v "JP2.TIF.ovr" | grep "_MS_" > `basename $f`/paths_ms.txt
        fi
    fi
done

# step 2: create the symlinks
for dataset in `find * -type d`; do
    if [ -f "$dataset/paths_panchro.txt" ] ; then
        create_links $dataset "paths_panchro.txt" ""
    fi
    if [ -f "$dataset/paths_ms.txt" ] ; then
        create_links $dataset "paths_ms.txt" "_color"
    fi
done
