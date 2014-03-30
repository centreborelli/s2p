#!/bin/bash

TMP_FILE="/tmp/list_srtm_tiles.txt"

function which_srtm_tiles ()
# $1: path to rpc.xml file
{
    rpc_file=$1

    lon=`grep FIRST_LON $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    lat=`grep FIRST_LAT $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    echo $lon $lat | srtm4_which_tile

    lon=`grep FIRST_LON $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    lat=`grep LAST_LAT $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    echo $lon $lat | srtm4_which_tile

    lon=`grep LAST_LON $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    lat=`grep FIRST_LAT $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    echo $lon $lat | srtm4_which_tile

    lon=`grep LAST_LON $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    lat=`grep LAST_LAT $rpc_file | cut -d ">" -f 2 | cut -d "<" -f 1`
    echo $lon $lat | srtm4_which_tile
}

touch $TMP_FILE
for f in *; do
    if [ -d $f ]; then
        if ls $f/dataset_* &> /dev/null; then
            # the dataset has subdatasets (multidate)
            for ff in $f/dataset_*/rpc??.xml; do
                which_srtm_tiles $ff >> $TMP_FILE
            done
        else
            # the dataset has no subdatasets
            for ff in $f/rpc??.xml; do
                which_srtm_tiles $ff >> $TMP_FILE
            done
        fi
    fi
done

# remove duplicate entries
sort -u $TMP_FILE
rm $TMP_FILE
