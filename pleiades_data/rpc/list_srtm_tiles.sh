#!/bin/bash

for f in *
do
    if [ -d "$f" ]
    then
        lon=`grep FIRST_LON $f/rpc01.xml | cut -d ">" -f 2 | cut -d "<" -f 1`
        lat=`grep FIRST_LAT $f/rpc01.xml | cut -d ">" -f 2 | cut -d "<" -f 1`
        echo $lon $lat | srtm4
    fi
done
