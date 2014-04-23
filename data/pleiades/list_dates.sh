#!/bin/bash

for f in *; do
    if [ -d $f ]; then
        if ls $f/dataset_* &> /dev/null; then
            # the dataset has subdatasets (multidate)
            for ff in $f/dataset_*; do
                echo $ff
                cat $ff/dim01.xml | grep "IMAGING_DATE" | cut -d '>' -f 2 | cut -d '<' -f 1
            done
        else
            # the dataset has no subdatasets
            echo $f
            cat $f/dim01.xml | grep "IMAGING_DATE" | cut -d '>' -f 2 | cut -d '<' -f 1
        fi
    fi
done
