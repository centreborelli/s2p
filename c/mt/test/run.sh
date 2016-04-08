#!/bin/sh
#
# test scripts runner

LOGFILE=test.log

# log file
_log_init() {
    rm -f $LOGFILE
    touch $LOGFILE
}

_log_clean() {
    rm -f $LOGFILE
}

_log() {
    echo "  $*"
    $* >> $LOGFILE 2>&1 
}

# run test scripts from the same folder
for SCRIPT in ${0%/*}/[0-9]*-*.sh; do
    echo - $SCRIPT
    . $SCRIPT
done;
