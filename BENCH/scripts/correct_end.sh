#!/bin/bash

# To check that croco has correct end
# => because with nvfortran return ieee_inexact even on success

# enforce strict
set -u

# vars
logfile=$(mktemp)

# run
"$@" 2>&1 | tee ${logfile}

# check
grep -q 'MAIN: DONE' ${logfile}
status=$?

# display error
if [[ ${status} != 0 ]]; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 1>&2
    echo "CROCO didn't reached end of execution !" 1>&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 1>&2
    tail -n 16 ${logfile}
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 1>&2
    echo "CROCO didn't reached end of execution !" 1>&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 1>&2
fi

# remove log
rm -f ${logfile}

# finish
exit $status
