#!/bin/bash

if [ $# -ne 3 ]; then echo "three command line arguments are needed; the path to the 1pe and multi_pe and the testcase name" ; exit;  fi

diff_lines=`diff <(ncdump $1) <(ncdump $2) | wc -l`
echo -e "\n========================================================================================="
if [ $diff_lines -gt 23 ]
then 
    echo -e "Failure for test '"$3"'\nComparing against "$1" and "$2" resulted in "$diff_lines" differences, non-reproducible">&2
    echo "========================================================================================="
    exit
else
    echo "Success for test "$3
    echo "========================================================================================="
fi

