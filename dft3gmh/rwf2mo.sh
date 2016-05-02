#!/bin/bash

usage="$(basename $0) [ G09 RWF File Name ] [ (optional) User-defined MO File Name ]"


if [[ $# -eq 0 ]]
then

    echo -e "\n${usage}\n"

    exit

elif [[ $# -ge 1 ]] && [[ $1 = "-h" ]]
then

    echo -e "\n${usage}\n"

    exit

elif [[ $# -eq 1 ]] && [[ $1 != "-h" ]]
then

    rwfname=$1

    moname="MO.deb"

elif [[ $# -eq 2 ]] && [[ $1 != "-h" ]]
then

    rwfname=$1

    moname=$2

else

    echo -e "\nExecution Error! Type \"$(basename $0) -h\" for more information\n\n "

    exit 9

fi

echo -e "\nRWF File = ${rwfname} ; MO File = ${moname}\n\n"

rwfdump "${rwfname}" "tmp.mo" 524R &&

~/bin/zm.gmx.tools/g_rwfdump2mo_d.py "tmp.mo" "${moname}" &&

rm -rf "tmp.mo" 




