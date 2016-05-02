#!/bin/bash

usage="$(basename $0) [ G09 RWF File Name ] [ (optional) User-defined sAO File Name ]"


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

    saoname="sao.deb"

elif [[ $# -eq 2 ]] && [[ $1 != "-h" ]]
then

    rwfname=$1

    saoname=$2

else

    echo -e "\nExecution Error! Type \"$(basename $0) -h\" for more information\n\n "

    exit 9

fi

echo -e "\nRWF File = ${rwfname} ; SAO File = ${saoname}\n\n"

rwfdump "${rwfname}" "tmp.sao" 514R &&

~/bin/zm.gmx.tools/g_rwfdump2sao_d.py "tmp.sao" "${saoname}" &&

rm -rf "tmp.sao" 




