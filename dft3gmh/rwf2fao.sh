#!/bin/bash

usage="$(basename $0) [ G09 RWF File Name ] [ (optional) User-defined FAO File Name ]"


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

    faoname="fao.deb"

elif [[ $# -eq 2 ]] && [[ $1 != "-h" ]]
then

    rwfname=$1

    faoname=$2

else

    echo -e "\nExecution Error! Type \"$(basename $0) -h\" for more information\n\n "

    exit 9

fi

echo -e "\nRWF File = ${rwfname} ; FAO File = ${faoname}\n\n"

rwfdump "${rwfname}" "tmp.fao" 536R &&

~/bin/zm.gmx.tools/g_rwfdump2fao_d.py "tmp.fao" "${faoname}" &&

rm -rf "tmp.fao" 




