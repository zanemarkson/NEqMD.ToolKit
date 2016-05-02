#!/bin/bash

usage="$(basename $0) [ G09 RWF File Name ] [ (optional) User-defined Dipole File Prefix ]"


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

    dipname="none"

elif [[ $# -eq 2 ]] && [[ $1 != "-h" ]]
then

    rwfname=$1

    dipname=$2

else

    echo -e "\nExecution Error! Type \"$(basename $0) -h\" for dipre information\n\n "

    exit 9

fi

echo -e "\nRWF File = ${rwfname} ; Dipole File Prefix = ${dipname}\n\n"

rwfdump "${rwfname}" "tmp.dip" 518R &&


if [[ "${dipname}" == "none" ]]
then
	~/bin/zm.gmx.tools/g_rwfdump2dip_d.py "tmp.dip" &&

	rm -rf "tmp.dip"

else	
	~/bin/zm.gmx.tools/g_rwfdump2dip_d.py "tmp.dip" "${dipname}" &&

	rm -rf "tmp.dip" 
fi




