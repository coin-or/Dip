#!/bin/bash

for i in `cat Dependencies | tr '\t' ' ' | tr -s ' '| cut -d ' ' -f 2`
do
    if [ `echo $i | cut -d '/' -f 5` == 'BuildTools' ]; then
	if [ `echo $i | cut -d '/' -f 6` != 'stable' ]; then
	    mkdir -p ThirdParty
	    svn co --non-interactive --trust-server-cert $i ThirdParty/`echo $i | cut -d '/' -f 7`
	fi
    else
	proj=`echo $i | cut -d '/' -f 5`
	svn co --non-interactive --trust-server-cert $i $proj
    fi
done
