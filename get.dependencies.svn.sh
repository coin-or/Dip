#!/bin/bash

for url in `cat Dependencies | tr '\t' ' ' | tr -s ' '| cut -d ' ' -f 2`
do
    if [ `echo $url | cut -d '/' -f 5` == 'BuildTools' ]; then
	if [ `echo $url | cut -d '/' -f 6` != 'stable' ]; then
	    mkdir -p ThirdParty
	    svn co --non-interactive --trust-server-cert $url ThirdParty/`echo $url | cut -d '/' -f 7`
	fi
    else
	svn_repo=`echo $url | cut -d '/' -f 5`
	if [ $svn_repo = "CHiPPS" ]; then
	    proj=`echo $url | cut -d '/' -f 6`
	elif [ $svn_repo = "Data" ]; then
	    proj=`echo $url | cut -d '/' -f 5-6`
	else
	    proj=`echo $url | cut -d '/' -f 5`
	fi
	svn co --non-interactive --trust-server-cert $url $proj
    fi
done
