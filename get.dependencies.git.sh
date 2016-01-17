#!/bin/bash

for i in `cat Dependencies | tr '\t' ' ' | tr -s ' '| cut -d ' ' -f 2`
do
    svn_repo=`echo $i | cut -d '/' -f 5`
    if [ $svn_repo = 'Data' ]; then
	svn co $i Data
    elif [ $svn_repo = 'BuildTools' ]; then
	if [ `echo $i | cut -d '/' -f 6` != 'stable' ]; then
	    proj=ThirdParty-`echo $i | cut -d '/' -f 7`
	    git clone https://github.com/coin-or-tools/$proj
	    if [ `echo $i | cut -d '/' -f 8` != 'trunk' ]; then
		cd $proj
		branch=`echo $i | cut -d '/' -f 8-9`
		git fetch
		git checkout $branch
		cd ..
	    fi
	fi
    else
	if [ $svn_repo = "CHiPPS" ]; then
	    git_repo=CHiPPS-`echo $i | cut -d '/' -f 6`
	    proj=`echo $i | cut -d '/' -f 6`
	    if [ `echo $i | cut -d '/' -f 7` = 'trunk' ]; then
		branch=master
	    else
		branch=`echo $i | cut -d '/' -f 7-8`
	    fi	
	else
	    git_repo=`echo $i | cut -d '/' -f 5`
	    proj=`echo $i | cut -d '/' -f 5`
	    if [ `echo $i | cut -d '/' -f 6` = 'trunk' ]; then
		branch=master
	    else
		branch=`echo $i | cut -d '/' -f 6-7`
	    fi	
	fi
	echo "Checking out" $proj "..."
	mkdir $proj
	cd $proj
	git init
	git remote add origin https://github.com/coin-or/$git_repo
	git config core.sparsecheckout true
	echo $proj/ >> .git/info/sparse-checkout
	git fetch
	git checkout $branch
	cd ..
    fi
done
