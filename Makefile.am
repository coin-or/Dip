#============================================================================#
# This file is part of the DIP Solver Framework.                             #
#                                                                            #
# DIP is distributed under the Eclipse Public License as part of the         #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Authors: Matthew Galati, SAS Institute Inc.                                #
#          Ted Ralphs, Lehigh University                                     #
#          Jiadong Wang, Lehigh University                                   #
#                                                                            #
# Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs     #
# All Rights Reserved.                                                       #
#============================================================================#

AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = doxydoc/all/doxygen.conf  \
	doxydoc/Dip/doxygen.conf    \
	doxydoc/CoinUtils/doxygen.conf \
	doxydoc/Osi/doxygen.conf       \
	doxydoc/Clp/doxygen.conf       \
	doxydoc/Cgl/doxygen.conf       \
	doxydoc/Cbc/doxygen.conf       \
	doxydoc/Alps/doxygen.conf 

########################################################################
#                          Subdirectories                              #
########################################################################
# subdirs is set by configure as the list of all subdirectories to recurse
SUBDIRS = $(subdirs)
SUBDIRS_DOXY = $(SUBDIRS) all

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd Dip; $(MAKE) test

test_cbc: all
	cd Dip; $(MAKE) test_cbc

test_symphony: all
	cd Dip; $(MAKE) test_symphony

test_cpx: all
	cd Dip; $(MAKE) test_cpx

test_grb: all
	cd Dip; $(MAKE) test_grb

tests: all
	for dir in $(SUBDIRS); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; $(MAKE) test) \
	  fi; \
	done

unitTests: tests

# Generate doxygen doc'n in subdirectories (except @PACKAGE_NAME@) if a doxydoc
# directory is present, then do the base, if present.

doxydoc:
	for dir in $(subdirs) ; do \
	  if test $$dir != @PACKAGE_NAME@ && test -r $$dir/doxydoc ; then \
	    (cd $$dir ; $(MAKE) doxydoc) \
	  fi ; \
	done ; \
	if test -r doxydoc/doxygen.conf ; then \
	  doxygen doxydoc/doxygen.conf ; \
	fi

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )

# DocInstallDir is defined in Makemain.inc and is specific to the package.
# For the short term, adopt the notion that we install only the package
# doxydoc.

install-doxydoc: doxydoc
	if test -r doxydoc/doxygen.conf ; then \
	  $(mkdir_p) $(DocInstallDir) ; \
	  cp -R doxydoc $(DocInstallDir) ; \
	fi

uninstall-doxydoc:
	rm -rf $(DocInstallDir)/doxydoc

clean-local: clean-doxydoc
# install-data-local: install-doxydoc
uninstall-local: uninstall-doxydoc

.PHONY: test tests unitTests doxydoc

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = coin_subdirs.txt

include BuildTools/Makemain.inc
