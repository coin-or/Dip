#============================================================================#
# This file is part of the DIP Solver Framework.                             #
#                                                                            #
# DIP is distributed under the Eclipse Public License as part of the         #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Author: Matthew Galati, SAS Institute Inc.                                 #
#                                                                            #
# Conceptual Design: Matthew Galati, SAS Institute Inc.                      #
#                    Ted Ralphs, Lehigh University                           #
#                                                                            #
# Copyright (C) 2002-2011, Lehigh University, Matthew Galati, Ted Ralphs     #
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

tests: all
	for dir in $(SUBDIRS); do \
	  if test -r $$dir/test/Makefile; then \
	    (cd $$dir; $(MAKE) test) \
	  fi; \
	done

doxydocd:
	(cd $(srcdir)/doxydoc/Dip; doxygen doxygen.conf)

doxydoc: 
	for dir in $(SUBDIRS_DOXY); do \
	  if test -r $(srcdir)/doxydoc/$$dir/doxygen.conf; then \
	    (cd $(srcdir)/doxydoc/$$dir; doxygen doxygen.conf) \
	  fi; \
	done

unitTests: tests

.PHONY: test tests unitTests doxydoc doxydocd

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = coin_subdirs.txt

include BuildTools/Makemain.inc
