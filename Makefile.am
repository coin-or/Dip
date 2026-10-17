#============================================================================#
# This file is part of the DIP Solver Framework.                             #
#                                                                            #
# DIP is distributed under the Eclipse Public License as part of the         #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Author: Matthew Galati,  SAS Institute Inc.                                #
#                                                                            #
# Conceptual Design: Matthew Galati, SAS Institute Inc.                      #
#                    Ted Ralphs, Lehigh University                           #
#                                                                            #
# Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs     #
# All Rights Reserved.                                                       #
#============================================================================#

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

########################################################################
#             Additional files to be included in tarball               #
########################################################################

pkgconfiglib_DATA = dip.pc

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

if COIN_HAS_CBC
test_cbc: all
	cd test; $(MAKE) test_cbc
endif

if COIN_HAS_SYMPHONY
test_symphony: all
	cd test; $(MAKE) test_symphony
endif

if COIN_HAS_CPX
test_cpx: all
	cd test; $(MAKE) test_cpx
endif

if COIN_HAS_GRB
test_grb: all
	cd test; $(MAKE) test_grb
endif

if COIN_HAS_XPR
test_xpr: all
	cd test; $(MAKE) test_xpr
endif

unitTest: test

# Doxygen documentation

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc uninstall-doxygen-docs

.PHONY: test unitTest doxydoc
