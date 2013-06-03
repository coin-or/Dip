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
# Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs     #
# All Rights Reserved.                                                       #
#============================================================================#

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles
EXTRA_DIST = 

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

.PHONY: test unitTest 

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = dip.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Dip

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs dip > $(addlibsdir)/dip_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libDecomp.lib @DIPLIB_LIBS_INSTALLED@" > $(addlibsdir)/dip_addlibs.txt
else
	echo -L@abs_lib_dir@ -lDecomp @DIPLIB_LIBS_INSTALLED@ > $(addlibsdir)/dip_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/dip_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
