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
#============================================================================#

## $Id: $

########################################################################
#                               libDippy                               #
########################################################################

DIPPYdir = $(pythondir)/dippy

DIPPYlibdir = $(pyexecdir)/dippy

DIPPY_PYTHON = __init__.py dipapi.py dippy.py

DIPPYlib_LTLIBRARIES = _dippy.la

# List all source files for this library, including headers
_dippy_la_SOURCES = \
	DippyDecompAlgo.cpp    \
	DippyDecompAlgo.h      \
	DippyDecompApp.cpp    \
	DippyDecompApp.h      \
	DippyPythonUtils.cpp    \
	DippyPythonUtils.h      \
	DippyDecompCut.h     \
	DippyModule.cpp   \
	DippySolve.cpp

########################################################################
#                            Additional flags                          #
########################################################################

# This is for libtool
_dippy_la_CPPFLAGS = $(PYTHON_CFLAGS) $(DIPLIB_CFLAGS)
_dippy_la_LIBADD = ../libDecomp.la $(CBCLIB_LIBS)
_dippy_la_LDFLAGS = $(LT_LDFLAGS) -module -avoid-version -export-symbols-regex initmyext -noundefined 

if DEPENDENCY_LINKING
_dippy_la_LIBADD += $(PYTHON_LIBS) $(DIPLIB_LIBS)
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I.. -I`$(CYGPATH_W) $(srcdir)` -I`$(CYGPATH_W) $(srcdir)/../` 

########################################################################
#                     Includes for LP Solver                           #
########################################################################
#because CGL depends on CLP, we need this in any case

if CLP_LP_SOLVER
  _dippy_la_CPPFLAGS += -D__DECOMP_LP_CLP__
endif

if CPLEX_LP_SOLVER
  _dippy_la_CPPFLAGS += -D__DECOMP_LP_CPX__
endif

########################################################################
#                     Includes for IP Solver                           #
########################################################################

if CBC_IP_SOLVER 
if COIN_HAS_CBC
  _dippy_la_CPPFLAGS += -D__DECOMP_IP_CBC__
endif
endif

if CPLEX_IP_SOLVER
  _dippy_la_CPPFLAGS += -D__DECOMP_IP_CPX__
endif

