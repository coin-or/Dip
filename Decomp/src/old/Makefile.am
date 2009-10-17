#============================================================================#
# This file is part of the DECOMP Solver Framework.                          #
#                                                                            #
# DECOMP is distributed under the Common Public License as part of the       #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)        #
#                                                                            #
# Copyright (C) 2002-2007, Lehigh University, Matthew Galati, and Ted Ralphs.#
# All Rights Reserved.                                                       #
#============================================================================#

## $Id: $

########################################################################
#                               libDecomp                              #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libDecomp.la

# List all source files for this library, including headers
libDecomp_la_SOURCES =

# DecompAlgo
libDecomp_la_SOURCES +=   \
	DecompAlgo.cpp    \
	DecompAlgo.h      \
	DecompAlgoC.cpp   \
	DecompAlgoC.h     \
	DecompAlgoD.cpp   \
	DecompAlgoD.h     \
	DecompAlgoPC.cpp  \
	DecompAlgoPC.h    \
	DecompAlgoRC.cpp  \
	DecompAlgoRC.h    \
	DecompCut.cpp     \
	DecompCut.h       \
	DecompCutOsi.h    \
	DecompCutPool.cpp \
	DecompCutPool.h  

# DecompApp
libDecomp_la_SOURCES +=  \
	DecompApp.cpp    \
	DecompApp.h      

# other 
libDecomp_la_SOURCES +=    \
	DecompConstants.h       \
	DecompStats.cpp         \
	DecompStats.h           \
	DecompDebug.cpp         \
	DecompModel.h           \
	DecompParam.h 	        \
	DecompPortable.h        \
	DecompTypes.h           \
	DecompVar.cpp           \
	DecompVar.h             \
	DecompVarPool.cpp       \
	DecompVarPool.h         \
	DecompBranch.cpp        \
	DecompMemPool.h         \
	DecompSolution.h        \
	DecompConstraintSet.cpp \
	DecompConstraintSet.h   \
	DecompWaitingCol.h	\
	DecompWaitingRow.h

# Util
libDecomp_la_SOURCES +=    \
	UtilGraphLib.cpp   \
	UtilGraphLib.h     \
	UtilHash.cpp       \
	UtilHash.h         \
	UtilMacros.cpp     \
	UtilMacros.h       \
	UtilParameters.cpp \
	UtilParameters.h   

#TODO - alps interface
if USE_ALPS_FRAMEWORK
libDecomp_la_SOURCES +=        \
	BcpsDecompModel.cpp    \
	BcpsDecompModel.h      \
	BcpsDecompSolution.cpp \
	BcpsDecompSolution.h   \
	BcpsDecompTreeNode.cpp \
	BcpsDecompTreeNode.h   \
	BcpsDecompNodeDesc.h   \
	UtilMacrosAlps.cpp     \
	UtilMacrosAlps.h       
endif

#TempFix
libDecomp_la_SOURCES +=            \
	OsiNullSolverInterface.hpp \
	OsiData.hpp

# This is for libtool (on Windows)
libDecomp_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.

LDADD=
ADDLIBS=

AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src`       \
	-I`$(CYGPATH_W) $(OSISRCDIR)/inc`

# CGL
if USE_CGL_CUTS
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir`
endif

#TODO: these should be optional
# CBC
# AM_CPPFLAGS += \
# 	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
# 	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc`

# ALPS/BCPS
if USE_ALPS_FRAMEWORK
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(ALPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(ALPSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BCPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BCPSOBJDIR)/inc`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                     Includes for LP Solver                           #
########################################################################
if CLP_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`   \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la    \
	$(CLPOBJDIR)/src/libClp.la
  ADDLIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
  AM_CPPFLAGS += -D__DECOMP_LP_CLP__
endif

if CPLEX_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
  AM_CPPFLAGS += -D__DECOMP_LP_CPX__
endif

if GLPK_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la 
  AM_CPPFLAGS += -D__DECOMP_LP_GLPK__
else 
endif

if OSL_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
  AM_CPPFLAGS += -D__DECOMP_LP_OSL__
endif

if SOPLEX_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
AM_CPPFLAGS += -D__DECOMP_LP_SPX__
endif

if XPRESS_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
  AM_CPPFLAGS += -D__DECOMP_LP_XPR__
endif


########################################################################
#                     Includes for IP Solver                           #
########################################################################
if CBC_IP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`   \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc`
  LDADD += $(OSIOBJDIR)/src/OsiCbc/libOsiCbc.la    \
	$(CBCOBJDIR)/src/libCbc.la
  ADDLIBS  += `cat $(CBCOBJDIR)/cbc_addlibs.txt`
  AM_CPPFLAGS += -D__DECOMP_IP_CBC__
endif

if CPLEX_IP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
  AM_CPPFLAGS += -D__DECOMP_IP_CPX__
endif

if GLPK_IP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la 
  AM_CPPFLAGS += -D__DECOMP_IP_GLPK__
else 
endif

if OSL_IP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
  AM_CPPFLAGS += -D__DECOMP_IP_OSL__
endif

if SOPLEX_IP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
AM_CPPFLAGS += -D__DECOMP_IP_SPX__
endif

if XPRESS_IP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
  AM_CPPFLAGS += -D__DECOMP_IP_XPR__
endif

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	../inc/config_decomp.h

#DecompAlgo
include_HEADERS +=        \
	DecompAlgo.h      \
	DecompAlgoC.h     \
	DecompAlgoD.h     \
	DecompAlgoPC.h    \
	DecompAlgoRC.h    \
	DecompCut.h       \
	DecompCutOsi.h    \
	DecompCutPool.h 

# DecompApp
include_HEADERS +=  \
	DecompApp.h      

# other 
include_HEADERS +=            \
	DecompConstants.h     \
	DecompStats.h         \
	DecompModel.h         \
	DecompParam.h 	      \
	DecompPortable.h      \
	DecompTypes.h         \
	DecompVar.h           \
	DecompVarPool.h       \
	DecompMemPool.h       \
	DecompSolution.h      \
	DecompConstraintSet.h \
	DecompWaitingCol.h    \
	DecompWaitingRow.h

# Util
include_HEADERS +=         \
	UtilGraphLib.h     \
	UtilHash.h         \
	UtilMacros.h       \
	UtilParameters.h  

if USE_ALPS_FRAMEWORK
include_HEADERS +=             \
	BcpsDecompModel.h      \
	BcpsDecompSolution.h   \
	BcpsDecompTreeNode.h   \
	BcpsDecompNodeDesc.h   \
	UtilMacrosAlps.h       
endif

# TempFix
include_HEADERS +=                  \
	OsiNullSolverInterface.hpp  \
	OsiData.hpp

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = DecompConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)
