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
	DecompAlgoCGL.cpp \
	DecompAlgoCGL.h   \
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
	Decomp.h                \
	DecompStats.cpp         \
	DecompStats.h           \
	DecompDebug.cpp         \
	DecompModel.cpp         \
	DecompModel.h           \
	DecompParam.h 	        \
	DecompVar.cpp           \
	DecompVar.h             \
	DecompVarPool.cpp       \
	DecompVarPool.h         \
	DecompBranch.cpp        \
	DecompMemPool.h         \
	DecompSolution.h        \
	DecompSolverResult.h    \
	DecompConstraintSet.cpp \
	DecompConstraintSet.h   \
	DecompWaitingCol.h	\
	DecompWaitingRow.h

# Util
libDecomp_la_SOURCES +=      \
	UtilGraphLib.cpp     \
	UtilGraphLib.h       \
	UtilHash.cpp         \
	UtilHash.h           \
	UtilMacros.cpp       \
	UtilMacros.h         \
	UtilMacrosDecomp.cpp \
	UtilMacrosDecomp.h   \
	UtilParameters.cpp   \
	UtilParameters.h     \
	UtilKnapsack.cpp     \
	UtilKnapsack.h       \
	UtilTimer.h

libDecomp_la_SOURCES +=        \
	AlpsDecompModel.cpp    \
	AlpsDecompModel.h      \
	AlpsDecompSolution.h   \
	AlpsDecompTreeNode.cpp \
	AlpsDecompTreeNode.h   \
	AlpsDecompNodeDesc.h   \
	AlpsDecompParam.h      \
	UtilMacrosAlps.cpp     \
	UtilMacrosAlps.h       

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
	-I`$(CYGPATH_W) $(OSISRCDIR)/inc`       \
	-DCOIN_ASSERT

# CGL
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

# ALPS/BCPS
AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(ALPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(ALPSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BCPSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BCPSOBJDIR)/inc`


# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc
#can't se just include <build-dir>/include/coin?

########################################################################
#                     Includes for LP Solver                           #
########################################################################
#because CGL depends on CLP, we need this in any case
AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	  -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`   \
	  -I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la    \
	  $(CLPOBJDIR)/src/libClp.la
ADDLIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`

if CLP_LP_SOLVER
  AM_CPPFLAGS += -D__DECOMP_LP_CLP__
endif

if CPLEX_LP_SOLVER
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
  AM_CPPFLAGS += -D__DECOMP_LP_CPX__
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

########################################################################
#                Headers that need to be installed                     #
########################################################################
# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
includecoindir = $(includedir)/coin

includecoin_HEADERS = \
	../inc/config_dip.h

#DecompAlgo
includecoin_HEADERS +=        \
	DecompAlgo.h      \
	DecompAlgoC.h     \
	DecompAlgoD.h     \
	DecompAlgoPC.h    \
	DecompAlgoRC.h    \
	DecompAlgoCGL.h   \
	DecompCut.h       \
	DecompCutOsi.h    \
	DecompCutPool.h 

# DecompApp
includecoin_HEADERS +=  \
	DecompApp.h      

# other 
includecoin_HEADERS +=            \
	Decomp.h              \
	DecompStats.h         \
	DecompModel.h         \
	DecompParam.h 	      \
	DecompVar.h           \
	DecompVarPool.h       \
	DecompMemPool.h       \
	DecompSolution.h      \
	DecompSolverResult.h  \
	DecompConstraintSet.h \
	DecompWaitingCol.h    \
	DecompWaitingRow.h

# Util
includecoin_HEADERS +=         \
	UtilGraphLib.h     \
	UtilHash.h         \
	UtilMacros.h       \
	UtilMacrosDecomp.h \
	UtilParameters.h   \
	UtilKnapsack.h     \
	UtilTimer.h

includecoin_HEADERS +=             \
	AlpsDecompModel.h      \
	AlpsDecompSolution.h   \
	AlpsDecompTreeNode.h   \
	AlpsDecompNodeDesc.h   \
	AlpsDecompParam.h      \
	UtilMacrosAlps.h       


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
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
