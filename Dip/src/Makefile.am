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

SUBDIRS = .
if HAVE_PYTHON
SUBDIRS += dippy
endif

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

# List all additionally required libraries
if DEPENDENCY_LINKING
libDecomp_la_LIBADD = $(DIPLIB_LIBS)
endif

# This is for libtool
libDecomp_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                            DIP binary                                #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = dip

# List all source files for this executable, including headers
dip_SOURCES = DecompMain.cpp

# Additional libraries
dip_LDADD = libDecomp.la $(DIPLIB_LIBS)
if COIN_HAS_HMETIS
dip_LDADD +=  $(HMETISLIB)
endif
dip_DEPENDENCIES = libDecomp.la $(DIPLIB_DEPENDENCIES)

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.

AM_CPPFLAGS = $(DIPLIB_CFLAGS)

if COIN_HAS_HMETIS
AM_CPPFLAGS += -I`$(CYGPATH_W) $(HMETISINCDIR)`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                     Includes for LP Solver                           #
########################################################################
#because CGL depends on CLP, we need this in any case

if CLP_LP_SOLVER
  AM_CPPFLAGS += -D__DECOMP_LP_CLP__
endif

if CPLEX_LP_SOLVER
  AM_CPPFLAGS += -D__DECOMP_LP_CPX__
endif

########################################################################
#                     Includes for IP Solver                           #
########################################################################

if CBC_IP_SOLVER 
if COIN_HAS_CBC
  AM_CPPFLAGS += -D__DECOMP_IP_CBC__
endif
endif

if CPLEX_IP_SOLVER
  AM_CPPFLAGS += -D__DECOMP_IP_CPX__
endif

########################################################################
#                Headers that need to be installed                     #
########################################################################
# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'

includecoindir = $(includedir)/coin

#DecompAlgo
includecoin_HEADERS =        \
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

########################################################################
#                     Installing manifest (MSVC++)                     #
########################################################################

if COIN_CXX_IS_CL
install-exec-hook:
	test -s dip.exe.manifest && \
	mt -manifest dip.exe.manifest -outputresource:dip.exe;\
	cp dip.exe $(bindir)/dip.exe
endif

#######################################################################
# Create the Config.h file that has all public defines and install it #
#######################################################################

install-exec-local:
	$(install_sh_DATA) config_dip.h $(DESTDIR)$(includecoindir)/DecompConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/DecompConfig.h
