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
# Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs     #
# All Rights Reserved.                                                       #
#============================================================================#

#############################################################################
#                       Names and other basic things                        #
#############################################################################
# autoconf version used by COIN
AC_PREREQ(2.59)

# name of project, version number, contact info
AC_INIT([Dip],[trunk],[dip@lists.coin-or.org])

# copyright for configure script
AC_COPYRIGHT([
//===========================================================================//
// This file is part of the DIP Solver Framework.                            //
//                                                                           //
// DIP is distributed under the Eclipse Public License as part of the        //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)       //
//                                                                           //
// Conceptual Design: Matthew Galati, SAS Institute Inc.                     //
//                    Ted Ralphs, Lehigh University                          //
//                                                                           //
// Copyright (C) 2002-2019, Lehigh University, Matthew Galati, Ted Ralphs    //
// All Rights Reserved.                                                      //
//===========================================================================//
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/DecompApp.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################
# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(Dip)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Dip)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

#AC_LIBTOOL_LANG_CXX_CONFIG
AC_LANG_PUSH(C++)

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# Get the EXEXT variable for CYGWIN
AC_EXEEXT

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_CHECK_PACKAGE(CoinUtils, [coinutils], [DipLib])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHECK_PACKAGE(Osi, [osi], [DipLib])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not available.])
fi
AC_COIN_CHECK_PACKAGE(Alps, [alps], [DipLib])
if test $coin_has_alps != yes ; then
  AC_MSG_ERROR([Required package Alps not available.])
fi

AC_COIN_CHECK_PACKAGE(Clp,  [osi-clp clp],  [DipLib])
AC_COIN_CHECK_PACKAGE(Sample,   [coindatasample])
AC_COIN_CHECK_PACKAGE(Cgl,  [cgl],      [DipLib])
AC_COIN_CHECK_PACKAGE(Cbc,  [osi-cbc cbc],  [DipLib])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available 

AC_COIN_CHECK_PACKAGE(Cpx,  [osi-cplex],  [DipLib])
AC_COIN_CHECK_PACKAGE(SYMPHONY, [osi-sym], [DipLib])
AC_COIN_CHECK_PACKAGE(Grb, [osi-gurobi], [DipLib])
AC_COIN_CHECK_PACKAGE(Xpr, [osi-xpress], [DipLib])

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library

AC_COIN_CHECK_PACKAGE(ASL,  [coinasl])
AC_COIN_CHECK_USER_LIBRARY([HMetis],[HMETIS],[hmetis.h],[HMETIS_PartRecursive])

#############################################################################
#                          OPENMP stuff                                     #
#############################################################################
AC_ARG_ENABLE([openmp],
[AC_HELP_STRING([--enable-openmp],
                [compile in OpenMP features])],
[case "${enableval}" in
   yes) dip_openmp=true
   ;;
   no)  dip_openmp=false
   ;;
   *) AC_MSG_ERROR(Bad value ${enableval} for --enable-openmp)
   ;;
esac],
[dip_openmp=auto])

if test $dip_openmp = auto; then 
   AC_MSG_CHECKING([whether compiler supports OpenMP])
   case "$CXX" in
     clang*) 
       AC_MSG_RESULT([no])
       dip_openmp=false
       ;;
     cl* | *\cl* | */cl*)
       AC_MSG_RESULT([yes])
       dip_openmp=true
       ;;
     g++*)
        dip_openmp=true
   	fopenmp_version=4.2
   	gcc_required_version=`echo ${fopenmp_version} | sed -e 's/\./0/g'`
   	get_gcc_version=`$CXX -dumpversion`
   	gcc_version=`echo ${get_gcc_version} | sed -e 's/\./0/g'`
   	major_gcc_version=`echo ${gcc_version:0:3}`
   	#AC_MSG_RESULT([$get_gcc_version])
   	if test ${major_gcc_version} -lt ${gcc_required_version}; then
       	   AC_MSG_RESULT([no])
           dip_openmp=false
        else
       	   AC_MSG_RESULT([yes])
           dip_openmp=true
	fi
        ;;
   esac
fi

AM_CONDITIONAL(OPENMP, [test $dip_openmp = true])
AC_SUBST(DIP_OPENMP_FLAGS)
AC_SUBST(DIP_OPENMP_LIBS)
if test $dip_openmp = true; then
  if test x$coin_cxx_is_cl = xtrue; then
    DIP_OPENMP_FLAGS=-openmp 
    DIP_OPENMP_LIBS=
  else
    DIP_OPENMP_FLAGS=-fopenmp
    DIP_OPENMP_LIBS=-lgomp
  fi
else
  DIP_OPENMP_FLAGS= 
  DIP_OPENMP_LIBS=
fi

#############################################################################
#                       Python Stuff (for Dippy)                            #
#############################################################################

AM_PATH_PYTHON([2.6],,[:])
if test "$PYTHON" != "" ; then
  PY_PREFIX=`$PYTHON -c 'import sys ; print sys.prefix'`
  echo $PY_PREFIX
  PY_PLATFORM=`$PYTHON -c 'import sys ; print sys.platform'`
  echo $PY_PLATFORM
  if test "$PY_PLATFORM" != "win32" ; then
     AC_COIN_CHECK_HEADER([$PY_PREFIX/include/python$PYTHON_VERSION/Python.h],,
			[PYTHON=])
  else
     PYTHON=""
  fi
fi		       
AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != "" ])
if test "$PY_PLATFORM" != "win32" ; then
  PYTHON_LIBS="-lpython$PYTHON_VERSION"
  PYTHON_CFLAGS="-I$PY_PREFIX/include/python$PYTHON_VERSION"
else
  PYTHON_LIBS="-lpython"
  PYTHON_CFLAGS="-I$PY_PREFIX\\include"
fi
AC_SUBST([PYTHON_LIBS])
AC_SUBST([PYTHON_CFLAGS])

#############################################################################
#                      Check for the LP solver                              #
#############################################################################

AC_MSG_CHECKING([that we have an LP solver])

found_lp_solver=false
if test $coin_has_clp = yes; then
   found_lp_solver=true
elif test $coin_has_cpx = yes; then
   found_lp_solver=true
elif test $coin_has_grb = yes; then
   found_lp_solver=true
elif test $coin_has_xpr = yes; then
   found_lp_solver=true
fi

if test $found_lp_solver = false; then 
  AC_MSG_ERROR(Couldn't find an LP solver!)
else
  AC_MSG_RESULT(yes)
fi

#############################################################################
#                      Check for the IP solver                              #
#############################################################################

AC_MSG_CHECKING([that we have an IP solver])

found_ip_solver=false
if test $coin_has_clp = yes; then
   found_ip_solver=true
elif test $coin_has_cpx = yes; then
   found_ip_solver=true
elif test $coin_has_grb = yes; then
   found_ip_solver=true
elif test $coin_has_xpr = yes; then
   found_ip_solver=true
fi

if test $found_ip_solver = false; then 
  AC_MSG_ERROR(Couldn't find an IP solver!)
else
  AC_MSG_RESULT(yes)
fi

#############################################################################
#                  Export symbols for available solvers                     #
#############################################################################

if test $coin_has_clp = yes ; then
  AC_DEFINE(DIP_HAS_CLP, [1], [Define to 1 if we have CLP])
fi
if test $coin_has_symphony = yes ; then
  AC_DEFINE(DIP_HAS_SYMPHONY, [1], [Define to 1 if we have SYMPHONY])
fi
if test $coin_has_cbc = yes ; then
  AC_DEFINE(DIP_HAS_CBC, [1], [Define to 1 if we have CBC])
fi
if test $coin_has_cpx = yes ; then
  AC_DEFINE(DIP_HAS_CPX, [1], [Define to 1 if we have CPLEX])
fi
if test $coin_has_grb = yes ; then
  AC_DEFINE(DIP_HAS_GRB, [1], [Define to 1 if we have Gurobi])
fi
if test $coin_has_xpr = yes ; then
  AC_DEFINE(DIP_HAS_XPR, [1], [Define to 1 if we have Xpress])
fi

##############################################################################
#                   VPATH links for unit test                                #
##############################################################################

AC_CONFIG_LINKS(test/MILP_DecompApp.cpp:examples/MILP/MILP_DecompApp.cpp
		test/MILP_DecompApp.h:examples/MILP/MILP_DecompApp.h
		test/MILP_Main.cpp:examples/MILP/MILP_Main.cpp
		test/MILP_Param.h:examples/MILP/MILP_Param.h)
                
AC_COIN_LINK_INPUT_CMD

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi Clp Cgl Cbc Alps)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Location of DIP for later use
AC_SUBST(DIPINCDIR)
AC_SUBST(DIPSRCDIR)
DipINCDIR=`cd $srcdir; pwd`/src
DipSRCDIR=`cd $srcdir; pwd`/src

# To be used in applications: default DIP defs
AC_SUBST(DipDEFS)
AC_SUBST(ARCH)
AC_SUBST(COIN_HAS_CLP) 

CLLCFLAGS=
case "$CXX" in
  cl* | *\cl* | */cl*)
     if test $dip_openmp = true; then
         CXXFLAGS="$CXXFLAGS -openmp"
     fi
  ;;
  g++* | *\g++* | */g++*)
     if test $dip_openmp = true; then
         CXXFLAGS="$CXXFLAGS -fopenmp"
     fi

     #Ignore warnings about unrecognized pragma lines and                                                                                    
     #deprecated string conversions                                                                                                          
     CXXFLAGS="$CXXFLAGS -Wno-unknown-pragmas -Wno-write-strings"

     #Ignore warnings about comparisons of signed and unsigned values and                                                                    
     #casting const pointers to ordinary pointers                                                                                            
     CXXFLAGS="$CXXFLAGS -Wno-sign-compare -Wno-cast-qual"
   ;;
esac

AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/dippy/Makefile
                 examples/Makefile.Examples
		 examples/TSP/Makefile
		 examples/TSP/tsp.parm
                 examples/TSP/data_extract.sh
                 examples/TSP/get_boost.sh
                 examples/TSP/get_concorde.sh
                 examples/SmallIP/Makefile
                 examples/SmallIP/smallip.parm
                 examples/SmallIP/README.txt
                 examples/SmallIP/test.sh
                 examples/MILP/Makefile
                 examples/MILP/milp.parm
                 examples/MILP/data_extract.sh
                 examples/MILP/README.txt
                 examples/MILPBlock/Makefile
                 examples/MILPBlock/milpblock.parm
                 examples/MILPBlock/milpblock.parm2
                 examples/MILPBlock/data_extract.sh
                 examples/MILPBlock/README.txt
                 examples/MILPBlock/test.sh
                 examples/AP3/Makefile
                 examples/AP3/ap3.parm
                 examples/GAP/Makefile
                 examples/GAP/gap.parm
                 examples/GAP/data_extract.sh
		 examples/GAP/get_pisinger.sh
		 examples/GAP/README.txt
                 examples/MAD/Makefile
                 examples/MAD/mad.parm
                 examples/MCF/Makefile
                 examples/MCF/mcf.parm
                 examples/MCF/data_extract.sh
		 examples/MCF/README.txt
                 examples/MMKP/Makefile
                 examples/MMKP/mmkp.parm
                 examples/MMKP/data_extract.sh
                 examples/MMKP/get_pisinger.sh
		 examples/MMKP/README.txt
                 examples/VRP/Makefile
                 examples/VRP/vrp.parm                 
                 examples/VRP/data_extract.sh                
                 examples/VRP/get_boost.sh
                 examples/VRP/get_concorde.sh
                 examples/VRP/get_cvrpsep.sh
                 examples/ATM/Makefile
                 examples/ATM/atm.parm
                 examples/ATM/data_extract.sh
		 examples/SDPUC/Makefile
		 examples/SDPUC/sdpuc.parm
		 examples/SDPUC/data_extract.sh
		 examples/SDPUC/README.txt
		 examples/ExternalSolvers/Knapsack/Pisinger/Makefile
		 examples/ExternalSolvers/Knapsack/Pisinger/get.pisinger
		 examples/ExternalSolvers/Knapsack/Pisinger/patch.pisinger
		 examples/ExternalSolvers/Knapsack/Pisinger/mcknap.h
                 examples/ExternalSolvers/VRP/get.cvrpsep
                 examples/ExternalSolvers/Graph/get.boost
                 examples/ExternalSolvers/TSP/get.concorde
                 test/Makefile
                 dip.pc
		 dip-uninstalled.pc
		 doxydoc/doxygen.conf])


# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_dip.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
