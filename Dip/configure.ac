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
# Copyright (C) 2002-2013, Lehigh University, Matthew Galati, Ted Ralphs     #
# All Rights Reserved.                                                       #
#============================================================================#

#############################################################################
#                       Names and other basic things                        #
#############################################################################
# autoconf version used by COIN
AC_PREREQ(2.59)

# name of project, version number, contact info
AC_INIT([Dip],[0.9.8],[dip@lists.coin-or.org])

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
// Copyright (C) 2002-2011, Lehigh University, Matthew Galati, Ted Ralphs    //
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
AC_COIN_PROJECTDIR_INIT(Dip,2:8:1)

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

AC_COIN_CHECK_PACKAGE(CoinDepend, [coinutils osi alps], [DipLib])
if test $coin_has_coindepend != yes ; then
  AC_MSG_ERROR([Required package CoinUtils, Osi, or Alps not available.])
fi

AC_COIN_CHECK_PACKAGE(Clp,  [osi-clp],  [DipLib])
AC_COIN_CHECK_PACKAGE(Sample,   [coindatasample])
AC_COIN_CHECK_PACKAGE(Cgl,  [cgl],      [DipLib])
AC_COIN_CHECK_PACKAGE(Cbc,  [osi-cbc],  [DipLib])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available 

AC_COIN_CHECK_PACKAGE(Cpx,  [osi-cplex],  [DipLib])

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
[dip_openmp=false])

if test $dip_openmp = true; then 
   AC_MSG_CHECKING([G++ version])
   fopenmp_version=4.2
   gcc_required_version=`echo ${fopenmp_version} | sed -e 's/\./0/g'`
   get_gcc_version=`$CXX -dumpversion`
   gcc_version=`echo ${get_gcc_version} | sed -e 's/\./0/g'`
   major_gcc_version=`echo ${gcc_version:0:3}`
   AC_MSG_RESULT([$get_gcc_version])
   if test ${major_gcc_version} -lt ${gcc_required_version}; then
       AC_MSG_ERROR(OpenMP feautures are supported by G++ 4.2.0 or higher.)
   DipDEFS="$DipDEFS -D_OPENMP"
   fi
fi      

AM_CONDITIONAL(OPENMP, [test $dip_openmp = true])



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

AC_MSG_CHECKING([for Dip default lp solver])
AC_ARG_WITH([lp-solver],
[AC_HELP_STRING([--with-lp-solver@<:@=lpsolver@:>@],
                [specify the LP solver in small
                letters(default lpsolver=clp)])],
		[dip_lp_solver=$withval],[dip_lp_solver=clp])

found_package=true
case $dip_lp_solver in 
   clp)
     AC_MSG_RESULT(Clp)
     if test $coin_has_clp != yes; then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cpx != yes; then 
         found_package=false
     fi
   ;; 
  *)
    AC_MSG_ERROR(Unknown LP solver!)
  ;;
esac

if test $found_package = false; then 
  AC_MSG_ERROR(Couldn't find the $dip_lp_solver package!)
fi

AM_CONDITIONAL(CLP_LP_SOLVER,[test $dip_lp_solver = clp &&
test $coin_has_clp != unavailable &&
test $coin_has_clp != skipping])
AM_CONDITIONAL(CPLEX_LP_SOLVER,[test $dip_lp_solver = cplex &&
test $coin_has_cpx != unavailable &&
test $coin_has_cpx != skipping])

#############################################################################
#                      Check for the IP solver                              #
#############################################################################

AC_MSG_CHECKING([for DIP default ip solver])
AC_ARG_WITH([ip-solver],
[AC_HELP_STRING([--with-ip-solver@<:@=ipsolver@:>@],
                [specify the IP solver in small
                letters(default ipsolver=cbc)])],
		[dip_ip_solver=$withval],[dip_ip_solver=cbc])

found_package=true
case $dip_ip_solver in 
   cbc)
     AC_MSG_RESULT(Cbc) 
     if test $coin_has_cbc != yes; then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cpx != yes; then
         found_package=false
     fi
   ;; 
  *)
    AC_MSG_ERROR(Unknown IP solver!)
  ;;
esac

if test $found_package = false; then 
  AC_MSG_ERROR(Couldn't find the $dip_ip_solver package!)
fi

AM_CONDITIONAL(CBC_IP_SOLVER,[test $dip_ip_solver = cbc &&
test $coin_has_cbc != unavailable &&
test $coin_has_cbc != skipping])
AM_CONDITIONAL(CPLEX_IP_SOLVER,[test $dip_ip_solver = cplex &&
test $coin_has_cpx != unavailable &&
test $coin_has_cpx != skipping])

##############################################################################
#                   VPATH links for unit test                                #
##############################################################################

AC_CONFIG_LINKS(test/MILP_DecompApp.cpp:examples/MILP/MILP_DecompApp.cpp
		test/MILP_DecompApp.h:examples/MILP/MILP_DecompApp.h
		test/MILP_Main.cpp:examples/MILP/MILP_Main.cpp
		test/MILP_Param.h:examples/MILP/MILP_Param.h)
                
AC_COIN_LINK_INPUT_CMD

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
		 dip-uninstalled.pc])


# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_dip.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
