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
# Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs     #
# All Rights Reserved.                                                       #
#============================================================================#

#############################################################################
#                       Names and other basic things                        #
#############################################################################
# autoconf version used by COIN
AC_PREREQ(2.59)

# name of project, version number, contact info
AC_INIT([Dip],[0.8],[ted@lehigh.edu])

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
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, Ted Ralphs    //
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
AC_COIN_PROJECTDIR_INIT

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

# Check if user wants to use readline library
AC_COIN_CHECK_GNU_READLINE
AM_CONDITIONAL(COIN_HAS_READLINE,
               [test $coin_has_readline = yes])

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(Data)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Alps)
AC_COIN_HAS_PROJECT(Dip)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available 

AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library

AC_COIN_HAS_ASL

#############################################################################
#                      Check for the LP solver                              #
#############################################################################
AM_CONDITIONAL(COIN_HAS_CLP,
               false)
AM_CONDITIONAL(COIN_HAS_CPLEX,
               false)

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
     AC_COIN_HAS_PROJECT(Clp)
     if [test $coin_has_clp = unavailable ||
         test $coin_has_clp = skipping]; then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cpx = false; then 
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
test $coin_has_cpx = true])

#############################################################################
#                      Check for the IP solver                              #
#############################################################################
AM_CONDITIONAL(COIN_HAS_CBC,
               false)

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
     AC_COIN_HAS_PROJECT(Cbc)
     if [test $coin_has_cbc = unavailable ||
         test $coin_has_cbc = skipping]; then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cpx = false; then 
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
test $coin_has_cpx = true])

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
AC_SUBST(DipINCDIR)
AC_SUBST(DipSRCDIR)
DipINCDIR=`cd $srcdir; pwd`/src
DipSRCDIR=`cd $srcdir; pwd`/src

# To be used in applications: default DIP defs
AC_SUBST(DipDEFS)
AC_SUBST(ARCH)
AC_SUBST(COIN_HAS_CLP) 

#Set CFLAGS to CXXFLAGS since we are using g++
CFLAGS=$CXXFLAGS
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(CLLCFLAGS)

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 examples/Makefile.Examples
		 examples/TSP/Makefile
		 examples/TSP/tsp.parm
                 examples/TSP/data_extract.sh
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
		 examples/ExternalSolvers/Knapsack/Pisinger/Makefile
		 examples/ExternalSolvers/Knapsack/Pisinger/get.pisinger
		 examples/ExternalSolvers/Knapsack/Pisinger/patch.pisinger
		 examples/ExternalSolvers/Knapsack/Pisinger/mcknap.h
                 examples/ExternalSolvers/VRP/get.cvrpsep
                 examples/ExternalSolvers/Graph/get.boost
                 examples/ExternalSolvers/TSP/get.concorde
                 test/Makefile
                 dip_addlibs.txt])


# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_dip.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
