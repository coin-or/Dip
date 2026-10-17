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

AC_INIT([Dip],[devel],[https://github.com/coin-or/Dip/issues/new],[coin-or-dip],[https://github.com/coin-or/Dip])

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

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (used in the examples)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                              COIN-OR components                           #
#############################################################################

m4_foreach_w([myvar],[CoinUtils Osi Cgl Alps],
  [ AC_COIN_CHK_PKG(myvar,[DipLib])
    if test $coin_has_[]m4_tolower(myvar) != yes ; then
      AC_MSG_ERROR([Required package myvar is not available.])
    fi
  ])

# LP and IP solvers. Each found solver defines DIP_HAS_<solver>. Dip uses
# the COIN-OR solvers both directly and through their Osi interfaces, so
# check for both (one .pc file per check, as some pkg-config versions do not
# accept several modules in one argument).
m4_foreach_w([myvar],[Clp Cbc SYMPHONY],
  [ AC_COIN_CHK_PKG(myvar,[DipLib])
    AC_COIN_CHK_PKG(Osi[]myvar,[DipLib],
      m4_case(myvar,[Clp],[osi-clp],[Cbc],[osi-cbc],[SYMPHONY],[osi-sym]))
    if test $coin_has_[]m4_tolower(myvar) = yes &&
       test $coin_has_osi[]m4_tolower(myvar) != yes ; then
      AC_MSG_ERROR([Found myvar, but not its Osi interface.])
    fi
  ])
AC_COIN_CHK_PKG(Cpx,[DipLib],[osi-cplex])
AC_COIN_CHK_PKG(Grb,[DipLib],[osi-gurobi])
AC_COIN_CHK_PKG(Xpr,[DipLib],[osi-xpress])

AC_MSG_CHECKING([that we have an LP solver])
if test $coin_has_clp = yes || test $coin_has_cpx = yes ||
   test $coin_has_grb = yes || test $coin_has_xpr = yes ; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Couldn't find an LP solver (Clp, CPLEX, Gurobi, or Xpress)!])
fi

# Instances for the tests
AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)

#############################################################################
#                            Other dependencies                             #
#############################################################################

# hMETIS (for partitioning into blocks); only used if requested with
# --with-hmetis, --with-hmetis-lflags, or --with-hmetis-cflags
AC_COIN_CHK_LIBHDR(HMetis,[DipLib],[-lhmetis],[],[],[],[],[default_skip])

# OpenMP (for solving the blocks in parallel); --disable-openmp turns it off
AC_OPENMP
if test -n "$OPENMP_CXXFLAGS" ; then
  DIPLIB_CFLAGS="$DIPLIB_CFLAGS $OPENMP_CXXFLAGS"
  DIPLIB_LFLAGS="$DIPLIB_LFLAGS $OPENMP_CXXFLAGS"
fi

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Doxygen documentation
AC_COIN_DOXYGEN([CoinUtils Osi Clp Cgl Cbc Alps])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile
                 dip.pc
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
                 examples/ExternalSolvers/TSP/get.concorde])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([src/config.h src/config_dip.h])

# Finalize libs and cflags for building DipLib.
AC_COIN_FINALIZE_FLAGS([DipLib])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
