#============================================================================#
# This file is part of the DECOMP Solver Framework.                          #
#                                                                            #
# DECOMP is distributed under the Common Public License as part of the       #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Author: Matthew Galati, SAS Institute Inc. (matthew.galati@sas.com)        #
#                                                                            #
# Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs.#
# All Rights Reserved.                                                       #
#============================================================================#

#############################################################################
#                       Names and other basic things                        #
#############################################################################
# autoconf version used by COIN
AC_PREREQ(2.59)

# name of project, version number, contact info
AC_INIT([Decomp],[0.8],[matthew.galati@sas.com])

# copyright for configure script
AC_COPYRIGHT([
//===========================================================================//
// This file is part of the Decomp Solver Framework.                         //
//                                                                           //
// Decomp is distributed under the Common Public License as part of the      //
// COIN-OR repository (http://www.coin-or.org).                              //
//                                                                           //
// Author: Matthew Galati, Lehigh University                                 //
//                                                                           //
// Copyright (C) 2002-2009, Lehigh University, Matthew Galati, and Ted Ralphs//
// All Rights Reserved.                                                      //
//===========================================================================//
])

# the version number of configure.in
AC_REVISION(0.8)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(configure.ac)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#     We only need automake to generate Makefiles for the distribution      #
#############################################################################

# Initialize automake
AC_COIN_CREATE_LIBTOOL

#############################################################################
#             Check which ThirdParty sources are available                  #
#############################################################################

AC_COIN_THIRDPARTY_SUBDIRS([Blas],[ThirdParty/Blas],[daxpy.f])
AC_COIN_THIRDPARTY_SUBDIRS([Lapack],[ThirdParty/Lapack],[LAPACK/SRC/dlarf.f])

#############################################################################
#                  Check which subprojects are there                        #
#############################################################################
AC_COIN_MAIN_SUBDIRS(CoinUtils)
AC_COIN_MAIN_SUBDIRS(Clp)
AC_COIN_MAIN_SUBDIRS(Osi)
AC_COIN_MAIN_SUBDIRS(Cgl)
AC_COIN_MAIN_SUBDIRS(Cbc)
AC_COIN_MAIN_SUBDIRS(Alps)
AC_COIN_MAIN_SUBDIRS(Decomp)
AC_COIN_MAIN_SUBDIRS(Data/Sample)

#############################################################################
#                            Third party solvers                            #
#############################################################################
# Check which third party solvers are available 
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################
# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
