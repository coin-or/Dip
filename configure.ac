#============================================================================#
# This file is part of the DIP Solver Framework.                             #
#                                                                            #
# DIP is distributed under the Eclipse Public License as part of the         #
# COIN-OR repository (http://www.coin-or.org).                               #
#                                                                            #
# Authors: Matthew Galati, SAS Institute Inc.                                #
#          Ted Ralphs, Lehigh University                                     #
#          Jiadong Wang, Lehigh University                                   #
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
#                  Check which subprojects are there                        #
#############################################################################

AC_COIN_MAIN_PACKAGEDIR(Blas,  ThirdParty,[daxpy.f])
AC_COIN_MAIN_PACKAGEDIR(Lapack,ThirdParty,[LAPACK/SRC/dlarf.f])
AC_COIN_MAIN_PACKAGEDIR(Sample,Data)
AC_COIN_MAIN_PACKAGEDIR(CoinUtils)
AC_COIN_MAIN_PACKAGEDIR(Osi)
AC_COIN_MAIN_PACKAGEDIR(Clp)
AC_COIN_MAIN_PACKAGEDIR(Cgl)
AC_COIN_MAIN_PACKAGEDIR(Cbc)
AC_COIN_MAIN_PACKAGEDIR(Alps)
AC_COIN_MAIN_PACKAGEDIR(SYMPHONY)
AC_COIN_MAIN_PACKAGEDIR(Dip)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi Clp Cgl Cbc Alps)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################
# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
