# DIP README

Welcome to the README for the Decomposition in Integer Programming (DIP)
Framework. DIP is distributed under the Eclipse Public License and is freely
redistributable. All source code and documentation is Copyright Matthew
Galati, Ted Ralphs, Lehigh University, and others. This README may be 
redistributed freely.

## Current Testing Status

[![Build Status](https://travis-ci.org/coin-or/Dip.svg?branch=master)](https://travis-ci.org/coin-or/Dip)

[![Build status](https://ci.appveyor.com/api/projects/status/eke5pr1qo9ywvs2h/branch/master?svg=true)](https://ci.appveyor.com/project/tkralphs/dip-qk8he/branch/master)

## Download and Install

[ ![Download](https://api.bintray.com/packages/coin-or/download/Dip/images/download.svg?version=0.92) ](https://bintray.com/coin-or/download/Dip/0.92/link)

Binaries for most platforms are available for download from [Bintray](https://bintray.com/coin-or/download/Dip).

For DipPy, see [Pypi](https://pypi.org/project/coinor.dippy/).

For installation, see the [INSTALL](https://github.com/coin-or/Dip/blob/master/INSTALL.md) file.

## Cite

[![DOI](https://zenodo.org/badge/23778922.svg)](https://zenodo.org/badge/latestdoi/23778922)

## DOCUMENTATION

For a quick start guide, please see the INSTALL file in this distribution.
Automatically generated documentation of the source code can be found here:

http://www.coin-or.org/Doxygen/Dip/

Further information can be found here:

http://github.com/coin-or/Dip

## SUPPORT

### List Serve

DIP users should either use the Dip mailing list or open an issue to ask questions to subscribe to the mailing list, 
go to http://list.coin-or.org/mailman/listinfo/dip

### Bug Reports

Bug reports should be reported on github at

https://github.com/Dip/issues/new

## SCREEN SHOTS

### Using DipPy in [Solver Studio](http://solverstudio.org)
![Solver Studio Pic 1](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/SolverStudioPic1.png)
### Choosing algorithm in [Solver Studio](http://solverstudio.org)
![Solver Studio Pic 3](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/SolverStudioPic3.png)
### Displaying search tree using [GrUMPy](https://github.com/coin-or/GrUMPy)
![GrUMPy Pic](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/GrUMPyPic.png)
### Editing DipPy Model with Eclipse and PyDev
![DipPy Pic](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/DipPyPic2.png)

## CHANGELOG

### Release 0.95.0
 * Python 3 support added

### Release 0.92.4
 * Update dependencies
 * Minor bug fixes for DipPy
 
### Release 0.92.3
 * Update dependencies
 * Minor bug fix
 * Add support for Appveyor and Travis
 * Install examples with DipPy

### Release 0.92.2
 * Get rid of pesky global variable DecompInf
 * Fix bugs in wedding planner example
 * Fix bugs in DipPy to allow returning no solutions, even when an exact 
   subprobem solver is used and to allow no branching candidates when 
   branching. 

### Release 0.91.6
 * Fixed bugs in Wedding Planner example.
 * Fixed bug in DipPy having to do with branching.

### Release 0.92.1
 * Fixed problem with dependency linking

### Release 0.91.5
 * Fixed bug with SYMPHONY when not all solutions are accepted.

### Release 0.92.0
 * Substantially re-designed internals
 * Renamed classes, functions, and parameters more intuitively
 * Eliminated unnecessary "MILP" parameter section and joined it to "DECOMP,"
   as well as making "DECOMP" the default parameter section name. 
 * Changed parameter setting mechanism to make it possible to pass parameters
   directly to solvers using native names.
 * Added interface to Gurobi
 * Added ability to select solver at run-time rather than compile-time.
 * In DipPy, the user can now return a status in the subproblem solve to
   indicate whether the subproblem was solved exactly. Previously, DipPy
   solved the subproblem to optimality internally whenever no solution was
   returned, which is unnecessary if the user's subproblem solver is exact. It
   also means that the user was previously required to provide a full
   description of the subproblem.

### Release 0.91.4
 * Fixed bugs in examples
 * Updates to dependencies
 * Samll bug fixes

### Release 0.91.4
 * Fixed bugs in examples
 * Updates to dependencies
 * Samll bug fixes

### Release 0.91.3
 * Fixes for correctly producing Doxygen documentation

### Release 0.91.2
 * Fixed issue with master only variable when solving master as an integer 
   program.
 * Added ability to generate multiple columns per iteration with SYMPHONY and 
   Cbc. 

### Release 0.91.1
 * Updating dependencies.
 * Fix for dependency linking
 * Fix to installation with {{{DESTDIR}}}

### Release 0.91.0
 * Multiple parallel modes added
   * Solution of individual subproblems can be parallelized
   * Multiple subproblems can be solved simultaneously
   * Multiple algorithms can be executed in parallel (branch and cut plus
        one or more decomposition-based algorithms)
 * Warm starting for solutions of subproblems is supported when using
   SYMPHONY as the subproblem solver. 
 * Unbounded feasible regions now supported.
 * Explicit treatment of master-only variables.

### Release 0.9.12:

 * Fixed long-standing issues with stand-alone apps

 * Small some bug fixes 

### Release 0.9.11:

 * Added some new DipPy examples

 * Small bug fixes 

### Release 0.9.10:

 * Fixes to stand-alone app examples

 * Fixes to Visual studio files for examples to support property pages 

##New Stable Version 0.91:

 * Changes to the DipPy callback interface to make it more user friendly

 * Changes to the application interface

 * Planning for other changes to the internal algorithm 

### Release 0.9.9:

 * Fixes to DipPy build and examples

 * Fixes to allow CGL cuts to be generated from within DipPy branch and price. 

### Release 0.9.8:

 * Fixes to DipPy examples

 * Support for dependency linking 

### Release 0.9.7:

 * Fixes to DipPy examples

 * Support for dependency linking 

### Release 0.9.6:

 * Fixes to allow proper installation of DipPy on Mac OS X 

### Release 0.9.5:

 * Small fixes to DipPy

### Release 0.9.4:

 * Fixes to parallel subproblem solution mode with OpenMP

### Release 0.9.3:

 * More updates to build system

### Release 0.9.2:

 * Updates to build system

### Release 0.9.1:

 * Fixes to Python installation

 * Fix to DipPy

### Release 0.9.0:

 * DIP now includes DipPy, a Python-based modeling language.

 * DIP is now a complete generic MILP solver, capable of automatically
   detecting block structure and applying a decomposition method.

 * DIP and DipPy build out of the box in Windows, OSX, and Linux with a
   combination of the autotools and a Python setup script.

 * There is now support for solving the subproblems in parallel when there is
   block structure using OpenMP.

 * Numerous bug fixes and improvements. 

### Release 0.82.2:

 * Fixes to build system

 * Updates to dependencies

### Release 0.82.1:

 * Fixes to build system

 * Updates to dependencies

### Release 0.82.0:

 * Support for MSVC++ version 10 added.

 * Support for BuildTools version 0.7 to incorporate recent enhancements,
   including proper library versioning in Linux, prohibiting installation of
   private headers, etc.

 * Updated externals to new stable versions of dependent projects.

 * Minor bug fixes. 

