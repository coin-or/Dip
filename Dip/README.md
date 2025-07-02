# Dip 0.95

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Dip?sort=semver)](https://github.com/coin-or/Dip/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

DIP (Decomposition in Integer Programming) is a framework for
implementing a wide variety of decomposition-based algorithms for solving
integer programs, including full generic column generation in which the user
need only indicate whic constraints to relax and the framework takes care of
all other aspects of the column-generation. The problem can be specified using
[DiPPy](https://pypi.org/project/coinor.dippy/), an extension of the
Python-based modeling language [PuLP](https://github.com/coin-or/PuLP) (see
[examples](Dip/src/dippy/examples)). If desired, custom subroutines for solving
the subproblem, branching, heuristics, etc. can be implemented in pure Python
and called automatically from C++ during execution of the algorithm.

Dip makes it easy to implement and compare different variety of
decomposition-based algorithm while keeping as many details of the algorithmic
implem,entation constant as possible, allowing for rigorous empirical
comparisons.

Using [DiPPy](https://pypi.org/project/coinor.dippy/), one can even implement
a full column-generation algorithm within an Excel spreadsheet using [Solver
Studio](http://solverstudio.org). See screen shots below.

Dip is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/eclipse-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org).

The Dip development site is https://github.com/coin-or/Dip.

## CITE

Code: [![DOI](https://zenodo.org/badge/23778922.svg)](https://zenodo.org/badge/latestdoi/23778922)

Paper: http://dx.doi.org/10.1007/s10107-005-0606-3

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/Dip/actions/workflows/windows-ci.yml/badge.svg?branch=stable/0.95)](https://github.com/coin-or/Dip/actions/workflows/windows-ci.yml?query=branch%3Astable/0.95)

[![Linux and MacOS Builds](https://github.com/coin-or/Dip/actions/workflows/linux-ci.yml/badge.svg?branch=stable/0.95)](https://github.com/coin-or/Dip/actions/workflows/linux-ci.yml?query=branch%3Astable/0.95)

## DOWNLOAD

What follows is a quick start guide for obtaining or building
Dip on common platforms. More detailed information is
available [here](https://coin-or.github.io/user_introduction.html).

### Docker image

There is a Docker image that provides Dip, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/Dip/releases). Older binaries
are archived as part of Dip
[here](https://www.coin-or.org/download/binary/Dip).

Due to license incompatibilities, pre-compiled binaries lack some 
functionality. If binaries are not available for your platform for the latest 
version and you would like to request them to be built and posted, feel free 
to let us know on the mailing list. 

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Dip from the
 [releases](https://github.com/coin-or/Dip/releases) page,
 * Cloning this repository from [Github](https://github.com/coin-or/Dip), or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

### Dependencies

Dip has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are aeasy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

These quick start instructions assume you are in a bash shell. 

### Using `coinbrew`

To download and build Dip from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Dip@0.95
./coinbrew build Dip
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Dip
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

## Quickstart with DIP/DipPy

If you are on Linux or OS X, DIP can be installed from source using
[coinbrew](https://github.com/coin-or/coinbrew) (see below). In Windows, there
are pre-built binary wheels for DipPy. After build and install of DIP, if you
are installing DipPy, make sure to set `PKG_CONFIG_PATH` to point to the
directory where the `.pc` files are installed (by default, this is in the
`lib/pkgconfig` directory in the installation directory). You may also need to
set either `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (OS X) to point to
the directory where the libraries are installed. Then simply do
```
pip install coinor.dippy
```

After installation, you can try running the examples to see if trhings are
working. These should all work out of the box. For a listing of examples and
some brief instructions, do
```
python -m coinor.dippy.examples
```
To get help for an example, do, e.g.,
```
python -m coinor.dippy.example.cflp --help
```
or just run
```
python -m coinor.dippy.examples.cflp
```

## Scrren Shots
### Using [DiPPy](https://pypi.org/project/coinor.dippy/) in [Solver Studio](http://solverstudio.org)
![Solver Studio Pic 1](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/SolverStudioPic1.png)
### Choosing algorithm in [Solver Studio](http://solverstudio.org)
![Solver Studio Pic 3](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/SolverStudioPic3.png)
### Displaying search tree using [GrUMPy](https://github.com/coin-or/GrUMPy)
![GrUMPy Pic](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/GrUMPyPic.png)
### Editing DipPy Model with Eclipse and PyDev
![DipPy Pic](https://raw.githubusercontent.com/coin-or/Dip/master/Dip/images/DipPyPic2.png)


## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxydoc` 

in the build directory. If Dip was built via `coinbrew`, then the build
directory will be `./build/Dip/0.95` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Dip/Doxygen).

## Additional Documentation

 * C++ [examples](Dip/examples)
 * [DiPPy](https://pypi.org/project/coinor.dippy/) (Python) [examples](Dip/src/dippy/examples)
 * Slides providing a brief introduction are [here](https://coral.ise.lehigh.edu/~ted/files/talks/DecompIFORS14.pdf)
 * Slides providing a much longer introduction are [here](https://coral.ise.lehigh.edu/~ted/files/talks/DecompCSIRO11.pdf)
 * For the methodological basis, see
   * [This](http://dx.doi.org/10.1007/s10107-005-0606-3) paper (pre-print [here](http://coral.ie.lehigh.edu/~ted/files/papers/DECOMP04.pdf)); or
   * [This](http://coral.ie.lehigh.edu/~ted/files/papers/DECOMP04.pdf) more detailed paper; or
   * [This](http://coral.ie.lehigh.edu/~ted/files/papers/MatthewGalatiDissertation09.pdf) dissertation (much more detailed)
 * http://coral.ie.lehigh.edu/~ted/files/papers/MatthewGalatiDissertation09.pdf


## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [Dip Discussion forum](https://github.com/coin-or/Dip/discussions)
 * [Report a bug](https://github.com/coin-or/Dip/issues/new)
 * [Doxygen generated documentation](http://coin-or.github.io/Dip/Doxygen)

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
   subproblem solver is used and to allow no branching candidates when 
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
 * Eliminated unnecessary MILP parameter section and joined it to DECOMP,
   as well as making DECOMP the default parameter section name. 
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

