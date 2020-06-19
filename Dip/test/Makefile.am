# Copyright (C) 2009 Lehigh University and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License

## $Id: Makefile.am 507 2006-06-03 03:20:16Z andreasw $

# Author:  Ted Ralphs           Lehigh    2009-10-18

AUTOMAKE_OPTIONS = foreign

TARGETS=

if COIN_HAS_CBC
TARGETS+=test_cbc
endif
if COIN_HAS_SYMPHONY
TARGETS+=test_symphony
endif
if COIN_HAS_CPX
TARGETS+=test_cpx
endif
if COIN_HAS_GRB
TARGETS+=test_grb
endif
if COIN_HAS_XPR
TARGETS+=test_xpr
endif

test: ../src/dip$(EXEEXT) $(TARGETS)

test_cbc:

	../src/dip$(EXEEXT) --BlockFileFormat Pair  --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.block --BestKnownLB 508.29976   --BestKnownUB 508.29976

	../src/dip$(EXEEXT) --BlockFileFormat List --Instance  `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.block --BestKnownLB 59704.02009   --BestKnownUB 59704.02009

	../src/dip$(EXEEXT) --BlockFileFormat Pair --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.block --BestKnownLB 11   --BestKnownUB 11

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.dec --BestKnownLB -88   --BestKnownUB -88

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --doCut 1 --doPriceCut 0 --BestKnownLB -88   --BestKnownUB -88

test_symphony:

	../src/dip$(EXEEXT) --BlockFileFormat Pair  --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.block --DecompIPSolver SYMPHONY --BestKnownLB 508.29976   --BestKnownUB 508.29976

	../src/dip$(EXEEXT) --BlockFileFormat List --Instance  `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.block --DecompIPSolver SYMPHONY --BestKnownLB 59704.02009   --BestKnownUB 59704.02009

	../src/dip$(EXEEXT) --BlockFileFormat Pair --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.block --DecompIPSolver SYMPHONY --BestKnownLB 11   --BestKnownUB 11

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.dec --DecompIPSolver SYMPHONY --BestKnownLB -88   --BestKnownUB -88

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --doCut 1 --doPriceCut 0 --DecompIPSolver SYMPHONY --BestKnownLB -88   --BestKnownUB -88

test_cpx:

	../src/dip$(EXEEXT) --BlockFileFormat Pair  --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.block --DecompIPSolver CPLEX --DecompLPSolver CPLEX --BestKnownLB 508.29976   --BestKnownUB 508.29976

	../src/dip$(EXEEXT) --BlockFileFormat List --Instance  `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.block --DecompIPSolver CPLEX --DecompLPSolver CPLEX --BestKnownLB 59704.02009   --BestKnownUB 59704.02009

	../src/dip$(EXEEXT) --BlockFileFormat Pair --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.block --DecompIPSolver CPLEX --DecompLPSolver CPLEX --BestKnownLB 11   --BestKnownUB 11

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.dec --DecompIPSolver CPLEX --DecompLPSolver CPLEX --BestKnownLB -88   --BestKnownUB -88

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --doCut 1 --doPriceCut 0 --DecompIPSolver CPLEX --DecompLPSolver CPLEX --BestKnownLB -88   --BestKnownUB -88

test_grb:

	../src/dip$(EXEEXT) --BlockFileFormat Pair  --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.block --DecompIPSolver Gurobi --DecompLPSolver Gurobi --BestKnownLB 508.29976   --BestKnownUB 508.29976

	../src/dip$(EXEEXT) --BlockFileFormat List --Instance  `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.block --DecompIPSolver Gurobi --DecompLPSolver Gurobi --BestKnownLB 59704.02009   --BestKnownUB 59704.02009

	../src/dip$(EXEEXT) --BlockFileFormat Pair --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.block --DecompIPSolver Gurobi --DecompLPSolver Gurobi --BestKnownLB 11   --BestKnownUB 11

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.dec --DecompIPSolver Gurobi --DecompLPSolver Gurobi --BestKnownLB -88   --BestKnownUB -88

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --doCut 1 --doPriceCut 0 --DecompIPSolver Gurobi --DecompLPSolver Gurobi --BestKnownLB -88   --BestKnownUB -88

test_xpr:

	../src/dip$(EXEEXT) --BlockFileFormat Pair  --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/retail3.block --DecompIPSolver Xpress --DecompLPSolver Xpress --BestKnownLB 508.29976   --BestKnownUB 508.29976

	../src/dip$(EXEEXT) --BlockFileFormat List --Instance  `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/atm_5_10_1.block --DecompIPSolver Xpress --DecompLPSolver Xpress --BestKnownLB 59704.02009   --BestKnownUB 59704.02009

	../src/dip$(EXEEXT) --BlockFileFormat Pair --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.mps --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/wedding_16.block --DecompIPSolver Xpress --DecompLPSolver Xpress --BestKnownLB 11   --BestKnownUB 11

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --BlockFile `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.dec --DecompIPSolver Xpress --DecompLPSolver Xpress --BestKnownLB -88   --BestKnownUB -88

	../src/dip$(EXEEXT) --Instance `$(CYGPATH_W) $(SAMPLE_DATA)`/block_milp.lp --doCut 1 --doPriceCut 0 --DecompIPSolver Xpress --DecompLPSolver Xpress --BestKnownLB -88   --BestKnownUB -88


.PHONY: test

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

#DISTCLEANFILES = flugpl.log flugpl.mps
