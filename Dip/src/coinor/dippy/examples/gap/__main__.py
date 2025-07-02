#!/usr/bin/env python

# Generalized Assignment Problem
# argument should be a problem file, see Dip/examples/GAP_Instance.cpp for format
# for an e.g. see gap0512-2.dat included in this directory

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import sys

from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize

try:
    from src.dippy import Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import Solve
    from coinor.dippy.examples.gen_func import *

from .gap_func import *

args = parseArgs()

# parse data file
prob = formulate(args.module)

dippyOpts = addDippyOpts(args)

dippyOpts['TolZero'] = '%s' % tol

Solve(prob, dippyOpts)

for m in prob.MACHINES:
    print() 
    print("Machine %d assigned tasks" %m, end=' ')
    for t in prob.TASKS:
        v = prob.assignVars[m][t].varValue
        if v:
            print("%d" %t, end=' ')
print()
            
if prob.display_mode != 'off':
    if (prob.Tree.attr['display'] == 'pygame') or (prob.Tree.attr['display'] == 'xdot'):
        prob.Tree.display()

