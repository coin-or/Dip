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
from os.path import dirname
from inspect import getfile
import coinor.dippy.examples.gap

from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize

try:
    from src.dippy import Solve
except ImportError:
    from coinor.dippy import Solve

from .gap_func import *
    
if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or len(sys.argv) > 2:
        print('Usage: coke <module_name>')
        print('       module_name : Python module containing instance data')
        print('                     For example file, check directory')
        print('                    ', dirname(getfile(coinor.dippy.examples.gap)))
        exit()
    else:
        module_name = sys.argv[1]
else:
    module_name = 'coinor.dippy.examples.coke.gap0515-2'

# parse data file
prob = formulate(module_name)

Solve(prob, {
    'TolZero': '%s' % tol,
    'doPriceCut': '1',
#    'logLevel': '3', 
})

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

