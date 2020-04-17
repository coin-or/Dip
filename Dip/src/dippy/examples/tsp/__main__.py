from __future__ import print_function
from builtins import range
from pulp import *
import sys
from os.path import dirname
from inspect import getfile
import coinor.dippy.examples.tsp

try:
    from src.dippy import Solve
except ImportError:
    from coinor.dippy import Solve

from .tsp_func import *
    
if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help' or len(sys.argv) > 2:
        print('Usage: tsp <module_name>')
        print('       module_name : Python module containing instance data')
        print('                     For example file, check directory')
        print('                    ', dirname(getfile(coinor.dippy.examples.tsp)))
        exit()
    else:
        module_name = sys.argv[1]
else:
    module_name = 'coinor.dippy.examples.tsp.tsp_data'

prob = formulate(module_name)

prob.generate_cuts = generate_cuts
prob.is_solution_feasible = is_solution_feasible

Solve(prob, {'doCut': '1'})

# print solution
for arc, var in list(prob.arc_vars.items()):
    if var.varValue:
        print(arc, var.varValue)
