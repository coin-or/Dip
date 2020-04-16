from __future__ import print_function
from builtins import range
from pulp import *
import sys

try:
    from src.dippy import Solve
except ImportError:
    from coinor.dippy import Solve

from .tsp_func import *
    
if len(sys.argv) > 1:
    prob = formulate(sys.argv[1])
else:
    prob = formulate('coinor.dippy.examples.tsp.tsp_data')

prob.generate_cuts = generate_cuts
prob.is_solution_feasible = is_solution_feasible

Solve(prob, {'doCut': '1'})

# print solution
for arc, var in list(prob.arc_vars.items()):
    if var.varValue:
        print(arc, var.varValue)
