from __future__ import print_function
from builtins import range
from pulp import *

try:
    from src.dippy import Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import Solve
    from coinor.dippy.examples.gen_func import *

from .tsp_func import *
    
args = parseArgs()
    
prob = formulate(args.module)

prob.generate_cuts = generate_cuts
prob.is_solution_feasible = is_solution_feasible

dippyOpts = addDippyOpts(args)

Solve(prob, dippyOpts)

# print solution
for arc, var in list(prob.arc_vars.items()):
    if var.varValue:
        print(arc, var.varValue)
