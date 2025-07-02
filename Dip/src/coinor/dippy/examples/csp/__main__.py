from __future__ import print_function
from builtins import range
from pulp import *
try:
    from src.dippy import Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import Solve
    from coinor.dippy.examples.gen_func import *

from .cutting_stock_func import *
    
args = parseArgs()

prob = formulate(args.module)

if args.useCustomSolver:
    prob.relaxed_solver = solve_subproblem

dippyOpts = addDippyOpts(args)
    
Solve(prob, dippyOpts)

for i, var in list(prob.useVars.items()):
    if var.varValue:
        print("Use", i, var.varValue)

for pat in prob.PATTERNS:
    for i in prob.ITEMS:
        if prob.cutVars[(pat, i)].varValue:
            print("Pat", pat, "item", i, \
                  prob.cutVars[(pat, i)].varValue)

##for (pat, w), var in cutVars.items():
##    if var.varValue:
##        print "Pat", pat, "item", w, var.varValue


