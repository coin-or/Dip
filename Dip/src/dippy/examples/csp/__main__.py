from __future__ import print_function
from builtins import range
from pulp import *

try:
    from src.dippy import Solve
except ImportError:
    from coinor.dippy import Solve

from .cutting_stock_func import *
    
if len(sys.argv) > 1:
    prob = formulate(sys.argv[1])
else:
    prob = formulate('coinor.dippy.examples.cflp.facility_ex2')

prob.relaxed_solver = solve_subproblem

Solve(prob, {
    'generateCuts': '1', 
    'doPriceCut':'1', 
    'SolveRelaxAsIp': '1', \
        # use default IP to solve subproblems
})

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


