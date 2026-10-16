from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from past.utils import old_div

from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize, LpAffineExpression

try:
    from src.dippy import DipSolStatOptimal, Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipSolStatOptimal, Solve
    from coinor.dippy.examples.gen_func import *

from .facility_location_func import *

args = parseArgs()

prob = formulate(args)

if debug_print_lp:
    prob.writeLP('facility_main.lp')
    for n, i in enumerate(LOCATIONS):
        prob.writeRelaxed(n, 'facility_relax%s.lp' % i);

if args.useCustomSolver:
    prob.relaxed_solver = solve_subproblem
if args.customInitRule == 'initOneEach':
    prob.init_vars = init_one_each
if args.customInitRule == 'initFirstFit':
    prob.init_vars = init_first_fit
if args.useCustomCuts:
    prob.generate_cuts = generate_weight_cuts
if args.rootHeuristic:
    prob.heuristics = heuristics
    prob.root_heuristic = True
if args.nodeHeuristic:
    prob.heuristics = heuristics
    prob.node_heuristic = True

dippyOpts = addDippyOpts(args)
    
dippyOpts['TolZero'] = '%s' % tol

Solve(prob, dippyOpts)

if prob.display_mode != 'off':
    numNodes = len(prob.Tree.get_node_list())
    if prob.Tree.attr['display'] == 'svg':
        prob.Tree.write_as_svg(filename = "facility_node%d" % (numNodes + 1), 
                               prevfile = "facility_node%d" % numNodes)
    prob.Tree.display()

# print solution
print("Optimal solution found!") 
print("************************************")
for i in prob.LOCATIONS:
    if prob.use_vars[i].varValue > 0:
        print("Location ", i, " is assigned: ", end=' ')
        print([j for j in prob.PRODUCTS if prob.assign_vars[(i, j)].varValue > 0])
print("************************************")
print()
