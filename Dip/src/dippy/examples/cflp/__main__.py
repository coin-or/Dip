#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from past.utils import old_div
import sys

from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize, LpAffineExpression

try:
    from src.dippy import DipSolStatOptimal, Solve
except ImportError:
    from coinor.dippy import DipSolStatOptimal, Solve

from .facility_location import *
    
tol = pow(pow(2, -24), old_div(2.0, 3.0))

#display_mode = 'xdot'
#layout = 'dot'

if len(sys.argv) > 1:
    prob = formulate(sys.argv[1])
else:
    prob = formulate('coinor.dippy.examples.cflp.facility_ex2')

if debug_print_lp:
    prob.writeLP('facility_main.lp')
    for n, i in enumerate(LOCATIONS):
        prob.writeRelaxed(n, 'facility_relax%s.lp' % i);

#prob.writeFull('facility.lp', 'facility.dec')

#prob.relaxed_solver = solve_subproblem
#prob.init_vars = init_one_each
#prob.init_vars = init_first_fit
#prob.generate_cuts = generate_weight_cuts
#prob.heuristics = heuristics
#prob.root_heuristic = True
#prob.node_heuristic = True

dippyOpts = {}
algo = 'Cut'
if len(sys.argv) > 2:
    algo = sys.argv[2]
if algo == 'PriceCut':
    dippyOpts['doPriceCut'] = '1'
    dippyOpts['CutCGL'] = '1'
elif algo == 'Price':
    dippyOpts['doPriceCut'] = '1'
    dippyOpts['CutCGL'] = '0'
elif algo == 'Direct':
    dippyOpts['doDirect'] = '1'
    dippyOpts['doCut'] = '1'
else:
    dippyOpts['doCut'] = '1'
    
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
