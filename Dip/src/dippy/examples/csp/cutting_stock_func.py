from __future__ import print_function
from builtins import range
from pulp import *
import importlib as ilib

try:
    from src.dippy import DipProblem, DipSolStatOptimal
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem, DipSolStatOptimal
    from coinor.dippy.examples.gen_func import *

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Solve a cutting stock problem.')
    parser.add_argument('--module', '-m', metavar = 'module name', 
                        help='name of the Python module from which to import data',
                        default = 'coinor.dippy.examples.csp.cutting_stock_data')
    parser.add_argument('--useCustomSolver', action='store_true', 
                        help='enable custom subproblem solver')

    addDippyArgs(parser)
    
    args = parser.parse_args()

    return(args)

def cross(i1, i2):
    r = []
    for a in i1:
        for b in i2:
            r.append((a, b))
    return r

def formulate(module_name):

    prob = DipProblem("Python", LpMinimize)

    m = ilib.import_module(module_name)
    
    # create variables
    useVars = LpVariable.dicts("Use", m.PATTERNS, 0, 1, LpBinary)
    prob.useVars = useVars

    cutVars = LpVariable.dicts("Cut", m.CUTS, 0, 10, LpInteger)
    prob.cutVars = cutVars

    # objective
    prob += lpSum(useVars[p] for p in m.PATTERNS), "min"

    # Meet demand
    for i in m.ITEMS:
        prob += lpSum(cutVars[(p, i)] for p in m.PATTERNS) \
                >= m.demand[i]

    # Ordering patterns
    for i, p in enumerate(m.PATTERNS):
        if p != m.PATTERNS[-1]:
            prob += useVars[p] >= useVars[m.PATTERNS[i+1]]

    for p in m.PATTERNS:
        prob.relaxation[p] += \
                lpSum(m.length[i] * cutVars[(p, i)] for i in m.ITEMS) \
                <= m.total_length[p] * useVars[p]

    prob.ITEMS = m.ITEMS
    prob.PATTERNS = m.PATTERNS
    
    return prob
        
def solve_subproblem(prob, keySub, redCosts, target):
    # get items with negative reduced cost
    item_idx = [i for i in ITEMS \
                if redCosts[cutVars[(keySub, i)]] < 0]
    vars = [cutVars[(keySub, i)] for i in item_idx]
    obj = [-redCosts[cutVars[(keySub, i)]] for i in item_idx]
    weights = [length[i] for i in item_idx]

    z, solution = kp(obj, weights, total_length[p])
    
    total_weight = sum(w * solution[i] \
                       for i, w in enumerate(weights))
    assert total_weight <= total_length[p]

    # add in reduced cost of useVars
    var_values = [(v, solution[i]) \
                  for i, v in enumerate(vars) \
                  if solution[i] > 0]
    var_values.append((useVars[keySub], 1))

    return DipSolStatOptimal, [var_values]

def kp(obj, weights, capacity):
    assert len(obj) == len(weights)
    n = len(obj)

    if n == 0:
        return 0, []

    if capacity == 0:
        return 0, [0 for i in range(n)]
    
    n = len(obj)

    # Don't include item
    zbest, solbest = kp(obj, weights, capacity - 1)
    # Check all items for inclusion
    for i in range(n):
        if weights[i] <= capacity:
            zyes, solyes = kp(obj, weights, \
                              capacity - weights[i])
            zyes += obj[i]
            solyes[i] += 1
            if zbest > zyes:
                zbest = zyes
                solbest = solyes

    return zbest, solbest

