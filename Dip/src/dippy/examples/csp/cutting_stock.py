from __future__ import print_function
from builtins import range
#!/usr/bin/env python

from pulp import *

try:
    import src.dippy as dippy
    from src.dippy import DipSolStatOptimal
except ImportError:
    import coinor.dippy as dippy
    from coinor.dippy import DipSolStatOptimal

length = {
"9cm": 9,
"7cm": 7,
"5cm": 5
}

ITEMS = list(length.keys())

demand = {
"9cm": 3,
"7cm": 2,
"5cm": 2
}

total_patterns = sum(demand[i] for i in ITEMS)

total_length = {}
for p in range(total_patterns):
    total_length[p] = 20
PATTERNS = list(total_length.keys())

def cross(i1, i2):
    r = []
    for a in i1:
        for b in i2:
            r.append((a, b))
    return r

CUTS = cross(PATTERNS, ITEMS)

prob = dippy.DipProblem("Python", LpMinimize)

# create variables

useVars = LpVariable.dicts("Use", PATTERNS, 0, 1, LpBinary)
prob.useVars = useVars

cutVars = LpVariable.dicts("Cut", CUTS, 0, 10, LpInteger)
prob.cutVars = cutVars

# objective
prob += lpSum(useVars[p] for p in PATTERNS), "min"

# Meet demand
for i in ITEMS:
    prob += lpSum(cutVars[(p, i)] for p in PATTERNS) \
            >= demand[i]

# Ordering patterns
for i, p in enumerate(PATTERNS):
    if p != PATTERNS[-1]:
        prob += useVars[p] >= useVars[PATTERNS[i+1]]

for p in PATTERNS:
    prob.relaxation[p] += \
    lpSum(length[i] * cutVars[(p, i)] for i in ITEMS) \
    <= total_length[p] * useVars[p]

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

prob.relaxed_solver = solve_subproblem

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

dippy.Solve(prob, {
    'generateCuts': '1', 
    'doPriceCut':'1', 
    'SolveRelaxAsIp': '1', \
        # use default IP to solve subproblems
})

for i, var in list(useVars.items()):
    if var.varValue:
        print("Use", i, var.varValue)

for pat in PATTERNS:
    for i in ITEMS:
        if cutVars[(pat, i)].varValue:
            print("Pat", pat, "item", i, \
                  cutVars[(pat, i)].varValue)

##for (pat, w), var in cutVars.items():
##    if var.varValue:
##        print "Pat", pat, "item", w, var.varValue


