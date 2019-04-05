#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMinimize, LpAffineExpression

try:
    import src.dippy as dippy
except ImportError:
    import coinor.dippy as dippy

CC = 1.3
BIG_M = 1e10

MINE_SUPPLY = {
    "M1": 25.8,
    "M2": 728,
    "M3": 1456,
    "M4": 49,
    "M5": 36.9,
    "M6": 1100,
}
MINES = list(MINE_SUPPLY.keys())
MINES.sort()

LOCATIONS = ["L1", "L2", "L3", "L4", "L5", "L6"]

SIZE_COSTS = {
    0: 0,
    75: 4.4,
    150: 7.4,
    225: 10.5,
    300: 13.5,
    375: 16.5,
    450: 19.6,
}
SIZES = list(SIZE_COSTS.keys())
SIZES.sort()

CUSTOMER_DEMAND = {
    "C1": 83,
    "C2": 5.5,
    "C3": 6.975,
    "C4": 5.5,
    "C5": 720.75,
    "C6": 5.5,
}
CUSTOMERS = list(CUSTOMER_DEMAND.keys())
CUSTOMERS.sort()

MINE_TRANS_DATA = """
        L1      L2      L3      L4      L5      L6 
M1      231737  46813   79337   195845  103445  45186
M2      179622  267996  117602  200298  128184  49046
M3      45170   93159   156241  218655  103802  119616
M4      149925  254305  76423   123534  151784  104081
M5      152301  205126  24321   66187   195559  88979
M6      223934  132391  51004   122329  222927  54357
"""

CUST_TRANS_DATA = """
        L1      L2      L3      L4      L5      L6 
C1      6736    42658   70414   45170   184679  111569
C2      217266  227190  249640  203029  153531  117487
C3      35936   28768   126316  2498    130317  74034
C4      73446   52077   108368  75011   49827   62850
C5      174664  177461  151589  153300  59916   135162
C6      186302  189099  147026  164938  149836  286307
"""

def read_table(data, coerce, transpose=False):
    lines = data.splitlines()
    headings = lines[1].split()
    result = {}
    for row in lines[2:]:
        items = row.split()
        for i, item in enumerate(items[1:]):
            if transpose: key = (headings[i], items[0])
            else: key = (items[0], headings[i])
            result[key] = coerce(item)
    return result

MINE_TRANS = read_table(MINE_TRANS_DATA, int)
for key in MINE_TRANS:
    MINE_TRANS[key] = MINE_TRANS[key] * CC

CUST_TRANS = read_table(CUST_TRANS_DATA, int, \
                        transpose=True)

ARC_COSTS = dict(MINE_TRANS)
ARC_COSTS.update(CUST_TRANS)
ARCS = list(ARC_COSTS.keys())

def cross(i1, i2):
    r = []
    for a in i1:
        for b in i2:
            r.append((a, b))
    return r

LOC_SIZES = cross(LOCATIONS, SIZES)

prob = dippy.DipProblem("Coke", LpMinimize)

# create variables
buildVars = LpVariable.dicts("Build", LOC_SIZES, None, \
                             None, LpBinary)
prob.buildVars = buildVars

# create arcs
flowVars = LpVariable.dicts("Arcs", ARCS)
for a in ARCS:
    flowVars[a].bounds(0, BIG_M)

prob.SIZES = SIZES

# objective
prob += 1e6 * lpSum(buildVars[(l, s)] * SIZE_COSTS[s] \
                    for (l, s) in LOC_SIZES) + \
              lpSum(flowVars[(s, d)] * ARC_COSTS[(s, d)] \
                    for (s, d) in ARCS), "min"

# plant availability
for loc in LOCATIONS:
    prob += lpSum(flowVars[(loc, i)] for i in CUSTOMERS) \
            <= lpSum(buildVars[(loc, s)] *s for s in SIZES)

# one size
for loc in LOCATIONS:
    prob += lpSum(buildVars[(loc, s)] for s in SIZES) == 1

# conserve flow (mines)
# flows are in terms of tonnes of coke
for m in MINES:
    prob += lpSum(flowVars[(m, j)] for j in LOCATIONS) <= \
            old_div(MINE_SUPPLY[m],CC)

for loc in LOCATIONS:
    prob += lpSum(flowVars[(m, loc)] for m in MINES) - \
            lpSum(flowVars[(loc, c)] for c in CUSTOMERS) \
            >= 0

for c in CUSTOMERS:
    prob += lpSum(flowVars[(loc, c)] \
                  for loc in LOCATIONS) \
            >= CUSTOMER_DEMAND[c]

def do_branch(prob, sol):
    tol = 1e-10
    SIZES = prob.SIZES
    buildVars = prob.buildVars
    for loc in LOCATIONS:
        sol_size = sum(sol[buildVars[loc, size]] * \
                       size for size in SIZES)
        # smallest index and size larger than sol_size
        if abs(sol_size - SIZES[-1]) < tol:
            continue
        i, bigger = [(i, s) for i, s in enumerate(SIZES) \
                     if s > sol_size][0]
        if i == 0:
            smaller = 0
        else:
            smaller = SIZES[i-1]
        if bigger - sol_size > tol and sol_size - smaller > tol:
            down_branch_ub = dict([(buildVars[loc, SIZES[j]], 0) \
                              for j in range(i, len(SIZES))])
            up_branch_ub = dict([(buildVars[loc, SIZES[j]], 0) \
                            for j in range(0, i)])

            return ({}, down_branch_ub, {}, up_branch_ub)
        
prob.branch_method = do_branch

dippy.Solve(prob, {
    'CutCGL': '0',
})

def print_table(rows, cols, fn):
    print("\t", "\t".join(cols))
    for r in rows:
        print(r,"\t", "\t".join(str(fn(r,c)) for c in cols))

def print_var_table(rows, cols, var, fn=lambda x: x):
    print_table(rows, cols, lambda x, y:
                fn(var[(x,y)].varValue))

for (l, s) in LOC_SIZES:
    if buildVars[(l,s)].varValue > 0:
        print("Build %s %s (%s)" % \
              (l, s, buildVars[(l,s)].varValue))
print()

print_var_table(MINES, LOCATIONS, flowVars, \
                fn=lambda x: CC*x)
print()
print_var_table(LOCATIONS, CUSTOMERS, flowVars)

