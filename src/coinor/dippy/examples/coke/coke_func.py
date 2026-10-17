from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import argparse
from pulp import LpVariable, lpSum, LpBinary, LpStatusOptimal

try:
    from src.dippy import DipProblem, Solve
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem, Solve
    from coinor.dippy.examples.gen_func import *

class CokeProb(object):
    def __init__(self, supply, demand, LOCATIONS, build_costs,
                 conversion_factor, transport_costs):
        self.MINES     = list(supply.keys())
        self.MINES.sort()
        self.CUSTOMERS = list(demand.keys())
        self.CUSTOMERS.sort()
        self.LOCATIONS = LOCATIONS
        self.SIZES     = list(build_costs.keys())
        self.SIZES.sort()
        self.ARCS      = list(transport_costs.keys())
        self.conversion_factor = conversion_factor
        self.supply            = supply
        self.demand            = demand
        self.build_costs       = build_costs
        self.transport_costs   = transport_costs

def parseArgs():
    
    parser = argparse.ArgumentParser(
        description='Solve a steel manufacturing problem.')
    parser.add_argument('--module', '-m', metavar = 'module name', 
                        help='name of the Python module from which to import data',
                        default = 'coinor.dippy.examples.coke.coke_data')
    parser.add_argument('--advancedBranch', action='store_true', 
                        help='enable advanced branching')

    addDippyArgs(parser)

    args = parser.parse_args()

    return args

def formulate(cp):

    prob = DipProblem("Coke", display_mode = 'xdot', display_interval = None)

    # create variables
    LOC_SIZES = [(l, s) for l in cp.LOCATIONS
                        for s in cp.SIZES]
    buildVars = LpVariable.dicts("Build", LOC_SIZES, cat=LpBinary)

    # create arcs
    flowVars = LpVariable.dicts("Arcs", cp.ARCS)
    BIG_M = max(sum(cp.supply.values()), sum(cp.demand.values()))
    for a in cp.ARCS:
        flowVars[a].bounds(0, BIG_M)

    # objective
    prob += 1e6 * lpSum(buildVars[(l, s)] * cp.build_costs[s] \
                        for (l, s) in LOC_SIZES) + \
                  lpSum(flowVars[(s, d)] * cp.transport_costs[(s, d)] \
                        for (s, d) in cp.ARCS), "min"

    # plant availability - assumes that SIZES are numeric,
    # which they should be
    for loc in cp.LOCATIONS:
        prob += lpSum(flowVars[(loc, i)] for i in cp.CUSTOMERS) \
             <= lpSum(buildVars[(loc, s)] * s for s in cp.SIZES)

    # one size
    for loc in cp.LOCATIONS:
        prob += lpSum(buildVars[(loc, s)] for s in cp.SIZES) == 1

    # conserve flow (mines)
    # flows are in terms of tonnes of coke
    for m in cp.MINES:
        prob += lpSum(flowVars[(m, j)] for j in cp.LOCATIONS) \
             <= cp.supply[m]

    # conserve flow (locations)
    # convert from coal to coke
    for loc in cp.LOCATIONS:
        prob += lpSum(flowVars[(m, loc)] for m in cp.MINES) - \
                cp.conversion_factor * \
                lpSum(flowVars[(loc, c)] for c in cp.CUSTOMERS) \
             >= 0

    for c in cp.CUSTOMERS:
        prob += lpSum(flowVars[(loc, c)] for loc in cp.LOCATIONS) \
             >= cp.demand[c]

    prob.cp = cp
    prob.buildVars = buildVars
    prob.flowVars = flowVars
    
    return prob
  
def do_branch(prob, sol):
    tol       = prob.tol
    LOCATIONS = prob.cp.LOCATIONS
    SIZES     = prob.cp.SIZES
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

        down_lbs = {}
        down_ubs = {}
        up_lbs = {}
        up_ubs = {}
        if bigger - sol_size > tol and sol_size - smaller > tol:
            down_ubs = dict([(buildVars[loc, SIZES[j]], 0) 
                             for j in range(i, len(SIZES))])
            up_ubs = dict([(buildVars[loc, SIZES[j]], 0) 
                           for j in range(0, i)])

            return down_lbs, down_ubs, up_lbs, up_ubs
        
def solve(prob, args):

    if args.advancedBranch:
        prob.branch_method = do_branch

    if args.algo != 'Cut':
        print("This application only works with the cutting plane method")
        
    dippyOpts = addDippyOpts(args)

    status, message, primals, duals = Solve(prob, dippyOpts)
    
    if status == LpStatusOptimal:
        return dict((var, var.value()) for var in prob.variables())
    else:
        return None

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

def print_table(rows, cols, fn):
    print("\t", "\t".join(cols))
    for r in rows:
        print(r,"\t", "\t".join(str(fn(r,c)) for c in cols))

def print_var_table(rows, cols, var, fn=lambda x: x):
    print_table(rows, cols, lambda x, y:
                fn(var[(x,y)].varValue))
