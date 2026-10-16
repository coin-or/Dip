# Generalized Assignment Problem
# argument should be a problem file, see Dip/examples/GAP_Instance.cpp for format
# for an e.g. see gap0512-2.dat included in this directory

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import importlib as ilib
import argparse
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize

try:
    from src.dippy import DipProblem, DipSolStatOptimal
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem, DipSolStatOptimal
    from coinor.dippy.examples.gen_func import *

debug_print = False

tol = pow(pow(2, -24), old_div(2.0, 3.0))

def parseArgs():
    
    parser = argparse.ArgumentParser(
        description='Solve a generalized assignment problem.')
    parser.add_argument('--module', '-m', metavar = 'module name', 
                        help='name of the Python module from which to import data',
                        default = 'coinor.dippy.examples.gap.gap0515-2')

    addDippyArgs(parser)

    args = parser.parse_args()

    return args

def formulate(module_name):

    m = ilib.import_module(module_name)

    lines = m.gap_data.splitlines()
    line = lines[1].split() #first line is blank
    NUM_MACHINES = int(line[0])
    NUM_TASKS = int(line[1])
    MACHINES = list(range(NUM_MACHINES))
    TASKS = list(range(NUM_TASKS))
    MACHINES_TASKS = [(m, t) for m in MACHINES for t in TASKS]

    COSTS = []
    line_num = 2
    for m in MACHINES:
        line = lines[line_num].split()
        assert len(line) == NUM_TASKS
        COSTS.append([int(f) for f in line])
        line_num += 1

    RESOURCE_USE = []
    for m in MACHINES:
        line = lines[line_num].split()
        assert len(line) == NUM_TASKS
        RESOURCE_USE.append([int(f) for f in line])
        line_num += 1

    line = lines[line_num].split()
    assert len(line) == NUM_MACHINES
    CAPACITIES = [int(f) for f in line]

    assignVars = []
    for m in MACHINES:
        v = []
        for t in TASKS:
            v.append(LpVariable("M%dT%d" % (m, t), cat=LpBinary))
        assignVars.append(v)

    prob = DipProblem("GAP")

    # objective
    prob += lpSum(assignVars[m][t] * COSTS[m][t] for m, t in MACHINES_TASKS), "min"

    # machine capacity (knapsacks, relaxation)
    for m in MACHINES:
        prob.relaxation[m] += lpSum(assignVars[m][t] * RESOURCE_USE[m][t] for t in TASKS) <= CAPACITIES[m]

        # assignment
    for t in TASKS:
        prob += lpSum(assignVars[m][t] for m in MACHINES) == 1

    prob.assignVars = assignVars
    prob.MACHINES = MACHINES
    prob.TASKS = TASKS

    return prob

def solve_subproblem(prob, machine, redCosts, target):
    if debug_print:
        print("solve_subproblem...")
        print(redCosts)
   
    # get tasks which have negative reduced costs
    task_idx = [t for t in TASKS if redCosts[assignVars[machine][t]] < 0]
    var = [assignVars[machine][t] for t in task_idx]
    obj = [-redCosts[assignVars[machine][t]] for t in task_idx]
    weights = [RESOURCE_USE[machine][t] for t in task_idx]

    z, solution = knapsack01(obj, weights, CAPACITIES[machine])

    # Get the reduced cost of the knapsack solution and waste
    if debug_print:
        print([(v, redCosts[v]) for v in var])
        print(obj)
        print("z, solution =", z, solution)
        print("rc", -z)

    if z < -tol: # Zero solution is optimal
        if debug_print:
            print("Zero solution is optimal")
        return DipSolStatOptimal, [{}]

    var_values = dict([(var[i], 1) for i in solution])

    if debug_print:
        rcCheck = 0.0
        for v in list(var_values.keys()):
            rcCheck += redCosts[v] * var_values[v]
        print("Checking rc calc", -z, rcCheck) 
        print(var_values)
        
    return DipSolStatOptimal, [var_values]

def knapsack01(obj, weights, capacity):
    """ 0/1 knapsack solver, maximizes profit. weights and capacity integer """
        
    debug_subproblem = False
    
    assert len(obj) == len(weights)
    n = len(obj)
    if n == 0:
        return 0, []

    if (debug_subproblem):
        relaxation = LpProblem('relaxation', LpMaximize)
        relax_vars = [str(i) for i in range(n)]
        var_dict   = LpVariable.dicts("", relax_vars, 0, 1, LpBinary)
        relaxation += lpSum(var_dict[str(i)] * weights[i] for i in range(n)) <= capacity
        relaxation += lpSum(var_dict[str(i)] * obj[i] for i in range(n))
        relaxation.solve()
        relax_obj = value(relaxation.objective)

        solution =  [i for i in range(n) if var_dict[str(i)].varValue > tol ]

        print(relax_obj, solution)


    c = [[0]*(capacity+1) for i in range(n)]
    added = [[False]*(capacity+1) for i in range(n)]
    # c [items, remaining capacity]
    # important: this code assumes strictly positive objective values
    for i in range(n):
        for j in range(capacity+1):
            if (weights[i] > j):
                c[i][j] = c[i-1][j]
            else:
                c_add = obj[i] + c[i-1][j-weights[i]]
                if c_add > c[i-1][j]:
                    c[i][j] = c_add
                    added[i][j] = True
                else:
                    c[i][j] = c[i-1][j]

    # backtrack to find solution
    i = n-1
    j = capacity

    solution = []
    while i >= 0 and j >= 0:
        if added[i][j]:
            solution.append(i)
            j -= weights[i]
        i -= 1
        
    return c[n-1][capacity], solution

