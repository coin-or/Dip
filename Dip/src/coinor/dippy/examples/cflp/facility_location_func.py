from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from past.utils import old_div
import importlib as ilib
import argparse

from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize, LpAffineExpression

try:
    from src.dippy import DipProblem, DipSolStatOptimal
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem, DipSolStatOptimal
    from coinor.dippy.examples.gen_func import *

tol = pow(pow(2, -24), old_div(2.0, 3.0))

debug_print = False

debug_print_lp = False

def parseArgs():
    
    parser = argparse.ArgumentParser(
        description='Solve a facility location problem.')
    parser.add_argument('--module', '-m', metavar = 'module name', 
                        help='name of the Python module from which to import data')
    parser.add_argument('--useCustomSolver', action='store_true', 
                        help='enable custom subproblem solver')
    parser.add_argument('--customInitRule', choices=['initOneEach', 'initFirstFit'],
                        default = 'Default', 
                        help='enable custom method of initializing columns')
    parser.add_argument('--useCustomCuts', action='store_true', 
                        help='enable custom cut generation')
    parser.add_argument('--rootHeuristic', '-r', action='store_true', 
                        help='enable root heuristic')
    parser.add_argument('--nodeHeuristic', '-n', action='store_true', 
                        help='enable node heuristic')
    parser.add_argument('--numLocations', '-L', type=int, default = 5, metavar = 'L',
                        help = 'number of locations for instance specified on command line')
    parser.add_argument('--numProducts', '-P', type=int, default = 5, metavar = 'P',
                        help = 'number of products for instance specified on command line')
    parser.add_argument('--demands', '-D', nargs='*', default = [7, 5, 3, 2, 2],
                        help = 'list of demands for instances specified on command line')
    parser.add_argument('--fixedCosts', '-F', nargs='*', default = None,
                        help = 'list of fixed costs for instance specified on command line')
    parser.add_argument('--capacity', '-C', type=int, default = 8, metavar = 'C',
                        help = 'capacity of facilities for instance specified on command line')
    
    addDippyArgs(parser)

    args = parser.parse_args()

    if len(args.demands) != args.numProducts:
        raise ('Number of demands specified is not equal to number of products!')
    
    if args.fixedCosts != None and len(args.fixedCosts) != args.numLocations:
        raise ('Number of fixed costs specified is not equal to number of locations!')
    
    return(args)

def formulate(args):

    FIXED_COST = None
    ASSIGNMENTS = None
    ASSIGNMENT_COSTS = None

    if args.module:
        m = ilib.import_module(args.module)
        LOCATIONS = m.LOCATIONS
        PRODUCTS = m.PRODUCTS
        DEMAND = m.DEMAND
        CAPACITY = m.CAPACITY
        if hasattr(m, 'FIXED_COST'):
            FIXED_COST = m.FIXED_COST
        if hasattr(m, 'ASSIGNMENTS'):
            ASSIGNMENTS = m.ASSIGNMENTS
        if hasattr(m, 'ASSIGNMENT_COSTS'):
            ASSIGNMENT_COSTS = m.ASSIGNMENT_COSTS
    else:
        LOCATIONS = range(args.numLocations)
        PRODUCTS = range(args.numProducts)
        DEMAND = [int(i) for i in args.demands]
        CAPACITY = args.capacity
        if args.fixedCosts != None:
            FIXED_COST = [int(i) for i in args.fixedCosts]
        
    if FIXED_COST == None:
        FIXED_COST = {i:1 for i in LOCATIONS}
    if ASSIGNMENTS == None:
        ASSIGNMENTS = [(i, j) for i in LOCATIONS for j in PRODUCTS]
    if ASSIGNMENT_COSTS == None:
        ASSIGNMENT_COSTS = {i:0 for i in ASSIGNMENTS}

    prob = DipProblem("Facility Location")

    assign_vars = LpVariable.dicts("x", ASSIGNMENTS, 0, 1, LpBinary)
    use_vars    = LpVariable.dicts("y", LOCATIONS, 0, 1, LpBinary)

    prob += (lpSum(use_vars[i] * FIXED_COST[i] for i in LOCATIONS) +
             lpSum(assign_vars[j] * ASSIGNMENT_COSTS[j] for j in ASSIGNMENTS), 
             "min")

    # assignment constraints
    for j in PRODUCTS:
        prob += lpSum(assign_vars[(i, j)] for i in LOCATIONS) == 1

    # Aggregate capacity constraints
    for i in LOCATIONS:
        prob.relaxation[i] += lpSum(assign_vars[(i, j)] * DEMAND[j]
                                    for j in PRODUCTS) <= CAPACITY * use_vars[i]

    # Disaggregate capacity constraints
    for i, j in ASSIGNMENTS:
        prob.relaxation[i] += assign_vars[(i, j)] <= use_vars[i]

    prob.LOCATIONS = LOCATIONS
    prob.PRODUCTS = PRODUCTS
    prob.assign_vars = assign_vars
    prob.use_vars = use_vars

    return prob
        
def solve_subproblem(prob, key, redCosts, target):
    if debug_print:
        print("solve_subproblem...")
        print("reduced costs:")
        print(redCosts)
        print("target value:", target)
   
    loc = key

    # Calculate effective objective coefficient of products

    avars = [assign_vars[(loc, j)] for j in PRODUCTS]
    obj = [max(-redCosts[assign_vars[(loc, j)]], 0) for j in PRODUCTS]
    weights = [DEMAND[j] for j in PRODUCTS]
   
    # Use 0-1 KP to max. total effective value of products at location
    z, solution = knapsack01(obj, weights, CAPACITY)
   
    # Get the reduced cost of the knapsack solution
    if debug_print:
        print([(v, redCosts[v]) for v in avars])
        print(obj)
        print("z, solution =", z, solution)
        print("redCosts[use_vars[loc]] =", redCosts[use_vars[loc]])
        print("Fixed cost, rc", FIXED_COST[loc], redCosts[use_vars[loc]] - z)

    if redCosts[use_vars[loc]] > z + tol: # ... or an empty location is "useful"
        if debug_print:
            print("Zero solution is optimal")
        return DipSolStatOptimal, [{}]

    var_values = dict([(avars[i], 1) for i in solution])
    var_values[use_vars[loc]] = 1

    if debug_print:
        rcCheck = 0.0
        for v in list(var_values.keys()):
            rcCheck += redCosts[v] * var_values[v]
        print("Checking rc calc", redCosts[use_vars[loc]] - z, rcCheck) 
        print(var_values)

    return DipSolStatOptimal, [var_values]

def knapsack01(obj, weights, capacity):
    """ 0/1 knapsack solver, maximizes profit. weights and capacity integer """
        
    debug_subproblem = False
    
    assert len(obj) == len(weights)
    n = len(obj)
    if n == 0:
        return 0, []

    if debug_subproblem:
        relaxation = LpProblem('relaxation', LpMaximize)
        relax_vars = [str(i) for i in range(n)]
        var_dict   = LpVariable.dicts("", relax_vars, 0, 1, LpBinary)
        relaxation += (lpSum(var_dict[str(i)] * weights[i] for i in range(n)) 
                       <= capacity)
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

def generate_weight_cuts(prob, sol):
    # Define mu and T for each knapsack
    mu = {}
    S = {}
    for i in LOCATIONS:
        mu[i] = CAPACITY
        S[i] = []
    
    # Use current assign_var values to assign items to locations
    assigning = True 
    while assigning:
        bestValue = 0
        bestAssign = None
        for i in LOCATIONS:
            for j in PRODUCTS:
                if j not in S[i]: # If this product is not in the subset
                    if (sol[assign_vars[(i, j)]] > bestValue) \
                    and (DEMAND[j] <= mu[i]):
                        # The assignment variable for this product is closer
                        # to 1 than any other product checked, and "fits" in
                        # this location's remaining space
                        bestValue = sol[assign_vars[(i, j)]]
                        bestAssign = (i, j)
        # Make the best assignment found across all products and locactions
        if bestAssign:
            (i,j) = bestAssign
            mu[i] -= DEMAND[j] # Decrease spare CAPACITY at this location
            S[i].append(j) # Assign this product to this location's set
        else:
            assigning = False # Didn't find anything to assign - stop

    # Generate the weight cuts from the sets found above
    new_cuts = []
    for i in LOCATIONS:
        if len(S[i]) > 0: # If an item assigned to this location
            con = LpAffineExpression() # Start a new constraint
            con += sum(DEMAND[j] * assign_vars[(i, j)] 
                            for j in S[i])
            con += sum(max(0, DEMAND[j] - mu[i]) *
                            assign_vars[(i, j)] for j in PRODUCTS 
                            if j not in S[i])
            new_cuts.append(con <= CAPACITY - mu[i])

    # Return the set of cuts we created to DIP
    if len(new_cuts) > 0:
        return new_cuts

def first_fit_heuristic():
    # Sort the items in descending weight order
    productReqs = [(DEMAND[j],j) for j in PRODUCTS]
    productReqs.sort(reverse=True)

    # Add items to locations, fitting in as much
    # as possible at each location.
    allLocations = []
    while len(productReqs) > 0:
        waste = CAPACITY
        currentLocation = []
        j = 0
        while j < len(productReqs):
            # Can we fit this product?
            if productReqs[j][0] <= waste:
                currentLocation.append(productReqs[j][1]) # index
                waste -= productReqs[j][0] # requirement
                productReqs.pop(j)
            else:
                # Try to fit next item
                j += 1
        allLocations.append((currentLocation, waste))

    # Return a list of tuples: ([products],waste)
    return allLocations

def first_fit():
    # Use generic first-fit heuristic that is shared
    # between heuristics and initial variable generation
    allLocations = first_fit_heuristic()

    # Convert to decision variable values
    sol = {}
    for i in LOCATIONS:
        for j in PRODUCTS:
            sol[assign_vars[(i, j)]] = 0
        sol[use_vars[i]] = 0

    index = 0
    for loc in allLocations:
        i = LOCATIONS[index]
        sol[use_vars[i]] = 1
        for j in loc[0]:
            sol[assign_vars[(i, j)]] = 1
        index += 1
        
    return sol

def frac_fit(xhat):
    # Initialise solution
    sol = {}
    waste = {}
    for i in LOCATIONS:
        for j in PRODUCTS: sol[assign_vars[(i, j)]] = 0
        sol[use_vars[i]] = 0
        waste[i] = 0
        
    # Get the list of non-zero fractional assignments
    fracAssigns = [ (xhat[assign_vars[(i, j)]], (i, j))
                     for i in LOCATIONS for j in PRODUCTS
                     if xhat[assign_vars[(i, j)]] > tol ]
    fracAssigns.sort()

    # Track which products and locations have been used
    notAllocated = dict((j,True) for j in PRODUCTS)
    notUsed      = dict((i,True) for i in LOCATIONS)
    while len(fracAssigns) > 0:
        fracX = fracAssigns.pop() # Get best frac. assignment left
        (i,j) = fracX[1]
        if notAllocated[j]:
            if notUsed[i]: # Create a new location if needed
                notUsed[i] = False
                sol[use_vars[i]] = 1
                waste[i] = CAPACITY
            if DEMAND[j] <= waste[i]: # Space left?
                sol[assign_vars[(i, j)]] = 1
                notAllocated[j] = False
                waste[i] -= DEMAND[j]
    
    # Allocate the remaining products
    unallocated = [(DEMAND[j],j) for j in PRODUCTS
                                      if notAllocated[j]]
    unallocated.sort(reverse=True)
    unused = [i for i in LOCATIONS if notUsed[i]]
    while len(unallocated) > 0:
        waste = CAPACITY
        index = 0
        loc = unused.pop()
        while index < len(unallocated):
            (j_req, j) = unallocated[index]
            if j_req <= waste:
                unallocated.pop(index)
                sol[assign_vars[(loc, j)]] = 1
                waste -= j_req
            else: index += 1
        sol[use_vars[loc]] = 1

    return sol

def heuristics(prob, xhat, cost):
    sols = []
    if prob.root_heuristic:
        prob.root_heuristic = False # Don't run twice
        sol = first_fit()
        sols.append(sol)
    if prob.node_heuristic:
        sol = frac_fit(xhat)
        sols.append(sol)

    if len(sols) > 0:
        return sols
      
def init_first_fit(prob):

    locations = first_fit_heuristic()
    bvs = []
    index = 0
    for loc in locations:
        i = LOCATIONS[index]
        if debug_print:
            print([(assign_vars[(i, j)], 1) for j in loc[0]])
        var_values = dict([(assign_vars[(i, j)], 1) for j in loc[0]])
        var_values[use_vars[i]] = 1
        dv = (loc[1], var_values)
        bvs.append((i, dv))
        index += 1
    if debug_print:
        print(bvs)
    return bvs

def init_one_each(prob):
    bvs = []
    if debug_print:
        print("LOCATIONS =", LOCATIONS)
    for index, loc in enumerate(LOCATIONS):
        lc = [PRODUCTS[index]]
        waste = CAPACITY - DEMAND[PRODUCTS[index]]
        var_values = dict([(assign_vars[(loc, j)], 1) for j in lc])
        var_values[use_vars[loc]] = 1

        dv = (waste, var_values)
        bvs.append((loc, dv))
    if debug_print:
        print(bvs)
    return bvs

