import coinor.pulp as pulp
from pulp import *
import coinor.dippy as dippy

from math import floor, ceil

tol = pow(pow(2, -24), 2.0 / 3.0)

from facility_ex1 import REQUIREMENT, PRODUCTS, LOCATIONS, CAPACITY
##from facility_test1 import REQUIREMENT, PRODUCTS, LOCATIONS, CAPACITY
##from facility_test2 import REQUIREMENT, PRODUCTS, LOCATIONS, CAPACITY
##from facility_test6 import REQUIREMENT, PRODUCTS, LOCATIONS, CAPACITY

prob = dippy.DipProblem("Facility_Location")

assign_vars = LpVariable.dicts("AtLocation",
              [(i, j) for i in LOCATIONS
                      for j in PRODUCTS],
              0, 1, LpBinary)
use_vars    = LpVariable.dicts("UseLocation",
              LOCATIONS, 0, 1, LpBinary)
waste_vars  = LpVariable.dicts("Waste",
              LOCATIONS, 0, CAPACITY)

# objective: minimise waste
prob += lpSum(waste_vars[i] for i in LOCATIONS), "min"

# assignment constraints
for j in PRODUCTS:
    prob += lpSum(assign_vars[(i, j)] for i in LOCATIONS) == 1

# Aggregate capacity constraints
for i in LOCATIONS:
    prob.relaxation[i] += lpSum(assign_vars[(i, j)] * REQUIREMENT[j]
                                for j in PRODUCTS) + waste_vars[i] \
                                    == CAPACITY * use_vars[i]

# Disaggregate capacity constraints
for i in LOCATIONS:
    for j in PRODUCTS:
        prob.relaxation[i] += assign_vars[(i, j)] <= use_vars[i]
        
# Ordering constraints
for index, location in enumerate(LOCATIONS):
    if index > 0:
        prob += use_vars[LOCATIONS[index-1]] >= use_vars[location]
        
# Anti-symmetry branches
def choose_antisymmetry_branch(prob, sol):
    num_locations = sum(sol[use_vars[i]] for i in LOCATIONS)
    up   = ceil(num_locations)  # Round up to next nearest integer 
    down = floor(num_locations) # Round down
    if  (up - num_locations   > tol) \
    and (num_locations - down > tol): # Is fractional?
        # Down branch: provide upper bounds, lower bounds are default
        down_branch_ub = [(use_vars[LOCATIONS[n]], 0) for
                                n in range(int(down), len(LOCATIONS))]
        # Up branch: provide lower bounds, upper bounds are default
        up_branch_lb = [(use_vars[LOCATIONS[n]], 1) for
                                n in range(0, int(up))]
        # Return the advanced branch to DIP
        return ([], down_branch_ub, up_branch_lb, [])

prob.branch_method = choose_antisymmetry_branch

def solve_subproblem(prob, index, redCosts, convexDual):
   loc = index

   # Calculate effective objective coefficient of products
   effs = {}
   for j in PRODUCTS:
       effs[j] = redCosts[assign_vars[(loc, j)]] \
               - redCosts[waste_vars[loc]] * REQUIREMENT[j]

   avars = [assign_vars[(loc, j)] for j in PRODUCTS]
   obj = [-effs[j] for j in PRODUCTS]
   weights = [REQUIREMENT[j] for j in PRODUCTS]
   
   # Use 0-1 KP to max. total effective value of products at location
   z, solution = knapsack01(obj, weights, CAPACITY)
   
   # Get the reduced cost of the knapsack solution and waste
   rc = redCosts[use_vars[loc]] -z + \
        redCosts[waste_vars[loc]] * CAPACITY
   waste = CAPACITY - sum(weights[i] for i in solution)
   rc += redCosts[waste_vars[loc]] * waste

   # Return the solution if the reduced cost is low enough
   if rc < -tol: # The location is used and "useful"
       if rc - convexDual < -tol:
           var_values = [(avars[i], 1) for i in solution]
           var_values.append((use_vars[loc], 1))
           var_values.append((waste_vars[loc], waste))

           dv = dippy.DecompVar(var_values, rc - convexDual, waste)
           return [dv]

   elif -convexDual < -tol: # An empty location is "useful"
           var_values = []

           dv = dippy.DecompVar(var_values, -convexDual, 0.0)
           return [dv]
       
   return []

prob.relaxed_solver = solve_subproblem

def knapsack01(obj, weights, capacity):
    """ 0/1 knapsack solver, maximizes profit. weights and capacity integer """
    assert len(obj) == len(weights)
    n = len(obj)

    if n == 0:
        return 0, []

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

def first_fit_heuristic():
    # Sort the items in descending weight order
    productReqs = [(REQUIREMENT[j],j) for j in PRODUCTS]
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
    
	
def first_fit(prob):
    locations = first_fit_heuristic()
    bvs = []
    index = 0
    for loc in locations:
        i = LOCATIONS[index]
        var_values = [(assign_vars[(i, j)], 1) for j in loc[0]]
        var_values.append((use_vars[i], 1))
        var_values.append((waste_vars[i], loc[1]))
        dv = dippy.DecompVar(var_values, None, loc[1])
        bvs.append((i, dv))
        index += 1
    return bvs

def one_each(prob):
   bvs = []
   for index, loc in enumerate(LOCATIONS):
       lc = [PRODUCTS[index]]
       waste = CAPACITY - REQUIREMENT[PRODUCTS[index]]
       var_values = [(assign_vars[(loc, j)], 1) for j in lc]
       var_values.append((use_vars[loc], 1))
       var_values.append((waste_vars[loc], waste))

       dv = dippy.DecompVar(var_values, None, waste)
       bvs.append((loc, dv))
   return bvs

prob.init_vars = first_fit
##prob.init_vars = one_each

prob.writeLP('facility_main.lp')
for n, i in enumerate(LOCATIONS):
    prob.writeRelaxed(n, 'facility_relax%s.lp' % i);

dippy.Solve(prob, {
    'TolZero': '%s' % tol,
    'doPriceCut': '1',
    'generateInitVars': '1', })

# print solution
for i in LOCATIONS:
    if use_vars[i].varValue > tol:
        print "Location ", i, \
              " produces ", \
              [j for j in PRODUCTS
               if assign_vars[(i, j)].varValue > tol]
