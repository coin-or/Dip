DEBUGGING = False

import sys

from pulp import *

if DEBUGGING:

    try:
        import path
    except ImportError:
        pass
    import dippy
else:
    import coinor.dippy as dippy

from math import floor, ceil

class BinPackProb:
    def __init__(self, ITEMS, volume, capacity):
        self.ITEMS = ITEMS
        self.volume = volume
        self.BINS = range(len(ITEMS)) # Create 1 bin for each item, indices start at 0
        self.capacity = capacity
    
def formulate(bpp):
    prob = dippy.DipProblem("Bin Packing",
                            display_mode = 'xdot',
#                           layout = 'bak',
                            display_interval = None,
                            )

    assign_vars = LpVariable.dicts("x",
                                   [(i, j) for i in bpp.ITEMS
                                    for j in bpp.BINS],
                                   cat=LpBinary)
    use_vars    = LpVariable.dicts("y", bpp.BINS, cat=LpBinary)
    waste_vars  = LpVariable.dicts("w", bpp.BINS, 0, None)

    prob += lpSum(waste_vars[j] for j in bpp.BINS), "min_waste"

    for j in bpp.BINS:
        prob += lpSum(bpp.volume[i] * assign_vars[i, j] for i in bpp.ITEMS) \
                + waste_vars[j] == bpp.capacity * use_vars[j]

    for i in bpp.ITEMS:
        prob += lpSum(assign_vars[i, j] for j in bpp.BINS) == 1

    for i in bpp.ITEMS:
        for j in bpp.BINS:
            prob += assign_vars[i, j] <= use_vars[j]

    for n in range(0, len(bpp.BINS) - 1):
        prob += use_vars[bpp.BINS[n]] >= use_vars[bpp.BINS[n + 1]]

    # Attach the problem data and variable dictionaries to the DipProblem 
    prob.bpp         = bpp
    prob.assign_vars = assign_vars
    prob.use_vars    = use_vars
    prob.waste_vars  = waste_vars

    return prob

def my_branch(prob, sol):
    bounds = None
   
    bounds = symmetry(prob, sol)
  
    if bounds is None:
        bounds = most_frac_use(prob, sol)
    
    if bounds is None:
        bounds = most_frac_assign(prob, sol)
    
    return bounds

def my_heuristics(prob, xhat, cost):
#  print "Heuristics..."
    sol = None
  
    if prob.is_root_node:
        prob.is_root_node = False
        if prob.root_heuristic:
            sol = first_fit(prob)
    else:
        if prob.node_heuristic:
            sol = frac_fit(prob, xhat)
  
    if sol is not None:
        return [sol]

def solve(prob):

#    prob.branch_method = my_branch
    prob.heuristics = my_heuristics
    prob.is_root_node = True
    prob.root_heuristic = True
    prob.node_heuristic = True
  
    dippyOpts = {'doPriceCut' : '1',
                 'CutCGL': '0',
#                'SolveMasterAsIp': '0'
#                'generateInitVars': '1',
#                 'LogDebugLevel': 5,
#                'LogDumpModel': 5,
                 }

    status, message, primals, duals = dippy.Solve(prob, dippyOpts)
  
    if status == LpStatusOptimal:
        return dict((var, var.value()) for var in prob.variables())
    else:
        return None

def most_frac_use(prob, sol):
    # Get the attached data and variable dicts
    bpp         = prob.bpp
    use_vars    = prob.use_vars
    tol         = prob.tol
  
    most   = float('-inf')
    bin = None
    for j in bpp.BINS:
        alpha = sol[use_vars[j]]
        up   = ceil(alpha)  # Round up to next nearest integer 
        down = floor(alpha) # Round down
        frac = min(up - alpha, alpha - down)
        if frac > tol: # Is fractional?
            if frac > most:
                most = frac
                bin = j
    
    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}
    if bin is not None:
#        print bin, sol[use_vars[bin]]
        down_ubs[use_vars[bin]] = 0.0
        up_lbs[use_vars[bin]] = 1.0
    
        return down_lbs, down_ubs, up_lbs, up_ubs

def most_frac_assign(prob, sol):
    # Get the attached data and variable dicts
    bpp         = prob.bpp
    assign_vars = prob.assign_vars
    tol         = prob.tol
  
    most   = float('-inf')
    assign = None
    for i in bpp.ITEMS:
        for j in bpp.BINS:
            up   = ceil(sol[assign_vars[i, j]])  # Round up to next nearest integer 
            down = floor(sol[assign_vars[i, j]]) # Round down
            frac = min(up - sol[assign_vars[i, j]], sol[assign_vars[i, j]] - down)
            if frac > tol: # Is fractional?
                if frac > most:
                    most = frac
                    assign = (i, j)
    
    down_lbs = {}
    down_ubs = {}
    up_lbs = {}
    up_ubs = {}
    if assign is not None:
#    print assign, sol[assign_vars[assign]]
        down_ubs[assign_vars[assign]] = 0.0
        up_lbs[assign_vars[assign]] = 1.0
    
        return down_lbs, down_ubs, up_lbs, up_ubs

def symmetry(prob, sol):
    # Get the attached data and variable dicts
    bpp      = prob.bpp
    use_vars = prob.use_vars
    tol      = prob.tol
  
    alpha = sum(sol[use_vars[j]] for j in bpp.BINS)
#      print "# bins =", alpha
    up   = int(ceil(alpha))  # Round up to next nearest integer 
    down = int(floor(alpha)) # Round down
    frac = min(up - alpha, alpha - down)
    if frac > tol: # Is fractional?
#    print "Symmetry branch"
    
        down_lbs = {}
        down_ubs = {}
        up_lbs = {}
        up_ubs = {}
        for n in range(up - 1, len(bpp.BINS)):
            down_ubs[use_vars[bpp.BINS[n]]] = 0.0
#           print down_ubs
        for n in range(up):
            up_lbs[use_vars[bpp.BINS[n]]] = 1.0
#           print up_lbs

        return down_lbs, down_ubs, up_lbs, up_ubs
  
def fit(prob, order):
    bpp         = prob.bpp
    use_vars    = prob.use_vars
    assign_vars = prob.assign_vars
    waste_vars  = prob.waste_vars
    tol         = prob.tol

    sol = {}

    for j in bpp.BINS:
        sol[use_vars[j]] = 1.0
        for i in bpp.ITEMS:
            sol[assign_vars[i, j]] = 0.0
        sol[waste_vars[j]] = bpp.capacity

    assigned = {}
    for i in bpp.ITEMS:
        assigned[i] = False

    for (i, j) in order:
        if (not assigned[i]) and (bpp.volume[i] <= sol[waste_vars[j]]):
            assigned[i] = True
            sol[assign_vars[i, j]] = 1.0
            sol[waste_vars[j]] -= bpp.volume[i]
  
    for j in bpp.BINS:
        if sol[waste_vars[j]] > bpp.capacity - tol:
            sol[use_vars[j]] = 0.0
            sol[waste_vars[j]] = 0.0
      
#   print sol
  
    return sol

import operator

def first_fit(prob):
#  print "first fit..."
  
    bpp = prob.bpp
    sorted_volume = sorted(bpp.volume.iteritems(), key=operator.itemgetter(1), reverse=True)
    sorted_ITEMS = [i for (i, v) in sorted_volume]
  
#   print sorted_ITEMS

    order = [(i, j) for i in sorted_ITEMS for j in bpp.BINS]
  
#   print order
  
    sol = fit(prob, order)
#   print sol
    return sol

def frac_fit(prob, xhat):
#    print "frac fit..."
  
    bpp = prob.bpp
    assign_vars = prob.assign_vars
  
    assign = dict(((i, j), xhat[assign_vars[i, j]]) for i in bpp.ITEMS for j in bpp.BINS)
#   print assign
  
    sorted_assign = sorted(assign.iteritems(), key=operator.itemgetter(1), reverse=True)
    order = [(i, j) for ((i, j), x) in sorted_assign]
   
#   print order
  
    sol = fit(prob, order)
#   print sol
    return sol
         
