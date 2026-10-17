from __future__ import division
from __future__ import absolute_import
from builtins import range
from past.utils import old_div
import unittest
from pulp import *
import coinor.dippy as dippy

from .dippy_tests import DippyTestCase


class TestFacilityProblem(DippyTestCase):

    def setUp(self):
        """
        sets up the coke problem as a unittest
        """
        self.prob, self.branch, self.cuts = create_facility_problem()
        self.tol = pow(pow(2, -24), old_div(2.0, 3.0))

    def test_pulp_solve(self):
        """
        tests that pulp can solve the problem
        """
        status = self.prob.solve()
        self.assertEqual(status, LpStatusOptimal)
        self.assertAlmostEqual(self.prob.objective.value(), 212.0)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)

    #~ def test_dippy_branch(self):
        #~ """
        #~ tests that the custom branch can solve the problem
        #~ """
        #~ self.prob.branch_method = self.branch
        #~ dippy.Solve(self.prob, {
            #~ 'TolZero': '%s' % self.tol,
            #~ })
        #~ self.assertAlmostEqual(self.prob.objective.value(), 212.0)
        #~ self.variable_feasibility_test(self.prob)
        #~ self.constraint_feasibility_test(self.prob)

    def test_dippy_cuts(self):
        """
        tests that the custom branch can solve the problem
        """
        self.prob.generate_cuts = self.cuts
        dippy.Solve(self.prob, {
            'TolZero': '%s' % self.tol,
            })
        self.assertAlmostEqual(self.prob.objective.value(), 212.0)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)

    def test_dippy_branch_cuts(self):
        """
        tests that the custom branch can solve the problem
        """
        self.prob.branch_method = self.branch
        self.prob.generate_cuts = self.cuts
        dippy.Solve(self.prob, {
            'TolZero': '%s' % self.tol,
            })
        self.assertAlmostEqual(self.prob.objective.value(), 212.0)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)



def create_facility_problem():
    """
    creates and returns the facility problem
    """

    from math import floor, ceil

    tol = pow(pow(2, -24), old_div(2.0, 3.0))

    # The requirements for the products

    REQUIREMENT = {
        1 : 66,
        2 : 4,
        3 : 85,
        4 : 93,
        5 : 68,
        6 : 76,
        7 : 74,
        8 : 39,
        9 : 66,
        10 : 17,
    }
    # Set of all products
    PRODUCTS = list(REQUIREMENT.keys())
    PRODUCTS.sort()

    # Set of all locations

    LOCATIONS = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    LOCATIONS.sort()

    # The capacity of the facilities

    CAPACITY = 100

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

    # aggregate CAPACITY constraints
    for i in LOCATIONS:
        prob += lpSum(assign_vars[(i, j)] * REQUIREMENT[j]
                      for j in PRODUCTS) + waste_vars[i] == \
                CAPACITY * use_vars[i]

    # disaggregate CAPACITY constraints
    for i in LOCATIONS:
        for j in PRODUCTS:
            prob += assign_vars[(i, j)] <= use_vars[i]

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
            down_branch_ub = dict([(use_vars[LOCATIONS[n]], 0)
                                   for n in range(int(down), len(LOCATIONS))])
            # Up branch: provide lower bounds, upper bounds are default
            up_branch_lb = dict([(use_vars[LOCATIONS[n]], 1)
                                 for n in range(0, int(up))])
            # Return the advanced branch to DIP
            return {}, down_branch_ub, up_branch_lb, {}

    def generate_weight_cuts(prob, sol):
    ##    print "In generate_weight_cuts, sol = ", sol

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
                        and (REQUIREMENT[j] <= mu[i]):
                            # The assignment variable for this product is closer
                            # to 1 than any other product checked, and "fits" in
                            # this location's remaining space
                            bestValue = sol[assign_vars[(i, j)]]
                            bestAssign = (i, j)
            # Make the best assignment found across all products and locactions
            if bestAssign:
                (i,j) = bestAssign
                mu[i] -= REQUIREMENT[j] # Decrease spare CAPACITY at this location
                S[i].append(j) # Assign this product to this location's set
            else:
                assigning = False # Didn't find anything to assign - stop

        # Generate the weight cuts from the sets found above
        new_cuts = []
        for i in LOCATIONS:
            if len(S[i]) > 0: # If an item assigned to this location
                con = LpAffineExpression() # Start a new constraint
                con += sum(REQUIREMENT[j] * assign_vars[(i, j)]
                                for j in S[i])
                con += sum(max(0, REQUIREMENT[j] - mu[i]) *
                                assign_vars[(i, j)] for j in PRODUCTS
                                if j not in S[i])
                new_cuts.append(con <= CAPACITY - mu[i])

        ##    print new_cuts
        # Return the set of cuts we created to DIP
        return new_cuts

        return sols


    return prob, choose_antisymmetry_branch, generate_weight_cuts


if __name__ == '__main__':
    unittest.main()
