from __future__ import absolute_import
from builtins import range
import unittest
from pulp import *
import coinor.dippy as dippy

from .dippy_tests import DippyTestCase


class TestCuttingStockProblem(DippyTestCase):


    def test_pulp_solve(self):
        """
        tests that pulp can solve the problem
        """
        self.prob, self.relaxation = create_cutting_stock_problem()
        status = self.prob.solve()
        self.assertEqual(status, LpStatusOptimal)
        self.assertAlmostEqual(self.prob.objective.value(), 2.0)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)

    def test_dippy_relaxation(self):
        """
        tests that the custom branch can solve the problem
        """
        self.prob, self.relaxation = create_cutting_stock_problem(doRelaxed=True)
        self.prob.relaxed_solver = self.relaxation
        dippy.Solve(self.prob, {
            'doPriceCut':1,
            'CutCGL': 0,
        })
        self.assertAlmostEqual(self.prob.objective.value(), 2.0)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)


def create_cutting_stock_problem(doRelaxed=False):
    """
    creates and returns the cutting_stock problem
    """

    length = {
    "9cm": 9,
    "7cm": 7,
    "5cm": 5
    }

    ITEMS = list(length.keys())

    demand = {
    "9cm": 1,
    "7cm": 1,
    "5cm": 4
    }

    total_patterns = sum(demand[i] for i in ITEMS)

    total_length = {}
    for p in range(total_patterns):
        total_length[p] = 20
    PATTERNS = list(total_length.keys())

    CUTS = [(p, i) for p in PATTERNS
                   for i in ITEMS]

    prob = dippy.DipProblem("Sponge_Rolls", LpMinimize)

    # create variables

    useVars = LpVariable.dicts("Use", PATTERNS, 0, 1, LpBinary)

    cutVars = LpVariable.dicts("Cut", CUTS, 0, 10, LpInteger)

    # objective
    prob += sum(useVars[p] for p in PATTERNS), "min"

    # Meet demand
    for i in ITEMS:
        prob += sum(cutVars[(p, i)] for p in PATTERNS) \
                >= demand[i]

    # Ordering patterns
    for i, p in enumerate(PATTERNS):
        if p != PATTERNS[-1]:
             prob += useVars[p] >= useVars[PATTERNS[i+1]]

    # Cut patterns
    for p in PATTERNS:
        if doRelaxed:
            prob.relaxation[p] += sum(length[i] *
                                      cutVars[(p, i)] for i in ITEMS) \
                                      <= total_length[p] * useVars[p]
        else:
            prob += sum(length[i] * cutVars[(p, i)] for i in ITEMS) \
                 <= total_length[p] * useVars[p]


    def relaxed_solver(prob, patt, redCosts, convexDual):
    ##    print patt, "in", PATTERNS
    ##    print "redCosts =", redCosts
    ##    print "convexDual =", convexDual
        # get items with negative reduced cost
        item_idx = [i for i in ITEMS \
                    if redCosts[cutVars[(patt, i)]] < 0]
        vars = [cutVars[(patt, i)] for i in item_idx]
        obj = [-redCosts[cutVars[(patt, i)]] for i in item_idx]
        weights = [length[i] for i in item_idx]

    ##    print "Using knapsack heuristic"
    ##    print "item_idx =", item_idx
    ##    print "obj =", obj
    ##    print "weights =", weights
    ##    print "capacity =", total_length[patt]
        z, solution = kp(obj, weights, total_length[patt])
    ##    print "Number of items = ", len(item_idx)
    ##    for i in range(len(item_idx)):
    ##        print "Item ", item_idx[i], " has profit ", obj[i], " and weight ", weights[i]
    ##    print "Knapsack has capacity ", total_length[patt]
    ##    print "Value = ", z
    ##    for i in range(len(item_idx)):
    ##        print "Included[", item_idx[i], "] = ", solution[i]

        total_weight = sum(w * solution[i] \
                           for i, w in enumerate(weights))
        assert total_weight <= total_length[patt]

        # add in reduced cost of useVars
        totalCut = sum(solution)
    ##    print "z, redCosts[useVars[", patt, "]], convexDual", z, redCosts[useVars[patt]], convexDual
        if totalCut > 0:
            rc = -z + redCosts[useVars[patt]] - convexDual
        else:
            rc = -convexDual

    ##    print "rc =", rc
    ##    sys.stdout.flush()
        if rc < 0: # Using this pattern
            var_values = dict([(v, solution[i]) \
                               for i, v in enumerate(vars) \
                               if solution[i] > 0])
            if totalCut > 0:
                var_values[useVars[patt]] = 1
                cost = 1
            else:
                cost = 0

            var_tuple = (cost, rc, var_values)
            return [var_tuple]

        return []

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
                if zyes > zbest:
                    zbest = zyes
                    solbest = solyes

        return zbest, solbest

    return prob, relaxed_solver


if __name__ == '__main__':
    unittest.main()
