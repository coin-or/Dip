from __future__ import absolute_import
import unittest
from pulp import *
import coinor.dippy as dippy

from .dippy_tests import DippyTestCase


class TestCokeProblem(DippyTestCase):

    def setUp(self):
        """
        sets up the coke problem as a unittest
        """
        self.prob, self.do_branch = create_coke_problem()

    def test_pulp_solve(self):
        """
        tests that pulp can solve the problem
        """
        status = self.prob.solve()
        self.assertEqual(status, LpStatusOptimal)
        self.assertAlmostEqual(self.prob.objective.value(), 191078860.725, 2)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)

    def test_dippy_branch(self):
        """
        tests that the custom branch can solve the problem
        """
        self.prob.branch_method = self.do_branch
        dippy.Solve(self.prob, {})
        self.assertAlmostEqual(self.prob.objective.value(), 191078860.725, 2)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)



def create_coke_problem():
    """
    creates and returns the coke problem
    """

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

    MINE_TRANSPORT_DATA = """
            L1      L2      L3      L4      L5      L6
    M1      231737  46813   79337   195845  103445  45186
    M2      179622  267996  117602  200298  128184  49046
    M3      45170   93159   156241  218655  103802  119616
    M4      149925  254305  76423   123534  151784  104081
    M5      152301  205126  24321   66187   195559  88979
    M6      223934  132391  51004   122329  222927  54357
    """

    CUST_TRANSPORT_DATA = """
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

    MINE_TRANSPORT = read_table(MINE_TRANSPORT_DATA, int)
    CUST_TRANSPORT = read_table(CUST_TRANSPORT_DATA, int, \
                            transpose=True)

    #add both parts of the network together
    ARC_COSTS = MINE_TRANSPORT.copy()
    ARC_COSTS.update(CUST_TRANSPORT)
    ARCS = list(ARC_COSTS.keys())

    LOCATIONS_SIZES = [(l, s) for l in LOCATIONS
                        for s in SIZES]

    prob = dippy.DipProblem("Coke", LpMinimize)

    # create variables
    buildVars = LpVariable.dicts("Build", LOCATIONS_SIZES,
                                 cat = LpBinary)

    # create arcs
    flowVars = LpVariable.dicts("Arcs", ARCS, lowBound=0,
                                upBound = BIG_M)

    # objective
    prob += 1e6 * sum(buildVars[(l, s)] * SIZE_COSTS[s]
                      for (l, s) in LOCATIONS_SIZES) + \
            sum(flowVars[(s, d)] * ARC_COSTS[(s, d)]
                for (s, d) in ARCS), \
            "cost_of_building_and_transport"

    # plant availability
    for loc in LOCATIONS:
        prob += sum(flowVars[(loc, c)] for c in CUSTOMERS) \
                <= sum(buildVars[(loc, s)] * s
                       for s in SIZES), \
                       "Plant_%s_availability"%loc

    # one size
    for loc in LOCATIONS:
        prob += sum(buildVars[(loc, s)] for s in SIZES) == 1, \
                "Plant_%s_size"%loc

    # conserve flow (mines)
    # flows are in terms of tonnes of coke
    for m in MINES:
        prob += sum(flowVars[(m, j)] for j in LOCATIONS) \
                <= MINE_SUPPLY[m], "Supply_mine_%s"%m

    for loc in LOCATIONS:
        prob += sum(flowVars[(m, loc)] for m in MINES) - \
                CC * sum(flowVars[(loc, c)] for c in CUSTOMERS) \
                == 0, "Conserve_flow_location_%s"%loc

    for c in CUSTOMERS:
        prob += sum(flowVars[(loc, c)]  for loc in LOCATIONS) \
                >= CUSTOMER_DEMAND[c], "Demand_cust_%s"%c

    def do_branch(prob, sol):
        tol = 1e-10
        for loc in LOCATIONS:
            sol_size = sum(sol[buildVars[loc, size]] * \
                           size for size in SIZES)
            # determine if solsize is bigger than the largest
            # or smaller tham the smallest
            if (abs(sol_size - max(SIZES)) < tol
                    or abs(sol_size - min(SIZES))) < tol:
                continue
            #find the first one bigger or equal to sol_size
            bigger = min([s for s in SIZES
                         if s >= sol_size - tol])
            if bigger - sol_size > tol:
                down_branch_ub = dict([(buildVars[loc, s], 0)
                                       for s in SIZES
                                       if s <= sol_size])
                up_branch_ub = dict([(buildVars[loc, s], 0)
                                     for s in SIZES
                                     if s > sol_size])

                return {}, down_branch_ub, {}, up_branch_ub

    return prob, do_branch


if __name__ == '__main__':
    unittest.main()
