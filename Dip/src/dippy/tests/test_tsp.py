from __future__ import absolute_import
from builtins import range
import unittest
from pulp import *
import coinor.dippy as dippy

from .dippy_tests import DippyTestCase


class TestTSPProblem(DippyTestCase):

    def setUp(self):
        self.prob, self.generate_cuts, self.is_solution_feasible = \
                    create_tsp_problem()
        self.tol = 1e-6


    def test_pulp_solve(self):
        """
        tests that pulp can solve the problem
        """
        status = self.prob.solve()
        self.assertEqual(status, LpStatusOptimal)
        self.assertAlmostEqual(self.prob.objective.value(), 14.24, 2)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)
        # would like this in here but not sure how to do it
        #~ self.assertFalse(self.is_solution_feasible())

    def test_dippy_solve(self):
        """
        tests that dippy can solve the problem
        """
        self.prob.generate_cuts = self.generate_cuts
        self.prob.is_solution_feasible = self.is_solution_feasible
        sol, duals = dippy.Solve(self.prob, {})
        self.assertAlmostEqual(self.prob.objective.value(), 18.39, 2)
        self.variable_feasibility_test(self.prob)
        self.constraint_feasibility_test(self.prob)
        # would like this in here but not sure how to do it
        self.assertTrue(self.is_solution_feasible(self.prob, sol, self.tol))

def create_tsp_problem(doRelaxed=False):
    """
    creates and returns the tsp problem
    """

    from math import sqrt

    # 2d Euclidean TSP with extremely simple cut generation

    # x,y coords of cities
    CITY_LOCS = [(0, 2), (0, 4), (1, 2), (1, 4), \
                 (4, 1), (4, 4), (4, 5),  (5, 0), \
                 (5, 2), (5, 5)]
    CITIES = list(range(len(CITY_LOCS)))

    ARCS = [] # list of arcs (no duplicates)
    ARC_COSTS = {} # distance

    # for each city, list of arcs into/out of
    CITY_ARCS = [[] for i in CITIES]

    # use 2d euclidean distance
    def dist(x1, y1, x2, y2):
        return sqrt((x1-x2)**2 + (y1-y2)**2)

    # construct list of arcs
    for i in CITIES:
        i_x, i_y = CITY_LOCS[i]
        for j in CITIES[i+1:]:
            j_x, j_y = CITY_LOCS[j]
            ARC_COSTS[(i,j)] = dist(i_x, i_y, j_x, j_y)
            ARCS.append((i, j))
            CITY_ARCS[i].append((i, j))
            CITY_ARCS[j].append((i, j))

    prob = dippy.DipProblem()

    arc_vars = LpVariable.dicts("UseArc", ARCS, 0, 1, LpBinary)

    # objective
    prob += sum(ARC_COSTS[x] * arc_vars[x] for x in ARCS)

    # degree constraints
    for city in CITIES:
        prob += sum(arc_vars[x] for x in CITY_ARCS[city]) == 2

    # generate subtour elimination constraints

    # dictionary for symmetric arcs, can be
    # accessed using (i, j) and (j, i)
    symmetric = {}
    for i in CITIES:
        for j in CITIES:
            if i < j:
                symmetric[(i, j)] = (i, j)
                symmetric[(j, i)] = (i, j)

    def generate_cuts(prob, sol):
        cons = []
        not_connected = set(CITIES)

        while not_connected:
     #       print "not_connected", [n for n in not_connected]
            start = not_connected.pop()
            nodes, arcs = get_subtour(sol, start)
            if len(nodes) == len(arcs) and \
               len(nodes) < len(CITIES):
                cons.append( sum(arc_vars[a] for a in arcs) \
                             <= len(arcs) - 1 )
     #       print "nodes", [n for n in nodes]
     #       print "arcs", [a for a in arcs]
            nodes.remove(start)
            not_connected -= nodes

    ##    print "# cons = ", len(cons)
    ##    print "cons = ", [ con for con in cons ]
        return cons

    def is_solution_feasible(prob, sol, tol):
        nodes, arcs = get_subtour(sol, 0)
    #    print "Checking: # nodes = ", len(nodes), \
    #          ", # cities = ", len(CITIES)
    #    print "nodes", [n for n in nodes]

        return len(nodes) == len(arcs) and \
               len(nodes) == len(CITIES)

    def get_subtour(sol, node):
        # returns: list of nodes and arcs
        # in subtour containing node
        nodes = set([node])
        arcs = set([])
        not_processed = set(CITIES)
        to_process = set([node])

        tol = 1e-6
        one = 1 - tol

        while to_process:
            c = to_process.pop()
            not_processed.remove(c)
            new_arcs = [ symmetric[(c, i)] \
                         for i in not_processed \
                         if sol[ \
                             arc_vars[symmetric[(c, i)]]]
                         > one]
            new_nodes = [ i for i in not_processed \
                          if symmetric[(i, c)] in new_arcs ]
     #       print "new_nodes", [n for n in new_nodes]
     #       print "new_arcs", [a for a in new_arcs]
            arcs |= set(new_arcs)
            nodes |= set(new_nodes)
            to_process |= set(new_nodes)
     #       print "not_processed", [n for n in not_processed]
     #       print "nodes", [n for n in nodes]
     #       print "arcs", [a for a in arcs]

        return nodes, arcs

    return prob, generate_cuts, is_solution_feasible


if __name__ == '__main__':
    unittest.main()
