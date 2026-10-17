"""
This is the testing framework used for the rest of the tests
"""
from __future__ import absolute_import
import unittest
from pulp import *
import coinor.dippy as dippy

class DippyTestCase(unittest.TestCase):

    def variable_feasibility_test(self, prob, tol = 1e-5):
        """
        tests that the problem is feasible for its class of variables
        """

        for var in prob.variables():
            self.assertTrue(var.valid(tol))
            self.assertLess(var.value(), var.upBound + tol)
            self.assertGreater(var.value(), var.lowBound - tol)
            if var.cat == LpInteger:
                self.assertAlmostEqual(int(var.value() + tol), var.value())

    def constraint_feasibility_test(self, prob, tol = 1e-5):
        """
        tests that the problem is feasible for its class of variables
        """

        for constraint_name, constraint in list(prob.constraints.items()):
            self.assertTrue(constraint.valid(tol))

if __name__ == '__main__':
    from .test_coke import *
    from .test_facility import *
    from .test_cutting_stock import *
    from .test_tsp import *
    unittest.main()
