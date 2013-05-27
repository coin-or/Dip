"""
How to use this file:

from _dippy import *

_Solve = Solve
def Solve(prob, params=None):
  # Set up prob here, must have DipAPI as a base class

  # call the Solve method from _dippy
  try:
    solList, dualList = _Solve(prob, processed)
    # solList  is a list of (col_name, value) pairs
    # dualList is a list of (row_name, value) pairs
  except:
    print "Error returned from _dippy"
    raise
"""

class DipAPIError(Exception):
  """
  DipAPI Exception
  """
  pass

class DipAPI(object):

  def getObjective(self):
    """
    Return objective as a dictionary with variables as keys
    and (non-zero) coefficients as values
    """
    raise DipAPIError("Bad function definition, DipAPI.getObjective must be overwritten")

  def getRows(self, problem=None):
    """
    Return constraints as a list of dictionaries with variables as keys
    and (non-zero) coefficients as values. Constraints also have
    getName, getLb and getUb methods

    problem = None implies the master problem, otherwise problem
    is a subproblem
    """
    raise DipAPIError("Bad function definition, DipAPI.getRows must be overwritten")

  def getCols(self, problem=None):
    """
    Returns a list of variables. Variables have getName, getLb,
    getUb and isInteger methods

    problem = None implies the master problem, otherwise problem
    is a subproblem
    """
    raise DipAPIError("Bad function definition, DipAPI.getCols must be overwritten")

  def getMasterAsTuple(self):
    """
    Returns all the master problem data as a tuple of other
    "data gathering" functions
    """
    return (self.getObjective(),
            self.getRows(),
            self.getCols())

  def getRelaxAsTuple(self, problem):
    """
    Returns all the subproblem constraints and variables
    """
    return (self.getRows(problem),
            self.getCols(problem))

  def getRelaxsAsDict(self):
    """
    Returns the relaxation subproblems as a dictionary with keys as
    defined by the user and values as subproblems
    """
    raise DipAPIError("Bad function definition, DipAPI.getRelaxsAsDict must be overwritten")

  def chooseBranchSet(self, xhat):
    """
    Finds the best branch for a fractional solution

    Inputs:
    xhat (list of (variable, value) tuples) = list of solution values for all variables

    Output:
    down_lb, down_ub, up_lb, up_ub (tuple of (variable, value) dictionaries) =
    lower and upper bounds for down branch, lower and upper bounds for up branch
    """
    raise DipAPIError("Bad function definition, DipAPI.chooseBranchSet must be overwritten")

  def postProcessNode(self, output):
    """
    Returns information from the node that has just been processed.

    Inputs:
    output (list of (parameter, value) tuples) = list of output values from the node
    """
    raise DipAPIError("Bad function definition, DipAPI.postProcess must be overwritten")

  def solveRelaxed(self, key, redCostX, convexDual):
    """
    Returns solutions to the whichBlock relaxed subproblem

    Inputs:  
    key (Python Object) = key of relaxed subproblem to be solved
    redCostX (list of (variable, value) tuples) = list of reduced costs for all variables
    convexDual (float) = dual for convexity constraint for this relaxed subproblem

    Output:
    varList (list of (cost, reduced cost, list of (variable, value) dictionaries)) =
    solution for this relaxed subproblem expressed as a cost, reduced cost and
    dictionary of non-zero values for variables
    """
    raise DipAPIError("Bad function definition, DipAPI.solveRelaxed must be overwritten")

  def isUserFeasible(self, sol, tol):
    """
    Lets the user decide if an integer solution is really feasible

    Inputs:
    sol (list of (variable, value) tuples) = list of solution values for all variables
    tol (double) = zero tolerance

    Outputs:
    (boolean) = false if not feasible (generate cuts) or true if feasible
    """
    raise DipAPIError("Bad function definition, DipAPI.isUserFeasible must be overwritten")

  def generateCuts(self, node):
    """
    Lets the user generate cuts to remove fractional "pieces" of the node solution

    Inputs:
    node (list of (string, object) tuples) = list of node properties

    Output:
    cutList (list of LpConstraints) =
    cuts for this fractional solution expressed as a list LpConstraints,
    i.e., a dictionary with LpVariables as keys and (non-zero) coefficients
    as values with getName, getLb and getUb bound methods
    """
    raise DipAPIError("Bad function definition, DipAPI.generateCuts must be overwritten")

  def solveHeuristics(self, xhat, costX):
    """
    Lets the user generate (heuristic) solutions from a fractional solution

    Inputs:
    xhat  (list of (variable, value) tuples) = list of solution values for all variables
    costX (list of (variable, value) tuples) = list of costs for all variables

    Outputs:
    solList (list of (variable, value) dictionaries) =
    solutions found from this fractional solution expressed as a
    dictionary of non-zero values for variables
    """
    raise DipAPIError("Bad function definition, DipAPI.solveHeuristics must be overwritten")

  def generateInitVars(self):
    """
    Returns initial solutions to relaxed subproblems

    Inputs:
    None

    Output:
    varList (list of (subproblem key, (cost, list of (variable, value) dictionaries))) =
    initial solutions for the relaxed subproblems expressed as a cost and
    dictionary of non-zero values for variables
    """
    raise DipAPIError("Bad function definition, DipAPI.generateInitVars must be overwritten")
