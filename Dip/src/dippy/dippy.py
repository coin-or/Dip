from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from past.builtins import basestring
from builtins import object
import pulp

import sys

try:
    import path.py
except ImportError:
    pass

from .dipapi import DipAPI

from ._dippy import *

gimpy_installed = True
try:
    from src.gimpy import BinaryTree
except ImportError:
    try:
        from coinor.gimpy import BinaryTree
    except ImportError:
        gimpy_installed = False

if gimpy_installed:
    try:
        from coinor.grumpy import BBTree
    except ImportError:
        grumpy_installed = False
    else:
        grumpy_installed = True
else:
    grumpy_installed = False

class DipError(Exception):
    """
    Dip Exception
    """
	
# DIP solver status
#enum DecompSolverStatus {
#   DecompSolStatError,
#   DecompSolStatOptimal,
#   DecompSolStatFeasible,
#   DecompSolStatInfeasible,
#   DecompSolStatNoSolution
#};
DipSolStatError      = 0
DipSolStatOptimal    = 1
DipSolStatFeasible   = 2
DipSolStatInfeasible = 3
DipSolStatNoSolution = 4
DipStatus = {
             DipSolStatError:      "Error",
             DipSolStatOptimal:    "Optimal",
             DipSolStatFeasible:   "Feasible",
             DipSolStatInfeasible: "Infeasible",
             DipSolStatNoSolution: "No solution"
            }

_Solve = Solve
def Solve(prob, params=None):
    """
    Solve a DipProblem instance, returning a solution object

    @param prob A DipProblem instance to solve
    @param params A dictionary of parameters to pass to DIP
    """

    # params is a dictionary, keys are strings, values are
    # strings or dictionaries

    # if value is a dictionary, key is a section, and items
    # of dictionary are names/values

    # if value is a string, then section is NULL and key is name
    # for these parameters we also assign them to the 'DECOMP'
    # section as a convenience

    # the dictionary is converted into a dictionary of
    # strings indexed by (section, name) tuples

    processed = {}
    if params is None:
        params = {}

    if (prob.branch_method == None) and (prob.display_mode == 'off'):
        params['pyBranchMethod'] = '0'
    if (prob.post_process_node == None) and (prob.display_mode == 'off'):
        params['pyPostProcessNode'] = '0'
    if (prob.post_process_branch == None) and (prob.display_mode == 'off'):
        params['pyPostProcessBranch'] = '0'
    if prob.relaxed_solver == None:
        params['pyRelaxedSolver'] = '0'
    if prob.is_solution_feasible == None:
        params['pyIsSolutionFeasible'] = '0'
    
    if (prob.generate_cuts == None) and (prob.generate_cuts_from_node == None):
        params['pyGenerateCuts'] = '0'

    if prob.generate_cuts != None:
        prob.gen_cuts = True
    else:
        prob.gen_cuts = False
    if prob.generate_cuts_from_node != None:
        prob.gen_cuts_node = True
    else:
        prob.gen_cuts_node = False
    
    if prob.heuristics == None:
        params['pyHeuristics'] = '0'
    if prob.init_vars == None:
        params['pyInitVars'] = '0'
    if prob.is_solution_feasible == None:
        params['pyIsSolutionFeasible'] = '0'

    for key, value in list(params.items()):
        valid_types = (basestring, int, float)
        if not isinstance(key, basestring):
            raise DipError("Bad key in parameter dictionary, expecting string")
        if isinstance(value, dict):
            section = key
            for name, param_val in list(value.items()):
                if not isinstance(param_val, valid_types):
                    raise DipError("Bad value '%s' in parameter dictionary, expecting string or number" % param_val)
                processed[(section, name)] = str(param_val)
        elif isinstance(value, valid_types):
            # add this parameter to both the 'None' section and the 'DECOMP' section
            processed[(None, key)] = str(value)
            processed[('DECOMP', key)] = str(value)
        else:
            raise DipError("Bad value '%s' in parameter dictionary, expecting string" % value)

    # DIP only solves minimisation problems
    if prob.sense == pulp.LpMaximize:
        raise DipError("DIP assumes a minimize objective, but DipProblem has "+ 
                       "maximize objective.\n" +
                       "Use prob.sense = pulp.LpMinimize and prob.objective " +
                       "*= -1 to remedy")
    
    # DIP only allows non-negative variables. This is difficult
    # to transform automatically, so request re-formulation
    for v in prob.variables():
        if v.lowBound < 0:
            raise DipError("Variable %s has negative lower bound, please " +
                           "re-formulate using sum of non-negative variables" 
                           % v.name)
        
    # call the Solve method from _dippy
    try:
        status, message, solList, dualList = _Solve(prob, processed)
        # solList  is a list of (col_name, value) pairs
        # dualList is a list of (row_name, value) pairs

        if prob.display_mode == 'svg' and gimpy_installed:
            if prob.display_interval is not None:
                prob.Tree.write_as_svg(filename = "%s_%d" % (prob.svg_prefix, 
                                                             prob.last_svg + 1),
                                       prevfile = "%s_%d" % (prob.svg_prefix, 
                                                             prob.last_svg))
                prob.last_svg += 1
            
    except Exception as ex:
        print("Error returned from _dippy")
        print(ex)
        raise

    if solList is None:
        solution = None
    else:
        solDict = dict(solList)
        setVars = set(prob.variables())
        setSolVars = set(solDict.keys())
        diff = setVars.symmetric_difference(setSolVars)
        if len(diff) > 0:
            raise DipError("Solution and variable list don't match in " +
                           "dippy.Solve")

        solution = solDict
        for v in prob.variables():
            v.varValue = solution[v]
        
    if dualList is None:
        duals = None
    else:
        dualDict = dict([(c.getName(), v) for (c, v) in dualList])
        setCons = set(prob.constraints)
        setDualCons = set(dualDict.keys())
        diff = setCons.symmetric_difference(setDualCons)
        if len(diff) > 0:
            raise DipError("Duals and constraint list don't match in dippy.Solve")

        duals = dualDict

    # return status, message, solution and duals
    return status, message, solution, duals

def createBranchLabel(lbs, ubs):
    maxLabelWidth = 20
      
    both = set(lbs.keys()) & set(ubs.keys())
    lbOnly = set(lbs.keys()) - both
    ubOnly = set(ubs.keys()) - both
                        
    bothStr = ''
    first = True
    currWidth = len(bothStr)
    for l in both:
        addStr = str(int(lbs[l])) + '<=' + \
                 ''.join(x for x in str(l) 
                         if x not in '()_') \
                 + '<=' + str(int(ubs[l]))
        if first:
            newWidth = len(addStr) + currWidth
        else:
            newWidth = len(addStr) + currWidth + 2
        if newWidth <= maxLabelWidth:
            if first:
                first = False
            else:
                bothStr += ', '
                currWidth += 2
            bothStr += addStr
            currWidth += len(addStr)
        else:
            bothStr += '\n' + addStr
            currWidth = len(addStr)
#            bothStr += ' ...'
#            break

    lbOnlyStr = ''
    first = True
    currWidth = len(lbOnlyStr)
    for l in lbOnly:
        addStr = str(int(lbs[l])) + '<=' + \
                 ''.join(x for x in str(l) 
                         if x not in '()_')
        if first:
            newWidth = len(addStr) + currWidth
        else:
            newWidth = len(addStr) + currWidth + 2
        if newWidth <= maxLabelWidth:
            if first:
                first = False
            else:
                lbOnlyStr += ', '
                currWidth += 2
            lbOnlyStr += addStr
            currWidth += len(addStr)
        else:
            lbOnlyStr += '\n' + addStr
            currWidth = len(addStr)
#            lbOnlyStr += ' ...'
#            break

    ubOnlyStr = ''
    first = True
    currWidth = len(ubOnlyStr)
    for l in ubOnly:
        addStr = ''.join(x for x in str(l) 
                         if x not in '()_') \
                 + '<=' + str(int(ubs[l]))
        if first:
            newWidth = len(addStr) + currWidth
        else:
            newWidth = len(addStr) + currWidth + 2
        if newWidth <= maxLabelWidth:
            if first:
                first = False
            else:
                ubOnlyStr += ', '
                currWidth += 2
            ubOnlyStr += addStr
            currWidth += len(addStr)
        else:
            ubOnlyStr += '\n' + addStr
            currWidth = len(addStr)
#            ubOnlyStr += ' ...'
#            break

    labelStr = ''
    first = True
    if len(bothStr) > 0:
        labelStr += bothStr
        first = False
    if len(lbOnlyStr) > 0:
        if first:
            first = False
        else:
            labelStr += '\n'
        labelStr += lbOnlyStr
    if len(ubOnlyStr) > 0:
        if first:
            first = False
        else:
            labelStr += '\n'
        labelStr += ubOnlyStr
        
    return labelStr

import string
def asCplexName(name):
    #to remove illegal characters from the names
    trans = str.maketrans("-+[] ->/","________")
    
    return str(name).translate(trans)

class DipProblem(pulp.LpProblem, DipAPI):

    def __init__(self, *args, **kwargs):
        # callback functions can be passed to class constructor as keyword 
        # arguments
        self.branch_method = kwargs.pop('branch_method', None)
        self.post_process_branch = kwargs.pop('post_process_branch', None)
        self.post_process_node = kwargs.pop('post_process_node', None)
        self.relaxed_solver = kwargs.pop('relaxed_solver', None)
        self.is_solution_feasible = kwargs.pop('is_solution_feasible', None)
        self.generate_cuts = kwargs.pop('generate_cuts', None)
        self.generate_cuts_from_node = kwargs.pop('generate_cuts_from_node', 
                                                  None)
        self.heuristics = kwargs.pop('heuristics', None)
        self.init_vars = kwargs.pop('init_vars', None)
        self.display_mode = kwargs.pop('display_mode', 'off')
        self.display_interval = kwargs.pop('display_interval', 1)
        self.layout = kwargs.pop('layout', 'dot')
        self.svg_prefix = kwargs.pop('svg_prefix', 'tree')
        
        if self.display_mode != 'off':
            if not gimpy_installed:
                print("GiMPy not installed. Display mode set to 'off'")
                self.display_mode = 'off'
            else:
                if grumpy_installed:
                    self.Tree = BBTree()
                else:
                    if self.layout == 'bak':
                        print("GrUMPy not installed. Display mode set to 'off'")
                        self.display_mode = 'off'
                    else:
                        self.Tree = BinaryTree()
        if self.display_mode != 'off':
            self.Tree.set_display_mode(self.display_mode)
            self.Tree.set_layout(self.layout)

        super(DipProblem, self).__init__(*args, **kwargs)
        self._subproblem = []
        self.relaxation = RelaxationCollection(self)

    def deepcopy(self):
        # callback functions can be passed to class constructor as keyword 
        # arguments
        dipcopy = DipProblem(name = self.name, sense = self.sense)
        dipcopy.branch_method = self.branch_method
        dipcopy.is_solution_feasible = self.is_solution_feasible
        dipcopy.generate_cuts = self.generate_cuts
        dipcopy.heuristics = self.heuristics
        dipcopy.init_vars = self.init_vars

        # This code is taken from pulp.py and needs to be coordinated
        # with pulp.py to avoid errors
        if dipcopy.objective != None:
            dipcopy.objective = self.objective.copy()
        dipcopy.constraints = {}
        for k,v in self.constraints.items():
            dipcopy.constraints[k] = v.copy()
        dipcopy.sos1 = self.sos1.copy()
        dipcopy.sos2 = self.sos2.copy()

        dipcopy._subproblem = self._subproblem[:]
        for k in list(self.relaxation.keys()):
            dipcopy.relaxation[k] = self.relaxation[k].copy()

        return dipcopy
    
    def variables(self):
        """
        Returns a list of the problem variables
        Overrides LpProblem.variables()
        
        Inputs:
            - none
        
        Returns:
            - A list of the problem variables
        """
        variables = {}
        if self.objective:
            variables.update(self.objective)
        for c in self.constraints.values():
            variables.update(c)
        for p in sorted(self.relaxation.keys()):
            for c in self.relaxation[p].constraints.values():
                variables.update(c)
        variables = list(variables)
        variables = sorted(variables, key=lambda variable: variable.name)
        
        return variables

    def getObjective(self):
        """
        Return objective as a dictionary with LpVariables as keys
        and (non-zero) coefficients as values
        """
        return self.objective

    def getRows(self, problem=None):
        """
        Return constraints as a list of dictionaries with LpVariables as keys
        and (non-zero) coefficients as values. Constraints also have
        getName, getLb and getUb methods (i.e., a LpConstraint)

        problem = None implies the master problem, otherwise problem
        is a subproblem
        """

        if problem is None:
            problem = self

        for n, c in problem.constraints.items():
            if c.name == None:
                c.name = n
        constraints = list(problem.constraints.values())

        return constraints

    def getCols(self, problem=None):
        """
        Returns a list of variables. Variables have getName, getLb,
        getUb and isInteger methods

        problem = None implies the master problem, otherwise problem
        is a subproblem
        """

        if problem is None:
            variables = self.variables()
        else:
            variables = {}
            for c in problem.constraints.values():
                variables.update(c)
            variables = list(variables)
            variables = sorted(variables, key=lambda variable: variable.name)

        return variables

    def getRelaxsAsDict(self):
        """
        Returns the relaxation subproblems as a dictionary with keys as
        defined by the user and values as subproblems
        """
        return self.relaxation.dict

    def writeFull(self, instancefile, blockfile, mip = True):
        f = open(instancefile, "w")
        b = open(blockfile, "w")
        f.write("\\* "+self.name+" *\\\n")
        if self.sense == 1:
            f.write("Minimize\n")
        else:
            f.write("Maximize\n")
        wasNone, dummyVar = self.fixObjective()
        objName = self.objective.name
        if not objName: objName = "OBJ"
        f.write(self.objective.asCplexLpAffineExpression(objName, constant = 0))
        f.write("Subject To\n")
        b.write("NBLOCKS\n")
        b.write("%i\n" % len(self.relaxation.dict))
        for k in self.constraints:
            f.write(self.constraints[k].asCplexLpConstraint(k))
        blockId = 0
        for r in self.relaxation.dict:
            rname = asCplexName(str(r))
            b.write("BLOCK %d\n" % blockId) 
            for k in self.relaxation.dict[r].constraints:
                f.write(self.relaxation.dict[r].constraints[k].asCplexLpConstraint(str(k)+'_'+rname))
                b.write(str(k)+'_'+rname+'\n')
            blockId += 1
        vs = list(self.variables())
        # check if any names are longer than 100 characters
        long_names = [v.name for v in vs if len(v.name) > 100]
        if long_names:
            raise PulpError('Variable names too long for Lp format\n' 
                                + str(long_names))
        # check for repeated names
        repeated_names = {}
        for v in vs:
            repeated_names[v.name] = repeated_names.get(v.name, 0) + 1
        repeated_names = [(key, value) for key, value in list(repeated_names.items())
                            if value >= 2]
        if repeated_names:
            raise PulpError('Repeated variable names in Lp format\n' 
                                + str(repeated_names))
        # Bounds on non-"positive" variables
        # Note: XPRESS and CPLEX do not interpret integer variables without 
        # explicit bounds
        if mip:
            vg = [v for v in vs if not (v.isPositive() and v.cat == pulp.LpContinuous) \
                and not v.isBinary()]
        else:
            vg = [v for v in vs if not v.isPositive()]
        if vg:
            f.write("Bounds\n")
            for v in vg:
                f.write("%s\n" % v.asCplexLpVariable())
        # Integer non-binary variables
        if mip:
            vg = [v for v in vs if v.cat == pulp.LpInteger and not v.isBinary()]
            if vg:
                f.write("Generals\n")
                for v in vg: f.write("%s\n" % v.name)
            # Binary variables
            vg = [v for v in vs if v.isBinary()]
            if vg:
                f.write("Binaries\n")
                for v in vg: f.write("%s\n" % v.name)
        f.write("End\n")
        f.close()
        self.restoreObjective(wasNone, dummyVar)

    def writeRelaxed(self, block, filename, mip = True):
        """
        Write the given block into a .lp file.
        
        This function writes the specifications (NO objective function,
        constraints, variables) of the defined Lp problem to a file.
        
        Inputs:
            - block -- the key to the block to write
            - filename -- the name of the file to be created.          
                
        Side Effects:
            - The file is created.
        """
        relaxation = self.relaxation[block]
        f = open(filename, "w")
        f.write("\\* "+relaxation.name+" *\\\n")
        f.write("Subject To\n")
        ks = list(relaxation.constraints.keys())
        ks.sort()
        for k in ks:
            f.write(relaxation.constraints[k].asCplexLpConstraint(k))
        vs = relaxation.variables()
        # check for repeated names
        relaxation.checkDuplicateVars()
        # Bounds on non-"positive" variables
        # Note: XPRESS and CPLEX do not interpret integer variables without 
        # explicit bounds
        if mip:
            vg = [v for v in vs if not (v.isPositive() and \
                                        v.cat == pulp.LpContinuous) \
                and not v.isBinary()]
        else:
            vg = [v for v in vs if not v.isPositive()]
        if vg:
            f.write("Bounds\n")
            for v in vg:
                f.write("%s\n" % v.asCplexLpVariable())
        # Integer non-binary variables
        if mip:
            vg = [v for v in vs if v.cat == pulp.LpInteger and \
                                   not v.isBinary()]
            if vg:
                f.write("Generals\n")
                for v in vg: f.write("%s\n" % v.name)
            # Binary variables
            vg = [v for v in vs if v.isBinary()]
            if vg:
                f.write("Binaries\n")
                for v in vg: f.write("%s\n" % v.name)
        f.write("End\n")
        f.close()
        
    def chooseBranchSet(self, xhat):
        """
        Finds the best branch for a fractional solution

        Inputs:
        xhat (list of (LpVariable, value) tuples) = list of solution values for all variables

        Output:
        down_lb, down_ub, up_lb, up_ub (tuple of (LpVariable, value) dictionaries) =
        lower and upper bounds for down branch, lower and upper bounds for up branch
        """
        try:
            if self.branch_method is None:
                return None
                        
            xhatDict = dict(xhat)
            setVars = set(self.variables())
            setXhatVars = set(xhatDict.keys())
            diff = setVars.symmetric_difference(setXhatVars)
            if len(diff) > 0:
                raise DipError("Solution and variable list don't match in chooseBranchSet")
    
            branch_sets = self.branch_method(self, xhatDict)

            if branch_sets is None:
                return None
                            
            if (branch_sets[0] or branch_sets[1]) and (branch_sets[2] or branch_sets[3]):
                return branch_sets
            else:
                raise DipError("Invalid bounds returned from user-specified branch_method")
        except Exception as ex:
            errorStr = "Error in chooseBranchSet\n%s" % ex
            raise DipError(errorStr)
          
    def decipherNode(self, output):
        outputDict = dict(output)
        if "xhat" in list(outputDict.keys()):
            xhat = outputDict["xhat"]
            outputDict["xhat"] = dict(xhat)
        if "bounds" in list(outputDict.keys()):
            bounds = outputDict["bounds"]
            outputDict["bounds"] = dict(bounds)
            
        return outputDict

    def postProcessNode(self, node):
        """
        Returns information from the node that has just been processed.

        Inputs:
        output (list of (parameter, value) tuples) = list of output values 
        from the node
        """
        try:        
            nodeDict = self.decipherNode(node)
            
            if gimpy_installed:
                nodeInd = nodeDict["nodeIndex"]
                parentInd = nodeDict["parentIndex"]
                nodeQuality = nodeDict["nodeQuality"]
                branchedDir = nodeDict["branchedDir"]
                nodeStatus = nodeDict["nodeStatus"]
                if branchedDir == -1:
                    branch_direction = 'L'
                else:
                    branch_direction = 'R'
                    
                if nodeStatus == 'Infeasible':
                    status = 'I'
                    BAKstatus = 'infeasible'
                    color = 'orange'
                elif nodeStatus == 'Candidate':
                    status = 'C'
                    BAKstatus = 'candidate'
                    color = 'yellow'
                elif nodeStatus == 'Solution':
                    status = 'S'
                    BAKstatus = 'integer'
                    color = 'lightblue'
                else:
                    status = 'P'
                    BAKstatus = 'fathomed'
                    color = 'red'
                    
                if nodeStatus != 'Infeasible':
                    label = status + ": " + "%.1f"%nodeQuality
                else:
                    label = 'I'

                numNodes = len(self.Tree.get_node_list())
                if parentInd == -1:
                    if self.layout == 'bak':
                        self.Tree.AddOrUpdateNode(nodeInd, parentInd,
                                                  branch_direction, BAKstatus, 
                                                  nodeQuality, None, None)
                    else:
                        self.Tree.add_root(nodeInd, label = label, 
                                           status = 'C', obj = nodeQuality, 
                                           color = color, style = 'filled', 
                                           fillcolor = color)
                    if self.Tree.attr['display'] == 'svg':
                        if self.display_interval is not None:
                            if numNodes % self.display_interval in [0, 1]:
                                self.Tree.write_as_svg(filename = "%s_0" 
                                                       % self.svg_prefix, 
                                                       nextfile = "%s_1" 
                                                       % self.svg_prefix, 
                                                       highlight = nodeInd)
                                self.last_svg = 0
                    numNodes += 1
                else:
                    if branch_direction == 'L':
                        n = self.Tree.get_left_child(parentInd)
                    else:
                        n = self.Tree.get_right_child(parentInd)
                    edge_label = self.Tree.get_edge_attr(parentInd, n, 'label')
                    self.Tree.del_node(n)
                    if self.layout == 'bak':
                        self.Tree.AddOrUpdateNode(nodeInd, parentInd, 
                                                  branch_direction, 'branched', 
                                                  nodeQuality, None, None)
                    elif branch_direction == 'L':
                        self.Tree.add_left_child(nodeInd, parentInd, 
                                                 label = label, 
                                                 status = status, 
                                                 obj = nodeQuality, 
                                                 color = color, 
                                                 style = 'filled', 
                                                 fillcolor = color)
                    else:
                        self.Tree.add_right_child(nodeInd, parentInd, 
                                                  label = label, 
                                                  status = status, 
                                                  obj = nodeQuality, 
                                                  color = color, 
                                                  style = 'filled', 
                                                  fillcolor = color)
                    if edge_label is not None:
                        self.Tree.set_edge_attr(parentInd, nodeInd, 
                                                'label', edge_label)
                    if self.Tree.attr['display'] == 'svg':
                        if self.display_interval is not None:
                            if numNodes % self.display_interval in [0, 1]:
                                self.Tree.write_as_svg(filename = "%s_%d" 
                                                       % (self.svg_prefix, 
                                                          self.last_svg + 1), 
                                                       prevfile = "%s_%d" 
                                                       % (self.svg_prefix, 
                                                          self.last_svg), 
                                                       nextfile = "%s_%d" 
                                                       % (self.svg_prefix, 
                                                          self.last_svg + 2), 
                                                       highlight = nodeInd)
                                self.last_svg += 1
                if self.display_interval is not None:
                    if numNodes % self.display_interval in [0, 1]:
                        self.Tree.display()

            if self.post_process_node is not None:
                self.post_process_node(self, nodeDict)
          
        except Exception as ex:
            errorStr = "Error in postProcessNode\n%s" % ex
            raise DipError(errorStr)

    def postProcessBranch(self, branchInfo):
        """
        Returns information from the node that has just been processed.

        Inputs:
        output (list of (parameter, value) tuples) describing branching decision
        """
        try:
          
            outputDict = dict(branchInfo)

            if gimpy_installed:
                nodeInd = outputDict['nodeIndex']
                nodeQuality = outputDict['nodeQuality']
                numNodes = len(self.Tree.get_node_list())
                for n in outputDict:
                    if n == 'pDownUB':
                        if self.layout == 'bak':
                            self.Tree.AddOrUpdateNode(-numNodes, 
                                                      nodeInd, 'L', 
                                                      'candidate', 
                                                      nodeQuality, None, None)
                        else:
                            self.Tree.add_left_child(-numNodes, 
                                                     nodeInd, 
                                                     label = 'C', 
                                                     status = 'C', 
                                                     obj = nodeQuality, 
                                                     color = 'yellow', 
                                                     style = 'filled', 
                                                     fillcolor = 'yellow')
                        if 'pDownLB' in outputDict:
                            lbs = outputDict['pDownLB']
                        else:
                            lbs = {}
                        ubs = outputDict['pDownUB']
                        labelStr = createBranchLabel(lbs, ubs)
                        self.Tree.set_edge_attr(nodeInd, 
                                                -numNodes, 
                                                'label', labelStr)
                        numNodes += 1

                    elif n == 'pUpLB':
                        if self.layout == 'bak':
                            self.Tree.AddOrUpdateNode(-numNodes, 
                                                      nodeInd, 'R', 
                                                      'candidate', 
                                                      nodeQuality, None, None)
                        else:
                            self.Tree.add_right_child(-numNodes, 
                                                      nodeInd, 
                                                      label = 'C', 
                                                      status = 'C', 
                                                      obj = nodeQuality, 
                                                      color = 'yellow', 
                                                      style = 'filled', 
                                                      fillcolor = 'yellow')
                        if 'pUpUB' in outputDict:
                            ubs = outputDict['pUpUB']
                        else:
                            ubs = {}
                        lbs = outputDict['pUpLB']
                        labelStr = createBranchLabel(lbs, ubs)
                        self.Tree.set_edge_attr(nodeInd, 
                                                -numNodes, 
                                                'label', labelStr)
                        numNodes += 1

                if self.Tree.get_node_attr(nodeInd, 'color') == 'yellow':
                    self.Tree.set_node_attr(nodeInd, 'color', 'green')
                    self.Tree.set_node_attr(nodeInd, 'fillcolor', 'green')
        
            if self.post_process_branch is not None:
                self.post_process_branch(self, outputDict)
  
        except Exception as ex:
            errorStr = "Error in postProcessBranch\n%s" % ex
            raise DipError(errorStr)

    def solveRelaxed(self, key, redCostX, target):
        """
        Returns solutions to the whichBlock relaxed subproblem
    
        Inputs:  
        key (Python Object) = key of relaxed subproblem to be solved
        redCostX (list of (variable, value) tuples) = list of reduced costs for all variables
        target (float) = any total reduced cost less than the target is "good" (results in a negative cost column)
    
        Output:
        status (integer) = status of the relaxation solve, will be one of
          DipSolStatOptimal    = no better columns can be found
          DipSolStatFeasible   = better columns can be found, but just use these for now
          DipSolStatNoSolution = use any columns returned, but also use DIP's default column finder
        varList (list of (cost, reduced cost, list of (variable, value) dictionaries)) =
        solution for this relaxed subproblem expressed as a cost, reduced cost and
        dictionary of non-zero values for variables
        """
        try:
          
            # transform redCostX into a dictionary
            redCostDict = dict(redCostX)
            setVars = set(self.variables())
            setRedCostVars = set(redCostDict.keys())
            diff = setVars.symmetric_difference(setRedCostVars)
            if len(diff) > 0:
                print(diff)
                raise DipError("Reduced cost and variable list don't match in",
                               "solveRelaxed")
    
            status, dvs = self.relaxed_solver(self, key, redCostDict, target)
    
            if len(dvs) > 0:
                dvs_with_costs = []
                for var in dvs:
                    if isinstance(var, dict):
                        cost = sum(self.objective[i]*var[i] for i in var
                                   if i in self.objective)
                        red_cost = sum(redCostDict[i]*var[i] for i in var
                                       if i in redCostDict)
                        dvs_with_costs.append((cost, red_cost, var))
                    else:
                        return status, dvs
                return status, dvs_with_costs
            else:
                return status, dvs
  
        except Exception as ex:
            errorStr = "Error in solveRelaxed\n%s" % ex
            raise DipError(errorStr)

    def isUserFeasible(self, sol, tol):
        """
        Lets the user decide if an integer solution is really feasible

        Inputs:
        sol (list of (LpVariable, value) tuples) = list of solution
        values for all variables tol (double) = zero tolerance

        Outputs:
        (boolean) = false if not feasible (generate cuts) or true if feasible
        """
        try:
              
            solDict = dict(sol)
            setVars = set(self.variables())
            setSolVars = set(solDict.keys())
            diff = setVars.symmetric_difference(setSolVars)
            if len(diff) > 0:
                raise DipError("Solution and variable list don't match in isUserFeasible")
    
            return self.is_solution_feasible(self, solDict, tol)
          
        except Exception as ex:
            errorStr = "Error in isUserFeasible\n%s" % ex
            raise DipError(errorStr)

    def generateCuts(self, node):
        """
        Lets the user generate cuts to remove fractional "pieces" of xhat

        Inputs:
        node (list of (string, object) tuples) = list of node properties

        Output:
        cutList (list of LpConstraints) =
        cuts for this fractional solution expressed as a list LpConstraints,
        i.e., a dictionary with LpVariables as keys and (non-zero) coefficients
        as values with getName, getLb and getUb bound methods
        """
        try:
    
            nodeDict = self.decipherNode(node)
            xhatDict = nodeDict["xhat"]
            setVars = set(self.variables())
            setXhatVars = set(xhatDict.keys())
            diff = setVars.symmetric_difference(setXhatVars)
            if len(diff) > 0:
                raise DipError("Solution and variable list don't match in generateCuts")
    
            # Generate a list of cuts as LpConstraints
            if self.gen_cuts:
                cuts = self.generate_cuts(self, xhatDict)
            else:
                cuts = None
    
            if self.gen_cuts_node:
                moreCuts = self.generate_cuts_from_node(self, nodeDict)
                if moreCuts is not None:
                    if cuts is None:
                        cuts = moreCuts
                else:
                    cuts.extend(moreCuts)
                    
            if cuts is not None:
                if len(cuts) > 0:
                    return cuts
                else:
                    print("Empty cut list in generateCuts, returning None")

        except Exception as ex:
            errorStr = "Error in generateCuts\n%s" % ex
            raise DipError(errorStr)

    def solveHeuristics(self, xhat, costX):
        """
        Lets the user generate (heuristic) solutions from a fractional solution

        Inputs:
        xhat  (list of (LpVariable, value) tuples) = list of solution values for all variables
        costX (list of (LpVariable, value) tuples) = list of costs for all variables

        Outputs:
        solList (list of (LpVariable, value) dictionaries) =
        solutions found from this fractional solution expressed as a
        dictionary of non-zero values for variables
        """
        try:
            
            # transform xhat into a dictionary
            xhatDict = dict(xhat)
            setVars = set(self.variables())
            setXhatVars = set(xhatDict.keys())
            diff = setVars.symmetric_difference(setXhatVars)
            if len(diff) > 0:
                raise DipError("Solution and variable list don't match in solveHeuristics")
    
            # transform costs into a dictionary
            costDict = dict(costX)
            setCostVars = set(costDict.keys())
            diff = setVars.symmetric_difference(setCostVars)
            if len(diff) > 0:
                raise DipError("Cost and variable list don't match in solveHeuristics")
    
            sols = self.heuristics(self, xhatDict, costDict)
            if sols is not None:
                if len(sols) > 0:
                    return sols
                else:
                    print("Empty solution list in solveHeuristics, returning None")

        except Exception as ex:
            errorStr = "Error in solveHeuristics\n%s" % ex
            raise DipError(errorStr)

    def generateInitVars(self):
        """
        Returns initial solutions to relaxed subproblems

        Inputs:
        None

        Output:
        varList (list of (subproblem key, (cost, (LpVariable, value) dictionaries))) =
        initial solutions for the relaxed subproblems expressed as a cost and
        dictionary of non-zero values for variables
        """
        try:
          
            bvs = self.init_vars(self)
            if bvs is not None:
                if len(bvs) > 0:
                    return bvs
                else:
                    print("Empty variable list in generateInitVars, returning None")

        except Exception as ex:
            errorStr = "Error in generateInitVars\n%s" % ex
            raise DipError(errorStr)
  

class RelaxationCollection(object):
    """
    A simple defaultdict for holding relaxation problems
    """
    PROBLEM_CLASS = pulp.LpProblem


    def __init__(self, parent):
        self.parent = parent
        self.dict = {}

    def __getitem__(self, name):
        if name not in self.dict:
            self.dict[name] = self.PROBLEM_CLASS()
        return self.dict[name]

    def __setitem__(self, name, value):
        self.dict[name] = value

    def keys(self):
        return list(self.dict.keys())

    def values(self):
        return list(self.dict.values())
