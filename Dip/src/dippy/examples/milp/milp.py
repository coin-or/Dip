'''
Created on Dec 29, 2013

@author: ted
'''
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div

#!/usr/bin/env python

import sys
import random
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize

try:
    import path
except ImportError:
    pass
        
try:
    import src.dippy as dippy
except ImportError:
    import coinor.dippy as dippy

def GenerateRandomBlock(VARIABLES, CONSTRAINTS, density = 0.2,
                        maxObjCoeff = 10, maxConsCoeff = 10, 
                        tightness = 2, rand_seed = 2):
    random.seed(rand_seed)
    OBJ = dict((i, random.randint(1, maxObjCoeff)) for i in VARIABLES)
    MAT = dict(((i, j), random.randint(1, maxConsCoeff) 
                        if random.random() <= density else 0)
                for j in CONSTRAINTS for i in VARIABLES)
    RHS = dict((i, random.randint(int(numVars*density*maxConsCoeff/tightness),
                                  int(numVars*density*maxConsCoeff/1.5)))
                                  for i in CONSTRAINTS)
    return OBJ, MAT, RHS

#display_mode = 'xdot'
#layout = 'dot'

tol = pow(pow(2, -24), old_div(2.0, 3.0))

prob = dippy.DipProblem("MILP")

numBlocks = 10
numBlockVars = [10]
numBlockCons = [5]
numLinkingCons = 10
numVars = sum(numBlockVars)
numCons = sum(numBlockCons) + numLinkingCons

VARIABLES = dict(((i, j), 0) for i in range(numBlocks) 
                 for j in range(numBlockVars[i]))  

CONSTRAINTS = []
for k in range(numBlocks):
    CONSTRAINTS.append(["C"+str(k)+"_"+str(j) for j in range(numCons)])
CONSTRAINTS.append(["C"+str(numBlocks)+"_"+str(j) for j in range(numCons)])

#Generate random MILP
var = LpVariable.dicts("x", VARIABLES, 0, 1, LpBinary)
numCons = len(CONSTRAINTS)
numVars = len(VARIABLES)

OBJ, MAT, RHS = GenerateRandomBlock(VARIABLES, CONSTRAINTS[numBlocks])

prob += -lpSum([OBJ[i]*var[i] for i in var]), "Objective"

#Linking constraints
for j in CONSTRAINTS[numBlocks]:
    prob += lpSum([MAT[i, j]*var[i] for i in var]) <= RHS[j], j

#Blocks    
for k in range(numBlocks):
    OBJ, MAT, RHS = GenerateRandomBlock([(k, i) for i in range(numBlockVars[k])], 
                                        CONSTRAINTS[k])
    for j in CONSTRAINTS[k]:
        prob.relaxation[k] += (lpSum([MAT[(k, i), j]*var[k, i] for i in range(numBlockVars[k])])
                               <=RHS[j], j)

dippy.Solve(prob,{
    'TolZero': '%s' % tol,
    'doDirect': '1',
    'CutCGL': '1',
    'SolveMasterAsIp': '0',
    'generateInitVars': '1',
    'LogDebugLevel': '3',
    'LogLevel': '4',
    'LogDumpModel': 5,
    'ALPS' :
    {'msgLevel' : 3}}
)

if prob.display_mode != 'off':
    numNodes = len(prob.Tree.get_node_list())
    if prob.Tree.attr['display'] == 'svg':
        prob.Tree.write_as_svg(filename = "facility_node%d" % (numNodes + 1), 
                               prevfile = "facility_node%d" % numNodes)
    prob.Tree.display()

