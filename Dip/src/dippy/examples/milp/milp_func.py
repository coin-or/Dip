'''
Created on Dec 29, 2013

@author: ted
'''
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import argparse
import random
from pulp import LpVariable, LpBinary, lpSum, value, LpProblem, LpMaximize

try:
    from src.dippy import DipProblem
    from src.dippy.examples.gen_func import *
except ImportError:
    from coinor.dippy import DipProblem
    from coinor.dippy.examples.gen_func import *

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Generate and solve a random block-structured MILP.')
    parser.add_argument('--randomSeed', '-r', metavar = 'R', type=int,
                        help='a random seed', default = 2)
    parser.add_argument('--numBlocks', '-b', metavar = 'NB', type=int,
                        help='number of blocks', default = 2) 
    parser.add_argument('--numVarsPerBlock', '-v', metavar = 'NV', type=int,
                        help='number of variables per block', default = 10)
    parser.add_argument('--numConsPerBlock', '-c', metavar = 'NC', type=int,
                        help='number of constraints per block', default = 5)
    parser.add_argument('--numLinkingCons', '-l', metavar = 'NL', type=int,
                        help='number of linking constraints', default = 10)
    parser.add_argument('--tightness', '-t', metavar = 'T', type=int,
                        help='how tight the constraints should be (higher is tighter)',
                        default = 2)
    parser.add_argument('--density', '-d', metavar = 'D', type=int,
                        help='density of the constraint matrix', default = 0.2)
    parser.add_argument('--maxConsCoeff', '-m', metavar = 'M', type=int,
                        help='maximum size of a constraint coefficient',
                        default = 10)

    addDippyArgs(parser)
    
    args = parser.parse_args()

    return(args)
    
def GenerateRandomBlock(VARIABLES, CONSTRAINTS, density = 0.2,
                        maxObjCoeff = 10, maxConsCoeff = 10, 
                        tightness = 2, rand_seed = 2):
    random.seed(rand_seed)
    numCons = len(CONSTRAINTS)
    numVars = len(VARIABLES)
    OBJ = dict((i, random.randint(1, maxObjCoeff)) for i in VARIABLES)
    MAT = dict(((i, j), random.randint(1, maxConsCoeff) 
                        if random.random() <= density else 0)
                for j in CONSTRAINTS for i in VARIABLES)
    RHS = dict((i, random.randint(int(numVars*density*maxConsCoeff/tightness),
                                  int(numVars*density*maxConsCoeff/1.5)))
                                  for i in CONSTRAINTS)
    return OBJ, MAT, RHS

def formulate(args):

    prob = DipProblem("MILP")

    numBlocks = args.numBlocks
    numBlockVars = [args.numVarsPerBlock for i in range(numBlocks)]
    numBlockCons = [args.numConsPerBlock for j in range(numBlocks)]
    numVars = sum(numBlockVars)
    numCons = sum(numBlockCons) + args.numLinkingCons

    VARIABLES = dict(((i, j), 0) for i in range(numBlocks) 
                 for j in range(numBlockVars[i]))  

    CONSTRAINTS = []
    for k in range(numBlocks):
        CONSTRAINTS.append(["C"+str(k)+"_"+str(j) for j in range(numCons)])
    CONSTRAINTS.append(["C"+str(numBlocks)+"_"+str(j) for j in range(numCons)])

    #Generate random MILP
    var = LpVariable.dicts("x", VARIABLES, 0, 1, LpBinary)

    OBJ, MAT, RHS = GenerateRandomBlock(VARIABLES, CONSTRAINTS[numBlocks],
                                        rand_seed = args.randomSeed, tightness = args.tightness,
                                        density = args.density)

    prob += -lpSum([OBJ[i]*var[i] for i in var]), "Objective"

    #Linking constraints
    for j in CONSTRAINTS[numBlocks]:
        prob += lpSum([MAT[i, j]*var[i] for i in var]) <= RHS[j], j

    #Blocks    
    for k in range(numBlocks):
        OBJ, MAT, RHS = GenerateRandomBlock([(k, i) for i in range(numBlockVars[k])], 
                                        CONSTRAINTS[k])
        for j in CONSTRAINTS[k]:
            prob.relaxation[k] += (lpSum([MAT[(k, i), j]*var[k, i] for i in
                                          range(numBlockVars[k])]) <=RHS[j], j)
    return prob

