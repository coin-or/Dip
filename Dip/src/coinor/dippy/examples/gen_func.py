import argparse

def addDippyArgs(parser):
    parser.add_argument('--CGLCuts', action='store_true',
                        help='enable CGL cuts')
    parser.add_argument('--algo', choices=['Cut', 'PriceCut', 'Price'],
                        help='which algorithm to use', default = 'Cut') 
    parser.add_argument('--solveMasterAsMip', action='store_true',
                        help='enable solving the master an an MIP')
    parser.add_argument('--nodeLogInterval', type=int, metavar = 'N',
                        help='node logging interval') 
    parser.add_argument('--nodeLimit', type=int, metavar = 'L',
                        help='node limit')
    parser.add_argument('--enforceBranchInSubproblem', action='store_true',
                        help='enforce branching constraints in subproblem')

def addDippyOpts(args):

    dippyOpts = {}

    if args.algo == 'PriceCut':
        dippyOpts['doPriceCut'] = '1'
        dippyOpts['CutCGL'] = '1'
    elif args.algo == 'Price':
        dippyOpts['doPriceCut'] = '1'
        dippyOpts['CutCGL'] = '0'
    else:
        dippyOpts['doCut'] = '1'

    if args.CGLCuts:
        dippyOpts['CutCGL'] = '1'
    else:
        dippyOpts['CutCGL'] = '0'

    if args.solveMasterAsMip:
        dippyOpts['SolveMasterAsMip'] = 1
    if args.enforceBranchInSubproblem:
        dippyOpts['enforceBranchInSubproblem'] = 1
        
    if args.nodeLimit or args.nodeLogInterval:
        dippyOpts['ALPS'] = {}
        if args.nodeLimit:
            dippyOpts['ALPS']['nodeLimit'] = args.nodeLimit
        if args.nodeLogInterval:
            dippyOpts['ALPS']['nodeLogInterval'] = args.nodeLogInterval

    return dippyOpts
