#!/usr/bin/python


import array
import os
import commands
import time
import math
import sys
import string
import fileinput
import re
import glob

# ---------------------------------------------------------------------------
# IMPORTANT PARAMETERS
# ---------------------------------------------------------------------------
#todo - try using os.path.join
HOME            = os.environ['HOME'] + "/"
COIN_DIR        = HOME + "COIN/coin-Decomp/"
COIN_DATA_DIR   = COIN_DIR + "Decomp/data/"
RUN_DIR         = HOME + "running/decomp/"

INST_SET        = "tsplib"
DATA_DIR        = COIN_DATA_DIR + "TSP/TSPLIB/"
#BEST_BOUND_FILE = DATA_DIR + "ap3.opt"
INST_LIST       = DATA_DIR + "tsplib.list"

OUTPUT_DIR1     = RUN_DIR + "tsp/tspc/"
OUTPUT_DIR2     = RUN_DIR + "tsp/tspp/"

INFTY=10e10
ABS_GAPTOL=0.1
EPS_TIME = 0.01
#TIME_LIMIT = 1800
delim = "  "
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------



def find_float(arr0,st0,fl0):
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			fl0 = float(re.search('-*[0-9]+\.[0-9]+',line).group())
			return 1,fl0
	return -1,fl0

def find_float_1(arr0,st0,fl0):
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			find = re.search('-*[0-9]+\.[0-9]+',line)
			if (find>=0):
				l = re.findall('-*[0-9]+\.[0-9]+',line)
				if (len(l)<3):
					return -1,fl0
				else:
					fl0 = float(l[2])
					return 1,fl0
			else:
				return -1,fl0
	return -1,fl0

def find_float_e(arr0,st0,fl0):
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			fl0 = float(re.search('-*[0-9]+\.[0-9]*e*\+*[0-9]*',line).group())
			return 1,fl0
	return -1,fl0

def find_int(arr0,st0,in0):
	in0=0
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			in0 = int(re.search('-*[0-9]+$',line).group())
			return 1,in0
	return -1,in0

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

for i in range(0,len(sys.argv)):
	if (sys.argv[i]=='-d'):
		if (i==len(sys.argv)):
			print "Missing argument to '-d'"
			print "usage: python report_decomp.py [-d <path to dir>]"
			sys.exit(0)
		else:
			OUTPUT_DIR=sys.argv[i+1]
			break

a = []
fl0 = 0.0
in0 = 0

flist=open(INST_LIST,'r')
a=flist.read().split()
flist.close()
a.sort()

print "### Instance set:", INST_SET



#print "##%9s"%"Instance",delim,"%12s"%"CPX Obj",delim,"%7s"%"miptime",delim,"%12s"%"tkegraph",delim,"%8s"%"Time",delim,"%6s"%"binit","%8s"%"cost-gap",delim,"%8s"%"time-LP",delim,"%8s"%"time-GS",delim,"%9s"%"time-cost",delim,"%6s"%"edges",delim,"%6s"%"constr"

#TODO - get cplex times... 

#check for errors first
print "### Reading from directory:", OUTPUT_DIR1
for instance in a:
        filewcard = OUTPUT_DIR1 + instance + ".condor.err.*"
        filename = glob.glob(filewcard)
        statinfo = os.stat(filename[0])
        if (statinfo.st_size):
                print "ERROR %10s"%instance

print "### Reading from directory:", OUTPUT_DIR1
for instance in a:
        filewcard = OUTPUT_DIR1 + instance + ".condor.err.*"
        filename = glob.glob(filewcard)[0]
        statinfo = os.stat(filename)
        if (statinfo.st_size):
                print "ERROR %10s" % instance
                fil        = open(filename, 'r')
                whole_file = fil.read().split('\n')
                fil.close()
                print whole_file



# print "### Reading from directory:", OUTPUT_DIR2
# for instance in a:
#         filewcard = OUTPUT_DIR2 + instance + ".condor.err.*"
#         filename = glob.glob(filewcard)[0]
#         statinfo = os.stat(filename)
#         if (statinfo.st_size):
#                 print "ERROR %10s" % instance
#                 fil        = open(filename, 'r')
#                 whole_file = fil.read().split('\n')
#                 fil.close()
#                 print whole_file













# for instance in a:
# 	print "%11s"%instance,

# 	fil=open(BEST_BOUND_FILE,'r')
# 	whole_file=fil.read().split('\n')
# 	fil.close()

# 	best_ub = INFTY
# 	find,best_ub=find_float(whole_file,instance,best_ub)
# 	print delim, "%12.1f"%best_ub,

# #TODO: check if file err file is size 0 


# #make this a function
# 	filename=OUTPUT_DIR1+"/"+INST_SET+"/"+instance+".condor.out"
# 	fil=open(filename,'r')
# 	whole_file=fil.read().split('\n')
# 	fil.close()
	
# 	bestObjValue = INFTY
# 	find,bestObjValue=find_float(whole_file,'flow cost =',bestObjValue)
# 	if (find<0 or bestObjValue >= INFTY):
# 		print delim, "%12s"%"NF",
# 	else:
# 		print delim, "%12.1f"%bestObjValue,

# 	solTime = INFTY
# 	find,solTime=find_float(whole_file,'time taken =',solTime)
# 	if (find<0 or solTime >= INFTY):
# 		print delim, "%8s"%"NF",
# 	else:
# 		if (solTime<EPS_TIME):
# 			solTime = EPS_TIME
# 		print delim, "%8.2f"%solTime,


# 	n = INFTY
# 	find,n = find_int(whole_file,'binary-search iterations =',n)
# 	if (find<0):
# 		print delim, "%6s"%"NF",
# 	else:
# 		print delim, "%6d"%n,

                

# 	"""
# 	flowCongestion = INFTY
# 	find,flowCongestion =find_float(whole_file,'flow congestion =',flowCongestion)
# 	if (find<0 or flowCongestion >= INFTY):
# 		print delim, "%8s"%"NF",
# 	else:
# 		print delim, "%8.2f"%flowCongestion,
# 	"""

# 	if (bestObjValue<INFTY and best_ub<INFTY):
# 		costFactor = abs(bestObjValue-best_ub)/bestObjValue*100,
# 		print delim, "%8.2f"%costFactor,
# 	else:
# 		print delim, "%8s"%"NF",

# 	timeInLP = INFTY
# 	find,timeInLP =find_float(whole_file,'time in LP solves =',timeInLP)
# 	if (find<0 or timeInLP >= INFTY):
# 		print delim, "%8s"%"NF",
# 	else:
# 		print delim, "%8.2f"%timeInLP,

# 	timeInGS = INFTY
# 	find,timeInGS =find_float(whole_file,'time in golden-section =',timeInGS)
# 	if (find<0 or timeInGS >= INFTY):
# 		print delim, "%8s"%"NF",
# 	else:
# 		print delim, "%8.2f"%timeInGS,

# 	timeInUF = INFTY
# 	find,timeInUF =find_float(whole_file,'time in update-cost =',timeInUF)
# 	if (find<0 or timeInUF >= INFTY):
# 		print delim, "%9s"%"NF",
# 	else:
# 		print delim, "%9.2f"%timeInUF,

# 	n = INFTY
# 	find,n = find_int(whole_file,'Number of edges              =',n)
# 	if (find<0):
# 		print delim, "%6s"%"NF",
# 	else:
# 		print delim, "%6d"%n,
	
# 	n = INFTY
# 	find,n = find_int(whole_file,'Number of mutual constraints =',n)
# 	if (find<0):
# 		print delim, "%6s"%"NF",
# 	else:
# 		print delim, "%6d"%n,

# 	print ''

