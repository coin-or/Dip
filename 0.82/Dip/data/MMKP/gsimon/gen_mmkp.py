#!/usr/bin/python

import commands
import os
from math import *
from   optparse import OptionParser
from xml.dom.minidom import parse
from random import *
import sys

def generate_uniform(n, min_v, max_v, seed):
    weights = []
    r = Random()
    for i in range(n):
	weights.append(r.uniform(min_v, max_v))
    return weights

def generate_linear(n, min_v, max_v):
    weights = []
    delta = max_v - min_v
    for i in range(n):
	weights.append(min_v + (delta*(i)) / float(n))
    return weights

def generate_fixed(value_sting):
    weights = []
    values = value_string.split()
    for value in values:
	weights.append(float(value))
    return weights


def generate_inv_strongly_correlated(profits, delta):
    weights = []
#    print profits
    max_profit = max(profits)
    for i in range(len(profits)):
	weights.append(max_profit - profits[i]/delta )
    return weights

def generate_inv_weakly_correlated(profits, delta):
    weights = []
    max_profit = max(profits)
    r = Random()

    for i in range(len(profits)):
	weights.append( r.uniform(0, max_profit- profits[i]/delta))
    return weights


def generate_strongly_correlated(profits, delta):
    weights = []
#    print profits
    max_profit = max(profits)
    for i in range(len(profits)):
	weights.append(profits[i] + max_profit/delta )
    return weights

def generate_weakly_correlated(profits, delta):
    weights = []
    max_profit = max(profits)
    r = Random()

    for i in range(len(profits)):
	weights.append( r.uniform(max(0, profits[i] - max_profit/delta), profits[i] + max_profit/delta))
    return weights

#####################################################
## Parse the command line options.
parser = OptionParser()
parser.add_option("--def_file", dest="def_file", type="string", default="",
		help="instance definition file" )


(options, args) = parser.parse_args()
out_file = sys.stdout
# Construct xml parser and read the definition file. 
dom = parse(options.def_file)
mmkps = dom.getElementsByTagName("mmkp")

for mmkp in mmkps:
#{
    assert mmkp.attributes.has_key('name')
    assert mmkp.attributes.has_key('classes')
    assert mmkp.attributes.has_key('dimensions')
    assert mmkp.attributes.has_key('instances')

    serie_name = mmkp.attributes['name'].value    
    num_classes = int(mmkp.attributes['classes'].value)
    num_dimensions= int(mmkp.attributes['dimensions'].value)
    num_instances = int(mmkp.attributes['instances'].value) 
    num_constraint_levels = int(mmkp.attributes['constraint_levels'].value) 

    if mmkp.attributes.has_key('p_value'):
	p_value = mmkp.attributes['p_value'].value
    else:
	p_value = "float" 
    
    if mmkp.attributes.has_key('w_value'):
	w_value = mmkp.attributes['w_value'].value
    else:
	w_value = "float" 

    if mmkp.attributes.has_key('c_value'):
	c_value = mmkp.attributes['c_value'].value
    else:
	c_value = "float" 

    # inst_type could be "generic" or "wsn"
    if mmkp.attributes.has_key('inst_type'):
	inst_type = mmkp.attributes['inst_type'].value
    else:
	inst_type = "generic" 
    
    if inst_type != "generic":
	print "Cannot process definition files for non-generic instances."
	exit()

    for inst_id in range(num_instances):
    #{
	
	instance_txt = [] 
	instance_txt.append("# MMKP instances generated from definition file: " + options.def_file + '\n')
        instance_txt.append("# serie_name\t" + serie_name + "\n")
	instance_txt.append("# num_instances\t" + repr(num_instances) +"\n")

	instance_txt.append('instance_id\t' +repr(inst_id+1) + '\n')
	instance_txt.append('num_class\t' +repr(num_classes)+ '\n')
	instance_txt.append('num_dimension\t' +repr(num_dimensions)+ '\n') 
	
	# Process each class in an mmkp
	classes = mmkp.getElementsByTagName('class')
	sum_min_weight = []
	sum_max_weight = []
	for d in range(num_dimensions):
	    sum_min_weight.append(0)
	    sum_max_weight.append(0)
	
	c_id = 0
	for one_class in classes: 
	#{
	    num_item = int(one_class.attributes['items'].value) 
	    
	    if one_class.attributes.has_key('repeat'):
		c_repeat = int(one_class.attributes['repeat'].value)  
	    else: 
		c_repeat = 1
	   
	    for repeat_id in range(c_repeat):
	    #{
		c_id = c_id + 1
		instance_txt.append('# class_id\t' + repr(c_id) + '\n')
		instance_txt.append('num_item\t' + repr(num_item)+ '\n')
		# Process profits
		profit_function = one_class.getElementsByTagName('profit')[0]
		method = profit_function.attributes['method'].value
		if method == 'uniform':
		#{
		    parameters = profit_function.getElementsByTagName('parameters')[0]
		    min_v = float(parameters.attributes['min'].value)
		    max_v = float(parameters.attributes['max'].value)
		    profit_values = generate_uniform(num_item, min_v, max_v, c_id)
		#}
		elif method == 'linear':
		#{
		    parameters = profit_function.getElementsByTagName('parameters')[0]
		    min_v = float(parameters.attributes['min'].value)
		    max_v = float(parameters.attributes['max'].value)
		    profit_values = generate_linear(num_item, min_v, max_v)
		#}
		elif method == 'fixed':
		#{
		    parameters = profit_function.getElementsByTagName('parameters')[0]
		    value_string= (parameters.attributes['values'].value)
		    profit_values = generate_fixed(value_string)
		#}

		# Process each type of weight
		weight_functions = one_class.getElementsByTagName('weight')	
		weight_value_matrix=[]
		
		w_id = 0
		for weight_function in weight_functions:
		#{
		    method = weight_function.attributes['method'].value
		    parameters = weight_function.getElementsByTagName('parameters')[0]
		    if weight_function.attributes.has_key('repeat'):
			w_repeat = int(weight_function.attributes['repeat'].value)  
		    else: 
			w_repeat = 1
		   
		    for weight_repeat_id in range(w_repeat):
		    #{

			w_id = w_id + 1
			if method == 'uniform':
			    min_v = float(parameters.attributes['min'].value)
			    max_v = float(parameters.attributes['max'].value)
			    weight_values = generate_uniform(num_item, min_v, max_v, c_id+w_id)

			elif method == 'strongly_correlated':
			    delta = float(parameters.attributes['delta'].value)
			    weight_values = generate_strongly_correlated(profit_values, delta)

			elif method == 'weakly_correlated':
			    delta = float(parameters.attributes['delta'].value)
			    weight_values = generate_weakly_correlated(profit_values, delta)
		    	elif method == 'inv_strongly_correlated':
			    delta = float(parameters.attributes['delta'].value)
			    weight_values = generate_inv_strongly_correlated(profit_values, delta)
			elif method == 'inv_weakly_correlated':
			    delta = float(parameters.attributes['delta'].value)
			    weight_values = generate_inv_weakly_correlated(profit_values, delta)

			weight_value_matrix.append(weight_values)
			sum_min_weight[w_id-1] = sum_min_weight[w_id-1] + min(weight_values)
			sum_max_weight[w_id-1] = sum_max_weight[w_id-1] + max(weight_values)
#		print sum_min_weight, min(weight_values)
#		print sum_max_weight, max(weight_values)
#		print weight_values
		    #}
		#}
		
		for i in range(len(profit_values)):
		    if p_value == "int":
			instance_txt.append("%d\t" % ceil(profit_values[i]))
		    else:
			instance_txt.append("%.2f\t" % profit_values[i])
		    for j in range(num_dimensions):
			if w_value == "int":    
			    instance_txt.append("%d\t" % ceil(weight_value_matrix[j][i]))
			else: 
			    instance_txt.append("%.2f\t" % weight_value_matrix[j][i])
		    instance_txt.append('\n')

	    #}
	#}
	
	constraint_txt = []
	for level in range(num_constraint_levels):
	#{
	    out_file_name = serie_name + '_' + repr(level) +'_' + repr(inst_id) + '.data'
	    out_file = open(out_file_name, "w")
	    for line in instance_txt:
		out_file.write(line)
	    
	    out_file.write("c:\t")
	    for d in range(num_dimensions):
		constraint = sum_min_weight[d] + (float(level+1)/float(num_constraint_levels+1)) * (sum_max_weight[d] - sum_min_weight[d])
		if c_value == "int":
		    out_file.write("%d\t" % ceil(constraint))
		else: 
		    out_file.write("%.2f\t" % constraint)
	    out_file.write('\n')
	    out_file.close()
	#}	
    #}
    
#}
