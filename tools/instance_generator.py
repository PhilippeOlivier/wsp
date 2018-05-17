# -*- coding: utf-8 -*-


################################################################################
# 
# This script generates 'num_instances' instances of size 'num_items'. Every
# instance shares the same item weights, but has a different cost matrix.
#
# Note: A necessary condition is that p+q+r <= 1.
#
################################################################################


################################## PARAMETERS ##################################

min_item_weight = 1
max_item_weight = 8
num_items = 50

p = 0.25 # Probability of c<0
q = 0.25 # Probability of c>0
r = 0.25 # Probability of conflict
conflict = 9999

min_cost = -5
max_cost = 5
num_instances = 1

################################################################################


import io
import os
import random
import sys


# Generates the cost matrices
def generate_matrices():
    matrices = []
    for _c in range(num_instances):
        matrix = [[0 for _ in range(num_items)] for _ in range(num_items)]
        for i in range(num_items):
            for j in range(0, i):
                x = random.random()
                if (x <= p):
                    cost = random.randint(min_cost, -1)
                    matrix[i][j] = cost
                    matrix[j][i] = cost
                elif (x <= p+q):
                    cost = random.randint(1, max_cost)
                    matrix[i][j] = cost
                    matrix[j][i] = cost
                elif (x <= p+q+r):
                    matrix[i][j] = conflict
                    matrix[j][i] = conflict
        matrices.append(matrix)
    return matrices


# Generates .wsp instance #n
def generate_instance(n):
    with io.open(os.path.join(dir_path+dir_name+"-"+str(n)+".wsp"), "w") as f:
        f.write(str(num_items)+"\n")

        for i in range(num_items):
            f.write(str(weights[i])+",")
        f.write("\n")

        for i in range(num_items):
            for j in range(num_items):
                f.write(str(matrices[n][i][j]) + ",")
            f.write("\n")
            
        f.close()

        
print("Working...")

weights = [random.randint(min_item_weight, max_item_weight)
           for _ in range(num_items)]
weights.sort(reverse=True)
matrices = generate_matrices()

num_bins = int(round(float(sum(weights))/10, 0))
print("Instance needs " + str(num_bins) + " bins.")

# Create the directory
dir_name = str(sys.argv[1])
dir_path = "./../data/"+dir_name+"/"
os.mkdir(dir_path)

# Create the instances and the instances.txt file
with io.open(os.path.join(dir_path+"instances.txt"), "w") as f:
    for i in range(0, num_instances):
        generate_instance(i)
        f.write(dir_name+"-"+str(i)+".wsp\n")
    f.close()

# Create the parameters.txt file
with io.open(os.path.join(dir_path+"parameters.txt"), "w") as f:
    f.write("num_items="+str(num_items)+"\n")
    f.write("num_bins="+str(num_bins)+"\n")
    f.write("total_weight="+str(sum(weights))+"\n")
    mean_load = float(sum(weights))/float(num_bins)
    f.write("mean_load="+str(mean_load)+"\n")
    f.write("min_item_weight="+str(min_item_weight)+"\n")
    f.write("max_item_weight="+str(max_item_weight)+"\n")
    f.write("min_cost="+str(min_cost)+"\n")
    f.write("max_cost="+str(max_cost)+"\n")
    f.write("p="+str(p)+"\n")
    f.write("q="+str(q)+"\n")
    f.write("r="+str(r)+"\n")
    
print("Finished generating instances.")
