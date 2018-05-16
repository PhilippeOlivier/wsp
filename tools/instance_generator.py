# -*- coding: utf-8 -*-


################################################################################
# 
# This script generates instances of size 'g' (the number of groups) for
# every p-q combination in the collection of relation matrices 'r'.
#
# Three distinct (small, random, and large) group layouts are generated. Each
# one of the layouts will be matched with every relation matrix. A total of 3*r
# instances will thus be generated.
#
# g : Number of groups.
# p : Probability that a pair of groups are "definitely apart".
# q : Probability that a pair of groups are either "rather together",
#     "rather apart", or "indifferent" (chosen randomly).
#
# Note: A necessary condition is that p+q <= 1.
#p, q, r: odds that c<0, c>0, conflict
#

# file layout:
# num_items
# total_weight
# num_bins
# w1,w2,w3,wn
# matrix

# instances are named .wsp. Also, no pickle

################################################################################


################################## PARAMETERS ##################################

min_item_weight = 1
max_item_weight = 8
num_items = 25

p = 0.25 # Probability of c<0
q = 0.25 # Probability of c>0
r = 0.25 # Probability of conflict
conflict = 9999

min_cost = -5
max_cost = 5
num_matrices = 5
convert_to_wsp = True

# add avg load per table when fixed number of tables

################################################################################


import datetime
import io
import math
import os.path
import pickle
import random
import sys


# Generates the cost matrices.
def generate_matrices():
    matrices = []
    for _c in range(num_matrices):
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
                else:
                    pass

        if (convert_to_wsp):
            # 0 = indifferent
            # 1 = rather apart
            # 2 = definitely apart
            # 3 = rather together
            for i in range(num_items):
                for j in range(0, i):
                    if (matrix[i][j] == conflict):
                        matrix[i][j] = 2
                        matrix[j][i] = 2
                    elif (matrix[i][j] < 0):
                        matrix[i][j] = 3
                        matrix[j][i] = 3
                    elif (matrix[i][j] > 0):
                        matrix[i][j] = 1
                        matrix[j][i] = 1

        matrices.append(matrix)

    return matrices


# Generates a .wpn instance and its corresponding pickle file.
def generate_instance(n):
    names = [["I"+str(j+1)+"-"+str(i+1)
              for i in range(weights[j])] for j in range(num_items)]

    # Create .wpn file
    with io.open(os.path.join(dir_path+dir_name+"-"+str(n)+".wpn"),
                 "w",
                 newline="\r\n") as f:
        
        f.write("THIS_FILE_IS_FOR_EXCLUSIVE_USE_OF_THE_WEDDING_SEAT" +
                "_PLANNER_www.weddingSeatPlanner.com_DO_NOT_EDIT_" +
                "Pat.no_X4993X33\n")
        
        f.write("0,"+str(num_items)+"\n")

        for i in range(num_items):
            for j in range(len(names[i])):
                f.write(str(names[i][j])+",")
            f.write("\n")
        
        for i in range(num_items):
            f.write(str(names[i][0]))
            if len(names[i]) > 1:
                f.write("+"+str(len(names[i])-1))
            f.write("\n")

        for i in range(num_items):
            for j in range(num_items):
                f.write(str(matrices[n][i][j]) + ",")
            f.write("\n")
            
        f.close()
    
    # pickle file
    pickle.dump([weights, num_bins, matrices[n]],
                open(os.path.join(dir_path+dir_name+"-"+str(n)+".pickle"),
                "wb"))

        
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
    for i in range(0, num_matrices):
        generate_instance(i)
        f.write(dir_name+"-"+str(i)+".wpn\n")
    f.close()

# Create the parameters.txt file
with io.open(os.path.join(dir_path+"parameters.txt"), "w") as f:
    f.write("num_items="+str(num_items)+"\n")
    f.write("num_bins="+str(num_bins)+"\n")
    f.write("total_weight="+str(sum(weights))+"\n")

    mean_load = float(sum(weights))/float(num_bins)
    f.write("mean_load="+str(mean_load)+"\n")
    minimum_feasible_deviation = num_bins * \
                                 (math.ceil(mean_load)-mean_load) * \
                                 (mean_load-math.floor(mean_load))
    f.write("minimum_feasible_deviation="+str(minimum_feasible_deviation)+"\n")

    f.write("min_item_weight="+str(min_item_weight)+"\n")
    f.write("max_item_weight="+str(max_item_weight)+"\n")
    f.write("min_cost="+str(min_cost)+"\n")
    f.write("max_cost="+str(max_cost)+"\n")
    f.write("wsp="+str(convert_to_wsp)+"\n")

    f.write("p="+str(p)+"\n")
    f.write("q="+str(q)+"\n")
    f.write("r="+str(r)+"\n")
    
print("Finished generating instances.")
