################################################################################
# 
# This script generates a figure showing the Pareto set for a batch of results.
# This script takes as its only argument the name of a results directory found
# in /wsp/results.
#
# Example:
# python3 visualize.py my_dataset_directory_name
#
# The figure is saved as figure.pdf in the directory passed as an argument,
# along with file average_times.txt, showing the average time for each model.
#
################################################################################


from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import operator
import os
import re
import sys


### FUNCTIONS ##################################################################


# Gets the results of a model for each slice of the Pareto set.
def get_results(model):
    results = []
    for i in range(1, d_max+1):
        results.append(get_slice(model, i-1, i))
    return results


# Gets the results for a slice, averaged for all instances.
# [d_min, d_max, actual deviation, costs, time]
def get_slice(model, d_min, d_max):
    model_dir = model
    if (model == "lb"):
        model_dir = "ipb"
    num_instances = 0
    num_feasible = 0
    averaged_slice = [d_min, d_max, 0, 0, 0]
    for instance in os.listdir(results_dir_path+"/"+model_dir):
        if (re.search(dataset+".+dev"+str(d_max)+".csv", instance)):
            instance_slice = []
            if (model == "cp"):
                instance_slice = get_slice_cp(instance)
            elif (model == "ipa"):
                instance_slice = get_slice_ipa(instance)
            elif (model == "ipb"):
                instance_slice = get_slice_ipb(instance)
            elif (model == "lb"):
                instance_slice = get_slice_lb(instance)
            num_instances += 1
            averaged_slice[4] += instance_slice[2]
            if (instance_slice[0] != None):
                averaged_slice[2] += instance_slice[0]
                averaged_slice[3] += instance_slice[1]
                num_feasible += 1
                
    if (num_feasible == 0):
        averaged_slice[2] = None
        averaged_slice[3] = None
    else:
        averaged_slice[2] /= float(num_feasible)
        averaged_slice[3] /= float(num_feasible)
    averaged_slice[4] /= num_instances
        
    return averaged_slice


# Gets the results for a single instance/slice of model CP.
# [deviation, costs, time]
def get_slice_cp(instance):
    with open(results_dir_path+"/cp/"+instance) as f:
        lines = [line.rstrip('\n') for line in f]

    if (len(lines) == 2 or len(lines) == 3):
        return [None, None, float(lines[-1])]

    try:
        instance_slice = [float(i) for i in lines[-3].split(',')]
        instance_slice[0], instance_slice[1] = instance_slice[1], instance_slice[0]
        instance_slice.append(float(lines[-1]))
    except ValueError:
        print("ERROR: unbound constrained variable in "+instance)
        exit(-1)
    
    return instance_slice


# Gets the results for a single instance/slice of model IPA.
# [deviation, costs, time]
def get_slice_ipa(instance):
    with open(results_dir_path+"/ipa/"+instance) as f:
        lines = [line.rstrip('\n') for line in f]
        
    if (re.search("No solution exists.", lines[0])):
        return [None, None, float(lines[1])]
    
    instance_slice = [float(i) for i in lines[3].split(',')] # TODO: Replace by lines[2] once IPA is fixed
    instance_slice[0], instance_slice[1] = instance_slice[1], instance_slice[0]
    
    return instance_slice


# Gets the results for a single instance/slice of model IPB.
# [deviation, costs, time]
def get_slice_ipb(instance):
    with open(results_dir_path+"/ipb/"+instance) as f:
        lines = [line.rstrip('\n') for line in f]
        
    if (re.search("No solution exists.", lines[0])):
        return [None, None, float(lines[1].split(',')[5])]
    
    instance_slice = [float(i) for i in lines[2].split(',')][2:]
    instance_slice[0], instance_slice[1] = instance_slice[1], instance_slice[0]
    
    return instance_slice


# Gets the results for a single instance/slice of model IPB (lower bound).
# [deviation, costs, time]
def get_slice_lb(instance):
    with open(results_dir_path+"/ipb/"+instance) as f:
        lines = [line.rstrip('\n') for line in f]
        
    if (re.search("No solution exists.", lines[0])):
        return [None, None, 0]
    
    instance_slice = [float(i) for i in lines[2].split(',')][:2]
    instance_slice[0], instance_slice[1] = instance_slice[1], instance_slice[0]
    instance_slice[1] = int(instance_slice[1]) # Round LB to integer
    instance_slice.append(0)
    
    return instance_slice

    
### MAIN #######################################################################


# Sanity check
results_dir_name = sys.argv[1]
results_dir_path = "../results/"+results_dir_name
if (not os.path.exists("../results/"+results_dir_name)):
    print("ERROR: Directory /results/"+results_dir_name+" not found.")
    exit(-1)

# Get parameters
with open(results_dir_path+"/parameters.txt") as f:
    params = [param.rstrip('\n') for param in f]
dataset = re.search("dataset=(.+)", params[0]).group(1)
norm = int(re.search("norm=(.+)", params[1]).group(1))
timelimit = int(re.search("timelimit=(.+)", params[2]).group(1))
d_max = int(re.search("d_max=(.+)", params[3]).group(1))

# Get results for each model
models = {"cp": None,
          "ipa": None,
          "ipb": None,
          "lb": None}

for model in models:
    if model in [sd for sd in os.listdir(results_dir_path) \
             if os.path.isdir(os.path.join(results_dir_path, sd))]:
        models[model] = get_results(model)
        if (model == "ipb"):
            models["lb"] = get_results("lb")
            
# Make the results suitable for a Pareto set
for model in models:
    if (models[model] is not None):
        # Ensure the lists are sorted
        sorted(models[model], key=operator.itemgetter(0))

        # Ensure solutions get increasingly better
        best = 999999
        for i in models[model]:
            # No feasible solution found yet
            if (i[3] is None) and \
               (best == 999999):
                pass
            # Current solution is at least as good as the previous
            elif (i[3] is not None):
                i[3] = min(best, i[3])
                best = i[3]
            # Infeasible interval
            elif (i[3] is None):
                i[3] = best

# Gather plotting information
plots = {"cp": None,
         "ipa": None,
         "ipb": None,
         "lb": None}

for model in models:
    if (models[model] is not None):
        plots[model] = dict()
        plots[model]["x"] = [row[1] for row in models[model]]
        plots[model]["y"] = [row[3] for row in models[model]]
        if (model == "cp"):
            plots[model]["color"] = "red"
            plots[model]["label"] = "CP"
            plots[model]["marker"] = "o"
        elif (model == "ipa"):
            plots[model]["color"] = "green"
            plots[model]["label"] = "IP(A)"
            plots[model]["marker"] = "+"
        elif (model == "ipb"):
            plots[model]["color"] = "blue"
            plots[model]["label"] = "IP(B)"
            plots[model]["marker"] = "x"
        elif (model == "lb"):
            plots[model]["color"] = "black"
            plots[model]["label"] = "L.B."
            plots[model]["linestyle"] = "-"

# Plot the results
fig = plt.figure(figsize = [9.6,4.8])
plt.ylabel('Objective')
plt.xlabel('Deviation')
plt.xticks([i for i in range(0, d_max+1, 1)])
plt.xlim([0, d_max+1])

for model in models:
    if (models[model] is not None and model != "lb"):
        plots[model]["id"] = plt.scatter(plots[model]["x"],
                                         plots[model]["y"],
                                         color=plots[model]["color"],
                                         label=plots[model]["label"],
                                         marker=plots[model]["marker"])

    if (models[model] is not None and model == "lb"):
        plots[model]["id"], = plt.plot(plots[model]["x"],
                                       plots[model]["y"],
                                       color=plots[model]["color"],
                                       label=plots[model]["label"],
                                       linestyle=plots[model]["linestyle"])

# Draw legend
legend_ids = []
legend_labels = []
for model in models:
    if (models[model] is not None):
        legend_ids.append(plots[model]["id"])
        legend_labels.append(plots[model]["label"])
plt.legend(legend_ids, legend_labels, bbox_to_anchor=(0.05, 0.33), loc=2, borderaxespad=0.)

# Save figure
pp = PdfPages(results_dir_path+"/figure.pdf")
pp.savefig(fig, bbox_inches='tight')
pp.close()

# Generate the times.txt file
open(results_dir_path+"/average_times.txt", 'w').close()
for model in models:
    if (models[model] is not None and model != "lb"):
        with open(results_dir_path+"/average_times.txt", 'a') as f:
            f.write(model+
                    "="+
                    str(sum([row[4] for row in models[model]])/len(models[model]))+
                    '\n')
