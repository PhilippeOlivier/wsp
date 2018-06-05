#!/bin/bash


################################################################################
# 
# This script automates testing for a model on a dataset. The four parameters
# below must be specified. This script takes as its only argument the name of a
# dataset directory found in /wsp/data.
#
# Example:
# ./test.sh my_dataset
#
# The results of these tests are saved in /wsp/results/my_dataset.
#
################################################################################


### PARAMETERS #################################################################

model=ipb     # {ipa, ipb}
norm=1        # {0, 1, 2, 3}
timelimit=5
d_max=3

################################################################################


# Sanity checks
if [ $model != "ipa" ] && [ $model != "ipb" ]; then
    echo ERROR: Unknown model.
    exit 1
elif [[ $norm -lt 0 || $norm -gt 3 ]]; then
    echo ERROR: Invalid norm.
    exit 1
elif [[ $timelimit -lt 0 ]]; then
    echo ERROR: Negative time limit.
    exit 1
elif [[ $d_max -lt 0 ]]; then
    echo ERROR: Invalid d_max.
    exit 1
fi

# The dataset name is provided as an argument
dataset=$1
dataset_dir=../data/$dataset
if [ ! -d "$dataset_dir" ]; then
    echo ERROR: Dataset not found.
    exit 1
fi

# Extract the number of bins from the parameters.txt file of the dataset
bins="$(grep -Po '^num_bins=\K\S+' $dataset_dir/parameters.txt)"

# Store all the instance names in an array
mapfile -t instances < $dataset_dir/instances.txt

# Create the results directory
results_dir=../results/$dataset-norm$norm-tl$timelimit-dmax$d_max/$model
if [ -d "$results_dir" ]; then
    echo ERROR: Directory $results_dir already exists.
    exit 1
fi
mkdir -p $results_dir

# Proceed with all tests
for instance in "${instances[@]}"
do
    instance_basename=`echo "$instance" | cut -d'.' -f1`
    for deviation in $( seq 1 $d_max )
    do
	echo Solving $instance \($((deviation-1))-$deviation\)
	if [ "$model" == "ipa" ]; then
	    ./../models/ipa/ipa -file $dataset_dir/$instance \
				-bins $bins \
				-norm $norm \
				-timelimit $timelimit \
				-dmin $((deviation-1)) \
				-dmax $deviation \
				> $results_dir/$instance_basename\_$deviation.csv
	elif [ "$model" == "ipb" ]; then
	    ./../models/ipb/ipb -file $dataset_dir/$instance \
				-bins $bins \
				-norm $norm \
				-timelimit $timelimit \
				-dmin $((deviation-1)) \
				-dmax $deviation \
				-cp \
				> $results_dir/$instance_basename\_$deviation.csv
	fi
    done
done

echo Done testing model $model with norm $norm, time limit $timelimit, and \
     d_max $d_max on dataset $dataset.
