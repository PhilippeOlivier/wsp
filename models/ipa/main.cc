/*
 * This is a modified version of the IP model introduced by Lewis and Carroll
 * in their 2016 paper: "Creating Seating Plans: A Practical Application".
 * 
 * In this modified model the deviation can be bounded, which is needed to
 * construct the Pareto set of solutions. This model also consider four
 * different load balancing norms.
 */


#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <stdlib.h>
#define RC_EPS 1.0e-6
#define EPSILON_2 0.0001
#define EPSILON_1 0.01
#define CONFLICT 9999 // IloIntMax causes overflow
#define L0_NORM 0
#define L1_NORM 1
#define L2_NORM 2
#define Li_NORM 3
ILOSTLBEGIN


ILOMIPINFOCALLBACK3(TimeLimitCallback,
		    IloCplex, cplex,
		    IloNum, time_start,
		    IloNum, time_limit) {
    if (cplex.getCplexTime()-time_start > time_limit) abort();
}


void Help();


int main(int argc, char* argv[]) {
    char* filename;

    IloInt d_min = 0;
    IloInt d_max = IloIntMax;
    IloInt norm = -1;
    IloInt num_bins = -1;
    IloNum time_limit = IloNumMax;
    
    if (argc == 0) {
	Help();
	exit(0);
    }
    for (int i=1; i<argc; i++) {
	if (!strcmp(argv[i], "-bins")) {
	    i++;
	    num_bins = atoi(argv[i]);
	}
	else if (!strcmp(argv[i], "-file")) {
	    i++;
	    filename = argv[i];
	}
	else if (!strcmp(argv[i], "-dmin")) {
	    i++;
	    d_min = atoi(argv[i]);
	}
	else if (!strcmp(argv[i], "-dmax")) {
	    i++;
	    d_max = atoi(argv[i]);
	}
	else if (!strcmp(argv[i], "-norm")) {
	    i++;
	    norm = atoi(argv[i]);
	}
	else if (!strcmp(argv[i], "-timelimit")) {
	    i++;
	    time_limit = atof(argv[i]);
	}

    }
    if ((num_bins <= 0) ||
	(d_min < 0) ||
	(d_max < d_min) ||
	(norm < 0 || norm > 3) ||
	(time_limit < 0)) {
	Help();
	exit(0);
    }

    IloEnv env;
    IloModel model(env);
    IloCplex solver(model);
    IloObjective objective = IloAdd(model, IloMinimize(env));
    IloNum time_start = solver.getCplexTime();
    solver.setOut(env.getNullStream()); // Supress Cplex output
    
    // Load data from a .wsp instance
    ifstream f(filename, ios::in);
    int n;
    char c;

    IloInt num_items;
    f >> num_items;

    IloIntArray weights(env);
    for (int i=0; i<num_items; i++) {
	f >> n;
	weights.add(n);
	f >> c;
    }
    
    IloArray<IloIntArray> costs(env);
    for (int i=0; i<num_items; i++) {
	costs.add(IloIntArray(env));
	for (int j=0; j<num_items; j++) {
	    f >> n;
	    costs[i].add(n);
	    f >> c;
	}
    }

    IloNum mean_load = IloSum(weights)/(float)num_bins;
    IloInt min_load = 0;
    IloInt max_load = IloIntMax;

    // Optimize bin load bounds
    if (d_max < IloIntMax) {
	if (norm == L0_NORM) {
	    if (d_max == 0) {
		min_load = (int)floor(mean_load);
		max_load = (int)ceil(mean_load);
	    }
	    else {
		min_load = 0;
		max_load = IloIntMax;
	    }
	}
	else if (norm == L1_NORM) {
	    min_load = (int)max((double)0,
				(double)ceil(mean_load-(float)d_max/2));
	    max_load = (int)floor(mean_load+(float)d_max/2);
	}
	else if (norm == L2_NORM) {
	    min_load = (int)max((double)0,
				(double)ceil(mean_load-
					     sqrt(d_max*
						  (num_bins-1)/num_bins)));
	    max_load = (int)floor(mean_load+sqrt(d_max*
						 (num_bins-1)/num_bins));
	}
	else if (norm == Li_NORM) {
	    min_load = (int)max((double)0,
				(double)ceil(mean_load-d_max));
	    max_load = (int)floor(mean_load+d_max);
	}
    }

    // Add a conflict for pairwise items of large combined weights
    for (int i=0; i<num_items-1; i++) {
	for (int j=i+1; j<num_items; j++) {
	    if (weights[i]+weights[j] > max_load) {
		costs[i][j] = CONFLICT;
		costs[j][i] = CONFLICT;
	    }
	}
    }


// TODO: remove the try catch once this works.
    try {

	/*
	 * Barebones model
	 */
	
	// Decision variables (9), constraints (20)
	// x[i][k] == 1 if item i is assigned to bin k
	IloArray<IloNumVarArray> x = IloArray<IloNumVarArray>(env, num_items);
	for (int i=0; i<num_items; i++) {
	    x[i] = IloNumVarArray(env, num_bins, 0, 1, ILOINT);
	}

	
	// Constraints (11)
	// An item is packed into a single bin
	for (int i=0; i<num_items; i++) {
	    model.add(IloSum(x[i]) == 1);
	}

	// Constraints (12)
	// Conflicting items are packed into separate bins
	for (int i=0; i<num_items-1; i++) {
	    for (int j=i+1; j<num_items; j++) {
		if (costs[i][j] == CONFLICT) {
		    for (int k=0; k<num_bins; k++) {
			model.add(x[i][k]+x[j][k] <= 1);
		    }
		}
	    }
	}

	// Constraints (13), (14), (21)
	// The load of every bin is between l and u
	for (int k=0; k<num_bins; k++) {
	    IloNumVarArray bin_load(env);
	    for (int i=0; i<num_items; i++) {
		bin_load.add(x[i][k]);
	    }
	    model.add(IloScalProd(bin_load, weights) >= min_load);
	    model.add(IloScalProd(bin_load, weights) <= max_load);
	    bin_load.end();
	}

	// Constraints (19)
	// Symmetry breaking
	for (int i=0; i<num_items; i++) {
	    for (int k=i+1; k<num_bins; k++) {
		model.add(x[i][k] == 0);
	    }
	}

	// Constraints (15), (16)
	// Auxiliary variable o[k] denotes the deviation of bin k
	IloNumVarArray o = IloNumVarArray(env, num_bins, 0, IloIntMax, ILOFLOAT);
	IloArray<IloNumExpr> w(env, num_bins);
	for (int k=0; k<num_bins; k++) {
	    w[k] = IloNumExpr(env);
	    for (int i=0; i<num_items; i++) {
		w[k] += x[i][k]*weights[i];
	    }
	    model.add(w[k]-mean_load <= o[k]);
	    model.add(mean_load-w[k] <= o[k]);
	}

	/*
	 * Load balancing norms
	 */
	
	IloNumExpr global_deviation(env, 0);
	IloNumExpr focc(env, 0);
	// TODO: Constraints (X), (X)....
	if (norm == L0_NORM) {
	    IloInt mean_load_floor = floor(mean_load);
	    IloInt mean_load_ceil = ceil(mean_load);
	    IloArray<IloIntVarArray> b = IloArray<IloIntVarArray>(env, num_bins);
	    for (int k=0; k<num_bins; k++) {
		b[k] = IloIntVarArray(env, 20, 0, 1); //TODO: replace 20 by something more sensible, representing the maximum possible deviation
		model.add(IloSum(b[k]) == 1);
	    }
	    for (int k1=0; k1<num_bins; k1++) {
		for (int k2=0; k2<20; k2++) { //TODO: 20
		    model.add(k2-o[k1]+1 <= b[k1][k2]);
		}
	    }

	    for (int k1=0; k1<num_bins; k1++) {
		for (int k2=0; k2<20; k2++) { //TODO:20
		    focc += b[k1][k2]*k2;
		}
	    }

	    for (int k=0; k<num_bins; k++) {
		global_deviation += b[k][mean_load_floor];
		global_deviation += b[k][mean_load_ceil];
	    }
	    model.add(global_deviation >= d_min);
	    model.add(global_deviation <= d_max);
	}
	
	// Constraints (17), (18)
	else if (norm == L1_NORM) {
	    global_deviation = IloSum(o);
	    model.add(global_deviation >= d_min);
	    model.add(global_deviation <= d_max);
	}

	// TODO: Constraints (X), (Y), ...
	else if (norm == L2_NORM) {
	    IloNum o_min = 0;
	    IloNum o_max = sqrt(d_max*(num_bins-1)/num_bins);
	    
	    IloArray<IloIntVarArray> y = IloArray<IloIntVarArray>(env, num_bins);
	    for (int k=0; k<num_bins; k++) {
	    	y[k] = IloIntVarArray(env, num_bins, 0, IloIntMax);
	    }

	    for (int k1=0; k1<num_bins; k1++) {
		for (int k2=0; k2<num_bins; k2++) {
		    model.add(y[k1][k2] >= o_min*o[k2]+o[k1]*o_min-o_min*o_min);
		    model.add(y[k1][k2] >= o_max*o[k2]+o[k1]*o_max-o_max*o_max);
		    model.add(y[k1][k2] <= o_max*o[k2]+o[k1]*o_min-o_max*o_min);
		    model.add(y[k1][k2] <= o[k1]*o_max+o_min*o[k2]-o_min*o_max);
		    global_deviation += y[k1][k2];
		}
	    }
	    model.add(global_deviation >= d_min); //TODO: still some problems here, try vsmall with dmax 50 and dmax 35, see the problem
	    model.add(global_deviation <= d_max);
	}

	// TODO: Constraints (X), (Y), ...
	else if (norm == Li_NORM) {
	    IloNumVar y = IloNumVar(env, d_min, d_max, ILOFLOAT);
	    for (int k=0; k<num_bins; k++) {
		model.add(o[k] <= y);
	    }
	    global_deviation = y;
	}
	
	// Objective (10)
	IloExprArray sum_costs(env);
	for (int k=0; k<num_bins; k++) {
	    for (int i=0; i<num_items-1; i++) {
		for (int j=i+1; j<num_items; j++) {
		    sum_costs.add(x[i][k] *
				  x[j][k] *
				  costs[i][j]);
		}
	    }
	}
	objective.setExpr(IloSum(sum_costs)+EPSILON_2*IloSum(o)/*+EPSILON_1*global_deviation*/);
	if (norm == L0_NORM) {
	    objective.setExpr(IloSum(sum_costs)+EPSILON_2*IloSum(o)+EPSILON_2*focc+EPSILON_1*global_deviation);
	}
	sum_costs.end();

	// Solve
	solver.use(TimeLimitCallback(env,
				     solver,
				     solver.getCplexTime(),
				     time_limit));
	
	if (!solver.solve()) {
	    std::cout << "No solution exists." << std::endl;
	    std::cout << solver.getCplexTime()-time_start << std::endl;
	    exit(0);
	}

        //TEMP
	cout <<"focc: " <<solver.getValue(focc) << endl;
	
	/*
	 * #bins,norm,dmin,dmax
	 * load1,...,loadk
	 * costs,deviation,time
	 */
	std::cout << num_bins
		  << ","
		  << norm
		  << ","
		  << d_min
		  << ","
		  << d_max
		  << std::endl;

	for (int k=0; k<num_bins; k++) {
	    std::cout << solver.getValue(w[k]);
	    (k != num_bins-1) ? std::cout << "," : std::cout << "\n";
	}

	std::cout << solver.getObjValue()-EPSILON_2*solver.getValue(IloSum(o))-EPSILON_1*solver.getValue(global_deviation) //TODO: round?
		  << ","
		  << solver.getValue(global_deviation)
		  << ","
		  << solver.getCplexTime()-time_start
		  << std::endl;

	

	// TEMP
	float mandev = 0;
	for (int k=0; k<num_bins; k++) {
	    float devk = solver.getValue(o[k]);
	    cout << devk << endl;
	    if (norm == 1) mandev += devk;
	    if (norm == 2) mandev += devk*devk;
	    if (norm == 3) {
		if (devk > mandev) mandev = devk;
	    }
	}
	cout << "Manual dev: " << mandev << endl;


	
    
    }
    
    catch (IloException& e) {
	cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
	cerr << "Unknown exception caught" << endl;
    }
	
    env.end();
    return 0;
}


void Help() {
    std::cout << "At a minimum, a filename, a number of bins, and a norm must "
	      << "be specified."
	      << std::endl;
    std::cout << "./ipa [options]"
	      << std::endl;
    std::cout << "-bins [number of bins]"
	      << std::endl;
    std::cout << "-file [myinstance.wsp]"
	      << std::endl;
    std::cout << "-dmin [minimum cumulative deviation (default unbounded)]"
	      << std::endl;
    std::cout << "-dmax [maximum cumulative deviation (default unbounded)]"
	      << std::endl;
    std::cout << "-norm [1 for L1-deviation, 2 for L2-deviation, 3 for "
	      << "Li-deviation]"
	      << std::endl;
    std::cout << "-timelimit [cutoff in seconds (default unlimited)]"
	      << std::endl;
    std::cout << "Example: ./ipa -file myinstance.wsp -bins 22 -norm 2"
	      << std::endl;
}
