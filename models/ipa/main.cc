/*
 * This is a modified version of the IP model introduced by Lewis and Carroll
 * in their 2016 paper: "Creating Seating Plans: A Practical Application".
 * 
 * In this modified model the deviation can be bounded, which is needed to
 * construct the Pareto set of solutions. This model also considers multiple
 * load balancing strategies.
 */


#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <stdlib.h>
#define RC_EPS 1.0e-6
#define RC_EPS2 1.0e-3
#define CONFLICT 9999 // IloIntMax causes overflow
#define L0_DEVIATION 0
#define L1_DEVIATION 1
#define L2_DEVIATION 2
#define Li_DEVIATION 3
// #define MINMAX 4
// #define MAXMIN 5
// #define SPREAD 6


ILOSTLBEGIN


ILOMIPINFOCALLBACK3(TimeLimitCallback,
		    IloCplex, cplex,
		    IloNum, time_start,
		    IloNum, time_limit) {
  if (cplex.getCplexTime()-time_start > time_limit) abort();
}


void Help();


int main(int argc, char* argv[]) {
  char* filename = NULL;

  IloInt d_min = 0; // Lower bound on cumulative deviation
  IloInt d_max = IloIntMax; // Upper bound on cumulative deviation
  IloInt strategy = -1; // Load balancing strategy
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
    else if (!strcmp(argv[i], "-strategy")) {
      i++;
      strategy = atoi(argv[i]);
    }
    else if (!strcmp(argv[i], "-timelimit")) {
      i++;
      time_limit = atof(argv[i]);
    }

  }
  if ((num_bins <= 0) ||
      (d_min < 0) ||
      (d_max < d_min) ||
      (strategy < 0 || strategy > 3) ||
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
    
  // Load data from an instance
  ifstream f(filename, ios::in);
  IloInt n;
  char c;

  IloInt num_items;
  f >> num_items;

  IloIntArray weights(env);
  for (IloInt i=0; i<num_items; i++) {
    f >> n;
    weights.add(n);
    f >> c;
  }
    
  IloArray<IloIntArray> costs(env);
  for (IloInt i=0; i<num_items; i++) {
    costs.add(IloIntArray(env));
    for (IloInt j=0; j<num_items; j++) {
      f >> n;
      costs[i].add(n);
      f >> c;
    }
  }

  IloNum mean_load = (IloNum)IloSum(weights)/num_bins;

  // Optimize bin load bounds (these bounds can be inferred algebraically)
  IloInt min_load = 0;
  IloInt max_load = IloIntMax;
  if (d_max < IloIntMax) {
    if (strategy == L0_DEVIATION) {
      if (d_max == 0) {
	// We allow this floor/ceil leeway in case the mean is fractional
	min_load = floor(mean_load);
	max_load = ceil(mean_load);
      }
      else {
	min_load = 0;
	max_load = IloIntMax;
      }
    }
    else if (strategy == L1_DEVIATION) {
      min_load = max((IloInt)0,
		     (IloInt)ceil(mean_load-(IloNum)d_max/2));
      max_load = floor(mean_load+(IloNum)d_max/2);
    }
    else if (strategy == L2_DEVIATION) {
      min_load = max((IloInt)0,
		     (IloInt)ceil(mean_load-
				  sqrt(d_max*(IloNum)(num_bins-1)/num_bins)));
      max_load = (IloInt)floor(mean_load+sqrt(d_max*
					      (IloNum)(num_bins-1)/num_bins));
    }
    else if (strategy == Li_DEVIATION) {
      min_load = max((IloInt)0,
		     (IloInt)ceil(mean_load-d_max));
      max_load = (IloInt)floor(mean_load+d_max);
    }
  }

  // Optimization: Add conflicts for pairwise items of large combined weights
  for (IloInt i=0; i<num_items-1; i++) {
    for (IloInt j=i+1; j<num_items; j++) {
      if (weights[i]+weights[j] > max_load) {
	costs[i][j] = CONFLICT;
	costs[j][i] = CONFLICT;
      }
    }
  }

  try {

    /***************************************************************************
     * Barebones model
     **************************************************************************/
	
    // Decision variables: x[i][k]==1 if item i is assigned to bin k, 0 otherwise
    IloArray<IloNumVarArray> x = IloArray<IloNumVarArray>(env, num_items);
    for (IloInt i=0; i<num_items; i++) {
      x[i] = IloNumVarArray(env, num_bins, 0, 1, ILOINT); // Integrality constraints
    }

    // Constraints: An item is packed into a single bin
    for (IloInt i=0; i<num_items; i++) {
      model.add(IloSum(x[i]) == 1);
    }

    // Constraints: Conflicting items are packed into separate bins
    for (IloInt i=0; i<num_items-1; i++) {
      for (IloInt j=i+1; j<num_items; j++) {
	if (costs[i][j] == CONFLICT) {
	  for (IloInt k=0; k<num_bins; k++) {
	    model.add(x[i][k]+x[j][k] <= 1);
	  }
	}
      }
    }

    // Constraints: The load of every bin is between the minimum/maximum loads
    for (IloInt k=0; k<num_bins; k++) {
      IloNumVarArray bin_load(env);
      for (IloInt i=0; i<num_items; i++) {
	bin_load.add(x[i][k]);
      }
      model.add(IloScalProd(bin_load, weights) >= min_load);
      model.add(IloScalProd(bin_load, weights) <= max_load);
    }

    // Constraints: Symmetry breaking
    for (IloInt i=0; i<num_items; i++) {
      for (IloInt k=i+1; k<num_bins; k++) {
	model.add(x[i][k] == 0);
      }
    }

    // Auxiliary variables o denote the deviation of the bins
    // i.e., bin k has a deviation of o[k]
    IloNumVarArray o = IloNumVarArray(env, num_bins, 0, IloIntMax, ILOFLOAT);
    IloArray<IloNumExpr> w(env, num_bins);
    for (IloInt k=0; k<num_bins; k++) {
      w[k] = IloNumExpr(env);
      for (IloInt i=0; i<num_items; i++) {
	w[k] += x[i][k]*weights[i];
      }
      model.add(w[k]-mean_load <= o[k]);
      model.add(mean_load-w[k] <= o[k]);
    }

    // Objective: Sum of the costs
    IloExprArray sum_costs(env);
    for (IloInt k=0; k<num_bins; k++) {
      for (IloInt i=0; i<num_items-1; i++) {
	for (IloInt j=i+1; j<num_items; j++) {
	  sum_costs.add(x[i][k] *
			x[j][k] *
			costs[i][j]);
	}
      }
    }

    // Auxiliary variable: Global deviation
    IloNumExpr global_deviation(env);
    
    /***************************************************************************
     * L0_DEVIATION
     **************************************************************************/

    if (strategy == L0_DEVIATION) {
      IloInt mean_load_floor = floor(mean_load);
      IloInt mean_load_ceil = ceil(mean_load);
      IloIntVarArray higher_than_floor(env, num_bins, 0, 1);
      IloIntVarArray lower_than_ceil(env, num_bins, 0, 1);
      IloIntVarArray balanced_bins(env, num_bins, 0, 1);
      for (IloInt k=0; k<num_bins; k++) {
	model.add(higher_than_floor[k] == (w[k] >= mean_load_floor));
	model.add(lower_than_ceil[k] == (w[k] <= mean_load_ceil));
	model.add(balanced_bins[k] == (higher_than_floor[k]+lower_than_ceil[k] >= 2));
      }

      global_deviation = num_bins-IloSum(balanced_bins);
      model.add(global_deviation >= d_min);
      model.add(global_deviation <= d_max);
      
      objective.setExpr(IloSum(sum_costs) + RC_EPS2*global_deviation);
    }
    
    /***************************************************************************
     * L1_DEVIATION
     **************************************************************************/
    
    if (strategy == L1_DEVIATION) {
      global_deviation = IloSum(o);
      model.add(global_deviation >= d_min);
      model.add(global_deviation <= d_max);
      
      objective.setExpr(IloSum(sum_costs));
    }

    /***************************************************************************
     * L2_DEVIATION
     **************************************************************************/
    
    if (strategy == L2_DEVIATION) {
      cout << "ERROR: This model does not support L2-deviation.\n"
      	   << "With this strategy, the problem becomes a mixed-integer "
      	   << "quadratically constrained and non-convex problem, which CPLEX "
	   << "does not handle."
      	   << endl;
      exit(0);
    }

    /***************************************************************************
     * Li_DEVIATION
     **************************************************************************/
    
    if (strategy == Li_DEVIATION) {
      IloNumVar y = IloNumVar(env, d_min, d_max, ILOFLOAT);
      for (IloInt k=0; k<num_bins; k++) {
      	model.add(o[k] <= y);
      }
      global_deviation = y;
      
      objective.setExpr(IloSum(sum_costs) + RC_EPS2*global_deviation);
    }

    /***************************************************************************
     * Solving and output
     **************************************************************************/
    
    solver.use(TimeLimitCallback(env,
				 solver,
				 solver.getCplexTime(),
				 time_limit));
	
    if (!solver.solve()) {
      cout << "No solution exists." << endl;
      cout << solver.getCplexTime()-time_start << endl;
      exit(0);
    }
    
    /*
     * #bins,strategy,dmin,dmax
     * load1,...,loadk
     * costs,deviation,time
     */
    cout << num_bins
	 << ","
	 << strategy
	 << ","
	 << d_min
	 << ","
	 << d_max
	 << endl;

    for (IloInt k=0; k<num_bins; k++) {
      cout << solver.getValue(w[k]);
      (k != num_bins-1) ? cout << "," : cout << "\n";
    }

    cout << round(solver.getObjValue()) // round() to remove the RC_EPS2 artifact
	 << ","
	 << solver.getValue(global_deviation)
	 << ","
	 << solver.getCplexTime()-time_start
	 << endl;
    
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
  cout << "At a minimum, a filename, a number of bins, and a strategy must "
       << "be specified."
       << endl;
  cout << "./ipa [options]"
       << endl;
  cout << "-bins [number of bins]"
       << endl;
  cout << "-file [myinstance.wsp]"
       << endl;
  cout << "-dmin [minimum cumulative deviation (default unbounded)]\n"
       << "      (only L1- and L2-deviation take this into account)"
       << endl;
  cout << "-dmax [maximum cumulative deviation (default unbounded)]"
       << endl;
  cout << "-strategy [0 for L0-deviation,\n"
       << "           1 for L1-deviation,\n"
       << "           2 for L2-deviation,\n"
       << "           3 for Li-deviation]"
       << endl;
  cout << "-timelimit [cutoff in seconds (default unlimited)]"
       << endl;
  cout << "Example: ./ipa -file myinstance.wsp -bins 22 -strategy 1"
       << endl;
}
