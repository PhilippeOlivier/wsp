#include <iostream>
#include <stack>
#include <vector>
#include "node.h"
#include "problem.h"


ILOSTLBEGIN


void Help();


int main(int argc, char* argv[]) {
  char* filename = nullptr;
    
  IloInt d_min = 0; // Lower bound on cumulative deviation
  IloInt d_max = IloIntMax; // Upper bound on cumulative deviation
  IloInt norm = -1; // Load balancing norm
  IloInt num_bins = -1;
  IloInt subproblem_type = -1;
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
    else if (!strcmp(argv[i], "-ip")) {
      subproblem_type = IP_SUBPROBLEM;
    }
    else if (!strcmp(argv[i], "-cp1")) {
      subproblem_type = CP_SUBPROBLEM1;
    }
    else if (!strcmp(argv[i], "-cp2")) {
      subproblem_type = CP_SUBPROBLEM2;
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
      (subproblem_type == -1) ||
      (time_limit < 0)) {
    Help();
    exit(0);
  }

  chrono::time_point<chrono::high_resolution_clock> time_start = chrono::high_resolution_clock::now();
  
  Problem *problem = new Problem(filename,
				 num_bins,
				 norm,
				 subproblem_type,
				 d_min,
				 d_max,
				 time_limit);

  // Generate an initial upper bound
  Problem *initial_upper_bound = new Problem(*problem);
  std::vector<IloInt> optimal_columns = initial_upper_bound->GetOptimalColumns();
  IloInt lower_bound = initial_upper_bound->GetLowerBound();
  IloNum lower_bound_deviation = initial_upper_bound->GetLowerBoundDeviation();
  IloInt upper_bound = initial_upper_bound->GetUpperBound();
  IloNum upper_bound_deviation = initial_upper_bound->GetUpperBoundDeviation();
  delete initial_upper_bound;

  // Root node
  std::stack<Node*> stack;
  Node *node = new Node();
  node->Solve(*problem);
  stack.push(node);
  
  // Branch and price
  IloBool solution_is_optimal = IloTrue;
  while (!stack.empty()) {
    if (chrono::duration<double>(chrono::high_resolution_clock::now()-time_start).count() >= time_limit) {
      solution_is_optimal = IloFalse;
      break;
    }
    node = stack.top();
    stack.pop();
    
    // 1.A: Node is feasible
    if (node->IsFeasible()) {

      // 2.A: Node is integral (update upper bound?)
      if (node->IsIntegral()) {
    	if (node->GetObjectiveValue() < upper_bound) { // A better integral solution has been found
    	  optimal_columns = problem->GetOptimalColumns();
    	  upper_bound = node->GetObjectiveValue();
    	  upper_bound_deviation = node->GetDeviation();
    	}
      }
      
      // 2.B: Node is not integral (prune or split?)
      else {
	if (node->GetObjectiveValue() < upper_bound) {
	  Node *left_node = new Node(*node, 0);
	  left_node->Solve(*problem);
	  stack.push(left_node);
	  Node *right_node = new Node(*node, 1);
	  right_node->Solve(*problem);
	  stack.push(right_node);
	}
      }
    }
    
    // 1.B: Node is infeasible
    else {
      ;
    }

    /* 
     * At this point the node is either:
     * - Infeasible, or integral/fractional but worse than the upper bound (prune)
     * - Fractional and better than the upper bound (split)
     * - Integral and better than the upper bound (update)
     * In every case, the current node is spent and can be discarded.
     */
    delete node;
  }

  delete problem;

  /*
   * #bins,norm,dmin,dmax,SPtype
   * load1,...,loadk
   * LBcosts,LBdeviation,UBcosts,UBdeviation,time
   */
  cout << num_bins
       << ","
       << norm
       << ","
       << d_min
       << ","
       << d_max
       << ","
       << subproblem_type
       << endl;

  IloBool first_bin = IloTrue;
  for (IloInt i=0; i<(IloInt)optimal_columns.size(); i++) {
    if (first_bin == IloFalse) {
      cout << ",";
    }
    cout << optimal_columns[i];
    first_bin = IloFalse;
  }
  cout << endl;
  
  cout << (solution_is_optimal ? upper_bound : lower_bound)
       << ","
       << (solution_is_optimal ? upper_bound_deviation : lower_bound_deviation)
       << ","
       << upper_bound
       << ","
       << upper_bound_deviation
       << ","
       << (chrono::duration<double>(chrono::high_resolution_clock::now()-time_start).count())
       << endl;
  
  return 0;
}


void Help() {
  cout << "At a minimum, a filename, a number of bins, a norm, "
       << "and a subproblem type must be specified."
       << endl;
  cout << "./ipb [options]"
       << endl;
  cout << "-bins [number of bins]"
       << endl;
  cout << "-file [myinstance.wsp]"
       << endl;
  cout << "-dmin [minimum cumulative deviation (default unbounded)]\n"
       << "      (Li-deviation does not take this into account)"
       << endl;
  cout << "-dmax [maximum cumulative deviation (default unbounded)]"
       << endl;
  cout << "-norm [0 for L0-deviation,\n"
       << "       1 for L1-deviation,\n"
       << "       2 for L2-deviation,\n"
       << "       3 for Li-deviation]"
       << endl;
  cout << "-ip [for an IP-based subproblem]"
       << endl;
  cout << "-cp1 [for a CP-based subproblem with 0-1 variables]"
       << endl;
  cout << "-cp2 [for a CP-based subproblem with a more natural CP model]"
       << endl;
  cout << "-timelimit [cutoff in seconds (default unlimited)]"
       << endl;
  cout << "Example: ./ipb -file myinstance.wsp -bins 22 -norm 2 -cp1"
       << endl;
}
