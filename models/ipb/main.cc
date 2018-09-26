#include <iostream>
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
    
  Problem *problem = new Problem(filename,
				 num_bins,
				 norm,
				 subproblem_type,
				 d_min,
				 d_max,
				 time_limit);

  // TODO: branch and price here
  problem->PrintSolution();
    
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
