#include <iostream>
#include "problem.h"
ILOSTLBEGIN


void Help();


int main(int argc, char* argv[]) {
    char* filename;
    IloInt norm;
    IloInt subproblem_type = -1;
    IloNum time_limit = IloNumMax;
    IloInt num_bins;
    IloInt d_min = 0;
    IloInt d_max = IloIntMax;
    
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
	else if (!strcmp(argv[i], "-cp")) {
	    subproblem_type = CP_SUBPROBLEM;
	}
	else if (!strcmp(argv[i], "-ip")) {
	    subproblem_type = IP_SUBPROBLEM;
	}
	else if (!strcmp(argv[i], "-timelimit")) {
	    i++;
	    time_limit = atof(argv[i]);
	}

    }
    if ((num_bins <= 0) ||
	(d_min < 0) ||
	(d_max < d_min) ||
	(norm != 1 && norm != 2 && norm != 3) ||
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

    problem->DisplaySolution();
    
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
    cout << "-dmin [minimum cumulative deviation (default unbounded)]"
	 << endl;
    cout << "-dmax [maximum cumulative deviation (default unbounded)]"
	 << endl;
    cout << "-norm [1 for L1-deviation, 2 for L2-deviation, 3 for Li-deviation]"
	 << endl;
    cout << "-cp [for a CP-based subproblem]"
	 << endl;
    cout << "-ip [for an IP-based subproblem]"
	 << endl;
    cout << "-timelimit [cutoff in seconds (default unlimited)]"
	 << endl;
    cout << "Example: ./ipb -file myinstance.wsp -bins 22 -norm 2"
	 << endl;
}
