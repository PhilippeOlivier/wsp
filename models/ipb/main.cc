#include <iostream>
#include "problem.h"
ILOSTLBEGIN


void Help();


int main(int argc, char* argv[]) {
    char* filename;
    IloInt norm;
    IloNum time_limit = IloNumMax;
    IloInt num_bins;
    IloInt min_deviation = 0;
    IloInt max_deviation = IloIntMax;
    
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
	else if (!strcmp(argv[i], "-maxdev")) {
	    i++;
	    max_deviation = atoi(argv[i]);
	}
	else if (!strcmp(argv[i], "-mindev")) {
	    i++;
	    min_deviation = atoi(argv[i]);
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
	(min_deviation < 0) ||
	(max_deviation < min_deviation) ||
	(norm != 1 && norm != 2 && norm != 3) ||
	(time_limit < 0)) {
	Help();
	exit(0);
    }

    
    Problem *bestproblem = new Problem(filename,
				       num_bins,
				       norm,
				       min_deviation,
				       max_deviation,
				       time_limit);

    bestproblem->DisplaySolution();
    
    return 0;
}


void Help() {
    cout << "At a minimum, a filename, a number of bins, and a norm must be specified."
	 << endl;
    cout << "./ipb [options]"
	 << endl;
    cout << "-bins [number of bins]"
	 << endl;
    cout << "-file [filename.wpn]"
	 << endl;
    cout << "-mindev [minimum cumulative deviation (default unbounded)]"
	 << endl;
    cout << "-maxdev [maximum cumulative deviation (default unbounded)]"
	 << endl;
    cout << "-norm [1 for L1-deviation, 2 for L2-deviation, 3 for Li-deviation]"
	 << endl;
    cout << "-timelimit [cutoff in seconds (default unlimited)]"
	 << endl;
    cout << "Example: ./ipb -file myinstance.wpn -bins 22 -norm 2"
	 << endl;
}
