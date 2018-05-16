#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>
#include <iomanip>
#include <iostream>
#include "problem.h"
#include <stdlib.h>
#include <vector>
ILOSTLBEGIN


void Help();


int main(int argc, char* argv[]) {
    char* filename;
    IloInt norm = 1;
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
    if ((argc < 2) ||
	(num_bins <= 0) ||
	(min_deviation < 0) ||
	(max_deviation < min_deviation) ||
	(norm != 1 && norm != 2) ||
	(time_limit < 0)) {
	Help();
	exit(0);
    }

    
    Problem *best = new Problem(filename,
				num_bins,
				norm,
				min_deviation,
				max_deviation,
				time_limit);

    best->DisplaySolution();
    
    return 0;
}


void Help() {
    cout << "At a minimum, a filename and a number of bins must be specified."
	 << endl;
    cout << "./wsp [options]"
	 << endl;
    cout << "-bins [number of bins]"
	 << endl;
    cout << "-file [filename.wpn]"
	 << endl;
    cout << "-mindev [minimum cumulative deviation (default unbounded)]"
	 << endl;
    cout << "-maxdev [maximum cumulative deviation (default unbounded)]"
	 << endl;
    cout << "-norm [1 for L1 or 2 for L2 (default L1-norm)]"
	 << endl;
    cout << "-timelimit [cutoff in seconds (default unlimited)]"
	 << endl;
    cout << "Example: ./wsp -file myinstance.wpn -bins 22 -maxdev 7"
	 << endl;
}
