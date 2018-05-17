#ifndef PROBLEM_H_
#define PROBLEM_H_


#include <algorithm>
#include <chrono>
#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>
#define RC_EPS 1.0e-6
#define CONFLICT 9999 // IloIntMax causes overflow


class Problem {
public:
    Problem(char* filename,
	    IloInt num_bins,
	    IloInt norm,
	    IloInt d_min,
	    IloInt d_max,
	    IloNum time_limit);
    ~Problem();
    void DisplaySolution();

private:
    void InitializeVariables();
    void LoadData(char* filename);
    void OptimizeBinLoadBounds();
    void GenerateInitialColumns();

    IloNum ComputePatternCost(IloNumArray pattern);
    IloNum ComputePatternDeviation(IloNumArray pattern);
    void SolveRelaxationIp();
    void SolveRelaxationCp();
    void SolveIntegrality();
    
    // Parameters
    IloInt d_min_;
    IloInt d_max_;
    IloInt norm_;
    IloNum time_limit_;
    std::chrono::time_point<std::chrono::high_resolution_clock> time_start_;
    
    // Problem data
    IloInt num_bins_;
    IloInt num_items_;
    IloInt min_load_;
    IloInt max_load_;
    IloNum mean_load_;
    IloIntArray weights_;
    IloArray<IloIntArray> costs_;

    // Master problem and column generation variables
    IloEnv env_;
    IloModel master_problem_;
    IloCplex master_solver_;
    IloObjective master_objective_;

    IloRangeArray x_;
    IloRange zeta_;
    IloRange gamma_;
    IloRange delta_;

    IloNumVarArray columns_;
    IloArray<IloNumArray> patterns_;
    IloNumArray pattern_deviations_;

    // Output variables
    IloNum lower_bound_;
    IloNum lower_bound_deviation_;
    IloNum upper_bound_;
    IloNum upper_bound_deviation_;
};


#endif  // PROBLEM_H_
