#ifndef PROBLEM_H_
#define PROBLEM_H_


#include <chrono>
#include <algorithm>
#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>
#define RC_EPS 1.0e-6
#define CONFLICT 9999 // IloIntMax causes overflow


class Problem {
public:
    Problem(char* filename,
	    IloInt num_bins,
	    IloInt norm,
	    IloInt min_deviation,
	    IloInt max_deviation,
	    IloNum time_limit);
    ~Problem();

    IloNum GetLowerBound();
    IloNum GetUpperBound();

    bool IsFeasible();
    bool IsIntegral();

    IloInt GetMinDeviation();
    IloInt GetMaxDeviation();
    IloInt GetNorm();

    void DisplaySolution();

private:
    void InitializeVariables();
    void LoadData(char* filename);
    void OptimizeBinLoads();
    void GenerateInitialColumns();

    void SolveRelaxation();
    void SolveIntegrality();
    
    void SolveSubproblem();
    IloNum ComputePatternCost(IloNumArray pattern);
    IloNum ComputePatternDeviation(IloNumArray pattern);
    void BranchingHeuristic();

    // Parameters
    IloInt min_deviation_;
    IloInt max_deviation_;
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
    //IloNumExprArg master_deviation_;
    IloRangeArray x_;
    IloRange gd;
    IloRange gamma_;
    IloRange delta_;
    IloNumVarArray columns_;
    IloArray<IloNumArray> patterns_;
    IloNumArray pattern_deviations_;

    // Branch and price variables
    IloBoolArray columns_fixed_at_0_; // Current partially integral solution
    IloBoolArray columns_fixed_at_1_; // Idem
    IloInt branching_column_;
    IloNum cp_solution_;
    IloNum cp_solution_deviation_;
    IloNum lower_bound_;
    IloNum lower_bound_deviation_;
    IloNum upper_bound_;
    IloNum upper_bound_deviation_;
    int fractional_columns_;

    // Flags
    bool relaxed_solution_is_feasible_;
    bool relaxed_solution_is_integral_;
};


#endif  // PROBLEM_H_
