#ifndef PROBLEM_H_
#define PROBLEM_H_


#include <algorithm>
#include <chrono>
#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>
#include <vector>


#define RC_EPS 1.0e-6
#define CONFLICT 9999 // IloIntMax causes overflow
#define IP_SUBPROBLEM 1
#define CP_SUBPROBLEM1 2
#define CP_SUBPROBLEM2 3
#define L0_NORM 0
#define L1_NORM 1
#define L2_NORM 2
#define Li_NORM 3


ILOSTLBEGIN


class Problem {
public:
  Problem(char* filename,
	  IloInt num_bins,
	  IloInt norm,
	  IloInt subproblem_type,
	  IloInt d_min,
	  IloInt d_max,
	  IloNum time_limit);
  Problem(Problem& problem);
  ~Problem();

  void SolveWithPartiallySetColumns(std::vector<IloInt> columns_set_at_0,
				    std::vector<IloInt> columns_set_at_1);
  IloInt GetFractionalColumn(); // Returns -1 if no column is fractional
  IloInt IsFeasible();
  IloInt GetLowerBound();
  IloNum GetLowerBoundDeviation();
  IloInt GetUpperBound();
  IloNum GetUpperBoundDeviation();
  std::vector<IloInt> GetOptimalColumns();

private:
  void LoadData(char* filename);
  void OptimizeBinLoadBounds();
  void GenerateInitialColumns();

  IloNum ComputePatternCost(IloNumArray pattern);
  IloNum ComputePatternDeviation(IloNumArray pattern);
  void UpdateConstraints(); // Must be updated after new columns are introduced
  IloBool AddNewColumn(IloNumArray new_pattern); // Returns IloTrue if the new pattern is added as a column
  void UpdateLowerBounds();

  void UnsetAllColumns();
  void SetColumns(std::vector<IloInt> to_be_set_at_0, std::vector<IloInt> to_be_set_at_1);
  
  void SolveRelaxation();
  void SolveRelaxationIp();
  void SolveRelaxationCp1();
  void SolveRelaxationCp2();
  void SolveIntegrality();
    
  // Parameters
  IloInt d_min_;
  IloInt d_max_;
  IloInt norm_;
  IloInt subproblem_type_;
  IloNum time_limit_;
  chrono::time_point<chrono::high_resolution_clock> time_start_;
    
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
  IloNumVarArray columns_;
  IloArray<IloNumArray> patterns_;
  IloIntArray pattern_costs_;
  IloNumArray pattern_deviations_;
  
  // Common variables for all norms
  IloNumExprArray bin_deviations_;
  IloRangeArray x_;
  IloRange zeta_;
  IloRange gamma_;
  IloRange delta_;

  // Branch and price variables
  IloConstraintArray constrain_column_at_0_;
  IloConstraintArray constrain_column_at_1_;
  IloBool is_feasible_;
  
  // Solution
  IloInt lower_bound_;
  IloNum lower_bound_deviation_;
  IloInt upper_bound_;
  IloNum upper_bound_deviation_;
};


#endif  // PROBLEM_H_
