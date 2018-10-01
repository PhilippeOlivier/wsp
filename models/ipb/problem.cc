#include "problem.h"


ILOSTLBEGIN


Problem::Problem(char* filename,
		 IloInt num_bins,
		 IloInt norm,
		 IloInt subproblem_type,
		 IloInt d_min,
		 IloInt d_max,
		 IloNum time_limit) {
  time_start_ = chrono::high_resolution_clock::now();
    
  master_problem_ = IloModel(env_);
  master_solver_ = IloCplex(master_problem_);
  master_solver_.setOut(env_.getNullStream()); // Supress Cplex output
  master_objective_ = IloAdd(master_problem_, IloMinimize(env_));

  weights_ = IloIntArray(env_);
  costs_ = IloArray<IloIntArray>(env_);
  columns_ = IloNumVarArray(env_);
  patterns_ = IloArray<IloNumArray>(env_);
  pattern_costs_ = IloIntArray(env_);
  pattern_deviations_ = IloNumArray(env_);

  LoadData(filename);
    
  num_bins_ = num_bins;
  mean_load_ = IloSum(weights_)/(IloNum)num_bins_;
  norm_ = norm;
  subproblem_type_ = subproblem_type;
  d_min_ = d_min;
  d_max_ = d_max;
  time_limit_ = time_limit;

  min_load_ = 0;
  max_load_ = IloSum(weights_); // The highest possible maximum load
    
  if (d_max_ < IloSum(weights_)) { // If d_max is reasonably bounded
    OptimizeBinLoadBounds();
  }

  // Optimization: Add conflicts for pairwise items of large combined weights
  for (IloInt i=0; i<num_items_-1; i++) {
    for (IloInt j=i+1; j<num_items_; j++) {
      if (weights_[i]+weights_[j] > max_load_) {
	costs_[i][j] = CONFLICT;
	costs_[j][i] = CONFLICT;
      }
    }
  }

  constrain_column_at_0_ = IloConstraintArray(env_);
  constrain_column_at_1_ = IloConstraintArray(env_);

  GenerateInitialColumns();
  SolveRelaxation();
}


/******************************************************************************/


Problem::Problem(Problem& problem) {
  master_problem_ = IloModel(env_);
  master_solver_ = IloCplex(master_problem_);
  master_solver_.setOut(env_.getNullStream()); // Supress Cplex output
  master_objective_ = IloAdd(master_problem_, IloMinimize(env_));

  norm_ = problem.norm_;
  num_items_ = problem.num_items_;
  num_bins_ = problem.num_bins_;
  weights_ = problem.weights_;
  costs_ = problem.costs_;
  min_load_ = problem.min_load_;
  max_load_ = problem.max_load_;
  mean_load_ = problem.mean_load_;
  d_min_ = problem.d_min_;
  d_max_ = problem.d_max_;
  patterns_ = problem.patterns_;
  pattern_costs_ = problem.pattern_costs_;
  pattern_deviations_ = problem.pattern_deviations_;
  
  x_ = IloAdd(master_problem_, IloRangeArray(env_,
					     num_items_,
					     1,
					     1));
  
  columns_ = IloNumVarArray(env_);
  for (IloInt i=0; i<patterns_.getSize(); i++) {
    IloNum pattern_cost = ComputePatternCost(patterns_[i]);
    columns_.add(IloNumVar(master_objective_(pattern_cost)+x_(patterns_[i])));
  }

  zeta_ = IloAdd(master_problem_, IloRange(env_, 0, 0));
  if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
    delta_ = IloAdd(master_problem_, IloRange(env_, 0, 0));
    gamma_ = IloAdd(master_problem_, IloRange(env_, 0, 0));
  }
  
  lower_bound_ = problem.lower_bound_;
  lower_bound_deviation_ = problem.lower_bound_deviation_;
  SolveIntegrality();
}


/******************************************************************************/


Problem::~Problem() {
  env_.end();
}


/******************************************************************************/


void Problem::SolveWithPartiallySetColumns(std::vector<IloInt> columns_set_at_0,
					   std::vector<IloInt> columns_set_at_1) {
  UnsetAllColumns();
  SetColumns(columns_set_at_0, columns_set_at_1);
  SolveRelaxation();
}


/******************************************************************************/


IloInt Problem::GetFractionalColumn() {
  // Heuristic: Find the lowest cost fractional column
  IloInt index = -1;
  IloNum cost = IloIntMax;
  for (IloInt i=0; i<columns_.getSize(); i++) {
    if (master_solver_.getValue(columns_[i]) >= 0+RC_EPS &&
	master_solver_.getValue(columns_[i]) <= 1-RC_EPS) {
      if (cost > pattern_costs_[i]) {
      	index = i;
      	cost = pattern_costs_[i];
      }
    }
  }
  return index;
}


/******************************************************************************/


IloInt Problem::IsFeasible() {
  return is_feasible_;
}


/******************************************************************************/


IloInt Problem::GetLowerBound() {
  return lower_bound_;
}


/******************************************************************************/


IloNum Problem::GetLowerBoundDeviation() {
  return lower_bound_deviation_;
}


/******************************************************************************/


IloInt Problem::GetUpperBound() {
  return upper_bound_;
}


/******************************************************************************/


IloNum Problem::GetUpperBoundDeviation() {
  return upper_bound_deviation_;
}


/******************************************************************************/


std::vector<IloInt> Problem::GetOptimalColumns() {
  std::vector<IloInt> optimal_columns;
  for (IloInt i=0; i<columns_.getSize(); i++) {
    // IloRound() because some columns are almost at 0 but not quite
    for (IloInt j=0; j<IloRound(master_solver_.getValue(columns_[i])); j++) {
      optimal_columns.push_back(IloScalProd(patterns_[i], weights_));
    }
  }
  return optimal_columns;
}


/******************************************************************************/


// Load data from an instance
void Problem::LoadData(char* filename) {
  ifstream f(filename, ios::in);
  IloInt n;
  char c;

  f >> num_items_;
  for (IloInt i=0; i<num_items_; i++) {
    f >> n;
    weights_.add(n);
    f >> c;
  }

  for (IloInt i=0; i<num_items_; i++) {
    costs_.add(IloIntArray(env_));
    for (IloInt j=0; j<num_items_; j++) {
      f >> n;
      costs_[i].add(n);
      f >> c;
    }
  }
}


/******************************************************************************/


// Optimize bin load bounds (these bounds can be inferred algebraically)
void Problem::OptimizeBinLoadBounds() {
  if (norm_ == L0_NORM) {
    if (d_max_ == 0) {
      // We allow this floor/ceil leeway in case the mean is fractional
      min_load_ = (IloInt)IloFloor(mean_load_);
      max_load_ = (IloInt)IloCeil(mean_load_);
    }
  }
  else if (norm_ == L1_NORM) {
    min_load_ = IloMax((IloInt)0,
		       (IloInt)IloCeil(mean_load_-(IloNum)d_max_/2));
    max_load_ = (IloInt)IloFloor(mean_load_+(IloNum)d_max_/2);
  }
  else if (norm_ == L2_NORM) {
    min_load_ = IloMax((IloInt)0,
		       (IloInt)IloCeil(mean_load_-sqrt(d_max_*(IloNum)(num_bins_-1)/num_bins_)));
    max_load_ = (IloInt)IloFloor(mean_load_+sqrt(d_max_*(IloNum)(num_bins_-1)/num_bins_));
  }
  else if (norm_ == Li_NORM) {
    min_load_ = IloMax((IloInt)0,
		       (IloInt)IloCeil(mean_load_-d_max_));
    max_load_ = (IloInt)IloFloor(mean_load_+d_max_);
  }
}


/******************************************************************************/


// Generate the initial columns with a compact CP model (uses all bins)
void Problem::GenerateInitialColumns() {
  IloModel cp_model(env_);
  IloIntVarArray initial_columns(env_, num_items_, 0, num_bins_-1);

  // Constraints: Conflicting items are packed into different bins
  for (IloInt i=0; i<num_items_-1; i++) {
    for (IloInt j=i+1; j<num_items_; j++) {
      if (costs_[i][j] >= CONFLICT) {
	cp_model.add(initial_columns[i] != initial_columns[j]);
      }
    }
  }
    
  IloIntVarArray bin_loads = IloIntVarArray(env_,
					    num_bins_,
					    min_load_,
					    max_load_);

  // Constraints: Bin loads are within bounds
  cp_model.add(IloPack(env_,
		       bin_loads,
		       initial_columns,
		       weights_,
		       IloIntExpr(env_, num_bins_)));

  // Constraints: Symmetry breaking
  for (IloInt k=0; k<num_bins_-1; k++) {
    cp_model.add(bin_loads[k] >= bin_loads[k+1]);
  }

  // Auxiliary variables: Bin deviations
  bin_deviations_ = IloNumExprArray(env_);
  for (IloInt i=0; i<num_bins_; i++) {
    bin_deviations_.add(IloAbs(bin_loads[i]-mean_load_));
  }

  // Constraints: Cumulative deviation is within d_min and d_max bounds
  if (norm_ == L0_NORM) {
    IloIntVarArray load_not_in_mean = IloIntVarArray(env_,
						     num_bins_,
						     0,
						     1);
    for (IloInt k=0; k<num_bins_; k++) {
      cp_model.add(load_not_in_mean[k] == (
		     (bin_loads[k] < (IloInt)IloFloor(mean_load_)) ||
		     (bin_loads[k] > (IloInt)IloCeil(mean_load_))));
    }
    cp_model.add(IloSum(load_not_in_mean) >= d_min_);
    cp_model.add(IloSum(load_not_in_mean) <= d_max_);
  }
  else if (norm_ == L1_NORM) {
    cp_model.add(IloSum(bin_deviations_) >= d_min_);
    cp_model.add(IloSum(bin_deviations_) <= d_max_);
  }
  else if (norm_ == L2_NORM) {
    cp_model.add(IloScalProd(bin_deviations_, bin_deviations_) >= d_min_);
    cp_model.add(IloScalProd(bin_deviations_, bin_deviations_) <= d_max_);
  }
  else if (norm_ == Li_NORM) {
    // d_min is not considered since the master problem does not consider it
    cp_model.add(IloMax(bin_deviations_) <= d_max_);
  }

  // Random-fit decreasing search strategy (works well with decreasing weights)
  IloSearchPhase phase(env_,
		       IloSelectSmallest(IloVarIndex(env_,
						     initial_columns)),
		       IloSelectRandomValue(env_));    
  IloCP cp_solver(cp_model);
  cp_solver.setSearchPhases(phase);
  cp_solver.setParameter(IloCP::TimeLimit, time_limit_);
  cp_solver.setOut(env_.getNullStream()); // Supress Cplex output
  cp_solver.setWarning(env_.getNullStream()); // Supress Cplex warnings
    
  if (!cp_solver.solve()) {
    /*
     * No solution exists.
     * #bins,norm,dmin,dmax,SPtype,time
     */
    cout << "No solution exists." << endl;
    cout << num_bins_
	 << ","
	 << norm_
	 << ","
	 << d_min_
	 << ","
	 << d_max_
	 << ","
	 << subproblem_type_
	 << ","
	 << (chrono::duration<double>(chrono::high_resolution_clock::now()-time_start_).count())
	 << endl;
    exit(0);
  }
    
  // Constraints (master problem): Items are packed in exactly one bin
  x_ = IloAdd(master_problem_, IloRangeArray(env_,
					     num_items_,
					     1,
					     1));

  // Constraints (master problem): Exactly num_bins_ bins are used
  zeta_ = IloAdd(master_problem_, IloRange(env_,
					   num_bins_,
					   IloSum(columns_),
					   num_bins_));
  
  // Add the new columns to the model
  for (IloInt i=0; i<num_bins_; i++) {
    IloNumArray new_pattern(env_, num_items_);
    for (IloInt j=0; j<num_items_; j++) {
      if (cp_solver.getValue(initial_columns[j] == i)) {
	new_pattern[j] = 1;
      }
    }
    patterns_.add(new_pattern);
    IloNum pattern_cost = ComputePatternCost(new_pattern);
    IloNum pattern_deviation = ComputePatternDeviation(new_pattern);
    pattern_deviations_.add(pattern_deviation);
    pattern_costs_.add(pattern_cost);
    columns_.add(IloNumVar(master_objective_(pattern_cost)+x_(new_pattern)));
    master_problem_.add(columns_[i] >= 0); // Ensures the IP subproblem doesn't use negative columns
    constrain_column_at_0_.add(IloConstraint(columns_[columns_.getSize()-1] == 0));
    constrain_column_at_1_.add(IloConstraint(columns_[columns_.getSize()-1] >= 1)); // >= since an empty bin can be used many times
  }
  
  
  // Constraints (master problem, except for Li_NORM): Cumulative deviation is within d_min and d_max bounds
  if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
    // Dummy constraints must be introduced for UpdateConstraints() to work properly
    delta_ = IloAdd(master_problem_, IloRange(env_, 0, 0));
    gamma_ = IloAdd(master_problem_, IloRange(env_, 0, 0));
  }
  UpdateConstraints();
}


/******************************************************************************/
    
  
void Problem::UpdateConstraints() {
  // Update the constraints by taking into account newly added columns
  master_problem_.remove(zeta_);
  zeta_ = IloAdd(master_problem_, IloRange(env_,
					   num_bins_,
					   IloSum(columns_),
					   num_bins_));

  if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
    master_problem_.remove(delta_);
    delta_ = IloAdd(master_problem_, IloRange(env_,
					      0,
					      IloScalProd(columns_,
							  pattern_deviations_),
					      d_max_));
    master_problem_.remove(gamma_);
    gamma_ = IloAdd(master_problem_, IloRange(env_,
					      d_min_,
					      IloScalProd(columns_,
							  pattern_deviations_),
					      IloInfinity));
  }
}


/******************************************************************************/


IloBool Problem::AddNewColumn(IloNumArray new_pattern) {
  // Ensure that patterns are really integral
  for (IloInt j=0; j<num_items_; j++) {
    if ((new_pattern[j] >= 0-RC_EPS) && (new_pattern[j] <= 0+RC_EPS)) {
      new_pattern[j] = 0;
    }
    else if ((new_pattern[j] >= 1-RC_EPS) &&
	     (new_pattern[j] <= 1+RC_EPS)) {
      new_pattern[j] = 1;
    }
  }
		    
  // Prevent an existing pattern from being recreated
  IloBool already_exists = IloFalse;
  for (IloInt j=0; j<patterns_.getSize(); j++) {
    IloInt similarity = 0;
    for (IloInt k=0; k<num_items_; k++) {
      if (new_pattern[k] == patterns_[j][k]) {
	similarity++;
      }
    }
    if (similarity == num_items_) {
      already_exists = IloTrue;
      break;
    }
  }
  if (already_exists == IloFalse) {
    patterns_.add(new_pattern);
    IloNum pattern_cost = ComputePatternCost(new_pattern);
    IloNum pattern_deviation = ComputePatternDeviation(new_pattern);
    pattern_deviations_.add(pattern_deviation);
    pattern_costs_.add(pattern_cost);
    columns_.add(IloNumVar(master_objective_(pattern_cost)+x_(new_pattern)));
    master_problem_.add(columns_[columns_.getSize()-1] >= 0); // Ensures the IP subproblem doesn't use negative columns
    constrain_column_at_0_.add(IloConstraint(columns_[columns_.getSize()-1] == 0));
    constrain_column_at_1_.add(IloConstraint(columns_[columns_.getSize()-1] == 1));
    return IloTrue;
  }
  return IloFalse;
}


/******************************************************************************/


void Problem::UpdateLowerBounds() {
  // We can round up since the optimal value is always integral
  lower_bound_ = (IloInt)IloCeil(master_solver_.getObjValue());
  
  if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
    lower_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_,
								 pattern_deviations_));
  }
  else if (norm_ == Li_NORM) {
    lower_bound_deviation_ = -1;
    for (IloInt i=0; i<columns_.getSize(); i++) {
      if (master_solver_.getValue(columns_[i]) >= 1-RC_EPS) {
	lower_bound_deviation_ = IloMax(lower_bound_deviation_,
					pattern_deviations_[i]);
      }
    }
  }
}


/******************************************************************************/


void Problem::UnsetAllColumns() {
  for (IloInt i=0; i<columns_.getSize(); i++) {
    master_problem_.remove(constrain_column_at_0_[i]);
    master_problem_.remove(constrain_column_at_1_[i]);
  }
}


/******************************************************************************/


void Problem::SetColumns(std::vector<IloInt> to_be_set_at_0, std::vector<IloInt> to_be_set_at_1) {
  for (IloInt i=0; i<(IloInt)to_be_set_at_0.size(); i++) {
    master_problem_.add(constrain_column_at_0_[to_be_set_at_0[i]]);
  }
  for (IloInt i=0; i<(IloInt)to_be_set_at_1.size(); i++) {
    master_problem_.add(constrain_column_at_1_[to_be_set_at_1[i]]);
  }
}


/******************************************************************************/


void Problem::SolveRelaxation() {
  if (subproblem_type_ == IP_SUBPROBLEM) {
    SolveRelaxationIp();
  }
  else if (subproblem_type_ == CP_SUBPROBLEM1) {
    SolveRelaxationCp1();
  }
  else if (subproblem_type_ == CP_SUBPROBLEM2) {
    SolveRelaxationCp2();
  }
}


/******************************************************************************/


void Problem::SolveRelaxationIp() {
  while (IloTrue) {
    UpdateConstraints();
    master_solver_.solve();
    is_feasible_ = IloTrue;
    if (!master_solver_.isPrimalFeasible()) {
      is_feasible_ = IloFalse;
      break;
    }
    UpdateLowerBounds();
      
    if ((chrono::duration<double>(chrono::high_resolution_clock::now()-time_start_).count()) > time_limit_) break;

    IloModel sub_problem(env_);
    IloObjective sub_objective = IloAdd(sub_problem, IloMaximize(env_));
    IloNumVarArray z(env_, num_items_, 0, 1, ILOINT);
	
    // Auxiliary variable: Beta (absolute deviation)
    IloNumVar beta(env_);
    if (norm_ == L0_NORM) {
      IloInt mean_load_floor = (IloInt)IloFloor(mean_load_);
      IloInt mean_load_ceil = (IloInt)IloCeil(mean_load_);
      IloIntVar higher_than_floor(env_, 0, 1);
      IloIntVar lower_than_ceil(env_, 0, 1);
      IloIntVar unbalanced_bin(env_, 0, 1);
      IloInt M1 = IloAbs(max_load_-mean_load_floor);
      IloInt M2 = IloAbs(max_load_-mean_load_ceil);

      // RC_EPS is needed in these constraints in case the mean is integral
      sub_problem.add(IloScalProd(z, weights_) >= mean_load_floor-M1*(1-higher_than_floor));
      sub_problem.add(mean_load_floor >= IloScalProd(z, weights_)-M1*(higher_than_floor)+RC_EPS);
      sub_problem.add(IloScalProd(z, weights_) <= mean_load_ceil+M2*(1-lower_than_ceil));
      sub_problem.add(mean_load_ceil <= IloScalProd(z, weights_)+M2*(lower_than_ceil)-RC_EPS);
      sub_problem.add(higher_than_floor+lower_than_ceil-2*(1-unbalanced_bin) <= 1);
      sub_problem.add(higher_than_floor+lower_than_ceil-2*(1-unbalanced_bin) >= 0);
      sub_problem.add(beta == unbalanced_bin);
    }
    else {
      sub_problem.add((IloScalProd(z, weights_)-mean_load_) <= beta);
      sub_problem.add(mean_load_-(IloScalProd(z, weights_)) <= beta);
    }

    // Constraints (Li_NORM): beta must not be greater than d_max
    if (norm_ == Li_NORM) {
      sub_problem.add(beta <= d_max_);
    }

    // Constraints: Bin loads
    sub_problem.add(IloScalProd(weights_, z) >= min_load_);
    sub_problem.add(IloScalProd(weights_, z) <= max_load_);
	
    // Constraints: Conflicting items
    for (IloInt i=0; i<num_items_-1; i++) {
      for (IloInt j=i+1; j<num_items_; j++) {
	if (costs_[i][j] == CONFLICT) {
	  sub_problem.add(z[i]+z[j] <= 1);
	}
      }
    }
	
    // Dual values
    IloNumArray y(env_, num_items_);
    master_solver_.getDuals(y, x_);
    IloNum zeta = master_solver_.getDual(zeta_);
    IloNum delta = 0;
    IloNum gamma = 0;
    if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
      delta = master_solver_.getDual(delta_);
      gamma = master_solver_.getDual(gamma_);
    }

    IloNumExpr obj1(env_);
    obj1 = IloScalProd(z, y);

    IloNumExpr obj2(env_);
    obj2 += zeta;

    IloNumExpr obj3(env_);
    if (norm_ == L0_NORM || norm_ == L1_NORM) {
      obj3 = beta*(gamma+delta);
    }
    else if (norm_ == L2_NORM) {
      obj3 = IloPower(beta, 2)*(gamma+delta);
    }

    // There are two definitions of obj3, for the two objectives
    if ((gamma+delta) <= 0) {
      obj3 = obj3;
    }
    if ((gamma+delta) > 0) {
      obj3 = -obj3;
    }
	
    IloNumExpr obj4(env_);
    for (IloInt i=0; i<num_items_-1; i++) {
      for (IloInt j=i+1; j<num_items_; j++) {
	obj4 += z[i]*z[j]*costs_[i][j];
      }
    }

    // Objective
    if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
      sub_objective.setExpr(obj1+obj2+obj3-obj4);
    }
    else if (norm_ == Li_NORM) {
      sub_objective.setExpr(obj1+obj2-obj4);
    }

    IloCplex sub_solver(sub_problem);
    sub_solver.setOut(env_.getNullStream()); // Supress Cplex output
    sub_solver.solve();

    IloInt new_columns_added = 0;
    for (IloInt i=0; i<IloMin(sub_solver.getSolnPoolNsolns(), num_bins_); i++) {
      if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM ?
	  (sub_solver.getValue(obj1, i)+sub_solver.getValue(obj2, i)+sub_solver.getValue(obj3, i)
	   >
	   sub_solver.getValue(obj4, i)) : // Else, Li_NORM
	  (sub_solver.getValue(obj1, i)+sub_solver.getValue(obj2, i)
	   >
	   sub_solver.getValue(obj4, i))) {

	IloNumArray new_pattern(env_, num_items_);
	sub_solver.getValues(new_pattern, z, i);

	if (AddNewColumn(new_pattern)) {
	  new_columns_added++;
	}
	
	if (new_columns_added == num_bins_) {
	  break;
	}
      }
    }

    // If no new column is added, no constraints were violated
    if (new_columns_added == 0) break;
  }
}


/******************************************************************************/


void Problem::SolveRelaxationCp1() {
  while (IloTrue) {
    UpdateConstraints();
    master_solver_.solve();
    is_feasible_ = IloTrue;
    if (!master_solver_.isPrimalFeasible()) {
      is_feasible_ = IloFalse;
      break;
    }
    UpdateLowerBounds();
      
    if ((chrono::duration<double>(chrono::high_resolution_clock::now()-time_start_).count()) > time_limit_) break;

    IloModel sub_problem(env_);
    IloIntVarArray z(env_, num_items_, 0, 1);

    // Auxiliary variable: Beta (absolute deviation)
    IloNumVar beta(env_);
    if (norm_ == L0_NORM) {
      sub_problem.add((IloScalProd(z, weights_) <= IloFloor(mean_load_)-RC_EPS ||
		       IloScalProd(z, weights_) >= IloCeil(mean_load_)+RC_EPS) == beta);
    }
    else {
      sub_problem.add(IloAbs(mean_load_-IloScalProd(z, weights_)) == beta);
    }

    // Constraints (Li_NORM): beta must not be greater than d_max
    if (norm_ == Li_NORM) {
      sub_problem.add(beta <= d_max_);
    }
      
    // Constraints: Bin loads
    sub_problem.add(IloScalProd(weights_, z) >= min_load_);
    sub_problem.add(IloScalProd(weights_, z) <= max_load_);
	
    // Constraints: Conflicting items
    for (IloInt i=0; i<num_items_-1; i++) {
      for (IloInt j=i+1; j<num_items_; j++) {
	if (costs_[i][j] == CONFLICT) {
	  sub_problem.add(z[i]+z[j] <= 1);
	}
      }
    }

    // Dual values
    IloNumArray y(env_, num_items_);
    master_solver_.getDuals(y, x_);
    IloNum zeta = master_solver_.getDual(zeta_);
    IloNum delta = 0;
    IloNum gamma = 0;
    if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
      delta = master_solver_.getDual(delta_);
      gamma = master_solver_.getDual(gamma_);
    }

    IloNumExpr obj1(env_);
    obj1 = IloScalProd(z, y);

    IloNumExpr obj2(env_);
    obj2 += zeta;

    IloNumExpr obj3(env_);
    if (norm_ == L0_NORM || norm_ == L1_NORM) {
      obj3 = beta*(gamma+delta);
    }
    else if (norm_ == L2_NORM) {
      obj3 = IloPower(beta, 2)*(gamma+delta);
    }

    // There are two definitions of obj3, for the two objectives
    if ((gamma+delta) <= 0) {
      obj3 = obj3;
    }
    if ((gamma+delta) > 0) {
      obj3 = -obj3;
    }
	
    IloNumExpr obj4(env_);
    for (IloInt i=0; i<num_items_-1; i++) {
      for (IloInt j=i+1; j<num_items_; j++) {
	obj4 += z[i]*z[j]*costs_[i][j];
      }
    }

    // The objective is moved into the constraints
    if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
      sub_problem.add(obj1+obj2+obj3 >= obj4+RC_EPS);
    }
    else if (norm_ == Li_NORM) {
      sub_problem.add(obj1+obj2 >= obj4+RC_EPS);
    }

    IloCP sub_solver(sub_problem);
    sub_solver.setOut(env_.getNullStream()); // Supress Cplex output
    sub_solver.startNewSearch();
    
    IloInt new_columns_added = 0;
    while (sub_solver.next()) {
      IloNumArray new_pattern(env_, num_items_);
      sub_solver.getValues(z, new_pattern);
      if (AddNewColumn(new_pattern)) {
	new_columns_added++;
      }
      if (new_columns_added == num_bins_) {
	break;
      }
    }

    // If no new column is added, no constraints were violated
    if (new_columns_added == 0) {
      sub_solver.endSearch();
      break;
    }
  }
}


/******************************************************************************/


void Problem::SolveRelaxationCp2() {
  while (IloTrue) {
    UpdateConstraints();
    master_solver_.solve();
    is_feasible_ = IloTrue;
    if (!master_solver_.isPrimalFeasible()) {
      is_feasible_ = IloFalse;
      break;
    }
    UpdateLowerBounds();
      
    if ((chrono::duration<double>(chrono::high_resolution_clock::now()-time_start_).count()) > time_limit_) break;

    IloModel sub_problem(env_);
    IloInt num_items_dummy = num_items_+1;
    IloInt z_size = 0;
    IloInt cumulative_weight = 0;
    // Ensure z_size is as small as possible
    for (IloInt i=num_items_-1; i>=0; i--) {
      cumulative_weight += weights_[i];
      z_size++;
      if (cumulative_weight >= max_load_) {
	break;
      }
    }
    IloIntVarArray z(env_, z_size, 0, num_items_dummy-1);

    // Weights and costs arrays with a weightless, costless, and relationless dummy item added
    IloIntArray weights_dummy(env_, num_items_dummy);
    for (IloInt i=0; i<num_items_dummy-1; i++) {
      weights_dummy[i] = weights_[i];
    }
    weights_dummy[num_items_dummy-1] = 0; // Dummy item
	
    IloArray<IloIntArray> costs_dummy(env_, num_items_dummy);
    for (IloInt i=0; i<num_items_; i++) {
      costs_dummy[i] = IloIntArray(env_, num_items_dummy);
      for (IloInt j=0; j<num_items_; j++) {
	costs_dummy[i][j] = costs_[i][j];
      }
    }
    costs_dummy[num_items_dummy-1] = IloIntArray(env_, num_items_dummy);
    for (IloInt i=0; i<num_items_dummy; i++) {
      costs_dummy[i][num_items_dummy-1] = 0;
      costs_dummy[num_items_dummy-1][i] = 0;
    }

    // Auxiliary variable: Beta (absolute deviation)
    IloIntExpr bin_weight(env_);
    for (IloInt i=0; i<z_size; i++) {
      bin_weight += weights_dummy[z[i]];
    }
    IloNumVar beta(env_);
    if (norm_ == L0_NORM) {
      sub_problem.add((bin_weight <= IloFloor(mean_load_)-RC_EPS ||
		       bin_weight >= IloCeil(mean_load_)+RC_EPS) == beta);
    }
    else {
      sub_problem.add(IloAbs(mean_load_-bin_weight) == beta);
    }

    // Constraints (Li_NORM): beta must not be greater than d_max
    if (norm_ == Li_NORM) {
      sub_problem.add(beta <= d_max_);
    }
	
    // Constraints: Bin loads
    sub_problem.add(bin_weight >= min_load_);
    sub_problem.add(bin_weight <= max_load_);

    // Constraints: Items can be packed only once (except the dummy item)
    IloIntVarArray cards = IloIntVarArray(env_, num_items_, 0, 1);
    cards.add(IloIntVar(env_, 0, num_items_));
    IloIntArray values = IloIntArray(env_, num_items_dummy);
    for (IloInt i=0; i<num_items_dummy; i++) {
      values[i] = i;
    }
    sub_problem.add(IloDistribute(env_, cards, values, z));

    // Symmetry breaking
    for (IloInt i=0; i<z_size-1; i++) {
      sub_problem.add(z[i] <= z[i+1]);
    }
	
    // Constraints: Conflicting items
    // Table constraint: the list of forbidden assignments
    IloIntTupleSet table = IloIntTupleSet(env_, 2);
    for (IloInt i=0; i<num_items_-1; i++) {
      for (IloInt j=i+1; j<num_items_; j++) {
	if (costs_[i][j] == CONFLICT) {
	  table.add(IloIntArray(env_, 2, i, j));
	  table.add(IloIntArray(env_, 2, j, i));
	}
      }
    }
	
    // Make sure pairwise variables are NOT in the table of forbidden assignments
    for (IloInt i=0; i<z_size-1; i++) {
      for (IloInt j=i+1; j<z_size; j++) {
	sub_problem.add(IloForbiddenAssignments(env_, z[i], z[j], table));
      }
    }

    // Dual values
    IloNumArray y_dummy(env_, num_items_);
    master_solver_.getDuals(y_dummy, x_);
    y_dummy.add(0); // Dummy item
    IloNum zeta = master_solver_.getDual(zeta_);
    IloNum delta = 0;
    IloNum gamma = 0;
    if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
      delta = master_solver_.getDual(delta_);
      gamma = master_solver_.getDual(gamma_);
    }

    IloNumExpr obj1(env_);
    for (IloInt i=0; i<z_size; i++) {
      obj1 += y_dummy[z[i]];
    }

    IloNumExpr obj2(env_);
    obj2 += zeta;

    IloNumExpr obj3(env_);
    if (norm_ == L0_NORM || norm_ == L1_NORM) {
      obj3 = beta*(gamma+delta);
    }
    else if (norm_ == L2_NORM) {
      obj3 = IloPower(beta, 2)*(gamma+delta);
    }

    // There are two definitions of obj3, for the two objectives
    if ((gamma+delta) <= 0) {
      obj3 = obj3;
    }
    if ((gamma+delta) > 0) {
      obj3 = -obj3;
    }

    IloIntVarArray together = IloIntVarArray(env_, num_items_, 0, 1);
    together.add(IloIntVar(env_, 0, num_items_)); // Dummy item
    sub_problem.add(IloDistribute(env_, together, z));

    IloNumExpr obj4(env_);
    for (IloInt i=0; i<num_items_dummy; i++) {
      obj4 += IloScalProd(together, costs_dummy[i])*together[i];
    }

    // The objective is moved into the constraints
    if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
      sub_problem.add(obj1+obj2+obj3 >= obj4+RC_EPS);
    }
    else if (norm_ == Li_NORM) {
      sub_problem.add(obj1+obj2 >= obj4+RC_EPS);
    }

    IloCP sub_solver(sub_problem);
    sub_solver.setOut(env_.getNullStream()); // Supress Cplex output
    sub_solver.startNewSearch();
      
    IloInt new_columns_added = 0;
    while (sub_solver.next()) {
      IloNumArray new_pattern(env_, num_items_);
      for (IloInt i=0; i<num_items_; i++) {
	new_pattern[i] = 0;
      }
      for (IloInt i=0; i<z_size; i++) {
	IloInt item = sub_solver.getValue(z[i]);
	if (item != (num_items_dummy-1)) {
	  new_pattern[item] = 1;
	}
      }
      if (AddNewColumn(new_pattern)) {
	new_columns_added++;
      }
      if (new_columns_added == num_bins_) {
	break;
      }
    }

    // If no new column is added, no constraints were violated
    if (new_columns_added == 0) {
      sub_solver.endSearch();
      break;
    }
  }
}


/******************************************************************************/


void Problem::SolveIntegrality() {
  master_problem_.add(IloConversion(env_, columns_, ILOINT));
  UpdateConstraints();
  master_solver_.solve();
  
  upper_bound_ = (IloInt)master_solver_.getObjValue();
  
  if (norm_ == L0_NORM || norm_ == L1_NORM || norm_ == L2_NORM) {
    upper_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_,
								 pattern_deviations_));
  }
  else if (norm_ == Li_NORM) {
    upper_bound_deviation_ = -1;
    for (IloInt i=0; i<columns_.getSize(); i++) {
      if (master_solver_.getValue(columns_[i]) >= 1-RC_EPS) {
	upper_bound_deviation_ = IloMax(lower_bound_deviation_,
					pattern_deviations_[i]);
      }
    }
  }
}


/******************************************************************************/


IloNum Problem::ComputePatternCost(IloNumArray pattern) {
  IloNum cost = 0;
  for (IloInt i=0; i<num_items_; i++) {
    for (IloInt j=i+1; j<num_items_; j++) {
      cost += pattern[i]*pattern[j]*costs_[i][j];
    }
  }
  
  return cost;
}


/******************************************************************************/


IloNum Problem::ComputePatternDeviation(IloNumArray pattern) {
  IloNum weight = IloScalProd(weights_, pattern);
  IloNum deviation = IloAbs(mean_load_-weight);

  if (norm_ == L0_NORM) {
    if (weight >= IloFloor(mean_load_) &&
	weight <= IloCeil(mean_load_)) {
      deviation = 0;
    }
    else {
      deviation = 1;
    }
  }
  else if (norm_ == L1_NORM || norm_ == Li_NORM) {
    deviation = deviation;
  }
  else if (norm_ == L2_NORM) {
    deviation = deviation*deviation;
  }
  
  return deviation;
}
