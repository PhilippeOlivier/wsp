#include "problem.h"
ILOSTLBEGIN


// Root node initialization
Problem::Problem(char* filename,
		 IloInt num_bins,
		 IloInt norm,
		 IloInt min_deviation,
		 IloInt max_deviation,
		 bool pareto,
		 IloNum time_limit) {
    time_start_ = std::chrono::high_resolution_clock::now();
    InitializeVariables();
    LoadData(filename);
    
    num_bins_ = num_bins;
    mean_load_ = IloSum(weights_)/(float)num_bins_;
    min_load_ = 0;
    max_load_ = IloIntMax;
    norm_ = norm;
    min_deviation_ = min_deviation;
    max_deviation_ = max_deviation;
    pareto_ = pareto;
    time_limit_ = time_limit;

    // Bin load bounds are optimized if possible
    if (max_deviation_ < IloIntMax) {
	OptimizeBinLoads();
    }

    // Add a conflict for pairwise items of large combined weights
    for (int i=0; i<num_items_-1; i++) {
	for (int j=i+1; j<num_items_; j++) {
	    if (weights_[i] + weights_[j] > max_load_) {
		costs_[i][j] = CONFLICT;
		costs_[j][i] = CONFLICT;
	    }
	}
    }

    GenerateInitialColumns();
    SolveRelaxation();
    SolveIntegrality();
}


Problem::~Problem() {
    env_.end();
}


void Problem::InitializeVariables() {
    master_problem_ = IloModel(env_);
    master_solver_ = IloCplex(master_problem_);
    master_objective_ = IloAdd(master_problem_, IloMinimize(env_));
    master_solver_.setOut(env_.getNullStream()); // Supress Cplex output

    weights_ = IloIntArray(env_);
    costs_ = IloArray<IloIntArray>(env_);
    columns_ = IloNumVarArray(env_);
    patterns_ = IloArray<IloNumArray>(env_);
    pattern_deviations_ = IloNumArray(env_);

    columns_fixed_at_0_ = IloBoolArray(env_);
    columns_fixed_at_1_ = IloBoolArray(env_);
}


// Load data from a .wpn instance of the weddingseatplanner.com website
void Problem::LoadData(char* filename) {
    ifstream f(filename, ios::in);
    string line;
    int n;
    char c;

    f >> line >> n >> c >> num_items_;
    for (int i=0; i<num_items_; i++) {
	f >> line;
	weights_.add(count(line.begin(), line.end(), ','));
    }
    for (int i=0; i<num_items_; i++) {
	f >> line;
    }
    for (int i=0; i<num_items_; i++) {
	costs_.add(IloIntArray(env_));
	for (int j=0; j<num_items_; j++) {
	    f >> n;
	    if (n == 1) {
		costs_[i].add(1);
	    }
	    else if (n == 2) {
		costs_[i].add(CONFLICT);
	    }
	    else if (n == 3) {
		costs_[i].add(-1);
	    }
	    else {
		costs_[i].add(0);
	    }
	    f >> c;
	}
    }
}


void Problem::OptimizeBinLoads() {
    if (norm_ == 1) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-(float)max_deviation_/2));
	max_load_ = (int)floor(mean_load_+(float)max_deviation_/2);
    }
    else if (norm_ == 2) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-sqrt(max_deviation_)));
	max_load_ = (int)floor(mean_load_+sqrt(max_deviation_));
    }
}


void Problem::GenerateInitialColumns() {
    try {
// TODO: put a time limit here maybe
	
	// Generate the initial columns with a compact CP model (uses all bins)
	IloModel cp_model(env_);
	IloIntVarArray initial_columns(env_, num_items_, 0, num_bins_-1);
	for (int i=0; i<num_items_-1; i++) {
	    for (int j=i+1; j<num_items_; j++) {
		if (costs_[i][j] >= CONFLICT) {
		    cp_model.add(initial_columns[i] != initial_columns[j]);
		}
	    }
	}
	IloIntVarArray bin_loads = IloIntVarArray(env_,
						  num_bins_,
						  min_load_,
						  max_load_);
	cp_model.add(IloPack(env_,
			     bin_loads,
			     initial_columns,
			     weights_,
			     IloIntExpr(env_, num_bins_)));
	
	IloNumExprArray bin_deviations = IloNumExprArray(env_);
	if (norm_ == 1) {
	    for (int i=0; i<num_bins_; i++) {
		bin_deviations.add(IloAbs(bin_loads[i] - mean_load_));
	    }
	    cp_model.add(IloSum(bin_deviations) >= min_deviation_);
	    cp_model.add(IloSum(bin_deviations) <= max_deviation_);
	}
	else if (norm_ == 2) {
	    for (int i=0; i<num_bins_; i++) {
		bin_deviations.add(IloAbs(bin_loads[i]-mean_load_)*IloAbs(bin_loads[i]-mean_load_));
	    }
	    cp_model.add(IloSum(bin_deviations) >= min_deviation_);
	    cp_model.add(IloSum(bin_deviations) <= max_deviation_);
	}
	
	IloCP cp_solver(cp_model);
	cp_solver.setParameter(IloCP::TimeLimit, time_limit_);
	cp_solver.setOut(env_.getNullStream()); // Supress Cplex output
	if (!cp_solver.solve()) {
	    cout << "No solution exists." << endl;
	    cout << "Time elapsed (s): "
		 << (std::chrono::duration<double>(std::chrono::high_resolution_clock::now()
						   - time_start_).count())
		 << endl;
	    exit(0);
	}

	x_ = IloAdd(master_problem_, IloRangeArray(env_,
						   num_items_,
						   1,
						   1));

    	gamma_ = IloAdd(master_problem_, IloRange(env_,
						  min_deviation_,
						  IloScalProd(columns_, pattern_deviations_),
						  IloInfinity));

	delta_ = IloAdd(master_problem_, IloRange(env_,
						  0,
						  IloScalProd(columns_, pattern_deviations_),
						  max_deviation_));

	// Add the new columns to the model
	cp_solution_ = 0;
	cp_solution_deviation_ = 0;
	for (int i=0; i<num_bins_; i++) {
	    IloNumArray new_pattern(env_, num_items_);
	    for (int j=0; j<num_items_; j++) {
		if (cp_solver.getValue(initial_columns[j] == i)) {
		    new_pattern[j] = 1;
		}
	    }
	    patterns_.add(new_pattern);
	    IloNum pattern_cost = ComputePatternCost(new_pattern);
	    IloNum pattern_deviation = ComputePatternDeviation(new_pattern);
	    pattern_deviations_.add(pattern_deviation);
	    if (pareto_) {
		columns_.add(IloNumVar(master_objective_(pattern_cost) +
				       x_(new_pattern)));
		cp_solution_ += pattern_cost;
	    }
	    else {
		columns_.add(IloNumVar(master_objective_(pattern_cost+pattern_deviation)
				       + x_(new_pattern)));
		cp_solution_ += pattern_cost + pattern_deviation;
	    }
	    columns_fixed_at_0_.add(IloFalse);
	    columns_fixed_at_1_.add(IloFalse);
	    cp_solution_deviation_ += pattern_deviation;
	}
	upper_bound_ = cp_solution_;
	upper_bound_deviation_ = cp_solution_deviation_;
    }

    catch(IloException& e) {
	env_.out() << " ERROR: " << e << endl;
    }
}


void Problem::SolveRelaxation() {
    SolveSubproblem();

    // Select the next column to branch on
    if (relaxed_solution_is_feasible_) {
	BranchingHeuristic();
    }
}


void Problem::SolveSubproblem() {
    try {
	IloModel sub_problem(env_);
	IloCplex sub_solver(sub_problem);
	IloObjective sub_objective = IloAdd(sub_problem, IloMaximize(env_));
	sub_solver.setOut(env_.getNullStream()); // Supress Cplex output
	IloNumVarArray z(env_, num_items_, 0, 1, ILOINT);

	// Constraints: Bin loads
	sub_problem.add(IloScalProd(weights_, z) >= min_load_);
	sub_problem.add(IloScalProd(weights_, z) <= max_load_);
	
	// Constraints: Conflicting items
	for (int i=0; i<num_items_-1; i++) {
	    for (int j=i+1; j<num_items_; j++) {
		if (costs_[i][j] == CONFLICT) {
		    sub_problem.add(z[i]+z[j] <= 1);
		}
	    }
	}
	
	for (;;) {
	    relaxed_solution_is_feasible_ = true;
	    master_problem_.remove(delta_);
	    delta_ = IloAdd(master_problem_, IloRange(env_,
						      0,
						      IloScalProd(columns_, pattern_deviations_),
						      max_deviation_));
	    master_problem_.remove(gamma_);
	    gamma_ = IloAdd(master_problem_, IloRange(env_,
						      min_deviation_,
						      IloScalProd(columns_, pattern_deviations_),
						      IloInfinity));
	    master_solver_.solve();
	    if (!master_solver_.isPrimalFeasible()) {
		relaxed_solution_is_feasible_ = false;
		break;
	    }
	    lower_bound_ = master_solver_.getObjValue();
	    //cout << lower_bound_ << endl;
	    lower_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_, pattern_deviations_));
	    fractional_columns_ = 0;
	    for (int i=0; i<columns_.getSize(); i++) {
		IloNum current_column = master_solver_.getValue(columns_[i]);
		if ((current_column > 0+RC_EPS) &&
		    (current_column < 1-RC_EPS)) {
		    fractional_columns_++;
		}
	    }

	    if ((std::chrono::duration<double>(std::chrono::high_resolution_clock::now()
	    				       - time_start_).count()) > time_limit_) break;
	    
	    
	    IloNumArray price(env_, num_items_); // Duals y*
	    master_solver_.getDuals(price, x_);

	    IloNumExpr obj1(env_);
	    for (int i=0; i<num_items_-1; i++) {
		for (int j=i+1; j<num_items_; j++) {
		    obj1 += z[i]*z[j]*costs_[i][j];
		}
	    }

	    IloNumVar t(env_);
	    sub_problem.add(t >= (mean_load_ - IloScalProd(z, weights_)));
	    sub_problem.add(t >= -(mean_load_ - IloScalProd(z, weights_)));

	    IloNum gamma_dual = master_solver_.getDual(gamma_);
	    IloNum delta_dual = master_solver_.getDual(delta_);

	    IloNumExpr obj2(env_);
	    if (norm_ == 1) {
		if ((gamma_dual+delta_dual) > 0) {
		    obj2 = -t*(gamma_dual+delta_dual);
		}
		if ((gamma_dual+delta_dual) < 0) {
		    obj2 = t*(gamma_dual+delta_dual);
		}
	    }
	    else if (norm_ == 2) {
		if ((gamma_dual+delta_dual) > 0) {
		    obj2 = -t*t*(gamma_dual+delta_dual);
		}
		if ((gamma_dual+delta_dual) < 0) {
		    obj2 = t*t*(gamma_dual+delta_dual);
		}
	    }
	    
	    IloNumExpr obj3(env_);
	    obj3 = IloScalProd(z, price);

	    if (pareto_) {
		sub_objective.setExpr(obj3+obj2-obj1);
	    }
	    else {
		sub_objective.setExpr(obj3-obj2-obj1); // Not sure that this is correct
	    }
	    sub_solver.solve();

	    bool new_column_added = false;
	    for (int i=0; i<sub_solver.getSolnPoolNsolns(); i++) {
		if ((pareto_) ?
		    (sub_solver.getValue(obj1, i) <=
		     sub_solver.getValue(obj3, i)+sub_solver.getValue(obj2, i)-RC_EPS) :
		    (sub_solver.getValue(obj1, i)+sub_solver.getValue(obj2, i) <=
		     sub_solver.getValue(obj3, i)-RC_EPS)) {
		    IloNumArray new_pattern(env_, num_items_);
		    sub_solver.getValues(new_pattern, z, i);

		    // Ensure that patterns are really integral
		    for (int j=0; j<num_items_; j++) {
			if ((new_pattern[j] >= 0-RC_EPS) &&
			    (new_pattern[j] <= 0+RC_EPS)) {
			    new_pattern[j] = 0;
			}
			else if ((new_pattern[j] >= 1-RC_EPS) &&
				 (new_pattern[j] <= 1+RC_EPS)) {
			    new_pattern[j] = 1;
			}
		    }
		    
		    // Prevent an existing pattern from being recreated
		    bool already_exists = false;
		    for (int j=0; j<patterns_.getSize(); j++) {
			int similarity = 0;
			for (int k=0; k<num_items_; k++) {
			    if (new_pattern[k] == patterns_[j][k]) {
				similarity++;
			    }
			}
			if (similarity == num_items_) {
			    already_exists = true;
			    break;
			}
		    }
		    if (already_exists == false) {
			new_column_added = true;
			patterns_.add(new_pattern);
			IloNum pattern_cost = ComputePatternCost(new_pattern);
			IloNum pattern_deviation = ComputePatternDeviation(new_pattern);
			pattern_deviations_.add(pattern_deviation);
			if (pareto_) {
			    columns_.add(IloNumVar(master_objective_(pattern_cost) +
						   x_(new_pattern)));
			}
			else {
			    columns_.add(IloNumVar(master_objective_(pattern_cost +
								     pattern_deviation)
						   + x_(new_pattern)));
			}
			columns_fixed_at_0_.add(IloFalse);
			columns_fixed_at_1_.add(IloFalse);
		    }
		}
	    }
	    
	    // If no new column is added, no constraints were violated
	    if (new_column_added == false) break;
	}
    }

    catch (IloException& e) {
	cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
	cerr << "Unknown exception caught" << endl;
    }
}


IloNum Problem::ComputePatternCost(IloNumArray pattern) {
    IloNum cost = 0;
    for (int i=0; i<num_items_; i++) {
	for (int j=i+1; j<num_items_; j++) {
	    cost += pattern[i]*pattern[j]*costs_[i][j];
	}
    }
    return cost;
}


IloNum Problem::ComputePatternDeviation(IloNumArray pattern) {
    IloNum deviation = 0;
    IloNum weight = IloScalProd(weights_, pattern);
    if (norm_ == 1) {
	deviation = IloAbs(mean_load_ - weight);
    }
    else if (norm_ == 2) {
	deviation = IloAbs(mean_load_ - weight);
	deviation = deviation*deviation;
    }
    return deviation;
}


void Problem::BranchingHeuristic() {
    relaxed_solution_is_integral_ = true;
    for (int i=0; i<columns_.getSize(); i++) {
	IloNum current_column = master_solver_.getValue(columns_[i]);
	if ((columns_fixed_at_0_[i] == IloFalse) &&
	    (columns_fixed_at_1_[i] == IloFalse) &&
	    (current_column > 0+RC_EPS) &&
	    (current_column < 1-RC_EPS)) {
	    branching_column_ = i;
	    relaxed_solution_is_integral_ = false;
	    break;
	}
    }
}


void Problem::SolveIntegrality() {
    try {
    master_problem_.add(IloConversion(env_, columns_, ILOINT));
    master_problem_.add(IloSum(columns_) == num_bins_);

    master_problem_.remove(delta_);
    delta_ = IloAdd(master_problem_, IloRange(env_,
					      0,
					      IloScalProd(columns_, pattern_deviations_),
					      max_deviation_));
    master_problem_.remove(gamma_);
    gamma_ = IloAdd(master_problem_, IloRange(env_,
					      min_deviation_,
					      IloScalProd(columns_, pattern_deviations_),
					      IloInfinity));
    
    master_solver_.solve();
    if (master_solver_.isPrimalFeasible()) {
	upper_bound_ = master_solver_.getObjValue();
	upper_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_, pattern_deviations_));
    }
    else {
	upper_bound_ = IloNumMax;
	upper_bound_deviation_ = IloNumMax;
    }

    }    catch (IloException& e) {
	cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
	cerr << "Unknown exception caught" << endl;
    }
}


IloNum Problem::GetLowerBound() {
    return lower_bound_;
}


IloNum Problem::GetUpperBound() {
    return upper_bound_;
}


bool Problem::IsFeasible() {
    return relaxed_solution_is_feasible_;
}


bool Problem::IsIntegral() {
    return relaxed_solution_is_integral_;
}


IloInt Problem::GetMinDeviation() {
    return min_deviation_;
}


IloInt Problem::GetMaxDeviation() {
    return max_deviation_;
}


IloInt Problem::GetNorm() {
    return norm_;
}


bool Problem::GetPareto() {
    return pareto_;
}


void Problem::DisplaySolution() {
    /*cout << "Relaxed solution is "
	 << ((relaxed_solution_is_feasible_) ? "feasible" : "INFEASIBLE")
	 << endl;
    
    cout << "Bin load LB: "
	 << min_load_
	 << " (auto)"
	 << endl;
    cout << "Bin load UB: "
	 << max_load_
	 << " (auto)"
	 << endl;
    */
    cout << "Deviation min/max: "
	 << min_deviation_
	 << "/";
    if (max_deviation_ == IloIntMax) {
	cout << "unbounded";
    }
    else {
	cout << max_deviation_;
    }
    cout << " with L"
	 << ((norm_ == 1) ? "1" : "2")
	 << "-norm"
	 << endl;
/*
    cout << "Initial CP solution: "
	 << cp_solution_
	 << " (dev. "
	 << cp_solution_deviation_
	 << ")"
	 << endl;
*/  
    if (pareto_) {
	cout << "Pareto LB: ";
    }
    else {
	cout << "Non-Pareto LB: ";
    }
    cout << lower_bound_
	 << " (dev. "
	 << lower_bound_deviation_
	 << ")"
	 << endl;
/*
    cout << "Columns fractional/total: "
	 << fractional_columns_
	 << "/"
	 << columns_.getSize()
	 << endl;
*/
    if (pareto_) {
	cout << "Integral solution: "
	     << upper_bound_
	     << " (dev. "
	     << upper_bound_deviation_
	     << ")"
	     << endl;
    }

    cout << "Time elapsed (s): "
	 << (std::chrono::duration<double>(std::chrono::high_resolution_clock::now()
					   - time_start_).count())
	 << endl;
    
    
    // num fractional columns (make a function for this)
    // have a boolean array representing the fractional columns?
}
