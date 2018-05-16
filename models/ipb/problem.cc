#include "problem.h"
ILOSTLBEGIN


// Root node initialization
Problem::Problem(char* filename,
		 IloInt num_bins,
		 IloInt norm,
		 IloInt min_deviation,
		 IloInt max_deviation,
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
    time_limit_ = time_limit;

    // Bin load bounds are optimized if possible
    if (max_deviation_ < IloIntMax) {
	OptimizeBinLoadBounds();
    }

    // Add a conflict for pairwise items of large combined weights
    for (int i=0; i<num_items_-1; i++) {
	for (int j=i+1; j<num_items_; j++) {
	    if (weights_[i]+weights_[j] > max_load_) {
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
}


// Load data from a .wpn instance of the weddingseatplanner.com website
// TODO: change instance, accept all costs
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


void Problem::OptimizeBinLoadBounds() {
    if (norm_ == 1) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-(float)max_deviation_/2));
	max_load_ = (int)floor(mean_load_+(float)max_deviation_/2);
    }
    else if (norm_ == 2) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-
					  sqrt(max_deviation_*
					       (num_bins_-1)/num_bins_)));
	max_load_ = (int)floor(mean_load_+sqrt(max_deviation_*
					       (num_bins_-1)/num_bins_));
    }
    else if (norm_ == 3) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-max_deviation_));
	max_load_ = (int)floor(mean_load_-max_deviation_);
    }
}


// Generate the initial columns with a compact CP model (uses all bins)
void Problem::GenerateInitialColumns() {
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
	    bin_deviations.add(IloAbs(mean_load_-bin_loads[i]));
	}
	cp_model.add(IloSum(bin_deviations) >= min_deviation_);
	cp_model.add(IloSum(bin_deviations) <= max_deviation_);
    }
    else if (norm_ == 2) {
	for (int i=0; i<num_bins_; i++) {
	    bin_deviations.add(IloAbs(bin_loads[i]-mean_load_)*
			       IloAbs(bin_loads[i]-mean_load_));
	}
	cp_model.add(IloSum(bin_deviations) >= min_deviation_);
	cp_model.add(IloSum(bin_deviations) <= max_deviation_);
    }
    else if (norm_ == 3) {
	for (int i=0; i<num_bins_; i++) {
	    bin_deviations.add(IloAbs(bin_loads[i]-mean_load_));
	}
	cp_model.add(IloMax(bin_deviations) >= min_deviation_);
	cp_model.add(IloMax(bin_deviations) <= max_deviation_);
    }
	
    IloCP cp_solver(cp_model);
    cp_solver.setParameter(IloCP::TimeLimit, time_limit_);
    cp_solver.setOut(env_.getNullStream()); // Supress Cplex output
    if (!cp_solver.solve()) {
	cout << "No solution exists." << endl;
	cout << (std::chrono::duration<double>(std::chrono::high_resolution_clock::now()
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
    zeta_ = IloAdd(master_problem_, IloRange(env_,
					     num_bins_,
					     IloSum(columns_),
					     num_bins_));

    // Add the new columns to the model
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
	columns_.add(IloNumVar(master_objective_(pattern_cost) +
			       x_(new_pattern)));
    }
}


void Problem::SolveRelaxation() {
    SolveSubproblem();
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
	    master_problem_.remove(zeta_);
	    zeta_ = IloAdd(master_problem_, IloRange(env_,
						     num_bins_,
						     IloSum(columns_),
						     num_bins_));
	    master_solver_.solve();
	    if (!master_solver_.isPrimalFeasible()) {
		relaxed_solution_is_feasible_ = false;
		break;
	    }
	    lower_bound_ = master_solver_.getObjValue();
	    //cout << lower_bound_ << endl;
	    lower_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_, pattern_deviations_));

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
	    IloNum zeta_dual = master_solver_.getDual(zeta_);

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
	    obj3 = IloScalProd(z, price)+zeta_dual;

	    sub_objective.setExpr(obj3+obj2-obj1);
	    sub_solver.solve();

	    bool new_column_added = false;
	    for (int i=0; i<sub_solver.getSolnPoolNsolns(); i++) {
		if ((sub_solver.getValue(obj1, i) <=
		     sub_solver.getValue(obj3, i)+sub_solver.getValue(obj2, i)-RC_EPS)) {
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
			columns_.add(IloNumVar(master_objective_(pattern_cost) +
					       x_(new_pattern)));
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
    else if (norm_ == 3) {
	deviation = IloAbs(mean_load_ - weight);
    }
    return deviation;
}


void Problem::SolveIntegrality() {
    try {

    master_problem_.add(IloConversion(env_, columns_, ILOINT));


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
    
    	    master_problem_.remove(zeta_);
	    zeta_ = IloAdd(master_problem_, IloRange(env_,
						     num_bins_,
						     IloSum(columns_),
						     num_bins_));
    master_solver_.solve();
    // THIS IS ALWAYS feasible since we start with a feasible solution with the CP model
//    if (master_solver_.isPrimalFeasible()) {
	upper_bound_ = master_solver_.getObjValue();
	upper_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_, pattern_deviations_));
/*    }
    else {
	upper_bound_ = IloNumMax;
	upper_bound_deviation_ = IloNumMax;
	}*/

    }    catch (IloException& e) {
	cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
	cerr << "Unknown exception caught" << endl;
    }
}


void Problem::DisplaySolution() {
/*
 * LBcosts,LBdeviation,UBcosts,UBdeviation,time
 *
 */
    std::cout << lower_bound_
	      << ","
	      << lower_bound_deviation_
	      << ","
	      << upper_bound_
	      << ","
	      << upper_bound_deviation_
	      << ","
	      << (std::chrono::duration<double>(std::chrono::high_resolution_clock::now()
						- time_start_).count())
	      << std::endl;
}
