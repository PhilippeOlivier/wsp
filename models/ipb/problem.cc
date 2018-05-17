#include "problem.h"
ILOSTLBEGIN


Problem::Problem(char* filename,
		 IloInt num_bins,
		 IloInt norm,
		 IloInt d_min,
		 IloInt d_max,
		 IloNum time_limit) {
    time_start_ = std::chrono::high_resolution_clock::now();
    InitializeVariables();
    LoadData(filename);
    
    num_bins_ = num_bins;
    mean_load_ = IloSum(weights_)/(float)num_bins_;
    norm_ = norm;
    d_min_ = d_min;
    d_max_ = d_max;
    time_limit_ = time_limit;

    if (d_max_ < IloIntMax) {
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
    SolveRelaxationIp();
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

    min_load_ = 0;
    max_load_ = IloIntMax;
    weights_ = IloIntArray(env_);
    costs_ = IloArray<IloIntArray>(env_);
    columns_ = IloNumVarArray(env_);
    patterns_ = IloArray<IloNumArray>(env_);
    pattern_deviations_ = IloNumArray(env_);
}


// Load data from a .wsp instance
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
			     (double)ceil(mean_load_-(float)d_max_/2));
	max_load_ = (int)floor(mean_load_+(float)d_max_/2);
    }
    else if (norm_ == 2) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-
					  sqrt(d_max_*
					       (num_bins_-1)/num_bins_)));
	max_load_ = (int)floor(mean_load_+sqrt(d_max_*
					       (num_bins_-1)/num_bins_));
    }
    else if (norm_ == 3) {
	min_load_ = (int)max((double)0,
			     (double)ceil(mean_load_-d_max_));
	max_load_ = (int)floor(mean_load_-d_max_);
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
	    bin_deviations.add(IloAbs(bin_loads[i]-mean_load_));
	}
	cp_model.add(IloSum(bin_deviations) >= d_min_);
	cp_model.add(IloSum(bin_deviations) <= d_max_);
    }
    else if (norm_ == 2) {
	for (int i=0; i<num_bins_; i++) {
	    bin_deviations.add(IloAbs(bin_loads[i]-mean_load_)*
			       IloAbs(bin_loads[i]-mean_load_));
	}
	cp_model.add(IloSum(bin_deviations) >= d_min_);
	cp_model.add(IloSum(bin_deviations) <= d_max_);
    }
    else if (norm_ == 3) {
	for (int i=0; i<num_bins_; i++) {
	    bin_deviations.add(IloAbs(bin_loads[i]-mean_load_));
	}
	cp_model.add(IloMax(bin_deviations) >= d_min_);
	cp_model.add(IloMax(bin_deviations) <= d_max_);
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

    // Master problem constraints
    x_ = IloAdd(master_problem_, IloRangeArray(env_,
					       num_items_,
					       1,
					       1));
    
    zeta_ = IloAdd(master_problem_, IloRange(env_,
					     num_bins_,
					     IloSum(columns_),
					     num_bins_));
    
    gamma_ = IloAdd(master_problem_, IloRange(env_,
					      d_min_,
					      IloScalProd(columns_, pattern_deviations_),
					      IloInfinity));

    delta_ = IloAdd(master_problem_, IloRange(env_,
					      0,
					      IloScalProd(columns_, pattern_deviations_),
					      d_max_));


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


void Problem::SolveRelaxationIp() {
    while (true) {
	// Update the constraints by taking into account newly added columns
	master_problem_.remove(zeta_);
	zeta_ = IloAdd(master_problem_, IloRange(env_,
						 num_bins_,
						 IloSum(columns_),
						 num_bins_));
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

	master_solver_.solve();
	lower_bound_ = master_solver_.getObjValue();
	lower_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_,
								     pattern_deviations_));

	if ((std::chrono::duration<double>(std::chrono::high_resolution_clock::now()
					   - time_start_).count()) > time_limit_) break;

	IloModel sub_problem(env_);
	IloCplex sub_solver(sub_problem);
	IloObjective sub_objective = IloAdd(sub_problem, IloMaximize(env_));
	sub_solver.setOut(env_.getNullStream()); // Supress Cplex output
	IloNumVarArray z(env_, num_items_, 0, 1, ILOINT);
	
	// Auxiliary variable: Beta (deviation)
	IloNumVar beta(env_);
	sub_problem.add((IloScalProd(z, weights_)-mean_load_) <= beta);
	sub_problem.add(mean_load_-(IloScalProd(z, weights_)) <= beta);

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
	
	// Dual values
	IloNumArray y(env_, num_items_);
	master_solver_.getDuals(y, x_);
	IloNum zeta = master_solver_.getDual(zeta_);
	IloNum gamma = master_solver_.getDual(gamma_);
	IloNum delta = master_solver_.getDual(delta_);

	//IloNumExpr obj1(IloScalProd(z, y));
	IloNumExpr obj1(env_);
	obj1 = IloScalProd(z, y);

        //IloNumExpr obj2(env_, zeta);
	IloNumExpr obj2(env_);
	obj2 += zeta;

	IloNumExpr obj3(env_);
	if (norm_ == 1) {
	    obj3 = beta*(gamma+delta);
	}
	else if (norm_ == 2) {
	    obj3 = beta*beta*(gamma+delta);
	}
	else if (norm_ == 3) {
	    obj3 = beta*(gamma+delta);
	}
	
	IloNumExpr obj4(env_);
	for (int i=0; i<num_items_-1; i++) {
	    for (int j=i+1; j<num_items_; j++) {
		obj4 += z[i]*z[j]*costs_[i][j];
	    }
	}

	// Objectives
	if ((gamma+delta) <= 0) {
	    sub_objective.setExpr(obj1+obj2+obj3-obj4);
	}
	else if ((gamma+delta) > 0) {
	    sub_objective.setExpr(obj1+obj2-obj3-obj4);
	}

	sub_solver.solve();

	bool new_column_added = false;
	for (int i=0; i<sub_solver.getSolnPoolNsolns(); i++) {
	    if ((sub_solver.getValue(obj1, i)+
	    	 sub_solver.getValue(obj2, i)+
	    	 sub_solver.getValue(obj3, i)
	    	 >
	    	 sub_solver.getValue(obj4, i))) {

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
	deviation = IloAbs(mean_load_-weight);
    }
    else if (norm_ == 2) {
	deviation = IloAbs(mean_load_-weight);
	deviation = deviation*deviation;
    }
    else if (norm_ == 3) {
	deviation = IloAbs(mean_load_-weight);
    }
    return deviation;
}


void Problem::SolveIntegrality() {
    master_problem_.add(IloConversion(env_, columns_, ILOINT));

    master_problem_.remove(zeta_);
    zeta_ = IloAdd(master_problem_, IloRange(env_,
					     num_bins_,
					     IloSum(columns_),
					     num_bins_));
    
    master_problem_.remove(delta_);
    delta_ = IloAdd(master_problem_, IloRange(env_,
					      0,
					      IloScalProd(columns_, pattern_deviations_),
					      d_max_));
    
    master_problem_.remove(gamma_);
    gamma_ = IloAdd(master_problem_, IloRange(env_,
					      d_min_,
					      IloScalProd(columns_, pattern_deviations_),
					      IloInfinity));

    master_solver_.solve();

    upper_bound_ = master_solver_.getObjValue();
    upper_bound_deviation_ = master_solver_.getValue(IloScalProd(columns_, pattern_deviations_));
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
