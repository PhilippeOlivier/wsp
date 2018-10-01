#include "node.h"


Node::Node() {
  env_ = IloEnv();
  depth_ = 0;
}


/******************************************************************************/


Node::Node(Node& node, IloInt side) {
  env_ = IloEnv();
  depth_ = node.GetDepth()+1;
  columns_set_at_0_ = node.GetColumnsSetAt0();
  columns_set_at_1_ = node.GetColumnsSetAt1();
  if (side == 0) {
    columns_set_at_0_.push_back(node.split_on_column_);
  }
  else if (side == 1) {
    columns_set_at_1_.push_back(node.split_on_column_);
  }
}


/******************************************************************************/


Node::~Node() {
  env_.end();
}


/******************************************************************************/


void Node::Solve(Problem& problem) {
  problem.SolveWithPartiallySetColumns(columns_set_at_0_, columns_set_at_1_);
  is_feasible_ = problem.IsFeasible();
  if (is_feasible_) {
    IloInt next_column_to_split = problem.GetFractionalColumn();
    if (next_column_to_split == -1) { // Solution is integral
      is_integral_ = IloTrue;
    }
    else {
      is_integral_ = IloFalse;
      split_on_column_ = next_column_to_split;
    }
    objective_value_ = problem.GetLowerBound();
    deviation_ = problem.GetLowerBoundDeviation();
  }
}


/******************************************************************************/


void Node::Print() {
  cout << "Depth: " << GetDepth() << endl;
  cout << "Feasible: " << IsFeasible() << endl;
  cout << "Integral: " << IsIntegral() << endl;
  cout << "Objective: " << GetObjectiveValue() << endl;
  cout << "Deviation: " << GetDeviation() << endl;
  cout << "Split on: " << GetSplitOnColumn() << endl;
  cout << "Columns set at 0: ";
  for (IloInt i=0; i<(IloInt)columns_set_at_0_.size(); i++) {
    cout << columns_set_at_0_[i] << " ";
  }
  cout << endl;
  cout << "Columns set at 1: ";
  for (IloInt i=0; i<(IloInt)columns_set_at_1_.size(); i++) {
    cout << columns_set_at_1_[i] << " ";
  }
  cout << endl;
}

  
/******************************************************************************/


IloInt Node::GetDepth() {
  return depth_;
}


/******************************************************************************/


std::vector<IloInt> Node::GetColumnsSetAt0() {
  return columns_set_at_0_;
}


/******************************************************************************/


std::vector<IloInt> Node::GetColumnsSetAt1() {
  return columns_set_at_1_;
}


/******************************************************************************/


IloBool Node::IsFeasible() {
  return is_feasible_;
}


/******************************************************************************/


IloBool Node::IsIntegral() {
  return is_integral_;
}


/******************************************************************************/


IloInt Node::GetObjectiveValue() {
  return objective_value_;
}


/******************************************************************************/


IloNum Node::GetDeviation() {
  return deviation_;
}


/******************************************************************************/


IloInt Node::GetSplitOnColumn() {
  return split_on_column_;
}
