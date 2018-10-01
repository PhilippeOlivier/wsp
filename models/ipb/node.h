#ifndef NODE_H_
#define NODE_H_


#include <ilcplex/ilocplex.h>
#include <vector>
#include "problem.h"


class Node {
public:
  Node();
  Node(Node& node, IloInt side);
  ~Node();

  void Solve(Problem& problem);
  void Print();
  
  IloInt GetDepth();
  std::vector<IloInt> GetColumnsSetAt0();
  std::vector<IloInt> GetColumnsSetAt1();

  IloBool IsFeasible();
  IloBool IsIntegral();
  IloInt GetObjectiveValue();
  IloNum GetDeviation();
  IloInt GetSplitOnColumn();

private:
  IloEnv env_;
  IloInt depth_;
  std::vector<IloInt> columns_set_at_0_;
  std::vector<IloInt> columns_set_at_1_;
  
  IloBool is_feasible_; // Is the node feasible?
  IloBool is_integral_; // Is the solution integral?
  IloInt objective_value_;
  IloNum deviation_;
  IloInt split_on_column_; // Split on this column next
};


#endif  // PROBLEM_H_
