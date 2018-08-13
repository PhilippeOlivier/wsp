#ifndef NODE_H_
#define NODE_H_


#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h> //remove?


class Node {
public:
  Node(Problem problem);
  Node(Problem problem,
       Node node,
       IloInt fix_at,
       IloInt column_)
  ~Node();

  IloInt GetDepth();
  IloBool IsIntegral();
  IloNum GetObjectiveValue();
  IloInt GetNumFractionalColumns();
    
private:
  IloInt depth_;
  IloInt num_fractional_cols_;
  IloArray<IloInt> columns_set_at_0_;
  IloArray<IloInt> columns_set_at_1_;
  IloNum objective_value_; // if num frac cols == 0 then this is integral
  // when checking the LR, if round(LR) (rounded towards 0) is equal to the best integral sol found yet, then prune
  // also, have a copy constructor in Problem, so that I can find an initial integral solution to start the branch and price with
};


#endif  // PROBLEM_H_
