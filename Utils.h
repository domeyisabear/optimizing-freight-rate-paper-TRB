#pragma once
#include <vector>
#include "gurobi_c++.h"
#include "Parameters.h"

GRBVar1d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d);
GRBVar2d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d1, int d2);

vector<double> getOptimalSolution(GRBVar1d& var);
vector<vector<double>> getOptimalSolution(GRBVar2d& var);
vector<double> getOptimalDualSolution(vector<GRBConstr>& constrs);

vector<double> getDemand(vector<double>& rateList);
vector<double> getDemandDerivative(vector<double>& rateList);

vector<vector<DemandLineApproxParam>> getLineApprox(vector<double>& rateMin, vector<double>& rateMax);

void saveResult(vector<double> result, string fpath);