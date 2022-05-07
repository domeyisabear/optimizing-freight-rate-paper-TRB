#pragma once
#include "Parameters.h"
#include <gurobi_c++.h>

class RelaxMasterProblem_exf
{
private:
	vector<double> rateMax;
	vector<double> rateMin;

	//GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	GRBVar1d r_od;
	GRBVar2d R_w_od, g_w_od, q_w_od_h;

	double scale_param = 1e4;
public:
	RelaxMasterProblem_exf(vector<double>& rate_min, vector<double>& rate_max);
	RelaxMasterSolution getOptimalValue();
	~RelaxMasterProblem_exf();
};

