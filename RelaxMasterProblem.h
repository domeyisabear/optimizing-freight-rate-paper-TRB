#pragma once
#include "Parameters.h"
#include "Utils.h"
#include "SubProblem.h"
#include "RelaxSubProblem.h"

class RelaxMasterProblem
{
private:
	vector<double> rateMax;
	vector<double> rateMin;

	vector<RelaxSubProblem*> SPList;

	GRBModel model = GRBModel(env);
	GRBVar Z_;  // scaled variable for Z
	GRBVar1d r_od_; // scaled variable for r_od

	RelaxSubDualSolution evalSolution(RelaxMasterSolution& rms);

	double scale_param = 1e4; // scale parameter used to avoid numerical issues in QCP
public:
	RelaxMasterProblem(vector<double>& rate_min, vector<double>& rate_max);
	RelaxMasterSolution getOptimalValue();
	~RelaxMasterProblem();
};

