#pragma once
#include "Parameters.h"
#include "RelaxSubProblem.h"

class RangeReduction
{
private:
	vector<double> rateMax;
	vector<double> rateMin;

	vector<RelaxSubProblem*> SPList;

	//GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar1d r_od_;
	GRBVar Z_;

	double out_UB;

	RelaxSubDualSolution evalSolution(RelaxMasterSolution& rms);

	double scale_param = 1e4;
public:
	RangeReduction(vector<double>& rate_min, vector<double>& rate_max, double out_ub);
	double getNewBound(int demandInd, bool isLB);
	~RangeReduction();
};

