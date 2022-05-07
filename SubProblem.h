#pragma once
#include "Parameters.h"

class SubProblem
{
private:
	//GRBEnv env=GRBEnv();
	GRBModel model=GRBModel(env);

	GRBVar1d g_od;
	GRBVar1d q_od_h;
	vector<GRBConstr> lambda,phi,pi;

	int scenarioInd;

	void optimize(vector<double>& rateList, vector<double>& demandList);

public:
	SubSolution getOptimalValue(vector<double>& rateList, vector<double>& demandList);
	SubDualSolution getOptimalDualValue(vector<double>& rateList, vector<double>& demandList);
	SubProblem(int sind);
	~SubProblem();
};

