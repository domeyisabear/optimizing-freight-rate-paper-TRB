#pragma once
#include "Parameters.h"
class RelaxSubProblem
{
private:
	int scenarioIndex;
	vector<double> rateMax;
	vector<double> rateMin;

	//GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	// dual vars
	GRBVar1d alpha_od,phi_i,pi_od;
	GRBVar2d beta_od_k, lambda_od_k;

	vector<vector<DemandLineApproxParam>> lineApx;

	vector<GRBConstr> R_od,g_od,q_od_h;

	void optimize(RelaxMasterSolution& rms);

public:
	RelaxSubProblem(int sind,vector<double>& rate_min,vector<double>& rate_max);
	RelaxSubDualSolution getDualSolution(RelaxMasterSolution& rms);
	SubSolution getPrimalSolution(RelaxMasterSolution& rms);
	~RelaxSubProblem();
};

