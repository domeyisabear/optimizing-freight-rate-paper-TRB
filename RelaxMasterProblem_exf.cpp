#include "RelaxMasterProblem_exf.h"
#include "Utils.h"
#include <gurobi_c++.h>


RelaxMasterProblem_exf::RelaxMasterProblem_exf(vector<double>& rate_min, vector<double>& rate_max)
{
	rateMin = rate_min;
	rateMax = rate_max;

	r_od = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "r_od", param->demandNum);
	R_w_od = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "R", param->scenarioNum, param->demandNum);
	g_w_od = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "g", param->scenarioNum, param->demandNum);
	q_w_od_h = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "q", param->scenarioNum, param->containerPathNum);

	vector<vector<DemandLineApproxParam>> lineApx = getLineApprox(rateMin, rateMax);
	for (int w = 0; w < param->scenarioNum; w++)
	{
		for (int i = 0; i < param->demandNum; i++)
		{
			//cout << i << endl;
			// constraint (8)
			model.addConstr(R_w_od[w][i] - rateMax[i] / scale_param * g_w_od[w][i] <= 0);
			// constraints (9)(10) or (12)(13)
			double ratio = param->demandList[i].ratioScenario[w];
			for (int k = 0; k < lineApx[i].size(); k++)
			{
				double a = lineApx[i][k].a;
				double b = lineApx[i][k].b;
				GRBQuadExpr expr = R_w_od[w][i] - rateMin[i] / scale_param * g_w_od[w][i] - ratio * (a*r_od[i] + b / scale_param)*(r_od[i] - rateMin[i] / scale_param)*scale_param;

				model.addQConstr(expr <= 0);
				model.addConstr(g_w_od[w][i] <= ratio * (a*r_od[i] * scale_param + b));
			}
			// constraint (4)
			GRBLinExpr expr = 0;
			vector<int> pathInd = param->demandList[i].containerPathIndex;
			for (auto ind : pathInd)
			{
				expr += q_w_od_h[w][ind];
			}
			model.addConstr(g_w_od[w][i] - expr == 0);
		}
		// constraint (3)
		for (int i = 0; i < param->shippingLegNum; i++)
		{
			vector<int> pathInd = param->shippingLegList[i].pathIndexList;
			GRBLinExpr expr = 0;
			for (auto ind : pathInd)
			{
				expr += q_w_od_h[w][ind];
			}
			double cap = param->shippingLegList[i].capacityScenario[w];
			model.addConstr(expr <= cap);
		}
	}
	// constraint (5)
	for (int i = 0; i < param->demandNum; i++)
	{
		model.addConstr(r_od[i] <= rateMax[i] / scale_param);
		model.addConstr(r_od[i] >= rateMin[i] / scale_param);
	}

	GRBLinExpr obj = 0;
	for (int w = 0; w < param->scenarioNum; w++)
	{
		GRBLinExpr cost = 0;
		for (int i = 0; i < param->containerPathNum; i++)
		{
			cost += q_w_od_h[w][i] * param->containerPathList[i].handleCost;
		}
		GRBLinExpr revenue = 0;
		for (int i = 0; i < param->demandNum; i++)
		{
			revenue += scale_param * R_w_od[w][i];
		}
		obj += param->scenarioProbList[w] * (cost - revenue);
	}
	model.setObjective(obj, GRB_MINIMIZE);
	//model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_NumericFocus, 3);
	model.set(GRB_IntParam_ScaleFlag, 3);
	model.set(GRB_DoubleParam_ObjScale, 1000);
	model.set(GRB_IntParam_BarHomogeneous, 1);
}

RelaxMasterSolution RelaxMasterProblem_exf::getOptimalValue()
{
	model.optimize();
	//cout << model.get(GRB_IntAttr_Status) << endl;
	double objVal = model.get(GRB_DoubleAttr_ObjVal);
	RelaxMasterSolution rms;
	rms.r_od = getOptimalSolution(r_od);
	for (int i = 0; i < param->demandNum; i++)
		rms.r_od[i] *= scale_param;
	rms.objVal = objVal;
	rms.LB = objVal;
	return rms;
}

RelaxMasterProblem_exf::~RelaxMasterProblem_exf()
{
}
