#include "RelaxSubProblem.h"
#include "Utils.h"

RelaxSubProblem::RelaxSubProblem(int sind, vector<double>& rate_min, vector<double>& rate_max)
{
	scenarioIndex = sind;
	rateMin = rate_min;
	rateMax = rate_max;

	lineApx = getLineApprox(rateMin, rateMax);

	alpha_od = addVars(model, -GRB_INFINITY, 0, GRB_CONTINUOUS, "alpha", param->demandNum);
	beta_od_k= vector<vector<GRBVar>>(param->demandNum);
	lambda_od_k = vector<vector<GRBVar>>(param->demandNum);
	phi_i = addVars(model, -GRB_INFINITY, 0, GRB_CONTINUOUS, "phi", param->shippingLegNum);
	pi_od = addVars(model, -GRB_INFINITY, GRB_INFINITY, GRB_CONTINUOUS, "pi", param->demandNum);
	for(int i=0;i<param->demandNum;i++)
	{
		vector<GRBVar> tmp(lineApx[i].size());
		for (int k = 0; k < lineApx[i].size(); k++)
			tmp[k] = model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, "beta");
		beta_od_k[i] = tmp;

		vector<GRBVar> tmp2(lineApx[i].size());
		for (int k = 0; k < lineApx[i].size(); k++)
			tmp2[k] = model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, "lambda");
		lambda_od_k[i] = tmp2;
	}

	// constraints
	// constraint (30)
	for(int i=0;i<param->demandNum;i++)
	{
		GRBLinExpr expr = 0;
		for(int k=0;k<lineApx[i].size();k++)
		{
			expr += beta_od_k[i][k];
		}
		R_od.push_back(model.addConstr(alpha_od[i] + expr == -1));
	}
	// constraint (31)
	for(int i=0;i<param->demandNum;i++)
	{
		GRBLinExpr expr1 = 0;
		GRBLinExpr expr2 = 0;
		for (int k = 0; k < lineApx[i].size(); k++)
		{
			expr1 += beta_od_k[i][k];
			expr2 += lambda_od_k[i][k];
		}
		g_od.push_back(model.addConstr(-rateMax[i] * alpha_od[i] - rateMin[i] * expr1 + expr2 + pi_od[i] <= 0));
	}
	// constraint (32)
	for(int i=0;i<param->containerPathNum;i++)
	{
		GRBLinExpr expr = 0;
		for(int j=0;j<param->containerPathList[i].shippingLegIndex.size();j++)
		{
			expr += phi_i[param->containerPathList[i].shippingLegIndex[j]];
		}
		q_od_h.push_back(model.addConstr(expr - pi_od[param->containerPathList[i].demandIndex] <= param->containerPathList[i].handleCost));
	}

	GRBLinExpr obj = 0;
	model.setObjective(obj, GRB_MAXIMIZE);

	model.set(GRB_IntParam_Threads, threadNumPerSolver);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_Method, 0);
	model.update();
}


void RelaxSubProblem::optimize(RelaxMasterSolution& rms)
{
	GRBLinExpr expr = 0;
	for (int i = 0; i < param->demandNum; i++)
	{
		double ratio = param->demandList[i].ratioScenario[scenarioIndex];
		for (int k = 0; k < lineApx[i].size(); k++)
		{
			double a = lineApx[i][k].a;
			double b = lineApx[i][k].b;
			expr += ratio * (a*rms.r_od[i] + b)*(rms.r_od[i] - rateMin[i])*beta_od_k[i][k] + ratio * (a*rms.r_od[i] + b)*lambda_od_k[i][k];
		}
	}
	for (int i = 0; i < param->shippingLegNum; i++)
	{
		expr += param->shippingLegList[i].capacityScenario[scenarioIndex] * phi_i[i];
	}
	model.setObjective(expr, GRB_MAXIMIZE);
	model.optimize();
	//cout << model.get(GRB_IntAttr_Status) << endl;
}

RelaxSubDualSolution RelaxSubProblem::getDualSolution(RelaxMasterSolution& rms)
{
	optimize(rms);

	RelaxSubDualSolution rsds;

	//if(model.get(GRB_IntAttr_Status)!=2)
	//{
	//	rsds.objVal = -MAX_NUM;
	//	return rsds;
	//}

	vector<vector<double>> beta_val(param->demandNum);
	vector<vector<double>> lambda_val(param->demandNum);
	vector<double> phi_val = getOptimalSolution(phi_i);
	for(int i=0;i<param->demandNum;i++)
	{
		beta_val[i] = getOptimalSolution(beta_od_k[i]);
		lambda_val[i] = getOptimalSolution(lambda_od_k[i]);
	}

	rsds.eta_0 = 0;
	rsds.eta_od_1 = vector<double>(param->demandNum);
	rsds.eta_od_2 = vector<double>(param->demandNum);
	for (int i = 0; i < param->shippingLegNum; i++)
	{
		rsds.eta_0 += param->shippingLegList[i].capacityScenario[scenarioIndex] * phi_val[i];
	}
	for(int i=0;i<param->demandNum;i++)
	{
		double ratio = param->demandList[i].ratioScenario[scenarioIndex];
		for(int k=0;k<lineApx[i].size();k++)
		{
			double a = lineApx[i][k].a;
			double b = lineApx[i][k].b;

			rsds.eta_0 += ratio * b*(lambda_val[i][k] - rateMin[i] * beta_val[i][k]);
			rsds.eta_od_1[i] += ratio * (b*beta_val[i][k] - a * rateMin[i] * beta_val[i][k] + a * lambda_val[i][k]);
			rsds.eta_od_2[i] += ratio * a*beta_val[i][k];
		}
	}

	rsds.objVal= model.get(GRB_DoubleAttr_ObjVal);

	return rsds;
}

SubSolution RelaxSubProblem::getPrimalSolution(RelaxMasterSolution& rms)
{
	optimize(rms);
	SubSolution ss;

	//if (model.get(GRB_IntAttr_Status) != 2)
	//{
	//	ss.objVal = -MAX_NUM;
	//	return ss;
	//}

	ss.g_od = getOptimalDualSolution(g_od);
	ss.q_od_h = getOptimalDualSolution(q_od_h);
	ss.R_od = getOptimalDualSolution(R_od);
	ss.objVal= model.get(GRB_DoubleAttr_ObjVal);

	return ss;
}


RelaxSubProblem::~RelaxSubProblem()
{
	
}
