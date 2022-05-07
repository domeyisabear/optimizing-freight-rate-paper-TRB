#include "SubProblem.h"
#include "Utils.h"

SubProblem::SubProblem(int sind)
{
	scenarioInd = sind;

	// variables
	g_od = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "g_od", param->demandNum);
	q_od_h = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "q_od_h", param->containerPathNum);

	// constraint (2)
	for (int i = 0; i < param->demandNum; i++)
	{
		lambda.push_back(model.addConstr(g_od[i] <= 0));
	}

	// constraint (4)
	for (int i = 0; i < param->demandNum; i++)
	{
		GRBLinExpr expr = 0;
		vector<int> pathInd = param->demandList[i].containerPathIndex;
		for (auto ind : pathInd)
		{
			expr += q_od_h[ind];
		}
		pi.push_back(model.addConstr(g_od[i] - expr == 0));
	}

	// constraint (3)
	for (int i = 0; i < param->shippingLegNum; i++)
	{
		vector<int> pathInd = param->shippingLegList[i].pathIndexList;
		GRBLinExpr expr = 0;
		for (auto ind : pathInd)
		{
			expr += q_od_h[ind];
		}
		double cap = param->shippingLegList[i].capacityScenario[sind];
		phi.push_back(model.addConstr(expr <= cap));
	}

	GRBLinExpr obj = 0;
	model.setObjective(obj, GRB_MINIMIZE);

	model.set(GRB_IntParam_Method, 1);
	model.set(GRB_IntParam_Threads, threadNumPerSolver);
	model.set(GRB_IntParam_OutputFlag, 0);
	model.update();
}

void SubProblem::optimize(vector<double>& rateList, vector<double>& demandList)
{
	// constraint (2)
	for (int i = 0; i < param->demandNum; i++)
	{
		double demandSize = param->demandList[i].ratioScenario[scenarioInd] * demandList[i];
		lambda[i].set(GRB_DoubleAttr_RHS, demandSize);
	}

	// obj
	GRBLinExpr cost = 0;
	for (int i = 0; i < param->containerPathNum; i++)
	{
		cost += q_od_h[i] * param->containerPathList[i].handleCost;
	}
	GRBLinExpr revenue = 0;
	for (int i = 0; i < param->demandNum; i++)
	{
		revenue += g_od[i] * rateList[i];
	}

	model.setObjective(cost - revenue, GRB_MINIMIZE);

	model.optimize();
}

SubSolution SubProblem::getOptimalValue(vector<double>& rateList, vector<double>& demandList)
{
	optimize(rateList, demandList);

	SubSolution ss;
	ss.g_od = getOptimalSolution(g_od);
	ss.q_od_h = getOptimalSolution(q_od_h);
	ss.objVal = model.get(GRB_DoubleAttr_ObjVal);

	return ss;
}

SubDualSolution SubProblem::getOptimalDualValue(vector<double>& rateList, vector<double>& demandList)
{
	optimize(rateList, demandList);

	SubDualSolution sds;
	sds.phi = getOptimalDualSolution(phi);
	sds.lambda = getOptimalDualSolution(lambda);
	sds.pi = getOptimalDualSolution(pi);
	sds.objVal= model.get(GRB_DoubleAttr_ObjVal);

	return sds;	
}

SubProblem::~SubProblem()
{

}
