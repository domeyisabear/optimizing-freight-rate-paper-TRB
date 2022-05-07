#include "RangeReduction.h"
#include "Utils.h"

RelaxSubDualSolution RangeReduction::evalSolution(RelaxMasterSolution& rms)
{
	RelaxSubDualSolution rsds_result;
	rsds_result.objVal = 0;
	rsds_result.eta_0 = 0;
	rsds_result.eta_od_1 = vector<double>(param->demandNum);
	rsds_result.eta_od_2 = vector<double>(param->demandNum);
	for (int j = 0; j < param->scenarioNum; j++)
	{
		RelaxSubDualSolution rsds = SPList[j]->getDualSolution(rms);
		rsds_result.objVal += param->scenarioProbList[j] * rsds.objVal;
		rsds_result.eta_0 += param->scenarioProbList[j] * rsds.eta_0;
		for (int i = 0; i < param->demandNum; i++)
		{
			rsds_result.eta_od_1[i] += param->scenarioProbList[j] * rsds.eta_od_1[i];
			rsds_result.eta_od_2[i] += param->scenarioProbList[j] * rsds.eta_od_2[i];
		}
	}
	return rsds_result;
}

RangeReduction::RangeReduction(vector<double>& rate_min, vector<double>& rate_max, double out_ub) :
	rateMin(rate_min),
	rateMax(rate_max),
	out_UB(out_ub)
{
	// initialize relax subproblem
	SPList = vector<RelaxSubProblem*>(param->scenarioNum);
	for (int i = 0; i < param->scenarioNum; i++)
	{
		SPList[i] = new RelaxSubProblem(i, rateMin, rateMax);
	}

	//vector<future<void>> future;
	//for (int i = 0; i < param->scenarioNum; i++)
	//{
	//	future.push_back(
	//		param->executor->commit(
	//			[](int w, vector<RelaxSubProblem*>* u, vector<double>& rateMin, vector<double>& rateMax)
	//	{
	//		(*u)[w] = new RelaxSubProblem(w, rateMin, rateMax);
	//	},
	//			i,
	//		&SPList,
	//		rateMin,
	//		rateMax
	//		)
	//	);
	//}
	//for (int i = 0; i < param->scenarioNum; i++)
	//	future[i].get();

	Z_ = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "Z");
	r_od_ = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "r_od", param->demandNum);

	// constraints
	for (int i = 0; i < param->demandNum; i++)
	{
		model.addConstr(scale_param*r_od_[i] <= rateMax[i]);
		model.addConstr(scale_param*r_od_[i] >= rateMin[i]);
	}
	model.addConstr(scale_param*Z_ <= out_UB);

	model.update();
	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_NumericFocus, 3);
	model.set(GRB_IntParam_ScaleFlag, 3);
	model.set(GRB_DoubleParam_ObjScale, 1000);
	model.set(GRB_IntParam_BarHomogeneous, 1);
	model.set(GRB_IntParam_Threads, threadNumPerSolver);
}

double RangeReduction::getNewBound(int demandInd, bool isLB)
{
	if (isLB)
	{
		GRBLinExpr obj = scale_param*r_od_[demandInd];
		model.setObjective(obj, GRB_MINIMIZE);
	}
	else
	{
		GRBLinExpr obj = scale_param*r_od_[demandInd];
		model.setObjective(obj, GRB_MAXIMIZE);
	}

	double objVal;
	try {
		model.optimize();
		vector<double> val_r_od = getOptimalSolution(r_od_);
		for (int i = 0; i < param->demandNum; i++)
			val_r_od[i] *= scale_param;

		RelaxMasterSolution rms;
		rms.r_od = val_r_od;
		RelaxSubDualSolution rsds = evalSolution(rms);

		int iteration = 1;

		while (rsds.objVal > out_UB + 1 && iteration<=30)
		{
			//add cut
			GRBQuadExpr expr = rsds.eta_0 / scale_param;
			for (int i = 0; i < param->demandNum; i++)
			{
				expr += rsds.eta_od_1[i] * r_od_[i] + rsds.eta_od_2[i] * scale_param * r_od_[i] * r_od_[i];
			}
			model.addQConstr(Z_ >= expr);

			model.optimize();
			//cout << model.get(GRB_IntAttr_Status) << endl;

			val_r_od = getOptimalSolution(r_od_);
			for (int i = 0; i < param->demandNum; i++)
				val_r_od[i] *= scale_param;

			rms.r_od = val_r_od;
			rsds = evalSolution(rms);

			iteration++;
		}
		objVal = model.get(GRB_DoubleAttr_ObjVal);
	}
	catch(...)
	{
		objVal = isLB ? rateMin[demandInd] : rateMax[demandInd];
	}
	 
	if (isLB)
	{
		if (objVal > rateMin[demandInd] + 1)
		{
			cout << "LB improve! index: " << demandInd << " from " << rateMin[demandInd] << " to " << objVal << endl;
			return objVal;
		}
		else
			return rateMin[demandInd];
	}
	else
	{
		if (objVal < rateMax[demandInd] - 1)
		{
			cout << "UB improve! index: " << demandInd << " from " << rateMax[demandInd] << " to " << objVal << endl;
			return objVal;
		}
		else
			return rateMax[demandInd];
	}
}

RangeReduction::~RangeReduction()
{
	for (int i = 0; i < param->scenarioNum; i++)
		delete SPList[i];
}
