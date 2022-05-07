#include "RelaxMasterProblem.h"
#include <iostream>

RelaxMasterProblem::RelaxMasterProblem(vector<double>& rate_min, vector<double>& rate_max):
rateMax(rate_max),rateMin(rate_min)
{
	// initialize relax subproblem
	SPList = vector<RelaxSubProblem*>(param->scenarioNum);

	vector<future<void>> future;
	for (int i = 0; i < param->scenarioNum; i++)
	{
		future.push_back(
			param->executor->commit(
				[](int w, vector<RelaxSubProblem*>* u, vector<double>& rateMin, vector<double>& rateMax)
		{
			(*u)[w] = new RelaxSubProblem(w, rateMin, rateMax);
		},
			i,
			&SPList,
			rateMin,
			rateMax
			)
		);
	}
	for (int i = 0; i < param->scenarioNum; i++)
		future[i].get();

	Z_ = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "Z");
	r_od_ = addVars(model, 0, GRB_INFINITY, GRB_CONTINUOUS, "r_od", param->demandNum);

	// constraints
	for (int i = 0; i < param->demandNum; i++)
	{
		model.addConstr(scale_param * r_od_[i] <= rateMax[i]);
		model.addConstr(scale_param * r_od_[i] >= rateMin[i]);
	}
	RelaxMasterSolution rms;
	rms.r_od = rateMin;
	RelaxSubDualSolution rsds = evalSolution(rms);
	GRBQuadExpr expr = rsds.eta_0 / scale_param;
	for(int i=0;i<param->demandNum;i++)
	{
		expr += rsds.eta_od_1[i] * r_od_[i] + rsds.eta_od_2[i] * scale_param * r_od_[i] * r_od_[i];
	}
	model.addQConstr(Z_ >= expr);

	model.update();

	model.set(GRB_IntParam_OutputFlag, 0);
	model.set(GRB_IntParam_NumericFocus, 3);
	model.set(GRB_IntParam_ScaleFlag, 3);
	model.set(GRB_DoubleParam_ObjScale, 1000);
	model.set(GRB_IntParam_BarHomogeneous, 1);
}

RelaxSubDualSolution RelaxMasterProblem::evalSolution(RelaxMasterSolution& rms)
{
	vector<future<RelaxSubDualSolution>> futures;
	for (int j = 0; j < param->scenarioNum; j++)
	{
		futures.push_back(param->executor->commit(
			[](RelaxSubProblem* sp, RelaxMasterSolution& rms)
		{
			return sp->getDualSolution(rms);
		},
			SPList[j],
			rms
			));
	}
	
	RelaxSubDualSolution rsds_result;
	rsds_result.objVal = 0;
	rsds_result.eta_0 = 0;
	rsds_result.eta_od_1 = vector<double>(param->demandNum);
	rsds_result.eta_od_2 = vector<double>(param->demandNum);
	for (int j = 0; j < param->scenarioNum; j++)
	{
		RelaxSubDualSolution rsds = futures[j].get();
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

RelaxMasterSolution RelaxMasterProblem::getOptimalValue()
{
	double LB = -MAX_NUM;
	double UB = MAX_NUM;
	double eps = (UB - LB) / fabs(LB + 1e-8);

	// obj
	GRBLinExpr obj = Z_ * scale_param;
	model.setObjective(obj, GRB_MINIMIZE);

	RelaxMasterSolution incum_rms;

	int iteration = 0;

	while (eps > in_eps && iteration<=30)
	{
		model.optimize();
		LB = model.get(GRB_DoubleAttr_ObjVal);

		//cout << model.get(GRB_IntAttr_Status) << endl;

		RelaxMasterSolution rms;
		rms.r_od = getOptimalSolution(r_od_);
		for (int i = 0; i < param->demandNum; i++)
			rms.r_od[i] *= scale_param;
		RelaxSubDualSolution rsds = evalSolution(rms);

		//add cut
		GRBQuadExpr expr = rsds.eta_0/scale_param;
		for (int i = 0; i < param->demandNum; i++)
		{
			expr += rsds.eta_od_1[i] * r_od_[i] + rsds.eta_od_2[i] * scale_param * r_od_[i] * r_od_[i];
		}
		model.addQConstr(Z_ >= expr);

		if (UB > rsds.objVal)
		{
			UB = rsds.objVal;
			incum_rms = rms;
		}

		// improve primal solution ///////////////////////////////////////////
		bool flag = false;
		{
			int primal_iter_num = 30;
			double primal_coeff = 0.5;
			RelaxMasterSolution rms = incum_rms;
			double primal_val = UB;
			int no_improve_iter = 0;

			for (int n = 0; n < primal_iter_num; n++)
			{
				// step size
				double sum = 0;
				for (int i = 0; i < param->demandNum; i++)
					sum += (rsds.eta_od_1[i] + 2 * rsds.eta_od_2[i] * rms.r_od[i])*(rsds.eta_od_1[i] + 2 * rsds.eta_od_2[i] * rms.r_od[i]);
				double step_size = primal_coeff * (primal_val - LB) / sum;
				//cout << step_size << endl;

				// update primal solution
				for (int i = 0; i < param->demandNum; i++)
				{
					rms.r_od[i] += -(rsds.eta_od_1[i] + 2 * rsds.eta_od_2[i] * rms.r_od[i])*step_size;
					if (rms.r_od[i] < rateMin[i])
						rms.r_od[i] = rateMin[i];
					else if (rms.r_od[i] > rateMax[i])
						rms.r_od[i] = rateMax[i];
				}

				// add cut
				rsds = evalSolution(rms);
				GRBQuadExpr expr = rsds.eta_0/scale_param;
				for (int i = 0; i < param->demandNum; i++)
				{
					expr += rsds.eta_od_1[i] * r_od_[i] + rsds.eta_od_2[i] * scale_param * r_od_[i] * r_od_[i];
				}
				model.addQConstr(Z_ >= expr);

				//cout << rsds.objVal << endl;
				primal_val = rsds.objVal;
				if (rsds.objVal < UB)
				{
					UB = rsds.objVal;
					incum_rms = rms;
					flag = true;
				}
				else
					no_improve_iter++;

				if (no_improve_iter >= 5)
				{
					primal_coeff *= 0.8;
					no_improve_iter = 0;
				}
			}
		}

		iteration++;
		eps = (UB - LB) / fabs(LB + 1e-8);
		//std::cout << "UB: " << UB << ", LB: " << LB << ", eps: " << eps * 100 << "%" << endl;
	}
	incum_rms.objVal = UB;
	incum_rms.LB = UB;
	return incum_rms;
}

RelaxMasterProblem::~RelaxMasterProblem()
{
	for (int i = 0; i < param->scenarioNum; i++)
	{
		delete SPList[i];
	}
}
