#pragma once
#include "Parameters.h"
#include "SubProblem.h"
#include "Utils.h"

class MasterProblem
{
private:
	vector<SubProblem*> v_sp=vector<SubProblem*>(param->scenarioNum);
public:
	double getOptimalValue(vector<double>& rateList)
	{
		vector<double> demandList = getDemand(rateList);

		vector<future<double>> results;
		for (int i = 0; i < param->scenarioNum; i++)
		{
			results.push_back(param->executor->commit(
				[](SubProblem* subp, vector<double>& rate, vector<double>& demand)
				{
					return subp->getOptimalValue(rate, demand).objVal;
				}, 
				v_sp[i], 
				rateList, 
				demandList));
		}
		double val = 0;
		for(int i=0;i<param->scenarioNum;i++)
		{
			double ss = results[i].get();
			val += ss*param->scenarioProbList[i];
		}
		return val;
	}

	vector<double> getDerivative(vector<double>& rateList)
	{
		vector<double> demandList = getDemand(rateList);
		vector<future<SubDualSolution>> futures;
		for (int i = 0; i < param->scenarioNum; i++)
		{
			futures.push_back(param->executor->commit(
				[](SubProblem* subp, vector<double>& rate, vector<double>& demand)
			{
				return subp->getOptimalDualValue(rate, demand);
			},
				v_sp[i],
				rateList,
				demandList));
		}

		vector<double> demand_der = getDemandDerivative(rateList);
		vector<double> result(param->demandNum);
		for (int w = 0; w < param->scenarioNum; w++)
		{
			SubDualSolution sds = futures[w].get();
			for(int i=0;i<param->demandNum;i++)
			{
				result[i] += param->scenarioProbList[w] * param->demandList[i].ratioScenario[w] * sds.lambda[i] * demand_der[i];
			}
		}

		return result;
	}

	MasterProblem()
	{
		vector<future<void>> future;
		for (int i = 0; i < param->scenarioNum; i++)
		{
			future.push_back(
				param->executor->commit(
				[](int w, vector<SubProblem*>* v)
				{
					SubProblem* subp = new SubProblem(w);
					(*v)[w]=subp;
				},
			i,
			&v_sp
			)
			);
		}
		for (int i = 0; i < param->scenarioNum; i++)
			future[i].get();
	}
	~MasterProblem()
	{
		for (int i = 0; i < param->scenarioNum; i++)
		{
			delete v_sp[i];
		}
	}
};
