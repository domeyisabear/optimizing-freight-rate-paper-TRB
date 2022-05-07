#include "Utils.h"
#include <fstream>
#include <iomanip>

GRBVar1d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d)
{
	vector<GRBVar> vars(d);
	for (int i = 0; i < d; i++)
	{
		string str = v_name + "_" + to_string(i);
		vars[i] = model.addVar(lb, ub, 0, vtype, str);
	}
	return vars;
}

GRBVar2d addVars(GRBModel& model, double lb, double ub, char vtype, std::string v_name, int d1, int d2)
{
	GRBVar2d vars(d1);
	for (int i = 0; i < d1; i++)
	{
		vars[i] = GRBVar1d(d2);
		for (int j = 0; j < d2; j++)
		{
			string str = v_name + "_" + to_string(i) + "_" + to_string(j);
			vars[i][j] = model.addVar(lb, ub, 0, vtype, str);
		}
	}
	return vars;
}

vector<double> getOptimalSolution(GRBVar1d& var)
{
	int len = var.size();
	vector<double> val = vector<double>(len);
	for (int i = 0; i < len; i++)
	{
		val[i] = var[i].get(GRB_DoubleAttr_X);
	}

	return val;
}

vector<vector<double>> getOptimalSolution(GRBVar2d& var)
{
	int len = var.size();
	vector<vector<double>> val = vector<vector<double>>(len);
	for (int i = 0; i < len; i++)
	{
		val[i] = vector<double>(var[i].size());
		for (int j = 0; j < var[i].size(); j++)
		{
			val[i][j] = var[i][j].get(GRB_DoubleAttr_X);
		}
	}

	return val;
}

vector<double> getOptimalDualSolution(vector<GRBConstr>& constrs)
{
	int len = constrs.size();
	vector<double> re = vector<double>(len);
	for (int i = 0; i < len; i++)
		re[i] = constrs[i].get(GRB_DoubleAttr_Pi);
	return re;
}

///////////////////////////////////////////////////////////
////// demand function: D(r)=param1-param2*r
vector<double> getDemand(vector<double>& rateList)
{
	vector<double> demandList(rateList.size());
	for (int i = 0; i < rateList.size(); i++)
	{
		double param1 = param->demandList[i].demandParam[0];
		double param2 = param->demandList[i].demandParam[1];
		double val = param1 - param2 * rateList[i];
		demandList[i] = val > 0 ? val : 0;
	}
	return demandList;
}

vector<double> getDemandDerivative(vector<double>& rateList)
{
	vector<double> result(rateList.size());
	for (int i = 0; i < rateList.size(); i++)
	{
		double param1 = param->demandList[i].demandParam[0];
		double param2 = param->demandList[i].demandParam[1];
		double val = param1 - param2 * rateList[i];
		result[i] = val > 0 ? -param2 : 0;
	}
	return result;
}

vector<vector<DemandLineApproxParam>> getLineApprox(vector<double>& rateMin, vector<double>& rateMax)
{
	vector<vector<DemandLineApproxParam>> result;
	vector<double> demandMax = getDemand(rateMin);
	vector<double> demandMin = getDemand(rateMax);

	for(int i=0;i<param->demandNum;i++)
	{
		vector<DemandLineApproxParam> item;
		DemandLineApproxParam tmp;
		tmp.a = (demandMin[i] - demandMax[i]) / (rateMax[i] - rateMin[i]);
		tmp.b = demandMin[i] - tmp.a*rateMax[i];
		item.push_back(tmp);
		result.push_back(item);
	}
	return result;
}

///////////////////////////////////////////////////////////////

void saveResult(vector<double> result, string fpath)
{
	ofstream outputfile(fpath);

	for(int i=0;i<result.size();i++)
	{
		outputfile << setprecision(16) << result[i] << "\n";
	}

	outputfile.close();
}


