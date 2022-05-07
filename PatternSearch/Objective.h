#pragma once
#include <vector>
#include <malloc.h>
#include "../Parameters.h"
#include "../MasterProblem.h"

extern double tolerance;
extern bool newInitialPoint;
using namespace std;

class Objective
{
public:
	int nvars = 0;
	std::vector<double>  initpnt;
	MasterProblem* mp;
	std::vector<double>  LB;
	std::vector<double>  UB;
	void initpt(double *&x)
	{
		double grace = 1;
		x = new double[nvars];
		for (int i = 0; i < nvars; i++)
		{
			x[i] = initpnt[i];
		}
	}
	void fcn(double * x, double & f, int & flag)
	{
		vector<double> rate(param->demandNum);
		for (int i = 0; i < nvars; i++)
			rate[i] = x[i];
		flag = 1;
		f = mp->getOptimalValue(rate);
		//cout << f << endl;
		return;
	}

	void low_bound(double *&x)
	{
		x = new double[nvars];
		for (int i = 0; i < nvars; i++)
		{
			x[i] = LB[i];
		}
	}

	void up_bound(double *&x)
	{
		x = new double[nvars];
		for (int i = 0; i < nvars; i++)
		{
			x[i] = UB[i];
		}
	}

	Objective() {};

	Objective(vector<double>& pnt, MasterProblem& t_mp, vector<double>& lower_bound, vector<double>& upper_bound) :initpnt(pnt), LB(lower_bound), UB(upper_bound)
	{
		mp = &t_mp;
		nvars = initpnt.size();
	}
	~Objective(){}
};

