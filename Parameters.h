#pragma once
#include <gurobi_c++.h>
#include<vector>
#include<set>
#include<cmath>
#include "threadpool.h"

#define Round(x) ((int)(x<0?x-0.5:x+0.5))
#define MAX_NUM 1e10
#define GRBVar3d vector<vector<vector<GRBVar>>>
#define GRBVar2d vector<vector<GRBVar>>
#define GRBVar1d vector<GRBVar>

using namespace std;

extern int timeLimit;
extern int iterationLim;
extern int threadPoolSize;
extern int threadNumPerSolver;
extern int segmentNum;
extern GRBEnv env;
extern double out_eps;
extern double in_eps;
extern bool local_search_used;
extern double local_search_rate;
extern bool range_reduction_used;
extern bool benders_used;

struct SubSolution
{
	double objVal;
	vector<double> g_od;
	vector<double> q_od_h;
	vector<double> R_od;
	vector<double> T_od;
};

struct SubDualSolution
{
	double objVal;
	vector<double> lambda;
	vector<double> phi;
	vector<double> pi;
};

struct RelaxMasterSolution
{
	vector<double> r_od;
	vector<double> g_od;
	double objVal;
	double LB;
};

struct RelaxSubDualSolution
{
	double objVal;
	vector<double> eta_od_2;
	vector<double> eta_od_1;
	double eta_0;
};

struct DemandLineApproxParam
{
	double a;
	double b;
};

struct ProblemNode
{
	vector<double> rateMin;
	vector<double> rateMax;
	double objVal;
	double LB;
	RelaxMasterSolution rms;
	bool operator < (const ProblemNode &a) const
	{
		return LB < a.LB;
	}
};

class Parameters
{
private:
	struct Demand
	{
		int originPort;
		int destinationPort;
		double rateMin;
		double rateMax;
		vector<double> demandParam;
		vector<double> ratioScenario;
		vector<int> containerPathIndex;
	};

	struct ContainerPath
	{
		int demandIndex;
		double handleCost;
		vector<int> shippingLegIndex;
	};

	struct ShippingLeg
	{
		int route;
		int tail;
		int head;
		vector<int> pathIndexList;
		vector<double> capacityScenario;
	};
	string workDir;

public:
	vector<Demand> demandList;
	vector<ContainerPath> containerPathList;
	vector<ShippingLeg> shippingLegList;
	vector<double> scenarioProbList;

	int demandNum;
	int containerPathNum;
	int shippingLegNum;
	int scenarioNum;

	threadpool* executor;

	string getWorkDir();
	
	Parameters(string fdir);
	~Parameters();
};

extern Parameters* param;