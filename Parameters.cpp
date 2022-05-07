#include "Parameters.h"
#include "ReadCsv.h"
#include "threadpool.h"
#include <thread>

int timeLimit = 86400;
int threadPoolSize = 12;
int threadNumPerSolver = 1;
int iterationLim = MAX_NUM;
int segmentNum = 1;
double out_eps = 0.0001;
double in_eps = 0.0001;
bool local_search_used = true;
double local_search_rate = 0.3;
bool range_reduction_used = true;
bool benders_used = true;

Parameters *param = NULL;

GRBEnv env = GRBEnv();

string Parameters::getWorkDir()
{
	return workDir;
}

Parameters::Parameters(string data_dir)
{
	workDir = data_dir;

	// read data
	vector<vector<double>> demandData = ReadCsv::read(workDir + "od_demand.csv");
	vector<vector<double>> pathData = ReadCsv::read(workDir + "path.csv");
	vector<vector<double>> legData = ReadCsv::read(workDir + "leg.csv");
	vector<vector<double>> ratioData = ReadCsv::read(workDir + "demand_ratio.csv");
	vector<vector<double>> capData = ReadCsv::read(workDir + "capacity.csv");
	vector<vector<double>> probData = ReadCsv::read(workDir + "probability.csv");

	if (demandData.size() == 0 || pathData.size() == 0 || legData.size() == 0 || ratioData.size() == 0 || capData.size() == 0 || probData.size() == 0)
		throw("no data record or incomplete data record!");

	// initialize
	demandNum = demandData.size();
	containerPathNum = pathData.size();
	shippingLegNum = legData.size();
	scenarioNum = probData[0].size();

	// demand list
	for(int i=0;i<demandNum;i++)
	{
		Demand d;
		d.originPort = Round(demandData[i][1])-1;
		d.destinationPort = Round(demandData[i][2])-1;
		d.rateMin = demandData[i][3];
		d.rateMax = demandData[i][4];
		vector<double> param;
		for(int j=5;j<demandData[i].size();j++)
		{
			param.push_back(demandData[i][j]);
		}
		d.demandParam = param;
		vector<double> ratio;
		for(int j=1;j<ratioData[i].size();j++)
		{
			ratio.push_back(ratioData[i][j]);
		}
		d.ratioScenario = ratio;
		d.containerPathIndex = vector<int>();
		demandList.push_back(d);
	}

	// container path list
	for(int i=0;i<containerPathNum;i++)
	{
		ContainerPath p;
		p.demandIndex = Round(pathData[i][1]) - 1;
		p.handleCost = pathData[i][2];
		containerPathList.push_back(p);
		demandList[p.demandIndex].containerPathIndex.push_back(i);
	}
	
	// shipping leg
	for(int i=0;i<shippingLegNum;i++)
	{
		ShippingLeg l;
		l.route = Round(legData[i][1]) - 1;
		l.tail = Round(legData[i][2]) - 1;
		l.head = Round(legData[i][3]) - 1;

		vector<int> path;
		for(int j=4;j<legData[i].size();j++)
		{
			int pInd = Round(legData[i][j]) - 1;
			path.push_back(pInd);
			containerPathList[pInd].shippingLegIndex.push_back(i);
		}
		l.pathIndexList = path;

		vector<double> cap;
		for(int j=1;j<capData[i].size();j++)
		{
			cap.push_back(capData[i][j]);
		}
		l.capacityScenario = cap;
		shippingLegList.push_back(l);
	}

	// probability
	scenarioProbList = probData[0];

	executor = new threadpool(threadPoolSize);
}

Parameters::~Parameters()
{
	delete executor;
}
