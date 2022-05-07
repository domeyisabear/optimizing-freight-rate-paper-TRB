#include "gurobi_c++.h"
#include "Parameters.h"
#include <iomanip>
#include "RelaxMasterProblem.h"
#include "BnC.h"
#include <iostream>
#include "RelaxMasterProblem_exf.h"

class InputParser {
public:
	InputParser(int &argc, char **argv) {
		for (int i = 1; i < argc; ++i)
			this->tokens.push_back(std::string(argv[i]));
	}
	/// @author iain
	const std::string& getCmdOption(const std::string &option) const {
		std::vector<std::string>::const_iterator itr;
		itr = std::find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
			return *itr;
		}
		static const std::string empty_string("");
		return empty_string;
	}
	/// @author iain
	bool cmdOptionExists(const std::string &option) const {
		return std::find(this->tokens.begin(), this->tokens.end(), option)
			!= this->tokens.end();
	}
private:
	std::vector <std::string> tokens;
};

int main(int argc, char *argv[])
{
	/// parameters: 
	/// -b: benders decomposition used (no argument)
	/// -e: out epsilon
	/// -f: data directory
	/// -h: thread num
	/// -i: iteration limit
	/// -l: local search (no argument)
	/// -p: local search rate
	/// -r: range reduction (no argument)
	/// -t: time limit
	
	InputParser input(argc, argv);

	if (input.cmdOptionExists("-b"))
		benders_used = true;
	else
		benders_used = false;

	if (input.cmdOptionExists("-l"))
		local_search_used = true;
	else
		local_search_used = false;

	if (input.cmdOptionExists("-r"))
		range_reduction_used = true;
	else
		range_reduction_used = false;

	const string &eps = input.getCmdOption("-e");
	if (!eps.empty())
		out_eps = stod(eps);

	string fn = input.getCmdOption("-f");
	if(fn.empty())
	{
		cout << "data path is needed!" << endl;
		//fn = R"(C:\Users\wang_\Google Drive\my papers\slot pricing\data\small_50_linear\)";
		return 1;
	}

	const string &hn = input.getCmdOption("-h");
	if (!hn.empty())
		threadPoolSize = stoi(hn);

	const string &it = input.getCmdOption("-i");
	if (!it.empty())
		iterationLim = stoi(it);

	const string &tl = input.getCmdOption("-t");
	if (!tl.empty())
		timeLimit = stoi(tl);

	const string &lr = input.getCmdOption("-p");
	if (!lr.empty())
		local_search_rate = stod(lr);

	cout << "---------------------------------------------------" << endl;
	cout << "Parameter summary:" << endl;
	cout << "data path: " << fn << endl;
	cout << "eps: " << out_eps << endl;
	cout << "thread pool size: " << threadPoolSize << endl;
	cout << "iteration limit: " << iterationLim << endl;
	cout << "time limit: " << timeLimit << endl;
	cout << "local search: " << local_search_used << endl;
	cout << "local search rate: " << local_search_rate << endl;
	cout << "range reduction: " << range_reduction_used << endl;
	cout << "benders: " << benders_used << endl;
	cout << "---------------------------------------------------" << endl;

	param = new Parameters(fn);
	clock_t start = clock();

	//vector<double> rateMin = vector<double>(param->demandNum);
	//vector<double> rateMax = vector<double>(param->demandNum);
	//for (int i = 0; i < param->demandNum; i++)
	//{
	//	rateMin[i] = param->demandList[i].rateMin;
	//	rateMax[i] = param->demandList[i].rateMax;
	//}
	//RelaxMasterProblem rmp = RelaxMasterProblem(rateMin, rateMax);
	//RelaxMasterSolution rms = rmp.getOptimalValue();
	//cout << rms.objVal << endl;

	BnC p = BnC();
	RelaxMasterSolution rms=p.getOptimalValue();
	cout << rms.objVal << endl;

	MasterProblem mp = MasterProblem();
	cout << mp.getOptimalValue(rms.r_od) << endl;

	string fname = "r_od";
	if (local_search_used)
		fname += "_localsearch";
	if (range_reduction_used)
		fname += "_rangereduction";
	if (benders_used)
		fname += "_benders";
	fname += ".csv";	
	saveResult(rms.r_od, param->getWorkDir() + fname);

	double time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
	cout << "time:" << time_elps << endl;
	delete param;
	return 0;
}

