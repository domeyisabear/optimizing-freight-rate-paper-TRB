#include "BnC.h"
#include "RelaxMasterProblem.h"
#include <algorithm>
#include "RangeReduction.h"
#include "PatternSearch/Objective.h"
#include "PatternSearch/EdHJSearch.h"
#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "RelaxMasterProblem_exf.h"

BnC::BnC()
{}

RelaxMasterSolution BnC::getOptimalValue()
{
	clock_t start = clock();

	// initialize
	double UB = MAX_NUM;
	double LB = -MAX_NUM;

	string fname = "UB_LB_BnB_quad";
	if (local_search_used)
		fname += "_localsearch";
	if (range_reduction_used)
		fname += "_rangereduction";
	if (benders_used)
		fname += "_benders";
	fname += ".csv";

	ofstream outputfile(param->getWorkDir() + fname);
	outputfile << setprecision(16) << "UB,LB,gap,time\n";

	ProblemNode root;
	root.rateMin = vector<double>(param->demandNum);
	root.rateMax = vector<double>(param->demandNum);
	for (int i = 0; i < param->demandNum; i++)
	{
		root.rateMin[i] = param->demandList[i].rateMin;
		root.rateMax[i] = param->demandList[i].rateMax;
	}

	RelaxMasterSolution rms;
	if(benders_used)
	{
		RelaxMasterProblem* rmp = new RelaxMasterProblem(root.rateMin, root.rateMax);
		rms = rmp->getOptimalValue();
		delete rmp;
	}
	else
	{
		RelaxMasterProblem_exf* rmp = new RelaxMasterProblem_exf(root.rateMin, root.rateMax);
		rms = rmp->getOptimalValue();
		delete rmp;
	}
	
	root.LB = rms.LB;
	root.objVal = rms.objVal;
	root.rms = rms;

	RelaxMasterSolution incum_rms = rms;
	LB = rms.LB;
	UB = mp.getOptimalValue(rms.r_od);
	double rela_gap = (UB - LB) / fabs(LB + 0.000000001);
	double time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;

	double range_reduct_time = 0;
	double local_search_time = 0;
	double qcp_time = 0;

	cout << "UB: " << UB << ", LB:" << LB << ", gap: " << rela_gap * 100 << "%"<< endl;
	outputfile << UB << "," << LB << "," << rela_gap << "," << time_elps << "\n";
	set<ProblemNode> probset;
	probset.insert(root);

	bool improve1 = false;
	bool improve2 = false;
	bool improve3 = false;
	bool improve4 = false;

	int iteration = 1;

	try
	{
		while (rela_gap > out_eps && time_elps <= timeLimit && iteration <= iterationLim && (!probset.empty()))
		{
			// set the problem node with the lowest objVal as the current node
			set<ProblemNode>::iterator nodeInd = min_element(begin(probset), end(probset));
			ProblemNode currentNode = *nodeInd;

			// update LB
			LB = currentNode.LB > LB ? currentNode.LB : LB;

			// determine the primal solution of current node
			vector<SubSolution> subsolution_list;
			vector<future<SubSolution>> futures;
			for (int j = 0; j < param->scenarioNum; j++)
			{
				futures.push_back(param->executor->commit(
					[](int sind, ProblemNode& node)
				{
					RelaxSubProblem sp(sind, node.rateMin, node.rateMax);
					return sp.getPrimalSolution(node.rms);
				},
					j,
					currentNode
					));
			}
			for (int j = 0; j < param->scenarioNum; j++)
			{
				subsolution_list.push_back(futures[j].get());
			}

			// determine the od to branch
			vector<double> demand_list = getDemand(currentNode.rms.r_od);
			vector<double> obj_diff(param->demandNum);
			double objVal = 0;
			for (int i = 0; i < param->demandNum; i++)
			{
				double revise_obj = 0;
				double origin_obj = 0;
				for (int w = 0; w < param->scenarioNum; w++)
				{
					double actualDemand = demand_list[i] * param->demandList[i].ratioScenario[w];
					double demand_exceed = subsolution_list[w].g_od[i] - actualDemand;

					double remain_exceed = demand_exceed >= 0 ? demand_exceed : 0;
					for (int p = param->demandList[i].containerPathIndex.size() - 1; p >= 0; p--)
					{
						int pathInd = param->demandList[i].containerPathIndex[p];
						double q_od_h = subsolution_list[w].q_od_h[pathInd];
						double reduction = (remain_exceed > q_od_h ? q_od_h : remain_exceed);

						double revise_q = q_od_h - reduction;
						revise_obj += param->scenarioProbList[w] * (param->containerPathList[pathInd].handleCost - currentNode.rms.r_od[i])*revise_q;
						origin_obj += param->scenarioProbList[w] * param->containerPathList[pathInd].handleCost * q_od_h;

						remain_exceed -= reduction;
					}

					double R_od = subsolution_list[w].R_od[i];
					origin_obj -= param->scenarioProbList[w] * R_od;
				}
				obj_diff[i] = revise_obj - origin_obj;
			}

			int splitInd = -1;
			double max_diff = 0;
			double sum_diff = 0;
			for (int i = 0; i < param->demandNum; i++)
			{
				if (obj_diff[i] > max_diff)
				{
					splitInd = i;
					max_diff = obj_diff[i];
				}
				sum_diff += obj_diff[i];
			}

			// update UB
			if (UB - 1e-3 > currentNode.objVal + sum_diff)
			{
				UB = currentNode.objVal + sum_diff;
				incum_rms = currentNode.rms;
				improve1 = true;
			}

			//// range reduction
			clock_t rr_start = clock();
			if (range_reduction_used)
			{
				vector<future<void>> futures;
				for (int i = 0; i < param->demandNum; i++)
				{
					futures.push_back(param->executor->commit(
						[](int ind, ProblemNode* node, double UB)
					{
						RangeReduction rr(node->rateMin, node->rateMax, UB);
						node->rateMax[ind] = rr.getNewBound(ind, false);
						node->rateMin[ind] = rr.getNewBound(ind, true);					
					},
						i, &currentNode, UB
						));
				}

				for (int i = 0; i < param->demandNum; i++)
				{
					futures[i].get();
				}
			}
			// check feasible
			bool is_feas = true;
			for(int i=0;i<param->demandNum;i++)
			{
				if(currentNode.rateMin[i]>currentNode.rateMax[i]-1e-3)
				{
					is_feas = false;
					break;
				}
			}
			if(!is_feas)
			{
				probset.erase(nodeInd);
				time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
				iteration++;
				rela_gap = (UB - LB) / fabs(LB + 0.000000001);
				cout << "UB: " << UB << ", LB:" << LB << ", gap: " << rela_gap * 100 << "%, iteration: " << iteration << ", time elapsed: " << time_elps << endl;
				outputfile << UB << "," << LB << "," << rela_gap << "," << time_elps << "\n";
				continue;
			}
			range_reduct_time += (double)(clock() - rr_start) / CLOCKS_PER_SEC;

			// split current node
			ProblemNode pn1 = currentNode;
			ProblemNode pn2 = currentNode;

			pn1.rateMax[splitInd] = currentNode.rms.r_od[splitInd];
			pn2.rateMin[splitInd] = currentNode.rms.r_od[splitInd];

			// evaluate two nodes
			clock_t qcp_start = clock();
			// pn1
			if(benders_used)
			{
				RelaxMasterProblem* rmp1 = new RelaxMasterProblem(pn1.rateMin, pn1.rateMax);
				pn1.rms = rmp1->getOptimalValue();
				delete rmp1;
			}
			else
			{
				RelaxMasterProblem_exf* rmp1 = new RelaxMasterProblem_exf(pn1.rateMin, pn1.rateMax);
				pn1.rms = rmp1->getOptimalValue();
				delete rmp1;
			}
			pn1.objVal = pn1.rms.objVal;
			pn1.LB = pn1.rms.LB;
			// pn2
			if(benders_used)
			{
				RelaxMasterProblem* rmp2 = new RelaxMasterProblem(pn2.rateMin, pn2.rateMax);
				pn2.rms = rmp2->getOptimalValue();
				delete rmp2;
			}
			else
			{
				RelaxMasterProblem_exf* rmp2 = new RelaxMasterProblem_exf(pn2.rateMin, pn2.rateMax);
				pn2.rms = rmp2->getOptimalValue();
				delete rmp2;
			}
			pn2.objVal = pn2.rms.objVal;
			pn2.LB = pn2.rms.LB;
			qcp_time += (double)(clock() - qcp_start) / CLOCKS_PER_SEC;

			// update UB
			double s1 = mp.getOptimalValue(pn1.rms.r_od);
			if (s1 < UB-1e-3)
			{
				UB = s1;
				incum_rms = pn1.rms;
				improve2 = true;
			}

			double s2 = mp.getOptimalValue(pn2.rms.r_od);
			if (s2 < UB - 1e-3)
			{
				UB = s2;
				incum_rms = pn2.rms;
				improve3 = true;
			}

			// local search
			clock_t ls_start = clock();
			if (local_search_used)
			{
				//srand((unsigned)time(NULL));
				double rv = rand() / double(RAND_MAX);
				//cout << rv << endl;
				if (rv < local_search_rate)
				{
					cout << "local search ...  ";
					maxCalls = 2000;
					initialStepLength = 50;
					Objective obj;
					if (improve1 || improve2 || improve3 || improve4)
						obj = Objective(incum_rms.r_od, mp, root.rateMin, root.rateMax);
					else
						obj = Objective((s1 < s2 ? pn1.rms.r_od : pn2.rms.r_od), mp, root.rateMin, root.rateMax);
					Vector<double> Edminimum(param->demandNum);
					EdHJSearch EdHJsearcher(obj);
					double EdMinimum;

					EdHJsearcher.ExploratoryMoves();
					EdHJsearcher.GetMinPoint(Edminimum);
					EdHJsearcher.GetMinVal(EdMinimum);

					improve1 = false;
					improve2 = false;
					improve3 = false;
					improve4 = false;

					if (EdMinimum < UB - 1e-3)
					{
						UB = EdMinimum;
						for (int i = 0; i < param->demandNum; i++)
							incum_rms.r_od[i] = Edminimum[i];
						improve4 = true;
					}
					cout << " complete!" << endl;
				}
			}
			local_search_time += (double)(clock() - ls_start) / CLOCKS_PER_SEC;

			time_elps = (double)(clock() - start) / CLOCKS_PER_SEC;
			iteration++;
			rela_gap = (UB - LB) / fabs(LB + 0.000000001);
			cout << "UB: " << UB << ", LB:" << LB << ", gap: " << rela_gap * 100 << "%, iteration: " << iteration << ", time elapsed: " << time_elps << endl;
			outputfile << UB << "," << LB << "," << rela_gap << "," << time_elps << "\n";

			// remove current node
			probset.erase(nodeInd);

			// add subproblem to set
			if (pn1.LB < UB)
				probset.insert(pn1);
			if (pn2.LB < UB)
				probset.insert(pn2);
			outputfile.flush();

		}
		cout << "total time : " << time_elps << endl;
		cout << "range reduction time : " << range_reduct_time << endl;
		cout << "local search time : " << local_search_time << endl;
		cout << "qcp time : " << qcp_time << endl;
	}
	catch(...)
	{
		cout << "Exception occured. Algorithm stops!" << endl;
	}

	incum_rms.objVal = UB;
	incum_rms.LB = LB;

	outputfile.close();
	return incum_rms;

}
