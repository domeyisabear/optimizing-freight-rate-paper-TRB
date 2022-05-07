#pragma once
#include "Parameters.h"
#include "RelaxSubProblem.h"
#include "MasterProblem.h"

class BnC
{
private:
	MasterProblem mp = MasterProblem();
public:
	BnC();
	RelaxMasterSolution getOptimalValue();
	~BnC()
	{};
};

