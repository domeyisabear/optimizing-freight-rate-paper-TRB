/* PatternSearch.cc
 * class implementation of pattern search optimization basis
 * Liz Dolan, The College of William and Mary, 1999
 *
 *
 *  The author of this software is Elizabeth  D. Dolan.
 *  Permission to use, copy, modify, and distribute this software
 *  for any purpose without fee is hereby granted, provided that
 *  this entire notice is included in all copies of any software
 *  which is or includes a copy or modification of this software
 *  and in all copies of the supporting documentation for such
 *  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT
 *  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, THE AUTHOR
 *  OFFERS NO REPRESENTATION OR WARRANTY OF ANY KIND
 *  CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS
 *  FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#include "PatternSearch.h"
#include <iostream>
using namespace std;
//constructors & destructor

int maxCalls = 2000;
double initialStepLength = 100;

PatternSearch::PatternSearch(Objective obj)
{
	objective = obj;
	variables = obj.nvars;
	patternLength = 0;
	functionCalls = 0;
	pattern = NULL;
	double * startPoint = NULL;
	objective.initpt(startPoint);
	minPoint = new Vector<double>(variables, startPoint);
	double * lowPnt = NULL;
	objective.low_bound(lowPnt);
	lowBound = new Vector<double>(variables, lowPnt);
	double * upPnt = NULL;
	objective.up_bound(upPnt);
	upBound = new Vector<double>(variables, upPnt);
	int flag = 1;
	latticeStepLength = initialStepLength;
	objective.fcn(startPoint, minValue, flag);
	
	if (flag != 1)
		std::cerr << "\nError signal in objective function at starting point.\n";
}//constructor initializes private data members to user defined starting values

PatternSearch::PatternSearch(const PatternSearch & Original)
{
	objective = Original.objective;
	variables = Original.GetVarNo();
	Original.GetPattern(pattern);
	Original.GetPatternLength(patternLength);
	functionCalls = Original.functionCalls;
	minPoint = new Vector<double>(*(Original.minPoint));
	latticeStepLength = Original.latticeStepLength;
	minValue = Original.minValue;
	lowBound = new Vector<double>(*(Original.lowBound));
	upBound = new Vector<double>(*(Original.upBound));
}//deep copy constructor

PatternSearch::~PatternSearch()
{
	delete pattern;
	pattern = NULL;
	delete minPoint;
	minPoint = NULL;
	delete lowBound;
	lowBound = NULL;
	delete upBound;
	upBound = NULL;
	//note: matrix and vector classes have their own destructors
}//destructor

//other search initialization routines

void PatternSearch::CopySearch(const PatternSearch & Original)
{
	int dimension = Original.GetVarNo();
	CleanSlate(dimension);
	Original.GetPattern(pattern);
	Original.GetPatternLength(patternLength);
	functionCalls = Original.functionCalls;
	minPoint = new Vector<double>(*(Original.minPoint));
	latticeStepLength = Original.latticeStepLength;
	minValue = Original.minValue;
	lowBound = new Vector<double>(*(Original.lowBound));
	upBound = new Vector<double>(*(Original.upBound));
}//deep copy  

void PatternSearch::CleanSlate(int dimensions)
//reinitialize all values
{
	variables = dimensions;
	functionCalls = 0;
	delete pattern;
	patternLength = 0;
	delete minPoint;
	pattern = NULL;
	minPoint = NULL;
	double * startPoint = NULL;
	objective.initpt(startPoint);
	minPoint = new Vector<double>(variables, startPoint);
	double * lowPnt = NULL;
	objective.low_bound(lowPnt);
	lowBound = new Vector<double>(variables, lowPnt);
	double * upPnt = NULL;
	objective.up_bound(upPnt);
	upBound = new Vector<double>(variables, upPnt);
	int flag = 1;
	latticeStepLength = initialStepLength;
	objective.fcn(startPoint, minValue, flag);
	if (flag != 1)
		std::cerr << "\nError signal in objective function at starting point.\n";
}//reinitialize all search values to appropriate user-defined start values

//algorithmic routines

void PatternSearch::NextPoint(int index, const Vector<double> & currentPoint, Vector<double> & nextPoint)
{
	//To get the next point I add the currentPoint to the
	//product of the pattern vector at index and lattice step length
	if (pattern != NULL && patternLength > index)
		nextPoint = currentPoint + (latticeStepLength * (*pattern).col(index));
}

void PatternSearch::ReplaceMinimum(const Vector<double> & newPoint, double newValue)
{
	(*minPoint) = newPoint;
	minValue = newValue;
}//UpdateMinPoint

void PatternSearch::ScalePattern(double scalar)
{
	latticeStepLength = latticeStepLength * scalar;
}//ScalePattern scale delta step length by scalar

bool PatternSearch::Stop()
{
	if (maxCalls > -1)
		if (functionCalls >= maxCalls)
			return true;
	if (latticeStepLength < stoppingStepLength)
		return true;
	else
		return false;
}//Stop based on cap on function evaluations or lower limit on trial step length 

void PatternSearch::fcnCall(int n, double *x, double & f, int & flag)
{
	objective.fcn(x, f, flag);
	functionCalls++;
}//record calls to objective function

//pattern-altering functions

void PatternSearch::InitializePattern(int patternSize, const Matrix<double> *pat)
{
	if (pattern != NULL && pattern != pat) delete pattern;
	pattern = new Matrix<double>((*pat));
	patternLength = patternSize;
}//InitializePattern

void PatternSearch::ReadPatternFile(std::istream & fp)
//Reads the pattern in column by column.
//The first digit read must be the number of columns.
{
	if (patternLength != 0)
	{
		delete pattern;
		pattern = NULL;
	}//if
	if (!fp) return;
	fp >> patternLength;  //the length of the pattern must precede the pattern
	pattern = new Matrix<double>(variables, patternLength);
	for (int i = 0; i < patternLength; i++)
	{
		for (int j = 0; j < variables; j++)
		{
			fp >> (*pattern)[j][i];
		}//inner for
	}//outer for
}//ReadPatternFile

//query functions

void PatternSearch::GetMinPoint(Vector<double> & minimum) const
{
	minimum = (*minPoint);
}//GetMinPoint

void PatternSearch::GetMinVal(double & value) const
{
	value = minValue;
}//GetMinVal

void PatternSearch::GetPattern(Matrix<double> * &pat) const
//user should pass just a pointer, without preallocated memory
{
	pat = new Matrix<double>((*pattern));
}//GetPattern

void PatternSearch::GetPatternLength(int & pattern) const
{
	pattern = patternLength;
}//GetPatternLength - trial steps in pattern (# of columns)

int PatternSearch::GetVarNo() const
{
	return variables;
}//GetVarNo = dimension of the problem

double PatternSearch::getLowBound(int dimension)
{
	return (*lowBound)[dimension];
}

double PatternSearch::getUpBound(int dimension)
{
	return (*upBound)[dimension];
}
















