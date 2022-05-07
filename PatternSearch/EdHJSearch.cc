/*  EdHJSearch.cc
 *  implementation of the EdHJSearch class to find a minimal objective function
 *  solution.
 *  Includes a minor modification to the basic Hooke and Jeeves strategy to
 *  avoid making pattern steps directly after contractions (which mostly cover
 *  the same ground that was already covered in the search step preceding
 *  the contraction.
 *  For a good description of the Hooke and Jeeves search algorithm
 *  I recommend Non-Linear Optimization Techniques by Box, Davies,
 *  and Swann, 1969.
 *  Liz Dolan, The College of William and Mary, 1999
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

#include "EdHJSearch.h"
using namespace std;

EdHJSearch::EdHJSearch(Objective& obj): PatternSearch(obj)
{
	step = initialStepLength;
	factor = 0.5;
}

EdHJSearch::EdHJSearch(const EdHJSearch & Original) : PatternSearch(Original)
{
	step = Original.step;
	factor = Original.factor;
}

EdHJSearch::~EdHJSearch()
{
}

void EdHJSearch::CopySearch(const EdHJSearch & Original)
{
	PatternSearch::CopySearch(Original);
	step = Original.step;
	factor = Original.factor;
}

void EdHJSearch::CleanSlate(int dimensions)
//reinitialize all values
{
	PatternSearch::CleanSlate(dimensions);
	step = initialStepLength;
}

void EdHJSearch::ExploratoryMoves()
{
	int dimens = GetVarNo();
	Vector<double> currentPoint(dimens);
	Vector<double> lastImprovingPoint(dimens); //last base point
	Vector<double> storage(dimens); //for intermediate storage to reduce rounding error
	Vector<double> direction(dimens, 0.0); //direction of pattern extending step
	double value;                   //objective function value
	double positiveValue;           //obj.fun. value in the positive step 
	double negativeValue;           // " for the negative step 
	double lastImprovingValue;      //obj.fun.value of last base point
	int success = 0;
	GetMinPoint(currentPoint);      //initialize to the user initial point
	GetMinPoint(lastImprovingPoint);
	GetMinPoint(storage);
	GetMinVal(value);
	GetMinVal(lastImprovingValue);
	bool foundImprove = false;
	bool contracted = false;
	do
	{
		for (int iteration = 0; iteration < dimens; iteration++)
		{
			double up = getUpBound(iteration);
			if (currentPoint[iteration] + step > up)
				currentPoint[iteration] = up;
			else
				currentPoint[iteration] += step; /////////////////// update in one direction
			fcnCall(GetVarNo(), currentPoint.begin(), positiveValue, success);

			if (success != 1)
			{
				positiveValue = value + 1.0;
				//if the call returned unsuccessfully, set positiveValue
					//to a value that will not be improving	    
			}
			if (positiveValue < value)
			{
				value = positiveValue;
				foundImprove = true;
				//continue search in other dimensions from here
			}//if positive is better

			if (!foundImprove)
			{
				currentPoint = storage;
				double low = getLowBound(iteration);
				if (currentPoint[iteration] - step < low)
					currentPoint[iteration] = low;
				else
					currentPoint[iteration] -= (step);   //////////////// update in one direction
				fcnCall(GetVarNo(), currentPoint.begin(), negativeValue, success);
				if (success != 1)
				{
					negativeValue = value + 1.0;
					//same kludge as in positive case
				}

				if (negativeValue < value)
				{
					value = negativeValue;
					foundImprove = true;
					//continue search in other dimensions from here
				}//if negative direction is better
			}//if we need to check the negative
			if (!foundImprove)
			{//reset to original position
				currentPoint = storage;
			}//if neither direction gave improvement
			else
			{
				storage = currentPoint;
			}
			foundImprove = false; //reset for next iteration
		}//for
		direction = currentPoint - lastImprovingPoint;
		//direction now holds the extended pattern step vector
		if (value < lastImprovingValue)
		{
			//check whether the "new" point is within factor*step of the old
			if (isnear(lastImprovingPoint, currentPoint, factor*step))
			{
				currentPoint = lastImprovingPoint;
				value = lastImprovingValue;
				storage = currentPoint;
			}
			else
			{
				lastImprovingValue = value;
				ReplaceMinimum(currentPoint, value);
				lastImprovingPoint = currentPoint;
				if (!contracted)   //my personal modification to the algorithm
				{
					//take the pattern extending step and find its value
					currentPoint = direction + currentPoint;
					for(int i=0;i<dimens;i++)
					{
						double up = getUpBound(i);
						double low = getLowBound(i);
						if (currentPoint[i] < low)
							currentPoint[i] = low;
						if (currentPoint[i] > up)
							currentPoint[i] = up;
					}
					storage = currentPoint;
					fcnCall(dimens, currentPoint.begin(), value, success);
				}
			}
			contracted = false;
		}
		else
		{
			if (isnear(currentPoint, lastImprovingPoint, factor * step))
			{
				step = step * factor;
				contracted = true;
			}
			else
			{
				//this case can only occur after an unsuccessful
				//search about a pattern-step-located point
				//move back to the point that was improving from the
				//search about the last base point
				GetMinPoint(currentPoint);
				GetMinVal(value);
				storage = currentPoint;
				contracted = false;
			}
		}
	} while (!Stop(step));//while we haven't stopped()
}//ExploratoryMoves


bool EdHJSearch::Stop(double pace)
//makes certain search has not exceeded maxCalls or 
//is stepping at less than stoppingStepLength
{
	if (maxCalls > -1)
		if (GetFunctionCalls() > maxCalls)
			return true;
	if (pace < stoppingStepLength)
		return true;
	return false;
}//Stop

double EdHJSearch::GetStepLength()
{
	return step;
}





