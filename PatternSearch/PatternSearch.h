/* PatternSearch.h
 * declarations of the PatternSearch base class member 
 * functions and data for optimization
 * Liz Dolan College of William and Mary 1999
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



#ifndef _PatternSearch_
#define _PatternSearch_

#include <istream>
#include <math.h>
#include <stdlib.h>
#include "Objective.h"  //has declaration of function to be minimized, starting point, n
#include "vec.h"        //vector and matrix classes
#include "cmat.h"
#include "subscrpt.h"

//default definitions; should be defined by user in objective.h, but just in case...
extern int maxCalls;    //means program will only terminate based on grid refinement
extern double initialStepLength;

#ifndef stoppingStepLength
#define stoppingStepLength  sqrt(fabs(3.0 * (4.0/3.0 -1.0) -1.0))
#endif

using namespace std;
typedef char file[32];

class PatternSearch
{
 public:
   
 //constructors & destructor

  PatternSearch(Objective obj);
  //default constructor

  PatternSearch(const PatternSearch & Original);
  //deep copy constructor

  virtual ~PatternSearch();
  //destructor

 //other search initialization routines

  virtual void CopySearch(const PatternSearch & Original);
  //deep copier of PatternSearch class objects 

  void CleanSlate(int dimensions);
  //reinitialize all values

 //algorithmic routines

  virtual void ExploratoryMoves() = 0;
  //must be implemented as the pattern search algorithm
 
  virtual void NextPoint(int index, const Vector<double> & currentPoint, Vector<double> & nextPoint);
  //calculates the next trial point by adding the pattern col. at index to the current vector
  //returns the prospect in nextPoint

  virtual void ReplaceMinimum(const Vector<double> & newPoint, double  newValue); 
  //replaces the minimizer & the minimum objective function value
 
  virtual void ScalePattern(double scalar);
  //scale lattice step length by scalar

  virtual bool Stop();  
  //gives default stopping criteria based on maxCalls and stoppingStepLength

  virtual void fcnCall(int n, double *x, double & f, int & flag);
  //indirection of function call for purposes of keeping an accurate 
  //tally of the number of function calls

//pattern-altering functions

  void InitializePattern(int patternSize, const Matrix<double> * pat);
  //deletes any existing pattern and replaces it with the one pointed to by pat
 
  void ReadPatternFile(std::istream & fp);
  //may also pass cin as the input stream as desired
  //input first the pattern length and then the values of each trial vector
  //(i.e. input the pattern by column)

 //query functions

  int GetFunctionCalls() const { return functionCalls; }; 
  //number of objective function evaluations

  void GetMinPoint(Vector<double> & minimum) const; 
  //minimum should be of the correct dimension

  void GetMinVal(double & value) const; 
  //best objective function value found thus far

  void GetPattern(Matrix<double>* &pat) const; 
  //user should pass just a pointer, without preallocated memory
  //deep copy of the pattern to a Matrix pointer 
  //points to a newly allocated chunk of memory upon return

  void GetPatternLength(int & pattern) const; 
  //returns the number of "columns" of the pattern matrix

  virtual double GetStepLength() const {return latticeStepLength;};
  //measure of grid refinement, a.k.a. trial step length, delta

  int GetVarNo() const;  
  //returns the number of dimensions

  double getLowBound(int dimension);

  double getUpBound(int dimension);

  Objective objective;
 private:

  Matrix<double> * pattern;    //a pointer to a pattern matrix, with each trial direction a column
  int variables;               //dimension of the problem
  Vector<double> * minPoint;   //the minimizer
  int patternLength ;          //number of "columns" in pattern matrix, i.e. trial points
  double minValue;             //best objective function value calculated thus far
  double latticeStepLength;    //density of underlying lattice
  long functionCalls;          //tally of the number of function calls
  Vector<double> * lowBound;
  Vector<double> * upBound;
};

#endif



