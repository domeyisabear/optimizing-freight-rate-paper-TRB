/*  EdHJSearch.h 
 *  header file for an edited Hooke and Jeeves search based 
 *  on the class PatternSearch.
 *  Includes a minor modification to the basic Hooke and Jeeves strategy to
 *  avoid making pattern steps directly after pattern contractions (which mostly cover
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
 *  ANY EXPRESS OR IMPLIED WARRANTY.  
 */                                                                 

#ifndef _EdHJSearch_
#define _EdHJSearch_

#include "PatternSearch.h"

using namespace std;

class EdHJSearch: public PatternSearch
{

  public :
  EdHJSearch(Objective& obj); 
  //the dimensions of the search space are required for the constructor
  EdHJSearch(const EdHJSearch & Original); 
  //deep copy constructor
  ~EdHJSearch();
  //destructor

  void CopySearch(const EdHJSearch & Original);
  //deep search copy
  void ExploratoryMoves();
  //searches in Hooke and Jeeves pattern for a better solution to the objective function
  //cancels pattern extending steps after mesh refinements
  bool Stop(double pace);
  //makes certain search has not already exceeded maxCalls or 
  //is stepping at less than stoppingStepLength
  double GetStepLength();
  //overrides the PatternSearch version with an acurate length
  void CleanSlate(int dimensions);
   //reinitialize all values
 private:
   double step;
   double factor; //factor by which the step is reduced
};//class EdHJSearch


#endif










