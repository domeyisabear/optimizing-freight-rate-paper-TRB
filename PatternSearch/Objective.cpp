#include "Objective.h"

/* I defined my stoppingStepLength to be tolerance, which I must define here.
 * Note that the value chosen is machine/compiler dependent as approximately
 *   the square root of machine epsilon
 */
double tolerance = sqrt(fabs(3.0 * (4.0 / 3.0 - 1.0) - 1.0));
/* Initial trial step length, reflecting the refinement of the search lattice
 */
double initialStep = .25;