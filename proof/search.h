/* Function to prove a lower bound on the approximation ratio of a
 * given rounding.
 */
#ifndef __SEARCH_H
#define __SEARCH_H

#include "configurations.h"
#include "alpha.h"

/* Attempt to prove that the approximatio ratio inside the cube
 * initbox is at least target_ratio when using the given rounding.
 * Proceed to depth at most max_depth.
 *
 * Returns:
 * +1 if success
 *  0 if inconclusive
 * -1 if counter-example found
 */
int search(const Conf& initbox, const Pairing &rounding, 
		   double target_ratio, int max_depth);

#endif
