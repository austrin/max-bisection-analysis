/* Implementation of the standard bivariate gaussian CDF.
 * Based on the Matlab implementation by Genz, converted to C++.
 */
#ifndef __BVNL_H
#define __BVNL_H

// Pr[X <= dh and Y <= dk] where X and Y are r-correlated gaussians.
// The two functions provide a lower and upper bound respectively.
double bvnl_down(double dh, double dk, double r);
double bvnl_up(double dh, double dk, double r);

#endif
