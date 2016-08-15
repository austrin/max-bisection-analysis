/* Functions related to the CDF of standard gaussian random
 * variables and vectors.
 *
 * All functions correspond to the functions in the paper with the
 * same names.
 */
#ifndef __GAUSSIAN_H
#define __GAUSSIAN_H

#include "intervals.h"

I Phi(const I& x);
I Phi_inv(const I& x);
I Gamma(const I& q1, const I& q2, const I& rho);
I Lambda(const I& r1, const I& r2, const I& trho);

#endif
