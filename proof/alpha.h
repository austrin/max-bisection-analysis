/* Functions for evaluating the approximation ratio alpha on a given
 * configuration, both using the basic linear rounding algorithm and
 * the pairing rounding algorithm.
 */
#ifndef __ALPHA_H
#define __ALPHA_H

#include "intervals.h"
#include "configurations.h"


// Representation of a pairing rounding algorithm
// The function f is f(x) = A*max(0, x-B)
class Pairing {
 public:
	const double c, A, B;

	Pairing(double c, double A, double B);
	
	// Evaluate f at x.  Requires x to be non-negative.
	I f(const I& x) const;
	
	// Compute the range of possible r for a given (range of) mu.
	// Corresponds to the weak bound in Lemma 5.1
	I r(const I& mu) const;
};

// Evaluate alpha at a configuration and range of r-values
I alpha(const Conf& X, const I& r1, const I& r2);

// Evaluate a configuration for the linear rounding algorithm
I alpha_linear(const Conf& X, double c);

// Evaluate a configuration using the pairing rounding algorithm
I alpha_pairing(const Conf &X, const Pairing &P);


#endif
